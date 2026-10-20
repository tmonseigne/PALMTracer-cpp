"""Vérifie les deux formats de rapport et l'échappement des diagnostics."""

import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from junit_to_markdown import parse_report, render_report, source_files, test_title, title_case


class JunitReportTests(unittest.TestCase):
	"""Contrôles du convertisseur sans exécution des algorithmes C++."""

	# --------------------------------------------------
	def test_google_test_statuses(self):
		"""Distingue succès, échec, erreur, cas ignoré et cas non exécuté."""
		root = ET.fromstring('''<testsuites time="1.5" timestamp="2026-10-14T12:00:00">
            <testsuite name="Suite">
                <testcase classname="Suite" name="Value/OneDegree" time="0.1"/>
                <testcase classname="Suite" name="Failure"><failure message="Different"/></testcase>
                <testcase classname="Suite" name="Error"><error>Exception</error></testcase>
                <testcase classname="Suite" name="Skipped"><skipped/></testcase>
                <testcase classname="Suite" name="Disabled" status="notrun"/>
            </testsuite></testsuites>''')
		report = parse_report(root, Path("local/gtest.xml"))
		self.assertEqual([case.status for case in report.cases], ["Reussi", "Echec", "Erreur", "Ignore", "Non execute"])
		self.assertEqual(report.cases[0].name, "Value/OneDegree")
		self.assertEqual(report.cases[0].suite, "Suite")
		self.assertEqual(report.duration, 1.5)

	# --------------------------------------------------
	def test_ctest_names_and_output(self):
		"""Sépare le nom complet CTest et conserve la sortie d'un cas en échec."""
		root = ET.fromstring('''<testsuite time="0.2">
            <testcase name="Reference/PSFTests.Gradient/OneDegree" classname="Reference/PSFTests.Gradient/OneDegree" time="0.2">
                <failure>Difference</failure><system-out>Valeur observee</system-out>
            </testcase></testsuite>''')
		report = parse_report(root, Path("linux/ctest-report.xml"))
		self.assertEqual(report.cases[0].suite, "Reference/PSFTests")
		self.assertEqual(report.cases[0].name, "Gradient/OneDegree")
		self.assertIn("Valeur observee", report.cases[0].details)

	# --------------------------------------------------
	def test_html_diagnostics_are_escaped(self):
		"""Le contenu d'une sortie de test ne devient jamais du HTML exécutable."""
		root = ET.fromstring('''<testsuite><testcase name="Value&lt;double&gt;">
            <failure>&lt;script&gt;unexpected&lt;/script&gt;</failure></testcase></testsuite>''')
		content = render_report([parse_report(root, Path("local/test.xml"))])
		self.assertIn("Value&lt;double&gt;", content)
		self.assertIn("&lt;script&gt;", content)
		self.assertNotIn("<script>", content)
		self.assertIn("<details open>", content)

	# --------------------------------------------------
	def test_invalid_format(self):
		"""Refuse un XML qui ne contient pas un rapport de tests."""
		with self.assertRaises(ValueError):
			parse_report(ET.fromstring("<unknown/>"), Path("test.xml"))

	# --------------------------------------------------
	def test_summary_and_parameter_groups(self):
		"""Affiche le bilan Python et replie les variantes réussies dans leur fichier."""
		root = ET.fromstring('''<testsuites time="80.2" timestamp="2026-09-14T16:58:56">
            <testsuite name="Modes/PSFTests">
                <testcase classname="Modes/PSFTests" name="Gradient/OneDegree"/>
                <testcase classname="Modes/PSFTests" name="Gradient/Zero"/>
            </testsuite></testsuites>''')
		report = parse_report(root, Path("windows/ctest.xml"), {"PSFTests": "test_sp_PSF.cpp"})
		content = render_report([report])
		self.assertIn("2 tests collected", content)
		self.assertIn("2 passed ✅, 0 failed ❌, 0 skipped ⏭️", content)
		self.assertIn("in 0:01:20s on 14/09/2026 at 16:58:56", content)
		self.assertEqual(content.count('title="test_sp_PSF.cpp">SP PSF</strong>'), 1)
		self.assertIn("Modes / Gradient", content)
		self.assertIn("One Degree", content)
		self.assertNotIn("[OneDegree]", content)
		self.assertNotIn("Modes/PSFTests.Gradient", content)
		self.assertIn('(2 cases)</small>', content)
		self.assertIn('<td colspan="3"><details><summary>', content)

	# --------------------------------------------------
	def test_source_file_attribute(self):
		"""Privilégie le fichier indiqué par GoogleTest et normalise les chemins Windows."""
		root = ET.fromstring('<testsuite><testcase name="Value" file="C:\\repo\\test_math.cpp"/></testsuite>')
		self.assertEqual(parse_report(root, Path("test.xml")).cases[0].file, "test_math.cpp")

	# --------------------------------------------------
	def test_readable_titles(self):
		"""Sépare les mots sans répéter les classes communes aux tests du fichier."""
		self.assertEqual(title_case("sp_PSF"), "SP PSF")
		self.assertEqual(test_title("Positions/PSFTests", "IsotropicValue"), "Positions / Isotropic Value")
		self.assertEqual(test_title("Shapes/ATrousImageTests", "ReconstructsInput"), "Shapes / Reconstructs Input")
		self.assertEqual(test_title("PSFTests", "IsotropicValue"), "Isotropic Value")

	# --------------------------------------------------
	def test_distinct_classes_with_same_method(self):
		"""Conserve une distinction lorsque deux classes du même fichier testent une méthode homonyme."""
		root = ET.fromstring('''<testsuite>
            <testcase classname="TrackAnalysisTests" name="Fit/Planar2D"/>
            <testcase classname="TrackFitTests" name="Fit/Power"/></testsuite>''')
		files = {"TrackAnalysisTests": "test_tr_Tracks.cpp", "TrackFitTests": "test_tr_Tracks.cpp"}
		content = render_report([parse_report(root, Path("test.xml"), files)])
		self.assertIn("Track Analysis / Fit / Planar 2D", content)
		self.assertIn("Track Fit / Fit / Power", content)

	# --------------------------------------------------
	def test_source_suite_mapping(self):
		"""Retrouve les suites paramétrées sans ajouter les sources à l'entrée de Doxygen."""
		files = source_files(Path(__file__).resolve().parents[2] / "Tests_CPU")
		self.assertEqual(files["PalmLocalizationTests"], "test_Palm.cpp")
		self.assertEqual(files["TrackAnalysisTests"], "test_tr_Tracks.cpp")

	# --------------------------------------------------
	def test_failed_variants_are_expanded(self):
		"""Ouvre automatiquement les variantes lorsqu'un diagnostic doit être consulté."""
		root = ET.fromstring('''<testsuite name="Suite">
            <testcase name="Value/One"><failure>Unexpected</failure></testcase>
            <testcase name="Value/Two"><skipped/></testcase></testsuite>''')
		content = render_report([parse_report(root, Path("test.xml"))])
		self.assertIn('<td colspan="3"><details open><summary>', content)
		self.assertIn("0 passed ✅, 1 failed ❌, 1 skipped ⏭️", content)

	# --------------------------------------------------
	def test_invalid_duration(self):
		"""Refuse les durées incohérentes au lieu de produire un bilan trompeur."""
		for duration in ("-1", "nan", "inf"):
			with self.subTest(duration=duration), self.assertRaises(ValueError):
				parse_report(ET.fromstring(f'<testsuite time="{duration}"/>'), Path("test.xml"))


if __name__ == "__main__":
	unittest.main()
