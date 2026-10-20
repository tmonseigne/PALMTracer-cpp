"""Vérifie les compteurs et le rendu du rapport de couverture sans lancer de C++."""

import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from coverage_to_markdown import parse_coverage, render_report


class CoverageReportTests(unittest.TestCase):
	"""Jeu de sources temporaire pour tester aussi le filtrage des chemins."""

	# --------------------------------------------------
	def setUp(self):
		"""Prépare un cœur CPU, un prototype GPU et un fichier hors périmètre."""
		self.directory = tempfile.TemporaryDirectory()
		self.addCleanup(self.directory.cleanup)
		self.root = Path(self.directory.name).resolve()
		(self.root / "lib" / "cuda").mkdir(parents=True)
		for name in ("lib/example.cpp", "lib/cuda/example.cpp", "test.cpp"):
			(self.root / name).write_text("first\nif (x < 3 && y > 2) {}\nlast\n", encoding="utf-8")

	# --------------------------------------------------
	def test_counts_deduplicate_classes(self):
		"""Une ligne dupliquée reste couverte si une classe ou un module l'a exécutée."""
		xml = ET.fromstring('''<coverage><packages><package><classes>
			<class filename="lib/example.cpp"><lines><line number="1" hits="0"/><line number="2" hits="0"/></lines></class>
			<class filename="lib/example.cpp"><lines><line number="1" hits="4"/><line number="3" hits="1"/></lines></class>
			</classes></package></packages></coverage>''')
		files = parse_coverage(xml, self.root)
		self.assertEqual(len(files), 1)
		self.assertEqual(files[0].lines, {1: 4, 2: 0, 3: 1})
		self.assertEqual(files[0].covered, 2)

	# --------------------------------------------------
	def test_scope_and_source_roots(self):
		"""Les chemins absolus ou relatifs sont acceptés sans inclure GPU, tests ou fichiers externes."""
		for filename in ("lib/example.cpp", str(self.root / "lib/example.cpp"), "example.cpp"):
			with self.subTest(filename=filename):
				xml = ET.Element("coverage")
				ET.SubElement(ET.SubElement(xml, "sources"), "source").text = str(self.root / "lib")
				classes = ET.SubElement(ET.SubElement(ET.SubElement(xml, "packages"), "package"), "classes")
				for name in (filename, "lib/cuda/example.cpp", "test.cpp", "../../elsewhere.cpp"):
					ET.SubElement(ET.SubElement(ET.SubElement(classes, "class", filename=name), "lines"), "line", number="1", hits="1")
				self.assertEqual([file.path for file in parse_coverage(xml, self.root)], [self.root / "lib/example.cpp"])

	# --------------------------------------------------
	def test_invalid_line_counts(self):
		"""Des compteurs invalides ne doivent pas donner un pourcentage trompeur."""
		for number, hits in (("0", "1"), ("1", "-1"), ("bad", "1")):
			with self.subTest(number=number, hits=hits), self.assertRaises(ValueError):
				xml = ET.fromstring(f'''<coverage><packages><package><classes><class filename="lib/example.cpp">
					<lines><line number="{number}" hits="{hits}"/></lines></class></classes></package></packages></coverage>''')
				parse_coverage(xml, self.root)

	# --------------------------------------------------
	def test_drive_source_root(self):
		"""Accepte la racine C: exportée par OpenCppCoverage, sans la traiter comme un chemin relatif."""
		xml = ET.Element("coverage")
		base = self.root.anchor.rstrip("\\/") or "/"
		ET.SubElement(ET.SubElement(xml, "sources"), "source").text = base
		classes = ET.SubElement(ET.SubElement(ET.SubElement(xml, "packages"), "package"), "classes")
		filename = (self.root / "lib/example.cpp").relative_to(self.root.anchor).as_posix()
		ET.SubElement(ET.SubElement(ET.SubElement(classes, "class", filename=filename), "lines"), "line", number="1", hits="1")
		self.assertEqual([file.path for file in parse_coverage(xml, self.root)], [self.root / "lib/example.cpp"])

	# --------------------------------------------------
	def test_render_escapes_source(self):
		"""Le bilan compte les lignes uniques et protège le HTML du code et des libellés."""
		xml = ET.fromstring('''<coverage><packages><package><classes><class filename="lib/example.cpp">
			<lines><line number="1" hits="1"/><line number="2" hits="0"/></lines></class></classes></package></packages></coverage>''')
		report = render_report(parse_coverage(xml, self.root), self.root, "Debug <MSVC>", "date")
		self.assertIn("50.00% line coverage", report)
		self.assertIn("1/2 lines covered", report)
		self.assertIn("Debug &lt;MSVC&gt;", report)
		self.assertIn("x &lt; 3 &amp;&amp; y &gt; 2", report)
		self.assertIn('<details class="coverage-lines">', report)

	# --------------------------------------------------
	def test_empty_report(self):
		"""L'absence de ligne instrumentée n'est pas présentée comme une couverture parfaite."""
		self.assertEqual(parse_coverage(ET.fromstring("<coverage/>"), self.root), [])
		self.assertIn("N/A line coverage", render_report([], self.root, "empty", "date"))


if __name__ == "__main__":
	unittest.main()
