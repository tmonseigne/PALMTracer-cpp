"""Convertit les resultats XML GoogleTest ou CTest en rapport Markdown pour Doxygen."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta
from html import escape
import math
import re
from pathlib import Path
import xml.etree.ElementTree as ET


@dataclass(frozen=True)
class TestCase:
    """Resultat d'un cas individuel, avec son nom parametre et son diagnostic."""

    suite: str
    name: str
    status: str
    duration: float
    details: str
    file: str = ""


@dataclass(frozen=True)
class TestReport:
    """Execution provenant d'un fichier XML, sans resultat deduit d'une autre execution."""

    source: Path
    timestamp: str
    duration: float
    cases: tuple[TestCase, ...]


def duration_seconds(value: str) -> float:
    """Valide une duree en secondes pour eviter les valeurs negatives ou non finies."""
    result = float(value)
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"Duree XML invalide : {value}")
    return result


def source_files(directory: Path) -> dict[str, str]:
    """Associe les suites GoogleTest a leur fichier, meme avec un ancien XML sans attribut file."""
    suites = {}
    for path in sorted(directory.glob("test_*.cpp")):
        content = path.read_text(encoding="utf-8-sig")
        for suite in re.findall(r"\bTEST(?:_P|_F)?\s*\(\s*(\w+)\s*,", content):
            suites[suite] = path.name
    return suites


def parse_report(root: ET.Element, source: Path, files: dict[str, str] | None = None) -> TestReport:
    """
    Extrait les cas d'un rapport GoogleTest ou JUnit CTest.

    :param root: Racine du document XML.
    :param source: Chemin utilise pour identifier l'execution dans le rapport.
    :param files: Correspondance facultative entre suites et fichiers sources.
    :return: Execution avec un statut explicite pour chaque cas.
    """
    if root.tag not in {"testsuite", "testsuites"}:
        raise ValueError(f"Format de rapport XML inconnu : {source}")
    parents = {child: parent for parent in root.iter() for child in parent}
    cases = []
    for element in root.iter("testcase"):
        name = element.attrib["name"]
        suite = element.get("classname", "")
        # CTest duplique le nom complet dans classname ; GoogleTest separe classe et methode.
        if "." in name and (not suite or suite == name):
            suite, name = name.split(".", 1)
        if not suite:
            suite = parents[element].get("name", "Tests")
        failure = element.find("failure")
        error = element.find("error")
        skipped = element.find("skipped")
        details = ""
        if error is not None or failure is not None:
            status = "Erreur" if error is not None else "Echec"
            diagnostic = error if error is not None else failure
            details = diagnostic.text or diagnostic.get("message", "")
            output = element.findtext("system-out", "")
            if output:
                details += "\n\n" + output
        elif skipped is not None:
            status = "Ignore"
            details = skipped.text or skipped.get("message", "")
        elif element.get("status") in {"notrun", "disabled"} or element.get("result") == "suppressed":
            status = "Non execute"
        else:
            status = "Reussi"
        filename = element.get("file", "").replace("\\", "/").rsplit("/", 1)[-1]
        if not filename:
            filename = (files or {}).get(suite.rsplit("/", 1)[-1], suite)
        cases.append(TestCase(suite, name, status, duration_seconds(element.get("time", "0")), details.strip(), filename))
    total_time = duration_seconds(root.get("time", str(sum(case.duration for case in cases))))
    return TestReport(source, root.get("timestamp", "Non fourni"), total_time, tuple(cases))


def status_label(status: str) -> str:
    """Produit un statut visible et comprehensible sans dependance a la couleur."""
    icons = {"Reussi": ("✅", "Passed"), "Echec": ("❌", "Failed"), "Erreur": ("❌", "Error"),
             "Ignore": ("⏭️", "Skipped"), "Non execute": ("⏭️", "Not run")}
    icon, label = icons[status]
    return f'<span class="test-status" title="{label}">{icon} {label}</span>'


def title_case(value: str) -> str:
    """Separe snake_case et CamelCase en conservant les acronymes scientifiques."""
    value = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1 \2", value)
    value = re.sub(r"([a-z0-9])([A-Z])", r"\1 \2", value).replace("_", " ")
    value = re.sub(r"([a-z])([0-9])", r"\1 \2", value)
    value = re.sub(r"\b([23]) D\b", r"\1D", value)
    acronyms = {"sp", "tr", "loc", "psf", "fft", "msd", "cpu", "gpu", "xy", "xyz", "dft", "id"}
    return " ".join(word.upper() if word.lower() in acronyms else word[:1].upper() + word[1:] for word in value.split())


def test_title(suite: str, name: str, include_suite: bool = False) -> str:
    """Conserve le prefixe de parametrage mais masque la classe, sauf en cas d'ambiguite."""
    parts = suite.split("/")
    prefix = [title_case(part) for part in parts[:-1]]
    if include_suite:
        prefix.append(title_case(re.sub(r"Tests?$", "", parts[-1])))
    return " / ".join([*prefix, title_case(name)])


def render_rows(cases: list[TestCase], title: str = "", parameters_only: bool = False) -> list[str]:
    """Construit les lignes et conserve les diagnostics dans des sections repliables."""
    lines = []
    for case in cases:
        name, separator, parameter = case.name.partition("/")
        label = title or test_title(case.suite, name)
        if separator:
            parameter = title_case(parameter)
            label = parameter if parameters_only else f"{label} / {parameter}"
        lines.append(f'<tr><td>{escape(label)}</td>'
                     f'<td>{status_label(case.status)}</td><td>{case.duration:.3f} s</td></tr>')
        if case.details:
            lines.append(f'<tr><td colspan="3"><details open><summary>Diagnostic</summary>'
                         f'<pre>{escape(case.details)}</pre></details></td></tr>')
    return lines


def render_report(reports: list[TestReport]) -> str:
    """
    Produit une page Doxygen avec bilan et tableaux regroupes par fichier et fonction.

    :param reports: Executions a presenter independamment, notamment par plateforme.
    :return: Markdown UTF-8, avec tableaux HTML et diagnostics echappes.
    """
    lines = ["# Test report {#test_report}", "", r"\htmlonly", '<div class="test-report">']
    header = '<thead><tr><th>Test / Parameters</th><th>Status</th><th>Duration</th></tr></thead>'
    for report in reports:
        counts = Counter(case.status for case in report.cases)
        label = f"{report.source.parent.name}/{report.source.name}"
        failed = counts['Echec'] + counts['Erreur']
        skipped = counts['Ignore'] + counts['Non execute']
        duration = str(timedelta(seconds=int(report.duration)))
        try:
            timestamp = datetime.fromisoformat(report.timestamp.replace("Z", "+00:00"))
            date = timestamp.strftime("on %d/%m/%Y at %H:%M:%S %z").rstrip()
        except ValueError:
            date = f"(timestamp: {escape(report.timestamp)})"
        lines.append('<section class="test-run">')
        if len(reports) > 1:
            lines.append(f'<h2>{escape(label)}</h2>')
        lines.extend([
            f'<p class="test-summary"><strong>{len(report.cases)} tests collected</strong>, '
            f'{counts["Reussi"]} passed ✅, {failed} failed ❌, {skipped} skipped ⏭️ '
            f'in {duration}s {date}</p>',
        ])
        files = defaultdict(list)
        for case in report.cases:
            files[case.file or case.suite].append(case)
        for filename, cases in files.items():
            failed = sum(case.status in {"Echec", "Erreur"} for case in cases)
            passed = sum(case.status == "Reussi" for case in cases)
            lines.extend([
                '<details class="test-file" open><summary>'
                f'<strong title="{escape(filename, quote=True)}">{escape(title_case(Path(filename).stem.removeprefix("test_")))}</strong>'
                f'<span>{passed}/{len(cases)} passed'
                f'{" · ❌ " + str(failed) + " failed" if failed else ""}</span></summary>',
                '<div class="test-table-scroll"><table class="test-results">', header, '<tbody>',
            ])
            groups = defaultdict(list)
            for case in cases:
                groups[(case.suite, case.name.partition("/")[0])].append(case)
            titles = Counter(test_title(suite, name) for suite, name in groups)
            for (suite, name), variants in groups.items():
                title = test_title(suite, name)
                if titles[title] > 1:
                    title = test_title(suite, name, include_suite=True)
                if len(variants) == 1:
                    lines.extend(render_rows(variants, title))
                    continue
                expanded = ' open' if any(case.status in {"Echec", "Erreur"} for case in variants) else ''
                passed = sum(case.status == "Reussi" for case in variants)
                lines.extend([
                    f'<tr class="test-group"><td colspan="3"><details{expanded}><summary>'
                    f'<span class="test-group-name">{escape(title)} '
                    f'<small>({len(variants)} cases)</small></span>'
                    f'<span>✅ {passed}/{len(variants)}</span>'
                    f'<span>{sum(case.duration for case in variants):.3f} s</span></summary>',
                    '<table class="test-results test-variants">', header, '<tbody>',
                    *render_rows(variants, parameters_only=True), '</tbody></table></details></td></tr>',
                ])
            lines.extend(["</tbody></table></div></details>"])
        lines.append('</section>')
    lines.extend(['</div>', r"\endhtmlonly", ""])
    return "\n".join(lines)


def main() -> None:
    """Genere le rapport depuis un ou plusieurs XML, sans lancer les tests."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, help="Rapports XML GoogleTest ou CTest/JUnit.")
    parser.add_argument("--output", type=Path, default=Path("Docs/reports/test_report.md"))
    parser.add_argument("--source-dir", type=Path, default=Path(__file__).resolve().parents[2] / "Tests_CPU")
    args = parser.parse_args()
    files = source_files(args.source_dir)
    reports = [parse_report(ET.parse(path).getroot(), path, files) for path in args.inputs]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(render_report(reports))
    print(f"Rapport genere : {args.output} ({len(reports)} execution(s)).")


if __name__ == "__main__":
    main()
