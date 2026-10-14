"""Convertit un export Cobertura en rapport de couverture des lignes pour Doxygen."""

from __future__ import annotations

import argparse
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import datetime
from html import escape
from pathlib import Path


@dataclass(frozen=True)
class FileCoverage:
	"""
	Lignes instrumentées d'un fichier, dédupliquées entre classes et modules.

	:ivar path: Chemin absolu du fichier source du cœur CPU.
	:ivar lines: Nombre d'exécutions par numéro de ligne instrumentée.
	"""

	path: Path
	lines: dict[int, int]

	# --------------------------------------------------
	@property
	def covered(self) -> int:
		"""Nombre de lignes exécutées au moins une fois."""
		return sum(hits > 0 for hits in self.lines.values())


# --------------------------------------------------
def source_path(filename: str, root: Path, sources: list[str]) -> Path | None:
	"""Résout les chemins Cobertura sans lire de fichier extérieur au cœur CPU."""
	name = Path(filename.replace("\\", "/"))
	candidates = [root / name]
	for source in sources:
		base = source.replace("\\", "/")
		# OpenCppCoverage exporte parfois la racine du disque comme C: au lieu de C:/.
		if len(base) == 2 and base[1] == ":":
			base += "/"
		candidates.append(root / base / name)
	for candidate in candidates:
		path = candidate.resolve()
		if path.is_relative_to(root / "lib") and not path.is_relative_to(root / "lib" / "cuda") and path.is_file():
			return path
	return None


# --------------------------------------------------
def parse_coverage(xml: ET.Element, root: Path) -> list[FileCoverage]:
	"""Lit les vrais compteurs par ligne, sans faire confiance aux pourcentages arrondis."""
	root = root.resolve()
	sources = [node.text for node in xml.findall("./sources/source") if node.text]
	files: dict[Path, dict[int, int]] = {}
	for cls in xml.findall("./packages/package/classes/class"):
		path = source_path(cls.get("filename", ""), root, sources)
		if path is None:
			continue
		lines = files.setdefault(path, {})
		for node in cls.findall("./lines/line"):
			number, hits = int(node.attrib["number"]), int(node.attrib["hits"])
			if number <= 0 or hits < 0:
				raise ValueError("Line numbers must be positive and hit counts nonnegative")
			lines[number] = max(lines.get(number, 0), hits)
	return [FileCoverage(path, lines) for path, lines in sorted(files.items()) if lines]


# --------------------------------------------------
def render_report(files: list[FileCoverage], root: Path, label: str, timestamp: str) -> str:
	"""Produit le bilan et les tableaux de lignes non couvertes, sans seuil de réussite."""
	covered, total = sum(file.covered for file in files), sum(len(file.lines) for file in files)
	summary_rate = f"{100 * covered / total:.2f}%" if total else "N/A"
	lines = ["# Coverage report {#coverage_report}", "", r"\htmlonly", '<div class="test-report coverage-report">',
			 (f'<p class="test-summary"><strong>{summary_rate} line coverage</strong> · {covered:,}/{total:,} lines covered · '
			  f'{total - covered:,} uncovered · {len(files)} files<br>{escape(label)} · {escape(timestamp)}</p>'),
			 ('<p>CPU library only (<code>lib/</code>, excluding <code>lib/cuda/</code>). '
			  'Counts refer to instrumented executable lines, not all source lines. Branch coverage is not measured. No target threshold.</p>'),
			 '<div class="test-table-scroll"><table class="test-results coverage-results">',
			 '<thead><tr><th>File / Uncovered Lines</th><th>Line Coverage</th><th>Covered / Total</th></tr></thead><tbody>']
	for file in files:
		name = file.path.relative_to(root).as_posix()
		rate = 100 * file.covered / len(file.lines)
		missing = sorted(number for number, hits in file.lines.items() if hits == 0)
		source = file.path.read_text(encoding="utf-8-sig").splitlines()
		lines.append(f'<tr><td><strong>{escape(name)}</strong>')
		if missing:
			lines.append(f'<details class="coverage-lines"><summary>{len(missing)} uncovered lines</summary><pre>')
			for number in missing:
				text = source[number - 1] if number <= len(source) else "[source line unavailable]"
				lines.append(f'<span><b>{number:5}</b>  {escape(text)}</span>')
			lines.append('</pre></details>')
		lines.append(f'</td><td><span>{rate:.2f}%</span><progress max="100" value="{rate:.4f}" '
					 f'aria-label="{escape(name, quote=True)} line coverage">{rate:.2f}%</progress></td>'
					 f'<td>{file.covered:,} / {len(file.lines):,}</td></tr>')
	lines.extend(['</tbody></table></div>', '</div>', r"\endhtmlonly", ""])
	return "\n".join(lines)


# --------------------------------------------------
def main():
	"""Génère la page depuis un export existant, sans lancer d'analyse ni de tests."""
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("input", type=Path, help="Export Cobertura XML, notamment OpenCppCoverage.")
	parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
	parser.add_argument("--output", type=Path, default=Path("Docs/reports/coverage_report.md"))
	parser.add_argument("--label", default="Local run", help="Plateforme, compilateur et configuration mesures.")
	args = parser.parse_args()
	root = args.root.resolve()
	xml = ET.parse(args.input).getroot()
	files = parse_coverage(xml, root)
	if not files:
		parser.error("No instrumented CPU library source found in this report")
	try:
		timestamp = datetime.fromtimestamp(float(xml.attrib["timestamp"])).astimezone().strftime("%d/%m/%Y at %H:%M:%S %z")
	except (KeyError, ValueError, OverflowError, OSError):
		timestamp = "Timestamp unavailable"
	args.output.parent.mkdir(parents=True, exist_ok=True)
	args.output.write_text(render_report(files, root, args.label, timestamp), encoding="utf-8", newline="\n")
	print(f"Coverage report generated: {args.output} ({len(files)} files)")


if __name__ == "__main__":
	main()
