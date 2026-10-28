"""Research report export (Milestone 4): Markdown report and CSV of sites."""
from __future__ import annotations

import csv
import io
from datetime import datetime, timezone

from .models import CutSite, EnzymeInfo, Fragment

BIOPYTHON_DOCS = "https://biopython.org/docs/latest/api/Bio.Restriction.html"

LIMITATIONS = [
    "Predictions are in silico only; they have not been confirmed by a laboratory digest.",
    "Enzymes come from the REBASE snapshot bundled with the installed Biopython "
    "version, so the list is neither exhaustive nor guaranteed to be current.",
    "Methylation sensitivity (Dam/Dcm/CpG), star activity, buffer conditions, "
    "partial digests and DNA topology effects are not modelled.",
    "Only the selected region is scanned: recognition sites that straddle the region "
    "boundary are missed, and enzymes that cut outside their recognition site "
    "(e.g. Type IIS) may report cuts near the edges incorrectly.",
    "The virtual gel is an empirical model of agarose migration, not a measurement; "
    "real band positions depend on buffer, voltage, DNA topology and loading.",
    "Cut positions are top-strand cuts (position of the first base after the cut, 1-based); "
    "fragment sizes ignore overhang lengths.",
    "Ambiguous bases (N, R, Y, ...) are matched using IUPAC rules and can inflate or hide sites.",
    "For research and education only; not for clinical or diagnostic use.",
]


def sites_to_csv(sites: list[CutSite]) -> str:
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["enzyme", "cut_position_bp", "recognition_site", "cut_pattern"])
    for s in sites:
        writer.writerow([s.enzyme, s.position, s.site, s.cut_pattern])
    return buf.getvalue()


def build_markdown_report(*, record_id: str, description: str, seq_len: int,
                          region: tuple[int, int], gc_percent: float,
                          circular: bool, enzyme_set_label: str,
                          enzymes_scanned: int, sites: list[CutSite],
                          infos: dict[str, EnzymeInfo],
                          digest_enzymes: list[str],
                          fragments: list[Fragment],
                          biopython_version: str) -> str:
    start, end = region
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        f"# Restriction Analysis Report — {record_id}",
        "",
        f"_Generated {now} by Restriction Enzyme Explorer (Biopython {biopython_version})_",
        "",
        "## Sequence",
        f"- **Identifier:** {record_id}",
        f"- **Description:** {description or '—'}",
        f"- **Total length:** {seq_len:,} bp",
        f"- **Selected region:** {start:,}–{end:,} ({end - start + 1:,} bp)",
        f"- **GC content (region):** {gc_percent:.1f}%",
        f"- **Topology assumed:** {'circular' if circular else 'linear'}",
        "",
        "## Analysis settings",
        f"- **Enzyme set:** {enzyme_set_label} ({enzymes_scanned} enzymes scanned)",
        f"- **Enzymes with at least one site:** {len({s.enzyme for s in sites})}",
        f"- **Total cut sites:** {len(sites)}",
        "",
        "## Detected sites",
    ]
    if sites:
        lines += ["| Enzyme | Cut position (bp) | Recognition site | Cut pattern |",
                  "|---|---|---|---|"]
        lines += [f"| {s.enzyme} | {s.position:,} | `{s.site}` | `{s.cut_pattern}` |"
                  for s in sites]
    else:
        lines.append("No sites detected with the selected enzymes in this region.")

    if digest_enzymes:
        lines += ["", f"## Simulated digest ({', '.join(digest_enzymes)})",
                  "| Fragment | Ends | Coordinates | Length (bp) |", "|---|---|---|---|"]
        lines += [f"| {i} | {f.left} \u2013 {f.right} | {f.start:,}\u2013{f.end:,} | {f.length:,} |"
                  for i, f in enumerate(fragments, 1)]

    lines += ["", "## Enzyme references",
              "| Enzyme | Site | Cut pattern | Ends | Temp (°C) | Commercial | Supplier codes |",
              "|---|---|---|---|---|---|---|"]
    for name in sorted({s.enzyme for s in sites}):
        i = infos[name]
        ends = i.overhang + (f" ({i.overhang_seq})" if i.overhang_seq else "")
        lines.append(
            f"| {i.name} | `{i.site}` | `{i.cut_pattern}` | {ends} | "
            f"{i.opt_temp if i.opt_temp is not None else '—'} | "
            f"{'yes' if i.commercial else 'no'} | {', '.join(i.suppliers) or '—'} |")

    lines += [
        "", "### Data sources",
        f"- Biopython `Bio.Restriction` documentation: {BIOPYTHON_DOCS}",
        "- Enzyme data: REBASE (Roberts RJ et al., Nucleic Acids Res. 2015), via Biopython.",
        "- Cock PJA et al. Biopython. Bioinformatics 25(11):1422–1423 (2009).",
        "", "## Limitations",
    ]
    lines += [f"- {text}" for text in LIMITATIONS]
    return "\n".join(lines) + "\n"
