"""IUPAC-aware motif search (Milestone 1: motif detection). No third-party imports."""
from __future__ import annotations

import re
from dataclasses import dataclass

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S",
    "W": "W", "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}


@dataclass(frozen=True)
class MotifHit:
    strand: str     # "+" or "-"
    start: int      # 1-based, forward-strand coordinates
    end: int
    matched: str    # forward-strand text that matched


def clean_motif(motif: str) -> str:
    motif = "".join(motif.split()).upper()
    bad = sorted(set(motif) - set(IUPAC))
    if not motif:
        raise ValueError("Motif is empty.")
    if bad:
        raise ValueError(f"Motif has non-IUPAC DNA characters: {', '.join(bad)}")
    return motif


def reverse_complement(motif: str) -> str:
    return "".join(COMPLEMENT[c] for c in reversed(motif))


def _regex(motif: str) -> re.Pattern:
    body = "".join(c if len(IUPAC[c]) == 1 else f"[{IUPAC[c]}]" for c in motif)
    return re.compile(f"(?=({body}))")  # lookahead => overlapping matches


def find_motif(seq: str, motif: str, start: int = 1, end: int | None = None,
               both_strands: bool = True) -> list[MotifHit]:
    """Find all (overlapping) occurrences of an IUPAC motif inside seq[start..end]."""
    motif = clean_motif(motif)
    end = len(seq) if end is None else end
    region = seq[start - 1:end]
    hits: list[MotifHit] = []
    strands = [("+", motif)]
    rc = reverse_complement(motif)
    if both_strands and rc != motif:  # palindromes would just duplicate hits
        strands.append(("-", rc))
    for strand, pattern in strands:
        for m in _regex(pattern).finditer(region):
            s = m.start() + start
            hits.append(MotifHit(strand, s, s + len(motif) - 1, m.group(1)))
    return sorted(hits, key=lambda h: (h.start, h.strand))
