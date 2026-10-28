"""Fragment-size calculation for a digest (top-strand cuts). No third-party imports."""
from __future__ import annotations

from typing import Iterable

from .models import Fragment


def compute_fragments(cut_positions: Iterable[int], start: int, end: int,
                      circular: bool = False) -> list[Fragment]:
    """Fragments of seq[start..end] given Biopython-style cut positions.

    A cut position p means the cut lies between bases p-1 and p.
    """
    length = end - start + 1
    if not circular:
        cuts = sorted({p for p in cut_positions if start < p <= end})
        starts = [start] + cuts
        ends = [c - 1 for c in cuts] + [end]
        return [Fragment(s, e, e - s + 1) for s, e in zip(starts, ends)]

    cuts = sorted({((p - start) % length) + start for p in cut_positions})
    if not cuts:
        return [Fragment(start, end, length)]
    n = len(cuts)
    frags = []
    for i, c in enumerate(cuts):
        nxt = cuts[(i + 1) % n]
        size = ((nxt - c) % length) or length
        frag_end = ((nxt - 1 - start) % length) + start
        frags.append(Fragment(c, frag_end, size))
    return frags


def digest_fragments(cuts, start: int, end: int, circular: bool = False) -> list[Fragment]:
    """Fragments of a digest with the enzyme(s) at each end.

    ``cuts`` are CutSite-like objects (``.position``, ``.enzyme``). Fragments are
    returned in sequence order.
    """
    from dataclasses import replace

    length = end - start + 1

    def norm(p: int) -> int:
        return ((p - start) % length) + start if circular else p

    names: dict[int, set] = {}
    for c in cuts:
        names.setdefault(norm(c.position), set()).add(c.enzyme)

    def label(pos: int, default: str) -> str:
        return " / ".join(sorted(names[pos])) if pos in names else default

    out = []
    for f in compute_fragments([c.position for c in cuts], start, end, circular):
        left = label(f.start, "start")
        right_pos = norm(f.end + 1) if circular else f.end + 1
        right = label(right_pos, "end")
        out.append(replace(f, left=left, right=right))
    return out


def fragment_sequence(seq: str, fragment: Fragment) -> str:
    """Top-strand sequence of a fragment (handles circular wrap-around)."""
    if fragment.start <= fragment.end:
        return seq[fragment.start - 1:fragment.end]
    return seq[fragment.start - 1:] + seq[:fragment.end]


def format_coordinates(fragment: Fragment) -> str:
    return f"{fragment.start:,}\u2013{fragment.end:,}"
