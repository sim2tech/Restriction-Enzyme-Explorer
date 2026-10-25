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
