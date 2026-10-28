"""Plain data containers shared across the package (no third-party imports)."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CutSite:
    """One predicted top-strand cut.

    ``position`` follows the Biopython convention: the 1-based position of the
    first base *after* the cut, in the coordinates of the full sequence.
    """

    enzyme: str
    position: int
    site: str          # recognition sequence, e.g. "GAATTC"
    cut_pattern: str   # e.g. "G^AATT_C"
    # Position of the bottom-strand cut in the same convention (palindromic
    # enzymes only; None when it cannot be stated reliably).
    bottom_position: int | None = None


@dataclass(frozen=True)
class Fragment:
    start: int         # 1-based, inclusive (may be > end for a circular wrap-around)
    end: int           # 1-based, inclusive
    length: int
    left: str = ""     # enzyme(s) cutting at the left end, or "start"
    right: str = ""    # enzyme(s) cutting at the right end, or "end"


@dataclass(frozen=True)
class EnzymeInfo:
    name: str
    site: str
    cut_pattern: str
    size: int                  # recognition-site length in bp
    overhang: str              # "blunt", "5' overhang", "3' overhang", "unknown"
    overhang_seq: str
    suppliers: tuple
    commercial: bool
    opt_temp: object = None    # optimal incubation temperature (°C) if documented
    inact_temp: object = None  # heat-inactivation temperature (°C) if documented
    isoschizomers: tuple = ()  # other enzymes recognising the same site


@dataclass(frozen=True)
class Feature:
    """An annotated feature (gene, promoter, origin ...), 1-based inclusive."""

    name: str
    type: str
    start: int
    end: int
    strand: int = 0                 # +1, -1 or 0 (unstranded)
    parts: tuple = ()               # ((start, end), ...) for joined/wrapping features

    @property
    def segments(self) -> tuple:
        return self.parts or ((self.start, self.end),)

    @property
    def length(self) -> int:
        return sum(e - s + 1 for s, e in self.segments)
