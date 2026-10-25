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


@dataclass(frozen=True)
class Fragment:
    start: int         # 1-based, inclusive
    end: int           # 1-based, inclusive
    length: int


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
