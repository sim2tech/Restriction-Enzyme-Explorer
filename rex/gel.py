"""Virtual agarose-gel model (pure Python): band positions and lane layout.

Migration is an empirical model, not a physical simulation: distance is linear in
log(size) across the resolving range of the chosen agarose percentage and
compresses smoothly outside it.
"""
from __future__ import annotations

import math

# Approximate useful separation range (bp) per agarose %.
GEL_RANGES = {0.7: (800, 12000), 1.0: (500, 10000), 1.5: (200, 3000),
              2.0: (100, 2000), 3.0: (50, 1000)}
LADDERS = {
    "None": [],
    "1 kb ladder (500 to 10,000 bp)": [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500],
    "100 bp ladder (100 to 1,500 bp)": [1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
}
MAJOR_TICKS = [100, 500, 1000, 5000, 10000]
MINOR_TICKS = [200, 300, 400, 600, 700, 800, 900, 2000, 3000, 4000, 6000, 7000, 8000, 9000,
               20000]
BAND_THRESHOLD_BP = 50   # smaller fragments are hard to see on a real gel


def mobility_fraction(bp: float, percent: float = 1.0) -> float:
    """Fraction (0 = well, 1 = gel end) of the run distance travelled by `bp`."""
    lo, hi = GEL_RANGES.get(percent, GEL_RANGES[1.0])
    bp = max(1.0, float(bp))
    t = (math.log(hi) - math.log(bp)) / (math.log(hi) - math.log(lo))
    if t < 0:
        return 0.12 * math.exp(1.2 * t)
    if t > 1:
        return 0.88 + 0.10 * (1 - math.exp(-1.5 * (t - 1)))
    return 0.12 + 0.76 * t


def migration_mm(bp: float, percent: float = 1.0, length_mm: float = 100.0) -> float:
    return mobility_fraction(bp, percent) * length_mm


def band_thickness_mm(bp: float, total_bp: float) -> float:
    """Band thickness grows with the DNA mass it holds (equimolar fragments)."""
    share = min(1.0, max(0.0, bp / total_bp)) if total_bp else 0.0
    return 0.35 + 1.5 * share ** 0.6


def axis_ticks(percent: float, length_mm: float = 100.0):
    """(bp, y_mm, is_major) for the size axis."""
    out = []
    for bp in sorted(MAJOR_TICKS + MINOR_TICKS):
        out.append((bp, migration_mm(bp, percent, length_mm), bp in MAJOR_TICKS))
    return out


def lane_sizes(fragments) -> list[int]:
    return sorted((f.length for f in fragments), reverse=True)
