"""Pure geometry and layout helpers for the circular map and sequence view."""
from __future__ import annotations

import math

FEATURE_COLORS = {
    "CDS": "#4F6D84", "gene": "#4F6D84", "promoter": "#2F8F83", "terminator": "#C86B6B",
    "rep_origin": "#E0B83A", "RBS": "#E9A35B", "protein_bind": "#D98B3A",
    "regulatory": "#D98B3A", "primer_bind": "#8E6BBF", "misc_feature": "#9AA3A8",
    "tag": "#E9A35B", "sig_peptide": "#6FA8A0",
}
DEFAULT_FEATURE_COLOR = "#8C8F94"


def feature_color(ftype: str) -> str:
    return FEATURE_COLORS.get(ftype, DEFAULT_FEATURE_COLOR)


def text_on(hex_colour: str) -> str:
    """Readable text colour (not pure black/white) for a given background."""
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (1, 3, 5))
    return "#2B2926" if 0.299 * r + 0.587 * g + 0.114 * b > 150 else "#F2EDE3"


def bp_to_angle(pos: float, length: int) -> float:
    """Radians clockwise from 12 o'clock for the left edge of base `pos` (1-based)."""
    return 2 * math.pi * (pos - 1) / length


def polar(r: float, theta: float) -> tuple[float, float]:
    return r * math.sin(theta), r * math.cos(theta)


def arc_points(a0: float, a1: float, r: float) -> list[tuple[float, float]]:
    steps = max(2, int(abs(a1 - a0) / math.radians(2)) + 2)
    return [polar(r, a0 + (a1 - a0) * i / (steps - 1)) for i in range(steps)]


def feature_polygon(a0: float, a1: float, r0: float, r1: float, strand: int = 0):
    """Closed outline (xs, ys) of a thick arc from angle a0 to a1; pointed at the strand end."""
    span = a1 - a0
    head = min(math.radians(6), span * 0.45) if strand else 0.0
    rm = (r0 + r1) / 2
    if strand < 0:
        body = a0 + head
        pts = [polar(rm, a0)] + arc_points(body, a1, r1) + arc_points(body, a1, r0)[::-1]
    else:
        body = a1 - head
        pts = arc_points(a0, body, r1) + ([polar(rm, a1)] if strand > 0 else []) \
            + arc_points(a0, body, r0)[::-1]
    pts.append(pts[0])
    return [p[0] for p in pts], [p[1] for p in pts]


def _bins(feature, length: int, n: int = 360) -> set[int]:
    out: set[int] = set()
    for s, e in feature.segments:
        b0, b1 = int((s - 1) * n / length), int((e - 1) * n / length)
        if b1 < b0:  # wraps the origin within one part
            b1 += n
        out.update(b % n for b in range(b0, b1 + 1))
    return out


def assign_lanes(features, length: int) -> list[int]:
    """Lane index per feature so overlapping features never share a lane."""
    lanes: list[set[int]] = []
    result = [0] * len(features)
    order = sorted(range(len(features)), key=lambda i: features[i].start)
    for i in order:
        covered = _bins(features[i], length)
        for k, used in enumerate(lanes):
            if not used & covered:
                used |= covered
                result[i] = k
                break
        else:
            lanes.append(set(covered))
            result[i] = len(lanes) - 1
    return result


def spread_1d(values: list[float], spacing: float, lo: float, hi: float) -> list[float]:
    """Move values (kept in order, high to low) so neighbours are >= spacing apart."""
    n = len(values)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: -values[i])
    y = [min(max(values[i], lo), hi) for i in order]
    for k in range(1, n):
        y[k] = min(y[k], y[k - 1] - spacing)
    if y[-1] < lo:                      # ran out of room at the bottom: pack upwards
        y[-1] = lo
        for k in range(n - 2, -1, -1):
            y[k] = max(y[k], y[k + 1] + spacing)
    out = [0.0] * n
    for k, i in enumerate(order):
        out[i] = y[k]
    return out


def nice_step(length: int, target: int = 8) -> int:
    for step in (10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000,
                 100000, 250000, 500000, 1000000):
        if length / step <= target:
            return step
    return 1000000


HIDDEN_FEATURE_TYPES = {"misc_binding", "old_sequence", "misc_difference", "repeat_region",
                        "mobile_element", "mat_peptide", "conflict", "variation", "STS", "gap"}
FEATURE_PRIORITY = {"CDS": 0, "rep_origin": 0, "promoter": 1, "terminator": 1, "RBS": 1,
                    "regulatory": 2, "protein_bind": 2, "sig_peptide": 3, "misc_feature": 3}
MAX_LANES = 3


def split_lanes(features, length: int, max_lanes: int = MAX_LANES):
    """Pack features into at most max_lanes lanes, most important first.

    Returns ([(feature, lane), ...] in original order, number of features dropped).
    """
    lanes: list[set[int]] = []
    placed: dict[int, int] = {}
    order = sorted(range(len(features)),
                   key=lambda i: (FEATURE_PRIORITY.get(features[i].type, 4), features[i].start))
    for i in order:
        covered = _bins(features[i], length)
        for k, used in enumerate(lanes):
            if not used & covered:
                used |= covered
                placed[i] = k
                break
        else:
            if len(lanes) < max_lanes:
                lanes.append(set(covered))
                placed[i] = len(lanes) - 1
    kept = [(features[i], placed[i]) for i in range(len(features)) if i in placed]
    return kept, len(features) - len(kept)


def text_angle(theta: float) -> float:
    """Rotation (degrees, clockwise) that keeps text along the ring and upright."""
    a = math.degrees(theta) % 360
    if a <= 90:
        return a
    return a - 180 if a <= 270 else a - 360


def bezier(p0, p1, p2, n: int = 14):
    """Points along a quadratic Bezier curve from p0 to p2 with control point p1."""
    pts = []
    for i in range(n + 1):
        t = i / n
        pts.append(((1 - t) ** 2 * p0[0] + 2 * (1 - t) * t * p1[0] + t * t * p2[0],
                    (1 - t) ** 2 * p0[1] + 2 * (1 - t) * t * p1[1] + t * t * p2[1]))
    return pts
