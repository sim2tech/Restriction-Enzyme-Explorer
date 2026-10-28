"""Circular plasmid map: feature arcs, curved leader lines, labelled cut sites (Plotly)."""
from __future__ import annotations

import math
from collections import Counter

import plotly.graph_objects as go

from .geometry import (MAX_LANES, bezier, bp_to_angle, feature_color, feature_polygon,
                       nice_step, polar, split_lanes, spread_1d, text_angle, text_on)
from .models import CutSite, Feature
from .plotting import FONT, PINK, THEMES

RING_R = 0.80
LANE_W = 0.095
LANE_GAP = 0.018
LABEL_X = 1.50
MAX_LABELS_PER_SIDE = 28
CHAR_UNITS = 0.026       # rough width of one 11px character in plot units
MIN_ARC = 0.03           # radians; single-base features are widened to this


def _generic(f: Feature) -> bool:
    return f.name == f.type


def build_circular_map(seq_id: str, length: int, topology: str, gc: float,
                       sites: list[CutSite], features: list[Feature], *,
                       highlight: str | None = None, labelled: set | None = None,
                       theme: str = "light", show_features: bool = True,
                       cards: dict | None = None) -> go.Figure:
    c = THEMES.get(theme, THEMES["light"])
    fig = go.Figure()
    kept, _dropped = split_lanes(list(features), length, MAX_LANES) if show_features else ([], 0)
    n_lanes = (max(lane for _, lane in kept) + 1) if kept else 0
    r_out = RING_R + LANE_GAP + n_lanes * (LANE_W + LANE_GAP) + 0.02
    annotations: list[dict] = []

    # Backbone ring.
    ring = [polar(RING_R, 2 * math.pi * i / 240) for i in range(241)]
    fig.add_trace(go.Scatter(x=[p[0] for p in ring], y=[p[1] for p in ring], mode="lines",
                             line=dict(color=c["muted"], width=2), hoverinfo="skip",
                             showlegend=False))
    # Coordinate ticks with labels that follow the ring.
    step = nice_step(length)
    tick_x, tick_y = [], []
    for pos in [1] + list(range(step, length + 1, step)):
        th = bp_to_angle(pos, length)
        a, b = polar(RING_R - 0.03, th), polar(RING_R, th)
        tick_x += [a[0], b[0], None]
        tick_y += [a[1], b[1], None]
        q = polar(RING_R - 0.10, th)
        annotations.append(dict(x=q[0], y=q[1], xref="x", yref="y", text=f"{pos:,}",
                                showarrow=False, textangle=text_angle(th),
                                font=dict(family=FONT, size=10, color=c["muted"])))
    fig.add_trace(go.Scatter(x=tick_x, y=tick_y, mode="lines", hoverinfo="skip",
                             line=dict(color=c["muted"], width=1), showlegend=False))

    # Feature arcs. Names that fit are written on the arc; the rest go in the label columns.
    right_labels, left_labels, seen = [], [], set()
    for f, lane in kept:
        r0 = RING_R + LANE_GAP + lane * (LANE_W + LANE_GAP)
        r1 = r0 + LANE_W
        rm = (r0 + r1) / 2
        colour = feature_color(f.type)
        segs = f.segments
        hover = (f"<b>{f.name}</b><br>{f.type}"
                 f"{' · complement' if f.strand < 0 else ''}<br>"
                 f"{f.start:,}\u2013{f.end:,} ({f.length:,} bp)")
        for k, (s, e) in enumerate(segs):
            a0, a1 = bp_to_angle(s, length), bp_to_angle(e + 1, length)
            if a1 < a0:
                a1 += 2 * math.pi
            if a1 - a0 < MIN_ARC:
                mid_a = (a0 + a1) / 2
                a0, a1 = mid_a - MIN_ARC / 2, mid_a + MIN_ARC / 2
            strand = f.strand if ((f.strand > 0 and k == len(segs) - 1) or
                                   (f.strand < 0 and k == 0)) else 0
            xs, ys = feature_polygon(a0, a1, r0, r1, strand)
            fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", fill="toself", fillcolor=colour,
                                     line=dict(color=colour, width=1), text=hover,
                                     hoveron="fills", hoverinfo="text", showlegend=False))
        if _generic(f) or f.name in seen:
            continue
        seen.add(f.name)
        s0, e0 = max(segs, key=lambda t: (t[1] - t[0]) % length)
        a0, a1 = bp_to_angle(s0, length), bp_to_angle(e0 + 1, length)
        if a1 < a0:
            a1 += 2 * math.pi
        mid = (a0 + a1) / 2
        if rm * (a1 - a0) >= CHAR_UNITS * len(f.name) + 0.05:
            x, y = polar(rm, mid)
            annotations.append(dict(x=x, y=y, xref="x", yref="y", text=f.name, showarrow=False,
                                    textangle=text_angle(mid),
                                    font=dict(family=FONT, size=11, color=text_on(colour))))
        else:
            ax, ay = polar(r1 + 0.005, mid)
            entry = dict(y=ay, anchor=(ax, ay), mid=mid,
                         text=f'<span style="color:{colour}">\u25A0</span> {f.name}', hl=False)
            (right_labels if ax >= 0 else left_labels).append(entry)

    # Cut sites: short ticks on the ring; labelled ones also get a name and a leader curve.
    counts = Counter(s.enzyme for s in sites)
    by_pos: dict[int, list[CutSite]] = {}
    for s in sites:
        by_pos.setdefault(s.position, []).append(s)
    tick_x, tick_y, hv_x, hv_y, hv_t, hl_x, hl_y = [], [], [], [], [], [], []
    cd_x, cd_y, cd_c = [], [], []     # hover points that carry an enzyme card
    for pos, group in by_pos.items():
        th = bp_to_angle(pos, length)
        a, b = polar(RING_R - 0.035, th), polar(RING_R + 0.035, th)
        tick_x += [a[0], b[0], None]
        tick_y += [a[1], b[1], None]
        names = sorted({s.enzyme for s in group})
        hx, hy = polar(RING_R, th)
        n = counts[names[0]]
        shown = [n_ for n_ in names if labelled is None or n_ in labelled or n_ == highlight]
        card = (cards or {}).get((shown[0], pos)) if shown else None
        if card:
            cd_x.append(hx); cd_y.append(hy); cd_c.append(card)
        else:
            hv_x.append(hx); hv_y.append(hy)
            hv_t.append(f"<b>{' · '.join(names)}</b><br>cuts after {pos - 1:,}<br>"
                        f"{'unique cutter' if n == 1 else f'{n} sites'}<br>{group[0].cut_pattern}")
        if not shown:
            continue
        ax, ay = polar(r_out, th)
        is_hl = highlight in shown
        if is_hl:
            hl_x.append(ax); hl_y.append(ay)
        entry = dict(y=ay, anchor=(ax, ay), mid=th, text=" \u00B7 ".join(shown) + f" ({pos:,})",
                     hl=is_hl, base=polar(RING_R + 0.035, th), card=card)
        (right_labels if ax >= 0 else left_labels).append(entry)
    cut_colour = c["cut"]
    fig.add_trace(go.Scatter(x=tick_x, y=tick_y, mode="lines", hoverinfo="skip",
                             line=dict(color=cut_colour, width=2), showlegend=False))
    fig.add_trace(go.Scatter(x=hv_x, y=hv_y, mode="markers", text=hv_t, hoverinfo="text",
                             marker=dict(size=12, color="rgba(0,0,0,0)"), showlegend=False))

    # Label columns with curved leaders.
    spacing = 0.085
    for labels, sign, anchor in ((right_labels, 1, "middle right"), (left_labels, -1, "middle left")):
        labels = labels[:MAX_LABELS_PER_SIDE]
        if not labels:
            continue
        ys = spread_1d([e["y"] for e in labels], spacing, -1.45, 1.45)
        lead = {False: ([], []), True: ([], [])}
        txt = {False: ([], [], []), True: ([], [], [])}
        stem_x, stem_y = [], []
        for e, y in zip(labels, ys):
            ax, ay = e["anchor"]
            k = (r_out + 0.32) / max(r_out, 1e-6)
            end = (sign * (LABEL_X - 0.04), y)
            for px, py in bezier((ax, ay), (ax * k, ay * k), end):
                lead[e["hl"]][0].append(px); lead[e["hl"]][1].append(py)
            lead[e["hl"]][0].append(None); lead[e["hl"]][1].append(None)
            txt[e["hl"]][0].append(sign * LABEL_X); txt[e["hl"]][1].append(y)
            txt[e["hl"]][2].append(e["text"])
            if e.get("card"):
                for k_ in range(max(1, int(CHAR_UNITS * 1.1 * len(e["text"]) / 0.12))):
                    cd_x.append(sign * (LABEL_X + 0.12 * k_)); cd_y.append(y); cd_c.append(e["card"])
            if "base" in e:
                stem_x += [e["base"][0], ax, None]; stem_y += [e["base"][1], ay, None]
        fig.add_trace(go.Scatter(x=stem_x, y=stem_y, mode="lines", hoverinfo="skip",
                                 line=dict(color=c["stem"], width=0.8), showlegend=False))
        for is_hl in (False, True):
            if not txt[is_hl][0]:
                continue
            fig.add_trace(go.Scatter(x=lead[is_hl][0], y=lead[is_hl][1], mode="lines",
                                     hoverinfo="skip", showlegend=False,
                                     line=dict(color=c["pink_text"] if is_hl else c["stem"],
                                               width=2 if is_hl else 1)))
            fig.add_trace(go.Scatter(x=txt[is_hl][0], y=txt[is_hl][1], mode="text",
                                     text=txt[is_hl][2], textposition=anchor, hoverinfo="skip",
                                     showlegend=False,
                                     textfont=dict(family=FONT, size=15 if is_hl else 13,
                                                   color=c["pink_text"] if is_hl else c["ink"])))
    if hl_x:
        for size, alpha in ((40, 0.10), (26, 0.20), (15, 0.36)):
            fig.add_trace(go.Scatter(x=hl_x, y=hl_y, mode="markers", hoverinfo="skip",
                                     showlegend=False,
                                     marker=dict(size=size, color=f"rgba({PINK},{alpha})",
                                                 line=dict(width=0))))

    if cd_x:   # carried to the page, which shows the card next to the pointer
        fig.add_trace(go.Scatter(x=cd_x, y=cd_y, mode="markers", customdata=cd_c,
                                 hoverinfo="none", showlegend=False,
                                 marker=dict(size=22, color="rgba(0,0,0,0)")))
    fig.add_trace(go.Scatter(x=[0], y=[0.10], mode="text", text=[seq_id], hoverinfo="skip",
                             textfont=dict(family=FONT, size=20, color=c["ink"]),
                             showlegend=False))
    fig.add_trace(go.Scatter(x=[0], y=[-0.06], mode="text", hoverinfo="skip", showlegend=False,
                             text=[f"{length:,} bp · {topology} · GC {gc:.1f}%"],
                             textfont=dict(family=FONT, size=13, color=c["muted"])))

    fig.update_layout(
        showlegend=False, dragmode="pan", height=660, annotations=annotations,
        paper_bgcolor=c["paper"], plot_bgcolor=c["paper"],
        margin=dict(l=8, r=8, t=8, b=8),
        xaxis=dict(visible=False, range=[-3.0, 3.0]),
        yaxis=dict(visible=False, range=[-1.55, 1.55], scaleanchor="x", scaleratio=1),
        font=dict(family=FONT, color=c["ink"]),
        hoverlabel=dict(bgcolor=c["hover_bg"], bordercolor=c["hover_border"],
                        font=dict(family=FONT, size=13, color=c["ink"])),
    )
    return fig
