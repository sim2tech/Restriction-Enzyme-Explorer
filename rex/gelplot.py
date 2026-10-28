"""Virtual agarose gel figure, laid out like NEBcutter's gel view (Plotly)."""
from __future__ import annotations

import plotly.graph_objects as go

from .gel import axis_ticks, band_thickness_mm, migration_mm
from .plotting import FONT

PALETTES = {
    "Black on paper": dict(bg="#F2EDE3", band="#1F1D1A", text="#2B2926", axis="#2B2926",
                           well="#B9B09E"),
    "White on black": dict(bg="#141414", band="#EDE9DF", text="#D8D6D2", axis="#BDB8AD",
                           well="#3A3A3A"),
    "Orange on black": dict(bg="#141210", band="#FF8A1F", text="#D8D6D2", axis="#BDB8AD",
                            well="#3A3A3A"),
}
LANE_W = 11.0
LANE_STEP = 15.0
AXIS_X = 0.0


def _rgba(hex_colour: str, alpha: float) -> str:
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (1, 3, 5))
    return f"rgba({r},{g},{b},{alpha:.2f})"


def build_gel(lanes: list[tuple[str, list[int]]], percent: float = 1.0,
              palette: str = "Black on paper", length_mm: float = 100.0,
              equal_mass_lanes: tuple = ("Marker",)) -> go.Figure:
    """lanes = [(title, [fragment sizes in bp]), ...]; lane 0 is usually the marker."""
    pal = PALETTES.get(palette, PALETTES["Black on paper"])
    fig = go.Figure()
    shapes = []

    # Size axis with major (labelled) and minor ticks, like the NEBcutter gel.
    shapes.append(dict(type="line", x0=AXIS_X, x1=AXIS_X, y0=0, y1=length_mm,
                       line=dict(color=pal["axis"], width=1.5)))
    tx, ty, lx, ly, lt = [], [], [], [], []
    for bp, y, major in axis_ticks(percent, length_mm):
        if not 0 < y < length_mm:
            continue
        size = 2.4 if major else 1.3
        tx += [AXIS_X - size, AXIS_X, None]
        ty += [y, y, None]
        if major:
            lx.append(AXIS_X - 3.0); ly.append(y); lt.append(f"{bp:,}")
    fig.add_trace(go.Scatter(x=tx, y=ty, mode="lines", hoverinfo="skip", showlegend=False,
                             line=dict(color=pal["axis"], width=1)))
    fig.add_trace(go.Scatter(x=lx, y=ly, mode="text", text=lt, textposition="middle left",
                             hoverinfo="skip", showlegend=False,
                             textfont=dict(family=FONT, size=12, color=pal["text"])))

    hx, hy, ht = [], [], []
    titles_x, titles_t = [], []
    for i, (title, sizes) in enumerate(lanes):
        cx = 10.0 + LANE_STEP * i + LANE_W / 2
        titles_x.append(cx); titles_t.append(title)
        shapes.append(dict(type="rect", x0=cx - LANE_W / 2, x1=cx + LANE_W / 2, y0=-3, y1=0,
                           line=dict(color=pal["well"], width=1.5), fillcolor="rgba(0,0,0,0)"))
        total = float(sum(sizes)) or 1.0
        for bp in sizes:
            y = migration_mm(bp, percent, length_mm)
            if title in equal_mass_lanes:
                th, alpha = 0.7, 0.85
            else:
                th = band_thickness_mm(bp, total)
                alpha = 0.55 + 0.45 * min(1.0, (bp / total) * 2.5)
            shapes.append(dict(type="rect", x0=cx - LANE_W / 2 + 0.6, x1=cx + LANE_W / 2 - 0.6,
                               y0=y - th / 2, y1=y + th / 2, line=dict(width=0),
                               fillcolor=_rgba(pal["band"], alpha)))
            hx.append(cx); hy.append(y); ht.append(f"<b>{title}</b><br>{bp:,} bp")
    fig.add_trace(go.Scatter(x=hx, y=hy, mode="markers", text=ht, hoverinfo="text",
                             marker=dict(size=18, color="rgba(0,0,0,0)"), showlegend=False))
    fig.add_trace(go.Scatter(x=titles_x, y=[-6.5] * len(titles_x), mode="text", text=titles_t,
                             hoverinfo="skip", showlegend=False,
                             textfont=dict(family=FONT, size=13, color=pal["text"])))

    width = 10.0 + LANE_STEP * max(1, len(lanes)) + 2
    fig.update_layout(
        shapes=shapes, showlegend=False, dragmode=False, height=560,
        paper_bgcolor=pal["bg"], plot_bgcolor=pal["bg"], margin=dict(l=8, r=8, t=8, b=8),
        xaxis=dict(visible=False, range=[-12, width], fixedrange=True),
        yaxis=dict(visible=False, range=[length_mm + 6, -10], fixedrange=True),
        font=dict(family=FONT, color=pal["text"]), hovermode="closest",
        hoverlabel=dict(font=dict(family=FONT, size=13)),
    )
    return fig
