"""Restriction map (Milestone 3), rendered with Plotly.

The sequence is drawn as a flat orange or green band on the page background.
Long sequences wrap into rows of at most MAX_BP_PER_ROW bases, or a single row
can be scrolled along the sequence. A selected enzyme gets a pink gradient
glow around its cut sites.
"""
from __future__ import annotations

import plotly.graph_objects as go

from .models import CutSite

MAX_BP_PER_ROW = 1200

# Colours seen on screen. The dark values are what you actually see in dark mode:
# the page is colour-inverted by CSS and the map is inverted back, so the map's
# own colours are used as written here.
THEMES = {
    "light": dict(paper="#F2EDE3", ink="#2B2926", ink_dim="rgba(43,41,38,0.38)",
                  stem="rgba(43,41,38,0.45)", muted="#7A7164", pink_text="#D81B7A",
                  hover_bg="#E8E1D3", hover_border="#CFC6B4", cut="#C0504D"),
    "dark": dict(paper="#17181A", ink="#E6E4E0", ink_dim="rgba(230,228,224,0.38)",
                 stem="rgba(230,228,224,0.45)", muted="#8B8883", pink_text="#FF6FB1",
                 hover_bg="#26292C", hover_border="#34373B", cut="#E36A63"),
}
PINK = "255,61,154"
PINK_CORE = "#FFB3D6"
FONT = "Arial Narrow, Liberation Sans Narrow, Helvetica Neue Condensed, Arial, sans-serif"

BANDS = {"orange": "#FF8A1F", "green": "#34D17A"}
# (half-width in bp at 1,200 bp/row, alpha) layers for the pink highlight.
PINK_LAYERS = [(14, 0.05), (9.5, 0.07), (6.5, 0.10), (4, 0.14), (2, 0.20)]

ROW_H = 3.9          # vertical distance between rows (plot units)
PX_PER_UNIT = 52     # pixels per plot unit


def plan_rows(region: tuple[int, int], bp_per_row: int = MAX_BP_PER_ROW
              ) -> list[tuple[int, int]]:
    """Split region (1-based, inclusive) into consecutive rows of <= bp_per_row."""
    start, end = region
    bp_per_row = max(1, min(bp_per_row, MAX_BP_PER_ROW))
    return [(s, min(s + bp_per_row - 1, end)) for s in range(start, end + 1, bp_per_row)]


def _tick_step(width: int) -> int:
    for step in (10, 20, 25, 50, 100, 200, 250, 500, 1000):
        if width / step <= 8:
            return step
    return 1000


def _pts():
    return dict(x=[], y=[], text=[], custom=[])


def build_map(rows: list[tuple[int, int]], sites: list[CutSite], *,
              bp_per_row: int = MAX_BP_PER_ROW, highlight: str | None = None,
              band: str = "orange", theme: str = "light",
              labelled: set | None = None, cards: dict | None = None) -> go.Figure:
    """Draw the given rows (each a (start, end) bp range) with their cut sites.

    ``labelled`` limits which enzymes get a name and stem; every other cut is still
    drawn as a notch on the band (None labels all of them).
    """
    c = THEMES.get(theme, THEMES["light"])
    colour = BANDS.get(band, BANDS["orange"])
    width = max(1, min(bp_per_row, MAX_BP_PER_ROW))
    scale = width / MAX_BP_PER_ROW
    step = _tick_step(width)
    dim = highlight is not None

    shapes: list[dict] = []
    stem_x, stem_y = [], []
    hl_stem_x, hl_stem_y = [], []
    notch_x, notch_y = [], []
    tick_x, tick_y = [], []
    tick_lx, tick_ly, tick_lt = [], [], []
    title_y, title_t = [], []
    reg, hl = _pts(), _pts()

    for i, (rs, re_) in enumerate(rows):
        y0 = -i * ROW_H
        length = re_ - rs + 1

        # The sequence itself: one flat band.
        shapes.append(dict(type="rect", x0=0, x1=length, y0=y0 - 0.16, y1=y0 + 0.16,
                           fillcolor=colour, line=dict(width=0), layer="below"))

        # Coordinate ticks under the band.
        for t in range(-(-rs // step) * step, re_ + 1, step):
            x = t - rs + 0.5
            tick_x += [x, x, None]
            tick_y += [y0 - 0.28, y0 - 0.40, None]
            tick_lx.append(x)
            tick_ly.append(y0 - 0.62)
            tick_lt.append(f"{t:,}")
        title_y.append(y0 - 1.05)
        title_t.append(f"bp {rs:,}\u2013{re_:,}")

        # Cut sites that belong to this row.
        row_sites = sorted((s for s in sites if rs <= s.position <= re_),
                           key=lambda s: (s.position, s.enzyme))
        k = -1
        for s in row_sites:
            x = s.position - rs
            notch_x += [x, x, None]
            notch_y += [y0 - 0.22, y0 + 0.22, None]
            if labelled is not None and s.enzyme not in labelled and s.enzyme != highlight:
                continue
            k += 1
            ly = y0 + 0.62 + 0.5 * (k % 4)
            custom = [s.enzyme, s.site, s.cut_pattern, s.position,
                      (cards or {}).get((s.enzyme, s.position), "")]
            if s.enzyme == highlight:
                hl_stem_x += [x, x, None]
                hl_stem_y += [y0 + 0.22, ly, None]
                bucket = hl
                for hw, alpha in PINK_LAYERS:  # soft pink column behind the cut
                    h = hw * scale
                    shapes.append(dict(type="rect", x0=x - h, x1=x + h,
                                       y0=y0 - 0.75, y1=ly + 0.35,
                                       fillcolor=f"rgba({PINK},{alpha})",
                                       line=dict(width=0), layer="below"))
            else:
                stem_x += [x, x, None]
                stem_y += [y0 + 0.22, ly - 0.08, None]
                bucket = reg
            bucket["x"].append(x)
            bucket["y"].append(ly)
            bucket["text"].append(s.enzyme)
            bucket["custom"].append(custom)

    fig = go.Figure()
    hover = ("<b>%{customdata[0]}</b><br>cuts before bp %{customdata[3]:,}"
             "<br>site %{customdata[1]}<br>%{customdata[2]}<extra></extra>")
    # With cards the page shows a rich card on hover, so Plotly's own tooltip is switched off.
    hov = {"hoverinfo": "none"} if cards else {"hovertemplate": hover}

    fig.add_trace(go.Scatter(x=tick_x, y=tick_y, mode="lines", hoverinfo="skip",
                             line=dict(color=c["muted"], width=1), showlegend=False))
    fig.add_trace(go.Scatter(x=tick_lx, y=tick_ly, mode="text", text=tick_lt,
                             textfont=dict(family=FONT, size=11, color=c["muted"]),
                             hoverinfo="skip", showlegend=False))
    fig.add_trace(go.Scatter(x=[0] * len(title_t), y=title_y, mode="text", text=title_t,
                             textposition="middle right",
                             textfont=dict(family=FONT, size=12, color=c["muted"]),
                             hoverinfo="skip", showlegend=False))
    # Dark notch across the band at every cut.
    fig.add_trace(go.Scatter(x=notch_x, y=notch_y, mode="lines", hoverinfo="skip",
                             line=dict(color=c["ink"], width=2.5), showlegend=False))
    fig.add_trace(go.Scatter(x=stem_x, y=stem_y, mode="lines", hoverinfo="skip",
                             line=dict(color=c["stem"], width=1), showlegend=False))
    fig.add_trace(go.Scatter(
        x=reg["x"], y=reg["y"], mode="markers+text", text=reg["text"],
        textposition="top center", customdata=reg["custom"], **hov,
        textfont=dict(family=FONT, size=13, color=c["ink_dim"] if dim else c["ink"]),
        marker=dict(size=7, color=c["ink_dim"] if dim else c["ink"], line=dict(width=0)),
        showlegend=False))

    if hl["x"]:
        fig.add_trace(go.Scatter(x=hl_stem_x, y=hl_stem_y, mode="lines", hoverinfo="skip",
                                 line=dict(color=f"rgba({PINK},0.25)", width=10),
                                 showlegend=False))
        fig.add_trace(go.Scatter(x=hl_stem_x, y=hl_stem_y, mode="lines", hoverinfo="skip",
                                 line=dict(color=c["pink_text"], width=2), showlegend=False))
        for size, alpha in ((44, 0.10), (30, 0.20), (19, 0.36)):  # radial halo
            fig.add_trace(go.Scatter(
                x=hl["x"], y=hl["y"], mode="markers", hoverinfo="skip", showlegend=False,
                marker=dict(size=size, color=f"rgba({PINK},{alpha})", line=dict(width=0))))
        fig.add_trace(go.Scatter(
            x=hl["x"], y=hl["y"], mode="markers+text", text=hl["text"],
            textposition="top center", customdata=hl["custom"], **hov,
            textfont=dict(family=FONT, size=15, color=c["pink_text"]),
            marker=dict(size=9, color=PINK_CORE, line=dict(color=c["pink_text"], width=2)),
            showlegend=False))

    n = max(1, len(rows))
    ymax = 2.9
    ymin = -(n - 1) * ROW_H - 1.6
    fig.update_layout(
        shapes=shapes, showlegend=False, dragmode="pan",
        paper_bgcolor=c["paper"], plot_bgcolor=c["paper"],
        margin=dict(l=8, r=8, t=8, b=8),
        height=int((ymax - ymin) * PX_PER_UNIT),
        xaxis=dict(visible=False, range=[-0.015 * width, 1.015 * width]),   # zoom/pan along x
        yaxis=dict(visible=False, range=[ymin, ymax], fixedrange=True),
        font=dict(family=FONT, color=c["ink"]),
        hoverlabel=dict(bgcolor=c["hover_bg"], bordercolor=c["hover_border"],
                        font=dict(family=FONT, size=13, color=c["ink"])),
    )
    return fig
