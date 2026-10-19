"""Interactive linear restriction map (Milestone 3) built with Plotly."""
from __future__ import annotations

import plotly.graph_objects as go
from plotly.colors import qualitative

from .models import CutSite


def build_map(seq_id: str, seq_len: int, region: tuple[int, int],
              sites: list[CutSite]) -> go.Figure:
    start, end = region
    fig = go.Figure()

    # Full-sequence backbone and highlighted analysis region.
    fig.add_shape(type="rect", x0=1, x1=seq_len, y0=-0.12, y1=0.12,
                  fillcolor="#cbd5e1", line_width=0, layer="below")
    fig.add_shape(type="rect", x0=start, x1=end, y0=-0.2, y1=0.2,
                  fillcolor="rgba(37,99,235,0.35)",
                  line=dict(color="#2563eb", width=1), layer="below")

    # Stagger labels over four heights so neighbouring names rarely overlap.
    ordered = sorted(sites, key=lambda s: s.position)
    level_of = {id(s): 1 + (i % 4) for i, s in enumerate(ordered)}

    names = sorted({s.enzyme for s in sites})
    colour = {n: qualitative.D3[i % len(qualitative.D3)] for i, n in enumerate(names)}

    for name in names:
        group = [s for s in ordered if s.enzyme == name]
        stem_x, stem_y = [], []
        for s in group:
            stem_x += [s.position, s.position, None]
            stem_y += [0, level_of[id(s)] - 0.15, None]
        fig.add_trace(go.Scatter(x=stem_x, y=stem_y, mode="lines",
                                 line=dict(color=colour[name], width=1),
                                 hoverinfo="skip", showlegend=False,
                                 legendgroup=name))
        fig.add_trace(go.Scatter(
            x=[s.position for s in group],
            y=[level_of[id(s)] for s in group],
            mode="markers+text",
            text=[name] * len(group),
            textposition="top center",
            marker=dict(size=9, color=colour[name]),
            name=f"{name} ({len(group)})",
            legendgroup=name,
            customdata=[[s.site, s.cut_pattern] for s in group],
            hovertemplate=(f"<b>{name}</b><br>cut before bp %{{x}}"
                           "<br>site %{customdata[0]}"
                           "<br>pattern %{customdata[1]}<extra></extra>"),
        ))

    fig.update_layout(
        title=f"{seq_id} — {len(sites)} cut site(s) in bp {start:,}–{end:,}",
        height=420,
        margin=dict(l=20, r=20, t=60, b=20),
        yaxis=dict(visible=False, range=[-1, 5.6]),
        xaxis=dict(title="Position (bp)", range=[1, seq_len + 1],
                   rangeslider=dict(visible=True)),
        legend=dict(orientation="h", yanchor="bottom", y=-0.6),
    )
    return fig
