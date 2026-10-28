"""Sequence view (pure Python -> HTML): both strands, cut sites, recognition sites, features.

Each base sits in a fixed-width cell so the two strands, cut marks and feature bars
line up even with a proportional font.
"""
from __future__ import annotations

import html

from .geometry import feature_color, text_on
from .motif import COMPLEMENT, find_motif

CELL_EM = 0.72


def complement(seq: str) -> str:
    return "".join(COMPLEMENT.get(c, c) for c in seq)


def _cells(seq: str, first: int, classes: dict[int, str]) -> str:
    out = []
    for i, base in enumerate(seq):
        cls = classes.get(first + i, "")
        out.append(f'<i class="b {cls}">{base}</i>' if cls else f'<i class="b">{base}</i>')
    return "".join(out)


def _feature_rows(features, rs: int, re_: int):
    """Greedy lane packing of feature segments that overlap [rs, re_]."""
    segs = []
    for f in features:
        n = len(f.segments)
        for k, (s, e) in enumerate(f.segments):
            if e < rs or s > re_:
                continue
            arrow_here = (f.strand > 0 and k == n - 1) or (f.strand < 0 and k == 0)
            segs.append((max(s, rs), min(e, re_), f, arrow_here))
    segs.sort(key=lambda t: t[0])
    lanes: list[int] = []   # last occupied position per lane
    placed = []
    for s, e, f, arrow in segs:
        for k, last in enumerate(lanes):
            if s > last:
                lanes[k] = e
                placed.append((k, s, e, f, arrow))
                break
        else:
            lanes.append(e)
            placed.append((len(lanes) - 1, s, e, f, arrow))
    return placed, len(lanes)


def build_sequence_html(seq: str, ws: int, we: int, sites, features=(), per_row: int = 60,
                        max_enzymes_highlighted: int = 40) -> str:
    """HTML for bases ws..we (1-based, inclusive) with both strands and annotations."""
    top_cut = {s.position: s for s in sites if ws <= s.position <= we}
    bot_cut = {s.bottom_position: s for s in sites
               if s.bottom_position is not None and ws <= s.bottom_position <= we}

    # Recognition sites for enzymes that cut near this window.
    rec: dict[int, str] = {}
    near = {(s.enzyme, s.site) for s in sites if ws - 40 <= s.position <= we + 40}
    for _, site in sorted(near)[:max_enzymes_highlighted]:
        for hit in find_motif(seq, site, ws, we, both_strands=True):
            for pos in range(hit.start, hit.end + 1):
                rec[pos] = "rs"

    blocks = []
    for rs in range(ws, we + 1, per_row):
        re_ = min(rs + per_row - 1, we)
        chunk = seq[rs - 1:re_]
        top_cls = {p: ("rs " if p in rec else "") + ("ct" if p in top_cut else "")
                   for p in range(rs, re_ + 1) if p in rec or p in top_cut}
        bot_cls = {p: "cb" for p in range(rs, re_ + 1) if p in bot_cut}

        labels = []
        in_block = sorted((p, s) for p, s in top_cut.items() if rs <= p <= re_)
        for k, (p, s) in enumerate(in_block):
            labels.append(f'<span class="sv-l" style="left:{(p - rs) * CELL_EM:.2f}em;'
                          f'top:{(k % 2) * 1.15:.2f}em">{html.escape(s.enzyme)}</span>')
        lab_h = 2.4 if len(in_block) > 1 else (1.25 if in_block else 0.25)

        placed, nlanes = _feature_rows(features, rs, re_)
        bars = []
        for lane, s, e, f, arrow in placed:
            colour = feature_color(f.type)
            width_em = (e - s + 1) * CELL_EM
            glyph = ("\u25B6" if f.strand > 0 else "\u25C0") if (arrow and f.strand) else ""
            label = f"{f.name} {glyph}" if f.strand > 0 else f"{glyph} {f.name}"
            title = f"{f.name} ({f.type}) {f.segments[0][0]:,}-{f.segments[-1][1]:,}"
            bars.append(
                f'<span class="sv-bar" title="{html.escape(title)}" style="left:'
                f'{(s - rs) * CELL_EM:.2f}em;top:{lane * 1.4:.2f}em;width:{width_em:.2f}em;'
                f'background:{colour};color:{text_on(colour)}">{html.escape(label.strip())}</span>')

        blocks.append(
            f'<div class="sv-block"><div class="sv-labs" style="height:{lab_h}em">'
            f'{"".join(labels)}</div>'
            f'<div class="sv-row"><span class="sv-n">{rs:,}</span>'
            f'<span class="sv-s">{_cells(chunk, rs, top_cls)}</span>'
            f'<span class="sv-n e">{re_:,}</span></div>'
            f'<div class="sv-row"><span class="sv-n"></span>'
            f'<span class="sv-s">{_cells(complement(chunk), rs, bot_cls)}</span></div>'
            f'<div class="sv-feats" style="height:{max(nlanes, 0) * 1.4}em">'
            f'{"".join(bars)}</div></div>')
    return '<div class="sv">' + "".join(blocks) + "</div>"


def snippet_html(seq: str, site, flank: int = 8) -> str:
    """Tiny two-strand view around one cut, for the enzyme detail card."""
    p = site.position
    ws, we = max(1, p - flank), min(len(seq), p + flank - 1)
    chunk = seq[ws - 1:we]
    top = {}
    for hit in find_motif(seq, site.site, ws, we, both_strands=True):
        for pos in range(hit.start, hit.end + 1):
            top[pos] = "rs"
    top[p] = (top.get(p, "") + " ct").strip()
    bot = {}
    if site.bottom_position is not None and ws <= site.bottom_position <= we:
        bot[site.bottom_position] = "cb"
    return ('<div class="sv rex-snip"><div class="sv-row"><span class="sv-n">5\u2032</span>'
            f'<span class="sv-s">{_cells(chunk, ws, top)}</span><span class="sv-n e">3\u2032</span></div>'
            '<div class="sv-row"><span class="sv-n">3\u2032</span>'
            f'<span class="sv-s">{_cells(complement(chunk), ws, bot)}</span>'
            '<span class="sv-n e">5\u2032</span></div></div>')


def enzyme_card_html(seq: str, name: str, info, enzyme_sites, total_cuts: int) -> str:
    """Flat detail panel for one enzyme: stats, both-strand snippet, pattern, isoschizomers."""
    first = enzyme_sites[0]
    unique = total_cuts == 1
    after = ", ".join(f"{s.position - 1:,}" for s in enzyme_sites[:3]) \
        + (" \u2026" if len(enzyme_sites) > 3 else "")
    same = ", ".join(info.isoschizomers[:6]) + (" \u2026" if len(info.isoschizomers) > 6 else "")
    parts = [
        '<div class="rex-card">',
        f'<div class="rc-head"><span class="rc-name">{html.escape(name)}</span>'
        f'<span class="rc-tag">{"unique cutter" if unique else f"{total_cuts} sites"}</span></div>',
        '<div class="rc-stats">'
        f'<div><b>{after}</b><span>cuts after</span></div>'
        f'<div><b>{html.escape(info.overhang)}</b><span>ends</span></div>'
        f'<div><b>{total_cuts}</b><span>{"site \u00B7 unique" if unique else "sites"}</span></div></div>',
        snippet_html(seq, first),
        f'<div class="rc-line">top cut {first.position - 1:,} | {first.position:,}'
        f' \u00B7 <code>{html.escape(info.cut_pattern)}</code></div>',
    ]
    if same:
        parts.append(f'<div class="rc-line">same site: {html.escape(same)}</div>')
    parts.append("</div>")
    return "".join(parts)
