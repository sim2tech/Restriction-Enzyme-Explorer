"""Restriction Enzyme Explorer — Streamlit app.

Run with:  streamlit run app.py
"""
from __future__ import annotations

import html
import math
from pathlib import Path

import Bio
import streamlit as st

from rex.analysis import COMMON_ENZYMES, describe_enzyme, find_sites, list_enzymes
from rex.fasta import gc_percent, parse_fasta
from rex.fragments import compute_fragments
from rex.motif import find_motif
from rex.plotting import MAX_BP_PER_ROW, build_map, plan_rows
from rex.report import build_markdown_report, sites_to_csv

EXAMPLE = Path(__file__).parent / "data" / "example.fasta"
ROWS_PER_PAGE = 6
MAX_LIST = 80
PLOT_CONFIG = {"displaylogo": False, "modeBarButtons": [["toImage"]],
               "toImageButtonOptions": {"format": "png", "filename": "restriction_map",
                                        "scale": 2}}

st.set_page_config(page_title="Restriction Enzyme Explorer", layout="wide")

# Type scale (1.25 ratio): 13 / 16 / 20 / 31 / 49 px, all in Arial Narrow Regular.
STYLE = """
<style>
:root{
  --paper:#F2EDE3; --paper-2:#E8E1D3; --field:#EFE9DC; --ink:#2B2926; --ink-2:#5A554C;
  --rule:#CFC6B4; --accent:#29465B; --accent-2:#4A6880; --code:#E3DBCB;
  --pk1:#FF2E93; --pk2:#FFB3D6; --pk-rgb:255,61,154;
  --s-1:.8125rem; --s0:1rem; --s1:1.25rem; --s2:1.953rem; --s3:3.052rem;
  --font:'Arial Narrow','Arial Narrow Regular','Liberation Sans Narrow','Nimbus Sans Narrow',
         'Helvetica Neue Condensed',Arial,sans-serif;
}
html{font-synthesis:none;}
html, body, p, label, li, h1, h2, h3, h4, dl, dt, dd, td, th, code, button, input, textarea,
[data-testid="stMarkdownContainer"] *, [data-baseweb="select"] *, [data-baseweb="tab"],
[data-baseweb="input"] *{font-family:var(--font) !important; font-weight:400;}
b, strong, h1, h2, h3, h4{font-weight:400 !important;}
html, body, .stApp, [data-testid="stAppViewContainer"], [data-testid="stMain"]{
  background:var(--paper) !important; color:var(--ink); font-size:var(--s0);}
[data-testid="stHeader"]{background:transparent;}
[data-testid="stSidebar"]{background:var(--paper-2) !important; border-right:1px solid var(--rule);}
h1,h2,h3,h4{color:var(--ink); letter-spacing:-.005em;}
[data-testid="stSidebar"] h2{font-size:var(--s1) !important; margin-top:1.4rem;
  padding-bottom:.3rem; border-bottom:1px solid var(--rule);}
.rex-title{font-size:var(--s3); line-height:1.05; letter-spacing:-.015em; margin:.2rem 0 .6rem;}
.rex-dek{font-size:var(--s1); line-height:1.45; color:var(--ink-2); max-width:40rem; margin:0;}
.rex-record{font-size:var(--s0); color:var(--ink-2); margin:0 0 1rem;}
.rex-record b{color:var(--ink);}
.rex-stats{display:flex; flex-wrap:wrap; margin:1.4rem 0 1.2rem; padding:0;
  border-top:2px solid var(--accent); border-bottom:1px solid var(--rule);}
.rex-stat{padding:.8rem 2.2rem .8rem 0; margin-right:2.2rem; border-right:1px solid var(--rule);}
.rex-stat:last-child{border-right:0; margin-right:0;}
.rex-stat dt{font-size:var(--s-1); letter-spacing:.06em; text-transform:uppercase; color:var(--ink-2);}
.rex-stat dd{margin:0; font-size:var(--s2); line-height:1.15;}
.rex-panel-title{font-size:var(--s1) !important; margin:0 0 .2rem !important;}
.rex-note{font-size:var(--s0); color:var(--ink-2); line-height:1.45; margin:0 0 .5rem;}
.rex-note b{color:var(--ink);}
button[data-baseweb="tab"]{font-size:var(--s0);}
button[data-baseweb="tab"][aria-selected="true"]{color:var(--accent);}
[data-baseweb="tab-highlight"]{background:var(--accent) !important;}
[data-baseweb="select"] > div, [data-baseweb="input"] > div,
[data-testid="stFileUploaderDropzone"], [data-baseweb="base-input"]{
  background:var(--field) !important; border-color:var(--rule) !important;}
code{background:var(--code) !important; color:var(--ink) !important;}
.stButton > button, .stDownloadButton > button, [data-testid="stPopover"] button{
  font-size:var(--s0); border-radius:3px; background:transparent; color:var(--ink);
  border:1px solid var(--accent);}
.stButton > button:hover, .stDownloadButton > button:hover, [data-testid="stPopover"] button:hover{
  background:var(--accent); color:var(--paper); border-color:var(--accent);}
[data-testid="stPopoverBody"], div[data-baseweb="popover"] > div{
  background:var(--paper) !important; border:1px solid var(--rule);}
/* Enzyme list: quiet rows, selected one gets a pink gradient outline + glow. */
[class*="st-key-enzpanel"] .stButton > button{
  justify-content:flex-start; text-align:left; border:1px solid var(--rule);
  font-size:var(--s0); padding:.4rem .7rem;}
[class*="st-key-enzpanel"] .stButton > button:hover{
  background:var(--paper-2); color:var(--ink); border-color:var(--accent);}
button[kind="primary"], button[data-testid="stBaseButton-primary"]{
  color:var(--ink) !important; border:2px solid transparent !important;
  background:linear-gradient(var(--paper),var(--paper)) padding-box,
             linear-gradient(135deg,var(--pk1) 0%,var(--pk2) 50%,var(--pk1) 100%) border-box !important;
  box-shadow:0 0 16px rgba(var(--pk-rgb),.45), inset 0 0 12px rgba(var(--pk-rgb),.16);}
.rex-tablewrap{max-height:520px; overflow:auto; border-top:2px solid var(--accent); margin:.4rem 0 1rem;}
.rex-table{width:100%; border-collapse:collapse; font-size:var(--s0); background:transparent;}
.rex-table th{position:sticky; top:0; background:var(--paper) !important; text-align:left;
  font-size:var(--s-1); letter-spacing:.06em; text-transform:uppercase; color:var(--ink-2) !important;
  padding:.55rem .9rem .45rem 0 !important; border:0 !important; border-bottom:1px solid var(--rule) !important;}
.rex-table td{color:var(--ink) !important; padding:.4rem .9rem .4rem 0 !important; background:transparent !important;
  border:0 !important; border-bottom:1px solid var(--rule) !important;}
.rex-fragstrip{display:flex; gap:2px; height:46px; margin:.4rem 0 1.2rem;}
.rex-frag{display:flex; align-items:center; justify-content:center; min-width:3px;
  overflow:hidden; white-space:nowrap; color:var(--paper); font-size:var(--s0);}
</style>
"""
st.markdown(STYLE, unsafe_allow_html=True)

# Dark mode: the whole page is colour-inverted with CSS so every widget follows, and the
# map is inverted back so it keeps its own colours (invert + hue-rotate(180deg) applied
# twice is the identity). The variables below are pre-inverted so the page lands on a
# neutral charcoal: page #17181A, sidebar #1F2124, text #E6E4E0, accent #8FB0C9.
dark = st.session_state.get("theme") == "Dark"
if dark:
    st.markdown("""
<style>
:root{--paper:#E6E7E9; --paper-2:#DCDEE1; --field:#D4D7DA; --ink:#1D1B17; --ink-2:#625F59;
  --rule:#C6C9CD; --accent:#395A73; --accent-2:#567A95; --code:#C6CCD2;
  --pk1:#D9498B; --pk2:#73274A; --pk-rgb:217,73,139;}
html{filter:invert(1) hue-rotate(180deg);}
[data-testid="stPlotlyChart"]{filter:invert(1) hue-rotate(180deg);}
</style>""", unsafe_allow_html=True)


@st.cache_data
def _enzyme_names(commercial_only: bool) -> list[str]:
    return list_enzymes(commercial_only)


@st.cache_data(show_spinner="Scanning sequence…")
def _scan(seq: str, names: tuple, start: int, end: int, circular: bool):
    return find_sites(seq, list(names), start, end, circular)


@st.cache_data
def _describe(name: str):
    return describe_enzyme(name)


def _toggle(name: str) -> None:
    st.session_state["hl"] = None if st.session_state.get("hl") == name else name


def _table(rows: list, max_rows: int = 1000) -> None:
    """Themed HTML table (follows light/dark and the Arial Narrow type)."""
    if not rows:
        return
    cols = list(rows[0])

    def cell(v) -> str:
        return html.escape(f"{v:,}" if isinstance(v, int) and not isinstance(v, bool) else str(v))

    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td>{cell(r[c])}</td>" for c in cols) + "</tr>"
                   for r in rows[:max_rows])
    st.markdown(f'<div class="rex-tablewrap"><table class="rex-table"><thead><tr>{head}'
                f'</tr></thead><tbody>{body}</tbody></table></div>', unsafe_allow_html=True)
    if len(rows) > max_rows:
        st.caption(f"Showing the first {max_rows:,} of {len(rows):,} rows.")


def _clear() -> None:
    st.session_state["hl"] = None


def _enzyme_window(label: str):
    """Open/close window for the enzyme list (popover; expander on old Streamlit)."""
    if hasattr(st, "popover"):
        try:
            return st.popover(label, use_container_width=True)
        except TypeError:
            return st.popover(label)
    return st.expander(label)


def _panel():
    """Scrollable container for the enzyme list (falls back on older Streamlit)."""
    for kwargs in ({"height": 520, "key": "enzpanel"}, {"height": 520}, {}):
        try:
            return st.container(**kwargs)
        except TypeError:
            continue
    return st.container()


def _fragment_strip(fragments) -> str:
    palette = ["#29465B", "#4A6880", "#3A566C", "#5E7A90"]
    total = sum(f.length for f in fragments) or 1
    parts = []
    for i, f in enumerate(fragments):
        label = f"{f.length:,} bp" if f.length / total >= 0.08 else ""
        parts.append(f'<div class="rex-frag" style="flex:{f.length} 1 0;'
                     f'background:{palette[i % len(palette)]}" '
                     f'title="{f.start:,}-{f.end:,} ({f.length:,} bp)">{label}</div>')
    return '<div class="rex-fragstrip">' + "".join(parts) + "</div>"


st.markdown(
    '<h1 class="rex-title">Restriction Enzyme Explorer</h1>'
    '<p class="rex-dek">Load a DNA sequence, choose a region, and see where known '
    'restriction enzymes cut it.</p>', unsafe_allow_html=True)

# ----------------------------------------------------------------- sidebar
with st.sidebar:
    st.header("Sequence")
    upload = st.file_uploader("FASTA file", type=["fa", "fasta", "fna", "txt"])
    if upload is not None:
        text = upload.getvalue().decode("utf-8-sig", errors="replace")
    else:
        text = EXAMPLE.read_text(encoding="utf-8")
        st.caption("No file uploaded. Showing the bundled synthetic example.")

try:
    records = parse_fasta(text)
except ValueError as err:
    st.error(str(err))
    st.stop()

with st.sidebar:
    record = records[0]
    if len(records) > 1:
        record = st.selectbox("Record", records,
                              format_func=lambda r: f"{r.id} ({len(r.sequence):,} bp)")
    seq = record.sequence
    seq_len = len(seq)
    if seq_len < 2:
        st.error("Sequence is too short to analyse.")
        st.stop()

    st.header("Region")
    start, end = st.slider("Region (bp)", 1, seq_len, (1, seq_len))
    whole = (start, end) == (1, seq_len)
    circular = st.checkbox("Treat as circular (whole sequence only)",
                           value=False, disabled=not whole) and whole

    st.header("Enzymes")
    mode = st.radio("Enzyme set", ["Common cloning enzymes",
                                   "All commercial enzymes", "Pick my own"])
    available = _enzyme_names(False)
    if mode == "Common cloning enzymes":
        names = [e for e in COMMON_ENZYMES if e in available]
    elif mode == "All commercial enzymes":
        names = _enzyme_names(True)
    else:
        names = st.multiselect("Enzymes", available, default=["EcoRI", "BamHI"])
    cut_filter = st.selectbox("Show enzymes that cut…",
                              ["any number of times", "at most 3 times",
                               "exactly once (unique cutters)"])
    max_cuts = {"any number of times": None, "at most 3 times": 3,
                "exactly once (unique cutters)": 1}[cut_filter]

    st.header("Display")
    layout = st.radio("Layout", ["Wrapped rows", "Single strip"],
                      help="Wrapped rows: up to 1,200 bases per row, the rest "
                           "continues below. Single strip: one row you scroll along.")
    bp_per_row = st.select_slider(
        "Bases per row" if layout == "Wrapped rows" else "Bases in view",
        options=[300, 600, 900, MAX_BP_PER_ROW], value=MAX_BP_PER_ROW)
    band = st.radio("Band colour", ["Orange", "Green"], horizontal=True)
    st.radio("Theme", ["Light", "Dark"], horizontal=True, key="theme")

if not names:
    st.warning("Select at least one enzyme.")
    st.stop()

# ---------------------------------------------------------------- analysis
region_seq = seq[start - 1:end]
all_sites = _scan(seq, tuple(names), start, end, circular)
counts: dict[str, int] = {}
site_of: dict[str, str] = {}
for s in all_sites:
    counts[s.enzyme] = counts.get(s.enzyme, 0) + 1
    site_of.setdefault(s.enzyme, s.site)
sites = [s for s in all_sites if max_cuts is None or counts[s.enzyme] <= max_cuts]
cutters = sorted({s.enzyme for s in sites})

hl = st.session_state.get("hl")
if hl not in cutters:
    hl = None

stats = [("Sequence", f"{seq_len:,} bp"), ("Region", f"{end - start + 1:,} bp"),
         ("GC (region)", f"{gc_percent(region_seq):.1f}%"),
         ("Enzymes cutting", f"{len(cutters)}"), ("Cut sites", f"{len(sites)}")]
st.markdown(
    '<dl class="rex-stats">' + "".join(
        f'<div class="rex-stat"><dt>{k}</dt><dd>{v}</dd></div>' for k, v in stats)
    + "</dl>"
    f'<p class="rex-record"><b>{html.escape(record.id)}</b> '
    f'{html.escape(record.description)}</p>', unsafe_allow_html=True)

tab_map, tab_sites, tab_digest, tab_motif, tab_report = st.tabs(
    ["Map", "Sites", "Digest", "Motif search", "Report"])

# --------------------------------------------------------------------- map
with tab_map:
    head_l, head_r = st.columns([3, 1.3], gap="medium")
    with head_l:
        if hl:
            st.markdown(f'<p class="rex-note">Highlighting <b>{html.escape(hl)}</b>.</p>',
                        unsafe_allow_html=True)
            st.button("Clear highlight", key="clear_hl", on_click=_clear)
        else:
            st.markdown('<p class="rex-note">Open the enzyme window on the right to '
                        'light up one enzyme\'s cuts.</p>', unsafe_allow_html=True)
    with head_r:
        with _enzyme_window(f"Enzymes in view ({len(cutters)})"):
            st.markdown('<p class="rex-note">Click one to highlight its cuts on the map. '
                        'Click again to clear.</p>', unsafe_allow_html=True)
            if cutters:
                with _panel():
                    for name in cutters[:MAX_LIST]:
                        st.button(f"{name} · {site_of[name]} · {counts[name]}×",
                                  key=f"enz_{name}", on_click=_toggle, args=(name,),
                                  type="primary" if name == hl else "secondary",
                                  use_container_width=True)
                if len(cutters) > MAX_LIST:
                    st.caption(f"Showing the first {MAX_LIST} of {len(cutters)}. "
                               "Use the cut filter to narrow the list.")
            else:
                st.caption("No enzymes cut in this region.")
            if hl:
                info = _describe(hl)
                where = [s.position for s in sites if s.enzyme == hl]
                shown = ", ".join(f"{p:,}" for p in where[:12]) + (" …" if len(where) > 12 else "")
                st.markdown(f'<p class="rex-note"><b>{html.escape(hl)}</b><br>'
                            f'<code>{html.escape(info.cut_pattern)}</code> · {info.overhang}<br>'
                            f'cuts before bp {shown}</p>', unsafe_allow_html=True)

    rows_all = plan_rows((start, end), bp_per_row)
    if layout == "Wrapped rows":
        pages = math.ceil(len(rows_all) / ROWS_PER_PAGE)
        page = 1
        if pages > 1:
            page = int(st.number_input(f"Rows page (1 to {pages})", 1, pages, 1))
        visible = rows_all[(page - 1) * ROWS_PER_PAGE: page * ROWS_PER_PAGE]
    else:
        visible = [rows_all[0]]
        if end - start + 1 > bp_per_row:
            ws = st.slider("Scroll along the sequence (window starts at bp)",
                           start, end - bp_per_row + 1, start,
                           step=max(1, bp_per_row // 20))
            visible = [(ws, ws + bp_per_row - 1)]

    fig = build_map(visible, sites, bp_per_row=bp_per_row, highlight=hl,
                    band=band.lower(), theme="dark" if dark else "light")
    st.plotly_chart(fig, use_container_width=True, config=PLOT_CONFIG, theme=None)
    st.markdown('<p class="rex-note">Each band is the sequence; a dark notch marks a cut. '
                'Hover a marker for its recognition site. The camera icon saves a PNG.</p>',
                unsafe_allow_html=True)
    if not sites:
        st.caption("No cut sites for these enzymes in this region.")

# ------------------------------------------------------------------- sites
with tab_sites:
    if sites:
        st.markdown(f'<p class="rex-dek">{len(sites)} cut site(s) from {len(cutters)} '
                    f'enzyme(s) in bp {start:,} to {end:,}.</p>', unsafe_allow_html=True)
        rows = []
        for s in sites:
            info = _describe(s.enzyme)
            rows.append({"Enzyme": s.enzyme, "Cut position (bp)": s.position,
                         "Recognition site": s.site, "Cut pattern": s.cut_pattern,
                         "Ends": info.overhang, "Cuts in region": counts[s.enzyme]})
        _table(rows)
        st.download_button("Download sites (CSV)", sites_to_csv(sites),
                           file_name=f"{record.id}_sites.csv", mime="text/csv")
    else:
        st.info("Nothing to show.")

# ------------------------------------------------------------------ digest
digest_enzymes: list[str] = []
fragments = []
with tab_digest:
    if cutters:
        digest_enzymes = st.multiselect("Enzymes in the digest", cutters,
                                        default=cutters[:1])
        cuts = [s.position for s in sites if s.enzyme in digest_enzymes]
        fragments = compute_fragments(cuts, start, end, circular)
        st.markdown(_fragment_strip(fragments), unsafe_allow_html=True)
        _table([{"Fragment": i, "Start": f.start, "End": f.end,
                       "Length (bp)": f.length}
                      for i, f in enumerate(fragments, 1)])
    else:
        st.info("No cutting enzymes to digest with.")

# ------------------------------------------------------------------- motif
with tab_motif:
    m_left, m_right = st.columns([1, 2], gap="large")
    with m_left:
        motif = st.text_input("Motif (IUPAC codes allowed)", placeholder="GAATTC or GRCGYC")
        both = st.checkbox("Search both strands", value=True)
    with m_right:
        if motif:
            try:
                hits = find_motif(seq, motif, start, end, both)
            except ValueError as err:
                st.error(str(err))
            else:
                st.markdown(f'<p class="rex-dek">{len(hits)} match(es) in bp '
                            f'{start:,} to {end:,}.</p>', unsafe_allow_html=True)
                if hits:
                    _table([{"Strand": h.strand, "Start": h.start, "End": h.end,
                                   "Matched": h.matched} for h in hits[:1000]])
        else:
            st.caption("Type a motif to search the selected region.")

# ------------------------------------------------------------------ report
with tab_report:
    infos = {n: _describe(n) for n in cutters}
    report = build_markdown_report(
        record_id=record.id, description=record.description, seq_len=seq_len,
        region=(start, end), gc_percent=gc_percent(region_seq), circular=circular,
        enzyme_set_label=mode, enzymes_scanned=len(names), sites=sites,
        infos=infos, digest_enzymes=digest_enzymes, fragments=fragments,
        biopython_version=Bio.__version__)
    st.download_button("Download report (Markdown)", report,
                       file_name=f"{record.id}_restriction_report.md",
                       mime="text/markdown")
    st.markdown(report)
