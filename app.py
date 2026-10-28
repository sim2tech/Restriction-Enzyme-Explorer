"""Restriction Enzyme Explorer — Streamlit app.

Run with:  streamlit run app.py
"""
from __future__ import annotations

import html
import json
import math
from pathlib import Path

import Bio
import streamlit as st

from rex.analysis import COMMON_ENZYMES, describe_enzyme, find_sites, list_enzymes
from rex.fasta import gc_percent
from rex.geometry import HIDDEN_FEATURE_TYPES, MAX_LANES, split_lanes
from rex.genbank import parse_any
from rex.circular import build_circular_map
from rex.fragments import digest_fragments, format_coordinates, fragment_sequence
from rex.gel import BAND_THRESHOLD_BP, GEL_RANGES, LADDERS, lane_sizes
from rex.gelplot import PALETTES, build_gel
from rex.motif import find_motif
from rex.plotting import MAX_BP_PER_ROW, build_map, plan_rows
from rex.report import build_markdown_report, sites_to_csv
from rex.seqview import build_sequence_html, enzyme_card_html

EXAMPLE = Path(__file__).parent / "data" / "pBR322.gb"
ROWS_PER_PAGE = 6
MAX_LIST = 80
INTERACTION_HTML = (
    "<script>(function () { try { var w = window.parent; "
    "if (w.__rexInstalled) return; var s = w.document.createElement('script'); "
    "s.textContent = " + json.dumps("\n(function () {\n  if (window.__rexInstalled) return;\n  window.__rexInstalled = true;\n  var tip = document.createElement('div');\n  tip.id = 'rex-tip';\n  document.body.appendChild(tip);\n  function place(e) {\n    var w = tip.offsetWidth, h = tip.offsetHeight;\n    var x = e.clientX + 16, y = e.clientY + 16;\n    if (x + w > window.innerWidth - 8) x = e.clientX - w - 16;\n    if (y + h > window.innerHeight - 8) y = window.innerHeight - h - 8;\n    tip.style.left = Math.max(8, x) + 'px';\n    tip.style.top = Math.max(8, y) + 'px';\n  }\n  function hook(gd) {\n    if (gd.__rexWheel !== true) {\n      gd.__rexWheel = true;\n      // Plain scrolling must move the page; only pinch (ctrl + wheel) may zoom the map.\n      gd.addEventListener('wheel', function (e) {\n        if (!(e.ctrlKey || e.metaKey)) e.stopImmediatePropagation();\n      }, {capture: true, passive: true});\n      gd.addEventListener('mousemove', function (e) {\n        if (tip.style.display === 'block') place(e);\n      });\n    }\n    if (gd.__rexHover !== true && typeof gd.on === 'function') {\n      gd.__rexHover = true;\n      gd.on('plotly_hover', function (ev) {\n        var pt = ev.points && ev.points[0];\n        if (!pt) return;\n        var cd = pt.customdata;\n        if (Array.isArray(cd)) cd = cd[cd.length - 1];\n        if (typeof cd === 'string' && cd.indexOf('rex-card') >= 0) {\n          tip.innerHTML = cd;\n          tip.style.display = 'block';\n          if (ev.event) place(ev.event);\n        }\n      });\n      gd.on('plotly_unhover', function () { tip.style.display = 'none'; });\n    }\n  }\n  function scan() { document.querySelectorAll('.js-plotly-plot').forEach(hook); }\n  scan();\n  new MutationObserver(scan).observe(document.body, {childList: true, subtree: true});\n})();\n") + "; "
    "w.document.head.appendChild(s); } catch (e) {} })();</script>"
)
MAX_NAMED = 40          # most enzymes named on a map at once
SEQ_LINES = 10          # lines per page in the sequence view
MAP_CONFIG = {"displaylogo": False, "scrollZoom": True,
              "modeBarButtons": [["zoomIn2d", "zoomOut2d", "resetScale2d", "toImage"]],
              "toImageButtonOptions": {"format": "png", "filename": "restriction_map", "scale": 2}}
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
.rex-card{border-top:2px solid var(--accent); border-bottom:1px solid var(--rule);
  padding:.7rem 0 .6rem; margin:.2rem 0 1rem; max-width:36rem;}
.rc-head{display:flex; justify-content:space-between; align-items:baseline; margin-bottom:.5rem;}
.rc-name{font-size:var(--s1);}
.rc-tag{font-size:var(--s-1); letter-spacing:.06em; text-transform:uppercase; color:var(--ink-2);}
.rc-stats{display:flex; margin-bottom:.6rem;}
.rc-stats div{padding:0 1.4rem 0 0; margin-right:1.4rem; border-right:1px solid var(--rule);}
.rc-stats div:last-child{border-right:0;}
.rc-stats b{display:block; font-size:var(--s1);}
.rc-stats span{font-size:var(--s-1); color:var(--ink-2);}
.rc-line{font-size:var(--s-1); color:var(--ink-2); margin-top:.35rem;}
.rex-snip{margin:.2rem 0 .3rem;}
#rex-tip{position:fixed; z-index:1000000; display:none; pointer-events:none; max-width:27rem;
  background:var(--paper); color:var(--ink); border:1px solid var(--rule); padding:.1rem 1rem .3rem;
  box-shadow:0 8px 28px rgba(0,0,0,.22); font-family:var(--font);}
#rex-tip .rex-card{margin:.4rem 0 .2rem; border-bottom:0; max-width:none;}
:root{--refilter:none;}
.sv{overflow-x:auto; padding-bottom:.4rem;}
.sv-block{margin:0 0 1.1rem;}
.sv-labs{position:relative; margin-left:3.6em;}
.sv-l{position:absolute; transform:translateX(-50%); white-space:nowrap; font-size:var(--s-1); color:var(--ink-2);}
.sv-row{display:flex; align-items:center; white-space:nowrap; line-height:1.4;}
.sv-n{display:inline-block; width:3.6em; text-align:right; padding-right:.7em; color:var(--ink-2); font-size:var(--s-1);}
.sv-n.e{text-align:left; padding-left:.7em;}
.sv-s{display:inline-block;}
.sv .b{display:inline-block; width:.72em; text-align:center; font-style:normal;}
.sv-row + .sv-row .b{color:var(--ink-2);}
.sv .rs{background:var(--code);}
.sv .ct, .sv .cb{border-left:2px solid var(--pk1);}
.sv-feats{position:relative; margin:.25em 0 0 3.6em;}
.sv-bar{position:absolute; height:1.25em; line-height:1.25em; font-size:var(--s-1); text-align:center;
  overflow:hidden; text-overflow:ellipsis; white-space:nowrap; border-radius:2px; filter:var(--refilter);}
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
  --pk1:#D9498B; --pk2:#73274A; --pk-rgb:217,73,139; --refilter:invert(1) hue-rotate(180deg);}
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


st.markdown(
    '<h1 class="rex-title">Restriction Enzyme Explorer</h1>'
    '<p class="rex-dek">Load a DNA sequence, choose a region, and see where known '
    'restriction enzymes cut it.</p>', unsafe_allow_html=True)

# ----------------------------------------------------------------- sidebar
with st.sidebar:
    st.header("Sequence")
    upload = st.file_uploader("FASTA or GenBank file",
                              type=["fa", "fasta", "fna", "txt", "gb", "gbk", "genbank"])
    if upload is not None:
        text = upload.getvalue().decode("utf-8-sig", errors="replace")
    else:
        text = EXAMPLE.read_text(encoding="utf-8")
        st.caption("No file uploaded. Showing pBR322 (GenBank J01749).")

try:
    records = parse_any(text)
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
                           value=(record.topology == "circular"), disabled=not whole) and whole

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

    st.header("Theme")
    st.radio("Theme", ["Light", "Dark"], horizontal=True, key="theme",
             label_visibility="collapsed")

if not names:
    st.warning("Select at least one enzyme.")
    st.stop()

# Map display settings live under the map (widgets below); read their saved values here.
present_types = sorted({f.type for f in record.features})
default_types = [t for t in present_types if t not in HIDDEN_FEATURE_TYPES]
if st.session_state.get("_rec_key") != (record.id, seq_len):
    st.session_state["_rec_key"] = (record.id, seq_len)
    st.session_state["feature_types"] = default_types
st.session_state["feature_types"] = [t for t in st.session_state.get("feature_types", default_types)
                                     if t in present_types]
feature_types = st.session_state["feature_types"]
layout = st.session_state.get("layout", "Wrapped rows")
bp_per_row = st.session_state.get("bp_per_row", MAX_BP_PER_ROW)
band = st.session_state.get("band", "Orange")
label_policy = st.session_state.get("label_policy", "Up to 3 cuts")
hover_cards = st.session_state.get("hover_cards", True)

# ---------------------------------------------------------------- analysis
region_seq = seq[start - 1:end]
all_features = list(record.features)
features = [f for f in all_features if f.type in feature_types]
all_sites = _scan(seq, tuple(names), start, end, circular)
counts: dict[str, int] = {}
site_of: dict[str, str] = {}
first_pos: dict[str, int] = {}
for s in all_sites:
    counts[s.enzyme] = counts.get(s.enzyme, 0) + 1
    site_of.setdefault(s.enzyme, s.site)
    first_pos.setdefault(s.enzyme, s.position)
sites = [s for s in all_sites if max_cuts is None or counts[s.enzyme] <= max_cuts]
cutters = sorted({s.enzyme for s in sites})

hl = st.session_state.get("hl")
if hl not in cutters:
    hl = None

# Which enzymes get a name on the maps (the rest stay visible as notches).
limit = {"Unique cutters (1 cut)": 1, "Up to 3 cuts": 3, "All enzymes": 10 ** 9}[label_policy]
candidates = sorted((n for n in cutters if counts[n] <= limit),
                    key=lambda n: (counts[n], first_pos[n]))
labelled = set(candidates[:MAX_NAMED])
names_capped = len(candidates) > MAX_NAMED
named_sites = [s for s in sites if s.enzyme in labelled or s.enzyme == hl]
cards = ({(s.enzyme, s.position): enzyme_card_html(seq, s.enzyme, _describe(s.enzyme), [s],
                                                   counts[s.enzyme]) for s in named_sites}
         if hover_cards else None)

stats = [("Sequence", f"{seq_len:,} bp"), ("Region", f"{end - start + 1:,} bp"),
         ("GC (region)", f"{gc_percent(region_seq):.1f}%"),
         ("Enzymes cutting", f"{len(cutters)}"), ("Cut sites", f"{len(sites)}")]
if all_features:
    stats.append(("Features", f"{len(all_features)}"))
st.markdown(
    '<dl class="rex-stats">' + "".join(
        f'<div class="rex-stat"><dt>{k}</dt><dd>{v}</dd></div>' for k, v in stats)
    + "</dl>"
    f'<p class="rex-record"><b>{html.escape(record.id)}</b> '
    f'{html.escape(record.description)}</p>', unsafe_allow_html=True)

(tab_map, tab_seq, tab_feat, tab_sites, tab_digest, tab_motif, tab_report) = st.tabs(
    ["Map", "Sequence", "Features", "Sites", "GEP", "Motif search", "Report"])

# --------------------------------------------------------------------- map
with tab_map:
    head_l, head_m, head_r = st.columns([1.3, 2, 1.3], gap="medium")
    with head_l:
        view = st.radio("View", ["Linear map", "Circular map"], horizontal=True,
                        index=1 if record.topology == "circular" else 0)
    with head_m:
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
        hl_sites = [s for s in sites if s.enzyme == hl]
        st.markdown(enzyme_card_html(seq, hl, _describe(hl), hl_sites, counts[hl]),
                    unsafe_allow_html=True)

    theme_name = "dark" if dark else "light"
    if view == "Circular map":
        fig = build_circular_map(record.id, seq_len, "circular" if record.topology == "circular"
                                 else "linear sequence", gc_percent(seq), sites, features,
                                 highlight=hl, labelled=labelled, theme=theme_name,
                                 show_features=bool(features), cards=cards)
        st.plotly_chart(fig, use_container_width=True, config=MAP_CONFIG, theme=None)
        note = ("Arcs are features (arrows show strand); red ticks are cuts. Hover an arc or a "
                "tick for details.")
        _, dropped = split_lanes(features, seq_len, MAX_LANES)
        if dropped:
            note += (f" {dropped} overlapping feature(s) are left off the ring to keep it "
                     "readable; see the Features tab.")
        if not whole:
            note += " The circle shows the whole sequence; cuts are those found in your region."
    else:
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
                        band=band.lower(), theme=theme_name, labelled=labelled, cards=cards)
        st.plotly_chart(fig, use_container_width=True, config=MAP_CONFIG, theme=None)
        note = ("Each band is the sequence; a dark notch marks a cut. Hover a marker for its "
                "recognition site. The camera icon saves a PNG.")
    note += (" Pinch (or Ctrl + scroll) over the map to zoom, drag to pan, double-click to "
             "reset. Plain scrolling always moves the page.")
    if names_capped:
        note += f" Naming {MAX_NAMED} of {len(candidates)} eligible enzymes; the rest are notches."
    st.markdown(f'<p class="rex-note">{note}</p>', unsafe_allow_html=True)
    if not sites:
        st.caption("No cut sites for these enzymes in this region.")

    st.markdown('<h3 class="rex-panel-title" style="margin-top:1.2rem !important">Display</h3>',
                unsafe_allow_html=True)
    linear_view = view == "Linear map"
    d1, d2, d3 = st.columns(3, gap="large")
    with d1:
        st.radio("Layout", ["Wrapped rows", "Single strip"], key="layout",
                 disabled=not linear_view,
                 help="Wrapped rows: up to 1,200 bases per row, the rest continues below. "
                      "Single strip: one row you scroll along.")
    with d2:
        st.select_slider("Bases per row (or in view)",
                         options=[300, 600, 900, MAX_BP_PER_ROW], value=MAX_BP_PER_ROW,
                         key="bp_per_row", disabled=not linear_view)
    with d3:
        st.radio("Band colour", ["Orange", "Green"], horizontal=True, key="band",
                 disabled=not linear_view)
    e1, e2, e3 = st.columns(3, gap="large")
    with e1:
        st.selectbox("Name these enzymes on maps",
                     ["Unique cutters (1 cut)", "Up to 3 cuts", "All enzymes"], index=1,
                     key="label_policy",
                     help="Naming every enzyme crowds the map. Other cuts stay visible as "
                          "notches or ticks.")
    with e2:
        st.multiselect("Feature types on the maps", present_types, key="feature_types",
                       disabled=not present_types,
                       help="Noisy types such as misc_binding are hidden by default. The "
                            "Features tab always lists everything.")
    with e3:
        st.checkbox("Enzyme card on hover", value=True, key="hover_cards",
                    help="Show the enzyme card next to the pointer when you hover a cut or "
                         "its name. Turn off for plain tooltips.")

# ---------------------------------------------------------------- sequence
with tab_seq:
    s1, s2, s3 = st.columns([1, 1, 2], gap="large")
    with s1:
        per_line = st.select_slider("Bases per line", options=[30, 60, 90], value=60)
    with s2:
        jump = int(st.number_input("Go to bp", 1, seq_len, 1))
    with s3:
        st.markdown('<p class="rex-note">Top strand runs 5\u2032 to 3\u2032, bottom strand 3\u2032 '
                    'to 5\u2032. A pink bar marks each cut on its strand and a tint marks the '
                    'recognition site. Features sit under the sequence.</p>',
                    unsafe_allow_html=True)
    ws = ((jump - 1) // per_line) * per_line + 1
    we = min(seq_len, ws + per_line * SEQ_LINES - 1)
    st.markdown(build_sequence_html(seq, ws, we, named_sites, features, per_line),
                unsafe_allow_html=True)
    st.caption(f"Showing bp {ws:,} to {we:,} of {seq_len:,}. Change Go to bp to move along "
               "the sequence.")

# ---------------------------------------------------------------- features
with tab_feat:
    if all_features:
        st.markdown(f'<p class="rex-dek">{len(all_features)} annotated feature(s) in '
                    f'{html.escape(record.id)}.</p>', unsafe_allow_html=True)
        _table([{"Name": f.name, "Type": f.type,
                 "Location": " + ".join(f"{a:,}\u2013{b:,}" for a, b in f.segments),
                 "Strand": {1: "+", -1: "\u2212", 0: ""}[f.strand],
                 "Length (bp)": f.length}
                for f in sorted(all_features, key=lambda f: f.start)])
    else:
        st.info("This sequence has no annotations. Upload a GenBank file (.gb or .gbk) to "
                "see genes, promoters, origins and other features here and on the maps.")

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
    if not cutters:
        st.info("No cutting enzymes to digest with.")
    else:
        digest_enzymes = st.multiselect("Enzymes in the digest", cutters, default=cutters[:1])
        d1, d2, d3 = st.columns(3, gap="large")
        with d1:
            percent = st.select_slider("Agarose gel", options=list(GEL_RANGES), value=1.0,
                                       format_func=lambda v: f"{v}%")
        with d2:
            marker = st.selectbox("Marker", list(LADDERS), index=1)
        with d3:
            palette = st.radio("Gel colours", list(PALETTES), index=1 if dark else 0)
        o1, o2 = st.columns(2, gap="large")
        with o1:
            show_uncut = st.checkbox("Include an uncut lane", value=True)
        with o2:
            per_enzyme = st.checkbox("Add a lane for each single-enzyme digest", value=False,
                                     disabled=len(digest_enzymes) < 2)

        if not digest_enzymes:
            st.caption("Pick at least one enzyme to run a digest.")
        else:
            dig_cuts = [s for s in sites if s.enzyme in digest_enzymes]
            fragments = digest_fragments(dig_cuts, start, end, circular)

            lanes = []
            if LADDERS[marker]:
                lanes.append(("Marker", LADDERS[marker]))
            if show_uncut:
                lanes.append(("Uncut", [end - start + 1]))
            if per_enzyme and len(digest_enzymes) > 1:
                for e in digest_enzymes[:5]:
                    single = digest_fragments([s for s in dig_cuts if s.enzyme == e],
                                              start, end, circular)
                    lanes.append((e, lane_sizes(single)))
            lanes.append((" + ".join(digest_enzymes) if len(digest_enzymes) < 3 else "Digest",
                          lane_sizes(fragments)))

            gel_col, table_col = st.columns([1.25, 1.75], gap="large")
            with gel_col:
                st.plotly_chart(build_gel(lanes, percent, palette), use_container_width=True,
                                config=PLOT_CONFIG, theme=None)
            with table_col:
                sizes = [f.length for f in fragments]
                st.markdown(
                    f'<p class="rex-dek">{len(fragments)} fragment(s) \u00B7 largest '
                    f'{max(sizes):,} bp \u00B7 smallest {min(sizes):,} bp</p>',
                    unsafe_allow_html=True)
                _table([{"#": i, "Ends": f"{f.left} \u2013 {f.right}",
                         "Coordinates": format_coordinates(f), "Length (bp)": f.length}
                        for i, f in enumerate(fragments, 1)])
                tiny = [f for f in fragments if f.length < BAND_THRESHOLD_BP]
                if tiny:
                    st.caption(f"{len(tiny)} fragment(s) under {BAND_THRESHOLD_BP} bp would be "
                               "hard to see on a real gel.")
                if len(set(sizes)) < len(sizes):
                    st.caption("Fragments of the same size run as one band.")
                csv_rows = ["fragment,left_end,right_end,start,end,length_bp"] + [
                    f"{i},{f.left},{f.right},{f.start},{f.end},{f.length}"
                    for i, f in enumerate(fragments, 1)]
                fasta = "".join(
                    f">{record.id}_fragment_{i} {f.start}-{f.end} {f.length} bp {f.left}-{f.right}\n"
                    + "\n".join(fragment_sequence(seq, f)[k:k + 70]
                                for k in range(0, f.length, 70)) + "\n"
                    for i, f in enumerate(fragments, 1))
                b1, b2 = st.columns(2)
                with b1:
                    st.download_button("Fragments (CSV)", "\n".join(csv_rows) + "\n",
                                       file_name=f"{record.id}_fragments.csv", mime="text/csv")
                with b2:
                    st.download_button("Fragments (FASTA)", fasta,
                                       file_name=f"{record.id}_fragments.fasta",
                                       mime="text/plain")
            st.markdown('<p class="rex-note">The gel is an empirical model of agarose '
                        'migration, not a measurement. Uncut circular DNA would run differently '
                        'on a real gel (supercoiled forms). Band thickness follows DNA mass.</p>',
                        unsafe_allow_html=True)

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

# Pinch-to-zoom without hijacking page scroll, and the hover enzyme card (see INTERACTION_HTML).
try:
    st.components.v1.html(INTERACTION_HTML, height=0)
except Exception:
    pass
