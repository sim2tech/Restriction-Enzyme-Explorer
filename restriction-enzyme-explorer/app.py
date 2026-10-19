"""Restriction Enzyme Explorer — Streamlit app.

Run with:  streamlit run app.py
"""
from __future__ import annotations

from pathlib import Path

import Bio
import streamlit as st

from rex.analysis import COMMON_ENZYMES, describe_enzyme, find_sites, list_enzymes
from rex.fasta import gc_percent, parse_fasta
from rex.fragments import compute_fragments
from rex.motif import find_motif
from rex.plotting import build_map
from rex.report import build_markdown_report, sites_to_csv

EXAMPLE = Path(__file__).parent / "data" / "example.fasta"

st.set_page_config(page_title="Restriction Enzyme Explorer", page_icon="🧬",
                   layout="wide")


@st.cache_data
def _enzyme_names(commercial_only: bool) -> list[str]:
    return list_enzymes(commercial_only)


@st.cache_data(show_spinner="Scanning sequence…")
def _scan(seq: str, names: tuple, start: int, end: int, circular: bool):
    return find_sites(seq, list(names), start, end, circular)


@st.cache_data
def _describe(name: str):
    return describe_enzyme(name)


st.title("🧬 Restriction Enzyme Explorer")
st.caption("Upload a FASTA file, choose a region, and map known restriction "
           "sites using Biopython's REBASE-derived enzyme database.")

# ----------------------------------------------------------------- sidebar
with st.sidebar:
    st.header("1 · Sequence")
    upload = st.file_uploader("FASTA file", type=["fa", "fasta", "fna", "txt"])
    if upload is not None:
        text = upload.getvalue().decode("utf-8-sig", errors="replace")
    else:
        text = EXAMPLE.read_text()
        st.info("No file uploaded — showing the bundled synthetic example.")

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

    st.header("2 · Region")
    start, end = st.slider("Region (bp)", 1, seq_len, (1, seq_len))
    whole = (start, end) == (1, seq_len)
    circular = st.checkbox("Treat as circular (whole sequence only)",
                           value=False, disabled=not whole) and whole

    st.header("3 · Enzymes")
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

if not names:
    st.warning("Select at least one enzyme.")
    st.stop()

# ---------------------------------------------------------------- analysis
region_seq = seq[start - 1:end]
all_sites = _scan(seq, tuple(names), start, end, circular)
counts: dict[str, int] = {}
for s in all_sites:
    counts[s.enzyme] = counts.get(s.enzyme, 0) + 1
sites = [s for s in all_sites if max_cuts is None or counts[s.enzyme] <= max_cuts]
cutters = sorted({s.enzyme for s in sites})

c1, c2, c3, c4 = st.columns(4)
c1.metric("Sequence", f"{seq_len:,} bp")
c2.metric("Region", f"{end - start + 1:,} bp")
c3.metric("GC (region)", f"{gc_percent(region_seq):.1f}%")
c4.metric("Enzymes cutting / sites", f"{len(cutters)} / {len(sites)}")
st.markdown(f"**{record.id}** {record.description}")

tab_map, tab_sites, tab_digest, tab_motif, tab_report = st.tabs(
    ["Map", "Sites", "Digest", "Motif search", "Report"])

with tab_map:
    if sites:
        st.plotly_chart(build_map(record.id, seq_len, (start, end), sites),
                        use_container_width=True)
        st.caption("Blue band = selected region. Click legend entries to hide "
                   "enzymes; drag the lower slider to zoom.")
    else:
        st.info("No sites found for these enzymes in this region.")

with tab_sites:
    if sites:
        rows = []
        for s in sites:
            info = _describe(s.enzyme)
            rows.append({"Enzyme": s.enzyme, "Cut position (bp)": s.position,
                         "Recognition site": s.site, "Cut pattern": s.cut_pattern,
                         "Ends": info.overhang,
                         "Cuts in region": counts[s.enzyme]})
        st.dataframe(rows, use_container_width=True, hide_index=True)
        st.download_button("Download sites (CSV)", sites_to_csv(sites),
                           file_name=f"{record.id}_sites.csv", mime="text/csv")
    else:
        st.info("Nothing to show.")

digest_enzymes: list[str] = []
fragments = []
with tab_digest:
    st.write("Simulate a digest with one or more of the enzymes that cut.")
    if cutters:
        digest_enzymes = st.multiselect("Enzymes in digest", cutters,
                                        default=cutters[:1])
        cuts = [s.position for s in sites if s.enzyme in digest_enzymes]
        fragments = compute_fragments(cuts, start, end, circular)
        st.dataframe([{"Fragment": i, "Start": f.start, "End": f.end,
                       "Length (bp)": f.length}
                      for i, f in enumerate(fragments, 1)],
                     use_container_width=True, hide_index=True)
    else:
        st.info("No cutting enzymes to digest with.")

with tab_motif:
    motif = st.text_input("Custom motif (IUPAC codes allowed, e.g. GAATTC or GRCGYC)")
    both = st.checkbox("Search both strands", value=True)
    if motif:
        try:
            hits = find_motif(seq, motif, start, end, both)
        except ValueError as err:
            st.error(str(err))
        else:
            st.write(f"{len(hits)} match(es) in bp {start:,}–{end:,}")
            if hits:
                st.dataframe([{"Strand": h.strand, "Start": h.start, "End": h.end,
                               "Matched": h.matched} for h in hits[:1000]],
                             use_container_width=True, hide_index=True)

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
    with st.expander("Preview", expanded=True):
        st.markdown(report)
