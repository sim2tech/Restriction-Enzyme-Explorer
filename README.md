# 🧬 Restriction Enzyme Explorer

An interactive bioinformatics tool for finding and visualising restriction-enzyme
sites in DNA sequences. Upload a FASTA file, select a region, scan it against
Biopython's REBASE-derived enzyme database, and export a research report.

## Features

- **FASTA or GenBank input**: GenBank files bring topology and feature annotations (genes, promoters, origins). pBR322 (GenBank J01749) is bundled as the demo
- **Region selection** with a slider; sites are reported in full-sequence coordinates
- **Real enzyme database** via `Bio.Restriction`: recognition sites, cut patterns, overhang type, optimum temperature, supplier codes
- **Linear map**: wrapped rows of up to 1,200 bp or a scrollable strip, flat orange or green band, a notch at every cut
- **Circular plasmid map**: feature arcs with strand arrows and names written on the arcs, ticks that follow the ring, red cut ticks with curved leader lines and hover details; feature types can be filtered
- **Enzyme detail card**: pick an enzyme to see its cut position, ends, both-strand snippet around the cut, cut pattern and enzymes with the same site
- **Zoom and pan**: pinch (or Ctrl + scroll) over a map to zoom, drag to pan, double-click to reset; plain scrolling still moves the page
- **Hover cards**: hover a cut or an enzyme name on either map to see its card (position, ends, both-strand snippet, pattern, same-site enzymes)
- **Display controls** sit under the map; only the theme switch is in the sidebar
- **Sequence view**: both strands, cut marks on each strand, recognition-site tint, features under the sequence
- **GEP tab** (digest and gel): multi-enzyme digests with fragment ends, coordinates (including circular wrap-around) and sizes; fragments download as CSV or FASTA
- **Virtual agarose gel**: log-scale size axis, marker ladders, agarose percentage, per-enzyme lanes, black-on-paper, white-on-black or orange-on-black
- **Enzyme window**: an open/close popup lists every enzyme in view; click one to light up its cuts with a pink gradient glow
- **Custom motif search** with IUPAC ambiguity codes, both strands
- **Report export**: Markdown report (sequence ID, region, sites, fragment ends, enzyme references, limitations) and CSV of sites
- **Light and dark theme** and **PNG export** of the maps and gel

## Quick start

```bash
git clone https://github.com/sim2tech/Restriction-Enzyme-Explorer.git
cd Restriction-Enzyme-Explorer
python -m venv .venv
```

Windows (PowerShell):

```powershell
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python -m streamlit run app.py
```

macOS / Linux:

```bash
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

Run the tests with `pytest`. The `.streamlit/config.toml` file sets the warm paper theme and blue-grey accent (text is set in Arial Narrow, so it looks best on Windows/Office machines), so keep it next to `app.py`.

## Project layout

```
app.py              Streamlit interface
rex/fasta.py        FASTA parsing and GC content
rex/genbank.py      GenBank reader (topology, features)
rex/motif.py        IUPAC motif search
rex/analysis.py     Biopython enzyme lookup and site finding
rex/fragments.py    Digest fragment sizes
rex/plotting.py     Linear restriction map
rex/circular.py     Circular plasmid map
rex/seqview.py      Sequence view (both strands, cuts, features)
rex/gel.py          Virtual gel model
rex/gelplot.py      Virtual gel figure
rex/geometry.py     Layout helpers for the circular map
.streamlit/config.toml  Theme (warm paper, dark-gray text)
rex/report.py       Markdown/CSV export
data/pBR322.gb      pBR322 complete sequence with annotations (NCBI GenBank J01749)
data/example.fasta  Synthetic 1,200 bp demo sequence
tests/              Unit tests (tests/data holds a synthetic test plasmid)
```

## Limitations

- Cut position = 1-based position of the first base after the top-strand cut (Biopython convention).
- Only the selected region is scanned, so sites straddling its boundary are missed.
- Bottom-strand cut marks are shown for palindromic enzymes only.
- The gel is an empirical model of agarose migration, not a measurement; uncut circular DNA would run differently on a real gel.
- Methylation sensitivity is not modelled.
- Features appear on the circular map and sequence view; the linear map does not draw them yet.
- Methylation sensitivity, star activity, buffers and partial digests are not modelled.
- Predictions are in silico; confirm in the lab. Not for clinical use.

## Roadmap

- [x] GenBank input with feature annotations
- [x] Circular map, sequence view and virtual gel
- [ ] Features on the linear map
- [ ] Double-digest buffer compatibility
- [ ] Methylation-sensitivity flags
- [ ] Primer design and PCR product preview
- [ ] Docker image / Streamlit Community Cloud deployment

## Credits

Cock PJA et al., *Biopython*, Bioinformatics 25(11):1422 (2009). Enzyme data from REBASE
(Roberts RJ et al., Nucleic Acids Res. 2015). Docs: https://biopython.org/docs/latest/api/Bio.Restriction.html
