# Restriction Enzyme Explorer

An interactive bioinformatics tool for finding and visualising restriction-enzyme
sites in DNA sequences. Upload a FASTA file, select a region, scan it against
Biopython's REBASE-derived enzyme database, and export a research report.

## Features

- **FASTA upload** (multi-record supported, IUPAC DNA validated; bundled synthetic example included)
- **Region selection** with a slider; sites are reported in full-sequence coordinates
- **Real enzyme database** via `Bio.Restriction`: recognition sites, cut patterns, overhang type, optimum temperature, supplier codes
- **Restriction map**: the sequence is a flat orange or green band with a dark notch at every cut
- **Two layouts**: wrapped rows of up to 1,200 bp with the rest continuing below, or a single strip you scroll along
- **Enzyme window**: an open/close popup on the right lists every enzyme in view; click one to light up its cuts with a pink gradient glow
- **Digest simulator**: fragment sizes and a proportional fragment bar for any combination of enzymes, linear or circular
- **Custom motif search** with IUPAC ambiguity codes, both strands
- **Report export**: Markdown report (sequence ID, region, sites, enzyme references, limitations) and CSV of sites
- **Light and dark theme** (Display > Theme in the sidebar)
- **PNG export** of the map from the camera icon on the plot

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
rex/motif.py        IUPAC motif search
rex/analysis.py     Biopython enzyme lookup and site finding
rex/fragments.py    Digest fragment sizes
rex/plotting.py     Restriction map
.streamlit/config.toml  Theme (warm paper, dark-gray text)
rex/report.py       Markdown/CSV export
data/example.fasta  Synthetic demo sequence (not a real gene)
tests/              Unit tests
```

## Limitations

- Cut position = 1-based position of the first base after the top-strand cut (Biopython convention).
- Only the selected region is scanned, so sites straddling its boundary are missed.
- Methylation sensitivity, star activity, buffers and partial digests are not modelled.
- only available for DNA, not open to RNA or protein models
- Predictions are in silico; confirm in the lab. Not for clinical use.

## Roadmap

- [ ] GenBank input with feature annotations (genes/CDS on the map)
- [ ] Double-digest buffer compatibility
- [ ] Methylation-sensitivity flags
- [ ] Virtual agarose gel image
- [ ] Docker image / Streamlit Community Cloud deployment


## Credits

Cock PJA et al., *Biopython*, Bioinformatics 25(11):1422 (2009). Enzyme data from REBASE
(Roberts RJ et al., Nucleic Acids Res. 2015). Docs: https://biopython.org/docs/latest/api/Bio.Restriction.html
