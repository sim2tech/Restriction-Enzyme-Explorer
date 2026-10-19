# Restriction Enzyme Explorer

An interactive bioinformatics tool for finding and visualising restriction-enzyme
sites in DNA sequences. Upload a FASTA file, select a region, scan it against
Biopython's REBASE-derived enzyme database, and export a research report.

## Features

- **FASTA upload** (multi-record supported, IUPAC DNA validated; bundled synthetic example included)
- **Region selection** with a slider; sites are reported in full-sequence coordinates
- **Real enzyme database** via `Bio.Restriction`: recognition sites, cut patterns, overhang type, optimum temperature, supplier codes
- **Interactive map** (Plotly): labelled cut sites, highlighted region, zoom range slider, clickable legend
- **Digest simulator**: fragment sizes for any combination of enzymes, linear or circular
- **Custom motif search** with IUPAC ambiguity codes, both strands
- **Report export**: Markdown report (sequence ID, region, sites, enzyme references, limitations) and CSV of sites

## Quick start

```bash
git clone https://github.com/sim2tech/Restriction-Enzyme-Explorer
cd restriction-enzyme-explorer
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

Run the tests with `pytest`.

## Project layout

```
app.py              Streamlit interface
rex/fasta.py        FASTA parsing and GC content
rex/motif.py        IUPAC motif search
rex/analysis.py     Biopython enzyme lookup and site finding
rex/fragments.py    Digest fragment sizes
rex/plotting.py     Interactive restriction map
rex/report.py       Markdown/CSV export
data/example.fasta  Synthetic demo sequence (not a real gene)
tests/              Unit tests
```

## Conventions and limitations

- Cut position = 1-based position of the first base after the top-strand cut (Biopython convention).
- Only the selected region is scanned, so sites straddling its boundary are missed.
- Methylation sensitivity, star activity, buffers and partial digests are not modelled.
- Predictions are in silico; confirm in the lab. Not for clinical use.

## Roadmap

- [ ] GenBank input with feature annotations (genes/CDS on the map)
- [ ] Double-digest buffer compatibility
- [ ] Methylation-sensitivity flags
- [ ] Virtual agarose gel image
- [ ] Docker image / Streamlit Community Cloud deployment

      ## Limitations

This tool is for research and education. Its results are **in silico predictions**, not lab confirmation.

- **Enzyme database:** enzymes come from the REBASE data bundled with the installed Biopython version. The list is not exhaustive and may not be fully current.
- **Region scanning:** only the selected region is searched. Recognition sites that straddle the region boundary are missed.
- **Cut positions:** a cut position is the 1-based position of the first base after the top-strand cut (Biopython convention). Fragment sizes ignore overhang lengths.
- **Enzymes that cut away from their site:** enzymes that cut outside their recognition sequence (e.g. Type IIS) may report cuts near the region edges inaccurately.
- **Not modelled:** methylation sensitivity (Dam/Dcm/CpG), star activity, buffer conditions, partial digests and DNA supercoiling.
- **Ambiguous bases:** N, R, Y and other IUPAC codes are matched by IUPAC rules and can add or hide sites.
- **Input:** DNA FASTA only (IUPAC alphabet). RNA, protein and GenBank files are not supported yet.
- **Not for clinical or diagnostic use.** Always confirm digests experimentally.

## Credits

Cock PJA et al., *Biopython*, Bioinformatics 25(11):1422 (2009). Enzyme data from REBASE
(Roberts RJ et al., Nucleic Acids Res. 2015). Docs: https://biopython.org/docs/latest/api/Bio.Restriction.html
