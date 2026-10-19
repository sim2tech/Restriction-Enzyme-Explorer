<img width="1476" height="745" alt="image" src="https://github.com/user-attachments/assets/2e0c4efe-8008-414c-b44b-cdda6afd3f3c" /># Restriction Enzyme Explorer

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
git clone <https://github.com/sim2tech/Restriction-Enzyme-Explorer>
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

## Credits

Cock PJA et al., *Biopython*, Bioinformatics 25(11):1422 (2009). Enzyme data from REBASE
(Roberts RJ et al., Nucleic Acids Res. 2015). Docs: https://biopython.org/docs/latest/api/Bio.Restriction.html
