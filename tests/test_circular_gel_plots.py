from pathlib import Path

import pytest

pytest.importorskip("plotly")

from rex.circular import build_circular_map  # noqa: E402
from rex.gelplot import PALETTES, build_gel  # noqa: E402
from rex.genbank import parse_any  # noqa: E402
from rex.models import CutSite  # noqa: E402

EXAMPLE = Path(__file__).parent / "data" / "synthetic_plasmid.gb"


def test_circular_map_builds_with_and_without_features():
    rec = parse_any(EXAMPLE.read_text())[0]
    sites = [CutSite("EcoRI", 244, "GAATTC", "G^AATT_C")]
    with_f = build_circular_map(rec.id, len(rec.sequence), "circular", 49.0, sites,
                                list(rec.features), highlight="EcoRI")
    without = build_circular_map(rec.id, len(rec.sequence), "circular", 49.0, sites, [],
                                 show_features=False)
    assert len(with_f.data) > len(without.data)


def test_gel_builds_for_every_palette():
    for palette in PALETTES:
        fig = build_gel([("Marker", [10000, 1000, 500]), ("Digest", [1500, 900])], 1.0, palette)
        assert len(fig.layout.shapes) > 5


def test_cards_ride_along_on_hover_points_instead_of_plotly_tooltips():
    from rex.plotting import build_map, plan_rows
    rec = parse_any(EXAMPLE.read_text())[0]
    site = CutSite("EcoRI", 244, "GAATTC", "G^AATT_C", 248)
    cards = {("EcoRI", 244): '<div class="rex-card">EcoRI</div>'}
    circ = build_circular_map(rec.id, len(rec.sequence), "circular", 49.0, [site],
                              list(rec.features), labelled={"EcoRI"}, cards=cards)
    assert any(t.customdata is not None and any("rex-card" in str(c) for c in t.customdata)
               for t in circ.data if t.mode == "markers")
    lin = build_map(plan_rows((1, 600)), [site], labelled={"EcoRI"}, cards=cards)
    carded = [t for t in lin.data if t.customdata is not None
              and any("rex-card" in str(c[-1]) for c in t.customdata)]
    assert carded and all(t.hoverinfo == "none" for t in carded)
