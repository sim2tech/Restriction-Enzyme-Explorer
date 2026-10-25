import pytest

pytest.importorskip("plotly")

from rex.models import CutSite  # noqa: E402
from rex.plotting import build_map, plan_rows  # noqa: E402


def test_plan_rows_wraps_at_1200():
    assert plan_rows((1, 2500)) == [(1, 1200), (1201, 2400), (2401, 2500)]
    assert plan_rows((1, 900)) == [(1, 900)]
    assert plan_rows((1, 5000), 5000)[0] == (1, 1200)  # capped at 1,200 per row


def test_build_map_handles_highlight_and_empty():
    sites = [CutSite("EcoRI", 302, "GAATTC", "G^AATT_C")]
    plain = build_map(plan_rows((1, 2500)), sites)
    lit = build_map(plan_rows((1, 2500)), sites, highlight="EcoRI", band="green")
    assert len(lit.data) > len(plain.data)       # pink halo traces added
    assert len(lit.layout.shapes) > len(plain.layout.shapes)
    build_map([(1, 300)], [], bp_per_row=300)    # no sites: still draws the band


def test_dark_theme_uses_dark_background():
    fig = build_map([(1, 300)], [], bp_per_row=300, theme="dark")
    assert fig.layout.paper_bgcolor != build_map([(1, 300)], [], bp_per_row=300).layout.paper_bgcolor
