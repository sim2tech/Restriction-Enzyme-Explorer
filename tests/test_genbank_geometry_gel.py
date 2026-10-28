from pathlib import Path

from rex.fragments import digest_fragments, format_coordinates, fragment_sequence
from rex.gel import axis_ticks, band_thickness_mm, migration_mm
from rex.genbank import parse_any, parse_location
from rex.geometry import assign_lanes, feature_polygon, nice_step, spread_1d
from rex.models import CutSite, Feature
from rex.seqview import build_sequence_html, complement

EXAMPLE = Path(__file__).parent / "data" / "synthetic_plasmid.gb"


def test_parse_location_variants():
    assert parse_location("101..120") == (1, [(101, 120)])
    assert parse_location("complement(1600..2460)") == (-1, [(1600, 2460)])
    assert parse_location("join(2761..2800,1..100)") == (1, [(2761, 2800), (1, 100)])
    assert parse_location("<1..>50") == (1, [(1, 50)])
    assert parse_location("garbage") is None


def test_example_genbank_loads():
    rec = parse_any(EXAMPLE.read_text())[0]
    assert rec.topology == "circular" and len(rec.sequence) == 2800
    names = {f.name: f for f in rec.features}
    assert names["T7 promoter"].start == 101 and names["AmpR"].strand == -1
    assert names["backbone"].segments == ((2761, 2800), (1, 100))
    assert rec.sequence[100:120] == "TAATACGACTCACTATAGGG"


def test_fasta_still_detected():
    assert parse_any(">x\nACGT\n")[0].sequence == "ACGT"


def test_lanes_do_not_overlap():
    feats = [Feature("a", "CDS", 1, 100), Feature("b", "CDS", 50, 150), Feature("c", "CDS", 200, 300)]
    lanes = assign_lanes(feats, 1000)
    assert lanes[0] != lanes[1] and lanes[2] == 0


def test_spread_1d_keeps_spacing_and_bounds():
    out = spread_1d([0.5, 0.5, 0.5, 0.49], 0.1, -1, 1)
    srt = sorted(out, reverse=True)
    assert all(a - b >= 0.1 - 1e-9 for a, b in zip(srt, srt[1:]))
    assert all(-1 <= v <= 1 for v in out)
    packed = spread_1d([-0.99, -0.99, -0.99], 0.1, -1, 1)
    assert min(packed) >= -1 and max(packed) <= 1


def test_polygon_closed_and_nice_step():
    xs, ys = feature_polygon(0.2, 1.0, 0.8, 0.9, strand=1)
    assert (xs[0], ys[0]) == (xs[-1], ys[-1])
    assert nice_step(2800) == 500 and nice_step(100) == 25


def test_digest_ends_and_wraparound():
    cuts = [CutSite("EcoRI", 20, "GAATTC", "G^AATT_C"), CutSite("BamHI", 70, "GGATCC", "G^GATC_C")]
    lin = digest_fragments(cuts, 1, 100, circular=False)
    assert [(f.left, f.right, f.length) for f in lin] == [
        ("start", "EcoRI", 19), ("EcoRI", "BamHI", 50), ("BamHI", "end", 31)]
    circ = digest_fragments(cuts, 1, 100, circular=True)
    assert sorted(f.length for f in circ) == [50, 50]
    wrap = [f for f in circ if f.start > f.end][0]
    assert (wrap.left, wrap.right) == ("BamHI", "EcoRI") and format_coordinates(wrap) == "70\u201319"
    seq = "".join("ACGT"[i % 4] for i in range(100))
    assert len(fragment_sequence(seq, wrap)) == 50


def test_gel_monotonic_and_ticks():
    sizes = [100, 300, 500, 1000, 3000, 5000, 10000, 20000]
    ys = [migration_mm(b, 1.0) for b in sizes]
    assert ys == sorted(ys, reverse=True) and all(0 < y < 100 for y in ys)
    # 300 bp runs near the end of a low-percentage gel but sits mid-gel (resolved) on 2%.
    assert migration_mm(300, 0.7) > 88 and 30 < migration_mm(300, 2.0) < 80
    assert band_thickness_mm(2000, 2500) > band_thickness_mm(100, 2500)
    assert any(major for _, _, major in axis_ticks(1.0))


def test_sequence_html_marks_cuts_both_strands_and_features():
    seq = "AAAAGAATTCAAAA"
    site = CutSite("EcoRI", 6, "GAATTC", "G^AATT_C", bottom_position=10)
    html_out = build_sequence_html(seq, 1, len(seq), [site],
                                   [Feature("x", "CDS", 3, 8, 1)], per_row=60)
    assert "EcoRI" in html_out and 'class="b ct"' in html_out.replace("rs ", "")
    assert "cb" in html_out and "sv-bar" in html_out
    assert complement("ACGT") == "TGCA"


PBR322 = Path(__file__).parent.parent / "data" / "pBR322.gb"


def test_pbr322_default_loads():
    rec = parse_any(PBR322.read_text())[0]
    assert (rec.id, len(rec.sequence), rec.topology) == ("SYNPBR322", 4361, "circular")
    bla = [f for f in rec.features if f.name == "bla" and f.type == "CDS"][0]
    assert (bla.start, bla.end, bla.strand) == (3293, 4153, -1)
    assert not any(f.name == "note" for f in rec.features)


def test_default_feature_filter_and_lane_cap():
    from rex.geometry import HIDDEN_FEATURE_TYPES, MAX_LANES, split_lanes
    rec = parse_any(PBR322.read_text())[0]
    shown = [f for f in rec.features if f.type not in HIDDEN_FEATURE_TYPES]
    assert 0 < len(shown) < len(rec.features)
    kept, dropped = split_lanes(shown, len(rec.sequence))
    assert max(lane for _, lane in kept) < MAX_LANES and len(kept) + dropped == len(shown)
    overlapping = [Feature(str(i), "CDS", 10, 200) for i in range(6)]
    kept, dropped = split_lanes(overlapping, 1000, 3)
    assert len(kept) == 3 and dropped == 3


def test_text_angle_keeps_text_upright_and_bezier_ends():
    import math
    from rex.geometry import bezier, text_angle
    assert text_angle(0) == 0 and text_angle(math.pi / 2) == 90
    assert text_angle(math.pi) == 0 and abs(text_angle(math.radians(135)) + 45) < 1e-9
    assert all(-90 <= text_angle(math.radians(a)) <= 90 for a in range(0, 360, 7))
    pts = bezier((0, 0), (1, 2), (2, 0))
    assert pts[0] == (0, 0) and pts[-1] == (2, 0)


def test_enzyme_card_has_strands_cut_and_isoschizomers():
    from rex.models import EnzymeInfo
    from rex.seqview import enzyme_card_html
    seq = "AAAAAAAAGAATTCTTTTTTTT"
    site = CutSite("EcoRI", 10, "GAATTC", "G^AATT_C", bottom_position=14)
    info = EnzymeInfo("EcoRI", "GAATTC", "G^AATT_C", 6, "5' overhang", "AATT", (), True,
                      isoschizomers=("MunI-like",))
    out = enzyme_card_html(seq, "EcoRI", info, [site], 1)
    assert "unique cutter" in out and "cuts after" in out and "MunI-like" in out
    assert "G^AATT_C" in out and 'class="b ct' in out.replace("rs ", "") and "cb" in out
