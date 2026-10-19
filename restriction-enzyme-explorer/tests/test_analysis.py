import pytest

pytest.importorskip("Bio")

from rex.analysis import describe_enzyme, find_sites, list_enzymes  # noqa: E402


def test_ecori_position():
    sites = find_sites("AAGAATTCAA", ["EcoRI"])
    assert [(s.enzyme, s.position) for s in sites] == [("EcoRI", 4)]


def test_region_offsets_positions():
    sites = find_sites("AAGAATTCAA", ["EcoRI"], start=3, end=10)
    assert sites[0].position == 4


def test_region_excludes_outside_sites():
    seq = "GAATTC" + "A" * 20 + "GGATCC"
    assert find_sites(seq, ["EcoRI"], 10, len(seq)) == []
    assert [s.position for s in find_sites(seq, ["BamHI"], 10, len(seq))] == [28]


def test_enzyme_metadata():
    info = describe_enzyme("EcoRI")
    assert info.site == "GAATTC"
    assert info.overhang == "5' overhang"
    assert describe_enzyme("EcoRV").overhang == "blunt"
    assert "EcoRI" in list_enzymes(commercial_only=False)


def test_unknown_enzyme():
    with pytest.raises(ValueError):
        find_sites("ACGT", ["NotAnEnzyme"])
