from rex.fragments import compute_fragments
from rex.motif import find_motif


def test_motif_forward_and_reverse():
    hits = find_motif("AAGGATCCAA", "GGATCC")  # palindrome -> reported once
    assert [(h.strand, h.start, h.end) for h in hits] == [("+", 3, 8)]
    assert find_motif("GAAATTT", "AAAC") == []  # neither strand contains it
    hits = find_motif("AAAGTTT", "AAAC")  # rc(AAAC)=GTTT found on minus strand
    assert [(h.strand, h.start) for h in hits] == [("-", 4)]


def test_motif_ambiguity_and_overlap():
    assert len(find_motif("AAAA", "AA")) == 3
    assert len(find_motif("GACGTC", "GRCGYC")) == 1


def test_motif_respects_region_offset():
    hits = find_motif("GAATTCAAGAATTC", "GAATTC", start=5, end=14)
    assert [h.start for h in hits] == [9]


def test_linear_fragments():
    frags = compute_fragments([41], 1, 100)
    assert [(f.start, f.end, f.length) for f in frags] == [(1, 40, 40), (41, 100, 60)]
    assert sum(f.length for f in compute_fragments([10, 50, 90], 1, 100)) == 100


def test_circular_fragments():
    one = compute_fragments([30], 1, 100, circular=True)
    assert [f.length for f in one] == [100]
    two = compute_fragments([20, 70], 1, 100, circular=True)
    assert sorted(f.length for f in two) == [50, 50]
    assert sum(f.length for f in compute_fragments([20, 30, 90], 1, 100, True)) == 100
