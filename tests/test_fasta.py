import pytest

from rex.fasta import gc_percent, parse_fasta


def test_multiline_and_description():
    recs = parse_fasta(">seq1 my plasmid\nACGT\nacgt\n")
    assert recs[0].id == "seq1"
    assert recs[0].description == "my plasmid"
    assert recs[0].sequence == "ACGTACGT"


def test_multiple_records_and_headerless():
    assert [r.id for r in parse_fasta(">a\nAC\n>b\nGT\n")] == ["a", "b"]
    assert parse_fasta("ACGT")[0].sequence == "ACGT"


def test_rejects_protein_and_rna():
    with pytest.raises(ValueError):
        parse_fasta(">p\nMKLVE\n")
    with pytest.raises(ValueError, match="RNA"):
        parse_fasta(">r\nACGU\n")


def test_empty_input():
    with pytest.raises(ValueError):
        parse_fasta("\n\n")


def test_gc():
    assert gc_percent("GGCC") == 100.0
    assert gc_percent("ATGC") == 50.0
