"""Dependency-free FASTA parsing and basic sequence statistics."""
from __future__ import annotations

from dataclasses import dataclass

VALID_BASES = set("ACGTRYSWKMBDHVN")  # IUPAC DNA alphabet


@dataclass(frozen=True)
class FastaRecord:
    """A sequence record. GenBank input also fills topology and features."""

    id: str
    description: str
    sequence: str
    topology: str = "linear"
    features: tuple = ()


def clean_sequence(raw: str, label: str = "unnamed") -> str:
    """Upper-case, strip whitespace and validate an IUPAC DNA sequence."""
    seq = "".join(raw.split()).upper()
    if not seq:
        raise ValueError(f"Record '{label}' has no sequence.")
    invalid = sorted(set(seq) - VALID_BASES)
    if invalid:
        hint = " (RNA? convert U to T first)" if "U" in invalid else ""
        raise ValueError(
            f"Record '{label}' contains non-DNA characters: {', '.join(invalid)}{hint}")
    return seq


def _build_record(header: str, chunks: list[str]) -> FastaRecord:
    seq = clean_sequence("".join(chunks), header or "unnamed")
    parts = header.split(None, 1)
    return FastaRecord(
        id=parts[0] if parts else "unnamed",
        description=parts[1] if len(parts) > 1 else "",
        sequence=seq,
    )


def parse_fasta(text: str) -> list[FastaRecord]:
    """Parse FASTA text. A header-less sequence is accepted as one record."""
    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(_build_record(header, chunks))
            header, chunks = line[1:].strip(), []
        else:
            if header is None:
                header = "unnamed"
            chunks.append(line)
    if header is not None:
        records.append(_build_record(header, chunks))
    if not records:
        raise ValueError("No sequence found. Is this a FASTA file?")
    return records


def gc_percent(seq: str) -> float:
    """GC% over unambiguous A/C/G/T bases only."""
    counted = sum(seq.count(b) for b in "ACGT")
    if counted == 0:
        return 0.0
    return 100.0 * (seq.count("G") + seq.count("C")) / counted
