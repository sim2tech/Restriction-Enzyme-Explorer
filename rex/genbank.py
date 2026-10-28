"""Small GenBank reader (no third-party imports): sequence, topology and features."""
from __future__ import annotations

import re

from .fasta import FastaRecord, clean_sequence, parse_fasta
from .models import Feature

_SKIP_TYPES = {"source"}
_NAME_KEYS = ("label", "gene", "product", "standard_name", "note")


def _split_top_level(text: str) -> list[str]:
    parts, depth, cur = [], 0, ""
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    if cur:
        parts.append(cur)
    return parts


def parse_location(loc: str):
    """Return (strand, parts) for a GenBank location string, or None if unsupported."""
    s = "".join(loc.split())
    strand = 1
    if s.startswith("complement(") and s.endswith(")"):
        strand, s = -1, s[len("complement("):-1]
    if s.startswith(("join(", "order(")):
        pieces = _split_top_level(s[s.index("(") + 1:-1])
    else:
        pieces = [s]
    parts = []
    for piece in pieces:
        if piece.startswith("complement(") and piece.endswith(")"):
            strand, piece = -1, piece[len("complement("):-1]
        piece = piece.replace("<", "").replace(">", "")
        m = re.fullmatch(r"(\d+)\.\.(\d+)", piece) or re.fullmatch(r"(\d+)()", piece)
        if not m:
            return None
        a = int(m.group(1))
        b = int(m.group(2)) if m.group(2) else a
        parts.append((a, b))
    return (strand, parts) if parts else None


def _feature_name(ftype: str, quals: dict) -> str:
    for key in _NAME_KEYS:
        value = quals.get(key, "").strip()
        if not value:
            continue
        if key == "note":
            if len(value) <= 40:     # long notes are descriptions, not names
                return value
            continue
        return value if len(value) <= 32 else value[:30] + "…"
    return ftype


def _parse_features(lines: list[str], seq_len: int) -> tuple:
    raw, cur = [], None
    for line in lines:
        if not line.strip():
            continue
        if line.startswith("     ") and len(line) > 5 and line[5] != " ":
            cur = {"type": line[5:21].strip(), "loc": line[21:].strip(), "q": {}, "last": None}
            raw.append(cur)
        elif cur is not None:
            text = line.strip()
            if text.startswith("/"):
                key, _, value = text[1:].partition("=")
                cur["q"][key] = value.strip('"')
                cur["last"] = key
            elif cur["last"] is not None:
                cur["q"][cur["last"]] += " " + text.strip('"')
            else:
                cur["loc"] += text
    feats = []
    for r in raw:
        if r["type"] in _SKIP_TYPES:
            continue
        parsed = parse_location(r["loc"])
        if not parsed:
            continue
        strand, parts = parsed
        if any(a < 1 or b > seq_len or a > seq_len for a, b in parts):
            continue
        multi = len(parts) > 1
        feats.append(Feature(name=_feature_name(r["type"], r["q"]), type=r["type"],
                             start=parts[0][0], end=parts[-1][1], strand=strand,
                             parts=tuple(parts) if multi else ()))
    # A "gene" feature that duplicates a CDS/other feature span adds only clutter.
    spans = {(f.start, f.end, f.strand) for f in feats if f.type != "gene"}
    feats = [f for f in feats if f.type != "gene" or (f.start, f.end, f.strand) not in spans]
    return tuple(feats)


def _parse_entry(lines: list[str]) -> FastaRecord:
    locus = next((ln for ln in lines if ln.startswith("LOCUS")), "")
    tokens = locus.split()
    rec_id = tokens[1] if len(tokens) > 1 else "unnamed"
    topology = "circular" if "circular" in locus.lower() else "linear"

    description, mode, feat_lines, origin = "", None, [], []
    for ln in lines:
        if ln.startswith("DEFINITION"):
            description, mode = ln[10:].strip(), "def"
        elif ln.startswith("FEATURES"):
            mode = "feat"
        elif ln.startswith("ORIGIN"):
            mode = "origin"
        elif ln[:1].strip() and ln.split()[0].isupper():
            mode = None
        elif mode == "def":
            description += " " + ln.strip()
        elif mode == "feat":
            feat_lines.append(ln)
        elif mode == "origin":
            origin.append(re.sub(r"[\d\s]", "", ln))
    seq = clean_sequence("".join(origin), rec_id)
    return FastaRecord(id=rec_id, description=description.strip().rstrip("."),
                       sequence=seq, topology=topology,
                       features=_parse_features(feat_lines, len(seq)))


def parse_genbank(text: str) -> list[FastaRecord]:
    entries, cur = [], []
    for line in text.splitlines():
        if line.strip() == "//":
            if cur:
                entries.append(cur)
            cur = []
        else:
            cur.append(line)
    if any(ln.strip() for ln in cur):
        entries.append(cur)
    records = [_parse_entry(e) for e in entries if any(ln.startswith("LOCUS") for ln in e)]
    if not records:
        raise ValueError("No GenBank record found.")
    return records


def parse_any(text: str) -> list[FastaRecord]:
    """Parse FASTA or GenBank text, detected from the first characters."""
    return parse_genbank(text) if text.lstrip().startswith("LOCUS") else parse_fasta(text)
