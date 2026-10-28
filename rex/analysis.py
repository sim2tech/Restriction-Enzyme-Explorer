"""Restriction-site analysis backed by Biopython's Bio.Restriction (REBASE data)."""
from __future__ import annotations

from Bio import Restriction
from Bio.Seq import Seq

from .models import CutSite, EnzymeInfo
from .motif import reverse_complement

# A handful of enzymes most people recognise from cloning work.
COMMON_ENZYMES = [
    "EcoRI", "BamHI", "HindIII", "XhoI", "NotI", "SalI", "PstI", "KpnI",
    "SacI", "SmaI", "NcoI", "NdeI", "XbaI", "SpeI", "EcoRV", "BglII",
]


def list_enzymes(commercial_only: bool = True) -> list[str]:
    """Names of enzymes in Biopython's database (optionally commercially available only)."""
    batch = Restriction.CommOnly if commercial_only else Restriction.AllEnzymes
    return sorted(e.__name__ for e in batch)


def get_enzyme(name: str):
    """Look an enzyme up by name in Biopython's database."""
    try:
        return Restriction.AllEnzymes.get(name)
    except (ValueError, KeyError, AttributeError) as exc:
        raise ValueError(f"Unknown enzyme: {name}") from exc


def _isoschizomers(enz, name: str) -> tuple:
    try:
        return tuple(sorted(e.__name__ for e in enz.isoschizomers() if e.__name__ != name))
    except Exception:
        return ()


def describe_enzyme(name: str) -> EnzymeInfo:
    """Recognition site, cleavage pattern and documented properties of one enzyme."""
    enz = get_enzyme(name)
    if enz.is_blunt():
        overhang = "blunt"
    elif enz.is_5overhang():
        overhang = "5' overhang"
    elif enz.is_3overhang():
        overhang = "3' overhang"
    else:
        overhang = "unknown"
    try:
        pattern = enz.elucidate()
    except Exception:  # undefined cut positions for some enzymes
        pattern = "n/a"
    return EnzymeInfo(
        name=name,
        site=enz.site,
        cut_pattern=pattern,
        size=enz.size,
        overhang=overhang,
        overhang_seq=str(getattr(enz, "ovhgseq", "") or ""),
        suppliers=tuple(getattr(enz, "suppl", ()) or ()),
        commercial=bool(enz.is_comm()),
        opt_temp=getattr(enz, "opt_temp", None),
        inact_temp=getattr(enz, "inact_temp", None),
        isoschizomers=_isoschizomers(enz, name),
    )


def bottom_shift(pattern: str, site: str) -> int | None:
    """Offset of the bottom-strand cut from the top-strand cut, e.g. EcoRI G^AATT_C -> +4.

    Only stated for palindromic sites, where the pattern reads the same on both
    strands. Returns None when the pattern cannot be read reliably.
    """
    try:
        if reverse_complement(site) != site or "^" not in pattern or "_" not in pattern:
            return None if "^" in pattern and "_" in pattern else (0 if "^" in pattern else None)
        top = len([c for c in pattern[:pattern.index("^")] if c.isalpha()])
        bottom = len([c for c in pattern[:pattern.index("_")] if c.isalpha()])
        return bottom - top
    except Exception:
        return None


def find_sites(sequence: str, enzyme_names: list[str], start: int = 1,
               end: int | None = None, circular: bool = False) -> list[CutSite]:
    """Search seq[start..end] (1-based, inclusive) with the chosen enzymes.

    Returned positions are in coordinates of the *full* sequence. Only the
    extracted region is scanned, so sites that straddle the region boundary
    are not reported. ``circular`` is only honoured for the whole sequence.
    """
    end = len(sequence) if end is None else end
    if not 1 <= start <= end <= len(sequence):
        raise ValueError("Region must satisfy 1 <= start <= end <= sequence length.")
    region = Seq(sequence[start - 1:end])
    whole = start == 1 and end == len(sequence)
    linear = not (circular and whole)

    sites: list[CutSite] = []
    for name in enzyme_names:
        enz = get_enzyme(name)
        try:
            pattern = enz.elucidate()
        except Exception:
            pattern = "n/a"
        shift = bottom_shift(pattern, enz.site)
        for pos in enz.search(region, linear=linear):
            top = pos + start - 1
            sites.append(CutSite(name, top, enz.site, pattern,
                                 None if shift is None else top + shift))
    return sorted(sites, key=lambda s: (s.position, s.enzyme))
