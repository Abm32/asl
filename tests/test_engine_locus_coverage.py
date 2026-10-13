"""``asl`` must pin every locus the installed engine can call.

This is the invariant that silently broke when ``anukriti-pgx-core`` 0.11.0
widened its allele table from 18 to 33 alleles. ``asl``'s ``PLUS_STRAND`` and
``CPIC_ALLELE`` are hand-maintained, so the widening did not propagate: the
engine gained 13 loci that no panel preset could name, and that were reachable
only if a VCF happened to supply them directly. Nothing failed — the product
just quietly could not order part of the coverage it advertised.

The failure mode is one-directional and that is why it went unnoticed. A locus
``asl`` pins but the engine cannot call fails loudly at interpretation time. A
locus the engine calls but ``asl`` does not pin fails *nothing*; it simply never
appears in a panel, a report's coverage strip, or an importer column.

So the assertion below is deliberately on the direction that stays silent.
"""

from __future__ import annotations

import csv
from pathlib import Path

import anukriti_pgx_core
import pytest

from asl.evidence import _load_gene
from asl.interpret import CPIC_ALLELE, PLUS_STRAND

ALLELE_TABLE = (
    Path(anukriti_pgx_core.__file__).parent
    / "calling"
    / "data"
    / "dpyd"
    / "DPYD_alleles_anukriti_v2024.01.tsv"
)


def _engine_rsids() -> set[str]:
    """Every distinct rsID in the installed engine's DPYD allele table.

    Read from the installed package rather than a checkout, because the pin in
    ``pyproject.toml`` is what ships — a source tree sitting next to this one
    may be any version at all.
    """
    if not ALLELE_TABLE.exists():  # pragma: no cover - engine layout changed
        pytest.skip(f"engine allele table not found at {ALLELE_TABLE}")
    rsids: set[str] = set()
    with ALLELE_TABLE.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            # A compound haplotype row (HapB3) lists its loci semicolon-joined.
            rsids.update(r.strip() for r in row["rsid"].split(";") if r.strip())
    assert rsids, "engine allele table parsed to zero rsIDs"
    return rsids


def test_asl_pins_every_locus_the_engine_calls() -> None:
    missing = sorted(_engine_rsids() - set(PLUS_STRAND))
    assert not missing, (
        f"the engine calls {len(missing)} locus/loci that asl does not pin: "
        f"{missing}. They cannot be named in a panel preset, so part of the "
        f"engine's coverage is not orderable through the product. Add the "
        f"plus-strand REF/ALT from CPIC sequence_location.chromosomelocation "
        f"(never from cDNA notation, which is minus-strand for DPYD)."
    )


def test_plus_strand_and_cpic_allele_cover_the_same_loci() -> None:
    """Neither dict may gain a locus without the other.

    ``PLUS_STRAND`` without ``CPIC_ALLELE`` raises ``KeyError`` deep inside the
    CLI's off-panel and unstated-locus messages. ``CPIC_ALLELE`` without
    ``PLUS_STRAND`` puts a locus in the ``full`` panel that then refuses at
    interpretation time.
    """
    assert set(PLUS_STRAND) == set(CPIC_ALLELE)


def test_every_pinned_allele_name_has_evidence() -> None:
    """A pinned name with no CPIC evidence row would report an untraceable claim."""
    known = set(_load_gene("DPYD")["alleles"])
    orphans = sorted(
        f"{rsid} -> {name}"
        for rsid, name in CPIC_ALLELE.items()
        if name not in known
    )
    assert not orphans, (
        f"pinned allele names absent from the CPIC evidence table: {orphans}"
    )


@pytest.mark.parametrize("rsid,pair", sorted(PLUS_STRAND.items()))
def test_plus_strand_pairs_are_distinct_single_bases(
    rsid: str, pair: tuple[str, str]
) -> None:
    ref, alt = pair
    assert ref in "ACGT" and alt in "ACGT", f"{rsid}: not a single base"
    assert ref != alt, f"{rsid}: REF equals ALT"
