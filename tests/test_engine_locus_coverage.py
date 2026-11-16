"""``asl`` must pin every locus the installed engine can call.

This is the invariant that silently broke when ``anukriti-pgx-core`` 0.11.0
widened its allele table from 18 to 33 alleles. ``asl``'s ``PLUS_STRAND`` and
``CPIC_ALLELE`` are hand-maintained, so the widening did not propagate: the
engine gained 13 loci that no panel preset could name, and that were reachable
only if a VCF happened to supply them directly. Nothing failed — the product
just quietly could not order part of the coverage it advertised.

Both directions are now asserted, and the second one is here because the
assumption written in this docstring turned out to be false.

A locus the engine calls but ``asl`` does not pin fails *nothing*; it simply
never appears in a panel, a report's coverage strip, or an importer column.

A locus ``asl`` pins but the engine cannot call was assumed to "fail loudly
at interpretation time". It does **not**, except in the degenerate case
where it is the only locus assessed — see
``test_engine_calls_every_locus_asl_pins``. On any real panel the engine
drops the unknown locus silently and a homozygous no-function patient reads
Normal Metabolizer on a standard dose.

So neither direction is safe to leave untested.
"""

from __future__ import annotations

import csv
from pathlib import Path

import anukriti_pgx_core
import pytest

from asl.evidence import _load_gene
from asl.interpret import CPIC_ALLELE, EVIDENCE_GAP, GRCH38_POS, PLUS_STRAND

ALLELE_TABLE = (
    Path(anukriti_pgx_core.__file__).parent
    / "calling"
    / "data"
    / "dpyd"
    / "DPYD_alleles_anukriti_v2024.01.tsv"
)


def _engine_alts(rsid: str) -> set[str]:
    """Every ALT token the installed engine's table defines at one locus.

    Paired positionally with the rsIDs on compound rows, the same way the
    engine's own matcher pairs them.
    """
    if not ALLELE_TABLE.exists():  # pragma: no cover - engine layout changed
        pytest.skip(f"engine allele table not found at {ALLELE_TABLE}")
    alts: set[str] = set()
    with ALLELE_TABLE.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rsids = [r.strip() for r in row["rsid"].split(";") if r.strip()]
            row_alts = [a.strip() for a in row["alt"].split(";") if a.strip()]
            if len(rsids) != len(row_alts):
                continue
            for r, a in zip(rsids, row_alts):
                if r == rsid:
                    alts.add(a)
    return alts


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
        if name not in known and name not in EVIDENCE_GAP
    )
    assert not orphans, (
        f"pinned allele names absent from the CPIC evidence table: {orphans}. "
        f"Either add the CPIC evidence row, or declare the gap in "
        f"asl.interpret.EVIDENCE_GAP with the reason."
    )
    # The declared gap must stay honest in the other direction too: a name
    # listed there that HAS evidence is a stale exemption, and exemptions
    # that outlive their reason are how a gap becomes permanent.
    stale = sorted(name for name in EVIDENCE_GAP if name in known)
    assert not stale, f"EVIDENCE_GAP lists names that now have evidence: {stale}"


@pytest.mark.parametrize("rsid,pair", sorted(PLUS_STRAND.items()))
def test_plus_strand_pairs_are_distinct_single_bases(
    rsid: str, pair: tuple[str, str]
) -> None:
    ref, alt = pair
    assert ref != alt, f"{rsid}: REF equals ALT"
    assert ref and alt, f"{rsid}: empty REF or ALT"
    assert set(ref) <= set("ACGT") and set(alt) <= set("ACGT"), (
        f"{rsid}: REF/ALT contain something that is not a nucleotide"
    )
    if len(ref) == len(alt):
        # A substitution. Still exactly one base each, as before.
        assert len(ref) == 1, f"{rsid}: equal-length pair that is not a single base"
        return
    # An indel. DPYD *3 and *7 are deletions, so a single-base pair cannot
    # express them, and the engine deliberately matches them by canonical
    # token rather than by nucleotide because the VCF spelling depends on the
    # caller's normalisation. The invariant that matters is therefore not
    # "one base" but "this canonicalises to the token the engine's allele
    # table pins" — which is what is asserted here, against the engine
    # itself rather than against a literal repeated in this file.
    from anukriti_pgx_core.calling.base import canonical_indel_token

    token = canonical_indel_token(ref, alt)
    assert token is not None, f"{rsid}: unequal lengths but not a recognisable indel"
    assert token in _engine_alts(rsid), (
        f"{rsid}: spelling {ref}>{alt} canonicalises to {token}, which the "
        f"engine's allele table does not define at this locus "
        f"({sorted(_engine_alts(rsid))})"
    )


def test_every_locus_has_a_grch38_position() -> None:
    """A position per locus, as data rather than a comment.

    The generator in synthetic_genomes matches 1000 Genomes records by
    position and REF/ALT, because the GRCh37 phase 3 release carries an empty
    ID column and rsID matching silently finds nothing. asl stays the source
    of truth for where each locus is. The positions come from CPIC; 13 of the
    32 were independently corroborated against the 1000 Genomes GRCh38
    high-coverage panel, the other 19 are absent from it (CPIC only).
    """
    assert set(GRCH38_POS) == set(PLUS_STRAND)
    # DPYD spans roughly chr1:97.08-97.89 Mb on GRCh38. The window is
    # deliberately coarse: this test checks internal consistency, and CPIC is
    # the single authority for positions, so no second source (e.g. dbSNP)
    # is imported here. The rs3918290 anchor below catches a whole-table
    # assembly relabel.
    for rsid, pos in GRCH38_POS.items():
        assert isinstance(pos, int), rsid
        assert 97_070_000 < pos < 97_890_000, (rsid, pos)
    # Positions are unique: two loci at one coordinate would make the
    # position-plus-alleles match ambiguous.
    assert len(set(GRCH38_POS.values())) == len(GRCH38_POS)
    # GRCh38, not GRCh37: rs3918290 is at 97,450,058 on GRCh38 and
    # 97,915,614 on GRCh37. 13 of the 32 loci were corroborated by
    # position-and-alleles matching against the 1000 Genomes GRCh38
    # high-coverage panel; the other 19 are absent there, CPIC only.
    assert GRCH38_POS["rs3918290"] == 97_450_058


def test_engine_calls_every_locus_asl_pins() -> None:
    """The other direction, which this file previously assumed was safe.

    The module docstring above says a locus ``asl`` pins but the engine
    cannot call "fails loudly at interpretation time". **That is only true
    when it is the sole locus assessed.** Measured against the pinned
    engine, a four-locus panel carrying one unknown positive is not loud at
    all:

        rs3918290  0/0   rs55886062 0/0   rs67376798 0/0   rs75017182 0/0
        rs72549303 1/1   <- homozygous DPYD*3, an engine that lacks the row

        -> *1/*1, Normal Metabolizer, STANDARD, refusal ''

    A homozygous no-function patient reported Normal Metabolizer on a full
    fluoropyrimidine dose, with nothing flagged. The ``NOT_ASSESSED``
    refusal fires only when ``sites_assessed`` is zero, so any real panel
    masks it.

    So pinning a locus ahead of the engine is not a harmless
    forward-declaration — it is a false negative waiting for a lab to order
    it. This asserts the ordering: ``asl`` may not advertise a locus its
    own pinned engine cannot interpret.
    """
    unsupported = sorted(set(PLUS_STRAND) - _engine_rsids())
    assert not unsupported, (
        f"asl pins {len(unsupported)} locus/loci the pinned engine cannot "
        f"call: {unsupported}. A panel preset could order them, a lab could "
        f"report one positive, and the engine would silently drop it — the "
        f"patient reads Normal Metabolizer. Bump anukriti-pgx-core first, "
        f"then add the locus."
    )
