"""CPIC's mandated HapB3 tag-SNP disclosure, and the alleles that reach it.

Why this module exists
----------------------

``HapB3`` is CPIC's most frequent actionable DPYD allele in this package's target
population: 0.019658 in CPIC's Central/South Asian frequencies, which is 77.5% of
what the four-variant EMA/CPIC guideline panel can find and 60.5% of every
actionable DPYD allele CPIC publishes. It is defined by two variants:

* ``c.1129-5923C>G`` (``rs75017182``), intron 10 — the **causal** variant, which
  creates a cryptic splice site (Nie 2017, PMID 28295243);
* ``c.1236G>A`` (``rs56038477``), exon 11 — a **benign synonymous** variant
  (p.Glu412=) that most assays genotype instead, being cheaper to type.

They are not in perfect linkage disequilibrium: 14 of 6265 tag carriers in
*All of Us* WGS lacked the causal variant (LD 0.9985, 0.223%; Turner/Haidar 2024,
PMID 38129972), replicated at 1 of 46 in a Spanish cohort (Int J Mol Sci
2025;26:8136). So a tag-only positive can be a **false positive**, and the cost
of acting on one is a 50% fluoropyrimidine dose reduction in curative-intent
treatment.

CPIC therefore requires the basis to be stated in the report:

    "if only c.1236G>A is tested or results are available for, it should be
    clearly stated in the test report that 'decreased function' was inferred by
    detecting the exonic tag SNP, and disclose that in rare cases, the causal
    decreased function variant c.1129-5923C>G may not be present despite having
    this tag SNP."

Until ``pgx-core`` 0.9.0 this package **could not** comply: the engine called
HapB3 from the tag alone and did not report which locus it had matched, so there
was nothing to base the disclosure on.

The invariant these tests defend
--------------------------------

The disclosure is **additive**. It never changes the phenotype, the clinical
action or the evidence level — those stay CPIC's, unchanged, in every basis.
That is what keeps this out of the failure class of the ``*6`` incident, where a
defensible judgement was written into a pinned table as if it were CPIC's.
"""

from __future__ import annotations

import pytest

from asl.interpret import CPIC_ALLELE, PLUS_STRAND, hapb3_disclosure
from asl.models import (
    AssayResult,
    Drug,
    Regimen,
    TestStatus,
    VariantCall,
    Zygosity,
)
from asl.pipeline import build
from asl.report import render

CAUSAL = "rs75017182"
TAG = "rs56038477"
C2279 = "rs112766203"

# The four variants the EMA/CPIC guideline panel tests. Note that it reaches
# HapB3 via the TAG, which is the whole problem.
GUIDELINE_PANEL = ("rs3918290", "rs55886062", "rs67376798", TAG)


def _assay(
    found: dict[str, Zygosity],
    tested: tuple[str, ...],
    *,
    regimen: Regimen = Regimen.SYSTEMIC,
    drug: Drug = Drug.FLUOROURACIL,
    panel: str = "test panel",
) -> AssayResult:
    calls = []
    for rsid in CPIC_ALLELE:
        if rsid in found:
            status, zygosity = TestStatus.TESTED, found[rsid]
        elif rsid in tested:
            status, zygosity = TestStatus.TESTED, Zygosity.HOM_REF
        else:
            status, zygosity = TestStatus.NOT_TESTED, None
        calls.append(
            VariantCall(
                rsid=rsid,
                status=status,
                zygosity=zygosity,
                allele_label=CPIC_ALLELE[rsid],
            )
        )
    return AssayResult(
        gene="DPYD",
        panel_name=panel,
        calls=tuple(calls),
        regimen=regimen,
        population="Central/South Asian",
        drug=drug,
    )


# ---------------------------------------------------------------------------
# 1. The disclosure fires exactly when CPIC requires it, and not otherwise
# ---------------------------------------------------------------------------


def test_tag_only_emits_cpics_required_disclosure() -> None:
    """The ordinary Indian case: a PCR panel that types only the tag SNP."""
    interp = build(
        _assay(
            {TAG: Zygosity.HET},
            GUIDELINE_PANEL,
            panel="EMA/CPIC 4-variant DPYD PCR panel (tag SNP for HapB3)",
        )
    )
    assert interp.hapb3_basis == "tag_only"

    rule_ids = {f.rule_id for f in interp.flags}
    assert "P5_DPYD_HAPB3_TAG_SNP_BASIS" in rule_ids

    text = render(interp)
    # CPIC's own words must appear, not a paraphrase.
    assert "inferred by detecting the exonic tag SNP" in text
    assert "may not be present despite having this tag SNP" in text
    assert "rs75017182" in text


def test_causal_variant_observed_needs_no_disclosure() -> None:
    """Nothing is inferred, so there is nothing to qualify.

    Emitting the caveat here would be noise, and a report that caveats
    everything trains its reader to skip caveats.
    """
    interp = build(
        _assay({CAUSAL: Zygosity.HET}, (*GUIDELINE_PANEL, CAUSAL))
    )
    assert interp.hapb3_basis == "causal"
    assert not any(f.rule_id.startswith("P5_DPYD_HAPB3") for f in interp.flags)


def test_full_haplotype_needs_no_disclosure() -> None:
    interp = build(
        _assay(
            {CAUSAL: Zygosity.HET, TAG: Zygosity.HET},
            (*GUIDELINE_PANEL, CAUSAL),
        )
    )
    assert interp.hapb3_basis == "haplotype"
    assert not any(f.rule_id.startswith("P5_DPYD_HAPB3") for f in interp.flags)


def test_ld_breakdown_gets_its_own_stronger_flag() -> None:
    """Tag found, causal tested and ABSENT — the documented discordant genotype.

    This is materially different from "the causal variant was never looked at":
    here the inference the recommendation rests on is actively contradicted at
    the causal locus, so it gets a distinct rule ID and says so.
    """
    interp = build(
        _assay({TAG: Zygosity.HET}, (*GUIDELINE_PANEL, CAUSAL))
    )
    assert interp.hapb3_basis == "tag_without_causal"

    rule_ids = {f.rule_id for f in interp.flags}
    assert "P5_DPYD_HAPB3_TAG_WITHOUT_CAUSAL" in rule_ids
    assert "P5_DPYD_HAPB3_TAG_SNP_BASIS" not in rule_ids

    text = render(interp)
    assert "WAS tested" in text and "absent" in text
    assert "FALSE POSITIVE" in text


def test_no_hapb3_allele_means_no_disclosure_and_no_basis() -> None:
    interp = build(
        _assay({"rs3918290": Zygosity.HET}, GUIDELINE_PANEL)
    )
    assert interp.hapb3_basis == ""
    assert not any(f.rule_id.startswith("P5_DPYD_HAPB3") for f in interp.flags)


# ---------------------------------------------------------------------------
# 2. The disclosure is ADDITIVE. It must never move a CPIC number.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "found,tested,expected_basis",
    [
        ({TAG: Zygosity.HET}, GUIDELINE_PANEL, "tag_only"),
        ({CAUSAL: Zygosity.HET}, (*GUIDELINE_PANEL, CAUSAL), "causal"),
        (
            {CAUSAL: Zygosity.HET, TAG: Zygosity.HET},
            (*GUIDELINE_PANEL, CAUSAL),
            "haplotype",
        ),
        ({TAG: Zygosity.HET}, (*GUIDELINE_PANEL, CAUSAL), "tag_without_causal"),
    ],
    ids=["tag_only", "causal", "haplotype", "tag_without_causal"],
)
def test_the_cpic_call_is_identical_whatever_the_basis(
    found, tested, expected_basis
) -> None:
    """One decreased-function chromosome is CPIC's Intermediate Metabolizer on a
    50% starting dose, however the laboratory came to know about it.

    If this fails, ``asl`` has started reclassifying rather than disclosing.
    """
    interp = build(_assay(found, tested))
    assert interp.hapb3_basis == expected_basis
    assert interp.phenotype == "Intermediate Metabolizer"
    assert interp.cpic_recommendation_action == "REDUCE_50PCT"


def test_disclosure_helper_returns_none_only_where_intended() -> None:
    assert hapb3_disclosure("") is None
    assert hapb3_disclosure("causal") is None
    assert hapb3_disclosure("haplotype") is None
    for basis in ("tag_only", "tag_without_causal"):
        result = hapb3_disclosure(basis)
        assert result is not None
        summary, detail = result
        assert summary and detail
        # CPIC's requirement is quoted in both.
        assert "exonic tag SNP" in detail


def test_no_basis_produces_an_all_clear() -> None:
    """The package-wide invariant, re-checked on this axis.

    A tag-only call is the one most likely to be a false positive, and it would
    be the most tempting place to soften the language into reassurance.
    """
    for basis_found, tested in [
        ({TAG: Zygosity.HET}, GUIDELINE_PANEL),
        ({TAG: Zygosity.HET}, (*GUIDELINE_PANEL, CAUSAL)),
        ({CAUSAL: Zygosity.HET}, (*GUIDELINE_PANEL, CAUSAL)),
    ]:
        text = render(build(_assay(basis_found, tested))).lower()
        for phrase in ("no risk", "all clear", "safe to proceed", "normal dose"):
            assert phrase not in text


# ---------------------------------------------------------------------------
# 3. The alleles that were previously uncallable
# ---------------------------------------------------------------------------


def test_c2279ct_is_interpretable() -> None:
    """``c.2279C>T`` is one of the three intermediate-metabolizer variants
    actually found in D-TORCH, the only Indian NGS DPYD study (Baskarane 2026,
    Front Pharmacol 17:1732128).

    Before ``pgx-core`` 0.9.0 the engine returned ``Indeterminate`` for it and
    ``asl`` refused the interpretation outright — for a real patient in the only
    Indian cohort that has looked. It is 18.6% of actionable Central/South Asian
    DPYD allele frequency, the second most common after HapB3.
    """
    interp = build(
        _assay({C2279: Zygosity.HET}, (*GUIDELINE_PANEL, CAUSAL, C2279))
    )
    assert interp.diplotype == "*1/c.2279C>T"
    assert interp.phenotype == "Intermediate Metabolizer"
    assert interp.cpic_recommendation_action == "REDUCE_50PCT"
    # It must carry provenance like any other found allele, not appear in the
    # diplotype alone.
    assert any(e.allele_label == "c.2279C>T" for e in interp.evidence)
    # No HapB3 basis: this allele has nothing to do with the tag SNP.
    assert interp.hapb3_basis == ""


def test_every_locus_has_a_pinned_strand_pair_and_a_cpic_label() -> None:
    """The two tables must not drift apart.

    A locus in ``CPIC_ALLELE`` but not ``PLUS_STRAND`` raises at interpretation
    time; the reverse is a locus that can never be labelled and so can never be
    reported.
    """
    assert set(PLUS_STRAND) == set(CPIC_ALLELE)


def test_the_two_hapb3_loci_map_to_different_cpic_alleles() -> None:
    """If both mapped to the same CPIC allele name the distinction would be
    cosmetic — the evidence section would show one row for two different
    findings."""
    assert CPIC_ALLELE[CAUSAL] != CPIC_ALLELE[TAG]
    assert CPIC_ALLELE[CAUSAL] == "c.1129-5923C>G"
    assert "HapB3" in CPIC_ALLELE[TAG]


# ---------------------------------------------------------------------------
# 4. The assay inventory must not contradict itself
# ---------------------------------------------------------------------------


def test_assay_inventory_does_not_claim_the_causal_variant_was_tested() -> None:
    """Found by walking a real patient.

    CPIC's allele name for the tag SNP is the compound
    ``"c.1129-5923C>G, c.1236G>A (HapB3)"`` — it *names the causal variant*.
    Printing that as what the assay tested, in a report that also lists
    ``c.1129-5923C>G`` under "absence is NOT established", tells the reader the
    same variant was both found and never looked at.

    The assay inventory therefore names the locus actually read, while evidence
    lookup keeps CPIC's compound name. Same distinction as the disclosure,
    applied to the inventory.
    """
    interp = build(
        _assay(
            {TAG: Zygosity.HET},
            GUIDELINE_PANEL,
            panel="EMA/CPIC 4-variant DPYD PCR panel",
        )
    )
    text = render(interp)

    inventory = text.split("[ASL-ASSAY-01]", 1)[1].split("[ASL-", 1)[0]
    found_line = next(
        ln for ln in inventory.splitlines() if "Tested & FOUND" in ln
    )
    # The tag is named; the causal variant is NOT claimed as tested.
    assert "c.1236G>A" in found_line
    assert "c.1129-5923C>G" not in found_line

    # And the causal variant is still correctly listed as not assayed.
    assert "Absence is NOT established for" in text
    untested = text.split("Absence is NOT established for", 1)[1]
    assert "c.1129-5923C>G" in untested


def test_causal_variant_tested_appears_in_the_inventory() -> None:
    """The converse: when the causal locus IS on the panel, say so."""
    interp = build(
        _assay(
            {CAUSAL: Zygosity.HET},
            (*GUIDELINE_PANEL, CAUSAL),
            panel="NGS DPYD panel A",  # deliberately does not echo a variant name
        )
    )
    text = render(interp)
    inventory = text.split("[ASL-ASSAY-01]", 1)[1].split("[ASL-", 1)[0]
    assert "c.1129-5923C>G" in inventory
    # It must not then also be listed among the loci that were never assayed.
    untested_labels = {label for label, _ in interp.uninformative_loci}
    assert "c.1129-5923C>G" not in untested_labels
    not_tested_rsids = {c.rsid for c in interp.assay.calls if c.rsid not in
                        {t.rsid for t in interp.assay.tested}}
    assert CAUSAL not in not_tested_rsids
