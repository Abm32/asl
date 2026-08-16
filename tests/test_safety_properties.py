"""The safety properties. These are the assertions that make asl trustworthy.

Each test here corresponds to a stated refusal in the architecture document. If
one of these fails, the package is not safe to put in front of a clinician,
regardless of whether everything else passes.
"""

from __future__ import annotations

import itertools

import pytest

from asl import (
    AssayResult,
    Regimen,
    TestStatus,
    VariantCall,
    Zygosity,
    build,
    render,
)
from asl.interpret import CPIC_ACTIONABLE, CPIC_ALLELE, InterpretationRefused

ALL_RSIDS = tuple(CPIC_ALLELE.keys())


def _tested(rsid: str, zyg: Zygosity) -> VariantCall:
    return VariantCall(
        rsid=rsid, status=TestStatus.TESTED, zygosity=zyg,
        allele_label=CPIC_ALLELE[rsid],
    )


def _untested(rsid: str) -> VariantCall:
    return VariantCall(
        rsid=rsid, status=TestStatus.NOT_TESTED, allele_label=CPIC_ALLELE[rsid]
    )


def _assay(genotypes: dict[str, Zygosity], **kw) -> AssayResult:
    tested = [_tested(r, genotypes.get(r, Zygosity.HOM_REF)) for r in CPIC_ACTIONABLE]
    untested = [_untested(r) for r in ALL_RSIDS if r not in CPIC_ACTIONABLE]
    kw.setdefault("regimen", Regimen.SYSTEMIC)
    kw.setdefault("population", "Central/South Asian")
    return AssayResult(
        gene="DPYD",
        panel_name="CPIC 4-variant DPYD panel",
        calls=tuple(tested + untested),
        **kw,
    )


# ---------------------------------------------------------------------------
# Refusal 4: no code path emits an all-clear.
# ---------------------------------------------------------------------------

FORBIDDEN = (
    "no risk", "low risk", "not at risk", "risk-free", "all clear",
    "all-clear", "normal risk", "no increased risk", "safe to proceed",
    "rules out", "ruled out", "excludes risk",
)

# Phrases in which a forbidden substring appears *negated*, and is therefore
# safe. Kept as an explicit allowlist rather than a looser regex so that adding
# a reassuring sentence to the report cannot accidentally slip past the check —
# a new phrasing has to be added here deliberately.
NEGATED_CONTEXTS = (
    "this is not an all-clear",
    "does not fully rule out",
    "not an all-clear",
)


def _scrub(text: str) -> str:
    """Remove known negated contexts before scanning for forbidden phrases."""
    out = text.lower()
    for phrase in NEGATED_CONTEXTS:
        out = out.replace(phrase, " ")
    return out


@pytest.mark.parametrize(
    "genotypes",
    [
        {},  # panel-negative — the case most likely to be read as reassurance
        {"rs3918290": Zygosity.HET},
        {"rs3918290": Zygosity.HOM_VAR},
        {"rs55886062": Zygosity.HOM_VAR},
        {"rs67376798": Zygosity.HET},
        {"rs56038477": Zygosity.HOM_VAR},
    ],
    ids=["panel_negative", "2A_het", "2A_hom", "13_hom", "c2846_het", "hapb3_hom"],
)
@pytest.mark.parametrize(
    "regimen", [Regimen.SYSTEMIC, Regimen.CHEMORADIATION, Regimen.UNKNOWN]
)
def test_no_input_produces_an_all_clear(genotypes, regimen) -> None:
    text = _scrub(render(build(_assay(genotypes, regimen=regimen))))
    for phrase in FORBIDDEN:
        assert phrase not in text, (
            f"report contains {phrase!r}. asl must never state the absence of "
            f"risk; a panel-negative result is a residual-risk statement."
        )


def test_exhaustive_actionable_combinations_never_emit_an_all_clear() -> None:
    """Every zygosity combination across the actionable panel."""
    zygs = (Zygosity.HOM_REF, Zygosity.HET, Zygosity.HOM_VAR)
    checked = 0
    for combo in itertools.product(zygs, repeat=len(CPIC_ACTIONABLE)):
        genotypes = dict(zip(CPIC_ACTIONABLE, combo))
        if sum({Zygosity.HOM_REF: 0, Zygosity.HET: 1, Zygosity.HOM_VAR: 2}[z]
               for z in combo) > 2:
            continue  # over-specified; refused separately
        text = _scrub(render(build(_assay(genotypes))))
        for phrase in FORBIDDEN:
            assert phrase not in text, f"{genotypes}: contains {phrase!r}"
        checked += 1
    assert checked >= 10, f"only {checked} combinations were reachable"


def test_every_report_carries_the_residual_risk_section() -> None:
    for genotypes in ({}, {"rs3918290": Zygosity.HET}):
        text = render(build(_assay(genotypes)))
        assert "RESIDUAL RISK — THIS IS NOT AN ALL-CLEAR" in text
        assert "PMID 29152729" in text  # the CPIC caveat is always cited


# ---------------------------------------------------------------------------
# Refusal 5: NOT_TESTED is never reported as absent.
# ---------------------------------------------------------------------------


def test_untested_loci_are_named_as_untested() -> None:
    text = render(build(_assay({})))
    assert "Absence is NOT established for" in text
    assert "An untested locus is not a negative result." in text
    # *6 is not on the CPIC four and must be named as untested.
    assert "c.2194G>A (*6)" in text


def test_a_tested_locus_requires_a_zygosity() -> None:
    with pytest.raises(ValueError, match="must record a zygosity"):
        VariantCall(rsid="rs3918290", status=TestStatus.TESTED)


def test_an_untested_locus_cannot_carry_a_zygosity() -> None:
    with pytest.raises(ValueError, match="cannot have a zygosity"):
        VariantCall(
            rsid="rs3918290",
            status=TestStatus.NOT_TESTED,
            zygosity=Zygosity.HOM_REF,
        )


def test_copies_is_undefined_for_an_untested_locus() -> None:
    call = _untested("rs1801160")
    with pytest.raises(ValueError, match="undefined for an untested locus"):
        _ = call.copies


# ---------------------------------------------------------------------------
# Refusal 6: no toxicity rate without a regimen.
# ---------------------------------------------------------------------------


def test_unknown_regimen_withholds_the_rate() -> None:
    """With no regimen, no rate may be attributed to this patient.

    The refusal text does mention both published rates, deliberately — it has to
    explain *why* it is refusing. What must not happen is a rate presented as
    applying to this patient, so the assertion targets the attribution line
    rather than the digits.
    """
    interp = build(_assay({}, regimen=Regimen.UNKNOWN))
    assert interp.residual.quoted_rate is None
    assert interp.residual.refusal_reason is not None

    text = render(interp)
    assert "RATE WITHHELD" in text
    assert "Regimen-matched rate" not in text
    assert "of patients without an actionable DPYD variant experienced" not in text


def test_systemic_and_chemoradiation_quote_different_rates() -> None:
    systemic = render(build(_assay({}, regimen=Regimen.SYSTEMIC)))
    chemorad = render(build(_assay({}, regimen=Regimen.CHEMORADIATION)))
    assert "22.7% (231/1018)" in systemic and "13.6%" not in systemic
    assert "13.6% (105/771)" in chemorad and "22.7%" not in chemorad


def test_the_attributable_fraction_is_never_quoted_as_sensitivity() -> None:
    """20-30% is an attributable fraction for early-onset events, not sensitivity.

    Conflating them overstates the panel roughly threefold, so the figure must
    not appear outside CPIC's own quoted context.
    """
    text = render(build(_assay({})))
    assert "3.5-21.6%" in text
    assert "20-30%" not in text
    assert "20–30%" not in text


# ---------------------------------------------------------------------------
# Refusal 1 and 3: no dose advice, no override of CPIC.
# ---------------------------------------------------------------------------


def test_report_states_it_does_not_recommend_doses() -> None:
    text = render(build(_assay({})))
    assert "asl does not recommend doses" in text
    assert "It does not recommend a dose" in text


def test_cpic_call_is_attributed_and_unmodified() -> None:
    interp = build(_assay({"rs3918290": Zygosity.HOM_VAR}))
    assert interp.phenotype == "Poor Metabolizer"
    assert interp.cpic_recommendation_action == "AVOID"
    text = render(interp)
    assert "CPIC CALL (unmodified, attributed to CPIC)" in text


def test_conflict_flag_states_the_cpic_call_stands() -> None:
    """*6 carries a conflict record; the flag must not reclassify it."""
    calls = [_tested(r, Zygosity.HOM_REF) for r in CPIC_ACTIONABLE]
    calls.append(_tested("rs1801160", Zygosity.HET))
    calls += [
        _untested(r) for r in ALL_RSIDS
        if r not in CPIC_ACTIONABLE and r != "rs1801160"
    ]
    assay = AssayResult(
        gene="DPYD", panel_name="panel + *6", calls=tuple(calls),
        regimen=Regimen.SYSTEMIC, population="Central/South Asian",
    )
    interp = build(assay)
    assert interp.phenotype == "Normal Metabolizer"  # unchanged by the flag
    text = render(interp)
    assert "P2_DPYD_STAR6_INVITRO_CLINICAL_CONFLICT" in text
    assert "The CPIC call" in text and "stands" in text
    # Both directions must be present.
    assert "Del Re" in text and "Baskarane" in text


# ---------------------------------------------------------------------------
# Refusal 10 and the over-specification guard.
# ---------------------------------------------------------------------------


def test_overspecified_genotype_is_refused() -> None:
    calls = [
        _tested("rs1801265", Zygosity.HOM_VAR),
        _tested("rs1801160", Zygosity.HET),
        _tested("rs2297595", Zygosity.HET),
    ]
    calls += [_tested(r, Zygosity.HOM_REF) for r in CPIC_ACTIONABLE]
    assay = AssayResult(
        gene="DPYD", panel_name="over-specified input", calls=tuple(calls),
        regimen=Regimen.SYSTEMIC, population="Central/South Asian",
    )
    with pytest.raises(InterpretationRefused, match="copies were detected"):
        build(assay)


# ---------------------------------------------------------------------------
# Assay-content refusal.
# ---------------------------------------------------------------------------


def test_unstated_panel_is_refused() -> None:
    assay = AssayResult(
        gene="DPYD", panel_name="   ",
        calls=(_tested("rs3918290", Zygosity.HOM_REF),),
    )
    with pytest.raises(InterpretationRefused, match="content is unstated"):
        build(assay)


def test_empty_assay_is_refused() -> None:
    assay = AssayResult(gene="DPYD", panel_name="named but empty", calls=())
    with pytest.raises(InterpretationRefused, match="no recorded loci"):
        build(assay)


# ---------------------------------------------------------------------------
# Population handling: omit, never substitute.
# ---------------------------------------------------------------------------


def test_star13_is_flagged_uninformative_in_south_asians() -> None:
    text = render(build(_assay({})))
    assert "LOCI THAT CANNOT RETURN A POSITIVE HERE" in text
    assert "c.1679T>G (*13)" in text
    assert "0.00000" in text


def test_no_population_means_no_frequency_claim() -> None:
    interp = build(_assay({"rs3918290": Zygosity.HET}, population=None))
    for rec in interp.evidence:
        assert rec.frequency_in_population is None
    assert "LOCI THAT CANNOT RETURN A POSITIVE HERE" not in render(interp)


def test_engine_version_is_recorded() -> None:
    interp = build(_assay({}))
    assert interp.engine_version.startswith("0.7.")
    assert interp.engine_version in render(interp)
