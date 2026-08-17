"""The safety properties. These are the assertions that make asl trustworthy.

Each test here corresponds to a stated refusal in the architecture document. If
one of these fails, the package is not safe to put in front of a clinician,
regardless of whether everything else passes.
"""

from __future__ import annotations

import itertools
import re

import pytest

from asl import (
    AssayResult,
    Drug,
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
    "regimen",
    [
        Regimen.SYSTEMIC,
        Regimen.CHEMORADIATION,
        Regimen.INDUCTION_TPF,
        Regimen.UNKNOWN,
    ],
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


# ---------------------------------------------------------------------------
# Refusal 6b: a regimen with no published comparator gets no borrowed rate.
#
# This is the Indian majority case. Induction TPF for head and neck cancer had
# no bucket of its own, so it mapped to SYSTEMIC and a patient on 5-FU
# 750 mg/m2/day D1-D5 plus docetaxel and cisplatin was shown 22.7% — measured in
# a Dutch, largely colorectal, mostly capecitabine cohort. Observed grade >=3
# toxicity on TPF in Indian routine practice is roughly threefold that
# (Patil 2016, J Community Support Oncol 14(10):412-419, n=58).
# ---------------------------------------------------------------------------


def test_induction_tpf_never_borrows_the_systemic_rate() -> None:
    """The regimen axis of refusal 11: no other group's rate, in either form."""
    interp = build(_assay({}, regimen=Regimen.INDUCTION_TPF))
    assert interp.residual.quoted_rate is None
    assert interp.residual.refusal_reason is not None

    text = render(interp)
    assert "RATE WITHHELD" in text
    assert "Regimen-matched rate" not in text
    assert "of patients without an actionable DPYD variant experienced" not in text


def test_induction_tpf_refusal_states_why_the_comparator_is_absent() -> None:
    """A refusal has to be actionable: it must name what is missing and why."""
    reason = build(_assay({}, regimen=Regimen.INDUCTION_TPF)).residual.refusal_reason
    assert reason is not None
    # Names the exposure that makes the published cohorts inapplicable.
    assert "750 mg/m2/day" in reason
    # Cites the observed Indian rates, and both source studies it is departing from.
    assert "PMID 30348537" in reason and "PMID 30361102" in reason
    assert "67.2%" in reason  # observed mucositis on TPF
    # And states that the observed figures are NOT a substitute comparator,
    # so this fix cannot be read as swapping one misattributed rate for another.
    assert "WHOLE-COHORT" in reason
    assert "not a substitute comparator" in reason


def test_every_regimen_either_quotes_a_matched_rate_or_refuses() -> None:
    """No regimen may fall through to a default rate.

    Asserted across the whole enum so that adding a member without adding either
    a rate or a refusal fails here rather than silently inheriting one.
    """
    for regimen in Regimen:
        rr = build(_assay({}, regimen=regimen)).residual
        if rr.quoted_rate is None:
            assert rr.refusal_reason, f"{regimen}: no rate and no refusal"
        else:
            assert rr.rate_source, f"{regimen}: rate without a source"
            assert regimen in (Regimen.SYSTEMIC, Regimen.CHEMORADIATION), (
                f"{regimen}: quoted a rate that was measured in another setting"
            )


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


def _assay_including(genotypes: dict[str, Zygosity], **kw) -> AssayResult:
    """An assay that TESTS every locus named in ``genotypes``, plus the CPIC four.

    The module-level ``_assay`` helper only tests ``CPIC_ACTIONABLE`` and marks
    everything else NOT_TESTED, which is right for the panel-negative properties
    but silently drops a genotype supplied for an off-panel locus.
    """
    tested_ids = tuple(dict.fromkeys(CPIC_ACTIONABLE + tuple(genotypes)))
    tested = [_tested(r, genotypes.get(r, Zygosity.HOM_REF)) for r in tested_ids]
    untested = [_untested(r) for r in ALL_RSIDS if r not in tested_ids]
    kw.setdefault("regimen", Regimen.SYSTEMIC)
    kw.setdefault("population", "Central/South Asian")
    return AssayResult(
        gene="DPYD",
        panel_name="CPIC four plus sequenced South Asian loci",
        calls=tuple(tested + untested),
        **kw,
    )


def test_report_lists_every_detected_allele_with_copy_counts() -> None:
    """pgx-core 0.7.3 §8 item 1: the diplotype must be auditable.

    A two-allele diplotype string asks the reader to trust that those two
    alleles are everything the engine found. Printing the per-allele copy counts
    lets them verify it instead. The over-specified case is refused upstream, so
    this line makes the remaining case checkable rather than merely correct.
    """
    interp = build(
        _assay_including({"rs1801265": Zygosity.HET, "rs1801160": Zygosity.HET})
    )
    assert dict(interp.allele_counts) == {"*6": 1, "*9A": 1}
    assert interp.copies_detected == 2

    flat = " ".join(render(interp).split())
    assert "Detected :" in flat
    assert "*6x1" in flat and "*9Ax1" in flat
    assert "2 variant copies total" in flat


def test_homozygous_call_reports_two_copies_of_one_allele() -> None:
    """A homozygote must not look like a single-copy finding.

    *9A is the most common DPYD variant in Central/South Asians (0.255), so
    homozygotes are expected rather than exotic, and the copy count is the only
    thing in the report that distinguishes *9A/*9A from a het.
    """
    interp = build(_assay_including({"rs1801265": Zygosity.HOM_VAR}))
    assert dict(interp.allele_counts) == {"*9A": 2}
    assert interp.copies_detected == 2
    flat = " ".join(render(interp).split())
    assert "*9Ax2" in flat


def test_zero_frequency_is_low_yield_not_impossibility() -> None:
    """A 0.00000 reference frequency must not be reported as "cannot occur".

    Chan 2024 (PMID 38886557) found *13 in a Tunisian patient with severe
    toxicity despite a 0% Middle Eastern reference frequency, and *2A in Chinese,
    Thai and Japanese patients despite a 0% East Asian reference frequency. A zero
    in a reference panel is an absence of observation, not proof of absence.

    Telling a clinician the variant "cannot" be found would license ignoring a
    positive call — the same conversion of missing data into certainty this
    package prevents, running in the opposite direction.
    """
    text = render(build(_assay({})))
    # Report text is wrapped, so compare against a whitespace-normalised copy.
    flat = " ".join(text.split())
    assert "cannot return a positive" not in flat
    assert "cannot occur" not in flat
    # It must still say the yield is low, and say why the zero is weak evidence.
    assert "very low expected yield" in flat.lower()
    assert "not proof the variant is absent" in flat
    assert "must still be acted on" in flat


def test_star13_is_flagged_uninformative_in_south_asians() -> None:
    text = render(build(_assay({})))
    assert "LOCI WITH VERY LOW EXPECTED YIELD HERE" in text
    assert "c.1679T>G (*13)" in text
    assert "0.00000" in text


def test_no_population_means_no_frequency_claim() -> None:
    interp = build(_assay({"rs3918290": Zygosity.HET}, population=None))
    for rec in interp.evidence:
        assert rec.frequency_in_population is None
    assert "LOCI WITH VERY LOW EXPECTED YIELD HERE" not in render(interp)


def test_engine_version_is_recorded() -> None:
    """The report must name the engine that produced the call, whatever it is.

    Asserting a specific version here would make this a chore to edit at every
    engine bump, and the point is not *which* version — it is that a call can be
    attributed to one. `pgx-core` 0.8.0 shipped with its own version string out of
    sync for exactly that reason.
    """
    interp = build(_assay({}))
    assert re.fullmatch(r"\d+\.\d+\.\d+", interp.engine_version)
    assert interp.engine_version in render(interp)


def test_an_unknown_population_key_is_refused() -> None:
    """A typo must not silently degrade into a population-blind report.

    ``"South Asian"`` is not a CPIC key — the key is ``"Central/South Asian"``.
    Before this guard the typo suppressed every frequency line and the
    uninformative-locus section while the report still printed the population in
    its header, so it looked population-aware and was not.
    """
    with pytest.raises(InterpretationRefused, match="not a CPIC population key"):
        build(_assay({}, population="South Asian"))

    # The correct key still works, and still produces frequency claims.
    interp = build(_assay({"rs3918290": Zygosity.HET}))
    assert any(r.frequency_in_population is not None for r in interp.evidence)


# ---------------------------------------------------------------------------
# Preconditions that would otherwise yield a confidently wrong report.
# ---------------------------------------------------------------------------


def test_an_assay_with_nothing_tested_is_refused() -> None:
    """Zero observations must not become "Normal Metabolizer".

    The engine is given no variants and correctly returns ``*1/*1``. Rendering
    that is the exact conversion of "we did not look" into "this patient is
    normal" that the package exists to prevent, so the pipeline refuses before
    the engine is called.
    """
    assay = AssayResult(
        gene="DPYD",
        panel_name="requisition raised, sample never run",
        calls=tuple(_untested(r) for r in ALL_RSIDS),
        regimen=Regimen.SYSTEMIC,
        population="Central/South Asian",
    )
    with pytest.raises(InterpretationRefused, match="NOT_TESTED"):
        build(assay)


def test_another_genes_assay_is_refused() -> None:
    """DPYD loci are the only thing the DPYD caller may be given."""
    assay = AssayResult(
        gene="CYP2D6",
        panel_name="CYP2D6 panel routed to the wrong interpreter",
        calls=(_tested("rs1801265", Zygosity.HET),),
        regimen=Regimen.SYSTEMIC,
        population="Central/South Asian",
    )
    with pytest.raises(InterpretationRefused, match="interprets DPYD only"):
        build(assay)


def test_a_found_allele_without_a_cpic_label_is_refused() -> None:
    """An allele that reaches the diplotype must also reach the evidence section.

    Without a label there is no provenance record, no frequency and no conflict
    flag, so the allele would appear in the CPIC call and nowhere else — a
    finding with its evidence stripped off.
    """
    calls = [
        VariantCall(
            rsid="rs1801265", status=TestStatus.TESTED, zygosity=Zygosity.HOM_VAR,
            allele_label=None,
        )
    ]
    calls += [_tested(r, Zygosity.HOM_REF) for r in CPIC_ACTIONABLE]
    assay = AssayResult(
        gene="DPYD", panel_name="unlabelled result", calls=tuple(calls),
        regimen=Regimen.SYSTEMIC, population="Central/South Asian",
    )
    with pytest.raises(InterpretationRefused, match="no CPIC allele label"):
        build(assay)


def test_a_duplicated_locus_is_refused() -> None:
    """Two readings of one position must be resolved upstream, not by us.

    The engine is keyed by rsID, so the second observation overwrote the first
    and the losing genotype vanished without appearing anywhere in the report.
    """
    with pytest.raises(ValueError, match="appears twice"):
        AssayResult(
            gene="DPYD",
            panel_name="duplicated locus",
            calls=(
                _tested("rs1801265", Zygosity.HET),
                _tested("rs1801265", Zygosity.HOM_VAR),
            ),
            regimen=Regimen.SYSTEMIC,
        )

    # A TESTED/NOT_TESTED pair for the same locus is the same defect.
    with pytest.raises(ValueError, match="appears twice"):
        AssayResult(
            gene="DPYD",
            panel_name="tested and untested at once",
            calls=(_tested("rs1801265", Zygosity.HET), _untested("rs1801265")),
            regimen=Regimen.SYSTEMIC,
        )


def test_an_offpanel_genotype_is_refused_rather_than_dropped() -> None:
    """A spreadsheet column outside the chosen panel must not be discarded.

    ``cpic4`` does not cover ``rs1801160`` (*6). Before this guard a CSV row
    recording ``*6`` heterozygous, interpreted against ``cpic4``, marked the
    locus NOT_TESTED and reported ``*1/*1`` — the observed variant disappeared
    entirely.
    """
    from asl.cli import _all_negative
    from asl.cli import _assay as cli_assay

    with pytest.raises(InterpretationRefused, match="does not cover"):
        cli_assay(
            "cpic4",
            _all_negative("cpic4") | {"rs1801160": Zygosity.HET},
            Regimen.SYSTEMIC,
            "Central/South Asian",
        )

    # The same genotype against a panel that does cover the locus is reported.
    interp = build(
        cli_assay(
            "cpic4_plus_sas",
            _all_negative("cpic4_plus_sas") | {"rs1801160": Zygosity.HET},
            Regimen.SYSTEMIC,
            "Central/South Asian",
        )
    )
    assert "*6" in interp.diplotype
    assert any("*6" in rec.allele_label for rec in interp.evidence)


def test_a_blank_csv_cell_is_a_no_call_not_a_negative() -> None:
    """The defect this guard closes: blank cell reported as "tested & negative".

    A spreadsheet where the *9A and *6 columns exist but are empty means the
    assay covered those loci and returned no call. Before this fix an absent or
    blank genotype defaulted to ``hom_ref``, so those loci appeared under
    "Tested & negative" and the patient was told absence had been established at
    positions nobody read.
    """
    from asl.cli import _assay as cli_assay
    from asl.cli import _row_genotypes

    row = {
        "rs3918290": "0", "rs55886062": "0", "rs67376798": "0", "rs56038477": "0",
        "rs1801265": "", "rs1801160": "", "rs1801159": "0", "rs2297595": "0",
    }
    genotypes = _row_genotypes(row)
    assert genotypes["rs1801265"] is None
    assert genotypes["rs1801160"] is None

    assay = cli_assay(
        "cpic4_plus_sas", genotypes, Regimen.SYSTEMIC, "Central/South Asian"
    )
    no_call = {"rs1801265", "rs1801160"}
    assert no_call.isdisjoint({c.rsid for c in assay.tested})
    assert no_call <= {c.rsid for c in assay.not_tested}

    # And the report says so, rather than listing them as negative.
    text = render(build(assay))
    assert "Absence is NOT established for" in text
    assert "c.85T>C (*9A)" in text.split("Absence is NOT established for")[1]
    # The panel name must not claim coverage the data did not deliver.
    assert "6 of 8 loci returned a call" in assay.panel_name


def test_a_panel_locus_with_no_column_at_all_is_refused() -> None:
    """An absent column is not a blank cell, and neither is a negative.

    A CSV that simply lacks the *6 column while ``cpic4_plus_sas`` is selected is
    a mismatch between the data and the claimed panel, so it is refused rather
    than filled in.
    """
    from asl.cli import _assay as cli_assay

    partial = _all_negative_except("cpic4_plus_sas", "rs1801160")
    with pytest.raises(InterpretationRefused, match="no entry for"):
        cli_assay("cpic4_plus_sas", partial, Regimen.SYSTEMIC, "Central/South Asian")


def _all_negative_except(panel: str, drop: str) -> dict[str, Zygosity | None]:
    from asl.cli import _all_negative

    genotypes = _all_negative(panel)
    genotypes.pop(drop)
    return genotypes


def test_quadrant_excludes_no_call_patients_rather_than_calling_them_negative(
    tmp_path,
) -> None:
    """The panel-negative/toxic cell must not be inflated by unread loci.

    That cell is the entire point of a site retrospective, and counting a no-call
    as a negative would bias it in the direction that flatters our own thesis.
    """
    import csv as _csv

    from asl.cli import cmd_quadrant

    path = tmp_path / "cohort.csv"
    cols = [
        "anonymised_id", "rs3918290", "rs55886062", "rs67376798", "rs56038477",
        "max_toxicity_grade",
    ]
    with path.open("w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        # Complete and negative, toxic.
        w.writerow(dict(zip(cols, ["P1", "0", "0", "0", "0", "3"])))
        # One locus never read, toxic — must be excluded, not counted negative.
        w.writerow(dict(zip(cols, ["P2", "0", "", "0", "0", "4"])))

    import argparse
    import io
    from contextlib import redirect_stdout

    buf = io.StringIO()
    with redirect_stdout(buf):
        cmd_quadrant(argparse.Namespace(csv=str(path), panel="cpic4"))
    out = buf.getvalue()

    assert "n = 1 patients" in out
    assert "1 excluded" in out
    assert "1/1 (100.0%) of panel-NEGATIVE" in out


# ---------------------------------------------------------------------------
# Refusal 11: a comparison group's rate is never attributed to a patient who
# is not in that group.
#
# The five tests below each pin a defect found by walking a *variant-positive*
# patient through the report, rather than the panel-negative patient every
# earlier test used. The panel-negative path was the one the package was
# designed around, so it was the one that was already correct.
# ---------------------------------------------------------------------------


def _residual_section(text: str) -> str:
    """The residual-risk block only, up to the panel-sensitivity line."""
    after = text.split("RESIDUAL RISK", 1)[1]
    return after.split("Panel sensitivity", 1)[0]


@pytest.mark.parametrize(
    "rsid",
    [
        "rs3918290",   # *2A, no function   -> Poor/Intermediate Metabolizer
        "rs55886062",  # *13, no function
        "rs67376798",  # c.2846A>T, decreased function
        "rs56038477",  # HapB3, decreased function
    ],
)
@pytest.mark.parametrize("zyg", [Zygosity.HET, Zygosity.HOM_VAR])
def test_wildtype_rate_is_withheld_from_a_carrier(rsid: str, zyg: Zygosity) -> None:
    """A carrier is not in the group Henricks and Lunenburg measured.

    Both published grade >=3 rates (22.7% systemic, 13.6% including
    chemoradiation) are rates among patients *without* an actionable DPYD
    variant. Printing one to a carrier attributes a comparison group's risk to
    someone outside it, and for a Poor Metabolizer it does so in the reassuring
    direction — the package's own central failure mode, aimed at the patient
    with the most to lose.
    """
    interp = build(_assay({rsid: zyg}))
    assert interp.residual.quoted_rate is None
    assert interp.residual.applies_to_this_patient is False
    section = _residual_section(render(interp))
    assert "22.7" not in section
    assert "13.6" not in section
    assert "DOES NOT APPLY TO THIS PATIENT" in section


def test_wildtype_rate_is_kept_for_a_normal_function_carrier() -> None:
    """`*9A` is CPIC *normal function*, so its carrier stays in the group.

    Henricks and Lunenburg defined their comparison group by the absence of an
    *actionable* (decreased/no function) variant, not by the absence of any
    variant. Withholding the rate from a `*9A` carrier would over-correct and
    strip the one number a panel-negative report exists to deliver.
    """
    calls = [
        _tested(r, Zygosity.HOM_REF) for r in CPIC_ACTIONABLE
    ] + [_tested("rs1801265", Zygosity.HET)] + [
        _untested(r)
        for r in ALL_RSIDS
        if r not in CPIC_ACTIONABLE and r != "rs1801265"
    ]
    assay = AssayResult(
        gene="DPYD",
        panel_name="CPIC 4 + *9A",
        calls=tuple(calls),
        regimen=Regimen.SYSTEMIC,
        population="Central/South Asian",
    )
    interp = build(assay)
    assert interp.residual.quoted_rate == "22.7% (231/1018)"
    assert interp.residual.applies_to_this_patient is True


def test_panel_negative_line_states_how_many_loci_returned_a_call() -> None:
    """"No variant was found" must be scoped to what was actually read.

    With three of the four CPIC loci returning no call, the earlier wording —
    "No variant was found among the loci tested" — read as a clean four-variant
    panel. One locus of reassurance was being presented as four.
    """
    calls = [_tested("rs3918290", Zygosity.HOM_REF)] + [
        _untested(r) for r in ALL_RSIDS if r != "rs3918290"
    ]
    assay = AssayResult(
        gene="DPYD",
        panel_name="CPIC 4-variant DPYD panel (1 of 4 loci returned a call)",
        calls=tuple(calls),
        regimen=Regimen.SYSTEMIC,
        population="Central/South Asian",
    )
    text = render(build(assay))
    assert "No variant was found at the 1 locus that returned a call." in text
    # Derived rather than hardcoded: the panel grows as pgx-core gains callable
    # alleles (13 -> 17 loci at engine 0.9.0), and a literal count here would
    # have to be edited every time — which is how a test starts asserting the
    # wrong thing confidently.
    assert f"{len(ALL_RSIDS) - 1} further loci were not assayed" in text
    assert "among the loci tested" not in text


@pytest.mark.parametrize("drug", list(Drug))
def test_cpic_action_always_names_the_drug_it_applies_to(drug: Drug) -> None:
    """A bare ``AVOID`` does not say what to avoid.

    CPIC publishes DPYD recommendations per drug and the engine is queried per
    drug, so the report must say which one it is displaying. The drug was
    hardcoded to capecitabine inside ``interpret`` and never surfaced, which
    made a 5-FU patient's report silently mislabelled.
    """
    assay = _assay({}, drug=drug)
    text = render(build(assay))
    assert f"(for {drug.value})" in text


def test_blank_clinical_action_is_refused_rather_than_printed() -> None:
    """The engine returns '' for an unknown pair by design; a report may not.

    ``action_for`` never raises — it is on the hot phenotype path — so an
    unrecognised phenotype/drug pair yields an empty string. Rendered, that is a
    "CPIC action :" line with nothing after it, which reads as "nothing to do"
    rather than "we could not determine what CPIC advises".
    """
    import asl.interpret as interpret_mod

    real = interpret_mod.action_for
    interpret_mod.action_for = lambda *_a, **_k: ""
    try:
        with pytest.raises(InterpretationRefused, match="no CPIC clinical action"):
            build(_assay({}))
    finally:
        interpret_mod.action_for = real


# ---------------------------------------------------------------------------
# Pending CPIC guideline changes. CPIC announced on 2026-07-09 that HapB3's
# allele value moves from 0.5 to 0.75, while explicitly leaving its database at
# 0.5 until the guideline publishes in fall 2026.
#
# This is the one class of change the --check drift guard cannot detect: the API
# and the pinned table agree, so nothing automated fires. It matters here more
# than anywhere else because HapB3 is 77.5% of the actionable allele frequency
# CPIC publishes for Central/South Asian — most Indian panel-positive patients
# are HapB3 carriers.
# ---------------------------------------------------------------------------


def test_hapb3_carrier_is_told_cpic_has_announced_a_change() -> None:
    interp = build(_assay({"rs56038477": Zygosity.HET}))
    ids = {f.rule_id for f in interp.flags}
    assert "P4_DPYD_HAPB3_CPIC_UPDATE_PENDING" in ids
    text = render(interp)
    assert "P4_DPYD_HAPB3_CPIC_UPDATE_PENDING" in text
    # Both values must appear, each labelled as what it is.
    assert "CPIC's value TODAY" in text
    assert "ANNOUNCED value (not applied)" in text


def test_the_announced_value_is_not_applied_to_the_phenotype() -> None:
    """The flag is notice, never a recalculation.

    Pre-applying a value from a pre-publication blog comment would be an
    undeclared deviation from the pinned source of truth — the same failure class
    as the 0.7.2 ``*6`` incident, pointed the other way. The engine's answer must
    be untouched.
    """
    interp = build(_assay({"rs56038477": Zygosity.HET}))
    assert interp.phenotype == "Intermediate Metabolizer"
    assert interp.cpic_recommendation_action == "REDUCE_50PCT"
    # And the evidence record still carries CPIC's published activity value.
    hapb3 = [e for e in interp.evidence if "HapB3" in e.allele_label]
    assert hapb3 and hapb3[0].activity_value == "0.5"


def test_no_pending_change_flag_for_an_allele_without_one() -> None:
    """The flag is per-allele, not a blanket disclaimer."""
    interp = build(_assay({"rs3918290": Zygosity.HET}))
    ids = {f.rule_id for f in interp.flags}
    assert "P4_DPYD_HAPB3_CPIC_UPDATE_PENDING" not in ids


def test_a_panel_negative_patient_gets_no_pending_change_flag() -> None:
    """No *pending-change* flag, which is what this test is about.

    This deliberately no longer asserts ``flags == ()``. It used to, and that was
    the wrong invariant: it encoded "a negative panel produces a clean report",
    which is the assumption P6_DPYD_CNV_NOT_INTERROGATED exists to break. A
    panel-negative SNV assay *must* carry the variant-class disclosure, because
    that is the case where the method's limit is least visible and most
    consequential. Asserting an empty flag tuple here would have made the fix
    look like a regression.
    """
    interp = build(_assay({}))
    ids = {f.rule_id for f in interp.flags}
    assert "P4_DPYD_HAPB3_CPIC_UPDATE_PENDING" not in ids
    assert not any(f.rule_id.startswith("P2_") for f in interp.flags)


# ---------------------------------------------------------------------------
# P6_DPYD_CNV_NOT_INTERROGATED — the variant class no locus-level refusal covers.
#
# Found by running the Delcuratolo 2026 patient (Front Pharmacol 17:1866621)
# through the pipeline: wild type on all four CPIC loci plus c.2194G>A, who then
# suffered grade 4 toxicity and septic shock from a heterozygous ~500 kb deletion
# of exons 5-20 that contained four of the five loci screened. asl reported
# "Tested & negative" for all of them, which is true of the bases and false of
# the chromosome — an SNV method reads only the surviving copy.
#
# This is not a NOT_TESTED case, which is why no existing refusal caught it: the
# loci WERE tested and DID return calls.
# ---------------------------------------------------------------------------


def test_an_snv_only_negative_panel_discloses_the_cnv_limit() -> None:
    """The all-negative report is the one that most needs the disclosure."""
    interp = build(_assay({}))
    ids = {f.rule_id for f in interp.flags}
    assert "P6_DPYD_CNV_NOT_INTERROGATED" in ids


def test_the_disclosure_survives_into_the_rendered_report() -> None:
    """A flag a clinician never sees is not a disclosure."""
    from asl.report import render

    text = render(build(_assay({})))
    assert "P6_DPYD_CNV_NOT_INTERROGATED" in text
    # The qualifier must sit on the line making the claim, not only in the flag.
    assert "does not mean two intact copies were observed" in text


def test_the_disclosure_is_emitted_for_a_carrier_too() -> None:
    """Finding one variant does not exclude a deletion on the other chromosome.

    The other flags are keyed to what was found. This one is keyed to what the
    method could not have seen, so it does not depend on the result.
    """
    interp = build(_assay({"rs3918290": Zygosity.HET}))
    ids = {f.rule_id for f in interp.flags}
    assert "P6_DPYD_CNV_NOT_INTERROGATED" in ids


def test_a_cnv_aware_method_gets_no_disclosure() -> None:
    """The flag states a real limit, so it must disappear when the limit does."""
    from dataclasses import replace

    from asl.models import AssayMethod

    base = _assay({})
    interp = build(replace(base, method=AssayMethod.CNV_AWARE))
    ids = {f.rule_id for f in interp.flags}
    assert "P6_DPYD_CNV_NOT_INTERROGATED" not in ids

    from asl.report import render

    assert "two intact copies were observed" not in render(interp)


def test_the_disclosure_attaches_no_probability_to_this_patient() -> None:
    """7% is the share of DPD-deficient patients, not a carrier frequency.

    Quoting it as this patient's risk would be the same denominator substitution
    ``residual`` refuses for the wild-type toxicity rate. It may appear only with
    its denominator stated.
    """
    from asl.report import render

    text = render(build(_assay({})))
    assert "7% (5/72) of DPD-deficient patients" in text
    assert "Denominator is patients already known to be DPD deficient" in text


def test_the_disclosure_recommends_no_test_and_changes_no_phenotype() -> None:
    """Recommending a test would move asl from informing to driving care."""
    interp = build(_assay({}))
    assert interp.phenotype == "Normal Metabolizer"
    assert interp.cpic_recommendation_action == "STANDARD"
    flag = next(
        f for f in interp.flags if f.rule_id == "P6_DPYD_CNV_NOT_INTERROGATED"
    )
    assert "no further testing is recommended" in flag.summary


def test_the_indian_deletion_rate_is_stated_as_unmeasured() -> None:
    """Unknown is not low, and the report must not let it read as low."""
    from asl.report import render

    text = render(build(_assay({})))
    assert "unmeasured, which is not the same as low" in text


# ---------------------------------------------------------------------------
# The shipped CSV template must actually work. Two defects found by running
# `asl interpret` against fixtures/cohort_template.csv — the file this package
# asks partner sites to fill in — rather than against a constructed dict.
# ---------------------------------------------------------------------------


def _template_path():
    from pathlib import Path

    return Path(__file__).resolve().parents[1] / "fixtures" / "cohort_template.csv"


def test_the_shipped_template_interprets_on_the_default_panel() -> None:
    """Every row of the template must produce a report, not a refusal.

    It did not. The template carries the four South Asian columns for sites that
    can fill them and leaves them blank for sites that cannot, and the off-panel
    check refused any *key* outside the panel — including a blank one. So the
    default invocation of the documented entry point refused all eight rows.

    A blank off-panel column observes nothing: it discards no genotype and widens
    no panel, so it is not the contradiction the refusal exists to catch.
    """
    import argparse
    import io
    from contextlib import redirect_stdout

    from asl.cli import cmd_interpret

    buf = io.StringIO()
    with redirect_stdout(buf):
        cmd_interpret(
            argparse.Namespace(
                csv=str(_template_path()),
                panel="cpic4",
                population="Central/South Asian",
                drug="capecitabine",
            )
        )
    out = buf.getvalue()
    assert "REFUSED" not in out
    assert out.count("DPYD PHARMACOGENOMIC INTERPRETATION") == 8


def test_a_real_off_panel_genotype_is_still_refused(tmp_path) -> None:
    """The relaxation above must not weaken the refusal it narrowed.

    A *called* genotype outside the panel is still contradictory and still a
    refusal — otherwise the panel silently widens.
    """
    import argparse
    import io
    from contextlib import redirect_stdout

    from asl.cli import cmd_interpret

    path = tmp_path / "off_panel.csv"
    path.write_text(
        "anonymised_id,rs3918290,rs55886062,rs67376798,rs56038477,rs1801265\n"
        "X1,0,0,0,0,1\n"
    )
    buf = io.StringIO()
    with redirect_stdout(buf):
        cmd_interpret(
            argparse.Namespace(
                csv=str(path),
                panel="cpic4",
                population="Central/South Asian",
                drug="capecitabine",
            )
        )
    out = buf.getvalue()
    assert "REFUSED" in out
    assert "does not cover" in out


def test_the_template_drug_column_is_read_not_defaulted() -> None:
    """The template calls it ``fluoropyrimidine_drug``; the API calls it ``drug``.

    Both are accepted. Falling back to the default instead would label the two
    5-FU patients in the template as capecitabine and say nothing about it.
    """
    import argparse
    import io
    from contextlib import redirect_stdout

    from asl.cli import cmd_interpret

    buf = io.StringIO()
    with redirect_stdout(buf):
        cmd_interpret(
            argparse.Namespace(
                csv=str(_template_path()),
                panel="cpic4",
                population="Central/South Asian",
                drug="capecitabine",
            )
        )
    out = buf.getvalue()
    assert out.count("(for fluorouracil)") == 2
    assert out.count("(for capecitabine)") == 6
