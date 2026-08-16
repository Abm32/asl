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
    interp = build(_assay({}))
    assert interp.engine_version.startswith("0.7.")
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
