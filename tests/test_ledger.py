"""Ledger invariants: no PHI, no unusable rows, append-only."""

from __future__ import annotations

import json

import pytest

from asl.ledger import (
    LedgerEntry,
    LedgerRefused,
    ToxicityEvent,
    append,
    from_mapping,
    quadrant,
    read_all,
)

PANEL = ("rs3918290", "rs55886062", "rs67376798", "rs56038477")


def _entry(**overrides) -> LedgerEntry:
    kwargs = {
        "anonymised_id": "MCC-0001",
        "cancer_type": "colorectal",
        "regimen": "CAPOX adjuvant",
        "drug": "capecitabine",
        "assay_method": "4-variant PCR panel",
        "rsids_tested": PANEL,
        "genotypes": {r: "0/0" for r in PANEL},
        "cpic_phenotype_at_dosing": "Normal Metabolizer",
        "engine_version": "0.7.3",
    }
    kwargs.update(overrides)
    return LedgerEntry(**kwargs)  # type: ignore[arg-type]


def test_a_valid_entry_records() -> None:
    e = _entry(toxicity=(ToxicityEvent("diarrhoea", 3, onset_cycle=2),))
    assert e.max_grade == 3
    assert e.as_dict()["max_toxicity_grade"] == 3


def test_rsids_tested_is_mandatory() -> None:
    """The field published Indian cohorts omit, and the reason they are unusable."""
    with pytest.raises(LedgerRefused, match="rsids_tested is empty"):
        _entry(rsids_tested=())


def test_a_genotype_cannot_come_from_an_untested_locus() -> None:
    with pytest.raises(LedgerRefused, match="not listed as tested"):
        _entry(genotypes={**{r: "0/0" for r in PANEL}, "rs1801160": "0/1"})


@pytest.mark.parametrize(
    "field_name", ["name", "mrn", "aadhaar", "phone", "dob", "address"]
)
def test_phi_field_names_are_refused_from_a_csv_row(field_name: str) -> None:
    """The realistic entry point: a site hands over a spreadsheet.

    ``LedgerEntry`` itself rejects unknown keywords with a TypeError, which is a
    stronger guarantee. But a CSV has whatever columns it has, so the denylist is
    enforced in ``from_mapping`` before anything is constructed.
    """
    row = {
        "anonymised_id": "MCC-0001",
        "cancer_type": "colorectal",
        "rsids_tested": ";".join(PANEL),
        field_name: "anything",
    }
    with pytest.raises(LedgerRefused, match="holds no PHI"):
        from_mapping(row)


def test_the_dataclass_itself_rejects_unknown_fields() -> None:
    with pytest.raises(TypeError):
        _entry(patient_name="Anita")  # type: ignore[arg-type]


def test_a_clean_csv_row_is_accepted() -> None:
    entry = from_mapping({
        "anonymised_id": "MCC-0007",
        "cancer_type": "colorectal",
        "regimen": "CAPOX adjuvant",
        "drug": "capecitabine",
        "assay_method": "4-variant PCR panel",
        "rsids_tested": ";".join(PANEL),
        "genotypes": "rs3918290=0/0;rs55886062=0/0",
        "cpic_phenotype_at_dosing": "Normal Metabolizer",
        "engine_version": "0.7.3",
        "toxicity": [{"term": "diarrhoea", "grade": 3, "onset_cycle": 2}],
    })
    assert entry.max_grade == 3
    assert entry.rsids_tested == PANEL


@pytest.mark.parametrize(
    "value",
    [
        "123456789012",          # Aadhaar-shaped
        "9876543210",            # mobile
        "patient@example.com",   # email
        "1974-03-11",            # full date
    ],
)
def test_phi_shaped_values_are_refused(value: str) -> None:
    with pytest.raises(LedgerRefused, match="holds no PHI"):
        _entry(anonymised_id=value)


def test_ctcae_grade_is_bounded() -> None:
    with pytest.raises(LedgerRefused, match="grade must be 0-5"):
        ToxicityEvent("diarrhoea", 7)


def test_append_is_append_only(tmp_path) -> None:
    path = tmp_path / "ledger.jsonl"
    append(_entry(anonymised_id="A-1"), path)
    append(_entry(anonymised_id="A-2"), path)
    rows = read_all(path)
    assert [r["anonymised_id"] for r in rows] == ["A-1", "A-2"]
    # A correction is a new row; nothing is rewritten.
    append(_entry(anonymised_id="A-1", cycles_received=8), path)
    assert len(read_all(path)) == 3


def test_quadrant_excludes_incomplete_panels(tmp_path) -> None:
    """A patient whose assay skipped a panel locus is not panel-negative."""
    entries = [
        _entry(anonymised_id="N-1").as_dict(),
        _entry(
            anonymised_id="P-1",
            genotypes={**{r: "0/0" for r in PANEL}, "rs3918290": "0/1"},
            toxicity=(ToxicityEvent("mucositis", 3),),
        ).as_dict(),
        _entry(
            anonymised_id="T-1", toxicity=(ToxicityEvent("diarrhoea", 4),)
        ).as_dict(),
        # Only two of four panel loci tested — must be excluded, not counted
        # as negative.
        _entry(
            anonymised_id="X-1",
            rsids_tested=("rs3918290", "rs55886062"),
            genotypes={"rs3918290": "0/0", "rs55886062": "0/0"},
            toxicity=(ToxicityEvent("diarrhoea", 3),),
        ).as_dict(),
    ]
    cells = quadrant(entries, PANEL)
    assert cells["negative_not_toxic"] == 1
    assert cells["positive_toxic"] == 1
    assert cells["negative_toxic"] == 1
    assert cells["excluded_incomplete_panel"] == 1


def test_entry_serialises_to_valid_json() -> None:
    payload = _entry(
        toxicity=(ToxicityEvent("HFS", 2, onset_cycle=1),),
        starting_dose_mg_m2=1000.0,
        bsa_m2=1.52,
        dose_modifications=("cycle 3: -25%",),
        cycles_received=6,
    ).as_dict()
    assert json.loads(json.dumps(payload))["max_toxicity_grade"] == 2


# ---------------------------------------------------------------------------
# The PHI scan reaches inside containers, and unreadable genotypes are a third
# state rather than a negative. Both defects were found by supplying the kind of
# row a site actually produces rather than the tidy row the tests above use.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "leak",
    [
        "mucositis, patient called on 9876543210",
        "grade 3 diarrhoea from 1978-04-12",
        "queried by onco@hospital.example",
        "referred under 123456789012",
    ],
)
def test_phi_inside_a_toxicity_term_is_refused(leak: str) -> None:
    """A free-text field nested in a tuple was never pattern-checked.

    ``as_dict`` renders ``toxicity`` as a list of dicts, and the original scan
    only pattern-checked top-level strings — so every value inside it passed
    untouched. The realistic leak is prose: a site types the phone number of the
    patient who rang about their mucositis into the event description.
    """
    with pytest.raises(LedgerRefused, match="toxicity"):
        _entry(toxicity=(ToxicityEvent(leak, 3),))


def test_phi_inside_dose_modifications_is_refused() -> None:
    """Same defect, in the other free-text container."""
    with pytest.raises(LedgerRefused, match="dose_modifications"):
        _entry(dose_modifications=("reduced 25% on 1978-04-12 per ward round",))


def test_phi_inside_genotypes_mapping_is_refused() -> None:
    """And in a nested mapping, which recursion has to reach too.

    The genotype value here is a legal token with PHI appended, because a bare
    illegal token is now refused earlier by the legibility check — and the point
    of this test is the PHI scan, not the token check.
    """
    with pytest.raises(LedgerRefused, match="genotypes"):
        _entry(genotypes={r: "0/0 read by 9876543210" for r in PANEL})


def test_phased_wildtype_is_read_as_negative_not_positive() -> None:
    """``0|0`` is legal VCF phased wild-type.

    The original check treated any token outside a four-item literal set as a
    variant call, so a phased wild-type genotype counted as panel-POSITIVE. That
    is the safe direction to be wrong for a single patient and the wrong
    direction for the cohort table: it moves patients out of the panel-negative
    cell the retrospective exists to measure.
    """
    entries = [
        _entry(
            genotypes={r: "0|0" for r in PANEL},
            toxicity=(ToxicityEvent("diarrhoea", 4),),
        ).as_dict()
    ]
    cells = quadrant(entries, PANEL)
    assert cells["negative_toxic"] == 1
    assert cells["positive_toxic"] == 0


def test_an_unreadable_genotype_is_excluded_not_counted_negative() -> None:
    """A blank or unrecognised token is neither a positive nor a negative.

    An empty string previously counted as panel-POSITIVE via the same literal
    set. Either reading is an assertion the data does not support, so the row is
    excluded and the exclusion is reported.

    The dicts are built by hand rather than through ``LedgerEntry``, because the
    ledger now refuses these tokens at entry. The ``quadrant`` guard is still
    needed and still tested: it reads JSONL that may have been written by an
    older version, or edited by hand, and defence at both layers is the point —
    the entry check stops new bad rows, this one stops old ones being miscounted.
    """
    base = _entry().as_dict()
    entries = [
        {**base, "anonymised_id": "U-1",
         "genotypes": {**{r: "0/0" for r in PANEL}, "rs3918290": ""},
         "max_toxicity_grade": 4},
        {**base, "anonymised_id": "U-2",
         "genotypes": {**{r: "0/0" for r in PANEL}, "rs3918290": "unknown"},
         "max_toxicity_grade": 0},
    ]
    cells = quadrant(entries, PANEL)
    assert cells["excluded_unreadable_genotype"] == 2
    assert cells["negative_toxic"] == 0
    assert cells["positive_toxic"] == 0


# ---------------------------------------------------------------------------
# Ledger adversarial pass. Until now the ledger had no adversarial attention at
# all — CONTRIBUTING.md names it as the open blind spot, and this is that pass.
#
# The method was the one that has worked every time: drive site-shaped data in
# and inspect what comes out, rather than re-reading the code. Eight defects.
# ---------------------------------------------------------------------------


def test_a_patient_entered_twice_is_counted_once() -> None:
    """The log is append-only, so a correction is a new row — and was double-counted.

    This is the severe one. The same ``anonymised_id`` recorded with grade 3 and
    later grade 0 produced **two** patients in **two contradictory cells** of the
    quadrant table. A site correcting an outcome is the ordinary reason for a
    duplicate, so this was not an edge case; it was the normal path.

    Later rows supersede earlier ones, and the count of superseded rows is
    reported rather than silently applied.
    """
    early = _entry(toxicity=(ToxicityEvent("diarrhoea", 3),)).as_dict()
    early["recorded_at"] = "2026-01-01T00:00:00+00:00"
    late = _entry(toxicity=(ToxicityEvent("diarrhoea", 0),)).as_dict()
    late["recorded_at"] = "2026-02-01T00:00:00+00:00"

    cells = quadrant([early, late], PANEL)
    assert cells["duplicates_superseded"] == 1
    # The later row wins: grade 0, so not toxic.
    assert cells["negative_not_toxic"] == 1
    assert cells["negative_toxic"] == 0
    assert sum(
        v for k, v in cells.items() if not k.startswith(("excluded", "duplicates"))
    ) == 1


def test_raw_log_can_still_be_audited_without_collapsing() -> None:
    """The append-only file is its own audit trail, so raw reading stays available."""
    early = _entry(toxicity=(ToxicityEvent("diarrhoea", 3),)).as_dict()
    early["recorded_at"] = "2026-01-01T00:00:00+00:00"
    late = _entry(toxicity=(ToxicityEvent("diarrhoea", 0),)).as_dict()
    late["recorded_at"] = "2026-02-01T00:00:00+00:00"

    cells = quadrant([early, late], PANEL, collapse_duplicates=False)
    assert cells["negative_toxic"] == 1
    assert cells["negative_not_toxic"] == 1
    assert cells["duplicates_superseded"] == 0


def test_an_unreadable_genotype_token_is_refused_at_entry() -> None:
    """An uninterpretable token is dropped by every later analysis silently.

    ``quadrant`` already excluded it, but the row was still written to the file,
    so a site could record a whole cohort that no analysis would ever count and
    never be told.
    """
    with pytest.raises(LedgerRefused, match="not a recognised token"):
        _entry(genotypes={**{r: "0/0" for r in PANEL}, "rs3918290": "banana"})


@pytest.mark.parametrize(
    "token", ["0/0", "0|0", "hom_ref", "wt", "0/1", "1|1", "het", "hom_var"]
)
def test_notation_variants_are_all_accepted(token: str) -> None:
    """Strict about legibility, permissive about notation.

    A site writes VCF, phased VCF or words depending on which instrument printed
    the report. Refusing a legal notation would push them to retype genotypes by
    hand, which is a worse failure than the one being prevented.
    """
    phenotype = (
        "Poor Metabolizer" if token in ("1|1", "hom_var") else "Normal Metabolizer"
    )
    entry = _entry(
        genotypes={**{r: "0/0" for r in PANEL}, "rs3918290": token},
        cpic_phenotype_at_dosing=phenotype,
    )
    assert entry.genotypes["rs3918290"] == token


def test_a_free_text_phenotype_is_refused() -> None:
    """Whether the phenotype predicted the outcome is the ledger's whole question.

    Free text in that column makes it unanalysable, which defeats the purpose of
    recording it at all.
    """
    with pytest.raises(LedgerRefused, match="not a CPIC metabolizer term"):
        _entry(cpic_phenotype_at_dosing="probably fine")


def test_a_blank_engine_version_is_refused() -> None:
    """A row that cannot be attributed to an engine version cannot be reanalysed.

    And the tables *do* change — CPIC has already announced a HapB3 revision. A
    row recorded under the old value has to be identifiable as such.
    """
    with pytest.raises(LedgerRefused, match="engine_version is empty"):
        _entry(engine_version="")


def test_a_blank_anonymised_id_is_refused() -> None:
    with pytest.raises(LedgerRefused, match="anonymised_id is empty"):
        _entry(anonymised_id="   ")


def test_a_genotype_that_contradicts_the_phenotype_is_refused() -> None:
    """A `*2A/*2A` homozygote cannot have been dosed as a Normal Metabolizer.

    CPIC calls that diplotype a Poor Metabolizer. If a row says otherwise the two
    columns came from different sources — different patients, or a hand-typed
    phenotype — and the row is evidence of nothing while polluting precisely the
    analysis the ledger exists for.
    """
    with pytest.raises(LedgerRefused, match="cannot both be true"):
        _entry(
            genotypes={**{r: "0/0" for r in PANEL}, "rs3918290": "1/1"},
            cpic_phenotype_at_dosing="Normal Metabolizer",
        )


def test_the_same_genotype_with_the_correct_phenotype_is_accepted() -> None:
    """The check catches a contradiction, not the genotype itself."""
    entry = _entry(
        genotypes={**{r: "0/0" for r in PANEL}, "rs3918290": "1/1"},
        cpic_phenotype_at_dosing="Poor Metabolizer",
    )
    assert entry.cpic_phenotype_at_dosing == "Poor Metabolizer"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"cycles_received": -5},
        {"cycles_received": 0},
        {"starting_dose_mg_m2": -1000.0},
        {"bsa_m2": 0.0},
        {"bsa_m2": 99.0},
        {"bsa_m2": 0.05},
    ],
)
def test_implausible_numerics_are_refused(kwargs: dict) -> None:
    """These are dose-intensity inputs, so a bad value corrupts the analysis quietly.

    Body surface area is a denominator: 99 m² does not fail loudly, it produces a
    dose intensity near zero and looks like extreme under-treatment.
    """
    with pytest.raises(LedgerRefused):
        _entry(**kwargs)


def test_toxicity_after_the_last_cycle_is_refused() -> None:
    """An event cannot occur in a cycle that never happened."""
    with pytest.raises(LedgerRefused, match="cycle"):
        _entry(
            toxicity=(ToxicityEvent("diarrhoea", 3, onset_cycle=99),),
            cycles_received=2,
        )


def test_onset_cycle_zero_is_refused() -> None:
    """Cycles are numbered from 1; a 0 sorts before every real event."""
    with pytest.raises(LedgerRefused, match="onset_cycle"):
        ToxicityEvent("diarrhoea", 3, onset_cycle=0)


def test_from_mapping_without_an_id_gives_a_named_refusal() -> None:
    """A bare KeyError told a site nothing about which column was missing."""
    with pytest.raises(LedgerRefused, match="no anonymised_id column"):
        from_mapping({"cancer_type": "colon", "rsids_tested": "rs3918290"})


def test_the_cli_quadrant_supersedes_a_duplicate_csv_row(tmp_path) -> None:
    """A spreadsheet gets edited by appending a corrected row below the original.

    The CLI reads a flat CSV rather than the ledger, so it has the duplicate
    problem independently and needs its own guard. Both layers report the count
    instead of applying it silently.
    """
    import argparse
    import io
    from contextlib import redirect_stdout

    from asl.cli import cmd_quadrant

    path = tmp_path / "dup.csv"
    path.write_text(
        "anonymised_id,rs3918290,rs55886062,rs67376798,rs56038477,"
        "max_toxicity_grade\n"
        "P1,0,0,0,0,4\n"      # original: severe toxicity
        "P1,0,0,0,0,0\n"      # correction appended below: none
    )
    buf = io.StringIO()
    with redirect_stdout(buf):
        cmd_quadrant(argparse.Namespace(table=str(path), panel="cpic4"))
    out = buf.getvalue()

    assert "n = 1 patients" in out
    assert "1 row(s) superseded" in out
    # The later row wins, so the single patient lands in the non-toxic cell.
    assert "0/1 (0.0%) of panel-NEGATIVE" in out
