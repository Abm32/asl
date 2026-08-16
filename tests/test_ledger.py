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
    kwargs = dict(
        anonymised_id="MCC-0001",
        cancer_type="colorectal",
        regimen="CAPOX adjuvant",
        drug="capecitabine",
        assay_method="4-variant PCR panel",
        rsids_tested=PANEL,
        genotypes={r: "0/0" for r in PANEL},
        cpic_phenotype_at_dosing="Normal Metabolizer",
        engine_version="0.7.3",
    )
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
