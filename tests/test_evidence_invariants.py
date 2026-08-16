"""Evidence-layer invariants.

Two rules from the architecture document, both of which exist because of the
``*6`` incident in ``anukriti-pgx-core`` 0.7.2 — a defensible clinical judgement
written straight into a pinned table with nothing recording why, which then
survived four releases as an undeclared deviation from CPIC:

1. The generated columns are copied from CPIC and never edited.
2. The hand-maintained columns never feed a phenotype.
"""

from __future__ import annotations

import json

import pytest

from asl import AssayResult, Regimen, TestStatus, VariantCall, Zygosity, build
from asl.evidence import (
    _DATA_DIR,
    EvidenceError,
    classify_basis,
    record_for,
    uninformative_in,
)
from asl.interpret import CPIC_ACTIONABLE, CPIC_ALLELE

GENES = ("DPYD", "NUDT15", "TPMT")


def _table(gene: str) -> dict:
    return json.loads((_DATA_DIR / f"{gene}_evidence.json").read_text())


@pytest.mark.parametrize("gene", GENES)
def test_evidence_table_exists_and_is_populated(gene: str) -> None:
    table = _table(gene)
    assert table["_gene"] == gene
    assert table["_allele_count"] == len(table["alleles"]) > 0
    assert "api.cpicpgx.org" in table["_source"]


@pytest.mark.parametrize("gene", GENES)
def test_every_allele_carries_the_generated_cpic_fields(gene: str) -> None:
    for label, row in _table(gene)["alleles"].items():
        for field in (
            "clinical_function",
            "evidence_strength",
            "findings_raw",
            "basis",
            "citations",
            "frequency",
        ):
            assert field in row, f"{gene} {label}: missing {field}"


@pytest.mark.parametrize("gene", GENES)
def test_basis_is_a_closed_vocabulary(gene: str) -> None:
    allowed = {"in_vitro", "ex_vivo", "in_vivo", "unclassified"}
    for label, row in _table(gene)["alleles"].items():
        assert row["basis"] in allowed, f"{gene} {label}: basis {row['basis']!r}"


def test_basis_classification_falls_back_rather_than_guessing() -> None:
    """An unrecognised findings string must not become in_vitro."""
    assert classify_basis("5-fluorouracil (in vitro)") == "in_vitro"
    assert classify_basis("Dihydrouracil (in vivo)") == "in_vivo"
    assert classify_basis("5-fluorouracil (ex vivo)") == "ex_vivo"
    assert classify_basis("something CPIC has not written before") == "unclassified"
    assert classify_basis(None) == "unclassified"
    assert classify_basis("") == "unclassified"


def test_star9a_basis_is_in_vivo_not_in_vitro() -> None:
    """The correction that mattered: *9A does not rest on Offer's in-vitro work.

    CPIC cites PMID 18452418 (He/Zhang 2008), a plasma DPD activity study in 142
    Chinese patients, at strength Strong. An earlier design document asserted all
    four South-Asian-common alleles rested on Offer 2013/2014 HEK293T assays.
    """
    rec = record_for("DPYD", "c.85T>C (*9A)", "Central/South Asian")
    assert rec.basis.value == "in_vivo"
    assert rec.evidence_strength == "Strong"
    assert "18452418" in rec.citations
    assert rec.clinical_function == "Normal function"


def test_south_asian_frequencies_match_cpic() -> None:
    """The arithmetic the problem statement rests on, read from the pinned table."""
    pop = "Central/South Asian"
    actionable = {
        "c.1905+1G>A (*2A)": 0.00508,
        "c.1679T>G (*13)": 0.0,
        "c.2846A>T": 0.00064,
        "c.1129-5923C>G, c.1236G>A (HapB3)": 0.01966,
    }
    common = {
        "c.85T>C (*9A)": 0.25526,
        "c.2194G>A (*6)": 0.09506,
        "c.1627A>G (*5)": 0.09373,
        "c.496A>G": 0.08517,
    }
    for label, expected in {**actionable, **common}.items():
        got = record_for("DPYD", label, pop).frequency_in_population
        assert got == pytest.approx(expected, abs=5e-5), f"{label}: {got}"

    a_sum = sum(
        record_for("DPYD", label, pop).frequency_in_population for label in actionable
    )
    c_sum = sum(
        record_for("DPYD", label, pop).frequency_in_population for label in common
    )
    assert a_sum == pytest.approx(0.02538, abs=1e-4)
    assert c_sum == pytest.approx(0.52921, abs=1e-4)
    # The four CPIC "Normal function" alleles are ~20x the actionable panel.
    assert c_sum / a_sum > 20


def test_star13_is_absent_in_south_asians_and_present_in_europeans() -> None:
    assert uninformative_in("DPYD", "c.1679T>G (*13)", "Central/South Asian")
    assert uninformative_in("DPYD", "c.1679T>G (*13)", "European") is None


def test_an_uncitable_allele_is_refused() -> None:
    with pytest.raises(EvidenceError, match="no CPIC evidence row"):
        record_for("DPYD", "*NOT_A_REAL_ALLELE")


# ---------------------------------------------------------------------------
# Rule 2: hand-maintained conflict metadata never changes a phenotype.
# ---------------------------------------------------------------------------


def _assay_with(rsid: str) -> AssayResult:
    calls = [
        VariantCall(
            rsid=r, status=TestStatus.TESTED, zygosity=Zygosity.HOM_REF,
            allele_label=CPIC_ALLELE[r],
        )
        for r in CPIC_ACTIONABLE
    ]
    calls.append(
        VariantCall(
            rsid=rsid, status=TestStatus.TESTED, zygosity=Zygosity.HET,
            allele_label=CPIC_ALLELE[rsid],
        )
    )
    return AssayResult(
        gene="DPYD", panel_name="conflict probe", calls=tuple(calls),
        regimen=Regimen.SYSTEMIC, population="Central/South Asian",
    )


@pytest.mark.parametrize(
    ("rsid", "label"),
    [("rs1801160", "*6"), ("rs2297595", "M166V"), ("rs1801265", "*9A")],
)
def test_alleles_with_conflict_records_still_call_normal_metabolizer(
    rsid: str, label: str
) -> None:
    """A conflict flag is display metadata. It must not move the phenotype.

    Each of these alleles has a hand-maintained conflict entry citing evidence
    against CPIC's Normal-function call. The phenotype must remain CPIC's.
    """
    interp = build(_assay_with(rsid))
    assert interp.phenotype == "Normal Metabolizer", (
        f"{label}: a conflict record changed the phenotype. Hand-maintained "
        f"metadata must never feed a phenotype call."
    )
    assert interp.cpic_recommendation_action == "STANDARD"
    assert any(f.rule_id.startswith(("P2_", "P3_")) for f in interp.flags)


def test_conflict_entries_state_both_directions() -> None:
    """An entry citing only evidence against CPIC is advocacy, not a record."""
    conflicts = json.loads(
        (_DATA_DIR / "DPYD_conflicts.json").read_text()
    )["conflicts"]
    assert conflicts, "no conflict entries found"
    for label, entry in conflicts.items():
        assert entry.get("rule_id"), f"{label}: no rule_id"
        assert entry.get("cpic_position"), f"{label}: CPIC's position is not stated"
        assert entry.get("contradicting_evidence"), f"{label}: no contrary evidence"
        assert entry.get("evidence_supporting_cpic"), (
            f"{label}: no evidence supporting CPIC. Every conflict record must "
            f"state both directions."
        )
        assert entry.get("resolution"), f"{label}: no resolution stated"


def test_haplotype_flag_carries_its_non_replication() -> None:
    """The Gentile claim must never travel without Kanai 2022's null."""
    flags = json.loads(
        (_DATA_DIR / "DPYD_conflicts.json").read_text()
    )["haplotype_flags"]
    for key, spec in flags.items():
        assert spec.get("clinical_non_replication"), f"{key}: no non-replication"
        assert spec.get("what_the_claim_is_not"), f"{key}: claim not bounded"
        assert "phase" in spec.get("phase_requirement", "").lower()
