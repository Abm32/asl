"""Load the CPIC-generated evidence tables and attach provenance per allele.

Reads only the pinned artifacts in ``data/evidence/``. No network access, which
is what lets this run in a hospital room with no connectivity.

Refuses when an allele row lacks citations or an evidence strength: a report
whose claims cannot be traced to a source is the failure this package exists to
prevent, so an untraceable row is an error rather than a blank field.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

from .models import Basis, EvidenceRecord

_DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "evidence"


class EvidenceError(RuntimeError):
    """Raised when the evidence layer cannot support a traceable claim."""


@lru_cache(maxsize=None)
def _load_gene(gene: str) -> dict:
    path = _DATA_DIR / f"{gene}_evidence.json"
    if not path.exists():
        raise EvidenceError(
            f"no pinned evidence table for {gene} at {path}. "
            f"Run: python scripts/regen_evidence_tables.py --gene {gene}"
        )
    data = json.loads(path.read_text())
    if not data.get("alleles"):
        raise EvidenceError(f"{path} contains no alleles")
    return data


@lru_cache(maxsize=None)
def _load_conflicts(gene: str) -> dict:
    path = _DATA_DIR / f"{gene}_conflicts.json"
    if not path.exists():
        return {"conflicts": {}, "haplotype_flags": {}}
    return json.loads(path.read_text())


def allele_names(gene: str) -> tuple[str, ...]:
    return tuple(_load_gene(gene)["alleles"].keys())


def frequency(gene: str, allele_label: str, population: str | None) -> float | None:
    """CPIC allele frequency in one population, or None if unavailable.

    Deliberately returns None rather than falling back to a European frequency.
    Substituting one population's frequency for another is exactly the error
    this package documents.
    """
    if population is None:
        return None
    row = _load_gene(gene)["alleles"].get(allele_label)
    if row is None:
        return None
    return (row.get("frequency") or {}).get(population)


def record_for(
    gene: str, allele_label: str, population: str | None = None
) -> EvidenceRecord:
    """Build the provenance record for one allele, entirely from CPIC data."""
    row = _load_gene(gene)["alleles"].get(allele_label)
    if row is None:
        raise EvidenceError(
            f"{gene}: no CPIC evidence row for allele {allele_label!r}. "
            f"asl will not report an allele it cannot cite."
        )
    if not row.get("citations"):
        raise EvidenceError(
            f"{gene} {allele_label}: CPIC row carries no citations; "
            f"asl will not make an uncitable claim."
        )
    if not row.get("evidence_strength"):
        raise EvidenceError(
            f"{gene} {allele_label}: CPIC row carries no evidence strength; "
            f"asl will not present a call whose strength is unknown."
        )

    return EvidenceRecord(
        allele_label=allele_label,
        clinical_function=row["clinical_function"],
        activity_value=row.get("activity_value"),
        evidence_strength=row.get("evidence_strength"),
        findings_raw=row.get("findings_raw"),
        basis=Basis(row.get("basis", "unclassified")),
        citations=tuple(row.get("citations") or []),
        frequency_in_population=frequency(gene, allele_label, population),
        population=population,
        biochemical_function=row.get("biochemical_function"),
    )


def conflict_for(gene: str, allele_label: str) -> dict | None:
    """Hand-maintained conflict entry for an allele, if one exists."""
    return _load_conflicts(gene).get("conflicts", {}).get(allele_label)


def haplotype_flags(gene: str) -> dict:
    return _load_conflicts(gene).get("haplotype_flags", {})


def uninformative_in(
    gene: str, allele_label: str, population: str | None
) -> str | None:
    """Reason this locus is uninformative in this population, or None.

    A variant at frequency 0.0 in the patient's population cannot be found in
    that patient, so testing it establishes nothing. Saying so is honest and
    costs nothing; not saying so lets a patient pay for a test that could not
    have returned a positive.
    """
    freq = frequency(gene, allele_label, population)
    if freq is None:
        return None
    if freq == 0.0:
        return (
            f"CPIC allele frequency in {population} is 0.00000 — "
            f"this locus cannot return a positive result in this population"
        )
    return None

# Ordered most specific first. CPIC `findings` strings look like
# "5-fluorouracil (in vitro)", "Dihydrouracil (in vivo)", "5-fluorouracil (ex vivo)".
_BASIS_PATTERNS: tuple[tuple[str, str], ...] = (
    (r"\bex\s*[-\s]?vivo\b", "ex_vivo"),
    (r"\bin\s*[-\s]?vitro\b", "in_vitro"),
    (r"\bin\s*[-\s]?vivo\b", "in_vivo"),
)


def classify_basis(findings_raw: str | None) -> str:
    """Classify CPIC's free-text ``findings`` into an assay basis.

    Returns ``in_vitro``, ``ex_vivo``, ``in_vivo``, or ``unclassified``.

    The fallback is deliberately ``unclassified`` rather than a guess. An
    unrecognised string means CPIC wrote something this function has not seen,
    and the honest report of that is "we do not know", not "probably in vitro".

    This is the canonical implementation; ``scripts/regen_evidence_tables.py``
    imports it so the pinned table and the library can never disagree.
    """
    if not findings_raw:
        return "unclassified"
    text = findings_raw.lower()
    for pattern, basis in _BASIS_PATTERNS:
        if re.search(pattern, text):
            return basis
    return "unclassified"
