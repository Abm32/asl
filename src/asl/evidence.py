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
from functools import cache
from pathlib import Path

from .models import Basis, EvidenceRecord

_DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "evidence"


class EvidenceError(RuntimeError):
    """Raised when the evidence layer cannot support a traceable claim."""


@cache
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


@cache
def _load_conflicts(gene: str) -> dict:
    path = _DATA_DIR / f"{gene}_conflicts.json"
    if not path.exists():
        return {"conflicts": {}, "haplotype_flags": {}}
    return json.loads(path.read_text())


def allele_names(gene: str) -> tuple[str, ...]:
    return tuple(_load_gene(gene)["alleles"].keys())


@cache
def population_keys(gene: str) -> frozenset[str]:
    """Every population key CPIC publishes a frequency under for this gene.

    Used to refuse an unrecognised population string. A typo such as
    ``"South Asian"`` for ``"Central/South Asian"`` would otherwise pass
    silently and turn every frequency lookup into ``None``, producing a report
    that looks population-aware and is not.
    """
    keys: set[str] = set()
    for row in _load_gene(gene)["alleles"].values():
        keys |= set((row.get("frequency") or {}).keys())
    if not keys:
        raise EvidenceError(
            f"{gene}: pinned evidence table publishes no population "
            f"frequencies; regenerate it before interpreting by population."
        )
    return frozenset(keys)


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


def pending_change_for(gene: str, allele_label: str) -> dict | None:
    """A change CPIC has announced but not yet published, if one exists.

    Distinct from ``conflict_for``. A conflict is the literature disagreeing with
    CPIC; this is **CPIC disagreeing with its own currently published value** and
    saying so in advance.

    It needs its own channel because the ``--check`` drift guard cannot see it.
    CPIC freezes its database until a guideline is published — for HapB3 it
    states explicitly that the allele value "will remain at AV=0.5 until then" —
    so the API and the pinned table agree, the drift check passes, and the fact
    that the recommendation is about to move is invisible to every automated
    guard we have. The announcement is the only signal, and until the guideline
    lands it has to be carried by hand.
    """
    return _load_conflicts(gene).get("pending_guideline_changes", {}).get(allele_label)


def uninformative_in(
    gene: str, allele_label: str, population: str | None
) -> str | None:
    """Reason this locus has very low prior yield in this population, or None.

    A variant at 0.00000 reference frequency in the patient's population has a
    very low prior probability of being found, so testing it adds little. Saying
    so is honest and costs nothing; not saying so lets a patient pay for a test
    that was nearly certain to be negative.

    **This is deliberately not phrased as "cannot occur."** A 0.00000 CPIC
    frequency is an *absence of observation in a reference panel*, not a proof of
    absence in the population. Chan 2024 (PMID 38886557) reported ``*13`` in a
    Tunisian patient with severe toxicity although its Middle Eastern reference
    frequency is 0%, and reported ``*2A`` in Chinese, Thai and Japanese patients
    although its East Asian reference frequency is 0%. Reference panels are small
    relative to India's population structure, so a zero there is weak evidence.
    Telling a clinician a variant "cannot" be found would license ignoring a
    positive call if the assay ever produced one — the same conversion of missing
    data into certainty that this package exists to prevent, inverted.
    """
    freq = frequency(gene, allele_label, population)
    if freq is None:
        return None
    if freq == 0.0:
        return (
            f"CPIC reports an allele frequency of 0.00000 in {population}, so a "
            f"positive result here is expected to be very rare and this locus "
            f"contributes little to this panel. NOTE: a 0.00000 reference "
            f"frequency is an absence of observation in a reference panel, not "
            f"proof the variant is absent from the population — this variant has "
            f"been reported in patients from populations where its reference "
            f"frequency is 0% (Chan 2024, PMID 38886557). A positive call at "
            f"this locus must still be acted on."
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
