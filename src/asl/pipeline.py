"""Assemble the full deterministic interpretation. Nothing here decides anything.

The order matters and mirrors the architecture document's data flow: the engine
is called first and its answer is fixed, then provenance, population context,
flags and residual risk are attached around it.
"""

from __future__ import annotations

from . import evidence, flags, interpret, residual
from .models import AssayResult, EvidenceRecord, Interpretation

GENE = "DPYD"


def build(assay: AssayResult) -> Interpretation:
    """Produce the record a report may render. Raises rather than guessing."""
    if not assay.panel_name.strip():
        raise interpret.InterpretationRefused(
            "cannot interpret a panel whose content is unstated: supply "
            "panel_name and the full list of loci the assay covers, including "
            "those it did not find. Which variants were tested is the "
            "distinction this package exists to preserve."
        )
    if not assay.calls:
        raise interpret.InterpretationRefused(
            "cannot interpret an assay with no recorded loci"
        )

    diplotype, phenotype, action, overspecified, copies = interpret.call_engine(assay)

    # Provenance for every allele actually found.
    records: list[EvidenceRecord] = []
    for call in assay.found:
        if call.allele_label:
            records.append(
                evidence.record_for(GENE, call.allele_label, assay.population)
            )

    # Loci that cannot return a positive in this population.
    uninformative: list[tuple[str, str]] = []
    for call in assay.tested:
        if not call.allele_label:
            continue
        reason = evidence.uninformative_in(GENE, call.allele_label, assay.population)
        if reason:
            uninformative.append((call.allele_label, reason))

    collected = list(flags.conflict_flags(GENE, assay))
    hap = flags.haplotype_flag(GENE, assay)
    if hap:
        collected.append(hap)

    return Interpretation(
        assay=assay,
        diplotype=diplotype,
        phenotype=phenotype,
        cpic_recommendation_action=action,
        evidence=tuple(records),
        flags=tuple(collected),
        residual=residual.compute(assay),
        uninformative_loci=tuple(uninformative),
        overspecified=overspecified,
        copies_detected=copies,
        engine_version=interpret.engine_version(),
    )
