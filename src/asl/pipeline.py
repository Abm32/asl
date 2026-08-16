"""Assemble the full deterministic interpretation. Nothing here decides anything.

The order matters and mirrors the architecture document's data flow: the engine
is called first and its answer is fixed, then provenance, population context,
flags and residual risk are attached around it.

Four preconditions are checked before the engine is called, each of which would
otherwise produce a confidently wrong report:

* **Gene.** This module interprets DPYD only. An assay for another gene is
  refused rather than run through the DPYD caller.
* **Something tested.** An assay whose every locus is ``NOT_TESTED`` would call
  ``*1/*1`` "Normal Metabolizer" from no observations at all.
* **Labelled findings.** A found allele with no CPIC label cannot be given
  provenance or a conflict flag, so it would appear in the diplotype and nowhere
  else in the report.
* **Population key.** An unrecognised population string silently disables every
  frequency claim and the uninformative-locus check.
"""

from __future__ import annotations

from . import evidence, flags, interpret, residual
from .models import AssayResult, EvidenceRecord, Interpretation

GENE = "DPYD"


def build(assay: AssayResult) -> Interpretation:
    """Produce the record a report may render. Raises rather than guessing."""
    if assay.gene != GENE:
        raise interpret.InterpretationRefused(
            f"asl interprets {GENE} only, but this assay is for "
            f"{assay.gene!r}. Running it through the {GENE} caller would "
            f"produce a {GENE} diplotype from another gene's loci."
        )
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
    if not assay.tested:
        raise interpret.InterpretationRefused(
            f"every one of the {len(assay.calls)} recorded loci is NOT_TESTED, "
            f"so nothing was observed. The engine would return *1/*1 Normal "
            f"Metabolizer from zero observations, which is the exact conversion "
            f"of 'we did not look' into 'this patient is normal' that asl "
            f"exists to prevent."
        )
    unlabelled = [c.rsid for c in assay.found if not c.allele_label]
    if unlabelled:
        raise interpret.InterpretationRefused(
            f"variant allele(s) found at {', '.join(unlabelled)} with no CPIC "
            f"allele label. Without a label asl cannot attach provenance, a "
            f"frequency or a conflict flag, so the allele would reach the "
            f"diplotype while appearing nowhere in the evidence section. "
            f"Supply allele_label (see interpret.CPIC_ALLELE)."
        )
    if assay.population is not None:
        known = evidence.population_keys(GENE)
        if assay.population not in known:
            raise interpret.InterpretationRefused(
                f"population {assay.population!r} is not a CPIC population key. "
                f"An unrecognised key silently suppresses every frequency claim "
                f"and the uninformative-locus check, so the report would look "
                f"population-aware while being population-blind. Known keys: "
                f"{', '.join(sorted(known))}."
            )

    diplotype, phenotype, action, overspecified, copies, counts = interpret.call_engine(
        assay
    )

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
    collected.extend(flags.pending_change_flags(GENE, assay))
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
        allele_counts=tuple(sorted(counts.items())),
        engine_version=interpret.engine_version(),
    )
