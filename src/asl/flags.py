"""Named uncertainties. Each is a flag alongside the CPIC call, never instead of it.

Two hard rules, both tested:

* A flag never alters a phenotype, an activity score, or a clinical action.
* A conflict is stated in **both** directions. An entry citing only the evidence
  against CPIC is advocacy, not a conflict record — and the platform's own ``U4``
  incident is what that rule exists to prevent.
"""

from __future__ import annotations

from . import evidence
from .models import AssayResult, Flag


def _fmt_studies(entries: list[dict]) -> str:
    parts = []
    for e in entries:
        bits = [e.get("citation", "?")]
        if e.get("pmid"):
            bits.append(f"PMID {e['pmid']}")
        if e.get("cohort_n"):
            bits.append(f"n={e['cohort_n']}")
        if e.get("cohort_ancestry"):
            bits.append(str(e["cohort_ancestry"]))
        head = "; ".join(bits)
        parts.append(f"      - {head}\n        {e.get('finding','')}")
        if e.get("caveats"):
            parts.append(f"        CAVEAT: {e['caveats']}")
    return "\n".join(parts)


def conflict_flags(gene: str, assay: AssayResult) -> tuple[Flag, ...]:
    """Emit a conflict flag for each found allele that has a conflict record."""
    out: list[Flag] = []
    for call in assay.found:
        label = call.allele_label
        if not label:
            continue
        entry = evidence.conflict_for(gene, label)
        if entry is None:
            continue

        against = entry.get("contradicting_evidence") or []
        supporting = entry.get("evidence_supporting_cpic") or []

        detail = (
            f"    CPIC position: {entry.get('cpic_position','')}\n"
            f"    Evidence pointing the other way:\n{_fmt_studies(against)}\n"
            f"    Evidence consistent with CPIC:\n{_fmt_studies(supporting)}\n"
            f"    Resolution: {entry.get('resolution','')}"
        )
        cites = tuple(
            e["pmid"] for e in (against + supporting) if e.get("pmid")
        )
        out.append(
            Flag(
                rule_id=entry["rule_id"],
                summary=(
                    f"{label}: the literature disagrees about this allele. "
                    f"The CPIC call above stands; asl does not reclassify it."
                ),
                detail=detail,
                citations=cites,
            )
        )
    return tuple(out)


def haplotype_flag(gene: str, assay: AssayResult) -> Flag | None:
    """Emit the Gentile-2016 co-occurrence flag when all three alleles are present.

    Always ``phase_unknown``: phase cannot be determined from unphased
    genotypes, so this is a flag and never a call. The clinical
    non-replication travels with it by construction — the claim is an ex-vivo
    degradation-rate finding, and reporting it without Kanai 2022's null would
    overstate it.
    """
    found_labels = {c.allele_label for c in assay.found}
    for spec in evidence.haplotype_flags(gene).values():
        required = set(spec.get("alleles_required") or [])
        if not required or not required.issubset(found_labels):
            continue

        detail = (
            f"    This patient carries all of: {', '.join(sorted(required))}.\n"
            f"    Claim: {spec.get('claim','')}\n"
            f"      Source: {spec.get('citation','')} "
            f"(PMID {spec.get('pmid','?')}, n={spec.get('cohort_n','?')}, "
            f"{spec.get('cohort_ancestry','?')})\n"
            f"    What the claim is NOT: {spec.get('what_the_claim_is_not','')}\n"
            f"    Clinical non-replication:\n"
            f"{_fmt_studies(spec.get('clinical_non_replication') or [])}\n"
            f"    Mechanistic support:\n"
            f"{_fmt_studies(spec.get('mechanistic_support') or [])}\n"
            f"    Phase: {spec.get('phase_requirement','')}"
        )
        return Flag(
            rule_id=spec["rule_id"],
            summary=(
                "Three co-occurring alleles reported to reduce 5-FU degradation "
                "ex vivo. Not replicated clinically. Phase unknown. This is a "
                "flag, not a finding of elevated risk."
            ),
            detail=detail,
            citations=(str(spec.get("pmid")),) if spec.get("pmid") else (),
            phase_unknown=True,
        )
    return None
