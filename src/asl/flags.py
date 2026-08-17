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


def variant_class_flag(gene: str, assay: AssayResult) -> Flag | None:
    """Disclose a variant class the assay method could not interrogate at all.

    Every other refusal in this package is keyed to a *locus*: it says the assay
    did not look here. This one is keyed to the *method*, and says the assay could
    not have seen this kind of variant anywhere in the gene, including at the loci
    it did read.

    That gap is why the flag is needed rather than covered by ``NOT_TESTED``. A
    locus inside a heterozygous deletion is not untested — it was assayed, it
    returned a call, and the call was ``hom_ref``, because an SNV method reads the
    surviving chromosome and cannot distinguish one reference copy from two. So
    the locus appears under "Tested & negative", which is true of the base and
    false of the chromosome, and no ``NOT_TESTED`` entry is generated to correct
    it. Without this flag the single most complete-looking report the package can
    produce — every panel locus tested, every one negative — is the one whose
    limitation is least visible.

    Deliberately does **not** attach a residual probability. The 7% figure
    (van Kuilenburg 2009, PMID 19296131) is the share of *already DPD-deficient*
    patients whose deficiency is due to a deletion; it is not a population carrier
    frequency, and printing it as a risk for this patient would be exactly the
    denominator substitution the residual module refuses elsewhere. It is quoted
    with its denominator stated, or not at all.

    Also does not recommend further testing. What to do about an unexcluded
    variant class is a clinical decision, and recommending a test would move this
    from informing clinical management to driving it.
    """
    spec = evidence.variant_class_limit(gene, assay.method.value)
    if spec is None:
        return None

    detail = (
        f"    Method: {assay.method.value} — detects "
        f"{spec.get('detects','')}.\n"
        f"    Cannot detect: {spec.get('cannot_detect','')}.\n"
        f"    Why: {spec.get('why','')}\n"
        f"    How often this matters:\n"
        f"{_fmt_studies(spec.get('frequency_evidence') or [])}\n"
        f"    Population: {spec.get('population_note','')}\n"
        f"    Detection: {spec.get('detection_methods','')}\n"
        f"    What asl does: {spec.get('what_asl_does','')}"
    )
    return Flag(
        rule_id=spec["rule_id"],
        summary=(
            "This method reads single bases, so it could not interrogate "
            "copy-number or structural variants anywhere in the gene — including "
            "at the loci it reported as negative. A locus inside a heterozygous "
            "deletion reads as homozygous reference, because only the surviving "
            "chromosome is amplified. A negative result here does not exclude "
            "this class of variant. No probability is attached, and no further "
            "testing is recommended: this states a limit of the method, not a "
            "finding about this patient."
        ),
        detail=detail,
        citations=tuple(
            e["pmid"]
            for e in (spec.get("frequency_evidence") or [])
            if e.get("pmid")
        ),
    )


def pending_change_flags(gene: str, assay: AssayResult) -> tuple[Flag, ...]:
    """Flag found alleles whose CPIC value is announced to change but unpublished.

    This is not a literature conflict — it is CPIC saying in advance that its own
    published value is about to move, while deliberately leaving the database at
    the old value until the guideline lands. The drift check therefore cannot see
    it: the API and the pinned table agree.

    The flag exists so the report says so, and it explicitly does **not**
    pre-apply the announced value. Dosing on a pre-publication blog comment would
    be an undeclared deviation from the pinned source of truth — the same failure
    class as the 0.7.2 ``*6`` incident, pointed the other way. What a clinician
    needs is CPIC's current recommendation plus notice that it is changing, and
    that is exactly what this emits.
    """
    out: list[Flag] = []
    for call in assay.found:
        label = call.allele_label
        if not label:
            continue
        spec = evidence.pending_change_for(gene, label)
        if spec is None:
            continue
        detail = (
            f"    Announced {spec.get('announced','')}, expected "
            f"{spec.get('expected_publication','')}.\n"
            f"    Source: {spec.get('source','')}\n"
            f"    CPIC's value TODAY (used above): "
            f"{spec.get('current_pinned_value','')}\n"
            f"    CPIC's ANNOUNCED value (not applied): "
            f"{spec.get('announced_value','')}\n"
            f"    Scope: {spec.get('scope_limit','')}\n"
            f"    Why this matters here: {spec.get('why_this_matters_here','')}\n"
            f"    What asl does: {spec.get('what_asl_does','')}\n"
            f"    CAVEAT: {spec.get('caveat','')}"
        )
        out.append(
            Flag(
                rule_id=spec["rule_id"],
                summary=(
                    f"{label}: CPIC has announced a revision to this allele's own "
                    f"value that is not yet published. The recommendation above is "
                    f"CPIC's current published one; asl does not pre-apply the "
                    f"announced change."
                ),
                detail=detail,
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
