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


#: Populations for which the European-derived four-variant panel is a poorer
#: instrument. Keys are CPIC population strings; "European" is deliberately
#: absent, and so is "Latino" — the systematic review below covers South Asian,
#: East Asian, Middle Eastern and African ancestry toxicity cases, and silence
#: about Latino patients is not the same as evidence of adequacy.
_NON_EUROPEAN_POPULATIONS = frozenset(
    {
        "Central/South Asian",
        "East Asian",
        "African American/Afro-Caribbean",
        "Sub-Saharan African",
    }
)

#: Alleles whose absence from an assay matters most outside European ancestry.
#: rsID -> (name, why it matters here).
_ANCESTRY_RELEVANT_LOCI = {
    "rs112766203": (
        "c.2279C>T",
        (
            "0.63% in GenomeIndia (n=9,768) and 0.000% in Europeans; AMP/ACMG "
            "2024 Tier 1, and absent from every commercial platform that "
            "consensus surveyed"
        ),
    ),
    "rs115232898": (
        "c.557A>G",
        (
            "AMP/ACMG 2024 Tier 1; 'relatively strong evidence' in African "
            "ancestry toxicity cases (Chan 2024) and not in the UK NHS panel"
        ),
    ),
    "rs146356975": (
        "c.868A>G",
        "AMP/ACMG 2024 Tier 1; reported in African ancestry",
    ),
}


def ancestry_panel_flag(gene: str, assay: AssayResult) -> Flag | None:
    """Disclose that a European-derived panel is a weaker test in this patient.

    The four variants that every guideline-conformant DPYD panel tests were
    selected on European data, and the evidence that they transfer is poor.
    Chan, Zhang & Pirmohamed (*Br J Cancer* 2024, PMID 38886557) reviewed 32
    studies, 53 variants and five ancestry groups: **of the four, only
    c.1905+1G>A also appears in South Asian, East Asian and Middle Eastern
    toxicity cases.** c.557A>G has relatively strong evidence in African
    ancestry and is not in the UK panel at all.

    The NHS Race and Health Observatory said it plainly in July 2024: the
    tested variants are "found primarily in people with White European
    ancestry" and the programme "has the potential to worsen health and race
    inequalities." That is a national health service describing its own
    programme, which is why this flag quotes it rather than arguing the point
    in our own voice.

    Fires only when the assay actually failed to read one of the loci that
    matter here. A panel that covered them has no gap to disclose, and a flag
    that fired anyway would train a reader to skip it.

    Does **not** change the phenotype, the diplotype or the CPIC action, and
    does not recommend a different test — which assay to run is a clinical and
    procurement decision, and recommending one would move this package from
    informing clinical management to directing it.
    """
    if gene != "DPYD" or assay.population not in _NON_EUROPEAN_POPULATIONS:
        return None

    tested_rsids = {c.rsid for c in assay.tested}
    missing = [
        (rsid, name, why)
        for rsid, (name, why) in _ANCESTRY_RELEVANT_LOCI.items()
        if rsid not in tested_rsids
    ]
    if not missing:
        return None

    lines = "\n".join(
        f"      {name} ({rsid}) — {why}" for rsid, name, why in missing
    )
    return Flag(
        rule_id="P7_DPYD_ANCESTRY_PANEL_MISMATCH",
        summary=(
            f"This panel was derived from European data, and this patient's "
            f"recorded population is {assay.population}. "
            f"{len(missing)} allele(s) that matter in non-European ancestry "
            f"were not assayed."
        ),
        detail=(
            "    Not assayed here, and relevant in this population:\n"
            f"{lines}\n"
            "    Chan TH, Zhang JE, Pirmohamed M. Br J Cancer 2024 (PMID\n"
            "    38886557), 32 studies / 53 variants / 5 ancestry groups: of\n"
            "    the four variants a guideline-conformant panel tests, ONLY\n"
            "    c.1905+1G>A also appears in South Asian, East Asian and\n"
            "    Middle Eastern toxicity cases.\n"
            "    NHS Race and Health Observatory, July 2024: the tested\n"
            "    variants are \"found primarily in people with White European\n"
            "    ancestry\" and the programme \"has the potential to worsen\n"
            "    health and race inequalities.\"\n"
            "    This does not change the CPIC call above. It states that a\n"
            "    negative result from this panel is weaker evidence in this\n"
            "    patient than the same result would be in a European one."
        ),
        citations=("38886557",),
    )


def monitoring_flag(gene: str, assay: AssayResult) -> Flag | None:
    """The label's monitoring instruction, which does not depend on genotype.

    Fires on every DPYD result, including — especially — the ones that look
    clean. The FDA capecitabine label instructs withholding or discontinuing
    "based on clinical assessment of the onset, duration, and severity of
    adverse reactions", with no genotype condition attached, and of its own
    variant table states: "This is not a complete listing of all DPYD variants
    that may result in DPD deficiency."

    The EMA's Direct Healthcare Professional Communication for capecitabine is
    blunter still: "Despite negative test results for DPD deficiency, severe
    toxicity may still occur."

    This is the one flag in the package that is unconditional, and that is the
    point. Every other disclosure here answers "what did this particular assay
    fail to establish?" This one answers "what remains true no matter what the
    assay found?" — and the answer is that genotype-guided dosing does not
    replace clinical monitoring. A report whose every caveat is contingent
    invites the reading that a sufficiently good assay would have no caveats.
    """
    if gene != "DPYD":
        return None
    return Flag(
        rule_id="P8_DPYD_MONITORING_UNCONDITIONAL",
        summary=(
            "Genotype-guided dosing does not replace clinical monitoring, "
            "whatever this result says."
        ),
        detail=(
            "    FDA capecitabine label: withhold or permanently discontinue\n"
            "    \"based on clinical assessment of the onset, duration, and\n"
            "    severity of adverse reactions\" — stated without any genotype\n"
            "    condition. Of its own DPYD variant table the label states:\n"
            "    \"This is not a complete listing of all DPYD variants that may\n"
            "    result in DPD deficiency.\"\n"
            "    EMA Direct Healthcare Professional Communication,\n"
            "    capecitabine: \"Despite negative test results for DPD\n"
            "    deficiency, severe toxicity may still occur.\"\n"
            "    NCCN Colon Cancer v5.2025: no FDA-authorised DPYD test is\n"
            "    available, and laboratory-developed tests vary in accuracy\n"
            "    and in which variants they detect.\n"
            "    asl does not tell you how to monitor. It states that the\n"
            "    obligation is not discharged by this result."
        ),
        citations=(),
    )
