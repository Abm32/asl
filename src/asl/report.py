"""Render the clinician-facing report. Deterministic; renders only the record.

Every line traces to something in the ``Interpretation``. There is no template
text asserting a clinical fact, and no code path that prints an all-clear.

Two engine artifacts are deliberately suppressed:

* the ``-1.0`` activity-score sentinel, which DPYD returns because it resolves
  through a named-diplotype table rather than a score sum. CPIC's per-allele
  activity values are printed instead.
* a truncated diplotype, which ``interpret`` refuses upstream rather than render.
"""

from __future__ import annotations

from .interpret import assayed_variant_name
from .models import AssayMethod, Interpretation, Regimen
from .residual import untested_statement

_W = 76


def _rule(char: str = "-") -> str:
    return char * _W


def render(interp: Interpretation) -> str:
    a = interp.assay
    out: list[str] = []

    out.append(_rule("="))
    out.append("DPYD PHARMACOGENOMIC INTERPRETATION")
    out.append(
        f"Panel: {a.panel_name}"
        + (f"   Regimen: {a.regimen.value}" if a.regimen is not Regimen.UNKNOWN else "")
        + (f"   Population: {a.population}" if a.population else "")
    )
    out.append(_rule("="))

    # --- CPIC call, attributed -------------------------------------------
    out.append("")
    out.append("CPIC CALL (unmodified, attributed to CPIC)")
    out.append(f"  Diplotype   : {interp.diplotype}")
    if interp.allele_counts:
        # pgx-core 0.7.3 §8 item 1: show every allele the engine detected, with
        # its copy count, so a reader can confirm the two named alleles account
        # for everything found rather than taking it on trust. Before 0.7.3 the
        # engine truncated silently; asl now refuses that case outright, and this
        # line makes the non-truncated case auditable too.
        detected = ", ".join(
            f"{allele}x{n}" for allele, n in interp.allele_counts
        )
        out.append(
            f"  Detected    : {detected}  "
            f"({interp.copies_detected} variant "
            f"cop{'y' if interp.copies_detected == 1 else 'ies'} total)"
        )
    out.append(f"  Phenotype   : {interp.phenotype}")
    out.append(f"  CPIC action : {interp.cpic_recommendation_action}  (for {a.drug.value})")
    out.append(
        "  Source      : CPIC DPYD guideline (Amstutz 2018, PMID 29152729) via "
        f"anukriti-pgx-core {interp.engine_version}"
    )
    out.append(
        "  NOTE: the dose recommendation above is CPIC's. asl does not "
        "recommend doses."
    )

    # --- tested vs not tested --------------------------------------------
    out.append("")
    out.append("[ASL-ASSAY-01]  WHAT WAS TESTED, AND WHAT WAS NOT")
    # Names the variant the assay actually read, not CPIC's compound allele
    # name. For the HapB3 tag those differ, and using the compound name here
    # would print "c.1129-5923C>G, c.1236G>A (HapB3)" as tested while the
    # untested list below correctly says c.1129-5923C>G was never assayed.
    tested = [assayed_variant_name(c.rsid, c.allele_label) for c in a.tested]
    neg = [
        assayed_variant_name(c.rsid, c.allele_label)
        for c in a.tested
        if c.copies == 0
    ]
    pos = [
        f"{assayed_variant_name(c.rsid, c.allele_label)} ({c.zygosity.value}, "
        f"{c.copies} cop{'y' if c.copies == 1 else 'ies'})"
        for c in a.found
    ]
    out.append(f"  Tested ({len(a.tested)}): {' | '.join(tested)}")
    out.append(f"  Tested & negative : {' | '.join(neg) if neg else '(none)'}")
    out.append(f"  Tested & FOUND    : {' | '.join(pos) if pos else '(none)'}")
    # "Negative" above is a statement about a base, not about a chromosome. On an
    # SNV-only method a locus inside a heterozygous deletion also reads negative,
    # so the qualifier belongs on the line that makes the claim rather than only
    # in a flag further down that a reader may not reach.
    if neg and a.method is AssayMethod.SNV_ONLY:
        out.extend(_wrap(
            "Negative above means no variant base was detected at that "
            "position. It does not mean two intact copies were observed: this "
            "method cannot detect a deletion of the locus, which reads as "
            "homozygous reference. See P6_DPYD_CNV_NOT_INTERROGATED below.",
            "  ",
        ))
    stmt = untested_statement(a)
    if stmt:
        out.append("")
        out.extend(_wrap(stmt, "  "))

    # --- uninformative loci ----------------------------------------------
    if interp.uninformative_loci:
        out.append("")
        out.append("[ASL-POP-01]  LOCI WITH VERY LOW EXPECTED YIELD HERE")
        for label, reason in interp.uninformative_loci:
            out.append(f"  {label}")
            out.extend(_wrap(reason, "    "))

    # --- evidence provenance ---------------------------------------------
    if interp.evidence:
        out.append("")
        out.append("[ASL-EVID-01]  EVIDENCE BASIS OF EACH FUNCTION CALL")
        for rec in interp.evidence:
            out.append(f"  {rec.allele_label}")
            out.append(
                f"    CPIC function : {rec.clinical_function}"
                + (f" (activity {rec.activity_value})" if rec.activity_value else "")
            )
            out.append(f"    Strength      : {rec.evidence_strength}")
            out.append(
                f"    Basis         : {rec.basis.value}"
                + (f"  [CPIC findings: {rec.findings_raw}]" if rec.findings_raw else "")
            )
            out.append(f"    Citations     : {', '.join('PMID ' + c for c in rec.citations)}")
            if rec.frequency_in_population is not None:
                out.append(
                    f"    Frequency     : {rec.frequency_in_population:.5f} "
                    f"in {rec.population} (CPIC)"
                )

    # --- flags ------------------------------------------------------------
    for flag in interp.flags:
        out.append("")
        out.append(f"[{flag.rule_id}]")
        out.extend(_wrap(flag.summary, "  "))
        out.append(flag.detail)
        if flag.phase_unknown:
            out.append("    phase_unknown = TRUE")

    # --- residual risk ----------------------------------------------------
    r = interp.residual
    out.append("")
    out.append("[ASL-RESID-01]  RESIDUAL RISK — THIS IS NOT AN ALL-CLEAR")
    if interp.is_panel_negative:
        n = len(a.tested)
        out.append(
            f"  No variant was found at the {n} "
            f"{'locus' if n == 1 else 'loci'} that returned a call."
        )
        if a.not_tested:
            m = len(a.not_tested)
            out.append(
                f"  {m} further {'locus was' if m == 1 else 'loci were'} not "
                f"assayed — see ASL-ASSAY-01 above. This is not a negative "
                f"result for {'it' if m == 1 else 'them'}."
            )
    out.append("")
    out.append(f'  CPIC states: "{r.cpic_caveat}."')
    out.append(f"    ({r.cpic_citation})")
    out.append("")
    if r.quoted_rate:
        out.extend(_wrap(
            f"Regimen-matched rate ({r.regimen.value}): {r.quoted_rate} of "
            f"patients without an actionable DPYD variant experienced grade >=3 "
            f"toxicity.",
            "  ",
        ))
        out.append(f"    Source: {r.rate_source}")
    else:
        out.append(
            "  RATE WITHHELD — DOES NOT APPLY TO THIS PATIENT:"
            if not r.applies_to_this_patient
            else "  RATE WITHHELD:"
        )
        out.extend(_wrap(r.refusal_reason or "", "    "))
    out.append("")
    out.extend(_wrap(f"Panel sensitivity: {r.panel_sensitivity}", "  "))

    # --- boundary ---------------------------------------------------------
    out.append("")
    out.append(_rule())
    out.append("WHAT THIS REPORT DOES NOT DO")
    for line in (
        "It does not recommend a dose. The CPIC recommendation above is CPIC's.",
        (
            "It does not override CPIC. Where the literature disagrees, the "
            "disagreement is flagged alongside CPIC's call, never instead of it."
        ),
        (
            "It does not state that this patient is at elevated risk. No such "
            "finding is established."
        ),
        "It does not treat an untested locus as a negative result.",
        (
            "It does not treat a negative result as evidence against variant "
            "classes the method could not detect."
        ),
    ):
        out.extend(_wrap(line, "  "))
    out.append(_rule())
    out.append(f"Generated {interp.generated_at.isoformat(timespec='seconds')}")
    out.append(
        "Research tool. Not a diagnostic device. Informs clinical management; "
        "does not drive it."
    )
    return "\n".join(out)


def _wrap(text: str, indent: str) -> list[str]:
    import textwrap

    return textwrap.wrap(
        text, width=_W, initial_indent=indent, subsequent_indent=indent
    ) or [indent.rstrip()]
