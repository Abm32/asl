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

from .models import Interpretation, Regimen, TestStatus
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
    out.append(f"  Phenotype   : {interp.phenotype}")
    out.append(f"  CPIC action : {interp.cpic_recommendation_action}")
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
    tested = [c.allele_label or c.rsid for c in a.tested]
    neg = [c.allele_label or c.rsid for c in a.tested if c.copies == 0]
    pos = [
        f"{c.allele_label or c.rsid} ({c.zygosity.value}, "
        f"{c.copies} cop{'y' if c.copies == 1 else 'ies'})"
        for c in a.found
    ]
    out.append(f"  Tested ({len(a.tested)}): {' | '.join(tested)}")
    out.append(f"  Tested & negative : {' | '.join(neg) if neg else '(none)'}")
    out.append(f"  Tested & FOUND    : {' | '.join(pos) if pos else '(none)'}")
    stmt = untested_statement(a)
    if stmt:
        out.append("")
        for line in _wrap(stmt, "  "):
            out.append(line)

    # --- uninformative loci ----------------------------------------------
    if interp.uninformative_loci:
        out.append("")
        out.append("[ASL-POP-01]  LOCI THAT CANNOT RETURN A POSITIVE HERE")
        for label, reason in interp.uninformative_loci:
            out.append(f"  {label}")
            for line in _wrap(reason, "    "):
                out.append(line)

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
        for line in _wrap(flag.summary, "  "):
            out.append(line)
        out.append(flag.detail)
        if flag.phase_unknown:
            out.append("    phase_unknown = TRUE")

    # --- residual risk ----------------------------------------------------
    r = interp.residual
    out.append("")
    out.append("[ASL-RESID-01]  RESIDUAL RISK — THIS IS NOT AN ALL-CLEAR")
    if interp.is_panel_negative:
        out.append("  No variant was found among the loci tested.")
    out.append("")
    out.append(f'  CPIC states: "{r.cpic_caveat}."')
    out.append(f"    ({r.cpic_citation})")
    out.append("")
    if r.quoted_rate:
        for line in _wrap(
            f"Regimen-matched rate ({r.regimen.value}): {r.quoted_rate} of "
            f"patients without an actionable DPYD variant experienced grade >=3 "
            f"toxicity.",
            "  ",
        ):
            out.append(line)
        out.append(f"    Source: {r.rate_source}")
    else:
        out.append("  RATE WITHHELD:")
        for line in _wrap(r.refusal_reason or "", "    "):
            out.append(line)
    out.append("")
    for line in _wrap(f"Panel sensitivity: {r.panel_sensitivity}", "  "):
        out.append(line)

    # --- boundary ---------------------------------------------------------
    out.append("")
    out.append(_rule())
    out.append("WHAT THIS REPORT DOES NOT DO")
    for line in (
        "It does not recommend a dose. The CPIC recommendation above is CPIC's.",
        "It does not override CPIC. Where the literature disagrees, the "
        "disagreement is flagged alongside CPIC's call, never instead of it.",
        "It does not state that this patient is at elevated risk. No such "
        "finding is established.",
        "It does not treat an untested locus as a negative result.",
    ):
        for ln in _wrap(line, "  "):
            out.append(ln)
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
