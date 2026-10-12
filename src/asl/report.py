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

from .interpret import CPIC_ALLELE, assayed_variant_name
from .models import AssayMethod, Interpretation, Regimen, TestStatus
from .residual import untested_statement

_W = 76

# Plain-language gloss for each CPIC DPYD phenotype, for the answer band.
#
# The gloss restates the phenotype; it never adds a claim. "Normal Metabolizer"
# carries no statement about the loci that were not read, which is why the band
# puts the coverage strip and the residual-risk pointer immediately beneath it
# rather than leaving the phenotype alone at the top of the page.
_PHENOTYPE_GLOSS: dict[str, str] = {
    "Normal Metabolizer": "normal DPD activity predicted from the loci that were read",
    "Intermediate Metabolizer": "REDUCED DPD activity",
    "Poor Metabolizer": "little or no DPD activity",
    "Indeterminate": "DPD activity could not be determined",
}

# Coverage-strip glyphs. ASCII on purpose: this report is read in a terminal, a
# printed PDF and an email body, and a box-drawing character that renders in one
# of those is a question mark in another.
_FOUND, _NEGATIVE, _UNASSAYED = "X", "#", "."


def _coverage_strip(interp: Interpretation) -> list[str]:
    """One glyph per locus: what the assay found, read, and never looked at.

    The strip is in a fixed locus order (``CPIC_ALLELE``), so the same position
    means the same variant on every patient's report and two reports can be laid
    side by side. The point it makes in one line is the one the whole package
    exists to make: how much of the gene a guideline-conformant panel leaves
    unread.
    """
    a = interp.assay
    by_rsid = {c.rsid: c for c in a.calls}
    ordered = [r for r in CPIC_ALLELE if r in by_rsid]
    ordered += [c.rsid for c in a.calls if c.rsid not in CPIC_ALLELE]

    glyphs = []
    for rsid in ordered:
        call = by_rsid[rsid]
        if call.status is TestStatus.NOT_TESTED:
            glyphs.append(_UNASSAYED)
        elif call.copies > 0:
            glyphs.append(_FOUND)
        else:
            glyphs.append(_NEGATIVE)
    found = glyphs.count(_FOUND)
    neg = glyphs.count(_NEGATIVE)
    un = glyphs.count(_UNASSAYED)
    return [
        f"  COVERAGE  [{''.join(glyphs)}]  {len(glyphs)} DPYD loci asl can call",
        (
            f"            {_FOUND} found ({found})   {_NEGATIVE} read, no "
            f"variant ({neg})   {_UNASSAYED} NOT assayed ({un})"
        ),
    ]


def _answer_band(interp: Interpretation) -> list[str]:
    """The top of the report: the answer, the dose, and what qualifies both.

    A report whose answer has to be assembled by a reader from a diplotype, an
    enum and a residual-risk paragraph three screens down is a report that gets
    skimmed to its first reassuring word. This band is the deliverable: it says
    what the result is, what CPIC advises in CPIC's own sentence, how much of
    the gene was actually read, and how many flags qualify the answer — before
    anything else on the page.
    """
    a = interp.assay
    out = [_rule("=")]
    gloss = _PHENOTYPE_GLOSS.get(interp.phenotype)
    out.append(
        f"  ANSWER    DPYD {interp.phenotype}"
        + (f" — {gloss}" if gloss else "")
    )
    out.append(
        f"            Diplotype {interp.diplotype}"
        + (
            f"   CPIC activity score {interp.cpic_activity_score_range}"
            if interp.cpic_activity_score_range
            else ""
        )
    )
    out.append("")
    # CPIC's sentence, quoted. See interpret.cpic_recommendation for why this is
    # a quotation and not a paraphrase.
    out.extend(_hang(
        f'CPIC, for {a.drug.value}: "{interp.cpic_recommendation_text}"',
        "  DOSE      ",
    ))
    cites = ", ".join("PMID " + c for c in interp.cpic_citations) or "PMID 29152729"
    out.extend(_wrap(
        f"CPIC DPYD guideline ({cites}), quoted. This is CPIC's "
        f"recommendation. asl does not recommend doses and has not modified "
        f"this one.",
        "            ",
    ))
    out.append("")
    out.extend(_coverage_strip(interp))
    out.append("")
    n = len(interp.flags)
    out.extend(_hang(
        (f"{n} flag{'' if n == 1 else 's'} below "
         f"qualif{'ies' if n == 1 else 'y'} this answer. " if n else "")
        + "ASL-RESID-01 states what this result does not exclude; it is part "
        "of the result, not an appendix to it.",
        "  READ ON   ",
    ))
    out.append(_rule("="))
    return out


def _signatory_block(signatory: str | None) -> list[str]:
    """Who authorised this report — or an explicit statement that nobody did.

    An unsigned report prints the refusal, never a blank line under the word
    "Signed", because a blank signature line on a printed page is read as an
    oversight in the printing rather than as an absent authorisation.

    The qualification requirement is India-specific and is stated rather than
    assumed: the Supreme Court's December 2017 ruling requires a laboratory
    report to be signed by a person holding a post-graduate qualification in
    pathology. asl cannot verify a qualification, so it records what it was
    given verbatim and says that it did not verify it.
    """
    out = ["", _rule(), "AUTHORISATION"]
    if signatory:
        out.extend(_wrap(f"Authorised by: {signatory}", "  "))
        out.extend(_wrap(
            "Recorded as supplied and not verified by asl. An Indian "
            "laboratory report must be signed by a person holding a "
            "post-graduate qualification in pathology (Supreme Court of India, "
            "December 2017).",
            "  ",
        ))
    else:
        out.extend(_wrap(
            "UNSIGNED. No signatory was supplied, so this output has not been "
            "authorised by anyone and is not a laboratory report. It is the "
            "deterministic interpretation of a genotype, for review by a "
            "qualified person who then takes responsibility for it.",
            "  ",
        ))
    return out


def _rule(char: str = "-") -> str:
    return char * _W


def render(interp: Interpretation, signatory: str | None = None) -> str:
    a = interp.assay
    out: list[str] = []

    out.extend(_answer_band(interp))
    out.append("")
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
    out.extend(_signatory_block(signatory))
    out.append(_rule())
    out.append(f"Generated {interp.generated_at.isoformat(timespec='seconds')}")
    out.append(
        "Research tool. Not a diagnostic device. Informs clinical management; "
        "does not drive it."
    )
    return "\n".join(out)


def _hang(text: str, label: str) -> list[str]:
    """Wrap ``text`` under a left-hand label, continuation lines aligned under it.

    ``_wrap`` indents every line identically, which repeats the label down the
    margin once a line wraps. The answer band is built of labelled rows, so it
    needs the hanging form.
    """
    import textwrap

    return textwrap.wrap(
        text, width=_W, initial_indent=label, subsequent_indent=" " * len(label)
    ) or [label.rstrip()]


def _wrap(text: str, indent: str) -> list[str]:
    import textwrap

    return textwrap.wrap(
        text, width=_W, initial_indent=indent, subsequent_indent=indent
    ) or [indent.rstrip()]
