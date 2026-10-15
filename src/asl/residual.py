"""Residual risk: what a panel-negative DPYD result does not exclude.

This module is the reason the package exists, and its central property is
negative: **there is no code path that returns "no risk" or "low risk".** That
is asserted by property test, not by convention.

Every quantitative claim here is published, and the qualitative claim is a CPIC
quotation rather than an assertion of ours — which keeps the novel-claim surface
minimal and is what supports an "informs clinical management" posture rather
than "drives".

The numbers, and why each is used where it is:

  Wild-type grade >=3 toxicity rate — the number that actually belongs in a
  panel-negative patient's report, because it answers "what is still true of me
  now that my test is negative?":
      systemic chemotherapy        231/1018 = 22.7%  Henricks 2018, PMID 30348537
      incl. chemoradiation         105/771  = 13.6%  Lunenburg 2018, PMID 30361102
  These are different populations at different fluoropyrimidine doses. Quoting
  one where the other applies misleads in both directions, so an unknown regimen
  produces a refusal rather than a default.

  Induction TPF gets NO rate, and that is the point. Both cohorts above are
  Dutch, largely colorectal, and predominantly capecitabine or lower-dose
  infusional 5-FU. Induction docetaxel/cisplatin/5-FU for head and neck cancer —
  5-FU 750 mg/m2/day continuous infusion D1-D5 plus two other cytotoxics, and
  one of the largest fluoropyrimidine exposures in India — is a different
  exposure with roughly threefold the observed toxicity (Patil 2016, J Community
  Support Oncol 14(10):412-419). No wild-type-stratified rate has been published
  for it, so this module refuses to quote one instead of reaching for the
  nearest available number. A missing comparator is stated as missing.

  Both are rates among patients **without** an actionable DPYD variant, so
  neither describes a carrier. A found decreased/no-function allele therefore
  also produces a refusal: attributing a comparison group's rate to a Poor
  Metabolizer understates their risk in the reassuring direction.

  Panel clinical sensitivity for grade >=3 toxicity: 3.5%-21.6% across 9 studies,
  Ontario HTA 2021, PMID 34484488. The HTA computed no pooled estimate, so the
  range is reported as a range.

  NOT USED HERE: the "20-30% of early-onset toxicity can be explained" figure.
  That is an attributable fraction for early-onset events, not a sensitivity,
  and conflating the two overstates the panel roughly threefold. It appears only
  inside the CPIC quotation below, in CPIC's own words and context.
"""

from __future__ import annotations

from . import evidence
from .interpret import CPIC_ALLELE, assayed_variant_name
from .models import AssayResult, Regimen, ResidualRisk

# Verbatim from the CPIC DPYD guideline (Amstutz et al., Clin Pharmacol Ther
# 2018;103(2):210-216). Quoted rather than paraphrased on purpose: the safety
# message is the guideline's own, so a clinician can check it against the source.
CPIC_CAVEAT = (
    "patients without a DPYD decreased/no function variant may still experience "
    "severe toxicity due to other genetic, environmental or other factors"
)
CPIC_CAVEAT_2 = (
    "a genetic test investigating only selected decreased/no function variants "
    "does not fully rule out DPD defects"
)
CPIC_CITATION = "CPIC DPYD guideline, Amstutz et al. 2018, PMID 29152729"

PANEL_SENSITIVITY = (
    "3.5-21.6% clinical sensitivity for grade >=3 fluoropyrimidine toxicity "
    "across 9 studies (Ontario Health Technology Assessment 2021, PMID 34484488; "
    "no pooled estimate was computed). Specificity 95-100%."
)

_WILDTYPE_RATES: dict[Regimen, tuple[str, str]] = {
    Regimen.SYSTEMIC: (
        "22.7% (231/1018)",
        (
            "Henricks et al. 2018, Lancet Oncol, PMID 30348537 — prospective, "
            "systemic fluoropyrimidine chemotherapy"
        ),
    ),
    Regimen.CHEMORADIATION: (
        "13.6% (105/771)",
        (
            "Lunenburg et al. 2018, Eur J Cancer, PMID 30361102 — cohort "
            "including chemoradiation at lower fluoropyrimidine doses"
        ),
    ),
}

_REFUSAL = (
    "No toxicity rate is quoted because the treatment regimen was not supplied. "
    "The published grade >=3 rate among patients without an actionable DPYD "
    "variant differs by setting (22.7% on systemic chemotherapy vs 13.6% in a "
    "cohort including chemoradiation), and quoting one where the other applies "
    "would mislead. Supply the regimen to obtain a rate."
)

# Induction TPF has no published wild-type-stratified grade >=3 rate, so it gets
# a refusal that states what IS observed on the regimen without dressing it up as
# a genotype-stratified comparator. The observed figures are whole-cohort rates,
# not rates among panel-negative patients, and the text says so — quoting them as
# if they were the wild-type comparator would repeat the error in the other
# direction.
_TPF_REFUSAL = (
    "No wild-type toxicity rate is quoted because none has been published for "
    "this regimen. Induction docetaxel/cisplatin/5-FU delivers 5-FU at "
    "750 mg/m2/day by continuous infusion on days 1-5 alongside two other "
    "cytotoxics, which is not the exposure in either cohort that produced the "
    "22.7% and 13.6% figures (Henricks 2018, PMID 30348537, and Lunenburg 2018, "
    "PMID 30361102 — both Dutch, largely colorectal, and predominantly "
    "capecitabine or lower-dose infusional 5-FU). Applying either to a patient "
    "on TPF would understate the risk: observed grade >=3 toxicity on TPF in "
    "Indian routine practice was neutropenia 56.9%, febrile neutropenia 20.7%, "
    "mucositis 67.2% and diarrhoea 74.1% (Patil 2016, J Community Support Oncol "
    "14(10):412-419, n=58). Those are WHOLE-COHORT rates, not rates among "
    "panel-negative patients, so they are not a substitute comparator either "
    "and no figure is quoted as one here. What applies is the CPIC caveat "
    "below, and the fact that this regimen's toxicity is high before genotype "
    "is considered at all."
)

_CARRIER_REFUSAL = (
    "No wild-type toxicity rate is quoted because a variant allele was found in "
    "this patient. The published grade >=3 rates are rates among patients "
    "WITHOUT an actionable DPYD variant (Henricks 2018, PMID 30348537; "
    "Lunenburg 2018, PMID 30361102). This patient is not in that group, so no "
    "figure from it is shown here — for a decreased- or no-function carrier it "
    "would understate the risk, and a number on the page is read even when the "
    "sentence around it says not to. What applies here is the CPIC "
    "recommendation shown above, together with the caveat below: the panel "
    "having found something does not mean it found everything."
)

# Loci whose variant alleles CPIC classifies as normal function, and which are
# therefore compatible with the wild-type comparison group. Read from the pinned
# CPIC evidence table rather than listed here, so the set cannot drift from it.
_NORMAL_FUNCTION = "normal function"


def actionable_loci(gene: str, loci: tuple[str, ...] | list[str]) -> tuple[str, ...]:
    """The subset of ``loci`` whose CPIC function is not normal.

    The quadrant tables divide patients into panel-positive and panel-negative
    and then print, underneath, two published rates whose own groups were
    defined by the absence of an *actionable* variant. A locus CPIC calls
    normal function cannot move a patient between those groups, so counting a
    call there as panel-positive contradicts the comparator printed beside it
    — and does so for the most common alleles in the gene, which is where it
    does the most damage: `*9A` alone is about a quarter of Central/South
    Asian chromosomes.

    This is the locus-level twin of ``_found_non_normal``, which already
    applies the same rule to found alleles. Both read the pinned CPIC table
    rather than naming alleles here, so neither can drift from it.

    An rsID with no citable CPIC row is kept, matching ``_found_non_normal``:
    when the function cannot be read, the safe direction is to let the locus
    count rather than silently drop it from the classification.
    """
    out: list[str] = []
    for rsid in loci:
        label = CPIC_ALLELE.get(rsid)
        if not label:
            out.append(rsid)
            continue
        try:
            fn = evidence.record_for(gene, label).clinical_function
        except evidence.EvidenceError:
            out.append(rsid)
            continue
        if (fn or "").strip().lower() != _NORMAL_FUNCTION:
            out.append(rsid)
    return tuple(out)


def _found_non_normal(assay: AssayResult) -> tuple[str, ...]:
    """Labels of found alleles CPIC does not call normal function.

    A found allele that CPIC calls *normal function* does not remove a patient
    from the wild-type comparison group — Henricks and Lunenburg both defined
    that group by the absence of an *actionable* (decreased/no function) variant,
    not by the absence of any variant at all. So a `*9A` carrier still belongs to
    the group the rate describes, and a `*2A` carrier does not.

    Falls back to treating an uncitable allele as non-normal, because the safe
    direction when the function call cannot be read is to withhold the rate.

    The *displayed* name is the variant the assay actually read, not CPIC's
    compound allele name: printing "Found: c.1129-5923C>G, c.1236G>A (HapB3)" for
    a tag-only positive asserts that the causal variant was observed, which is
    the exact inference CPIC requires to be qualified rather than stated. The
    function lookup still uses CPIC's name, because that is the row carrying it.
    """
    out: list[str] = []
    for call in assay.found:
        label = call.allele_label
        display = assayed_variant_name(call.rsid, label)
        if not label:
            out.append(call.rsid)
            continue
        try:
            fn = evidence.record_for(assay.gene, label).clinical_function
        except evidence.EvidenceError:
            out.append(display)
            continue
        if (fn or "").strip().lower() != _NORMAL_FUNCTION:
            out.append(display)
    return tuple(out)


def compute(assay: AssayResult) -> ResidualRisk:
    """Build the residual-risk record for an assay result.

    Returns a ``ResidualRisk`` in every case. There is no input for which this
    function reports the absence of risk — a negative panel yields a residual
    risk statement, and so does a positive one.

    A rate is quoted only when both conditions hold: the regimen is known, and
    the patient belongs to the group the rate was measured in.
    """
    non_normal = _found_non_normal(assay)
    if non_normal:
        return ResidualRisk(
            quoted_rate=None,
            rate_source=None,
            regimen=assay.regimen,
            panel_sensitivity=PANEL_SENSITIVITY,
            cpic_caveat=f"{CPIC_CAVEAT}; and {CPIC_CAVEAT_2}",
            cpic_citation=CPIC_CITATION,
            refusal_reason=(
                _CARRIER_REFUSAL
                + " Found: "
                + ", ".join(non_normal)
                + "."
            ),
            applies_to_this_patient=False,
        )

    rate, source = _WILDTYPE_RATES.get(assay.regimen, (None, None))

    if rate is None:
        refusal = (
            _TPF_REFUSAL if assay.regimen is Regimen.INDUCTION_TPF else _REFUSAL
        )
        return ResidualRisk(
            quoted_rate=None,
            rate_source=None,
            regimen=assay.regimen,
            panel_sensitivity=PANEL_SENSITIVITY,
            cpic_caveat=f"{CPIC_CAVEAT}; and {CPIC_CAVEAT_2}",
            cpic_citation=CPIC_CITATION,
            refusal_reason=refusal,
        )

    return ResidualRisk(
        quoted_rate=rate,
        rate_source=source,
        regimen=assay.regimen,
        panel_sensitivity=PANEL_SENSITIVITY,
        cpic_caveat=f"{CPIC_CAVEAT}; and {CPIC_CAVEAT_2}",
        cpic_citation=CPIC_CITATION,
        refusal_reason=None,
    )


def untested_statement(assay: AssayResult) -> str | None:
    """State plainly that untested loci establish nothing.

    Returns None when the assay tested everything it could have.
    """
    if not assay.not_tested:
        return None
    labels = ", ".join(
        c.allele_label or c.rsid for c in assay.not_tested
    )
    return (
        f"Absence is NOT established for: {labels}. "
        f"These variants were not assayed by '{assay.panel_name}'. "
        f"An untested locus is not a negative result."
    )
