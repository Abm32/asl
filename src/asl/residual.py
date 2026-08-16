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
  patient report, because it answers "what is still true of me now that my test
  is negative?":
      systemic chemotherapy        231/1018 = 22.7%  Henricks 2018, PMID 30348537
      incl. chemoradiation         105/771  = 13.6%  Lunenburg 2018, PMID 30361102
  These are different populations at different fluoropyrimidine doses. Quoting
  one where the other applies misleads in both directions, so an unknown regimen
  produces a refusal rather than a default.

  Panel clinical sensitivity for grade >=3 toxicity: 3.5%-21.6% across 9 studies,
  Ontario HTA 2021, PMID 34484488. The HTA computed no pooled estimate, so the
  range is reported as a range.

  NOT USED HERE: the "20-30% of early-onset toxicity can be explained" figure.
  That is an attributable fraction for early-onset events, not a sensitivity,
  and conflating the two overstates the panel roughly threefold. It appears only
  inside the CPIC quotation below, in CPIC's own words and context.
"""

from __future__ import annotations

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
        "Henricks et al. 2018, Lancet Oncol, PMID 30348537 — prospective, "
        "systemic fluoropyrimidine chemotherapy",
    ),
    Regimen.CHEMORADIATION: (
        "13.6% (105/771)",
        "Lunenburg et al. 2018, Eur J Cancer, PMID 30361102 — cohort including "
        "chemoradiation at lower fluoropyrimidine doses",
    ),
}

_REFUSAL = (
    "No toxicity rate is quoted because the treatment regimen was not supplied. "
    "The published grade >=3 rate among patients without an actionable DPYD "
    "variant differs by setting (22.7% on systemic chemotherapy vs 13.6% in a "
    "cohort including chemoradiation), and quoting one where the other applies "
    "would mislead. Supply the regimen to obtain a rate."
)


def compute(assay: AssayResult) -> ResidualRisk:
    """Build the residual-risk record for an assay result.

    Returns a ``ResidualRisk`` in every case. There is no input for which this
    function reports the absence of risk — a negative panel yields a residual
    risk statement, and so does a positive one.
    """
    rate, source = _WILDTYPE_RATES.get(assay.regimen, (None, None))

    return ResidualRisk(
        quoted_rate=rate,
        rate_source=source,
        regimen=assay.regimen,
        panel_sensitivity=PANEL_SENSITIVITY,
        cpic_caveat=f"{CPIC_CAVEAT}; and {CPIC_CAVEAT_2}",
        cpic_citation=CPIC_CITATION,
        refusal_reason=None if rate else _REFUSAL,
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
