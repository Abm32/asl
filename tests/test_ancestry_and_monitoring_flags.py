"""Two disclosures that do not depend on what the assay found.

Every other refusal in this package answers *"what did this assay fail to
establish?"* These two answer different questions:

``P7_DPYD_ANCESTRY_PANEL_MISMATCH`` — *is this panel the right instrument for
this patient?* The four variants a guideline-conformant DPYD panel tests were
selected on European data. Chan, Zhang & Pirmohamed (Br J Cancer 2024, PMID
38886557) reviewed 32 studies across five ancestry groups: **only
c.1905+1G>A also appears in South Asian, East Asian and Middle Eastern
toxicity cases.** The NHS Race and Health Observatory (July 2024) says the
tested variants are "found primarily in people with White European ancestry"
and the programme "has the potential to worsen health and race inequalities".

``P8_DPYD_MONITORING_UNCONDITIONAL`` — *what remains true whatever the assay
found?* The FDA capecitabine label instructs withholding or discontinuing on
clinical assessment with no genotype condition, and says of its own variant
table that it "is not a complete listing". The EMA DHPC: "Despite negative
test results for DPD deficiency, severe toxicity may still occur."

A report whose every caveat is contingent invites the reading that a
sufficiently good assay would have no caveats. P8 exists to deny that.
"""

from __future__ import annotations

import pytest

from asl.cli import PANELS, _all_negative, _assay
from asl.models import Drug, Regimen
from asl.pipeline import build

SAS = "Central/South Asian"


def _interp(panel="cpic4", population=SAS, genotypes=None):
    g = _all_negative(panel) | (genotypes or {})
    return build(_assay(panel, g, Regimen.SYSTEMIC, population, Drug.CAPECITABINE))


def _rule_ids(interp):
    return {f.rule_id for f in interp.flags}


def _flag(interp, rule_id):
    return next(f for f in interp.flags if f.rule_id == rule_id)


# ---------------------------------------------------------------------------
# P7 — the panel is a weaker instrument here
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "population",
    [
        "Central/South Asian",
        "East Asian",
        "African American/Afro-Caribbean",
        "Sub-Saharan African",
    ],
)
def test_a_european_panel_is_disclosed_as_weaker_outside_europe(population) -> None:
    interp = _interp(population=population)
    assert "P7_DPYD_ANCESTRY_PANEL_MISMATCH" in _rule_ids(interp)


def test_it_does_not_fire_for_a_european_patient() -> None:
    """The panel was built for this patient. There is no mismatch to disclose."""
    assert "P7_DPYD_ANCESTRY_PANEL_MISMATCH" not in _rule_ids(
        _interp(population="European")
    )


def test_a_panel_that_covers_the_loci_has_no_gap_to_disclose() -> None:
    """The point of the flag, demonstrated by its own silence.

    `amp_tier1` assays c.2279C>T, c.557A>G and c.868A>G, so for a South Asian
    patient there is nothing left to warn about — the gap was closed by
    choosing a better panel rather than by writing a caveat. A flag that fired
    anyway would train a reader to skip it, which would cost us the cases where
    it matters.
    """
    assert "rs112766203" in PANELS["amp_tier1"][1]
    assert "P7_DPYD_ANCESTRY_PANEL_MISMATCH" not in _rule_ids(
        _interp(panel="amp_tier1")
    )


def test_it_names_the_alleles_that_were_not_read_and_why_they_matter() -> None:
    # "Some alleles were missed" is not actionable. Which ones, and why they
    # matter in THIS population, is.
    flag = _flag(_interp(), "P7_DPYD_ANCESTRY_PANEL_MISMATCH")
    assert "c.2279C>T" in flag.detail
    assert "c.557A>G" in flag.detail
    assert "GenomeIndia" in flag.detail
    assert "38886557" in flag.citations


def test_it_does_not_change_the_call_or_recommend_another_test() -> None:
    """Disclosure, not direction.

    Which assay to run is a clinical and procurement decision. Recommending one
    would move this package from informing clinical management to directing it,
    which is a different regulatory class.
    """
    plain = _interp(population="European")
    flagged = _interp(population=SAS)
    assert flagged.phenotype == plain.phenotype
    assert flagged.cpic_recommendation_action == plain.cpic_recommendation_action
    assert flagged.diplotype == plain.diplotype

    text = _flag(flagged, "P7_DPYD_ANCESTRY_PANEL_MISMATCH").detail.lower()
    for directive in ("we recommend", "you should order", "order a", "repeat the test"):
        assert directive not in text


# ---------------------------------------------------------------------------
# P8 — what is true whatever the assay found
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("panel", sorted(PANELS))
@pytest.mark.parametrize("population", ["European", SAS])
def test_the_monitoring_obligation_is_disclosed_on_every_result(
    panel, population
) -> None:
    assert "P8_DPYD_MONITORING_UNCONDITIONAL" in _rule_ids(
        _interp(panel=panel, population=population)
    )


def test_it_fires_on_the_cleanest_possible_result() -> None:
    """Especially on the clean ones.

    An all-negative panel is the report most likely to be read as an all-clear,
    and it is the one where an unconditional monitoring statement does the most
    work.
    """
    interp = _interp()
    assert interp.is_panel_negative
    assert "P8_DPYD_MONITORING_UNCONDITIONAL" in _rule_ids(interp)


def test_it_quotes_the_regulators_rather_than_asserting() -> None:
    flag = _flag(_interp(), "P8_DPYD_MONITORING_UNCONDITIONAL")
    assert "not a complete listing" in flag.detail
    assert "severe toxicity may still occur" in flag.detail
    assert "FDA" in flag.detail and "EMA" in flag.detail


def test_it_does_not_tell_anyone_how_to_monitor() -> None:
    # Stating that an obligation is not discharged is informing. Specifying the
    # schedule, the test or the threshold would be directing care.
    detail = _flag(_interp(), "P8_DPYD_MONITORING_UNCONDITIONAL").detail.lower()
    for directive in ("every cycle", "weekly", "check ", "measure ", "order "):
        assert directive not in detail
