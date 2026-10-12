"""What the top of the page says, and whose voice it says it in.

``test_safety_properties`` pins what the report may never say. This file pins
what it must say, and where. The two are different obligations: a report can
satisfy every refusal in this package and still be useless, because the answer
is assembled on page two out of a diplotype, an enum and a paragraph of
residual risk that a clinician between patients will not reach.

The band exists because the output is the product. A result nobody reads to the
end has to be correct in its first six lines.
"""

from __future__ import annotations

import pytest

from asl.cli import _all_negative, _assay
from asl.models import Drug, Regimen, Zygosity
from asl.pipeline import build
from asl.report import render

POP = "Central/South Asian"


def _report(genotypes=None, panel="cpic4", **kw) -> str:
    g = _all_negative(panel) | (genotypes or {})
    assay = _assay(panel, g, kw.pop("regimen", Regimen.SYSTEMIC), POP,
                   kw.pop("drug", Drug.CAPECITABINE))
    return render(build(assay), **kw)


# ---------------------------------------------------------------------------
# The dose, in CPIC's voice.
# ---------------------------------------------------------------------------


def test_the_dose_is_cpics_sentence_quoted_not_an_enum() -> None:
    """``REDUCE_50PCT`` is auditable and unreadable; the sentence is both.

    A reader cannot check a closed-enum token against a guideline. They can
    check a quoted sentence against PMID 29152729 in a minute, which is the
    whole basis on which this report claims to inform rather than decide.
    """
    text = _report({"rs3918290": Zygosity.HET})
    assert (
        '"Reduce starting dose by 50%, then titrate to toxicity or therapeutic '
        'drug monitoring."' in " ".join(text.split())
    )
    assert "PMID 29152729" in text


def test_the_dose_is_attributed_to_cpic_in_the_same_breath() -> None:
    """Attribution is the regulatory posture, not a courtesy.

    CDSCO classifies medical-device software by the voice a claim is made in.
    The same sentence quoted from CPIC and attributed to CPIC is a guideline
    informing a clinician; written in asl's own voice it is asl driving a
    dosing decision, which is a different device class. So the attribution and
    the disclaimer travel with the sentence, in the band, rather than only in
    the boundary section at the foot of the page.
    """
    flat = " ".join(_report({"rs3918290": Zygosity.HET}).split())
    assert "This is CPIC's recommendation." in flat
    assert "asl does not recommend doses and has not modified this one." in flat


@pytest.mark.parametrize(
    ("genotypes", "expect"),
    [
        ({}, "Use label-recommended dosage and administration."),
        ({"rs3918290": Zygosity.HET}, "Reduce starting dose by 50%"),
        ({"rs3918290": Zygosity.HOM_VAR}, "Avoid capecitabine."),
    ],
    ids=["normal", "intermediate", "poor"],
)
def test_every_phenotype_reaches_the_band_with_cpics_own_wording(
    genotypes, expect
) -> None:
    assert expect in " ".join(_report(genotypes).split())


def test_the_answer_band_precedes_everything_else() -> None:
    text = _report({"rs3918290": Zygosity.HET})
    assert text.index("ANSWER") < text.index("DPYD PHARMACOGENOMIC INTERPRETATION")
    assert text.index("DOSE") < text.index("CPIC CALL")


# ---------------------------------------------------------------------------
# The coverage strip.
# ---------------------------------------------------------------------------


def test_the_strip_counts_found_read_and_unassayed_separately() -> None:
    """Three states, never two.

    Collapsing "read, no variant" and "not assayed" into a single negative is
    the defect this package exists to prevent, and the strip is the most
    compressed place it could happen. A guideline-conformant 4-variant panel
    reads 4 of the 19 loci asl can call; the strip is where that is visible
    without reading a paragraph.
    """
    text = _report({"rs3918290": Zygosity.HET})
    flat = " ".join(text.split())
    assert "X found (1)" in flat
    assert "# read, no variant (3)" in flat
    assert ". NOT assayed (15)" in flat
    assert "19 DPYD loci asl can call" in flat


def test_the_strip_is_in_a_fixed_locus_order_across_patients() -> None:
    """Position means the same variant on every report, so two can be compared."""
    import re

    def strip(text: str) -> str:
        return re.search(r"COVERAGE  \[([X#.]+)\]", text).group(1)

    a = strip(_report({"rs3918290": Zygosity.HET}))
    b = strip(_report({"rs67376798": Zygosity.HET}))
    assert len(a) == len(b) == 19
    assert a.index("X") != b.index("X")  # different locus, different position


def test_a_wider_panel_shows_more_read_loci_in_the_strip() -> None:
    narrow = " ".join(_report(panel="cpic4").split())
    wide = " ".join(_report(panel="amp_tier1").split())
    assert "# read, no variant (4)" in narrow
    assert "# read, no variant (8)" in wide


# ---------------------------------------------------------------------------
# Authorisation.
# ---------------------------------------------------------------------------


def test_an_unsigned_report_says_so_rather_than_leaving_a_blank_line() -> None:
    """A blank line under "Signed" reads as a printing fault, not as unsigned.

    asl has no signatory and cannot acquire one, so the honest output is the
    statement that nobody authorised this — which is also what keeps the output
    from being mistaken for a laboratory report it is not.
    """
    flat = " ".join(_report().split())
    assert "UNSIGNED." in flat
    assert "is not a laboratory report" in flat


def test_a_supplied_signatory_is_recorded_verbatim_and_marked_unverified() -> None:
    flat = " ".join(
        _report(signatory="Dr A Nair, MD (Pathology), KMC Reg. 12345").split()
    )
    assert "Authorised by: Dr A Nair, MD (Pathology), KMC Reg. 12345" in flat
    assert "not verified by asl" in flat
    assert "post-graduate qualification in pathology" in flat
    assert "UNSIGNED" not in flat
