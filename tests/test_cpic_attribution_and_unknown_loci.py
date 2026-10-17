"""Two things the payload must carry, and one it must not pretend to know.

Why this module exists
----------------------

**CPIC attribution.** ``Interpretation`` has carried ``cpic_recommendation_text``
and ``cpic_activity_score_range`` since the recommendation table landed, but the
API response did not copy them out. A caller therefore received the action code
(``cpic_action``) with nothing to attribute it to. An unattributed dose is this
package recommending rather than CPIC informing, which is a different regulatory
class and not a question of style — so the quoted sentence and its score range
travel with the code or the payload is not reportable.

**Unknown loci.** ``Diplotype.unknown_loci`` and the ``UNKNOWN_LOCUS_CARRIES_VARIANT``
refusal landed in ``pgx-core`` *after* the 0.11.0 release, and 0.11.0 is the exact
pin here. So the engine in use cannot answer the question at all. The failure to
avoid is ``getattr(result, "unknown_loci", ())`` — which hands a consumer an empty
list, indistinguishable from "the engine looked and found none". That is the same
silent, falsely reassuring collapse that ``require_hapb3_capable_engine`` exists
to stop, and the same class as the defect where the engine dropped unknown loci
without saying so.

Three states, kept distinct end to end:

* ``None`` — this engine cannot report unknown loci. Ignorance.
* ``()`` — it reported none. A negative result.
* non-empty — it reported these.

A hard refusal is deliberately *not* used: every released engine would fail it
and the endpoint would stop answering. The gap is disclosed instead, and it
closes by itself once the pin moves past 0.11.0.
"""

from __future__ import annotations

from asl import interpret
from asl.interpret import CPIC_ALLELE, unknown_loci_from
from asl.models import (
    AssayResult,
    Drug,
    Regimen,
    TestStatus,
    VariantCall,
    Zygosity,
)
from asl.pipeline import build

# *2A (rs3918290), a CPIC no-function allele — the least ambiguous way to reach
# a phenotype that carries a real dose recommendation.
NO_FUNCTION = "rs3918290"


def _assay(found: dict[str, Zygosity]) -> AssayResult:
    calls = [
        VariantCall(
            rsid=rsid,
            status=TestStatus.TESTED,
            zygosity=found.get(rsid, Zygosity.HOM_REF),
            allele_label=CPIC_ALLELE[rsid],
        )
        for rsid in CPIC_ALLELE
    ]
    return AssayResult(
        gene="DPYD",
        panel_name="test panel",
        calls=tuple(calls),
        regimen=Regimen.SYSTEMIC,
        population="Central/South Asian",
        drug=Drug.FLUOROURACIL,
    )


class _Stub:
    """Stands in for the engine's ``Diplotype`` class and one of its results."""

    def __init__(self, **kw: object) -> None:
        for key, value in kw.items():
            setattr(self, key, value)


def test_cpic_dose_sentence_travels_with_the_action_code() -> None:
    """The action code is never the only thing the payload says about the dose."""
    interp = build(_assay({NO_FUNCTION: Zygosity.HET}))

    assert interp.cpic_recommendation_action, "no action code"
    # The code alone is asl asserting a dose. The quoted sentence and the score
    # range are what make it CPIC's.
    assert interp.cpic_recommendation_text.strip(), (
        "action code present with no CPIC sentence to attribute it to"
    )
    assert interp.cpic_activity_score_range.strip(), (
        "CPIC dose sentence present with no activity-score range it applies to"
    )


def test_unknown_loci_is_none_when_the_engine_cannot_report_it(monkeypatch) -> None:
    """Ignorance must not serialise as a negative result.

    This is the state of the pinned 0.11.0 engine today.
    """
    monkeypatch.setattr(interpret, "Diplotype", _Stub)  # no `unknown_loci`
    assert unknown_loci_from(_Stub()) is None


def test_unknown_loci_is_empty_when_the_engine_reports_none(monkeypatch) -> None:
    class _Capable:
        unknown_loci: tuple[str, ...] = ()

    monkeypatch.setattr(interpret, "Diplotype", _Capable)
    assert unknown_loci_from(_Capable()) == ()


def test_unknown_loci_is_carried_through_when_the_engine_reports_some(
    monkeypatch,
) -> None:
    class _Capable:
        unknown_loci: tuple[str, ...] = ()

    monkeypatch.setattr(interpret, "Diplotype", _Capable)
    assert unknown_loci_from(_Stub(unknown_loci=("rs1801265",))) == ("rs1801265",)


def test_pipeline_matches_the_live_engine_s_capability() -> None:
    """The field tracks the installed engine, and is never quietly dropped.

    On the 0.11.0 pin this is ``None``; past it, a tuple. Both are correct —
    what would be wrong is the pipeline losing the field on the way out.
    """
    interp = build(_assay({NO_FUNCTION: Zygosity.HET}))

    if hasattr(interpret.Diplotype, "unknown_loci"):
        assert isinstance(interp.unknown_loci, tuple)
    else:
        assert interp.unknown_loci is None
