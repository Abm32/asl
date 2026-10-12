"""Did this result arrive in time to change the first dose?

The binding constraint on DPYD testing is not assay turnaround. Virk, Challa &
Sadri (Am J Clin Pathol 2026;166(3), PMID 42803377) name it directly — the
problem is "results returning after the first treatment decision" — and in
their cohort testing was ordered for 86.9% of patients while only **44.9%** had
a result before cycle 1. UPenn's prospective series found 57.4% (Tuteja 2025,
PMID 40773711).

A result landing on day 8 into a decision taken on day 3 is not a safety
intervention; it is a post-hoc audit document. The report should say which one
it is rather than print two dates in different places and leave the reader to
subtract.
"""

from __future__ import annotations

import dataclasses
from datetime import UTC, date, datetime, timedelta

import pytest

from asl.cli import _all_negative, _assay, _with_start_date
from asl.interpret import InterpretationRefused
from asl.models import Drug, Regimen
from asl.pipeline import build


def _today() -> date:
    """UTC, because the pipeline compares against ``datetime.now(UTC).date()``.

    Using the local date here would make these tests flake for anyone running
    them within a few hours of midnight in a non-UTC timezone — and the
    failure would look like a logic bug in the flag rather than a timezone
    mismatch in the test.
    """
    return datetime.now(UTC).date()


def _interp(start=None):
    assay = _assay(
        "cpic4", _all_negative("cpic4"), Regimen.SYSTEMIC, "European",
        Drug.CAPECITABINE,
    )
    assay = dataclasses.replace(assay, treatment_start_date=start)
    return build(assay)


def _p9(interp):
    return next(f for f in interp.flags if f.rule_id.startswith("P9_"))


def test_no_recorded_date_is_not_the_same_as_in_time() -> None:
    """The commonest case, and the one most easily mistaken for success.

    A programme that never records the planned start date cannot measure the
    gap it most needs to close, so silence here is reported as silence.
    """
    flag = _p9(_interp(None))
    assert flag.rule_id == "P9_DPYD_TIMELINESS_UNKNOWN"
    assert "cannot be determined" in flag.summary
    assert "not the same as arriving in time" in flag.detail
    assert "44.9%" in flag.detail


def test_a_result_that_beat_the_start_date_says_so_with_the_margin() -> None:
    flag = _p9(_interp(_today() + timedelta(days=5)))
    assert flag.rule_id == "P9_DPYD_RESULT_BEFORE_TREATMENT"
    assert "5 day(s) before" in flag.summary


def test_a_late_result_says_the_first_dose_was_decided_without_it() -> None:
    flag = _p9(_interp(_today() - timedelta(days=3)))
    assert flag.rule_id == "P9_DPYD_RESULT_AFTER_TREATMENT_START"
    assert "first dose was decided without it" in flag.summary
    # What it changes, stated: the recommendation now bears on later cycles,
    # and any toxicity already seen was seen at an unadjusted dose.
    assert "subsequent cycles" in flag.detail
    assert "unadjusted dose" in flag.detail


def test_a_result_on_the_start_date_counts_as_late() -> None:
    # The decision is taken on or before the morning of day one. A result
    # timestamped the same day cannot be assumed to have preceded it, and
    # assuming it did would flatter us in exactly the direction that matters.
    flag = _p9(_interp(_today()))
    assert flag.rule_id == "P9_DPYD_RESULT_AFTER_TREATMENT_START"


def test_it_never_advises_delaying_treatment() -> None:
    """Whether to wait for a pending result is a clinical judgement.

    The Swiss SGPPT position that five working days is tolerable rests on
    fluoropyrimidines not being used for oncologic emergencies — a premise this
    package cannot verify for any given patient.
    """
    detail = _p9(_interp(_today() - timedelta(days=3))).detail.lower()
    assert "does not advise delaying" in detail
    for directive in ("delay treatment", "withhold the first", "postpone"):
        assert directive not in detail.replace("does not advise delaying", "")


def test_the_flag_and_the_report_timestamp_cannot_disagree() -> None:
    # Computed against the same timestamp the Interpretation carries, so the
    # flag text and the report's own "Generated" line are always consistent.
    interp = _interp(_today() - timedelta(days=2))
    assert interp.generated_at.date().isoformat() in _p9(interp).detail


def test_exactly_one_timeliness_flag_is_ever_emitted() -> None:
    for start in (None, _today(), _today() + timedelta(days=9)):
        flags = [f for f in _interp(start).flags if f.rule_id.startswith("P9_")]
        assert len(flags) == 1, start


class TestParsingTheDate:
    def _assay(self):
        return _assay(
            "cpic4", _all_negative("cpic4"), Regimen.SYSTEMIC, "European",
            Drug.CAPECITABINE,
        )

    def test_an_absent_date_is_simply_absent(self) -> None:
        for raw in (None, "", "   "):
            assert _with_start_date(self._assay(), raw).treatment_start_date is None

    def test_an_iso_date_is_read(self) -> None:
        got = _with_start_date(self._assay(), "2026-10-17")
        assert got.treatment_start_date == date(2026, 10, 17)

    @pytest.mark.parametrize("raw", ["17/10/2026", "Oct 17 2026", "tomorrow", "2026-13-01"])
    def test_an_unparseable_date_is_refused_not_dropped(self, raw) -> None:
        """Silently ignoring it would be the worst of both worlds.

        The report would then say the timeliness of this result was unknown,
        when in fact it was supplied and thrown away — the same conflation of
        "unknown" with "not recorded" that this package exists to prevent,
        turned on its own input.
        """
        with pytest.raises(InterpretationRefused) as exc:
            _with_start_date(self._assay(), raw)
        assert "YYYY-MM-DD" in str(exc.value)
        assert "supplied and discarded" in str(exc.value)
