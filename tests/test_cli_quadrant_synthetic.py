"""The quadrant table must say when the cohort it counts was generated.

The panel-negative toxicity rate prints directly above two real published
comparators. Without a marker, a screenshot of that pair reads as a finding
about a real cohort, which is the one thing synthetic data must never support.
"""

from __future__ import annotations

import csv

from asl.cli import PANELS, main

COLUMNS = [
    "anonymised_id", "max_toxicity_grade", *PANELS["cpic4"][1],
]


def _table(path, ids):
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for i, pid in enumerate(ids):
            row = {r: "0" for r in PANELS["cpic4"][1]}
            row["anonymised_id"] = pid
            row["max_toxicity_grade"] = "3" if i % 2 else "1"
            w.writerow(row)
    return path


def test_an_all_synthetic_cohort_is_banner_marked(tmp_path, capsys):
    path = _table(tmp_path / "c.csv", ["SYNTH-00001", "SYNTH-00002"])
    assert main(["quadrant", "--table", str(path), "--panel", "cpic4"]) == 0
    out = capsys.readouterr().out
    assert "SYNTHETIC DATA" in out
    # Above the table, not below it, and above the published comparators.
    assert out.index("SYNTHETIC DATA") < out.index("Four-quadrant table")
    assert out.index("SYNTHETIC DATA") < out.index("Published comparators")


def test_a_real_cohort_gets_no_banner(tmp_path, capsys):
    path = _table(tmp_path / "c.csv", ["PT-1", "PT-2"])
    assert main(["quadrant", "--table", str(path), "--panel", "cpic4"]) == 0
    assert "SYNTHETIC DATA" not in capsys.readouterr().out


def test_one_real_patient_removes_the_banner(tmp_path, capsys):
    """A mixed table is not a synthetic table, and must not be labelled one."""
    path = _table(tmp_path / "c.csv", ["SYNTH-00001", "PT-2"])
    assert main(["quadrant", "--table", str(path), "--panel", "cpic4"]) == 0
    assert "SYNTHETIC DATA" not in capsys.readouterr().out
