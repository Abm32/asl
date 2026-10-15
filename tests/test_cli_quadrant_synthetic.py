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


def _wide_table(path, rows):
    """A cpic4_plus_sas table: the four actionable loci plus four CPIC
    normal-function ones."""
    loci = list(PANELS["cpic4_plus_sas"][1])
    cols = ["anonymised_id", "max_toxicity_grade", *loci]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for pid, grade, calls in rows:
            row = {r: "0" for r in loci}
            row.update(calls)
            row["anonymised_id"] = pid
            row["max_toxicity_grade"] = grade
            w.writerow(row)
    return path


def test_a_normal_function_call_is_not_panel_positive(tmp_path, capsys):
    """`*9A` is CPIC normal function and about a quarter of Central/South
    Asian chromosomes. Counting it as panel-positive pulls those patients out
    of the panel-negative denominator, which is the cell this table exists to
    measure — and the comparators printed underneath define their own group by
    the absence of an ACTIONABLE variant, so the table would contradict the
    rate it is being compared to, for the most common allele in the gene.
    """
    path = _wide_table(tmp_path / "w.csv", [
        ("P1", "4", {}),                       # wild-type, toxic
        ("P2", "1", {"rs1801265": "1"}),       # *9A het, not toxic
        ("P3", "4", {"rs1801265": "1"}),       # *9A het, toxic
    ])
    assert main([
        "quadrant", "--table", str(path), "--panel", "cpic4_plus_sas"
    ]) == 0
    out = capsys.readouterr().out
    assert "panel POSITIVE                   0           0" in out
    assert "panel NEGATIVE                   2           1" in out
    assert "2/3 (66.7%) of panel-NEGATIVE" in out


def test_an_actionable_call_is_still_panel_positive(tmp_path, capsys):
    """The guard must not swallow the alleles that do count. `*2A` is CPIC
    no function and has to land in the positive row."""
    path = _wide_table(tmp_path / "a.csv", [
        ("P1", "4", {"rs3918290": "1"}),       # *2A het, toxic
        ("P2", "1", {}),                       # wild-type, not toxic
    ])
    assert main([
        "quadrant", "--table", str(path), "--panel", "cpic4_plus_sas"
    ]) == 0
    out = capsys.readouterr().out
    assert "panel POSITIVE                   1           0" in out
    assert "panel NEGATIVE                   0           1" in out


def test_a_no_call_at_a_normal_function_locus_does_not_exclude(tmp_path, capsys):
    """A blank at a normal-function locus cannot change the classification,
    so discarding the patient costs a row and buys nothing."""
    path = _wide_table(tmp_path / "n.csv", [
        ("P1", "4", {"rs1801160": ""}),        # no call at *6, normal function
        ("P2", "1", {}),
    ])
    assert main([
        "quadrant", "--table", str(path), "--panel", "cpic4_plus_sas"
    ]) == 0
    out = capsys.readouterr().out
    assert "n = 2 patients" in out
    assert "panel NEGATIVE                   1           1" in out


def test_no_panel_label_sells_normal_function_loci_as_population_coverage():
    """The label is what the report attributes its conclusions to, and it is
    printed inside the sentence naming what was NOT assayed. A preset whose
    extra loci are all CPIC normal function must not describe them in terms a
    reader takes for wider coverage for their population.

    `cpic4_plus_sas` said "+ 4 alleles common in South Asians" until
    2026-10-14. Both halves were wrong: the four add no actionable coverage,
    and three of them are not South Asian-enriched (*9A is African-enriched,
    *5 is lowest in South Asians, c.496A>G is level with Europeans).

    `full` is deliberately exempt from the population clause: "All 32 DPYD
    loci" is a count, true as written, and claims nothing about ancestry.
    """
    from asl.residual import actionable_loci

    population_words = ("south asian", "indian", "common in", "ancestry")
    for key, (label, loci) in PANELS.items():
        if len(actionable_loci("DPYD", loci)) == len(loci):
            continue  # every locus can change a dose; nothing to disclaim
        lowered = label.lower()
        if any(w in lowered for w in population_words):
            assert "no added actionable coverage" in lowered, (
                f"panel {key!r} describes non-actionable loci in population "
                f"terms without saying they add no coverage: {label!r}"
            )


def test_the_sas_preset_label_states_it_adds_no_actionable_coverage():
    label, loci = PANELS["cpic4_plus_sas"]
    from asl.residual import actionable_loci

    assert len(actionable_loci("DPYD", loci)) == 4, (
        "the four extras must still be the CPIC-normal ones; if a later edit "
        "repoints this preset at actionable loci, rewrite the label too"
    )
    assert "no added actionable coverage" in label
    assert "common in South Asians" not in label
