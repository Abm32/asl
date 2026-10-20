"""A site's own panel must be expressible, and the variant list must stay live.

Two surfaces that exist for a partner laboratory rather than for a clinician.

``resolve_panel`` is the safety-relevant one. A panel determines what a report
may conclude, so a locus that cannot be placed must stop the run rather than
quietly shrink the panel underneath the conclusion.

``variants`` is the artifact a site maps its assay against. It is generated
from the pinned CPIC table at run time precisely so it cannot drift from the
engine, and the property worth pinning is the one a site gets wrong: a common
allele that changes no dose must not read as something worth testing for.
"""

from __future__ import annotations

import csv

import pytest

from asl.cli import PANELS, main, resolve_panel
from asl.interpret import InterpretationRefused
from asl.residual import actionable_loci


def test_a_site_panel_resolves_the_names_a_lab_report_actually_prints():
    # rsID, cDNA, bare star and prefixed star are four spellings of two loci.
    _, loci = resolve_panel("rs3918290,c.1679T>G,*9A,DPYD*6")
    assert loci == ("rs3918290", "rs55886062", "rs1801265", "rs1801160")


def test_the_hapb3_tag_and_causal_variants_stay_distinct_in_a_site_panel():
    # The whole point of the hand-pinned HapB3 synonyms: a site that genotypes
    # the exonic tag has not genotyped the causal variant, and a panel that
    # conflated them would report an inferred call as a direct one.
    _, tag = resolve_panel("HapB3")
    _, causal = resolve_panel("c.1129-5923C>G")
    assert tag == ("rs56038477",)
    assert causal == ("rs75017182",)


def test_a_locus_named_twice_is_carried_once():
    _, loci = resolve_panel("HapB3,c.1236G>A,rs56038477")
    assert loci == ("rs56038477",)


def test_an_unplaceable_locus_refuses_rather_than_shrinking_the_panel():
    with pytest.raises(InterpretationRefused) as exc:
        resolve_panel("rs3918290,DPYD*99")
    assert "DPYD*99" in str(exc.value)
    # The refusal must not leave the caller holding a usable 1-locus panel.
    assert "rs3918290" not in str(exc.value).split("Run `asl variants`")[0][:40]


def test_an_empty_panel_refuses():
    with pytest.raises(InterpretationRefused):
        resolve_panel(" , ")


def test_preset_names_still_win_over_locus_parsing():
    assert resolve_panel("cpic4") == PANELS["cpic4"]


def test_a_normal_function_locus_on_a_site_panel_cannot_make_a_patient_positive():
    """The quadrant split must not move with panel size.

    `*9A` is a quarter of Central/South Asian chromosomes and CPIC calls it
    normal function. A site that adds it to its assay must not thereby turn
    every second patient into a "panel positive", which would contradict the
    published comparators the quadrant table prints beneath itself.
    """
    assert actionable_loci("DPYD", ("rs1801265", "rs1801159", "rs2297595")) == ()
    assert actionable_loci("DPYD", ("rs3918290",)) == ("rs3918290",)


def _cohort(path, *, panel_loci, grades):
    cols = ["anonymised_id", "max_toxicity_grade", *panel_loci]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for i, grade in enumerate(grades):
            row = {r: "0" for r in panel_loci}
            row["anonymised_id"] = f"SYNTH-{i:05d}"
            row["max_toxicity_grade"] = str(grade)
            w.writerow(row)
    return path


def test_quadrant_accepts_a_site_panel_and_names_it(tmp_path, capsys):
    spec = "*2A,*13,c.2846A>T,HapB3,*9A,*5,*6"
    _, loci = resolve_panel(spec)
    path = _cohort(tmp_path / "c.csv", panel_loci=loci, grades=[3, 3, 1, 1])
    assert main(["quadrant", "--table", str(path), "--panel", spec]) == 0
    out = capsys.readouterr().out
    assert "site-defined panel, 7 loci" in out
    # All genotypes are hom_ref, so all four sit in the panel-NEGATIVE row.
    assert "2/4 (50.0%) of panel-NEGATIVE patients" in out


def test_a_bad_panel_prints_a_refusal_not_a_traceback(tmp_path, capsys):
    path = _cohort(tmp_path / "c.csv", panel_loci=("rs3918290",), grades=[1])
    bad = "rs3918290,NOTALOCUS"
    assert main(["quadrant", "--table", str(path), "--panel", bad]) == 1
    assert "REFUSED" in capsys.readouterr().err


def test_variants_lists_every_locus_the_engine_calls(capsys):
    from asl.interpret import CPIC_ALLELE

    assert main(["variants", "--format", "csv"]) == 0
    rows = list(csv.DictReader(capsys.readouterr().out.splitlines()))
    assert [r["rsid"] for r in rows] == list(CPIC_ALLELE)


def test_variants_marks_the_common_normal_function_trap(capsys):
    """`*9A` is the row a site is most likely to misread.

    It is the most frequent variant in the table for this population and it
    changes no dose. If `actionable` ever reads yes here, the table is telling
    a laboratory to buy an assay that cannot change a prescription.
    """
    assert main(["variants", "--format", "csv"]) == 0
    rows = {r["rsid"]: r for r in csv.DictReader(capsys.readouterr().out.splitlines())}

    star9a = rows["rs1801265"]
    assert star9a["actionable"] == "no"
    assert star9a["clinical_function"] == "Normal function"
    assert float(star9a["frequency"]) > 0.2

    star2a = rows["rs3918290"]
    assert star2a["actionable"] == "yes"
    assert float(star2a["frequency"]) < float(star9a["frequency"])
    assert star2a["citations_pmid"]


def test_every_variant_row_carries_a_citation_or_says_it_cannot(capsys):
    assert main(["variants", "--format", "csv"]) == 0
    for row in csv.DictReader(capsys.readouterr().out.splitlines()):
        assert row["citations_pmid"] or "UNCITABLE" in row["clinical_function"]
