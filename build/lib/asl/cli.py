"""Command line interface.

Designed for a clinic room: no network, no database, no configuration. Reads a
CSV or a few flags, prints a report.

    asl demo                          # the walkthrough patient
    asl demo --panel cpic4            # panel-negative on the CPIC four
    asl interpret --csv patients.csv   # a batch from a spreadsheet
    asl panels                        # what panel presets exist
    asl quadrant --csv patients.csv    # the four-quadrant table for a cohort

``quadrant`` is the one to run against a partner site's spreadsheet: it answers
"how many patients had grade 3+ toxicity while negative on every panel variant?"
which is the question a retrospective exists to settle.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from . import render
from .interpret import CPIC_ACTIONABLE, CPIC_ALLELE, InterpretationRefused
from .models import AssayResult, Drug, Regimen, TestStatus, VariantCall, Zygosity
from .pipeline import build

# Panel presets. A panel is a *list of loci*, and choosing one is choosing what
# the report may conclude — so presets are named and explicit rather than
# inferred from whichever columns happen to be present in a spreadsheet.
PANELS: dict[str, tuple[str, tuple[str, ...]]] = {
    "cpic4": (
        "CPIC 4-variant DPYD panel",
        CPIC_ACTIONABLE,
    ),
    "cpic4_plus_sas": (
        "CPIC 4-variant panel + 4 alleles common in South Asians",
        CPIC_ACTIONABLE + ("rs1801265", "rs1801160", "rs1801159", "rs2297595"),
    ),
    "full": (
        "All 13 DPYD loci anukriti-pgx-core calls",
        tuple(CPIC_ALLELE.keys()),
    ),
}

_ZYG = {
    "0": Zygosity.HOM_REF, "hom_ref": Zygosity.HOM_REF, "wt": Zygosity.HOM_REF,
    "1": Zygosity.HET, "het": Zygosity.HET,
    "2": Zygosity.HOM_VAR, "hom_var": Zygosity.HOM_VAR, "hom": Zygosity.HOM_VAR,
}


def _assay(
    panel: str,
    genotypes: dict[str, Zygosity | None],
    regimen: Regimen,
    population: str | None,
    drug: Drug = Drug.CAPECITABINE,
) -> AssayResult:
    """Build an AssayResult from a panel preset plus observed genotypes.

    ``genotypes`` must contain an entry for **every** locus the panel covers.
    ``None`` means the assay covered the locus but returned no call, which
    becomes ``NOT_TESTED`` — not a negative. Every locus outside the panel
    becomes ``NOT_TESTED`` too. That mapping is the whole point: it is what stops
    an unread variant being reported as absent.

    There is deliberately no default. An earlier version treated an absent key as
    ``hom_ref``, which made "the assay returned no call" and "the assay returned
    a negative" the same value — so a blank spreadsheet cell was reported as
    "tested & negative". That is the precise conflation this package exists to
    prevent, and a default is what allowed it.

    Refuses when an actual genotype is supplied for a locus the chosen panel does
    not cover. That combination is contradictory — the spreadsheet says the locus
    was read, the panel says it was never assayed — and resolving it either way
    loses information: marking the locus NOT_TESTED discards an observed variant,
    and marking it TESTED silently widens the panel the report attributes its
    conclusions to.

    A ``None`` for an off-panel locus is **not** that contradiction and is
    accepted. It means an empty column, which observes nothing: it discards no
    genotype and widens no panel. Refusing it made
    ``fixtures/cohort_template.csv`` — the layout this package asks sites to fill
    in — refuse every row under the default panel, because the template carries
    the four South Asian columns for sites that can fill them and leaves them
    blank for sites that cannot.
    """
    label, loci = PANELS[panel]
    off_panel = sorted(r for r in set(genotypes) - set(loci) if genotypes[r] is not None)
    if off_panel:
        raise InterpretationRefused(
            "genotype supplied for "
            + ", ".join(f"{r} ({CPIC_ALLELE[r]})" for r in off_panel)
            + f", which panel '{panel}' does not cover. asl will not drop an "
            "observed genotype, and it will not silently widen the panel its "
            "conclusions are attributed to. Choose a panel that covers these "
            "loci (try `asl panels`) or remove the columns."
        )
    unstated = sorted(set(loci) - set(genotypes))
    if unstated:
        raise InterpretationRefused(
            "no entry for "
            + ", ".join(f"{r} ({CPIC_ALLELE[r]})" for r in unstated)
            + f", which panel '{panel}' does cover. State a zygosity, or None "
            "if the assay returned no call at that locus. asl has no default: "
            "silently assuming a negative is how an unread locus becomes a "
            "clean result."
        )

    no_call = [r for r in loci if genotypes[r] is None]
    calls: list[VariantCall] = []
    for rsid in loci:
        zyg = genotypes[rsid]
        calls.append(
            VariantCall(
                rsid=rsid,
                status=TestStatus.NOT_TESTED if zyg is None else TestStatus.TESTED,
                zygosity=zyg,
                allele_label=CPIC_ALLELE[rsid],
            )
        )
    for rsid in CPIC_ALLELE:
        if rsid not in loci:
            calls.append(
                VariantCall(
                    rsid=rsid,
                    status=TestStatus.NOT_TESTED,
                    allele_label=CPIC_ALLELE[rsid],
                )
            )
    if no_call:
        # The panel name is what the report attributes its conclusions to, so it
        # must not claim coverage the data did not deliver.
        label = (
            f"{label} ({len(loci) - len(no_call)} of {len(loci)} loci returned "
            f"a call)"
        )
    return AssayResult(
        gene="DPYD",
        panel_name=label,
        calls=tuple(calls),
        regimen=regimen,
        population=population,
        drug=drug,
    )


def _all_negative(panel: str) -> dict[str, Zygosity | None]:
    """Every locus in a panel, explicitly tested and negative.

    Used by the demo and by tests. Written out rather than defaulted, because a
    default is what let a blank cell become a negative result.
    """
    return {rsid: Zygosity.HOM_REF for rsid in PANELS[panel][1]}


def cmd_panels(_args: argparse.Namespace) -> int:
    for key, (label, loci) in PANELS.items():
        print(f"{key:<16} {label}")
        print(f"{'':16} {len(loci)} loci: {', '.join(CPIC_ALLELE[r] for r in loci)}")
        print()
    return 0


def cmd_demo(args: argparse.Namespace) -> int:
    """The walkthrough patient: panel-negative, then the same patient sequenced."""
    pop = "Central/South Asian"
    drug = Drug(args.drug)
    if args.panel == "cpic4":
        print("PATIENT: 52F, Kerala, stage III colon adenocarcinoma, adjuvant CAPOX.")
        print("LAB: 4-variant DPYD PCR panel. All four negative.\n")
        assay = _assay("cpic4", _all_negative("cpic4"), Regimen.SYSTEMIC, pop, drug)
    else:
        print("PATIENT: same patient, with research sequencing added.")
        print("FOUND: *9A heterozygous, *6 heterozygous.\n")
        assay = _assay(
            "cpic4_plus_sas",
            _all_negative("cpic4_plus_sas")
            | {"rs1801265": Zygosity.HET, "rs1801160": Zygosity.HET},
            Regimen.SYSTEMIC,
            pop,
            drug,
        )
    print(render(build(assay)))
    return 0


def _row_genotypes(row: dict[str, str]) -> dict[str, Zygosity | None]:
    """Read the rsID columns of one CSV row.

    A column that is present but blank maps to ``None`` — the assay covered the
    locus and returned no call — which ``_assay`` renders as NOT_TESTED. A column
    that is absent entirely is not in the returned map at all, and ``_assay``
    refuses if the chosen panel covers it. Blank and absent are different facts
    and are kept different.
    """
    out: dict[str, Zygosity | None] = {}
    for key, value in row.items():
        k = (key or "").strip()
        if k not in CPIC_ALLELE:
            continue
        raw = (value or "").strip()
        if not raw:
            out[k] = None
            continue
        z = _ZYG.get(raw.lower())
        if z is None:
            raise ValueError(
                f"{k}: cannot read genotype {value!r}. Use 0/1/2, "
                f"hom_ref/het/hom_var, or wt/het/hom. Leave the cell blank "
                f"only if the assay returned no call at this locus."
            )
        out[k] = z
    return out


def _row_drug(row: dict[str, str], default: str) -> str:
    """Read the drug column, accepting either name the template has used.

    ``fixtures/cohort_template.csv`` calls it ``fluoropyrimidine_drug`` because
    that is unambiguous on a spreadsheet a site fills in by hand; the API and the
    flags call it ``drug``. Both are accepted rather than silently falling back to
    the default, because falling back would label a 5-FU patient's report with
    capecitabine and say nothing about it.
    """
    for key in ("drug", "fluoropyrimidine_drug"):
        value = (row.get(key) or "").strip()
        if value:
            return value.lower()
    return default


def cmd_interpret(args: argparse.Namespace) -> int:
    path = Path(args.csv)
    rows = list(csv.DictReader(path.open()))
    if not rows:
        print(f"{path}: no rows", file=sys.stderr)
        return 1
    for i, row in enumerate(rows, 1):
        pid = row.get("anonymised_id") or row.get("patient_id") or f"row{i}"
        print(f"\n{'#' * 76}\n# PATIENT {pid}\n{'#' * 76}")
        try:
            assay = _assay(
                args.panel,
                _row_genotypes(row),
                Regimen(row.get("regimen", "unknown").strip().lower() or "unknown"),
                row.get("population") or args.population,
                Drug(_row_drug(row, args.drug)),
            )
            print(render(build(assay)))
        except (InterpretationRefused, ValueError) as exc:
            print(f"REFUSED: {exc}")
    return 0


def cmd_quadrant(args: argparse.Namespace) -> int:
    """The four-quadrant table: panel status x toxicity.

    This is the question a site retrospective exists to answer. The cell that
    matters is panel-negative with grade >=3 toxicity: those are the patients a
    guideline-conformant test reassured and who then toxified anyway.

    A patient with any no-call among the panel loci is **excluded**, not counted
    as negative. Counting them as negative would inflate exactly the cell the
    analysis turns on, and in the direction that flatters our own thesis.
    """
    path = Path(args.csv)
    rows = list(csv.DictReader(path.open()))
    _, loci = PANELS[args.panel]

    # A spreadsheet can list the same patient twice — a corrected row appended
    # below the original is how a site actually edits a CSV. Counting both would
    # double-count one patient, and in the observed ledger case put them in two
    # contradictory cells. Later rows win, and the count is reported rather than
    # silently applied.
    deduped: dict[str, dict[str, str]] = {}
    unidentified = 0
    for row in rows:
        pid = (row.get("anonymised_id") or row.get("patient_id") or "").strip()
        if not pid:
            unidentified += 1
            continue
        deduped[pid] = row
    superseded = len(rows) - unidentified - len(deduped)
    rows = list(deduped.values())

    cells = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
    missing = 0
    incomplete = 0
    for row in rows:
        g = _row_genotypes(row)
        unread = [r for r in loci if g.get(r) is None]
        if unread:
            incomplete += 1
            continue
        positive = any(g[r] is not Zygosity.HOM_REF for r in loci)
        grade = (row.get("max_toxicity_grade") or "").strip()
        if not grade.isdigit():
            missing += 1
            continue
        cells[(positive, int(grade) >= 3)] += 1

    n = sum(cells.values())
    print(f"Four-quadrant table — panel: {PANELS[args.panel][0]}")
    print(f"n = {n} patients with a readable toxicity grade"
          + (f" ({missing} skipped, grade not recorded)" if missing else ""))
    if superseded:
        print(
            f"    {superseded} row(s) superseded: the same anonymised_id appeared "
            f"more than once, and the last occurrence was used. Counting both "
            f"would count one patient twice."
        )
    if unidentified:
        print(
            f"    {unidentified} row(s) skipped: no anonymised_id, so duplicates "
            f"cannot be detected for them."
        )
    if incomplete:
        print(
            f"    {incomplete} excluded: at least one panel locus returned no "
            f"call. Not counted as negative — an unread locus is not a clean "
            f"result, and counting it as one would inflate the panel-negative "
            f"cell this table exists to measure."
        )
    print()
    print(f"{'':22}{'grade >=3':>12}{'grade <3':>12}")
    print(f"{'panel POSITIVE':<22}{cells[(True, True)]:>12}{cells[(True, False)]:>12}")
    print(f"{'panel NEGATIVE':<22}{cells[(False, True)]:>12}{cells[(False, False)]:>12}")
    print()
    neg_total = cells[(False, True)] + cells[(False, False)]
    if neg_total:
        pct = 100.0 * cells[(False, True)] / neg_total
        print(
            f"{cells[(False, True)]}/{neg_total} ({pct:.1f}%) of panel-NEGATIVE "
            f"patients had grade >=3 toxicity."
        )
        print(
            "  Published comparators for patients without an actionable DPYD "
            "variant:\n"
            "    22.7% (231/1018) systemic chemotherapy   Henricks 2018, PMID 30348537\n"
            "    13.6% (105/771)  incl. chemoradiation    Lunenburg 2018, PMID 30361102"
        )
    print()
    print("NOTE: this is a descriptive table, not a hypothesis test. It does not")
    print("establish that any allele causes toxicity, and asl does not claim so.")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="asl", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("demo", help="run the walkthrough patient")
    p.add_argument("--panel", default="cpic4", choices=sorted(PANELS))
    p.add_argument("--drug", default="capecitabine", choices=[d.value for d in Drug])
    p.set_defaults(func=cmd_demo)

    p = sub.add_parser("panels", help="list panel presets")
    p.set_defaults(func=cmd_panels)

    p = sub.add_parser("interpret", help="interpret a CSV of patients")
    p.add_argument("--csv", required=True)
    p.add_argument("--panel", default="cpic4", choices=sorted(PANELS))
    p.add_argument("--population", default="Central/South Asian")
    p.add_argument("--drug", default="capecitabine", choices=[d.value for d in Drug])
    p.set_defaults(func=cmd_interpret)

    p = sub.add_parser("quadrant", help="four-quadrant table for a cohort CSV")
    p.add_argument("--csv", required=True)
    p.add_argument("--panel", default="cpic4", choices=sorted(PANELS))
    p.set_defaults(func=cmd_quadrant)

    args = ap.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
