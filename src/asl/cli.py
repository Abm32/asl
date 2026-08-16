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
from .models import AssayResult, Regimen, TestStatus, VariantCall, Zygosity
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
    genotypes: dict[str, Zygosity],
    regimen: Regimen,
    population: str | None,
) -> AssayResult:
    """Build an AssayResult from a panel preset plus observed genotypes.

    Every locus in the panel becomes a TESTED call; every locus outside it
    becomes NOT_TESTED. That mapping is the whole point: it is what stops an
    untested variant being reported as absent.
    """
    label, loci = PANELS[panel]
    calls: list[VariantCall] = []
    for rsid in loci:
        calls.append(
            VariantCall(
                rsid=rsid,
                status=TestStatus.TESTED,
                zygosity=genotypes.get(rsid, Zygosity.HOM_REF),
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
    return AssayResult(
        gene="DPYD",
        panel_name=label,
        calls=tuple(calls),
        regimen=regimen,
        population=population,
    )


def cmd_panels(_args: argparse.Namespace) -> int:
    for key, (label, loci) in PANELS.items():
        print(f"{key:<16} {label}")
        print(f"{'':16} {len(loci)} loci: {', '.join(CPIC_ALLELE[r] for r in loci)}")
        print()
    return 0


def cmd_demo(args: argparse.Namespace) -> int:
    """The walkthrough patient: panel-negative, then the same patient sequenced."""
    pop = "Central/South Asian"
    if args.panel == "cpic4":
        print("PATIENT: 52F, Kerala, stage III colon adenocarcinoma, adjuvant CAPOX.")
        print("LAB: 4-variant DPYD PCR panel. All four negative.\n")
        assay = _assay("cpic4", {}, Regimen.SYSTEMIC, pop)
    else:
        print("PATIENT: same patient, with research sequencing added.")
        print("FOUND: *9A heterozygous, *6 heterozygous.\n")
        assay = _assay(
            "cpic4_plus_sas",
            {"rs1801265": Zygosity.HET, "rs1801160": Zygosity.HET},
            Regimen.SYSTEMIC,
            pop,
        )
    print(render(build(assay)))
    return 0


def _row_genotypes(row: dict[str, str]) -> dict[str, Zygosity]:
    out: dict[str, Zygosity] = {}
    for key, value in row.items():
        k = (key or "").strip()
        if k in CPIC_ALLELE and (value or "").strip():
            z = _ZYG.get(value.strip().lower())
            if z is None:
                raise ValueError(
                    f"{k}: cannot read genotype {value!r}. Use 0/1/2, "
                    f"hom_ref/het/hom_var, or wt/het/hom."
                )
            out[k] = z
    return out


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
    """
    path = Path(args.csv)
    rows = list(csv.DictReader(path.open()))
    _, loci = PANELS[args.panel]

    cells = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
    missing = 0
    for row in rows:
        g = _row_genotypes(row)
        positive = any(g.get(r, Zygosity.HOM_REF) is not Zygosity.HOM_REF for r in loci)
        grade = (row.get("max_toxicity_grade") or "").strip()
        if not grade.isdigit():
            missing += 1
            continue
        cells[(positive, int(grade) >= 3)] += 1

    n = sum(cells.values())
    print(f"Four-quadrant table — panel: {PANELS[args.panel][0]}")
    print(f"n = {n} patients with a readable toxicity grade"
          + (f" ({missing} skipped, grade not recorded)" if missing else ""))
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
    p.set_defaults(func=cmd_demo)

    p = sub.add_parser("panels", help="list panel presets")
    p.set_defaults(func=cmd_panels)

    p = sub.add_parser("interpret", help="interpret a CSV of patients")
    p.add_argument("--csv", required=True)
    p.add_argument("--panel", default="cpic4", choices=sorted(PANELS))
    p.add_argument("--population", default="Central/South Asian")
    p.set_defaults(func=cmd_interpret)

    p = sub.add_parser("quadrant", help="four-quadrant table for a cohort CSV")
    p.add_argument("--csv", required=True)
    p.add_argument("--panel", default="cpic4", choices=sorted(PANELS))
    p.set_defaults(func=cmd_quadrant)

    args = ap.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
