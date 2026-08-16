"""Regenerate the per-allele evidence table verbatim from CPIC's live API.

Why this script exists
----------------------
An earlier design proposed hand-authoring an evidence table carrying the basis,
strength and citations behind each allele-function call. That was wrong: CPIC's
API already publishes almost all of it. Hand-authoring it would repeat the
failure class that ``anukriti-pgx-core`` 0.7.2 was released to fix — a
defensible-looking value written into a data file with nothing recording where
it came from.

So every field below is **copied from CPIC** and never edited:

    clinical_function      <- allele.clinicalfunctionalstatus
    biochemical_function   <- allele.functionalstatus
    activity_value         <- allele.activityvalue
    evidence_strength      <- allele.strength      (Definitive/Strong/Moderate/Limited)
    findings_raw           <- allele.findings      (encodes assay type + substrate)
    citations              <- allele.citations     (PMIDs)
    frequency              <- allele.frequency     (per CPIC population)

One field is *derived*, with an explicit fallback:

    basis  <- classified from findings_raw into in_vitro / ex_vivo / in_vivo

``findings_raw`` is free text, so classification can fail. When it does the
value is ``unclassified`` — never ``in_vitro``. Guessing the assay type would be
the same error in miniature.

Three fields cannot be generated because CPIC does not publish them, and they
live in ``src/asl/data/<GENE>_conflicts.json`` instead: ``cohort_n``,
``cohort_ancestry`` and ``conflict``. They are display metadata only and never
feed a phenotype.

Usage (from the asl repo root)::

    python scripts/regen_evidence_tables.py            # write
    python scripts/regen_evidence_tables.py --check    # exit 1 if drifted

``--check`` is what CI runs: regenerate in memory, diff against the committed
file, so divergence from CPIC fails a build instead of surviving to an audit.

Network: this script is the only code in the package that touches the network.
The library itself reads the pinned artifact and stays pure, which is what makes
a hospital demo work without connectivity.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from pathlib import Path

_API = "https://api.cpicpgx.org/v1"

# Genes asl interprets. DPYD first because that is the shipped workflow;
# NUDT15/TPMT are the planned second workflow (CPIC thiopurine 2025 update,
# PMID 41618934) and are generated now so the table is ready.
_GENES = ("DPYD", "NUDT15", "TPMT")

_OUT_DIR = Path("src/asl/data")


def _get(path: str) -> list[dict]:
    req = urllib.request.Request(
        f"{_API}/{path}", headers={"Accept": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=120) as fh:  # noqa: S310 (pinned host)
        return json.load(fh)


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from asl.evidence import classify_basis  # noqa: E402


def build_gene_table(gene: str) -> dict:
    """Build the evidence table for one gene, entirely from CPIC."""
    alleles = _get(f"allele?genesymbol=eq.{gene}")
    if not alleles:
        raise SystemExit(f"refusing to write an empty table: CPIC returned no {gene} alleles")

    rows: dict[str, dict] = {}
    for a in alleles:
        name = (a.get("name") or "").strip()
        if not name:
            continue
        findings_raw = a.get("findings")
        rows[name] = {
            "clinical_function": a.get("clinicalfunctionalstatus"),
            "biochemical_function": a.get("functionalstatus"),
            "activity_value": a.get("activityvalue"),
            "evidence_strength": a.get("strength"),
            "findings_raw": findings_raw,
            "basis": classify_basis(findings_raw),
            "citations": sorted(a.get("citations") or []),
            "frequency": a.get("frequency") or {},
        }

    return {
        "_gene": gene,
        "_source": (
            f"CPIC API allele table (api.cpicpgx.org/v1/allele?genesymbol=eq.{gene})"
        ),
        "_note": (
            "Every field except `basis` is copied verbatim from CPIC and must "
            "never be edited. `basis` is classified from `findings_raw` and "
            "falls back to 'unclassified' rather than guessing. Regenerate with "
            "scripts/regen_evidence_tables.py; CI runs it with --check so drift "
            "from CPIC fails a build."
        ),
        "_allele_count": len(rows),
        "alleles": dict(sorted(rows.items())),
    }


def render(table: dict) -> str:
    return json.dumps(table, indent=2, sort_keys=False) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--check",
        action="store_true",
        help="Do not write; exit 1 if a committed table differs from CPIC.",
    )
    ap.add_argument(
        "--gene",
        action="append",
        help="Limit to one gene (repeatable). Defaults to all.",
    )
    args = ap.parse_args(argv)

    genes = tuple(args.gene) if args.gene else _GENES
    _OUT_DIR.mkdir(parents=True, exist_ok=True)

    drifted: list[str] = []
    for gene in genes:
        table = build_gene_table(gene)
        text = render(table)
        path = _OUT_DIR / f"{gene}_evidence.json"

        if args.check:
            if not path.exists():
                drifted.append(f"{path} (missing)")
            elif path.read_text() != text:
                drifted.append(str(path))
        else:
            path.write_text(text)
            print(f"wrote {path} ({table['_allele_count']} alleles)")

    if args.check:
        if drifted:
            print(
                "evidence tables have drifted from CPIC's live API:\n  - "
                + "\n  - ".join(drifted)
                + "\n\nRun: python scripts/regen_evidence_tables.py",
                file=sys.stderr,
            )
            return 1
        print(f"evidence tables match CPIC ({len(genes)} gene(s)).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
