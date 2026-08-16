"""The only call site into ``anukriti-pgx-core``. Never overrides it.

``asl`` does not re-derive a phenotype, does not adjust an activity score, and
does not reclassify an allele. It calls the engine, treats the answer as final,
and attaches provenance and residual risk around it. Everything of value in this
package is in what surrounds the call, not in the call.

Two engine behaviours are surfaced rather than absorbed:

* **Over-specification.** If more than two variant-allele copies are detected the
  engine truncates the diplotype to two (``pgx-core`` 0.7.3 flags this via
  ``overspecified``). ``asl`` refuses rather than papering over it, because a
  diplotype that omits detected alleles understates the genotype.
* **The activity-score sentinel.** DPYD resolves through a named-diplotype table
  rather than an activity-score sum, so the engine returns ``-1.0``. That value
  must never reach a report; ``report`` prints CPIC's per-allele activity values
  instead.
"""

from __future__ import annotations

from anukriti_pgx_core import __version__ as ENGINE_VERSION
from anukriti_pgx_core.calling import DPYDCaller
from anukriti_pgx_core.phenotype.clinical_action import action_for
from anukriti_pgx_core.types import VCFVariant

from .models import AssayResult, TestStatus, Zygosity

# Plus-strand genomic REF/ALT per rsID, from CPIC
# sequence_location.chromosomelocation (g.<pos><REF>><ALT>, GRCh37).
#
# This is the same convention pgx-core 0.7.2+ uses in its allele table and the
# same convention a VCF and a dbSNP-forward assay report use, so the letters
# pass through unchanged. There is deliberately no translation layer here — the
# anukriti-api flip table that existed to compensate for the pre-0.7.2 engine
# inverted homozygous calls once the engine was fixed, and its removal is pinned
# by a test in that repo.
PLUS_STRAND: dict[str, tuple[str, str]] = {
    "rs3918290": ("C", "T"),   # *2A        g.97450058C>T
    "rs55886062": ("A", "C"),  # *13        g.97515787A>C
    "rs67376798": ("T", "A"),  # c.2846A>T  g.97082391T>A
    "rs56038477": ("C", "T"),  # HapB3 tag  g.97573863C>T
    "rs1801265": ("A", "G"),   # *9A        g.97883329A>G
    "rs1801160": ("C", "T"),   # *6         g.97305364C>T
    "rs1801159": ("T", "C"),   # *5         g.97515839T>C
    "rs2297595": ("T", "C"),   # M166V      g.97699535T>C
    "rs1801158": ("C", "T"),   # *4         g.97515865C>T
    "rs1801266": ("G", "A"),   # *8         g.97691776G>A
    "rs1801268": ("C", "A"),   # *10        g.97079071C>A
    "rs78060119": ("C", "A"),  # *12        g.97573943C>A
    "rs55971861": ("T", "G"),  # c.1906A>C  g.97382461T>G
}

# rsID -> CPIC allele name, for evidence lookup.
CPIC_ALLELE: dict[str, str] = {
    "rs3918290": "c.1905+1G>A (*2A)",
    "rs55886062": "c.1679T>G (*13)",
    "rs67376798": "c.2846A>T",
    "rs56038477": "c.1129-5923C>G, c.1236G>A (HapB3)",
    "rs1801265": "c.85T>C (*9A)",
    "rs1801160": "c.2194G>A (*6)",
    "rs1801159": "c.1627A>G (*5)",
    "rs2297595": "c.496A>G",
    "rs1801158": "c.1601G>A (*4)",
    "rs1801266": "c.703C>T (*8)",
    "rs1801268": "c.2983G>T (*10)",
    "rs78060119": "c.1156G>T (*12)",
    "rs55971861": "c.1906A>C",
}

# The four variants CPIC's guideline makes actionable. A guideline-conformant
# panel tests these; in CPIC's Central/South Asian reference data they sum to
# 0.02538 allele frequency, and *13 is 0.00000.
CPIC_ACTIONABLE = (
    "rs3918290",
    "rs55886062",
    "rs67376798",
    "rs56038477",
)

_GT = {Zygosity.HOM_REF: "0/0", Zygosity.HET: "0/1", Zygosity.HOM_VAR: "1/1"}


class InterpretationRefused(RuntimeError):
    """Raised when the engine's answer cannot be reported as-is."""


def call_engine(assay: AssayResult) -> tuple[str, str, str, bool, int, dict[str, int]]:
    """Run the deterministic engine.

    Returns ``(diplotype, phenotype, clinical_action, overspecified,
    copies_detected, allele_counts)``.

    ``allele_counts`` is the engine's per-allele copy count. It is returned
    alongside the diplotype so a report can show that the two named alleles
    account for everything detected, rather than asking the reader to trust it
    (``pgx-core`` 0.7.3, §8 item 1).

    Raises ``InterpretationRefused`` when more variant copies were detected than
    a diploid locus admits, because the engine's diplotype then omits alleles
    that were found and reporting it would understate the genotype.
    """
    variants: dict[str, VCFVariant] = {}
    for call in assay.tested:
        pair = PLUS_STRAND.get(call.rsid)
        if pair is None:
            raise InterpretationRefused(
                f"{call.rsid}: no pinned plus-strand REF/ALT. asl will not guess "
                f"a strand convention — add the CPIC chromosomelocation pair."
            )
        ref, alt = pair
        variants[call.rsid] = VCFVariant(
            ref=ref, alt=alt, genotype=_GT[call.zygosity]  # type: ignore[index]
        )

    result = DPYDCaller().call(variants)
    phenotype = (
        result.phenotype.phenotype
        if hasattr(result.phenotype, "phenotype")
        else str(result.phenotype)
    )
    action = action_for("DPYD", phenotype, "capecitabine")

    overspecified = bool(getattr(result, "overspecified", False))
    copies = int(getattr(result, "copies_detected", assay.total_variant_copies))
    counts: dict[str, int] = dict(getattr(result, "allele_counts", {}) or {})

    if overspecified:
        detected = ", ".join(sorted(result.alleles_detected))
        raise InterpretationRefused(
            f"{copies} variant-allele copies were detected at a diploid locus "
            f"({detected}), so the engine truncated the diplotype to "
            f"{result.diplotype!r} and it no longer represents everything found. "
            f"asl refuses to report a truncated diplotype. Resolve the input "
            f"(phasing, sample identity, or assay artifact) before interpreting."
        )

    return result.diplotype, phenotype, action, overspecified, copies, counts


def engine_version() -> str:
    return ENGINE_VERSION
