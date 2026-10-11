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
from anukriti_pgx_core.types import Diplotype, VCFVariant

from .models import AssayResult, Zygosity

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
    "rs75017182": ("G", "C"),  # HapB3 CAUSAL c.1129-5923C>G  g.97579893G>C
    "rs56038477": ("C", "T"),  # HapB3 tag  c.1236G>A         g.97573863C>T
    "rs112766203": ("G", "A"),  # c.2279C>T g.97305279G>A
    "rs55674432": ("C", "A"),  # c.2639G>T  g.97098616C>A
    "rs72549304": ("G", "A"),  # c.1475C>T  g.97549609G>A
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
#
# The two HapB3 loci map to DIFFERENT CPIC allele names, which is the whole
# point: rs75017182 is the causal variant CPIC publishes as an allele in its own
# right, and rs56038477 is the benign synonymous tag SNP whose decreased-function
# reading is CPIC's inference from the haplotype. See ``hapb3_disclosure``.
CPIC_ALLELE: dict[str, str] = {
    "rs3918290": "c.1905+1G>A (*2A)",
    "rs55886062": "c.1679T>G (*13)",
    "rs67376798": "c.2846A>T",
    "rs75017182": "c.1129-5923C>G",
    "rs56038477": "c.1129-5923C>G, c.1236G>A (HapB3)",
    "rs112766203": "c.2279C>T",
    "rs55674432": "c.2639G>T",
    "rs72549304": "c.1475C>T",
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

# rsID -> the name of the variant the assay ACTUALLY interrogates at that locus.
#
# For every locus except the HapB3 tag this is the same as ``CPIC_ALLELE``. The
# tag SNP is the exception, and it matters: CPIC's allele name for it is the
# compound ``"c.1129-5923C>G, c.1236G>A (HapB3)"``, which *names the causal
# variant*. Displaying that as what was tested, in a report that simultaneously
# lists ``c.1129-5923C>G`` under "absence is not established", contradicts
# itself — the reader is told the same variant was both found and not looked at.
#
# So evidence lookup keeps CPIC's compound name (it is the row carrying the
# function and the frequency), while *display of what was assayed* names the
# locus that was actually read. This is the same distinction the disclosure
# makes, applied to the assay inventory.
ASSAYED_VARIANT: dict[str, str] = dict(CPIC_ALLELE)
ASSAYED_VARIANT["rs56038477"] = "c.1236G>A (HapB3 exonic tag SNP)"


def assayed_variant_name(rsid: str, fallback: str | None = None) -> str:
    """What the assay read at ``rsid``, for display in the assay inventory."""
    return ASSAYED_VARIANT.get(rsid) or fallback or rsid

_GT = {Zygosity.HOM_REF: "0/0", Zygosity.HET: "0/1", Zygosity.HOM_VAR: "1/1"}


class InterpretationRefused(RuntimeError):
    """Raised when the engine's answer cannot be reported as-is."""


def require_hapb3_capable_engine() -> None:
    """Refuse to interpret on an engine that cannot report the HapB3 basis.

    ``hapb3_basis`` arrived in ``pgx-core`` 0.9.0. On an older engine the
    attribute is simply absent, and ``getattr(..., "")`` would quietly yield
    ``""`` — which this package reads as "no HapB3-family allele was called".

    That failure would be **silent and falsely reassuring in exactly the way
    this package exists to prevent**: a tag-only call would be reported as an
    Intermediate Metabolizer with CPIC's 50% dose reduction and *no disclosure*,
    because the report would believe there was nothing to disclose. The
    difference between "no HapB3 allele" and "a HapB3 allele inferred from a
    benign tag SNP" would vanish into the same empty string.

    An older engine also cannot call ``c.2279C>T``, ``c.2639G>T`` or
    ``c.1475C>T`` — 21.3% of actionable Central/South Asian DPYD allele
    frequency — and would return ``Indeterminate`` for them.

    So the version floor is checked once, loudly, rather than degrading.
    """
    if not hasattr(Diplotype, "hapb3_basis"):
        raise InterpretationRefused(
            f"anukriti-pgx-core {ENGINE_VERSION} does not report "
            f"Diplotype.hapb3_basis, which arrived in 0.9.0. Without it asl "
            f"cannot tell a tag-SNP-inferred HapB3 call from a directly "
            f"observed one, and CPIC requires that distinction to be stated in "
            f"the report. Refusing rather than emitting a report with the "
            f"disclosure silently missing. Upgrade to >=0.9.0."
        )


def call_engine(
    assay: AssayResult,
) -> tuple[str, str, str, bool, int, dict[str, int], str]:
    """Run the deterministic engine.

    Returns ``(diplotype, phenotype, clinical_action, overspecified,
    copies_detected, allele_counts, hapb3_basis)``.

    ``allele_counts`` is the engine's per-allele copy count. It is returned
    alongside the diplotype so a report can show that the two named alleles
    account for everything detected, rather than asking the reader to trust it
    (``pgx-core`` 0.7.3, §8 item 1).

    ``hapb3_basis`` is which locus produced a HapB3-family decreased-function
    call (``pgx-core`` 0.9.0). It is carried out of the engine unchanged so
    ``hapb3_disclosure`` can emit CPIC's mandated tag-SNP statement; ``asl`` does
    not compute it, because which locus was observed is the engine's knowledge,
    not this package's.

    Raises ``InterpretationRefused`` when more variant copies were detected than
    a diploid locus admits, because the engine's diplotype then omits alleles
    that were found and reporting it would understate the genotype.
    """
    require_hapb3_capable_engine()
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
    action = action_for("DPYD", phenotype, assay.drug.value)
    if not action or action == "INDETERMINATE":
        # A blank action reaching a report would render as an empty "CPIC action"
        # line, which reads as "nothing to do" rather than "we could not
        # determine what CPIC advises". The engine returns '' for an unknown
        # phenotype/drug pair by design (it never raises), so the check belongs
        # here.
        raise InterpretationRefused(
            f"the engine returned no CPIC clinical action for phenotype "
            f"{phenotype!r} with {assay.drug.value}. asl will not print a blank "
            f"or INDETERMINATE action, because an empty recommendation line "
            f"reads as 'nothing to do'. Resolve the phenotype or the drug pair."
        )

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

    return (
        result.diplotype,
        phenotype,
        action,
        overspecified,
        copies,
        counts,
        str(getattr(result, "hapb3_basis", "") or ""),
    )


# CPIC's own wording, quoted rather than paraphrased. The guideline requires this
# to appear in a test report whenever "decreased function" rests on the tag SNP.
_CPIC_TAG_DISCLOSURE = (
    '"decreased function" was inferred by detecting the exonic tag SNP '
    "c.1236G>A (rs56038477), and in rare cases the causal decreased function "
    "variant c.1129-5923C>G (rs75017182) may not be present despite having "
    "this tag SNP"
)


def hapb3_disclosure(basis: str) -> tuple[str, str] | None:
    """CPIC's mandated tag-SNP disclosure, as ``(summary, detail)``.

    Returns ``None`` when no disclosure is required — that is, when no
    HapB3-family allele was called, or when the call rests on the causal variant
    itself and there is nothing to qualify.

    Why this exists
    ---------------

    ``HapB3`` is CPIC's most frequent actionable DPYD allele in this population
    (0.019658 in CPIC's Central/South Asian frequencies, 60.5% of all actionable
    DPYD allele frequency). It is defined by two variants: the causal intronic
    ``c.1129-5923C>G`` (rs75017182), which creates a cryptic splice site, and the
    **benign synonymous** ``c.1236G>A`` (rs56038477), which most assays genotype
    instead because it is cheaper to type.

    They are not in perfect linkage disequilibrium — 14 of 6265 tag carriers in
    *All of Us* WGS lacked the causal variant (PMID 38129972) — so a tag-only
    result can be a false positive, and a false positive here means a 50%
    fluoropyrimidine dose reduction in a patient who did not need one, during
    curative-intent treatment.

    CPIC therefore requires the basis to be stated in the report. ``asl`` emits
    CPIC's sentence and does **not** alter the recommendation: the phenotype and
    the clinical action are CPIC's, unchanged, in every basis. What changes is
    that the reader is told what the call rests on.
    """
    if basis in ("", "causal", "haplotype"):
        return None

    if basis == "tag_without_causal":
        return (
            ("This decreased-function call rests on the tag SNP, and the causal "
            "variant was tested and NOT found."
            ),
            "    CPIC requires this to be disclosed: " + _CPIC_TAG_DISCLOSURE + ".\n"
            "    In this patient the causal variant c.1129-5923C>G (rs75017182) "
            "WAS tested\n"
            "    and was absent, while the tag SNP c.1236G>A (rs56038477) was "
            "found. That is\n"
            "    the documented linkage-disequilibrium breakdown (Turner/Haidar "
            "2024,\n"
            "    PMID 38129972: 14 of 6265 tag carriers in All of Us WGS lacked "
            "the causal\n"
            "    variant; replicated 1 of 46 in Int J Mol Sci 2025;26:8136), so "
            "this is the\n"
            "    genotype most likely to be a FALSE POSITIVE decreased-function "
            "call.\n"
            "    The recommendation above is CPIC's and is unchanged — asl does "
            "not\n"
            "    reclassify it. What this flag adds is that the inference the "
            "recommendation\n"
            "    rests on is contradicted at the causal locus in this patient, "
            "and that a\n"
            "    50% dose reduction in curative-intent treatment is the cost of "
            "acting on\n"
            "    a false positive. Confirmatory review of the causal locus is "
            "warranted\n"
            "    before reducing dose.",
        )

    # basis == "tag_only"
    return (
        ("This decreased-function call was inferred from the exonic tag SNP; the "
        "causal variant was not tested."
        ),
        "    CPIC requires this to be disclosed: " + _CPIC_TAG_DISCLOSURE + ".\n"
        "    The causal variant c.1129-5923C>G (rs75017182) was NOT among the "
        "loci this\n"
        "    assay tested, so its presence is inferred rather than observed. The "
        "two are\n"
        "    in high but imperfect linkage disequilibrium (LD 0.9985; 14 of 6265 "
        "tag\n"
        "    carriers in All of Us WGS lacked the causal variant — 0.223% —"
        " Turner/Haidar\n"
        "    2024, PMID 38129972; replicated 1 of 46 in Int J Mol Sci "
        "2025;26:8136).\n"
        "    The observed discordance was almost entirely in participants of "
        "European\n"
        "    ancestry and LD was perfect in the South Asian cohort "
        "(Suarez-Kurtz 2024,\n"
        "    PMID 38634417) — but no study has measured this linkage in an "
        "Indian cohort,\n"
        "    so the rate here is unmeasured rather than known to be zero.\n"
        "    The recommendation above is CPIC's and is unchanged. Testing the "
        "causal\n"
        "    variant directly resolves this; it is what CPIC now recommends and "
        "what\n"
        "    some commercial DPYD panels already genotype.",
    )


def engine_version() -> str:
    return ENGINE_VERSION
