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
from anukriti_pgx_core.phenotype.clinical_action import action_for, details_for_action
from anukriti_pgx_core.types import Diplotype, VCFVariant

from .models import AssayResult, Zygosity

# Plus-strand genomic REF/ALT per rsID, from CPIC
# sequence_location.chromosomelocation (g.<pos><REF>><ALT>). The positions are
# GRCh38, anchored by rs3918290 being chr1:97450058 on GRCh38 against 97915614
# on GRCh37. 13 of the 32 were independently corroborated by position-and-alleles
# matching against the 1000 Genomes GRCh38 high-coverage panel; the other 19 are
# absent from that panel (too rare for 3,202 samples) and so carry CPIC as their
# only source.
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
    "rs115232898": ("T", "C"),  # c.557A>G  g.97699474T>C
    "rs146356975": ("T", "C"),  # c.868A>G  g.97595149T>C
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
    # The remaining loci the engine calls. Added 2026-10-13: pgx-core 0.11.0
    # widened its allele table to 33 alleles across 32 sites, and until this
    # block existed `asl` pinned only the 19 above — so the other 13 could not
    # be named in a panel preset and were reachable only if a VCF happened to
    # supply them. Part of the engine's coverage was therefore not orderable
    # through the product.
    #
    # All thirteen are No function except c.1314T>G (Decreased function). Each
    # REF/ALT pair below was read from CPIC `sequence_location.chromosomelocation`
    # and re-verified against the live api.cpicpgx.org on 2026-10-13; the same
    # derivation reproduces all 19 pairs above unchanged, which is what makes
    # the new ones trustworthy rather than merely plausible.
    "rs72549310": ("G", "A"),   # c.61C>T    g.97883353G>A
    "rs72549308": ("T", "G"),   # c.601A>C   g.97699430T>G
    "rs72549307": ("T", "C"),   # c.632A>G   g.97699399T>C
    "rs183385770": ("C", "T"),  # c.1024G>A  g.97593322C>T
    "rs143154602": ("G", "A"),  # c.1057C>T  g.97593289G>A
    "rs186169810": ("A", "C"),  # c.1314T>G  g.97573785A>C
    "rs111858276": ("T", "C"),  # c.1484A>G  g.97549600T>C
    "rs59086055": ("G", "A"),   # c.1774C>T  g.97450190G>A
    "rs138616379": ("C", "T"),  # c.1775G>A  g.97450189C>T
    "rs145773863": ("C", "T"),  # c.1777G>A  g.97450187C>T
    "rs137999090": ("C", "T"),  # c.2021G>A  g.97373598C>T
    "rs141044036": ("T", "C"),  # c.2872A>G  g.97082365T>C
    "rs72547601": ("T", "C"),   # c.2933A>G  g.97079121T>C
}

#: GRCh38 chr1 position of each locus, as data rather than a comment.
#:
#: Consumers match 1000 Genomes records on position plus REF/ALT. The GRCh37
#: phase 3 release has an empty ID column, so matching on rsID there finds
#: nothing at all and does so silently.
GRCH38_POS: dict[str, int] = {
    "rs3918290": 97450058,
    "rs55886062": 97515787,
    "rs67376798": 97082391,
    "rs75017182": 97579893,
    "rs56038477": 97573863,
    "rs112766203": 97305279,
    "rs115232898": 97699474,
    "rs146356975": 97595149,
    "rs55674432": 97098616,
    "rs72549304": 97549609,
    "rs1801265": 97883329,
    "rs1801160": 97305364,
    "rs1801159": 97515839,
    "rs2297595": 97699535,
    "rs1801158": 97515865,
    "rs1801266": 97691776,
    "rs1801268": 97079071,
    "rs78060119": 97573943,
    "rs55971861": 97382461,
    "rs72549310": 97883353,
    "rs72549308": 97699430,
    "rs72549307": 97699399,
    "rs183385770": 97593322,
    "rs143154602": 97593289,
    "rs186169810": 97573785,
    "rs111858276": 97549600,
    "rs59086055": 97450190,
    "rs138616379": 97450189,
    "rs145773863": 97450187,
    "rs137999090": 97373598,
    "rs141044036": 97082365,
    "rs72547601": 97079121,
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
    "rs115232898": "c.557A>G",
    "rs146356975": "c.868A>G",
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
    # Same thirteen as the tail of PLUS_STRAND, in the same order. Every name
    # below is a key in the pinned CPIC evidence table, so each carries
    # citations and a function assignment; a name that was not would make
    # `evidence` refuse rather than report an untraceable claim.
    #
    # This dict's iteration order is the report's coverage-strip order, so the
    # thirteen are appended rather than interleaved: a locus keeps the column
    # it had before the widening.
    "rs72549310": "c.61C>T",
    "rs72549308": "c.601A>C",
    "rs72549307": "c.632A>G",
    "rs183385770": "c.1024G>A",
    "rs143154602": "c.1057C>T",
    "rs186169810": "c.1314T>G",
    "rs111858276": "c.1484A>G",
    "rs59086055": "c.1774C>T",
    "rs138616379": "c.1775G>A",
    "rs145773863": "c.1777G>A",
    "rs137999090": "c.2021G>A",
    "rs141044036": "c.2872A>G",
    "rs72547601": "c.2933A>G",
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


def unknown_loci_from(result: object) -> tuple[str, ...] | None:
    """The engine's unknown loci, or ``None`` when it cannot report them.

    ``Diplotype.unknown_loci`` landed in ``pgx-core`` after the 0.11.0 release,
    so on the pinned engine the attribute is absent. ``getattr(..., ())`` would
    yield an empty tuple, which a consumer reads as "the engine looked and found
    none" — the same silent, falsely reassuring collapse that
    ``require_hapb3_capable_engine`` exists to stop.

    So the three states are kept distinct and carried into the payload:
    ``None`` (this engine cannot report unknown loci), ``()`` (it reported
    none), and a non-empty tuple (it reported these). A hard refusal is not
    used here because every released engine would fail it; the gap is disclosed
    instead, and disappears on its own once the pin moves past 0.11.0.
    """
    if not hasattr(Diplotype, "unknown_loci"):
        return None
    return tuple(getattr(result, "unknown_loci", ()) or ())


def call_engine(
    assay: AssayResult,
) -> tuple[str, str, str, bool, int, dict[str, int], str, tuple[str, ...] | None]:
    """Run the deterministic engine.

    Returns ``(diplotype, phenotype, clinical_action, overspecified,
    copies_detected, allele_counts, hapb3_basis, unknown_loci)``.

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
        unknown_loci_from(result),
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


def cpic_recommendation(phenotype: str, drug: str) -> tuple[str, str, tuple[str, ...]]:
    """CPIC's own dosing sentence for this phenotype, with its activity range.

    Returns ``(recommendation, activity_score_range, citations)`` read verbatim
    from the engine's pinned CPIC clinical-action table.

    This exists so the report can print a dose recommendation **in CPIC's
    words**. Before it, the report's only dosing line was the closed-enum token
    ``REDUCE_50PCT`` — correct, traceable, and not a sentence anyone can act on
    or audit. A reader cannot check an enum against the guideline; they can
    check "Reduce starting dose by 50%, then titrate to toxicity or therapeutic
    drug monitoring" against PMID 29152729 in one minute.

    It is also the regulatory posture. CDSCO classifies software by the voice a
    claim is made in, not by where it sits on the page: a dose attributed to
    CPIC and quoted from CPIC is CPIC informing the clinician, where the same
    dose in asl's own words would be asl driving the decision. So the sentence
    is quoted, never paraphrased, and the attribution travels with it.

    Raises ``InterpretationRefused`` rather than returning an empty sentence,
    for the same reason ``call_engine`` refuses a blank action: a report with a
    dosing heading and nothing under it reads as "no dose change needed".
    """
    record = details_for_action("DPYD", phenotype, drug)
    text = str(record.get("recommendation") or "").strip()
    if not text:
        raise InterpretationRefused(
            f"the engine's CPIC clinical-action table carries no recommendation "
            f"sentence for phenotype {phenotype!r} with {drug}. asl will not "
            f"print a dosing section with nothing in it, nor substitute wording "
            f"of its own for CPIC's."
        )
    return (
        text,
        str(record.get("activity_score_range") or ""),
        tuple(str(c) for c in record.get("citations") or ()),
    )


def engine_version() -> str:
    return ENGINE_VERSION
