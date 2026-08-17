"""Types encoding the boundary between what was tested and what was concluded.

The single most important distinction in this package is between:

    NOT_TESTED  — the assay never looked for this variant
    ABSENT      — the assay looked and did not find it

Conflating those two is the defect this whole package exists to prevent. A
4-variant PCR panel that never assayed ``*6`` has established nothing about
``*6``, and a report that says "Normal Metabolizer" without saying so has
converted an absence of looking into an absence of risk.

Extraction and adjudication are split as *types*, following
``cohortfit/models.py``: an ``AssayResult`` records what a laboratory observed,
and an ``Interpretation`` records what can be concluded from it. Nothing in the
first may silently become the second.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import UTC, datetime
from enum import Enum


class Zygosity(str, Enum):
    """Observed zygosity at a tested locus."""

    HOM_REF = "hom_ref"
    HET = "het"
    HOM_VAR = "hom_var"


class TestStatus(str, Enum):
    """Whether a variant was assayed at all.

    ``NOT_TESTED`` is not a kind of negative result. It is the absence of a
    result, and it must survive all the way to the report.
    """

    TESTED = "tested"
    NOT_TESTED = "not_tested"


class Basis(str, Enum):
    """Assay basis behind an allele-function call, classified from CPIC."""

    IN_VITRO = "in_vitro"
    EX_VIVO = "ex_vivo"
    IN_VIVO = "in_vivo"
    UNCLASSIFIED = "unclassified"


class Regimen(str, Enum):
    """Fluoropyrimidine exposure setting.

    This exists because the published grade >=3 toxicity rate among patients
    without an actionable DPYD variant differs by setting — 22.7% on systemic
    chemotherapy (Henricks 2018) against 13.6% where chemoradiation is included
    at lower fluoropyrimidine doses (Lunenburg 2018). Quoting one where the
    other applies misleads in both directions, so ``residual`` refuses to quote
    a rate when this is ``UNKNOWN``.

    ``INDUCTION_TPF`` is separate from ``SYSTEMIC`` on purpose, and it is the
    setting that matters most in India. Induction docetaxel/cisplatin/5-FU for
    locally advanced head and neck cancer delivers 5-FU at 750 mg/m2/day as a
    continuous infusion on days 1-5 alongside two other cytotoxics — and head
    and neck cancer is one of India's largest cancer groups (oral cavity and
    pharynx alone: 198,438 incident cases in 2022, ICMR-NCDIR, PMID 36510887).
    Observed grade >=3 toxicity in Indian routine practice on this regimen is
    roughly threefold the Henricks figure: neutropenia 56.9%, febrile
    neutropenia 20.7%, mucositis 67.2%, diarrhoea 74.1% (Patil 2016,
    J Community Support Oncol 14(10):412-419, n=58 — "the toxicities differ
    substantially from the toxicities seen in trial settings").

    Before this member existed, TPF had no correct bucket and mapped to
    ``SYSTEMIC``, so an Indian head-and-neck patient was shown a 22.7% rate
    measured in a Dutch, largely colorectal, mostly capecitabine cohort. There
    is no published wild-type-stratified grade >=3 rate for this regimen, so
    ``residual`` **refuses to quote one** rather than substituting the nearest
    available number. The refusal is the correct output: the honest answer is
    that the comparator does not exist yet, and inventing a proxy for it is how
    a Dutch rate ends up in an Indian patient's report.
    """

    SYSTEMIC = "systemic"
    CHEMORADIATION = "chemoradiation"
    INDUCTION_TPF = "induction_tpf"
    UNKNOWN = "unknown"


class Drug(str, Enum):
    """The fluoropyrimidine the CPIC recommendation is being read for.

    CPIC publishes DPYD recommendations per drug, and ``asl`` displays CPIC's
    recommendation, so the report has to say which drug it is displaying one
    for. Without this the report showed a bare ``AVOID`` with no drug attached —
    correct for a fluoropyrimidine and meaningless, or worse, over-general, for
    anything else a reader might have in mind.

    The enum is closed to the two agents CPIC's DPYD guideline covers.
    ``asl.interpret`` refuses a drug for which the engine returns no action
    rather than printing a blank one.
    """

    CAPECITABINE = "capecitabine"
    FLUOROURACIL = "fluorouracil"


@dataclass(frozen=True)
class VariantCall:
    """One locus on one assay: what was looked for, and what was seen."""

    rsid: str
    status: TestStatus
    zygosity: Zygosity | None = None
    allele_label: str | None = None  # CPIC allele name, e.g. "c.2194G>A (*6)"

    def __post_init__(self) -> None:
        if self.status is TestStatus.TESTED and self.zygosity is None:
            raise ValueError(
                f"{self.rsid}: a TESTED locus must record a zygosity. "
                f"If the assay did not look, use TestStatus.NOT_TESTED — "
                f"an untested locus is not a negative result."
            )
        if self.status is TestStatus.NOT_TESTED and self.zygosity is not None:
            raise ValueError(
                f"{self.rsid}: a NOT_TESTED locus cannot have a zygosity."
            )

    @property
    def copies(self) -> int:
        """Variant-allele copies observed. Raises for an untested locus."""
        if self.status is TestStatus.NOT_TESTED:
            raise ValueError(
                f"{self.rsid}: copies is undefined for an untested locus"
            )
        return {
            Zygosity.HOM_REF: 0,
            Zygosity.HET: 1,
            Zygosity.HOM_VAR: 2,
        }[self.zygosity]  # type: ignore[index]


@dataclass(frozen=True)
class AssayResult:
    """What a laboratory actually did and observed. No conclusions.

    ``panel_name`` and the full ``calls`` list are both required: an assay whose
    variant content is unstated cannot be interpreted, and ``assay.parse``
    refuses it rather than assuming a panel.
    """

    gene: str
    panel_name: str
    calls: tuple[VariantCall, ...]
    regimen: Regimen = Regimen.UNKNOWN
    population: str | None = None  # CPIC population key, e.g. "Central/South Asian"
    drug: Drug = Drug.CAPECITABINE

    def __post_init__(self) -> None:
        # A repeated locus means two contradictory or redundant observations of
        # the same position. The engine is keyed by rsID, so a duplicate would
        # silently overwrite its predecessor and the losing genotype would
        # disappear without appearing anywhere in the report.
        seen: dict[str, VariantCall] = {}
        for call in self.calls:
            prior = seen.get(call.rsid)
            if prior is not None:
                raise ValueError(
                    f"{call.rsid} appears twice in this assay "
                    f"({prior.status.value}"
                    + (f"/{prior.zygosity.value}" if prior.zygosity else "")
                    + f" and {call.status.value}"
                    + (f"/{call.zygosity.value}" if call.zygosity else "")
                    + "). A locus has one observation. Resolve the duplicate "
                    "before interpreting — asl will not choose between two "
                    "readings of the same position."
                )
            seen[call.rsid] = call

    @property
    def tested(self) -> tuple[VariantCall, ...]:
        return tuple(c for c in self.calls if c.status is TestStatus.TESTED)

    @property
    def not_tested(self) -> tuple[VariantCall, ...]:
        return tuple(c for c in self.calls if c.status is TestStatus.NOT_TESTED)

    @property
    def found(self) -> tuple[VariantCall, ...]:
        """Tested loci carrying at least one variant copy."""
        return tuple(c for c in self.tested if c.copies > 0)

    @property
    def total_variant_copies(self) -> int:
        return sum(c.copies for c in self.tested)


@dataclass(frozen=True)
class EvidenceRecord:
    """The provenance behind one allele-function call, read from CPIC."""

    allele_label: str
    clinical_function: str
    activity_value: str | None
    evidence_strength: str | None
    findings_raw: str | None
    basis: Basis
    citations: tuple[str, ...]
    frequency_in_population: float | None
    population: str | None
    biochemical_function: str | None = None


@dataclass(frozen=True)
class Flag:
    """A named uncertainty. Never a phenotype change.

    ``rule_id`` is stable and greppable so a clinician or auditor can trace any
    line of a report back to the record that produced it.
    """

    rule_id: str
    summary: str
    detail: str
    citations: tuple[str, ...] = ()
    phase_unknown: bool = False


@dataclass(frozen=True)
class ResidualRisk:
    """What a result does *not* exclude.

    ``quoted_rate`` is ``None`` in two distinct cases, and both are refusals
    rather than missing values:

    * the regimen was not supplied, so no setting-matched rate applies;
    * **a variant was found**, so the wild-type rate does not describe this
      patient at all. Henricks 2018 and Lunenburg 2018 both report the grade >=3
      rate among patients *without* an actionable variant. Printing it to a
      carrier attributes a comparison group's risk to someone who is not in it,
      and for a Poor Metabolizer it understates the risk in the reassuring
      direction — the exact failure this package exists to prevent, applied to
      the patient with the most to lose.

    ``applies_to_this_patient`` records which case holds, so ``report`` can say
    why a rate is absent instead of leaving a reader to infer it.
    """

    quoted_rate: str | None
    rate_source: str | None
    regimen: Regimen
    panel_sensitivity: str
    cpic_caveat: str
    cpic_citation: str
    refusal_reason: str | None = None
    applies_to_this_patient: bool = True


@dataclass(frozen=True)
class Interpretation:
    """The full deterministic record. Everything a report may say.

    ``report`` and ``narrate`` may only render what is in here. If a claim is
    not present in this object it may not appear in the output.
    """

    assay: AssayResult
    diplotype: str
    phenotype: str
    cpic_recommendation_action: str
    evidence: tuple[EvidenceRecord, ...]
    flags: tuple[Flag, ...]
    residual: ResidualRisk
    uninformative_loci: tuple[tuple[str, str], ...] = ()  # (rsid, reason)
    overspecified: bool = False
    copies_detected: int = 0
    allele_counts: tuple[tuple[str, int], ...] = ()  # (allele, copies) from engine
    engine_version: str = ""
    generated_at: datetime = field(
        default_factory=lambda: datetime.now(UTC)
    )

    @property
    def is_panel_negative(self) -> bool:
        """True when no actionable variant was found among tested loci."""
        return len(self.assay.found) == 0
