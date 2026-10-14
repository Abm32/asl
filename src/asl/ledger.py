"""Append-only outcome capture. The part that compounds.

No genotype-to-toxicity registry is open to Indian sites — U-PGx/PREPARE is
EU-only, IGNITE and eMERGE are US-only, and India's PvPI has no genotype field
and handles individual case safety reports rather than cohorts. So the data that
could eventually correct an allele-function call for this population does not
exist, and cannot be joined from anywhere.

This module records it. Two design decisions carry the weight:

**``rsids_tested`` is mandatory.** It is the field published Indian studies omit,
and it is why Patil 2019's 1064 genotyped patients cannot answer the question
their genotypes were collected for: without knowing which loci were assayed, a
negative is uninterpretable. An entry without it is refused.

**No PHI, enforced rather than requested.** Field names are checked against a
denylist and free-text values are pattern-checked, because a site will paste a
spreadsheet column and the failure mode is silent.

Where the data goes, once there is enough of it: **publish first**, then submit a
CPIC inquiry. CPIC's reassessment route is real and its trigger is specific —
"an inquiry which cites **published literature** supportive of evidence against an
allele function assignment" (Tibben 2025, PMID 41175864), after which a staff
member circulates an interim evidence table to the whole Pharmacogene Curation
Expert Panel and a new assignment needs **70% consensus**. A JSONL file is not an
input to that process; a paper derived from it is. Then PharmVar for haplotype
registration, and ClinPGx, which accepts primary data.

There is a close precedent: in **October 2025, in response to a member inquiry**,
CPIC convened the SLCO1B1 panel to reassess "select alleles with **elevated
frequencies in underrepresented populations**" against emerging myotoxicity
evidence. That is this argument, for a different gene, already accepted as grounds
for convening a panel.

**One thing to be clear about, because it narrows the claim.** CPIC has already
*lowered* its evidence bar for DPYD specifically — the framework paper names DPYD
as its worked example of modifying "the threshold of evidence for clinical
actionability… in the setting of limited data", on explicit type-II-error grounds.
So the four alleles Indians carry are called *Normal function* against a threshold
deliberately tilted toward calling things actionable, with Strong evidence for
`*9A` and `*5`. CPIC has not neglected them. What CPIC has not had is **clinical
outcome data from this population**, which is the only thing left that could move
the call — and is exactly what this module collects.

That is also why changing our own table unilaterally is forbidden: it is precisely
what the `pgx-core` 0.7.2 release forbade.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from .residual import actionable_loci

# Field names that must never appear in a ledger entry. Checked on keys, so a
# site pasting a spreadsheet column called "patient_name" fails loudly.
_PHI_FIELDS = frozenset({
    "name", "patient_name", "first_name", "last_name", "surname",
    "mrn", "medical_record_number", "hospital_number", "uhid", "ip_number",
    "aadhaar", "aadhar", "pan", "phone", "mobile", "email", "address",
    "dob", "date_of_birth", "birth_date", "nok", "next_of_kin",
})

# Value patterns that look like identifiers regardless of the field name.
_PHI_PATTERNS = (
    (re.compile(r"\b\d{12}\b"), "a 12-digit number (Aadhaar-shaped)"),
    (re.compile(r"\b[6-9]\d{9}\b"), "a 10-digit mobile number"),
    (re.compile(r"[^@\s]+@[^@\s]+\.[a-z]{2,}", re.IGNORECASE), "an email address"),
    (re.compile(r"\b\d{4}-\d{2}-\d{2}\b"), "a full date (use onset_cycle instead)"),
)


class LedgerRefused(RuntimeError):
    """Raised when an entry cannot be recorded safely or usefully."""


@dataclass(frozen=True)
class ToxicityEvent:
    """One CTCAE-graded event. Onset is recorded as a cycle, never a date."""

    term: str
    grade: int
    onset_cycle: int | None = None

    def __post_init__(self) -> None:
        if not 0 <= self.grade <= 5:
            raise LedgerRefused(f"CTCAE grade must be 0-5, got {self.grade}")
        if self.onset_cycle is not None and self.onset_cycle < 1:
            raise LedgerRefused(
                f"onset_cycle must be 1 or greater, got {self.onset_cycle}. "
                f"Cycles are numbered from 1; 0 or a negative value is a data "
                f"entry error, and it would sort before every real event."
            )


# Recognised genotype tokens. The ledger stores what a site recorded, so it is
# deliberately permissive about *notation* — VCF, phased VCF, and word forms are
# all accepted — but not about tokens it cannot interpret at all.
WILD_TYPE_TOKENS = frozenset({
    "0", "0/0", "0|0", "hom_ref", "wt", "ref", "wild_type", "wild-type",
})
VARIANT_TOKENS = frozenset({
    "1", "2", "0/1", "1/0", "0|1", "1|0", "1/1", "1|1",
    "het", "hom_var", "hom", "var",
})
_READABLE_TOKENS = WILD_TYPE_TOKENS | VARIANT_TOKENS

# CPIC's closed set of DPYD metabolizer terms, plus the engine's refusal value.
# Free text here would make the column unanalysable, which defeats the point of
# recording it: the whole question the ledger exists to answer is whether the
# phenotype assigned at dosing predicted what happened.
CPIC_PHENOTYPES = frozenset({
    "Normal Metabolizer",
    "Intermediate Metabolizer",
    "Poor Metabolizer",
    "Indeterminate",
})

# Zygosity implied by a token, for the genotype/phenotype consistency check.
_COPIES = {
    **{t: 0 for t in WILD_TYPE_TOKENS},
    "0/1": 1, "1/0": 1, "0|1": 1, "1|0": 1, "het": 1, "1": 1,
    "1/1": 2, "1|1": 2, "hom_var": 2, "hom": 2, "2": 2, "var": 1,
}

# The four CPIC-actionable DPYD loci. A no-function homozygote at any of these
# cannot be a Normal Metabolizer, and that is what makes the check possible.
_ACTIONABLE = ("rs3918290", "rs55886062", "rs67376798", "rs56038477")
_NO_FUNCTION = ("rs3918290", "rs55886062")


@dataclass(frozen=True)
class LedgerEntry:
    """One patient's genotype and what happened to them.

    ``rsids_tested`` must list every locus the assay covered, including those
    that were negative. That is what makes the row usable.
    """

    anonymised_id: str
    cancer_type: str
    regimen: str
    drug: str
    assay_method: str
    rsids_tested: tuple[str, ...]
    genotypes: dict[str, str]
    cpic_phenotype_at_dosing: str
    engine_version: str
    toxicity: tuple[ToxicityEvent, ...] = ()
    starting_dose_mg_m2: float | None = None
    bsa_m2: float | None = None
    dose_modifications: tuple[str, ...] = ()
    cycles_received: int | None = None
    recorded_at: datetime = field(
        default_factory=lambda: datetime.now(UTC)
    )

    def __post_init__(self) -> None:
        # PHI first, before any other validation. A leak must be refused whether
        # or not the rest of the row is well-formed: the failure mode is a site
        # pasting a spreadsheet, and a row that is rejected for a bad genotype
        # token gives them no reason to remove the phone number before retrying.
        _reject_phi(self.as_dict())

        if not self.rsids_tested:
            raise LedgerRefused(
                "rsids_tested is empty. Record every locus the assay covered, "
                "including the negatives — without it a negative result is "
                "uninterpretable and the row cannot support any later analysis. "
                "This is the field published Indian cohorts omit."
            )
        unknown = set(self.genotypes) - set(self.rsids_tested)
        if unknown:
            raise LedgerRefused(
                f"genotypes reported for loci not listed as tested: "
                f"{sorted(unknown)}. A result cannot come from an untested locus."
            )
        if not str(self.anonymised_id).strip():
            raise LedgerRefused(
                "anonymised_id is empty. Every row needs a site-assigned "
                "identifier — without one, duplicates cannot be detected and a "
                "row cannot be corrected later by superseding it."
            )
        if not str(self.engine_version).strip():
            raise LedgerRefused(
                "engine_version is empty. The recorded phenotype is only "
                "meaningful against the engine and CPIC tables that produced it; "
                "a row that cannot be attributed to a version cannot be "
                "reanalysed when the tables change — and they do change (see "
                "P4_DPYD_HAPB3_CPIC_UPDATE_PENDING)."
            )

        # Genotype tokens must be interpretable. The ledger is permissive about
        # notation and strict about legibility: an uninterpretable token is
        # silently dropped from every later analysis, which is worse than a
        # refusal at entry because nobody finds out.
        for rsid, raw in self.genotypes.items():
            if str(raw).strip().lower() not in _READABLE_TOKENS:
                raise LedgerRefused(
                    f"{rsid}: genotype {raw!r} is not a recognised token. Use "
                    f"0/0, 0/1, 1/1 (phased forms accepted), or "
                    f"hom_ref/het/hom_var. An unreadable token is excluded from "
                    f"every later analysis, so recording it achieves nothing."
                )

        if self.cpic_phenotype_at_dosing not in CPIC_PHENOTYPES:
            raise LedgerRefused(
                f"cpic_phenotype_at_dosing {self.cpic_phenotype_at_dosing!r} is "
                f"not a CPIC metabolizer term. Use one of: "
                f"{', '.join(sorted(CPIC_PHENOTYPES))}. Free text makes the "
                f"column unanalysable, and whether the phenotype assigned at "
                f"dosing predicted the outcome is the entire question this "
                f"ledger exists to answer."
            )

        # A no-function homozygote cannot have been a Normal Metabolizer. This is
        # the one internal contradiction worth catching: it means the genotype and
        # the phenotype came from different patients, or the phenotype was typed
        # by hand. Either way the row is evidence of nothing and would pollute
        # exactly the analysis the ledger is for.
        for rsid in _NO_FUNCTION:
            token = str(self.genotypes.get(rsid, "")).strip().lower()
            if _COPIES.get(token) == 2 and (
                self.cpic_phenotype_at_dosing == "Normal Metabolizer"
            ):
                raise LedgerRefused(
                    f"{rsid} is recorded as homozygous variant (a no-function "
                    f"DPYD allele) while cpic_phenotype_at_dosing is "
                    f"'Normal Metabolizer'. Those cannot both be true — CPIC "
                    f"calls that diplotype a Poor Metabolizer. The genotype and "
                    f"the phenotype are from different sources; resolve which is "
                    f"wrong before recording the row."
                )

        for name, value, low in (
            ("cycles_received", self.cycles_received, 0),
            ("starting_dose_mg_m2", self.starting_dose_mg_m2, 0),
            ("bsa_m2", self.bsa_m2, 0),
        ):
            if value is not None and value <= low:
                raise LedgerRefused(
                    f"{name} must be greater than {low}, got {value}. A "
                    f"non-positive value is a data entry error and would enter "
                    f"any dose-intensity calculation as if it were real."
                )
        if self.bsa_m2 is not None and not 0.2 <= self.bsa_m2 <= 3.0:
            raise LedgerRefused(
                f"bsa_m2 {self.bsa_m2} is outside the physiologically plausible "
                f"range 0.2-3.0 m². Body surface area is a dose denominator, so "
                f"an implausible value silently corrupts dose intensity."
            )
        if self.cycles_received is not None:
            late = [
                t.onset_cycle
                for t in self.toxicity
                if t.onset_cycle is not None and t.onset_cycle > self.cycles_received
            ]
            if late:
                raise LedgerRefused(
                    f"toxicity recorded at cycle(s) {late} but only "
                    f"{self.cycles_received} cycle(s) were received. An event "
                    f"cannot occur in a cycle that never happened; this is "
                    f"usually a transcription slip between the two columns."
                )

    @property
    def max_grade(self) -> int:
        return max((t.grade for t in self.toxicity), default=0)

    def as_dict(self) -> dict:
        return {
            "anonymised_id": self.anonymised_id,
            "cancer_type": self.cancer_type,
            "regimen": self.regimen,
            "drug": self.drug,
            "assay_method": self.assay_method,
            "rsids_tested": list(self.rsids_tested),
            "genotypes": dict(self.genotypes),
            "cpic_phenotype_at_dosing": self.cpic_phenotype_at_dosing,
            "engine_version": self.engine_version,
            "toxicity": [
                {"term": t.term, "grade": t.grade, "onset_cycle": t.onset_cycle}
                for t in self.toxicity
            ],
            "starting_dose_mg_m2": self.starting_dose_mg_m2,
            "bsa_m2": self.bsa_m2,
            "dose_modifications": list(self.dose_modifications),
            "cycles_received": self.cycles_received,
            "max_toxicity_grade": self.max_grade,
            "recorded_at": self.recorded_at.isoformat(timespec="seconds"),
        }


def _reject_phi(payload: dict, _path: str = "") -> None:
    """Refuse anything that looks like patient-identifying information.

    Recurses into lists and nested mappings. A flat scan was not enough: the
    realistic leak is a free-text field that is itself a container —
    ``dose_modifications=["reduced 25% on 1978-04-12 per Dr Rao, ph 98..."]`` and
    a toxicity ``term`` — and those are exactly the fields a site fills in prose.
    A scan that only reads top-level strings passes them all.
    """
    for key, value in payload.items():
        where = f"{_path}{key}"
        if key.lower() in _PHI_FIELDS:
            raise LedgerRefused(
                f"field {where!r} is not accepted: the ledger holds no PHI. "
                f"Use an anonymised identifier assigned at the site."
            )
        _scan_value(where, value)


def _scan_value(where: str, value: object) -> None:
    """Pattern-check one value, recursing through containers."""
    if isinstance(value, str):
        for pattern, what in _PHI_PATTERNS:
            if pattern.search(value):
                raise LedgerRefused(
                    f"field {where!r} contains what looks like {what}. "
                    f"The ledger holds no PHI."
                )
    elif isinstance(value, dict):
        _reject_phi(value, f"{where}.")
    elif isinstance(value, (list, tuple)):
        for i, item in enumerate(value):
            _scan_value(f"{where}[{i}]", item)


def from_mapping(row: dict) -> LedgerEntry:
    """Build an entry from an arbitrary mapping, e.g. a CSV row.

    This is the realistic PHI entry point. ``LedgerEntry`` itself is a frozen
    dataclass, so a stray ``patient_name=`` keyword is a ``TypeError`` — but a
    site hands over a spreadsheet, and a spreadsheet has whatever columns it has.
    So the denylist is enforced here, before anything is constructed.
    """
    _reject_phi(row)

    rsids = row.get("rsids_tested") or ()
    if isinstance(rsids, str):
        rsids = tuple(r.strip() for r in rsids.split(";") if r.strip())

    genotypes = row.get("genotypes") or {}
    if isinstance(genotypes, str):
        genotypes = dict(
            part.split("=", 1) for part in genotypes.split(";") if "=" in part
        )

    toxicity = tuple(
        ToxicityEvent(
            term=t.get("term", "unspecified"),
            grade=int(t.get("grade", 0)),
            onset_cycle=t.get("onset_cycle"),
        )
        for t in (row.get("toxicity") or ())
    )

    if "anonymised_id" not in row or not str(row["anonymised_id"]).strip():
        raise LedgerRefused(
            "no anonymised_id column. Every row needs a site-assigned identifier "
            "so duplicates can be detected and a row can be superseded by a "
            "correction later. A bare KeyError here told a site nothing about "
            "which column was missing or why it matters."
        )

    return LedgerEntry(
        anonymised_id=str(row["anonymised_id"]),
        cancer_type=str(row.get("cancer_type", "")),
        regimen=str(row.get("regimen", "")),
        drug=str(row.get("drug", "")),
        assay_method=str(row.get("assay_method", "")),
        rsids_tested=tuple(rsids),
        genotypes=dict(genotypes),
        cpic_phenotype_at_dosing=str(row.get("cpic_phenotype_at_dosing", "")),
        engine_version=str(row.get("engine_version", "")),
        toxicity=toxicity,
        starting_dose_mg_m2=row.get("starting_dose_mg_m2"),
        bsa_m2=row.get("bsa_m2"),
        dose_modifications=tuple(row.get("dose_modifications") or ()),
        cycles_received=row.get("cycles_received"),
    )


def append(entry: LedgerEntry, path: Path) -> Path:
    """Append one entry as a JSON line. Append-only by construction.

    JSONL is chosen because a corrected row is a *new* row: nothing is rewritten,
    so the file is its own audit trail.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry.as_dict(), sort_keys=True) + "\n")
    return path


def read_all(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [
        json.loads(line) for line in path.read_text().splitlines() if line.strip()
    ]


def latest_per_patient(entries: list[dict]) -> list[dict]:
    """Collapse an append-only log to one row per patient: the last recorded.

    The file is append-only by design — a correction is a *new* row, so the file
    is its own audit trail. The consequence is that reading it back raw
    double-counts any patient who was entered twice, and a site correcting an
    outcome is the normal reason for that. In the observed failure the same
    ``anonymised_id`` appeared with grade 3 and then grade 0, and ``quadrant``
    counted them as two patients in two different cells.

    Later rows supersede earlier ones, ordered by ``recorded_at`` and falling back
    to file order. This is the only correct way to read the log for analysis, and
    ``quadrant`` now applies it by default.
    """
    by_id: dict[str, tuple[int, dict]] = {}
    for i, entry in enumerate(entries):
        pid = str(entry.get("anonymised_id", ""))
        key = (entry.get("recorded_at") or "", i)
        prior = by_id.get(pid)
        if prior is None or key > (prior[1].get("recorded_at") or "", prior[0]):
            by_id[pid] = (i, entry)
    return [entry for _, entry in sorted(by_id.values(), key=lambda p: p[0])]


def quadrant(
    entries: list[dict],
    panel_rsids: tuple[str, ...],
    collapse_duplicates: bool = True,
) -> dict[str, int]:
    """Four-quadrant counts: panel status x grade >=3 toxicity.

    Only counts entries where every panel locus was actually tested. A patient
    whose assay skipped a panel variant cannot be called panel-negative, and
    silently counting them as such is the error this package exists to prevent.

    A genotype token this function does not recognise is also **excluded**, not
    read as negative. The earlier implementation tested membership of a short
    literal set and treated everything else as a variant call, which made a
    blank string a *positive* and would have read ``"0|0"`` — legal VCF phased
    wild-type — as a variant too. Both directions are wrong, and unreadable is
    a third state rather than either.

    ``collapse_duplicates`` applies :func:`latest_per_patient` first, which is
    the default because the log is append-only: a corrected row is a new row, so
    a patient entered twice would otherwise be counted twice — and in the
    observed case, in two contradictory cells. ``duplicates_superseded`` reports
    how many rows were collapsed, so the number is visible rather than silently
    applied. Pass ``False`` only to audit the raw log.
    """
    superseded = 0
    if collapse_duplicates:
        collapsed = latest_per_patient(entries)
        superseded = len(entries) - len(collapsed)
        entries = collapsed

    cells = {
        "positive_toxic": 0,
        "positive_not_toxic": 0,
        "negative_toxic": 0,
        "negative_not_toxic": 0,
        "excluded_incomplete_panel": 0,
        "excluded_unreadable_genotype": 0,
        "duplicates_superseded": superseded,
    }
    # Computed once: the lookup reads the pinned CPIC table per locus.
    deciding = set(actionable_loci("DPYD", panel_rsids))

    for e in entries:
        tested = set(e.get("rsids_tested") or [])
        if not set(panel_rsids).issubset(tested):
            cells["excluded_incomplete_panel"] += 1
            continue
        gts = e.get("genotypes") or {}
        tokens = [str(gts.get(r, "")).strip().lower() for r in panel_rsids]
        if any(t not in _READABLE_TOKENS for t in tokens):
            cells["excluded_unreadable_genotype"] += 1
            continue
        # Only an actionable locus decides the split. A call at a locus CPIC
        # calls normal function leaves the patient in the group the published
        # comparators describe, so counting it as positive would move the most
        # common carriers out of the denominator those rates are compared to.
        positive = any(
            t in VARIANT_TOKENS
            for r, t in zip(panel_rsids, tokens, strict=True)
            if r in deciding
        )
        toxic = int(e.get("max_toxicity_grade") or 0) >= 3
        key = f"{'positive' if positive else 'negative'}_{'toxic' if toxic else 'not_toxic'}"
        cells[key] += 1
    return cells
