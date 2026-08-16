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

Where the data goes, once there is enough of it: publish, then submit to CPIC's
PCEP allele-function reassessment process (Tibben 2025, PMID 41175864 — it
requires published literature and a 70% consensus vote), to PharmVar for
haplotype registration, and to ClinPGx/PharmGKB, which accepts primary data.
That is the honest route to changing a call in the guideline rather than in our
own table — and changing it in our own table unilaterally is precisely what the
``pgx-core`` 0.7.2 release forbade.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

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
    (re.compile(r"[^@\s]+@[^@\s]+\.[a-z]{2,}", re.I), "an email address"),
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
        default_factory=lambda: datetime.now(timezone.utc)
    )

    def __post_init__(self) -> None:
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
        _reject_phi(self.as_dict())

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


def _reject_phi(payload: dict) -> None:
    """Refuse anything that looks like patient-identifying information."""
    for key, value in payload.items():
        if key.lower() in _PHI_FIELDS:
            raise LedgerRefused(
                f"field {key!r} is not accepted: the ledger holds no PHI. "
                f"Use an anonymised identifier assigned at the site."
            )
        if isinstance(value, str):
            for pattern, what in _PHI_PATTERNS:
                if pattern.search(value):
                    raise LedgerRefused(
                        f"field {key!r} contains what looks like {what}. "
                        f"The ledger holds no PHI."
                    )


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


def quadrant(entries: list[dict], panel_rsids: tuple[str, ...]) -> dict[str, int]:
    """Four-quadrant counts: panel status x grade >=3 toxicity.

    Only counts entries where every panel locus was actually tested. A patient
    whose assay skipped a panel variant cannot be called panel-negative, and
    silently counting them as such is the error this package exists to prevent.
    """
    cells = {
        "positive_toxic": 0,
        "positive_not_toxic": 0,
        "negative_toxic": 0,
        "negative_not_toxic": 0,
        "excluded_incomplete_panel": 0,
    }
    for e in entries:
        tested = set(e.get("rsids_tested") or [])
        if not set(panel_rsids).issubset(tested):
            cells["excluded_incomplete_panel"] += 1
            continue
        gts = e.get("genotypes") or {}
        positive = any(
            str(gts.get(r, "0")).strip() not in ("0", "0/0", "hom_ref", "")
            for r in panel_rsids
        )
        toxic = int(e.get("max_toxicity_grade") or 0) >= 3
        key = f"{'positive' if positive else 'negative'}_{'toxic' if toxic else 'not_toxic'}"
        cells[key] += 1
    return cells
