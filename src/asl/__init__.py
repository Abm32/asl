"""asl — Anukriti Oncology Safety Ledger.

Reports what a pharmacogenomic result can and cannot rule out for a patient's
population, with the provenance and strength of every allele-function call
attached.

The deterministic layer decides; this package explains. It never overrides
``anukriti-pgx-core``, never recommends a dose, and has no code path that
returns an all-clear.
"""

from __future__ import annotations

from .models import (
    AssayMethod,
    AssayResult,
    Basis,
    Drug,
    EvidenceRecord,
    Flag,
    Interpretation,
    Regimen,
    ResidualRisk,
    TestStatus,
    VariantCall,
    Zygosity,
)
from . import vocabulary
from .pipeline import build
from .report import render

__version__ = "0.1.0"

__all__ = [
    "AssayMethod",
    "AssayResult",
    "Basis",
    "Drug",
    "EvidenceRecord",
    "Flag",
    "Interpretation",
    "Regimen",
    "ResidualRisk",
    "TestStatus",
    "VariantCall",
    "Zygosity",
    "__version__",
    "build",
    "render",
    # Which phenotype words are real, for which gene, at which layer. Exported
    # because consumers outside this package need it: prahari answers
    # thiopurines on TPMT and NUDT15, whose vocabulary is not DPYD's.
    "vocabulary",
]
