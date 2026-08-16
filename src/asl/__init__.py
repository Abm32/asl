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
    AssayResult,
    Basis,
    EvidenceRecord,
    Flag,
    Interpretation,
    Regimen,
    ResidualRisk,
    TestStatus,
    VariantCall,
    Zygosity,
)
from .pipeline import build
from .report import render

__version__ = "0.1.0"

__all__ = [
    "AssayResult",
    "Basis",
    "EvidenceRecord",
    "Flag",
    "Interpretation",
    "Regimen",
    "ResidualRisk",
    "TestStatus",
    "VariantCall",
    "Zygosity",
    "build",
    "render",
    "__version__",
]
