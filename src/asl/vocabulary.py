"""Which phenotype words are real, for which gene, at which layer.

This module exists because a single flat set of phenotype strings is wrong in
four separate ways, each measured against the live CPIC API and recorded in
``data/CPIC_VOCABULARY.json`` under ``_measurements``:

1. CPIC's standardised ``term`` table defines **11** phenotype terms; its
   ``diplotype`` table emits **33** distinct values; **24** of those are absent
   from the standardised list.
2. There is **not one exact-case match** between the two. The term table writes
   ``Normal metabolizer``, the diplotype table ``Normal Metabolizer``. Matching
   literally on CPIC's own published spelling rejects every real CPIC result,
   so every comparison here casefolds.
3. The vocabulary is **per gene**, and genes do not even answer on the same
   dimension: DPYD has 3 terms, TPMT and NUDT15 5 each, CYP2C19 8; SLCO1B1
   answers in *Function* terms and G6PD in ``Deficient``/``Normal``/
   ``Variable``. A global set can only ever be right for the one gene it was
   written against — here, DPYD.
4. ``No Result`` — a gene that was not assayed — appears **251 times** as a
   ``recommendation.phenotypes`` key and **never** as a ``diplotype.generesult``.
   CPIC models it at the recommendation layer alone. So a gene has *two*
   vocabularies, and the recommendation one is strictly larger.

The practical consequence for a report: a term may be one CPIC uses without
being one CPIC has defined. ``Possible Intermediate Metabolizer`` — which both
TPMT and NUDT15 emit — is such a term. This module will not invent a gloss for
it. ``gloss()`` returns CPIC's own ``functionaldef`` or nothing at all, and
``has_standardised_definition()`` lets a caller say "CPIC publishes no
definition for this term" instead of quietly printing the phenotype bare, which
reads as though nothing were missing.
"""

from __future__ import annotations

import json
from functools import cache
from pathlib import Path

_DATA = Path(__file__).resolve().parent / "data" / "CPIC_VOCABULARY.json"

#: CPIC's own phenotype value for a gene the report did not establish. It is
#: not "normal", it is not our word for the gap, and it is valid only where a
#: recommendation is being looked up.
NOT_ESTABLISHED = "No Result"


class VocabularyError(RuntimeError):
    """Raised when the vocabulary snapshot is missing or unusable."""


@cache
def snapshot() -> dict:
    if not _DATA.exists():
        raise VocabularyError(
            f"no CPIC vocabulary snapshot at {_DATA}. "
            f"Run: python scripts/regen_cpic_vocabulary.py"
        )
    data = json.loads(_DATA.read_text())
    if not data.get("gene_phenotypes"):
        raise VocabularyError(f"{_DATA} carries no gene phenotypes")
    return data


def measurements() -> dict:
    """The four measurements this module's shape is derived from."""
    return dict(snapshot()["_measurements"])


def genes() -> tuple[str, ...]:
    return tuple(snapshot()["gene_phenotypes"])


@cache
def _folded(gene: str, layer: str) -> dict[str, str]:
    """casefolded value -> the spelling CPIC publishes for this gene."""
    data = snapshot()
    if layer == "diplotype":
        values = list(data["gene_phenotypes"].get(gene, ()))
        # Terms the engine may return that CPIC does not emit as a diplotype
        # result. Widened deliberately and recorded in the snapshot, rather
        # than mixed into the gene set as though CPIC had published them.
        values += list(data.get("engine_terms", ()))
    elif layer == "recommendation":
        values = list(data["recommendation_phenotypes"].get(gene, ()))
        values += list(data["gene_phenotypes"].get(gene, ()))
        values += list(data.get("engine_terms", ()))
    else:
        raise ValueError(
            f"unknown layer {layer!r}. Use 'diplotype' or 'recommendation'."
        )
    return {v.casefold(): v for v in values}


def phenotypes_for(gene: str, *, layer: str = "diplotype") -> frozenset[str]:
    """Every phenotype this gene may carry at this layer, CPIC's spelling.

    ``layer='diplotype'`` is what a genotype may resolve to.
    ``layer='recommendation'`` additionally admits ``No Result``, because a
    recommendation may be keyed on a gene nobody assayed.
    """
    return frozenset(_folded(gene, layer).values())


def canonical(gene: str, value: str, *, layer: str = "diplotype") -> str | None:
    """CPIC's spelling of `value` for this gene, or None if it is not one.

    Matching is casefolded and whitespace-tolerant. It has to be: CPIC's own
    standardised term list and its own result table disagree on capitalisation
    for every single term, so a literal comparison is guaranteed to fail.
    """
    if value is None:
        return None
    return _folded(gene, layer).get(" ".join(str(value).split()).casefold())


def is_recognised(gene: str, value: str, *, layer: str = "diplotype") -> bool:
    return canonical(gene, value, layer=layer) is not None


#: Which half of CPIC's term table a lookup means.
ALLELE_FUNCTION = "allele"
PHENOTYPE = "phenotype"


@cache
def _terms_folded(kind: str) -> dict[str, dict]:
    """casefolded term -> its row, within one half of the term table.

    Keyed by kind because three of CPIC's seventeen standardised terms are
    ambiguous without it. ``Normal function``, ``Decreased function`` and
    ``Increased function`` each appear twice: once under "Allele functional
    status: all genes" and once under "Phenotype: transporters", with
    different definitions. ``Normal function`` is "Fully functional/wild-type"
    as an allele and "Fully functional transporter function" as a phenotype.
    A lookup keyed on the word alone silently returns whichever row happened
    to be last, so an allele's function can be glossed with a transporter
    phenotype's definition and nothing says so.
    """
    if kind not in (ALLELE_FUNCTION, PHENOTYPE):
        raise ValueError(
            f"unknown kind {kind!r}. Use {ALLELE_FUNCTION!r} or {PHENOTYPE!r}."
        )
    wanted = "Allele functional status" if kind == ALLELE_FUNCTION else "Phenotype"
    return {
        t["term"].casefold(): t
        for t in snapshot()["terms"]
        if (t.get("category") or "").startswith(wanted)
    }


def ambiguous_terms() -> tuple[str, ...]:
    """Terms CPIC defines twice, with different meanings per category."""
    allele = set(_terms_folded(ALLELE_FUNCTION))
    phenotype = set(_terms_folded(PHENOTYPE))
    return tuple(
        sorted(_terms_folded(ALLELE_FUNCTION)[t]["term"] for t in allele & phenotype)
    )


def standardised_terms() -> tuple[dict, ...]:
    """CPIC's consensus term list (Caudle 2017, PMID 27441996), verbatim."""
    return tuple(snapshot()["terms"])


def _row(value: str, kind: str) -> dict:
    if value is None:
        return {}
    return _terms_folded(kind).get(" ".join(str(value).split()).casefold()) or {}


def has_standardised_definition(value: str, *, kind: str = PHENOTYPE) -> bool:
    """Did CPIC's consensus term list define this word at all?

    False for terms CPIC uses but never standardised — ``Indeterminate``,
    ``Possible Intermediate Metabolizer``, ``Likely Poor Metabolizer`` and
    twenty-one others. A report that knows this can say so.
    """
    return bool(_row(value, kind).get("functional_definition"))


def gloss(value: str, *, kind: str = PHENOTYPE) -> str:
    """CPIC's own functional definition of a term, or "".

    Never invented. An empty string means CPIC has not defined the term, which
    is a fact about the term and worth printing as one — see
    ``has_standardised_definition``. Pass ``kind=ALLELE_FUNCTION`` to define an
    allele's function rather than a phenotype; three terms differ between the
    two and returning the wrong one is silent.
    """
    return _row(value, kind).get("functional_definition", "")


def genetic_definition(value: str, *, kind: str = PHENOTYPE) -> str:
    """CPIC's own genotype-side definition of a term, or "".

    CPIC writes "N/A" where the term has no genotype-side definition — every
    allele-function row does. That is an absence, so it is returned as one
    rather than printed as the literal string "N/A".
    """
    definition = _row(value, kind).get("genetic_definition", "")
    return "" if definition.strip().upper() == "N/A" else definition


def undefined_terms_for(gene: str, *, layer: str = "diplotype") -> tuple[str, ...]:
    """This gene's phenotypes that CPIC uses but never standardised."""
    return tuple(
        sorted(
            v
            for v in phenotypes_for(gene, layer=layer)
            if not has_standardised_definition(v)
        )
    )
