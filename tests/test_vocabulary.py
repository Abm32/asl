"""The four measurements that forced this module's shape, as assertions.

These are not style tests. Each one encodes a defect that a single flat set of
phenotype strings has, and that this package had until 2026-10-25. The counts
are asserted against the snapshot's own recorded measurements rather than
against literals, so CPIC publishing a new term does not fail the suite —
``scripts/regen_cpic_vocabulary.py --check`` is what catches that, and it is
the right place for it.
"""

import pytest

from asl import vocabulary

# --- finding 1: the standardised list does not cover what CPIC emits -------


def test_cpic_emits_far_more_terms_than_it_standardised():
    m = vocabulary.measurements()
    assert m["distinct_emitted_generesults"] > m["standardised_phenotype_terms"]
    assert m["emitted_but_not_standardised"]


def test_the_terms_this_package_needs_most_are_among_the_undefined():
    """`Possible Intermediate Metabolizer` is emitted by both thiopurine genes
    and defined by nobody. So is `Indeterminate`."""
    undefined = set(vocabulary.measurements()["emitted_but_not_standardised"])
    assert "Possible Intermediate Metabolizer" in undefined
    assert "Indeterminate" in undefined


# --- finding 2: matching on CPIC's own spelling matches nothing ------------


def test_no_emitted_value_matches_a_standardised_term_exactly():
    """The single most important assertion here.

    CPIC's term table writes `Normal metabolizer`; its diplotype table writes
    `Normal Metabolizer`. A consumer that validates against CPIC's published
    vocabulary, spelled as CPIC publishes it, rejects every real CPIC result.
    """
    assert vocabulary.measurements()["exact_case_matches"] == []


def test_casefolded_lookup_recovers_the_spelling_cpic_actually_uses():
    assert vocabulary.canonical("DPYD", "normal metabolizer") == "Normal Metabolizer"
    assert vocabulary.canonical("DPYD", "  POOR   metabolizer ") == "Poor Metabolizer"
    assert vocabulary.gloss("Normal Metabolizer")  # the term table spells it lower


# --- finding 3: the vocabulary is per gene ---------------------------------


def test_dpyd_and_the_thiopurine_genes_do_not_share_a_vocabulary():
    dpyd = vocabulary.phenotypes_for("DPYD")
    tpmt = vocabulary.phenotypes_for("TPMT")
    assert "Possible Intermediate Metabolizer" in tpmt
    assert "Possible Intermediate Metabolizer" not in dpyd


def test_a_thiopurine_phenotype_is_not_a_dpyd_phenotype():
    assert vocabulary.is_recognised("TPMT", "Possible Intermediate Metabolizer")
    assert not vocabulary.is_recognised("DPYD", "Possible Intermediate Metabolizer")


def test_some_genes_do_not_answer_in_metabolizer_terms_at_all():
    """SLCO1B1 is a transporter and answers in Function terms. A set written
    against DPYD cannot be widened into correctness for it."""
    slco = vocabulary.phenotypes_for("SLCO1B1")
    assert any("Function" in v for v in slco)
    assert not any("Metabolizer" in v for v in slco)


# --- finding 4: "not established" lives at a different layer ---------------


def test_no_result_is_never_a_diplotype_outcome():
    assert vocabulary.measurements()["not_established_is_recommendation_layer_only"]
    assert not vocabulary.is_recognised("TPMT", vocabulary.NOT_ESTABLISHED)


def test_no_result_is_valid_where_a_recommendation_is_looked_up():
    assert vocabulary.is_recognised(
        "TPMT", vocabulary.NOT_ESTABLISHED, layer="recommendation"
    )
    assert vocabulary.NOT_ESTABLISHED in vocabulary.phenotypes_for(
        "TPMT", layer="recommendation"
    )


def test_the_recommendation_layer_is_a_superset_of_the_diplotype_layer():
    for gene in ("TPMT", "NUDT15"):
        assert vocabulary.phenotypes_for(gene) <= vocabulary.phenotypes_for(
            gene, layer="recommendation"
        )


def test_an_unknown_layer_is_refused_rather_than_defaulted():
    with pytest.raises(ValueError):
        vocabulary.phenotypes_for("DPYD", layer="whatever")


# --- the gloss never invents ----------------------------------------------


def test_the_gloss_is_cpics_own_words():
    gloss = vocabulary.gloss("Poor Metabolizer")
    assert gloss
    assert gloss in {t["functional_definition"] for t in vocabulary.standardised_terms()}


def test_a_term_cpic_never_defined_gets_no_gloss_rather_than_an_invented_one():
    assert vocabulary.gloss("Possible Intermediate Metabolizer") == ""
    assert not vocabulary.has_standardised_definition(
        "Possible Intermediate Metabolizer"
    )


def test_each_gene_can_name_the_terms_cpic_left_undefined():
    for gene in ("TPMT", "NUDT15"):
        assert "Possible Intermediate Metabolizer" in vocabulary.undefined_terms_for(
            gene
        )


def test_the_loinc_column_exists_and_is_empty():
    """CPIC's schema carries it; CPIC has not populated it. Recorded so its
    arrival is visible rather than assumed."""
    assert all(t["loinc"] == "" for t in vocabulary.standardised_terms())


def test_na_genetic_definitions_are_not_printed_as_na():
    """Every allele-function row carries CPIC's literal "N/A"."""
    assert (
        vocabulary.genetic_definition(
            "Normal function", kind=vocabulary.ALLELE_FUNCTION
        )
        == ""
    )


# --- finding 5, found by this package's own tests -------------------------


def test_three_standardised_terms_mean_two_different_things():
    """`Normal function` is an allele status AND a transporter phenotype, with
    different definitions. A lookup keyed on the word alone returns whichever
    row came last, and nothing says which one you got."""
    assert set(vocabulary.ambiguous_terms()) == {
        "Normal function",
        "Decreased function",
        "Increased function",
    }


def test_the_same_word_glosses_differently_by_kind():
    as_allele = vocabulary.gloss("Normal function", kind=vocabulary.ALLELE_FUNCTION)
    as_phenotype = vocabulary.gloss("Normal function", kind=vocabulary.PHENOTYPE)
    assert as_allele == "Fully functional/wild-type"
    assert "transporter" in as_phenotype
    assert as_allele != as_phenotype


def test_an_unknown_kind_is_refused_rather_than_defaulted():
    with pytest.raises(ValueError):
        vocabulary.gloss("Normal function", kind="whatever")


# --- the ledger uses it per gene ------------------------------------------


def _entry(**kw):
    from asl.ledger import LedgerEntry

    base = {
        "anonymised_id": "SITE-1",
        "cancer_type": "colorectal",
        "regimen": "CAPOX",
        "drug": "capecitabine",
        "assay_method": "qPCR",
        "rsids_tested": ("rs3918290",),
        "genotypes": {"rs3918290": "0/0"},
        "cpic_phenotype_at_dosing": "Normal Metabolizer",
        "engine_version": "0.11.0",
    }
    base.update(kw)
    return LedgerEntry(**base)


def test_the_ledger_still_accepts_a_dpyd_row_unchanged():
    assert _entry().gene == "DPYD"


def test_the_ledger_accepts_a_thiopurine_phenotype_for_a_thiopurine_gene():
    entry = _entry(
        gene="TPMT",
        drug="azathioprine",
        regimen="azathioprine",
        cpic_phenotype_at_dosing="Possible Intermediate Metabolizer",
    )
    assert entry.gene == "TPMT"


def test_the_ledger_refuses_that_same_phenotype_against_dpyd():
    from asl.ledger import LedgerRefused

    with pytest.raises(LedgerRefused, match="not a CPIC phenotype for DPYD"):
        _entry(cpic_phenotype_at_dosing="Possible Intermediate Metabolizer")


def test_the_ledger_accepts_cpics_other_capitalisation():
    """It has to: CPIC's own two tables disagree on it."""
    assert _entry(cpic_phenotype_at_dosing="normal metabolizer")
