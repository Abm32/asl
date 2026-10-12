"""The genotype table is where a wrong answer enters, so it is where asl refuses.

Everything downstream of ``importer`` is deterministic and audited. A defect
here is different in kind: it does not produce a refusal or a visible error, it
produces a *correct interpretation of the wrong genotype*, which looks exactly
like a correct report.

These tests pin the three ways that can happen.
"""

from __future__ import annotations

import zipfile

import pytest

from asl.importer import (
    LOCUS_SYNONYM,
    ImportRefused,
    locus_key,
    read_rows,
    read_zygosity,
    row_genotypes,
)
from asl.models import Zygosity

# ---------------------------------------------------------------------------
# 1. A cell that says a variant is present without saying how many copies.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "text", ["Detected", "POSITIVE", "present", "Mutant", "carrier", "Yes"]
)
def test_a_copy_unstated_cell_is_refused_not_read_as_heterozygous(text) -> None:
    """The one that would kill someone.

    A lab report reading "DPYD *2A: Detected" states presence and not dosage.
    One copy is an Intermediate Metabolizer at 50% of the starting dose; two
    copies is a Poor Metabolizer for whom CPIC contraindicates fluoropyrimidines
    outright. Heterozygous is the likelier of the two by a wide margin, which is
    precisely what makes guessing it attractive and wrong: the error is rare,
    silent, and lands on the only patient for whom the drug is lethal.

    So the word is refused by name. The refusal names the consequence too,
    because "cannot parse" invites someone to type 1 and move on.
    """
    with pytest.raises(ImportRefused) as exc:
        read_zygosity(text, "c.1905+1G>A (*2A)")
    assert "not how many copies" in str(exc.value)
    assert "Poor Metabolizer" in str(exc.value)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("Heterozygous", Zygosity.HET),
        ("het", Zygosity.HET),
        ("0/1", Zygosity.HET),
        ("1", Zygosity.HET),
        ("Wild type", Zygosity.HOM_REF),
        ("Not detected", Zygosity.HOM_REF),
        ("0", Zygosity.HOM_REF),
        ("Homozygous mutant", Zygosity.HOM_VAR),
        ("hom_var", Zygosity.HOM_VAR),
        ("2", Zygosity.HOM_VAR),
    ],
)
def test_the_words_a_laboratory_actually_prints_are_read(text, expected) -> None:
    assert read_zygosity(text) is expected


def test_a_blank_cell_is_a_no_call_and_not_a_negative() -> None:
    """The distinction the whole package is built on, at its entry point."""
    assert read_zygosity("") is None
    assert read_zygosity("   ") is None
    assert read_zygosity("0") is Zygosity.HOM_REF


# ---------------------------------------------------------------------------
# 2. Header synonyms, and the HapB3 trap.
# ---------------------------------------------------------------------------


def test_the_two_hapb3_loci_do_not_collapse_into_each_other() -> None:
    """CPIC's label for the TAG locus names the CAUSAL variant.

    ``CPIC_ALLELE["rs56038477"]`` is ``"c.1129-5923C>G, c.1236G>A (HapB3)"`` —
    the compound haplotype name. Deriving header synonyms from that string
    alone would map a column headed ``c.1129-5923C>G`` onto rs56038477, so a
    laboratory that genotyped the causal intronic variant would have its result
    recorded as a tag-SNP result.

    The report would then print CPIC's tag-SNP disclosure — "decreased function
    was inferred by detecting the exonic tag SNP" — about a patient in whom the
    causal variant was directly observed. That is a false caveat, and the
    mirror-image failure (a tag result recorded as causal) suppresses a caveat
    CPIC requires. The two are pinned apart.
    """
    assert LOCUS_SYNONYM[locus_key("c.1129-5923C>G")] == "rs75017182"
    assert LOCUS_SYNONYM[locus_key("c.1236G>A")] == "rs56038477"
    assert LOCUS_SYNONYM[locus_key("HapB3")] == "rs56038477"
    assert LOCUS_SYNONYM[locus_key("HapB3 causal")] == "rs75017182"


@pytest.mark.parametrize(
    ("header", "rsid"),
    [
        ("rs3918290", "rs3918290"),
        ("DPYD*2A", "rs3918290"),
        ("*2A", "rs3918290"),
        ("c.1905+1G>A", "rs3918290"),
        ("c.1905+1G>A (*2A)", "rs3918290"),
        ("c.2279C>T", "rs112766203"),
        ("c.2846A>T", "rs67376798"),
    ],
)
def test_a_header_is_resolved_however_the_site_wrote_it(header, rsid) -> None:
    assert LOCUS_SYNONYM[locus_key(header)] == rsid


def test_a_non_locus_column_is_ignored_rather_than_guessed_at() -> None:
    out = row_genotypes(
        {"anonymised_id": "MCC-001", "regimen": "systemic", "rs3918290": "0"}
    )
    assert out == {"rs3918290": Zygosity.HOM_REF}


def test_two_headers_naming_one_locus_are_refused() -> None:
    """Two readings of one position. Choosing between them is not asl's job.

    ``AssayResult`` already refuses a duplicated rsID, but by then both columns
    have been read and one has silently overwritten the other — the refusal
    fires on the survivor, not on the conflict. Catching it at the header means
    the message can name both columns, which is what the person fixing the
    spreadsheet needs.
    """
    with pytest.raises(ImportRefused) as exc:
        row_genotypes({"rs56038477": "0", "HapB3": "1"})
    assert "rs56038477" in str(exc.value)


# ---------------------------------------------------------------------------
# 3. XLSX, and the cell that is not there.
# ---------------------------------------------------------------------------


def _xlsx(path, header: list[str], rows: list[dict[str, str]]) -> None:
    """Write a minimal XLSX, OMITTING cells whose value is blank.

    This is what Excel itself does: an empty cell is absent from the row XML
    rather than present and empty. Writing the fixture the same way is the
    point of the fixture.
    """
    strings: list[str] = []

    def sid(text: str) -> int:
        if text not in strings:
            strings.append(text)
        return strings.index(text)

    body = []
    cells = "".join(
        f'<c r="{chr(65 + i)}1" t="s"><v>{sid(h)}</v></c>'
        for i, h in enumerate(header)
    )
    body.append(f'<row r="1">{cells}</row>')
    for n, row in enumerate(rows, start=2):
        cells = "".join(
            f'<c r="{chr(65 + i)}{n}" t="s"><v>{sid(row[h])}</v></c>'
            for i, h in enumerate(header)
            if row.get(h, "")
        )
        body.append(f'<row r="{n}">{cells}</row>')

    shared = (
        '<?xml version="1.0"?><sst xmlns="http://schemas.openxmlformats.org/'
        'spreadsheetml/2006/main">'
        + "".join(f"<si><t>{s}</t></si>" for s in strings)
        + "</sst>"
    )
    sheet = (
        '<?xml version="1.0"?><worksheet xmlns="http://schemas.openxmlformats.'
        'org/spreadsheetml/2006/main"><sheetData>'
        + "".join(body)
        + "</sheetData></worksheet>"
    )
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("[Content_Types].xml", "<Types/>")
        z.writestr("xl/sharedStrings.xml", shared)
        z.writestr("xl/worksheets/sheet1.xml", sheet)


def test_an_omitted_xlsx_cell_stays_blank_in_its_own_column(tmp_path) -> None:
    """The defect a naive XLSX reader ships with.

    XLSX omits an empty cell from the row entirely. Reading cells in document
    order therefore shifts every later value one column to the left — and here
    the columns are loci, so a genotype lands on the wrong variant. A
    heterozygous HapB3 becomes a heterozygous *2A: same report shape, different
    allele, and *2A is no-function where HapB3 is decreased-function.

    The reader decodes each cell's ``r`` attribute (``C2``) to a column index
    instead, so the gap stays where it belongs.
    """
    path = tmp_path / "lab.xlsx"
    _xlsx(
        path,
        ["anonymised_id", "rs3918290", "HapB3", "max_toxicity_grade"],
        [{"anonymised_id": "P1", "rs3918290": "", "HapB3": "Heterozygous",
          "max_toxicity_grade": "3"}],
    )
    (row,) = read_rows(path)
    assert row["rs3918290"] == ""
    assert row["HapB3"] == "Heterozygous"

    genotypes = row_genotypes(row)
    assert genotypes["rs3918290"] is None  # NOT_TESTED, not a negative
    assert genotypes["rs56038477"] is Zygosity.HET


def test_csv_and_xlsx_of_the_same_table_agree(tmp_path) -> None:
    header = ["anonymised_id", "rs3918290", "rs56038477"]
    rows = [{"anonymised_id": "P1", "rs3918290": "0", "rs56038477": "het"}]

    xlsx = tmp_path / "t.xlsx"
    _xlsx(xlsx, header, rows)
    csv_path = tmp_path / "t.csv"
    csv_path.write_text("anonymised_id,rs3918290,rs56038477\nP1,0,het\n")

    assert row_genotypes(read_rows(xlsx)[0]) == row_genotypes(read_rows(csv_path)[0])


def test_the_old_binary_excel_format_is_refused_by_name(tmp_path) -> None:
    """Handing a .xls to the CSV reader produces gibberish headers, not an error.

    Every locus would then be unstated, ``_assay`` would refuse, and the
    refusal would talk about panel coverage rather than about the file format —
    sending the person to fix the wrong thing.
    """
    path = tmp_path / "old.xls"
    path.write_bytes(b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1")
    with pytest.raises(ImportRefused) as exc:
        read_rows(path)
    assert "Save As" in str(exc.value)
