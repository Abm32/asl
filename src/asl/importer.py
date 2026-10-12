"""Read a genotype table the way a laboratory actually wrote it.

``asl`` interprets a DPYD result; this module is how that result gets in. It
exists because the artifact an Indian site hands over is almost never a VCF. It
is a 4-variant PCR report, transcribed into a spreadsheet by a person, with
column headers naming alleles rather than rsIDs and genotype cells reading
"Heterozygous" rather than ``1``.

Three jobs, and the third is the one that matters:

* **Locus synonyms.** ``rs56038477``, ``c.1236G>A``, ``HapB3`` and ``DPYD*2A``
  are all headers a real sheet has used. They are resolved to one rsID.
* **File formats.** CSV, TSV and XLSX. XLSX is read with the standard library
  (it is a zip of XML) rather than by adding a dependency.
* **Refusing an ambiguous genotype.** This is the safety surface. A cell
  reading "Detected", "Positive" or "Mutant" states that a variant is present
  and does **not** state how many copies. One copy of ``*2A`` is an
  Intermediate Metabolizer at 50% dose; two copies is a Poor Metabolizer for
  whom CPIC contraindicates the drug outright. Guessing heterozygous — the
  likelier of the two — would turn a contraindication into a dose reduction in
  the patient with the most to lose. So those words are refused by name, with
  the reason, rather than read as ``het``.

The same principle governs blank cells, and it is inherited from ``cli``: a
blank is ``NOT_TESTED``, never a negative. A header that is present with an
empty column means the assay covered nothing there; it does not mean zero.
"""

from __future__ import annotations

import csv
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree

from .interpret import ASSAYED_VARIANT, CPIC_ALLELE
from .models import Zygosity


def locus_key(text: str) -> str:
    """Normalise a header or cell for matching: casefold, strip non-alphanumerics.

    ``DPYD *2A``, ``dpyd*2a`` and ``DPYD_2A`` are the same column. Punctuation
    is dropped rather than mapped because a spreadsheet's punctuation is noise —
    except in cDNA notation, where ``c.1129-5923C>G`` and ``c.1236G>A`` differ
    in their digits, which survive.
    """
    return re.sub(r"[^a-z0-9]+", "", text.casefold())


def _locus_synonyms() -> dict[str, str]:
    """Build the header-to-rsID map from the pinned allele names.

    Derived rather than hand-written, so a locus added to ``CPIC_ALLELE`` is
    reachable by its CPIC name the moment it exists, and the two cannot drift.
    For each rsID this accepts the rsID itself, the full CPIC allele label, the
    cDNA part, and the star name if there is one.

    The HapB3 pair is the one place a hand-written entry is needed, and it is
    the place that most matters. CPIC's label for the tag locus rs56038477 is
    the compound ``"c.1129-5923C>G, c.1236G>A (HapB3)"``, which *names the
    causal variant* — so deriving synonyms from it alone would map a column
    headed ``c.1129-5923C>G`` onto the tag SNP, silently recording a causal
    genotype as a tag genotype and erasing the distinction the engine's
    ``hapb3_basis`` exists to carry. The two are pinned apart explicitly.
    """
    out: dict[str, str] = {}

    def add(text: str, rsid: str) -> None:
        k = locus_key(text)
        if not k:
            return
        prior = out.get(k)
        if prior is not None and prior != rsid:
            raise ValueError(
                f"header {text!r} is ambiguous between {prior} and {rsid}"
            )
        out[k] = rsid

    for rsid, label in CPIC_ALLELE.items():
        add(rsid, rsid)
        if rsid in _HAPB3:
            continue
        add(label, rsid)
        # "c.1905+1G>A (*2A)" -> "c.1905+1G>A" and "*2A", plus "DPYD*2A".
        cdna, _, star = label.partition(" (")
        add(cdna, rsid)
        if star:
            star = star.rstrip(")")
            add(star, rsid)
            add("DPYD" + star, rsid)

    for rsid, names in _HAPB3.items():
        for name in names:
            add(name, rsid)
        add(ASSAYED_VARIANT[rsid], rsid)
    return out


#: The HapB3 pair, pinned by hand. See ``_locus_synonyms``.
_HAPB3: dict[str, tuple[str, ...]] = {
    "rs75017182": (
        "c.1129-5923C>G",
        "HapB3 causal",
        "HapB3 intronic",
        "DPYD intron 10",
    ),
    "rs56038477": (
        "c.1236G>A",
        "HapB3",
        "DPYD HapB3",
        "HapB3 tag",
        "E412E",
        "p.Glu412=",
    ),
}

LOCUS_SYNONYM: dict[str, str] = _locus_synonyms()

#: Cell text that states a zygosity unambiguously.
_ZYGOSITY: dict[str, Zygosity] = {
    locus_key(text): zyg
    for text, zyg in (
        ("0", Zygosity.HOM_REF),
        ("hom_ref", Zygosity.HOM_REF),
        ("homref", Zygosity.HOM_REF),
        ("wt", Zygosity.HOM_REF),
        ("wild type", Zygosity.HOM_REF),
        ("wildtype", Zygosity.HOM_REF),
        ("normal", Zygosity.HOM_REF),
        ("negative", Zygosity.HOM_REF),
        ("not detected", Zygosity.HOM_REF),
        ("nd", Zygosity.HOM_REF),
        ("absent", Zygosity.HOM_REF),
        ("0/0", Zygosity.HOM_REF),
        ("ref/ref", Zygosity.HOM_REF),
        ("1", Zygosity.HET),
        ("het", Zygosity.HET),
        ("heterozygous", Zygosity.HET),
        ("heterozygote", Zygosity.HET),
        ("0/1", Zygosity.HET),
        ("1/0", Zygosity.HET),
        ("ref/var", Zygosity.HET),
        ("2", Zygosity.HOM_VAR),
        ("hom", Zygosity.HOM_VAR),
        ("hom_var", Zygosity.HOM_VAR),
        ("homvar", Zygosity.HOM_VAR),
        ("homozygous variant", Zygosity.HOM_VAR),
        ("homozygous mutant", Zygosity.HOM_VAR),
        ("homozygous", Zygosity.HOM_VAR),
        ("homozygote", Zygosity.HOM_VAR),
        ("1/1", Zygosity.HOM_VAR),
        ("var/var", Zygosity.HOM_VAR),
    )
}

#: Cell text asserting a variant is present without saying how many copies.
#: Refused by name rather than read as heterozygous — see the module docstring.
_COPY_UNSTATED: frozenset[str] = frozenset(
    locus_key(t)
    for t in (
        "detected",
        "positive",
        "present",
        "mutant",
        "variant",
        "carrier",
        "abnormal",
        "+",
        "yes",
    )
)


class ImportRefused(ValueError):
    """Raised when a cell or header cannot be read without guessing."""


def read_zygosity(text: str, locus: str = "") -> Zygosity | None:
    """Read one genotype cell. ``None`` for a blank — meaning NOT_TESTED.

    Raises ``ImportRefused`` for text that is either unrecognised or
    recognised-but-ambiguous. The two get different messages because they need
    different fixes: an unrecognised cell is a transcription the vocabulary has
    not seen, while an ambiguous one is a lab report that genuinely did not
    record the copy number and has to be gone back to.
    """
    raw = (text or "").strip()
    if not raw:
        return None
    k = locus_key(raw)
    zyg = _ZYGOSITY.get(k)
    if zyg is not None:
        return zyg
    where = f"{locus}: " if locus else ""
    if k in _COPY_UNSTATED:
        raise ImportRefused(
            f"{where}{raw!r} states that a variant is present but not how many "
            f"copies, and asl will not assume one. For a no-function allele "
            f"the difference is CPIC's entire recommendation: one copy is an "
            f"Intermediate Metabolizer at 50% of the starting dose, two copies "
            f"is a Poor Metabolizer for whom CPIC contraindicates the drug. "
            f"Guessing heterozygous — the likelier of the two — would turn a "
            f"contraindication into a dose reduction. Record het / hom_var "
            f"(or 1 / 2) from the laboratory's own report."
        )
    raise ImportRefused(
        f"{where}cannot read genotype {raw!r}. Use 0/1/2, "
        f"hom_ref/het/hom_var, wt/het/hom, or the words a report prints "
        f"(wild type, heterozygous, homozygous). Leave the cell blank only if "
        f"the assay returned no call at this locus — a blank is NOT_TESTED, "
        f"not a negative."
    )


def row_genotypes(row: dict[str, str]) -> dict[str, Zygosity | None]:
    """Read the locus columns of one table row, keyed by rsID.

    A column whose header is not a known DPYD locus is ignored — that is how
    ``anonymised_id``, ``regimen`` and a site's own free-text columns pass
    through. A column that *is* a locus and is blank maps to ``None``: the
    sheet carried the locus and the assay returned no call there.

    Two headers naming the same locus are refused rather than merged. A sheet
    with both ``rs56038477`` and ``HapB3`` columns has two readings of one
    position, and choosing between them is exactly what this package does not
    do.
    """
    out: dict[str, Zygosity | None] = {}
    seen: dict[str, str] = {}
    for header, value in row.items():
        rsid = LOCUS_SYNONYM.get(locus_key(header or ""))
        if rsid is None:
            continue
        prior = seen.get(rsid)
        if prior is not None:
            raise ImportRefused(
                f"columns {prior!r} and {header!r} both name {rsid} "
                f"({CPIC_ALLELE[rsid]}). That is two observations of one "
                f"position, and asl will not choose between them. Remove one."
            )
        seen[rsid] = header or ""
        out[rsid] = read_zygosity(value, CPIC_ALLELE[rsid])
    return out


# ---------------------------------------------------------------------------
# File formats
# ---------------------------------------------------------------------------

_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_REL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _column(ref: str) -> int:
    """Zero-based column index from an A1-style cell reference."""
    n = 0
    for ch in ref:
        if not ch.isalpha():
            break
        n = n * 26 + (ord(ch.upper()) - 64)
    return n - 1


def _xlsx_rows(path: Path) -> list[list[str]]:
    """Read the first worksheet of an XLSX as a rectangular list of strings.

    Written against the standard library because the alternative is a
    dependency for one file format, and the format is a zip of XML.

    The detail that has to be right is the **empty cell**. XLSX omits an empty
    cell from the row entirely, so reading cells in document order shifts every
    later value one column left — which here would move a genotype onto the
    wrong locus. Each cell's ``r`` attribute (``C4``) is decoded to a column
    index and written at that index, so an omitted cell stays blank where it
    belongs. That is the whole reason this is not four lines.
    """
    with zipfile.ZipFile(path) as z:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ElementTree.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall(f"{_NS}si"):
                shared.append("".join(t.text or "" for t in si.iter(f"{_NS}t")))

        sheets = [n for n in z.namelist() if n.startswith("xl/worksheets/sheet")]
        if not sheets:
            raise ImportRefused(f"{path}: no worksheet in this workbook")
        root = ElementTree.fromstring(z.read(min(sheets)))

        rows: list[list[str]] = []
        for row in root.iter(f"{_NS}row"):
            cells: list[str] = []
            for c in row.findall(f"{_NS}c"):
                idx = _column(c.get("r", "")) if c.get("r") else len(cells)
                if idx < 0:
                    idx = len(cells)
                while len(cells) <= idx:
                    cells.append("")
                if c.get("t") == "s":
                    v = c.find(f"{_NS}v")
                    i = int(v.text) if v is not None and v.text else -1
                    cells[idx] = shared[i] if 0 <= i < len(shared) else ""
                elif c.get("t") == "inlineStr":
                    cells[idx] = "".join(t.text or "" for t in c.iter(f"{_NS}t"))
                else:
                    v = c.find(f"{_NS}v")
                    cells[idx] = (v.text or "") if v is not None else ""
            rows.append(cells)
        return rows


def read_rows(path: str | Path) -> list[dict[str, str]]:
    """Read a genotype table into a list of header-keyed rows.

    Dispatches on the file extension: ``.xlsx``, ``.tsv``, anything else CSV.
    An ``.xls`` (the old binary format) is refused by name rather than handed
    to the CSV reader, which would read its binary header as a row of
    gibberish column names and then report every locus as unstated.
    """
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix == ".xls":
        raise ImportRefused(
            f"{p}: .xls is the pre-2007 binary Excel format, which asl does "
            f"not read. Open it and 'Save As' .xlsx or .csv — and when saving "
            f"as CSV check that blank cells stayed blank, because a blank is "
            f"how asl records 'the assay returned no call here'."
        )
    if suffix == ".xlsx":
        rows = _xlsx_rows(p)
        if not rows:
            return []
        header = [h.strip() for h in rows[0]]
        out: list[dict[str, str]] = []
        for cells in rows[1:]:
            if not any(c.strip() for c in cells):
                continue  # a trailing blank row Excel left behind
            padded = cells + [""] * (len(header) - len(cells))
            out.append(dict(zip(header, padded)))
        return out

    delimiter = "\t" if suffix == ".tsv" else ","
    with p.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


__all__ = [
    "LOCUS_SYNONYM",
    "ImportRefused",
    "locus_key",
    "read_rows",
    "read_zygosity",
    "row_genotypes",
]
