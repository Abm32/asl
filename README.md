# asl — Anukriti Oncology Safety Ledger

> **Status: research tool, pre-1.0, not clinically validated.** Not a diagnostic
> device. No regulatory clearance. The go/no-go site study described under
> [Status](#status) has not been done. Read
> [Honest limitations](#honest-limitations) before drawing conclusions from it.

**What it does.** Takes a DPYD pharmacogenomic result and reports what that
result *can and cannot rule out for this patient's population*, with the
provenance and evidence strength of every allele-function call attached. Then it
captures the toxicity outcome afterwards, so the calls can eventually be
corrected with evidence instead of asserted.

**What it does not do.** It does not recommend doses, does not override CPIC, and
has no code path that returns an all-clear. Those are enforced by tests, not by
convention.

```
pip install -e .
asl demo                                # the worked patient case
asl quadrant --csv fixtures/cohort_template.csv
```

Runs entirely offline. CPIC's data is pinned in `data/evidence/`, so a clinic
demonstration does not depend on connectivity.

### Where to start

| If you want | Read |
|---|---|
| the argument, in one table | [The problem, in one table](#the-problem-in-one-table) |
| a real report, without installing | [`out/demo_panel_negative.txt`](out/demo_panel_negative.txt), [`out/demo_with_sequencing.txt`](out/demo_with_sequencing.txt) |
| what it will not do, and why | [What it refuses to do](#what-it-refuses-to-do) |
| how the pieces fit | [Architecture](#architecture) |
| to run it on your own cohort | [Running it against a site's data](#running-it-against-a-sites-data) |
| to contribute | [`CONTRIBUTING.md`](CONTRIBUTING.md) |
| to report a wrong report | [`SECURITY.md`](SECURITY.md) — a misleading output is a reportable defect here |

### Repository layout

```
src/asl/
  models.py      the NOT_TESTED vs ABSENT distinction, as types
  interpret.py   the only call site into anukriti-pgx-core; never overrides it
  evidence.py    reads the pinned CPIC tables; refuses uncitable claims
  residual.py    what a result does not exclude; no "no risk" code path exists
  flags.py       named uncertainties (P2/P3 conflicts, P4 pending CPIC changes)
  pipeline.py    assembles the record; four preconditions before the engine runs
  report.py      renders the record and nothing else
  ledger.py      append-only outcome capture; recursive PHI refusal
  cli.py         asl demo | panels | interpret | quadrant
data/evidence/   generated from CPIC's API, never hand-edited (CI checks drift)
scripts/         regen_evidence_tables.py [--check]
tests/           140 tests; every refusal above has one
out/             worked outputs, regenerated on release
fixtures/        the CSV column layout to ask a site for
```

Sole dependency: `anukriti-pgx-core==0.7.3`, pinned exactly — a new engine
version carries CPIC table updates that need clinical regression review first.

---

## The problem, in one table

Computed from CPIC's own Central/South Asian allele frequencies under
Hardy-Weinberg — no Indian data required, and the command is
`asl demo` plus `tests/test_evidence_invariants.py`:

| | allele frequency (CPIC, Central/South Asian) |
|---|---|
| The 4 variants a guideline-conformant DPYD panel tests | **0.02538** |
| ⤷ of which `*13` alone | **0.00000** — cannot occur in this population |
| The 4 variants CPIC calls "Normal function" that Indians actually carry | **0.52921** (20.9×) |
| **⇒ patients reported "Normal Metabolizer"** | **95.0%** |
| **⇒ patients carrying ≥1 of those "Normal function" alleles** | **77.8%** |

The only Indian NGS cohort observed 71.1% (54/76, D-TORCH, AIIMS). An
expectation derived from the guideline's own reference tables lands within seven
points of what was measured.

So the test is not being run badly. **The interpretation layer converts "we found
nothing on our list" into "this patient is normal", and destroys the distinction
on the way.** That conversion is a software behaviour.

CPIC already says the test cannot do that:

> "patients without a DPYD decreased/no function variant may still experience
> severe toxicity due to other genetic, environmental or other factors"
>
> "a genetic test investigating only selected decreased/no function variants
> does not fully rule out DPD defects"
>
> — CPIC DPYD guideline, Amstutz et al. 2018, PMID 29152729

`asl` does not make a novel clinical claim. It declines to discard a caveat the
guideline publishes.

---

## What a report contains

Run `asl demo` to see it. Six sections, each with a greppable rule ID so any line
can be traced back to the record that produced it:

| Section | Rule ID | What it adds |
|---|---|---|
| CPIC call | — | diplotype, phenotype, CPIC's own action, attributed to CPIC |
| Tested vs not tested | `ASL-ASSAY-01` | **an untested locus is never reported as absent** |
| Uninformative loci | `ASL-POP-01` | `*13` is 0.00000 here — very low expected yield, **not** proof of absence |
| Evidence basis | `ASL-EVID-01` | strength, assay type, PMIDs, population frequency — per allele, read from CPIC |
| Named conflicts | `P2_*` / `P3_*` | where the literature disagrees, **stated in both directions** |
| Pending CPIC changes | `P4_*` | where **CPIC has announced a revision to its own value** but not published it |
| Residual risk | `ASL-RESID-01` | regimen-matched wild-type toxicity rate + panel sensitivity + the CPIC caveat |

### The one flag that no automated guard can catch

`P4_DPYD_HAPB3_CPIC_UPDATE_PENDING` exists because of a real gap. On 9 July 2026
CPIC announced that HapB3's allele value moves from **0.5 to 0.75** — a
heterozygote's recommendation going from 50% of standard dose to **75%, starting
cycle 1** — while stating that its database "will remain at AV=0.5 until then"
(full guideline expected fall 2026).

So the API returns 0.5, `data/evidence/` says 0.5, `--check` passes, and the fact
that the recommendation is about to move is invisible to every automated guard.
**In CPIC's own Central/South Asian data HapB3 is 77.5% of the entire actionable
allele frequency** (0.019658 of 0.025380) — an Indian DPYD panel is, to a first
approximation, a HapB3 test. So this lands on nearly every Indian patient who
tests positive.

`asl` reports CPIC's current value, flags the announced one alongside it, and
**does not pre-apply it**. Dosing on a pre-publication comment would be an
undeclared deviation from the pinned source of truth — the `*6` failure class,
pointed the other way.

Note also which direction this cuts: for the ~5% who test positive, a
European-derived activity score has been **over-restricting** them, while the
argument elsewhere in this README is that the panel **over-reassures** the
negatives. Both are the same defect — one number doing two jobs.

### The numbers it quotes, and why each

- **Wild-type grade ≥3 toxicity rate** — the number that answers "what is still
  true of me now that my test is negative?"
  - systemic chemotherapy: **22.7% (231/1018)** — Henricks 2018, PMID 30348537
  - including chemoradiation: **13.6% (105/771)** — Lunenburg 2018, PMID 30361102
  - These are different populations at different doses. **With no regimen
    supplied, `asl` refuses to quote a rate** rather than defaulting.
- **Panel clinical sensitivity: 3.5–21.6%** across 9 studies — Ontario HTA 2021,
  PMID 34484488. Reported as a range because the HTA computed no pooled estimate.
- **Deliberately not used:** the "20–30% of early-onset toxicity can be
  explained" figure. That is an attributable fraction, not a sensitivity, and
  conflating them overstates the panel roughly threefold. It appears only inside
  CPIC's own quoted context. Asserted by test.

---

## Architecture

```
CPIC live API ──(scripts/regen_evidence_tables.py --check, gated in CI)──┐
                                                                        ▼
anukriti-pgx-core ==0.7.3 ──▶ interpret ──▶ pipeline ──▶ report      data/evidence/
   (CPIC truth, decides)         │            │                      (pinned, offline)
                                 │            ├── evidence  (provenance per allele)
                                 │            ├── flags     (named conflicts)
                                 │            └── residual  (what is not excluded)
                                 │
                                 └── refuses: over-specified genotype,
                                     unstated panel content, uncitable allele
                            ledger ──▶ append-only JSONL ──▶ CPIC PCEP submission
```

**The evidence layer is generated, never hand-authored.** CPIC's API already
publishes `strength`, `findings` (which encodes assay type), `citations` and
per-population `frequency` per allele. `asl` copies all of it and edits none of
it; `--check` fails a build on drift. Only three fields are hand-maintained —
cohort size, cohort ancestry, and the conflict flag — because CPIC does not
publish them, and **none of the three may feed a phenotype** (asserted by test).

That discipline exists because of a real incident: `pgx-core` carried `*6` as
*Decreased function* against CPIC's *Normal function*. The clinical evidence
supported our value, but it was an **undeclared deviation from the pinned source
of truth**, and it survived four releases. It is now restored to CPIC's word,
with the disagreement recorded as a named conflict — which is what
`P2_DPYD_STAR6_INVITRO_CLINICAL_CONFLICT` is.

---

## What it refuses to do

Each is enforced by a test in `tests/test_safety_properties.py`.

1. **Does not recommend doses.** CPIC's recommendation is displayed as CPIC's.
   Crossing this line makes the software *drive* clinical management — Class B/C
   under CDSCO — rather than *inform* it.
2. **Does not invent a population-specific activity score.** No such values
   until the ledger supports them.
3. **Does not override CPIC.** Disagreements are named uncertainties *alongside*
   the call.
4. **Does not emit an all-clear.** No input produces one; checked exhaustively
   across every reachable zygosity combination on the actionable panel.
5. **Does not treat `NOT_TESTED` as absent.**
6. **Does not quote a toxicity rate without a regimen.**
7. **Does not claim a sub-population frequency.** CPIC's `Central/South Asian` is
   the finest granularity the data supports; a missing population omits the
   frequency rather than substituting European.
8. **Does not let an LLM touch a number.** (`narrate` is not yet implemented; when
   it is, it renders only what is in the record and is validated against it.)
9. **Does not accept PHI.** Field-name denylist plus value pattern checks on the
   CSV path, because that is how a spreadsheet column actually arrives.
10. **Does not resolve an over-specified genotype by dropping alleles.** If more
    than two variant copies are detected it refuses, because the engine's
    diplotype then omits alleles that were found.
11. **Does not quote a comparison group's rate to a patient outside it.** The
    22.7% / 13.6% figures are rates among patients *without* an actionable
    variant. A decreased/no-function carrier is not in that group, so no figure
    from it appears in their report — for a Poor Metabolizer it would understate
    the risk. A `*9A` carrier *does* stay in the group, because CPIC calls that
    allele normal function and the published cohorts were defined by the absence
    of an actionable variant, not of any variant.
12. **Does not display a recommendation without naming the drug it applies to.**
    A bare `AVOID` does not say what to avoid. `Drug` is a closed enum, and a
    phenotype/drug pair for which the engine returns no action is refused rather
    than rendered as an empty line.

---

## Running it against a site's data

`fixtures/cohort_template.csv` is the column layout to ask for. One row per
patient, no identifiers.

```
asl quadrant --csv their_data.csv
```

produces the four-quadrant table — panel status × toxicity — and the cell that
matters:

```
                         grade >=3    grade <3
panel POSITIVE                   2           0
panel NEGATIVE                   3           3

3/6 (50.0%) of panel-NEGATIVE patients had grade >=3 toxicity.
  Published comparators: 22.7% systemic (PMID 30348537),
                         13.6% incl. chemoradiation (PMID 30361102)
```

(That is the shipped 8-row template; two rows are excluded from the percentage
because they are panel-positive.)

**That is the question a retrospective exists to settle:** how many patients had
severe toxicity while testing negative on every panel variant? It is descriptive,
not a hypothesis test, and the output says so.

`asl.ledger.quadrant()` is stricter than the CLI: it **excludes** any patient
whose assay did not cover every panel locus, because such a patient cannot be
called panel-negative. Counting them as negative is the exact error this package
exists to prevent.

Both paths also **supersede duplicates**: a patient recorded twice is counted
once, using the later row, and the superseded count is printed. The log is
append-only by design, so a correction *is* a new row — before this was handled,
one patient appearing twice was counted as two, in two contradictory cells.

---

## Honest limitations

- **The Indian evidence base is thin and partly biased.** Two of the five
  published Indian DPYD studies tested *only* patients who had already developed
  grade 3+ toxicity (Patil 2016: 12 of 34; Sahu 2016: 28 of 506), which inflates
  any apparent genotype–toxicity association. Genotyped **without** toxicity
  selection: ~451 patients across two centres — Pavithran 2021 (n=375, Amrita
  Kochi, ASCO abstract only) and D-TORCH (n=76, AIIMS). **The methodologically
  stronger of the two found no association** (OR 0.71, p=0.612).
- **c.496A>G and `*6` are contested, not established.** Pavithran reports
  c.496A>G in 32 of 47 variant-positive patients with toxicity despite
  pre-emptive dose reduction; D-TORCH found it in 7/76 with no signal, in the
  same national population with a stronger method. `asl` reports both and claims
  no elevated risk.
- **The residual-risk numbers are not Indian.** Henricks and Lunenburg are Dutch
  cohorts. That is a limitation, and it is also why the argument holds without
  Indian data.
- **Plasma uracil phenotyping would be better in principle** — it measures enzyme
  activity and so captures causes a European-derived genotype panel misses. It is
  not proposed because a 955-patient, 17-centre Dutch study failed to validate the
  16 ng/mL threshold at all (de With 2022, PMID 35397172: OR 0.997, p=0.71), the
  assay degrades within ~47 minutes at room temperature (PMID 36104927), and
  nobody has characterised those kinetics at Indian ambient temperatures. **This
  remains a real gap.** CMC Vellore has a validated LC-MS/MS assay and published
  Indian data (DOI 10.2217/pme-2022-0042) and would be the site to test it.
- **No Indian mandate exists.** No CDSCO, ICMR, MoHFW or NCG document recommends
  pharmacogenomic testing before chemotherapy; DPYD testing costs ₹5,000–35,000
  and is not covered by PM-JAY or any state scheme. Adoption has to come from
  clinical usefulness, not policy.
- **Regulatory posture.** Under CDSCO's MDSW guidance (`CDSCO/MD/GD/MDSW/01/2026`,
  final 21 Jul 2026) this is *informs clinical management* — Class A, but still
  requiring registration under Chapter IIIB of MDR 2017. Not exempt. Under FDA's
  21st Century Cures §520(o)(1)(E) it likely qualifies as non-device CDS
  *because* it enables independent review of the basis, which is what the
  evidence-provenance display is for.
- **Not a diagnostic device.** Research tool. Informs clinical management; does
  not drive it.

---

## Status

| Phase | State |
|---|---|
| 0 · engine correctness (`pgx-core` 0.7.2 + 0.7.3) | **done**, tagged |
| 1 · generated evidence tables, CI drift check | **done** — 156 alleles across DPYD, NUDT15, TPMT |
| 2 · `assay` + `residual`, no-all-clear guarantee | **done** — 140 tests |
| 3 · `interpret` + `flags` + `report` + CLI | **done** |
| 4 · `ledger` + API + auth | ledger done + hardened; HTTP layer live in `anukriti-api` (authenticated) |
| 5 · `narrate` (the only LLM call) | not started |
| 6 · site retrospective | **the go/no-go** |

NUDT15/TPMT evidence tables are already generated for the second workflow
(CPIC thiopurine 2025 update, PMID 41618934), and `pgx-core` 0.8.0 adds
`NUDT15Caller`, but no NUDT15 interpretation path is wired here yet — and the
shape it should take is now narrower than first assumed.

**Why the DPYD argument does not transfer to NUDT15.** For DPYD, the activity
score is "validated in the Caucasian population, remains unverified in other
ethnic groups" (White 2021, PMID 34916829) — that published caveat is what the
residual-risk section rests on. CPIC says the *opposite* for thiopurines: "given
the robustness of the genotype–phenotype associations for both TPMT and NUDT15
across diverse populations… the recommendations are not limited to any specific
ancestry group." Nobody has shown CPIC wrong for NUDT15 in Indians, and the
variant that carries the Indian signal (`*3`, rs116855232) is the most-studied in
the gene.

What CPIC *does* publish for thiopurines is a single-gene coverage caveat with a
number on it: if only one of TPMT/NUDT15 is tested, "up to **10–15% of patients**
[have] actionable variants in the nontested gene" — which is the realistic Indian
scenario, since the validated low-cost ARMS-PCR assay reads NUDT15 c.415C>T and
TPMT`*3C` and nothing else. That, plus CPIC's own note that an untested rare
variant makes a patient "`*1` **by default**", is the honest basis for a
thiopurine workflow. Coverage, not reinterpretation.

## Where this sits in the platform

| Repo | Role |
|---|---|
| [`anukriti-pgx-core`](https://github.com/AnukritiAi-hq/anukriti-pgx-core) | the deterministic CPIC engine. Decides the phenotype. `asl` never overrides it. |
| **`asl`** (this repo) | what the phenotype leaves out: provenance, residual risk, outcome capture |
| [`anukriti-api`](https://github.com/AnukritiAi-hq/anukriti-api) | mounts `asl` at `POST /asl/dpyd/interpret`, on the **authenticated** surface because it takes a patient genotype |
| [`anukriti-main`](https://github.com/AnukritiAi-hq/anukriti-main) | frontend; `api.asl.*` client seam exists, no UI yet |
| [`anukriti_docs`](https://github.com/AnukritiAi-hq/anukriti_docs) | `ONCOLOGY_SOLUTION_AND_ARCHITECTURE_2026-08-16.md` — the full problem derivation, five research passes, and every defect round in this package |

`asl` is a sibling repo rather than a PyPI package, so `anukriti-api` imports it
via `PYTHONPATH` and degrades to `503` if it is absent.

## Citation

If this package or its argument informs published work, cite the analysis rather
than the code — the load-bearing claims are other people's:

- CPIC DPYD guideline: Amstutz U et al., *Clin Pharmacol Ther* 2018;103(2):210-216. PMID 29152729
- The Indian cohort: Baskarane, Divakar, … Batra, *Front Pharmacol* 2026;17:1732128. doi 10.3389/fphar.2026.1732128
- Residual-risk rates: Henricks LM et al., *Lancet Oncol* 2018. PMID 30348537; Lunenburg CATC et al., *Eur J Cancer* 2018. PMID 30361102
- Zero reference frequency ≠ absence: Chan SL, Zhang J, Pirmohamed M, *Br J Cancer* 2024;131:498-514. PMID 38886557

## Licence

Apache-2.0. See [`LICENSE`](LICENSE).
