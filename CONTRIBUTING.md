# Contributing to `asl`

Before anything else, read [`README.md`](README.md) — particularly *What it
refuses to do*. Those twelve refusals are the package. Code that removes one is
not an optimisation, and it will fail a test.

## The one rule everything else follows from

**The deterministic layer decides; this package explains.**

`anukriti-pgx-core` owns the phenotype call. `asl` calls it, treats the answer as
final, and attaches provenance and residual risk *around* it. It does not
re-derive a phenotype, adjust an activity score, or reclassify an allele.

This is not stylistic. The reason it is a hard rule is a real incident:
`pgx-core` once carried DPYD `*6` as *Decreased function* against CPIC's *Normal
function*. The clinical evidence supported our value — and it was still wrong to
ship, because it was an **undeclared deviation from the pinned source of truth**,
and it survived four releases undetected. If you believe CPIC is wrong, the
mechanism is a named conflict flag *alongside* CPIC's call, never instead of it.

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest            # 114 tests, all offline
ruff check src tests
```

No network access is needed to run the suite, and none should be introduced into
the library. CPIC's data is pinned in `data/evidence/`; only
`scripts/regen_evidence_tables.py` talks to the network.

## The evidence layer is generated, never hand-authored

`data/evidence/*_evidence.json` is produced by
`scripts/regen_evidence_tables.py` from CPIC's API. Do not hand-edit those files.
CI runs the script with `--check` and fails the build if the committed tables
have drifted from CPIC.

Exactly four things are hand-maintained, all in `*_conflicts.json`, because CPIC
does not publish them:

| field | why |
|---|---|
| `cohort_n`, `cohort_ancestry` | CPIC's citations don't carry them |
| conflict flags (`P2_*`, `P3_*`) | CPIC has no field for "the literature disagrees" |
| `pending_guideline_changes` (`P4_*`) | CPIC announces changes before publishing them, and freezes its database meanwhile — so no automated check can see them |

**None of those four may feed a phenotype.** Asserted by test.

## Testing

Every refusal in the README corresponds to a test in
`tests/test_safety_properties.py`. If you add a refusal, add the test; if you
relax one, the test must be updated deliberately and the reason written down.

Two conventions that are worth more than they look:

- **A new test must fail against the pre-fix code.** Verify it before you commit.
  A test that passes both ways documents nothing.
- **Drive patient-shaped inputs through the pipeline; don't only read the code.**
  Every defect round in this package's history was found this way, and the
  blind spots have been systematic:

  | round | what was varied | what was missed |
  |---|---|---|
  | 1 | input shape | the carrier — every case was panel-negative |
  | 2 | patient type | the shipped CSV template, which refused all 8 rows |
  | 3 | the entry point | *(open — the ledger and NUDT15 have had no adversarial pass)* |

  An adversarial pass inherits the blind spot of whatever example it is built
  around. Pick an example the last one didn't use.

## Style

- Black-compatible formatting, `ruff` clean, line length 88.
- Type hints on public functions; `from __future__ import annotations` at the top.
- Docstrings explain *why*, and especially why a refusal exists. A future
  contributor deleting a check they don't understand is the realistic failure
  mode, so the comment is the guard.

## Commits

Conventional Commits: `feat(scope):`, `fix:`, `test:`, `docs:`, `data(scope):`.
Imperative subject under ~70 characters; the body explains why, with concrete
sub-changes as a bullet list.

**No commit trailers.** No `Co-authored-by`, no `Generated-by`.

## What needs doing

See the *Status* table in the README. The honest short version: the API and auth
layer, the `narrate` step, and — the actual go/no-go — a site retrospective that
can arbitrate whether the contested alleles carry risk in this population. No
amount of code substitutes for that.
