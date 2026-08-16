# Security and safety policy

This repository has two kinds of risk, and the second one is the unusual one.

## 1. Reporting a vulnerability

Do not open a public issue. Email **security@anukriti.ai** with the affected
version, reproduction steps, and impact. Expect an acknowledgement within a few
working days.

## 2. Reporting a clinical-safety defect

**A wrong report is a security defect here, and it is the more likely one.**

`asl` exists because a pharmacogenomic result presented without its uncertainty
misleads clinicians. A bug in this package therefore has a route to patient harm
that does not involve an attacker at all — and the direction matters: a defect
that makes a result look *more* reassuring than the evidence supports is the
severe case.

Report these the same way, and please say plainly which of these it is:

| severity | example |
|---|---|
| **critical** | any output that reads as an all-clear; an untested locus presented as negative; a decreased/no-function carrier shown a rate or statement derived from wild-type patients |
| **high** | a dose or dose fraction that is not CPIC's, presented as if it were; a phenotype altered by anything other than the pinned engine; an activity value applied that CPIC has not published |
| **high** | PHI accepted by the ledger, or written to disk |
| **medium** | a citation that does not support the claim attached to it; a population frequency substituted from another population |
| **medium** | a refusal silently converted into a default |

Real examples that were found and fixed, to calibrate what counts:

- A `*2A/*2A` Poor Metabolizer was shown "22.7% of patients **without** an
  actionable DPYD variant experienced grade ≥3 toxicity." True sentence, wrong
  patient, and it erred toward reassurance for the patient with the most to lose.
- The report said "no variant was found among the loci tested" on a panel where
  three of four loci returned no call — one locus of reassurance rendered in the
  language of four.
- The PHI scan did not recurse into lists, so a phone number typed into a free-text
  toxicity description was accepted and written to the ledger.
- `ASL-POP-01` said a zero-frequency locus "cannot return a positive." That is an
  overclaim in the *opposite* direction — it would license dismissing a true
  positive as an assay error.

## Handling patient data

- **The library accepts no PHI.** The ledger enforces a field-name denylist plus
  value pattern checks, recursively, and refuses rather than sanitising.
- **The CLI is offline.** CPIC data is pinned in `data/evidence/`; nothing in
  `src/asl/` opens a socket. A clinic demonstration does not depend on, or leak
  to, the network.
- **Ledger files are gitignored** (`ledger*.jsonl`). If you ever find one in a
  commit, treat it as a data incident: rotate nothing, but purge the history and
  notify the site.
- **The HTTP surface is authenticated.** `anukriti-api` mounts `/asl/*` on the
  private surface precisely because it takes a patient genotype, and that is
  pinned by a test in that repo. Do not add an unauthenticated route to it.

## Scope limits, stated because they are safety-relevant

`asl` is a **research tool**. It is not a diagnostic device, has not been through
a site validation, and has no regulatory clearance. Under CDSCO's MDSW guidance
it is positioned as *informing* clinical management rather than driving it, which
is why it displays CPIC's recommendation as CPIC's and never generates a dose.

If you are considering clinical use: the go/no-go study described in the README's
*Status* table has not been done. Read *Honest limitations* first.

## Supported versions

Pre-1.0 and unreleased. Fixes land on `master`; there are no maintained release
branches yet. The `anukriti-pgx-core` pin is exact on purpose — a new engine
version carries CPIC table updates that require clinical regression review before
adoption.
