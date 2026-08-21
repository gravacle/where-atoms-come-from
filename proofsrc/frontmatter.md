# PROOF — WHERE ATOMS COME FROM — V002 — 2026-08-20

> **"We need a model to work on that represents the full project."** — the principal, 2026-08-20
>
> **"The model should be the overall representation of the proof."** — the principal, 2026-08-20

This document is the narration of the **UNIVERSAL RECORD MODEL** — `model/project_model.py`, the URM.
Every section below is one of the model's layers. Every claim is one model function, its validator
gate, its grounding, and its ledger row. **There is no step outside the model.** A statement the URM
does not carry is not in this proof.

**`PROOF_V001.md` is not superseded as mathematics.** Its Theorems A–D stand unchanged; what this
document changes is their siting. They are statements about the **DEF-A corner** — §4 — and they are
`FORMAL`. The program's proof is no longer those theorems; it is the model, and this is the model's
narration.

---

## 0.1 WHAT `PROOF` MEANS IN THIS DOCUMENT

> **The principal, 2026-08-20:** *"The PROOF must be something that any physicist anywhere in the
> world can run against their real world data and confirm that it works as asserted. That is a proof.
> Otherwise we're not proving anything other than that we can construct equations."*

A row is `PROVED` only when it names **(1)** the measurable quantity in units, **(2)** the predicted
value or relation, **(3)** what would falsify it, and **(4)** at least **two structurally different
real record surfaces** — different mechanism, never different parameters.

**Only requirement (4) is enforced by a tool.** `ledger/status.py` refuses the status unless the row's
grounding carries a `RECORDS VERIFIED:` line naming at least two distinct surfaces. Requirements
(1)–(3) are enforced by the registrar reading the row, and by nothing else. That asymmetry is stated
here because a bar half-enforced by a tool and half by a person is exactly the shape every guard in
this program had before it started refusing.

**One row in this program holds that status: `C-71`, and §3 states it in full.** Its measurement is
a *within-part comparison* — the shift of a programmed page's `V_t` distribution against the same
part's own erased population — so it needs no absolute baseline, and that is why it survived the
erratum below.

**`C-72` was `PROVED` and is now `PARTIAL`, on this document's own finding.** Writing §3 surfaced
that its registered prediction, *"the ratio is 1 for every data pattern"*, was an **algebraic identity
of scoring an erased cell as exactly zero charge** — the sealed lane's own source says so — while the
same model layer gives those cells a `±5 e` residual, under which the ratio is `0.968`–`0.980` and
depends on the pattern. That is `C-72`'s **own registered falsifier**, fired by its own model. A
replacement protocol was drafted and then **refuted by two independent refuters**: the ratio is not
translation-free, so it needs `V_t,neutral` to better than the tolerance it is bounding, which is not
a datasheet number. **What survives is computed and confirmed** — a data-independent floor with no `N`
in it, against an orientation surface that screens as `N^{-1/2}`, so the two encodings differ in kind
and their separation grows as `√N`. What does not survive is a runnable prediction. The repair is
`T-50`.

**`C-71`'s surface is MODELLED** — constants pinned to literature classes, patterns drawn from sealed
seeds. `RECORDS VERIFIED` names two structurally different real *mechanisms* whose standards were
scored inside the model; **no device was measured.** What is offered to the outside is the falsifier.

Every other claim here is `FORMAL`, `DEFINED`, `CANDIDATE`, `PARTIAL` or `OPEN`, and each block prints
its own. **A `FORMAL` result is real mathematics about the program's own stipulated definition and is
not a claim about the world** (`H-3`, `PARTIAL`, standing).

**This document does not claim that gravity has been derived.** §5 states what has been computed;
§6 makes the one comparison the program permits, with the conditions that are not earned named in
the same sentences as the result.

---

## 0.2 HOW TO READ A CLAIM BLOCK

Every claim is one block. The five cells are the claim's whole warrant:

| cell | what it is |
|---|---|
| **model** | the URM function that carries the claim. `none` means narration, and the scope cell says why |
| **gate** | the validator check that fires on it, as `file :: check name`. `none` means no gate exists yet |
| **grounding** | where the numbers come from: a sealed lane, a `D-25` provenance entry, or a pinned external source |
| **rows** | the rows it rests on, each with its **status** and its **carrier mark**. A row is a claim in `ledger/status_ledger.tsv` or a task in `ledger/plan.tsv`; a task prints its plan status and the tier `PLAN`, and is never carrier evidence |
| **scope** | the caveat that travels with the claim wherever it is quoted — what the claim does **not** say |

The carrier mark is the `T-9` audit's verdict (`LANE_T9_AUDIT/T9_carrier_audit.tsv`): `TWO-CARRIER`
means the result stands on two structurally different carriers; `SINGLE-CARRIER` means one;
`NOT-CARRIER-SHAPED` means the result is not the kind of thing a carrier carries.

**`UNAUDITED` means the row carries no mark at all**, and it is not a weaker `SINGLE-CARRIER`. The
audit's original scope was the then-live `FORMAL`, `PROVED` and `MEASURED` rows, and its T-52
extension brought every later in-scope row through the same audit. `DEFINED`, `PARTIAL`, `CANDIDATE`,
`OPEN` and `BLOCKED` rows remain unaudited unless an extension explicitly reaches them. **A block
resting on no `TWO-CARRIER` row opens its scope
cell with `SINGLE-CARRIER —`**, and an unaudited row never lifts that requirement.

**Of the 150 audited rows, 113 are `SINGLE-CARRIER`, 22 are `TWO-CARRIER`, 15 are
`NOT-CARRIER-SHAPED`.** That is the program's state, printed rather than described.

---

## 0.3 HOW TO CHECK THIS DOCUMENT

Each line runs from the repository root:

```bash
python3 replicate/check_proof.py      # this document's own gate — expect GATE PASSED
python3 model/validate_project.py     # base project + D-25 gates — expect 24 PASS, 0 FAIL
python3 model/validate_geometry.py    # geometry layer, then project chain — expect 33 + 24 PASS
python3 model/validate_urm.py         # four folded families, then both chains — expect 176 family / 233 full PASS
./replicate/reproduce.sh              # re-runs every lane script against its sealed output
```

**The fifth does not pass, and this document will not pretend otherwise.** The most recent completed
run in the repository records **46 differing scripts against 267 identical**
(`replicate/reproduce_full2.log`): thirteen nondeterministic scripts, unnormalised wall-clock lines,
one script that aborts, and stale seals. A fresh full run is in flight. **The debt is `T-35`, `TODO`,
and §8 states it.** The first four commands pass as written. The validator stack contains **24 base
project/D-25 gates, 33 geometry gates, and 176 folded-family gates: 233 model gates in the full
umbrella run**; the proof gate is separate.

`check_proof.py` **refuses** a block with a missing field, a row that is in neither the ledger nor the
plan, a row that is `WITHDRAWN` or `FAILED`, a status or carrier mark that has gone stale against the
record, an unmarked single-carrier block, a model function or gate name that does not exist, the word
`PROVED` without a `PROVED` row behind it, or a banned classical-gravity comparison outside §6. It
exits `2` rather than passing when it cannot read the ledger or the audit.

**What it does not reach is section prose.** R1–R9 and R11 parse `### P-` blocks only and close at the
next heading; only R10's `D-1` scan reads the whole file. **Every section header in this document,
§0 included, cites no ledger row and is enforced by nothing but that scan** — which is why the headers
below carry no claim that is not also in a block. Closing that gap is the gate's own standing debt.

**What this document does not have, stated up front.** {{STATS}} The blocks without a gate are
narration over sealed lane output, and each one says so in its scope cell. Closing that gap — every
claim in this proof gated by a check in `model/` — is the proof's own standing debt.
