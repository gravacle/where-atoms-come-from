# REPLICATE — how to check this program yourself

**Everything here is re-runnable. Nothing asks you to take a number on trust.**

```bash
git clone https://github.com/gravacle/where-atoms-come-from && cd where-atoms-come-from
./replicate/reproduce.sh --seals     # fast: verify every SHA-256 manifest
./replicate/reproduce.sh             # full: re-run every script, diff against sealed output
```

Requires **Python 3.9+** and **numpy** only. No network, no GPU, no data download. The heaviest
script needs a few minutes; `--quick` skips those.

`reproduce.sh` **exits non-zero** if any seal fails, any script's output differs from the sealed
`.txt` beside it, or the status grid does not regenerate byte-for-byte from its ledger.

---

## WHAT IS CLAIMED

Read in this order. Each is short and each links to the numbers behind it.

| | document | what it is |
|---|---|---|
| 1 | `CORE_FRAMEWORK_V001.md` | **the anchor** — the five clauses defining a record |
| 2 | `THIRD_TERM_V001.md` | Γ, defined **only by the function it performs** (R1–R3) |
| 2b | `GLOSSARY.md` | **every term, and whose it is** — BORROWED, RENAMING, or OURS. Read it before believing any claim is a discovery |
| 3 | `MODEL.md` | **the model** — records constructed from `(H,{L_k})` and nothing else |
| 4 | `STATUS_LEDGER_V001.md` | **every claim, with its status and its evidence** — generated, never typed |
| 4b | `THE_PLAN_V001.md` | **the road to a proven process** — every task with a checkable DONE WHEN. `./ledger/plan.py` |
| 5 | `PROCESS_V002.md` | **the process as it now stands** — every sentence cites a ledger row |
| 5b | `RECORD_FORMATION_V001.md` | what a formation process must deliver |
| 6 | `REGISTER_V001.md` | the append-only record, including **everything withdrawn** |

**Start with the ledger.** It is the whole program in one grid: 90+ rows, each with a status from a
**closed vocabulary**, the evidence, what blocks it, and **what would refute it**.

---

## HOW TO CHECK ONE CLAIM

Every row names its lane. Example — the record-count law, row **C-14**:

```bash
cd model && python3 count_law.py
```

You should see **22 PASS, 0 FAIL**, and a control showing the naive reading fails on `[3,3]`,
`[6,6]`, `[5,5]`. If you get anything else, the claim is wrong and we want to know.

**Is any of it an artefact of one carrier?** That is the question `T-9` exists for:

```bash
cd LANE_T9_CARRIERINDEP && python3 t9_sweep.py
```

**32 PASS, 0 FAIL** across `[[8,2,2]]` toric, `[[8,1,2]]` **non-manifold**, and `[[4,2,2]]` which is
not a lattice at all. **A row that survives only one carrier is NOT marked so in the ledger — this sentence was false when written and is corrected here, 2026-08-20, by external review.** The string `SINGLE-CARRIER` appears nowhere in `ledger/status_ledger.tsv`; roughly 24 of 162 rows cite two structurally different carriers. Plan task T-9 was marked DONE while satisfying neither branch of its own DONE_WHEN and is REOPENED.

---

## HOW TO READ A NUMBER HONESTLY

This program has caught itself reporting artefacts as results more than once. Three rules were
adopted because of it, and every script follows them:

- **A zero is only a measurement if a POSITIVE CONTROL beside it would have registered a non-zero.**
  Look for the control before believing any zero.
- **A check that cannot fail is not a check** (`D-8`). `check(..., True)` is decoration; a fit with no
  noise floor beside it is decoration.
- **The ensemble average `Tr(Rρ)` is blind to record formation** (`D-6`). A fair coin has mean zero
  and every flip is definite. Negative verdicts warranted only by an expectation value are void.

---

## WHAT IS **NOT** ESTABLISHED

Stated here so you do not have to find it: **empirical contact is zero.** No number in this program
came from a measurement, and **no distinguishing prediction has been identified** (`X-4`, `T-VI.3`,
`T-VI.4`, all BLOCKED). The results are exact statements about a defined object, verified by
computation. **They are not yet physics in the testable sense**, and the ledger says so on its face.

The road from here is the **`PF` series** in the ledger — the eight steps to a full proof, each with
a status you can check.
