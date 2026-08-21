# INTEGRATION_arrow — folding the ARROW layer into the URM (T-54/T-55)

Builder: the arrow-family lane, 2026-08-21, under the principal's directive that the URM
is the world model — the framework new observations are added INTO.  This document is the
registrar's integration spec.  The builder wrote ONLY:

- `model/arrow.py` — the ported/delegating machinery (new file)
- `model/checks_arrow.py` — the gate block, initially 26 checks; 27 after the
  adversarial custom-model refusal repair, all PASS (new file)
- this file

Nothing else was touched: not `project_model.py`, not any validator, lane, ledger,
register, or document.

## 1. The layer

The family PROOF_V002 §3 says is owed a LAYER: the arrow results F-17 (threshold at
weight d), F-18 (entanglement without record bits), F-19 (irreversibility from inside),
F-21 (fragment redundancy carries the threshold), with the PF-2 history that F-20 rests
beside.  Sealed sources: `LANE_F1_ARROW` (f1_arrow.txt, f1b_invariance.txt),
`LANE_PF2_DYNAMICAL` (pf2_history.txt), the T-9 battery rows for F-20/F-21
(`LANE_T9_AUDIT/T9_carrier_audit.tsv`).

Existing machinery reused, not duplicated: the sealed PF-2 bath IS
`record_model.Environment`'s default; `RecordModel.formation` was already proved
sealed-equal by `validate_formation.py`; **`RecordModel.redundancy` — called by no
validator before this fold-in — is now wired** (checks "F-21 RecordModel.redundancy
(WIRED)" and the weight-1 zero row both call it).

## 2. ProjectModel methods to add (the ARROW section)

Append to `ProjectModel`, after the GEOMETRY section, following its delegate idiom
(`import arrow as AW` inside each method).  Signatures and docstring text verbatim:

```python
# ---------------------------------------------------------------- ARROW (T-54/T-55; F-17..F-21)
# Every method below delegates to model/arrow.py -- machinery ported from the sealed
# arrow lanes or delegated to record_model, per that module's header.
def arrow_threshold(self, lam=0.8, weights=(1, 2), coupling=None):
    """F-17: THE ARROW CARRIES THE RECORD'S OWN THRESHOLD.  chi(O:B) swept over every
       observable of each listed weight in the mean-force state under the weight-d
       coupling -- sealed: 0.00000000 on all 24 weight-1 observables, 0.11448276 at
       weight 2 = d, closed form from Z_B(+-1) exact (LANE_F1_ARROW part 4, sealed).
       Owner: ORIGINAL.  Scope: SINGLE-CARRIER, toric-2x2 (sealed T-9 audit)."""
    import arrow as AW
    return AW.arrow_threshold(lam=lam, weights=weights, coupling=coupling)

def arrow_ledger(self, lam=0.8):
    """F-18: THE RECORD'S ARROW IS NOT AMBIENT DECOHERENCE.  The four-row coupling
       ledger (weight, I(S:B), chi(record:B)) -- sealed: the weight-1 coupling entangles
       I(S:B) = 0.04549256 yet transfers ZERO record bits (LANE_F1_ARROW/f1b part (c),
       sealed).  Owner: ORIGINAL.  Scope: SINGLE-CARRIER."""
    import arrow as AW
    return AW.arrow_ledger(lam=lam)

def arrow_invariance(self, n_unitaries=12, seed=5, lam=0.8):
    """F-19: IRREVERSIBILITY FROM INSIDE.  System-only unitaries cannot reduce I(S:B)
       (sealed 3.686e-14 over 12 unitaries; covariance check passes; chi about the FIXED
       label moves 1.145e-01) -- the copy is relocatable, never erasable, from inside
       (LANE_F1_ARROW/f1b parts (a,b), sealed).  Owner: BORROWED (textbook invariance);
       the relative-arrow reading for records ours (T-III.6)."""
    import arrow as AW
    return AW.arrow_invariance(n_unitaries=n_unitaries, seed=seed, lam=lam)

def arrow_history(self, times, coupling=None, lam=0.8, env=None, keep_states=()):
    """The arrow as a HISTORY (F-20-adjacent): chi(record:B)(t) from a product state
       with <record> exactly conserved -- the record is READ, not written; negative
       times are the sealed reversal control (LANE_PF2_DYNAMICAL, sealed; one
       eigendecomposition serves every time).  RecordModel.formation is the single-time
       twin (validate_formation.py).  Owner: ASSEMBLED.  F-20's mechanism TWO-CARRIER
       (toric-2x2; bouquet); this venue toric."""
    import arrow as AW
    return AW.arrow_history(times, coupling=coupling, lam=lam, env=env,
                            keep_states=keep_states)

def arrow_redundancy(self, coupling=None, lam=0.8, t=4.0, env=None):
    """F-21: REDUNDANCY CARRIES THE RECORD'S THRESHOLD.  Whole-bath and per-fragment
       chi through RecordModel.redundancy -- sealed: fragments 0.789366/0.048377/
       0.678602 under the weight-d coupling, EXACTLY ZERO under weight-1
       (LANE_PF2_DYNAMICAL parts 3-4, sealed).  Owner: ASSEMBLED (quantum-Darwinism-
       style apparatus, Zurek/Blume-Kohout; the threshold-in-fragments finding ours).
       Scope: fragment bits SINGLE-CARRIER toric-only -- the T-9 battery replicated
       only the whole-bath weight-1 null on [[8,1,2]]/[[4,2,2]]."""
    import arrow as AW
    return AW.arrow_redundancy(coupling=coupling, lam=lam, t=t, env=env)

def arrow_observation(self, env, coupling, record=None, model=None, lam=0.8, t=4.0,
                      tier="world", provenance=None):
    """OBSERVATION ENTRY for the arrow family (T-54/T-55): score a NEW bath/fragment
       observation through the family's own instruments -- I(S:B), whole-bath chi,
       per-fragment chi, and the verdicts (holds_record_bits, entangled_without_record
       = the F-18 class, redundant_fragments).  A custom RecordModel must bring its own
       explicit record; the toric default is used only with the default model.  D-25 AT
       THE GATE: world-tier baths
       require provenance; corner baths must self-declare 'DEF-A'.  Every outcome
       registers -- entangled-without-record is a RESULT, not a failure."""
    import arrow as AW
    return AW.score_bath_observation(env, coupling, record=record, model=model, lam=lam, t=t,
                                     tier=tier, provenance=provenance)
```

## 3. Where the checks chain in

`model/checks_arrow.py` exposes `run_arrow_checks(check)` in the `validate_geometry.py`
idiom (`check(name, cond, detail="")`).  The builder considered two integration options:

- **Recommended: a sibling `model/validate_arrow.py`** mirroring `validate_geometry.py`'s
  skeleton — its own `check()` counter, `run_arrow_checks(check)`, then the CHAIN section
  running `validate_project.py`, exit `0` iff both.  This keeps `validate_geometry.py`'s
  verified 31-gate count untouched (the T-46 verifier counted it) and keeps the ~108s
  arrow runtime out of the geometry gate.
- Or append `from checks_arrow import run_arrow_checks; run_arrow_checks(check)` inside
  `validate_geometry.py` before its summary — one validator, one conjunction, but the
  gate count and runtime both change.

**Registrar disposition:** the shared `model/validate_urm.py` umbrella was selected. It
runs ARROW as a separately counted 27-gate family, then the other three homed families,
geometry, and project/D-25. The added gate is the custom-model/default-record refusal
with an explicit matching-record positive control. Standalone remains available at
27 PASS, 0 FAIL.

## 4. The observation-entry story

How a NEW observation of this family's kind — a bath, or fragments of one, coming to
hold (or not hold) a record — enters the URM:

1. **The bath spec enters through the gate.** The observer supplies a
   `record_model.Environment` (qubit count, energies, beta — the observed environment)
   and the coupling the system meets it through.  `arrow_observation` REFUSES a
   world-tier bath without provenance (the real environment it models, constants'
   pinned sources) and a corner bath without the `DEF-A` self-declaration — D-25
   mirrored from `URM.surface`, so the exact idealisation can never silently pose as
   the world at this layer either.
2. **The family's own instruments score it.**  One evolution
   (`RecordModel.evolve`, shared eigendecomposition), then: I(S:B) — does it entangle
   at all; whole-bath chi — does it hold record bits (the F-17/F-18 discriminator);
   per-fragment chi and the redundancy count (F-21's instrument).
3. **Every outcome registers.**  `holds_record_bits` with redundant fragments is the
   F-21 class; `entangled_without_record` is the F-18 class (a result, not a failure);
   all-zeros with the identity-class coupling is the ledger's bottom row.  The check
   block demonstrates the whole path on a never-sealed 2-qubit bath: the gate refuses
   the undeclared entry, then scores the declared one, and the F-18 discriminator fires
   on the new bath — the threshold is structural, not a memorized state.
4. **A new LAW of this family** (say, a redundancy scaling law) would enter as a further
   layer method whose validator gates it the same way: sealed anchors + computed
   comparisons + probes beyond the gated range.

## 5. Numbers the gates hold (all reproduced, sealed sources named in the check block)

| Row | Sealed value | Reproduced |
|---|---|---|
| F-17 | weight-1: 0.00000000 (all 24); weight-2: 0.11448276 (252 swept, argmax ZZ[0,2]); closed form exact | exact at print precision (measured weight-1 max 3.9e-15) |
| F-18 | Ze row: I(S:B)=0.04549256, chi=0.00000000; Zbar 0.11448276/0.11448276; Zbar2 0.11448276/0; identity 0/0 | exact at print precision |
| F-19 | I(S:B) invariance 3.686e-14 (12 unitaries, seed 5); covariance PASS (<1e-8); fixed-label movement 1.145e-01 | 3.686e-14 and 1.145e-01 exact at print precision; covariance 3.9e-15 (blockwise route; lane 9.992e-16 — same statement, the lane's own 1e-8 bound is the gate) |
| PF-2 | chi(t)=0 / 0.40660635 / 0.81447230 / 0.97527192 / 0.78665760 / 0.90811968; <Zbar> exactly constant; reversal exact; I(S:B)=chi | all exact at print precision |
| F-21 | whole 0.90811968; fragments 0.789366/0.048377/0.678602; weight-1 all EXACTLY ZERO | exact (through RecordModel.redundancy — the wire) |

Probes beyond every gated range: lam=0.5 (chi == closed form at a never-sealed value);
t=±3.0 (reversal equality off the sealed grid); fragment pair {0,2} = 0.90454725
(data processing, never sealed); the never-sealed 2-qubit bath through the entry gate.

## 6. What the registrar should also know

- **Packaging economies, both identical maps** (module header, "PACKAGING NOTES"): the
  conditioning projector applied blockwise instead of via explicit `kron`; the history's
  joint Hamiltonian diagonalised on the real path when its imaginary part is exactly
  zero.  The F-21 cross-instrument gate (complex path in `RecordModel.redundancy` vs
  real path in `arrow_history`, agreement 4.8e-15) polices the second.
- The one sealed number visibly moved by the first economy is the machine-epsilon
  covariance bound (3.9e-15 vs 9.992e-16) — gated at the lane's own 1e-8 bound with the
  sealed print quoted.
- Scope honesty is in the docstrings: F-17/F-18/F-19/F-21 SINGLE-CARRIER per the sealed
  T-9 audit; F-21's fragment bits toric-only; F-20 TWO-CARRIER (toric-2x2; bouquet).
  Nothing in this fold-in upgrades a scope.
- Runtime: ~108s total, dominated by three unavoidable dim-2048 eigendecompositions
  (the history's shared one, and one inside each of the two `RecordModel.redundancy`
  calls — weight-d and weight-1 rows).
- `arrow.py` is importable standalone (verified from a foreign cwd), does no heavy work
  at import (carrier/model/env are lazy singletons), and returns DATA.
