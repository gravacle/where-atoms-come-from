# PLAN TO PROOF — WHERE ATOMS COME FROM — V001 — 2026-08-18

**THE PRINCIPAL, 2026-08-18:** *"this project isn't about building toy models. We want a proof."*

## 0. WHAT "PROOF" CAN MEAN HERE, STATED SO THE PLAN IS HONEST

Physics has no proof in the mathematical sense. It has two things, and a claim about a fundamental
process needs **both**:

1. **A THEOREM** — *given assumptions X, records must form this way.* Assumptions that are minimal
   and uncontroversial, and **not** "given a `Z_2` lattice."
2. **A PREDICTION THAT COULD HAVE FAILED** — a measurable consequence that does not follow from the
   existing account, which then survives testing.

**Everything this program has produced is measurement on one instance.** The plan converts that into
the two things above, and says where it will stop being able to.

## PHASE A — CONVERT MEASUREMENTS INTO THEOREMS (cheap, immediate, high yield)

**Several results were MEASURED that are actually PROVABLE.** Measuring a theorem on one lattice is
strictly weaker than proving it, and we did that four times.

| result | current status | what it should be |
|---|---|---|
| durability ⟹ unwritable | measured on 4 carriers × 2 groups | **general theorem.** `d<R>/dt = Tr(R·L[ρ])`; the coherent term gives `-i Tr([R,H]ρ) = 0`; the dissipator vanishes whenever `[L_k,R]=0`. **No lattice, no gauge group, no dimension enters.** A W-31 design agent noted it holds for NON-unitary jumps too — that generalisation is unverified and is part of this phase |
| reading requires a closed path | measured (cut-ring `0.000000`) | **standard gauge theory:** a gauge-invariant operator built from parallel transporters is supported on a closed path or terminates on charges. State and prove it |
| `capacity = cycle rank − 1` | proven for lattices (W-42) | **generalise:** the regions whose boundary avoids a given link form the kernel of ONE linear functional over GF(2). Holds for any complex, not just planar patches |
| `legibility = perimeter/2 − 1` | measured on `n×n` blocks (W-46) | prove from the rank of the boundary pairing; establish the class of regions it holds for |

**DELIVERABLE:** a theorems document with proofs and their exact hypotheses.
**WHAT IT BUYS:** the necessity results stop being facts about a `Z_2` lattice.
**FALSIFIER:** a proof attempt that fails reveals a hidden lattice assumption — which is itself the
most useful thing this phase can produce.

## PHASE B — UNIVERSALITY: WHAT SURVIVES CHANGING THE MODEL

Phase A covers what is provable. The rest must be **shown not to depend on the carrier.** Ordered by
how likely each is to break something:

1. **Non-abelian group.** `SU(2)`. Plaquettes no longer commute, so W-52's independence-at-`g²=0`
   and the sieve's operator ranking may both change.
2. **Non-Markovian bath.** Every result here assumes Lindblad. A structured bath with memory is
   where einselection stories most often fail.
3. **Dimension.** Everything is 2D planar. In 3D the boundary of a region is a SURFACE, storage and
   legibility scale differently, and W-46's exact pair must be recomputed.
4. **Continuum limit.** Whether any relation has a finite limit as the lattice spacing goes to zero.

**FALSIFIER, and it is the point:** if a result survives none of these, it was a lattice artifact and
must be withdrawn. **Expect withdrawals. That is the phase working.**

## PHASE C — EXTRACT A PREDICTION THAT COULD FAIL

**W-48 already killed the obvious candidate** — the decoherence half of the composed law reduces to
`exp(-ΓT)`, so taking it to interferometry data would fit a textbook result.

**The surviving candidate is the storage/legibility split:** storage scales with the bulk, external
legibility with the boundary, as an exact pair on the same object. **Standard decoherence says
nothing about it.** Its natural setting is engineered quantum memories, where the number of storable
states and the number of externally readable ones are both controlled quantities.

**THE OBLIGATION BEFORE TOUCHING ANY DATA:** establish that the split predicts something the existing
account of those systems does not. **If it does not, say so and stop** — that is the same discipline
that produced W-48, and it is what prevents fitting.

## PHASE D — CONTACT

Only reachable if C produces a prediction that survives. **Not before.**

---

## HONEST FEASIBILITY

**Phase A is doable now** and will make the program's four necessity results considerably stronger
than they are.
**Phase B is doable and expensive**, and its most likely product is withdrawals.
**Phase C is speculative.** It may terminate in "no distinguishing prediction exists," which would
mean the account is a re-derivation of known decoherence physics in gauge language.
**Phase D is contingent on C.**

> **The plan does not promise a proof. It promises to find out whether one is available, and to say
> so either way.**
