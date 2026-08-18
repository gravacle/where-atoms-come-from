# WHAT WE BUILT — WHERE ATOMS COME FROM — V001 — 2026-08-18

**The charter:** *find the process responsible for record formation and explain the roles of EM,
gravity, and alpha in that process.*

This document states what was built, what is measured, what is an identification rather than a
derivation, and what is not established. Every number is sealed in `REGISTER_V001.md` with its lane.

---

## 1. THE PROCESS

A **record** is a physical quantity that acquires a definite value and then keeps it. In this
construction it requires all of the following, each measured separately:

| | requirement | evidence |
|---|---|---|
| 1 | **`H` has an exactly degenerate eigenvalue.** A record must satisfy `[H,R]=0` and not be a function of the energy; if every eigenvalue is simple the commutant is abelian and no record exists | `dim commutant = Σ(multiplicity)²`, exact in every case; **0 of 400** random Hermitians are degenerate (W-60) |
| 2 | **The record is a closed loop.** Gauge invariance forces it, and reading it requires traversing the whole loop | T2, proved; cut-ring control **exactly `0.000000`** (W-37) |
| 3 | **It is selected, not chosen.** The dynamics picks which observable is the record | predictability sieve returns the full boundary at **28×**, nothing nominated (W-34) |
| 4 | **Durability requires `[H,R]=0` and `[L_k,R]=0`** | T1, proved for **any** GKSL generator with arbitrary jumps; and shown to need **no Markov assumption at all** (Phase A, Phase B) |
| 5 | **It becomes definite on every run**, invisibly to the ensemble | `E[⟨R⟩] = −0.015 ± 0.040` beside `E[⟨R⟩²] = 0.999`; 99.5% of runs definite (W-35) |
| 6 | **Capacity is exactly `m − 1`**; beyond it the carrier evicts one record rather than degrading all | T3, proved for arbitrary complexes, any prime field, any dimension; eviction in **0 of 793** baths protecting four (W-41, W-42) |
| 7 | **No partial legibility.** A region resolves the record or it does not | `0.000000` at 7 of 8 links, `1.000000` at 8 (W-57) |

---

## 2. THE THREE TERMS

### EM — the carrier and the record itself
Supplies the gauge field, the Gauss constraint, the boundary that is simultaneously a cycle and a
separator, and the holonomy that becomes the record. **The record observable is selected by the
dynamics, not nominated.**

### ALPHA — the dial, with three measured roles
- `1/(8g²)` — the number of times a record can be written and re-read before it decays (W-32)
- `9.6·g⁴·P/τ < 1` — the condition for a region to be read before it forgets (W-37b)
- **and the sharpest: alpha BREAKS the degeneracy records live in.** `H_magnetic` alone carries **65**
  records; **any** electric term drops it to **14**, immediately and not gradually (W-60)

### GRAVITY — supplying the exact degeneracy, via topology

**THE MEASUREMENT.** Same theory, same perturbation, two topologies:

```
                              unperturbed      splitting at eps = 1e-06
  TORUS   (genus 1)            4-fold                    4.9e-13
  DISK    (genus 0)            none                          --
  SYMMETRY-sourced (W-60)      4-fold                    2.0e-06   (LINEAR in eps)
```

> **Symmetry-sourced degeneracy splits linearly with the perturbation and dies. Topology-sourced
> degeneracy splits FOUR MILLION TIMES LESS at the same perturbation. The difference between the two
> carriers is GENUS — not dynamics, not symmetry, not a metric.**

**WHAT IS MEASURED:** the exact degeneracy that records require can be supplied by the **topology of
the space**, and when it is, it is protected in a way symmetry cannot match.

**WHAT IS AN IDENTIFICATION, NOT A DERIVATION:** calling that *gravity*. It satisfies the principal's
own criterion — *the same function in a different form* — and the function is real and load-bearing:
**supplying the exact degeneracy records live in.** But the topology here is the **lattice's**, not
spacetime's, and nothing in this program makes it respond to matter content. **Topological order is
established physics (Wen, Kitaev); what is ours is the chain W-60 → W-61 connecting it to the
record-formation requirement.**

**AND IT MEETS AN INDEPENDENT CONSTRAINT.** Harlow–Ooguri: quantum gravity admits **no exact global
symmetries** (proved in AdS/CFT). So a gravitating world **cannot** use global symmetry to protect
its records. **Topological degeneracy needs no symmetry — which is what such a world has left.** The
two results were obtained independently and point the same way.

---

## 3. THE OBSTRUCTIONS COHERE

Six obstructions reduce to four, and **three of the four are not dynamical** (W-59):

```
  T1  a conserved quantity is conserved          DYNAMICAL, near-tautological
  T2  gauge invariance forces closed loops       KINEMATIC        => implies W-57
  T3  capacity = m - 1                           COMBINATORIAL    => implies eviction
  --  spectral coincidence is non-generic        MEASURE-THEORETIC => W-55 and W-58 are ONE fact
```

W-57 follows from T2 (exact indicator, zero partial values). W-55's "no factorisation" and W-58's "no
global history" collapse together at the same coupling — kernel `4 → 1` while the sumset residual
jumps `3.0e-16 → 5.2e-02`.

---

## 4. WHAT IS NOT ESTABLISHED

- **Zero empirical contact.** Not one number came from a measurement.
- **A2 — background time.** 125 of 125 lanes evolve in an external `t`; W-58 built the first
  relational carrier, but the gravity-relevant case (a clock whose rate depends on what is there)
  collapsed the kernel. **Unresolved.**
- **A1 — factorisation.** Being repaired: the algebraic framework needs none, and given a region its
  boundary data is derived (W-56). But W-36, W-40 and W-43b still carry partial-trace numbers
  computed in a structure the system lacks.
- **Why the world's Hamiltonian has the topology it has.** W-61 moves the question from *"why these
  symmetries"* to *"why this topology"* — better posed, still unanswered.
- **The topology here does not respond to content.** ~~For the identification with gravity to be more
  than functional, it would have to.~~ **SECOND SENTENCE WITHDRAWN (`5c01e47`)** — it makes a classical
  measure the criterion (**D-1**). The live question is **G-6**: is `H₁` the *only* structure
  satisfying R1–R3?

---

## 5. THE ONE-LINE STATEMENT

> **A record is an exactly degenerate, gauge-invariant closed loop, selected by the dynamics, made
> definite on each run by monitoring, bounded in number by `area − 1`, and legible only in whole.
> EM supplies it, alpha breaks the degeneracy it needs, and the topology of the space supplies that
> degeneracy in the only form robust enough to survive perturbation.**
