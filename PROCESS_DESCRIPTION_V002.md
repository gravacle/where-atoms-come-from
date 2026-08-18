# THE PROCESS — DESCRIPTION AND TERMS — V002 — 2026-08-17

**PROPOSED, NOT ADOPTED.** Supersedes `PROCESS_DESCRIPTION_V001.md`, which is wrong in two places
and is kept for the record. Everything below is measured and sealed in `REGISTER_V001.md`.

---

## 1. THE PROCESS

**A record forms when five conditions hold at once. Each was measured separately; none was assumed.**

| # | condition | measured | lane |
|---|---|---|---|
| 1 | **A carrier with a boundary that is both a CYCLE and a SEPARATOR** | perimeter degrees all 2; removal → 5 components; discrete Stokes exact `0.00e+00` | W-27 |
| 2 | **A record observable `R` commuting with `H`** — else `d<R>/dt = 0` fails and the bath erases it | 4 carriers × 2 groups | W-30a |
| 3 | **An environment that CANNOT DISTURB `R` but is CORRELATED with it** — `[L,R] = 0`. This is a quantum non-demolition measurement | `Γ ~ g⁴` (shielded) vs `Γ` flat (unshielded) | W-32 |
| 4 | **Continuous monitoring**, which collapses each individual run to a definite value | `E[<R>²]: 0 → 0.945 → 0.999`; 99.5% of runs definite | W-35 |
| 5 | **Redundant copying**, so that many observers agree — *the one condition still missing* | `I(|F|=1)/I(all) = 0.047`: no redundancy | W-36 |

**Conditions 1–4 are met. Condition 5 is not, and the reason is precise (§4).**

## 2. THE DYNAMICS

**CARRIER.** A finite complex — vertices, links, plaquettes — with the field on links.
**CONSTRAINT.** A Gauss law at every vertex, imposed exactly, with **no coupling constant**.
**HAMILTONIAN.** `H = -(1/g²) Σ_plaquettes W - g² Σ_links Z`, magnetic and electric.
**ENVIRONMENT.** Lindblad coupling at rate `γ`. **Required** — a closed pure system has `S = 0`
forever and can pay for nothing.
**RECORD.** The rim Wilson loop. **Not nominated — SELECTED**: the predictability sieve, run over
the whole operator algebra with nothing named in advance, returns the full boundary and only the
full boundary, **28× slower than any other loop**, with single, double and triple loops all clustered
together. Move the bath and the pointer basis moves with it (W-34).

## 3. THE THREE TERMS

**EM — PRESENT, AND IT IS BOTH HALVES PLUS THE CONSTRAINT.** Electric flux on a CUT, magnetic flux
on a CYCLE, Gauss at the vertices: Maxwell's equations entire. EM supplies the carrier's structure
and the record-bearing observable itself.

**ALPHA — PRESENT, AND IT SETS HOW LONG A RECORD LASTS.** `ω = 4g²` is the rate the record can be
driven (first order); `Γ = 32g⁴` is the rate it decays (second order); so

> **`RATIO = 1/(8g²)` = the number of times a record can be written and re-read before it decays.**
> Slope `-1.000` across two decades. At `g² = 0.0005` that is **250 uses**, doubling every time the
> coupling halves. A pure number, read off the Lindbladian spectrum, not fitted. (W-32)

**GRAVITY — STILL NO TERM, AND V001's CLAIM ABOUT IT IS WITHDRAWN.** V001 said "without a responsive
carrier there are no records at all." **That is false.** Approximate records form and persist on a
static carrier, and approximate is what every real record is. What remains true is narrower: the
carrier is static, the gauge group here rotates the field while GR's constraints deform the geometry,
and **no measurement in this program yet requires the difference.**

## 4. THE OPEN OBSTRUCTION — LEGIBLE vs PROTECTED

**`[L,R] = 0` is what protects the record — and it is also what stops any local fragment of the
environment from carrying a copy.**

```
  env coupled to Z on the cut  (LOCAL)     I(|F|=1)/I(all) = 0.047   one delocalised copy
  env coupled to R itself   (NONLOCAL)     I(|F|=1) = 0.999 bits, FLAT   perfect redundancy
```

> **THE OBSERVABLE THAT SURVIVES IS THE ONE NOBODY CAN READ.** `R` is a product over 8 rim links —
> nonlocal — and a local environment cannot couple to a nonlocal observable.

**This is the same shape as W-30's obstruction, one level up: there, writable vs durable; here,
legible vs protected.**

## 5. WHAT THE INSTRUMENT MISSED, AND WHY IT MATTERS

**Every lane from W-28 to W-34 measured `<R> = Tr(R ρ)` — an ensemble average.** A fair coin has mean
zero and every flip is still definite. `E[<R>] = -0.015 ± 0.040` while `E[<R>²] = 0.945`: **the tie
is broken on every single run, and the value is a fair coin.** W-29c's "no bath selects a sector"
and W-33's "preserves but never amplifies" were both measuring the coin's fairness. **Selection was
never missing. The instrument could not see it.** (W-35)

**And record formation is not a property of the master equation at all** — the jump unravelling of
the identical generator gives `E[<R>²] = 0.000000`, with the same `ρ(t)`. **It is a fact about what
the environment does with what it learns.** That is why §1 condition 5 is not optional bookkeeping:
redundancy is what makes the answer observer-independent.

## 6. THE NEXT STEP

**In gauge theory a nonlocal observable is measured by transporting a probe charge around the loop —
Aharonov–Bohm. Local at every instant; the nonlocality is in the PATH, paid for in TIME proportional
to the loop's perimeter.** So a local environment *can* read the rim loop, and case B above is the
effective description of that.

**BUILD: the probe transported around the rim.** Then the closing question is a competition of
timescales — **readout time grows with the perimeter, while the record decays at `32g⁴`** — and
whether a region can be read before it forgets is a relation between the coupling and the size of
the region. **That relation, if it exists, is the next pure number in this program.**
