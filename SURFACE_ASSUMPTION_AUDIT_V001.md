# SURFACE ASSUMPTION AUDIT — V001 — 2026-08-18

**The principal:** *"Are we sure that we don't have any incorrect assumptions imported into the
surface that are a root cause of not being able to find a process from here?"*

**The audit was never done.** Every assumption in the five-step process, classified by whether it was
derived, imported and tested, imported and untested, or **imported after being measured to fail**.

## THE FIVE STEPS AND WHAT EACH ASSUMES

| # | assumption | where | status |
|---|---|---|---|
| **A1** | the Hilbert space **factorises** into system ⊗ environment | step 1, and every partial trace in the program | **IMPORTED AND MEASURED TO FAIL** |
| **A2** | there is a **background time parameter `t`** | steps 2–5; `dρ/dt`, `e^{−iHt}` | **IMPORTED, NEVER EXAMINED — 125 of 125 lanes** |
| **A3** | there is a **fixed self-adjoint `H`** generating evolution | steps 2–5 | imported; W-49 showed carriers can respond, so `H` fixed is a special case |
| **A4** | "durable" = **expectation value constant** | step 2 | imported; W-35 showed ensembles are blind to record formation |
| **A5** | individual **trajectories are physical** | step 4 | imported — this IS Gap 1 |
| **A6** | "the record" = **the slowest-decaying observable** | step 3 | a chosen criterion, never justified |
| **A7** | a **probe exists distinct from the record** and can move | step 5 | imported |
| **A8** | the bath is **Markovian** | steps 2–4 | imported — **CLEARED**: Phase B/B1 showed T1 needs no Markov assumption |

## A1 — WE MEASURED THIS FAILS AND KEPT USING IT

The register already carries it: *"In a gauge theory the physical (gauge-invariant) algebra does not
factorize into system and [environment]"*, logged as imported premise **IMP-2** (Donnelly–Wall,
Casini–Huerta–Rosabal), and the AA lane measured `physical + gauge = L` exactly. **W-37d broke on it
directly** — the `Z_3` measurement failed because the gauge field is not a tensor factor, and the fix
was to abandon partial traces for classical mutual information.

> **Every Holevo quantity, every partial trace, every "reduced state of the environment" in this
> program presupposes a factorisation the program itself measured does not exist.**

**And Gap 2 is this same assumption.** W-55 asked where factorisations come from and found: from the
modeller. **A1 and Gap 2 are one thing seen twice.**

## A2 — THE ONE NEVER EXAMINED, AND THE ROOT-CAUSE CANDIDATE

**125 of 125 lanes evolve in a background time parameter. Zero derive it.** Searching the entire
project for relational or emergent time returns nothing.

> **GRAVITY IS THE DYNAMICS OF SPACE-TIME STRUCTURE. IF `t` IS A FIXED EXTERNAL PARAMETER, THERE IS
> NOWHERE FOR GRAVITY TO ACT.**

**This explains all three failed gravity routes at once, and none of them touched it:**

| route | what it varied | what it left fixed |
|---|---|---|
| W-31 topology | which graph (space) | `t` |
| W-44 metric | link lengths, potentials, a clock **read against `t`** | `t` |
| W-45/W-50 markers | nothing dynamical at all | `t` |

**And it matches the paper's own statement of where gravity enters:** *metric/**proper-time** action*.
**We have no proper time.** We have one global background parameter shared by every part of the
carrier — which is exactly the structure general relativity does not have.

**W-44's clock build is the sharpest evidence.** It attached a clock to a probe — and the clock ran
against the same background `t` as everything else. **A clock that reads the background parameter
cannot register a difference in proper time, because there is only one.**

## WHAT THIS AUDIT DOES AND DOES NOT ESTABLISH

**ESTABLISHED:** A1 is imported after being measured to fail, and is load-bearing in every
information-theoretic quantity reported. A2 is imported, load-bearing everywhere, and was never
examined once in 125 lanes. A8 is cleared.

**NOT ESTABLISHED:** that removing A2 would produce a process. **This is a root-cause HYPOTHESIS, not
a result.** It is falsifiable in the useful direction: build a carrier with no background time — where
evolution is relational, one part read against another — and see whether anything gravity-shaped has
somewhere to act. If it still does not, A2 was not the obstruction.

> **THE HONEST ANSWER TO THE PRINCIPAL'S QUESTION IS YES. There are two imported assumptions at the
> surface. One we measured to fail and used anyway. The other we never looked at, and it is the one
> that would have to be wrong for gravity to have had anywhere to enter.**
