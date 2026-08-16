# WHERE ATOMS COME FROM — REGISTER V001

**PURPOSE.** Indexed by QUESTION, phrased the way someone would re-ask it. Before opening any
line of work, search this file for the question first.

Every row carries: the question as it would naturally be re-asked · the ruling · **where the
proof is** · and **the exact condition under which the row reopens**. A row with no reopen
condition is closed permanently.

**APPEND-ONLY.** Rows are never edited, only superseded by a later row citing the earlier one.
**This register adopts nothing and derives nothing.** It records where rulings live.

**POINTER RULE (custody §1):** no term appears in a row without a digest, a `file:line`, or a
named ruling behind it.

---

## HOW TO USE THIS REGISTER

Search by the question you are about to ask. If it is here, read the row and its reopen
condition before commissioning anything.

---

## W-01 — CAN A FORMATION CONDITION BE WRITTEN ON K1'S OWN DATA, WITHOUT A STIPULATED SPLIT? (S2; construction pair, build + adversarial audit. **AUDIT REFUTED THE BUILD AND WROTE THE CONDITION THE BUILD DECLARED UNWRITABLE.**)

**RULING: YES — WITH ZERO ADDITION AND NO IMPORTED SPLIT.**

**THE CONSTRUCTION.** Extend loop transport from a single fibre to the section space
`Gamma(L) = C^5` that the build itself derived:

```
(M_gamma s)(v) = W(gamma) · s(v)   for v on the loop
(M_gamma s)(v) = s(v)              otherwise
```

These are **unitary and non-scalar** on `Gamma(L)`. Then

```
<M_dF s, M_c s> = conj(W_F)·W_C·p_0  +  conj(W_F)·(p_1+p_2)  +  W_C·(p_3+p_4)
```

which vanishes **iff 0 lies in the convex hull of three unit-modulus coefficients** — verified
against brute-force simplex minimisation, 1369 grid points, 0 mismatches.

**IT FIRES ON S1'S OWN PUBLISHED CONNECTION.** With `W_F = -1`, `W_C = -i` (S1 §6) and
`p = (1/2, 0, 0, 1/4, 1/4)`, the overlap is `i/2 - i/2 = 0` **exactly**.

**PROPERTIES, ALL CHECKED:** gauge-invariant · correct trivial limit (no formation at
`W_F = W_C = 1`) · **the root can never fire**, independently reproduced · and it
**distinguishes curvature from flat holonomy**, which K1 exists to separate and which the
build's own rival construction could not see (its self-reported defect F1).

**THE ONLY LOAD-BEARING OBSTRUCTION IS DURABILITY.** Over `n` circuits the condition fires —
`0.0247` at `n = 42` — and then **recurs to `0.99994`**. K1 is finite and closed; this is a
**reversible write, not a durable record**, reproduced from first principles on five vertices
rather than inherited.

**WHAT THE BUILD CLAIMED, AND WHY IT WAS WRONG.** It returned a threefold impossibility —
rank, no time, no split. **Two are false and one is overstated.** Its rank theorem was proved
on `L_v0 = C`, an alternative its **own choice ledger had rejected**, and then applied to the
`Gamma(L) = C^5` it had chosen; `|<s, M_dF s>| = 0` on K1's own connection refutes it directly.
Its no-split theorem — resting on `dim = 5` being prime — is refuted: reduced states exist
(the build computed one), a canonical three-way direct-sum split exists, and `C^4 = C^2 (x) C^2`
on the non-root subspace. Only a weaker structural half survives: no *canonical* factorisation.
Its no-time claim is overstated — **circuit count is carrier-supplied discrete time.**
And its "minimum addition" question dissolves: **zero addition suffices.**

**A RESULT INSIDE THE CORRECTIONS THAT MATTERS ON ITS OWN.** A canonical three-way split
**exists on K1 and is derivable from its incidence** — `{v0} / {v1,v2} / {v3,v4}`: the root, the
filled triangle's other vertices, the unfilled triangle's other vertices. The
source/record/environment *shape* the predecessor introduced by imperative is **available here
without being imposed.**

**WHERE THE PROOF IS.** Build `S2_FORMATION_CONDITION_ON_K1_V001.md`
sha256 `248ce856efaef157c68e818dde589d0200bbc1dd9fd9fc1fcc8cdc7bc88734d9`.
Audit `S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md`
sha256 `0bea11bd4b7764f65c8d44cc8812d43bdcc569d34f51d0c46d82851a2efbd0d5`.
All 16 of the build's file:line pointers resolve; all 5 digests verify at bytes; alpha not
engaged. Numeric corrections carried: the build's §3.4 table is non-reproducible and its
corrected row values are `0.000000 / 0.816497 / 0.816497 / 1.000000`; the recurrence figure
`0.994373` is a window artefact and the true supremum is at least `0.999793`.

**GRADE: ADVERSARIALLY-CHECKED, NOT INDEPENDENTLY-CORROBORATED.** Build and audit share a model
lineage (custody §4). A failure mode common to that lineage passes through both.

**REOPENS IF:** a lineage-independent lane fails to reproduce the firing on S1's connection ·
the convex-hull criterion fails on a carrier where it should hold · or the three-way split's
derivability is shown to depend on K1's particular vertex count rather than on its incidence.

---

## W-02 — DOES A DIRECTED SYSTEM OVER K1 ESCAPE RECURRENCE, AND DOES IT PRODUCE THE NINE? (S3; construction pair. **AUDIT CONFIRMED-WITH-CORRECTIONS — every computable claim reproduced in a fresh independent implementation, no arithmetic error found anywhere.**)

**RULING: THE CROSSING EXISTS. THE LIMIT ESCAPES. SEVEN OF THE NINE HOLD OUTRIGHT; TWO HOLD IN THE LIMIT AND FAIL AT EVERY FINITE STAGE — WHICH IS WHY THE PREDECESSOR RECORDED ONE OF THEM FALSE.**

**THE SYSTEM.** `A_N = M_5(C) (x) M_2(C)^(x)N`, dimension `25·4^N`, embeddings
`iota_N(X) = X (x) I_2` — verified unital, multiplicative, *-preserving and isometric — with
`A_infinity = UHF(5·2^inf)`: simple, unital, infinite-dimensional. **Cost: one qubit per cell,
proved minimal. This is not zero addition.** The record slot is non-abelian by necessity: a
unitary write onto a reset carrier leaves *pure* slot states, and pure states of an abelian
algebra have overlap only 0 or 1.

**THE TRAP WAS REAL AND IS DISARMED BY COMPUTATION.** Repeated circuits of one loop span
**3 dimensions at N=1 and 3 at N=100.** Circuits grow no algebra and escape nothing.

**THE MECHANISM IS NOT WHAT IT LOOKS LIKE, AND THE AUDIT'S CORRECTION IS THE FINDING.**
The escape is **not** "recurrence is gone" — **the carrier's recurrence is untouched**,
`sup |Z_k| = 1`. It is **"the record grew."** And the monotonicity theorem carrying the escape
is a **tautology of any per-cell tensor system**: it holds identically on families that never
form at all. **All the content sits in the divergence `Sum (1 - z_n) = inf`.**

**AND THE ESCAPE IS CONDITIONAL IN TWO WAYS, NOT ONE.** The build flagged the **cell schedule**
— an adversarial schedule locked to the carrier's near-recurrences defeats the crossing
entirely. The audit found the second, flagged nowhere: **the support of the ready state.**

**THE CRITERION, CORRECTED AND PROVED, VERIFIED 10/10 ACROSS ALL FOUR SUPPORT CLASSES.**
With `u = conj(W_F)`, `v = W_C`, characters `chi_0 = uv`, `chi_F = u`, `chi_C = v`,
support `S = supp(p)`, and `G = <chi_a/chi_b : a,b in S>`:

```
FORMATION OCCURS  <=>  G != {1}
```

`|S|=3 -> G=<u,v>` · `S={0,C} -> <u>` · `S={0,F} -> <v>` · `S={F,C} -> <u/v>` ·
`|S|=1 -> never`, which **recovers W-01's "the root can never fire" as a special case.**
**Formation is a group-theoretic condition on where the record sits.**

**THE NINE.** P-2 persistence · P-3 thresholded non-return (**exactly, not asymptotically**) ·
P-4 recoverability · P-5 redundancy · P-6 asymptotic centrality · P-7 sector-hood ·
P-8 inductive compatibility — **all HOLD**, each with its computation.
**P-1 durability-with-irreversibility and P-9 orthogonal reduced supports HOLD IN THE LIMIT AND
FAIL AT EVERY FINITE STAGE**, exhibited with an explicit reversing unitary at N = 1..5.
**P-9 is resolved: `||omega_F^N - omega_C^N|| -> 2.000000000000`. The predecessor recorded it
false because it is a limit property and was evaluated at a finite stage.**

**TRANSPORT (Q5): YES, MULTIPLICATIVELY.** `Phi_{N+1} = Phi_N · Z_{N+1}`, and the difference
functional is compatible with the embeddings. **Firing is absorbing.** And the limit converts
W-01's **binary** convex-hull criterion into a **rate** `lambda = -0.766802`, defined everywhere
off the trivial connection.

**CORRECTIONS THAT MATTER.** A claimed unital embedding of `M_5 (x) M_2^{(x)N}` into
`M_5^{(x)(N+1)}` is **false for every N**, by the divisibility criterion the same section
states. And the universal claim that every non-trivial connection forms is **false** — four
exhibited families never form, one of them **on K1's own published ready state**.

**WHERE THE PROOF IS.** Build `S3_THE_CROSSING_V001.md`
sha256 `cbf1d79679ca2ecf3ee260e8a6467062e3a93260b2325cfc287c261e9b4469cb`;
audit `S3_THE_CROSSING_AUDIT_V001.md` per its sidecar. The build self-reported eleven defects
including a one-character digest error in its own custody block, caught by the pointer rule and
recorded rather than silently fixed.

**GRADE: ADVERSARIALLY-CHECKED, NOT INDEPENDENTLY-CORROBORATED** (custody §4).

**REOPENS IF:** a schedule condition is found that is carrier-intrinsic rather than adversarial ·
the ready-state support condition is shown to follow from the incidence rather than being free ·
the rate `lambda` fails to vary with the connection at S4 · or a lineage-independent lane fails
to reproduce the character-ratio criterion.

---

## ERRATUM AGAINST W-02 — THE RATE AS REGISTERED IS WRONG (registrar error, 2026-08-16)

W-02 records `lambda = -0.766802`. **That number is wrong under every reading and the error is
mine, not the build's.** It is S3's **N = 4000 finite-stage value of an oscillating sequence**,
transcribed into the register as if it were a converged rate. The corrected statement of record:
S3's headline connection `f = 2.0, c = 1.1` is **exactly resonant** (`-11f + 20c = 0`), so its
orbit is not dense in `T^2` and its average is taken over a **subtorus**; it converges to
**`-0.767014993`**, while the generic torus value is **`-0.767507880` = m(0.4+0.3x+0.3y)**.
S3 section 6(f)'s "rationally independent: orbit dense in T^2" is FALSE and section 4.3's
"Weyl equidistribution AGREES" is a false agreement. Every row of S3 section 5.7 is a subtorus
value. **No W-02 verdict is withdrawn** — the crossing, the UHF limit and the nine stand — but
the mechanism sentence and the rate are corrected here.

---

## W-03 — S4, THE MEASUREMENT, AND ITS AUDIT: **THE HEADLINE IS REFUTED. THE STRUCTURE IS NOT UNIVERSAL.**

**RULING: "THE TOPOLOGY IS INERT; FORMATION SEES LOOP INCIDENCE" IS NOT EARNED. THE FIRST HALF
IS TRUE BUT IS A THEOREM OF THE TRANSPORT'S DEFINITION, NOT A FINDING OF THE CONTROLS. THE
SECOND HALF IS FALSE.**

Build `S4_THE_MEASUREMENT_V001.md` (`4277f0ef…`). Audit: six independent lanes re-deriving from
scratch with their own code, refuters attacking each lane's load-bearing claim, a completeness
critic. **31 agents. Every reproducible number in the build reproduced** — worst agreement
`0.0e+00`, and thirteen of thirteen rows of the S3 corpus block verified for the first time.
**The arithmetic was never the problem.**

**WHY THE CONTROLS COULD NOT HAVE FAILED.** All four controls vary exactly **one** object: the
pushforward of the ready state onto the four vertex classes, `pi = (p_00,p_10,p_01,p_11)`.
Control 1 (the fill) is **vacuous** — `pi` is unchanged, so `lambda` is unchanged **by
identity**. And `d2` — the 2-cells, the entire topological content — **enters the formation
functional nowhere**. "Topology is inert" is therefore an **analytic fact about the
construction**, true before any carrier was built. The controls tested nothing.

**AND FORMATION DOES NOT SEE INCIDENCE.** `lambda_B` is a function of the **MULTISET** of the
four class weights — 24 of 24 permutations invariant, worst spread `2.4e-15`. **The incidence
labels are invisible.** Two exact refutations of the build's own key exhibits:
- A carrier **retaining the pinch** (filled triangle + unfilled square at `v0`) reproduces the
  full `0.0634` shift, `|diff| = 0.0e+00`. The pinch is **not necessary** for it.
- Moving weight **off** the pinch with the incidence **untouched** closes the gap continuously
  to exactly zero. The pinch is **not sufficient** either.
- **PINCH AND SPECTATOR ARE THE SAME OBJECT**, exactly: multiplying `Z_k` by
  `conj(u)^k conj(v)^k` leaves `|Z_k|` fixed and maps class `(a,b) -> (1-a,1-b)`, so
  `00 <-> 11` and `10 <-> 01` is an exact symmetry at **every** connection (max deviation
  `6.55e-15` over 2000 samples).

**THE MODALITY NAMED IN S4'S OWN BRIEF AND NEVER RUN: CHARGE.** Named again in the sealed
upstream S2 as "**this is the S4 knob**". The build never ran it; six audit lanes and their
refuters never noticed. The critic ran it: **THEOREM S4-1 FAILS.** Exponents `(1,0),(2,0),(3,0)`
give `|S| = 3` with `rank G = 1`, against S4-1's "rank 2 iff `|S| >= 3`"; and `q = (1,2,2,2,2)`
moves `lambda` from `-0.767508` to `-1.200555`. **The four-class taxonomy is a charge-1
statement.** Winding was also never non-vacuously run; the critic ran it and it **discriminates
the two schedules** — bearing directly on F1, the flag S4 declared undecidable.

**THE UNIVERSALITY QUESTION, ANSWERED EARLY AND NEGATIVELY.** Two of S4's three deepest results
are theorems about the **character lattice of an abelian group** — "`lambda_B` is a function of
`L` alone" is Pontryagin duality on `T^2`, and Theorem S4-1 is the four corners of a square in
`Hom(U(1)^2, U(1)) = Z^2`. **Neither has an SU(2) analogue.** Section 2's foundational sentence,
"every vertex phase of `s` cancels", is **false the moment the holonomy is not a scalar**:
six ready states with **identical class weights** give `|Z_1|` spread **0.4247** under SU(2),
against exactly `0.0` under U(1). **The construction's core is group-specific.**

**FURTHER CORRECTIONS SURVIVING ATTACK.** `B0b` is `log(4/9) = -0.8109302162163288` exactly, so
**nine of nine** carrier rates are exact, not eight; `(pi,pi)` is a **third strict saddle**, not
a local minimum (the build's own table contradicts `chi(T^2) = 0` at `1-2+2 = +1`); the
exceptional-value split is **527/314/213**, not 638/380/258; section 5 misquotes a reopen clause
**the register does not contain**. Lawton's theorem (1983) is missing from the IMPORT AUDIT.

**THE CARRIER AXIS IS A LIST, NOT A FAMILY.** No cellular collapse, quotient, or map connects
any two of the ten complexes. And the entire C*-algebraic layer of W-02 is **absent from S4** —
which silently moves the first factor from `M_5(C)` to `M_V(C)` as `V` runs 4 to 11, and never
says so.

**CORRECTED HEADLINE OF RECORD:** *under a U(1) connection at unit charge with two designated
loops, formation's rate is a function of exactly two things — the pushforward of the ready state
onto the four vertex classes, and (under the canonical clock) the relation lattice of the two
holonomies; the carrier's topology is inert **by construction rather than by experiment**, the
incidence labels are invisible off a dense measure-zero set, and whether any of it survives
charge, a non-abelian group, or a third loop was never tested.*

**GRADE: ADVERSARIALLY-CHECKED on the axes it ran; NOT INDEPENDENTLY-CORROBORATED, and not
adversarially checked at all on the axes it did not run.** Route to corroboration specified:
a lane outside this model lineage; the charge and three-loop runs with the four controls re-run
under each; publication of `d1`, `d2`, grid and seed conventions; this erratum entered.

**REOPENS ON:** the charge run · the SU(2) run · a third loop · a third schedule (only two exist
anywhere on the page, and the A/B dichotomy may be a knife-edge between two arbitrary choices).
