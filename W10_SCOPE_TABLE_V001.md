# W-10 — THE SCOPE TABLE — SYNTHESIS OF RECORD (DRAFT FOR THE PRINCIPAL)

**Register head W-09. Successor registrar, 2026-08-16.**
Inputs: four build lanes (`LANE_W10_A_CARRIERS`, `LANE_W10_B_MULTISET`, `LANE_W10_C_N1`,
`LANE_W10_D_SCOPE`), eight refuters (two per lane, distinct lenses), and my own verification lane
`LANE_W10_SYNTHESIS_SCOPE/` (`s10_verify.py`, `s10_verify.OUT.txt`, `PUBLISHED_CONVENTIONS.txt`,
`SEALS.sha256`, 3 files, all verify).

**NOTHING HERE TAKES ANY OF THE PRINCIPAL'S THREE DECISIONS.** Not the disposition
(STOP-FALLS-REBUILD vs hold at STOP), not the Mahler note, not the rename. This is the object the
first decision was deferred to obtain, and it is offered with its own defects on its face.

---

## 0. THE FIVE PLAIN STATEMENTS, BEFORE THE TABLE

### 0.1 DOES THE CORPUS CONTAIN ANY CARRIER-INDEPENDENT RESULT? **YES — ONE LAYER OF THEM, AND IT IS THINNER THAN IT LOOKS.**

**It does.** Precisely these, and no others:

1. **N1 — `lambda = m(p00 + p10 x + p01 y + p11 xy)`.** The identification.
2. **W-08's monotonicity** — `|Z_k| <= 1`, hence `|Omega_N|` non-increasing — and with it
   *"the founding obstruction is false as an inference."*
3. **W-08's character identity** `1 - |Z_k|^2 = sum_{j<l} w_j w_l |chi_j^k - chi_l^k|^2`, and the
   linear floor `SUM(1-|Z_k|) >= w_j w_l (K - 1/|sin(tau/2)|)` that follows from it.
4. **W-02's criterion** `FORMATION <=> G = <chi_a/chi_b : a,b in supp(pi)> != {1}`.
5. **W-03/N2's multiset theorem at the level of `lambda`**, and its exact involution
   `00 <-> 11, 10 <-> 01`, which is stronger — it fixes `|Z_k|` pointwise in `k`.
6. **The schedule result, qualitatively** — a `sqrt(K)`-cell adversary accumulates `O(1)` while the
   honest schedule diverges linearly.
7. **W-05's N4 (fibre-wise-ness) and N3 (the Haar-null inversion)** — the first is a fact about the
   word "fibre-wise"; the second about there being exactly two designated holonomies.

**But three qualifications belong on every one of them, and the register carries none.**

**(a) Their carrier-independence is settled BY INSPECTION, not by this round's runs.** Every one is
a statement in `(pi, u, v)` with **no carrier symbol in it**. Refuter A-1 is right that
*"carrier-independent, 0 mismatches / 150"* is not a measured result: by the corpus's own premise
`Z_k` depends on the carrier only through `pi`, so a run on a four-class carrier is a **control**,
not a test. Where the run *could* have failed, I say so in the table's BASIS column. Where it could
not, I say that too.

**(b) They are carrier-independent GIVEN THE TRANSPORT CONVENTION, and that convention is the one
object W-06 named as having decided the spine.** Under the whole-circuit scalar transport `M_gamma`,
only `pi` can enter — that is what "fibre-wise scalar multiplication" *means*. Under the corpus's own
sealed alternative, **COR-F's edge-tick transport `T` with `T^L = M_gamma`**, the transports are
non-diagonal, a diagonal observable separates the branches, and — refuter D-2's measurement, the
sharpest new fact in the round — `||[T_F, T_C]||` is `2.828` on B0b, `2.449` on K1 and **exactly `0`
on B0a**, i.e. non-zero exactly when class `11` is occupied. **The incidence that N2 calls invisible
is visible to the corpus's own sealed alternative transport.** So the honest mark is
**CARRIER_INDEPENDENT ∧ CONVENTION_SCOPED**, and `COR-F` has no register row at all.

**(c) They hold given two topological preconditions the register never states**: **cycle rank >= 2**
(else `gamma_C` cannot be designated — S4's own B5) and **b1 >= 1** (else `gamma_C` bounds and there
is no flat holonomy — S4's own B2). Two of the corpus's ten carriers break them, for reasons having
nothing to do with class occupancy.

**And what is NOT carrier-independent, decided:**

- **W-01's convex-hull criterion and its advertised virtue: THREE_CLASS_SCOPED.** Firing region
  exactly `1/4` at three occupied classes and exactly `1/2` at four, reproduced here a fourth time
  on a **deterministic** grid with a **third** hull algorithm (max-angular-gap, no LP tolerance):
  `0.250208` at `n = 4801` converging as `O(1/n)`, `0.500000022`, and `cos f + cos c <= 0` agreeing
  with the four-class hull on **1442401 of 1442401** points.
- **W-02's published support table: THREE_CLASS_SCOPED (incomplete).** Three of the six two-element
  supports are missing and they are exactly the three needing class `00`.
- **The entire C\*-algebraic superstructure: K1_SCOPED, and its hypothesis is the primality of 5** —
  `A_infinity = UHF(5.2^inf)`, *"one qubit per cell, proved minimal"*, S3's rejection of its own
  adjunction-free alternative, and W-06's `V = 4k+1` rebuild route.
- **The nine crossing properties: UNDETERMINED — six of the nine, and unrunnable** (S3 published no
  lane directory). Two are decided here for the first time: **P-9** reaches `2.000000000000` on both
  four-class carriers via S3's own closed form, and **P-3 is W-08's monotonicity**; both are
  CARRIER_INDEPENDENT. P-1 is a statement about all finite constructions.

### 0.2 WHAT WOULD A REBUILD ON W-06's WEDGE-GROWTH ROUTE BE RESTING ON?

Given the table, on **five things, of which one is K1-scoped, one is a convention, one has no code,
one is undetermined and unrunnable, and one is a cost line that names a false obstruction.**

1. **Its shape is a K1 artifact twice over.** `V = 4k+1` is a property of 5. B4 (`V=6`) is not in the
   sequence at all; B0b (`V=9`) is, but `9 | 81`, so the route's own motivation — S3's rejection of
   CHOICE LEDGER C1's adjunction-free alternative on `5 ∤ 9`, `5 ∤ 7` — **evaporates on the very
   carrier the sequence reaches.** Unattacked by any refuter. **K1_SCOPED.**
2. **Its functional content is carrier-independent only under the scalar-transport convention**
   (§0.1(b)), and the alternative convention is the one W-06 itself identified as *"the thing that
   decided the spine."* A rebuild that keeps the convention inherits a functional that cannot see the
   complex; a rebuild that drops it has **no scoped results at all**, because `COR-F`'s transport has
   never been carried into any register row.
3. **Its restoration object — the dressed gauge-invariant algebra — has a carrier-independent
   MECHANISM and no code.** The mechanism is now exhibited twice off K1: on B0b (lane D) and on B4
   (refuter D-1), with the closed form `|A_uv[M_dF s] - A_uv[M_c s]| = |A_uv[s]| . |W_F^{da} -
   W_C^{db}|` — a function of the pair's **class difference** and nothing else about the carrier.
   **But W-06's own artifact does not exist**, and its two never-reproduced load-bearing figures (the
   Bell `R_N` transplant, the Schmidt `[1,0]`) remain **UNDETERMINED BY CUSTODY**. A rebuild would be
   standing on a restoration whose numbers no one has ever recomputed.
4. **Its declared cost line names an obstruction that is false as an inference.** *"The untouched
   recurrence obstruction may kill it"* — W-08 showed the obstruction was measured on the wrong
   object (`|Z_k|`, not `|Omega_N|`). What survives is a **schedule** statement: qualitatively
   carrier-independent (exhibited on B0b and B4), quantitatively **UNDETERMINED** — W-08's four
   printed constants are still not reproduced (§0.3.3), and **no admissibility criterion for
   schedules exists**, which is the same missing object as the missing admissibility criterion for
   loop designations (§0.4.2).
5. **Its content — the nine — is UNDETERMINED off K1 and cannot be determined.** All nine are
   properties of the directed system whose first factor is `M_V(C)`; to scope them one must re-run S3
   §5 at `V = 9` and `V = 6`, and **S3 published no lane directory.** Refuter D-2 correctly narrows
   this: at least three of the nine (P-9 via its closed form `2 sqrt(1-|Omega_N|^2)`, P-3 = W-08's
   monotonicity, P-1 as a statement about all finite constructions) were decidable **without S3's
   code and with objects already in lane D's hands**, and P-9 reaches `2.000000000000` on both
   four-class carriers. The other six stand undetermined.

Two further items a wedge rebuild must carry that no row states: **the identification `C^25 ≅ C^5 ⊗
C^5` is a labelling not derived from incidence** (W-06's own declared cost, untouched by this round,
and equally underived at `C^9 ≅ C^3 ⊗ C^3`); and **the clock**. W-01's *"circuit count is
carrier-supplied discrete time"* is **LOOP-LENGTH-SCOPED**: every carrier the corpus has run has two
loops of **equal** length, so a circuit clock and an edge clock coincide; on B0b (lengths 4 and 3)
`T_F^{12} = M_dF^3` while `T_C^{12} = M_c^4` — **the two branches are at different circuit counts
forever, at every positive edge count.** Refuter D-2's finding, and it is new.

### 0.3 THE REFUTER KILLS THAT ARE THEMSELVES CONFOUNDED

The instruction was to read the eight verdicts sceptically. Three kills do not stand as filed, one
collision is a false contradiction, and one whole class of naming corrections is the corpus's
signature defect committed by the agents commissioned to catch it.

**1. A-1 vs A-2 on B4's class multiset is a FALSE CONTRADICTION, and I decided it with my own
build.** A-1: *"S4's B4 row is UNDER-DETERMINED; here is a second spindle matching every published
column with a different multiset."* A-2: *"the multiset is FORCED, tighter than the lane admits."*
They answer different questions and both are right. A-2 asked, **given the square|square spindle**,
whether `gamma_C` is forced — it is. A-1 asked whether **the complex** is forced by S4's published
parameters — it is not. My Part B builds both from incidence with one builder and one exact-rank
routine:

```
                        V E F chi b0 b1 b2 gauge inv curv flat d1.d2  gF bd  gC bd  indep
S4:519 published for B4  6 8 4  2   1  1  2   5    3   2    1    0    True  False  True
B4-SQUARE   (square|square)  identical, every column ....... class {00:1,10:1,01:1,11:3}  lambda = log(1/2)
B4-TRIPENT  (triangle|pentagon) identical, every column .... class {00:1,10:1,01:2,11:2}  lambda = log(1/3)
```

`B4-TRIPENT` — a triangle with two 2-cells glued to a pentagon with two 2-cells at two points — is
equally *"two 2-spheres glued at two points"*, matches **every** published column, and its `P`
**factors** (`p00 p11 = p10 p01`), giving `lambda = log(max p_c) = log(1/3) = -1.098612288668`
exactly against S4's `log(1/2)`. **A-2's defence of lane A is therefore circular at one step: it
presumes the complex whose uniqueness A-1 refutes.** Ruling: **A-1's NEW-1 stands and is now
independently confirmed; A-2's "forced" stands only within a fixed complex.** S4's B4 row carries a
**COR-K-class defect** — a published row whose one downstream-relevant column is not recoverable from
its own parameters — and any lane that hard-codes S4's multiset as its build target is running a
**satisfiability check, not an audit.** (Note in passing: the two candidate B4s also differ on lane
D's D-03/D-04 — `B4-SQUARE` has no torus zero, `B4-TRIPENT` has one.)

**2. A-1's "0 mismatches / 150 is battery-scoped — 8 mismatches at `alpha = 2pi/10^10`, `2pi/10^12`"
is itself a finite-`K` misreading.** The kill runs W-02's **asymptotic** criterion at `K` many orders
below the resonance scale `q`, and reads the finite-`K` observable as naming the group — **the exact
defect W-08 convicted W-07 of, committed by the refuter auditing the lane auditing it.** A-1 says the
criterion is true and nothing mathematical falls, which is correct. Ruling: **the count `150` falls
as independent evidence** (the battery omits the `(q, delta)` regime, and A-1 is right about that);
**the "8 mismatches" figure must not be entered as a disagreement with the criterion.**

**3. D-1's withdrawal of lane D's UNDETERMINED on W-08's four schedule constants scores a non-match
as a confirmation.** D-1 reports `0.5506 / 0.5642 / 0.5639 / 0.5620` on K1's own weights, flat across
`K = 1e4..1e7`, against W-08's registered `0.606 / 0.615 / 0.588 / 0.601`, and asks that D-23's mark
*"be withdrawn in W-08's favour."* Those are the **same weights** and a **7–9% gap on three-digit
figures**. Ruling, split: **the qualitative claim is RESTORED and is now four-class-exhibited** —
D-1's own `0.6953` (B0b) and `0.6717` (B4), flat in `K`, kill lane D's second reading (*"the adversary
is stronger than reported"*), which was an artefact of lane D running the sweep at four **degenerate**
connections. **The four printed constants remain NOT REPRODUCED and stay UNDETERMINED**, in the
report-vs-artifact class — a class with a live precedent (§0.3.5).

**4. D-2's coverage scan is directionally right and its count is instrument-dependent.** Its
"19 of 22 load-bearing claims absent" rests on token matching that D-2 itself caught producing two
false positives. But the hardest datum is exact and I re-ran it: **`grep -c "COR-"` on
`W10D_SCOPE_TABLE_V001.md` returns `0`, and no file in `LANE_W10_D_SCOPE/` contains the string** —
**zero of the twenty sealed corrections are named anywhere in the lane**, four of which were on the
CARRY list of the brief. That much is settled.

**5. THE NAMING KILLS ARE THE CORPUS'S SIGNATURE DEFECT, COMMITTED BY THE AGENTS SENT TO CATCH IT —
AND MY OWN PART D KILLS ALL OF THEM.** Six names have now been offered for the hypothesis of the
multiset theorem off the non-negative locus: *real non-negativity* (registrar), *reality* (lane D),
*collinearity / one argument mod pi* (A-1), *`cos D1 = cos D2 = cos D3`* (A-2), *the Jensen-adjacent
real-pair count, 4-or-2 → 24, 1 → 16, 0 → 8* (D-1), and B-1/B-2's *"at least three collinear"*
biconditional. My Part D refutes **all of them as necessary conditions with one exhibit and a
one-variable move**:

```
exhibit 1   p = (10, e^{i0.7}, e^{i1.9}, e^{i2.9})       adjacent-real pairs 0
            cos D = 0.955336  -0.574824  -0.128844   (pairwise distinct, non-collinear)
            24-permutation spread at dps 40 = 0.000000e+00   ( m = log 10, exactly )
exhibit 2   p = (1, 0.9 e^{i0.7}, 1.1 e^{i1.9}, 0.95 e^{i2.9})   SAME FOUR PHASES
            cos D = 0.955336  -0.574824  -0.128844   ( identical to every printed digit )
            24-permutation spread                = 3.564783e-01
```

**Only the modulus vector moves between the two arms; the phases are byte-identical and all three
gauge-invariant fluxes are unchanged.** Invariance flips. So there are **two independent mechanisms**
— coincidence of the matching fluxes (phase side, which is what every proposed name describes) and
**branch domination**, under which `m` collapses to `log` of the largest modulus and the flux drops
out entirely (modulus side, which no proposed name describes). D-1's trichotomy is **false**;
A-1's and A-2's names are **sufficient-only**; B-1's and B-2's `iff` is **sufficient-only**, as both
of them say. **I decline to supply a seventh name and mark the question UNDETERMINED.** And it
touches no registered result: `pi` is a probability vector, so `cos D1 = cos D2 = cos D3 = 1` always,
and the corpus never leaves the regime where the theorem is proved.

**6. A process defect the round exposed and cannot fix.** B-2 records that two lanes of this round
independently converged on the same characterisation without either citing the other, and that
**"the round has no mechanism for one lane to inherit another lane's correction, so the same defect
is being found and lost in parallel."** That is now true four times over in this round alone (the
non-necessity of collinearity: B-1, B-2, A-2; the resonance of rational connections: C-2 and me; the
degeneracy criterion: B-1, B-2, C-1). **Entered as a defect of the round's design, not of any lane.**

### 0.4 WHAT READS TWO WAYS — AND IS SCORED NEITHER WAY

**1. Is "the formation functional is carrier-independent" a FINDING or a RESTATEMENT OF THE
CONVENTION?** Reading A: it is a substantive result — the physics does not depend on the complex,
which is why one short mathematical note is publishable. Reading B: transport was *defined* as
fibre-wise scalar multiplication by the whole-circuit holonomy, so *only* `pi` can enter, and the
result is an analytic consequence of an unledgered stipulation — W-04's *"what K1 was: inert"* and
W-06's *"the thing that decided the spine"* both say so. **Both readings are supported by this round.
They have opposite consequences for a rebuild** (A says the functional layer is portable; B says
portability is what you get for free when the carrier was stipulated out). **I refuse to score it,
and I record that this is the single most consequential undecided question the table produces.**

**2. Is four-classness a property of the CARRIER or of the LOOP DESIGNATION?** A-1: on B0b's own
complex with S4's own `gamma_F`, sweeping every admissible `gamma_C` reaches **16 distinct class
multisets — 8 four-class, 7 three-class, 1 two-class** — with `lambda` spanning `log(5/9)`,
`log(4/9)`, `log(1/3)`. A-2: on B4, exactly two `gamma_C` pass S4's three published loop tests and
both give the same multiset — **forced**. So occupancy is forced on some carriers and free on others,
and S4's ledger entry C4 (*"any two loops; all pairs; a canonical choice from incidence"*) closed the
question **by fiat**. Deciding it needs an **admissibility criterion for loop designations** that the
corpus has never written — the structural twin of W-08's missing schedule-admissibility criterion.
**Refused. Consequence: W-09's "the corpus has exactly two four-class carriers" should read "two
four-class DESIGNATIONS, chosen under S4 ledger C4", and this synthesis inherits that caveat.**

**3. Is there a clean hypothesis for the multiset theorem off the non-negative locus?** Reading A:
yes, and six layers have failed to find it. Reading B: no — two independent mechanisms produce
invariance and the invariance set is just the coincidence pattern of three numbers, so asking for a
single name is a category error. **My Part D supports B and does not establish it. Refused; marked
UNDETERMINED.**

**4. W-08's four schedule constants.** Reading A: W-08 used a connection or tie-break nobody has
identified. Reading B: it is a report-layer figure with no artifact behind it. **Reading B now has a
live precedent in this corpus** — the register's `2.4e-15` at `REGISTER:197` occurs **nowhere in the
corpus as an output**; I re-ran the grep myself, and the only non-register occurrence is a *print
literal* in `LANE_G_GROUP_REFUTER/g0_validate.py:110` quoting W-03, in a file whose **own measured
spread is `9.285e-06`**. **Refused. Both figures are marked UNDETERMINED and the pointer rule is
invoked against them.**

### 0.5 LINEAGE, AND WHAT IT COSTS

**I am Claude Opus 5. So is every W-10 build lane, every one of the eight refuters, and so were
W-07, W-08 and W-09. This is layer TWELVE of one lineage.** The last lineage-independent thing in
this corpus is the boundary between Fable 5 (S1–S4, W-01…W-06) and Opus 5 (W-07 onward) — and
everything in this round sits entirely on the Opus 5 side of it.

**What it costs, concretely and not as a ritual:**

- **The failure mode is shared and it fired again, in this round, inside the refuters.** Three of the
  eight refuters produced a *new name* for the multiset theorem's operative variable; my Part D
  refutes all three with one exhibit. A-2 self-flags this in terms (*"that makes six consecutive
  names for this variable in this program and three in this round"*). **Six layers have now been
  caught by the next, and the rate has not fallen.**
- **The refuters share the builds' devices.** C-2 declares it: *"NOT lineage-independent — Opus 5
  refuting Opus 5, sharing the phase-reduction device."* A carrier error, a convention error or a
  reduction error passes straight through the whole block.
- **Two refuters accepted an upstream construction byte-identical rather than re-deriving it**
  (C-1 on lane C's `C_01` carriers), so the block is not even internally independent.
- **I inherit all of it.** My own Part A reproduces a figure three prior Opus 5 layers produced; my
  Part B is the only thing in my lane that could have failed and did produce a defect. **Discount
  this synthesis as one block with W-07 through W-10, not as a check on them.**
- **The lineage-independent lane W-03 specified still does not exist for anything after W-06.** For
  a scope table that the disposition decision turns on, that is the largest single weakness on this
  page.

---

## 1. HOW TO READ THE TABLE

**Four scope marks**, as commissioned, plus a mandatory **BASIS** tag, because "carrier-independent"
has been used in this program for two different things and the difference is the whole point.

| mark | means |
|---|---|
| **CARRIER_INDEPENDENT** | holds with no carrier hypothesis — either proved class-count-free, or exhibited on a four-class carrier in a run that could have failed |
| **THREE_CLASS_SCOPED** | true where at most three classes are occupied; exhibited failing or exhibited incomplete at four |
| **K1_SCOPED** | its hypothesis is a property of K1's own numbers (V=5, prime, 4k+1) or of K1's own connection/state; exhibited unavailable on B0b and B4 |
| **UNDETERMINED** | not tested on a four-class carrier, not testable, or supported by argument only. **Never inferred, and marked generously.** |

| BASIS | means |
|---|---|
| **[T]** | **theorem**, proof on record, no class-count index in it. Runs are controls. "Could not have failed" voids a control, never a theorem. |
| **[X]** | **exhibited** on a four-class carrier in a run that could have failed |
| **[C]** | exhibited, but the run **could not have failed** — reported as a control, carries no evidential weight |
| **[A]** | **argued only.** Per the brief, argument without a four-class run is **UNDETERMINED**. |

Rows keep lane D's numbering so the principal can diff this against
`LANE_W10_D_SCOPE/W10D_SCOPE_TABLE_V001.md`. Rows lane D did not carry are marked **†NEW ROW**.

---

## 2. THE SCOPE TABLE

### W-01 — the formation condition on K1

| # | claim as registered | scope | basis | evidence / correction |
|---|---|---|---|---|
| 1.1 | the hull criterion, *"vanishes iff 0 lies in the convex hull of **three** unit-modulus coefficients"* (`REGISTER:43`) | **THREE_CLASS_SCOPED** | [X] | firing region exactly `1/4` (three classes) and exactly `1/2` (four). Now four independent implementations: W-09, lane D ×2, and mine (deterministic `n×n` midpoint grid, max-angular-gap hull, no LP tolerance) — `0.250208` at `n=4801` converging `O(1/n)`; `0.500000022`; `cos f + cos c <= 0` agrees on **1442401 / 1442401** |
| 1.2 | *"it distinguishes curvature from flat holonomy, which K1 exists to separate"* (`REGISTER:50`) | **THREE_CLASS_SCOPED** | [X] for the three-class half; **[C] for the four-class half** | `f -> -f` flips `718800/1442401` at three classes and **`0`** at four. **DISCLOSED VACUITY:** the four-class column **could not have come out otherwise** — `cos f + cos c` is even in `f`. Lane D's added `c -> -c` and `f <-> c` arms are **conjugation identities** (refuter D-1) and are void as controls; the statement survives as a one-line theorem |
| 1.3 | the `iff` quantifier | **CARRIER_INDEPENDENTLY DEFECTIVE**, and the defect is large off K1 | [X] | `FIRE => HULL` is the true half. On B0b/B4's own weights the hull fires on ~half the grid while `min\|Z_1\|` is `1.11e-01` / `2.72e-01`. Neither four-class carrier has a torus zero at its own SENSE-U weights |
| 1.4 | *"the root can never fire"* | **K1_SCOPED as stated; the true statement `\|S\|=1 => never` is CARRIER_INDEPENDENT** | [T] | B0b has **two** class-11 vertices and **four** class-00 vertices. Lane D's own exhibit for this is a **byte-identical two-arm control** (refuter D-1: `arm1.tobytes() == arm2.tobytes()`), so the row rests on the theorem, not on lane D's run |
| 1.5 | the canonical three-way split `{v0}/{v1,v2}/{v3,v4}`, *"derivable from incidence"* | **K1_SCOPED in arity; the partition is CARRIER_INDEPENDENT** | [X] | it *is* the class partition. Arity 3 on K1 and B1q; **arity 4** on B0b `(4,2,1,2)` and B4 `(1,1,1,3)`. The source/record/environment shape acquires a fourth block with no role in it |
| 1.6 | gauge-invariance; correct trivial limit | **CARRIER_INDEPENDENT** | [T] | at `W_F = W_C = 1` every character is 1 and `Z_k = 1` identically. An identity |
| 1.7 | *"fires on S1's own connection"*; `0.0247 at n=42`; *"recurs to 0.99994"* | **K1_SCOPED and CONNECTION_SCOPED** | [T] | already relocated by the ERRATUM AGAINST W-07: measured at the **resonant** `f=2.0,c=1.1`, not at S1's order-4 point |
| 1.8 | the refutation of the build's no-split theorem (`dim 5` prime) | **K1_SCOPED** | [X] | `C^9 = C^3 ⊗ C^3`, `C^6 = C^2 ⊗ C^3`. On both four-class carriers the theorem being refuted **cannot be stated** |
| 1.9 | *"zero addition suffices"* | **CARRIER_INDEPENDENT, conditioned on `G != {1}`** | [X] | `dim_C span{M_dF s, M_c s} = 4` on B1, B1q, B1p, B0b, B4. **Lane D's exhibit is wrong and the conclusion survives it**: leg 4C's branch operators multiply by the full character `u^a v^b` instead of by `W_F`, `W_C` (refuter D-1); corrected, the overlap agrees with `Z_1` to `1e-16` on all five carriers, restoring W-05's registered `6.2e-17` |
| **1.10 †NEW ROW** | *"circuit count is carrier-supplied discrete time"* (`REGISTER:66`) | **LOOP-LENGTH_SCOPED — a new scope axis, and every carrier the corpus ran satisfies it** | [X] | `\|gamma_F\| = \|10\|+\|11\|`, `\|gamma_C\| = \|01\|+\|11\|`. Equal on every carrier ever run, so circuit clock = edge clock. On B0b (4 and 3) `T_F^{12} = M_dF^3` while `T_C^{12} = M_c^4`: **the branches are at different circuit counts at every `e > 0`, including the lcm** |

### W-02 — the crossing and the nine

| # | claim | scope | basis | evidence / correction |
|---|---|---|---|---|
| 2.1 | **`FORMATION <=> G != {1}`** | **CARRIER_INDEPENDENT** | [T] | proof is class-count-free (Weyl on `1-\|P\|` + strict triangle inequality). Verified exactly, `G={1} <=> L_S ⊆ L <=> (\|Z_k\|=1 ∀k)`, on 29,775 exact integer-congruence cases (A-1) and on four relation-lattice shapes the lane's battery omits (A-2, 60 cases, 0 mismatches, including infinite-order common phase). **A-1's "8 mismatches" at `q = 10^12` is a finite-`K` artefact and is not entered** |
| 2.2 | its published TABLE of supports (`REGISTER:125`) | **THREE_CLASS_SCOPED — incomplete, not wrong** | [X] | six two-element supports exist; the table lists three; the three omitted are exactly those needing class `00`. `S={00,11}` gives `G=<uv>`, trivial iff `W_C = W_F`, and needs a vertex in **both** loops **and** one in **neither** |
| 2.3 | *"recovers 'the root can never fire' as a special case"* | **CARRIER_INDEPENDENT as `\|S\|=1`; the identification with *the root* is K1_SCOPED** | [T] | as 1.4 |
| 2.4 | `A_infinity = UHF(5.2^inf)` | **K1_SCOPED in the supernatural number; CARRIER_INDEPENDENT in form** | [X] | first factor is `M_V(C)`: `M_9` on B0b, `M_6` on B4 |
| 2.5 | *"one qubit per cell, proved minimal. This is not zero addition."* | **K1_SCOPED** | [X] | minimality rests on `5` prime and `5 ∤ 9`, `5 ∤ 7`. On B0b `9 \| 81`, so `M_9^{⊗N}` is a directed system with unital embeddings and **no adjunction at all**; on B4 `C^6` factorises outright |
| 2.6 | `Phi_{N+1} = Phi_N . Z_{N+1}`; firing is absorbing | **CARRIER_INDEPENDENT** | [X] | `Omega_N` to `N = 10^6` on both four-class carriers |
| 2.7 | `sup_k \|Z_k\| = 1`, *"the carrier's recurrence is untouched"* | **CARRIER_INDEPENDENT in form; CONNECTION_SCOPED in ATTAINED-vs-APPROACHED** | [X] | see 7.2 |
| 2.8 | the monotonicity theorem is *"a tautology of any per-cell tensor system"* | **CARRIER_INDEPENDENT** | [T] | `\|Z_k\| <= 1` by the triangle inequality on non-negative weights summing to 1. No incidence input |
| 2.9 | **the nine, P-1 … P-9** | **UNDETERMINED — seven of nine; two are now decided** | [A]/[X] | lane D's blanket *"requires S3's code"* is **over-broad** (refuter D-2). **P-9 is decided**: S3 §5.9's closed form `\|omega_F^N - omega_C^N\| = 2 sqrt(1-\|Omega_N\|^2)` needs only `\|Omega_N\|`, and reaches `2.000000000000` on **both** four-class carriers → **CARRIER_INDEPENDENT [X]**. **P-3 is decided**: it *is* W-08's monotonicity → **CARRIER_INDEPENDENT [T]**. P-1 is a statement about all finite constructions. **P-2, P-4, P-5, P-6, P-7, P-8 remain UNDETERMINED and unrunnable — S3 published no lane directory** |
| 2.10 | `lambda = -0.766802` | already **ERRATUM'd** | — | untouched |
| **2.11 †NEW ROW** | *"repeated circuits span 3 dimensions at N=1 and 3 at N=100 — the trap, disarmed by computation"* (`REGISTER:104-105`, S3:264-271 *"EXACTLY the three-way split algebra and nothing more"*) | **the "no growth" half CARRIER_INDEPENDENT; the clause "EXACTLY the class algebra and nothing more" is FALSE at four classes** | [X] | refuter D-2: span of `{M_dF^n, M_c^n}` is flat in `N` on all seven carriers, but the value is the occupied-class count, and on B0b/B4 the **pure powers span 3 while the algebra they GENERATE is 4** — because `(e10+e11)+(e00+e01) = (e01+e11)+(e00+e10)` is a relation that exists only when class 00 is occupied |
| **2.12 †NEW ROW** | *"the record slot is non-abelian by necessity"* (`REGISTER:101-102`) | **UNDETERMINED** | [A] | an argument about pure states of abelian algebras; never run at four classes by anyone, including this round |
| **2.13 †NEW ROW** | P-9's finite-stage exception clause (*"on W-01's firing locus `Omega_N = 0` exactly for all N >= n"*) | **VACUOUS on both four-class carriers** | [X] | neither has a torus zero at its own SENSE-U weights, so the clause has no instance there. `min\|Z_1\|` = `1/9` (B0b) and `sqrt(6)/9` (B4), both exact |

### ERRATUM AGAINST W-02 — the rate and the resonance

| # | claim | scope | basis | evidence / correction |
|---|---|---|---|---|
| E2.1 | `f=2.0, c=1.1` is exactly resonant (`-11f + 20c = 0`); its orbit is on a subtorus; `-0.767014993` vs generic `-0.767507880` | **CARRIER_INDEPENDENT (arithmetic of the connection)** | [T] | confirmed in exact rationals here (primitive relation `(11,-20)`), and reproduced independently by C-1 in an exact rational series and by C-2 by roots at dps 40 |
| **E2.2 †NEW ROW** | **the typing rule the corpus needs and does not have** | **CARRIER_INDEPENDENT — and it convicts three more published connections** | [T] | **if `f` and `c` are both rational then `mf + nc` is rational and `2 pi j` is irrational for `j != 0`, so `j = 0` and a nonzero relation always exists: EVERY RATIONAL CONNECTION IS EXACTLY RESONANT.** Exact census (my Part C): `f=2.0,c=1.1` → `(11,-20)`; `f=2.0,c=2.0` → `(1,-1)`; **`S4:973`'s `3.14159, 1.57080`, published in terms as *"they are generic"*, → `(157080, -314159)`**; lane D's hard-coded `1.3, 2.0` → `(20,-13)`; lane D's "golden" pair → `(1,1)` since `1/phi + 1/phi^2 = 1` exactly. **`S4:603`'s `f = 1.0, c = sqrt(2)` is the ONLY generic connection the corpus publishes — and it is the one S4 used to verify the entire `lambda` column.** S4:973's mislabelling sits **inside the paragraph correcting the corpus's first two mislabelled connections** |

### W-03 — the measurement and the multiset theorem

| # | claim | scope | basis | evidence / correction |
|---|---|---|---|---|
| 3.1 | **the multiset theorem, "24 of 24 permutations invariant"** | **CARRIER_INDEPENDENT at real non-negative weights — i.e. on every object the corpus contains** | [T] | proved pointwise from `\|a + b e^{it}\| = \|b + a e^{it}\|` for real non-negative `a,b`, plus the two Jensen pairings. **Its extension off the non-negative locus is UNDETERMINED and six names have failed** (§0.3.5). **The registered figure `2.4e-15` has NO ARTIFACT IN THE CORPUS** — verified by me: it occurs only in the register and as a print literal quoting the register; the sealed four-class runs measure `9.285e-06` and `3.775e-15`. **Pointer-rule breach; the figure is UNDETERMINED** |
| 3.2 | **the LEVEL at which it holds** (not registered anywhere) | **CARRIER_INDEPENDENT, OBSERVABLE_SCOPED. NEW and it survives.** | [T] | **2 of 24** permutations preserve `\|Z_k\|` pointwise (identity + the involution); **24 of 24** preserve `lambda`. The labels are invisible **asymptotically** and **visible at every finite `N`**. **Lane D's novelty claim falls** (`REGISTER:196-200` states both halves), the **substance does not**: what is new is the sharpness — no third permutation survives |
| 3.3 | the involution `00 <-> 11`, `10 <-> 01` | **CARRIER_INDEPENDENT, and stronger than 3.1** | [T] | fixes `\|Z_k\|` pointwise in `k`; the only non-identity permutation that does |
| 3.4 | *"the carrier's topology is inert"* | **CARRIER_INDEPENDENT — and UNTESTABLE BY CONSTRUCTION** | [C] | `d2` enters the functional nowhere. Lane D declared its own exhibit a zero-variable control; **refuter D-2 found lane D committed the same void anyway in `w10a`-style form elsewhere in the round** (see 5.1 note). The row is correct and carries no evidence |
| 3.5 | *"formation does not see incidence"* | **CARRIER_INDEPENDENT, with a wording correction** | [T] | the class **multiset** is exactly the incidence content that survives; what is invisible is the **labelling**. **And the correction has a limit**: under COR-F's transport the labels are visible again (§0.1(b)) |
| 3.6 | B0b `= log(4/9)` exactly; *"nine of nine rates are exact"* | **CONFIRMED and STRENGTHENED** | [T] | branch domination gives closed forms on both four-class carriers. Verified at `mp.dps = 30` with panels split at the exact crossings: B0b `log(4/9)`, B4 `log(1/2)`, SENSE C `log(1/4)`, quadrature-minus-closed-form **exactly 0** (A-1). **S4:599's "QUADRATURE ONLY" label was wrong.** *Novelty correction:* the branch-dominance argument is **S4:594-595's own** (its `0.2222 + 0.1111 cos y` **is** the exact `C = -2/9`, `D = -1/9`); what is new is the general `\|D\| <= \|C\|` form |
| 3.7 | the charge run (S4-1 fails at charge ≠ 1) and the SU(2) run | **UNDETERMINED w.r.t. carrier** | [A] | orthogonal axes. **And lane D's row 5.x restating Theorem S4-1 as CARRIER_INDEPENDENT/HIGH is MISMARKED** — `REGISTER:210-212` records that theorem **failing** under charge; the mark must read **CARRIER_INDEPENDENT AND CHARGE-1-SCOPED** |
| 3.8 | *"the carrier axis is a list, not a family"* | **CARRIER_INDEPENDENT (a fact about S4)** | [X] | 2 of 11 named carriers occupy four classes — **as designated by S4** (see §0.4.2) |
| **3.9 †NEW ROW** | the **exceptional-set structure** — `lambda_B` constant on a full-measure set, 19 distinct values, `527/314/213`, *"discontinuity everywhere"*, *"neither upper- nor lower-semicontinuous"* | **STRONGLY WEIGHT-SCOPED, and no register row states that** | [X] | refuter D-2 derived the **four-class** resonant closed form `lambda_res(m,n) = m(p01 + p00 z^m + p11 z^n + p10 z^{m+n})` — S4 publishes only the three-term form because K1 has `p00 = 0` — validated against S4's 16 published values to `4.8e-10`. At the corpus's own headline resonance `(11,20)` the departure is `+4.93e-04` on K1, `+1.17e-08` on B0b, `-4.84e-10` on B4; by `(41,53)` the four-class columns are at `1e-14` while K1 is at `1e-5`. **Mechanism: torus zeros (lane D's own D-04), never connected to this block** |
| **3.10 †NEW ROW** | `(pi,pi)` is a **third strict saddle** (Morse/landscape correction) | **CARRIER-DEPENDENT** | [X] | refuter D-2, with `sum (-1)^index = chi(T^2) = 0` as an unfudgeable check: every three-class arm gives (1 max, 3 saddles, 2 conical zeros); **B0b gives (1 max, 2 saddles, 1 min) — four critical points, and both conical zeros, which S4 calls *"this is W-01's convex-hull criterion"*, are GONE** |
| **3.11 †NEW ROW** | §2's *"every vertex phase of `s` cancels"* | **CONVENTION_SCOPED — its hypothesis is the transport convention, not the occupancy** | [C] | exhibited at `1.4e-16` over 200 random complex states at four classes and **declared a zero-variable control by the refuter who ran it**: it could not have failed under fibre-wise scalar transport, which is the point |

### W-04 — the errata

| # | claim | scope | basis | evidence / correction |
|---|---|---|---|---|
| 4.1 | **ERR-1** — the operative variable is SCALARITY, not commutativity | **CARRIER_INDEPENDENT (fibre-scoped)** | [T] | reproduced; no carrier or class occupancy enters |
| 4.2 | **ERR-2** — the SU(2) run was a three-way confound | **CARRIER_INDEPENDENT (a fact about the control)** | [T] | an argument about what varied |
| 4.3 | **ERR-3** — the gravitational term was dropped, unrecorded | **CARRIER_INDEPENDENT (corpus fact), re-verified** | [X] | `coupling constant` 0 artifacts; `backreaction`, `edge mode`, `plaquette` 1 each, all inside the errata |
| 4.4 | **ERR-4** — S1 is unaudited; *"smallest complex"* is false | **K1_SCOPED and CONFIRMED** | [X] | and **the counterexample carrier is itself three-class**, `{11:2, 10:1, 01:1}` |
| 4.5 | *"what K1 was: inert"* + four unledgered choices | **(i)–(iii) CARRIER_INDEPENDENT; (iv) primality K1_SCOPED** | [X] | and **(ii), loop transport as scalar multiplication, is the convention of §0.1(b) — it is the choice this table's headline is conditional on** |
| 4.6 | the residue (gauge non-factorisation ⇒ edge modes) | **UNDETERMINED** | [A] | W-06's IMP-2 records it leans on a theorem the register refuted two rows earlier |

### W-05 — the survivor list and the STOP legs

| # | claim | scope | basis | evidence / correction |
|---|---|---|---|---|
| 5.1 | **N1 — `lambda = m(P)`** | **CARRIER_INDEPENDENT** | [T] + [X] | nothing in the identification uses `p00 = 0`, three classes, or K1's incidence. Direct schedule-B simulation on B0b and B4 agrees with `m(P)` to `2.3e-07` / `3.6e-07`. **Publication conditions, all three new and all three required:** (i) `L = {0}` **plus** either (H2a) the strict zero-free inequality `hi.lo > 0`, under which Weyl alone suffices with **no Diophantine input**, or (H2b) an inhomogeneous Diophantine condition (lane C, confirmed by both its refuters — (H2a) is the genuinely new clause); (ii) Lawton 1983 is **not** the missing theorem, **Lind–Schmidt–Ward 1990 is**; (iii) **`S4:582`'s SENSE-C row and lane A's `w10a_4_lambda` rows are a ZERO-VARIABLE CONTROL**: `B0b C` and `B4 C` are the same byte-identical `pi = (1/4,1/4,1/4,1/4)` (input sha256 identical), printed as two carriers reproducing a `lambda` column. **8 printed rows, 6 independent data** (refuter A-2) |
| 5.2 | **N2 — the multiset theorem** | **CARRIER_INDEPENDENT with 3.1, 3.2, 3.3 attached** | [T] | **must not be published without the level qualification (3.2) and without the note that its registered figure has no artifact (3.1)** |
| 5.3 | **N3 — the Haar-null inversion** | **CARRIER_INDEPENDENT, given a two-loop designation** | [X] | firing/resonance measure is a statement about there being exactly two designated holonomies. **But lane D's supporting rows are noise**: at `(11,20)` the true departure on B0b is `+1.17e-08` where lane D printed `9.66e-07`, and at `c = -f` the true departure is **exactly zero** where lane D printed `5.80e-07`. Two of its four evidence rows are measurement noise; the two low-height resonances carry the row |
| 5.4 | **N4 — fibre-wise-ness** | **CARRIER_INDEPENDENT** | [T] | `dim A^G = V` for the fibre-wise group, settled **exactly** by solving the fixed-point equations at nine `(V,r)` pairs (refuter D-1) rather than by the Monte-Carlo projector lane D had already had to self-correct once |
| 5.5 | LEG ONE — `M_gamma` lies in the gauge group | **CARRIER_INDEPENDENT, and half-scoped** | [T] | fibre-wise on any carrier — a definition. Under the **full** gauge law it is a gauge transformation only if `theta` is constant on every edge |
| 5.6 | LEG TWO — the slot is already inside the carrier, `dim_C = 4` | **CARRIER_INDEPENDENT, conditioned on `G != {1}`** | [X] | on all five carriers, under both the correct branch operator and lane D's incorrect one |
| 5.7 | LEG THREE — the rejections rest on `5 ∤ 9`, `5 ∤ 7` | **K1_SCOPED, CONFIRMED** | [X] | |
| 5.8 | *"has a field: NO"* | **CARRIER_INDEPENDENT structurally; UNDETERMINED experimentally** | [A] | no action exists anywhere in the corpus; W-05's experimental headline is already withdrawn as confounded |
| 5.9 | the WITHDRAWN list X1–X5 | **CARRIER_INDEPENDENT** | — | withdrawals, not claims |
| 5.10 | the rediscovery ledger | **CARRIER_INDEPENDENT (bibliography)** | [A] | corrected once by W-06; not re-adjudicated |

### W-06 — the audit of the audits

| # | claim | scope | basis | evidence / correction |
|---|---|---|---|---|
| 6.1 | the dressed restoration — S3-0's broad form is FALSE, not vacuous | **MECHANISM CARRIER_INDEPENDENT; NUMBERS K1+state+connection scoped** | [X] | exhibited on B0b (lane D) **and now on B4** (refuter D-1: `V` 9→6, `E` 18→8, `chi` 0→2, different tree, different multiset; gauge invariance `1.11e-16`, diagonal separation `1.11e-16`, dressed separation `0.352750896307`, closed form to `1.112e-16` on all 36 pairs). **Two arms now, not one** |
| 6.2 | `4.45e-16`, `3sqrt3/10`, `2.221e-16`, `1000 of 4000`, the Bell `R_N` transplant, the Schmidt `[1,0]`, the 52-subgroup sweep | **UNDETERMINED BY CUSTODY** | [A] | W-06 has no artifact and no lane code. **The Bell transplant and the Schmidt `[1,0]` are load-bearing and have never been reproduced by anyone** |
| 6.3 | the wedge-growth route `V = 4k+1` | **K1_SCOPED** | [X] | a property of 5. B4 is not in the sequence; B0b is, but `9 \| 81` dissolves the motivation |
| 6.4 | `U(1)^V` is the 0-cochain group (leg a); the 52-sweep (leg b) | **leg (a) CARRIER_INDEPENDENT; leg (b) VOID as a control** | [T]/[C] | as W-07 ruled |
| 6.5 | IMP-1 — *"could not have failed"* voids a control, never a theorem | **CARRIER_INDEPENDENT (methodological)** | [T] | **applied throughout this table via the BASIS column, which is what it asks for** |
| 6.6 | *"the crossing contains no measurement"* | **UNDETERMINED** | [A] | no code; never reproduced |
| **6.7 †NEW ROW — THE LARGEST OMISSION IN LANE D's TABLE** | **COR-F / the edge-tick transport** — *"the thing that decided the spine was never imported… `T` with `T^3 = M_gamma`"* (`REGISTER:527-542`) | **CARRIER_INDEPENDENT AND UNSCOPED — it has no register row, and every CARRIER_INDEPENDENT verdict above is conditional on rejecting it** | [X] | built on B0b by refuter D-2: `\|\|T_F^4 - M_dF\|\| = 1.24e-16`, `T` non-diagonal, does not preserve `diag`, **diagonal-observable branch separation `0.244292` against `2.78e-17` under the scalar convention**, and `\|\|[T_F,T_C]\|\|_F` = `2.828427` (B0b), `2.449490` (K1), **`0.000000` (B0a)** — non-zero **exactly when class 11 is occupied**. Lane D's only nominal cover, its row 4.5, cites a leg that computes **vertex-count arithmetic and contains no transport** |
| **6.8 †NEW ROW** | the twenty sealed corrections `S3 COR-A…COR-L`, `S2 COR-A…COR-H` | **UNSCOPED — zero coverage** | — | `grep -c "COR-"` on the scope table returns **0**, and no file in the lane contains the string. **Four of them were on the CARRY list of the brief.** COR-B, COR-E, COR-F, COR-K are all live and all bear on rows above |

### W-07 and the ERRATUM against it

| # | claim | scope | basis | evidence / correction |
|---|---|---|---|---|
| 7.1 | `ord(rho)` as the operative variable | **FALLEN** (W-08) | — | not re-litigated |
| 7.2 | **ATTAINED vs APPROACHED** | **CARRIER_INDEPENDENT, CONNECTION_SCOPED** | [X] | identical figures on K1, B0b and **B4**: `1000/4000` at `W_F=-1,W_C=-i`; `2000/4000` at `W_F=W_C=-1`; `0/4000` at random and at irrational connections. Different `V`, `E`, `F`, tree, multiset — the cut is arithmetic of the connection |
| 7.3 | the dressed reconstruction reproducing W-06 | **K1_SCOPED (numbers)** | [T] | the closed form explains why: `\|A_uv\|` is a state-and-tree quantity |
| 7.4 | *"W-06 has no artifact and no lane code"* | **CARRIER_INDEPENDENT (custody fact)** | [X] | still true at the bytes |
| 7.5 | ERRATUM E-1 and E-2 | **CONNECTION_SCOPED / K1+state_SCOPED** | [T] | untouched, and E-1 is reinforced by E2.2 above |

### W-08 — the race

| # | claim | scope | basis | evidence / correction |
|---|---|---|---|---|
| 8.1 | `\|Z_k\| <= 1`, so `\|Omega_N\|` is monotone non-increasing | **CARRIER_INDEPENDENT** | [T] | triangle inequality, no incidence input. Its four-class run is a **corollary of 8.3 and is not independent evidence** (refuter A-1) |
| 8.2 | **the founding obstruction is FALSE AS AN INFERENCE** | **CARRIER_INDEPENDENT** | [T] | follows from 8.1 |
| 8.3 | the exact character identity | **CARRIER_INDEPENDENT and CLASS-COUNT-INDEPENDENT** | [T] | one-line proof for any number of classes: `1-\|Z\|^2 = \|sum w_j\|^2 - \|sum w_j chi_j\|^2 = sum_{j<l} w_j w_l \|chi_j - chi_l\|^2`. Verified as a **formal polynomial identity** for `n = 2..6` classes (A-2: 0 non-zero monomials) and in the **group ring** over 26 moduli including primes to 37 (A-1: 0 residuals). My Part E re-runs it exactly at four/three/two classes: **residual 0**, and is **DECLARED A CONTROL THAT COULD NOT HAVE FAILED** |
| 8.4 | the linear floor; `G = {1}` ⟺ no formation | **CARRIER_INDEPENDENT** | [T] | follows from 8.3 + Dirichlet, no Diophantine input |
| 8.5 | the measured decay densities `0.4919 / 0.4692 / 0.5295 / 0.4692` | **K1_SCOPED (measurements)** | [X] | reproduced independently at W-08's own `K = 1e7`: `0.491886` (order-4), `0.469183` (resonant) — refuter A-2 |
| 8.6 | **the schedule result** | **CARRIER_INDEPENDENT (qualitative)** | [X] | exhibited on B0b and B4 |
| 8.7 | **the four printed constants `0.606, 0.615, 0.588, 0.601` ("flat in K")** | **the K-INDEPENDENCE is CARRIER_INDEPENDENT [X]; the FOUR CONSTANTS remain UNDETERMINED** | [X]/[A] | lane D's second reading (*"the adversary is stronger than reported"*) is **killed**: its four sweep connections were **all degenerate** (`(1.3,2.0)`, `(2.0,1.1)`, the golden pair, order-4 — see E2.2). At genuinely generic connections refuter D-1 gets `0.5506–0.5642` (K1), `0.6953–0.7056` (B0b), `0.6717–0.6819` (B4), **flat across `K = 1e4..1e7`**. **But `0.55 != 0.606` on the same weights: the constants are still not reproduced, and D-1's withdrawal of the UNDETERMINED mark is a null scored as a confirmation (§0.3.3)** |
| 8.8 | the exponents `K^{-1/2}` at `d_eff = 2`, `K^{-1/3}` at `d_eff = 1` | **CARRIER_INDEPENDENT** | [T] | `d_eff` = rank of the relation lattice = 2 for every occupied set of size ≥ 3 **at charge 1**; at charge ≠ 1 the register's own counterexample gives rank 1 |
| 8.9 | the onset non-uniformity `K_0 ~ t^{-2/3}` | **UNDETERMINED** | [A] | never tested off K1 |
| 8.10 | *"the weakest load-bearing claim is `p00 = 0`"* | **SUPERSEDED by W-09**, and W-09 in turn qualified by §0.4.2 | [X] | |
| 8.11 | the isolation audit's finding (*"the commonest FATAL defect is ZERO variables moved"*) | **CARRIER_INDEPENDENT (methodological), and it fired FIVE more times in this round** | [X] | lane A's `w10a_4_lambda` SENSE-C pair (byte-identical `pi`, 8 rows/6 data); lane A's script-3 PART C2 arms (exact involution images, max dev `4.441e-16`); lane D's leg 3C (`arm1.tobytes() == arm2.tobytes()`); lane D's D-02 conjugation arms; refuter B-1's own B.4 sweep (self-flagged). **And the guard is mismatched to the defect**: lanes hash **input** vectors; output-diffing at `1e-12` turns "150 cases" into **56 distinct** and "120 cases" into **40 distinct** |
| **8.12 †NEW ROW** | **the FLOOR half** — *"the floor's entire contribution is `O(1)` and `K`-independent… THE FLOOR AND THE RATE ARE INDEPENDENT COORDINATES. THEY DO NOT RACE"* | **UNDETERMINED off K1** | [A] | half of W-08's own headline, and **no lane in this round ran it on a four-class carrier** |

### W-09

| # | claim | scope | basis | evidence / correction |
|---|---|---|---|---|
| 9.1 | firing region exactly `1/4` / exactly `1/2` | **REPRODUCED, now four times, on three algorithms and two sampling schemes** | [X] | mine is the deterministic-grid arm; `1442401/1442401` closed-form agreement |
| 9.2 | the operative variable is **all four classes occupied** | **CONFIRMED, and narrowed twice** | [X] | (i) what dies at four classes is **sign-sensitivity in each holonomy separately**, not role-distinction — the `f <-> c` exchange was always blind (a conjugation identity); (ii) *"needs a pinch AND a spectator"* is a statement about the **designation**, not the carrier (§0.4.2) |
| 9.3 | *"the corpus has never run a four-class carrier through anything"* | **TRUE WHEN WRITTEN; DISCHARGED BY THIS ROUND, with two caveats on the discharge** | [X] | (a) B4's row is **under-determined at source** (§0.3.1); (b) B0b's four-classness is **one of 16 reachable designations** on its own complex. The corpus has now run **two four-class DESIGNATIONS**, one of which is not recoverable from the row it was quoted from |

### PRECONDITIONS — †ALL NEW ROWS, and every CARRIER_INDEPENDENT verdict above is conditional on them

| # | precondition | scope | basis | evidence |
|---|---|---|---|---|
| P.1 | **cycle rank >= 2** — else `gamma_C` cannot be designated | **CARRIER_INDEPENDENT precondition, broken by S4's own B5** | [X] | S4:519 in terms: *"the formation datum does not exist on it — this is the one place in S4 where topology, and nothing else, decides an outcome."* B5 appears in no scope row |
| P.2 | **b1 >= 1** — else `gamma_C` bounds and there is no flat holonomy | **CARRIER_INDEPENDENT precondition, broken by S4's own B2** | [X] | B2's class multiset is **K1's exactly**, so it breaks W-01's advertised virtue for a reason having nothing to do with class occupancy |
| P.3 | **equal loop lengths** — else circuit clock ≠ edge clock, forever | **satisfied by every carrier ever run; broken on B0b** | [X] | see 1.10 |
| P.4 | **whole-circuit scalar transport** (`M_gamma`, not COR-F's `T`) | **the convention on which the entire CARRIER_INDEPENDENT column rests** | [X] | see 6.7 and §0.1(b) |
| P.5 | **the loop designation** — `gamma_F`, `gamma_C` chosen under S4 ledger C4 *"by fiat"* | **UNDETERMINED — no admissibility criterion exists** | [A] | 16 multisets reachable on B0b; forced on B4 and K1. See §0.4.2 |

---

## 3. CORRECTED WORDINGS OF RECORD

Where the register states a claim unqualified that needs a qualification, here is the exact
replacement text. **These are proposed corrections, not entered rulings** — entering them is a
registrar act and I am flagging, not patching.

**`REGISTER:43` (W-01's criterion).** For *"vanishes iff 0 lies in the convex hull of three
unit-modulus coefficients"* read:
> *`<M_dF s, M_c s>` vanishes **only if** 0 lies in the convex hull of the occupied characters, and
> the converse fails wherever the class weights are not realisable as the convex coefficients. On
> K1 exactly three characters are occupied because no vertex lies outside both loops; the firing
> region is then exactly `1/4` (Wendel 1962). Where all four classes are occupied `uv` is determined
> by `u` and `v`, Wendel does not apply, and the region is exactly `1/2`, with closed form
> `cos f + cos c <= 0`.*

**`REGISTER:50` (W-01's advertised virtue).** For *"it distinguishes curvature from flat holonomy,
which K1 exists to separate"* read:
> *It is sensitive to the **sign** of each holonomy separately — never to their **roles**, since the
> exchange `f <-> c` and the conjugation `c -> -c` are blind at every class count by an identity of
> complex conjugation. The sign-sensitivity is a **three-class** property: `f -> -f` flips the
> verdict on `~50%` of the parameter square at three occupied classes and on **exactly none** at
> four, because `cos f + cos c` is even in `f`.*

**`REGISTER:125` (W-02's support table).** Append:
> *Six two-element supports exist. The three above are the three K1 can realise. The three omitted —
> `{00,10} -> <u>`, `{00,01} -> <v>`, `{00,11} -> <uv>` — each require a vertex in **neither** loop;
> `{00,11}` requires one in **neither** and one in **both**, and is trivial exactly when
> `W_C = W_F`.*

**`REGISTER:196-197` (the multiset theorem).** For *"a function of the MULTISET of the four class
weights — 24 of 24 permutations invariant, worst spread `2.4e-15`"* read:
> *…a function of the multiset of the class weights, at the level of `lambda` and at no finer level:
> **24 of 24** permutations preserve `lambda`; **exactly 2 of 24** — the identity and the involution
> — preserve `|Z_k|` pointwise, so the labels are invisible asymptotically and visible at every
> finite `N`. The hypothesis of record is **real non-negative** coefficients, which every
> pushforward `pi` satisfies; the extension to complex coefficients is **open**, six proposed names
> have failed, and two independent mechanisms (flux coincidence and branch domination) are known to
> produce invariance. **The figure `2.4e-15` has no artifact in this corpus** and is withdrawn under
> the pointer rule; the sealed four-class runs on disk measure `9.285e-06` and `3.775e-15`.*

**`REGISTER:456-460` (N1).** Append to the publication line:
> *…and it is publishable only with (i) `L = {0}` plus either the strict zero-free inequality
> `hi.lo > 0` (under which Weyl suffices with no Diophantine input) or an inhomogeneous Diophantine
> condition; (ii) **Lind–Schmidt–Ward 1990**, not Lawton 1983, named as the governing theorem;
> (iii) the two topological preconditions (cycle rank `>= 2`, `b1 >= 1`); and (iv) the statement
> that its carrier-independence is a consequence of the whole-circuit scalar-transport convention
> and does not survive COR-F's edge-tick transport.*

**`REGISTER:620-624` (the wedge route).** Append:
> *`V = 4k+1` is a property of the number 5. The sequence's own motivation — S3's rejection of
> CHOICE LEDGER C1 on `5 ∤ 9`, `5 ∤ 7` — **evaporates at `V = 9`, where `9 | 81`** and the
> adjunction-free directed system exists outright. The cost line's *"untouched recurrence
> obstruction"* names an obstruction W-08 showed **false as an inference**; what survives is a
> schedule statement with an exponent and no admissibility criterion.*

**`REGISTER:818-823` (the schedule constants).** For *"accumulates `0.606, 0.615, 0.588, 0.601` nats
at `K = 10^4..10^7`"* read:
> *…accumulates `O(1)` nats, flat in `K` across three decades, at every generic connection and on
> every carrier tested (K1 `~0.56`, B0b `~0.70`, B4 `~0.68`). **The four constants as printed have
> not been reproduced by any later lane** and the connection or tie-break producing them is not
> recorded.*

**`REGISTER:920` (W-09's headline).** For *"K1 has only one of them, and so does every carrier the
corpus ever ran"* read:
> *…every carrier-**and-designation** the corpus ever ran. Class occupancy is a property of the
> designated loop pair, not of the complex: on B0b's own complex with S4's own `gamma_F`, 16
> distinct class multisets are reachable by admissible `gamma_C`. On K1, B1q and B4 the designation
> is forced; on B0b it is not, and S4 ledger C4 closed the choice by fiat.*

**PROPOSED ERRATUM AGAINST S4 (new).**
> *`S4:973`'s "the truncated decimals `3.14159` and `1.57080` are not `pi` and `pi/2`; **they are
> generic**" is false. Both are rational, so the connection is exactly resonant with primitive
> relation `(157080, -314159)`. **Every rational `(f,c)` is exactly resonant**, because `mf + nc` is
> then rational and `2 pi j` is irrational for `j != 0`. Of the connections the corpus publishes,
> exactly one is generic: `S4:603`'s `f = 1.0, c = sqrt(2)` — the connection S4 used to verify the
> entire `lambda` column. The numerical conclusions are unaffected to any attainable precision
> (Boyd–Lawton), but the typing is wrong, and it is wrong inside the paragraph that corrects the
> corpus's first two mislabelled connections.*

**PROPOSED CORRECTION AGAINST S4's B4 ROW (new, COR-K class).**
> *`S4:519, :580`'s B4 row is not recoverable from its own parameters. A second complex — a triangle
> with two 2-cells glued at two points to a pentagon with two 2-cells — reproduces every published
> column (`V=6 E=8 F=4 chi=2 b=(1,1,2)`, gauge 5, inv 3, curv 2, flat 1, `d1.d2 = 0`, `gamma_F`
> bounds, `gamma_C` does not, independent) with class multiset `{00:1,10:1,01:2,11:2}` and
> `lambda(U) = log(1/3)` against the published `log(1/2)`. Publish `d1` and `d2`, or the row's one
> downstream-relevant column is a free parameter.*

---

## 4. NEW DEFECTS OF RECORD FROM THIS ROUND

| # | severity | against | defect |
|---|---|---|---|
| N-1 | **high** | S4, and lane A's headline | S4's B4 row is under-determined by its own parameters; a second admissible spindle gives a different multiset and `lambda = log(1/3)`. Independently confirmed here (Part B). Lane A hard-coded S4's multiset as its build target, so its match **could not have failed** |
| N-2 | **high** | W-09, the round's framing, and this synthesis | class occupancy is a property of the **loop designation**, not the carrier; 16 multisets reachable on B0b; no admissibility criterion for designations exists |
| N-3 | **high** | the register | **COR-F / the edge-tick transport has no scope row and no register row**, and every CARRIER_INDEPENDENT verdict is conditional on rejecting it. `\|\|[T_F,T_C]\|\|` is non-zero **exactly when class 11 is occupied** |
| N-4 | **high** | S4:973, lane D, and the corpus's typing discipline | **every rational `(f,c)` is exactly resonant.** Three further published/used connections mislabelled "generic"; lane D hard-coded one in **eight legs** |
| N-5 | **medium** | the register (pointer rule) | `2.4e-15` at `REGISTER:197` has no artifact in the corpus. So do W-06's `3sqrt3/10`-adjacent figures, and W-08's four schedule constants are not reproduced |
| N-6 | **medium** | five lanes of this round | **five more zero-variable controls**, including two inside blocks written to correct earlier zero-variable controls, and an arms-diff guard that hashes **inputs** while the collapses are in the **outputs** (150 cases → 56 distinct; 120 → 40) |
| N-7 | **medium** | six consecutive layers | the multiset theorem's hypothesis has now been misnamed six times, three of them **this round, by refuters commissioned to catch misnaming**. Two mechanisms exist; no name covers both. **UNDETERMINED, and I decline to add a seventh** |
| N-8 | **medium** | lane D's method statement | zero of twenty sealed corrections named anywhere in the lane; `grep -c "COR-"` = 0, re-verified |
| N-9 | **medium** | the round's design | no mechanism exists for one lane to inherit another lane's correction; the same defect was found and lost in parallel at least four times |
| N-10 | **low-medium** | lane A, lane C, lane D (novelty) | four novelty overclaims contradicted by the corpus's own bytes or the lanes' own scripts: branch dominance is S4:594-595's argument; the SENSE-C factorisation is S4:598's; three of four "never-audited" rows were reproduced in `LANE_R_MAPS_REFUTER` in the W-03 round; `\|S\|=1 => never` is `REGISTER:126` verbatim |
| N-11 | **low** | lane C, latent | `C_04`'s phase engine forms `k * A_num` in int64 and reaches 88% of the ceiling at `K = 1e7`; **it wraps silently at `K ~ 1.13e7`**, which is the first thing any successor raising `K` will do |
| N-12 | **mine, recorded not patched** | this lane | (i) my Part A four-class `f -> -f` column **could not have failed** and is printed as an identity, disclosed **before** the number; (ii) my Part E is a control that could not have failed and is declared as one; (iii) my first Jensen evaluation of `B4-TRIPENT` returned `-1.098611627631` against the exact `log(1/3) = -1.098612288668` — a `6.6e-07` quadrature error caused by the log singularity the factorisation creates. I left both figures on the page rather than replacing the quadrature number, because it is the round's own quadrature warning firing on the round's own synthesis; (iv) my "independent" column first used the wrong test (image of `d2` rather than rank of the two cycle vectors) and returned `False` on both spindles — caught because it contradicted S4's published column, i.e. by the artifact and not by my reading |

---

## 5. WHAT THIS SYNTHESIS DOES NOT ESTABLISH — STATED BEFORE THE LEDGER, NOT AFTER

1. **It does not decide the disposition, the Mahler note, or the rename.**
2. **Thirteen rows carry an UNDETERMINED mark — one of them covering six of the nine crossing
   properties — and none was guessed.** The largest blocks are those six crossing properties
   (unrunnable — S3 published no lane directory), everything W-06 computed (no artifact), W-08's
   floor half off K1, and the multiset theorem's hypothesis off the non-negative locus.
3. **S1 and S4 remain unaudited, and this round made that worse, not better.** The four-class
   evidence rests on two designations of two complexes, one of which (B4) **is not recoverable from
   the row it was quoted from** and the other of which (B0b) is one of 16 admissible choices.
   **S4's carrier code is not in the repo and no adversary has ever read it.**
4. **The central CARRIER_INDEPENDENT column is conditional on a convention (P.4) that the corpus's
   own sealed audit flagged as unledgered, and that no register row scopes.** If the principal wants
   one number from this page it is that one.
5. **This is layer twelve of one lineage.** Every lane, every refuter and this synthesis are Opus 5.
   Six consecutive layers have been caught by the next and the rate has not fallen. **Discount this
   one.**

---

## 6. THE ONE-LINE LEDGER

> **The corpus does contain carrier-independent results — the whole formation-functional layer, and
> they are theorems, not measurements. But their carrier-independence is a consequence of the
> transport convention W-06 named as having decided the spine, and it holds under two topological
> preconditions the register never states. What is three-class-scoped is the corpus's most-quoted
> sentence and the table under it. What is K1-scoped is the entire C\*-algebraic superstructure,
> whose hypothesis is the primality of 5 — which is also the hypothesis of the rebuild route. And
> the four-class evidence this round was commissioned to produce rests on two loop DESIGNATIONS, one
> of which is not recoverable from the S4 row it was quoted from.**
