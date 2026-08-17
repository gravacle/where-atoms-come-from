# W-08 — SYNTHESIS: **THE RACE IS DECIDED. THE DECAY OUTRUNS THE FLOOR, AND THE TWO NEVER COMPETE. W-07 SURVIVES ITS HEADLINE FIGURE AND LOSES ITS HEADLINE SENTENCE.**

Material for register entry W-08. Four build lanes (M1 IDENTIFICATION, M2 FLOOR, M3 ZEROSET,
M4 REFUTE-W07), eight refuters, one isolation audit over 64 comparisons. This page is the
synthesis and it commits. **It takes none of the principal's three decisions.**

Synthesis lane code, conventions, outputs and seals: `LANE_W08_SYNTHESIS/`, `SEALS.sha256`,
20 files, all verify. Ten scripts, `S_01`…`S_10`. Double precision is the default and is said
so; every precision-sensitive claim below is re-checked in exact rationals, exact Gaussian
rationals, or 60-digit decimal, and the check is named where it is used.

---

## 0. LINEAGE, DECLARED BEFORE ANY VERDICT

W-07 is Claude Opus 5. M1–M4 are Opus 5. All eight refuters are Opus 5. **This synthesis is
Opus 5 and is layer NINE of a single lineage.** Per `feedback_read_sealed_corrections_first`,
same-lineage adversaries share failure modes, and this corpus has now watched that happen three
times in one sub-chain: W-07 named `ord(rho)`, M4 named `min_k ||k theta||`, M4's refuter 2
named `(q, delta)` and then declared that its own name is stated at fixed `K` and fixed
`tol/amp` and so names the effect at a point of a grid rather than across it. **Discount this
page as one block with W-07 and M1–M4, not as an independent ninth check.** The
lineage-independent lane W-03 specified still does not exist for anything in W-08.

My own likeliest failure mode is the corpus's: I name an operative variable below (§1, §5) and
nobody has attacked it.

---

## 1. THE RACE — **ANSWERED: THE DECAY OUTRUNS THE FLOOR, ON EVERY CONNECTION OFF `G = {1}`, WITH NO DIOPHANTINE HYPOTHESIS. AND THE INTERESTING PART IS THAT THEY ARE NOT RACING.**

W-07's REOPENS clause asks: *does the record's decay rate outrun the near-return floor?*
**Yes.** Not narrowly, and not by a Diophantine accident. Three steps, none of which needs a
hypothesis on the arithmetic of the connection.

**(R1) THE RECORD CANNOT RECOVER, AND THE FOUNDING OBSTRUCTION'S INFERENCE IS DRAWN FROM THE
WRONG OBJECT.** `|Z_k| = |p11 chi_11^k + p10 chi_10^k + p01 chi_01^k| <= p11+p10+p01 = 1` by the
triangle inequality on non-negative weights summing to one. Hence `|Omega_N| = prod_{k<=N}|Z_k|`
is **monotone non-increasing**. Checked rather than asserted (`S_03_schedule_and_monotone.OUT.txt`,
S3a): over `k <= 10^6`, `max_k(|Z_k| - 1) = -1.446e-09` at `f=2.0, c=1.1` and `-1.591e-06` on the
Schmidt pair, with **0** exceedances in either.

A near-return at cell `k` means **that circuit writes nothing**. It does not un-write the previous
`N-1`. The founding obstruction's load-bearing inference — *"It fires, and then it un-fires. A
reversible write, not a record"* (`S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md:447`; repeated at
`S3_THE_CROSSING_V001.md:429`) — is drawn from the **single-cell** observable `|Z_k|` and never
from the product `Omega_N`, which is the object durability is a property of. **That is the
defect, and it is a defect of the OBSERVABLE, not of the measurement point.** W-07 looked for the
fault in *where* the obstruction was measured; the fault is in *what* was measured. The corpus
already half-owned this: S3-audit **COR-C** struck `|Omega_N| <= e^{lambda N}` and replaced it
with *"the monotonicity and `(1/N)log|Omega_N| -> lambda` that §4.3–4.4 actually prove"* — the
monotonicity was in the record and was never used against the recurrence sentence.

**(R2) THE DECAY IS LINEAR IN `K`, WITH AN EXPLICIT CONSTANT, UNCONDITIONALLY.** M2-1's identity,
re-verified here from the characters with no import from M2 (`S_06_theorem_checks.OUT.txt`, S6a):

```
|Z_k|^2  =  1  -  sum_{j<l} w_j w_l | chi_j^k - chi_l^k |^2
worst residual over 3 states x 6 root-of-unity orders x q^3 character triples x k<=6 : 4.441e-16
```

(both M2 refuters obtained **exactly 0** in `Fraction` arithmetic over 768 and 336 cases; the
identity is algebraic, not numerical). Hence for any pair with `chi_j != chi_l`,
`SUM_{k<=K}(1-|Z_k|) >= w_j w_l (K - 1/|sin(tau/2)|)` — **linear in `K`**, constant `0.12` on
RS-G (`= w11 w10`). Measured densities at `K = 10^7`, RS-G, one connection moving and nothing
else (`S_02_the_race.OUT.txt`): `0.491886` (order-4 ATTAINED), `0.469183` (resonant),
`0.529471` (W07GEN, rank 1), `0.469189` (Schmidt), `0.469189` (random seed 20260816). Every
one exceeds the proved `0.12`; none is near zero.

**(R3) THE FLOOR'S TOTAL CONTRIBUTION IS A FIXED, MINUTE SHARE OF THE DECAY BUDGET, INDEPENDENT
OF `K`.** Schmidt connection, RS-G, `K = 10^7` (`S_03…OUT.txt`, S3c). Total budget
`SUM(1-|Z_k|) = 4691887.4124` nats. The **entire** contribution of **every** cell in the
near-return band:

| threshold `eps` | # cells | `#/K` | M2's unfitted `0.838820*eps` | band SUM | **band share of budget** |
|---|---|---|---|---|---|
| `1e-1` | 861917 | 8.619e-02 | 8.388e-02 | 43491.06 | **9.269e-03** |
| `1e-3` | 8393 | 8.393e-04 | 8.388e-04 | 4.1986 | **8.949e-07** |
| `1e-6` | 8 | 8.000e-07 | 8.388e-07 | 0.0000 | **9.168e-13** |

The share scales like `0.8388 eps^2 / 2c` and **does not depend on `K` at all**. The `#/K` column
reproduces M2-4's `N(eps)/K = 0.838820*eps` — a constant M2 derived from `det M = 0.036` **with
no fitting** — to three digits at every threshold across six decades, from an implementation
sharing no code with M2. That is the best-supported quantitative law in W-08 and I add a
sixth-decade confirmation of it.

**THE SHARPEST FORM OF THE ANSWER IS IN M2's OWN LEG D4, WHICH THE LANE DID NOT READ AS AN ANSWER
TO THIS QUESTION** (`LANE_W08_M2_FLOOR/m2_d_attained.OUT.txt`, D4). Approaching S1's order-4 point
along `(alpha,beta) = (1/2 + t alpha*, 3/4 + t beta*)`, as `t` runs `1e-1 -> 1e-6` the floor
`F(1e6)` falls `2.189e-06 -> 3.226e-11` — **five orders deeper** — while the asymptotic rate does
**not move at all**: it is `-0.767507880` in every single row.

> **THE FLOOR AND THE RATE ARE INDEPENDENT COORDINATES. THEY DO NOT RACE. The floor can be driven
> arbitrarily deep without moving the rate by anything, because the floor is a property of ONE
> cell and the rate is a property of ALL of them.** Asking whether the decay outruns the floor is
> asking whether a linear quantity outruns a per-term one; the answer is yes, and the reason it is
> yes is that the question compares incommensurables.

**WHAT ACTUALLY DEFEATS DURABILITY — TWO THINGS, AND RECURRENCE IS NEITHER.**

1. **`G = {1}`.** M1-T4 / M2-2: `|Omega_N| -> 0` **iff** `G = <chi_a/chi_b : a,b in supp(pi)> != {1}`.
   Proved by Weyl on the **continuous** function `1-|P|` plus the strict triangle inequality; no
   Diophantine input. This is S3-audit **COR-B**'s four non-forming families, re-derived as a
   criterion. It survived all four refuters who attacked it, including the two degenerate
   sub-cases M1 self-flagged (M1-R1's `R1_07(A)`: all three cases have `G != {1}`, densities
   `2.7e-1..3.3e-1`).
2. **THE SCHEDULE — and this is where the founding obstruction genuinely survives.** Reproduced
   independently (`S_03…OUT.txt`, S3b): on the Schmidt connection, an adversary writing only the
   `K^{-1/2}·K = sqrt(K)` cells of smallest `1-|Z_k|` writes `J = 100, 316, 1000, 3162` cells at
   `K = 1e4..1e7` and accumulates `0.6008, 0.5965, 0.5967, 0.5961` nats: `|Omega| ~ 0.55` forever,
   **with unboundedly many writes**. M2-refuter-2's figures (`0.607/0.596/0.597/0.596`) reproduce
   to three digits. Under the corpus's **own** registered mechanism — W-02's divergence of
   `SUM(1-z_n)` — durability fails on an APPROACHED connection too. **Durability is a property of
   the (connection, schedule) pair. The corpus has never stated a schedule stipulation and M2
   self-flagged that its schedule half has no theorem.**

**I SUPPLY THE MISSING SCHEDULE THEOREM, DERIVED FROM M2-4's OWN `N(eps)` LAW AND CONFIRMED**
(`S_09_adversary_density.OUT.txt`). Let `J_max(K)` be the largest number of cells an adversary can
write out of `k <= K` while keeping the whole record above `|Omega| >= 1/e`. From
`N(eps) ~ C_d K eps^{d/2}`, the `J`-th smallest gap is `~ (J/(C_d K))^{2/d}` and the budget is
`~ (d/(d+2)) J (J/(C_d K))^{2/d}`, giving at fixed budget

```
  d_eff = 2 :  J_max ~ K^{1/2},  admissible write density ~ K^{-1/2}
  d_eff = 1 :  J_max ~ K^{2/3},  admissible write density ~ K^{-1/3}
```

Measured, RS-G, budget 1 nat, one variable (`d_eff`) moving:

| connection | `J(1e4)` | `J(1e5)` | `J(1e6)` | `J(1e7)` | last-decade exponent | theory |
|---|---|---|---|---|---|---|
| `d_eff=2` Schmidt | 129 | 409 | 1294 | 4095 | **0.5003** | 0.500 |
| `d_eff=2` random seed 20260816 | 129 | 409 | 1294 | 4095 | **0.5003** | 0.500 |
| `d_eff=1` resonant `f=2.0,c=1.1` | 129 | 377 | 1598 | 7416 | **0.6666** | 0.667 |
| `d_eff=1` W07GEN `2pi.phi, 2pi.phi^2` | 465 | 2163 | 10043 | 46619 | **0.6667** | 0.667 |

**A RANK-1 CONNECTION IS MORE EXPOSED TO THE ADVERSARY THAN A RANK-2 ONE, NOT LESS** — and both
of the corpus's distinguished connections, plus W-07 leg D's row labelled "GENERIC (badly
approximable)", are rank 1. And the attained case is the limit of this: on S1's order-4
connection `1-|Z_k| = 0` **exactly** on `250000` of `10^6` cells, so `J_max = infinity` at **zero**
budget. *That, and only that, is what ATTAINED buys the adversary: a **fixed positive** write
density instead of one that must fall like `K^{-1/2}` or `K^{-1/3}`.* W-07's ATTAINED/APPROACHED
cut is real, it is a cut on the **schedule** axis and nowhere else, and it is a difference of
**exponent**, not the "absolute versus cosmetic" dichotomy M2-7 stated and M2-refuter-2 correctly
attacked as a yardstick artefact.

**THE ONE NON-UNIFORMITY, WHICH IS MINE AND IS NOT IN THE RECORD** (`S_07_onset.OUT.txt`).
`lambda < 0` at every non-trivial connection, but the ONSET diverges. On the rank-1 locus
`f = c = t`, the number of circuits needed to write **one nat** is `2, 10, 47, 216, 1000, 4642` as
`t` runs `1 -> 1e-5` (`K_0 ~ t^{-2/3}`), while `lambda = log 0.3 = -1.203972804` in every row; on a
rank-2 pair `(t, t sqrt2)` it is `2, 9, 41, 189, 879` with `lambda = m(P)` in every row.
**Durability holds for every non-trivial connection and holds uniformly for none.** The
obstruction that actually survives on K1 is not recurrence and not finiteness of the spectrum; it
is that the write rate is not bounded below over the connection space.

**COMMITTED ANSWER, ONE SENTENCE.** *Off the `G = {1}` locus and on the honest schedule
`k_n = n`, the record's decay outruns the near-return floor by a linear-versus-`O(1)` margin that
requires no Diophantine hypothesis and is unchanged by driving the floor five orders deeper; the
founding obstruction "a finite discrete spectrum is recurrent, therefore a reversible write and
not a durable record" is **FALSE as an inference**, because it reads a single-cell recurrence as a
recurrence of the record, and `|Omega_N|` is monotone; what survives of it is a **schedule**
statement with an exponent, `K^{-1/2}` at `d_eff = 2` and `K^{-1/3}` at `d_eff = 1`, and a
**non-uniformity** statement, `K_0 ~ t^{-2/3}`.*

---

## 2. RULING ON W-07, WITH THE DISCOUNT APPLIED TO ITS REFUTERS EXACTLY AS HARD

### 2.1 WHAT SURVIVES

**(W-07/S1) THE REPRODUCTION OF W-06 — SURVIVES, AND I ADD A THIRD INDEPENDENT CONFIRMATION.**
Rebuilding the dressing from S1's own edge transports and the seed `20260816`
(`S_08_w07_reproduction.OUT.txt`) gives `amp = 0.271776442688027` and
`D_1 = amp |W_F^{-1} - W_C| = 0.384349931183` — W-07's registered separation, to `7.85e-14`, from
code that shares nothing with W-07's. `1000 of 4000` reproduces exactly, in **exact rational
arithmetic** at 60 digits (`S_05_ord_rho_name.OUT.txt`, row 1: `count = 1000`, `min D = 0`
exactly). M4-refuter-1 re-ran all 15 W-07 seals and all its legs byte-identically. **W-06's
restoration of S3-0 from VACUOUS to FALSE stands, confirmed twice from outside its own lane.**

**(W-07/S6) THE 52-PARTITION SWEEP IS A VACUOUS CONTROL — SURVIVES, AND M4-5's ATTACK ON IT
FALLS.** M4-5 claimed the sweep is contingent because the 5-cycle `C5` gives a different answer.
It does not: `m4_f_sweep.py:84` passes the **numeral** `2` to `C5`, whose own invariant count is
`E - rank(d1) = 5 - 4 = 1`. At `C5`'s own target the sweep returns exactly **one** winner and it
**is** the discrete partition — M4-refuter-1 re-ran it, and on 178 random connected graphs the
unique winner is the discrete partition **178 of 178**. M4's own leg F2 (`invariants = E - k + 1`,
forced by connectedness) proves W-07's point more sharply than W-07 did, three blocks above M4's
verdict against it. **W-07's application of the vacuity disqualifier was correct and correctly
scoped to a CONTROL, never to a theorem.**

**(W-07) THE ATTAINED/APPROACHED DISTINCTION — SURVIVES AS A STATEMENT, RELOCATED AS A FINDING.**
It is true, it is one line (`|Z_{k_0}| = 1` forces all supported characters to agree at `k_0`,
hence at every multiple, hence blank density exactly `1/|H|`), and §1 above shows it is a cut on
the **schedule exponent** — the only axis on which it does any work. W-07 was right that nothing
in the register distinguished them and right that they are different obstructions. It was wrong
about which axis the difference lives on.

**(W-07) "NO DURABLE RECORD IS SHOWN … THE READING IS TWO-WAY … DISCOUNT THIS ONE."** W-07's own
limitation block is the most accurate paragraph in the row and W-08 vindicates it: the reading
**was** two-way, and W-08 resolves it in neither of W-07's two directions but in a third.

### 2.2 WHAT FALLS

**(D-W07-1) THE OPERATIVE VARIABLE `ord(rho)`, FINITE VERSUS INFINITE — FALLS, AND FALLS AT BOTH
EDGES, THOUGH NOT FOR M4's REASONS.** M4-1 killed the infinite edge (an irrational `theta` with
`ord = infinity` returning `1000 of 4000`), and both refuters confirmed the closed-form bound at
`m4_g G3b`. M4-refuter-2 then supplied the cleaner and stronger kill, which I reproduce in exact
`Fraction` arithmetic with 60-digit decimal sines (`S_05_ord_rho_name.OUT.txt`):

```
theta                                    ord(rho)     floor(K/ord)   ACTUAL count<1e-9    min D
S1 PUBLISHED  theta = -1/4                      4             1000                1000    0 exactly
theta = 1/4 + 1e-13  (EXACT rational)      10^13                0                1000    6.830e-13
theta = 1/1000 + 4e-16 (EXACT rational)  2.5x10^15              0                   4    6.830e-13
```

Row 2 has **finite** order `10^13`, so W-07 §3's own sharp form `1000 = 4000/4 is ord(rho)`,
i.e. `count = floor(K/ord)`, predicts **0** — and the observable W-07 **tabulates** returns
**1000**. Irrationality was never needed; "finite versus infinite" was never the axis. Rows 2 and
3 pin `min_k ||k theta||` to the last digit and the count moves by 250x, so M4's replacement name
`min_k ||k theta||` fails too. **`ord(rho)` survives only for the count of EXACT zeros; W-07's
table does not count exact zeros, it counts cells below `1e-9`.** The correct name at fixed
`(K, tol, amp)` is the pair `(q, delta)` — the denominator of the nearby rational and the offset —
with `count = min(floor(K/q), ceil(eps/(q|delta|)) - 1)`.

**(D-W07-2) THE HEADLINE SENTENCE — "THE RECURRENCE OBSTRUCTION WAS MEASURED AT `ord(rho) = 4`" —
IS FALSE, AND W-07's OWN LEG D TABLE CONTAINS THE CONTRADICTION.** This is the most consequential
finding in W-08 and it was made by M4's **refuter**, against a lane whose entire assignment was to
refute W-07 and which instead certified the sentence as *"Confirmed true"* in its own M4-9. I
confirm it from the corpus text and from my own arithmetic.

The register's W-01 durability figures — *"the condition fires — `0.0247` at `n = 42` — and then
recurs to `0.99994`"* (`REGISTER_V001.md:55`), the figures STOP-FALLS-REBUILD was issued against —
were measured at `f = 2.0, c = 1.1, p = (0.4,0.15,0.15,0.15,0.15)`. The corpus states its own test
point verbatim at `S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md:438-444` and at
`S3_THE_CROSSING_V001.md:423-427`. My independent recomputation (`S_01_where_measured.OUT.txt`):

```
f = 2.0, c = 1.1, RS-G :   |Z_42| = 0.024654       (register's "0.0247 at n=42")
                           record maxima at k = 1, 2, 4, 6, 63, 154, 377, 6723, 7100, 99023, 106123
                           |Z_377| = 0.999941230   (register's "recurs to 0.99994", COR-E's 0.999941)
                           |Z_106123| = 0.999999981 (COR-E's k<=200000 figure, exact match)
```

**That connection is exactly resonant (`-11f + 20c = 0`), of INFINITE order, and `sup|Z_k|` there
is APPROACHED, never attained** — `0` cells above `1-1e-12` over `k <= 10^6`, by W-07's **own** leg
D table. At S1's order-4 connection the same observable takes exactly three values
`{1/sqrt(10), 2/5, 1}` (verified in exact Gaussian rationals: `|Z_k|^2 = 1/10, 4/25, 1/10, 1`,
product `1/625`, `lambda = -(1/2)log 5 = -0.804718956217050`), and `0.024654` is **`2.916e-01`
away from the nearest attained value**. W-06's `1000 of 4000` is the only recurrence figure in the
corpus measured at the order-4 point. **W-07 §1's headline — that the obstruction the register
calls undented was tested at the single point where it cannot fail — is false of the register's
own headline figures. It is true of exactly one figure, W-06's, and W-06 is the lane with no
artifact.**

**(D-W07-3) "`3*sqrt(3)/10` DOES NOT REPRODUCE — A FACTOR `sqrt(3)` NEEDS AN ELEMENT OF ORDER 3
AND `Z_4` HAS NONE" — FALSE.** Exact, in `Fraction` (`S_06_theorem_checks.OUT.txt`, S6c). On S1's
published connection the pair `(v0, v3)` has `(dF,dC) = (-1,0)`, so `D_k = amp|(-1)^k - 1|` takes
values in `{0, 2 amp}`. With `|s|^2 = (3/4, 4/25, 0, 9/100, 0)` (sums to 1 exactly),
`amp^2 = 27/400` and `(2 amp)^2 = 27/100 = (3 sqrt3/10)^2` **exactly**. The `sqrt 3` comes from
the ready state's amplitudes, not from the group. **The register's W-07 row is wrong on this
sentence.** W-07's page is not: it wrote *"Either a different observable, a different connection,
or a DIFFERENT NORMALISATION produced it"*, and this is escape three of the three it listed. The
custody conclusion — W-06's computational content is flagged, not inherited — is untouched:
reachability is not reproduction.

**(D-W07-4) LEG D's ROW LABELLED "GENERIC (BADLY APPROXIMABLE)" IS RANK 1, AND W-07 CAUGHT THE
CONFOUND IN LEG E AND DID NOT PROPAGATE IT.** `w07_d_carrier_recur.py:19` builds the row from
`2pi phi` and `2pi phi^2`; `phi^2 - phi = 1`, so `c - f = 2pi` exactly and `W_C = W_F`, `chi_0 = 1`.
`w07_e_isolation.py:1-3` records exactly this self-correction for leg B — *"Recorded rather than
silently fixed"* — and leg D was not revisited. **The register's `1.4e-12`** (*"generically `|Z_k|`
still returns to within `1.4e-12`"*) **is `1 - 0.999999999998574` from that rank-1 row at RS-P**,
not from a generic connection. M2-8 rediscovered this without citing W-07's own record of it, and
then M2-refuter-1 showed M2-8's correction moves TWO variables (connection rank **and** ready
state, `m2_f_floorlaw.py:120`'s `break` ending the reproduction loop after RS-P) and that the
missing fourth cell **reverses** it: genuinely rank-2 pairs **on RS-P** give `1.5055e-12`,
`2.0191e-15`, `1.0392e-12`, `3.8612e-12` against W-07's `1.4256e-12`. **The label is wrong; the
number is typical; and M2-8's "six orders of magnitude too deep" falls with it.**

### 2.3 THE VERDICT ON W-07

**W-07 SURVIVES ITS FIRST REFUTATION AS A REPRODUCTION AND AS A CUSTODY FINDING. IT DOES NOT
SURVIVE AS A NAMING.** Its figures are right, its vacuity ruling is right, its ATTAINED/APPROACHED
observation is right and misplaced, its self-limitation is exemplary. Its **operative variable is
misnamed** — the failure mode it predicted for itself, in the field it predicted it in — and its
**headline sentence about where the corpus measured recurrence is false**, refuted by a document
the corpus published at `S3_THE_CROSSING_V001.md:423` and by W-07's own leg D table. Its
supersession clause reads *"W-07 is superseded if … `ord(rho)` is shown not to be the operative
variable."* **It has been, three times over, and W-07 is superseded on that clause.**

---

## 3. THE EIGHT REFUTERS, READ SKEPTICALLY. **THREE KILLS ARE THEMSELVES CONFOUNDED.**

Discount applied to audits exactly as to builds. The isolation audit found 27 confounds across
64 comparisons; I checked the ones that carry a kill.

**KILLS I CONFIRM BY INDEPENDENT COMPUTATION** — M1-R1's **D2**, the `(2584,1597)` window artefact
(§4, item 2, where I go further than the refuter); M1-R1's **D1**, the rank-1 locus (§4, item 3);
M1-R2's **A2 == A0** vacuous control; M4-R1's **where-measured** kill (§2.2); M4-R2's **D1/D2**,
which the isolation audit rightly calls the cleanest isolation in W-08 and which I reproduce in
exact arithmetic; M2-R2's **adversarial schedule** and its **D1-vs-D2 under-read** (I read
`m2_d_attained.OUT.txt` and D1's *"the attained case is not slower; it is the FASTEST of the four"*
does sit twelve lines above D2's `+0.309362514` at `n=2` and `+0.343185988` at `n=3`); M3-R2's
**BARB / B1q** counterexample to M3's headline; M4-R1's **M4-5 refutation**.

**KILLS THAT ARE THEMSELVES CONFOUNDED, AND MUST NOT BE CARRIED AS STATED:**

1. **M2-REFUTER-1's `F_VWA` ROW (its "fitted exponent `-2.999` at `dim H = 2`") IS A STEP FUNCTION
   FITTED AS A POWER LAW.** The isolation audit found this and it is right: leg B's `F_VWA` pair
   gives `F = 3.551e-08` at `K = 1e3, 1e4, 1e5, 1e6` and then `3.609e-23` at `1e7`. A five-point
   log-log slope through four identical values and one jump is the slope of the jump, not a
   scaling exponent — and the same refuter's own verdict elsewhere says *"a 5-point fit of a
   running minimum is not a `d_eff` classifier"*. **The number `-2.999` must not enter the
   register.** What survives is the **omission** finding: `F_VWA` is genuinely absent from
   `m2_f_floorlaw.py:33-38`'s hard-coded `CONN` list and from `m2_e_schedule.py:107` while the
   pickle those legs load contains it, and it is the one class contradicting both laws. The
   omission is a real defect; the replacement figure is not evidence.

2. **M2-REFUTER-1's "M2-4 IS ONE-SIDED" IS CORRECT; ITS "THE TWO CLUSTERS DO NOT OVERLAP IS FALSE"
   IS CORRECT AND CHANGES NOTHING LOAD-BEARING.** I record it as sustained but downgraded: the
   refuter's own resampling gives means `-1.010`/`-0.986` against theory `-1` and
   `-1.979`/`-1.932` against theory `-2`, i.e. **it confirms the law in the mean while refuting
   "fixed by `d_eff` alone"**. My `S_09` result is the stronger statement in the same direction and
   is derived from the `N(eps)` law rather than fitted.

3. **M1-REFUTER-1's D1 IS RIGHT AND ITS EXHIBIT IS MIS-STATED.** The mathematics is correct and I
   verified it exactly (§4, item 3). But the isolation audit's catch stands: `R1_05`'s two arms sit
   at different `N` (A1 reported at `N >= 1e3`, A2 only at `N = 10`) and the quoted "differ by 50.8
   at `N=10`" is an A2-versus-`m(Q)` dip — a Birkhoff average set against a closed-form Mahler
   measure — **which is exactly the statistic substitution R1 charges M1's COMPARISON 5 with.** The
   audit supplied the missing arm (A1 at `N = 10` is `-0.121090780`, true gap `51.90`) and the
   conclusion survives. **Carry the conclusion, not the exhibit.** And both of R1's "identical
   relation lattice `L = Z(1,1)`" points are exact rationals (denominators `10^60`, `10^221`),
   hence rank-2 lattices; the rationality is flagged for one and not the other.

4. **M3-REFUTER-2's "THE FRAMING INVERTS UNDER MEASURE" SWAPS YARDSTICKS MID-SENTENCE.** Its two
   numbers (connection-side firing region `1/4 -> 1/2`; state-side `1/4 -> 1/4`) are exact,
   independently derived, and not in dispute — the closed form `0 in conv{1,u,v,uv} <=>
   cos f + cos c <= 0` is a genuine new result. But M3's state-side finding was never about the
   **measure** of the firing set; it was about which **predicate** decides it (triangle inequality
   `-> w1+w4 <= w2+w3`), and that does change. The word "content" does different work on the two
   sides. **The inversion claim is not sustained; the closed form and the doubling are.**

5. **M2-REFUTER-2's "88.5% OF ATTAINED POINTS DEVIATE BY LESS THAN 4.8%" AND "60.5% DECAY SLOWER"
   ARE MEASURED UNDER A PRIOR THE LANE DOES NOT DECLARE AS SUCH** — counting measure on closures of
   exact order `<= 120`. The refuter *does* print the prior-sensitivity elsewhere (its own volumes
   move from `1/4, 1/2` to `0.28057, 0.55985` under a different prior) and does not carry the
   caveat to these two percentages. The **directional** finding survives outright, because it is
   an existence claim refuted by the lane's own table (`n=2`, `n=3`); the **percentages** should be
   quoted with their prior or not at all.

**AND ONE PATTERN THE ISOLATION AUDIT NAMED THAT I RATIFY AS THE MOST USEFUL SENTENCE IN W-08.**
The commonest FATAL defect across all twelve lanes is not "two variables moved" but **"zero
variables moved"** — five of eleven FATALs are controls that could not have failed, one of them
(`M1_03_three_connections.py:168`, `A2 == A0 == 10937044409`, difference exactly 0) a control in
which literally nothing varied and which is reported in F10 as a confirmation. **An isolation
ledger cannot detect any of these, because a ledger records what the author intended to vary.**
Three build lanes wrote ledgers and it did not stop them; five of the seven FATAL build confounds
sit in the arm the lane's own ledger names as decisive. **Treat "the ledger says one variable" as
evidence of nothing.**

---

## 4. CORRECTIONS OF RECORD FROM W-08, NUMBERED

1. **W-01's registered recurrence figures `0.0247 at n=42` and `0.99994` are measured at
   `f = 2.0, c = 1.1`, RS-G — the EXACTLY RESONANT connection of the erratum against W-02, of
   INFINITE order, where `sup|Z_k|` is APPROACHED.** Not at S1's order-4 connection.
   `REGISTER_V001.md:55` and the W-07 row §1/§4 both need the correction.
   (`S_01_where_measured.OUT.txt`; `S2…AUDIT:438-444`; `S3_THE_CROSSING_V001.md:423-427`.)
2. **`LANE_W08_M1_IDENTIFICATION/M1_07_lawton.py:50` publishes `+4.716e-09` for the `(2584,1597)`
   Boyd-Lawton row and the true value is `~+4.24e-08`, nine times larger.** The `np.roots` branch
   caps at degree 1200, so that row alone switches evaluator to a `2^23`-point quadrature which
   landed on an accidental zero-crossing. I go past the refuter: besides reproducing
   `2^23 -> +4.7157e-09`, `2^24 -> +4.0392e-08`, `2^25 -> +4.2415e-08`, `2^26 -> +4.2432e-08`,
   `2^27 -> +4.2432e-08`, **I ran the `np.roots` route at degree 4181 anyway and got `+4.2491e-08`**
   — a third route, agreeing with the converged quadrature to `5.9e-11`
   (`S_04_lawton_window.OUT.txt`). The `(610,377)` row `+5.2505e-07` is converged by both routes
   and stands. Boyd-Lawton accumulation is **not** overturned; the headline figure is.
   **COR-E defect class, inside the lane that cites COR-E.**
3. **M1-T4's confinement of the Diophantine fragility to `H = T^2` is false, on a codimension-one
   locus that is row 3 of M1's own 16-row table.** On `uv = 1` (i.e. `W_C = W_F`),
   `P|_H · z = 0.3z^2 + 0.4z + 0.3`. Clearing denominators, `3z^2 + 4z + 3` is palindromic (roots
   in reciprocal pairs) with discriminant `16 - 36 = -20 < 0` (roots complex conjugate); conjugate
   **and** reciprocal forces `|r|^2 = r·(1/r) = 1`, so **both roots lie exactly on `|z| = 1`**, by
   two lines over `Z` with no numerics (`S_06…OUT.txt`, S6b; numerical check `2.220e-16`). Hence
   `log|Z_k|` is unbounded below there and its Birkhoff average carries the same inhomogeneous
   Diophantine hypothesis. `m = log 0.3 = -1.203972804`. **T4's conclusion is untouched** — every
   pathology on this locus drives `lambda` **down**, and `1-|P|` stays continuous.
4. **M1_04's PART 3 "`lambda_H = m(P)` iff `H = T^2` iff `u,v` have no multiplicative relation" is
   false on a positive-measure set**, and M1's own nine-case table contains two counterexamples it
   read as generic: row 1 (`u=e^{-2i}, v=e^{1.1i}`) is the erratum's exactly-resonant point
   labelled *"full support, generic connection"*, and row 5 (`u=e^{0.7i}, v=e^{1.3i}`) satisfies
   `13(0.7) - 7(1.3) = 0` exactly. **Zero of the nine rows has `rank L = 0`.** The correct operative
   variable is the pair `(L, pi)`: on the open region `max(p10,p01,p11) > 1/2`, `log|P|` is
   continuous on `T^2`, plain Weyl suffices and **no** Diophantine hypothesis is needed at any rank.
5. **The register's W-07 row sentence "a factor `sqrt 3` needs an element of order 3 and `Z_4` has
   none" is false** — `(2 amp)^2 = 27/100` exactly at `|s|^2 = (3/4, 4/25, 0, 9/100, 0)` on the
   order-4 connection (§2.2 D-W07-3). The custody finding is unaffected.
6. **M2-9's `D_H = 1 - 2/pi = 0.363380227639` is wrong in the 12th digit**; `1 - 2/pi =
   0.3633802276324186`. My measured density on RS-P is `0.363380` at `K = 10^7` from four
   different connections (`S_02_the_race.OUT.txt`), all agreeing, so the corpus value is fine and
   only the printed constant is off by `6.58e-12`. **COR-K defect class, minor.**
7. **M2-8's "six orders of magnitude too deep" and M2's leg F4 fall**: `m2_f_floorlaw.py:120`'s
   `break` ends the "reproduction first" loop after RS-P, so the comparison moves ready state
   **and** connection rank; the missing fourth cell reverses the direction (§2.2 D-W07-4).
8. **M3's headline "it fails on ANY carrier with a spectator vertex" is refuted by the corpus's own
   control carrier `B1q`** (`S4_THE_MEASUREMENT_V001.md:519`, multiset at `:582`), which the lane
   **cited and printed** at `m3_4_threeclass_vs_fourclass.py:34-36`. The operative condition is **all
   four classes occupied**, which the lane's own isolation C5 proves and its headline replaces
   with a false weaker form.
9. **`m3_5_theoremA_both_sides.py:85-97` and `m3_2_fourclass.py:89-90` are controls that could not
   have failed** — the first restates the definition of a convex hull, the second prints `"yes"`
   from a conditional that only checks one direction. Both claims are true; neither was measured.
   Voids the controls, never the theorems.

---

## 5. WHAT READS TWO WAYS, AND IS NOT SCORED AS EITHER

1. **THE `RESONANT` ROW OF MY OWN `S_09` TABLE.** Its four-point fitted exponent is `0.591`
   against a theory of `0.667`; its **last-decade** exponent is `0.6666`. A window in which a
   primitive relation of height 31 has not yet resolved is exactly the COR-E class. **I state the
   exponent as a lower bound with the last-decade figure, and I do NOT score the row as
   confirmation of the `d_eff = 1` law.** The other three rows agree at four points and at the last
   decade both.
2. **THE ADVERSARIAL SCHEDULE.** `|Omega| ~ 0.55` forever with unboundedly many writes reads as
   *"durability is a schedule stipulation the corpus never made"* **or** as *"`K^{-1/2}` write
   density is not an admissible schedule and the stipulation is the right one"*. **W-08 does not
   distinguish them and does not score itself as having.** M2's own self-flag says the schedule
   half has no theorem; §1 supplies the exponent, not the admissibility criterion.
3. **THE NULL "`ord(rho) = infinity` GIVES 0 OF 4000, FIVE TIMES OUT OF FIVE" (W-07 leg E).** It
   reads as *"finite order is what does it"* **or** as *"the five draws were far from every low
   denominator, which is a probability-one event and was never a test"*. M4-6 said the second, M4's
   refuter 1 showed M4-6 slid a probabilistic near-certainty into a **logical** disqualifier (the
   failure set is non-empty **inside W-07's own table**: `theta = -1/4`; positive measure
   `>= 4.685e-06`), and M4's refuter 2 showed the two named-algebraic rows are deterministically
   safe by a continued-fraction bound. **The correct reading is the third: five-for-five is a
   density statement, not an isolation of `ord(rho)`, and IMP-1 does not apply.**
4. **`lambda` CAN BE SLOWER THAN GENERIC OFF THE DEGENERATE POINT.** M2's D2 gives `-0.458145366`
   at `n=2` and `-0.424321892` at `n=3` — attained, and `+40%`/`+45%` slower than
   `-0.767507880`. This reads as *"degeneracy can help the record persist"* **or** as *"`lambda`
   is simply a non-monotone function of `H` and direction has nothing to do with degeneracy"*. The
   second is supported by M2's D3 (same order, different direction, `-0.860504844` vs
   `-0.804718956`). **Neither is scored. What IS decided is that every one of these rates is
   strongly negative, so the race verdict of §1 does not depend on which reading is right.**

---

## 6. THE WEAKEST LOAD-BEARING CLAIM IN THE CORPUS

**It is not in W-08. It is `W-01`'s convex-hull criterion as the register states it, at
`REGISTER_V001.md:43`, and specifically the clause "three unit-modulus coefficients".**

The criterion `<M_dF s, M_c s> = 0 <=> 0 in conv{conj(W_F)W_C, conj(W_F), W_C}` is W-01's ruling,
is the thing "the audit refuted the build and wrote", is quoted forward by S3, S4, W-03 and every
W-08 lane, and is the sentence STOP-FALLS-REBUILD's rebuild would be built on. It is
**K1-scoped twice over and the register says so nowhere**:

- **"three"** is `p00 = 0`, an incidence fact about K1 (no vertex outside both loops). On any
  carrier with a spectator vertex there are **four** coefficients and the criterion becomes
  `cos f + cos c <= 0`. **I verified this myself rather than carrying the refuter's report**
  (`S_10_w01_offK1.OUT.txt`, 200000 draws, one variable moving — the number of occupied
  characters — with the same hull test, grid and evaluator in both arms): the closed form agrees
  with the four-character hull on **200000 of 200000** points, and the firing region moves from
  `0.250850` (exact `1/4`) to `0.501075` (exact `1/2`). **The firing region doubles.**
- **The quantifier is dropped.** M3-F3 is right that `iff` is false read literally; M3's refuter 1
  is right that the failing direction is HULL `=>` FIRE (`17154 of 46800`) while
  `Z_1 = 0 => 0 in conv` is unconditionally true and is the half W-01's own exhibit travels on.
  **No W-01 result is touched — and that is exactly why the defect has survived seven layers
  unexamined.**
- **W-01's own advertised property does not survive one added vertex.** *"It distinguishes
  curvature from flat holonomy, which K1 exists to separate"* — on K1-plus-a-spectator the
  criterion `cos f + cos c <= 0` is **separable** and invariant under `f -> -f` alone. My own
  run (`S_10_w01_offK1.OUT.txt`, same 200000 draws, only the character count moving): sending
  `f -> -f` alone changes the verdict at **100215 of 200000** points on K1 and at **0 of 200000**
  on the spectator carrier. Witness: `(pi/2, pi/2)` fires on K1, `(-pi/2, pi/2)` does not, and
  **both** fire on the spectator carrier. **The property the register advertises as the
  criterion's virtue is a coincidence of `p00 = 0`.**

Every W-08 lane rests on it, none of the four build lanes tested it off K1, and the two lanes that
went near it (M3 and its refuters) established the generalisation and then mis-stated its scope in
opposite directions. **It is load-bearing, it is quoted forward unqualified in the corpus's most
cited sentence, and it is the only load-bearing claim in the corpus that has never been checked on
a second carrier by a lane whose job it was.**

Runner-up, for the record: **M2-4's floor law stated two-sidedly.** `F(K) <= C_2 K^{-1}` is
Dirichlet, one-sided; the matching lower bound is badly-approximable, a measure-zero hypothesis;
and the lane's own `F_VWA` data breaks the lower half. That one is at least *known* to be weak,
which the W-01 sentence is not.

---

## 7. WHAT W-08 DOES NOT ESTABLISH — BEFORE THE REOPENS, NOT AFTER

**No durable record is constructed here.** §1 decides a race between two computed objects on K1;
it does not build anything. **The schedule half has an exponent and no admissibility criterion.**
**Nothing off K1 is decided**: BARB and K1S are refuter constructions, `B1q` and `B0a` are S4
controls, and W-08's four build lanes ran on K1 alone. **The dressed observable is still a
reconstruction and W-06's code still does not exist.** **And this is layer nine of one lineage,
with no refuter of this page.**

---

## 8. REOPENS

**W-08's race verdict (§1) REOPENS IF:**
- a connection with `G != {1}` is exhibited on which `SUM_{k<=K}(1-|Z_k|)` is sublinear in `K`
  — which would refute M2-1's identity, so this is a request to find an arithmetic error in
  `S_06_theorem_checks.OUT.txt` / the two refuters' exact-`Fraction` versions; **or**
- an *intrinsic* admissibility criterion for schedules is written down under which the corpus's
  `SUM(1-z_n)` test is recovered, and the `K^{-1/2}` / `K^{-1/3}` adversary is excluded — in which
  case the ATTAINED/APPROACHED cut of §1 changes character; **or**
- the monotonicity of `|Omega_N|` is shown not to be the right reading of "durable record" —
  i.e. someone states a durability observable on which a single-cell near-return **does** undo
  earlier writes. **That is the one move that would restore the founding obstruction, and nobody
  in this corpus has attempted it.**

**W-08's ruling on W-07 (§2) REOPENS IF:** W-06's lane code is produced and its recurrence figure
turns out to be measured somewhere other than the order-4 point (which would move D-W07-2 in
W-07's favour) · or the `(q, delta)` naming of §2.2 is shown not to survive varying `K` and
`tol/amp`, which its own author has already flagged as untested.

**THE ONE COMPUTATION THAT WOULD MOVE W-08 MOST, AND IT IS NOT A REBUILD:** run the W-01
convex-hull criterion, in the register's own words, on **`B1q` and `B0b`** — two carriers the
corpus already owns and has already published class multisets for (`S4:519`, `S4:582`) — and
report which of the two readings (`0 in conv`, or `w1+w4 <= w2+w3`) each obeys, with the
quantifier restored. It is one script against two existing objects and it decides §6.

**W-08 IS SUPERSEDED IF:** a lineage-independent lane fails to reproduce `|Z_42| = 0.024654` at
`f=2.0, c=1.1` · or fails to reproduce the band-share table of §1(R3) · or shows that
`|Omega_N|`'s monotonicity is not what the corpus means by a record.
