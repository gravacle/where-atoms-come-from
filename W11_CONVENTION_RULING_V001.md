# W-11 — THE CONVENTION RULING

**THE QUESTION (W-10 §0.4.1, "the single most consequential undecided question the table produces").**
Is *"the formation functional is carrier-independent"* a **FINDING**, or a restatement of the
**TRANSPORT CONVENTION**?

**RULING: READING B.**

**Ruled from `/Users/bgm/MB Work/where-atoms-come-from/LANE_W11_RULING/` — 7 scripts, 6 captured outputs,
`PUBLISHED_CONVENTIONS.txt` and `SEALS.sha256`; 15 files, 14 sealed entries, all verify. Seed `20260818`, shared with no
lane. Connection `f = 1.0, c = sqrt(2)` — `S4_THE_MEASUREMENT_V001.md:603`, the only generic pair
the corpus publishes. Nothing is imported from any W-11 lane.**

---

## 0. THE ANSWER IN ONE LINE, AND WHY IT IS NOT A JUDGEMENT CALL

`Z_n = <branch_F, branch_C> = <s, Q_n s>` where **`Q_n = (branch_F)^* (branch_C)`** — the
**relative branch operator at the tick the record is read**. Everything in the corpus's functional
layer is a statement about `Q`.

W-01 defines `(M_gamma s)(v) = W(gamma)·s(v)` on the loop and `s(v)` off it. Therefore

```
diag(M_F)_v = W_F^(1[v in gamma_F])          measured deviation: 0.00e+00
diag(M_C)_v = W_C^(1[v in gamma_C])          measured deviation: 0.00e+00
Q_k = diag( conj(W_F)^(k·a) · W_C^(k·b) )    at class (a,b);  1.13e-15 (K1), 1.56e-15 (B0b)
                                             over 500 connections, k = 1,2,3,7
```

**`M_gamma`'s diagonal *is* the incidence indicator.** The four classes are the joint level sets of
the two operators' diagonals, *by the definition of both*. So

```
Z_k = SUM_v (class character)_v |s_v|^2 = p00 + p10·u^k + p01·v^k + p11·(uv)^k
```

is reached **by substitution**. No lemma intervenes. No other outcome was available. That is what
"restatement" means, and it is prior to K1, prior to B0b, prior to any carrier at all.

**The corpus already owns half of this and never joined it up.** W-03 registered: *"'topology is
inert' is an **analytic fact about the construction**, true before any carrier was built. The
controls tested nothing"* (`REGISTER_V001.md:191`). W-06 registered the mechanism: *not scalar
multiplication — **FIBRE-WISE-NESS**"* (`:577`). W-10 registered that `M_gamma` *"is literally an
element of the gauge group `U(1)^V`"* (`:409`). Three rows, one fact, never assembled. §0.4.1 asks
whether the *rest* of carrier-independence — everything beyond topology — has the same status.
**It does, by the same identity, and this ruling assembles it.**

---

## 1. WHAT IS ACTUALLY BEING RULED — THE TEST MADE OPERATIONAL

A **restatement** is a statement whose truth follows from the definitions by substitution, with no
auxiliary lemma, and which therefore could not have come out otherwise. A **finding** requires a
step that could have failed. The question is decidable on that test, and I do not return
UNDECIDABLE.

The corpus states as one sentence three claims of **three different statuses**, and Reading A draws
all its force from the conflation:

| | claim | status |
|---|---|---|
| **(1)** | the functional depends on the state only through `pi` | **RESTATEMENT** — the identity above |
| **(2)** | the functional depends on the carrier only through `pi` | **RESTATEMENT** — the same identity |
| **(3)** | `lambda = m(p00 + p10 x + p01 y + p11 xy)`, a logarithmic Mahler measure | **THEOREM, with a hypothesis** |

Reading A's own wording is *"a substantive result — the physics does not depend on the complex,
**which is why** one short mathematical note (N1) is publishable."* **That "which is why" is a
non-sequitur, and it is the whole of Reading A's apparent strength.** (3) does not rest on (1) or
(2). Killing (1) and (2) does not touch (3). **See §4.**

---

## 2. RULING ON `T`'s CANONICITY — the registrar's own declared weak point

**`T` IS NOT CANONICAL AMONG ALL UNITARIES WITH `U^L = M_gamma`. `T` IS CANONICAL — EXACTLY — AMONG
OPERATORS ADMISSIBLE UNDER THE CORPUS'S OWN CLAUSE. AND THE FINDING DOES NOT NEED EITHER FACT.**

### 2.1 Not canonical in the wide family, and the wide family is large

`U^L = M_gamma` on the loop means `U = w·V` with `w^L = W` and `V^L = I`, so the family is
positive-dimensional. `{n : U^n is diagonal} = d·Z` with `d | L` — verified over 1200 draws per
loop per carrier: K1 `{1:122, 3:1078}` and `{1:126, 3:1074}`; B0b `{1:24, 2:123, 4:1053}` and
`{1:127, 3:1073}`; every observed `d` divides `L`, `max||U^L - M|| = 7.4e-15`.

### 2.2 Canonical under the corpus's own admissibility clause — and the clause is already written

**`S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md:657`, CHOICE LEDGER A1**, verified at the bytes,
states the corpus's own test for an admissible loop operator: *"(a) unitary, (b) uses no data
beyond fibres/edges/orientation/connection, (c) reduces to the build's `T_gamma` on `L_v0`."*

**Clause (c) is CIRCULAR and I strike it.** `T` maps `L_v0` into `L_v1`; requiring an operator to
"reduce to `T_gamma` on `L_v0`" requires it to be an endomorphism of each fibre — i.e. it
presupposes fibre-wise-ness, which is the property at issue. A1's own flag block says A1 is
*"the one place Construction A could be attacked and an S3 lane should attack it there first"*
(`:686`). No lane ever did.

**Clauses (a) and (b) survive, and with LOCALITY they settle it.** Locality is not imported: it is
`S1_CARRIER_K1_V001.md:52-53`, which *defines* parallel transport — *"Parallel transport along
`e : u -> v` is `z |-> U_e z`. Reverse traversal transports by `U_e^{-1}`."* — the corpus's own
definition of the thing `M_gamma` is named after, and one that `M_gamma` does not satisfy.

> **THEOREM (mine, four lines; `r4_uniqueness.py`).** On an `L`-cycle write the general local
> operator as `U[v,v] = g_v`, `U[v+1,v] = c_v U_v`, `U[v-1,v] = d_v conj(U_{v-1})`.
> 1. A closed walk of `L` steps has `(#forward - #backward)` in `{-L, 0, +L}`, so each **diagonal**
>    entry of `U^L` is `(prod c)·W + N_vv + (prod d)·W^{-1}`, with `N` constant in the connection.
> 2. `U^L = W·I` identically in the connection forces `prod c_v = 1`, `prod d_v = 0`, `N = 0`.
> 3. Column `v` has exactly those three entries, so unitarity gives
>    `|g_v|^2 + |c_v|^2 + |d_v|^2 = 1`, hence `|c_v| <= 1`.
> 4. `prod |c_v| = 1` with every `|c_v| <= 1` forces `|c_v| = 1`, hence `g_v = d_v = 0` for **all** `v`.
>
> **`U = Lam·T`** with `Lam` unimodular diagonal, `prod c_v = 1` — COR-F's `T` times an
> `(L-1)`-torus of **adjoined per-vertex phases**. Clause (b) supplies no vertex-indexed phase, so
> setting them to `1` gives **exactly `T`**. A pure backward shift delivers `W^{-1}`, not `W`, so
> the loop's own orientation (`S1:27`) fixes the direction.

Steps 3–4 are arithmetic. Steps 1–2 are verified numerically by a statistic that prints `O(1)`
when false: diagonal Fourier support off `{-L,0,+L}` is `4.4e-16 / 8.9e-16 / 3.6e-15` at
`L = 3,4,5`, with `A_{+L} = (prod c)·W` to `1.4e-15` and `A_{-L} = (prod d)·conj(W)` to `2.0e-15`;
and the root defect over an `(s,r)` sweep of the diagonal and reverse bands is zero **in exactly one
cell**, `s = r = 0`.

### 2.3 The two criteria that discriminate, one that does not, and one that is void

| criterion | source | `T` | `M_gamma` | `D = diag(W^{1/L})` |
|---|---|---|---|---|
| unitary | A1(a) | pass | pass | pass |
| carrier data only | A1(b) | pass — entries are single `U_e` | pass — entries are the loop product | **FAIL** — needs a lift of `W`; exhaustive `m in [-4,4]^3` finds **no** monomial equal to `W^{1/3}` |
| continuous in the connection | degree: `L·deg(h) = 1` is impossible | pass, max step `1.6e-03` | pass | **FAIL**, max step `3.0e+00` over one winding |
| gauge covariance | COR-J | `4.6e-15` | `9.4e-15` | `3.3e-15` — **discriminates nothing** |
| locality | `S1:52-53` | pass | **FAIL** — moves nothing | **FAIL** |

**COR-J does not do the work three lanes assign to it.** `T` is gauge-covariant, so `Z^T_n` is
gauge-**invariant** at every partial tick. *"The record must be gauge-invariant"* excludes neither
`T` nor `D`.

### 2.4 WHAT THE FINDING DEGRADES TO IF `T` IS NOT CANONICAL — **NOTHING**

`T`'s role is **exhibitory**: it shows the class-constant-diagonal condition is not automatic. That
needs one **admissible** witness, not a **canonical** one — and the witness is supplied by the
corpus's own sealed audit, graded CONFIRMED-WITH-CORRECTIONS at D1. Moreover the whole admissible
family is a witness: **400 draws from the full `(Lam·T) × (Lam·T)` family, 800 root checks all
passed, `0` pi-blind, smallest spread `6.0e-01` (K1) and `5.0e-01` (B0b).**

So the degradation is: from *"the corpus's canonical rival transport makes incidence visible"* to
*"**every** admissible edge tick makes incidence visible, and the corpus's operator is the one point
of the family at which it does not."* **That is logically stronger, because it quantifies over the
family and needs no canonical member.**

### 2.5 THE BRIEF'S DISJUNCTION IS NOT EXHAUSTIVE AND BOTH HORNS ARE FALSE AS WRITTEN

- **Horn 1** — *"if a different, equally natural edge tick restores invisibility, Reading B falls."*
  `D` restores invisibility (`5.6e-16` at every tick). **But `D` is not an edge tick.** It moves no
  fibre value (`|<e0, D e0>| = 1.000` against `0.000` for `T`), it is discontinuous in the
  connection, it fails A1(b) — and `diag(D)_v = w^(1[v in gamma])` is *the class indicator again*.
  **`D` is the corpus's own convention at a finer clock, not a rival to it.** Horn 1's antecedent is
  never met.
- **Horn 2** — *"if every unitary with `U^L = M_gamma` except the diagonal ones breaks invisibility."*
  **False.** Correlated non-diagonal pairs `U_F = L_F R`, `U_C = L_C R` (`R` block-diagonal on the
  four classes) are pi-blind: **200/200 draws with `U_F` non-diagonal, worst spread `6.3e-15` (K1),
  `5.8e-15` (B0b).** An independent sampler cannot find them at any sample size — the pi-blind set
  is the correlated locus, measure zero for independent draws. Every lane that reported "`0` of
  `N` random draws preserve invisibility" reported a control that could not have failed.

**The correct trichotomy: the pi-blind set is neither "the diagonal ones" nor "the canonical one".
It is exactly the set of pairs whose RELATIVE operator is class-constant diagonal — see §3.**

---

## 3. THE OPERATIVE VARIABLE, NAMED — the seventh attempt, with the guard run first

**NAMING GUARD.** `grep` over `REGISTER_V001.md` and all six sealed corpus artifacts:
`"fibre-wise"` — register **4**, corpus 0 (W-06's N4 correction, `:577`).
`"class-constant"`, `"class function"`, `"incidence class"`, `"advance"`, `"comparison time"`,
`"class-uniform"` — **0 everywhere.** I therefore attach a clause to an existing register term
rather than coin a seventh name (W-10 N-7).

> **THE OPERATIVE VARIABLE IS WHETHER THE *RELATIVE BRANCH OPERATOR*
> `Q_n = (branch_F^n)^* (branch_C^n)`, AT THE TICK THE RECORD IS READ, IS MULTIPLICATION BY A
> FUNCTION OF THE INCIDENCE CLASS — i.e. FIBRE-WISE (`REGISTER:577`) **AND** CLASS-CONSTANT.**
>
> **BICONDITIONAL.** For *any* pair of operators and any tick, `|Z_n|` is a function of `pi` alone
> **iff** `Q_n` is class-constant diagonal. Verified over **9000 (pair, tick) cells** — 4500 per
> carrier, five operator families including COR-F's own `T` and unrestricted Haar pairs —
> **AGREE 9000, DISAGREE 0**, with **3900 cells** on the class-constant side and **5100** on the
> other, so the test could have failed and did not.
>
> **The corpus's convention makes `Q_k` class-constant diagonal BY CONSTRUCTION**, its diagonal
> being literally `conj(W_F)^{ka} W_C^{kb}` at class `(a,b)`. **That sentence is Reading B.**

**Both coordinates are operative, and each is refuted alone.**
- Hold the tick at `n = 1`, move only the operator: `T` spreads `6.1e-01` (K1) / `8.1e-01` (B0b);
  `D` spreads `3.3e-16` / `2.2e-16`. **So "the comparison time" is not the operative variable.**
- Hold the operator at `T`, move only the tick: spread `7.5e-01` over `n <= 24`, and `5.6e-16`
  restricted to `n = 0 mod lcm(L_F, L_C)`. **So "the transport" is not the whole of it either.**
- **"DIAGONAL" is wrong in both directions.** A *diagonal* `L`-th root of `M_gamma` whose labels
  differ **within** a class is fibre-wise at every tick and separates the arms at `6.75e-01` (K1),
  `3.07e-01` (B0b), one variable moved. And the correlated pairs above are non-diagonal and blind.
- **"Fibre-wise and LOOP-constant" is too strong.** A *class*-constant `A_F` that is not
  loop-constant (different phases on classes `10` and `11`) is still pi-blind — the weaker
  condition is the right one.

`Q` **appears in no register row.** Every result in the corpus's functional layer is a statement
about it.

---

## 4. RULING ON N1 — the result proposed for publication

**N1 SURVIVES READING B INTACT. IT NEEDS TWO STATED HYPOTHESES, NEITHER OF WHICH IT CURRENTLY
CARRIES, AND ONE SENTENCE OF ITS FRAMING MUST BE WITHDRAWN.**

**Why N1 is not a restatement.** The Mahler step is a genuine theorem about `(pi, characters)`, not
about `M_gamma`. It survives *changing the operator* inside the fibre-wise class: under the
fibre-wise root `D` — a finer clock, fractional winding, neither branch ever "closed" —
`Z^D_n = p00 + p10·u^{n/L_F} + p01·v^{n/L_C} + p11·(...)` holds at every tick to `2.5e-16` (K1) and
`2.2e-16` (B0b), and the rate is `m(pi)` to `3.9e-07` and `8.0e-06` at `N = 2e5`. A restatement of
`M_gamma` would not survive replacing `M_gamma`.

**Why it needs a hypothesis.** `m(P)` is the *torus* average. Measured on K1, `pi = (0,.3,.3,.4)`,
`N = 2e5`, one variable moving:

| connection | rate | `|rate - m(P)|` |
|---|---|---|
| `f = 1.0, c = sqrt(2)` — generic, `S4:603` | `-0.767504139` | `3.7e-06` |
| `f = 2.0, c = 1.1` — resonant, `-11f + 20c = 0` | `-0.767026255` | **`4.8e-04`** |
| `f = pi, c = 3pi/2` — order 4, S1 §6's own connection | `-0.804718956` | **`3.7e-02`** |

`m(P) = -0.767507880357` by Jensen at `2^20`, cross-checked by a Jensen-free 2-D quadrature
(`3000^2` grid; agreeing to `7e-8` at this `pi`, where `p00 = 0` puts a zero of `|P|` on the torus
and slows the quadrature, and to `1e-12` at the other two `pi` tested) — the register's value is right. But the corpus's **ERRATUM AGAINST W-02** (`REGISTER:162-175`) is exactly this
correction already, and **N1 as registered does not carry it.**

> **N1, AS IT SHOULD BE STATED.**
>
> *Let a finite complex carry a rank-one Hermitian bundle with a `U(1)` connection, let `gamma_F`
> and `gamma_C` be two designated loops with holonomies `W_F`, `W_C`, and let `pi = (p00,p10,p01,p11)`
> be the pushforward of `|s|^2` onto the four incidence classes `(1[v in gamma_F], 1[v in gamma_C])`.*
>
> ***HYPOTHESIS 1 (the convention, and it is a stipulation, not a theorem).** The two branches are
> generated by operators whose relative operator `Q_n` is multiplication by a function of the
> incidence class — in particular by the whole-circuit holonomy operators `M_gamma`. This is
> exactly the condition under which the left-hand side below depends on `s` only through `pi`. The
> corpus's own sealed audit exhibits an admissible alternative — COR-F's edge tick `T`,
> `S3_THE_CROSSING_AUDIT_V001.md:160-209, :794` — under which it does not.*
>
> ***HYPOTHESIS 2 (equidistribution).** `(conj(W_F), W_C)` generates a dense subgroup of `T^2`.*
>
> *Then `lim_N (1/N) SUM_{k<=N} log |<M_F^k s, M_C^k s>| = m(p00 + p10 x + p01 y + p11 xy)`, the
> logarithmic Mahler measure. If instead the pair generates a proper closed subgroup `H < T^2`, the
> limit is the average of `log|P|` over `H`, which differs from `m(P)` in general — by `4.8e-04` at
> `f = 2.0, c = 1.1` and by `3.7e-02` at `f = pi, c = 3pi/2`.*

**WITHDRAW from N1's framing:** *"the physics does not depend on the complex, which is why one short
mathematical note is publishable."* The note is publishable. **That is not why.** What is
publishable is the identification of a branch-comparison decay rate with a logarithmic Mahler
measure, and thence with the entropy of an algebraic `Z^2`-action. What is **not** publishable as a
discovery is the blindness to the complex, which is the hypothesis, restated.

**PUBLISH N1. Under Hypotheses 1 and 2, stated. Six to eight pages, one theorem.**

---

## 5. WHAT READS TWO WAYS, AND IS SCORED NEITHER WAY

**(i) THE TRIVIAL-CONNECTION CONTACT POINT.** Pre-registered at `FOUNDING_DESIGN_V001.md:117-118`
and `S2:583` (*"No formation at trivial connection is the known trivial answer"*). Measured at
`a = 0`: under `M_gamma`, `min_n |Z_n| = 1.000000000000` exactly on both carriers; under `T`,
`0.375665` (K1) and `0.568562` (B0b) — **`T` fires at zero field.** So the contact point *does*
discriminate.

*Reading (α):* the criterion's admitted set **is** the fibre-wise category, and it excludes every
dynamics that moves amplitude — including the corpus's **own ledgered alternative**, S2 audit A2's
*"a real parameter `t` with a Hamiltonian"*. A criterion whose admitted set is exactly the
stipulated class restates the stipulation; it does not derive it.
*Reading (β):* an externally pre-registered physical requirement that provably forces the class is
a substantive selection, and the convention is motivated rather than arbitrary.

**I refuse to score it, and my verdict does not depend on which is right — because MOTIVATING A
STIPULATION DOES NOT CONVERT ITS ANALYTIC CONSEQUENCES INTO FINDINGS.** Even granting (β) in full,
`Z_k = SUM over classes` remains reached by substitution from the operator's form. (β) would buy the
convention a *justification*, not its consequences a *discovery*. That distinction is the whole
ruling and I state it here rather than let it be found.

**(ii) THE HOLONOMY READING.** If the corpus's object is *by nature* a closed-loop invariant, then
fibre-wise-ness is forced by what the object is, invisibility is a theorem about a genuinely defined
object, and mid-path comparison is a category error — so `T` is not a rival at all. **No lane closed
this and neither do I.** What the record shows is that the corpus is **not entitled to it for free**:
`CHOICE LEDGER A2` (`S2 audit :658`) warrants the **circuit** clock by citing that *"edge count is
carrier-supplied combinatorics (S1 `:16-22`)"* — and `S1:16-22` **is the edge list**; and `S3`'s
`CHOICE LEDGER C3` (`:981`) cites *"audit COR-F"*, whose text (`S2 audit :639-641`) reads *"edge
traversals **/** circuits"* — **a disjunction cited to close one disjunct**, with the edge unit absent
from its alternatives column. **Twice, the warrant offered for the circuit clock is a fact about
edges.** All four quotations verified at the bytes. Again: my verdict survives either reading,
because under the holonomy reading `Z_k = SUM over classes` is *still* substitution.

**(iii) THE LAYER COUNT.** W-10 says layer **TWELVE**. My brief says **THIRTEEN**. All five
cross-refutations say **FOURTEEN**. Three numbers for one round, in a program whose central
discipline is applied *by number*. Recorded; not adjudicated.

---

## 6. WHAT I STRIKE

**FROM THE REGISTRAR (`LANE_W11_CONVENTION_TEST`, seals 9/9 OK).**
1. **Its conclusion sentence.** *"Invisibility holds exactly where both branch operators are
   DIAGONAL"* is wrong in both directions (§3). The verdict it carries is right; the sentence
   carrying it is not.
2. **Leg B1 and every circuit row of leg C3 are COULD-NOT-HAVE-FAILED CONTROLS**, and the lane
   reports them as *"the invisibility theorem, working exactly as registered."* They are identities.
   I score them as zero evidence in both directions — including against the registrar.
3. **Leg B is confounded on K1 and the lane's own leg B3 says so** (`1.97e-15`): both loops have
   length 3, so the circuit convention **is** the edge convention sampled every third tick. B1 vs B2
   on K1 moves the sampling times, not the transport. **Leg C3 on B0b is the load-bearing arm** and
   it is clean — `lcm(4,3) = 12`, and *no* `n <= 20000` puts both branches at the same circuit count.
4. **`PUBLISHED_CONVENTIONS.txt` states no connection at all.** The connection used is
   `f = 57/25, c = 2 + sqrt(2)`, recoverable only from `w11_b_decisive.py:31`, and it is not the
   corpus's published generic pair. Everything reproduces at `S4:603` with my own states and seed,
   so nothing numerical moves — but a lane commissioned to settle a *convention* question published
   no convention for the parameter its effect is most sensitive to.
5. **Leg D clause 2** (*"the edge rate is not `m(P)/3`"*) is well posed on K1 and **ill-posed on
   B0b**, where the loops have different lengths and there is no single `L` to rescale by. Clause 1
   — the edge rate is **state-dependent** while the circuit rate is not — is the strong clause,
   reproduces on **12** states rather than 3 (spread `3.7e-01` K1, `2.1e-02` B0b, arms diffed at
   `min||s_i - s_j|| = 0.76 / 0.82`, `max|pi_i - pi_0| = 2.2e-16`), and is not a units question.

**FROM THE FIVE CROSS-REFUTATIONS** (all seals verify: 26, 33, 31, 16, 10 files).
6. **`LANE_W11_R_UNIQUENESS`'s criterion (b) is a ZERO-VARIABLE CONTROL WITH HARD-CODED VERDICT
   STRINGS.** Confirmed at the bytes: `t3_canonicity.py:33-34` prints the literals `root: CHANGED`
   and `(= a cube root of unity)` beside computed values that its own sealed
   `t3_canonicity.OUT.txt:13-14` records as `0.000000` and `1.000000-0.000000j`. `root_op` takes
   `np.angle(W)`, and `W` is unchanged by a `2pi` shift of any `a_e`, so the two arms are identical
   **by construction**. `LANE_W11_R_UNIQUENESS_CROSS` found this; **I confirm it independently and
   uphold the strike.** This is the **seventh** zero-variable control of record, in the leg whose
   brief was canonicity.
7. **`LANE_W11_R_STEELMAN_CROSS`'s leg E3 is not exhaustive.** It claims the gauge-covariant root
   variety on a 3-cycle is `3^3 = 27` branch choices and enumerates all 729 pairs. `Lam·T` is
   gauge-covariant (`3.5e-15`), a root (`1.9e-15`), continuous, and **commutes with `T` in 0 of 200
   draws** — so it is not among the 27, and the covariant root variety is **positive-dimensional**.
   Its *"EXACTLY ONE passes"* is not established. Its verdict is unharmed (every member of `Lam·T`
   is pi-visible), but the one leg it offers as **exact and exhaustive** is neither.
8. **`LANE_W11_R_STEELMAN_CROSS`'s name for the operative variable, "THE COMPARISON TIME", falls.**
   §3: hold the tick at `n = 1`, move only the operator, and the answer flips. Its argument is
   scoped to the root variety read at loop closure, where the transport axis is trivially inert.
9. **`LANE_W11_R_CLOCK_CROSS`'s name, "the branch operator must be FIBRE-WISE and CONSTANT ON ITS
   LOOP", falls twice.** It is **not necessary** (200/200 non-diagonal correlated pairs are blind at
   `6.3e-15`) and it is **too strong** (class-constant-but-not-loop-constant is already blind).
   Lane M's `C1` and lane I's `D2` have the right object; three of the five lanes do not.
10. **`LANE_W11_R_CLOCK`'s published connection is not the one its headline leg uses.** Confirmed at
    the bytes: `w11c_1_lattice.py:98` uses the *secondary* pair on K1 and `:112` an **undeclared**
    `uniform(0, 2pi)` draw on B0b, while its conventions file declares `(1.0, sqrt(2))` primary.
11. **THE ROUND'S CONSENSUS THAT "THE CORPUS HAS NO ADMISSIBILITY CRITERION FOR TICKS/DYNAMICS" IS
    FALSE, AND THIS IS THE SHARPEST CORRECTION I MAKE TO THE REFUTERS.** Three of five cross-lanes
    file it as a missing object, a "fourth" alongside W-08's schedules and W-09's designations.
    **`CHOICE LEDGER A1` at `:657` states one**, and I ruled with it (§2.2). What the corpus lacks is
    not the criterion but its **application**: A1's alternatives column lists *"extend by zero; act
    only at the root"* and never lists the edge tick, and A1's third clause is circular. **The defect
    is a ledger closed on an incomplete alternatives column and a circular why — not an absent rule.**

**A CHARGE I CHECKED AND DO NOT UPHOLD.** I suspected `LANE_W11_R_UNIQUENESS_CROSS`'s
`3 × edge = -2.756952` of being unreproducible, since it matches no registrar figure while its
neighbours match exactly. It is in its own sealed `x4_scope_and_defects.OUT.txt:22`, labelled
`state B2` — its own state, not the registrar's. **The figure reproduces. The charge fails.**

---

## 7. WHAT SURVIVES

1. **THE VERDICT: READING B**, on the identity of §0, which needs no rival transport, no canonical
   tick, no carrier, and no numerics.
2. **THE REGISTRAR'S LEG A, ENTIRE, AND NOW IN EXACT ARITHMETIC.** COR-F's sealed exhibit reproduces
   in a fourth implementation (`||T*T - I|| = 0.00e+00`, `W_C = -0.029200+0.999574j`,
   `diag(T rho T*) = [0.15 0.15 0.15 0.40 0.15]`). `T` is unitary (`3.5e-16`), gauge-covariant
   (`4.8e-15`), and `T^L = M_gamma` on **both** carriers (`7.2e-16`) — and **exactly, in Gaussian-
   rational `Fraction` arithmetic with zero rounding, at `L = 3` and `L = 4`.**
3. **THE REGISTRAR'S LEGS B2, C3 AND D-CLAUSE-1**, reproduced from independent code, with **my** seed,
   **my** 40-state and 12-state families, and **the corpus's own published generic connection** rather
   than the registrar's. K1 edge spread `7.5e-01` against circuit `5.6e-16`; B0b edge `8.9e-01`
   against circuit `4.4e-16`; both return to `<= 5.6e-16` exactly on `n = 0 mod lcm(L_F, L_C)`.
   **The registrar's picture is not an artefact of its unpublished connection or of its three
   hand-built states.**
4. **THE BICONDITIONAL** (§3) — 9000 cells, 0 disagreements, both outcomes in quantity. This is the
   complete answer to §0.4.1 and it is one line of algebra.
5. **THE UNIQUENESS THEOREM** (§2.2) — the admissible local family is exactly `Lam·T`, and **all of
   it is pi-visible** (400 draws, 800 root checks, 0 blind, smallest spread `5.0e-01`).
6. **THE DIVISOR STRUCTURE** `{n : U^n diagonal} = d·Z`, `d | L`, reproduced independently
   (`LANE_W11_R_CLOCK_CROSS`'s theorem; 1200 draws per loop per carrier, every observed `d` divides
   `L`).
7. **N1's VALUE AND ITS CROSS-CHECKS.** `m(0.4 + 0.3x + 0.3y) = -0.767507880357` by Jensen at `2^20`,
   cross-checked by a Jensen-free 2-D quadrature; W-10 N-3's `log(4/9) = -0.810930216216` reproduced
   exactly; the order-4 rate is `-0.804718956`, which is COR-K's number and equals `-(1/4)log 25`.
8. **THE DOCUMENTARY FINDINGS, ALL VERIFIED AT THE BYTES BY ME:** `S1:52-53`'s edge-wise definition
   of parallel transport, which the corpus's own `M_gamma` violates; `S1:16-22` as the edge list;
   `S2 audit :657` (A1) and `:658` (A2); `S2 audit :639-641` (COR-F's disjunction); `S3:981` (C3);
   `S3 audit :794` (COR-F) and `:798` (COR-J, *"Add it to the CHOICE LEDGER"* — never added);
   `S4:575` and `S4:603`; `REGISTER:409`, `:577`.
9. **EVERY W-11 LANE REPRODUCES AT THE BYTES.** All eleven `SEALS.sha256` verify: 9, 16, 16, 18, 16,
   12, 26, 33, 31, 16, 10 files, **0 failures**. **The arithmetic was never the problem in this round
   either.** Every defect above is a defect of naming, scope, control design, or reporting.

---

## 8. NEW DEFECTS OF RECORD

- **W11-1.** The registrar's conclusion sentence names **diagonality**, which is neither necessary
  nor sufficient. §6.1, §3.
- **W11-2.** `Q_n = (branch_F^n)^* branch_C^n` — the object every functional-layer result is about —
  **has no register row and no name in the corpus.** The class-constant-diagonal biconditional
  should be entered against W-10's REOPENS clause on COR-F.
- **W11-3.** **`CHOICE LEDGER A1` (`S2 audit :657`) is closed on an incomplete alternatives column
  and a circular why-clause**, and it is the entry §0.4.1 should have named — not A2. Its clause (c)
  excludes the edge tick only by presupposing fibre-wise-ness. Its own flag `G2` (`:686`) said an S3
  lane should attack it first; none did, through nine subsequent rounds.
- **W11-4.** **The seventh zero-variable control of record**, with hard-coded verdict strings beside
  numbers that refute them (§6.6). W-08's isolation audit named this the commonest fatal defect and
  said a ledger cannot catch it; it is now committed inside a lane whose brief quoted that warning.
- **W11-5.** **The corpus's `M_gamma` violates `S1:52-53`'s own definition of parallel transport** —
  it moves no fibre value — and no sealed artifact or register row notices that the operator the
  entire functional layer runs on is excluded by the carrier page's definition of the thing it is
  named after. This is why `S1:52-53` cannot be used to exclude `D` without also excluding the spine.
- **W11-6.** **N1 as registered carries neither of its hypotheses** (§4), although the corpus's own
  ERRATUM AGAINST W-02 is exactly the equidistribution correction.
- **W11-7.** **`COR-J` still has no register row** — one occurrence in `REGISTER_V001.md`, inside
  W-05's narrative; *"gauge covariance"* occurs zero times — and this round establishes that it is
  **not** the criterion doing the work: `T` is gauge-covariant and `Z^T` gauge-invariant, so COR-J
  excludes nothing here. Two lanes assign it decisive force it does not have.
- **W11-8.** **Three different layer numbers for one round** (§5.iii).

---

## 9. SELF-FLAGS

- **LINEAGE, AND WHAT IT COSTS.** I am Opus 5. So is the registrar, so are all five refuters, so are
  all five cross-refuters, so is every row from W-07. The last lineage-independent boundary in this
  corpus is Fable 5 → Opus 5, and **everything in W-07 through W-11 sits on one side of it**. My
  agreement with all ten lanes on the verdict is therefore **worth less than my disagreements with
  them on the evidence** — agreement across one lineage is the cheapest thing this program produces.
  **Discount this ruling as one block with W-07 through W-11.** The likeliest failure mode is the
  program's recurring one, and I have now committed it once myself (below).
- **TWO SELF-DEFECTS, IN THE SEALED CODE, NOT PATCHED OUT.** (i) My first uniqueness search was a
  numerical-gradient descent that returned *"0 of 24 random starts converged"* — a null about the
  **search**, not the variety — and an Adam retry that **diverged** to loss `~9e+16`. I caught it
  only because I ran the positive control (COR-F's own `T`, loss `1.0e-30`) **first**; without that
  rule I would have reported the null as a classification. (ii) My first Fourier test asserted a
  support claim about the whole matrix that is true only of its **diagonal**, and printed `6.13`
  against its own sentence. Both are documented in `r4_uniqueness.py`'s header. A third, smaller
  one: `r3_canonicity.py` printed `||U^L - M||` for the `M_gamma` row, where the quantity is
  meaningless; the label is corrected and the correction is recorded here.
- **MY UNIQUENESS THEOREM RESTS ON A *STRICT* READING OF A1(b) THAT IS MINE, NOT THE CORPUS'S.**
  Read strictly ("no adjoined constants whatever"), `T` is exactly canonical. Read loosely, the
  `(L-1)`-torus is admissible too. **Nothing in the verdict moves** — the whole torus is pi-visible
  (`0/400`) — but a reader who thinks canonicity **is** the question should read §2.2 as resting on
  my reading of one clause. This is the place to attack me.
- **THE ONE THING I DID NOT DO, AND IT IS THE SAME GAP THE ROUND INHERITED.** I did not run a
  dynamics that is **not** an `L`-th root of `M_gamma` — in particular not the corpus's own ledgered
  alternative, S2 audit A2's *"a real parameter `t` with a Hamiltonian"* (the magnetic Laplacian) at
  generic connections. Every "rival" tested by every lane in this round, mine included, satisfies
  `U^L = M_gamma` — **the round explored the root variety of the very operator whose form is the
  question.** My verdict does not need it (the identity of §0 needs no rival at all), but the
  *strength* claims of §2.4 are scoped to that variety.
- **RANK ONE, `U(1)`, TWO LOOPS, TWO CARRIERS.** Nothing here touches charge (W-03's THEOREM S4-1
  failure), `SU(2)`, a third loop, or the C\*-algebraic layer. At higher rank "fibre-wise" and
  "diagonal" come apart and §3's identification would have to be redone.
- **MY B0b IS A RECONSTRUCTION**, as every lane's is: S4's carrier code is not in the repo (W-10
  N-10). I re-derived the walks, the incidence and the class multiset in my own code and checked the
  multiset against `S4:575`, but the grid index functions are the same shape every lane uses.
  **Narrowed, not removed** — the defect W-10 recorded against two of its own refuters.
- **I HAVE NOW MISNAMED NOTHING, WHICH IS EXACTLY WHAT THE PREVIOUS SIX LAYERS BELIEVED.** My name
  (§3) is a two-clause extension of `REGISTER:577`, guarded by `grep` before use. The first thing a
  refuter should attack is whether `Q` is the right object at all, and specifically whether my
  "class-constant diagonal" collapses at higher rank or on a carrier with a class of size `>= 2` in
  **both** loops — `LANE_W11_R_MATH_CROSS`'s `SHARE2` is that carrier and I did **not** rebuild it.

---

## 10. WHAT THIS RULING DOES AND DOES NOT DECIDE

**IT DOES NOT TAKE ANY OF THE PRINCIPAL'S THREE DECISIONS.** It answers the question the second
deferral was for, and nothing else.

**WHAT IT HANDS THE DECISION.**
- **Accept STOP-FALLS-REBUILD, or hold at STOP.** Reading B says the functional layer's portability
  is **what you get for free when the carrier was stipulated out** — so a rebuild inherits no
  portable physics from it, only a correct piece of mathematics. W-10 posed this as the axis with
  *"opposite consequences for a rebuild"*, and the ruling comes down on the side that gives a rebuild
  **less**, not more, to carry forward. **That bears on the first decision; it does not take it.**
- **Publish the Mahler note.** **YES — and Reading B does not withdraw it.** N1's publishable content
  is untouched, because the Mahler identification is a theorem about `(pi, characters)` that survives
  changing the operator inside the fibre-wise class. What Reading B withdraws is the claim that N1's
  blindness to the complex is a **discovery about physics**. Publish it with Hypotheses 1 and 2
  stated as in §4 — Hypothesis 1 in particular, which converts the weakest sentence in the corpus's
  framing into the paper's own honest scope line.
- **Rename the project.** Untouched by this ruling.

**REOPENS IF:** a dynamics that is **not** an `L`-th root of `M_gamma` — the corpus's own ledgered
Hamiltonian alternative is the obvious one — is run at generic connections and found pi-blind ·
or an admissibility criterion is written from the carrier's data that admits `M_gamma` and excludes
COR-F's `T` **without** presupposing fibre-wise-ness (A1 clause (c) does not qualify) · or a
**lineage-independent** lane reads this ruling.

