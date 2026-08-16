# S2 AUDIT — ADVERSARIAL, DEFAULT-REFUTE — V001 — 2026-08-16

**OVERALL: REFUTED.**

Target: `S2_FORMATION_CONDITION_ON_K1_V001.md`
sha256 `248ce856efaef157c68e818dde589d0200bbc1dd9fd9fc1fcc8cdc7bc88734d9` — verified at bytes,
both sidecars agree.
Carrier: `S1_CARRIER_K1_V001.md` sha256 `3eb70375bfd0900e4dd56cae294fa31b3b6e19cf6634853501fab5ffcebd92ac`.
Custody: `CUSTODY_V001.md` sha256 `6629ae3b5aedbafda4e56b2a743b84a888949c59f988296056a2a7abeb20aa49`.
Design: `FOUNDING_DESIGN_V001.md` sha256 `ca25b79c76531d909e75fcb58163ad6456b10f086d59638f47224210d19b13e9`.
Inheritance: `INHERITED_FROM_THE_PREDECESSOR.md` sha256 `96c1d305c7a4eeceab8dc78f971832eb0c2e34980a72d1c12c5ece5cd4ad0079`.

This is an audit of **mathematics**, not of citations. Every number below was recomputed from the
definitions in S1 at my own hand. Testimony carried no weight (custody §4, `:44–53`).

---

## 0. THE VERDICT, STATED FIRST

The build's headline is: *a formation condition of the predecessor's character cannot be written
on K1, and the obstruction is threefold and proved* (target `:22–28`).

**Two of the three obstructions are false, and the third is overstated.** The false ones are
`O1` (rank) and `O3` (split) — the two the build called its sharpest results. I exhibit, on
**K1's own published connection** (S1 `:114–115`), two evolutions of a common ready state that
become **exactly orthogonal**, with **zero addition** to the carrier.

```
GRADE BY DIMENSION
  D1  state space          CONFIRMED-WITH-CORRECTIONS
  D2  evolutions           REFUTED      (undeclared change of state space between Q1 and Q2)
  D3  orthogonality        REFUTED      (O1 false, O3 false, Thm 3's exhibit non-reproducible)
  D4  minimum addition     REFUTED      (not minimal — zero addition suffices)
  D5  split test           CONFIRMED-WITH-CORRECTIONS
  D6  custody              CONFIRMED-WITH-CORRECTIONS
OVERALL                    REFUTED
```

**What survives, and it is not nothing.** S1's counts re-derive exactly. The generator argument of
Q2 is correct. Every number in Q4's Laplacian analysis reproduces to six digits, including the
root bound `0.621268` at `(pi, pi)` — that is real work and it stands. The pointer discipline is
the best I have audited in this corpus: all sixteen `file:line` pointers resolve correctly. The
FLAG BLOCK is honest and F2 anticipates part of COR-B.

**What is refuted is the page's three governing claims**, and the refutation is constructive: the
condition the build declared unwritable is written below, in §4, on K1's own data.

---

## 1. THE UNIVERSE SWEPT, DECLARED

Custody §3 (`:34–42`) requires the universe declared and the count reported.

```
Corpus swept: /Users/bgm/MB Work/where-atoms-come-from/*.md      hits = 6 files
  CUSTODY_V001.md · FOUNDING_DESIGN_V001.md · INHERITED_FROM_THE_PREDECESSOR.md
  PUBLICATION_V001.md · REGISTER_V001.md · S1_CARRIER_K1_V001.md
  + the target S2_FORMATION_CONDITION_ON_K1_V001.md
Opened in full: 6 of 6, plus the target.  Not opened: PUBLICATION_V001.md, REGISTER_V001.md
  (neither is cited by the target and neither carries a governing clause bearing on it).
Digests recomputed: 5 of 5 cited digests. All 5 match. Sidecars: 2 of 2 match.
Predecessor corpus: NOT swept. Set aside (INHERITED :6-9). No finding here depends on it.
```

---

## 2. D1 — THE STATE SPACE — **CONFIRMED-WITH-CORRECTIONS**

### 2.1 What reproduces

`Gamma(L) = C^5` as a direct sum over vertices, with `<s,t> = sum_v conj(s(v)) t(v)`: correct, and
it is the inner product the fibres carry (S1 `:48`) summed with counting measure. Gauge acts
fibrewise and unitarily: correct.

Both saturated counts recomputed by ranking the real Jacobian of the gauge orbit:

```
VECTORS   dim_R Gamma(L) = 10 ; generic orbit rank = 5 ; invariants = 5     CONFIRMED
          (and the orbit rank drops to 4 when one component vanishes — the action is not free
           there. The build says "freely where all s(v) != 0" and is right to say it.)
RAYS      dim_R CP^4 = 8 ; gauge through U(1)^5/U(1) = 4 ; invariants = 4    CONFIRMED
```

S1's own count re-derived independently: the gauge generator matrix `M = d1^T` (6x5) has
**rank 4**, so `6 - 4 = 2` invariants, and `w_F = (1,1,1,0,0,0)`, `w_C = (0,0,0,1,1,1)` both lie in
its left null space. The build's warning that the three effective dimensions differ (5 / 4 / 4) is
correct and worth the space it takes.

### 2.2 The correction — COR-E

Target `:108–110` states as a **derived result**:

> *the complete gauge-invariant content of a normalised state on K1 is a probability distribution
> over its five vertices, and nothing else.* All phase information in a state is pure gauge. Only
> the connection's phases (§3) are physical.

**The first sentence is true of a state considered alone. The third sentence is false**, and S2 is
precisely the stage at which it matters, because a formation condition is by construction a
statement about a state **and** a connection together.

Recomputed on the joint configuration space:

```
joint dim_R  =  6 (connection)  +  10 (state)   =  16
gauge orbit rank at a generic point             =   5        (recomputed, not assumed)
JOINT GAUGE INVARIANTS                          =  11
build's exhibited invariants: 5 moduli + 2 (W_F, W_C)  =  7        MISSING: 4
```

The eleven are exhibited: gauge-fix `g_v = exp(-i arg s(v))` where `s(v) != 0`, which makes the
section real and positive and leaves **no** residual gauge freedom. Then all five moduli `|s(v)|`
and all **six** edge phases `a_e + arg s(u) - arg s(v)` are invariants. Verified invariant under
random gauge to `1e-9`. `5 + 6 = 11 = 16 - 5`.

**Consequence.** A nowhere-vanishing section is a *trivialisation*, and against it the connection
has six meaningful phases, not two. The build's own §1.4 conclusion foreclosed a space of four
further gauge-invariant, connection-dependent quantities before Q3 ever looked for a condition to
write. This is not the reason Q3 fails — §4 below writes a condition using only the seven the
build did exhibit — but the count as stated is wrong and it narrowed the search.

---

## 3. D2 — THE EVOLUTIONS — **REFUTED**

### 3.1 What reproduces, and it is correct

```
rank d1 = 4 , nullity 2                                    CONFIRMED (computed)
free rank of the based-loop group = E - V + 1 = 2          CONFIRMED
d1 @ d2 = 0  (d^2 = 0 on this complex)                     CONFIRMED
chi = 0, b1 = 1                                            CONFIRMED
pi_1(K1, v0) = Z, not Z^2 — the face kills the filled class CONFIRMED
```

So every closed transport at `v0` is `W_F^m W_C^n`. The build's claim that the two triangles are
not a selection but a generating set is **correct**, and its structural note that one is a
curvature and one is flat is correct and is the best paragraph on the page.

**And the build did not conflate parallel transport with time evolution.** The brief instructed me
to check for exactly that. It did not happen: §2.2 explicitly refuses the Laplacian *because*
`exp(-i t Delta_A)` needs a `t` that K1 does not supply, and says so in order not to "smuggle O2
in unannounced." That refusal is correct and to the build's credit.

### 3.2 The refutation — COR-A's mechanism

**Q1 derives the state space to be `Gamma(L) = C^5`. Q2 then places the evolutions on
`L_v0 = C`. The change of space is never declared, and it is the load-bearing move of the entire
page.**

The build's own CHOICE LEDGER convicts it. Entry C2 (target `:593`) reads:

> C2 | States = sections of the line bundle, `(+)_v L_v` | *alternatives:* states at the root only
> (`L_v0`, dim 1); … | **closed** — root-only states are the dim-1 special case, covered by
> Theorem 1(c)

The ledger considers `L_v0` as an **alternative and rejects it**. Theorem 1 is then proved *on the
rejected alternative* and applied as a global obstruction to the *chosen* space. That is backwards.
`L_v0` is not "the dim-1 special case covered by Theorem 1(c)"; it is the only space on which
Theorem 1 has content.

The build half-notices at `:292–294`: *"Theorem 1 forbids orthogonality within a fibre. It says
nothing about `Gamma(L) = C^5`."* Having said it, it asks only whether the connection can carve
orthogonal **subspaces**, and never asks the obvious next question — **do the two transports
themselves act on `Gamma(L)`?** They do, canonically, and §4 does it.

### 3.3 The second correction — COR-F

O2 is stated at `:26` as `K1 carries no time parameter`. What is proved at `:283–288` is that K1
carries no *continuous* time. K1 does carry a discrete one: **the number of edges traversed**.
Path length is combinatorial data of the carrier's own cells (S1 `:16–22`), not an import. The
build's sentence *"transport is path-indexed, not time-indexed"* treats those as exclusive; a path
index **is** a discrete time index, monotone and ordered. §4.4 below solves the gate's *"at a time
solving an equation"* clause with `n` = circuit count and no added parameter.

---

## 4. D3 — THE ORTHOGONALITY RESULT — **REFUTED** — AND MY INDEPENDENT Q3

### 4.1 Theorem 1 is true, trivial, and mis-scoped

Reproduced at the build's own sample size:

```
200,000 independent draws of (a_1..a_6, z):
max | |<T_F z, T_C z>| - |z|^2 |  =  1.0658e-14
```

matching the build's `1.07e-14`. As a statement about `L_v0 = C` it is correct and it is one line:
`U(1)` is closed under multiplication, and `C` has no non-zero orthogonal pairs. Fine.

**What is false is O1 as a governing claim.** Target `:25` and `:218–221`:

> `O1 RANK rank-one fibres => transport is a scalar => rays are fixed => never orthogonal`
> … *"A gate that fires when two evolutions of a common ready state become orthogonal requires
> those evolutions to move the ray and to move it into a complement that exists. On K1 neither
> holds: nothing moves the ray, and no complement exists."*

Rank-one fibres make transport a scalar **on each fibre**. The collection of those fibre-scalars is
a **diagonal, non-scalar** operator on the state space Q1 derived — and diagonal non-scalar
unitaries move rays and produce orthogonality in abundance.

### 4.2 CONSTRUCTION A — the object the build missed

**Definition, on K1's own data and nothing else.** For a closed edge path `gamma` based at `v0`,
let `V(gamma)` be its vertex set and `W(gamma) in U(1)` its holonomy. Define

```
M_gamma : Gamma(L) -> Gamma(L)
  (M_gamma s)(v)  =  W(gamma) . s(v)      for v in V(gamma)
  (M_gamma s)(v)  =  s(v)                 otherwise
```

In words: **transport each of the loop's own vertices once around the loop; leave the rest alone.**
This is parallel transport — the *same* operation the build used — applied at every vertex the loop
passes through rather than only at the root. It uses fibres, edges, orientation and the connection.
No time, no rank, no split, no weight, no measure.

Well-defined because `U(1)` is abelian, so the holonomy based at any vertex of `gamma` equals
`W(gamma)`. Verified by explicit composition:

```
hol(dF) based at v0, v1, v2  =  (-1+0j), (-1+0j), (-1+0j)      all equal W_F
hol(c)  based at v0, v3, v4  =  (-0-1j), (-0-1j), (-0-1j)      all equal W_C
```

Applied to K1's two loops, with `V(dF) = {v0,v1,v2}` and `V(c) = {v0,v3,v4}` (S1 `:19–25`):

```
M_dF = diag( W_F, W_F, W_F,  1,   1  )
M_c  = diag( W_C,  1,   1,  W_C, W_C )
```

Both **unitary**. Neither **scalar** — verified. And therefore Theorem 1 does not apply to them.

### 4.3 THEOREM A — the formation condition, written

For a normalised section `s` with probability vector `p_v = |s(v)|^2` — which is, by the build's
own Q1, the complete gauge-invariant content of the ready state:

```
<M_dF s, M_c s>  =  conj(W_F) W_C . p_0  +  conj(W_F) . (p_1 + p_2)  +  W_C . (p_3 + p_4)
```

Verified against direct matrix action to `1e-10`. Every ingredient is gauge-invariant (the two
holonomies by S1 `:80–82`, the moduli by §2.1), so the whole is. Checked under random gauge: exact.

> **THEOREM A (firing criterion).** There exists a normalised `s` with `<M_dF s, M_c s> = 0` if and
> only if `0` lies in the convex hull of the three unit-modulus numbers
> `{ e^{i(c-f)}, e^{-if}, e^{ic} }`, equivalently iff the three angles `{0, c, c+f} mod 2pi` have
> every consecutive gap `<= pi`. Here `f = a1+a2+a3`, `c = a4+a5+a6`.

**Proof.** The three coefficients are unit-modulus and the three weights `p_0`, `p_1+p_2`,
`p_3+p_4` range over exactly the 2-simplex as `s` ranges over normalised sections. A convex
combination of three unit vectors vanishes iff `0` is in their convex hull; for points on a circle
that holds iff no open half-plane contains all three, i.e. iff the maximal arc gap is at most
`pi`. ∎

**Checked, not asserted:** the convex-hull criterion against brute-force minimisation of
`|<M_dF s, M_c s>|` over a 120-step triangulation of the simplex, at **1369** grid points of
`(f, c)` — **0 mismatches**.

### 4.4 The exhibit — on S1's own published connection, exact

S1's worked instance (S1 `:114–115`) is `a_1=a_2=a_3=pi/3`, `a_4=a_5=a_6=pi/2`, giving
`W_F = -1`, `W_C = -i`. Take the ready state `p = (1/2, 0, 0, 1/4, 1/4)`:

```
<M_dF s, M_c s>  =  (i)(1/2)  +  (-1)(0)  +  (-i)(1/4 + 1/4)  =  i/2 - i/2  =  0     EXACTLY
computed: |<M_dF s, M_c s>| = 1.27e-16
gauge-invariant under 5 random gauge transformations, max deviation < 1e-12
```

**Two evolutions of a common ready state, both parallel transport around K1's own two loops,
become exactly orthogonal on K1's own published connection.** That is the predecessor's gate minus
its time clause, written on the carrier's own data with nothing added.

Sharper still, on the same instance:

```
|<s, M_dF s>|  =  |-1/2 + 1/4 + 1/4|  =  0
```

A **single** loop transport carries the ready state to a state orthogonal to it. The build's
sentence *"nothing moves the ray"* (`:219–220`) is false on K1's own worked connection.

### 4.5 The four properties the build required, checked against Construction A

```
TRIVIAL-CONNECTION LIMIT   W_F = W_C = 1 -> all three coefficients = 1 -> <.,.> = 1 for every
                           ready state. NO FORMATION.  Correct trivial answer.        PASSES
CRITERION C4 (build's own)  the verdict is a non-constant function of (W_F, W_C) — the firing
                           region is a proper subset of the torus (Thm A).            PASSES
GAUGE INVARIANCE           every factor is a gauge invariant; checked numerically.     PASSES
NO SPLIT                   no tensor factorisation, no partial trace, no subsystem. Only the
                           vertex-indexed direct sum (K1's) and the two loop vertex-sets (K1's).
                                                                                       PASSES
```

**The root can never fire.** For `p = delta_{v0}`, `<M_dF s, M_c s> = conj(W_F) W_C`, of modulus
`1` for every connection. K1's own designated origin is the wrong place to start — reproducing the
build's §4.2 root result **by a completely independent route**. That agreement is the one
cross-corroboration in this audit and it is worth recording: two unrelated constructions both find
that `v0`'s joint membership in both loops disqualifies it as a ready state.

### 4.6 Construction A closes the build's F1 — its own named gap for S3

The build's F1 (`:618–623`) says Construction 1 is symmetric in `W_F`/`W_C`, cannot tell a
curvature from a flat holonomy, and that *"that is the gap S3 should close first."*

`|<M_dF s, M_c s>|` is invariant under `f <-> c` **iff** `p_1 + p_2 = p_3 + p_4`. Exhibited:

```
ready state p = (p_v0, p_v1+p_v2, p_v3+p_v4) = (0.30, 0.45, 0.25)
  at (f, c) = (1.221730, 2.591814):  |<M_dF s, M_c s>| = 0.001828     FIRES
  at (c, f) = (2.591814, 1.221730):  |<M_dF s, M_c s>| = 0.378633     DOES NOT FIRE

max asymmetry over 400 random (f,c),  q = r = 0.35   :  0.000000
max asymmetry over 400 random (f,c),  q=0.45, r=0.25 :  0.398988
```

For any ready state that does not weight the two triangles equally, the condition **distinguishes
the curvature from the flat holonomy**. F1 is closed, at S2, at zero cost.

### 4.7 THEOREM 2 / O3 — **REFUTED**

Target `:369–373`:

> `Gamma(L) = C^5` admits NO non-trivial tensor factorisation. Hence on K1 there is no subsystem,
> no environment, no partial trace, and no reduced state — not as an omission, but as an
> impossibility.

The arithmetic is right and the conclusion does not follow. `5` is prime, so `C^5` has no
non-trivial tensor factorisation **of the whole space**. Everything after "Hence" is false.

**(a) Reduced states exist, and the build computed one.** Take the vertex algebra
`A = C^5 subset M_5` — the diagonal, which is exactly K1's own vertex indexing. The pinching
`E(X) = diag(X)` is a conditional expectation: idempotent, positive, unital, trace-preserving —
all four verified. The reduced state of `psi` on `A` is:

```
reduced state on A          =  [0.431595 0.025621 0.052237 0.450998 0.039548]
build's own Q1 invariant p_v =  [0.431595 0.025621 0.052237 0.450998 0.039548]      IDENTICAL
```

**The build's §1.4 "derived result" — that the gauge-invariant content of a state is the
probability vector — *is* a reduced state.** Q1 computed the object Q3 declares impossible. The
build's supporting sentence *"direct sums do not admit partial traces"* (`:392`) is false: the
pinching onto a block-diagonal subalgebra is a partial trace in the algebraic sense, and it is
canonical.

**(b) A three-way split exists, and it is K1's own.** Verified as a conditional expectation
(idempotent, positive, trace 1):

```
Gamma(L)  =  C_{v0}  (+)  C^2_{v1,v2}  (+)  C^2_{v3,v4}          1 + 2 + 2 = 5
algebra M_1 (+) M_2 (+) M_2 in M_5 ; commutant C (+) C (+) C ; block pinching is the expectation
```

Blocks two and three are **exactly the supports of the build's own `P_C` and `P_F`** — the row
norms it printed at `:352–353` are `[0,1,1,0,0]` and `[0,0,0,1,1]`, which I reproduce exactly. So
a three-way decomposition of the shape the predecessor stipulated is not merely available on K1,
it is the decomposition the build's own Construction 1 produces.

**(c) A genuine tensor factorisation exists on a K1-designated subspace.** The non-root vertices
carry two independent K1-supplied indices — *which triangle* (S1 `:24–25`) and *which step from
the root along the oriented loop* (S1 `:14`, `:19–21`):

```
v1 = (F,1)   v2 = (F,2)   v3 = (C,1)   v4 = (C,2)
=>  C^4 = C^2 (x) C^2      and      Gamma(L) = C  (+)  ( C^2 (x) C^2 )
partial trace over "position" of a random non-root state: [0.136977 0.863023], trace 1.0
```

The claim that a three-way split "needs `dim >= 8`" (`:383`) is true only for a *tensor* product of
three non-trivial factors. It is false for subsystems in the sense that carries partial traces and
reduced states — which is the sense the design's own inductive-limit template uses
(`FOUNDING_DESIGN_V001.md:62–66`, an inductive limit of finite **record algebras**).

**What survives** is the build's own F2, and only that: `C^5` admits no non-trivial tensor
factorisation of the whole space, and K1 designates no **canonical** bipartition of the whole
space. That is a real observation. It is not "impossible", and Q5's boast that the negative answer
comes "by arithmetic" rather than restraint (`:36–38`) does not stand.

### 4.8 THEOREM 3 is true; its exhibit is not — COR-D

The Lemma and Theorem 3 are correct. I reproduce the theorem itself exactly:

```
20,000 random connections, W_F=1 and W_C=1 each forced in half the draws
theorem-vs-computation mismatches:  0                                   CONFIRMED
P_F row norms [0,0,0,1,1] ; P_C row norms [0,1,1,0,0]                   CONFIRMED
```

**But three of the four rows of the exhibited table at `:335–340` are wrong.** Recomputed with
orthonormal bases from the null space, so that `max|<P_F,P_C>|` is the cosine of the smallest
principal angle — the only basis-independent meaning the quantity has:

```
   case              dim P_F  dim P_C   build printed   CORRECT VALUE
W_F!=1, W_C!=1          2        2        0.000000        0.000000     agrees
W_F =1, W_C!=1          3        2        0.803456        0.816497 = sqrt(2/3)
W_F!=1, W_C =1          2        3        0.430039        0.816497 = sqrt(2/3)
W_F =1, W_C =1          3        3        0.607040        1.000000  exactly
```

Two independent defects.

**First, an internal contradiction.** The build's own proof at `:327–328` reads: *"If
`W_F = W_C = 1`, both have dimension 3 and `3 + 3 = 6 > 5`, so they intersect non-trivially."* A
non-trivial intersection forces the largest principal cosine to be **exactly 1**. I confirm the
intersection is 1-dimensional (the globally flat section) and the principal cosines at the trivial
connection are `[1.0, 0.6667, 0.0]`. The build's own theorem forbids the number its own table
prints.

**Second, a broken symmetry.** The two middle cases are mirror images of one another under
exchanging the roles of the triangles — both are a 3-dimensional loop-parallel space against a
2-dimensional delta span, and both give `sqrt(2/3)`. The build printed `0.803456` and `0.430039`,
which are not equal to each other and not equal to anything. They are artefacts of non-orthonormal
spanning vectors. Custody's "exhibit, do not assert" is not satisfied by numbers that cannot be
recomputed.

The theorem stands; the exhibit must be replaced.

### 4.9 MY INDEPENDENT Q3 RESULT

> **A formation condition of the predecessor's character CAN be written on K1's own data, with
> zero addition to the carrier.** It is `<M_dF s, M_c s> = 0`: two parallel transports around K1's
> own two loops, applied to a common ready state in the state space Q1 itself derived, become
> orthogonal exactly on the region of the invariant torus given by Theorem A. It is
> connection-dependent, gauge-invariant, has the correct trivial-connection limit, fires exactly
> on S1's own published connection, distinguishes the curvature from the flat holonomy, and uses
> no split of any kind.
>
> **`O1` is false. `O3` is false. `O2` survives only in the weakened form "K1 supplies no
> *continuous* time" — the circuit count `n` is carrier-supplied and solves the gate's "at a time"
> clause.**
>
> **What genuinely obstructs S2 is neither rank nor split nor time. It is DURABILITY**, exactly as
> `FOUNDING_DESIGN_V001.md:51–66` promoted it. Exhibited on Construction A with `n` circuits, at
> `f = 2.0`, `c = 1.1`, `p = (0.4, 0.15, 0.15, 0.15, 0.15)`:
>
> ```
> |<M_dF^n s, M_c^n s>|  :  min over n <= 400  =  0.024654  at n = 42
>                           max over n >= 1    =  0.999941
> ```
>
> **It fires, and then it un-fires.** A reversible write, not a record — reached at zero cost
> rather than at the cost of one real parameter and an open Hamiltonian choice. That is the same
> conclusion the build reached in §4.2, and it is the only one of its three obstructions that was
> ever load-bearing.

---

## 5. D4 — THE MINIMUM ADDITION — **REFUTED**

### 5.1 Every number in §4.2 reproduces. This is the strongest part of the build.

```
degrees (v0..v4) = (4,2,2,2,2), sum 12 = 2E                              CONFIRMED
Delta_A Hermitian                                                        CONFIRMED
spectrum gauge-invariant, and a function of (f,c) alone                  CONFIRMED
min over the 241x241 torus grid of p_max at v0 = 0.621268 at (pi,pi)     CONFIRMED exactly
=> bound |A(t)| >= 2 p_max - 1 = 0.242536                                CONFIRMED
independent search, 40,000 random (f,c,t) at v0: min |A| = 0.245612      bound holds, tight
min p_max at v1 = 0.250000  -> bound vacuous                             CONFIRMED
build's firing point f=3.663319822, c=2.194049746, t=51.456642361:
   |A| = 1.0458e-06                                                      CONFIRMED
   spectrum [0.695897, 1.0, 1.802623, 3.0, 5.501479]                     CONFIRMED
   weights  [0.117930, 0.378913, 0.177096, 0.229476, 0.096584]           CONFIRMED
```

The inequality `|A(t)| >= 2 p_max - 1` is correct (reverse triangle inequality on weights summing
to one), and I confirm the degenerate-eigenvalue grouping does not change it. §4.3's rank-2 exhibit
is correct: `T_C (1,0) = (0,1)`, orthogonal; and rank-2-with-scalar-connection indeed still gives
overlap 1. The parameter arithmetic `6 x 4 = 24` against `5 x 4 = 20` is right (and for the record,
the correct residual count is `24 - 19 = 5`, since the constant central `U(1)` acts trivially — the
build wisely did not assert a number).

### 5.2 The refutation — COR-C

Target `:528`: *"The minimum that answers the brief's question is (i): one real parameter."*

**Zero parameters answer it.** Construction A (§4) writes a firing, connection-dependent,
gauge-invariant formation condition with:

- no real `t`,
- no Hamiltonian choice — so the build's **C3, which it left OPEN, does not arise at all**,
- no rank increase,
- no split,
- and a *"time"* that is the carrier's own circuit count.

And Construction A reaches the *same* qualitative verdict as addition (i) — fires, then recurs.
So the build paid one real parameter plus an open modelling choice for a conclusion available for
nothing. Addition (i) is not the minimum; it is not even necessary.

The build's own ranking table at `:520–526` should read: **CHEAPEST — add nothing.**

### 5.3 One numerical correction — COR-H

Target `:465`: `sup |<psi(0), psi(t')>| for t' > t* + 5 = 0.994373`. Over a longer window
(400,000 samples out to `t* + 5000`) I get **0.999793**. The build's figure understates the
recurrence, which strengthens rather than weakens its verdict, but the displayed number is a
window artefact and should be labelled as a lower bound on the supremum, not the supremum.

### 5.4 The cost claims, checked

- *"A split? No"* for addition (i) — **correct**. Orthogonality of two vectors in `C^5` uses no
  factorisation. I found none.
- *"Raising the rank to 2 hands you a bipartition for free"* — **correct**, `C^10 = C^2 (x) C^5`,
  and correctly qualified as non-canonical and gauge-dependent. This is the most careful analysis
  on the page.
- *"more vertices — NOT A FIX"* — **incorrect as reasoning**, though the conclusion is unaffected
  by anything here. The stated reason is that Theorem 1 is rank-based; but per §4.7, more vertices
  *would* change what Theorem 2 forbids (`C^6 = C^2 (x) C^3`), which the build itself concedes at
  F2. The row conflates the two obstructions.

---

## 6. D5 — THE SPLIT TEST — **CONFIRMED-WITH-CORRECTIONS**

I scanned the whole artifact for a smuggled decomposition: a preferred basis, a chosen subsystem, a
tensor factorisation, a marked vertex used as a system boundary.

**The literal answer survives. I found no imported source/record/environment split, no stipulated
subsystem, and no partial trace taken.** Specifically:

- `Gamma(L) = (+)_v L_v` — K1's own (S1 `:48`). Confirmed as a genuine direct sum.
- The root `v0` — K1's own (S1 `:14`). It is used as a **base point** for closed transports, which
  is what a root is for. It is **not** used as a system boundary. Checked: no result on the page
  partitions the state space at `v0`.
- The two loops — K1's own, and the build's proof that they *generate* is correct, so they are not
  a selection. Confirmed.
- Unit weights (C1) — declared, and the claim that every result is invariant under positive
  reweighting is **correct**, since §3.4's supports are disjoint for any `w_v > 0`.
- The vertex basis in `span{delta_v3, delta_v4}` — not a smuggle: those subspaces are
  characterised support-wise, which is basis-free.

**Three corrections.**

**(1) The justification is refuted, not the answer.** Target `:539–540`: *"By Theorem 2 no such
split exists on `Gamma(L) = C^5`. I did not decline to import one; there was no object to import it
onto."* Per §4.7 there are at least three such objects. The correct statement is the modest one:
*this page imported no split*, which is true and sufficient. The escalation to *impossible* is not.

**(2) A reduced state was computed and not declared as one.** Q1's probability vector is the
restriction of the state to the vertex algebra (§4.7(a), verified identical). Under custody's
disclosure standard that belongs in the Q5 table with its provenance — **K1's own**, since the
vertex algebra is K1's own indexing — not omitted because it was not called a partial trace.

**(3) A three-way decomposition was derived, used, and not listed.** Construction 1 produces
`{v0} / {v1,v2} / {v3,v4}`, and the build itself remarks that "the root lies in neither" sector
(`:356–358`). That is a three-way grouping of the vertex set of exactly the shape the predecessor
stipulated. It is **derived, not stipulated**, which is the important distinction and is to the
build's credit — but the Q5 table claims "no subsystem was designated" while the page's central
construction designates three. It should be listed as **K1's own, derived**.

---

## 7. D6 — CUSTODY — **CONFIRMED-WITH-CORRECTIONS**

### 7.1 What passes, and it is a high standard

**Pointer rule (custody §1, `:7–20`) — every pointer checked, every pointer resolves.** All
sixteen `file:line` references were opened and matched:

```
S1 :14 root · :19-21 edges · :24-25 faces · :40-41 b1=1 · :48-49 fibres · :51-54 connection
   :59-60 gauge · :66-71 count · :75-78 invariants · :89-92 f=da · :109 the point · :114-115 instance
FOUNDING_DESIGN :51-66 obstruction · :54-55 quotes · :93-98 S2 spec · :97 falsifier · :117-118 contact
CUSTODY :44-53 pairing · :55-68 alpha · :70-74 metric      INHERITED :6-9 set-aside · :11-16 salvage
ALL RESOLVE CORRECTLY.  Digests: 5 cited, 5 verified at bytes.
```

**Alpha (custody §5).** Not engaged. No coupling, no measured number, no target-driven selection
anywhere on the page. I6's claim is correct. My own Construction A likewise introduces no number
from outside; the only constants in it are `W_F`, `W_C` and a probability vector.

**Predecessor material (custody §7).** Cited by digest, never copied, never load-bearing. I2 and I3
are correctly graded, and §4.2 exhibits recurrence numerically rather than inheriting it — which is
the right way to use a premise whose transfer grade has not been issued.

**Self-report.** F1–F6 are honest. F2 anticipates part of COR-B and F6 correctly refuses the
"independently-corroborated" grade. That refusal is now discharged: this audit is
**adversarially-checked**, by the same model, and per custody §4 it may not be graded higher.

### 7.2 Corrections

**(a) Assertion in place of exhibit — custody's standard, §3.4's table.** Three of four printed
numbers are not reproducible and one contradicts the build's own proof (COR-D). This is the
clearest custody breach on the page: the build says "exhibited, all four cases" and what is
exhibited is not recomputable.

**(b) A governing claim beyond what is proved.** `:22` — *"It cannot, and the obstruction is
threefold and proved."* Two of the three are false and the third is scoped narrower than stated. A
headline is a governing clause.

**(c) A modelling choice absent from the CHOICE LEDGER.** The change of state space between Q1 and
Q2 (§3.2). C2 gestures at it and closes it **backwards**, which is worse than omitting it: it
records the rejected alternative as covered by a theorem that in fact only holds on the rejected
alternative.

**(d) C3 was left OPEN and did not need to be.** Per §5.2 the whole Hamiltonian question is an
artefact of choosing addition (i). With zero addition there is no Hamiltonian to choose. An OPEN
ledger entry that a better construction dissolves should be recorded as dissolved, not carried.

**(e) F5 is correctly bounded.** The "first in this lineage" claim is properly limited to the page.
No correction; recorded because I checked it and it holds as bounded.

---

## 8. CORRECTIONS, IN SEVERITY ORDER

**COR-A — `O1` is false. Theorem 1 is a theorem about `L_v0 = C`, promoted to a global obstruction
on the state space `Gamma(L) = C^5` that Q1 itself derived.** The loop transports extend
canonically to `Gamma(L)` as diagonal non-scalar unitaries `M_dF`, `M_c`, and on K1's own published
connection (S1 `:114–115`) with ready state `p = (1/2,0,0,1/4,1/4)` they give
`<M_dF s, M_c s> = 0` exactly. *"Nothing moves the ray"* is false: `|<s, M_dF s>| = 0` on the same
instance. Rewrite `:22–28`, `:25`, `:181–221`, and the falsifier verdict at `:672–675`.

**COR-B — `O3`/Theorem 2's conclusion is false.** Primality of 5 forbids only a tensor
factorisation *of the whole space*. Reduced states exist (the build computed one in Q1), a
canonical three-way direct-sum subsystem split exists on K1's own data, and a genuine tensor
factorisation `C^4 = C^2 (x) C^2` exists on K1's designated non-root subspace. Retain only F2's
weaker structural half. Rewrite `:27`, `:36–38`, `:369–392`, `:539–540`.

**COR-C — Q4's minimum is not minimal.** Zero addition suffices (Construction A), including the
"at a time solving an equation" clause via the carrier's own circuit count. This also dissolves the
OPEN ledger entry C3. Rewrite `:518–530`.

**COR-D — §3.4's exhibited table is non-reproducible and self-contradictory.** Correct values are
`0.000000 / 0.816497 / 0.816497 / 1.000000`; the trivial-connection row **must** be `1.000000`
because the build's own proof establishes a non-trivial intersection there. Theorem 3 itself
stands (0 mismatches in 20,000 draws, reproduced). Replace the table at `:335–340`.

**COR-E — Q1's §1.4 "derived result" is false for the joint system.** `(connection, state)` mod
gauge has `16 - 5 = 11` invariants; the page exhibits 7. State phases are pure gauge only when the
connection is ignored. Rewrite `:108–110`.

**COR-F — `O2` is overstated.** K1 supplies no *continuous* time; it does supply a discrete one
(edge traversals / circuits), which is enough to write the gate's time clause. Rewrite `:26`,
`:283–288`.

**COR-G — CHOICE LEDGER C2 is closed backwards**, and the Q5 disclosure table omits two
decompositions actually used: the reduced state on the vertex algebra, and the derived three-way
grouping `{v0}/{v1,v2}/{v3,v4}`. Both are K1's own; both should be listed. Rewrite `:593`,
`:544–551`.

**COR-H — §4.2's recurrence figure `0.994373` is a window artefact.** Over `t* + 5000` the
supremum is `0.999793`. Label it a lower bound. Rewrite `:465`.

---

## CHOICE LEDGER — MINE

| # | Choice | Alternatives | Why | Status |
|---|---|---|---|---|
| A1 | Extend loop transport to `Gamma(L)` by acting on **every vertex of the loop** and leaving off-loop vertices fixed | extend by zero off the loop; act only at the root (the build's choice) | it is the unique extension that is (a) unitary, (b) uses no data beyond fibres/edges/orientation/connection, (c) reduces to the build's `T_gamma` on `L_v0`. Extension by zero was computed and gives a connection-**independent** verdict (`z_0 = 0`), failing the build's own C4 | **closed** — the alternative was computed and rejected on the build's own criterion |
| A2 | Time = number of circuits `n` of the loop | a real parameter `t` with a Hamiltonian (the build's addition (i)); no time at all | edge count is carrier-supplied combinatorics (S1 `:16–22`); a real `t` is not. `n` is monotone and ordered, which is all the gate's clause requires | **closed** — and it strictly dominates addition (i), which needs `t` *and* an open Hamiltonian choice |
| A3 | Unit weights in `<s,t>`, inherited from the build's C1 | any positive `w_v` | I audit the build's object, and independently: Theorem A's coefficients are unit-modulus regardless, and reweighting only reparametrises the simplex, so the firing criterion is unchanged | **closed** — verified: the convex-hull criterion is weight-independent |
| A4 | Principal angles (SVD of the Gram of orthonormal bases) as the meaning of `max\|<P_F,P_C>\|` | any spanning-set Gram | it is the only basis-independent reading; the build's numbers are not reproducible under any orthonormal convention | **closed** |
| A5 | Ready states reported as probability vectors `p` | full complex sections | `p` is the build's own Q1 complete invariant for rays, so Theorem A is stated in the carrier's own coordinates | **closed** — checked against direct matrix action, agreement to `1e-10` |

---

## IMPORT AUDIT — MINE

| # | Notion | Source | Defined here? | Survives without it? |
|---|---|---|---|---|
| J1 | K1 in full | `S1_CARRIER_K1_V001.md` sha256 `3eb70375…d92ac` | restated where used | n/a — the subject |
| J2 | The target artifact | sha256 `248ce856…734d9`, verified at bytes | quoted by line | n/a — the subject |
| J3 | Conditional expectation, block pinching, algebraic subsystem | standard operator algebra | not defined here; **exhibited numerically** (idempotence, positivity, unitality, trace preservation all checked) | **YES** — COR-B also stands on the direct-sum split and the `C^2 (x) C^2` factorisation alone, both elementary linear algebra |
| J4 | Numerical range / convex hull of points on a circle | standard | Theorem A is proved from scratch above | **YES** |
| J5 | The predecessor's gate wording | the brief for this audit, restating predecessor material | quoted, **not** premise | **YES** — it is the target of the test. Carries no transfer grade (`INHERITED :11–16`); nothing here is founded on it |
| J6 | Alpha | not used | — | **YES** — no coupling, no measured number, no target-driven selection in this audit or in Construction A |

---

## FLAG BLOCK — DEFECTS IN MY OWN AUDIT

**G1 — Construction A is transport, not dynamics, and I must not let COR-A obscure that.** `M_dF`
and `M_c` are parallel transport operators. They defeat `O1` and they solve the gate's *orthogonal*
clause and, via A2, its *at a time* clause. They do **not** supply a Hamiltonian flow, and the
build's refusal to call transport a dynamics (§2.2) remains correct and remains to its credit. What
I refute is the claim that rank forbids the condition — not the claim that K1 lacks dynamics.

**G2 — A1 is a choice and Theorem A depends on it.** Extension by zero, also defensible, gives a
connection-independent verdict and would **not** write a formation condition. I computed and
declared the alternative, but a different extension convention is the one place Construction A
could be attacked, and an S3 lane should attack it there first.

**G3 — Theorem A's firing criterion is symmetric in `f, c`; only the *fixed-ready-state* condition
is asymmetric.** The existential statement "some ready state fires" has a `f <-> c` symmetric
region (the arc-gap multiset for `{0, c, c+f}` is symmetric). So §4.6 closes F1 for a *designated*
ready state, not for the existence question. That is a real closure but a narrower one than a
casual reading of §4.6 would suggest, and I state it here rather than let it be found.

**G4 — the recurrence exhibit in §4.9 is numerical over `n <= 400`, not a proof.** The
almost-periodicity of `conj(W_F^n) W_C^n p_0 + conj(W_F^n) q + W_C^n r` on the closure of the
subgroup generated by `(W_F, W_C)` in the 2-torus is standard, but I did not prove it here; I
exhibited it. The qualitative verdict (fires, then un-fires) is robust; the value `n = 42` is not.

**G5 — COR-B(c)'s `C^4 = C^2 (x) C^2` uses an ordering of the two non-root vertices of each
triangle.** That ordering is K1's own (edge orientation from the root, S1 `:19–21`), but it is a
labelling, and a reader who rejects it still has COR-B(a) and COR-B(b), which use no ordering at
all. The refutation of Theorem 2 does not depend on this item.

**G6 — same-model lineage, custody §4 (`:44–53`).** Brief-writer, builder and auditor are the same
model. A failure mode shared across those roles passes through all three invisibly. Everything here
is graded **adversarially-checked**, never **independently-corroborated**. In particular COR-A
found a gap the build half-saw and walked past (`:292–294`) — which is exactly the shape of failure
custody §4 predicts, and its recurrence in a second pass by the same model cannot be excluded.

**G7 — I did not re-audit S1.** I re-derived every S1 count I used (`chi`, `b1`, `rank d1`, the
`6-4=2` invariant count, `d^2 = 0`) and all reproduce, but S1 was not the target and no verdict
here should be read as clearing it.

---

## 9. WHAT THIS AUDIT DELIVERS

**Delivers:** all five cited digests verified at bytes; S1's topology and gauge counts
independently re-derived; Theorem 1 reproduced at 200,000 samples and correctly scoped;
Theorem 3 reproduced at 20,000 samples with a corrected exhibit; every number in Q4's Laplacian
analysis reproduced to six digits; the joint invariant count the build did not do; **Construction
A**, which writes the formation condition the build declared unwritable, with zero addition, on
K1's own published connection, exactly, and which closes the build's own F1; and the identification
of **durability, not rank and not split**, as the only obstruction that was ever load-bearing.

**Does not deliver:** any crossing (S3); any dynamics on K1; a proof that Construction A's
recurrence is unavoidable; independent lineage.

**The build's falsifier was "if it cannot be stated without a stipulated split, say so plainly"
(`FOUNDING_DESIGN_V001.md:97`). It can be stated without one — and, contrary to the build, it can
also be stated without rank, without a split's impossibility, and without adding anything at all.
What it cannot yet do is last.**

## 10. CUSTODY

Audited under `CUSTODY_V001.md` sha256
`6629ae3b5aedbafda4e56b2a743b84a888949c59f988296056a2a7abeb20aa49`. Default verdict REFUTED;
every dimension graded against recomputation, not testimony. Every claim exhibited on this page.
Predecessor material cited by digest only, never as foundation. Grade: **adversarially-checked**,
never **independently-corroborated** (custody §4, G6). No git action taken. Sealed on creation.
