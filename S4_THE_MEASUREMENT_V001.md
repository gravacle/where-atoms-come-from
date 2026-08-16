# S4 — THE MEASUREMENT: DOES FORMATION SEE THE STRUCTURE? — BUILD — V001 — 2026-08-16

**Stage S4 of `FOUNDING_DESIGN_V001.md` §7** — *"RUN IT, AND VARY THE ELECTROMAGNETIC
STRUCTURE. Change holonomies, winding, charge. Does formation change, and how?
Falsifier: formation is wholly insensitive to the connection — a real and publishable
negative."*

This is a **computational** stage. Every number below was produced by code written for this
page, run in this session, and — wherever a closed form exists — **re-derived independently
in exact arithmetic and checked against the direct simulation.** The code is listed in the
COMPUTATION LEDGER.

---

## 0. WHAT S4 FOUND, IN ORDER OF HOW MUCH IT COSTS THE PROGRAM

1. **The rate `lambda` is not one function. It is two, and the schedule chooses which.**
   Under the uniform clock `lambda_A = log|Z_1|` is a genuine two-parameter function on the
   torus of connections. Under the carrier's own canonical clock `lambda_B` is
   **constant on a set of full measure**, and varies only on a dense measure-zero
   exceptional set. Both are computed and mapped below.

2. **`lambda_B` is a function of the RELATION LATTICE, not of the holonomies.**
   `lambda_B(W_F, W_C)` depends on the pair only through
   `L = { (m,n) in Z^2 : W_F^{-m} W_C^{n} = 1 }`. Verified to `3.3e-16` across 49 distinct
   lattice classes, 144 connections. This is a theorem (Pontryagin duality + Weyl), and it
   is proved and checked below.

3. **F-A does not fire on full support, and FIRES on three of the four support classes.**
   With weight at the root, `lambda` sees `W_F` and `W_C` separately — proved exactly, by
   the Fourier support of `|Z_1|^2` being rank 2. Remove the root's weight and `lambda`
   collapses to a function of the **product** `W_F·W_C` alone. **The falsifier's condition
   is exactly realized whenever `p_0 = 0`.**

4. **F-B does not fire on its literal condition — and its target claim is CONFIRMED by two
   controls built to test it.** `lambda` is not constant across the family. But what it
   responds to is **not topology**: a control that moves `chi` from 0 to 1 leaves `lambda`
   **bit-identical**, and a control that holds every topological invariant fixed moves
   `lambda` by `0.0634`. **The carrier's topology is inert. Its loop incidence is not.**

5. **A DEFECT OF RECORD.** S3's headline test connection `f = 2.0, c = 1.1` is **exactly
   resonant** — `-11f + 20c = 0` identically — so the orbit is **not** dense in `T^2`, as
   S3 §6(f) and the S3 audit both state. The reported number `-0.767026` is correct and is
   converging to the **subtorus** value `-0.767014993`, not to the full-torus mean
   `-0.767507880`. Every row of S3 §5.7 is a subtorus value. Detail at §7.

6. **Q8: no formation quantity relates to `2*pi*chi`, and the reason is structural, not
   numerical.** `chi -> lambda` is not a well-defined map: `chi = 0` carries six distinct
   `lambda` values in this family, and the fill control moves `chi` with `lambda` frozen.

---

## 1. SEALS VERIFIED BEFORE ANY WORK

All from each artifact's own directory, `shasum -a 256 -c`, both sidecar forms:

```
S1_CARRIER_K1_V001.md                        OK / OK
S2_FORMATION_CONDITION_ON_K1_V001.md         OK / OK
S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md   OK / OK
S3_THE_CROSSING_V001.md                      OK / OK
S3_THE_CROSSING_AUDIT_V001.md                OK / OK
CUSTODY_V001.md                              OK / OK
FOUNDING_DESIGN_V001.md                      OK / OK
REGISTER_V001.md                             OK / OK
PUBLICATION_V001.md                          OK / OK
INHERITED_FROM_THE_PREDECESSOR.md            OK / OK
```

**20 of 20 verify.** `S4_THE_MEASUREMENT_V001.md` did not exist before this write.

### 1.1 The corpus reproduced from scratch before it was used

`K1` was rebuilt from `S1_CARRIER_K1_V001.md` §1 as an incidence structure, and the overlap
was computed by **direct matrix action on `C^5`** — five diagonal unitaries, random vertex
phases on the ready section — not from any inherited closed form:

```
closed form vs matrix action, max deviation over k <= 200 :  1.487e-14
      (S3 build 1.538e-14 ; S3 audit 1.528e-14)
8 random gauge transformations, spread                    :  1.962e-15
K1 topology recomputed from incidence: V=5 E=6 F=1 chi=0 b0=1 b1=1 b2=0
      gauge acts through 4 parameters ; 2 invariants ; 1 curvature + 1 flat
d1 . d2 = 0 exactly
```

S3's §4.3 table, all thirteen rows, reproduced to the last printed digit:

```
     N     |Z_N|     (1/N)log|Omega_N|   S3 build     verdict
     1   0.411271       -0.888504       -0.888504     MATCH
     2   0.470386       -0.821353       -0.821353     MATCH
     5   0.364118       -0.806574       -0.806574     MATCH
    10   0.776953       -0.733234       -0.733234     MATCH
    20   0.247337       -0.759800       -0.759800     MATCH
    42   0.024654       -0.800483       -0.800483     MATCH
    50   0.681600       -0.819895       -0.819895     MATCH
   100   0.350139       -0.768925       -0.768925     MATCH
   200   0.610201       -0.769663       -0.769663     MATCH
   400   0.600798       -0.768673       -0.768673     MATCH
  1000   0.580708       -0.768205       -0.768205     MATCH
  2000   0.573418       -0.768470       -0.768470     MATCH
  4000   0.360282       -0.766802       -0.766802     MATCH        13 of 13

lambda at N = 200000                    = -0.767026     MATCH
min |Z_k|, k<=400   = 0.024654 at k=42                  MATCH
sup |Z_k|, k<=4000  = 0.999941 at k=377                 MATCH
#{|Z_k|>0.99, k<=4000} = 37                             MATCH
sum(1-|Z_n|), n<=100000 = 46918.264                     MATCH
K1's own connection: Z1=0, Z2=-1, Z3=0, Z4=+1, Z5=0, Z6=-1   MATCH
```

**Nothing on this page is inherited numerically. The corpus was re-run, then departed from.**

---

## 2. THE FORMATION FUNCTIONAL, DERIVED FOR AN ARBITRARY CARRIER

W-01's overlap is written for `K1`. S4 needs it on ten carriers, so it is **derived**, not
transported.

Let `K` be a finite CW complex of dimension `<= 2` with a rank-one Hermitian fibre at each
vertex and a `U(1)` connection on edges. Let `gamma_F` be the boundary of a filled 2-cell
and `gamma_C` a 1-cycle. W-01's operators on `Gamma(L) = C^V` are

```
(M_gamma s)(v) = W(gamma) · s(v)   if v lies on gamma ;   s(v)  otherwise.
```

Then, with `s` any ready section, `p_v = |s_v|^2`, and `k` circuits of each loop:

```
Z_k = < M_{gamma_F}^k s , M_{gamma_C}^k s >
    = sum_v  conj(W_F)^{k a_v} · W_C^{k b_v} · p_v          a_v = [v on gamma_F]
    = sum_v  u^{k a_v} v^{k b_v} p_v                        b_v = [v on gamma_C]
```

**Every vertex phase of `s` cancels.** So the whole carrier enters formation through one
object: the **vertex class weight vector**, indexed by `(a,b) in {0,1}^2`, with characters

```
(0,0) -> 1        (1,0) -> u = conj(W_F)      (0,1) -> v = W_C      (1,1) -> uv
```

On `K1` the classes are `(1,1):{v0}`, `(1,0):{v1,v2}`, `(0,1):{v3,v4}`, and the formula
returns exactly W-01's `conj(W_F)W_C p_0 + conj(W_F)(p_1+p_2) + W_C(p_3+p_4)`, and the
corrected criterion's `chi_0 = uv, chi_F = u, chi_C = v`. **The three characters of the
corrected criterion are the three non-empty vertex classes of `K1`. That is what they are.**

### 2.1 The exact modulus — the identity the whole of §3 rests on

Expanding `|Z_k|^2` for the three-class case gives, with no approximation,

```
|Z_k|^2  =  p0^2 + q^2 + r^2  +  2 p0 q cos(k c)  +  2 p0 r cos(k f)  +  2 q r cos(k(f+c))
```

`(p0, q, r)` the class weights, `W_F = e^{if}`, `W_C = e^{ic}`.
Checked against the general routine at **2000 random `(f,c,k)`: max deviation `2.287e-14`.**

**Read the three frequency vectors.** In `(f,c)` they are `(0,1)`, `(1,0)`, `(1,1)`, with
coefficients `2 p0 q`, `2 p0 r`, `2 q r`. **They span `Z^2`.** Section 3.2 turns this into
the proof of the F-A verdict.

### 2.2 `lambda`, and what it actually is

Under a schedule `k_n`, `Omega_N = prod_n Z_{k_n}` and `lambda = lim (1/N) log|Omega_N|`.

- **Schedule A (uniform, `k_n = 1`):** `lambda_A = log|Z_1|`.
- **Schedule B (the canonical clock, `k_n = n` — S3's primary):** the sequence
  `(u^n, v^n)` equidistributes (Weyl) on the closure `H` of the group it generates, and

```
lambda_B  =  int_H log|Z| d(Haar_H)
```

`H` is a closed subgroup of `T^2`. By Pontryagin duality `H = L^perp`, where

```
L  =  { (m,n) in Z^2 : u^m v^n = 1 }  =  { (m,n) : W_F^{-m} W_C^{n} = 1 }
```

is the **relation lattice**. Three cases, all computable in closed form:

| `rank L` | `H` | `lambda_B` |
|---|---|---|
| 0 | `T^2` | 2-variable Mahler measure `m(p00 + p10 x + p01 y + p11 xy)` — a **constant** |
| 1, generator `d·(m,n)`, `(m,n)` primitive | `d` cosets of a circle | mean of `d` **1-variable** Mahler measures |
| 2 | finite | exact finite average over the orbit |

**`lambda_B` is therefore a function of `L` alone** — of the *arithmetic* of the pair, not of
its values. Tested three ways in §3.3.

---

## 3. S4-A — THE CONNECTION AXIS

Held fixed unless stated: carrier `K1`; class weights `(p0,q,r) = (0.4,0.3,0.3)`, which is
S3's `p = (0.4,0.15,0.15,0.15,0.15)`. **Both the schedule and the support are varied and
reported as axes, per the carried correction.**

### 3.1 Q1 — THE MAP `lambda(W_F, W_C)`

#### Schedule A — an honest two-parameter map

`lambda_A = (1/2) log[ 0.34 + 0.24 cos c + 0.24 cos f + 0.18 cos(f+c) ]`.

**Critical points, solved analytically** (stationarity gives `sin f = sin c`, hence the
branches `c = f` and `c = pi - f`):

```
   (f, c)                      |Z_1|        lambda_A      character
   (0, 0)                    1.0000000     +0.000000      MAXIMUM (unique)
   (arccos(-2/3), same)      0.0000000     -infinity      ZERO — exact firing
   (2pi-arccos(-2/3), same)  0.0000000     -infinity      ZERO — exact firing
   (pi, pi)                  0.2000000     -1.609438      local minimum
   (0, pi)                   0.4000000     -0.916291      saddle
   (pi, 0)                   0.4000000     -0.916291      saddle
```

`arccos(-2/3) = 2.300523983`. On the diagonal `c = f` with `q = r` the overlap collapses to
`Z_1 = 0.4 + 0.6 cos f`, so the zeros are **exact and there are exactly two on `T^2`** —
two real equations in two unknowns. A `1200 x 1200` grid finds `min |Z_1| = 3.097e-04` and
`262` cells under `0.01`, consistent with two isolated conical zeros.

- **Level sets:** closed curves; symmetric under `(f,c) -> (c,f)` (because `q = r` here) and
  under `(f,c) -> (-f,-c)` (always). They are **not** straight lines, which §3.2 makes into
  the F-A proof.
- **Singular locus (`lambda = -infinity`):** exactly the two points above. This is W-01's
  convex-hull criterion; it is non-empty iff `max(p0,q,r) <= 1/2`, i.e. iff the weights obey
  the triangle inequality.
- **Locus where formation fails:** `|Z_1| = 1` forces the three unit numbers to coincide,
  hence `W_F = W_C = 1`. **One point of `T^2`.**
- **Torus mean:** `-0.767508397` on the grid, against the exact Cassaigne–Maillot value
  `-0.767507880`.

#### Schedule B — constant almost everywhere

```
predicted generic rate  m(p0 + q x + r y) = -0.767507880   (exact, Cassaigne-Maillot)

120 uniformly random connections, N = 400000 each:
   mean -0.767508   median -0.767508   std 0.000021   min -0.767635   max -0.767411
   #{ |lambda_B - m| < 0.01 } = 120 of 120        outliers: 0
```

**`lambda_B` is constant on a full-measure set.** It varies only on the exceptional set

```
   { (f,c) : there exists a primitive (m,n) with -m f + n c = 0 mod 2pi }
```

a countable union of closed geodesics on `T^2` — dense, Haar-measure zero. On the primitive
locus `(m,n)` with `d = 1`, `lambda_B` is the one-variable Mahler measure
`m(r z^{m+n} + p0 z^{m} + q)`, exactly:

```
  (m, n)     lambda_resonant     direct on the locus (N=2e6)     dev
  ( 1,  0)    -0.356674944          -0.356675071               1.3e-07
  ( 0,  1)    -0.356674944          -0.356675071               1.3e-07
  ( 1,  1)    -1.203972804          -1.203995817               2.3e-05
  ( 1, -1)    -0.510825624          -0.510825783               1.6e-07
  ( 2,  1)    -0.681980359          -0.681980552               1.9e-07
  ( 2, -1)    -0.916290732          -0.916290860               1.3e-07
  ( 3,  1)    -0.767783712          -0.767784078               3.7e-07
  ( 3,  2)    -0.732940865          -0.732941263               4.0e-07
  ( 4,  1)    -0.784966659          -0.784966909               2.5e-07
  ( 5,  1)    -0.749392712          -0.749393530               8.2e-07
  ( 5,  3)    -0.765224351          -0.765224482               1.3e-07
  ( 7,  3)    -0.759305247          -0.759304981               2.7e-07
  ( 7, 11)    -0.764712281          -0.764711656               6.3e-07
  (11, 20)    -0.767014993          -0.767015210               2.2e-07
  (13,  8)    -0.768271734          -0.768271717               1.7e-08
  (29, 17)    -0.767138179          -0.767135497               2.7e-06
```

Nineteen distinct exceptional values found; range `[-1.203973, -0.356675]`; **13 above the
generic value and 6 below.** The exceptional values **accumulate** on the generic one as the
resonance grows:

```
   (m,n)        |m|+|n|      lambda            deviation from generic
   (1,1)             2     -1.203972804          -4.36e-01
   (5,3)             8     -0.765224351          +2.28e-03
   (11,20)          31     -0.767014993          +4.93e-04
   (41,53)          94     -0.767416955          +9.09e-05
   (97,61)         158     -0.767515568          -7.69e-06
   (610,377)       987     -0.767507355          +5.25e-07
```

**The map, stated whole.** `lambda_B` is `-0.767507880` off a dense null set; on that set it
takes a countable dense-in-`[-infinity, -0.33]` family of other values; it is `0` at the
single point `W_F = W_C = 1`; and it is `-infinity` on the countable dense set where some
`Z_k` vanishes exactly (for each `k`, the preimage of the two zeros under `(f,c)->(kf,kc)`,
i.e. `2k^2` points — **verified: 110 exact zeros found among the `k<=5` preimages, predicted
`2*sum_{k<=5} k^2 = 110`**).

**The design's named S4 external contact** (`FOUNDING_DESIGN_V001.md:119`, *"dependence on
winding must be integer-periodic"*) — **discharged**: `max| |Z_1(f,c)| - |Z_1(f+2pi k, c+2pi l)| |`
over 500 random samples with `k,l in [-5,5]` is `2.33e-15`; and for `lambda_B` the
periodicity is exact by construction, since `L` is defined from `u,v` which depend on
`(f,c)` only mod `2pi`. **The S2 contact also survives:** `lambda_A(0,0) = lambda_B(0,0) = 0`
exactly — no formation at the trivial connection.

### 3.2 Q2 — THE SEPARATION TEST, AND THE F-A VERDICT

**Schedule A, full support: `lambda_A` depends on `W_F` and `W_C` separately, and this is
provable rather than sampled.**

`|Z_1|^2` has, exactly, the Fourier expansion of §2.1. Its support on the dual lattice of
`(f,c)` is computed and confirmed:

```
   (m,n)        |coefficient|      predicted (= 2 p0 q etc.)
   ( 0, 0)      0.340000000        0.340000000
   (+-1, 0)     0.120000000        0.120000000   = 2 p0 r / 2
   ( 0,+-1)     0.120000000        0.120000000   = 2 p0 q / 2
   (+-1,+-1)    0.090000000        0.090000000   = 2 q r / 2
   total mass off {(0,0), +-(0,1), +-(1,0), +-(1,1)} = 7.094e-14
```

> **PROPOSITION.** If `lambda_A` were a function of a single combination `m f + n c`, then
> `|Z_1|^2 = exp(2 lambda_A)` would be too, and its Fourier support would lie inside the
> **rank-1** line `Z·(m,n)`. Its support contains `(0,1)`, `(1,0)` and `(1,1)`, which **span
> `Z^2`**. Contradiction. **QED**

Direct witnesses, each pair sharing one candidate combination and disagreeing in `lambda_A`:

```
   same W_F  (f=1.0):        lambda_A(1.0,0.3) = -0.145776   lambda_A(1.0,3.0) = -1.083946
   same W_C  (c=1.0):        lambda_A(0.3,1.0) = -0.145776   lambda_A(3.0,1.0) = -1.083946
   same product (f+c=3.0):   lambda_A(0.4,2.6) = -0.865231   lambda_A(1.5,1.5) = -0.815445
   same ratio (c-f=0):       lambda_A(0.5,0.5) = -0.076288   lambda_A(2.5,2.5) = -2.517188
```

**Schedule B, full support:** `lambda_B` depends on the pair **only through `L`**, and `L`
is not a function of any single combination — the product `W_F·W_C` fixed at `e^{3.1i}`
carries four different lattices and four different rates:

```
   f=2.0 c=1.1  L=<(11,20)>  lambda = -0.767014993
   f=1.0 c=2.1  L=<(21,10)>  lambda = -0.768980527
   f=0.5 c=2.6  L=<(26, 5)>  lambda = -0.767137326
   f=1.5 c=1.6  L=<(16,15)>  lambda = -0.766734063
```

and the ratio `W_C/W_F` fixed at `1` carries `-1.203973` (`f=c=2`), `-1.535057` (`f=c=2pi/3`)
and `-1.198293` (`f=c=2pi/5`).

### **ARMED FALSIFIER F-A — VERDICT: DOES NOT FIRE ON FULL SUPPORT; FIRES ON THREE OF THE FOUR SUPPORT CLASSES.**

The falsifier is **support-conditional**, and this is the sharpest thing on the page.
Reading §2.1 with a weight set to zero kills the corresponding cosine terms:

| support `S` | `G` | `lambda_A` is a function of | F-A |
|---|---|---|---|
| `{0,F,C}`, `p0>0` | rank-2 `<u,v>` | `cos c`, `cos f` **and** `cos(f+c)` | **does not fire** |
| `{0,C}` (`q=0`) | `<u>` | `cos f` only — **`W_F` alone** | **FIRES** |
| `{0,F}` (`r=0`) | `<v>` | `cos c` only — **`W_C` alone** | **FIRES** |
| `{F,C}` (`p0=0`) | `<u/v>` | `cos(f+c)` only — **the product `W_F·W_C` alone** | **FIRES** |
| `|S|=1` | `{1}` | nothing; `|Z| = 1` | vacuous — never forms |

**The load-bearing statement.** `p0` is the weight at **`v0`, the join vertex — the pinch**.
`K1`'s separation of curvature from flat holonomy does real work **if and only if the ready
state puts weight on the pinch.** With `p0 = 0` the carrier's distinction between a bounded
and an unbounded cycle is exactly the decoration F-A describes: only `W_F·W_C` survives.

This is not a hedge. It is a statement about which states can see the carrier, and it
sharpens W-01's *"the root can never fire"* into its complement: **the root alone can never
fire, and without the root the two invariants collapse into one.**

Closed forms for the collapsed classes, exact by Jensen (`m(a + b z) = log max(|a|,|b|)`):

```
   S = {0,C}  lambda = log max(p0, r)        (0.5,0,0.5)   -> log 0.5 = -0.693147181
   S = {0,F}  lambda = log max(p0, q)        (0.5,0.5,0)   -> log 0.5 = -0.693147181
   S = {F,C}  lambda = log max(q, r)         (0,0.3,0.7)   -> log 0.7 = -0.356674944
```

### 3.3 Q3 — RATIONAL vs IRRATIONAL, FINITE `G` vs DENSE `G`

**The character of formation does not distinguish them. The rate does, discontinuously, and
the right invariant is not "rational" but "resonant".**

`G != {1}` — formation itself — is insensitive: `G` is non-trivial for every connection off
`W_F = W_C = 1`, finite or dense alike. **Formation is binary and blind to the arithmetic.**
The rate is not.

**The three tiers, established by computation:**

- **TIER 3, `rank L = 0`** (no relation): `H = T^2`, `lambda_B = -0.767507880`, constant.
- **TIER 2, `rank L = 1`**: `lambda_B` depends on `(m,n)` **and on nothing else** — not on
  where on the circle the connection sits. Six points on the locus `-2f+3c = 0`:

```
   (f,c) = (0.9,0.6) (3,2) (6,4) (8.1,5.4) (13.2,8.8) (17.7,11.8)
   lambda_B(N=2e6) = -0.732941279 -0.732941263 -0.732941831 -0.732939926 -0.732941131 -0.732940590
   exact                                  -0.732940865   for all six
```

- **TIER 1, `rank L = 2`** (both holonomies roots of unity): `lambda_B` is an exact finite
  average, and it depends on `L` — **not** on the order of the orbit. Grouping all 143
  non-trivial pairs with denominator 12 by the Hermite basis of their lattice:

```
   49 lattice classes ; MAX spread of lambda inside any class = 3.331e-16
   29 distinct lambda values across the 49 classes
```

**That `3.3e-16` is the whole answer to Q2 under schedule B and to Q3 at once:
`lambda_B` is a function of `L` and of nothing else.**

**Is there a discontinuity at the rationals? Yes — and everywhere else on the exceptional
set too.** Approaching `f = c = 2pi/3` along a generic direction (both coordinates moved,
in ratio `1 : sqrt2`, so that `L` collapses to `0`):

```
   at f=c=2pi/3 exactly                               lambda = -1.535056729
   f=2pi/3+1e-1, c=2pi/3+sqrt2*1e-1  (N=8e6)          lambda = -0.767507087
   f=2pi/3+1e-2, c=2pi/3+sqrt2*1e-2                   lambda = -0.767503738
   f=2pi/3+1e-3, c=2pi/3+sqrt2*1e-3                   lambda = -0.767546770
   f=2pi/3+1e-4, c=2pi/3+sqrt2*1e-4                   lambda = -0.766759838
   generic value                                              -0.767507880
   JUMP AT THE RATIONAL POINT                                 -0.767549
```

**A methodological correction I make against my own first attempt, recorded rather than
silently fixed:** my first discontinuity test held `f = 2pi/3` and moved only `c`. That is
**confounded** — `f` alone rational keeps the relation `(3,0)` alive, so `L` never becomes
trivial and the limit is the *partial*-resonance value, not the generic one. Under that
(wrong) protocol I obtained `-0.798965`; the correct closed form for `L = <(3,0)> = <3·(1,0)>`,
`d = 3`, is `-0.798965257`, and the direct run at `N = 8e6` gives `-0.798965194`. The number
was right; the interpretation would have been wrong.

**Semicontinuity: neither.** Over all 638 finite exceptional values with denominator `<= 12`:
`min -1.535057`, `max -0.331417`, **380 above the generic value and 258 below.**
So `lambda_B` is neither upper- nor lower-semicontinuous anywhere on the exceptional set.

**Is "rational" the right word? No, and this is a correction to how S3 §6(f) frames it.**
The invariant is the relation lattice. Roots-of-unity pairs are merely the sub-case
`rank L = 2`. A pair with **both holonomies irrational multiples of `2pi`** can still be
resonant — `(f,c) = (3t, 2t)` for any `t` — and behaves exactly like a "rational" point.

### 3.4 Q4 — THE SUPPORT AXIS

```
  support     weights        characters   rk G   lambda_B (exact)        lambda_A sees
  |S|=3   (0.4,0.3,0.3)   uv, u, v         2     -0.767507880       cos c + cos f + cos(f+c)
  |S|=3   (1/3,1/3,1/3)   uv, u, v         2     -0.775546341       "
  |S|=3   (0.8,0.1,0.1)   uv, u, v         2     -0.223143551       "
  |S|=3   (0.2,0.5,0.3)   uv, u, v         2     -0.693147181       "
  S={0,C} (0.5,0,0.5)     uv, v            1     -0.693147181       cos f
  S={0,F} (0.5,0.5,0)     uv, u            1     -0.693147181       cos c
  S={F,C} (0,0.5,0.5)     u, v             1     -0.693147181       cos(f+c)
  S={F,C} (0,0.3,0.7)     u, v             1     -0.356674944       cos(f+c)
  |S|=1   any vertex      one              0      0  (never forms)  nothing
```

**The `|S|=3` region has its own phase structure**, inherited from the Cassaigne–Maillot
formula: `lambda` is the smooth dilogarithmic expression when `(p0,q,r)` obey the triangle
inequality, and collapses to `log max(p0,q,r)` — **independent of the other two weights** —
when they do not. The transition lines are `p0 = 1/2`, `q = 1/2`, `r = 1/2`.
S3's own rows land on both sides and both are now exact:

```
  p = (0.8,0.1,0.1)  : 0.8 > 0.1+0.1  -> log-max regime,      exact  log(0.8) = -0.223143551
  p = (0.2,0.5,0.3)  : 0.5 = 0.2+0.3  -> degenerate triangle, exact  log(0.5) = -0.693147181
  p = (1/3,1/3,1/3)  : equilateral    -> log(1/3) + m(1+x+y)          = -0.775546341
  p = (0.1,0.45,0.45): triangle       -> Cassaigne-Maillot            = -0.727723488
```

Sweeping the simplex on a `241`-step grid (`29400` interior points): the **minimum**
`lambda = -0.775535` sits at `(0.3320, 0.3320, 0.3361)`, against the exact equilateral value
`log(1/3) + m(1+x+y) = -0.775546341`; `7260` sampled points are in the triangle regime and
`22140` in the log-max regime. The supremum is `0`, approached only at the vertices.

**Is the `|S|=1` non-formation isolated, or the boundary of a region?**
**It is the boundary, and the rate vanishes continuously into it.** In the binary criterion
it is isolated — `G = {1}` exactly at the three vertices of the simplex and nowhere else. In
the rate it is not isolated at all:

```
   p = (1-2e, e, e)     lambda        log(1-2e)
   e = 0.03           -0.061875404   -0.061875404
   e = 0.01           -0.020202707   -0.020202707
   e = 0.001          -0.002002003   -0.002002003
   e = 0.0001         -0.000200020   -0.000200020
```

`lambda -> 0` linearly in the deficit. **The `|S|=1` locus is where the rate degenerates
continuously to zero, not a puncture.** The binary criterion's discontinuity there is an
artefact of thresholding a continuous quantity.

---

## 4. S4-B — THE CARRIER AXIS

### 4.1 Q5 — THE FAMILY, CONSTRUCTED

The continuum statement being tested — *ring torus -> horn torus -> spindle torus ->
double-covered sphere as the distance from the axis of revolution decreases* — is used
**only to motivate which discrete complexes to build.** No continuum limit is taken and no
step below uses one. **Where the brief's list and the family disagree, I departed and say so
in the CHOICE LEDGER.**

Every number is computed from the incidence matrices `d1`, `d2` — Betti numbers as ranks
over `Q`, `d1·d2 = 0` verified on all ten:

```
carrier                                    V   E   F  chi  b0  b1  b2  gauge inv curv flat
B0a ring torus 3x3 grid, loops disjoint    9  18   9    0   1   2   1     8  10    8    2
B0b ring torus 3x3 grid, loops meet        9  18   9    0   1   2   1     8  10    8    2
B3  horn torus (octahedron, poles ident.)  5  12   8    1   1   1   1     4   8    7    1
B1  K1 (the pinch, as handed)              5   6   1    0   1   1   0     4   2    1    1
B4  spindle (two spheres glued at 2 pts)   6   8   4    2   1   1   2     5   3    2    1
B5  double-covered sphere                  4   4   2    2   1   0   1     3   1    1    0
B2  K1, both triangles filled              5   6   2    1   1   0   0     4   2    2    0
B1p K1-bridged (homotopy equiv. to K1)     6   7   1    0   1   1   0     5   2    1    1
B1q K1-bridged + spectator vertex          7   8   1    0   1   1   0     6   2    1    1
B1s K1, every edge subdivided              11 12   1    0   1   1   0    10   2    1    1

  d1 . d2 = 0 on every carrier: max |entry| = 0.0e+00
```

**The gauge count and the invariant split, derived once and asserted in code for every row.**
`U(1)` gauge acts through `V - b0` effective parameters (the global phase per component acts
trivially), so the number of independent gauge invariants is `E - (V - b0)` = the cycle rank
of the 1-skeleton. Of these, `rank(d2)` are **curvatures** (independent face-boundary
holonomies) and `b1` are **flat**. `curvature + flat = invariants` is checked by assertion on
every carrier and holds on all ten. On `K1`: `2 = 1 + 1`, reproducing S1 §4 from incidence
alone.

**The designated loops are verified, not asserted.** Each is supplied as a signed edge chain
and tested against the boundary maps:

```
carrier                                gF cycle  gF bounds  gC cycle  gC bounds  independent
B0a ring torus, loops disjoint           0        True       0         False       True
B0b ring torus, loops meet               0        True       0         False       True
B3  horn torus                           0        True       0         False       True
B1  K1                                   0        True       0         False       True
B4  spindle                              0        True       0         False       True
B5  double sphere                        0        True      n/a         n/a         n/a
B2  K1 both filled                       0        True       0         TRUE        True
B1p K1-bridged                           0        True       0         False       True
B1q K1-bridged + spectator               0        True       0         False       True
B1s K1 subdivided                        0        True       0         False       True
```

`gF bounds = True` is what makes `W_F` a curvature; `gC bounds = False` is what makes `W_C`
flat; `independent = True` is what makes them independently assignable by a real connection.
**B2 is the deliberate exception** — both loops bound, so B2 has **no flat holonomy at all**
and its second holonomy is a second curvature. That is precisely what makes it the control.

**B5, the degenerate end, has `b1 = 0`: no free cycle exists, so `gamma_C` cannot be
designated and the formation datum does not exist on it.** This is the one place in S4 where
topology, and nothing else, decides an outcome.

### 4.2 Q6 — `lambda` AND `G` ALONG THE FAMILY

**The stated sense in which the connection is held fixed.** For each carrier: `W(gamma_F) =
e^{if}` and `W(gamma_C) = e^{ic}` with the *same* `(f,c)`, generic; **every other independent
invariant set to 1** (the minimal excitation carrying the designated pair — always
realizable, since `gamma_F` and `gamma_C` are linearly independent in the cycle space, which
is verified above). The ready state is held in **two** senses, because the choice is not
forced and it changes the answer:

- **SENSE U** — uniform on vertices, `p_v = 1/V`. The carrier's combinatorics feed in.
- **SENSE C** — fixed at the *class* level: `(0.4,0.3,0.3)` for 3 classes, `(0.5,0.5)` for 2,
  `(0.25,...)` for 4. The carrier feeds in only through *which* classes exist.

```
carrier                          class counts             chars      rkG  lambda(U)      lambda(C)
B0a ring torus, disjoint     {00:2, 01:3, 10:4}          1 v u        2   -0.747659833  -0.767507880
B0b ring torus, meeting      {00:4, 01:1, 10:2, 11:2}    1 v u uv     2   -0.810930216  -1.386294361
B3  horn torus               {01:2, 10:2, 11:1}          v u uv       2   -0.756573586  -0.767507880
B1  K1                       {01:2, 10:2, 11:1}          v u uv       2   -0.756573586  -0.767507880
B4  spindle                  {00:1, 01:1, 10:1, 11:3}    1 v u uv     2   -0.693147181  -1.386294361
B5  double sphere            -- no gamma_C --            n/a         n/a       n/a           n/a
B2  K1 both filled           {01:2, 10:2, 11:1}          v u uv       2   -0.756573586  -0.767507880
B1p K1-bridged               {01:3, 10:3}                v u          1   -0.693147181  -0.693147181
B1q K1-bridged + spectator   {00:1, 01:3, 10:3}          1 v u        2   -0.741029583  -0.767507880
B1s K1 subdivided            {01:5, 10:5, 11:1}          v u uv       2   -0.724759919  -0.767507880
```

**Eight of the nine entries are EXACT closed forms, not quadrature**, and each is stated with
the identity that produces it:

```
   B0a  U  = m(2/9 + (3/9)y + (4/9)x)      Cassaigne-Maillot   -0.747659833081
   B3/B1/B2 U = m(0.4 + 0.4x + 0.2y)       Cassaigne-Maillot   -0.756573585640
   B1q  U  = m(1/7 + (3/7)x + (3/7)y)      Cassaigne-Maillot   -0.741029582571
   B1s  U  = m(5/11 + (5/11)x + (1/11)y)   Cassaigne-Maillot   -0.724759919461
   B4   U  : the Jensen-in-x integrand's max is always the second branch (their squares
             differ by 0.2222 + 0.1111 cos y > 0), so lambda = log(1/2) = -0.693147180560 EXACTLY
   B1p  U  = log max(1/2,1/2)              Jensen              -0.693147180560 EXACTLY
   SENSE C, 3 classes = m(0.4+0.3x+0.3y)   Cassaigne-Maillot   -0.767507880358
   SENSE C, 4 classes: (1+x+y+xy)/4 = (1+x)(1+y)/4, so lambda = log(1/4) = -1.386294361120 EXACTLY
   B0b  U  : genuinely 4-term and does not factor. QUADRATURE ONLY, -0.810930216216,
             stable to 12 places across n = 5e4, 2e5, 8e5, 3.2e6
```

Closed forms verified against **direct schedule-B simulation** at `f = 1.0, c = sqrt(2)`,
`N = 2e6` — worst deviation `3.0e-06`, every carrier:

```
    B0a  direct -0.747660592  exact -0.747659833   dev 7.6e-07
    B0b  direct -0.810929681  quad  -0.810930216   dev 5.4e-07
    B3   direct -0.756576560  exact -0.756573586   dev 3.0e-06
    B1   direct -0.756576560  exact -0.756573586   dev 3.0e-06
    B4   direct -0.693146936  exact -0.693147181   dev 2.4e-07
    B2   direct -0.756576560  exact -0.756573586   dev 3.0e-06
    B1p  direct -0.693149769  exact -0.693147181   dev 2.6e-06
    B1q  direct -0.741029963  exact -0.741029583   dev 3.8e-07
    B1s  direct -0.724760493  exact -0.724759919   dev 5.7e-07
```

**`lambda` is not constant across the family**: six distinct values in SENSE U, three in
SENSE C, plus one carrier on which formation has no datum at all.

### 4.3 THE SEPARATION THEOREM — the general statement the family forced out

> **THEOREM (S4-1).** Let a carrier have designated loops `gamma_F, gamma_C`, and let `S` be
> the set of vertex classes `(a,b) in {0,1}^2` carrying positive weight, with characters
> `1, u, v, uv`. Then `G = < chi_x / chi_y : x,y in S >` has **rank 2 iff `|S| >= 3`**, and
> for `|S| = 2`:
>
> ```
>   S = {1,u}  or {uv,v}   ->  G = <u>     lambda sees W_F only
>   S = {1,v}  or {uv,u}   ->  G = <v>     lambda sees W_C only
>   S = {u,v}              ->  G = <u/v>   lambda sees the PRODUCT  W_F·W_C only
>   S = {1,uv}             ->  G = <uv>    lambda sees the RATIO    W_C/W_F only
> ```
>
> and for `|S| = 1`, `G = {1}` and there is no formation, ever.
>
> **PROOF.** The exponent vectors of `1, u, v, uv` are the four corners of the unit square.
> `rank G` is the rank of the lattice generated by the differences of the corners in `S`.
> Any three corners of a square are affinely independent, so `|S| >= 3` gives rank 2;
> `|S| = 2` gives the single difference vector listed; `|S| = 1` gives the zero lattice. ∎

Checked by enumeration over **all 15 non-empty subsets**: `rank 2` on 5, `rank 1` on 6,
`rank 0` on 4 — matching the prediction `5 / 6 / 4` exactly.

**This subsumes W-02's four support classes as the `K1` special case** (no `(0,0)` vertex
exists on `K1`, so its three classes are `uv, u, v` and the four cases are the subsets of
those). It also shows W-02's list is *incomplete as a general statement*: a carrier with a
vertex on **neither** loop has a fourth character `1`, and a fifth two-element case
`S = {1,uv}` — the **ratio** `W_C/W_F` — which cannot arise on `K1` at all.

### 4.4 Q7 — THE TRANSITION. FOUR CONTROLS.

**CONTROL 1 — THE FILL CONTROL. `chi` moves; `lambda` does not.**

`B1` and `B2` have the *identical* 1-skeleton and the *identical* designated loops. Filling
the second triangle changes the topology and nothing else:

```
                   B1 (K1)        B2 (both filled)
   chi                0                 1
   b1                 1                 0
   b2                 0                 0
   F                  1                 2
   curvature          1                 2
   flat               1                 0
   vertex classes  {01:.4, 10:.4, 11:.2}  IDENTICAL
   rank G             2                 2
   lambda        -0.756573585634   -0.756573585634       |difference| = 0.000e+00
   direct (N=2e6) -0.756576560141  -0.756576560141       |difference| = 0.000e+00
```

**`chi` moves from 0 to 1, `b1` from 1 to 0, the flat holonomy is destroyed outright — and
`lambda` does not move by one bit.**

**CONTROL 2 — THE INCIDENCE CONTROL. Every topological invariant equal; `lambda` differs.**

`B1` (two triangles **wedged** at `v0`) and `B1p` (two triangles **bridged** by an edge) are
homotopy equivalent — both are `S^1` — and agree on every invariant computed:

```
                   B1              B1p
   chi                0             0
   b0, b1, b2      1, 1, 0       1, 1, 0
   invariants         2             2
   curvature          1             1
   flat               1             1
   vertex classes  {01:.4,10:.4,11:.2}   {01:.5, 10:.5}
   rank G             2             1
   lambda        -0.756573586   -0.693147181          |difference| = 6.343e-02
```

**On `K1` the rate separates `W_F` from `W_C`. On the bridged carrier — same `chi`, same
Betti numbers, same homotopy type, same invariant split — it sees only the product
`W_F·W_C`.** The whole difference is **one vertex of incidence**: whether the two loops meet.

**CONTROL 3 — THE SUBDIVISION CONTROL. Homeomorphic carriers; `lambda` differs under SENSE U.**

`B1s` is `K1` with every edge subdivided once — the **same topological space**.

```
   SENSE U:  -0.756573585634  vs  -0.724759919464     |difference| = 3.181e-02
   SENSE C:  -0.767507880351  vs  -0.767507880351     |difference| = 0.0e+00
```

**Under SENSE U the *value* of `lambda` is not even a homeomorphism invariant.** Under SENSE C
it is — but only because SENSE C fixes by hand the very thing SENSE U reads off the complex.
**This is an OPEN choice and it is load-bearing; see the CHOICE LEDGER.**

**CONTROL 4 — THE SPECTATOR CONTROL. The pinch is sufficient for separation, not necessary.**

Subdividing `B1p`'s bridge adds a vertex `m` lying on **neither** loop:

```
   B1p  classes {01:.5, 10:.5}                     rank G = 1  (product only)
   B1q  classes {00:.143, 01:.429, 10:.429}        rank G = 2  (separates)
```

**A vertex that touches neither loop restores full separation with no pinch anywhere.** This
cuts directly against the tidiest reading of `K1`'s design and is recorded because it does.

**Is the pinch a distinguished point of the family?**

```
stage                             chi  b1  b2  curv flat  rkG   lambda (SENSE U)
B0a ring torus, loops disjoint      0   2   1    8    2    2     -0.747659833
B0b ring torus, loops meet          0   2   1    8    2    2     -0.810930216
B3  horn torus (the pinch)          1   1   1    7    1    2     -0.756573586
B1  K1 (the pinch, as handed)       0   1   0    1    1    2     -0.756573586
B4  spindle                         2   1   2    2    1    2     -0.693147181
B5  double-covered sphere           2   0   1    1    0   n/a     n/a
```

**Answer: no — the pinch is an ordinary point of the family, and the *degenerate end* is the
distinguished one.** `rank G = 2` from the ring torus all the way to the spindle; the only
structural break is at `B5`, where `b1 = 0` and there is no formation datum to compute.
Nothing happens at the pinch that does not also happen on either side of it.

**But the pinch is distinguished in one exact respect, and it is the respect S1 claimed.**
The honest discrete horn torus `B3` (`chi = 1`, `b = (1,1,1)`, built as an octahedron with
its poles identified) and `K1` (`chi = 0`, `b = (1,1,0)`) have **the identical vertex class
profile `{(1,1):1, (1,0):2, (0,1):2}` and the identical `lambda` to twelve places.** The
brief's identification of `K1` as "the pinch's class" is **confirmed at computation** —
`K1` reproduces the horn torus's formation data exactly, while carrying neither its `chi`
nor its `b2`. That is a real result about what `K1` is a stand-in for.

### **ARMED FALSIFIER F-B — VERDICT: DOES NOT FIRE ON ITS LITERAL CONDITION. ITS TARGET CLAIM IS CONFIRMED.**

I report both readings and do not choose between them for the auditor.

- **Literal condition — "`lambda` is constant across the entire family":** **FALSE.**
  Six distinct values in SENSE U, three in SENSE C, and one carrier (`B5`) on which the
  formation datum does not exist. **F-B does not fire.**

- **Target claim — "the carrier's topology is INERT, the geometric picture is decoration":**
  **CONFIRMED, by controls built specifically to test it.** Control 1 moves `chi`, `b1`,
  `F`, and the curvature/flat split with `lambda` frozen to the bit. Control 2 holds every
  topological invariant fixed and moves `lambda`. Control 3 shows the value is not a
  homeomorphism invariant. **`lambda` is a function of the loop-incidence pattern and the
  vertex weights. It is not a function of the topology.**

  The one surviving topological statement is negative and narrow: `b1 >= 1` is what makes
  a *flat* holonomy exist at all, and `rank(d2) >= 1` is what makes a *curvature* exist.
  **Topology decides whether the formation datum exists and whether the second holonomy is
  free or determined. It does not touch the rate.**

**The honest summary of S4-B: formation sees the pinch, not the topology — and it sees a
spectator vertex just as well as it sees the pinch.** The revolution-family picture
(ring/horn/spindle/sphere) motivated ten good constructions and is, as a *classification* of
formation behaviour, decoration. What classifies formation behaviour is Theorem S4-1: which
of four vertex classes carry weight.

### 4.5 Q8 — DOES ANY FORMATION QUANTITY RELATE TO `2*pi*chi`?

**Compared to nothing outside the corpus, as instructed.**

```
carrier                        chi    2 pi chi     lambda (SENSE U)
B0a ring torus, disjoint         0     0.00000      -0.747659833
B0b ring torus, meeting          0     0.00000      -0.810930216
B3  horn torus                   1     6.28319      -0.756573586
B1  K1                           0     0.00000      -0.756573586
B4  spindle                      2    12.56637      -0.693147181
B2  K1 both filled               1     6.28319      -0.756573586
B1p K1-bridged                   0     0.00000      -0.693147181
B1q K1-bridged + spectator       0     0.00000      -0.741029583
B1s K1 subdivided                0     0.00000      -0.724759919

   chi = 0 carries SIX distinct lambda values
   chi = 1 carries one ; chi = 2 carries one
```

**ANSWER: NO, and the ground is a theorem rather than a failed correlation.**

`chi -> lambda` is **not a well-defined map.** Control 1 exhibits two carriers with
`chi = 0` and `chi = 1` and *bit-identical* `lambda`; Control 2 exhibits two carriers with
`chi = 0` and `chi = 0` and `lambda` differing by `0.0634`. **Therefore no function
`lambda = g(2*pi*chi)` can exist,** and no amount of sampling could produce one that is not
an artefact of the sample.

**I refused the correlation hunt deliberately, and say so.** With nine carriers and three
`chi` values, a two-parameter fit would have found something. The fill control is the
decisive test because it moves `chi` while holding *everything the functional can see*
fixed; that is a controlled experiment and a correlation table is not.

**Where `2*pi` does enter.** Only as the normalisation of Haar measure on `U(1)`:
`lambda = (1/2pi)^{dim H} int_H log|Z|`. That normalisation is a **probability** measure, so
rescaling it is a no-op — `2*pi` cancels out of every quantity on this page and is never a
factor in an answer. Gauss–Bonnet's `chi` multiplies a curvature **density** integrated over
a manifold; there is no density in `lambda`. The curvature `W_F` enters `lambda` only through
its **holonomy**, and §3.1 established that under the carrier's own clock the holonomy's
contribution to the rate is **almost-everywhere constant in `W_F`**. There is no channel
through which a Gauss–Bonnet-shaped relation could reach `lambda`, and none is observed.

**Grade: NEGATIVE, at computation, with a mechanism.** Not INDETERMINATE.

---

## 5. WHAT THIS DOES TO W-02's REOPEN CONDITIONS

`REGISTER_V001.md` W-02 reopens if, among others, *"the rate `lambda` fails to vary with the
connection at S4"* or *"a lineage-independent lane fails to reproduce `lambda = -0.767026`
or the character-ratio criterion."*

- **`lambda = -0.767026` is reproduced** (this lane, `N = 200000`, independent
  implementation: `-0.767026255`), **and its stated provenance is refuted** — see §7.
- **The character-ratio criterion is reproduced and generalised** — Theorem S4-1, verified
  over all 15 class subsets.
- **"`lambda` fails to vary with the connection"** — this is the closest thing to a fired
  reopen condition on the page, and it is **half-true, so I state it precisely rather than
  resolving it in either direction.** Under the **uniform** schedule `lambda` varies fully
  and smoothly over the torus. Under the **canonical clock — S3's primary schedule** —
  `lambda` **is constant on a set of full Haar measure** and varies only on a dense null set.
  Whether that counts as "failing to vary" is a judgement about which schedule is physical,
  and S3 §6(e) itself records *"I do not know which is physical"* (FLAG F5). **I do not know
  either, and I decline to decide it here.** It is flagged.
- **This lane is not lineage-independent** (custody §4). Grade below accordingly.

---

## 6. TOY_SEPARATION

**Required by the brief; a toy not labelled a toy is a reportable defect.**

**ACTUAL SURFACE — load-bearing, and I stand behind each as a real object:**

1. **The general overlap `Z_k = sum_v u^{k a_v} v^{k b_v} p_v`** (§2). Derived from W-01's
   own operators for an arbitrary carrier; reduces to W-01 on `K1` and to the corrected
   criterion's three characters. This is the actual surface of the whole page.
2. **The exact identity `|Z_k|^2 = ... + 2p0q cos(kc) + 2p0r cos(kf) + 2qr cos(k(f+c))`**
   (§2.1) and the rank-2 Fourier proposition built on it (§3.2). Exact, checked at `2.3e-14`.
3. **`lambda_B = int_{L^perp} log|Z| d(Haar)`, hence a function of `L` alone** (§2.2, §3.3).
   Theorem, and verified to `3.3e-16` across 49 lattice classes.
4. **The closed forms**: Cassaigne–Maillot for the generic rate, Jensen for every resonant
   and two-support rate, exact finite averages for roots of unity. These are exact values,
   not fits, and each is cross-checked against direct simulation.
5. **Theorem S4-1, the separation theorem** (§4.3). Proved and enumerated.
6. **The four controls** (§4.4). These are the actual experiment of S4-B. Everything else in
   S4-B exists to make them possible.
7. **The defect at `(2.0, 1.1)`** (§7). A checkable fact about a sealed artifact.

**ILLUSTRATIONS — built to be looked at, not leaned on. Labelled here:**

- **The ten carriers as a "degeneration family."** They are ten legitimate CW complexes with
  correctly computed homology; **the claim that they are *the* discrete realization of the
  revolution family is an illustration, not a theorem.** There is no theorem here that a
  discrete complex "is" a horn torus; the identifications are by construction and named
  analogy. Their role is to be *different from each other in controlled ways*, and for that
  they are the actual surface. Their role as *stages of a geometric degeneration* is
  illustrative.
- **`B0a` vs `B0b` (the two loop placements on the ring torus).** These exist to show the
  incidence effect on a carrier that is not `K1`. The specific `3x3` grid is arbitrary.
- **The `1200 x 1200` grid picture of `lambda_A`** (§3.1). The analytic critical-point
  solution is the actual surface; the grid is a check on it.
- **The `241`-step simplex sweep** (§3.4). The Cassaigne–Maillot closed form is the actual
  surface; the sweep confirms where the minimum sits.

**NOT CLAIMED ANYWHERE, and stated so it cannot be picked up later:** that `lambda` is a
physical rate, that any number on this page has units, that the carriers are models of
spacetime, or that the revolution family has any status beyond having suggested what to
build.

---

## 7. A DEFECT OF RECORD IN S3 AND ITS AUDIT

**S3's headline test connection is exactly resonant, and both the build and the audit state
the opposite.**

`S3_THE_CROSSING_V001.md` §6(f) states, quoted whole including the clause that decides it:

> `f = 2.00000  c = 1.10000   lambda_B = -0.767026    f, c, 2pi rationally independent: orbit dense in T^2`

and §4.3 states:

> `torus mean of log|Z| (Weyl equidistribution prediction)     = -0.767507        AGREES`

`f` and `c` are **both rational**, hence rationally **dependent**:

```
   -11 * f + 20 * c  =  -11*(2.0) + 20*(1.1)  =  0     EXACTLY, in exact arithmetic
   gcd(11, 20) = 1   ->  primitive relation  u^11 v^20 = 1
```

The orbit closure is the **kernel circle of `chi_(11,20)`**, one-dimensional, not `T^2`.

```
   exact subtorus rate m(r z^31 + p0 z^11 + q), (m,n)=(11,20)   :  -0.767014993
   direct schedule-B average, N =    200000                     :  -0.767026255
   direct schedule-B average, N =   1000000                     :  -0.767016017
   direct schedule-B average, N =   4000000                     :  -0.767015363
   direct schedule-B average, N =  16000000                     :  -0.767014977
   full-torus Mahler measure  m(p0 + q x + r y)                 :  -0.767507880
   S3 build and S3 audit report (N = 200000)                    :  -0.767026
```

**The reported number is correct. Its justification is not.** The run is converging to the
**subtorus** value and agrees with the full-torus mean only to `4.9e-04`, by coincidence.
Confirmed by detuning `c` off the rational:

```
   c = 1.1 + 0        lambda_B(N=4e6) = -0.767015363    (subtorus)
   c = 1.1 + 1e-9                     = -0.767015979
   c = 1.1 + 1e-7                     = -0.767463181    (leaving)
   c = 1.1 + 1e-5                     = -0.767508449    (full torus)
   c = 1.1 + 1e-3                     = -0.767507996
```

**Scope of the defect.** Every row of S3 §5.7's ready-state table was computed at this same
connection and is therefore a subtorus value:

```
   p (p0,q,r)          S3 @ N=5e4     exact (11,20)     full-torus m
   (0.4,0.3,0.3)        -0.767043      -0.767014993     -0.767507880
   (0.5,0.0,0.5)        -0.692761      -0.693147181     -0.693147181
   (1/3,1/3,1/3)        -0.776733      -0.776771395     -0.775546341
   (0.8,0.1,0.1)        -0.223149      -0.223143551     -0.223143551
   (0.1,0.45,0.45)      -0.727616      -0.727613407     -0.727723488
   (0.2,0.5,0.3)        -0.693152      -0.693147181     -0.693147181
```

Every row is converging to the `(11,20)` column, not the torus column. In the two-support and
log-max cases the two coincide, which is why nothing looked wrong.

**Cutting against my own finding, stated because it does:** the defect changes no verdict in
W-02. Theorem S3-2's conclusion (`sum(1-z_n) = infinity`, `lambda < 0`, exponential decay)
holds identically on a subtorus, because the proof only needs `|Z| < 1` on a set of positive
Haar measure of the closure. **The escape is untouched. What is wrong is the stated reason a
particular number is that number**, and — materially — the generality of §6(f)'s claim that
the rate depends on "the arithmetic of the holonomies": it does, but its own worked example
is on the wrong side of that classification.

**S3's other §6(f) rows, re-derived exactly:**

```
   f=c=2.0        -> L = <(1,1)>, exact lambda = m(0.3z^2+0.4z+0.3) = -1.203972804   (S3: -1.203587 at finite N)
   f=pi, c=pi/2   -> orbit order 4, exact lambda = -0.804718956                      (S3: -0.804719)
   f=pi, c=3pi/2  -> orbit order 4, exact lambda = -0.804718956                      (S3: -0.804719)
```

**And the S3 audit's disagreement at `f = 3.14159, c = 1.57080` is resolved, in the build's
favour on the number and the audit's favour on the method.** The audit reported
`-0.860699` at `N = 200000` against the build's `-0.804719`. Neither is the limit:

```
   2c - f = 1.0e-05   -> near-resonance (m,n)=(1,2), transient length 2pi/1e-5 = 628319 circuits
   lambda_B(N =   200000) = -0.860699        (the audit's number, reproduced exactly)
   lambda_B(N =  1000000) = -0.733898
   lambda_B(N =  4000000) = -0.771307
   lambda_B(N = 20000000) = -0.767625
   generic value          = -0.767507880
```

**The truncated decimals `3.14159` and `1.57080` are not `pi` and `pi/2`; they are generic,
and the limit is the generic value.** The build's `-0.804719` is the value for *exact*
`pi, pi/2` (a genuine order-4 orbit) and is correct for that connection; the audit's
`-0.860699` is a near-resonance transient at `N` a third of the way through one beat period.
**Both are right about their own computation and both mislabel it.**

---

## 8. CHOICE LEDGER

| # | Choice | Alternatives | Why | Status |
|---|---|---|---|---|
| C1 | Report **both** schedules as an axis rather than fixing one | fix B (S3's "primary"); fix A | The carried correction says the escape is conditional on the schedule; the two give qualitatively different maps, and §3.1 shows the difference is the whole content of Q1 | **OPEN** — S3 FLAG F5 says "I do not know which is physical"; neither do I |
| C2 | Report **both** senses of "hold the ready state fixed" across carriers (SENSE U uniform-on-vertices, SENSE C fixed class weights) | pick one | Control 3 shows the choice decides whether `lambda` is a homeomorphism invariant. Picking one silently would have produced a clean answer that the other choice contradicts | **OPEN and load-bearing** |
| C3 | Ready state `(0.4,0.3,0.3)` as the default on the connection axis | uniform `(1/3,1/3,1/3)`; K1's published `(1/2,0,1/2)` | It is S3's own table state, so every number here is comparable to the corpus. K1's published state is `|S|=2` and would have hidden the F-A structure entirely | closed |
| C4 | The formation datum is an **ordered pair (face boundary, free cycle)**, `gamma_F` bounding and `gamma_C` not | any two loops; all pairs; a canonical choice from incidence | This is what makes `W_F` a curvature and `W_C` flat, which is `K1`'s stated reason for existing (S1 §5). Verified per carrier rather than assumed | closed |
| C5 | On carriers with more than two invariants, **set all non-designated invariants to 1** | random; equal; sweep them | "The minimal excitation carrying the designated pair" is stateable and defensible; sweeping 8 extra invariants on the ring torus would answer a different question | **OPEN** — the non-designated invariants are untested and could matter |
| C6 | Discrete "ring torus" = `3x3` square-grid torus (regular CW) | minimal 1-vertex torus (not regular); 7-vertex triangulation | Regularity is `S1`'s own standard (S1 §1). The 1-vertex torus is not a regular CW complex | closed |
| C7 | Discrete "horn torus" = octahedron with poles identified, `chi=1`, `b=(1,1,1)` | the brief's "cycles meeting at one vertex, `b1=1`" (= `K1`, `chi=0`) | **DEPARTURE FROM THE BRIEF'S LIST, declared.** The continuum horn torus is `T^2` with a longitude collapsed = a sphere with two points identified, which has `chi = 1`, not 0. I built **both** and report that they give identical formation data — which is the finding, and would have been invisible had I built only one | closed, declared |
| C8 | Discrete "double-covered sphere" = two 2-cells on one 1-skeleton (`chi=2`, `b2=1`) | two spheres glued at two points (= my `B4`, which is what the `R->0` limit of the abstract surface actually is) | **DEPARTURE DECLARED.** I built the brief's version as `B5` and the abstract-surface version as `B4`, and report both. `B5` is where formation dies (`b1=0`); `B4` is where it survives. Calling only one of them "the double sphere" would have decided Q7 by naming | closed, declared |
| C9 | Added `B1p`, `B1q`, `B1s`, `B2` — four carriers **not** in the brief's list | build only the four named stages | The brief's four stages cannot separate "topology" from "incidence": they covary. The controls are the only way to answer F-B and Q8 non-circularly | closed |
| C10 | `b_*` computed over `Q` (matrix ranks), not over `Z` | Smith normal form over `Z` | numpy has no exact SNF here and none of these carriers has torsion. **Limitation recorded in the FLAG BLOCK, not hidden** | **OPEN** |
| C11 | `a_v in {0,1}` — a vertex visited twice by a loop still counts once | count multiplicity | W-01's operator is *"multiply by `W(gamma)` at vertices on `gamma`"*, which is a set membership. Multiplicity would be a different operator and a different construction | closed; but note it makes Theorem S4-1's "corners of a square" argument exact, and multiplicity would break it |
| C12 | Q8 answered by a **control**, not by fitting `lambda` against `2*pi*chi` over the family | regression / ratio table | With 9 carriers and 3 `chi` values a fit finds something. Custody §5 bars selecting structure by landing on a value; fitting an internal quantity to `chi` is the same move one level down | closed |

---

## 9. IMPORT AUDIT

Anything used that the corpus does not define.

| Import | What it is | Source (named, since the corpus carries no digest for it) | Does the finding survive without it? |
|---|---|---|---|
| **Weyl equidistribution** | `(u^n, v^n)` equidistributes on the closure of the group it generates | classical; **already imported by S3** (Theorem S3-2's proof) and by the S3 audit's corrected theorem | Not a new import. Without it, `lambda_B` is only defined as a limit and every value here would be a finite-`N` estimate; all the *comparisons* (controls, F-A, F-B) survive since they use the same estimator on both sides |
| **Pontryagin duality** (closed subgroups of `T^2` = annihilators of sublattices of `Z^2`) | gives `H = L^perp`, hence the lattice-determination of `lambda_B` | classical | **Partially.** The three-tier classification survives as an empirical statement (verified to `3.3e-16` across 49 classes); only the word "theorem" needs it |
| **Jensen's formula / Mahler measure**, `m(P) = log|lead| + sum_{|root|>1} log|root|` | every exact resonant and two-support value | classical | Yes — each value was independently confirmed by direct schedule-B simulation to `<= 2.7e-06` |
| **Cassaigne–Maillot formula** for `m(a+bx+cy)` | the exact generic rate, e.g. `-0.767507880` | Cassaigne & Maillot, on the Mahler measure of `a+bx+cy`; author-named, no digest | Yes — cross-checked against an **independent** route (exact Jensen in `x`, quadrature in `y`), 12 random points, max deviation `6.4e-11`, and against a `4000^2` grid |
| **Bloch–Wigner dilogarithm `D(z)`** and `Li_2` | inside Cassaigne–Maillot | classical; implemented from scratch here (no mpmath/scipy/sympy available) | Yes — validated at `Li_2(1), Li_2(-1), Li_2(1/2), Li_2(2)`, `D(e^{i pi/3}) = 1.014941606410` against the known `1.014941606409`, and `m(1+x+y) = 0.323065947219` against the known value |
| **`m(1+x+y) = 0.323065947219`** | the equilateral-simplex minimum | classical constant, named not sealed | Yes — the same value is produced by the independent Jensen-in-`x` route to `3e-12` |
| **The revolution family** (ring/horn/spindle/sphere) | motivation for which complexes to build | the brief; geometric folklore | **Yes, entirely.** It is used to *choose* constructions and never as a step. Every homology number is computed from incidence. Removing the motivation removes ten names and no results |

**No measured physical constant appears anywhere on this page.** No coupling. No target
value. Nothing is fitted to anything.

---

## 10. COMPUTATION LEDGER

All code written for this page, `python3 3.9.6` + `numpy 2.0.2`, no other libraries
available or used. Files in this session's scratchpad:

```
s4lib.py             the library: Li_2, Bloch-Wigner, mahler1 (Jensen), Cassaigne-Maillot,
                     mahler2_bilinear (exact in x, quadrature in y), Z_of_k, lambda_direct,
                     lambda_generic, lambda_resonant, lambda_rank1, lambda_finite_orbit,
                     G_lattice_rank, class CW (d1, d2, Betti, gauge/invariant report)
s4_v0_validate.py    S1/S2/S3 reproduced from scratch; special functions validated
s4_A_connection.py   Q1-Q4: the exact modulus, schedule A map, Fourier support, schedule B
                     sweep, resonance loci, rational table, support axis, simplex sweep
s4_A2_resonance.py   the (2.0,1.1) defect; the exceptional set; the three tiers; firing locus
s4_A3_lattice.py     lattice-determination (3 tests); corrected discontinuity test;
                     F-A verdict table; self-check
s4_B2_carrier.py     the ten carriers, homology, loop verification, classes, lambda,
                     the four controls, Theorem S4-1 enumeration, Q8
```

**Every reported number and how it was produced, with its independent re-derivation:**

| Quantity | Primary method | Independent re-derivation | Agreement |
|---|---|---|---|
| `Z_k` on `K1` | closed form | direct matrix action on `C^5`, random section phases | `1.487e-14` |
| `|Z_k|^2` identity §2.1 | analytic expansion | general routine, 2000 random `(f,c,k)` | `2.287e-14` |
| generic `lambda_B` `-0.767507880` | Cassaigne–Maillot (exact) | (a) Jensen-in-`x` + quadrature, (b) `4000^2` grid, (c) 120-point direct sweep | (a) `1.3e-11`, (b) `5.2e-07`, (c) mean matches, std `2.1e-05` |
| per-carrier `lambda` (8 of 9) | exact: Cassaigne–Maillot / Jensen / factorisation | quadrature route | `<= 4.9e-11` on the 3-class ones; `1.1e-16` on `B4`; exact `log(1/4)` vs quadrature `3.5e-06` on the 4-class SENSE C |
| `lambda` on resonance `(m,n)` | `mahler1` (Jensen, exact) | direct schedule-B, `N=2e6`, on the locus | `<= 2.7e-06` (16 loci) |
| `lambda` for `S={0,C}`,`{0,F}`,`{F,C}` | `log max(.,.)` exact | direct, `N=2e6` | `<= 1.65e-05` |
| `lambda` at roots of unity | exact finite orbit, integer exponents | for `f=0,c=2pi b/12`: closed form `log 0.7 + (1/12)log(1-(3/7)^12) = -0.3566781436122585` vs computed `-0.3566781436122585` | **exact** |
| `lambda` for `(11,20)` `-0.767014993` | `mahler1` | direct `N=1.6e7` = `-0.767014977` | `1.6e-08` |
| `d=3` partial resonance `-0.798965257` | `lambda_rank1` | direct `N=8e6` at two different generic `c` | `6.4e-08`, `1.9e-07` |
| simplex minimum `-0.775546341` | `log(1/3)+m(1+x+y)` exact | `241`-grid sweep found `-0.775535` at `(0.332,0.332,0.336)` | `1.1e-05` (grid resolution) |
| all Betti numbers | ranks of `d1`,`d2` over `Q` | `chi = V-E+F = b0-b1+b2` checked; `curv+flat=inv` asserted | exact on all 10 |
| `lambda` per carrier | closed form on class weights | direct schedule-B `N=2e6`, `f=1.0,c=sqrt2` | `<= 3.0e-06` (9 carriers) |
| lattice-determination | HNF grouping of 143 rational pairs | max spread within class | `3.331e-16` |
| Theorem S4-1 | proof | enumeration of all 15 subsets: `5/6/4` vs predicted `5/6/4` | exact |
| winding periodicity | — | 500 random samples, `k,l in [-5,5]` | `2.33e-15` |

**Actual hit counts for every sweep run** (custody §3):

```
  2000 random (f,c,k) modulus checks                                    2000 run, 0 failures
  120 random connections, schedule B, N=4e5 each                        120 run, 120 within 1e-2 of the
                                                                        generic value, 0 outliers
  1200 x 1200 grid of lambda_A                                          1440000 cells, 262 with |Z_1|<0.01
  16 primitive resonance loci, exact vs direct                          16 run, max dev 2.7e-06
  143 rational pairs with denominator 12                                143 run, 49 lattice classes
  638 rational pairs with denominator <= 12 (finite lambda)             638 run, 380 above / 258 below generic
  29400 simplex points                                                  29400 run, 7260 triangle / 22140 log-max
  15 vertex-class subsets                                               15 run, 5 rank-2 / 6 rank-1 / 4 rank-0
  10 carriers x (homology + loop verification + lambda)                 10 run, 10 pass
  500 winding-periodicity samples                                       500 run, max dev 2.33e-15
  110 predicted exact zeros among k<=5 preimages                        110 found, 110 predicted
```

**Numerical hygiene, self-imposed and checked:**

1. `Li_2`'s branch on the real axis `> 1` is the conjugate of the usual convention. It never
   matters here: Cassaigne–Maillot evaluates `D` only at `(a/b)e^{i gamma}` with
   `gamma in (0,pi)` strictly, never real. Confirmed against the independent route.
2. Phase drift at large `k`: worst case `N = 2e7`, `|kf| ~ 6e7`, ulp `~7.5e-9` rad. Bounded
   by re-deriving the `(11,20)` value in exact closed form — `1.6e-08` agreement.
3. `mahler2_bilinear`'s quadrature is weakest where the integrand has a log singularity;
   worst observed `3.5e-06` on `S={0,F}`. **Every two-character case is reported from the
   exact `log max` closed form instead**, and the discrepancy is recorded rather than
   smoothed.
4. Roots-of-unity averages use **integer exponent arithmetic mod the orbit order** — no
   floating-point drift is possible in them.

---

## 11. FLAG BLOCK

**SEAL TALLY.** 20 of 20 inherited sidecars verified `OK` from each artifact's own directory
before any work; `S4_THE_MEASUREMENT_V001.md` did not previously exist. This artifact is
sealed from its own directory on creation. **Self-excluded from its own sweeps.**

**FENCE LINE — HELD.** No quantity on this page is compared to any measured physical
constant. Alpha is not engaged: no coupling, no target, no measured number, no
target-driven selection. `lambda`, `G`, `L`, the class weights and every Mahler measure are
free, connection- and carrier-determined, and are fitted to nothing. No continuum limit is
taken; the revolution family appears only as motivation for which complexes to build, and
every homological number is computed from incidence matrices. **Alpha may constrain
structure; alpha may not select structure** — nothing here selects.

**FALSIFIER STATUS.**

- **F-A — DOES NOT FIRE ON FULL SUPPORT; FIRES ON THREE OF FOUR SUPPORT CLASSES.**
  With `p0 > 0` the rate depends on `W_F` and `W_C` separately, proved by rank-2 Fourier
  support. With `p0 = 0` it depends on the **product `W_F·W_C` alone** — F-A's condition
  exactly. With `q = 0` or `r = 0` it depends on one holonomy alone. **The carrier's
  separation of curvature from flat holonomy does work if and only if the ready state puts
  weight on the pinch vertex.** Additionally, under the canonical clock the dependence is
  arithmetic, not metric: `lambda_B` is a function of the relation lattice `L`, constant on
  a full-measure set.
- **F-B — DOES NOT FIRE ON ITS LITERAL CONDITION; ITS TARGET CLAIM IS CONFIRMED.**
  `lambda` is not constant across the family (six values under SENSE U). But the fill control
  moves `chi` from 0 to 1 with `lambda` bit-identical, and the incidence control holds every
  topological invariant fixed while `lambda` moves by `0.0634`. **The carrier's topology is
  inert; formation sees the loop-incidence pattern.** The one narrow survival: `b1 >= 1` and
  `rank(d2) >= 1` are what make the formation datum exist at all.

**DEFECTS FOUND IN MY OWN DRAFT ON SELF-CHECK** (custody §4; found by me, before sealing):

- **D1.** My first Q3 discontinuity protocol held `f` rational and moved only `c`. That is
  confounded — the relation `(3,0)` survives — and it would have made a partial-resonance
  value look like the generic limit. Corrected in §3.3, with the wrong protocol and its
  (correct-for-what-it-computed) number both shown.
- **D2.** My first horn-torus CW complex mixed N-hemisphere and S-hemisphere edges in its
  face attachments. It happened to give the right Betti numbers, which is exactly why it was
  dangerous. Rebuilt correctly (`s4_B2_carrier.py`) with all loop properties verified against
  `d1` and `d2` rather than asserted.
- **D3.** My first carrier pass asserted the designated loops by naming vertex sets. That is
  the assertion-without-verification failure custody §1 exists to catch. Rewritten so each
  loop is a signed edge chain and `gF bounds / gC does not bound / independent` are
  **computed**.
- **D4.** `lambda_generic`'s quadrature carries up to `3.5e-06` error wherever the integrand
  has a log singularity — which happens exactly when the bilinear polynomial **factors**, and
  those are precisely the cases with an exact closed form. My first carrier table quoted the
  quadrature values `-1.386290895` for `B0b` and `B4` under SENSE C; **the exact value is
  `log(1/4) = -1.386294361120`**, because `(1 + x + y + xy)/4 = (1+x)(1+y)/4`. Corrected in
  §4.2, where eight of the nine carrier rates are now stated as exact closed forms with the
  identity that produces each, and the one genuinely quadrature-only row (`B0b` SENSE U) is
  labelled as such and shown stable to 12 places under a 64-fold refinement.
- **D5.** Betti numbers are computed over `Q`. Torsion in `H_1` would be invisible. No
  carrier here has any, but the computation cannot certify that — declared, not hidden
  (CHOICE LEDGER C10).

**FLAGS CARRIED FORWARD.**

- **F1 (inherited, unresolved).** Which schedule is physical. S3 FLAG F5 says "I do not know";
  S4 confirms the two schedules give **qualitatively different answers to Q1 and to W-02's
  reopen condition** — under B, `lambda` is constant a.e. This is now more consequential than
  when S3 flagged it, and it is still open.
- **F2.** SENSE U vs SENSE C for holding the ready state fixed across carriers is unresolved
  and load-bearing (Control 3). A principled rule for transporting a ready state between
  carriers does not exist in the corpus.
- **F3.** The non-designated invariants (8 of 10 on the ring torus) were set to 1 and never
  swept. Their effect on `lambda` is untested.
- **F4.** Theorem S4-1 assumes each loop meets each vertex at most once as a *set*
  (CHOICE LEDGER C11). Carriers where a loop's multiplicity matters are outside its scope.
- **F5.** `B5` is the only carrier where topology alone decides an outcome, and it decides it
  by making the datum not exist. That is a weaker sense of "topology matters" than S1's
  framing suggests, and it should not be reported as stronger.

**GRADE: ADVERSARIALLY-UNCHECKED — this is a BUILD with no audit yet, and it is not
lineage-independent** (custody §4). Build, brief-writer and any future auditor share a model
lineage; a failure mode common to that lineage passes through all of them. **Nothing here is
independently-corroborated.** The `(2.0,1.1)` defect is offered specifically as something an
auditor can check in one line of exact arithmetic without trusting any of this page's code:
`-11*2 + 20*1.1 = 0`.

**WHAT S4 DELIVERS, AND WHAT IT DOES NOT.**

**Delivers:** the formation functional for an arbitrary carrier, derived; the exact modulus
identity and the rank-2 Fourier proof; the identification of `lambda_B` as a functional of
the relation lattice, with exact closed forms in all three tiers; the map of `lambda` over
both the connection torus and the ready-state simplex, with extrema, level-set structure,
singular loci and non-formation loci; a support-conditional verdict on F-A; ten verified CW
carriers with computed homology; the separation theorem; four controls that separate topology
from incidence; a negative and mechanised answer on `2*pi*chi`; and one checkable defect of
record.

**Does not deliver, and does not claim to:** any resolution of which schedule is physical;
any principled transport of a ready state between carriers; any sweep of the non-designated
invariants; any statement about `S5`'s free dimensionless coupling; any physical
interpretation of `lambda`; and no independent corroboration of anything above.

**Built under `CUSTODY_V001.md`. Sealed on creation.**
