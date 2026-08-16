# S2 — THE FORMATION CONDITION ON K1's OWN DATA — V001 — 2026-08-16

Stage S2 of `FOUNDING_DESIGN_V001.md`
sha256 `ca25b79c76531d909e75fcb58163ad6456b10f086d59638f47224210d19b13e9` (§7, lines 93–98).
Carrier: `S1_CARRIER_K1_V001.md`
sha256 `3eb70375bfd0900e4dd56cae294fa31b3b6e19cf6634853501fab5ffcebd92ac`.
Custody: `CUSTODY_V001.md`
sha256 `6629ae3b5aedbafda4e56b2a743b84a888949c59f988296056a2a7abeb20aa49`.

**This is a construction attempt, not an archaeology.** Everything below is exhibited on this
page. Where the answer is an impossibility, it is proved, not hedged (custody §6, `:70–74`).

---

## 0. WHAT THIS STAGE RETURNS, STATED FIRST

The stage was asked whether a formation condition of the predecessor's character — *two
source-conditioned evolutions from a common ready state become orthogonal, at a time solving an
equation* — can be written on K1's own data.

**It cannot, and the obstruction is threefold and proved.** Three independent absences, each
fatal on its own:

```
O1  RANK   rank-one fibres => transport is a scalar => rays are fixed => never orthogonal
O2  TIME   K1 carries no time parameter => "at a time solving an equation" is not expressible
O3  SPLIT  dim Gamma(L) = 5, and 5 is PRIME => no tensor factorisation exists at all
```

**And one thing was nevertheless built.** A genuinely connection-dependent, gauge-invariant
orthogonality does exist on K1's own data — not between two evolved states, but between the two
**loop-parallel sectors**, and it fires exactly when both invariants are non-trivial
(Construction 1, §3.4). It has the correct trivial-connection limit, which discharges S2's named
external contact point (§6).

**O3 is the sharpest result on this page.** The split the predecessor introduced by imperative is
not merely absent from K1. On K1's state space it is **impossible** — and therefore Q5 is answered
in the negative not by my restraint but by arithmetic.

---

## 1. Q1 — THE STATE SPACE, DERIVED

### 1.1 The object

K1 places `L_v = C` at each of the five vertices, with `<z,w> = conj(z) w`, rank one
(S1 `:48–49`). A **state on K1** is a section of that line bundle: an assignment of one vector to
each fibre.

```
Gamma(L)  =  L_v0 (+) L_v1 (+) L_v2 (+) L_v3 (+) L_v4   =   C^5
```

**This is a DIRECT SUM, not a tensor product.** The distinction is not decorative; it is the
whole of §5. Nothing here is a subsystem of anything.

### 1.2 The inner product

```
<s, t>  =  sum over v of  conj(s(v)) t(v)
```

the direct sum of the five fibre inner products. K1 supplies the fibre inner products (S1 `:48`)
and the vertex set; the sum with **unit weights** is the counting measure on vertices. K1 carries
no volume form, so no other weighting is derivable from it (CHOICE LEDGER C1).

### 1.3 What gauge does to it

Gauge is `g_v in U(1)` at each vertex (S1 `:59–60`), acting **fibrewise**:

```
(g . s)(v)  =  g_v s(v)
```

This is unitary for the inner product of §1.2, so it preserves norms and orthogonality — which is
why every orthogonality verdict on this page is gauge-invariant, and this is checked numerically
at §3.4.

### 1.4 The dimension count, displayed

Two counts, and they differ from the connection's count for a reason worth stating:

```
STATES (vectors)
  dim_R Gamma(L)                   = 10        (5 complex fibres)
  gauge group U(1)^5, dim_R        =  5
  global phase acts NON-trivially on sections (it rescales s by a common phase)
  => gauge acts through 5 effective parameters, freely where all s(v) != 0
  10 - 5                           =  5 real invariants
  and exactly five exist:  ( |s(v0)|, |s(v1)|, |s(v2)|, |s(v3)|, |s(v4)| )

STATES (rays — the physical states)
  dim_R P(Gamma(L)) = dim_R CP^4   =  8
  gauge acts through U(1)^5 / U(1)_global = 4 effective parameters
  8 - 4                            =  4 real invariants
  and exactly four exist:  the probability vector p_v = |s(v)|^2 / ||s||^2,
  a point of the 4-simplex  { p_v >= 0, sum p_v = 1 }
```

Both counts are saturated: the exhibited invariants are complete, because the moduli determine
the gauge orbit (given equal moduli, choose `g_v` to match the phases vertex by vertex).

**Contrast with S1, and the place one could slip.** On the *connection*, the global phase acts
trivially, giving `6 - 4 = 2` (S1 `:66–71`). On *sections* it does not. The two counts use
different effective gauge dimensions — 5 for vectors, 4 for rays, 4 for the connection — and
conflating them would silently produce the wrong invariant count.

**Derived result:** *the complete gauge-invariant content of a normalised state on K1 is a
probability distribution over its five vertices, and nothing else.* All phase information in a
state is pure gauge. Only the connection's phases (§3) are physical.

---

## 2. Q2 — THE CANDIDATE EVOLUTIONS, DERIVED

### 2.1 The two closed transports

Transport along `e: u -> v` is `z |-> U_e z`; reverse traversal by `U_e^{-1}` (S1 `:51–54`).
From the root `r = v0` (S1 `:14`), following the orientations of S1 `:19–21`:

```
FILLED triangle,  e1 . e2 . e3 :
  z at v0  ->  U_1 z at v1  ->  U_2 U_1 z at v2  ->  U_3 U_2 U_1 z at v0
  T_F z = exp(i(a1+a2+a3)) z = W_F z

UNFILLED triangle,  e4 . e5 . e6 :
  z at v0  ->  U_4 z at v3  ->  U_5 U_4 z at v4  ->  U_6 U_5 U_4 z at v0
  T_C z = exp(i(a4+a5+a6)) z = W_C z
```

`T_F` and `T_C` are endomorphisms of the **same** one-dimensional space `L_v0`, and their
coefficients are precisely K1's two invariants (S1 `:75–78`). Re-derived here at my own hand: the
gauge action on the six edge parameters has rank 4, its invariant subspace is 2-dimensional, and
`(a1+a2+a3)` and `(a4+a5+a6)` are annihilated by it and independent.

**Are these the natural candidates for the gate's "two evolutions"? Yes** — and they are the only
pair the carrier privileges. Reason: the based-loop group of the 1-skeleton is free of rank

```
E - V + 1  =  6 - 5 + 1  =  2       (confirmed: rank of the boundary matrix d1 is 4, nullity 2)
```

on exactly these two generators, and since U(1) is abelian the holonomy representation factors
through `Z^2 -> U(1)`, sending the two generators to `W_F` and `W_C`. So **every** closed
transport at `v0` is `W_F^m W_C^n` for integers `m, n`. The two triangles are not a selection
among many; they generate everything.

A structural note K1 forces: `pi_1(K1, v0) = Z`, not `Z^2` — the face kills the filled
triangle's class (S1 `:40–41`, `b1 = 1`). The loop `e1.e2.e3` is null-homotopic in K1 yet has
holonomy `W_F != 1` in general. That is exactly what "curvature" means here, and it is why the
two candidates are genuinely different in kind: one is a curvature, one is flat.

### 2.2 The other candidate the carrier supplies, and why I did not take it

**The covariant graph Laplacian.** K1's edges and connection assemble, with no further input,
into a Hermitian operator on the *whole* state space:

```
H_{v,u} = U_e   for e : u -> v ,     H_{u,v} = conj(U_e) ,      D = diag(deg v)
Delta_A = D - H          degrees (v0..v4) = (4, 2, 2, 2, 2),  sum = 12 = 2E
```

Hermiticity verified. Gauge acts by conjugation `Delta_A -> G Delta_A G^{-1}`, so its spectrum is
gauge-invariant and depends only on `(W_F, W_C)`.

**I did not take it as the gate's evolution, and the reason is exact: it is not an evolution
until a time is supplied, and K1 has no time.** `Delta_A` is an operator, not a dynamics;
`exp(-i t Delta_A)` requires `t in R`, which is not among K1's cells, fibres, or phases. Taking
it here would have smuggled O2 in unannounced. It is instead priced honestly in Q4 (§4.2), where
it turns out to be the *cheapest* addition — and to fail for a different, already-proved reason.

**Not taken, with reasons:** open transports `L_v0 -> L_v1` etc. (they land in different fibres,
so no inner product between them exists without an added identification); higher charge
`U_e |-> exp(i q a_e)` (still rank one, so Theorem 1 applies verbatim — this is the S4 knob and
it does not touch S2).

---

## 3. Q3 — CAN AN ORTHOGONALITY CONDITION BE WRITTEN?

### 3.1 THEOREM 1 (Ray rigidity of a rank-one connection) — the obstruction O1

> **Let `gamma` be any edge path in K1 from `u` to `v`, traversing edges in either direction.
> Then parallel transport `P(gamma) : L_u -> L_v` is multiplication by a single element of
> `U(1)`. Consequently:**
>
> **(a)** for any two paths `gamma, gamma'` from `v0` to `v0` and any `z in L_v0`,
>   `|<P(gamma) z, P(gamma') z>| = |z|^2` ;
> **(b)** in particular `|<T_F z, T_C z>| = |<W_F z, W_C z>| = |z|^2`, which is non-zero for
>   every `z != 0` and **every** value of the connection ;
> **(c)** more deeply, no non-zero `z in L_v` has *any* non-zero orthogonal partner in `L_v`.

**Proof.** Each forward edge transport is multiplication by `U_e = exp(i a_e) in U(1)`; each
reverse traversal by `U_e^{-1} in U(1)`. `U(1)` is closed under multiplication and inversion, so
the composite is `P(gamma) = exp(i phi)` for a single real `phi` — the signed sum of the `a_e`
along `gamma`. Then

```
<e^{i phi} z, e^{i phi'} z>  =  conj(e^{i phi} z) (e^{i phi'} z)  =  e^{i(phi' - phi)} |z|^2
|<...>| = |z|^2 .
```

For (c): `L_v = C`, so `z^perp = { w in C : conj(z) w = 0 } = {0}` whenever `z != 0`. ∎

**Verified rather than asserted.** 200,000 independent random draws of the six phases
`a_1..a_6` and of `z in C`:

```
max over 200,000 samples of   | |<T_F z, T_C z>| - |z|^2 |   =   1.07e-14
```

— machine zero. The refutation attempt failed. Ray check on a random connection:
`(T_F z)/z = -0.97927 - 0.20256i` and `(T_C z)/z = -0.27950 - 0.96014i`, both of modulus 1.

**The precise obstruction, stated exactly.** `dim_C L_v = 1`. Transport therefore lands in the
scalars `C^x`, and the scalars act trivially on rays. **The ray `[z]` is a fixed point of the
entire holonomy group of K1.** A gate that fires when two evolutions of a common ready state
become orthogonal requires those evolutions to *move the ray* and to move it into a complement
that exists. On K1 neither holds: nothing moves the ray, and no complement exists. This is not a
statement about the two triangles, about `b1 = 1`, or about the face — it is a statement about
rank, and it survives every choice of connection, charge, root, and path.

### 3.2 The readings that "fire" but are disqualified, and why

Custody demands the distinction be made rather than blurred, so each candidate reading is tested
against one criterion I state and own:

> **CRITERION (mine, declared — CHOICE LEDGER C4).** A formation condition must be a
> **non-constant function of the connection's gauge invariants** `(W_F, W_C)`. A condition that
> returns the same verdict at the trivial connection as at a generic one is not a formation
> condition — it is a fact about the carrier's shape, true before any connection exists.

This criterion is *entailed* by S2's named external contact point, "the trivial-connection limit
must give the known trivial answer" (`FOUNDING_DESIGN_V001.md:117–118`).

**(i) The pair `(W_F z, W_C z)` in `C^2`.** Explicitly named in the brief, so explicitly tested.
Form `v_F = (W_F z, 0)` and `v_C = (0, W_C z)`. Then `<v_F, v_C> = 0`. Computed:

```
TRIVIAL connection (all a_e = 0) :  <v_F, v_C> = 0
generic random connection        :  <v_F, v_C> = 0
```

**Disqualified twice over.** First, the orthogonality comes entirely from the direct-sum index —
from having *written down two slots* — and not from the connection; it fires identically at the
trivial connection, failing the criterion. Second, the object it lives on is **not K1's**: `C^2`
here is `L_v0 (+) L_v0`, a duplication of one fibre indexed by the two loop classes. K1 supplies
the index set; the direct sum over it is supplied by me. This reading is bookkeeping wearing the
costume of a result.

**(ii) Orthogonality of the two loops as chains.** In `C_1(K1;R)` with the cells orthonormal,
`z_F = e1+e2+e3` and `z_C = e4+e5+e6` have disjoint support, so `<z_F, z_C> = 0`. This *is* on
K1's own data and it *is* true — it is the statement that the two triangles meet only at `v0`.
But it holds for every connection and indeed before any connection exists. **Disqualified by the
criterion.** It is geometry, not formation.

**(iii) Orthogonality of histories.** A "history" on K1 is an edge path; a decoherence-functional
reading would ask `D(gamma, gamma') = conj(W(gamma)) W(gamma') x (weights)` to vanish off the
diagonal. Without weights, `|D| = 1` always, by Theorem 1(a). Vanishing requires a **measure on
paths** — amplitudes attached to the two loop classes. **K1 supplies no path measure**: it has
cells, fibres, and phases, and no weights. Added, not carrier-supplied.

**(iv) Reduced states.** Requires a partial trace, which requires a tensor factorisation.
**Impossible on K1 — see Theorem 2, §3.5.** Not "not supplied": impossible.

**(v) Maximal discrepancy — the closest K1 genuinely comes.** The one connection-dependent
gauge-invariant scalar available at `v0` is the relative holonomy `W_F^{-1} W_C`, and

```
|| T_F z - T_C z ||^2  =  |W_F - W_C|^2 |z|^2  =  4 sin^2( (f - c)/2 ) |z|^2 ,
      f = a1+a2+a3 ,  c = a4+a5+a6
```

maximal when `f - c = pi`, i.e. `W_F^{-1} W_C = -1`. This is a genuine condition, on K1's own
data, connection-dependent, gauge-invariant, and trivial at the trivial connection. **It is not
orthogonality and it is not offered as orthogonality.** It is antipodality on a circle. In rank
one the state space of a fibre is a single ray and "as far apart as possible" is a relative phase
of `pi`, whose overlap is still `|z|^2`. The gap between antipodal and orthogonal *is* the gap
between `U(1)` and rank two, and naming it as orthogonality would be exactly the error this stage
exists to avoid.

### 3.3 The second absence, independent of rank — O2

Even granting orthogonality, the gate's clause *"at a time solving an equation rather than
chosen"* requires a time variable. K1's transport is **path-indexed, not time-indexed**: `P(gamma)`
depends on which edges are traversed, not on how long anything took. There is no `t`, no metric,
no ordering beyond edge orientation, and therefore no equation for a `t` to solve. **O2 is
independent of O1:** it would survive unchanged if the fibres were rank 100.

### 3.4 CONSTRUCTION 1 — what *does* fire on K1's own data

Theorem 1 forbids orthogonality *within a fibre*. It says nothing about `Gamma(L) = C^5`, where
orthogonality is abundant. So the honest question is whether the connection can carve orthogonal
subspaces out of `Gamma(L)`. **It can.**

**Definition (on K1's own data).** For a closed edge path `gamma` at `v0`, let

```
P(gamma)  =  { s in Gamma(L) :  s(v) = U_e s(u)  for every edge e : u -> v of gamma }
```

— the sections **parallel along `gamma`**. This uses only fibres, edges, and the connection.

**LEMMA.** Let `W(gamma)` be the holonomy and `V(gamma)` the vertex set of `gamma`. Then

```
W(gamma) != 1  =>  every s in P(gamma) vanishes on V(gamma);  dim P(gamma) = 5 - |V(gamma)|
W(gamma)  = 1  =>  dim P(gamma) = 5 - |V(gamma)| + 1
```

**Proof.** Going once around, `s(v0) = W(gamma) s(v0)`, so `(1 - W) s(v0) = 0`. If `W != 1` then
`s(v0) = 0`, and each subsequent vertex on `gamma` carries a unit-modulus multiple of `s(v0)`,
hence also `0`; vertices off `gamma` are unconstrained. If `W = 1`, `s(v0)` is free and determines
`s` along `gamma`; vertices off `gamma` remain free. ∎

For K1: `V(partial F) = {v0,v1,v2}` and `V(c) = {v0,v3,v4}`, each of size 3, meeting exactly in
the root, and covering all five vertices.

> **THEOREM 3 (the two loop-parallel sectors).**
> ```
> P(partial F)  PERPENDICULAR TO  P(c)      <=>      W_F != 1  AND  W_C != 1
> ```

**Proof.**
*(<=)* Both holonomies non-trivial: by the Lemma, `P(partial F) = span{ delta_v3, delta_v4 }` and
`P(c) = span{ delta_v1, delta_v2 }`. Disjoint supports, so orthogonal.
*(=>) by contraposition.* If `W_F = W_C = 1`, both have dimension 3 and `3 + 3 = 6 > 5`, so they
intersect non-trivially and cannot be orthogonal. If exactly one is trivial, say `W_F = 1` and
`W_C != 1`, then `P(c) = span{delta_v1, delta_v2}` while `s in P(partial F)` with `s(v0) = z != 0`
has `s(v1) = U_1 z != 0`, giving `<delta_v1, s> = U_1 z != 0`. ∎

**Exhibited, all four cases** (the first row is S1's own worked instance, S1 `:114–115`):

```
   case               W_F                W_C          dim P_F  dim P_C  max|<P_F,P_C>|  ORTHOGONAL
W_F!=1, W_C!=1   -1.0000+0.0000i   -0.0000-1.0000i       2        2        0.000000       True
W_F =1, W_C!=1   +1.0000+0.0000i   -0.0000-1.0000i       3        2        0.803456       False
W_F!=1, W_C =1   -1.0000+0.0000i   +1.0000-0.0000i       2        3        0.430039       False
W_F =1, W_C =1   +1.0000+0.0000i   +1.0000-0.0000i       3        3        0.607040       False
```

Cross-checked against the theorem on 20,000 further random connections (with `W_F = 1` and
`W_C = 1` each forced in half the draws): **0 mismatches**.

**Gauge invariance, checked not assumed.** Three random gauge transformations
`a_e -> a_e + theta_target - theta_source` applied to the generic case: dimensions `(2,2)`
unchanged, `max|<P_F,P_C>| = 0.0` exactly, `(W_F, W_C)` unchanged.

**Supports, and a detail K1 chose rather than I:**

```
P_F row-norms over (v0..v4) = [0, 0, 0, 1, 1]      P_F lives on {v3, v4}
P_C row-norms over (v0..v4) = [0, 1, 1, 0, 0]      P_C lives on {v1, v2}
```

Each sector lives on the *other* triangle's free vertices, and **the root `v0` lies in neither**:
`2 + 2 = 4 < 5`, with `span{delta_v0}` the orthogonal complement of both. The join vertex, the one
place the two triangles meet, is excluded from both sectors by the mechanism itself.

**What Construction 1 is, stated precisely and without inflation.** It is a genuine,
gauge-invariant, connection-dependent orthogonality of two subspaces of K1's own state space,
determined by K1's own two loops, with the correct trivial-connection limit. **It is not the
predecessor's gate.** There is no time, no evolution, and no common ready state; and the
mechanism is *annihilation*, not separation — non-trivial holonomy destroys sections on its own
loop, and the two survivors are orthogonal because what remains of them is disjointly supported.
Orthogonality here arises from destruction, not from dynamics. That is a real result about K1 and
it is offered as exactly that.

### 3.5 THEOREM 2 (No factorisation) — the obstruction O3, and the sharpest result here

> **`Gamma(L) = C^5` admits NO non-trivial tensor factorisation. Hence on K1 there is no
> subsystem, no environment, no partial trace, and no reduced state — not as an omission, but as
> an impossibility.**

**Proof.** A Hilbert space of dimension `n` factorises as `H_A (x) H_B` only with
`dim H_A . dim H_B = n`. Enumerating the divisors of 5:

```
factorisations 5 = a . b :   (1,5), (5,1)
non-trivial (a >= 2 and b >= 2) :   NONE.       5 is prime.
```

A three-way split with all factors non-trivial requires `dim >= 2.2.2 = 8 > 5`, so the
predecessor's **source / designated-record / environment** factorisation is excluded *a fortiori*.
∎

**Why this is the sharpest thing on the page.** The brief asks whether a formation condition can
be written "without importing a split". Theorem 2 upgrades the answer from *I did not need one*
to *one cannot be had*. The predecessor's Hilbert space is written by an imperative — *"Let the
source, designated record subsystem, and required environment form one closed Hilbert space"* —
and on K1 that imperative has no legal object. `Gamma(L)` is a direct sum over vertices; direct
sums do not admit partial traces.

**And the honest limit of the theorem, stated here rather than left to be found.** The primality
of 5 is a property of *this* carrier, not a structural law. A complex with 6 vertices would give
`C^6 = C^2 (x) C^3`, and a factorisation would become available — non-canonically, but available.
So the durable half of O3 is the weaker, structural half: **K1 supplies no canonical
factorisation, because `Gamma(L)` is a direct sum indexed by vertices and nothing in K1 groups
those vertices into subsystems.** The primality makes it, for this carrier, not merely
non-canonical but non-existent. Both halves are stated; only the first generalises. Flagged at
FLAG BLOCK F2.

---

## 4. Q4 — THE MINIMUM ADDITION, AND WHAT IT COSTS

Nothing below is added to K1. Each is described and priced.

### 4.1 What is minimally required, stated as a requirement before it is shopped for

A condition of the gate's character needs operators `T, T'` on a common space `H` and a
`z != 0` in `H` with `<T z, T' z> = 0`. By Theorem 1 this needs **both**:

```
(R1)  dim H >= 2                       — otherwise no complement exists
(R2)  T, T' not both scalars           — otherwise the ray is fixed
```

Rank one fails R1 and R2 simultaneously. Note that **rank 2 alone is not enough**: a rank-2
bundle with a scalar `U(1) x I` connection still gives `|<T_F z, T_C z>| = 1.0` — computed. The
connection must be non-abelian in its action, not merely the fibre wider.

### 4.2 Addition (i) — a TIME PARAMETER only. The cheapest, and it fails elsewhere

**What is added:** one real number `t`, plus a choice of Hamiltonian. The operator itself is
**K1's own** (§2.2): `Delta_A = D - H`, built from edges and connection with nothing imported.
Only `t` comes from outside — K1 has no time (O2).

**Does it write the condition? Partly, and the details are unexpected.** Take the ready state to
be a vertex delta and ask for `|<psi(0), exp(-i t Delta_A) psi(0)>| = 0`.

**From the root `v0`, it can NEVER fire — and this is exact, not numerical.** For
`A(t) = sum_j p_j exp(-i lambda_j t)` with spectral weights `p_j` summing to 1, the triangle
inequality gives `|A(t)| >= 2 p_max - 1`. Scanning the full invariant torus `(f, c)` on a 241x241
grid:

```
min over (f,c) of  p_max  =  0.621268   at (f,c) = (pi, pi)
=>  |<psi(0), psi(t)>|  >=  2(0.621268) - 1  =  0.242536  >  0   for ALL t and ALL connections
```

**The root ready state can never become orthogonal to itself.** Its degree is 4 of the total 12
(§2.2), which concentrates spectral weight above one half and closes the polygon inequality. K1's
own designated origin (S1 `:14`) is precisely the wrong place to start a gate.

**From a non-root vertex it does fire.** With ready state `delta_v1` (degree 2), `min p_max`
falls to `0.25`, the bound goes vacuous, and joint search over `(f, c, t)` — three unknowns
against two real equations, so solutions are expected on a curve — reaches

```
|<psi(0), psi(t)>| = 1.045e-06   at  f = 3.663319822, c = 2.194049746, t = 51.456642361
spectrum = [0.695897, 1.0, 1.802623, 3.0, 5.501479]
weights  = [0.117930, 0.378913, 0.177096, 0.229476, 0.096584]   p_max = 0.3789 <= 1/2
W_F = -0.866960 - 0.498378i        W_C = -0.583680 + 0.811984i
```

**Cost, itemised.**
- One real parameter `t`, from **outside** K1.
- A choice of Hamiltonian: `D - H` versus `H` versus a weighted variant. K1 does not select one
  (CHOICE LEDGER C3, **OPEN**).
- **A split? No.** Orthogonality in `C^5` is between two *vectors*; no factorisation is used or
  needed. This route writes a formation-type condition with **no split whatsoever**.
- **But it walks straight into the proved obstruction.** Measured on the same run:
  ```
  sup |<psi(0), psi(t')>|  for t' > t* + 5   =   0.994373
  ```
  The state returns to 99.4% of itself. Five eigenvalues means an almost-periodic return
  amplitude. This is exactly *"a finite discrete spectrum is recurrent"* and *"the one-cell
  operator therefore provides a reversible write, not a durable record"*
  (`FOUNDING_DESIGN_V001.md:54–55`, quoting the predecessor; **grade: adopted-as-premise-of-the-
  design, not independently re-derived from the predecessor corpus** — see IMPORT AUDIT I3).
  Orthogonality fires and then un-fires. It is a write, not a record.

**Verdict on (i): the cheapest addition, and it reproduces the predecessor's failure on K1 rather
than escaping it.** That is informative: the failure is not an artefact of the predecessor's
split. It is the finiteness, and it is waiting for S3.

### 4.3 Addition (ii) — a SECOND FIBRE DIRECTION. Rank 2 with a non-abelian structure group

**What is added:** `L_v = C^2` and structure group `SU(2)` or `U(2)`. **Entirely from outside** —
S1 fixes rank one (S1 `:48–49`).

**It writes the condition immediately.** Exhibited:

```
T_F = I,   T_C = [[0, -1], [1, 0]] in SU(2),   ready state z = (1,0)
T_F z = (1, 0)        T_C z = (0, 1)        <T_F z, T_C z> = 0        ORTHOGONAL
```

R1 and R2 are both met and Theorem 1 no longer applies.

**Cost, itemised — and the cost is not where one would look for it.**
- Rank is a stipulation, not a derivation. K1 gives no second direction at a vertex: there is no
  spare index anywhere in five vertices, six phases, one face.
- The invariant count is destroyed and must be rebuilt: `6 x dim U(2) = 24` edge parameters
  against `5 x 4 = 20` gauge parameters; holonomies no longer commute, and the invariants become
  conjugacy-class data (traces of words in `T_F, T_C`), not two angles. S1's clean `6 - 4 = 2`
  does not survive.
- **A split? Yes — and this is the non-obvious cost.** `Gamma(L)` becomes `C^10`, and
  ```
  factorisations of 10 :  (1,10), (2,5), (5,2), (10,1)     — 10 is COMPOSITE
  ```
  `C^10 = C^2 (x) C^5` = (internal) (x) (positional). **Raising the rank to 2 hands you a
  bipartition for free**, and Theorem 2's protection is gone the moment the fibre widens.
- **The precise honest qualification.** That factorisation is **not canonical**: identifying each
  `L_v` with a fixed `C^2` requires a trivialisation, i.e. a gauge choice, and a gauge
  transformation `g_v in U(2)` changes the factorisation. A gauge-invariant trivialisation exists
  only when the connection is flat — which is the case the whole project is not interested in. So
  rank 2 supplies a bipartition that is *available but gauge-dependent*. It is **not** a
  source/record/environment split: it is internal-versus-positional, it has two factors and not
  three, and nothing in it designates a record.

**Verdict on (ii): the minimum addition that makes the gate's condition literally writable is
rank 2 with a non-scalar connection.** It does not drag in a source/record/environment split. It
does drag in a gauge-dependent internal/positional bipartition, which is a weaker but real cost,
and it should be declared under custody §5 if S3 takes this route.

### 4.4 Ranked, since the brief asks for the minimum

```
CHEAPEST  (i)  add t only            — no split at all; fails durability by recurrence (proved)
NEXT      (ii) add rank 2 + SU(2)    — writes the gate exactly; costs a gauge-dependent bipartition
NOT A FIX      more vertices         — rank-one no-go is rank-based, not size-based; Theorem 1 is
                                       untouched by enlarging the complex
NOT A FIX      higher charge q       — still rank one; Theorem 1 applies verbatim
```

**The minimum that answers the brief's question is (i): one real parameter.** The minimum that
makes the *predecessor's gate as written* expressible is (i) **and** (ii) together — time for the
"when", rank for the "orthogonal".

---

## 5. Q5 — THE SPLIT TEST

**Plainly: NO. Nothing in this attempt required a source/record/environment factorisation, a
preferred subsystem, a system/environment tensor split, or a chosen decomposition of that kind.**

And the answer is stronger than compliance. By **Theorem 2** (§3.5) no such split *exists* on
`Gamma(L) = C^5`. I did not decline to import one; there was no object to import it onto.

**Full disclosure of every decomposition and every preference actually used, with provenance:**

| Used | What it is | Whose? |
|---|---|---|
| `Gamma(L) = (+)_v L_v` | the vertex-indexed **direct sum** | **K1's** — S1 `:48` places one fibre per vertex; this is the definition of a section, not a split. No partial trace is definable over it. |
| the root `v0` | base point for closed transports | **K1's** — S1 `:14`, "Root: r = v0" |
| the two loops `partial F`, `c` | the two sector labels in Construction 1 | **K1's** — S1 `:24–25`; and by §2.1 they *generate* all closed transports, so they are not a selection |
| filled vs unfilled | which loop is a curvature | **K1's** — S1 `:24–25`, `:40–41` |
| unit weights in `<s,t>` | counting measure on vertices | **mine**, declared — C1. K1 carries no volume form. Every result on this page is invariant under reweighting (the supports in §3.4 are disjoint regardless of positive weights). |
| the criterion of §3.2 | "must depend on `(W_F, W_C)`" | **mine**, declared — C4, and entailed by `FOUNDING_DESIGN_V001.md:117–118` |

No subsystem was designated. No environment was posited. No state was traced out. **This is, so
far as this page can establish, the first construction in this lineage to state a
formation-type condition without a stipulated split** — Construction 1 (§3.4) is stated entirely
on `(W_F, W_C)` and K1's loops, and addition (i) (§4.2) writes a firing orthogonality condition
with no factorisation anywhere.

The claim is bounded exactly as follows: *this page* uses no split, and *this carrier* admits
none. Whether a crossing at S3 can be built without one is not settled here, and §4.3 shows the
most natural strengthening already costs a bipartition.

---

## 6. THE EXTERNAL CONTACT POINT, DISCHARGED

`FOUNDING_DESIGN_V001.md:117–118` names S2's contact as *"the trivial-connection limit must give
the known trivial answer."*

**Discharged, and it is a real test rather than a formality — two candidate conditions were
disqualified by it** (§3.2 (i) and (ii), both of which fire at the trivial connection and were
rejected for exactly that). The surviving objects behave correctly:

```
TRIVIAL CONNECTION  (all a_e = 0, so W_F = W_C = 1) :
  Construction 1 :  dim P_F = dim P_C = 3 in C^5;  3+3 > 5;  NOT orthogonal.   No formation. ✓
  Discrepancy    :  |W_F - W_C| = 0.  No separation of any kind.                            ✓
  Addition (i)   :  Delta_A becomes the ordinary graph Laplacian; the root bound
                    |A(t)| >= 0.2425 still holds; nothing fires from the root.               ✓
  Theorem 1      :  |<T_F z, T_C z>| = |z|^2, unchanged — as it must be, being rank-based.   ✓
```

**No formation at trivial connection is the known trivial answer, and Construction 1 returns
exactly it.**

---

## CHOICE LEDGER

| # | Choice | Alternatives | Why | Status |
|---|---|---|---|---|
| C1 | Inner product on `Gamma(L)` = direct sum with **unit** weights | any positive vertex weights `w_v`; degree weighting | K1 carries no volume form (S1 has no metric data), so counting measure is the only weighting derivable without addition | **closed** — every result is invariant under positive reweighting; §3.4's supports are disjoint for any `w_v > 0` |
| C2 | States = sections of the line bundle, `(+)_v L_v` | states at the root only (`L_v0`, dim 1); states on edges; states on the face | the fibres are placed at vertices by S1 `:48`; a section is the unique object using all of them | **closed** — root-only states are the dim-1 special case, covered by Theorem 1(c) |
| C3 | Hamiltonian for addition (i) = `Delta_A = D - H` | magnetic adjacency `H`; normalised Laplacian `D^{-1/2} H D^{-1/2}`; any positive function of these | `D - H` is the discrete `d*d`, matching S1's `f = da` idiom (S1 `:89–92`) | **OPEN** — K1 does not select one. The recurrence verdict is robust across all of them (any finite Hermitian operator on `C^5` has 5 eigenvalues); the *value* of `t*` is not. Flagged for S3. |
| C4 | Criterion: a formation condition must depend non-constantly on `(W_F, W_C)` | accept any vanishing inner product as formation | otherwise readings (i) and (ii) of §3.2 "succeed" while carrying zero information about the connection | **closed** — entailed by `FOUNDING_DESIGN_V001.md:117–118` |
| C5 | Ready state for addition (i) taken as a vertex delta, then arbitrary | superposition ready states only | vertex deltas are the states K1's own basis supplies; arbitrary states were also searched (§4.2) | **closed** — both reported; the root bound of §4.2 holds for the root delta specifically and is stated as such |
| C6 | Construction 1 uses **loop-parallelism**, not transport, to build sectors | transported states (blocked by Theorem 1); eigenspaces of `Delta_A` (needs C3) | it is the only connection-dependent subspace construction using no addition whatsoever | **closed** |

---

## IMPORT AUDIT

Every notion used that is not defined on this page, with whether the finding survives without it.

| # | Notion | Source | Defined here? | Survives without it? |
|---|---|---|---|---|
| I1 | K1: cells, incidence, orientation, root, fibres, connection, gauge, `W_F`, `W_C` | `S1_CARRIER_K1_V001.md` sha256 `3eb70375…d92ac`, `:13–25`, `:48–54`, `:59–64`, `:75–78` | restated in full at §1–§2 | n/a — this is the subject |
| I2 | The predecessor's gate: "two source-conditioned evolutions from a common ready state become orthogonal, at a time solving an equation"; and the imperative *"Let the source, designated record subsystem, and required environment form one closed Hilbert space"* | the brief for this stage, restating predecessor material | quoted, **not** used as premise | **YES.** It is the *target being tested*, never a foundation. Theorems 1–3 are proved from K1 alone and do not cite it. **Carries no transfer grade** — per `INHERITED_FROM_THE_PREDECESSOR.md:11–16` (sha256 `96c1d305…d0079`) no predecessor asset is premise until it appears in the salvage assessment with a grade. Used here only as a description of what to test for. |
| I3 | "a finite discrete spectrum is recurrent"; "the one-cell operator therefore provides a reversible write, not a durable record" | `FOUNDING_DESIGN_V001.md:54–55` (sha256 `ca25b79c…13e9`), itself quoting the predecessor | no | **YES.** §4.2 does not rely on the quotation: recurrence is **exhibited numerically on K1** (`sup|overlap| = 0.994373`) and follows from `dim = 5` finite. The quotation names the phenomenon; it does not establish it here. Grade: adopted-as-design-premise, not independently re-derived from the predecessor corpus. |
| I4 | "the trivial-connection limit must give the known trivial answer" | `FOUNDING_DESIGN_V001.md:117–118` | no | **PARTLY.** Theorems 1, 2, 3 stand without it. Criterion C4 and the disqualifications of §3.2 (i)–(ii) are motivated by it; without it those readings would be merely uninformative rather than disqualified. |
| I5 | Alpha | not used | — | **YES.** No coupling, no numeric target, and no measured number enters this page. Custody §5 (`:55–68`) is not engaged: nothing here was selected by where it lands. |
| I6 | Standard linear algebra: `U(1)`, `SU(2)`, Hermitian spectral theorem, `dim(A(x)B) = dim A . dim B`, triangle inequality, rank–nullity | mathematics | no | **YES** — ambient, not project-specific. |

---

## FLAG BLOCK — defects found in my own draft on self-check

**F1 — Construction 1 does not see the face, and that is a real weakness.** Theorem 3's criterion
is symmetric in `W_F` and `W_C`. It cannot tell a curvature from a flat holonomy, so K1's
distinguishing feature — one filled triangle, one unfilled, `b1 = 1` (S1 `:24–25`, `:40–41`) —
is **invisible** to it. The carrier was built to separate the two (S1 `:109`) and this
construction does not use the separation. Any S3 condition that still cannot tell `W_F` from
`W_C` is not using the carrier it was given.

**F2 — Theorem 2's strongest form is carrier-specific, not structural.** It leans on 5 being
prime. A 6-vertex complex would admit `C^6 = C^2 (x) C^3`. Stated in place at §3.5; repeated here
so it cannot be quoted away. The generalising half is the weaker one: no *canonical*
factorisation, because a direct sum over vertices designates no subsystem.

**F3 — C3 is OPEN and the numbers in §4.2 depend on it.** The specific `t*`, the value
`0.621268`, and the bound `0.242536` are computed for `Delta_A = D - H`. The *qualitative*
verdicts — root can never fire; non-root can; everything recurs — are robust (any Hermitian
operator on `C^5` has five eigenvalues, and the `2 p_max - 1` bound is Hamiltonian-independent
given the weights), but the numbers are not portable to another choice.

**F4 — §4.2's root bound is proved over a 241x241 grid, not analytically.** `min p_max = 0.621268`
is a numerical minimum over the invariant torus, not a closed-form one. The *inequality*
`|A(t)| >= 2 p_max - 1` is exact; the *value* of `min p_max` is grid-bounded. The margin
(`0.6213` versus the critical `0.5`) is wide and the function is smooth in `(f,c)`, so the
conclusion is not delicate — but it is not a proof, and is labelled accordingly here and at §4.2.

**F5 — "first in this lineage" (§5) is a claim about this page, not a survey.** I did not sweep
the predecessor's 5,512 archive files to confirm no earlier split-free formation condition exists;
that corpus is set aside (`INHERITED_FROM_THE_PREDECESSOR.md:6–9`). The defensible claim is the
bounded one already stated: this page uses no split, and this carrier admits none.

**F6 — Custody §4 pairing not satisfied in lineage.** This build was self-checked, not audited by
an independent lane. Per custody §4 (`:44–53`) the findings here are graded
**adversarially-checked** at best, never **independently-corroborated** — and in fact the
adversarial pass was performed by the same model that built. The 200,000-sample refutation of
Theorem 1 and the 20,000-sample cross-check of Theorem 3 are the strongest independence available
here, and they are machine checks, not lineage independence.

---

## 7. WHAT S2 DELIVERS, AND WHAT IT DOES NOT

**Delivers:** the state space derived with both dimension counts saturated; the two candidate
evolutions derived and shown to *generate* all closed transports; **Theorem 1**, the rank-one
no-go, proved and machine-refuted 200,000 times; **Theorem 2**, no factorisation, from the
primality of 5; **Theorem 3 / Construction 1**, a connection-dependent gauge-invariant
orthogonality on K1's own data with the correct trivial limit; the minimum addition priced two
ways, including the exact result that K1's own root can never fire; and a negative answer to the
split test that is arithmetic rather than restraint.

**Does not deliver, and does not claim to:** the predecessor's gate (it cannot be written here,
and that is the finding); any crossing (S3); any durability — §4.2 shows the cheapest addition
recurs to 99.4% and is a reversible write, which is precisely the obstruction S3 must answer
(`FOUNDING_DESIGN_V001.md:51–66`); any dependence on the face, which F1 records as the gap S3
should close first.

**S2's falsifier was "if it cannot be stated without a stipulated split, say so plainly"**
(`FOUNDING_DESIGN_V001.md:97`). It can be stated without one — Construction 1 does, and
addition (i) does — and the gate *as the predecessor wrote it* cannot be stated at all, for three
independent reasons, each proved. Both halves are said plainly.

## 8. CUSTODY

Built under `CUSTODY_V001.md` sha256 `6629ae3b5aedbafda4e56b2a743b84a888949c59f988296056a2a7abeb20aa49`.
Every claim exhibited on this page; every computation reproducible from the definitions given
here. Predecessor material cited by digest only, never copied, never as foundation (§IMPORT AUDIT
I2, I3). No git action taken. Sealed on creation.
