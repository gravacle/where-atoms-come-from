# S3 — THE CROSSING — V001 — 2026-08-16

**This is a construction.** Stage S3 of `FOUNDING_DESIGN_V001.md` — *build the crossing, against the
nine requirements, answering the recurrence obstruction.* Built under `CUSTODY_V001.md`.

Everything below is exhibited: object stated, computation run, result displayed. Every number on
this page was produced by the code described beside it and is reproducible from the definitions
given here alone.

---

## 0. THE VERDICT, STATED FIRST

**THE LIMIT ESCAPES RECURRENCE, AND THE ESCAPE IS CONDITIONAL.**

1. **No carrier can hold the record — proved in general, not for K1 only** (§2.4). Every loop
   transport on a rank-one Hermitian fibre bundle with a `U(1)` connection is a **diagonal**
   unitary; the gauge-invariant algebra of *any* such complex is exactly the diagonal; and a
   diagonal unitary leaves the diagonal of a density matrix fixed. So the two branch states are
   **indistinguishable by every gauge-invariant carrier observable, at every circuit count, on
   every complex, however large.** Growing the carrier cannot help. **This is why a crossing is
   needed at all**, and it is now a theorem rather than an obstruction inherited from elsewhere.

2. **The directed system is forced to be multiplicative, not cellular** (§3). Repeated circuits
   generate a **3-dimensional** algebra at `N = 1` and the same 3-dimensional algebra at `N = 100`
   — the trap, disarmed by computation. Cellular growth of K1 admits **no unital connecting
   homomorphism** (5 ∤ 9, 5 ∤ 7). What grows is a **record slot per cell**:
   `A_N = M_5(C) (x) M_2(C)^{(x)N}`, `A_infinity = UHF(5 * 2^infinity)`.

3. **Recurrence fails in the limit, and it fails by monotonicity rather than by hypothesis**
   (§4). The one-cell overlap still returns — to `0.999941` at circuit 377, exactly as W-01
   recorded. The **record functional** `Omega_N = prod_{n<=N} Z_n` is **strictly decreasing at
   every step**, hence injective, hence returns to no earlier value ever. It falls as
   `e^(lambda N)` with `lambda = -0.766802` on the audit's own test connection.

4. **All nine hold in the limit** (§5). **P-1 and P-9 hold in the limit and FAIL at every finite
   stage** — which is precisely why the predecessor recorded P-9 false: it is a limit property and
   was evaluated at a finite stage. P-9 additionally holds *at finite stage* exactly on W-01's
   firing locus, and does so on K1's own published connection at the **first cell**.

5. **The formation condition transports, multiplicatively, and firing is absorbing** (§6).
   `Omega_{N+1} = Omega_N * Z_{N+1}`. Once fired, fired forever.

6. **The price, stated rather than hidden** (§7, FLAG BLOCK). The escape requires the cell
   schedule to sample the carrier's clock with divergent deficiency. **An adversarial schedule
   locked to the carrier's near-recurrence times defeats the construction entirely** — exhibited
   at §4.6 with explicit circuit numbers. And the construction is **not zero-addition**: it costs
   exactly one qubit per cell, which is proved minimal (§2.5) but is not nothing. S2 was free;
   S3 is not.

---

## 1. POINTERS AND DIGESTS

Custody §1: no term in a governing clause without a digest, a `file:line`, or a definition on this
page. Digests computed in this directory at build time with `shasum -a 256`:

```
S1_CARRIER_K1_V001.md                       3eb70375bfd0900e4dd56cae294fa31b3b6e19cf6634853501fab5ffcebd92ac
S2_FORMATION_CONDITION_ON_K1_V001.md        248ce856efaef157c68e818dde589d0200bbc1dd9fd9fc1fcc8cdc7bc88734d9
S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md  0bea11bd4b7764f65c8d44cc8812d43bdcc569d34f51d0c46d82851a2efbd0d5
REGISTER_V001.md                            95ecd65d1d1b4e0440d29a69e6ec1cb36ada8059d3e88948ff2177d26f62000b
CUSTODY_V001.md                             6629ae3b5aedbafda4e56b2a743b84a888949c59f988296056a2a7abeb20aa49
FOUNDING_DESIGN_V001.md                     ca25b79c76531d909e75fcb58163ad6456b10f086d59638f47224210d19b13e9
INHERITED_FROM_THE_PREDECESSOR.md           96c1d305c7a4eeceab8dc78f971832eb0c2e34980a72d1c12c5ece5cd4ad0079
```

**Terms defined on this page, in the order they are first used:** *record algebra* (§2.1),
*gauge-invariant subalgebra* (§2.2), *loop algebra* (§2.3), *record slot* (§2.5), *directed system*
(§3.1), *unital \*-homomorphism* (§3.1), *direct limit* (§3.1), *UHF algebra* (§3.4), *write map*
(§3.3), *schedule* (§3.5), *deficiency* (§4.5), *record functional* (§4.2), *disjoint states*
(§5.7), *asymptotically central* (§5.2), *sector* (§5.7).

**Terms taken from elsewhere by digest, never copied and never load-bearing:** the nine
requirements P-1…P-9 (`INHERITED_FROM_THE_PREDECESSOR.md`, cited as a specification only, per
custody §7); W-01's formation condition and the recurrence figure (`REGISTER_V001.md`, and both are
**re-derived from scratch on this page**, §4.1). The inductive-limit template
(`FOUNDING_DESIGN_V001.md:61–65`) is cited as a lead, and what this page builds departs from it in
one respect stated at §3.2.

---

## 2. Q1 — THE RECORD ALGEBRA ON K1, DERIVED

### 2.1 The observable algebra

The state space is `Gamma(L) = (+)_{v} L_v = C^5`: five rank-one Hermitian fibres, S1 `:48–49`,
assembled into a direct sum by the S2 build's Q1 and adopted by W-01 (`REGISTER_V001.md`, "State
space `Gamma(L) = C^5`"). The algebra of observables on a finite-dimensional Hilbert space `H` is
`B(H)`, all bounded operators; here

```
B(Gamma(L))  =  M_5(C)          dim_C  =  25
```

There is no choice in this: `M_5(C)` is the full operator algebra of `C^5`.

### 2.2 The gauge-invariant subalgebra, with the counts

Gauge acts at vertices, `g_v = exp(i theta_v)` (S1 `:59–63`). On sections, `(g s)(v) = g_v s(v)`;
on observables, `X -> g X g*`. Entrywise:

```
(g X g*)_{uv}  =  exp( i (theta_u - theta_v) ) X_{uv}
```

so `X_{uv}` carries the character `theta -> exp(i(theta_u - theta_v))` of `U(1)^5`. That character
is trivial iff `u = v`. The fixed-point algebra is therefore the **diagonal**:

```
A_inv  =  D_5  =  C^5        dim_C  =  5
```

**Exhibited, exactly.** The average of `X -> g X g*` over the finite subgroups `Z_2^5` (32 elements)
and `Z_3^5` (243 elements) of `U(1)^5` — both large enough that every non-trivial character above
averages to zero — applied to a random complex `X`:

```
|| E(X) - diag(X) ||_F   over Z_2^5 :   1.0728e-15
|| E(X) - diag(X) ||_F   over Z_3^5 :   1.0568e-14
```

**The counts, displayed and saturated:**

```
dim_C M_5(C)                                              =  25
matrix units carrying the TRIVIAL character (u = v)       =   5
distinct NON-trivial characters exp(i(theta_u-theta_v))   =  20
matrix units carrying them                                =  20
                                                    25   =   5  +  20
```

Nothing is hidden: the 25 dimensions split exactly into 5 invariant and 20 non-invariant, and the
5 are the diagonal. For comparison, S1's connection count is `6 - 4 = 2` (S1 `:66–70`), also
saturated. These are different counts of different objects and neither constrains the other.

### 2.3 The loop algebra — the three-way split, DERIVED here

W-01 (`REGISTER_V001.md`) extends loop transport to `Gamma(L)` as

```
M_dF  =  diag( W_F, W_F, W_F,  1,   1  )        V(dF) = {v0,v1,v2}
M_c   =  diag( W_C,  1,   1,  W_C, W_C )        V(c)  = {v0,v3,v4}
```

Both are diagonal with gauge-invariant entries, so **both lie in `A_inv`**. Ask what they generate.
The unital \*-algebra `C*(M_dF, M_c)` is computed by taking all words up to length 8 in
`{I, M_dF, M_c, M_dF*, M_c*}` and measuring the rank of their span:

```
  connection                        W_F              W_C            dim C*(M_dF, M_c)
  generic  (f=2.0, c=1.1)      -0.4161+0.9093j   0.4536+0.8912j            3
  S1 published instance        -1.0000+0.0000j  -0.0000-1.0000j            3
  W_F = 1                       1.0000+0.0000j   0.4536+0.8912j            2
  trivial connection            1.0000+0.0000j   1.0000+0.0000j            1
```

**Generic value 3, and the three blocks are read straight off the generators:**

```
M_dF values by vertex :  [ W_F  W_F  W_F   1    1  ]
M_c  values by vertex :  [ W_C   1    1   W_C  W_C ]
   v0 -> (W_F, W_C)      v1,v2 -> (W_F, 1)      v3,v4 -> (1, W_C)
```

The pair separates the vertex set into exactly `{v0} / {v1,v2} / {v3,v4}`. So

```
C*(M_dF, M_c)  =  C.P_0  (+)  C.P_F  (+)  C.P_C  =  C^3
```

**This is W-01's canonical three-way split, obtained here as the algebra generated by K1's own two
loop transports.** W-01 derived it from incidence; it is re-derived here from the transports, by an
independent route. **A caveat W-01 does not carry:** the dimension is 3 only for `W_F != 1` and
`W_C != 1`; it drops to 2 and then 1 as the connection trivialises. The split is
**connection-dependent, not a fixed structure of the complex.** Recorded at FLAG F8.

### 2.4 THEOREM S3-0 — THE CARRIER CANNOT HOLD THE RECORD (general, not K1-specific)

> **Theorem.** Let `K` be any finite complex with a rank-one Hermitian fibre at each vertex and a
> `U(1)` connection on edges. Let `M_gamma` be the loop transport of §2.3 for any closed path
> `gamma`. Then for every state `rho` on `Gamma(L)` and every `n`,
> `diag( M_gamma^n rho M_gamma^{-n} ) = diag(rho)`. Consequently, for every **gauge-invariant**
> observable `X` and any two loops `gamma, delta`,
> `<M_gamma^n s, X M_gamma^n s> = <M_delta^n s, X M_delta^n s>`.
>
> **The two branch states are indistinguishable by every gauge-invariant carrier observable, at
> every circuit count, on every such complex, however many cells it has.**

**Proof.** `M_gamma` is diagonal, with entries of modulus 1. Hence
`(M rho M*)_{vv} = |M_{vv}|^2 rho_{vv} = rho_{vv}`. By §2.2 the gauge-invariant algebra of
`Gamma(L)` is exactly the diagonal for any vertex set, so a gauge-invariant `X` sees only
`diag(rho)`, which is loop-transport-invariant. Both branches start from the same `s`, so both give
`diag(s s*)`. ∎

**Exhibited on K1**, three random ready states and random circuit counts:

```
n = 11 :  max_v | |x(v)|^2 - |y(v)|^2 |  =  8.882e-16
n = 31 :  max_v | |x(v)|^2 - |y(v)|^2 |  =  5.829e-16
n =  8 :  max_v | |x(v)|^2 - |y(v)|^2 |  =  2.776e-16
```

**This theorem is the reason S3 exists.** It says the crossing is not optional and not a matter of
making the carrier bigger. It also disposes, in one line, of every construction whose record
algebra is the gauge-invariant algebra of a complex — including the two the brief names as
candidates (more triangles at the root; a growing sequence of complexes). §3.2 returns to them.

### 2.5 Where the record lives, and why the slot has dimension 2

The record cannot be a subalgebra of the carrier's gauge-invariant algebra (§2.4), and it must not
be a subalgebra of a *finite* algebra at all (the recurrence obstruction, reproduced at §4.1). So
the record lives in an algebra the carrier does not contain. Call its per-cell factor `R`, the
**record slot**.

**What `R` must do.** A cell's write must leave the two branches in two slot states whose overlap
is the cell's own comparison value `Z` (that value cannot be improved on — see §3.3, a unitary
cannot manufacture orthogonality). If the write is unitary and returns the carrier to its ready
state, then the slot states are **pure** states of `R`.

**Hence `R` is non-abelian.** The pure states of an abelian algebra `C^k` are its `k` characters,
which are pairwise orthogonal: two of them have overlap `0` or `1`, never a general `Z`. An abelian
slot can therefore hold no comparison value other than *fully formed* or *not formed at all*.

**Hence `dim R >= 4`, i.e. `R = M_2(C)` at the minimum**, since two unit vectors span at most two
dimensions and `M_2(C)` is the operator algebra of `C^2`. And **2 is K1's own number**: it is the
count of independent based loops on K1 — one from `b1 = 1` and one from the single face
(S1 `:33–44`) — and it is the count of gauge invariants, `6 - 4 = 2` (S1 `:66–77`). The slot has
exactly one complex dimension per loop the carrier can compare.

```
                                                          dim_C
carrier, full observable algebra        M_5(C)              25
carrier, gauge-invariant subalgebra     D_5 = C^5            5
carrier, loop algebra (the split)       C^3                  3      (connection-dependent)
record slot, one cell                   R = M_2(C)           4
record slot, classical pointer part     C^2 subset M_2       2
```

**The record lives in `R^{(x)N}` and its limit, and in no subalgebra of `M_5(C)`.** That statement
is the crossing.

---

## 3. Q2 — THE DIRECTED SYSTEM, BUILT

### 3.1 Definitions, given here

A **directed system** of C\*-algebras is a sequence `A_0, A_1, A_2, ...` with connecting maps
`iota_N : A_N -> A_{N+1}`. Each `iota_N` must be a **unital \*-homomorphism**: linear,
`iota(XY) = iota(X) iota(Y)`, `iota(X*) = iota(X)*`, and `iota(1) = 1`. The **direct limit**
`A_infinity` is the completion of the union of the images in the norm inherited along the system;
it is a C\*-algebra, and it is **not finite-dimensional** whenever `dim A_N -> infinity`.

### 3.2 What grows, and why the alternatives were not taken

**Taken: a fresh record slot per cell.**

```
A_N  =  M_5(C)  (x)  M_2(C)^{(x)N}          iota_N(X)  =  X (x) I_2
```

**Not taken — repeated circuits of one loop.** The brief's warning is correct and here is the
computation. The span of `{ M_dF^n, M_c^n : n <= N }` inside `M_5(C)`:

```
  N =    1 :  dim = 3          N =   10 :  dim = 3
  N =    2 :  dim = 3          N =   25 :  dim = 3
  N =    3 :  dim = 3          N =  100 :  dim = 3
  N =    5 :  dim = 3
```

**Constant at 3 forever.** Powers of block-constant functions are block-constant, so circuits
generate exactly the three-way split algebra `C^3` of §2.3 and nothing more. Circuits supply
**discrete time** (audit COR-F) and **no algebra**. They are used here as the clock (§3.5) and never
as the growth.

**Not taken — more unfilled triangles at the root (`b1` grows), or a chain
`K^(1) subset K^(2) subset ...` adding cells.** Two independent computations kill these.

*(i) Branch-blindness, which is decisive and growth-independent.* By THEOREM S3-0 (§2.4), the
gauge-invariant algebra of **any** complex is the diagonal, and the diagonal cannot distinguish the
two branches at any stage or in the limit. A bigger complex has a bigger blind algebra. No amount
of cellular growth changes this.

*(ii) No unital connecting map exists along the carrier's own growth.* A unital \*-homomorphism
`M_m -> M_n` exists **iff `m | n`** (a unital representation of `M_m` on `C^n` is a direct sum of
`k` copies of the defining one, so `n = km`). The two growths the carrier suggests give:

```
wedge of full K1 copies at the root : |V_N| = 1 + 4N  ->  d_N = 5, 9, 13, 17, ...
      5 | 9 ?  False        9 | 13 ?  False       13 | 17 ?  False
wedge of unfilled triangles at root : |V_N| = 3 + 2N  ->  d_N = 5, 7, 9, 11, ...
      5 | 7 ?  False        7 | 9  ?  False        9 | 11 ?  False
```

The only linear \*-preserving map available is extension by zero, `X -> X (+) 0`, which is **not
unital** — and P-8 then destroys the system outright:

```
omega^{N+1}( I_N (+) 0 )      =  omega^N(I)  =  1
omega^{N+1}( 0   (+) I_new )  =  1 - 1       =  0
```

so the stage-`N+1` state must assign zero weight to every new cell. **P-8 forces the process to
stop.** Cellular growth of the carrier is not a directed system in the sense the requirements need.

*This is where the construction departs from the inherited template.* The template is *"cell-indexed
finite record algebras"* (`FOUNDING_DESIGN_V001.md:61–65`). Cell-indexing survives — `N` counts
cells. What does not survive is the idea that the cell's algebra is the **carrier's** algebra: on
this carrier it cannot be, by §2.4. The cell contributes a **record slot**, not a fibre.

**Not taken — tensor powers of the carrier algebra, `M_5^{(x)N`.** `5 | 25`, so this one *is* a
legitimate directed system, with limit `UHF(5^infinity)`. It is rejected as **not minimal**: it
adds 24 complex dimensions per cell where 3 are needed (`M_2` beyond the trivial slot), and §2.5
proves 4 is the floor. Custody's "exhibit the object, add nothing unforced" selects the minimum.
`M_5 (x) M_2^{(x)N}` embeds unitally in `M_5^{(x)(N+1)}` for `N >= 1`, so nothing is lost.

### 3.3 The write map — the durability map, exhibited as a unitary

**Definition.** Fix a ready section `s` (normalised) and a circuit count `k`. Put

```
x = M_dF^k s        y = M_c^k s        Z = <x, y>
```

and define the slot pair in `R = C^2`, with `z = |Z|`, `xi = arg Z`, `alpha = (1/2) arccos z`:

```
phi^F  =  ( cos alpha ,  sin alpha )
phi^C  =  e^{i xi} ( cos alpha , -sin alpha )              <phi^F, phi^C>  =  e^{i xi} cos 2 alpha  =  Z
```

**The write map `V`** is the unitary on `Gamma(L) (x) R = C^10` determined by

```
V ( x (x) |0> )  =  s (x) phi^F
V ( y (x) |0> )  =  s (x) phi^C
```

on the two-dimensional span, extended arbitrarily to the orthogonal complement.

**It exists because the Gram matrices agree** — that is exactly the isometry condition, and it is
exactly the statement that the write **relocates** the branch distinction and does not manufacture
it. In words: *the write moves the comparison out of the reusable finite carrier into a fresh slot,
and resets the carrier to ready.*

**Exhibited** at `f = 2.0, c = 1.1, p = (0.4, 0.15, 0.15, 0.15, 0.15), k = 1`:

```
Gram(domain) - Gram(codomain)              =  2.220e-16       (the isometry condition)
|| V* V - I_10 ||                          =  2.617e-15       V IS UNITARY
|| V( x (x) |0> ) - s (x) phi^F ||         =  3.054e-16
|| V( y (x) |0> ) - s (x) phi^C ||         =  4.594e-16
<phi^F, phi^C>  =  0.259878772772 - 0.318757783880j
Z_1             =  0.259878772772 - 0.318757783880j      agree to 5.55e-17
|| s ||         =  1.000000                              carrier RESET in both branches
```

The extension off the branch span is a free choice; two different extensions differ by
`|| V - V' || = 4.1688` and give the identical record (`|| V'(x (x) |0>) - s (x) phi^F || =
3.246e-16`). Nothing computed on this page depends on it. Recorded at C5 and FLAG F2.

### 3.4 The system, and that its maps are what they must be

```
A_N  =  M_5(C) (x) M_2(C)^{(x)N}                 iota_N(X) = X (x) I_2

    N :   0       1       2       3        4         5          6
dim :   25     100     400    1600     6400     25600     102400          =  25 * 4^N
```

```
unital    :  || iota(I_20) - I_40 ||                 =  0.000e+00
*-hom     :  || iota(XY) - iota(X) iota(Y) ||        =  0.000e+00
adjoints  :  || iota(X*) - iota(X)* ||               =  0.000e+00
isometric :  |  ||iota(X)||_op - ||X||_op  |         =  5.329e-15
injective :  X (x) I = 0  =>  X = 0                  (immediate)
```

**The limit.** `M_5 (x) M_2^{(x)N} = M_{5 . 2^N}` and `5.2^N | 5.2^{N+1}`, so the system is a UHF
system. A **UHF algebra** is a direct limit of full matrix algebras along unital embeddings; it is
determined by its supernatural number.

```
A_infinity  =  UHF( 5 * 2^infinity )
```

unital, simple, with a unique tracial state, and **infinite-dimensional** — `dim A_N = 25 . 4^N`
diverges. The record subalgebra is `R_infinity = UHF(2^infinity)`, the CAR algebra.

### 3.5 The schedule, and the cell

A **cell** is one formation event. Cell `n` consists of: the carrier, in ready state `s`, runs
`k_n` circuits of each of K1's two loops, producing the branch pair `(M_dF^{k_n} s, M_c^{k_n} s)`;
then the write `V_n` of §3.3 transfers the comparison into slot `n` and resets the carrier.

A **schedule** is the sequence `(k_n)` of circuit counts. Two are displayed throughout:

- **Schedule B, the canonical clock: `k_n = n`.** Cell `n` reads the carrier at carrier-time `n`.
  Circuit count is the carrier's own discrete time (audit COR-F), so this is the schedule the
  carrier supplies. **This is the primary schedule** and it is the one that confronts the
  recurrence head-on, because its cell values `Z_n` are exactly the recurrent quantity.
- **Schedule A, uniform: `k_n = 1`.** Every cell reads one circuit. `Z_n = Z_1` for all `n`.

Every theorem below is stated for a general schedule; the two are used to show nothing hinges on
the choice, and §4.6 shows what does.

---

## 4. Q3 — DOES THE LIMIT ESCAPE RECURRENCE?

### 4.1 The frozen carrier: the recurrence, reproduced from scratch

The one-cell overlap, from W-01, verified here against direct matrix action on `C^5`:

```
Z_k  =  e^{ik(c-f)} p_0  +  e^{-ikf} (p_1+p_2)  +  e^{ikc} (p_3+p_4)

closed form vs matrix action, max deviation over k <= 200 :   1.538e-14
```

At the audit's own test point `f = 2.0, c = 1.1, p = (0.4, 0.15, 0.15, 0.15, 0.15)`:

```
min |Z_k| over k <= 400   =  0.024654   at k = 42        W-01 records 0.0247 at n = 42     MATCH
sup |Z_k| over k <= 4000  =  0.999941   at k = 377       W-01 records 0.99994              MATCH
|Z_k| > 0.99  at  37  of the first 4000 circuit counts
```

**It fires, it un-fires, and it keeps coming back.** Reproduced here from the definitions, not
inherited.

### 4.2 The record functional

**Definition.** The stage-`N` branch vectors and the **record functional**:

```
Xi_F^N  =  s (x) phi_1^F (x) phi_2^F (x) ... (x) phi_N^F
Xi_C^N  =  s (x) phi_1^C (x) phi_2^C (x) ... (x) phi_N^C

Omega_N  =  < Xi_F^N , Xi_C^N >  =  prod_{n=1}^{N} Z_{k_n}
```

The product form is immediate: the slot factors are independent and `<s,s> = 1`.

**There is no conflict with unitarity.** The two branches are two *different* unitary histories of
the *same* initial state `s (x) |0>^N` — branch F applies `M_dF` at each cell, branch C applies
`M_c`. Distinct unitaries applied to a common state may produce any overlap. Nothing here is
non-unitary, and no environment has been introduced.

### 4.3 THE COMPUTATION — schedule B, `k_n = n`

```
     N        |Z_N|         |Omega_N|            (1/N) log|Omega_N|
     1      0.411271     4.112706e-01              -0.888504
     2      0.470386     1.934559e-01              -0.821353
     5      0.364118     1.772339e-02              -0.806574
    10      0.776953     6.540411e-04              -0.733234
    20      0.247337     2.514545e-07              -0.759800
    42      0.024654     2.505486e-15              -0.800483
    50      0.681600     1.571127e-18              -0.819895
   100      0.350139     4.036647e-34              -0.768925
   200      0.610201     1.405872e-67              -0.769663
   400      0.600798    2.936168e-134              -0.768673
  1000      0.580708      < 1e-300                 -0.768205
  2000      0.573418      < 1e-300                 -0.768470
  4000      0.360282      < 1e-300                 -0.766802
```

```
Lyapunov rate      lambda = lim (1/N) log|Omega_N|          =  -0.766802   (N = 4000)
                                                            =  -0.767026   (N = 200000)
torus mean of log|Z| (Weyl equidistribution prediction)     =  -0.767507        AGREES

|Omega_N| non-increasing at every step  :  True
|Omega_N| STRICTLY decreasing at every step :  True
```

**Read the second column against the third.** At `N = 10` the one-cell overlap is `0.777` — nearly
un-fired. At `N = 50` it is `0.682`. The single-cell quantity wanders back up to within `6e-5` of
`1`. **The record functional does none of that.** It falls monotonically through `1e-134` by
`N = 400` and never rises.

### 4.4 THE THEOREMS

> **THEOREM S3-1 (monotone, and the exact criterion).** Let `z_n = |Z_{k_n}| in [0,1]`. Then
> `|Omega_N| = prod_{n<=N} z_n` is non-increasing in `N`, and
>
> ```
> lim_N |Omega_N|  >  0      iff      every z_n > 0   and   sum_n (1 - z_n)  <  infinity.
> ```
>
> **Proof.** `|Z| <= 1` by the triangle inequality on three unit-modulus numbers with weights
> summing to 1, so each factor lies in `[0,1]` and the product is non-increasing. The criterion is
> the standard convergence test for infinite products: for `a_n in [0,1)`, `prod (1-a_n)` has a
> non-zero limit iff `sum a_n < infinity`; a single zero factor sends the product to `0`. ∎

> **THEOREM S3-2 (the canonical clock always escapes).** Take `k_n = n` and a ready state with
> `p_0, p_1+p_2, p_3+p_4 > 0`. If `(W_F, W_C) != (1,1)` then `sum_n (1 - |Z_n|) = infinity` and
> `(1/N) log|Omega_N| -> lambda < 0`. Hence `|Omega_N| -> 0` exponentially.
>
> **Proof.** Write `u = conj(W_F)`, `v = W_C`, so `Z_n = (uv)^n p_0 + u^n q + v^n r` with
> `q = p_1+p_2`, `r = p_3+p_4`. The sequence `(u^n, v^n)` runs over the cyclic group generated by
> `(u,v)` in `T^2` and equidistributes (Weyl) on its closure `T`, a compact group with Haar measure
> `mu`. `|Z| = 1` requires the three unit vectors `uv, u, v` to coincide on the support of `p`;
> with `p` fully supported that forces `u = v = 1`, the identity of `T`. If `T` is infinite,
> `mu({e}) = 0`; if `T` is finite of order `m >= 2`, `mu({e}) = 1/m < 1`. In both cases `|Z| < 1`
> on a set of positive measure, so `lambda = int_T log|Z| dmu < 0` (the logarithmic singularity at
> zeros of `Z` is integrable, and `lambda = -infinity` there is permitted and means immediate exact
> formation). `T` is trivial iff `u = v = 1` iff `W_F = W_C = 1`. Divergence of
> `sum (1 - |Z_n|)` follows because `1 - |Z_n|` is bounded below by a positive constant on a set of
> positive density. ∎

> **THEOREM S3-3 (recurrence is gone, not merely unproved).** If `z_n < 1` for all `n`, then
> `|Omega_N|` is **strictly decreasing**, hence **injective as a function of `N`**. The record
> functional therefore returns to no value it has ever taken, at any later stage, ever. In
> particular
>
> ```
> sup_{M >= N} |Omega_M|  =  |Omega_N|          for every N.
> ```
>
> **Proof.** Immediate from strict monotonicity. ∎

**This is the answer the brief demands.** It is not "the theorem's hypothesis is unmet". The
returning behaviour is *gone*, and the reason is structural: the record functional is a **product**
over cells, and a product of factors in `[0,1)` cannot rise. **A recurrent factor cannot undo a
monotone product.**

**And the hypothesis is also unmet, for the record.** The inherited theorem binds closed systems
with finitely many unitary degrees of freedom. `dim A_N = 25 . 4^N -> infinity`, and
`A_infinity = UHF(5 . 2^infinity)` is infinite-dimensional. The system is **closed** — no
environment, no bath, no open dynamics has been added. **The escape is by unbounded cell count
alone.**

### 4.5 On K1's own published connection — exact, and the sharpest exhibit on this page

S1 `:114–115`: `W_F = -1`, `W_C = -i`. Ready state `p = (1/2, 0, 0, 1/4, 1/4)` (W-01's).
Then `Z_k = conj(W_F^k) W_C^k / 2 + W_C^k / 2`, evaluated in exact arithmetic:

```
  k :   Z_k                |Z_k|      Omega_k
  1 :   i/2 - i/2  =  0      0        0        EXACT
  2 :          -1            1        0
  3 :             0          0        0
  4 :          +1            1        0
  5 :             0          0        0
  6 :          -1            1        0
cross-check by direct matrix action, k = 1..8, max deviation :  8.951e-16
```

**Read it.** `|Z_2| = 1` **exactly** — the one-cell quantity returns *completely*, not nearly, at the
second circuit. Total recurrence. And `Omega_k = 0` for every `k >= 1`: **the record is exact and
permanent from the first cell.** On K1's own published connection the crossing succeeds
immediately, and the carrier's complete recurrence at `k = 2` is simply irrelevant to it.

### 4.6 THE PRICE — an adversarial schedule defeats the construction

Theorem S3-1 is an *iff*, so the escape is conditional and the condition is checkable. The
**deficiency** of cell `n` is `1 - z_n`. If a schedule is locked to the carrier's near-recurrence
times, the deficiencies are summable and the product does not vanish.

Record-breaking near-recurrences of the frozen carrier at `f = 2.0, c = 1.1`, searched to
`k = 200000`:

```
   k =      6    |Z_k| = 0.972113557    deficiency = 2.789e-02
   k =     63    |Z_k| = 0.979059564    deficiency = 2.094e-02
   k =    154    |Z_k| = 0.994783559    deficiency = 5.216e-03
   k =    377    |Z_k| = 0.999941230    deficiency = 5.877e-05
   k =   6723    |Z_k| = 0.999948938    deficiency = 5.106e-05
   k =   7100    |Z_k| = 0.999999729    deficiency = 2.708e-07
   k =  99023    |Z_k| = 0.999999855    deficiency = 1.450e-07
   k = 106123    |Z_k| = 0.999999981    deficiency = 1.949e-08

sum of deficiencies along this thinned schedule (last 40 terms) =  1.533095   (converging)
product along it                                                =  0.117102   bounded away from 0
```

Because the orbit is dense, the deficiencies can be driven below `2^-j`, so a schedule with
`sum (1 - z_n) < infinity` **always exists**, and along it there is **no escape at all**.

**Stated plainly: the crossing is not a property of the algebra alone. It is a joint property of
the algebra and the cell schedule.** The canonical clock `k_n = n` always escapes (Theorem S3-2);
so does the uniform schedule `k_n = 1`, since `z_n = |Z_1| < 1` is constant. What fails is a process
whose cells are timed to the carrier's own returns. This is a genuine structural condition the
crossing imposes on any process that would use it, and it is a result, not a caveat.

For the canonical schedule the divergence is linear and unmistakable:

```
sum_{n <=     10} (1 - |Z_n|)  =       4.850
sum_{n <=    100} (1 - |Z_n|)  =      47.240
sum_{n <=   1000} (1 - |Z_n|)  =     469.347
sum_{n <=  10000} (1 - |Z_n|)  =    4692.037
sum_{n <= 100000} (1 - |Z_n|)  =   46918.264
```

### 4.7 Gauge invariance of the whole construction

Every quantity above is a function of `(W_F, W_C, p)`, which are gauge invariants (S1 `:76–82`;
W-01). Checked directly by transforming the ready section by random vertex phases and recomputing
the whole 30-cell product:

```
|Omega_30| = 1.192688e-10 ;  max deviation over 8 random gauge transformations = 5.147e-25
```

---

## 5. Q4 — THE NINE, ONE AT A TIME

Test point throughout: `f = 2.0`, `c = 1.1`, `p = (0.4, 0.15, 0.15, 0.15, 0.15)`, schedule
`k_n = n`, `lambda = -0.7668`.

### 5.1 P-8 — INDUCTIVE COMPATIBILITY — **HOLDS**

*Requirement:* `omega^{N+1}( iota_N(A) ) = omega^N(A)`.

*Demonstration.* `omega^N(X) = <Xi^N, X Xi^N>` for the stage-`N` branch vector. Since
`Xi^{N+1} = Xi^N (x) phi_{N+1}` and `iota_N(X) = X (x) I_2`,

```
omega^{N+1}(X (x) I)  =  <Xi^N (x) phi, (X (x) I)(Xi^N (x) phi)>  =  <Xi^N, X Xi^N> . ||phi||^2  =  omega^N(X)
```

*Computed*, random Hermitian `X` at each stage, both branches:

```
N = 1 :  dim A_N =   100   max | omega^{N+1}(X (x) I) - omega^N(X) |  =  2.312e-16
N = 2 :  dim A_N =   400                                             =  4.965e-16
N = 3 :  dim A_N =  1600                                             =  2.221e-15
N = 4 :  dim A_N =  6400                                             =  2.221e-15
```

**HOLDS exactly, and it is the cheapest of the nine** — a product state along a tensor system is
automatically compatible. It is also the property that killed cellular growth (§3.2).

### 5.2 P-6 — ASYMPTOTIC CENTRALITY — **HOLDS, with the inherited constant**

*Requirement:* `|| [M_N, O] ||  <=  2m ||O|| / N`.

*Definition.* An element is **asymptotically central** if its commutator with every fixed local
observable tends to zero. Take the mean pointer `M_N = (1/N) sum_{n=1}^{N} sigma_x^{(n)}`.

*Proof.* For `O in A_m` (carrier plus the first `m` slots), `sigma_x^{(n)}` commutes with `O` for
`n > m`. So `[M_N, O] = (1/N) sum_{n<=m} [sigma_x^{(n)}, O]`, a sum of `m` commutators each of norm
at most `2 ||sigma_x|| ||O|| = 2||O||`. ∎

*Computed*, `||O|| = 1`:

```
   m    N     ||[M_N,O]||      bound 2m||O||/N     ratio
   1    2       0.432914          1.000000        0.4329
   1    4       0.216457          0.500000        0.4329
   1    6       0.144305          0.333333        0.4329
   1    8       0.108228          0.250000        0.4329
   1   10       0.086583          0.200000        0.4329
   2    2       0.464446          2.000000        0.2322
   2    4       0.232223          1.000000        0.2322
   2    6       0.154815          0.666667        0.2322
   2    8       0.116112          0.500000        0.2322
   2   10       0.092889          0.400000        0.2322
```

The inherited bound holds with room. `1/N` decay confirmed to four figures.

### 5.3 P-1 — DURABILITY AND IRREVERSIBILITY JOINTLY — **HOLDS IN THE LIMIT; FAILS AT EVERY FINITE STAGE**

*Durability* is §4.3–4.4: `|Omega_N| <= e^{lambda N}`, monotone, no return.

*Irreversibility.* In the limit the two record states are **disjoint** (§5.7). If `u` is a unitary
in `A_infinity` then `omega . Ad(u)` is quasi-equivalent to `omega` — the GNS representation of
`omega . Ad(u*)` is `pi_omega` composed with an inner automorphism, unitarily implemented by
`pi_omega(u)`. Disjoint states are not quasi-equivalent. **Therefore no unitary in `A_infinity`
carries `omega_F` to `omega_C`.** Every finite product of cell writes lies in some `A_N subset
A_infinity`, so no dynamics assembled from the construction's own moves can undo the record.

*And at finite `N` it is reversible* — exhibited, because this is the sharp half of the result:

```
N = 1  dim =  10   ||U*U - I|| = 2.88e-15   || U Xi_F - Xi_C || = 2.88e-16    U in A_N : YES
N = 2  dim =  20   ||U*U - I|| = 5.04e-15                       = 5.05e-16    U in A_N : YES
N = 3  dim =  40   ||U*U - I|| = 7.97e-15                       = 6.38e-16    U in A_N : YES
N = 4  dim =  80   ||U*U - I|| = 1.44e-14                       = 4.62e-16    U in A_N : YES
N = 5  dim = 160   ||U*U - I|| = 2.30e-14                       = 8.71e-16    U in A_N : YES
```

Any two unit vectors in a finite-dimensional space are related by a unitary of that space. **So
irreversibility is not available at any finite stage, on any construction whatever.** It is a limit
property or it is nothing. **GRADE: HOLDS in the limit. FAILS at every finite `N`, necessarily.**

### 5.4 P-3 — THRESHOLDED NON-RETURN — **HOLDS**

Given a threshold `eps`, the first stage below it, and — by strict monotonicity (Theorem S3-3) —
the supremum over all later stages, which equals the value at the crossing:

```
   eps        first N with |Omega_N| < eps     |Omega_N| there      sup_{M>=N}|Omega_M|
  1e-1                  3                       7.613e-02          = |Omega_N|
  1e-3                  9                       8.418e-04          = |Omega_N|
  1e-6                 20                       2.515e-07          = |Omega_N|
  1e-12                37                       5.203e-13          = |Omega_N|
  1e-30                90                       3.368e-31          = |Omega_N|
  1e-100              300                       4.031e-101         = |Omega_N|
asymptotic estimate  ceil(log eps / lambda)  gives  4, 10, 19, 37, 91, 301   —  agrees to +/-1
```

**Non-return is exact here, not asymptotic**: monotonicity makes `sup_{M>=N}` equal to `|Omega_N|`
identically. Contrast the frozen carrier, whose `sup_{k>=42} |Z_k| = 0.999941`.

### 5.5 P-2 — PERSISTENCE UNDER LATER CELLS — **HOLDS**

The cell-`n` write `V_n` acts on `Gamma(L) (x) R_n` and therefore lies in `M_5 (x) M_2^{(n)}`. It
commutes with every observable of slot `m != n`. Computed for cell 2 against `sigma_x` on slot 1,
inside `M_5 (x) M_2 (x) M_2`:

```
|| [ V_(cell 2) , sigma_x on slot 1 ] ||  =  0.000e+00
```

A record written at cell `m` is never touched again — its slot state is exactly `phi_m^{F/C}` at
every later stage. Combined with P-8, `omega^M |_{A_m} = omega^m` for all `M >= m`. **The record is
not merely undisturbed, it is stage-independent.**

### 5.6 P-4 — RECOVERABILITY — **HOLDS**

The pointer `M_N = (1/N) sum sigma_x^{(n)}` takes value `+s̄_N` in branch F and `-s̄_N` in branch C,
with `s̄_N = (1/N) sum_{n<=N} sqrt(1 - |Z_n|^2)`. Fluctuation in either branch is
`(1/N) sqrt(sum |Z_n|^2) <= 1/sqrt(N)`.

```
     N     omega_F(M_N)   omega_C(M_N)     gap        s.d. of M_N in a branch
     1       +0.911513      -0.911513    1.823027          0.411271
     5       +0.882705      -0.882705    1.765409          0.208487
    10       +0.804660      -0.804660    1.609320          0.175377
    50       +0.797975      -0.797975    1.595949          0.081230
   200       +0.788639      -0.788639    1.577279          0.041121
  1000       +0.785987      -0.785987    1.571973          0.018434
  4000       +0.785995      -0.785995    1.571990          0.009219
```

Signal-to-noise `~ sqrt(N)`. Reading the sign of `M_N` recovers the branch with error
`O(e^{-c N})`. **The record is readable, by a single observable, from the record algebra alone.**

### 5.7 P-7 — SECTOR-HOOD — **HOLDS**

*Definitions.* Two states are **disjoint** if no sub-representation of one GNS representation is
unitarily equivalent to a sub-representation of the other; a **sector** is a quasi-equivalence
class. Disjointness is what makes a record classical: the limit pointer becomes a central element
taking a definite value on each sector.

*Elementary proof, self-contained.* By P-6 the pointer `M_N` is asymptotically central; by P-4 its
fluctuation vanishes as `1/sqrt(N)`. So in the limit `M_infinity` is a **sharp central observable**
with value `+s̄` on `omega_F` and `-s̄` on `omega_C`, where

```
s̄  =  lim (1/N) sum_{n<=N} sqrt(1 - |Z_n|^2)  =  0.785995   (N = 4000)   0.785940   (N = 50000)
```

A central element taking two different values distinguishes two sectors. Since `s̄ > 0` whenever the
connection is non-trivial, the two branches lie in **different sectors of the limit**.

`s̄` is a function of the ready state, and is *not* a constant of the construction (all rows at
`N = 50000`):

```
   p = (0.4000, 0.3000, 0.3000)   s̄ = 0.785940      lambda = -0.767043
   p = (0.5000, 0.0000, 0.5000)   s̄ = 0.636630      lambda = -0.692761
   p = (0.3333, 0.3333, 0.3333)   s̄ = 0.790162      lambda = -0.776733
   p = (0.8000, 0.1000, 0.1000)   s̄ = 0.560484      lambda = -0.223149
   p = (0.1000, 0.4500, 0.4500)   s̄ = 0.730358      lambda = -0.727616
   p = (0.2000, 0.5000, 0.3000)   s̄ = 0.760015      lambda = -0.693152
```

(The first value is near `pi/4 = 0.785398`. It is a coincidence of one ready state, it is not
stable across the table, and **nothing on this page uses it.** Recorded so that it cannot be
picked up later as if it meant something.)

*Second, independent route, cited but not load-bearing.* The classical criterion for product states
on an infinite tensor product (von Neumann 1939; Powers 1967; Bures 1969) makes two pure product
states quasi-equivalent iff `sum_n (1 - |<psi_n, phi_n>|) < infinity`. Here that sum is
`46918.264` at `N = 1e5` and grows linearly (§4.6), so the states are **disjoint** by that route
too. This is external mathematics named by author, not by digest, and the elementary route above
does not depend on it.

### 5.8 P-5 — REDUNDANCY — **HOLDS**

Every subset of slots of positive density carries the whole record. Both the product and the
pointer estimate are computed over sub-collections of the first 4000 slots:

```
   subset             |S|      prod_{n in S} |Z_n|      pointer estimate
   every 2nd slot     2000        < 1e-300                 +0.786057
   every 10th          400        1.931e-133               +0.786590
   every 100th          40        3.966e-13                +0.796081
   random 5%           200        7.955e-66                +0.783661
   ---------------------------------------------------------------------
   full record        4000                                 +0.785995
```

Forty slots out of four thousand already give `|Omega_S| = 4e-13` and the pointer to 1.3%. **The
record is not localised in any slot and destroying any density-zero set of slots destroys nothing.**

### 5.9 P-9 — ORTHOGONAL REDUCED SUPPORTS — **HOLDS IN THE LIMIT; FAILS AT EVERY FINITE STAGE EXCEPT ON W-01's FIRING LOCUS**

**This was the sharpest target: recorded false of record in the predecessor and produced by
nothing.**

*The reduced states.* Restricting `Xi_F^N`, `Xi_C^N` to the record algebra `R_N` gives two pure
states with overlap `Omega_N`. Their supports are the rank-one projections onto those vectors;
those are orthogonal iff `Omega_N = 0`. The state-norm distance is
`|| omega_F^N - omega_C^N || = 2 sqrt(1 - |Omega_N|^2)`, maximal at 2.

```
     N        |Omega_N|          || omega_F^N - omega_C^N ||
     0      1.000000e+00            0.000000000000
     1      4.112706e-01            1.823026602854
     2      1.934559e-01            1.962217957798
     5      1.772339e-02            1.999685856739
    10      6.540411e-04            1.999999572230
    20      2.514545e-07            2.000000000000
    42      2.505486e-15            2.000000000000
   100      4.036647e-34            2.000000000000
```

*In the limit:* the states are disjoint (§5.7), so their central supports in the enveloping von
Neumann algebra are **orthogonal central projections** — orthogonal reduced supports, in the
strongest available sense. **HOLDS.**

*At finite `N`, generically:* `Omega_N != 0`, so the supports are **not** orthogonal. **FAILS.**

*At finite `N` on W-01's firing locus:* if `Z_{k_n} = 0` for some `n` — exactly W-01's convex-hull
criterion, exactly the condition that `0` lies in the convex hull of
`{e^{i(c-f)}, e^{-if}, e^{ic}}` — then `Omega_N = 0` **exactly** for all `N >= n`. On K1's own
published connection this occurs at the **first cell** (§4.5). **HOLDS exactly, at finite stage,
on K1's own connection.**

**Why the predecessor recorded it false.** P-9 is a limit property. Evaluated at any finite stage
of any construction it is false generically, because two pure states of a finite-dimensional space
with non-zero overlap have non-orthogonal supports and no unitary construction can change that. The
predecessor's `false` was a correct finite-stage evaluation of a property that only exists in the
limit. **The crossing is precisely the passage that makes it true.**

### 5.10 THE NINE, SUMMARISED

```
P-1  durability + irreversibility     HOLDS in the limit  ·  FAILS at every finite N (necessarily)   §5.3
P-2  persistence under later cells    HOLDS                                                          §5.5
P-3  thresholded non-return           HOLDS  (exact, by monotonicity, not asymptotic)                §5.4
P-4  recoverability                   HOLDS  (single pointer, SNR ~ sqrt N)                          §5.6
P-5  redundancy                       HOLDS  (every positive-density subset suffices)                §5.8
P-6  asymptotic centrality            HOLDS  (inherited bound 2m||O||/N, with room)                  §5.2
P-7  sector-hood                      HOLDS  (two sectors, label +/- s̄, elementary proof)            §5.7
P-8  inductive compatibility          HOLDS  (exact; and it is what kills cellular growth)           §5.1
P-9  orthogonal reduced supports      HOLDS in the limit  ·  HOLDS at finite N on W-01's firing
                                      locus  ·  FAILS at finite N generically                        §5.9
```

**Nine of nine in the limit. Two of them — P-1 and P-9 — provably unavailable at any finite
stage.** That is not a weakness of this construction; §5.3 shows it is a theorem about all finite
constructions. It is the precise content of the phrase *"a missing crossing"*.

---

## 6. Q5 — DOES THE FORMATION CONDITION TRANSPORT?

**Yes, multiplicatively, and firing is absorbing.**

**(a) The stage functional.** W-01's condition on K1 is `Z = 0`. Its stage-`N` counterpart is the
record functional `Phi_N = Omega_N = <Xi_F^N, Xi_C^N>`, which is the *same object* — the overlap of
the two branch evolutions of a common ready state — evaluated on the stage-`N` algebra. At `N = 1`
with `k_1 = 1` it is exactly `Z_1`, W-01's quantity.

**(b) The compatibility.** By construction

```
Phi_{N+1}  =  Phi_N  .  Z_{k_{N+1}}
```

so the family `(Phi_N)` is determined by its predecessors and one new cell. And by P-8 the
comparison functional itself is compatible with the system's maps:

```
( omega_F^{N+1} - omega_C^{N+1} ) . iota_N   =   omega_F^N - omega_C^N        (max deviation 2.2e-15, §5.1)
```

**So the formation condition is a compatible family of functionals on the directed system and
defines a functional on `A_infinity`.** Firing at stage `N` is a statement about `A_N`, and it
remains that same statement about `A_N` when read at stage `N+1` and in the limit.

**(c) Firing is absorbing.** `Phi_N = 0` implies `Phi_M = 0` for all `M >= N`, immediately from
(b). **Once formed, formed forever.** That is exactly the durability the predecessor's gate wanted
and could not deliver.

**(d) What the limit does to W-01's criterion.** W-01's convex-hull criterion is the condition for
**exact firing in a single cell**, and it cuts out a proper subset of the invariant torus. In the
limit that binary criterion becomes a **rate**:

```
Phi_infinity  =  0     for every connection with (W_F, W_C) != (1,1)        [Theorem S3-2]
rate          =  lambda  =  int over the orbit closure of log|Z|
```

The convex-hull locus is where `lambda = -infinity`. Everywhere else off the trivial connection
`lambda` is finite and negative. **The limit converts a measure-zero exact-firing criterion into a
formation rate defined on the whole torus minus one point.**

**(e) The trivial-connection limit, preserved — and one honest complication.** At `W_F = W_C = 1`,
`Z_n = 1` for all `n`, `Omega_N = 1` for all `N`: **no formation, ever.** The design's S2 external
contact (`FOUNDING_DESIGN_V001.md:117–118`) survives the crossing intact. And at every *fixed*
stage the limit is continuous:

```
     (f, c)          |Z_1|      |Omega_10|     |Omega_100|      lambda_B     lambda_A = log|Z_1|
  (0.0000, 0.0000)  1.000000   1.0000e+00    1.0000e+00        +0.000000        +0.000000
  (0.0010, 0.0006)  1.000000   9.9992e-01    9.3556e-01        -0.760309        -0.000000
  (0.0100, 0.0060)  0.999980   9.9245e-01    1.0686e-03        -0.765479        -0.000020
  (0.0500, 0.0300)  0.999508   8.2631e-01    9.4865e-35        -0.765183        -0.000492
  (0.2000, 0.1000)  0.992964   4.8838e-02    8.2888e-31        -0.681963        -0.007061
  (1.0000, 0.6000)  0.813939   5.3873e-04    3.6970e-34        -0.765229        -0.205869
  (2.0000, 1.1000)  0.411271   6.5404e-04    4.0366e-34        -0.767026        -0.888504
```

`|Omega_N| -> 1` as the connection trivialises, for each fixed `N`. **But the asymptotic rate
`lambda_B` under the canonical clock is discontinuous at the trivial connection** — for any dense
orbit it equals the orbit-closure average however small `(f,c)` are, because the carrier's clock
eventually explores the whole closure. Under the uniform schedule the rate is `log|Z_1|` and *is*
continuous. **The two schedules are distinguished exactly at the trivial connection, and nowhere
else.** Recorded at FLAG F5; I do not know which is physical.

**(f) The rate depends on the arithmetic of the holonomies, not only their values** — a result that
belongs to S4 and is recorded here because it was computed here:

```
  f = 2.00000  c = 1.10000   lambda_B = -0.767026    f, c, 2pi rationally independent: orbit dense in T^2
  f = 2.00000  c = 2.00000   lambda_B = -1.203587    f = c: orbit confined to the diagonal subtorus
  f = 3.14159  c = 1.57080   lambda_B = -0.804719    both rational multiples of pi: finite orbit
  f = 3.14159  c = 4.71239   lambda_B = -0.804719    finite orbit, order 4
```

And winding periodicity — the design's S4 contact (`:119`), checked early because it was free:

```
|Z_1(f, c)| = 0.411270593796    |Z_1(f + 2pi, c - 4pi)| = 0.411270593796    difference 1.67e-16
```

---

## 7. WHAT S3 DELIVERS, AND WHAT IT DOES NOT

**Delivers.** The record algebra on K1 with its counts, derived and saturated. A theorem that no
carrier can hold the record, general to all such complexes. A directed system, exhibited with its
objects, its maps checked to be unital \*-homomorphisms, and its limit identified. A durability map
exhibited as an explicit unitary. A proof and a computation that recurrence is gone — by
monotonicity, not by hypothesis. All nine requirements graded with demonstrations, including the
two that provably cannot hold at finite stage. The formation condition's transport, its absorbing
character, and its conversion from criterion to rate.

**Does not deliver, and does not claim to.** Any variation of the electromagnetic structure beyond
the four rate values computed in passing at §6(f) — that is S4. Any dimensionless coupling slot —
that is S5. Any account of *why* the process should supply cells at all, or on what schedule; §4.6
shows the schedule is load-bearing and this page does not derive it. Any connection to the
predecessor's magnitude demands (`FOUNDING_DESIGN_V001.md:77–84` already says so).

**S3's falsifier was "no construction satisfies the nine, or none escapes recurrence."** A
construction satisfying all nine in the limit is exhibited, and it escapes recurrence. The
falsifier is not triggered. **The design's S3 external contact — "the recurrence no-go, answered
rather than checked" — is discharged at §4.4:** the no-go's hypothesis is shown unmet *and* the
returning behaviour is shown to be absent, by strict monotonicity of an injective sequence.

---

## CHOICE LEDGER

| # | Choice | Alternatives | Why | Status |
|---|---|---|---|---|
| C1 | Growth by a **record slot per cell** (tensor system) | (a) repeated circuits; (b) more unfilled triangles at the root; (c) a chain of complexes; (d) tensor powers `M_5^{(x)N}` | (a) computed: `dim = 3` at `N = 1` and at `N = 100`, no growth at all. (b),(c) killed twice over — by THEOREM S3-0 (the gauge-invariant algebra of any complex is branch-blind) and by the unital-embedding arithmetic (`5 ∤ 9`, `5 ∤ 7`), with P-8 then forcing the process to stop. (d) is legitimate but not minimal: 24 extra complex dimensions per cell where §2.5 proves 4 is the floor | **closed** — all four alternatives computed, §3.2 |
| C2 | Record slot `R = M_2(C)` | `C^2` (abelian); `C^3` (the three-way split); `M_3`; `M_5` | pure states of an abelian algebra have overlap `0` or `1` only, so an abelian slot cannot hold a general comparison value; two unit vectors span two dimensions, so `M_2` is the exact minimum. And `2` is K1's own count twice over — loops (`b1 + faces = 1 + 1`) and invariants (`6 - 4`) | **closed** — §2.5 |
| C3 | Schedule `k_n = n`, the carrier's own circuit clock | `k_n = 1`; any other schedule | circuit count is carrier-supplied discrete time (audit COR-F); `k_n = n` is the schedule that makes the cell values *equal to the recurrent quantity*, so the escape is tested against the obstruction rather than around it | **closed, with a live condition** — Theorem S3-1 is stated for all schedules and §4.6 exhibits one that defeats the construction |
| C4 | Display parameters `f = 2.0, c = 1.1, p = (0.4,0.15,0.15,0.15,0.15)` | any | they are the audit's own §4.9 test point, so the frozen-carrier numbers here are directly comparable with W-01's `0.0247 at n=42` and `0.99994` — both reproduced. Every governing statement is a theorem, not these numbers | **closed** |
| C5 | Slot pair `phi^F = (cos a, sin a)`, `phi^C = e^{i xi}(cos a, -sin a)` | any pair with overlap `Z` | symmetric under branch exchange, which makes the pointer `sigma_x` take manifestly opposite values `+/- sqrt(1-z^2)` and gives P-7 a signed sector label at no cost | **closed** — nothing else on the page depends on it |
| C6 | Unitary extension of the write off the branch span | any | free; two extensions differing by `||V - V'|| = 4.17` give identical records | **closed** — verified extension-independent, §3.3 |

---

## IMPORT AUDIT

| # | Notion | Source | Defined on this page? | Survives without the import? |
|---|---|---|---|---|
| I1 | The nine requirements P-1…P-9 | `INHERITED_FROM_THE_PREDECESSOR.md` digest `96c1d3…` — a **specification**, custody §7 | graded here, each with its own demonstration | the construction stands; the grading would have no target |
| I2 | W-01's formation condition and its recurrence figure | `REGISTER_V001.md` digest `95ecd6…` | **re-derived from scratch**, §4.1: closed form checked against matrix action to `1.538e-14`; `0.024654 at k = 42` and `0.999941` reproduced | **yes** — it is recomputed, not inherited |
| I3 | K1's topology, gauge count, invariants, worked instance | `S1_CARRIER_K1_V001.md` digest `3eb703…`, lines cited inline | recomputed at §2.1–2.3 and §4.5 | yes |
| I4 | The canonical three-way split | W-01 (register), derived there from incidence | **independently re-derived** here as `C*(M_dF, M_c) = C^3`, §2.3, with a caveat W-01 lacks | yes |
| I5 | The inductive-limit template | `FOUNDING_DESIGN_V001.md:61–65` digest `ca25b7…` | used as a lead; **departed from** at §3.2 (the cell contributes a record slot, not a fibre) | yes |
| I6 | Direct limit, UHF, disjointness, quasi-equivalence, asymptotic centrality | standard C\*-algebra theory | all defined on this page, §3.1, §3.4, §5.2, §5.7 | n/a |
| I7 | Product-state quasi-equivalence criterion (von Neumann; Powers; Bures) | external mathematics, named by author, **not by digest** | used only as a **second, non-load-bearing** route in §5.7; the elementary mean-observable proof is displayed and independent | yes |
| I8 | **Alpha** | — | **NOT ENGAGED.** No coupling, no measured number, no target-driven selection anywhere on this page. `lambda` and `s̄` are free, connection- and state-determined, and are fitted to nothing | n/a |

---

## FLAG BLOCK — DEFECTS IN MY OWN CONSTRUCTION

**F1 — The write map is gauge-covariant, not gauge-invariant.** `V` is built from the ready section
`s`, which transforms under gauge. The correct statement is `V_{g.data} = (g (x) I) V (g (x) I)*`:
the write is a covariant *function of the data*, like a connection, not an invariant operator.
Every recorded quantity is gauge-invariant (checked, `5.147e-25`, §4.7), but the operator is not,
and I have not proved the covariance relation in general — only checked the invariance of the
outputs.

**F2 — The crossing is not unique.** The unitary extension off the branch span is arbitrary
(§3.3), and the choice of `M_2` basis is arbitrary (C5). What is pinned down is the record content;
the map that writes it is a family.

**F3 — The escape is conditional and the condition is not derived.** §4.6 exhibits a schedule
locked to the carrier's near-recurrence times along which the construction fails completely. This
page proves the canonical clock escapes; it does **not** derive that a physical process must use
the canonical clock. That is a real gap and it is where an adversary should attack first.

**F4 — This is not a zero-addition construction. S2 was; S3 is not.** The cost is exactly one
qubit per cell. §2.5 proves that is the floor given a unitary write and a reset carrier, and §2.4
proves the carrier cannot pay it. But "minimum, and proved minimum" is not "free", and the brief's
standard for S2 was free.

**F5 — The asymptotic rate is discontinuous at the trivial connection under the canonical clock**
(§6(e)) and continuous under the uniform one. I do not know which is physical, and the page does
not decide it.

**F6 — The record slots carry no gauge action and are not cells of K1.** They are additions to the
carrier, justified as the minimum a write needs, not derived from the complex. A stronger
construction would find them inside a carrier; THEOREM S3-0 says they cannot be found inside *this*
kind of carrier's invariant algebra, which is an argument that they must be added, not a derivation
of what is added.

**F7 — Disjointness is proved by one self-contained route and corroborated by one external one.**
The elementary mean-observable argument (§5.7) is displayed in full. The product-state criterion is
cited by author name, not by digest, and I did not verify its statement against a primary source.
Nothing on this page rests on it.

**F8 — I found a caveat missing from W-01 and record it against the register, not against
myself.** The canonical three-way split has dimension 3 only for `W_F != 1` and `W_C != 1`; it
collapses to 2 and then to 1 as the connection trivialises (§2.3). W-01 states the split's
derivability without this connection-dependence. This does not affect W-01's ruling, which is about
a generic connection, but the register row should carry it.

**F9 — This is an UNAUDITED BUILD.** Custody §4 pairs every build with a default-refute audit that
re-derives at bytes. No audit has been run. **Grade: unaudited build. Not adversarially-checked,
and certainly not independently-corroborated** — and per custody §4 it could not reach the latter
grade even after an audit, since builder and auditor share a model lineage.

**F11 — A digest in this artifact's own custody block was transcribed wrongly in the first draft,
and the pointer rule caught it.** §9 carried `6629ae5b...` for `CUSTODY_V001.md` where the true
digest is `6629ae3b...` — one character, in the one place a reader would trust without checking. It
was found by re-extracting every hex token from the finished page and matching it against a fresh
`shasum -a 256` run, which is now the procedure and not an afterthought. **All seven cited digests
verify at bytes; no other hex token on the page is a digest.** Recorded because custody §1 exists
for exactly this failure and a page that fixed it silently would be worth less than one that
reports it.

**F10 — Floating-point floor.** `|Omega_N|` underflows below `N ~ 750`; all values past that are
reported as cumulative logs and exponentiated, never as raw products. The `< 1e-300` entries in
§4.3 are honest underflow, and the `(1/N)log|Omega_N|` column is the quantity actually computed.

---

## 8. FOR THE REGISTER

Question, phrased as it would be re-asked:

> **W-02 — CAN A CROSSING BE BUILT ON K1 THAT ESCAPES RECURRENCE AND SATISFIES THE NINE?**
> **RULING: YES, IN THE LIMIT, AT THE COST OF ONE QUBIT PER CELL, AND CONDITIONAL ON THE CELL
> SCHEDULE.** Proof on this page. **Reopens if:** a schedule condition weaker than
> `sum (1 - |Z_{k_n}|) = infinity` is shown necessary or sufficient · a record slot smaller than
> `M_2(C)` is exhibited · THEOREM S3-0 is shown to fail on a complex with higher-rank fibres or a
> non-abelian structure group · or a lineage-independent lane fails to reproduce `lambda` or the
> monotonicity of `|Omega_N|`.

---

## 9. CUSTODY

Built under `CUSTODY_V001.md` digest `6629ae3b5aedbafda4e56b2a743b84a888949c59f988296056a2a7abeb20aa49`
(recomputed at §1, and re-verified against a fresh `shasum -a 256` of every cited file after drafting —
see FLAG F11). Pointer rule observed: every governing term carries a digest, a `file:line`, or
a definition in §1's index. Exhibit-not-assert observed: every number above was produced by the
computation described beside it, and the two prior-work figures it reproduces (`0.0247 at n = 42`,
`0.99994`) were recomputed rather than copied. Alpha not engaged (I8). Predecessor material cited by
digest, never copied, never foundational (I1, I5). CHOICE LEDGER, IMPORT AUDIT and FLAG BLOCK above.
Sealed on creation to `.sha256` and `.seal.sha256`.
