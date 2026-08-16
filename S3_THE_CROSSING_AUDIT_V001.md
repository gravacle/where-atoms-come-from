# S3 — THE CROSSING — ADVERSARIAL AUDIT — V001 — 2026-08-16

**Default verdict REFUTED.** Testimony carries no weight; every number below was recomputed from
the definitions rather than read off the build. Custody `CUSTODY_V001.md` §4.

**TARGET.** `S3_THE_CROSSING_V001.md`
sha256 `cbf1d79679ca2ecf3ee260e8a6467062e3a93260b2325cfc287c261e9b4469cb`
(recomputed at bytes in this directory; both sidecars `shasum -a 256 -c` **OK**).

---

## 0. THE VERDICT, STATED FIRST

**OVERALL: CONFIRMED-WITH-CORRECTIONS.** The default is not sustained.

The construction survives. Every computable claim on the build's page reproduces — the algebra
counts, the loop-algebra dimensions, the circuit-span collapse, the embedding checks, the frozen
carrier's `0.024654 at k = 42` and `0.999941 at k = 377`, the whole `Omega_N` table to six figures,
`lambda = -0.766802` at `N = 4000` and `-0.767026` at `N = 200000`, the torus mean `-0.767508`,
the exact `K1`-connection column, the deficiency sums to three decimals, and every one of the nine
demonstration tables. That is an unusually clean numerical record and it is stated first because
the corrections below are about **what the page says its numbers mean**, not about the numbers.

**But four statements on the page are false**, and one of them is false by the page's own criterion
three paragraphs earlier:

- **§3.2's** `M_5 (x) M_2^{(x)N}` **embeds unitally in** `M_5^{(x)(N+1)}` **for N >= 1** is FALSE for
  every `N >= 1`. §3.2 itself states the test — a unital \*-homomorphism `M_m -> M_n` exists iff
  `m | n` — and `5.2^N | 5^{N+1}` requires `2^N | 5^N`. **This is the pattern the S2 audit refuted
  its build for: a claim asserted on a structure the same page's own ruling excludes.**
- **§6(d) and §0 item 3's** *"`Phi_infinity = 0` for every connection with `(W_F, W_C) != (1,1)`"*
  is FALSE. I exhibit a family of non-trivial connections on which `|Omega_N| = 1` for every `N` —
  **no formation, ever** — and one member of that family uses `K1`'s own published ready state.
  The page's Theorem S3-2 carries the hypothesis that makes it true; the verdict block, §4.4's
  summary and §6(d) drop it. **This is a second conditionality on the escape and the FLAG BLOCK
  does not carry it** (F3 flags only the schedule).
- **§5.3's** `|Omega_N| <= e^{lambda N}`, offered as P-1's durability half, is FALSE for **63% of
  `N <= 200000`**, first at `N = 4`.
- **§5.4's** `sup_{k>=42} |Z_k| = 0.999941` is a window artefact stated as an equality. Over
  `k <= 200000` the supremum is `0.999999981` at `k = 106123`, and the true supremum is `1`. **The
  register already carries this exact defect class as a correction against the S2 build (COR-H).**

**MY INDEPENDENT Q3 RESULT is at §4.** In short: the escape is real, it is not the escape the page
advertises, and it is conditional in **two** ways, not one.

---

## 1. POINTERS AND DIGESTS

Recomputed in this directory at audit time with `shasum -a 256`. All seven digests the build cites
verify at bytes, and no other hex token on its page is a digest:

```
S3_THE_CROSSING_V001.md                     cbf1d79679ca2ecf3ee260e8a6467062e3a93260b2325cfc287c261e9b4469cb
S1_CARRIER_K1_V001.md                       3eb70375bfd0900e4dd56cae294fa31b3b6e19cf6634853501fab5ffcebd92ac
S2_FORMATION_CONDITION_ON_K1_V001.md        248ce856efaef157c68e818dde589d0200bbc1dd9fd9fc1fcc8cdc7bc88734d9
S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md  0bea11bd4b7764f65c8d44cc8812d43bdcc569d34f51d0c46d82851a2efbd0d5
REGISTER_V001.md                            95ecd65d1d1b4e0440d29a69e6ec1cb36ada8059d3e88948ff2177d26f62000b
CUSTODY_V001.md                             6629ae3b5aedbafda4e56b2a743b84a888949c59f988296056a2a7abeb20aa49
FOUNDING_DESIGN_V001.md                     ca25b79c76531d909e75fcb58163ad6456b10f086d59638f47224210d19b13e9
INHERITED_FROM_THE_PREDECESSOR.md           96c1d305c7a4eeceab8dc78f971832eb0c2e34980a72d1c12c5ece5cd4ad0079
```

**Pointer sweep — universe declared, hit counts reported.** Pattern
`(S1|FOUNDING_DESIGN|REGISTER|CUSTODY|INHERITED|S2…)…:[0-9]+` over the build's page, newlines
normalised: **12 hits, 10 distinct** (`FOUNDING_DESIGN:61–65` occurs three times), plus one bare
`` `:119` `` with no filename (§6(f)) that the pattern does not catch. **All eleven distinct
pointers were opened and all eleven land on content.** The build claims "all 10 file:line pointers
resolve"; the count is 11 once the filename-less one is included.

```
S1 :33–44   b1 = 1, the one surviving cycle                  RESOLVES
S1 :48–49   L_v = C, rank one, Hermitian                     RESOLVES
S1 :59–63   gauge at vertices, g_v = exp(i theta_v)          RESOLVES
S1 :66–70   6 - 4 = 2                                        RESOLVES
S1 :66–77   the two invariants                               RESOLVES
S1 :76–82   W_F, W_C gauge invariant                         RESOLVES
S1 :114–115 W_F = -1, W_C = -i                               RESOLVES
FOUNDING_DESIGN :61–65  the inductive-limit template  (x3)   RESOLVES
FOUNDING_DESIGN :77–84  honest scope                         RESOLVES
FOUNDING_DESIGN :117–118 S2/S3 external contact              RESOLVES
FOUNDING_DESIGN :119 (cited bare, as ":119")                 RESOLVES — but the pointer is
                                                             filename-less (COR-L)
```

Named rulings: *"audit COR-F"* resolves to `S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md:639`;
*"the audit's own §4.9 test point"* resolves to that file `:424` and `:440`, where
`f = 2.0, c = 1.1, p = (0.4, 0.15, 0.15, 0.15, 0.15)`, `min 0.024654 at n = 42`,
`max 0.999941` appear verbatim. Both CONFIRMED.

**Sidecars.** All nine `.md` files in the directory verify against their `.sha256`. Two sidecar
*content* conventions coexist in the directory (absolute path in five files, bare filename in four);
S3 followed the bare form, matching `REGISTER_V001.md` and both S2 artifacts. Not S3's introduction;
recorded so custody §2 can normalise additively.

**Terms defined on this page:** *record functional*, *deficiency*, *support class*, *the corrected
criterion* (§4.4), *edge-transport operator* (§2.4).

---

## 2. D1 — THE ALGEBRA — **CONFIRMED-WITH-CORRECTIONS**

### 2.1 The counts, recomputed

Group averaging of `X -> g X g*` over finite subgroups of `U(1)^5`, random complex `X`, and the
rank of the averaging projector acting on all 25 matrix units:

```
|| E(X) - diag(X) ||_F  over Z_2^5  (32 elements)   =  1.1061e-15
|| E(X) - diag(X) ||_F  over Z_3^5  (243 elements)  =  6.2687e-15
|| E(X) - diag(X) ||_F  over Z_5^5  (3125 elements) =  7.6159e-14     [my own third check]
rank of the averaging projector on M_5(C)          =  5              EXACT
dim_C M_5(C) = 25 ;  trivial characters (u=v) = 5 ;  distinct non-trivial = 20 ;  5 + 20 = 25
```

**CONFIRMED.** `dim = 25`, `A_inv = D_5 = C^5`, `dim = 5`, counts saturated. The 20 non-trivial
characters are pairwise distinct as characters of `U(1)^5` — checked, not assumed.

### 2.2 The loop algebra

Rank of the span of all words up to length 6 in `{I, M_dF, M_c, M_dF*, M_c*}`:

```
  connection                    W_F                 W_C            dim C*(M_dF,M_c)   build
  generic (f=2.0, c=1.1)   -0.4161+0.9093j     0.4536+0.8912j            3              3
  S1 published instance    -1.0000+0.0000j    -0.0000-1.0000j            3              3
  W_F = 1                   1.0000+0.0000j     0.4536+0.8912j            2              2
  trivial connection        1.0000+0.0000j     1.0000+0.0000j            1              1
```

**CONFIRMED**, including F8's caveat: the split is connection-dependent. The independent
re-derivation of `{v0}/{v1,v2}/{v3,v4}` from the transports rather than from incidence is real and
is a genuine second route to W-01's result.

### 2.3 The record slot

`R` non-abelian: correct. Pure states of `C^k` are its `k` characters, pairwise orthogonal, so
overlap `in {0,1}`. `dim R >= 4`: correct, and I strengthen the build's argument — every
finite-dimensional C\*-algebra is a direct sum of matrix algebras, so its dimension is a sum of
squares; `3 = 1+1+1` is the only decomposition of 3 and is abelian, and `4 = 4` or `1+1+1+1`, so
`M_2(C)` is the **unique** minimal non-abelian choice, not merely a minimal one. The build's
conclusion is right and under-argued in its favour.

The write-map hypotheses that make the floor binding are correctly stated and I verify the key one:
a unitary write with a reset carrier forces `<phi^F, phi^C> = Z` exactly, so orthogonal slot states
require `Z = 0`.

```
Z_1              =  0.259878772772 - 0.318757783880j      (build: identical to 12 places)
<phi^F, phi^C>   =  0.259878772772 - 0.318757783880j      difference  5.55e-17
Gram(domain) - Gram(codomain), Frobenius                  7.850e-17
```

**The unitary `V` exists — I built one independently** (matched Gram–Schmidt bases plus an arbitrary
unitary between the 8-dimensional complements):

```
|| V* V - I_10 ||_F                      =  1.888e-15
|| V( x (x) |0> ) - s (x) phi^F ||       =  2.459e-16
|| V( y (x) |0> ) - s (x) phi^C ||       =  3.640e-16
|| s ||                                  =  1.000000        carrier reset in both branches
extension-independence: a second extension gives  || V'( x (x) |0> ) - s (x) phi^F || = 2.625e-16
```

**CONFIRMED.**

### 2.4 THEOREM S3-0 — correct in the body, over-broad in the verdict

The theorem as stated at §2.4 is **correct and correctly scoped**: it says *"let `M_gamma` be the
loop transport of §2.3"*, and for those operators the proof is airtight — they are diagonal
unimodular, so `(M rho M*)_{vv} = rho_{vv}`, and §2.2's fixed algebra is the diagonal.

**§0 item 1 drops the scope** and asserts: *"Every loop transport on a rank-one Hermitian fibre
bundle with a `U(1)` connection is a **diagonal** unitary."* That is false. Define the
**edge-transport operator** `T` for the unfilled cycle `v0 -> v3 -> v4 -> v0` by moving each fibre
value one edge along the loop, `(Ts)(v3) = U_{e4} s(v0)`, `(Ts)(v4) = U_{e5} s(v3)`,
`(Ts)(v0) = U_{e6} s(v4)`, identity off the loop. With `a_4, a_5, a_6 = 0.7, 1.3, -0.4`:

```
|| T*T - I ||             =  0.00e+00        T is unitary
T diagonal?                  False
T^3                       =  diag( W_C, 1, 1, W_C, W_C )      with W_C = exp(i(a4+a5+a6))
                                                              = -0.029200 + 0.999574j
diag( T rho T* )          =  [0.15  0.15  0.15  0.40  0.15]
diag( rho )               =  [0.40  0.15  0.15  0.15  0.15]     NOT PRESERVED
```

`T` is a *bona fide* parallel transport around the same closed cycle of the same complex, its cube
is exactly W-01's scalar operator, and it is not diagonal — so a diagonal (gauge-invariant)
observable **does** separate branches built from it. The theorem is a theorem about W-01's
scalar-multiplication convention for `M_gamma`, not about "every loop transport".

This does not damage the construction — the branches S3 actually uses are W-01's — but §0 item 1's
*"proved in general, not for `K1` only"* and §2.4's *"disposes, in one line, of every construction
whose record algebra is the gauge-invariant algebra of a complex"* are broader than what was proved.
**COR-F.**

### 2.5 An undeclared premise

The step §2.4 → §2.5 → *"this is why a crossing is needed at all"* requires the premise **"the
record must be gauge-invariant."** That premise appears nowhere in the CHOICE LEDGER or the IMPORT
AUDIT. It is load-bearing and it is applied asymmetrically: on the carrier it annihilates every
candidate, and on the added slots it is **vacuous**, because those slots carry no gauge action at
all (the build's own F6). The conclusion "the record must be added" is therefore part theorem and
part bookkeeping convention. **COR-J.** (Not fatal: the underlying fact — `x` and `y` have equal
moduli at every vertex, so only off-diagonal observables can separate them — is true and is the
real content.)

**GRADE D1: CONFIRMED-WITH-CORRECTIONS** (COR-F, COR-J).

---

## 3. D2 — THE DIRECTED SYSTEM — **CONFIRMED-WITH-CORRECTIONS**

### 3.1 Does anything actually grow? Yes — and the trap is genuinely disarmed

The brief's specific risk is repeated circuits of one loop. Recomputed rank of
`span{ M_dF^n, M_c^n : n <= N }` inside `M_5(C)`:

```
  N =   1 : 3      N =   5 : 3      N =  25 : 3
  N =   2 : 3      N =  10 : 3      N = 100 : 3
  N =   3 : 3
  unital *-algebra generated by all powers and products (words <= 6) : 3
```

**CONFIRMED — constant at 3, forever.** And the reason is structural, not numerical: in the
block coordinates `(a_0, a_F, a_C)` of §2.3, `M_dF^n = (W_F^n, W_F^n, 1)` and
`M_c^n = (W_C^n, 1, W_C^n)` are block-constant, and powers of block-constant functions are
block-constant. Circuits supply time and no algebra. **The build did not take the trap.**

### 3.2 The system that was taken

```
unital     || iota(I_20) - I_40 ||               =  0.000e+00
*-hom      || iota(AB) - iota(A) iota(B) ||      =  0.000e+00
adjoints   || iota(A*) - iota(A)* ||             =  0.000e+00
isometric  | ||iota(A)||_op - ||A||_op |         =  1.776e-15
dim A_N = 25 . 4^N :  25, 100, 400, 1600, 6400, 25600, 102400
5 . 2^N | 5 . 2^{N+1} for N = 0..7                            True
```

**CONFIRMED.** `iota_N(X) = X (x) I_2` is a unital, multiplicative, \*-preserving, isometric,
injective embedding; `A_N = M_{5.2^N}`; the system is UHF; `dim A_N -> infinity`. The limit is
`UHF(5 . 2^infinity)` and is infinite-dimensional. **Something genuinely grows.**

The divisibility obstruction to cellular growth also reproduces exactly:

```
wedge of full K1 copies at the root : d_N = 5, 9, 13, 17    5|9 False   9|13 False   13|17 False
wedge of unfilled triangles         : d_N = 5, 7,  9, 11    5|7 False   7|9  False    9|11 False
```

### 3.3 COR-A — the false claim, and it is false by §3.2's own criterion

§3.2 closes its rejection of alternative (d) with:

> *"`M_5 (x) M_2^{(x)N}` embeds unitally in `M_5^{(x)(N+1)}` for `N >= 1`, so nothing is lost."*

Apply the criterion the same section states — *"a unital \*-homomorphism `M_m -> M_n` exists **iff
`m | n`**"*:

```
   N     m = 5.2^N     n = 5^{N+1}     m | n ?     equivalent test  2^N | 5^N
   1          10              25        False                       False
   2          20             125        False                       False
   3          40             625        False                       False
   4          80            3125        False                       False
   5         160           15625        False                       False
   6         320           78125        False                       False
```

`5 . 2^N | 5^{N+1}` reduces to `2^N | 5^N`, which fails for every `N >= 1` because
`gcd(2,5) = 1`. It also fails in the limit: `5 . 2^infinity` does not divide `5^infinity`, so
`UHF(5 . 2^infinity)` does not embed unitally in `UHF(5^infinity)` either. **The claim is false at
every reading, and it is the page's own test that refutes it.**

Severity: the sentence is a reassurance appended to a rejection that stands on minimality (§2.5's
floor of 4), so nothing downstream collapses. But it is a false theorem, stated in the section that
supplies the tool to refute it, and carried into CHOICE LEDGER C1(d). **COR-A.** The correct
statement is: *the two systems are incomparable — neither `UHF(5.2^infinity)` nor `UHF(5^infinity)`
embeds unitally in the other.*

### 3.4 A minor labelling defect

The build's §3.2 says extension by zero is *"the only linear \*-preserving map available."* There
are many linear \*-preserving maps between those algebras; what is meant, and what is true, is
*the only injective unital-on-its-image \*-homomorphism up to unitary equivalence*. Recorded at
COR-L.

**GRADE D2: CONFIRMED-WITH-CORRECTIONS** (COR-A, COR-L). The system is real, the maps are what the
build says, growth is genuine, and the trap was avoided by computation.

---

## 4. D3 — ESCAPE FROM RECURRENCE — **CONFIRMED-WITH-CORRECTIONS** (heaviest)

### 4.1 Everything numerical reproduces

Independent implementation, closed form checked against direct matrix action on `C^5` with random
vertex phases on the ready section:

```
closed form vs matrix action, max deviation over k <= 200   :  1.528e-14      (build 1.538e-14)
min |Z_k| over k <= 400    =  0.024654  at k =  42                            (build identical)
sup |Z_k| over k <= 4000   =  0.999941  at k = 377                            (build identical)
#{ |Z_k| > 0.99, k <= 4000 } = 37                                             (build 37)
```

The record functional, `k_n = n`, `f = 2.0`, `c = 1.1`, `p = (0.4, 0.15, 0.15, 0.15, 0.15)`:

```
     N     |Z_N|        |Omega_N|          (1/N)log|Omega_N|     build      verdict
     1   0.411271     4.112706e-01            -0.888504       -0.888504     MATCH
     2   0.470386     1.934559e-01            -0.821353       -0.821353     MATCH
     5   0.364118     1.772339e-02            -0.806574       -0.806574     MATCH
    10   0.776953     6.540411e-04            -0.733234       -0.733234     MATCH
    20   0.247337     2.514545e-07            -0.759800       -0.759800     MATCH
    42   0.024654     2.505486e-15            -0.800483       -0.800483     MATCH
    50   0.681600     1.571127e-18            -0.819895       -0.819895     MATCH
   100   0.350139     4.036647e-34            -0.768925       -0.768925     MATCH
   200   0.610201     1.405872e-67            -0.769663       -0.769663     MATCH
   400   0.600798    2.936168e-134            -0.768673       -0.768673     MATCH
  1000   0.580708        < 1e-300             -0.768205       -0.768205     MATCH
  2000   0.573418        < 1e-300             -0.768470       -0.768470     MATCH
  4000   0.360282        < 1e-300             -0.766802       -0.766802     MATCH

lambda at N = 200000                             = -0.767026    (build -0.767026)   MATCH
torus mean of log|Z| (independent 4000^2 grid)   = -0.767508    (build -0.767507)   MATCH
|Omega_N| strictly decreasing for all N <= 200000 : True   (max |Z_n| = 0.999999981)
```

Deficiency sums, `K1`'s exact published column, the near-recurrence records, and the gauge check all
reproduce:

```
sum_{n<=10} (1-|Z_n|) = 4.850   n<=100 : 47.240   n<=1000 : 469.347
n<=10000 : 4692.037             n<=100000 : 46918.264                    ALL MATCH
K1's connection (W_F=-1, W_C=-i, p=(1/2,0,0,1/4,1/4)) :
   Z_1 = 0, Z_2 = -1, Z_3 = 0, Z_4 = +1, Z_5 = 0, Z_6 = -1 ; matrix-action dev <= 9.97e-16
record near-recurrences k = 6, 63, 154, 377, 6723, 7100, 99023, 106123  ALL MATCH to 9 figures
|Omega_30| = 1.192688e-10  (build identical) ; 8 random gauge transformations, spread 2.585e-25
```

**On the numbers, the build is clean.** I found no arithmetic error anywhere in §4.

### 4.2 COR-B — the escape is conditional in TWO ways, and the page carries only one

The build's Theorem S3-2 is stated with the hypothesis *"a ready state with `p_0`, `p_1+p_2`,
`p_3+p_4 > 0`"*. **That hypothesis is doing real work, and three other places on the page drop it:**

- §0 item 3 (the verdict);
- §4.4's closing *"The escape is by unbounded cell count alone"*;
- **§6(d)**: *"`Phi_infinity = 0` for every connection with `(W_F, W_C) != (1,1)`"*.

Counterexamples, computed over `k <= 50000` with `k_n = n`:

```
  W_F = 1,   W_C = e^{i},     p = (1/2, 0, 0, 1/4, 1/4)   min|Z_k| = 1.000000000000   |Omega| = 1
  W_C = 1,   W_F = e^{i},     p = (1/2, 1/2, 0, 0, 0)     min|Z_k| = 1.000000000000   |Omega| = 1
  W_F W_C = 1 (f=1, c=-1),    p = (0, 1/2, 1/2, 0, 0)     min|Z_k| = 1.000000000000   |Omega| = 1
  any connection,             p = (1, 0, 0, 0, 0)         min|Z_k| = 1.000000000000   |Omega| = 1
```

All have `(W_F, W_C) != (1,1)`. **None ever forms, at any `N`, under any schedule.** And the first
row uses **`K1`'s own published ready state** — the very state §4.5 calls "the sharpest exhibit on
this page".

**This is a new restriction the build does not carry anywhere, including its FLAG BLOCK.** F3 flags
the *schedule*; nothing flags the *support of the ready state*. **COR-B.**

### 4.3 COR-C, COR-D, COR-E — three further false statements in §4–§5

**COR-C.** §5.3 gives P-1's durability half as `|Omega_N| <= e^{lambda N}`. Tested against the
page's own `lambda` values:

```
  lambda = -0.766802  (page's §4.3 value) :  1719 of N <= 200000 violate; first at N = 4
  lambda = -0.767026  (page's N=200000 value) : 127042 of 200000 violate — 63.5%
  witness N = 10 : |Omega_10| = 6.540411e-04   >   e^{-0.766802 x 10} = 4.675426e-04
```

`(1/N)log|Omega_N|` oscillates about `lambda` and exceeds it infinitely often, so no such uniform
bound exists. The **true** durability statement — monotonicity plus `(1/N)log|Omega_N| -> lambda` —
is on the page at §4.3–4.4 and is correct. The inequality is a false decoration on a true property.

**COR-D.** §0 item 3 states strict decrease and injectivity **unconditionally**:
*"strictly decreasing at every step, hence injective, hence returns to no earlier value ever."*
Theorem S3-3 carries the hypothesis `z_n < 1 for all n`. On **`K1`'s own published connection**
(§4.5) `z_2 = 1` exactly and `Omega_N = 0` for every `N >= 1` — constant, hence **not** strictly
decreasing and **not** injective. The page's headline and the page's showcase contradict each other.
(The property that matters, `sup_{M>=N}|Omega_M| = |Omega_N|`, survives; injectivity does not.)

**COR-E.** §5.4 contrasts the record with *"the frozen carrier, whose `sup_{k>=42} |Z_k| =
0.999941`."* Recomputed:

```
sup_{42 <= k <=   4000} |Z_k| = 0.999941230  at k =    377
sup_{42 <= k <= 200000} |Z_k| = 0.999999981  at k = 106123
true sup over k >= 42 = 1, not attained (the orbit closure of (u^k, v^k) contains the identity)
```

A windowed maximum stated as a supremum. **The register already carries this exact defect against
the S2 build** — *"the recurrence figure `0.994373` is a window artefact… label it a lower bound"*
(`REGISTER_V001.md` `95ecd6…`). The build repeated it while quoting the row that corrects it. Note
this cuts the build's own way: the true supremum is `1`, so the contrast with the record is
*stronger* than claimed, not weaker.

### 4.4 A CONSTRUCTION FROM THE AUDIT — the corrected criterion

Custody §6: a clean result is a construction. Theorem S3-2 covers one of four support classes. Here
is the statement that covers all of them.

Put `u = conj(W_F)`, `v = W_C`, and let `p` have three block weights `p_0` (root), `p_F = p_1+p_2`,
`p_C = p_3+p_4` with `chi_0 = uv`, `chi_F = u`, `chi_C = v`. Let `S = { a : p_a > 0 }` be the
**support class**, and let

```
G  =  the closed subgroup of the unit circle generated by  { chi_a / chi_b  :  a, b in S }.
```

> **THEOREM (corrected S3-2).** Along the canonical clock `k_n = n`:
> `|Z_n| = 1 for every n` **iff** `G = {1}`, in which case `|Omega_N| = 1` for all `N` and there is
> **no formation, ever**. If `G != {1}` then `sum_n (1 - |Z_n|) = infinity` and
> `(1/N)log|Omega_N| -> lambda = int_G log|Z| dmu_G < 0`, so `|Omega_N| -> 0` exponentially.
>
> **Proof.** `|Z_n| = 1` with weights summing to 1 forces the unit numbers `{chi_a^n : a in S}` to
> coincide, i.e. `(chi_a/chi_b)^n = 1` for all `a,b in S`, i.e. `chi_a/chi_b` lies in the `n`-th
> roots of unity for all `n`, i.e. each ratio is `1`. So `|Z_n| = 1` for all `n` iff `G = {1}`. If
> `G != {1}`, `(chi_a^n)` equidistributes on the closure `T` of the group it generates, `|Z| < 1`
> off a proper closed subgroup of `T`, hence on a set of positive Haar measure and positive
> density; the rest is the build's own S3-1 and S3-2 argument. ∎

Reading off the four classes:

```
  S = {0, F, C}   G = <u, v>      trivial iff  W_F = W_C = 1        (this is the build's S3-2)
  S = {0, C}      G = <u>         trivial iff  W_F = 1
  S = {0, F}      G = <v>         trivial iff  W_C = 1
  S = {F, C}      G = <u/v>       trivial iff  W_F . W_C = 1
  |S| = 1         G = {1}         NEVER forms   (recovers W-01's "the root can never fire")
```

**Verified, 10 cases, all four classes, `k <= 50000`:**

```
  case                          min |Z_k|         |Omega_50000|     prediction
  |S|=3  generic             0.023011609829      0.0000000000      G non-trivial   -> forms
  |S|=3  (W_F,W_C)=(1,1)     1.000000000000      1.0000000000      G = {1}         -> never
  S={0,C}, W_F = 1           1.000000000000      1.0000000000      G = {1}         -> never
  S={0,C}, W_F = e^{i}       0.000015072177      0.0000000000      G non-trivial   -> forms
  S={0,F}, W_C = 1           1.000000000000      1.0000000000      G = {1}         -> never
  S={0,F}, W_C = e^{i 0.3}   0.000006985336      0.0000000000      G non-trivial   -> forms
  S={F,C}, W_F W_C = 1       1.000000000000      1.0000000000      G = {1}         -> never
  S={F,C}, W_F W_C != 1      0.000043013526      0.0000000000      G non-trivial   -> forms
  |S|=1  root only           1.000000000000      1.0000000000      G = {1}         -> never
  |S|=1  face block only     1.000000000000      1.0000000000      G = {1}         -> never
```

10 of 10. This is the statement §6(d) should carry, and it makes W-01's "the root can never fire"
a corollary rather than a separate fact.

### 4.5 MY INDEPENDENT Q3 RESULT

> **THE LIMIT DOES ESCAPE THE INHERITED NO-GO. THE ESCAPE IS REAL, WEAKER THAN ADVERTISED, AND
> CONDITIONAL IN TWO WAYS RATHER THAN ONE.**
>
> **(a) The hypothesis is genuinely unmet — CONFIRMED.** `A_N = M_5 (x) M_2^{(x)N}` with
> `iota_N(X) = X (x) I_2` is a directed system of unital \*-homomorphisms (verified to `0.0` /
> `1.8e-15`), `dim A_N = 25 . 4^N -> infinity`, and `A_infinity = UHF(5 . 2^infinity)` is
> infinite-dimensional. The inherited theorem binds constructions whose durable record is a closed
> system with **finitely many** unitary degrees of freedom. This one is not in that class. That much
> is not rhetoric.
>
> **(b) But the mechanism is not "recurrence is gone" — it is "the record grew."** The carrier's
> recurrence is **entirely untouched**: `sup_{k >= 42} |Z_k| = 1` (§4.3, COR-E), and slot by slot
> the branch overlap `Z_{k_n}` recurs exactly as W-01 recorded. What does not return is the
> **product** `Omega_N = prod Z_{k_n}`, and its monotonicity is a **tautology of the tensor
> construction**: `|Omega_{N+1}| = |Omega_N| . z_{N+1}` with `z <= 1`, for any per-cell tensor
> system whatever, formation or no formation. THEOREM S3-3 therefore carries almost no information —
> it holds identically on the never-forming families of §4.2, where `|Omega_N| = 1` for all `N`.
> **All the content is in `sum (1 - z_n) = infinity`,** i.e. in Theorem S3-2 and its corrected form.
>
> **(c) So the honest statement of the escape is:** *given a non-trivial group `G` (§4.4) and a
> schedule whose deficiencies diverge, the record functional falls exponentially at rate
> `lambda = int_G log|Z| dmu_G` and never returns.* Both conditions are necessary. The page carries
> the schedule condition (F3) and does not carry the support condition (COR-B).
>
> **(d) "Closed" is a stretch, though not a lie.** §4.4's *"The system is **closed** — no
> environment, no bath, no open dynamics has been added"* describes an unbounded supply of fresh
> ancillas prepared in a fixed pure state, which is the standard model of a bath. The claim is
> defensible only in the specific sense that **nothing is traced out** and every state stays pure —
> and that saving distinction is not stated. What the page should say, and does say elsewhere (F4),
> is that the escape costs one qubit per cell.
>
> **(e) The no-go is answered as a statement about a monotone functional, not as a dynamical
> statement.** No one-parameter dynamics on `A_infinity` is exhibited whose recurrence could be
> tested; "time" is the direction of the directed system, which cannot decrease by construction.
> Under the canonical clock cell number equals carrier circuit count, which is the strongest form
> available here and is a fair answer to the design's S3 contact — but it is an answer at the level
> of *the record's index*, not of *a flow*.
>
> **(f) I found no escape the build missed, and no error in its escape.** The build took the right
> system (the trap is genuinely disarmed by computation), and every number in §4 reproduces. The
> corrections are to what it claims the system shows.

**GRADE D3: CONFIRMED-WITH-CORRECTIONS** (COR-B, COR-C, COR-D, COR-E).

---

## 5. D4 — THE NINE — **CONFIRMED-WITH-CORRECTIONS**

Every HOLDS was attacked. Test point as the build's.

### 5.1 P-8 — verified **as stated**, not as paraphrased

Requirement: `omega^{N+1}( iota_N(A) ) = omega^N(A)`. Random Hermitian `A in A_N` at each stage,
both branches, `iota_N(A) = A (x) I_2`:

```
N = 1 :  dim A_N =  100    max | omega^{N+1}(A (x) I) - omega^N(A) |  =  2.483e-16
N = 2 :  dim A_N =  400                                              =  4.965e-16
N = 3 :  dim A_N = 1600                                              =  3.545e-16
N = 4 :  dim A_N = 6400                                              =  4.628e-16
```

**HOLDS.** The identity is exact for product vectors along a tensor system, and the demonstration
is correct. CONFIRMED.

### 5.2 P-6 — verified **as stated**, with the exact inequality

Requirement: `|| [M_N, O] || <= 2m ||O|| / N` for `O in A_m`. With `M_N = (1/N) sum_{n<=N}
sigma_x^{(n)}` and my own random `O` (operator norm 1), computed as operator norms of the actual
commutators in `M_5 (x) M_2^{(x)N}`:

```
   m   N    ||[M_N,O]||     bound 2m||O||/N     ratio     holds
   1   2      0.696949         1.000000        0.6969     True
   1   4      0.418960         0.500000        0.8379     True
   1   6      0.248826         0.333333        0.7465     True
   1   8      0.161969         0.250000        0.6479     True
   2   2      1.131163         2.000000        0.5656     True
   2   4      0.576268         1.000000        0.5763     True
   2   6      0.411070         0.666667        0.6166     True
   2   8      0.267037         0.500000        0.5341     True
```

My ratios vary because I drew a fresh `O` at each row; for a **fixed** `O` the ratio is exactly
constant in `N`, since `[M_N,O] = (1/N) sum_{n<=m}[sigma_x^{(n)},O]`, which is why the build's
column reads `0.4329` and `0.2322` throughout. That constancy is a correct consequence, not a
coincidence. The proof is valid: `sigma_x^{(n)}` commutes with `A_m` for `n > m`, and each of the
`m` surviving commutators has norm at most `2||sigma_x|| ||O||`. **HOLDS.** CONFIRMED.

### 5.3 P-9 — the sharpest target: **produced, not asserted**

The build does exhibit orthogonal reduced supports rather than claiming them, and the exhibit is
correct. The reduced states on `R_N` are the pure states of the product vectors, with overlap
`Omega_N`; supports are the rank-one projections onto them; orthogonal iff `Omega_N = 0`; and

```
     N    |Omega_N|          || omega_F^N - omega_C^N || = 2 sqrt(1 - |Omega_N|^2)
     0   1.000000e+00              0.000000000000
     1   4.112706e-01              1.823026602854
     2   1.934559e-01              1.962217957798
     5   1.772339e-02              1.999685856739
    10   6.540411e-04              1.999999572230
    20   2.514545e-07              2.000000000000
    42   2.505486e-15              2.000000000000
   100   4.036647e-34              2.000000000000
```

reproduces to twelve places. The finite-stage failure is correctly identified as necessary (two pure
states of a finite-dimensional space with non-zero overlap cannot have orthogonal supports), and the
finite-stage **success** on `K1`'s own published connection at the first cell is exact and I
verified it in closed form (`Z_1 = 0` exactly). The limit statement rests on disjointness, i.e. on
P-7 — see the next item. **HOLDS**, with the P-7 caveat. CONFIRMED.

### 5.4 P-7 — COR-H: the "self-contained elementary proof" is a sketch, and it names an object that
does not exist

§5.7 argues: *"in the limit `M_infinity` is a **sharp central observable** with value `+s̄` on
`omega_F` and `-s̄` on `omega_C`."*

**There is no such element of `A_infinity`.** `M_N` has spectrum spread over `[-1,1]` at every `N`
(its extreme eigenvalues are `+1` and `-1`), so for every `N` and every scalar `c`,
`|| M_N - c.1 || >= 1`. The sequence `(M_N)` has no norm limit and converges to a scalar in no
algebra containing it. The correct object lives in the **GNS weak closures**, and getting it there
needs two steps the page does not write:

1. From `omega_F((M_N - s̄_N)^2) -> 0` one gets `(pi_F(M_N) - s̄)Omega_F -> 0` **only on the cyclic
   vector**. Extending to a dense set uses asymptotic centrality:
   `(pi_F(M_N) - s̄) pi_F(a) Omega_F = pi_F(a)(pi_F(M_N) - s̄)Omega_F + [pi_F(M_N), pi_F(a)]Omega_F`,
   both terms `-> 0` by P-4 and P-6; with `||M_N|| <= 1` this gives `pi_F(M_N) -> s̄ . 1` strongly.
2. If `omega_F` and `omega_C` were quasi-equivalent there would be a **normal** \*-isomorphism
   `theta` with `theta . pi_F = pi_C`; normality carries the strong limit across, giving
   `theta(s̄ . 1) = -s̄ . 1`, hence `s̄ = 0`, contradicting `s̄ = 0.785995 > 0`.

I supply both steps here, so **P-7 HOLDS** — but as written on the build's page the elementary route
is a sketch and the rigorous load falls on I7, the product-state criterion the page explicitly calls
*"cited but not load-bearing"*. That inversion is the defect. **COR-H.**

Numerically, `s̄` and `lambda` reproduce exactly across ready states (`N = 50000`):

```
  p = (0.4000, 0.3000, 0.3000)   s̄ = 0.785940   lambda = -0.767043      build identical
  p = (0.5000, 0.0000, 0.5000)   s̄ = 0.636630   lambda = -0.692761      build identical
  p = (0.3333, 0.3333, 0.3333)   s̄ = 0.790162   lambda = -0.776733      build identical
  p = (0.8000, 0.1000, 0.1000)   s̄ = 0.560484   lambda = -0.223149      build identical
  p = (0.1000, 0.4500, 0.4500)   s̄ = 0.730358   lambda = -0.727616      build identical
  p = (0.2000, 0.5000, 0.3000)   s̄ = 0.760015   lambda = -0.693152      build identical
```

The build's handling of the `s̄ ~ pi/4` near-coincidence — recording it, showing it is not stable
across the table, and declaring nothing uses it — is exactly right under custody §5 and is noted in
the build's favour.

### 5.5 P-1 — COR-C, already stated

Irreversibility in the limit is correct and its argument is sound: the GNS representation of
`omega . Ad(u)` is `Ad(pi_omega(u)) . pi_omega`, unitarily equivalent to `pi_omega`, hence
quasi-equivalent, hence not disjoint — so no unitary of `A_infinity` carries `omega_F` to `omega_C`.
I checked that this argument is self-contained on the page and it is.

Finite-stage reversibility is correct and trivially so: any two unit vectors of a
finite-dimensional space are related by a unitary of that space. The exhibited table is a
demonstration of a triviality, which the build says.

The **durability** half is stated as `|Omega_N| <= e^{lambda N}` and is false (COR-C). **HOLDS**
after replacing that inequality with the monotonicity and the Lyapunov limit that §4.3–4.4 actually
prove.

### 5.6 P-2, P-3, P-4, P-5 — all CONFIRMED

```
P-2  || [ V_(cell 2), sigma_x on slot 1 ] ||  = 0.000e+00 ; slot m state is stage-independent
P-3  eps      first N   |Omega_N|        sup_{M>=N}|Omega_M| == |Omega_N| ?   ceil(log eps/lambda)
     1e-1        3      7.613e-02                 True                              4
     1e-3        9      8.418e-04                 True                             10
     1e-6       20      2.515e-07                 True                             19
     1e-12      37      5.203e-13                 True                             37
     1e-30      90      3.368e-31                 True                             91
     1e-100    300      4.031e-101                True                            301
P-4  N=1 gap 1.823027 sd 0.411271 · N=200 gap 1.577279 sd 0.041121 · N=4000 gap 1.571990 sd 0.009219
     (identical to the build at every row; omega_F(M_N) = (1/N) sum sqrt(1-|Z_n|^2) verified)
P-5  every 2nd  |S|=2000  prod < 1e-300   pointer +0.786057
     every 10th |S|= 400  prod 1.931e-133 pointer +0.786590
     every 100th|S|=  40  prod 3.966e-13  pointer +0.796081        (build identical)
     random 5%  |S|= 200  prod 2.124e-67  pointer +0.792498        (my own draw)
     full 4000                            pointer +0.785995
```

P-3's non-return is exact, not asymptotic, and I confirmed `sup_{M>=N}|Omega_M| = |Omega_N|`
directly for every threshold row rather than inferring it. All four **HOLD**.

**GRADE D4: CONFIRMED-WITH-CORRECTIONS** (COR-C, COR-H). Nine of nine survive the attack in the
limit; the two finite-stage failures are correctly identified as necessary; two of the nine carry
defective supporting statements.

---

## 6. D5 — CONDITION TRANSPORT — **CONFIRMED-WITH-CORRECTIONS**

**(a) The multiplicative law, computed.** `Phi_{N+1} = Phi_N . Z_{N+1}` verified elementwise for
`N = 1..99`: **True**, to machine precision. `Phi_1 = Z_1` at `k_1 = 1`: verified to `5.55e-17`.

**(b) The compatibility, computed.** `(omega_F^{N+1} - omega_C^{N+1}) . iota_N = omega_F^N -
omega_C^N` is P-8 applied to each branch and subtracted; my §5.1 figures (`2.5e-16` … `4.6e-16`)
verify it. **CONFIRMED.**

**(c) Firing is absorbing.** `Phi_N = 0 => Phi_M = 0` for `M >= N` follows immediately from (a).
**CONFIRMED**, and it is the one place where the crossing delivers something W-01's condition alone
did not.

**(d) Criterion-to-rate — FALSE AS STATED (COR-B).** See §4.2 and the corrected criterion at §4.4.
The correct statement is `Phi_infinity = 0` iff `G != {1}`, not iff `(W_F,W_C) != (1,1)`.

**(e) The trivial-connection limit — CONFIRMED, every row.**

```
     (f, c)          |Z_1|      |Omega_10|     |Omega_100|      lambda_B      lambda_A
  (0.0000, 0.0000)  1.000000   1.0000e+00    1.0000e+00       +0.000000     +0.000000
  (0.0010, 0.0006)  1.000000   9.9992e-01    9.3556e-01       -0.760309     -0.000000
  (0.0100, 0.0060)  0.999980   9.9245e-01    1.0686e-03       -0.765479     -0.000020
  (0.0500, 0.0300)  0.999508   8.2631e-01    9.4865e-35       -0.765183     -0.000492
  (0.2000, 0.1000)  0.992964   4.8838e-02    8.2888e-31       -0.681963     -0.007061
  (1.0000, 0.6000)  0.813939   5.3873e-04    3.6970e-34       -0.765229     -0.205869
  (2.0000, 1.1000)  0.411271   6.5404e-04    4.0366e-34       -0.767026     -0.888504
```

Identical to the build at every entry. The design's S2 external contact survives: `|Omega_N| -> 1`
at each fixed `N` as the connection trivialises. FLAG F5's discontinuity of `lambda_B` is real and
correctly diagnosed — for any dense orbit the rate is the orbit-closure average however small
`(f,c)` is.

**(f) COR-K — two rows are not reproducible from their displayed parameters.**

```
  f = 3.14159, c = 1.57080  ->  lambda_B = -0.860699        page reports -0.804719
  f = 3.14159, c = 4.71239  ->  lambda_B = -0.831362        page reports -0.804719
  f = pi,      c = pi/2     ->  lambda_B = -0.804719        (orbit order 4)   MATCHES the page
  f = pi,      c = 3pi/2    ->  lambda_B = -0.804719        (orbit order 4)   MATCHES the page
```

The reported value is the exact order-4 orbit value, `(2 log 0.316228 + log 0.4 + log 1)/4 =
-0.804719`, so the computation behind it is right and the *displayed parameters* are five-decimal
truncations that do not produce it — at `N = 200000` the truncation drifts by `0.53` radians. The
page's own standard is *"every number… is reproducible from the definitions given here alone."*
A reader who types `3.14159` gets a different number. **COR-K.** (The rows labelled dense-in-`T^2`
and `f = c` reproduce exactly: `-0.767026` and `-1.203587`.) Winding periodicity reproduces:
`|Z_1(f,c)| = |Z_1(f+2pi, c-4pi)| = 0.411270593796`, difference `1.11e-16`.

**GRADE D5: CONFIRMED-WITH-CORRECTIONS** (COR-B, COR-K).

---

## 7. D6 — CUSTODY AND SMUGGLING — **CONFIRMED-WITH-CORRECTIONS**

**What is clean.** All seven cited digests verify at bytes; the sealed target verifies against both
sidecars; all eleven file:line pointers resolve and land on content; the two named rulings resolve;
CHOICE LEDGER, IMPORT AUDIT and FLAG BLOCK are present and substantive; W-01's figures are genuinely
recomputed rather than copied (I reproduced both from the definitions); predecessor material is
cited by digest and never copied. **Alpha is not engaged** — I swept the page for coupling slots,
inserted measured constants and target-driven selection and found none; `lambda` and `s̄` are free
functions of connection and ready state, fitted to nothing, and the one numerical near-coincidence
(`s̄ ~ pi/4`) is disclosed, shown unstable and explicitly quarantined. That is the correct handling
under custody §5. **F11's self-report of a one-character digest error is the behaviour custody §1
exists to produce and is noted in the build's favour.**

**COR-G — the pointer rule is breached on the specification itself.** The build grades nine
requirements whose **exact statements do not exist anywhere in the declared corpus**:

```
sweep: "P-6" | "P-9" | "asymptotic centrality" | "2m||O||" over *.md in this directory
  FOUNDING_DESIGN_V001.md:69-75  — the nine by NAME and in order; no labels, no formulas
  INHERITED_FROM_THE_PREDECESSOR.md:22-26 — "enumerated and quoted in the predecessor's register
                                            at Q-1179"; no path, no digest, no text
  S3_THE_CROSSING_V001.md        — 8 hits, all the build's own
  total hits elsewhere: 0
```

§5.2 calls `2m ||O|| / N` *"the inherited constant"* and §5.10 *"inherited bound"*; IMPORT AUDIT I1
gives the source as `INHERITED_FROM_THE_PREDECESSOR.md` digest `96c1d3…`. **That file contains no
such bound, no such equation, and no labels P-1…P-9.** The labels and order are reconstructable from
`FOUNDING_DESIGN_V001.md:69–75`, which the build does not cite for them; the **formulas** are
reconstructable from nothing on disk. Under custody §1 they entered a governing clause from a brief
that is not in the corpus — the predecessor's dominant failure, named in the custody file's own
"why". The right repair is cheap: the build proves both statements itself, so it should either point
I1 at `FOUNDING_DESIGN_V001.md:69–75` for the names and declare the formulas **authored here**, or
put Q-1179's text in the corpus with a digest.

**COR-I — imported operator-algebra results presented as though the page had produced them.**
IMPORT AUDIT I6 says direct limit, UHF, disjointness, quasi-equivalence and asymptotic centrality
are *"all defined on this page"*. The **definitions** are. The **theorems used** are not, and are
neither proved nor attributed:

```
  §3.4  "a UHF algebra ... is determined by its supernatural number"     — Glimm's theorem
  §3.4  "unital, simple, with a unique tracial state"                    — asserted
  §3.4  "R_infinity = UHF(2^infinity), the CAR algebra"                  — asserted
  §5.7  the product-state quasi-equivalence criterion                    — attributed (I7), good
```

None of the first three is load-bearing — the only property the argument needs is
infinite-dimensionality, which the page proves from `dim A_N = 25.4^N`. But the I6 row is inaccurate
as written, and this is exactly the risk the brief names for this stage.

**COR-J — the undeclared premise** (§2.5 above): *the record must be gauge-invariant*. Load-bearing
for the passage from THEOREM S3-0 to "a crossing is needed", absent from the CHOICE LEDGER, and
vacuous on the slots it licenses.

**One further structural observation, not a correction.** The write map `V` is defined **by its
required input–output pairs**; existence follows from the Gram equality and I verified it
independently. But `V_n` depends on the cell index, on the circuit count and on the connection
through `Z`, and **no fixed interaction generating the family `(V_n)` is exhibited anywhere.** The
object the design names as the thing to be built — *the durability map* — is therefore delivered as
a family of unitaries specified by what they must do, together with a proof that such unitaries
exist and that they cannot manufacture the record they carry. That is a real construction and the
page is honest that the write "relocates" rather than creates; but a reader should not take §3.3 as
having produced a dynamical law. F2 gestures at this; it deserves to be stated at the strength above.

**GRADE D6: CONFIRMED-WITH-CORRECTIONS** (COR-G, COR-I, COR-J, COR-L).

---

## 8. CORRECTIONS, IN SEVERITY ORDER

| # | Where | Correction |
|---|---|---|
| **COR-A** | §3.2, C1(d) | *"`M_5 (x) M_2^{(x)N}` embeds unitally in `M_5^{(x)(N+1)}` for `N >= 1`"* is **FALSE for every `N >= 1`**, by the criterion §3.2 states itself: `5.2^N | 5^{N+1}` requires `2^N | 5^N`. Also false in the limit. Replace with: *the two UHF systems are incomparable; neither embeds unitally in the other.* Nothing downstream depends on it; the rejection of (d) stands on minimality. |
| **COR-B** | §0 item 3, §4.4, §6(d) | *"`Phi_infinity = 0` for every connection with `(W_F,W_C) != (1,1)`"* is **FALSE**. Four exhibited families of non-trivial connections never form, one of them on `K1`'s own published ready state. Theorem S3-2's positivity hypothesis must be carried. **This is a second conditionality on the escape — the support of the ready state — and the FLAG BLOCK carries only the schedule (F3).** Replace with the corrected criterion of §4.4 above: formation iff `G != {1}`. |
| **COR-C** | §5.3 (P-1) | `|Omega_N| <= e^{lambda N}` is **false for 63.5% of `N <= 200000`**, first at `N = 4`; witness `N = 10`, `6.540411e-04 > 4.675426e-04`. Replace with the monotonicity and `(1/N)log|Omega_N| -> lambda` that §4.3–4.4 actually prove. |
| **COR-D** | §0 item 3 | Strict decrease and injectivity are stated **unconditionally** and are contradicted by the page's own §4.5: on `K1`'s published connection `z_2 = 1` and `Omega_N = 0` for all `N >= 1` — constant, not injective. Carry Theorem S3-3's hypothesis `z_n < 1`, and note that non-return survives in the form that matters. |
| **COR-E** | §5.4 | `sup_{k>=42}|Z_k| = 0.999941` is a **window artefact stated as an equality**; over `k <= 200000` it is `0.999999981` at `k = 106123` and the true supremum is `1`. Label it a lower bound. **The register already carries this defect class against the S2 build (COR-H) in the very row this page cites.** |
| **COR-F** | §0 item 1, §2.4 remark | *"Every loop transport … is a diagonal unitary"* is false; the edge-transport operator `T` of §2.4 above is a unitary transport around the same cycle, is not diagonal, does not preserve `diag(rho)`, and has `T^3` equal to W-01's scalar operator. THEOREM S3-0's body is correctly scoped; the verdict block and the "disposes in one line" remark are not. |
| **COR-G** | I1, §5.2, §5.10 | **Pointer-rule breach (custody §1).** The exact statements of the nine — `2m||O||/N`, `omega^{N+1}(iota_N(A)) = omega^N(A)`, the labels P-1…P-9 — appear in **no file in the declared corpus**; the cited digest `96c1d3…` contains none of them. Point I1 at `FOUNDING_DESIGN_V001.md:69–75` for the names and order, and declare the formulas **authored on the page** (the page proves both), or place Q-1179's text in the corpus with a digest. |
| **COR-H** | §5.7 (P-7) | The *"elementary, self-contained"* disjointness proof names `M_infinity`, *"a sharp central observable"*, which **is not an element of `A_infinity`**: `||M_N - c.1|| >= 1` for every `N` and every scalar `c`. Two steps are missing (strong convergence on `pi(A)Omega` from variance + asymptotic centrality; normality of the quasi-equivalence isomorphism). Supplied at §5.4 above, so **P-7 holds** — but as written the load falls on I7, which the page calls non-load-bearing. |
| **COR-I** | I6, §3.4 | Glimm's classification of UHF algebras by supernatural number, simplicity, unique trace, and `CAR ~ UHF(2^infinity)` are **imported theorems asserted without proof or attribution**, while I6 says they are "all defined on this page". Definitions are; theorems are not. None load-bearing — infinite-dimensionality is proved from `dim A_N`. |
| **COR-J** | §2.4→§2.5, CHOICE LEDGER | The premise **"the record must be gauge-invariant"** is undeclared, load-bearing for "a crossing is needed at all", and applied asymmetrically (fatal inside the complex, vacuous on the added slots, F6). Add it to the CHOICE LEDGER. |
| **COR-K** | §6(f) rows 3–4 | Not reproducible from the displayed parameters: `f = 3.14159, c = 1.57080` gives `-0.860699` and `f = 3.14159, c = 4.71239` gives `-0.831362`, not the reported `-0.804719` (the exact `(pi, pi/2)` / `(pi, 3pi/2)` order-4 value). Display the exact parameters or the computed values, not one of each. |
| **COR-L** | several | Minor, collected: (i) **F10**'s underflow point is wrong — `|Omega_N|` first underflows float64 at `N = 968`, not "`N ~ 750`" (`|Omega_750| = 2.5e-251`). (ii) §4.6's *"last 40 terms"* — there are **11** record-breakers to `k = 200000`; `1.533095` and `0.117102` are the totals over all 11 and reproduce exactly, but the label is wrong. (iii) §3.3's `||V - V'|| = 4.1688` **cannot be an operator norm** (unitaries differ by at most 2); the norm is unlabelled while §3.4 labels `||.||_op` explicitly. (iv) §5.3's table column headed `dim` is the **Hilbert-space** dimension `5.2^N` while §5.1's is the **algebra** dimension `25.4^N`. (v) §3.2's *"only linear \*-preserving map"* should read *only injective \*-homomorphism up to unitary equivalence*. (vi) §6(f)'s pointer `` `:119` `` carries no filename. |

---

## 9. GRADES

```
D1  the algebra (Q1)                CONFIRMED-WITH-CORRECTIONS     COR-F, COR-J
D2  the directed system (Q2)        CONFIRMED-WITH-CORRECTIONS     COR-A, COR-L
D3  escape from recurrence (Q3)     CONFIRMED-WITH-CORRECTIONS     COR-B, COR-C, COR-D, COR-E
D4  the nine (Q4)                   CONFIRMED-WITH-CORRECTIONS     COR-C, COR-H
D5  condition transport (Q5)        CONFIRMED-WITH-CORRECTIONS     COR-B, COR-K
D6  custody and smuggling           CONFIRMED-WITH-CORRECTIONS     COR-G, COR-I, COR-J, COR-L

OVERALL                             CONFIRMED-WITH-CORRECTIONS
```

**Why not REFUTED.** The default was tested and did not hold. The system is a genuine directed
system of unital \*-homomorphisms with an infinite-dimensional limit; the trap the brief named
(repeated circuits) was avoided by computation, not by assertion; the nine survive attack in the
limit with the two finite-stage failures correctly identified as necessary rather than as this
construction's weakness; and every numerical claim on 59 kB of page reproduced independently, most
of them to twelve significant figures. **Why not CONFIRMED.** Four statements are false, one of them
refuted by the page's own criterion three paragraphs earlier; a second conditionality on the escape
is unflagged; and the specification the page grades itself against has no resolvable pointer.

**GRADE OF THIS AUDIT: ADVERSARIALLY-CHECKED, NOT INDEPENDENTLY-CORROBORATED.** Custody §4 — build
and audit share a model lineage. A failure mode common to that lineage passes through both.

---

## 10. FOR THE REGISTER

> **W-02 — CAN A CROSSING BE BUILT ON `K1` THAT ESCAPES RECURRENCE AND SATISFIES THE NINE?**
> **RULING: YES, IN THE LIMIT, AT THE COST OF ONE QUBIT PER CELL, AND CONDITIONAL ON *TWO* THINGS —
> THE CELL SCHEDULE AND THE SUPPORT OF THE READY STATE.**
> `A_N = M_5(C) (x) M_2(C)^{(x)N}`, `iota_N(X) = X (x) I_2`, `A_infinity = UHF(5 . 2^infinity)`.
> The record functional `Omega_N = prod_n Z_{k_n}` is non-increasing, and falls exponentially at
> rate `lambda = int_G log|Z| dmu_G` **iff** the group `G` of §4.4 of this audit is non-trivial and
> the schedule's deficiencies diverge. Build `S3_THE_CROSSING_V001.md` sha256 `cbf1d7…`; audit this
> page. **Corrections carried: COR-A … COR-L above**, of which COR-A (a false unital-embedding
> claim, refuted by the build's own criterion) and COR-B (the missing ready-state condition) are
> material. **The escape is real; the recurrence of the one-cell quantity is untouched and its true
> supremum is `1`, not `0.999941`.**
> **REOPENS IF:** a lineage-independent lane fails to reproduce `lambda = -0.767026` or the
> monotonicity of `|Omega_N|` · the corrected criterion of §4.4 is shown to miss a support class ·
> a record slot smaller than `M_2(C)` is exhibited · THEOREM S3-0 is shown to fail for W-01's own
> transports on some complex · or a fixed interaction generating the write family `(V_n)` is
> exhibited, which would upgrade §3.3 from an existence proof to a dynamical law.

---

## CHOICE LEDGER

| # | Choice | Alternatives | Why | Status |
|---|---|---|---|---|
| A1 | Recompute every number from the definitions in a fresh implementation, sharing no code with the build | spot-check; trust the build's tables | custody §4: testimony carries no weight. The build's page is 59 kB of numbers and the only way to grade it is to produce them again | **closed** — every table in §2–§6 above is my own output |
| A2 | Attack the **hypotheses** of the build's theorems rather than its arithmetic, once the arithmetic reproduced | keep hunting for numerical error | after §4.1 reproduced to twelve figures the remaining risk was scope, not accuracy. That is where COR-A…COR-F were found | **closed** |
| A3 | Supply the missing steps in P-7's proof rather than fail it | grade P-7 REFUTED on the sketch | custody §6: a construction beats a hedge, and the argument is genuinely two lines from complete. The defect is recorded, the requirement holds | **closed** — §5.4 |
| A4 | State a **corrected criterion** covering all four support classes rather than only exhibit the counterexample | exhibit the counterexample and stop | the metric is constructions completed. The counterexample alone would leave the page's §6(d) unrepairable; the criterion repairs it | **closed** — §4.4, verified 10/10 |
| A5 | Grade CONFIRMED-WITH-CORRECTIONS, not REFUTED | REFUTED on COR-A + COR-B | the S2 precedent for REFUTED was a theorem whose *conclusion* failed on the space actually chosen. Here the false claims are (a) an appended reassurance whose section stands without it and (b) a dropped hypothesis that the page's own theorem statement carries. The construction itself survives every attack I could mount | **closed, with the reasoning displayed** |
| A6 | Use `f = 2.0, c = 1.1, p = (0.4, 0.15, 0.15, 0.15, 0.15)` for comparison and my own parameters for attack | the build's parameters only | comparability against the build and the S2 audit `:440` requires the first; independence requires the second. §4.2, §4.4 and §2.4 use parameters the build never touched | **closed** |

---

## IMPORT AUDIT

| # | Notion | Source | Defined here? | Survives without the import? |
|---|---|---|---|---|
| J1 | `K1`'s topology, gauge action, invariants, published instance | `S1_CARRIER_K1_V001.md` `3eb703…`, lines cited inline and opened | recomputed at §2.1–2.2 and §4.1 | yes |
| J2 | W-01's loop operators, formation condition, recurrence figures | `REGISTER_V001.md` `95ecd6…`; audit `0bea11…` `:424`, `:440`, `:639` | **recomputed from the definitions**, §4.1; `0.024654 at k=42` and `0.999941 at k=377` reproduced independently | yes |
| J3 | The nine, by name and order | `FOUNDING_DESIGN_V001.md:69–75` | used as the grading target only | the grading would have no target — see COR-G |
| J4 | The exact statements `2m||O||/N` and `omega^{N+1}(iota_N(A)) = omega^N(A)` | **no resolvable source in the corpus** — this is COR-G | verified as written on the build's page; both are proved there | n/a — the defect is the pointer, not the mathematics |
| J5 | Unital \*-homomorphism `M_m -> M_n` exists iff `m \| n` | standard representation theory of matrix algebras; the build states the one-line reason | used at §3.3 to refute the build's own §3.2 claim; the arithmetic `2^N \| 5^N` is elementary and displayed | yes |
| J6 | GNS, quasi-equivalence, disjointness, normal isomorphisms | standard C\*-algebra theory, named as such | used at §5.4 to **complete** the build's sketch; every step displayed | the completion would be unavailable; the defect COR-H would stand as a bare finding |
| J7 | Weyl equidistribution on a compact abelian group | standard; used in the corrected criterion §4.4 | the criterion's *"iff `G = {1}`"* half is elementary and self-contained; only the rate needs equidistribution | the iff survives; the rate formula would not |
| J8 | **Alpha** | — | **NOT ENGAGED.** Swept for coupling slots, inserted constants and target-driven selection: none found on the build's page or introduced here. The `s̄ ~ pi/4` near-coincidence is disclosed by the build, shown unstable, and used by nothing | n/a |

---

## FLAG BLOCK — DEFECTS IN MY OWN AUDIT

**G1 — Shared lineage, stated first.** Custody §4. I am the same model as the builder. Every finding
here is **adversarially-checked**, never **independently-corroborated**, and a failure mode common
to the lineage passes through both roles invisibly. COR-A is the kind of error that survives such a
pairing and it survived until an explicit divisibility table was written out; there may be others of
that kind that I did not think to tabulate.

**G2 — My reproductions are numerical, at one test point, in one language.** §4.1's agreement is
double-precision numpy against double-precision numpy. The exact column at §4.5 I verified in closed
form by hand as well as numerically; the rest I did not. A lane in different arithmetic could still
find a disagreement in the sixth figure of `lambda`.

**G3 — The corrected criterion of §4.4 is proved for the canonical clock only.** For a general
schedule `(k_n)` the "never forms" half still holds when `G = {1}`, but the exponential half needs
the schedule to equidistribute, which §4.6 of the build shows can fail. I did not state the general
schedule version and it is not obvious.

**G4 — I did not verify the build's §5.2 and §5.5 numbers, only its inequalities.** The specific
figures `0.432914`, `0.2322`, and the `0.000e+00` commutator depend on a random `O` and a specific
basis choice I cannot reconstruct; I verified the **bound** and the **1/N law** with my own operators
instead. If those particular figures were fabricated I would not have caught it — though the
constancy of the build's ratio column is the signature of a correct fixed-`O` computation.

**G5 — I graded P-9 HOLDS partly on P-7, which I had to repair.** If my repair at §5.4 is wrong,
P-7 and P-9's limit halves both fall, and the grade of D4 falls with them. The repair is standard
and I displayed every step, but it is mine and not the build's.

**G6 — COR-K may be a display convention rather than an error.** If the build computed with `np.pi`
and printed five decimals, the mathematics is right and only the page is unreproducible. I graded it
as a reproducibility defect on the page's own stated standard rather than as a numerical error, and
a reader may reasonably think that harsh.

**G7 — I did not open the predecessor's Q-1179.** It is outside this corpus and
`INHERITED_FROM_THE_PREDECESSOR.md` gives no path or digest for it. COR-G is therefore a statement
about the *declared corpus*, not a claim that the nine are misquoted. They may be quoted perfectly;
nothing on disk lets me check.

**G8 — This audit is itself unaudited**, and per custody §4 could not reach
independently-corroborated even if it were.

---

## 11. CUSTODY

Built under `CUSTODY_V001.md` digest
`6629ae3b5aedbafda4e56b2a743b84a888949c59f988296056a2a7abeb20aa49` (recomputed at §1). Pointer rule
observed: every governing term carries a digest, an opened `file:line`, or a definition on this page.
Exhibit-not-assert observed: every number above is my own output from the definitions, and where I
could not reproduce a figure I said which and why (G4). Sweep discipline observed: the universe for
the pointer sweep and the nine-requirements sweep is declared and the hit counts reported (§1, §7).
Predecessor material cited by digest, never copied, never foundational. Alpha not engaged (J8).
No git action taken. Sealed on creation to `.sha256` and `.seal.sha256`.
