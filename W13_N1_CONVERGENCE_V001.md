# W-13 — THE N1 CONVERGENCE RULING

**Register head W-12. Registrar's ruling, 2026-08-17. Lane: `LANE_W13_RULING/`, 11 files, `SEALS.sha256`.**

Four build lanes (Z, C, R, L) and four refutations were commissioned to settle whether

```
lambda  =  m(P),        P(x,y) = p00 + p10 x + p01 y + p11 xy,
lambda  =  lim_N (1/N) SUM_{k<=N} log|P(u^k, v^k)|,   u = conj(W_F), v = W_C
```

is true as stated and under exactly what hypothesis. This document rules. It is not a summary of
the eight; it re-derives the small number of facts on which they **disagree**, from code that
imports nothing from any of them, and it decides against three of them on the central question.

---

## 0. THE RULING, BEFORE THE EVIDENCE

> **N1 IS TRUE AT A NAMED CONNECTION, AND THE NAMED CONNECTION IS THE CORPUS'S OWN.**
>
> The sharp hypothesis is **H2 + (D1) + (D2)**, stated in §1. It is not merely an "almost every"
> statement. **(D1) and (D2) both follow from a published theorem — Baker's, on INHOMOGENEOUS
> linear forms in logarithms — whenever the two loop ANGLES `f` and `c` are algebraic real
> numbers with irrational ratio and the ready state has algebraic weights.** `f = 1.0,
> c = sqrt(2)` — W-10 N-4's *"the ONLY generic connection the corpus publishes"*, and the pair
> S4 used to verify its entire `lambda` column — is such a pair. **N1 holds there,
> unconditionally.**
>
> **Lane C's central ruling — *"N1 IS NOT VERIFIED AT ANY NAMED CONNECTION IN THIS CORPUS AND
> CANNOT BE WITH PRESENT MATHEMATICS"* — is FALSE. Lane L's §2.5 — *"NO effective Diophantine
> result covers `u = e^{-i}`, `v = e^{i sqrt2}` … Baker requires algebraic arguments"* — is
> FALSE, and its refuter confirmed the error rather than catching it. Lane R's R-7 is RIGHT,
> and its refuter was right to call it the round's genuinely new content.**
>
> The error in C and L is one inference: *u and v are transcendental, therefore Baker does not
> apply.* **The transcendence of `u` is not an obstruction — it is the hypothesis.** What enters
> the linear form is the **angle** `f`, not `e^{if}`; if `f` is a non-zero algebraic real then
> `u = e^{-if}` is necessarily transcendental **and** the constant term of the form is algebraic
> and non-zero, which is precisely the case Baker's inhomogeneous theorem covers.
>
> **N1 is nevertheless FALSE at four of the five connections this corpus has published**, all
> four for the already-registered reason that they are exactly resonant (W-10 N-4), and it is
> **false on a comeager set of connections** — which is Lebesgue-null. Measure and category
> disagree, and no counterexample can have algebraic angles.
>
> **The identification `lambda = m(P)` is NOT NEW.** For almost every starting phase it is
> Birkhoff plus Jensen, and in the form the corpus needs it is the standard first step of
> Herman's 1983 method for bounding Lyapunov exponents of one-frequency cocycles from below.
> This is the **third** time this program has re-derived a named result (after Hepp/Bell and
> after Cassaigne–Maillot). What is not in the literature is the **assembly at a single phase**,
> and the note must prove it rather than cite it.

---

## 1. THE SHARP HYPOTHESIS, IN PUBLISHABLE FORM

Write `alpha = arg(u)/2pi = -f/2pi`, `beta = arg(v)/2pi = c/2pi`, `theta = (alpha,beta) in T^2`,
`Z(P) = {P = 0} cap T^2`, `A_N = (1/N) SUM_{k<=N} log|Z_k|`.

### THEOREM W13-1 (SUFFICIENCY). The three hypotheses, and no more.

> Let `pi` be a probability vector and suppose `Z(P)` is **finite** (equivalently: `pi` is not on
> one of the three codimension-2 curve strata of lane Z's Z1). Suppose
>
> * **(H2)** `L(theta) := {(m,n) in Z^2 : m alpha + n beta in Z} = {0}`;
> * **(D1)** there are `c1, tau1 > 0` with `|| m alpha + n beta || >= c1 max(|m|,|n|)^{-tau1}` for
>   all `(m,n) != (0,0)`;
> * **(D2)** there are `c2, tau2 > 0` with `dist((u^k,v^k), Z(P)) >= c2 k^{-tau2}` for all `k >= 1`.
>
> Then `A_N -> m(P)`.

**Basis: PROVED, and the proof is on disk twice already.** It is `M1_08_THEOREMS.txt` T2(c)
(sealed W-08, 2026-08-16), re-derived independently as lane R's R-2 and lane C's C4. **Lane C's
version is the correct one and supersedes the other two**: M1_08 and R-2 both count orbit points
in a metric ball, which needs the *isotropic* discrepancy and costs a square root in dimension 2
(the L-refuter found this and is right); lane C's dyadic **sup-norm square annuli** are
differences of two axis-aligned boxes, so the ordinary star discrepancy suffices and no extra
import is needed. Adopt C4's form. `(D1)` feeds the discrepancy through Erdős–Turán–Koksma
(`D*_N = O(N^{-1/(tau1+1)} (log N)^2)`); `(D2)` truncates the shell sum at `j_max = O(log N)`.

**(D1) and (D2) are independent of each other and neither follows from H2.** (D2) is the one that
decides; it is *inhomogeneous* — a condition on the orbit relative to `Z(P)`, and `Z(P)` is
determined by the **ready state**. So **the hypothesis that decides N1 is not a condition on the
connection alone.** That is lane R's R-6, and priority belongs to M1_08 T2(c), which named the
condition as inhomogeneous relative to `Z(P)` one round earlier. `L-8` counted the consequence:
the word *"inhomogeneous"* occurs **0 times in `REGISTER_V001.md`**. Recounted by me: still 0.

### THEOREM W13-2 (THE LICENCE). Algebraic angles, irrational ratio.

> Let `pi` have **algebraic** entries with `Z(P)` finite, and let `f, c` be **algebraic real
> numbers with `f/c` irrational**. Then (H2), (D1) and (D2) all hold, and therefore
> `lambda = m(P)`.

**Proof.**

*(H2).* If `m f + n c = 2 pi q` with `q != 0` the left side is algebraic and the right is
transcendental, so `q = 0`; then `m f + n c = 0` forces `f/c in Q`. Hence for algebraic `f, c`,
**H2 is exactly `f/c irrational`** — a one-line criterion the corpus does not have.

*(D2).* `X01` (exact, in the field `Q(i sqrt5)` which happens to hold both readings of S4's SENSE
C): on `T^2`, `P = 0` forces `cos(arg x) = (p01^2+p11^2-p00^2-p10^2)/(2(p00 p10 - p01 p11))`, a
rational function of `pi`; so `x0` and `y0 = -A(x0)/B(x0)` are **algebraic**. At K1's registered
`pi = (0, 3/10, 3/10, 2/5)`: `cos s0 = -2/3` exactly, `x0 = (-2 + i sqrt5)/3`, minimal polynomial
`3z^2 + 4z + 3`, `y0 = conj(x0)`, and `P(x0,y0) = 0` **exactly in rational arithmetic**. Then the
first-coordinate distance is `|Lambda|/2pi` with

```
        Lambda  =  -k f  -+ s0  +  2 pi n            (s0 = arg x0)
      i Lambda  =  (-i k f)   -+  Log(x0)   +   n (2 pi i)
                   \______/       \______/        \______/
                    beta_0        log of an        log of 1
                   ALGEBRAIC      ALGEBRAIC        (a period)
```

`beta_0 = -i k f` is algebraic because **`f` is algebraic**, and non-zero for `k >= 1`, with
height polynomial in `k`. `Log(x0)` and `2 pi i` are **Q-linearly independent** because
`s0/2pi` is irrational — **Niven's theorem**, since `cos s0 = -2/3` is not in `{0, +-1/2, +-1}`.
Baker's inhomogeneous theorem then gives `|Lambda| > C H^{-K}` with `H` polynomial in `k`, i.e.
**(D2)**. The same form with `f` replaced by `c` covers the second coordinate; both zeros are
covered because `Z(P)` is finite; `dist >= (2/pi)|Lambda|/2pi` on the circle.

*(D1).* `|| m alpha + n beta || = |(-m f + n c) - 2 pi q| / 2pi`, and
`i Lambda = i(-m f + n c) - q(2 pi i)` — one logarithm, `beta_0` algebraic and non-zero, heights
polynomial in `max(|m|,|n|)` (and `|q| <= max(|m|,|n|)` since the form must be small). Same
theorem. **Sanity check on the SHAPE, and it is a check that could have failed:** at `n = 0` the
bound is `|m f - 2 pi q| > C H^{-K}`, which for `f = 1` is exactly the **effective irrationality
measure of `pi`** — a documented, effective, polynomial bound (Mahler 1953, exponent 42;
Salikhov 2008, 7.6063…). The form is of a shape for which effective bounds are known to exist. □

**Basis: PROVED, modulo one IMPORT quoted by shape and not read at the source** — Baker's
inhomogeneous bound (§4). Flagged here exactly as lane R flagged it (its D-4) and lane L flagged
Lawton (its custody note).

**Scope, stated so it cannot be over-read.** The licensed set is **countable** (algebraic angles),
so this is not a measure-theoretic improvement on §2. It is a **different kind** of statement:
it covers the connections one can *write down*. Every real connection anyone in this corpus has
ever named is either a rational multiple of `pi` in both angles (resonant, N1 false) or an
algebraic pair (licensed). There is no third named case on the page.

### THEOREM W13-3 (THE COUNTEREXAMPLE). H2 is not sufficient, and the failure set is comeager.

> There is a dense `G_delta` set `G` in `T^2` on which **H2 holds** and
> `liminf_N A_N = -infinity`, while `limsup_N A_N <= m(P)` always.

**Basis: PROVED.** Three independent proofs now exist (lane C's C6, the Z-refuter's R4, the
C-refuter's transposition of C6 from `theta` to the starting point `x_0`), and one explicit
ladder (lane R's R-3, rebuilt at 900 digits by its refuter). I verified the **density step**
myself, which is the only computable ingredient, at the corpus's own baseline (`X03` leg C, four
decades of `k`, perturbation rebuilt at each `k` from the fixed baseline `f=1, c=sqrt2`, moved
along `(+d,+d)` so `(uv)^k != 1` and the exact-hit relation Theorem Z4 forbids is never formed):

```
       k     |theta - base|   dist(k theta, Z)      log|Z_k|      log|Z_k|/k
     100          0.002816          2.048e-131      -300.0000       -3.000000
    1000         0.0004789         5.204e-1304     -3000.0000       -3.000000
   10000          8.443e-6        5.825e-13030    -30000.0000       -3.000000
  100000          1.396e-6         1.8e-130289   -300000.0000       -3.000000
```

The perturbation needed falls like `1/(2k)` while the depth is held at exactly 3 nats per unit
`k`. The local expansion carrying it was validated against direct high-precision evaluation
(`|alpha+beta| = 0.4` exactly; relative deviation `1.64e-42` at `d = 1e-24`).

**And a structural consequence nobody stated.** Theorem W13-2 says every algebraic-angle pair
with irrational ratio satisfies (D2). Therefore **`G` contains no algebraic-angle connection.**
Every counterexample must have a transcendental loop angle with prescribed approximation
properties — *constructed*, never merely *named*. That is why the corpus's own published pair is
safe and why the corpus's numerics could never have found the failure.

### THE THREE-WAY STATEMENT OF RECORD

| the connection | verdict on `lambda = m(P)` | authority |
|---|---|---|
| **algebraic angles `f, c`, `f/c` irrational** (incl. `f=1, c=sqrt2`) | **TRUE** | THEOREM W13-2 (Baker) |
| **Lebesgue-almost every connection** | **TRUE** | Borel–Cantelli (lane C's C5; lane L's 2.4) |
| **a comeager set of connections** | **FALSE**, `liminf = -infinity` | THEOREM W13-3 (Baire) |
| **both angles rational multiples of `pi`** (all four resonant published pairs) | **FALSE**, limit is a subtorus or finite-orbit average | erratum against W-02; W-10 N-4 |
| **`pi` on a curve stratum** (S1's own registered ready state) | **OPEN in general; the object is a one-variable Sudler product and the classical literature applies** | §5 |

---

## 2. A.E. CONNECTION VERSUS NAMED CONNECTION — AND WHICH N1 NEEDS

**N1 needs the named statement.** Lane C is right that the a.e. statement is available and clean,
and right that it transfers to every carrier and designation through W-12's Haar-pushforward
theorem. It is **wrong** that this is all N1 needs and all it can have.

The reason is in the register already. **W-12 Corollary 1**: nothing structural constrains
`(W_F, W_C)`; the connection is a free point of `T^2`. **W-14**: the functional sees exactly
three things — `pi`, `u`, `v`. So a statement about almost every connection is a statement about
a *randomly drawn* connection, and this construction has no measure on connections and no
mechanism that draws one. Every number this corpus has ever computed was computed at a
**named** connection, and **S4 verified its entire `lambda` column at exactly one of them**
(`S4:603`, `f = 1.0, c = sqrt(2)`, worst deviation `3.0e-06` across nine carriers). If N1 is
published as an a.e. statement only, the corpus's own verification table is outside its scope.
That is lane L's §2.4 read correctly — and lane L then denied the corpus the one thing that fixes
it.

**The measurement, six decades, one variable (the arithmetic type of `theta`), at K1's registered
`pi`** (`X03` leg A/A2; phases advanced in exact integers mod 2^64, surrogate-rational error at
`N = 1e8` bounded by `2.7e-12` turns and printed):

```
connection                                        1e3        1e4        1e5        1e6        1e7        1e8
f=1, c=sqrt2   ALGEBRAIC ANGLES  [S4:603]    +2.02e-03  +4.54e-06  -3.87e-05  -3.27e-06  -3.53e-07  +7.07e-09
f=sqrt3, c=sqrt5   ALGEBRAIC ANGLES          +2.35e-04  +3.76e-04  -6.76e-06  +4.05e-07  -2.72e-07
f=2^(1/3), c=5^(1/5)  ALGEBRAIC ANGLES       -1.22e-03  +1.80e-04  -4.97e-05  +3.14e-06  +1.71e-07
f=1, c=e       MIXED                         -5.49e-03  -6.03e-04  -5.68e-05  +9.62e-07  +1.82e-07
Haar draw (seed 20260817)                    -4.93e-04  -2.72e-04  +1.04e-05  +1.23e-06  -2.43e-08
f=2.0, c=1.1   RESONANT CONTROL              -6.97e-04  +4.73e-04  +4.78e-04  +4.92e-04  +4.93e-04  +4.93e-04
f=pi, c=pi/2   ORDER-4 CONTROL               -3.72e-02  -3.72e-02  -3.72e-02  -3.72e-02  -3.72e-02  -3.72e-02
```

Deviations from `m(P) = -0.767507880357785` (Jensen midpoint at `2^24`, tabulated over five node
decades; Cassaigne–Maillot closed form agrees to `2.6e-12`). **The two controls could have failed
and did not**: the resonant arm sits flat at the erratum-against-W-02 subtorus value
`-0.767015004` (register: `-0.767014993`), and the order-4 arm sits at `-0.804718956` = `-(1/2)log 5`
exactly. **The convergent arms are not evidence for the theorem** — lane C's C7 is right that no
simulation of any size can see the failure, so agreement is worth nothing on its own. They are
here because a *disagreement* would have refuted Theorem W13-2, and none appeared over six
decades.

**Lane C's recommended register text — *"never 'for the connection (f,c) = …', for any named pair
whatsoever"* — is withdrawn.** The replacement is §8.

---

## 3. WHAT DOES *NOT* CHANGE

* **`lambda = m(P)` is FALSE at four of the five published connections.** All four are exactly
  resonant (W-10 N-4: every rational `(f,c)` is). Gaps: `3.7e-02` at S1's own connection,
  `4.9e-04` at S3/S4's headline. Already of record; reconfirmed here at a sixth decade.
* **N1 is false at S1's published connection *with* S1's published ready state** for a second
  reason: `Z_k = 0` at every odd `k`, so `A_N = -infinity` for every `N` (M1_08 T2(e)).
* **Durability is untouched.** M1_08 T4 (`|Omega_N| -> 0` iff `G != {1}`) runs on Weyl applied to
  the **continuous** function `1 - |P|` with no Diophantine input. Nothing in this round bears on
  it. A reader who takes Theorem W13-3 as damaging durability has misread it.
* **The EMPTY stratum needs nothing.** Where `Z(P)` is empty, `log|P|` is continuous and bounded,
  Weyl alone gives convergence under H2, and `m(P) = log` of the dominant weight (lane Z's Z2/Z4(i),
  which is the elementary branch of Cassaigne–Maillot). Correct, and it is **three quarters of the
  state simplex** — I re-derived the fraction by exact lattice counting rather than Monte Carlo:
  `0.285337 / 0.268194 / 0.259234` at denominators `60 / 120 / 240`, converging to `1/4` singular
  from above.

---

## 4. THE CITATION

**Boyd–Lawton does not license N1, and that is now established twice over.** Lawton's theorem
(Lawton, *J. Number Theory* **16** (1983) 356–362, proving Boyd's conjecture) concerns
`m_1(P^{(r)})` — Haar averages of `log|P|` over **connected one-parameter subgroups**, i.e.
integrals. N1's `lambda` is a **discrete Birkhoff sum along a single Z-orbit**. The field draws
this distinction unprompted: Lind–Schmidt–Verbitskiy, *Entropy and growth rate of periodic points
of algebraic Z^d-actions*, Contemp. Math. **532** (2010), §9 — Lawton's result applies to
sequences of **compact connected** subgroups *"and so the diophantine issues disappear"*. Lane L
verified that quotation at the bytes; its refuter verified it a second time. **W-08's M1_08 T2(d)
is correct and is now correct on the record of the field.**

### WHAT DOES LICENSE IT

> **A. BAKER, "Linear forms in the logarithms of algebraic numbers. III", *Mathematika* **14**
> (1967) 220–228** — the **inhomogeneous** case: for algebraic `beta_0 != 0` and algebraic
> `alpha_j` whose logarithms are Q-linearly independent,
> `|beta_0 + beta_1 log alpha_1 + … + beta_n log alpha_n| > C H^{-K}` with `H` the maximum height
> of the `beta_j` and `C, K` effectively computable. **Refined with usable constants in Baker IV,
> *Mathematika* **15** (1968) 204–216.**
>
> **This is the citation the corpus is missing, and it is not the citation any of the four lanes
> named.** Lane L named Baker but for the wrong reason and then declared him inapplicable; lane R
> named him correctly by theorem-shape and did not identify which paper. **CUSTODY FLAG: I did not
> read Baker III or IV.** Both statements above are second-hand (encyclopaedic sources read
> 2026-08-17). The note must either read them or state the dependence.

**Supporting, each in its exact role:**

* **Erdős–Turán–Koksma** (Kuipers–Niederreiter, *Uniform Distribution of Sequences*, Ch. 2):
  (D1) `->` polynomial star-discrepancy decay. *Not read at the source by me either.*
* **Lawton 1983 Theorem 1 / E. Dobrowolski, "A note on Lawton's theorem", *Canad. Math. Bull.*
  **60** (2017) 484–489, Thm 1.4** — the small-value bound `mu_l{|P| <= v} <= C (v/h)^{1/sum(k_i-1)}`,
  **uniform in `pi`**. Lane L is right that this exists and right that Lawton must be cited twice.
  **But its refuter is right that it is OPTIONAL here**: Theorem W13-1 fixes `pi` before any limit,
  so a `pi`-dependent local constant is admissible, and lane Z's exact local expansion is sharper
  by four to nine orders at the shell scales the proof visits. Cite it as the route that makes the
  argument **uniform over the simplex**, which is a real gain and not the one lane L claimed.
* **Lind–Schmidt–Verbitskiy 2010** — the structurally identical problem, the naming of the
  difficulty as Diophantine, and the sentence *"such estimates … do not appear to be available"*.
* **V. Dimitrov, "Convergence to the Mahler measure and the distribution of periodic points for
  algebraic Noetherian Z^d-actions", arXiv:1611.04664** — **verified independently by me at the
  abstract**: it closes LSV's open problem for averages over torsion points, for all non-zero
  integer Laurent polynomials. **It does not touch a single Z-orbit at an infinite-order point**,
  so lane L's headline survives — but lane L's evidence block, which rests on LSV's 2010 *"it is
  not known"*, is **nine years stale**, and the paper was in lane L's own fetch queue. The
  L-refuter found this; it is the single best finding in the four refutations and it is
  **SUSTAINED**.
* **Cassaigne–Maillot, "Hauteur des hypersurfaces et fonctions zêta d'Igusa", *J. Number Theory*
  **83** (2000) 226–255** — the closed form for `m(a+bx+cy)`; the elementary branch **is** lane Z's
  Theorem Z2. Lane L's citation is correct at **2000**; **the L-refuter's "(1997)" is wrong** and
  I correct it here. The register cites Cassaigne eighteen times through S4 and **no register row
  carries the reference.**
* **Boyd is cited half.** The Boyd height `mu(r)` comes from Boyd, *J. Number Theory* **13** (1981)
  116–121, not from Boyd, *Canad. Math. Bull.* **24** (1981) 453–469. Both entries are needed.

### AND THE PLAIN STATEMENT

**No published theorem states Theorem W13-2 or Theorem W13-1. The note must prove them.** The
proof is three quoted theorems plus about one page of dyadic shell counting; it is not new
mathematics and must not be presented as such. This is the same verdict lane L reached, reached
for a different reason, and with the opposite consequence for the corpus's own connection.

---

## 5. NOVELTY — RULED, AND IT GOES AGAINST THE CORPUS FOR THE THIRD TIME

**"The rate is a logarithmic Mahler measure" is not a new identification.** It is the normal state
of affairs in this corner of mathematical physics: dimer free energy (Kenyon–Okounkov–Sheffield
2006), spanning-tree entropy (Lyons 2005), entropy of algebraic `Z^d`-actions (Lind–Schmidt–Ward,
*Invent. Math.* **101** (1990) 593–629). Lane L established this and it stands.

**And it is worse than lane L said, in the direction that matters.** `Omega_N = prod_{k<=N} Z_k`
is the transfer product of a scalar (`GL(1,C)`) quasi-periodic cocycle over the rotation of `T^2`
by `theta`, and `lambda` is its **Lyapunov exponent**. For almost every starting phase, the
Lyapunov exponent of such a cocycle equals `INT log|P| = m(P)` by Birkhoff plus Jensen — and using
**Jensen's formula and the Mahler measure to evaluate the Lyapunov exponent of a one-frequency
cocycle is the opening move of M. R. Herman, "Une méthode pour minorer les exposants de Lyapounov
…", Comment. Math. Helv. 58 (1983) 453–502.** *That is forty-three years old, it is cited in every
survey of one-frequency Schrödinger operators, and the word "Lyapunov" occurs twice in this entire
sealed corpus, in no register row and no lane directory* (lane L's L-12, recounted by its refuter,
recounted again by me).

**The d = 1 shadow is a named 65-year-old problem.** At S1's **own registered ready state**
`pi = (0,0,1/2,1/2)` — a CURVE state — `P = (1/2) y (1+x)`, so `|Z_k| = |cos(pi k alpha)|`
exactly and `SUM log|Z_k|` is a **Sudler product** (Sudler 1964; Erdős–Szekeres 1959; Lubinsky,
*J. Number Theory* **76** (1999) 217–247; Aistleitner–Technau–Zafeiropoulos). Measured, five
decades, `alpha = sqrt2 - 1`: `-0.68164 / -0.69196 / -0.69182 / -0.69317 / -0.69314` against
`m(P) = -log 2 = -0.69315`. **The L-refuter's correction to the identity is right and I confirm it
exactly**: the corpus's object is `prod|1 + q^k|`, not Sudler's `prod|1 - q^k|`, and the exact
relation is `prod_{k<=N}|1+z^k| = P_N(2a)/P_N(a)` — verified to `5.8e-12` at `N = 1e5`. So the
corpus's *other* registered state is a classical object about which the classical literature has
theorems, and no lane looked them up.

> **RULING ON NOVELTY: the identification is NOT new; the note must say so in its first
> paragraph.** What is not in the literature, as far as three searches by three parties can tell,
> is (i) this particular branch-comparison rate being *this* Mahler measure, and (ii) the
> **single-phase** assembly of Theorem W13-2. Both are modest. **A null literature search in this
> program has a bad record** — W-04 found Bell, W-06 found the corrected bibliography was *worse*,
> lane L found Cassaigne, the L-refuter found Dimitrov, and I found Herman. Score the novelty claim
> as **weak and leaning negative**, and write the note so that being scooped costs it nothing.

---

## 6. THE FOUR REFUTATIONS, READ AS SCEPTICALLY AS THE LANES

All eight documents' seals verify; I re-ran every `SEALS.sha256` in the eight directories.

### SUSTAINED (checked by me, from my own code or at the bytes)

* **The C-refuter's FATAL isolation finding.** `W13C_04_spacing_and_start.py:169` reads
  `if (s1, s2) == (0, 0) and k0 < K1 <= k0+n: lz[K1-k0-1] = logpv`. **Confirmed at the bytes.**
  The mpmath dip value is injected into exactly one arm, selected by a literal test on that arm's
  own label, at a `theta` engineered so the dip exists at `x_0 = 0` and nowhere else. Without the
  splice that arm reads `log(1.06e-15)/1e4 + m(P) = -0.770956`, i.e. the same as the other four to
  `3e-3`, not `-4.767319`. **C8's exhibit is a control that could not have failed, and C8's
  headline — *"the tied starting point … is THE WHOLE OF THE FAILURE"* — is false**, because the
  same Baire argument transposed from `theta` to `x_0` gives a comeager set of starting points at
  **every** H2 connection. **What survives of C8 is the narrow and correct point that Birkhoff's
  theorem gives a.e. STARTING POINT, `{0}` is null, and any citation of Birkhoff for N1 is a
  category error.** That point is worth keeping.
* **The C-refuter's misnaming finding.** The C-lane's `theta` is a rational plus `10^{-M} sqrt2`
  in each coordinate — a **quadratic irrational**, badly approximable, the *opposite* of Liouville.
  I rebuilt the same family (`X03` leg B) and confirm the depths `-1150.7 / -4604.6 / -18420.5 /
  -39999.9` nats and the arithmetic identity `A_{k1} - m(P) = log|Z_{k1}|/k_1`. **It is one dip
  whose depth is free, not a divergence**, and four rows of an arithmetic identity are not
  evidence. This is the **seventh** misnaming in this program.
* **The Z-refuter's B0b transcription finding.** `S4:575` writes `{00:4, 01:1, 10:2, 11:2}`, i.e.
  `pi = (4/9, 2/9, 1/9, 2/9)`; lane Z coded `(4/9, 2/9, 2/9, 1/9)`. **Confirmed, with the
  consequence measured**: the transposed state **factors** (`9P = (2+x)(2+y)`), so W-10 N-3's own
  non-factoring certificate `p00 p11 = 8/81 != 2/81 = p10 p01` is **false of it** (both
  cross-products are `4/81`). Blind to the transposition: stratum (EMPTY both), `min|P|`
  (`1/9` both), `m(P)` (`-0.810930216216` both) and every Jensen-derived figure. **Not blind:
  `|P|(x,y)` itself, max difference `0.222222` over `2e5` torus points** — hence `Z_k`, hence
  lane Z's only EMPTY-stratum convergence exhibit.
* **The R-refuter's ND-1, and it is the best mathematical finding of the four.** A rank-1 relation
  lattice need not be primitively generated. With `L = Z.(d,d)` the closure is
  `{(x,y) : (xy)^d = 1}`, which **contains** `{xy = 1}` — and at K1's registered `pi` the two zeros
  lie exactly on `{xy = 1}` (X01, exact). **So the closure meets `Z(P)` for every `d`, not only
  `d = 1`.** Reproduced from my own code, one variable moved (`d`), winding ratio held at `(1,1)`:
  `d = 1: -1.203973305`, `d = 2: -0.891400674`, `d = 3: -0.823368270`, `d = 5: -0.787137397`,
  against R-5's prediction `log(0.3) = -1.203972804` on all four. **R-5's "the ONE singular line
  is `c = f`" is wrong; the singular resonant set is `{(c-f)/2pi rational}`, a countable dense
  family.** N1's scope is unchanged — the set is still Haar-null — but the *description* of where
  N1 fails must be corrected before publication.
* **The L-refuter's Dimitrov finding** (§4) and its **Sudler ratio identity** (§5).
* **The Z-refuter's labelling finding (R-NEW-4), and it is sharper than it was stated.**
  `S4:566` fixes SENSE C only as *"(0.4,0.3,0.3) for 3 classes"*. Exact, all three readings
  (`X01`): the registered `(0,3/10,3/10,2/5)` gives `cos s0 = -2/3` and `x0 y0 = 1`; the reading
  `(0,2/5,3/10,3/10)` gives `cos s0 = -1/9`, minimal polynomial `9z^2+2z+9`, and
  `x0 y0 = (22 - 7 i sqrt5)/27 != 1`. **`m(P)` is identical in all three** (multiset), so Theorem
  Z4's anti-diagonal, and with it lane C's collapse of the 2-D inhomogeneous problem to the 1-D
  irrationality of `z_1`, are properties of a **labelling**, not of the registered multiset.
  **Note what this does NOT touch: in all three readings `x0` and `y0` are algebraic of degree 2,
  so Theorem W13-2 is labelling-robust.** The licence does not depend on the ambiguity; the
  lanes' collapse arguments do.

### OVERTURNED

* **The L-refuter's confirmation of L-4** (*"N1 is unlicensed at f = 1.0, c = sqrt(2)… survives
  even after my RL4(C) reduction, because `arg(uv)/2pi = (sqrt2-1)/2pi` is still a linear form in
  1, sqrt2 and pi"*). **A linear form in `1`, `sqrt2` and `pi` with integer coefficients is exactly
  what Baker's inhomogeneous theorem bounds** — `beta_0 = i(-m + n sqrt2)` algebraic and non-zero,
  one logarithm `2 pi i`. The refuter reached the right object and drew the opposite conclusion
  from it. **Both lane L and its refuter are overturned on this point.**
* **Lane C's C9 ruling** (*"cannot be with present mathematics"*), and the recommended register
  sentence built on it.
* **The C-refuter's own count correction is SUSTAINED but its framing is not.** It is right that
  `S4:973`'s `f = 3.14159, c = 1.57080` converges to `m(P)` to `~1e-08` despite lying in the
  exceptional set (a relation of order `5e5` is a Boyd–Lawton-regime subtorus). That is a
  genuinely good catch against an unrun row. But it does not rehabilitate the connection: it is
  still resonant, still outside H2, and its agreement is an accident of relation size.

### NOT SCORED EITHER WAY

* **The Z-refuter's headline-inversion finding (R-NEW-1).** It says lane Z's *"the zero set is
  empty on 3/4 of the simplex"* inverts M3-2(d)'s *"the firing region is exactly where the
  analytically interesting rate lives."* Both halves are true statements and the disagreement is
  about which one is the headline. **That is a framing question, not a finding**, and I decline to
  score it. What I *do* rule: the 3/4 is Lebesgue measure on a space the corpus adopts no measure
  on (R-NEW-5, and it is right), and the corpus's own SENSE-U carrier column is mostly singular —
  my recount from `S4:575` as written, nine runnable rows: **TWO 6, EMPTY 2, CURVE 1**, i.e.
  **7 of 9 singular against a Lebesgue 1/4.** But calling that column "the corpus's only
  unchosen sample" is itself a stretch: the *carriers* were chosen. **Reads two ways; §7.**

---

## 7. WHAT FALLS THAT IS NOT N1 — AND IT IS THE SECOND RESULT PROPOSED FOR PUBLICATION

**N2, THE MULTISET THEOREM, INHERITS N1's HYPOTHESIS EXACTLY.** `REGISTER:196-197` states
*"`lambda_B` is a function of the MULTISET of the four class weights — 24 of 24 permutations
invariant, worst spread 2.4e-15"*, and `W10_SCOPE_TABLE:3.2` restates *"24 of 24 permutations
preserve `lambda`"* with no connection qualification. Measured by me, **one variable moving (the
permutation), same connection, same estimator, same `N = 2e6`**, on K1's registered weights (12
distinct orderings, since two weights coincide):

```
f=1, c=sqrt2   GENERIC          spread 3.670e-06     (Birkhoff noise at N=2e6; all -> m(P))
f=2.0, c=1.1   RESONANT         spread 7.247e-04     6 distinct values
f=pi, c=pi/2   ORDER 4          spread 5.579e-02     2 distinct values, -0.804718956 / -0.860504844
```

**The invariance is a theorem about `m(P)`, not about `lambda`.** `m(P)` is `S_4`-symmetric
because Jensen makes `|A|` depend on the unordered pair `{p00,p10}`, `|B|` on `{p01,p11}`, and
`x <-> y` swaps `p10 <-> p01` — and `(12),(34),(23)` generate `S_4`. So **24-of-24 holds exactly
where `lambda = m(P)`, i.e. exactly under Theorem W13-1's hypotheses.** What survives at *every*
connection is only W-03's own involution `00<->11, 10<->01`, which I confirm holds to
`0.000e+00` at all three connections above.

**CUSTODY, AND IT IS THE PROGRAM'S NAMED FAILURE MODE AGAIN.** This is **not** my discovery. It is
sealed on disk in `LANE_W10_B_MULTISET/b4_involution_labels.OUT.txt` leg E (2026-08-16) —
*"SIGHTED, AND THIS IS THE ONE THAT MATTERS — the RATE ITSELF at a degenerate connection … W-03's
'the incidence labels are invisible' is therefore a statement about the GENERIC connection, and
the corpus computed almost nothing at a generic connection"* — and it was attacked and
**strengthened** by `LANE_W10_B_MULTISET_REFUTE_2/r4`. It is in **no register row**, and W-10's
own scope table states the unqualified opposite one row away. `grep -n "B-11" REGISTER_V001.md`
returns nothing. **Two of the four results this corpus proposes to publish were entangled, the
entanglement was found ten hours before this round opened, and it was under-read.** My only
addition is the co-extensivity: N1 and N2 are one result with one hypothesis.

---

## 8. THE REGISTER TEXT FOR N1

Replace *"under H2"* with:

> **N1.** `lambda = m(p00 + p10 x + p01 y + p11 xy)` under **H1** (W-11: the relative branch
> operator is class-constant diagonal — a stipulation, with COR-F an admissible alternative under
> which it fails) and the **three arithmetic hypotheses of W-13 Theorem 1**: `L(theta) = {0}`; a
> polynomial homogeneous Diophantine bound on `(alpha,beta)`; and a polynomial **inhomogeneous**
> bound relative to `Z(P)`, which is a joint condition on the connection **and the ready state**.
> The third is the one that decides, and it is not implied by the first two (W-13 Theorem 3;
> M1_08 T2(b)).
>
> **It holds unconditionally whenever the two loop angles `f, c` are algebraic reals with `f/c`
> irrational and `pi` has algebraic entries and finite `Z(P)` — by Baker's theorem on
> inhomogeneous linear forms in logarithms (Baker III, Mathematika 14 (1967) 220–228). In
> particular it holds at `f = 1.0, c = sqrt(2)` (S4:603), the only generic connection this corpus
> publishes and the pair S4 used to verify its whole `lambda` column.** It holds for
> Lebesgue-almost every connection, on every carrier and every designation (W-12 + W-13 Theorem 2
> of lane C). **It FAILS on a comeager set of connections and at all four published resonant
> pairs.** No counterexample has algebraic angles. **Boyd–Lawton does not license it and no
> published theorem states it; the note proves it.** The identification of a decay rate with a
> logarithmic Mahler measure is standard (Lind–Schmidt–Ward 1990; Kenyon–Okounkov–Sheffield 2006;
> Herman 1983 for the Lyapunov-exponent form) and **no novelty is claimed for it.**
>
> **N2 carries the same hypothesis**: 24-of-24 permutation invariance holds exactly where
> `lambda = m(P)`. At every connection, only W-03's involution survives.

---

## 9. WHAT READS TWO WAYS, AND IS SCORED NEITHER WAY

1. **MEASURE VERSUS CATEGORY.** The failure set is Lebesgue-null **and** comeager. *Reading A*:
   the a.e. statement is the right one and the comeager failure set is a curiosity, as it is
   throughout Diophantine analysis. *Reading B*: a result whose hypothesis fails on a
   topologically generic point is a result about an idealisation. **Both supported. Not scored —
   and for once it does not matter**, because Theorem W13-2 decides the corpus's own connection
   outright and neither reading touches it.
2. **WHETHER THE CORPUS LIVES IN THE HARD REGIME BY NECESSITY.** Lebesgue on the state simplex
   says singular on `1/4` (exact, re-derived); the SENSE-U carrier column says singular in 7 of 9.
   *Reading A*: the corpus's sample is unrepresentative because its carriers were chosen.
   *Reading B*: Lebesgue on the simplex is a measure nothing in the construction supplies, so the
   carrier column is the better sample. **Both supported. Not scored.** W-15's discriminator
   applies and convicts the pair: `hard regime / easy regime` is a **coined** dichotomy, and the
   third thing it hides is that the two regimes need *different theorems*, not a verdict.
3. **THE STATUS OF A COUNTABLE LICENCE.** Theorem W13-2 covers a **countable** set of connections.
   *Reading A*: that is exactly the set of connections anyone can name, so it is the physically
   meaningful statement. *Reading B*: a countable set is Haar-null, so by the corpus's own N3
   logic it is a measure-zero curiosity and the a.e. statement is the only robust one. **Both
   supported. Not scored.** They are not in conflict — both statements are true and they answer
   different questions — but which one belongs in the abstract is a judgement the note's author
   must make, not the registrar.
4. **WHY THE CORPUS'S NUMERICS ALWAYS AGREED.** *Reading A*: because N1 is true at the
   connections it computed at (which Theorem W13-2 now says it is). *Reading B*: because no
   simulation can see the failure (lane C's C7, which is right). **Both are true at once**, which
   is why agreement was never evidence and is still not.

---

## 10. CORRECTIONS AGAINST MY OWN BRIEF, AND MY OWN DEFECTS

**B-1. THE BRIEF'S BICONDITIONAL IS THREE-CLASS-SCOPED, AND I WROTE IT UNQUALIFIED.** The brief
for this round asserted *"P has a torus zero IF AND ONLY IF W-01's criterion fires. THE SINGULAR
CASE IS THE FORMATION CASE."* The L-refuter's counterexample is right and I verified it exactly:
`pi = (1/2, 3/10, 1/10, 1/10)` has `max = 1/2 <= 1/2` (hull fires) and `(S1-S2)(D1-D2) = 3/25 > 0`
(**no torus zero**). Exhaustive exact sweep: **0 disagreements on the `p00 = 0` face at
denominators 12/18/24; 52 / 200 / 500 disagreements out of 455 / 1330 / 2925 on the full
four-class simplex (11.4% / 15.0% / 17.1%).** The surviving direction is `torus zero => hull`.
This is the same accident W-09 convicted W-01's advertised virtue of, committed by the registrar
in the brief that quotes W-09.

**B-2. THE BRIEF'S `min|P| = 2.0e-04 on a 2048^2 grid` IS A WINDOW ARTEFACT AND I STATED IT AS A
FACT.** `min_{T^2}|P| = 0` **exactly** at K1's registered `pi` (X01, exact). Lane Z caught this
and is right; it is COR-E's own defect class, committed in the brief that cites COR-E.

**B-3. MY OWN VOID CONTROL, RECORDED NOT PATCHED.** `X02`'s last arm was built to violate (D2) at
`k = 1e4` and **does not**: I engineered only the **first** coordinate while `dist` is the sup over
both, so the second coordinate holds the distance up. Its row is indistinguishable from the good
arms and proves nothing. A correct two-coordinate violator is `X03` leg B and does crash. **This
is precisely the class W-08's isolation audit names as the commonest fatal defect, committed by
the registrar in the lane written to adjudicate it.** Both the void arm and the working one are in
the sealed output.

**B-4. MY OWN MISLABEL, RECORDED NOT PATCHED.** `X02`'s first version labelled the reference decay
of `min_{k<=N} dist` as `N^{-1}`; for a **point** target in `T^2` the equidistributed rate is
`N^{-1/2}`. The numbers were unaffected; the label was wrong. M1_06 and lane C both state
`N^{-1/2}` and I did not read my own carry-forward carefully enough.

**B-5. NOT READ AT THE SOURCE.** Baker III/IV, Erdős–Turán–Koksma, Lawton 1983, Dobrowolski 2017,
Cassaigne–Maillot, Herman 1983, Lubinsky 1999, Dimitrov (abstract only). All are quoted by shape
or from secondary sources. **The note cannot be written from this ruling alone.**

**B-6. WHAT I DID NOT TEST.** Whether Dimitrov's sub-Liouville estimate extends off torsion points
— the L-refuter names this as the most consequential thing it left open and it is still open, and
if it did extend it would bear directly on (D2). Whether Theorem W13-2 survives on the **curve**
strata, where `Z(P)` is infinite: it does not, as stated, and the reduction to a one-variable
Sudler product (§5) is the right route there but I did not carry it to a theorem.

---

## 11. LINEAGE

**I am Claude Opus 5. So are lanes Z, C, R and L, all four refutations, and W-07 through W-12.**
This is **layer fifteen of one block**, and the last lineage-independent boundary in this corpus
is still Fable 5 → Opus 5 at W-07. Every disagreement adjudicated above was adjudicated inside the
block; the finding that overturns three of the eight documents was produced by the same model that
produced them. **Discount this row as one block with W-07 through W-13, not as an independent
check.**

The block's named failure mode is **misnaming the operative variable**, and I ran the guard: I
grepped the register before naming anything, and the name I attach is not new — it is
*inhomogeneous Diophantine condition relative to `Z(P)`*, which M1_08 T2(c) wrote on 2026-08-16
and which the register has never contained. The block's second failure mode is **under-reading**,
and this ruling found two more instances of it (B-11's multiset finding; Dimitrov in lane L's own
fetch queue) and committed at least one itself (B-4). **The rate has not fallen. Discount this
layer too.**

**NOTHING HERE TAKES ANY OF THE PRINCIPAL'S DECISIONS.** It answers the question the round was
opened for and it attaches an answer to the second standing decision: **the Mahler note is
publishable, with a named-connection theorem it did not know it had, under a hypothesis the
register does not currently carry, with a novelty claim that must be dropped, and with N2 folded
into N1 rather than published beside it.**

---

## APPENDIX — WHERE THE PROOF IS

`LANE_W13_RULING/`, `SEALS.sha256`:

* `X00_CONVENTIONS.txt` — what was read before any code was written; conventions; isolation rule.
* `X01_exact_geometry.py/.OUT.txt` — the zero set exactly in `Q(i sqrt5)` for all three readings of
  SENSE C; minimal polynomials; `P(x0,y0) = 0` exactly; the anti-diagonal decided per labelling;
  `m(P)` over five node decades with the trend; Cassaigne–Maillot as a second route; the conical
  local shape over six decades.
* `X02_baker_forms.py/.OUT.txt` — the three linear forms written out with every ingredient
  certified; (D2) measured over seven decades on seven arms; (D1) measured over five decades of
  `M`; H2 detection; the 2^64 / 2^128 / exact-integer phase cross-check; **my void control,
  recorded**.
* `X03_ladder_and_counterexample.py/.OUT.txt` — six decades of `A_N` on seven connections at two
  labellings, with the two resonant controls; the two-coordinate one-dip construction at
  L = 500…17372; the divergence-density ladder over four decades with the local expansion
  validated to `1.6e-42`.
* `X04_audit_of_the_four.py/.OUT.txt` — N2 under all permutations at three connections; the B0b
  transcription with both factorisation certificates and the measured `|P|` spread; the
  `d = 1,2,3,5` closure ladder; the Sudler ratio identity over five decades; S1's registered state
  as a Sudler product; lane C's splice arithmetic.
* `X05_brief_and_census.py/.OUT.txt` — the brief's biconditional swept exactly at three
  denominators; the SENSE-U carrier census from `S4:575` as written; the `1/4` measure by exact
  lattice counting.
