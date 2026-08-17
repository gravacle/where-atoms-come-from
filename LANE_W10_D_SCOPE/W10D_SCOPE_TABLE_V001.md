# W-10 LANE D — THE SCOPE TABLE — V001 — 2026-08-16

**THE QUESTION THE PRINCIPAL DEFERRED A DISPOSITION TO GET ANSWERED:** which registered results
are K1-scoped and which are carrier-independent?

**METHOD.** Every load-bearing claim in `REGISTER_V001.md` — W-01 through W-09, both errata, and
W-05's N1–N4 survivor list — was read, and then **either exhibited on a four-class carrier,
exhibited failing there, or marked UNDETERMINED with the precise reason.** No row is marked
CARRIER_INDEPENDENT on the strength of reading alone.

**CARRIERS.** `B0b` (ring torus 3×3 grid, loops meet) and `B4` (spindle) — the corpus's own two
four-class rows, class multisets quoted verbatim from `S4_THE_MEASUREMENT_V001.md:576, :580`.
Controls: `B1` = K1, `B1q` (spectator), `B1p` (two classes).

**LANE FILES.** `PUBLISHED_CONVENTIONS.txt`, `w10d_1_hull.py` … `w10d_7_corpus.py` with captured
outputs, `SEALS.sha256`.

---

## 0. THE FOUR SCOPE VERDICTS, DEFINED

| verdict | means |
|---|---|
| **CARRIER_INDEPENDENT** | exhibited holding on a four-class carrier (or proved from the class-weight vector alone, with the proof on this page) |
| **THREE_CLASS_SCOPED** | true wherever at most three classes are occupied; **exhibited failing, or exhibited incomplete, at four** |
| **K1_SCOPED** | its hypothesis is a property of K1's particular numbers (V = 5, prime, 4k+1) or of K1's particular connection/state; **exhibited unavailable on B0b and B4** |
| **UNDETERMINED** | not tested here, or not testable because the artifact/code does not exist. **Never inferred.** |

---

## 1. THE HEADLINE, IN FOUR SENTENCES

**The formation functional itself is carrier-independent.** Every claim that is a statement about
the class-weight 4-vector — N1's Mahler identification, N2's multiset theorem, W-02's
character-ratio criterion, W-08's monotonicity and its exact character identity, the schedule
result, N3's Haar-null inversion, N4's fibre-wise-ness, W-06's dressed restoration *as a
mechanism* — **holds verbatim on both four-class carriers, and eleven of those were exhibited
there for the first time on this page.**

**What is K1-scoped is not the physics — it is the arithmetic of the number 5.** S2's
no-factorisation theorem, S3's rejection of its own adjunction-free alternative, W-02's *"cost:
one qubit per cell, proved minimal"*, and W-06's `V = 4k+1` wedge-growth rebuild route are all
statements whose hypothesis is *five vertices being prime*. **On B0b `C^9 = C^3 ⊗ C^3` and
`9 | 81`; on B4 `C^6 = C^2 ⊗ C^3`. Every one of those four arguments is unavailable on both
four-class carriers the corpus owns.**

**What is three-class-scoped is W-01 — the corpus's most-quoted row — and W-02's table.**
W-09's ruling is reproduced here independently and narrowed. W-02's criterion is
carrier-independent but **its published table enumerates three of the six two-element supports,
and the three it omits are exactly the three requiring class 00, which is empty on K1.**

**And one scope defect is not about carriers at all and is new here: the multiset theorem holds
at the level of `lambda` and at no finer level.** Two of 24 permutations preserve `|Z_k|`
pointwise; twenty-four of 24 preserve `lambda`. **The incidence labels are invisible
asymptotically and visible at every finite `N`.** N2 is offered for publication without that
qualification.

---

## 2. THE SCOPE TABLE

### W-01 — the formation condition on K1

| # | claim as registered | scope | evidence |
|---|---|---|---|
| 1.1 | the convex-hull criterion, *"vanishes iff 0 lies in the convex hull of **three** unit-modulus coefficients"* | **THREE_CLASS_SCOPED** | leg 1A: three occupied classes fire on **exactly 1/4** (Wendel 1962), four on **exactly 1/2** (`cos f + cos c ≤ 0`, agreeing with the hull on 200000/200000). Two structurally independent hull algorithms agree on every draw. W-09 reproduced from outside its own code. |
| 1.2 | *"it distinguishes curvature from flat holonomy, which K1 exists to separate"* | **THREE_CLASS_SCOPED, and narrower than W-09 states** | leg 1A: `f → −f` flips 99785/200000 on three classes and **0/200000** on four — W-09 confirmed. **New:** `c → −c` flips *the same 99785* on three classes and **0** on four, and the exchange `f ↔ c` flips **0 on three classes as well as on four**. So the criterion never distinguished the two *roles*; it was sensitive to each holonomy's **sign**, and that is what four classes destroys. |
| 1.3 | the `iff` quantifier | **CARRIER_INDEPENDENTLY DEFECTIVE, and the defect is large off K1** | leg 1C: on B0b's and B4's **own published weights** the hull fires on 0.4989 of the grid while `min |Z_1|` over the whole grid is `1.11e-01` and `2.72e-01`. The gap is half the parameter space, not a measure-zero technicality. leg 1D proves it exactly: `P` has **no torus zero** at either carrier's SENSE-U weights. |
| 1.4 | *"the root can never fire"* | **K1_SCOPED as stated; the true statement is CARRIER_INDEPENDENT** | leg 3B/3C: the theorem is `|S| = 1 ⇒ never`, a **class** fact. B0b has **two** class-11 vertices (no "the root") and **four** class-00 vertices; a state spread over all four class-00 vertices also never fires. Root and class coincide on K1 only because class 11 there has exactly one vertex. |
| 1.5 | the canonical **three-way** split `{v0}/{v1,v2}/{v3,v4}`, *"derivable from incidence"* | **K1_SCOPED (arity); the partition is CARRIER_INDEPENDENT** | leg 3C: it *is* the loop-membership class partition. Arity 3 on K1 (1,2,2) and B1q (1,3,3); **arity 4 on B0b (4,2,1,2) and B4 (1,1,1,3)**. The "source/record/environment shape" the row says is available without imposition acquires a fourth block with no role in that shape. |
| 1.6 | gauge-invariance of the criterion; correct trivial limit | **CARRIER_INDEPENDENT** | leg 7D: at `W_F = W_C = 1` every character is 1, so `Z_k = Σp = 1` identically on any carrier. An identity, not a test. |
| 1.7 | *"it fires on S1's own published connection"*, `0.0247 at n=42`, `recurs to 0.99994` | **K1_SCOPED and CONNECTION_SCOPED** | already relocated by the ERRATUM AGAINST W-07: those two figures were measured at the **resonant** `f=2.0, c=1.1`, not at S1's order-4 point. Untouched here. |
| 1.8 | the refutation of the build's **no-split** theorem (`dim = 5` prime ⇒ no factorisation) | **K1_SCOPED** | leg 4B: `C^9 = C^3 ⊗ C^3`, `C^6 = C^2 ⊗ C^3`. On both four-class carriers **the whole space factorises** and the theorem being refuted cannot even be stated. |
| 1.9 | *"zero addition suffices"* | **CARRIER_INDEPENDENT, conditioned on `G ≠ {1}`** | leg 4C: the compression to `span{M_dF s, M_c s}` has `dim_C = 4` on B1, B1q, B1p, B0b, B4 alike; it needs only two occupied classes with distinct characters. |

### W-02 — the crossing and the nine

| # | claim | scope | evidence |
|---|---|---|---|
| 2.1 | **the character-ratio criterion, FORMATION ⟺ `G ≠ {1}`** | **CARRIER_INDEPENDENT** | leg 3A/3B: all 15 supports enumerated; the two supports **K1 cannot realise** — `S = {00,11}` (`G = ⟨uv⟩`, sees only the product) and `S = {00,10}` (`G = ⟨u⟩`) — were run. Each has a non-trivial connection with `|Ω_N| = 1` exactly and a generic one with `|Ω_N| → 0`. |
| 2.2 | **its published TABLE** (`|S|=3 → ⟨u,v⟩`, `{0,C} → ⟨u⟩`, `{0,F} → ⟨v⟩`, `{F,C} → ⟨u/v⟩`) | **THREE_CLASS_SCOPED — incomplete, not wrong** | leg 3A: there are **six** two-element supports. The table lists three. **The three it omits are exactly `{00,10}`, `{00,01}`, `{00,11}` — the three requiring class 00, which is empty on K1.** No carrier the corpus ran could have revealed the omission. |
| 2.3 | *"recovers W-01's 'the root can never fire' as a special case"* | **CARRIER_INDEPENDENT as `|S| = 1`; the identification with *the root* is K1_SCOPED** | leg 3C, as 1.4. |
| 2.4 | `A_∞ = UHF(5·2^∞)`, simple, unital, infinite-dimensional | **K1_SCOPED in the supernatural number; CARRIER_INDEPENDENT in form** | leg 4B: the first factor is `M_V(C)`; on B0b it is `M_9`, on B4 `M_6`. W-03 already recorded that S4 silently moves this factor and never says so. |
| 2.5 | **"cost: one qubit per cell, proved minimal. This is not zero addition."** | **K1_SCOPED** | leg 4B: the minimality rests on `dim = 5` prime and on `5 ∤ 9, 5 ∤ 7`. On B0b, `9 \| 81`, so `M_9^{⊗N}` is a directed system with unital embeddings and **no adjunction at all** — S3's own CHOICE LEDGER C1 rejection ground has no analogue. On B4, `C^6` factorises outright. |
| 2.6 | `Φ_{N+1} = Φ_N · Z_{N+1}`; firing is absorbing | **CARRIER_INDEPENDENT** | legs 2B, 3B: `Ω_N = Π Z_k` computed to `N = 10^6` on both four-class carriers. |
| 2.7 | `sup_k |Z_k| = 1`, *"the carrier's recurrence is untouched"* | **CARRIER_INDEPENDENT in form; CONNECTION_SCOPED in ATTAINED-vs-APPROACHED** | leg 3E and leg 5E; see 7.2. |
| 2.8 | the monotonicity theorem is *"a tautology of any per-cell tensor system"* | **CARRIER_INDEPENDENT** | leg 3E: `\|Z_k\| ≤ 1` by the triangle inequality on non-negative weights summing to 1 — no incidence input. `max_k(\|Z_k\|−1) ≤ 0` over `k ≤ 10^6` on all four carriers at three connections. |
| 2.9 | **the nine, P-1 … P-9** | **UNDETERMINED** | not tested on any four-class carrier by anyone, including this lane. They are properties of the directed system, whose first factor is `M_V(C)`; **to determine their scope one must re-run S3 §5 with `V = 9` and `V = 6`, which requires S3's code — S3 published no lane directory.** |
| 2.10 | `λ = −0.766802` | already **ERRATUM'd**: a finite-stage value of an oscillating sequence at an exactly resonant connection | untouched here; leg 2C reproduces the resonance mechanism on B0b (`c = f` moves the rate by `3.45e-01`). |

### W-03 — the measurement, and the multiset theorem

| # | claim | scope | evidence |
|---|---|---|---|
| 3.1 | **the multiset / invisibility theorem, "24 of 24 permutations invariant"** | **CARRIER_INDEPENDENT — and its hypothesis is misnamed twice** | leg 2D. Reproduced at four classes: spread `4.4e-16` (B0b), `0.0e+00` (B4), `1.6e-11` (generic). **(a)** The registrar named the operative hypothesis *REAL NON-NEGATIVITY*. **It is REALITY alone.** Two arms with negative entries also give `\|G\| = 24`; the pointwise identity `\|a+be^{it}\| = \|b+ae^{it}\|` needs only `conj(a)b = conj(b)a`, i.e. `ab ∈ R`. **(b)** With **complex** coefficients `\|G\|` drops to exactly **8** — the Newton-polygon `D4` the registrar originally predicted. His prediction was right about the group and wrong about which hypothesis excludes it. |
| 3.2 | **the level at which it holds** — *not registered anywhere* | **CARRIER_INDEPENDENT, OBSERVABLE_SCOPED. NEW.** | legs 3D + 4A: **2 of 24** permutations preserve `\|Z_k\|` pointwise (identity and W-03's involution); **24 of 24** preserve `λ`. The incidence labels are invisible **asymptotically** and **visible at every finite `N`** in `\|Ω_N\|`. W-05 offers N2 for publication without this qualification. |
| 3.3 | the exact involution `00 ↔ 11`, `10 ↔ 01` | **CARRIER_INDEPENDENT, and stronger than the rest of 3.1** | leg 3D: it fixes `\|Z_k\|` **pointwise in k** (max dev `1.1e-16`), because it multiplies `Z_k` by a unit-modulus factor. It is the only permutation besides the identity that does. |
| 3.4 | *"the carrier's topology is inert"* | **CARRIER_INDEPENDENT — and UNTESTABLE BY CONSTRUCTION** | leg 4D. `d2` enters the functional nowhere; W-03 said so itself. **The four-class pair makes this vivid and the exhibit is declared a ZERO-VARIABLE CONTROL on this page:** B0b (`χ=0,b1=2,b2=1`) and B4 (`χ=2,b1=1,b2=2`) at SENSE-C feed the *same* weight vector to the *same* functional and of course agree to the last digit. That confirms nothing. **No experiment on any carrier can bear on this row.** |
| 3.5 | *"formation does not see incidence; the incidence labels are invisible"* | **CARRIER_INDEPENDENT, with a wording correction** | the class **multiset** is exactly the incidence content that survives, and it is what sets `λ` (B0b `log(4/9)`, B4 `log(1/2)`). What is invisible is the **labelling** of the four classes, not the incidence. |
| 3.6 | B0b `= log(4/9)` exactly; *"nine of nine carrier rates are exact"* | **CONFIRMED and STRENGTHENED** | leg 2A: Jensen's branches for B0b differ by a **constant** `9/81 > 0` at every `t`, so `m = log max(4/9,1/9) = log(4/9)` exactly (`\|diff\| 6.7e-16`); B4's differ by `(8+4cos t)/36 > 0`, so `m = log(1/2)` exactly. **Both four-class carriers have closed forms; S4's "QUADRATURE ONLY" label on B0b was wrong and W-03's correction is right.** |
| 3.7 | the charge run (`S4-1` fails at charge ≠ 1) and the SU(2) run | **UNDETERMINED w.r.t. carrier** | orthogonal axes; ERR-2 already convicts the SU(2) arm of a three-way confound. Not touched here. |
| 3.8 | *"the carrier axis is a list, not a family"* | **CARRIER_INDEPENDENT (a fact about S4)** | leg 7C: census of all 11 named carriers. **Two of eleven occupy four classes.** |

### W-04 — the errata

| # | claim | scope | evidence |
|---|---|---|---|
| 4.1 | **ERR-1**: the operative variable is **SCALARITY**, not commutativity | **CARRIER_INDEPENDENT (fibre-scoped)** | leg 6A reproduced: `T_F = I`, `T_C = diag(1,−1)`, `z = (1,1)/√2` gives overlap `0.0000000000000000` with commutator `0.0e+00` and `T_C` non-scalar. No carrier or class occupancy enters. |
| 4.2 | **ERR-2**: the SU(2) run was a three-way confound | **CARRIER_INDEPENDENT (a fact about the control)** | not re-run; it is an argument about what varied, and nothing on the carrier axis bears on it. |
| 4.3 | **ERR-3**: the gravitational term was dropped, unrecorded | **CARRIER_INDEPENDENT (corpus fact), and re-verified** | leg 7A: over the 14 `.md` artifacts now sealed, `grep -in gravit` returns the founding-design target line plus **only** the register's and the handoff's own quotations of ERR-3. `coupling constant`: 0 artifacts. `backreaction`, `edge mode`, `plaquette`: 1 each, all inside the errata themselves. |
| 4.4 | **ERR-4**: S1 is unaudited and *"smallest complex"* is false | **K1_SCOPED (it is a claim about K1) and CONFIRMED** | leg 7B: two triangles sharing an edge, one filled — `V=4 E=5 F=1 χ=0 b0=1 b1=1 b2=0`, **2 invariants**, `d1·d2 = 0`. K1's complete profile on one fewer vertex and edge. **AND ITS CLASS MULTISET IS `{11:2, 10:1, 01:1}` — three classes.** The carrier W-04 says S1 should have used is also three-class. |
| 4.5 | *"What K1 was: inert"* and the four unledgered choices | **(i)–(iii) CARRIER_INDEPENDENT; (iv) "five vertices being prime" K1_SCOPED** | leg 4B. |
| 4.6 | the residue (gauge non-factorisation ⇒ edge modes) | **UNDETERMINED** | W-06's IMP-2 already records that this leans on a theorem the register had refuted two rows earlier. Not adjudicated here. |

### W-05 — the survivor list, and the STOP legs

| # | claim | scope | evidence |
|---|---|---|---|
| 5.1 | **N1 — `λ = m(p00 + p10 x + p01 y + p11 xy)`** | **CARRIER_INDEPENDENT** | legs 2A + 2B. Jensen reduction agrees with S4's published nine to `≤ 3.5e-13`; **direct schedule-B simulation to `N = 10^6` at random generic `(f,c)` on B0b and B4 agrees with `m(P)` to `2.3e-07` and `3.6e-07`.** Nothing in the identification uses `p00 = 0`, three classes, or K1's incidence. |
| 5.2 | **N2 — the multiset theorem** | **CARRIER_INDEPENDENT, with 3.1(a),(b) and 3.2 attached** | as above. **It must not be published without the level qualification.** |
| 5.3 | **N3 — the Haar-null inversion** | **CARRIER_INDEPENDENT, given a two-loop designation** | leg 2C: 40 random connections on B0b at `N = 4·10^5` give sd `1.0e-05` and max deviation from `m(P)` of `4.7e-05`; the resonant points `c = f` and S1's order-4 point depart by `3.45e-01` and `7.89e-02`. The resonance set is Haar-null in `T^2` because there are exactly two designated holonomies — a statement about the loop designation, not about class occupancy. |
| 5.4 | **N4 — fibre-wise-ness** | **CARRIER_INDEPENDENT** | leg 6C: `dim_C A^G = V` for the fibre-wise group at `V = 5` and `V = 9`, ranks 1, 2, 3 (and `V=6, r=2`), by Haar projector spectrum with a clean gap. One-line proof on the page. **Its content is entirely the word "fibre-wise"; no carrier can bear on it.** |
| 5.5 | LEG ONE — `M_γ` lies in the gauge group | **CARRIER_INDEPENDENT, and half-scoped** | leg 6B: `M_γ` is fibre-wise on any carrier — a definition. Under the **full** gauge law it is a gauge transformation only if `θ` is constant on every edge; on B0b it would shift **8 of 18** edges. Same as on K1. |
| 5.6 | LEG TWO — the slot is already inside the carrier, `dim_C = 4` | **CARRIER_INDEPENDENT, conditioned on `G ≠ {1}`** | leg 4C: `dim span{M_dF s, M_c s} = 2` hence `dim_C = 4` on all five carriers tested, including both four-class ones and the two-class B1p. |
| 5.7 | LEG THREE — the rejections rest on `5 ∤ 9`, `5 ∤ 7` | **K1_SCOPED, and CONFIRMED** | leg 4B. |
| 5.8 | *"has a field: NO"* | **CARRIER_INDEPENDENT on the structural half; UNDETERMINED on the experimental half** | no action on 2-cells exists anywhere in the corpus (leg 7A: `coupling constant` 0 artifacts) — that half needs no carrier. W-05's own experimental headline is already withdrawn as confounded. |
| 5.9 | the WITHDRAWN list X1–X5 | **CARRIER_INDEPENDENT** | withdrawals, not claims. |
| 5.10 | the rediscovery ledger (Hepp / Bell / Zurek / SBS / AB / Donnelly-Wall) | **CARRIER_INDEPENDENT (bibliography)**, and already corrected once by W-06 | not re-adjudicated here. |

### W-06 — the audit of the audits

| # | claim | scope | evidence |
|---|---|---|---|
| 6.1 | **the dressed-algebra restoration — S3-0's broad form is FALSE, not vacuous** | **MECHANISM: CARRIER_INDEPENDENT. NUMBERS: K1 + state + connection scoped.** | leg 5, on a **reconstructed B0b validated against S4's published row** (`V=9 E=18 F=9 χ=0 b1=2 b2=1`, `d1·d2 = 0`, class multiset `{00:4,01:1,10:2,11:2}` — exact match). The tree-dressed `A_uv = conj(t_u) t_v` is invariant under the full gauge action to `1.7e-16`; S3's diagonal observables separate the branches by `5.6e-17`; **the dressed one separates by `0.1848`.** Closed form, deviation `5.6e-17`: `\|A_uv[M_dF s] − A_uv[M_c s]\| = \|A_uv[s]\| · \|W_F^{a_v−a_u} − W_C^{b_v−b_u}\|` — **so the dressed separation is a function of the pair's CLASS DIFFERENCE and nothing else about the carrier, and it is non-zero under exactly the same character condition as formation.** |
| 6.2 | the figures `4.45e-16`, `3√3/10`, `2.221e-16`, `1000 of 4000`, the Bell `R_N` transplant, the Schmidt `[1,0]`, the 52-subgroup sweep | **UNDETERMINED BY CUSTODY** | **W-06 has no artifact and no lane code** (W-07 §custody). Two figures reproduce from independent reconstructions (W-07, W-08); the rest cannot be scoped because there is nothing to scope. The Bell transplant and the Schmidt `[1,0]` are load-bearing and have **never** been reproduced by anyone. |
| 6.3 | **the wedge-growth rebuild route `V = 4k+1`** | **K1_SCOPED** | leg 4B: `V = 4k+1` is a property of 5. B4 has `V = 6` and is not in the sequence; B0b has `V = 9` which *is* `4k+1` **but `9 \| 81`, so the route's own motivation (S3's `5 ∤ 9` rejection) evaporates there.** The rebuild route is a K1 artifact twice over. |
| 6.4 | the gauge-group hypothesis is refuted; `U(1)^V` is the 0-cochain group | **CARRIER_INDEPENDENT (leg (a)); leg (b) already VOID as a control (W-07)** | the 0-cochain argument is intrinsic to a bundle over any CW base. Not re-run. |
| 6.5 | IMP-1: *"could not have failed"* voids a control, never a theorem | **CARRIER_INDEPENDENT (methodological), and applied against this lane's own 4D** | leg 4D declares its own topology exhibit a zero-variable control rather than reporting it as a confirmation. |
| 6.6 | *"the crossing contains no measurement"* (Schmidt `[1,0]`) | **UNDETERMINED** | no code; not reproduced by anyone; see 6.2. |

### W-07 and the ERRATUM against it

| # | claim | scope | evidence |
|---|---|---|---|
| 7.1 | `ord(ρ)` as the operative variable | **FALLEN** (W-08), superseded on W-07's own clause | not re-litigated. |
| 7.2 | **ATTAINED vs APPROACHED** | **CARRIER_INDEPENDENT, CONNECTION_SCOPED** | leg 5E on B0b: with `W_F=−1, W_C=−i` the dressed record is annihilated on **1000 of 4000** cells (`min 1.6e-18`); with `W_F=W_C=−1` on **2000 of 4000**; at a random connection and at `e^{i√2}, e^{i√3}` on **0 of 4000** (`min 1.1e-02`, `9.0e-03`). **The identical figure as on K1, on a four-class carrier with a different `V`, different `E`, different tree and different state.** The cut is arithmetic of the connection, not a property of any carrier. |
| 7.3 | the dressed reconstruction reproducing W-06 (`0.384349931183`, `1000 of 4000`) | **K1_SCOPED (numbers)** | leg 5D's closed form explains why: the value is `\|A_uv\|` times a character difference, and `\|A_uv\|` is a state-and-tree quantity. |
| 7.4 | *"W-06 has no artifact and no lane code"* | **CARRIER_INDEPENDENT (custody fact)** | still true at the bytes: no `LANE_W06_*` directory exists in the repo. |
| 7.5 | ERRATUM E-1 (where the recurrence figures were measured) and E-2 (`3√3/10` is reachable) | **CONNECTION_SCOPED / K1+state_SCOPED** | untouched. |

### W-08 — the race

| # | claim | scope | evidence |
|---|---|---|---|
| 8.1 | **`\|Z_k\| ≤ 1`, so `\|Ω_N\|` is MONOTONE NON-INCREASING** | **CARRIER_INDEPENDENT** | leg 3E: `max_k(\|Z_k\|−1)` over `k ≤ 10^6` at the generic, resonant and order-4 connections is `0.0e+00` on B1, B1q, B4 and `−1.1e-16` on B0b. |
| 8.2 | **the founding obstruction is FALSE AS AN INFERENCE** | **CARRIER_INDEPENDENT** | follows from 8.1 with no incidence input. |
| 8.3 | the exact character identity `\|Z_k\|² = 1 − Σ_{j<l} w_j w_l \|χ_j^k − χ_l^k\|²` | **CARRIER_INDEPENDENT — and it needs no computation at all** | leg 3E: one-line proof for any number of classes; checked at four classes over all `(a,b)` at `q = 7`, `k = 1..5`, max residual `3.3e-16`. |
| 8.4 | `SUM(1−\|Z_k\|) ≥ w_j w_l (K − 1/\|sin(τ/2)\|)`; `G = {1}` ⟺ no formation | **CARRIER_INDEPENDENT** | leg 3B: exhibited on the two supports K1 cannot realise. |
| 8.5 | the measured decay densities `0.4919 / 0.4692 / 0.5295 / 0.4692` | **K1_SCOPED (measurements)** | not reproduced off K1; the *bound* they are compared to is carrier-independent (8.4). |
| 8.6 | **the schedule result** — an adversary writing the `√K` smallest-`(1−\|Z_k\|)` cells accumulates `O(1)` while the honest schedule diverges linearly | **CARRIER_INDEPENDENT (qualitative)** | leg 4E: on B0b the adversary accumulates `0.78, 0.48, 0.15, 0.05` nats at `K = 10^4…10^7` while the honest `k_n = n` accumulates `810932` at `10^6`; B4 `0.73, 0.51, 0.16, 0.05` against `693148`. |
| 8.7 | **the four printed constants `0.606, 0.615, 0.588, 0.601` ("flat in K")** | **NOT REPRODUCED — UNDETERMINED** | leg 5F: on K1's own weights at four connections (generic, S3/S4 headline, golden-ratio, order-4) the accumulation **never stays flat through `K = 10^7`**; it decays (`0.57, 0.62, 0.36, 0.11` at the generic point). **Reads two ways and is scored as neither:** either W-08 used a connection/tie-break this lane has not identified, **or the adversary is stronger than W-08 reported.** W-08's qualitative ruling survives either way. |
| 8.8 | the exponents `K^{-1/2}` at `d_eff = 2`, `K^{-1/3}` at `d_eff = 1` | **CARRIER_INDEPENDENT** | leg 3A: `d_eff` = rank of the relation lattice = 2 for **every** occupied set of size ≥ 3, three- and four-class alike. |
| 8.9 | the onset non-uniformity `K_0 ~ t^{-2/3}` | **UNDETERMINED** | not tested here or anywhere off K1. |
| 8.10 | *"the weakest load-bearing claim is W-01's 'three unit-modulus coefficients'"*, i.e. `p00 = 0` | **SUPERSEDED by W-09** — the operative variable is four-class occupancy, not `p00` | leg 1A reproduces both halves independently: B1q has `p00 > 0` and reproduces K1 to the digit. |
| 8.11 | the isolation audit's finding (*"the commonest FATAL defect is ZERO variables moved"*) | **CARRIER_INDEPENDENT (methodological)** | applied against this lane at 4D and at 1A (B0b and B4 are **the same arm** for W-01's criterion and are reported as one). |

### W-09

| # | claim | scope | evidence |
|---|---|---|---|
| 9.1 | firing region **exactly 1/4** (three classes) / **exactly 1/2** (four) | **REPRODUCED INDEPENDENTLY** | leg 1A/1B: `0.250700` and `0.498925` on the same seed; Wendel for the first, `cos f + cos c ≤ 0` (200000/200000) and the measure-preserving `(f,c) → (π−f, π−c)` for the second. Two structurally independent hull algorithms, zero disagreements. |
| 9.2 | the operative variable is **all four classes occupied**, needing a pinch **and** a spectator | **CONFIRMED, and narrowed** | leg 1A: what dies at four classes is **sign-sensitivity in each holonomy separately**; the `f ↔ c` exchange was already blind on three classes. W-09's "curvature-aware / curvature-blind" wording overstates what the three-class case had. |
| 9.3 | *"the corpus has never run a four-class carrier through anything"* | **CONFIRMED AND EXTENDED** | leg 7C: **2 of 11** named carriers occupy four classes, and **W-04's own ERR-4 counterexample carrier is three-class too.** |

---

## 3. ISOLATION LEDGER

| leg | the one thing that moves | held fixed | arm diff verified |
|---|---|---|---|
| 1A–1D | **which characters the incidence occupies** | hull test (both algorithms), grid, seed 20260816, 200000 draws, evaluator, `(f,c)` sample | printed at the head of `w10d_1_hull.OUT.txt`; **4 distinct arms of 5 rows declared — B0b and B4 are the same arm and are reported as one** |
| 2A–2C | the class weight 4-vector | Jensen reduction, node count, `k`-range, seed | six distinct weight vectors printed |
| 2D | **the coefficient FIELD** (non-negative real → real → complex) | integrator, node count, the 24 permutations, tolerance `1e-10` | seven distinct coefficient vectors printed |
| 3A–3E | the occupied support `S` / the class weight vector | character map, `k`-range, code path | 15 supports enumerated exhaustively; no two identical |
| 4B–4E | `V` (5, 6, 7, 9, 11) and the weight vector | arithmetic routine, `K` sweep, `(f,c)` | **4D is DECLARED A ZERO-VARIABLE CONTROL on the page** rather than reported as a confirmation |
| 5 | the **carrier complex itself** (K1 → reconstructed B0b: `V` 5→9, `E` 6→18, `F` 1→9, tree, state) | dressing recipe, full gauge action, branch definition, `k`-range | reconstruction validated against S4's published row before use |
| 6 | `V` and the fibre rank `r` | Haar projector, sample count | 7 distinct `(V,r)` |
| 7 | — (byte-level and census legs, no arms) | — | n/a |

**Not an isolation:** legs 2B, 2C and 5E vary the connection as well, and say so; those rows are labelled CONNECTION_SCOPED, not carrier-scoped.

---

## 4. WHAT THIS LANE DOES NOT ESTABLISH — STATED BEFORE THE VERDICT

1. **Nine of the register's load-bearing claims are marked UNDETERMINED and none of them was
   guessed.** The largest block is **W-02's nine crossing properties**, which no lane has run on
   any carrier but K1 and which cannot be run without S3's code — **S3 published no lane
   directory.** The second is **everything W-06 computed**, which cannot be scoped because the
   artifact does not exist.
2. **`S1` and `S4` were never audited, and this lane did not audit them either.** It quotes S4's
   published class-multiset table verbatim and, for leg 5, *reconstructs* B0b and validates the
   reconstruction against that table. **If S4's table is wrong, every four-class row above is
   wrong with it, and no adversary has ever read S4's carrier code — which is not in the repo.**
3. **The four-class evidence rests on two carriers that are the same arm for the criterion that
   matters most.** B0b and B4 have different multisets and different topology, but for W-01's
   criterion they present the identical occupied set and the identical firing region. **A third
   four-class carrier with a different multiset would be a genuinely new arm and does not exist.**
4. **One W-08 figure did not reproduce (8.7) and the disagreement is unresolved.** It is recorded
   two ways and scored as neither.
5. **A scope table is not a disposition.** It says where the results live, not whether the program
   should continue. **Nothing here takes any of the principal's three decisions.**
6. **This is layer TEN of one lineage.** W-07, W-08, W-09 and this lane are all Opus 5. Three of
   the four named an operative variable and were corrected by the next. **Discount this one.**
   Its likeliest failure mode is the corpus's recurring one — misnaming the operative variable.
   The names offered here are: **for W-01, four-class occupancy (W-09's, verified, narrowed);
   for the multiset theorem, REALITY of the coefficients and the ASYMPTOTIC level; for the
   minimality/rebuild block, the primality of 5.**

---

## 5. THE ONE-LINE LEDGER

> **The formation functional is carrier-independent and always was — it never sees anything but
> the class-weight vector. What is K1-scoped is the C\*-algebraic superstructure built on top of
> it, and its hypothesis is not a topological fact but the primality of the number 5. What is
> three-class-scoped is the corpus's most-quoted sentence and the table under it. And the one
> result offered for publication that needs a qualification nobody has written is the multiset
> theorem, which is true of the rate and false of the record at every finite stage.**
