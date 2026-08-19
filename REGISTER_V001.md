# WHERE ATOMS COME FROM — REGISTER V001

**PURPOSE.** Indexed by QUESTION, phrased the way someone would re-ask it. Before opening any
line of work, search this file for the question first.

Every row carries: the question as it would naturally be re-asked · the ruling · **where the
proof is** · and **the exact condition under which the row reopens**. A row with no reopen
condition is closed permanently.

**APPEND-ONLY.** Rows are never edited, only superseded by a later row citing the earlier one.
**This register adopts nothing and derives nothing.** It records where rulings live.

**POINTER RULE (custody §1):** no term appears in a row without a digest, a `file:line`, or a
named ruling behind it.

---

## HOW TO USE THIS REGISTER

Search by the question you are about to ask. If it is here, read the row and its reopen
condition before commissioning anything.

---

## W-01 — CAN A FORMATION CONDITION BE WRITTEN ON K1'S OWN DATA, WITHOUT A STIPULATED SPLIT? (S2; construction pair, build + adversarial audit. **AUDIT REFUTED THE BUILD AND WROTE THE CONDITION THE BUILD DECLARED UNWRITABLE.**)

**RULING: YES — WITH ZERO ADDITION AND NO IMPORTED SPLIT.**

**THE CONSTRUCTION.** Extend loop transport from a single fibre to the section space
`Gamma(L) = C^5` that the build itself derived:

```
(M_gamma s)(v) = W(gamma) · s(v)   for v on the loop
(M_gamma s)(v) = s(v)              otherwise
```

These are **unitary and non-scalar** on `Gamma(L)`. Then

```
<M_dF s, M_c s> = conj(W_F)·W_C·p_0  +  conj(W_F)·(p_1+p_2)  +  W_C·(p_3+p_4)
```

which vanishes **iff 0 lies in the convex hull of three unit-modulus coefficients** — verified
against brute-force simplex minimisation, 1369 grid points, 0 mismatches.

**IT FIRES ON S1'S OWN PUBLISHED CONNECTION.** With `W_F = -1`, `W_C = -i` (S1 §6) and
`p = (1/2, 0, 0, 1/4, 1/4)`, the overlap is `i/2 - i/2 = 0` **exactly**.

**PROPERTIES, ALL CHECKED:** gauge-invariant · correct trivial limit (no formation at
`W_F = W_C = 1`) · **the root can never fire**, independently reproduced · and it
**distinguishes curvature from flat holonomy**, which K1 exists to separate and which the
build's own rival construction could not see (its self-reported defect F1).

**THE ONLY LOAD-BEARING OBSTRUCTION IS DURABILITY.** Over `n` circuits the condition fires —
`0.0247` at `n = 42` — and then **recurs to `0.99994`**. K1 is finite and closed; this is a
**reversible write, not a durable record**, reproduced from first principles on five vertices
rather than inherited.

**WHAT THE BUILD CLAIMED, AND WHY IT WAS WRONG.** It returned a threefold impossibility —
rank, no time, no split. **Two are false and one is overstated.** Its rank theorem was proved
on `L_v0 = C`, an alternative its **own choice ledger had rejected**, and then applied to the
`Gamma(L) = C^5` it had chosen; `|<s, M_dF s>| = 0` on K1's own connection refutes it directly.
Its no-split theorem — resting on `dim = 5` being prime — is refuted: reduced states exist
(the build computed one), a canonical three-way direct-sum split exists, and `C^4 = C^2 (x) C^2`
on the non-root subspace. Only a weaker structural half survives: no *canonical* factorisation.
Its no-time claim is overstated — **circuit count is carrier-supplied discrete time.**
And its "minimum addition" question dissolves: **zero addition suffices.**

**A RESULT INSIDE THE CORRECTIONS THAT MATTERS ON ITS OWN.** A canonical three-way split
**exists on K1 and is derivable from its incidence** — `{v0} / {v1,v2} / {v3,v4}`: the root, the
filled triangle's other vertices, the unfilled triangle's other vertices. The
source/record/environment *shape* the predecessor introduced by imperative is **available here
without being imposed.**

**WHERE THE PROOF IS.** Build `S2_FORMATION_CONDITION_ON_K1_V001.md`
sha256 `248ce856efaef157c68e818dde589d0200bbc1dd9fd9fc1fcc8cdc7bc88734d9`.
Audit `S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md`
sha256 `0bea11bd4b7764f65c8d44cc8812d43bdcc569d34f51d0c46d82851a2efbd0d5`.
All 16 of the build's file:line pointers resolve; all 5 digests verify at bytes; alpha not
engaged. Numeric corrections carried: the build's §3.4 table is non-reproducible and its
corrected row values are `0.000000 / 0.816497 / 0.816497 / 1.000000`; the recurrence figure
`0.994373` is a window artefact and the true supremum is at least `0.999793`.

**GRADE: ADVERSARIALLY-CHECKED, NOT INDEPENDENTLY-CORROBORATED.** Build and audit share a model
lineage (custody §4). A failure mode common to that lineage passes through both.

**REOPENS IF:** a lineage-independent lane fails to reproduce the firing on S1's connection ·
the convex-hull criterion fails on a carrier where it should hold · or the three-way split's
derivability is shown to depend on K1's particular vertex count rather than on its incidence.

---

## W-02 — DOES A DIRECTED SYSTEM OVER K1 ESCAPE RECURRENCE, AND DOES IT PRODUCE THE NINE? (S3; construction pair. **AUDIT CONFIRMED-WITH-CORRECTIONS — every computable claim reproduced in a fresh independent implementation, no arithmetic error found anywhere.**)

**RULING: THE CROSSING EXISTS. THE LIMIT ESCAPES. SEVEN OF THE NINE HOLD OUTRIGHT; TWO HOLD IN THE LIMIT AND FAIL AT EVERY FINITE STAGE — WHICH IS WHY THE PREDECESSOR RECORDED ONE OF THEM FALSE.**

**THE SYSTEM.** `A_N = M_5(C) (x) M_2(C)^(x)N`, dimension `25·4^N`, embeddings
`iota_N(X) = X (x) I_2` — verified unital, multiplicative, *-preserving and isometric — with
`A_infinity = UHF(5·2^inf)`: simple, unital, infinite-dimensional. **Cost: one qubit per cell,
proved minimal. This is not zero addition.** The record slot is non-abelian by necessity: a
unitary write onto a reset carrier leaves *pure* slot states, and pure states of an abelian
algebra have overlap only 0 or 1.

**THE TRAP WAS REAL AND IS DISARMED BY COMPUTATION.** Repeated circuits of one loop span
**3 dimensions at N=1 and 3 at N=100.** Circuits grow no algebra and escape nothing.

**THE MECHANISM IS NOT WHAT IT LOOKS LIKE, AND THE AUDIT'S CORRECTION IS THE FINDING.**
The escape is **not** "recurrence is gone" — **the carrier's recurrence is untouched**,
`sup |Z_k| = 1`. It is **"the record grew."** And the monotonicity theorem carrying the escape
is a **tautology of any per-cell tensor system**: it holds identically on families that never
form at all. **All the content sits in the divergence `Sum (1 - z_n) = inf`.**

**AND THE ESCAPE IS CONDITIONAL IN TWO WAYS, NOT ONE.** The build flagged the **cell schedule**
— an adversarial schedule locked to the carrier's near-recurrences defeats the crossing
entirely. The audit found the second, flagged nowhere: **the support of the ready state.**

**THE CRITERION, CORRECTED AND PROVED, VERIFIED 10/10 ACROSS ALL FOUR SUPPORT CLASSES.**
With `u = conj(W_F)`, `v = W_C`, characters `chi_0 = uv`, `chi_F = u`, `chi_C = v`,
support `S = supp(p)`, and `G = <chi_a/chi_b : a,b in S>`:

```
FORMATION OCCURS  <=>  G != {1}
```

`|S|=3 -> G=<u,v>` · `S={0,C} -> <u>` · `S={0,F} -> <v>` · `S={F,C} -> <u/v>` ·
`|S|=1 -> never`, which **recovers W-01's "the root can never fire" as a special case.**
**Formation is a group-theoretic condition on where the record sits.**

**THE NINE.** P-2 persistence · P-3 thresholded non-return (**exactly, not asymptotically**) ·
P-4 recoverability · P-5 redundancy · P-6 asymptotic centrality · P-7 sector-hood ·
P-8 inductive compatibility — **all HOLD**, each with its computation.
**P-1 durability-with-irreversibility and P-9 orthogonal reduced supports HOLD IN THE LIMIT AND
FAIL AT EVERY FINITE STAGE**, exhibited with an explicit reversing unitary at N = 1..5.
**P-9 is resolved: `||omega_F^N - omega_C^N|| -> 2.000000000000`. The predecessor recorded it
false because it is a limit property and was evaluated at a finite stage.**

**TRANSPORT (Q5): YES, MULTIPLICATIVELY.** `Phi_{N+1} = Phi_N · Z_{N+1}`, and the difference
functional is compatible with the embeddings. **Firing is absorbing.** And the limit converts
W-01's **binary** convex-hull criterion into a **rate** `lambda = -0.766802`, defined everywhere
off the trivial connection.

**CORRECTIONS THAT MATTER.** A claimed unital embedding of `M_5 (x) M_2^{(x)N}` into
`M_5^{(x)(N+1)}` is **false for every N**, by the divisibility criterion the same section
states. And the universal claim that every non-trivial connection forms is **false** — four
exhibited families never form, one of them **on K1's own published ready state**.

**WHERE THE PROOF IS.** Build `S3_THE_CROSSING_V001.md`
sha256 `cbf1d79679ca2ecf3ee260e8a6467062e3a93260b2325cfc287c261e9b4469cb`;
audit `S3_THE_CROSSING_AUDIT_V001.md` per its sidecar. The build self-reported eleven defects
including a one-character digest error in its own custody block, caught by the pointer rule and
recorded rather than silently fixed.

**GRADE: ADVERSARIALLY-CHECKED, NOT INDEPENDENTLY-CORROBORATED** (custody §4).

**REOPENS IF:** a schedule condition is found that is carrier-intrinsic rather than adversarial ·
the ready-state support condition is shown to follow from the incidence rather than being free ·
the rate `lambda` fails to vary with the connection at S4 · or a lineage-independent lane fails
to reproduce the character-ratio criterion.

---

## ERRATUM AGAINST W-02 — THE RATE AS REGISTERED IS WRONG (registrar error, 2026-08-16)

W-02 records `lambda = -0.766802`. **That number is wrong under every reading and the error is
mine, not the build's.** It is S3's **N = 4000 finite-stage value of an oscillating sequence**,
transcribed into the register as if it were a converged rate. The corrected statement of record:
S3's headline connection `f = 2.0, c = 1.1` is **exactly resonant** (`-11f + 20c = 0`), so its
orbit is not dense in `T^2` and its average is taken over a **subtorus**; it converges to
**`-0.767014993`**, while the generic torus value is **`-0.767507880` = m(0.4+0.3x+0.3y)**.
S3 section 6(f)'s "rationally independent: orbit dense in T^2" is FALSE and section 4.3's
"Weyl equidistribution AGREES" is a false agreement. Every row of S3 section 5.7 is a subtorus
value. **No W-02 verdict is withdrawn** — the crossing, the UHF limit and the nine stand — but
the mechanism sentence and the rate are corrected here.

---

## W-03 — S4, THE MEASUREMENT, AND ITS AUDIT: **THE HEADLINE IS REFUTED. THE STRUCTURE IS NOT UNIVERSAL.**

**RULING: "THE TOPOLOGY IS INERT; FORMATION SEES LOOP INCIDENCE" IS NOT EARNED. THE FIRST HALF
IS TRUE BUT IS A THEOREM OF THE TRANSPORT'S DEFINITION, NOT A FINDING OF THE CONTROLS. THE
SECOND HALF IS FALSE.**

Build `S4_THE_MEASUREMENT_V001.md` (`4277f0ef…`). Audit: six independent lanes re-deriving from
scratch with their own code, refuters attacking each lane's load-bearing claim, a completeness
critic. **31 agents. Every reproducible number in the build reproduced** — worst agreement
`0.0e+00`, and thirteen of thirteen rows of the S3 corpus block verified for the first time.
**The arithmetic was never the problem.**

**WHY THE CONTROLS COULD NOT HAVE FAILED.** All four controls vary exactly **one** object: the
pushforward of the ready state onto the four vertex classes, `pi = (p_00,p_10,p_01,p_11)`.
Control 1 (the fill) is **vacuous** — `pi` is unchanged, so `lambda` is unchanged **by
identity**. And `d2` — the 2-cells, the entire topological content — **enters the formation
functional nowhere**. "Topology is inert" is therefore an **analytic fact about the
construction**, true before any carrier was built. The controls tested nothing.

**AND FORMATION DOES NOT SEE INCIDENCE.** `lambda_B` is a function of the **MULTISET** of the
four class weights — 24 of 24 permutations invariant, worst spread `2.4e-15`. **The incidence
labels are invisible.** Two exact refutations of the build's own key exhibits:
- A carrier **retaining the pinch** (filled triangle + unfilled square at `v0`) reproduces the
  full `0.0634` shift, `|diff| = 0.0e+00`. The pinch is **not necessary** for it.
- Moving weight **off** the pinch with the incidence **untouched** closes the gap continuously
  to exactly zero. The pinch is **not sufficient** either.
- **PINCH AND SPECTATOR ARE THE SAME OBJECT**, exactly: multiplying `Z_k` by
  `conj(u)^k conj(v)^k` leaves `|Z_k|` fixed and maps class `(a,b) -> (1-a,1-b)`, so
  `00 <-> 11` and `10 <-> 01` is an exact symmetry at **every** connection (max deviation
  `6.55e-15` over 2000 samples).

**THE MODALITY NAMED IN S4'S OWN BRIEF AND NEVER RUN: CHARGE.** Named again in the sealed
upstream S2 as "**this is the S4 knob**". The build never ran it; six audit lanes and their
refuters never noticed. The critic ran it: **THEOREM S4-1 FAILS.** Exponents `(1,0),(2,0),(3,0)`
give `|S| = 3` with `rank G = 1`, against S4-1's "rank 2 iff `|S| >= 3`"; and `q = (1,2,2,2,2)`
moves `lambda` from `-0.767508` to `-1.200555`. **The four-class taxonomy is a charge-1
statement.** Winding was also never non-vacuously run; the critic ran it and it **discriminates
the two schedules** — bearing directly on F1, the flag S4 declared undecidable.

**THE UNIVERSALITY QUESTION, ANSWERED EARLY AND NEGATIVELY.** Two of S4's three deepest results
are theorems about the **character lattice of an abelian group** — "`lambda_B` is a function of
`L` alone" is Pontryagin duality on `T^2`, and Theorem S4-1 is the four corners of a square in
`Hom(U(1)^2, U(1)) = Z^2`. **Neither has an SU(2) analogue.** Section 2's foundational sentence,
"every vertex phase of `s` cancels", is **false the moment the holonomy is not a scalar**:
six ready states with **identical class weights** give `|Z_1|` spread **0.4247** under SU(2),
against exactly `0.0` under U(1). **The construction's core is group-specific.**

**FURTHER CORRECTIONS SURVIVING ATTACK.** `B0b` is `log(4/9) = -0.8109302162163288` exactly, so
**nine of nine** carrier rates are exact, not eight; `(pi,pi)` is a **third strict saddle**, not
a local minimum (the build's own table contradicts `chi(T^2) = 0` at `1-2+2 = +1`); the
exceptional-value split is **527/314/213**, not 638/380/258; section 5 misquotes a reopen clause
**the register does not contain**. Lawton's theorem (1983) is missing from the IMPORT AUDIT.

**THE CARRIER AXIS IS A LIST, NOT A FAMILY.** No cellular collapse, quotient, or map connects
any two of the ten complexes. And the entire C*-algebraic layer of W-02 is **absent from S4** —
which silently moves the first factor from `M_5(C)` to `M_V(C)` as `V` runs 4 to 11, and never
says so.

**CORRECTED HEADLINE OF RECORD:** *under a U(1) connection at unit charge with two designated
loops, formation's rate is a function of exactly two things — the pushforward of the ready state
onto the four vertex classes, and (under the canonical clock) the relation lattice of the two
holonomies; the carrier's topology is inert **by construction rather than by experiment**, the
incidence labels are invisible off a dense measure-zero set, and whether any of it survives
charge, a non-abelian group, or a third loop was never tested.*

**GRADE: ADVERSARIALLY-CHECKED on the axes it ran; NOT INDEPENDENTLY-CORROBORATED, and not
adversarially checked at all on the axes it did not run.** Route to corroboration specified:
a lane outside this model lineage; the charge and three-loop runs with the four controls re-run
under each; publication of `d1`, `d2`, grid and seed conventions; this erratum entered.

**REOPENS ON:** the charge run · the SU(2) run · a third loop · a third schedule (only two exist
anywhere on the page, and the A/B dichotomy may be a knife-edge between two arbitrary choices).

---

## W-04 — THE STEP-BACK: FOUR ERRATA, THE EXPLANANDUM ADJUDICATED, AND THE RESIDUE LOCATED

Nine lanes: archaeology at the bytes, two opposed positions argued at full strength, three judges
on distinct lenses, an explanandum lane with literature access, a synthesis required to commit to
one of four recommendations including STOP. **22/22 sidecars verified; 309,047 bytes swept.**

### ERR-1 — S2's CENTRAL CONDITION IS GLOSSED FALSELY, AND THE FALSE HALF PROPAGATED EVERYWHERE

S2 states its gate condition **correctly** at `:414-416` — **(R1)** `dim H >= 2`, **(R2)** `T, T'`
not both **scalars** — and correctly again at `:514`. It then glosses it at `:421` as *"The
connection must be non-abelian in its action"* and titles §4.3 *"Rank 2 with a non-abelian
structure group."* **THE GLOSS IS FALSE.** One line refutes it: `T_F = I`,
`T_C = diag(1,-1)`, `z = (1,1)/sqrt(2)` gives overlap **exactly 0** with commutator **exactly 0**.
**An abelian, commuting, non-scalar connection writes the gate.** The only control S2 offered for
the gloss tested the scalar case — a vacuous control. **The operative variable is SCALARITY, not
COMMUTATIVITY.** The false half propagated into S3, S4, W-03, and into every lane and judge of
this synthesis. **This is the single most consequential undetected error in the corpus.**

### ERR-2 — W-03's UNIVERSALITY KILL IS NOT ESTABLISHED (registrar error in reporting)

**SU(2) has no faithful rank-one representation.** The SU(2) run therefore changed the **fibre
rank**, the **transport's scalarity**, and the **group's commutativity** simultaneously, and the
resulting `|Z_1|` spread of `0.4247` was attributed to the group alone. **That is a three-way
confounded control assigned to one factor — structurally the identical defect W-03 convicted S4's
controls of committing.** W-03's own ruling reads "the controls vary exactly one object" and calls
that vacuity; its own control varied three and named one. **CORRECTED TEXT OF RECORD:** *the core
is specific to **scalar transport on a rank-one fibre**; the gauge group was never independently
varied.* The universality question is **OPEN**, not answered.

### ERR-3 — THE GRAVITATIONAL TERM WAS DROPPED BEFORE ANY CONSTRUCTION, UNRECORDED

`grep -in "gravit"` over all 11 sealed artifacts returns **EXACTLY ONE LINE** in 309,047 bytes:
`FOUNDING_DESIGN_V001.md:13`, inside the target question. **Zero** in S1, S2, S3, S4, their
audits, the register, custody. §1 poses a tri-partite target (gravitational, electromagnetic,
alpha); **§3 declares the object to be the durability map and the gravitational term is simply
gone.** No register row records the drop. Structurally there is also **no action, no coupling
constant, no Gauss law, no plaquette weight, no backreaction** anywhere — and *a U(1) connection
with no action is a fixed background phase assignment, not a field.* **The tri-partite target has
ZERO of three terms built, not one of three: the electromagnetic term is not banked either.**

### ERR-4 — S1 IS UNAUDITED AND ITS ONE STATED REASON IS FALSE

`S1:43-44` claims K1 is *"the smallest complex carrying one face and one independent cycle at the
same time."* **False at the bytes:** two triangles sharing an **edge**, one filled — `V=4, E=5,
F=1`, `chi=0`, `b1=1`, `rank(d2)=1`, **2 invariants** — is K1's complete invariant profile on one
fewer vertex and one fewer edge, **and it has no pinch.** S1 carries **no CHOICE LEDGER, no
IMPORT AUDIT, no FLAG BLOCK, no register row, and no paired default-refute audit**, against
custody §4. **The artifact fixing the carrier, the fibre rank, the gauge group and the
two-invariant split is the only one in the corpus no adversary ever read.** The pinch itself is
provenance *"Carrier chosen by the principal, 2026-08-16"* with no reason recorded; the word
"pinch" enters at S4, retrofitted.

### THE EXPLANANDUM: **PARTIAL-RESIDUE** — MOSTLY SOLVED AND UNRECOGNISED

**The nine crossing properties map one-for-one onto published criteria.** P-5 redundancy is
quantum Darwinism's `R_delta`. **P-9 is verbatim the Spectrum Broadcast Structure condition
`rho_k^i rho_k^j = 0`** (Korbicz et al.). P-7 sector-hood is **Hepp disjointness**. And the
construction meeting them is **the Coleman-Hepp model, with Bell's 1975 objection re-derived as a
discovery** — S3's "P-1 and P-9 hold in the limit and fail at every finite stage, with an explicit
reversing unitary" **is Bell's objection**. Corpus literature contact: Zurek **0**, Hepp **0**,
Bell **0**, Haag **0**, Araki **0**, "superselection" **0**, "decoherence" **1** (an aside).
**The corpus imports its mathematics by name (Mahler 12, Cassaigne 18, Glimm 2, Lawton 1) and
reinvents its physics anonymously.**

### THE RESIDUE — REAL, AND NOT WHAT THE PROGRAM WAS CHASING

**In a gauge theory the physical (gauge-invariant) algebra does not factorize into system and
environment. Every version of the decoherence / einselection / quantum-Darwinism account
presupposes a given tensor factorization — so on gauge-invariant ground it has no place to
stand — and nothing in the literature or in this corpus DERIVES which adjoined boundary degrees
of freedom constitute the record rather than stipulating them.**

**S2 Theorem 2's no-factorization result IS that non-factorization, and the one-qubit-per-cell
adjunction IS the edge-mode construction** (Donnelly-Wall; Donnelly-Giddings) — which exists
precisely because a gauge theory's Hilbert space must be extended to factorize. **The program
arrived there from a direction that field does not take, did not recognise it, and did not stay.**
It arrived by way of **5 being prime**, not by way of the constraint. Corpus occurrences of
"edge mode", "entangling", "Gauss law", "constraint" in the gauge sense: **0**.

### WHAT K1 WAS

**Inert.** The complex is invisible to the functional. Everything that mattered was decided in
**four unledgered choices laid on top of it**: (i) **rank-one fibres** (S1 §3, an inherited
predecessor convention) — which *entails* U(1) and entails S2's no-go; (ii) **loop transport as
scalar multiplication** rather than edge-by-edge parallel transport — the entire content of
Theorem S3-0, already refuted in scope by the corpus's own sealed **COR-F**, which exhibits a
bona fide non-diagonal unitary transport `T` around the same cycle with `T^3` exactly W-01's
operator; (iii) **no action on the 2-cells** — the complete and trivial explanation of "topology
is inert"; (iv) **five vertices being prime** — which forced the central architectural move.
**U(1) was never chosen.** It is the connected automorphism group of a one-dimensional Hermitian
space; the group follows by necessity from the fibre. The U(1)-versus-SU(2) question is malformed.

### RECOMMENDATION OF RECORD: **NARROW**

Publish the one real result — *formation is a function of exactly the weight multiset and the
relation lattice, blind to topology, to incidence labels and to the pinch, with pinch and
spectator exchanged by an exact symmetry* — **with the Hepp / Bell / Zurek / SBS /
Aharonov-Bohm / Lind-Schmidt-Ward placements it lacks.** Enter these errata. Run **exactly two
cheap decisive tests**: (1) **THE ACTION TEST** — write any action (`beta*Re(W_F)`) and sample the
connection; if `lambda` moves the program has a field for the first time, and if it does not the
construction is provably a statement about a fixed background phase; (2) **THE SCALARITY TEST** —
re-run the S2 gate and Theorem S3-0 at **rank two under U(1)xU(1)**, abelian by construction,
commutator exactly zero, non-scalar. **And rename the project.** Confidence HIGH on the diagnosis.

**IF THE SCALARITY TEST RETURNS "S3-0 FAILS", THE HONEST RECOMMENDATION CONVERTS TO STOP** — the
program's spine would then be an artefact of a convention nobody chose.

### THE FINDING THAT GOVERNS HOW EVERYTHING ABOVE IS WEIGHED

**THE AUDITS COMMITTED THE EXACT DEFECT THEY CONVICTED THE BUILDS OF, AT THEIR OWN HEADLINE, AND
NOBODY AUDITED THEM.** Verified twice at bytes (ERR-1, ERR-2). **The corpus's negatives are not
more reliable than its positives — they came off the same production line with the same failure
mode, one floor higher.** Discount an audit finding by the same factor as a build finding.

---

## W-05 — THE TWO DECISIVE TESTS: **NO FIELD. THE CROSSING WAS NEVER NECESSARY. RECOMMENDATION OF RECORD: STOP.**

Both tests W-04 commissioned ran, each verified on two independent lenses (isolation; mathematics).

### ISOLATION HELD ON NEITHER TEST — ALL FOUR VERIFIERS RETURNED "DOES NOT HOLD"

**The action test's positive half is destroyed twice, independently.** Its headline compared two
carriers "differing only in which triangle the 2-cell is glued to." **`K1^F` and `K1^C` are the
SAME 2-COMPLEX** — four explicit simplicial automorphisms carry `FACE_F` to `FACE_C`,
`max|P·D1 - D1·P| = 0.0`. And the entire `0.267` effect is **reproduced bit-for-bit with the face
never moved**, by the transposition `p10 <-> p01` of the class weights alone, **max deviation
`0.000e+00` over 400 random `(p,beta)`**. Its declared held-fixed background (weights 0.45/0.15)
is an unledgered import supplying **100% of the signal**; under every symmetric ready state in the
corpus the effect is **exactly zero**.

**The scalarity test's headline arm is a four-variable move**, driving the face holonomy from
S1's published `W_F = -1` to `+1` — **zeroing the curvature on the one carrier built to separate
curvature from flat holonomy** — with a ready state supported on `v0` alone, so the four non-root
fibres carry weight exactly zero and the computation reduces to W-04's abstract 2x2 exhibit.

**THIS MAKES FOUR CONSECUTIVE CONFOUNDED HEADLINES: S4's controls (caught by W-03), W-03's SU(2)
run (caught by W-04), and now BOTH W-05 tests — the last two committed by lanes that had READ the
conviction and opened with ISOLATION LEDGERS written to prevent it. That is not bad luck. It is a
property of the production line.** Anything from either test requiring an experiment is WITHDRAWN.

### HAS A FIELD: **NO** — ESTABLISHED, NOT MERELY UNESTABLISHED

For any fixed connection the rate is a deterministic function of `(f,c)` and does not depend on
`beta` under either schedule. **The action supplies a PRIOR over a fixed background, not a field:**
no backreaction, no Gauss law, no constraint, no equation of motion. The build's own flag conceded
it — *"What was built is an ENSEMBLE, not a field. The coupling runs one way."*
**And the affirmative claim is not merely unsupported, it is INVERTED:** a genuine topological
change (two 2-cells on one triangle, `chi 0->1`, `b2 0->1`) is **exactly absorbed by
`beta -> 2*beta`, deviation `0.000e+00`**, while two carriers sharing `chi` and `F` differ by
`0.349`. **The action did not make topology matter. It made the plaquette COUNT matter, and the
lane read one for the other.** The "two schedules, opposite verdicts" result is an **order-of-limits
artefact**: `Omega_1 = Z_1`, so A and B coincide exactly at `N=1` and differ by `C(beta)/N`.

### CROSSING NECESSARY: **NO. IT WAS NEVER NECESSARY.**

**LEG ONE — THEOREM S3-0 IS VACUOUS, WHICH IS WORSE THAN FALSE.** W-01's loop transport
`M_gamma` is **literally an element of the gauge group `U(1)^V`** —
`||M_gamma - gauge_op||_F = 0.000e+00` on S1's published connection and over 2000 random ones;
`max||M* X M - X||_F = 9.7e-16` over 2000 random gauge-invariant `X`. **So S3-0 is the identity
that gauge-invariant observables are invariant under gauge transformations.** S3's own proof says
it twice — *"M_gamma is diagonal"* and *"the gauge-invariant algebra is exactly the diagonal"* are
the same fact. **It could not have failed under any variation of rank, group, scalarity or
commutativity — which is exactly why it survived the test designed to break it. The program's
spine is a definition.** Shrink the gauge group to the global U(1) and an ordinary carrier
observable separates the branches at **2.0000000000000004**.

**LEG TWO — THE SLOT IS ALREADY INSIDE THE CARRIER.** S3 §2.5 specifies `R` as non-abelian,
holding both branches as PURE states with overlap the cell's comparison value `Z`, `dim_C >= 4`.
**The compression of `M_5(C)` to `span{M_dF s, M_c s}` supplies exactly that with nothing
adjoined:** `dim_C = 4`, non-abelian, both branches pure (norms `1.000000000000`), overlap
agreeing with `Z` to `6.2e-17` — **on GENERIC ready states**, not the root-delta state.
**S3 proved `dim R >= 4` by the two-non-parallel-vectors argument and then looked for the `C^2`
outside instead of inside.** What the crossing actually bought was ONE thing: **exemption from an
undeclared premise.** "The record must be gauge-invariant" (COR-J, sealed) is applied
**asymmetrically** — fatal inside the complex, waived by fiat outside it, since S3's own F6 says
*"The record slots carry no gauge action."* The same `M_2(C)` is excluded in the carrier and
admitted outside it, by a premise waived at the boundary.

**LEG THREE — RECURRENCE DOES NOT FORCE ADJUNCTION, BY S3'S OWN LEDGER.** The obstruction is real
and is NOT withdrawn: `K1` is finite and closed, `sup|Z_k| = 1`, so no finite stage holds a durable
record. **But that forces infinitely many CELLS, not an algebra outside the carrier.** S3's CHOICE
LEDGER C1 considered the adjunction-free system `M_5(C)^(x)N` and recorded it **"LEGITIMATE BUT NOT
MINIMAL: 24 extra complex dimensions per cell where §2.5 proves 4 is the floor."** It was rejected
on **a dimension count**, not on necessity. The two alternatives rejected on necessity were killed
by S3-0 (the tautology) and by `5 ∤ 9, 5 ∤ 7` (an accident of five vertices being prime).
**Strip the tautology and the prime and nothing forces the crossing at all.**

### THE CONVERSION, HANDLED HONESTLY IN BOTH DIRECTIONS

W-04's condition — *"if the scalarity test returns 'S3-0 FAILS', the recommendation converts to
STOP"* — **DID NOT TRIGGER ON ITS LETTER. S3-0 SURVIVED**, at rank two, abelian, non-scalar, and
non-abelian too. **The conversion is made anyway, on a ground W-04 did not anticipate and which is
strictly stronger: a theorem that FAILS is contentful. S3-0 is VACUOUS.** And the entailment fails
independently of whether S3-0 is true, because the slot is already in the carrier.

### **RECOMMENDATION OF RECORD: STOP — THE PROGRAM AS CONSTITUTED**

*"The object it was built to explain — a forced crossing out of the carrier — does not exist."*
**STOP is not "burn it": it is "the thing you were chasing is not there, and here is the real thing
you found while not finding it."**

### WHAT SURVIVES — THE PUBLISHABLE RESULT

**N1. THE RATE IS A LOGARITHMIC MAHLER MEASURE:**
`lambda = m(p00 + p10·x + p01·y + p11·xy)` — the per-cell decay rate of a branch comparison on a
gauge-theoretic carrier is exactly the Mahler measure of the polynomial whose coefficients are the
ready state's pushforward onto the four loop-membership classes. **The identification is new; it
inherits the entropy theory of algebraic `Z^d`-actions wholesale.** Six to eight pages, one theorem.

**N2. THE MULTISET / INVISIBILITY THEOREM** with its exact involution — *the one thing in the
corpus that survived every attack, including this round's.*

**N3. THE NULL, INVERTED AND CORRECTLY STATED:** the rate is invariant under every absolutely
continuous connection measure, so **no local Wilson action of any form at any finite coupling can
move it** — because the whole resonance structure lives on a Haar-null set. **Averaging does not
blur the resonance structure, it DELETES it.** `P(rank L > 0) = 0`.

**N4. THE VACUITY OF THE SPINE**, as a methodological result that generalises past this program:
**defining transport as scalar multiplication by the holonomy collapses transport into gauge, and
every "the carrier cannot hold it" conclusion downstream is then a tautology.**

**WITHDRAWN — DO NOT PUBLISH OR CITE:** X1 "the 2-cells enter under an action" · X2 "scalarity
governs the gate" (the gate is written at rank one, abelian, scalar, at exactly zero — W-01's own
ruling, reproduced at `3.5e-17`) · X3 T1-2's criterion · X4 "the slot bought rank two" ·
X5 "two physical clocks."

**REDISCOVERY, ALL UNCITED IN 309,047 SEALED BYTES:** P-9 = the SBS condition verbatim (Korbicz) ·
P-5 = quantum Darwinism's `R_delta`, P-2/P-3/P-6 the einselection cluster (Zurek) · P-7 = Hepp 1972
disjointness, the directed system = the Coleman-Hepp model · **"P-1 and P-9 hold in the limit and
fail at every finite stage with an explicit reversing unitary" IS Bell's 1975 objection to Hepp —
the corpus's most celebrated result, and it is fifty years old** (and correct: re-derived honestly
from five vertices, *"a good sign about the machinery and a bad sign about the literature search"*) ·
the flat holonomy carrying all formation weight at strong coupling = discrete Aharonov-Bohm 1959 ·
one-qubit-per-cell on a non-factorizing algebra = Donnelly-Wall edge modes ·
Mahler / Lind-Schmidt-Ward / Lawton 1983.

**THE LEDGER IN ONE LINE:** *the corpus imports its mathematics by name and reinvents its physics
anonymously, and the reinventions are competent — which is why the honest publication is a short
mathematical note with six citations, not a foundations paper.*

---

## W-06 — THE AUDIT OF THE AUDITS: **THREE IMPORTS FOUND. THE SPINE WAS KILLED BY NONE OF THEM. DISPOSITION: STOP-FALLS-REBUILD.**

Commissioned at the principal's direction — *"let's make sure that we aren't importing machinery
that isn't suitable for the record level to destroy observations."* Four legs, eight refuters.
**Answer: IN PART — and the single largest correction is not an import at all.**

### THE THREE IMPORTS

**IMP-1 — THE VACUOUS-CONTROL DISQUALIFIER, AND IT IS THE ONE THAT CARRIED STOP.**
*"It could not have failed under variation, therefore it is void"* entered this record **lawfully**
at REGISTER:191 against one of S4's **CONTROLS**, and was then carried across a category boundary
at :414-416 onto a **THEOREM**. **On a control the norm is sound. On a theorem it is incoherent —
a proved statement cannot fail; that is what "theorem" means.** Applied to Theorem S3-0 it produced
the word "vacuous", and that word carried the STOP. It entered a governing clause with **no digest,
no file:line, no named ruling** — which under custody §1 makes it **flagged, not inherited**.

**IMP-2 — THE GAUGE NON-FACTORIZATION PREMISE**, imported from Donnelly-Wall/Giddings and pinned to
S2's Theorem 2 — **a theorem this register had ALREADY REFUTED two rows earlier at REGISTER:63**,
and which W-04's own ERR-3 disqualifies. *A refuted theorem cannot BE a correct one.* **Registrar
error: both halves of that contradiction are in this register and are corrected here.**

**IMP-3 — IDENTIFICATION BY HEADLINE PARAPHRASE.** Decoherence and Bell attributions matched on
words rather than papers. See the bibliography corrections below.

### AND THE PROCEDURAL VOID

W-04 pre-registered: *"IF THE SCALARITY TEST RETURNS 'S3-0 FAILS', THE RECOMMENDATION CONVERTS TO
STOP."* **W-05 reports the test returned NEGATIVE — "S3-0 SURVIVED" — and converts anyway**, on a
post-hoc ground built from IMP-1. **Replacing a pre-registered falsifier that came back negative
with a post-hoc one is precisely the abuse the experimental-design norm exists to forbid —
committed in the act of importing that norm.** No registered reopen condition fired.

### **BUT THE THING THAT DECIDED THE SPINE WAS NEVER IMPORTED — AND THE CORPUS ALREADY KNEW**

It was an **unledgered stipulation the corpus made about itself**: that loop transport means
**multiplying by the whole-circuit holonomy** `M_gamma`, rather than moving **edge by edge** via `T`
with `T^3 = M_gamma`. S3's CHOICE LEDGER C3 ledgers only *which circuit schedule*; the edge tick
appears in **no alternatives column in any sealed artifact.**

**THE CORPUS'S OWN SEALED AUDIT FOUND THIS AND WROTE IT DOWN AS COR-F, AT
`S3_THE_CROSSING_AUDIT_V001.md:794`, BEFORE W-03, W-04 AND W-05 RAN.** W-04 cited it in one
subordinate clause. W-05 reached past it. `grep -n "COR-" REGISTER_V001.md` returns **two hits in
491 lines, neither in the W-02 row that records S3-0's ruling.** Two independent legs of this audit,
starting from unrelated suspicions, walked straight into it.

> **THE CHAIN WAS NOT UNDER-ADVERSARIAL. IT WAS UNDER-READ.**

**This is a registrar failure, not a lane failure: I registered W-02 without carrying COR-F into it.**

### THE GAUGE-GROUP HYPOTHESIS: REFUTED, AND THE QUESTION RETIRED

The registrar's suspicion — that local `U(1)^V` is a continuum import and global U(1) the
record-native choice — is **WRONG, twice over.**
**(a) `U(1)^V` is derivable from S1's own bytes:** S1 §5 writes `f = da` and verifies `dd = 0`,
deploying K1's cochain complex; the gauge parameter is that complex's **0-cochain**, and U(1)-valued
0-cochains on K1 **ARE** `U(1)^V` by definition. Intrinsic to a bundle over a CW base, not a
discretisation of anything. **(b) It is UNIQUELY correct on K1's own count:** sweeping all **52**
partition subgroups of `U(1)^5`, **exactly ONE** yields S1 §4's invariant parameter count
`2 = b1 + #faces = E - rank(d1)` — the full local group. Global U(1) gives 6; the class group 4.
**And the escape route W-05 advertised is itself the import:** "global U(1)" smuggles in a
cross-vertex phase frame the complex does not contain — the same *"Let"* that FOUNDING_DESIGN's S2
clause exists to catch. **Across all 52 subgroups the invariant separation is exactly zero in
precisely the 4 that CONTAIN the transports, and 1.32-2.00 in the other 48: S3-0 IS TRUE IF AND
ONLY IF IT IS VACUOUS, AND EVERY EXIT FROM VACUITY IS AN EXIT INTO FALSITY.**

### WHAT IS RESTORED — A QUESTION, NOT A HEADLINE

**THE QUESTION S3-0 WAS BUILT TO CLOSE IS RESTORED TO OPEN: can the carrier's own gauge-invariant
algebra hold the record? THE ANSWER: IT SEES THE RECORD AND IT DOES NOT HOLD IT.**
S3 cited the whole gauge law (S1:59-63) and **implemented half of it** — it computed the fixed
algebra with the **connection held fixed**, while S1:63 (`a_e -> a_e + theta_v - theta_u`) *is* the
connection law. Under the FULL action, a **Wilson-line-dressed observable built only from S1's own
edge transports** is gauge-invariant to `4.45e-16` and separates the branches at exactly
**`3*sqrt(3)/10`** where S3's own test returns zero. **So S3:187-188's "indistinguishable by EVERY
gauge-invariant carrier observable" is FALSE, not vacuous** — and on W-05's own criterion, *a
theorem that fails is contentful*. **S3-0's load-bearing BROAD form (S3:16-22, :206-207), the form
that killed CHOICE LEDGER C1's alternatives, is RESTORED FROM VACUOUS TO FALSE.**

**BUT IT RECURS.** Over `n <= 4000` cells the dressed separation returns to `2.221e-16` and falls
below `1e-9` on **1000 of 4000 cells**. **The recurrence obstruction is the one thing in this corpus
that no demolition and no restoration has dented.**

**N4's MECHANISM RESTORED CORRECTED:** not *scalar* multiplication — **FIBRE-WISE-NESS**. Any
fibre-wise unitary lies in the local gauge group at every rank (`dim A^G = 5` at rank 1, 2 and 3,
verified against genuinely non-scalar and non-abelian transports). W-05 named the wrong variable —
**the same error as ERR-1, recurring one floor up.**

### THE CORPUS IS SMALLER AFTER THIS AUDIT, NOT LARGER

- **THE BELL 1975 ATTRIBUTION IS STRUCK AND REPLACED BY HEPP 1972** — two years older; the register's
  own phrase "reversing unitary" occurs **0 times** in the ten non-register artifacts; Bell's abstract
  draws the **opposite** moral.
- **AND BELL'S ACTUAL OBJECTION TRANSPLANTS INTO S3'S SYSTEM AND WORKS.** An explicit
  Hermitian-unitary `R_N`, one factor per slot, has branch cross-expectation of modulus **exactly
  1.000000000000 at every `N = 1..9`** while `|Omega_N|` falls from `4.11e-01` to `8.42e-04`. That is
  Bell's `lambda_n = c` with `c = 1`. **The record is undone at every finite stage by an observable
  of Bell's own type. This is new, and it cuts against the corpus.**
- **THE AHARONOV-BOHM DEMOTION FALLS ON ITS STATED SUBJECT** — it demoted an observation the corpus
  does not contain, and S4:672 asserts its **negation**.
- **NEW, UNCONTESTED: THE CROSSING CONTAINS NO MEASUREMENT.** S3's write map produces Schmidt
  spectrum `[1, 0]` — entanglement entropy `<= 1.6e-15` — on every superposition; the carrier is the
  same pure state in both branches (trace distance `0.0`). **Carrier and record are never correlated.
  The "record" records which of two counterfactual transports was applied — a CONTROL SETTING, not a
  property of any system.**
- **P-5 = quantum Darwinism and P-9 = the SBS condition SURVIVE.** The leg that struck them tested
  Korbicz's **macrofraction-indexed, asymptotic** condition **at single-slot grain**. Both refuters
  caught it; **the audit recorded the error against itself.**
- **ZERO OF NINE crossing properties are novel, and the corrected bibliography is WORSE:** not four
  rediscoveries spanning 1972-2021 but substantially **ONE body of work** (Hepp, Coleman-Hepp, Bell,
  Lanford-Ruelle) rediscovered whole.

### WHAT STAYS DEAD

Theorem S3-0 as stated · "the carrier cannot hold the record, therefore a crossing is necessary"
(dead four times over) · any hope that a different gauge group restores anything · both of W-05's
own test headlines · "has a field: no" · W-03's kill of "formation sees loop incidence" ·
**all four of W-04's errata, attacked and unbroken** · the novelty of the nine.

### **DISPOSITION: STOP-FALLS-REBUILD**

**STOP's OBJECT is correctly dead and the restorations make it DEADER** — a carrier-internal
observable that holds the branch distinction is one more reason nothing forced a crossing.
**STOP's GROUND has flipped sign** — a vacuous spine is a dead end; a false spine with an exhibited
counterexample is a construction site. **STOP's DISPOSITION was issued on an unpointered import
after the pre-registered falsifier came back negative.** The live route: the carrier's own dressed
gauge-invariant algebra over a **wedge-growth sequence** (`V = 4k+1`; `k = 1,6,31,156,781` gives
`V = 5,25,125,625,3125`, every consecutive pair dividing — **so S3's rejection of C1(b)/(c) on
`5∤9` and `5∤7` is an accident of two arithmetic progressions, not a theorem**). Zero adjoined
dimensions, at S2's zero-addition standard. **Cost, stated: the identification `C^25 ≅ C^5 (x) C^5`
is a labelling not derived from incidence, and the untouched recurrence obstruction may kill it.**

### CHAIN RELIABILITY

**Reliable at arithmetic. Unreliable at disposition.** Every demolition's numbers survived four
independent audits. But the demolitions' failure mode is **MISNAMING THE OPERATIVE VARIABLE** —
S2 said commutativity where it was scalarity; W-03 said the group where it was three things;
W-05 said scalarity where it is fibre-wise-ness — **each caught only by the next level, the third
only now, inside the one result W-05 offered for publication.** W-04 wrote the antidote
(*"discount an audit finding by the same factor as a build finding"*), W-05 quoted it, applied it
downward, and **did not apply it to itself. The chain wrote its own discount and did not take it.**
**FIVE consecutive layers have each been caught by the next, and the rate has not fallen. Discount
this layer too.**

---

## W-07 — THE LINEAGE-INDEPENDENT LANE, RUN AT LAST: **THE RECURRENCE OBSTRUCTION WAS MEASURED AT `ord(rho) = 4`. W-06's RESTORATION SURVIVES. W-06's COMPUTATIONAL CONTENT IS NOT ON DISK.**

The lane W-03 specified and nobody ran. All 9 commits of this repo carry `Co-Authored-By: Claude
Fable 5`; this lane is **Claude Opus 5**, so custody §4's shared-lineage caveat does not apply
between it and the corpus it audits. Sealed corrections read first, in full, before any computation:
S3 audit **COR-A…COR-L** including **COR-F** at `S3_THE_CROSSING_AUDIT_V001.md:794`; S2 audit
**COR-A…COR-H**; **S1 entire**; the erratum against W-02.

**WHERE THE PROOF IS.** `W07_RECURRENCE_ISOLATION_V001.md`
sha256 `0f84fa564b1324837f2187070db2a2b5184ad13f9e3ff64f0abd8480a4a24397`. Lane code, conventions, seeds, grid
and outputs: `LANE_W07_RECURRENCE_ISOLATION/`, `SEALS.sha256`, 15 files, all OK.

### W-06 REPRODUCED FIRST, INDEPENDENTLY — TWO OF ITS THREE FIGURES

The dressed observable, rebuilt from S1's own edge transports under the **full** gauge action
(S1:63 on the connection *and* on the section — the half S3 did not implement), is gauge-invariant
at `3.600e-16` (W-06: `4.45e-16`) and separates the branches at `0.384349931183` where every
diagonal S3-gauge-invariant observable returns `6.939e-18`. **S3:187-188 is FALSE, and W-06's
restoration of S3-0's broad form from VACUOUS to FALSE is CONFIRMED from outside the lineage.**
And W-06's recurrence figure reproduces exactly: **`1000 of 4000` cells below `1e-9`, min `6.7e-19`.**

### THE FINDING: **`W_F = -1`, `W_C = -i`, `<W_F, W_C> = Z_4`. S1's PUBLISHED CONNECTION HAS FINITE ORDER 4.**

`1000 = 4000/4` is `ord(rho)`, the order of the branch ratio in `U(1)`. **Isolating that one variable
— carrier, ready state, observable, dressing, `k`-range and code path all held fixed — the effect
disappears five times out of five:**

```
ord(rho) = 4  (S1 published)      1000 of 4000 below 1e-9,  min 6.729e-19
ord(rho) = inf  sqrt2/sqrt3          0 of 4000,             min 1.567e-04
ord(rho) = inf  random seeds 1,2,3   0 of 4000,             min 1.7e-05 .. 3.1e-04
```

Scaling to `K = 1e7`: on the published connection the dressed record is annihilated **to exact zero
on exactly `K/4` cells at every `K`, forever**; off it the separation is **never zero at any `k`**
and the worst near-return floor falls like `~2pi/K`. **THE OPERATIVE VARIABLE IS `ord(rho)`, FINITE
VERSUS INFINITE** — not "the connection", which is too coarse to be a finding.

### AND THE SAME CUT ONE LEVEL DOWN, ON THE FOUNDING OBSTRUCTION

`FOUNDING_DESIGN §4` promotes *"a finite discrete spectrum is recurrent"* to the obstruction the
construction must answer **at the start**. On S1's published connection `sup_k |Z_k| = 1` is
**ATTAINED**, exactly and periodically — `500000` of `10^6` cells at the published ready state,
`250000` at a generic one. On a generic connection and on S3/S4's headline `f=2.0,c=1.1` alike it is
`0.9999999999986` and **`0` cells**: **APPROACHED, never reached.** **COR-E saw half of this and
labelled the number a lower bound. Nothing in this register distinguishes ATTAINED from APPROACHED,
and they are different obstructions** — one is absolute, the other is a Diophantine rate question
nobody has asked. **The corpus has two distinguished connections — S1's, of order 4, and S3/S4's,
exactly resonant (erratum against W-02) — and both are arithmetically degenerate in different ways.**

### CUSTODY: **W-06 HAS NO ARTIFACT AND NO LANE CODE**

`4.45e-16`, `3*sqrt(3)/10`, `2.221e-16`, `1000 of 4000`, `4.11e-01`, `8.42e-04`, the 52-subgroup
sweep, the Bell `R_N` transplant, the Schmidt `[1,0]`, the wedge sequence — **every figure W-06
produced traces to exactly one file: this register.** W-01/W-02 left sealed build+audit artifacts;
W-03/W-05 left ten lane directories. **W-06 left nothing.** By its own IMP-1 test — *"no digest, no
file:line, no named ruling makes it flagged, not inherited"* — **W-06's computational content is
flagged, not inherited, and the disposition of record rests on it.** Mitigating and recorded: two of
its three figures reproduce exactly from an independent reconstruction, which is not what one expects
of numbers that were never computed. **`3*sqrt(3)/10` does not reproduce** — a factor `sqrt(3)` needs
an element of order 3 and `Z_4` has none. Unresolvable while the code does not exist.

### THE 52-PARTITION SWEEP IS A VACUOUS CONTROL

Re-run and the figure is right — `Bell(5)=52`, distribution `{2:1, 3:10, 4:25, 5:15, 6:1}`, unique
winner the discrete partition. **And it could not have come out otherwise:** a `k`-block subgroup
acts through at most `k-1` parameters, so `invariants >= 7-k`, so `2` forces `k=5`. One line, no
sweep. **"Could not have failed" voids a CONTROL** — and this is a control. **W-06's leg (b) is
void; its leg (a) — `U(1)^V` is the 0-cochain group of the complex S1 §5 already deploys in
`f = da` — is an argument, not a sweep, and SURVIVES. The gauge-group hypothesis stays refuted.**

### WHAT THIS DOES NOT ESTABLISH — STATED BEFORE THE VERDICT, NOT AFTER

**No durable record is shown.** Generically `|Z_k|` still returns to within `1.4e-12` of complete
indistinguishability by `k = 10^6`, which may erase a record as thoroughly as exact recurrence.
**The reading is TWO-WAY** — either the founding obstruction is an artefact of a degenerate
connection, **or** it survives generically in quantitative form and only the register's phrasing is
defective. **This lane does not distinguish them and does not score itself as having.** The dressed
observable is a **reconstruction**, not W-06's code. **And this is layer six**: five layers have each
been caught by the next and the rate has not fallen. **Discount this one.** Its likeliest failure
mode is the corpus's recurring one — misnaming the operative variable. It is named `ord(rho)` here.

**GRADE: INDEPENDENTLY-DERIVED (Opus 5 against a Fable 5 corpus), NOT ADVERSARIALLY-CHECKED.** No
refuter has read it. That is the inverse of every other grade in this corpus, and the missing half
is the weaker one.

**NOTHING HERE TAKES ANY OF THE PRINCIPAL'S THREE DECISIONS.** It changes the **shape** of the first:
STOP-FALLS-REBUILD was issued with the recurrence obstruction standing as *"the one thing no
demolition and no restoration has dented"*, and its rebuild cost line reads *"the untouched
recurrence obstruction may kill it."* Both are measured at `ord(rho) = 4`.

**REOPENS / NEXT, AND IT IS ONE CALCULATION, NOT A REBUILD:** `lambda = m(p00 + p10 x + p01 y +
p11 xy)` is the entropy of an algebraic `Z^2`-action; the near-return depth of `|rho^k - 1|` is
governed by the same arithmetic. **Does the record's growth rate outrun the near-return floor?**
Both sides are already computed objects here, and it decides whether the founding obstruction
survives off the degenerate point. **W-07 is superseded if:** a refuter shows the dressed
reconstruction is not W-06's object · `ord(rho)` is shown not to be the operative variable ·
or W-06's lane code is produced and contradicts §2.

---

## ERRATUM AGAINST W-07 — TWO SENTENCES OF THAT ROW ARE FALSE (registrar error, 2026-08-16)

Both found by W-08's refuters, both re-verified by me from the corpus's own bytes before being
carried here (`LANE_W08_REGISTRAR_VERIFY/`, 5 files sealed).

**E-1. W-07's HEADLINE IS FALSE OF THE FIGURES IT NAMES.** W-07 §1 states that *"every recurrence
figure in this corpus was measured"* at S1's order-4 connection. **The corpus states its own test
point verbatim, twice, and it is a different connection:** `S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md:439-440`
and `S3_THE_CROSSING_V001.md:423-427` both give `f = 2.0, c = 1.1, p = (0.4,0.15,0.15,0.15,0.15)`.
Recomputed independently: `min|Z_k| = 0.024654 at k = 42` and `sup|Z_k| = 0.999941 at k = 377` —
**W-01's registered `0.0247 at n=42` and `0.99994`, the figures STOP-FALLS-REBUILD was issued
against.** That connection is the **exactly resonant** one of the erratum against W-02
(`-11f + 20c = 0`), of **INFINITE order**, where the supremum is APPROACHED and never attained.
At the order-4 point `|Z_k|` takes exactly three values `{1/sqrt(10), 2/5, 1}` and `0.024654` is
`0.2916` from the nearest. **W-07's claim is true of exactly one figure — W-06's `1000 of 4000` —
and W-06 is the lane with no artifact.** The ATTAINED/APPROACHED distinction is not withdrawn; the
claim about where the corpus measured is.

**E-2. "A FACTOR `sqrt(3)` NEEDS AN ELEMENT OF ORDER 3 AND `Z_4` HAS NONE" IS FALSE.** Exactly, in
rationals: with `|s_v|^2 = (3/4, 4/25, 0, 9/100, 0)` (sums to 1) the pair `(v0,v3)` has
`(dF,dC) = (-1,0)`, so `D_k = amp|(-1)^k - 1|`, `amp^2 = 27/400`, `(2·amp)^2 = 27/100 = (3√3/10)^2`
**exactly**. The `sqrt(3)` is in the ready state's amplitudes, not in the group. **`3√3/10` is
reachable on S1's own published connection.** The W-07 *page* listed "a different normalisation"
among three escapes; the *row* did not, and the row is what is corrected. **The custody finding is
untouched — reachability is not reproduction, and W-06's code still does not exist.**

---

## W-08 — THE RACE: **THE DECAY OUTRUNS THE FLOOR, AND THE TWO NEVER COMPETE. THE FOUNDING OBSTRUCTION IS FALSE AS AN INFERENCE. WHAT SURVIVES OF IT IS A SCHEDULE STATEMENT WITH AN EXPONENT.**

Commissioned on W-07's REOPENS clause. Four build lanes (identification, floor, zero set, and a
refuter aimed at W-07), eight refuters on distinct lenses, an isolation audit over **64 comparisons**,
and a synthesis required to commit. **14 agents, 0 errors.** Every finding carried below was
re-verified by the registrar from the corpus's own bytes before entry.

**WHERE THE PROOF IS.** `W08_THE_RACE_V001.md` sha256 `3a9fa0ebfe88b3f19562e6251132fd11faebddca62d30d295f693352bf133e79`. Lane code:
fourteen directories `LANE_W08_*`, each with `SEALS.sha256`, **all verify**.

### THE ANSWER: **YES, UNCONDITIONALLY — AND THE QUESTION COMPARES INCOMMENSURABLES**

**`|Z_k| <= 1` always**, by the triangle inequality on non-negative weights summing to one — so
**`|Omega_N| = prod_{k<=N} |Z_k|` is MONOTONE NON-INCREASING.** Checked, not asserted: over
`k <= 10^6`, `max_k(|Z_k| - 1) = -1.446e-09` at the resonant connection, `0` exceedances anywhere.
**A near-return means that circuit writes nothing. It does not un-write the previous `N-1`.**

> **THE FOUNDING OBSTRUCTION IS FALSE AS AN INFERENCE.** *"It fires, and then it un-fires. A
> reversible write, not a record"* (`S2..._AUDIT_V001.md:447`, repeated `S3_THE_CROSSING_V001.md:429`)
> is drawn from the **single-cell** observable `|Z_k|` and never from the product, **which is the
> object durability is a property of.** W-07 looked for the fault in *where* the obstruction was
> measured. **The fault is in *what* was measured.** The corpus already half-owned this: COR-C
> struck `|Omega_N| <= e^{lambda N}` and replaced it with *"the monotonicity … that §4.3-4.4
> actually prove"* — **the monotonicity was in the record and was never turned on the recurrence
> sentence.** Under-read again, and this time the correction was the registrar's own to make.

**THE DECAY IS LINEAR WITH AN EXPLICIT CONSTANT, WITH NO DIOPHANTINE HYPOTHESIS.** From the
character identity `|Z_k|^2 = 1 - sum_{j<l} w_j w_l |chi_j^k - chi_l^k|^2` (exact `0` residual in
`Fraction` arithmetic over 768 and 336 cases, two independent refuters): for any pair with
`chi_j != chi_l`, `SUM_{k<=K}(1-|Z_k|) >= w_j w_l (K - 1/|sin(tau/2)|)`. Measured densities at
`K = 10^7`, one variable moving: `0.4919` (order-4), `0.4692` (resonant), `0.5295` (rank 1),
`0.4692` (random) — every one far above the proved floor `0.12`.

**THE FLOOR'S ENTIRE CONTRIBUTION IS `O(1)` AND `K`-INDEPENDENT.** At `K = 10^7` the whole
near-return band below `eps = 1e-3` contributes `8.9e-07` of the decay budget; below `1e-6`,
`9.2e-13`. **And the floor can be driven five orders deeper without moving the rate at all** —
approaching the order-4 point, `F(10^6)` falls `2.189e-06 -> 3.226e-11` while `lambda` reads
`-0.767507880` in every row. **THE FLOOR AND THE RATE ARE INDEPENDENT COORDINATES. THEY DO NOT
RACE.** The floor is a property of ONE cell; the rate is a property of ALL of them.

### WHAT ACTUALLY DEFEATS DURABILITY — AND RECURRENCE IS NEITHER

**(1) `G = {1}`.** `|Omega_N| -> 0` **iff** `G = <chi_a/chi_b : a,b in supp(pi)> != {1}` — proved by
Weyl on the *continuous* function `1-|P|` plus the strict triangle inequality, no Diophantine
input. This is **COR-B's four non-forming families re-derived as a criterion**, and it survived all
four refuters who attacked it.

**(2) THE SCHEDULE — AND THIS IS WHERE THE FOUNDING OBSTRUCTION GENUINELY SURVIVES.** Verified by
the registrar: an adversary writing only the `sqrt(K)` cells of smallest `1-|Z_k|` accumulates
`0.606, 0.615, 0.588, 0.601` nats at `K = 10^4..10^7` — **`|Omega| ~ 0.55` forever, with
unboundedly many writes** — while the honest schedule `k_n = n` accumulates `7678 -> 7675061`.
**Durability is a property of the (connection, SCHEDULE) pair, and the corpus has never stated a
schedule stipulation.** The admissible adversarial write density is `K^{-1/2}` at `d_eff = 2` and
`K^{-1/3}` at `d_eff = 1` (measured last-decade exponents `0.5003, 0.5003, 0.6666, 0.6667` against
theory `0.500, 0.500, 0.667, 0.667`).

**AND THIS IS WHERE W-07's ATTAINED/APPROACHED CUT DOES ITS WORK — THE SCHEDULE AXIS, AND NO OTHER.**
At the order-4 point `1-|Z_k| = 0` exactly on `250000` of `10^6` cells, so the adversary has a
**fixed positive** write density at **zero** budget. Off it the density must fall like `K^{-1/2}`
or `K^{-1/3}`. **A difference of EXPONENT — not the "absolute versus cosmetic" dichotomy.**
**W-07 was right that the register distinguished neither, and wrong about which axis it lives on.**

**AND ONE NON-UNIFORMITY, NEW AND IN NO PRIOR ROW:** `lambda < 0` at every non-trivial connection,
but the **onset** diverges — the circuits needed to write one nat run `2, 10, 47, 216, 1000, 4642`
as `t -> 1e-5` on the rank-1 locus, with `lambda = log 0.3` in every row (`K_0 ~ t^{-2/3}`).
**Durability holds for every non-trivial connection and holds uniformly for none.**

### RULING ON W-07

**SURVIVES:** the reproduction of W-06 (a third independent confirmation, `D_1 = 0.384349931183`
to `7.85e-14` from code sharing nothing with W-07's; `1000 of 4000` exact in 60-digit rationals) ·
the **vacuity ruling on the 52-partition sweep**, whose attacker's counterexample was itself wrong
(the discrete partition is the unique winner on **178 of 178** random connected graphs) · the
ATTAINED/APPROACHED observation, relocated · and its own limitation block, which W-08 vindicates.
**FALLS:** `ord(rho)` **as the operative variable** — an exact rational `theta = 1/4 + 1e-13` has
**finite** order `10^13`, so W-07's own sharp form predicts `0` and the observable returns `1000`;
the replacement name `min_k||k theta||` fails too, and the correct name at fixed `(K, tol, amp)` is
the pair `(q, delta)`. **W-07 is superseded on its own clause** — *"superseded if `ord(rho)` is
shown not to be the operative variable."* Plus the two errata above. **The corpus's recurring
defect, committed by the row that named it.**

### THE WEAKEST LOAD-BEARING CLAIM IN THE CORPUS — AND IT IS NOT IN W-08

**`W-01`'s convex-hull criterion as this register states it at `:43`, and specifically the clause
"three unit-modulus coefficients".** Registrar-verified, 200000 draws, one variable moving (the
number of occupied characters, same hull test and grid in both arms):

- **"three" is `p00 = 0`, an incidence fact about K1** — no vertex lies outside both loops. With
  one spectator vertex there are **four** characters, the criterion acquires the closed form
  `cos f + cos c <= 0` (agreeing with the hull on **200000 of 200000**), and **the firing region
  doubles from `0.2507` to `0.4989`** — exactly `1/4` to `1/2`.
- **W-01's own advertised virtue does not survive one added vertex.** *"It distinguishes curvature
  from flat holonomy, which K1 exists to separate"*: sending `f -> -f` alone changes the verdict at
  **99785 of 200000** points on K1 and at **0 of 200000** on the spectator carrier. **The property
  the register advertises as the criterion's virtue is a coincidence of `p00 = 0`.**
- The quantifier is dropped: `iff` is false read literally; `Z_1 = 0 => 0 in conv` is the true half
  and is the one W-01's exhibit travels on.

**No W-01 result is touched — and that is exactly why it has survived seven layers unexamined.**
It is quoted forward by S3, S4, W-03 and every W-08 lane; it is the sentence a rebuild would be
built on; and it is the only load-bearing claim in the corpus never checked on a second carrier by
a lane whose job it was.

### THE ISOLATION AUDIT'S FINDING, WHICH GENERALISES PAST THIS PROGRAM

**The commonest FATAL defect across twelve lanes is not "two variables moved" — it is "ZERO
variables moved."** Five of eleven FATALs are controls that could not have failed, one of them a
control in which the two arms were byte-identical (`difference exactly 0`) and which was reported
as a confirmation. **An isolation ledger cannot detect any of these, because a ledger records what
the author intended to vary.** Three build lanes wrote ledgers; five of the seven FATAL build
confounds sit in the arm the lane's own ledger names as decisive. **Treat "the ledger says one
variable" as evidence of nothing — read the code.**

### WHAT READS TWO WAYS, AND IS SCORED AS NEITHER

**The adversarial schedule** reads as *"durability is a schedule stipulation the corpus never
made"* **or** as *"`K^{-1/2}` write density is not an admissible schedule and the stipulation is
right."* W-08 supplies the exponent, not the admissibility criterion, and does not distinguish
them. **W-07 leg E's five-for-five null** is a density statement, not an isolation of `ord(rho)`.
**And `lambda` can be ~40% slower at degenerate points than generic** — which reads as "degeneracy
helps the record persist" or as "`lambda` is simply non-monotone in `H`"; the second is better
supported and neither is scored. **Every one of those rates is strongly negative, so the race
verdict does not depend on which reading is right.**

### WHAT W-08 DOES NOT ESTABLISH

**No durable record is constructed.** §1 decides a race between two computed objects on K1; it
builds nothing. **The schedule half has an exponent and no admissibility criterion. Nothing off K1
is decided** — all four build lanes ran on K1 alone. **The dressed observable is still a
reconstruction and W-06's code still does not exist.** **And this is layer NINE of one lineage:**
W-07, M1-M4, all eight refuters and the synthesis are Opus 5. Three of them named an operative
variable and the next caught it each time. **Discount this row as one block with W-07, not as an
independent check.** The lineage-independent lane W-03 specified still does not exist for anything
in W-08.

**NOTHING HERE TAKES ANY OF THE PRINCIPAL'S THREE DECISIONS.** It bears on the first: the
disposition's cost line — *"the untouched recurrence obstruction may kill it"* — names an
obstruction that is **false as an inference** and survives only as a schedule statement the corpus
never wrote down.

**REOPENS IF:** a connection with `G != {1}` is exhibited on which `SUM(1-|Z_k|)` is sublinear ·
an **intrinsic** admissibility criterion for schedules is written down under which the corpus's own
`SUM(1-z_n)` test is recovered and the `K^{-1/2}` adversary excluded · **or someone states a
durability observable on which a single-cell near-return DOES undo earlier writes — the one move
that would restore the founding obstruction, and nobody in this corpus has attempted it.**
**NEXT, AND IT IS ONE SCRIPT AGAINST TWO OBJECTS THE CORPUS ALREADY OWNS:** run W-01's criterion,
in this register's own words, on carriers `B1q` and `B0b` (`S4_THE_MEASUREMENT_V001.md:519, :582`)
and report which reading each obeys, quantifier restored. It decides the weakest claim above.

---

## W-09 — W-01's CRITERION OFF K1: **THE OPERATIVE VARIABLE IS ALL FOUR CLASSES OCCUPIED, WHICH NEEDS A PINCH *AND* A SPECTATOR. K1 HAS ONLY ONE OF THEM, AND SO DOES EVERY CARRIER THE CORPUS EVER RAN.**

The computation W-08 named as next. **One script, against carriers the corpus already owns** —
class multisets quoted verbatim from `S4_THE_MEASUREMENT_V001.md:575-590`; nothing new constructed.
`LANE_W09_W01_OFF_K1/`, `SEALS.sha256`, 8 files, all verify. Isolation: hull test, grid, seed and
evaluator identical in every row; **the one thing that moves is which characters the incidence
occupies.**

### THE RESULT

| carrier (S4's own row) | occupied | characters | firing region | `f -> -f` flips |
|---|---|---|---|---|
| `B1` K1 as handed, `B2`, `B1s`, `B3` | `{01,10,11}` | `v u uv` | **exactly 1/4** | 99785/200000 |
| **`B1q` K1-bridged + SPECTATOR** | `{00,01,10}` | `1 v u` | **exactly 1/4** | 99785/200000 |
| `B1p` K1-bridged | `{01,10}` | `v u` | `0` — never fires | 0 |
| **`B0b` ring torus loops meet, `B4` spindle** | **all four** | `1 v u uv` | **exactly 1/2** | **0/200000** |

**BOTH VALUES ARE EXACT, NOT MEASURED.** Three occupied classes reduce, on dividing by one
character, to `{1, e^{i th1}, e^{i th2}}` with the two angles independent uniform — so
**Wendel's theorem (1962)** gives `P(0 in hull) = 1 - 2^{-2}(C(2,0)+C(2,1)) = 1/4`. With all four,
`uv` is determined by `u` and `v`, Wendel does **not** apply, and the closed form does:
`0 in conv{1,u,v,uv} <=> cos f + cos c <= 0` (agrees on **400000 of 400000**), which is exactly
`1/2` because `(f,c) -> (pi-f, pi-c)` preserves the measure and flips the sign.

### AND THE NAMING — W-08's SYNTHESIS NAMED THIS WRONG, AND SO DID I

W-08 §6 called the operative variable *"three unit-modulus coefficients"*, i.e. `p00 = 0`.
**`B1q` refutes that: it has `p00 > 0` and reproduces K1 to the last digit.** The reason is in this
register already — **W-03's involution, `00 <-> 11` and `10 <-> 01`, exchanges K1's occupied set
`{01,10,11}` with `B1q`'s `{00,01,10}` exactly.** K1 has a vertex in **both** loops and none in
neither; `B1q` has one in **neither** and none in both. *Pinch-only and spectator-only are the same
carrier for this criterion, by a symmetry the corpus proved at W-03 and never applied here.*

> **CORRECTED NAME OF RECORD: the operative variable is whether the incidence occupies ALL FOUR
> CLASSES — which requires a vertex in both loops AND a vertex in neither. Not the coefficient
> count, not `p00`, not the pinch.** Every three-class carrier fires on exactly `1/4` and is
> curvature-aware; the four-class carriers fire on exactly `1/2` and are **curvature-blind**.
> **This is the fourth consecutive layer to misname the operative variable, and the third time the
> correct name was already in an earlier row.**

### WHAT FALLS

**W-01's ADVERTISED VIRTUE IS A THREE-CLASS ACCIDENT.** *"It distinguishes curvature from flat
holonomy, which K1 exists to separate"* (`REGISTER_V001.md:47`). Reversing the **curvature alone**
and leaving the flat holonomy untouched changes the verdict at `99785/200000` on every three-class
carrier and at **`0` of `200000`** on every four-class one — `cos f + cos c` is even in `f`.
**The property the register advertises as the criterion's virtue disappears the moment a carrier
carries both a pinch and a spectator, and no carrier the corpus ever ran carries both.** The ten
carriers of `S4:519` include exactly two four-class rows, `B0b` and `B4`; neither was run against
W-01's criterion by any lane before this one.

**NO W-01 VERDICT IS WITHDRAWN.** The criterion is correct on K1, the firing on S1's connection is
correct, and the `FIRE => HULL` direction is unconditional. What falls is its **scope** and its
**stated reason for existing**, both quoted forward unqualified by S3, S4, W-03 and every W-08 lane.

### REGISTRAR DEFECT, RECORDED RATHER THAN PATCHED

`w09c_exact.py`'s hull helper first read `np.diff(..., 0)`, where the second positional argument is
the difference **order**, not the axis — so leg C2 printed `0.000000` and a meaningless agreement
count on its first run. Legs A and B pass `axis=0` by keyword and were never affected; their
figures stand unchanged. The defect and its correction are in the sealed file.

**NOTHING HERE TAKES ANY OF THE PRINCIPAL'S THREE DECISIONS.**

**REOPENS IF:** a four-class carrier is exhibited on which the firing region is not `1/2` · or the
`f -> -f` blindness is shown to depend on the uniform measure on `(f,c)` rather than on the parity
of `cos f + cos c` · or a lineage-independent lane fails to reproduce the two exact values.
**NEXT:** every registered result derived on K1 that is quoted forward without a carrier
qualification — W-01's criterion, W-02's character-ratio criterion, W-03's multiset theorem, N1's
polynomial — should be re-stated with the class-occupancy condition attached. **The corpus has
never run a four-class carrier through anything.**

---

## ERRATUM AGAINST W-09 — "FOUR-CLASS CARRIER" IS NOT A THING (registrar error, 2026-08-16)

W-09 concluded *"of S4's ten carriers exactly two are four-class"* and *"the corpus has never run a
four-class carrier through anything."* **Class occupancy is a property of the LOOP DESIGNATION, not
of the complex.** The four classes are defined by *(is `v` in `gamma_F`? is `v` in `gamma_C`?)*, so
changing `gamma_C` changes them. Exhibited: on **B0b's own complex with S4's own `gamma_F`**,
sweeping every admissible simple cycle for `gamma_C` reaches **16 distinct class multisets** — 8
four-class, 7 three-class, 1 two-class — with `lambda` spanning `log(5/9)`, `log(4/9)`, `log(1/3)`
(`LANE_W10_A_CARRIERS_REFUTE_1/r1_rebuild.OUT.txt:81`). **CORRECTED TEXT OF RECORD:** *two
four-class **designations**, chosen under S4's CHOICE LEDGER C4, which closed the loop-designation
question by fiat.* Nothing in W-09's arithmetic changes — the `1/4` and `1/2` regions and the
curvature-blindness are properties of class occupancy and stand. **What falls is the phrase "the
corpus has exactly two four-class carriers", and with it the implicit claim that occupancy is
forced.** Deciding it needs an **admissibility criterion for loop designations**, which the corpus
has never written — the structural twin of W-08's missing schedule-admissibility criterion.

---

## W-10 — THE SCOPE TABLE: **THE CORPUS DOES CONTAIN CARRIER-INDEPENDENT RESULTS. THEY ARE THEOREMS, NOT MEASUREMENTS, AND THEIR CARRIER-INDEPENDENCE IS A CONSEQUENCE OF THE TRANSPORT CONVENTION W-06 NAMED AS HAVING DECIDED THE SPINE.**

Commissioned at the principal's direction, who **deferred the disposition decision by one step to
obtain it**. Four build lanes, eight refuters, a synthesis with its own verification lane, and the
registrar's own exact re-check of everything that corrects a register row. **13 agents, 0 errors.**

**WHERE THE PROOF IS.** `W10_SCOPE_TABLE_V001.md` sha256 `3b84120c895ea1024ac3de08cb6c5834a230610bdd0ed84c62dcc6917f9a90b1` — full table, ~50 rows, each with a SCOPE mark and
a **BASIS** tag distinguishing *theorem* from *exhibited* from *could-not-have-failed* from
*argued-only*. Fourteen sealed lane directories `LANE_W10_*`, all verify.

### THE ANSWER, IN THE FORM THE DECISION NEEDS

**CARRIER_INDEPENDENT:** N1's identification `lambda = m(p00 + p10 x + p01 y + p11 xy)` ·
W-08's monotonicity and the falsity of the founding obstruction as an inference · W-08's character
identity and its linear floor · W-02's criterion `FORMATION <=> G != {1}` · W-03/N2's multiset
theorem and its exact involution · the schedule result qualitatively · N3 and N4.
**Every one is a statement in `(pi, u, v)` with no carrier symbol in it — so a four-class run is a
CONTROL, not a test.** The table says so in the BASIS column rather than scoring the runs as
evidence.

**AND THE QUALIFICATION THAT MATTERS MOST, WHICH NO REGISTER ROW CARRIES:** they are
carrier-independent **given the transport convention** — whole-circuit scalar multiplication by
`M_gamma`. Under the corpus's own sealed alternative, **COR-F's edge-tick transport `T` with
`T^L = M_gamma`**, the transports are non-diagonal and `||[T_F, T_C]||` is `2.828` on B0b, `2.449`
on K1 and **exactly `0` on B0a — non-zero precisely when class 11 is occupied.** **The incidence
that N2 calls invisible is visible to the corpus's own sealed alternative transport, and COR-F
still has no register row.** Honest mark: **CARRIER_INDEPENDENT ∧ CONVENTION_SCOPED.**

**THREE_CLASS_SCOPED:** W-01's convex-hull criterion and its advertised virtue (W-09, now
reproduced a fourth time on a deterministic grid with a third hull algorithm) · W-02's published
support table, which is **incomplete, not wrong** — six two-element supports exist, three are
listed, and the three omitted are exactly those needing class 00.

**K1_SCOPED:** **the entire C\*-algebraic superstructure, and its hypothesis is the primality of 5**
— `UHF(5·2^inf)`, *"one qubit per cell, proved minimal"*, S3's rejection of its own adjunction-free
alternative, **and W-06's `V = 4k+1` rebuild route.** On B0b, `9 | 81`, so `M_9^{(x)N}` is a directed
system with unital embeddings and **no adjunction at all** — the route's own motivation evaporates
on the very carrier its sequence reaches.

**UNDETERMINED — thirteen rows, none guessed.** Six of the nine crossing properties (**unrunnable:
S3 published no lane directory**) · everything W-06 computed (no artifact) · W-08's four schedule
constants (not reproduced) · the multiset theorem's hypothesis off the non-negative locus.
**P-9 and P-3 are decided here for the first time and are CARRIER_INDEPENDENT.**

### NEW DEFECTS OF RECORD, REGISTRAR-VERIFIED

**N-4 — EVERY RATIONAL CONNECTION IS EXACTLY RESONANT, AND IT CONVICTS THREE MORE PUBLISHED ROWS.**
One line: if `f, c` are rational then `mf + nc` is rational, `2 pi j` is irrational for `j != 0`, so
`j = 0`, and `(m,n) = (r q, -p s)` is a nonzero relation. **Exact census, verified in `Fraction`:**
`f=2.0,c=1.1 -> (11,-20)` (already of record) · `f=2.0,c=2.0 -> (1,-1)` ·
**`S4:973`'s `3.14159, 1.57080`, published in terms as *"they are generic"*, `-> (157080,-314159)`**
· W-10 lane D's own hard-coded `1.3, 2.0 -> (20,-13)`. **S4:973's mislabelling sits inside the
paragraph correcting the corpus's first two mislabelled connections.** **`S4:603`'s `f = 1.0,
c = sqrt(2)` is the ONLY generic connection the corpus publishes — and it is the one S4 used to
verify the entire `lambda` column**, which is why that column survives.

**N-1 — S4's B4 ROW IS UNDER-DETERMINED BY ITS OWN PARAMETERS.** A second spindle — a triangle with
two 2-cells glued to a pentagon with two 2-cells at two points — matches **every published column**
(`V=6 E=8 F=4 chi=2 b=(1,1,2)`, gauge 5, inv 3, curv 2, flat 1, `d1·d2=0`, `gamma_F` bounds,
`gamma_C` does not, independent) and delivers class multiset `{00:1,10:1,01:2,11:2}` with
`lambda = log(1/3)` against S4's `log(1/2)`. **COR-K defect class against S4. And the build lane
hard-coded S4's multiset as its target, so its match could not have failed** — caught by its own
refuter.

**N-3 — S4:599's *"genuinely 4-term and does not factor. QUADRATURE ONLY, -0.810930216216"* IS HALF
WRONG.** The non-factoring half is exact (`p00·p11 = 8/81 != 2/81 = p10·p01`). The other half falls:
`|a+be^{it}|^2 - |c+de^{it}|^2 = 5/27 + (4/27)cos t`, which would need `cos t = -5/4` to vanish, so
**one Jensen branch dominates everywhere** and `lambda = m(4/9 + (2/9)x) = log(4/9)` **exactly**.
Independent Jensen quadrature at `n = 2^22` returns `-0.810930216216329`; branches never cross.
**Non-factoring does not imply no closed form.**

**N-6 — FIVE MORE ZERO-VARIABLE CONTROLS**, two of them inside blocks written to correct earlier
zero-variable controls, and an arms-diff guard that hashes **inputs** while the collapse is in the
**outputs**. **N-7 — the multiset theorem's hypothesis has now been misnamed SIX times, three of
them this round by refuters commissioned to catch misnaming, and my own "real non-negativity" is
one of the six.** Two independent mechanisms produce the invariance — flux coincidence and branch
domination — and no proposed name covers both. **Marked UNDETERMINED; the synthesis declined to
supply a seventh name and so do I.** It touches no registered result: `pi` is a probability vector,
so the corpus never leaves the regime where the theorem is proved.

**N-10 — "S4 WAS NEVER AUDITED" IS TOO STRONG, AND I REPEATED IT.** B1, B1q and B0b had their full
S4 rows — `d1`, `d2^T`, invariants, class multisets, `lambda` by three methods — reproduced in
`LANE_R_MAPS_REFUTER/rm_1_validate.OUT.txt` during the **W-03** round. The correct statement is that
**S4's artifact carries no sealed audit document and its carrier code is not in the repo**; parts of
its table were checked in passing by a lane doing something else.

### WHAT READS TWO WAYS — AND THE SYNTHESIS REFUSED TO SCORE IT, CORRECTLY

**IS "THE FORMATION FUNCTIONAL IS CARRIER-INDEPENDENT" A FINDING OR A RESTATEMENT OF THE
CONVENTION?** *Reading A:* a substantive result — the physics does not depend on the complex, which
is why a short mathematical note is publishable. *Reading B:* transport was **defined** as fibre-wise
scalar multiplication by the whole-circuit holonomy, so only `pi` **can** enter, and the result is an
analytic consequence of an unledgered stipulation — W-04's *"what K1 was: inert"* and W-06's *"the
thing that decided the spine"* both say so. **Both readings are supported. They have opposite
consequences for a rebuild: A says the functional layer is portable; B says portability is what you
get for free when the carrier was stipulated out.** Not scored. **This is the single most
consequential undecided question the table produces.**

### LINEAGE

**Layer TWELVE of one lineage.** W-07 through W-10, every lane and every refuter, are Opus 5. The
last lineage-independent boundary in this corpus is Fable 5 (S1-S4, W-01…W-06) to Opus 5 (W-07 on),
and everything in this round sits on one side of it. Two refuters accepted an upstream construction
byte-identical rather than re-deriving it. **Six consecutive layers have each been caught by the
next and the rate has not fallen. Discount this row as one block with W-07 through W-09.**

**NOTHING HERE TAKES ANY OF THE PRINCIPAL'S THREE DECISIONS.**

**REOPENS IF:** an admissibility criterion for loop designations is written (it would decide N-2 and
the B0b sweep) · COR-F's edge-tick transport is carried into a register row, which would rescope
every CARRIER_INDEPENDENT verdict · S3's lane code is produced, which would make six of the nine
crossing properties runnable · or a **lineage-independent** lane reads this table.
**NEXT, IF ANYTHING:** the convention question above is decidable — run the whole formation
functional under COR-F's `T` instead of `M_gamma` on K1 and on B0b, and see whether the incidence
becomes visible. It is the one experiment that would tell reading A from reading B.

---

## W-11 — THE CONVENTION RULING: **READING B. CARRIER-INDEPENDENCE IS A RESTATEMENT OF THE TRANSPORT CONVENTION. N1 IS NOT, SURVIVES INTACT, AND IS PUBLISHABLE UNDER TWO HYPOTHESES IT DOES NOT CURRENTLY CARRY.**

Commissioned at the principal's direction, who deferred the disposition decision a **second** time
to obtain it. The registrar's test (4 legs), then eleven agents sent to refute it: a uniqueness
attack on `T`, a clock-correspondence attack, a mathematics attack, a full-strength steelman of
Reading A, a blind independent rebuild, one cross-refuter each, and a ruling. **11 agents, 0
errors.** Thirteen sealed lane directories, all verify. Every load-bearing claim re-verified by the
registrar, **including the one that cuts against the registrar's own test.**

**WHERE THE PROOF IS.** `W11_CONVENTION_RULING_V001.md` sha256 `4a71ca885ade912dc672d99c2e8e3b82e3041a3e159527e36658efdc5e3cf027`.
Test: `LANE_W11_CONVENTION_TEST/`. Registrar's check: `LANE_W11_REGISTRAR_VERIFY/`.

### THE RULING, AND IT IS NOT A JUDGEMENT CALL

`W-01` **defines** `(M_gamma s)(v) = W(gamma)·s(v)` on the loop and `s(v)` off it. So `M_gamma`'s
diagonal **is the incidence indicator**, the four classes are the joint level sets of the two
operators' diagonals **by the definition of both**, and

```
Z_k = SUM_v (class character)_v |s_v|^2 = p00 + p10 u^k + p01 v^k + p11 (uv)^k
```

**is reached by substitution. No lemma intervenes. No other outcome was available.** That is what
"restatement" means, and it is prior to K1, prior to B0b, prior to any carrier.
**The corpus owns three-quarters of this and never assembled it:** W-03's *"an analytic fact about
the construction, true before any carrier was built"* (`:191`), W-06's *"not scalar multiplication —
FIBRE-WISE-NESS"* (`:577`), W-10's *"`M_gamma` is literally an element of the gauge group"* (`:409`).

**AND THE SENTENCE CONFLATES THREE CLAIMS OF THREE DIFFERENT STATUSES.** (1) the functional depends
on the state only through `pi` — **restatement**; (2) on the carrier only through `pi` —
**restatement**; (3) `lambda = m(P)` — **theorem with a hypothesis**. Reading A's own wording is
*"the physics does not depend on the complex, **which is why** one short note is publishable."*
**That "which is why" is a non-sequitur and is the whole of Reading A's apparent strength.**

### THE OPERATIVE VARIABLE, NAMED ON THE SEVENTH ATTEMPT — WITH THE GUARD RUN FIRST

Six names have failed in this program, three of them from refuters commissioned to catch misnaming
(W-10 N-7). The ruling grepped the register first and **attached a clause to an existing term
rather than coining a seventh.**

> **`|Z_n|` is a function of `pi` alone IFF the RELATIVE BRANCH OPERATOR
> `Q_n = (branch_F^n)^* (branch_C^n)`, at the tick the record is read, is MULTIPLICATION BY A
> FUNCTION OF THE INCIDENCE CLASS — fibre-wise (`REGISTER:577`) AND class-constant.**

Registrar-verified over 600 cells across four operator families: **AGREE 600, DISAGREE 0**, with
**316 blind and 284 not**, so the test could have failed. **The corpus's convention makes `Q_k`
class-constant diagonal BY CONSTRUCTION** — its diagonal is literally `conj(W_F)^{ka} W_C^{kb}`.
**That sentence is Reading B.** `Q` appears in no register row, and every result in the functional
layer is a statement about it.

### `T`'s CANONICITY — THE REGISTRAR'S OWN DECLARED WEAK POINT, DISSOLVED RATHER THAN DEFENDED

`T` is **not** canonical among all unitaries with `U^L = M_gamma`; that family is
positive-dimensional. **It is canonical under the corpus's own admissibility clause** — S2 audit
CHOICE LEDGER A1 (`:657`), whose clause (c) the ruling **strikes as circular** (it presupposes
fibre-wise-ness, the property at issue) and whose clauses (a) and (b), with locality from S1's own
definition of parallel transport (`S1:52-53`), give a four-line theorem: `U = Lam·T`, COR-F's `T`
times an `(L-1)`-torus of per-vertex phases the carrier does not supply.

**AND THE FINDING NEEDS NEITHER FACT.** The load-bearing result is a **biconditional quantified
over the whole family**, with `U^L = M_gamma` **dropped entirely** — 40 gradient descents over
`U(3)xU(3)`, 18 free parameters, 38/40 converged, positive control `5.49e-12`, **0 non-fibre-wise
solutions**. So the degradation if `T` is not canonical is: from *"the corpus's canonical rival
makes incidence visible"* to *"**every** admissible edge tick does, and the corpus's operator is the
one point of the family at which it does not."* **Logically stronger.**

**AND BOTH HORNS OF THE BRIEF'S DISJUNCTION — MINE — ARE FALSE AS WRITTEN.** The diagonal `L`-th
root `D` restores invisibility, but `D` is not an edge tick: it moves no fibre value and its
diagonal is the class indicator again — **the corpus's own convention at a finer clock**. And
"only the diagonal ones are blind" is false: **correlated non-diagonal pairs are `pi`-blind, 200 of
200, registrar-verified.** The blind set is the **correlated locus**, measure zero for independent
draws — **so every lane in this program that reported "0 of N random draws preserve invisibility"
ran a control that could not have failed.**

### THE COST TO READING B, CARRIED RATHER THAN BURIED

**COR-F's `T` FIRES AT THE TRIVIAL CONNECTION.** `FOUNDING_DESIGN:117-118` and `S2:583`
pre-register *"no formation at trivial connection"*. Registrar-verified at `a = 0`: under
`M_gamma`, `min|Z_n| = 1.000000000000` **exactly** — contact point met; under `T`, `0.429`. At zero
field the transports are pure cyclic shifts and the branches differ because the **paths** differ.
**The edge convention fails the corpus's own pre-registered contact point.** So the corpus has a
*reason* for its convention. **The ruling's distinction, which decides the question and which I
adopt: MOTIVATING A STIPULATION DOES NOT CONVERT ITS ANALYTIC CONSEQUENCES INTO FINDINGS.** Granting
the contact point in full buys the convention a **justification**, not its consequences a
**discovery** — `Z_k = SUM over classes` is still substitution.

### N1 — THE RESULT PROPOSED FOR PUBLICATION — **SURVIVES INTACT**

**N1 is NOT a restatement.** It is a theorem about `(pi, characters)`, and it **survives replacing
`M_gamma`**: under the fibre-wise root `D` — a finer clock, fractional winding, neither branch ever
closed — the identification holds at every tick to `2.5e-16` and the rate is `m(pi)` to `3.9e-07`.
A restatement of `M_gamma` would not survive replacing `M_gamma`. **PUBLISH IT — under two
hypotheses it does not currently carry.** *H1: the relative branch operator is multiplication by a
function of the incidence class — **a stipulation, not a theorem**, and the corpus's own sealed
COR-F exhibits an admissible alternative under which it fails.* *H2: `(conj(W_F), W_C)` generates a
dense subgroup of `T^2`* — without which the limit is the average over the proper closed subgroup,
differing from `m(P)` by `4.8e-04` at the resonant connection and `3.7e-02` at S1's own.
**H2 IS THE ERRATUM AGAINST W-02 (`:162-175`), AND N1 AS REGISTERED DOES NOT CARRY IT.**
**WITHDRAW from N1's framing:** *"the physics does not depend on the complex, which is why one
short note is publishable."* The note is publishable. **That is not why.**

### WHAT READS TWO WAYS, SCORED NEITHER WAY

**THE HOLONOMY READING** — if the corpus's object is *by nature* a closed-loop invariant, `T` is not
a rival at all and mid-path comparison is a category error. **No lane closed it and neither does the
ruling.** What the record shows is that the corpus is **not entitled to it for free: twice, the
warrant offered for the circuit clock is a fact about edges.** CHOICE LEDGER A2 (`S2 audit :658`)
warrants the circuit clock by citing *"edge count is carrier-supplied combinatorics (S1 :16-22)"* —
and `S1:16-22` **is the edge list**; S3's C3 (`:981`) cites COR-F, whose text reads *"edge traversals
**/** circuits"* — **a disjunction cited to close one disjunct**, the edge unit absent from its
alternatives column. All four quotations verified at the bytes. **The verdict survives either
reading, because under the holonomy reading `Z_k = SUM over classes` is still substitution.**

### LINEAGE, AND A DISAGREEMENT LEFT UNADJUDICATED

Every lane, every cross-refuter, the ruling and the registrar are Opus 5. **W-10 says layer twelve;
the brief said thirteen; all five cross-refutations say fourteen. Three numbers for one round, in a
program whose central discipline is applied by number. Recorded, not adjudicated.**

**NOTHING HERE TAKES ANY OF THE PRINCIPAL'S THREE DECISIONS.** It answers the question the second
deferral was for, and it attaches an answer to the second decision: **the Mahler note is
publishable, under H1 and H2 stated, with one framing sentence withdrawn.**

**REOPENS IF:** the holonomy reading is closed in either direction · an admissible operator is
exhibited whose relative operator is class-constant diagonal but which is not fibre-wise, or vice
versa (it would break the biconditional) · or a **lineage-independent** lane reads this row.

---

## W-12 — THE TOPOLOGICAL-CONSTRAINT ROUTE: **CLOSED, ON EVERY CARRIER. THE PROPOSAL WAS THE REGISTRAR'S AND IT IS REFUTED — AND THE REFUTATION UPGRADES N3 FROM ARGUED TO PROVED AND FIXES N1's SECOND HYPOTHESIS AS PURELY ARITHMETIC.**

Commissioned by the principal, 2026-08-17, on an idea borrowed from **DNA supercoiling**: the
linking number `Lk` is an **integer** topological invariant, conserved under every continuous
deformation and changed only by cutting, and the physics is driven by the mismatch it enforces
(`Lk = Tw + Wr`). **The registrar's proposal:** N3 kills every **absolutely continuous** connection
measure — a Wilson action is one — because the resonance structure is Haar-null; a **topological**
constraint is **singular** and so escapes N3's own wording. `S1:102-104` states that Bianchi is
vacuous on K1 because `b2 = 0`; `S4:511-515` owns five carriers with `b2 >= 1` on which it is not,
and no lane ever imposed it. `LANE_W12_BIANCHI/`, 5 files, sealed, all verify.

### THE PROPOSAL FAILS, AND THE FIRST REASON IS THAT BIANCHI IS NOT A CONSTRAINT AT ALL

For a 2-cycle `z`, `SUM_F z_F f(F) = a(boundary of z) = a(0) = 0`. **With the connection given by
edge phases — which is how `S1 §3` defines it — this holds IDENTICALLY for every `a`.** Verified on
B0b over 2000 random connections: `max |SUM_F z_F f(F)| = 1.69e-14`. **Bianchi restricts which
CURVATURE ASSIGNMENTS are realizable; it does not restrict which connections exist.** And it never
pins a single designated holonomy — `W_F` is one face curvature among nine and one relation among
nine quantities leaves any one of them free. `(W_F, W_C)` covers `T^2`: 20000 draws fill **400 of
400** grid cells.

### AND THE SECOND REASON IS A THEOREM THAT CLOSES THE ROUTE EVERYWHERE

> **The map `phi : (R/2piZ)^E -> T^2`, `a |-> (<gamma_F,a>, <gamma_C,a>)`, is a continuous
> homomorphism of compact connected groups, so its image is a CONNECTED CLOSED subgroup of `T^2` —
> hence `{1}`, a circle, or `T^2` — of dimension equal to the R-rank of the `2 x E` incidence
> matrix. For any two DISTINCT designated loops that rank is 2, because a simple cycle has
> `0, +-1` coefficients and two simple cycles are R-dependent only if they share support.
> **So the image is all of `T^2` and `phi` pushes Haar to Haar.**

Verified: rank 2 on every corpus carrier and designation, including an adversarial pair sharing two
of three edges; and the pushforward is uniform on `T^2` to sampling error (`chi2/dof = 1.007, 0.918,
0.914`; **576 of 576** cells occupied on each).

**COROLLARY 1 — NO CARRIER AND NO DESIGNATION CAN CONSTRAIN `(W_F, W_C)`.** Leg A's null is
structural, not a property of B0b. **The topological-constraint route is closed for every complex,
and it cannot be reopened by choosing a different carrier.**

**COROLLARY 2 — N1's HYPOTHESIS H2 IS PURELY ARITHMETIC.** It is a condition on the individual
connection and is **never violable by the carrier or by the loop designation.** The one degenerate
designation that confines the pair — `gamma_C = +- gamma_F`, which CHOICE LEDGER C4's *"any two
loops"* permits — is disposed of twice over: at `gamma_C = +gamma_F`, `uv = 1` and `Z_k = 1` for all
`k`, so `G = {1}` and **W-02's own criterion excludes it**; at `gamma_C = -gamma_F` the pair is
confined to the diagonal circle and H2 genuinely fails, **but only classes 00 and 11 are occupied,
so `P = p00 + p11 xy` depends on the product alone and the subtorus average equals `m(P)` exactly**
— deterministic trapezoid at `2^20`, difference `5.55e-17`. **Whenever formation occurs at all, H2
cannot fail for topological or designational reasons.** This is a strengthening of N1: its second
hypothesis has no structural failure mode.

**COROLLARY 3 — N3 UPGRADES FROM ARGUED TO PROVED, AND FROM K1-SCOPED TO CARRIER-INDEPENDENT.**
W-10's scope table marked N3 carrier-independent on **basis [A], argued only**. Since `phi` pushes
Haar to Haar on every carrier, **any** absolutely continuous measure on connections pushes forward
absolutely continuous on `T^2`, where the resonant set is Haar-null. **So no local action of any
form at any finite coupling can move the rate, on any carrier, for any designation.** Basis [A] -> [T].

### WHAT THIS DOES AND DOES NOT SAY

It does **not** say the borrowed idea was empty: the diagnosis it produced — that the corpus's three
missing admissibility criteria are missing because there is no energy function — is untouched.
**What is now closed is the specific escape hatch N3's wording appeared to leave open.** The wall
that N3 describes is not a feature of K1 and cannot be walked around by changing the carrier.
**A null here reads TWO WAYS and is scored as neither:** the topological route is closed, **or** the
registrar's formulation of it was the wrong one and a different singular constraint — one on the
*state* rather than on the connection — is untested. Nothing here touched the state side.

### REGISTRAR DEFECT, RECORDED RATHER THAN PATCHED

Leg B's diagonal-circle average was first computed by 400000 Monte Carlo draws and printed a
difference of `1.10e-03` — about two sigma of MC noise on a bounded smooth periodic integrand, and
**not evidence of a gap**. Replaced with a deterministic trapezoid, which converges spectrally
there; the difference is `5.55e-17`. The first version and its correction are in the sealed script.

**NOTHING HERE TAKES ANY OF THE PRINCIPAL'S THREE DECISIONS.**

**REOPENS IF:** a singular constraint on the **ready state** rather than on the connection is
exhibited that moves the rate · or a carrier is exhibited whose designated loops are R-dependent
without collapsing the occupied classes to `{00,11}`.

---

## W-14 — THE STATE SIDE, AND A RELEVANCE TEST FOR BORROWED MACHINERY: **THE STATE SIDE IS OPEN AND IS A KNOB, NOT A CHANNEL. THE INVENTORY OF WHAT THIS CONSTRUCTION CAN SEE IS NOW COMPLETE.**

Commissioned by the principal, 2026-08-17, on W-12's own REOPENS clause — *"a singular constraint on
the READY STATE rather than on the connection is untested"* — together with a standing instruction
that borrowed ideas be checked for relevance **at the record level** and **tested, not assumed**.
`LANE_W14_STATE_SIDE/`, 7 files, sealed, all verify.

### THE STATE SIDE IS STRUCTURALLY UNLIKE THE CONNECTION SIDE, AND IS GENUINELY OPEN

W-12 closed the connection side by a theorem: `a |-> (W_F, W_C)` is onto `T^2` with Haar
pushforward, so `lambda` is **constant almost everywhere** in the connection and no absolutely
continuous measure moves it. **On the state side there is no such constant.** Over 4000 uniform
draws on the 3-simplex, `lambda` ranges `-1.331196` to `-0.027776` — **a spread of `1.303420` nats
against exactly `0` across the whole connection torus.** And averaging does not collapse, because
`m` is not affine: `E[m(P(pi))] = -0.673126` against `m(P(E[pi])) = -1.380508`. **There is no N3 for
states, and there cannot be one.**

### BUT IT IS A KNOB, NOT A CHANNEL — AND THAT IS WHAT CLOSES IT

A constraint helps only if it lets the functional **see** something it could not see before. With
`pi` pinned exactly on **B0b, a four-class carrier**, and 60 states differing in **both** within-class
weight and phase:

```
  k        min |Z_k|            max |Z_k|            spread
  1   0.38232295728697   0.38232295728697          3.33e-16
  2   0.98022550160520   0.98022550160520          5.55e-16
 20   0.33432356245583   0.33432356245583          2.78e-16
```

**Exactly blind.** The state reaches the functional only through `pi`. **So a singular constraint on
the ready state selects a `pi` and cannot make the construction sensitive to anything new.** The
"wrong side" possibility is closed — for a different reason than the connection side: not because
the rate cannot move, but because nothing else gets in.

### THE ONE STATE RULE THAT IS NOT A STIPULATION, AND EXACTLY WHAT IT BUYS

Every ready state in this corpus was **chosen**. The exception is **SENSE U**, the uniform state,
where `pi` is fixed by the carrier's own class **sizes**. Under it `lambda` is a function of the
carrier alone — so the carrier does enter. **All seven recomputed values match S4's published
column exactly.** And the channel is narrow:

> **THE CARRIER ENTERS ONLY THROUGH THE MULTISET OF CLASS SIZES, AND S4's OWN TABLE CONTAINS THE
> EXHIBIT UNREMARKED.** `B3` (horn torus), `B1` (K1) and `B2` (K1 both triangles filled) share the
> class multiset `{01:2, 10:2, 11:1}` and differ in **every topological invariant on the page** —
> `chi = 1, 0, 1`; `b1 = 1, 1, 0`; `b2 = 1, 0, 0` — and give the **identical** rate
> `-0.756573585640`, **spread `0.000e+00`**. Nothing about `chi`, `b1`, `b2`, the pinch, the faces
> or the 2-cells survives.

### **THE RELEVANCE TEST FOR BORROWED MACHINERY** — the methodological result of this row

The functional sees **exactly three things: `pi`, `u`, `v`.** Every imported idea acts on some
variable. **Before importing, name the variable and check the functional can see it.**

| the import acts on | status | authority |
|---|---|---|
| the **connection** | **CLOSED.** The map is onto `T^2` with Haar pushforward on every carrier, so no absolutely continuous measure moves the rate | W-12, proved |
| the **ready state** | **OPEN but a KNOB.** It selects a `pi`; it is not a channel and cannot add sensitivity | W-14, exhibited on a four-class carrier |
| the **transport convention** | **THE ONLY LIVE CHANNEL.** Changing it is what makes the incidence visible — at the cost of firing at the trivial connection | W-11 |
| **anything else** | the functional cannot see it at all | the definition of `Z_k` |

**COROLLARY, AND IT SORTS THE TWO BORROWS THIS PROGRAM HAS TRIED.** The only import that can change
what this construction sees is one that changes the transport convention or **adds a dynamics**.
The DNA/topological-constraint borrow acted on the **connection** and was therefore closed before it
was run (W-12) — the relevance test would have predicted it. The kinetic-proofreading borrow adds a
**dynamics**, and is the right kind of import; it remains untested and its literature contact in
this corpus is **zero** (Landauer 0, Hopfield 0, proofreading 0, dissipation 0, entropy production 0,
free energy 0, Bennett 0, detailed balance 0, ratchet 0).

### WHAT READS TWO WAYS

**The narrowness of the carrier channel** reads as *"the construction is blind to the complex, so it
is not about geometry at all"* **or** as *"class-size multiset is exactly the right coarse-graining
and the rest is gauge"*. **Not scored.** W-11's ruling bears on it — the blindness is a consequence
of the convention — but which reading is right is a question about what the object is for.

**NOTHING HERE TAKES ANY OF THE PRINCIPAL'S THREE DECISIONS.**

**REOPENS IF:** an import is exhibited that acts on none of `pi`, `u`, `v` and still moves the rate
(it would refute the inventory) · or the uniform-state rule is shown to be a stipulation after all.

---

## ERRATUM AGAINST W-14 — ITS EVIDENCE WAS A CONTROL THAT COULD NOT HAVE FAILED (registrar error, 2026-08-17)

Raised by the principal, who asked what *"we know it by theorem rather than by exhaustion"* meant.
The audit that question forced is against me.

**W-14 leg B pinned `pi` on B0b and reported 60 states seen identically to `3e-16` as the finding
that the state side is a knob rather than a channel. THAT RUN COULD NOT HAVE FAILED.** `W-01`
defines `M_gamma` diagonal, so `Z_k = SUM_v (class character)_v |s_v|^2` **by substitution**, and
the state cannot enter except through class sums. The run was a check on my implementation of
B0b's class map, not evidence for the claim — **and it is exactly the defect class W-08's isolation
audit named as the commonest fatal one in this program, committed in the row that registers the
test for it.** Under the program's own rule, *"could not have failed" voids a CONTROL, never a
THEOREM*: **leg B is void; the claim stands on the theorem below.**

**THE PROOF W-14 SHOULD HAVE CARRIED, IN TWO LINES.** `Z_k = <s, Q_k s> = SUM_{u,v} conj(s_u)
Q_k[u,v] s_v`.
**(a)** If `Q` has any off-diagonal entry the sum contains a cross-vertex product `conj(s_u) s_v`
with `u != v`, which `pi` does not determine — hold `|s|` fixed and rotate the phase of `s_v` alone
and the term moves. So `pi`-dependence **forces `Q` diagonal**.
**(b)** For diagonal `Q`, `Z_k = SUM_v d_v |s_v|^2`, a function of the class sums **iff** `d_v` is
constant on each class. **QED.**
This also supplies **W-11's biconditional**, which that row verified over 9000 cells rather than
proving. **Both rows now rest on the same two lines.**

**AND THE PROOF PREDICTS WHERE A RUN WOULD SEE SOMETHING, WHICH IS THE TEST LEG B SHOULD HAVE BEEN**
(`LANE_W14_STATE_SIDE/w14_d_proofs.py`, sealed): everything held fixed, **one** off-diagonal entry
added **within a class** — spread over 200 phase draws goes from `2.220e-16` to `1.200e-01`.
**A live failure mode, found exactly where (a) says it is.**

**THE INVENTORY, RE-GRADED.** `connection` — **theorem** (W-12). `ready state` — **theorem** (the two
lines above), not W-14's run. `transport convention` — **theorem** (the same two lines), not W-11's
9000 cells. `anything else` — **definitional**: `Z_k` contains `s`, `W_F`, `W_C` and no other symbol.
**The inventory is closed by proof. W-14's verdict is unchanged and its stated ground is replaced.**

**WHY THIS MATTERS BEYOND THE ROW.** Exhaustion over a sample is structurally blind to measure-zero
sets, and **the interesting sets in this problem are measure-zero** — the resonant connections are
Haar-null (N3), and the `pi`-blind operator pairs are a correlated locus that W-11 found only after
five lanes reported *"0 of N random draws"* and were sampling the wrong set entirely. **This program
has been defeated by exhaustion twice. A closure claim backed by sampling is worth nothing here,
and I made one.**

---

## W-15 — THE REGISTRAR'S ANALYTIC SHAPE, MEASURED: **THE PHRASE TIC IS INHERITED. THE HABIT OF CONVERTING QUESTIONS INTO TWO-VALUED OBJECTS IS THE REGISTRAR'S, AT 2.5x AND 6x THE PREDECESSOR'S RATE — AND IT IS THE DEFECT THE PROGRAM'S FOUNDING QUESTION HAD.**

Raised by the principal, 2026-08-17, from a prose complaint that escalated into a claim about
reasoning: *"You're locked to an x/y way of thinking about almost everything."* Measured on this
register, split at W-07, the Fable 5 / Opus 5 boundary. `LANE_W15_BINARY_AUDIT/`, sealed.

| construction | Fable 5 /1k words | Opus 5 /1k words | ratio |
|---|---|---|---|
| `rather than` | 1.50 | 1.53 | 1.02 |
| `X, not Y` | 3.60 | 2.51 | **0.70** |
| `versus` / `vs` | 0.15 | 0.22 | 1.44 |
| **biconditional (`iff`)** | 0.30 | 0.76 | **2.53** |
| **`two ways` / `two readings`** | 0.30 | 1.75 | **5.85** |

**THE PHRASE-LEVEL TIC IS HOUSE STYLE AND IS NOT THE REGISTRAR'S** — on `X, not Y` the registrar is
**below** the predecessor. **What is the registrar's is structural: biconditionals at 2.5x and
two-way readings at nearly 6x.** The defect is not *talking* in x/y. It is **converting questions
into x/y**, because a two-valued question is decidable by a computation and an open one is not.
Every headline the registrar wrote has that shape: attained/approached (W-07), decay/floor (W-08),
three-class/four-class (W-09), carrier-independent/K1-scoped (W-10), Reading A/Reading B (W-11),
closed/open (W-12), knob/channel (W-14). **Sixteen "reads two ways, scored neither" verdicts in
eight rows is a template, not a judgement** — it wears the form of restraint while performing a
conversion.

### AND THE PROGRAM'S FOUNDING QUESTION HAS THE SAME DEFECT, WHICH IS WHY THIS IS A FINDING

*"A forced crossing — does a durable record live INSIDE the carrier, or must it be adjoined
OUTSIDE?"* **That is a dichotomy, and W-05 established it was malformed:** *"THE SLOT IS ALREADY
INSIDE THE CARRIER … S3 proved `dim R >= 4` by the two-non-parallel-vectors argument and then looked
for the `C^2` OUTSIDE instead of inside"* (`REGISTER:415-422`). **The corpus was spent on an
inside/outside binary that was never real.** So the instrument this registrar has been auditing the
program with carries the program's own defect — **the W-04 pattern (*"the audits committed the exact
defect they convicted the builds of"*) recurring one floor up, structurally rather than locally.**

### THE DISCRIMINATOR, SINCE "USE FEWER DICHOTOMIES" IS NOT A RULE

> **DID THE MATHEMATICS HAND ME THE TWO BRANCHES, OR DID I NAME THEM?**
> `diagonal / non-diagonal` was forced and has a proof under it. `knob / channel` was coined — the
> words are the registrar's and nothing guarantees the options are exhaustive. **Every coined pair
> is a place to look for a third thing.**

**Corroborating evidence from the corpus's own quality ordering:** the two most useful artifacts
this registrar produced are the ones that are **not** two-valued — W-10's **four**-valued scope
table with its BASIS tags, and W-14's relevance test, a **three**-way inventory over `pi`, `u`, `v`.
Both reached further than any of the oppositions.

**REOPENS / APPLIES TO:** every coined pair in W-07 through W-14 is now flagged for a third-option
check. **The first three to check, in order: `knob / channel` (W-14), `finding / restatement`
(W-10 §0.4.1, which W-11 in fact split THREE ways and scored better for it), and
`theorem / exhaustion` (which omits the middle case that carries most real evidence — exhibited
with a live failure mode).**

---

## ERRATUM AGAINST W-11 — ITS HEADLINE NAMES THE WRONG OBJECT, AND TWO LANES SAID SO AND WERE OVERRIDDEN (registrar error, 2026-08-17)

W-11 registered: *"CARRIER-INDEPENDENCE IS A RESTATEMENT OF THE TRANSPORT CONVENTION."*
**Two W-11 lanes returned `REFUTED_AS_POSED` against that framing, independently, and the synthesis
folded them into a seven-to-one majority for READING_B.** `LANE_W11_R_BLIND` (a blind rebuild that
never opened the registrar's code) and `LANE_W11_R_CLOCK` both wrote: **there is no transport
convention to restate.**

**VERIFIED BY THE REGISTRAR** (`LANE_W16_ADVANCE_LATTICE/`, sealed):
`||T_F^3 - M_dF|| = 6.18e-16`, `||T_C^3 - M_c|| = 4.34e-16`. **`M_gamma` lies inside the cyclic group
`T` generates — the corpus's operator is a POWER of COR-F's.** There is **one** transport and **two
schedules on it**: both conventions are rays in one family `Y(mF,mC) = <T_F^mF s, T_C^mC s>`, the
corpus's being `(L_F k, L_C k)` and COR-F's `(n,n)`.

**AND THE OPERATIVE VARIABLE IS THE ADVANCE LATTICE.** Over a `26x26` lattice with 64 same-`pi`
states, one variable moving: **invisible cells `81`, advance sublattice `L_F Z x L_C Z` `81`,
SET EQUALITY `True`**, max spread off the lattice `9.085e-01`. **Invisibility is a property of
WHERE THE RECORD IS READ — whether a branch is mid-loop — not of a choice of transport.**

**W-11's VERDICT SURVIVES; ITS NAME DOES NOT.** Carrier-independence is still not a finding — it
holds exactly on the advance lattice and fails everywhere else, and the corpus reads only there.
**Corrected text of record:** *carrier-independence is a consequence of reading the record only on
the advance lattice `L_F Z x L_C Z`, where both branch operators are simultaneously diagonal.*
**This is the SEVENTH consecutive misnaming of the operative variable in this program, and the
first one where the correct name was supplied, in writing, by lanes the registrar overruled.**

---

## W-16 — A PROCESS AGAINST THE BINARY, PROPOSED FOR THE PRINCIPAL'S SEAL: **THE ESCAPE HATCH WAS NOT MISSING. IT WAS OUTVOTED.**

Commissioned by the principal, 2026-08-17: *"Is there a process more than a log note that can be
used to free ourselves from this tyranny of the binary?"* **A log note will not do it, and this
corpus already proves that** — W-08's isolation audit found three build lanes wrote isolation
ledgers and it did not stop them, with five of seven fatal confounds sitting in the arm the lane's
own ledger named as decisive. **A check that requires self-report gets self-reported.**

### THE FAILURE, MEASURED

Across all workflow journals: **46 agent results, 30 carrying a verdict. `REFUTED_AS_POSED`
returned 2 times (7%). `UNDECIDABLE_AS_STATED`, present in every schema, returned 0 times.**
The mechanism: **the brief poses the binary, which RECRUITS.** Eight of ten lanes answer inside the
frame; the two who reject it are a minority the brief created; the synthesis reads verdicts as
votes; the majority wins; the register records it. **Frame-level findings are destroyed by a step
the registrar wrote.**

**AND A SCHEMA DEFECT UNDERNEATH IT.** `REFUTED_AS_POSED` was an ALTERNATIVE to answering. A lane
that wants to say *"B, and the question is also wrong"* **could not**. Both objectors had to
surrender their answer to raise the objection. **Frame-refusal was made to compete with answering.**

### THE FOUR GATES — mechanical, executed, greppable

1. **THE FRAME LANE RUNS FIRST AND ALONE.** Its brief carries the QUESTION and no candidate
   answers. Its only deliverable is the enumeration of options, each tagged **DERIVED** (with the
   citation that forces it) or **COINED** (a registrar word). Answer lanes are commissioned against
   ITS enumeration. This removes the binarization from the one point where it occurs.
2. **`frame_objection` BECOMES A SEPARATE FIELD, NOT A VERDICT.** A lane answers *and* objects.
   Removing the competition is most of the fix.
3. **A FRAME OBJECTION BLOCKS SYNTHESIS AND IS NEVER A VOTE.** The round re-poses and re-runs, or
   the register headline records the question as unresolved.
4. **AT REGISTRATION, MECHANICALLY:** grep the round's journal for frame objections; **no register
   row may close a question over a frame objection it does not name and address.** This is the
   pointer rule (custody §1) applied to dissent, and its absence is detectable by `grep` and not by
   the registrar's honesty.

### WHAT THIS DOES NOT FIX

**The registrar still writes the frame lane's brief.** Keeping candidate answers out of it is as far
out of the registrar's hands as this goes without a different model — and the lineage-independent
lane W-03 specified still does not exist for anything after W-06. **Whether these gates generalise
past this program is untested.** And gate 1 is itself a coined structure; per W-15's own
discriminator it is a place to look for a third thing.

**ADOPTION IS THE PRINCIPAL'S ACT, NOT THE REGISTRAR'S.** `FOUNDING_DESIGN_V001.md` reads *"PROPOSED,
NOT ADOPTED. Governing only at the principal's seal"*, and no principal act is on disk. **This row
proposes an amendment to custody §1 and does not make one.**

**REOPENS IF:** a round runs under all four gates and still destroys a frame objection · or a frame
lane is shown to inherit the registrar's shape through its brief despite carrying no answers.

---

## W-13 — N1's CONVERGENCE, SETTLED: **TRUE AT THE CORPUS'S OWN PUBLISHED CONNECTION BY BAKER'S THEOREM. FALSE ON A COMEAGER SET AND AT ALL FOUR PUBLISHED RESONANT PAIRS. THE IDENTIFICATION IS NOT NEW — THE THIRD TIME THIS PROGRAM HAS RE-DERIVED A NAMED RESULT.**

The one open mathematical question in the corpus and the ground of the principal's second standing
decision. Four build lanes, four refuters, a ruling that re-derives what they disagree on and
**decides against three of them.** 9 agents, 0 errors. `W13_N1_CONVERGENCE_V001.md` sha256 `a0d131762263bc78b42a7cd2dec45565bb7486991bd29bbf6e4d9ce9bce84097`;
lanes `LANE_W13_*`; registrar's own check at `LANE_W13_REGISTRAR_VERIFY/`.

### THE SHARP HYPOTHESIS — THREE CONDITIONS, AND THE THIRD IS THE ONE THAT DECIDES

`(H2)` `L(theta) = {0}` · **(D1)** a polynomial **homogeneous** Diophantine bound on
`(alpha, beta)` · **(D2)** a polynomial **inhomogeneous** bound on `dist((u^k,v^k), Z(P))`.
**(D1) and (D2) are independent and neither follows from H2.** W-11's H2 alone is **not
sufficient**. And **(D2) is a joint condition on the connection AND the ready state**, because
`Z(P)` is fixed by the state — so the hypothesis that decides N1 is not a condition on the
connection at all. Priority belongs to `LANE_W08_M1_IDENTIFICATION/M1_08_THEOREMS.txt` T2(c),
sealed 2026-08-16, which named it and which the register never carried.

### **N1 IS TRUE AT `f = 1.0, c = sqrt(2)` — UNCONDITIONALLY**

**(D1) and (D2) both follow from Baker's theorem on INHOMOGENEOUS linear forms in logarithms
(Baker III, Mathematika 14 (1967) 220-228) whenever the two loop ANGLES are algebraic reals with
irrational ratio and the ready state has algebraic weights.** `f = 1.0, c = sqrt(2)` is such a pair
— **W-10 N-4's "the only generic connection the corpus publishes", and the pair S4:603 used to
verify its entire `lambda` column.**

**Two lanes ruled the opposite and were overturned.** Lane C: *"N1 is not verified at any named
connection in this corpus and cannot be with present mathematics."* Lane L: *"no effective
Diophantine result covers `u = e^{-i}`, `v = e^{i sqrt2}` … Baker requires algebraic arguments."*
**Both false, on one inference:** *u and v are transcendental, therefore Baker does not apply.*
**The transcendence of `u` is not an obstruction — it is the hypothesis.** What enters the linear
form is the **angle** `f`, not `e^{if}`. **Lane L's refuter confirmed the error instead of catching
it.**

### AND IT IS FALSE ON A COMEAGER SET — MEASURE AND CATEGORY DISAGREE

**THEOREM (Baire).** There is a dense `G_delta` in `T^2` on which **H2 holds** and
`liminf_N A_N = -infinity`. Three independent proofs, one explicit ladder, verified over four
decades: perturbation falling like `1/(2k)` while the depth is held at exactly **3 nats per unit
`k`** (`-3.0627, -3.0081, -3.0004, -3.0000`). **The failure set is Lebesgue-null AND comeager.**
**The corpus's entire N3/W-12 apparatus is measure-theoretic and has no vocabulary for the second.**
**And no counterexample can have algebraic angles** — so the corpus's own pair is safe, and the
corpus's numerics could never have found the failure.

**THE THREE-WAY STATEMENT OF RECORD:** algebraic angles with irrational ratio — **TRUE** (Baker) ·
Lebesgue-a.e. connection — **TRUE** (Borel-Cantelli) · a comeager set — **FALSE** (Baire) · both
angles rational multiples of `pi`, which is **all four published resonant pairs** — **FALSE** ·
`pi` on a curve stratum, **which is S1's own registered ready state** — **OPEN, and it is a
one-variable Sudler product where the classical literature applies.**

### NOVELTY: **NOT NEW, FOR THE THIRD TIME**

*"The rate is a logarithmic Mahler measure"* is the normal state of affairs in this corner —
Lind-Schmidt-Ward 1990, Kenyon-Okounkov-Sheffield 2006, Lyons 2005. **Worse: `Omega_N` is the
transfer product of a scalar quasi-periodic cocycle and `lambda` is its LYAPUNOV EXPONENT, and
using Jensen and the Mahler measure to evaluate it is the opening move of Herman 1983** (*Comment.
Math. Helv.* **58**, 453-502). **Forty-three years old, in every survey of one-frequency
Schrödinger operators, and the word "Lyapunov" occurs twice in the whole sealed corpus, in no
register row and no lane directory.** After Hepp/Bell (W-04) and Cassaigne-Maillot, this is the
third. **Score the novelty claim weak and leaning negative.**

**AND S1's OWN REGISTERED READY STATE IS A 65-YEAR-OLD NAMED PROBLEM.** `p = (1/2,0,0,1/4,1/4)`
gives `pi = (0,0,1/2,1/2)`, so `P = (1/2) y (1+x)` and **`|Z_k| = |cos(pi k alpha)|` exactly**
— registrar-verified to `3.33e-16` over `k <= 2e5`, with the average approaching
`m(P) = -log 2 = -0.693147181`. **That is a Sudler product** (Sudler 1964; Erdos-Szekeres 1959;
Lubinsky 1999). **No lane looked it up.**

---

## ERRATUM AGAINST W-10 — N2 IS STATED WITHOUT ITS HYPOTHESIS, AND N1 AND N2 ARE ONE RESULT (registrar error, 2026-08-17)

`REGISTER:196-197` and `W10_SCOPE_TABLE_V001.md` §3.2 state W-03/N2's multiset theorem —
*"24 of 24 permutations invariant"* — **with no connection qualification**, and W-10 marked it
**CARRIER_INDEPENDENT**. **Registrar-verified, one variable moving (the permutation), same
connection, same estimator, `N = 2e6`, K1's registered weights:**

```
f=1, c=sqrt2   GENERIC     spread 3.670e-06   (Birkhoff noise; all -> m(P))
f=2.0, c=1.1   RESONANT    spread 7.247e-04   6 distinct values
f=pi, c=pi/2   ORDER 4     spread 5.579e-02   2 distinct values
```

**The `S_4` invariance is a theorem about `m(P)`, not about `lambda`.** It holds **exactly where
`lambda = m(P)`, i.e. exactly under W-13's hypotheses.** What survives at *every* connection is only
**W-03's own involution `00<->11, 10<->01`** — registrar-verified at `0.000e+00`, `0.000e+00`,
`1.110e-16` at the three connections above. **N1 and N2 are one result with one hypothesis.**

**AND THE CUSTODY FAILURE IS THE PROGRAM'S NAMED ONE.** This was sealed on disk in
`LANE_W10_B_MULTISET/b4_involution_labels.OUT.txt` leg E on 2026-08-16 — *"W-03's 'the incidence
labels are invisible' is a statement about the GENERIC connection, and the corpus computed almost
nothing at a generic connection"* — attacked and **strengthened** by its own refuter, and entered
in **no register row**, while W-10's scope table states the unqualified opposite one row away.
**Two of the four results this corpus proposes to publish were entangled, the entanglement was
found ten hours before W-13 opened, and it was under-read. Again.**

**REOPENS IF:** a counterexample with algebraic angles is exhibited (it would refute W-13's Baker
argument) · or the Sudler literature settles the curve stratum, which would close the one line
still open.

---

## W-17 — THE FRAME INSTRUMENT, VALIDATED BY RETRODICTION: **IT WORKS AT TWO FIFTHS OF ITS SIZE. AND ITS ONE CERTIFIED CATCH IS THAT THE SENTENCE WHICH LICENSED THE PROGRAM'S CENTRAL STAGE IS AN INVALID INFERENCE THAT SIXTEEN REGISTER ROWS NEVER EXAMINED.**

Directed by the principal, 2026-08-17: a process, applied at route decisions, that challenges the
binary frame — testing not only whether there are more options but **whether the binary is false**.
Instrument: `FRAME_CHALLENGE_V001.md`. Validation: **retrodiction under reading cutoffs** over four
route decisions, each lane forbidden to read the corpus past the point its decision was taken, plus
a cutoff auditor and an instrument auditor told to argue the instrument is ceremony. 7 agents,
0 errors. `W17_INSTRUMENT_VALIDATION_V001.md`, sealed. Reduced instrument: `FRAME_CHALLENGE_V002.md`.

### THE SUBSTANTIVE FINDING — **T1/F5 AT R1, AND IT IS NEW AGAINST THE WHOLE CORPUS**

`FOUNDING_DESIGN_V001.md:63-65`: *"An inductive limit of finite objects is not finite — **which is
precisely how it escapes recurrence**."* That sentence is the founding design's **single stated
escape** from its own §4 obstruction, and `:99-102` sends stage S3 out *"with the inductive-limit
template as the starting point."* **It is the sentence that licensed the program's central stage.**

**IT INFERS NON-RECURRENCE FROM NON-FINITENESS, AND NON-FINITENESS IS NECESSARY, NOT SUFFICIENT.**
Registrar-verified, needing no complex, no gauge group and no fibre — only a commensurate spectrum:

```
   dim = 10, 100, 10000, 1000000     E_n = n     |A(t = 2 pi)| = 1.000000000000000
   contrast, incommensurate spectra: max|A| over t in [1,1e4] = 0.9476 / 0.9618 / 0.9592
```

**An infinite commensurate spectrum recurs exactly, at every dimension including the limit.**
Two lines of arithmetic. **And it was never run:** across all **1696 lines** of this register,
W-01 through W-16, occurrences of **`inductive limit`: 0. `quasi-local`: 0.** Sixteen rows, eleven
adversarial, **and not one examined the escape route that licensed S3.** W-08's *"the founding
obstruction is false as an inference"* attacks a **different** defect — single-cell `|Z_k|` versus
the product — and neither restates the other. **This attacks the ESCAPE, not the obstruction.**

**AND IT IS THE ONE CATCH NO LEAK REACHES**, which matters because the round was contaminated (below):
the instrument's worked examples cover F2 and F4 only; **F5 has no worked example anywhere.**

### THE INSTRUMENT IS PART FAULT MACHINE, AND THE CONTROL THAT PROVED IT WAS OMITTED FROM THE ROUND

The ruling ran the control the round did not: **point the instrument at frames known to be sound.**

| frame | F1 | F2 | F3 | F4 | F5 | fired |
|---|---|---|---|---|---|---|
| **PARITY** (evens / odds) | FIRE | **FIRE** | ok | ok | ok | **2/5** |
| **`diagonal / non-diagonal`** on `U(3)` | FIRE | **FIRE** | **FIRE** | indeterminate | ok | **3.5/5** |

**The second control is W-15's own named paradigm of a binary the mathematics handed over** —
*"`diagonal / non-diagonal` was forced and has a proof under it"* (`REGISTER:1509-1511`).
**F2 fires on evens versus odds**: `f(n) = n+1` carries the odd arm exactly onto the even arm,
symmetric difference `2` of `100000`. **F2 fires on diagonality**: for a 3-cycle `P` and any diagonal
`D`, `(P·D)^3` is diagonal — off-diagonal norm of `U^3` **`0.000e+00`** over 20000 draws, on a
positive-dimensional slab. **A map exists between almost any two arms. F3** voids a predicate whose
cell is measure-zero, and `mu(diagonal) = 0` under Haar while the predicate partitions `U(3)`
exactly — **and this corpus's own interesting sets are all measure-zero** (N3's resonances, W-11's
correlated locus). **F4's verdict is chosen by the lane**: two defensible presupposition sentences
give opposite answers on the same frame.

**AND F1 AS WRITTEN IS TRIVIALLY SATISFIABLE** — two arms always lie in some space. The missing
clause is **EXCLUSION**: exhibit a point in **neither** arm. Parity: `0 of 100000`. Diagonality:
`0 of 20000`. **Both controls clear under it, and F1 keeps every productive result of the round.**

### **RULING: ADOPT TWO TESTS, NOT FIVE.** `FRAME_CHALLENGE_V002.md`

**KEEP `T1` (the null option), run first** — it produced the round's only certified catch at a cost
of two lines. **KEEP `T2` (embed), amended with EXCLUSION + RELEVANCE + PROVENANCE.**
**DELETE `F2` and `F3`. DEMOTE `F4`** to a grep gate: keep the requirement to *write* the
presupposition, drop the requirement to *score* it.

### REGISTRAR DEFECT — **I CONTAMINATED MY OWN EXPERIMENT**

`FRAME_CHALLENGE_V001.md:54-56` carries **R1's own F4 answer as a worked example**, citing a register
row far past that lane's cutoff — and the brief **ordered every lane to read the instrument first.**
The R1 lane detected this itself and refused to certify its own retrodiction, writing:
*"Availability is demonstrated. Discoverability is not."* **A second leak: the session auto-memory
line carried "forced crossing does not exist" — R1's answer — at turn zero.** **Worked examples in a
frame instrument are leak vectors.** V002 carries none.

### WHAT READS TWO WAYS

**Can a frame lane's brief be kept clean by THIS registrar at all?** The round says no on the first
attempt and the fix is mechanical (delete the examples), **or** it says the contamination is
structural because the registrar writes the brief. **Not scored.** The lineage-independent lane
W-03 specified still does not exist for anything after W-06, and **this is layer fourteen of one
block.**

**ADOPTION IS THE PRINCIPAL'S ACT.** `FOUNDING_DESIGN_V001.md` reads *"PROPOSED, NOT ADOPTED.
Governing only at the principal's seal"*, and no principal act is on disk.

**REOPENS IF:** a sound frame is exhibited on which `T1` or amended `T2` fires · or the founding
escape inference is shown to have been tested somewhere in the corpus after all.

---

## W-19 — ATTEMPT TWO, FIRST COMPUTATION: **NO CARRIER OF THIS CLASS WORKS AT ANY SIZE. THE PLATEAU IS A GAUSS IDENTITY. AND THE RESIDUE IS PHASE-INDEXED — ZERO BITS ON THE ELECTRIC RECORD, ONE FULL BIT ON THE MAGNETIC, AND THE COUPLING DECIDES WHICH.**

Commissioned under R-1 as amended. **First row entered under the gates adopted 2026-08-17**: frame
objections reported first and never as votes, and a next step named or the row is returned.
5 agents, 0 errors. `W19_CARRIER_THRESHOLD_V001.md` sha256 `abf83e5b971e9d394abb847e090b6029a32151d5f6fe110e96350faf4b7ef7dd`; lanes `LANE_W19_*`; registrar's own
check at `LANE_W19_REGISTRAR_VERIFY/`.

### THE FRAME OBJECTIONS, FIRST AND SEPARATELY — **BOTH SUSTAINED**

**The question named a minimum without naming the class it minimises over**, and the class does all
the work: degree-2 allowed gives 9 · min degree 3 with multi-edges 9 · simple 12 · girth >= 4 15 ·
girth >= 6 21 · nested-from-both-endpoints unreachable — **every one a correct answer to the question
as written.** And *"report a THRESHOLD, not an impression"* — **the registrar's own instruction** —
rewarded a number and penalised the finding that the number is not the object.

**THE MEASUREMENT THAT SETTLES IT.** Hold carrier, state, coupling, algebra, fragment count and
`delta` fixed; move **only the partition** into disjoint fragments:

```
heawood (L=21)   R_delta on the chosen BFS-cut partition : 5 of 5
                 over 200 random equal-size disjoint partitions : mean 0.400, max 2, ZERO in 127
```

**`R_delta` as computed is a property of the partition, and the partition was chosen — correctly and
knowingly — to be the one the Gauss law guarantees.**

### THE NUMBER, DELIVERED AND MARKED — AND LANE A's PROOF IS FALSE

`L = 12` simple (`tri_chain12`, V=8, cubic, girth 3), `L = 9` multigraph. **Lane A claimed 21
(Heawood) and marked it PROVED**; the proof conflates *girth through `l` >= 6* with *global girth
>= 5*, and cycles avoiding `l` are unconstrained. Both counterexamples were run **through lane A's
own code at lane A's own coupling.** Exhaustive: no simple graph on `V <= 7` with min degree 3
reaches `d = 5`; `V >= 8` forces `L >= 12`; **2520 of 19355 labelled cubic graphs on 8 vertices
attain it.** General floor `L_min = ceil(3(P+2)/2)`.

### **BUT UNDER A CRITERION THAT CAN FAIL, NO CARRIER OF THIS CLASS WORKS — AT ANY `L`**

Pin the criterion: **P1** gauge-invariant system algebra · **P2** fragments pairwise disjoint ·
**P3** `|F| <= |E|/2` · **P4** four points · **P5** *the plateau must be able to fail.*
**The obstruction is structural, not computational, and two halves exhaust the class:**

**(1) WHERE THE GAUSS LAW CARRIES THE PLATEAU, P5 FAILS IDENTICALLY. REGISTRAR-VERIFIED:**
if `F` contains a `u`-`v` cut of `G - l`, then `Z_l · prod_cut Z = I` on the physical sector —
**`0.000e+00`** — so `I(S:F) = H(S)` **exactly, for every state.** Measured on three independent
Haar physical states sharing nothing but the constraint: `I/H(S) = 1.000000` on all three.
The ruling's own arms: ground state at `g^2 = 0.50`, ground state at `g^2 = 3.00` **where
`H_elec(S) = 0.001283` — an EMPTY record** — and two Haar states, **all giving `1.000000 x5`.**
**"Could not have failed" voids a control, and this was being read as evidence.**

**(2) WHERE P5 CAN BE SATISFIED, THE PLATEAU HAS NO GAUGE CONTENT.** On `d = 1` carriers a plateau
is carried by the state and does fail — `R_delta = 5` for a constructed electric GHZ, **`0` for the
ground state at every coupling and for Haar.** And what it certifies is **a classical repetition code
in the X basis, whose curves are identical on `L` bare qubits with no Gauss law and no plaquettes.**
**Zero variables moved with respect to gauge invariance.**

### **THE FINDING: THE BOUNDARY-ALGEBRA DEPENDENCE IS REAL AND IT IS PHASE-INDEXED**

One object moves — the algebra assigned to the system link. `EXT` = full matrix algebra in the
extended Hilbert space; gauge-invariant = `alg{X_l}`, the only non-trivial gauge-invariant
subalgebra of a single link.

| record | `EXT/H(S)` | gauge-invariant`/H` | **difference** |
|---|---|---|---|
| **ELECTRIC**, non-enclosing `F` | `1.000000000` | `1.000000000` | **`0.000000000` bits, exactly, on every fragment** |
| **MAGNETIC** (the exact `g^2 -> 0` ground state) | `1.000000000` | **`0.000000000`** | **`1.000000000` bits. Full swing.** |

**Small `g^2` is where the ground state is magnetic and the algebra disagreement is one full bit.
Large `g^2` is where it is electric and the disagreement is exactly zero.** So *which boundary
algebra* is **not a free convention — it is undetermined precisely in the phase the dynamics selects
at small coupling**, and Gauss-forced in the other. **This is W-04's residue, located, with a
Hamiltonian and a coupling underneath it for the first time — and it is the only thing in this round
that required attempt two's genuine dynamics to say.**
Consequence recorded: the algebra choice also moves the plateau count, so the answers now on the
table are **8, 9, 12, 15, 21 and unreachable — all correct under different unstated conventions.**

### NEXT STEP — NAMED, AS THE STANDING RULE REQUIRES

**ONE RUN. `Z_2` GAUGE + STATIC CHARGE ON `tri_chain12`, MAGNETIC SYSTEM REGION, CRITERION PINNED IN
WRITING BEFORE ANY STATE IS BUILT.** `S` = the triangle `{1,2,3}` — the smallest region on a simple
min-degree-3 carrier whose gauge-invariant algebra is **non-abelian and carries a Wilson loop**.
Nine environment links, four disjoint fragments **declared before the state is computed**.
**The isolated variable is the charge sector and nothing else**: `G_v = +1` everywhere versus
`G_v = -1` at two vertices, arms byte-identical apart from the sign pattern — **and the diff printed**.
Physical dimension `2^13 = 8192` with matter against `2^22` on Heawood, a **512x saving**, so there
is no ceiling. **And the question is INVERTED: do not ask where the plateau appears — ask where it
can fail.**

**REOPENS IF:** a carrier in this class is exhibited satisfying P1-P5 · or the phase-indexing of the
algebra disagreement fails to reproduce on a second carrier.

---

## C-1 — **A CLAIM, NOT A FINDING: THE RECORD IS GAUGE-VARIANT.** (registrar, 2026-08-17)

**A new row type, entered at the principal's direction** — *"if we engage in fitting we find it, we
slap our own hands and get on with it; we can't tiptoe around the process if we hope to learn
something."* **Detection changes the calculus:** the sham-target test caught a live doubt in one
script, so the cost of a wrong claim has fallen and the optimal boldness has risen. **A CLAIM row is
a commitment stated flat so it can be shot at.** It is not evidence and must never be cited as one.

> **C-1. THE RECORD IS GAUGE-VARIANT. That is why every gauge-invariant measurement this program
> has made has come back either FORCED or FLAT.**

**THE EVIDENCE, ALL ALREADY ON DISK.**
- **W-19:** the redundancy plateau in the gauge-invariant channel is a **Gauss identity** —
  `I(S:F) = H(S)` for every state, coupling and history. Forced; carries no information.
- **W-19:** its only live control lay **entirely in the gauge-variant part**, recorded there as a
  defect of the lane rather than as a result.
- **FT-1** (`LANE_FT_SHAM/`, sealed): across the whole coupling sweep the gauge-invariant channel
  swings **exactly `0.000000`**, while gauge-variant descriptions swing `0.021611` to `0.990639`.
  Invariance checked, not asserted: `Z` commutes with every Gauss operator at `0.000e+00`, `X` fails
  at `1.6e+01`.
- **Attempt one:** `Z_k` moved, but the crossing it motivated adjoined a record slot **carrying no
  gauge action at all** — S3's own flag F6.
- **AND THE CORPUS FLAGGED IT ON DAY ONE AND NEVER READ IT.** **COR-J**, sealed in the S3 audit:
  *the premise "the record must be gauge-invariant" is undeclared, load-bearing, and applied
  asymmetrically — fatal inside the complex, vacuous on the added slots.* **Attempt one built a
  gauge-variant record, noticed the inconsistency, and excused it. Attempt two has been measuring
  the invariant channel and finding nothing there.** This is the program's signature failure —
  a sealed correction, unread, for the whole life of the program.

**WHAT IT WOULD MEAN.** A gauge-variant record lives in the **edge modes** — the degrees of freedom
that must be added to factorize a gauge theory's Hilbert space (Donnelly-Wall). That is **W-04's
located residue**, it is the principal's **boundary**, and it is TERRITORY_MAP region **R-III**.
**Four objects that have been circling each other are one object.**

**THE FALSIFIER, AND IT IS CHEAP.** Exhibit a **gauge-invariant** quantity carrying record structure
that is **not a constraint identity**. Every lane that has measured an invariant channel has had the
opportunity to produce one and none has. **One counterexample kills C-1.**

**FITTING DISCLOSURE, PER THE STANDING GUARD.** C-1 was assembled **after** the measurements, from
findings gathered for other purposes. It is exactly the shape of a fitted story. It is entered as a
CLAIM for that reason — **so that the fitting tests can be pointed at it** rather than at a
conclusion already dressed as a result. **What the registrar got wrong by hedging:** the
`0.000000` appeared in the registrar's own sham output and was reported as a limitation of the test.
**It is the result.**

**NEXT STEP.** W-20's ledger arms are already pointed at this without having been designed for it:
its ARM 1 removes the gauge structure and its ARM 3 varies the boundary's formation. **Add the
falsifier as a standing arm in every subsequent round: measure a gauge-invariant channel and try to
make it carry something that is not forced.**

---

## W-22 — THE FRAME CHALLENGE ON TWO INHERITED BINARIES: **BOTH ARE TWO-POINT SAMPLES, AND WHAT THEY EXCLUDE IS A VALUE GETTING ITS START.**

Commissioned by the principal, 2026-08-17: *"let's question all imported assumptions — variant /
invariant — allow / require — we shouldn't be afraid to ask — could we be looking at a point where
certain values are just getting their start."* First application of `FRAME_CHALLENGE_V002` to
inherited machinery. `LANE_FC_BINARIES/`, sealed.

### **VARIANT / INVARIANT — EXCLUSION CLAUSE FIRES.**

On an **open** carrier — six bulk vertices constrained, one dangling link whose far end is boundary
and carries no Gauss operator — the boundary charge `Q = Gauss(6)` is verified **central** against
every bulk constraint at `0.000e+00`, with `Q^2 = I`, giving exactly **two superselection sectors**.
All twenty single-link operators then sort by two **computed** predicates into **three** categories:

```
 INVARIANT (bulk observables)   10 operators
 VARIANT (genuinely junk)        9 operators
 SECTOR-CHANGING                 1 operator   -- X_9, on the dangling link
```

**`X_9` lies in NEITHER named arm.** It fails to commute with the bulk constraint, so the binary
files it under *variant* and discards it — **and it is the unique operator that moves between
superselection sectors, i.e. that creates or destroys boundary charge.** The classification throws
away the one thing that makes a value.

**AND THAT IS THE PRINCIPAL'S QUESTION IN THE ARITHMETIC.** A superselection sector is a value that
is or is not there, with no superposition between. `X_9` is the operator that starts one. **What
variant/invariant cannot see is precisely a value getting its start.**
*Caveat of record: `X_9` is unique because the construction has one dangling link — the COUNT is a
property of the construction. The three-way SPLIT is not; it follows from `Q` being central.*

### **ALLOW / REQUIRE — FAILS T1 AND T2, AND NO COMPUTATION IS NEEDED.**

**T1, the null option:** both branches presuppose the value **already exists**, to be permitted or
forced. If it **originates** at the surface, neither obtains and the question is malformed there.
**T2, embed:** the modal axis is not even two-valued — **forbid** is missing — so allow/require is a
two-point sample of a three-point axis. **And origination is not on that axis at all**: it is
generative, not modal. The binary samples two of three points on one axis while the live option sits
on a different axis entirely.

**WHAT THIS COSTS A REGISTERED INFERENCE.** The predecessor found alpha **allow-side** — permitted,
forced by nothing — and `FOUNDING_DESIGN §2` read that as the warrant for treating alpha as an
INPUT: *"a process that leaves the coupling's slot open is exactly what a process should look like
when the coupling comes from outside."* **Under the third option that inference does not go
through.** A slot that looks open because nothing forces the value is **indistinguishable** from a
slot where the value is being made, to any test that asks only whether something is forced.
**Allow-side and origination have the same signature under the test that produced the finding.**

### WHAT THIS DOES AND DOES NOT SAY

It does **not** show that any value originates anywhere. It shows that **two classifications this
program inherited cannot represent the possibility**, and that one registered inference depends on
the missing option not existing. **A null under a classification that cannot express the alternative
is not evidence about the alternative** — the same form as the standing correction on classical-shaped
nulls.

**NEXT STEP.** The sector-changing operator is now a named object on a carrier we own. **Ask whether
it is dynamical**: does any Hamiltonian built from the carrier's own data move amplitude between
superselection sectors, or is the boundary charge frozen by construction? A frozen charge is a label;
a charge the dynamics can change is a value with a history. **That is the difference between a
bookkeeping sector and a value getting its start, and it is computable on this carrier.**

---

## ERRATUM AGAINST W-22 — "TWO SUPERSELECTION SECTORS" IS FALSE ON THE PHYSICAL SUBSPACE, AND THE DISCRIMINATOR THAT FOLLOWED IT WAS VOID (registrar error, 2026-08-17)

Found by a diagnostic the registrar printed and nearly read past. Both defects are the same defect.

**E-1. THE CARRIER WAS TOO SMALL TO CARRY THE CLAIM.** W-22 reported `Q^2 = I` and inferred **two**
superselection sectors. On the **physical** subspace there is **one**. With a single dangling link,
every link with both ends in the bulk contributes `Z` twice and cancels, while the dangling link
contributes once — so `prod_bulk G_v = Z_9 = Q` **exactly** (`0.000e+00`), and imposing bulk Gauss
`= +1` **forces** `Q = +1`. Measured: eigenvalues of `Q` on the physical subspace `= [1.0]`.
**W-22's three-way operator split stands** — it rests on `Q` being central, which is independently
verified — **but the sector count does not.**

**E-2. AND THE DISCRIMINATOR BUILT ON IT COULD NOT HAVE FAILED.** The registrar then tested whether
`H` is a function of `Q` alone — the property separating the gravitational structure (where the
Hamiltonian **is** the boundary term) from the electromagnetic one (where charge is a surface term
and energy a bulk density) — and reported **EM-shaped**. **With `Q` constant on the physical
subspace, `H` cannot be a function of it, and the verdict was forced.** The run printed
`distinct eigenvalues of Q on the physical subspace : 1` in its own output. **VOID. The EM-versus-
gravity question is UNTESTED, not answered.**

**E-3. AND THE CORRECTED STRUCTURE IS BETTER THAN THE ONE CLAIMED.** With **two** dangling links
(physical dimension 32): `Q_6` takes **both** values `[-1.0, +1.0]` while the total `Q_6·Q_7` remains
forced at `[1.0]`. **The boundary charge is not free in total — it is free to be DISTRIBUTED, with
only its sum fixed by the bulk.** That is a different and richer object than a global label, and a
single dangling link cannot exhibit it.

**RECORDED WITH IT — THE INERTNESS RESULT STANDS AND IS NOT AFFECTED.** `[H, Q] = 0.000e+00` at
every coupling: the dangling link lies on **no cycle**, so no plaquette term can contain it, and in
**pure** gauge theory nothing in the Hamiltonian can move boundary charge. **The sector is frozen;
it is a label, not a process. Only matter can unfreeze it** — adding a term containing the
sector-changing operator gives `||[H,Q]|| = 4.480e+01`.

**AND A CORRECTION TO THE REGISTRAR'S LANGUAGE, AT THE PRINCIPAL'S DIRECTION.** The registrar wrote
that a negative result would "kill it". **It would not.** Under the charter (R-3) a failed test kills
**one assembly**, never the program — and this one did not even do that. **It killed one carrier that
was too small to show the structure**, which is the same failure as the theta graph at T1: the
question was sound, the construction could not carry it.

**NEXT STEP, AND IT IS NOW CHEAP AND SPECIFIED.** Re-run the EM-versus-gravity discriminator on the
**two-dangling-link** carrier where `Q` genuinely varies: does `H` restricted to the physical
subspace vary **within** a fixed `Q` sector? If it does not, the Hamiltonian is a boundary term and
the structure is gravity-shaped. If it does, energy is a bulk density and the structure is
EM-shaped. **The test is live for the first time, and the carrier exists.**

---

## W-24 — **OUR BOUNDARY WAS NEVER A BOUNDARY. IT WAS A FRINGE — AND A FRINGE CANNOT CARRY ANYTHING.**

Commissioned by the principal, 2026-08-17, after the observation *"lot of pairs"*. `LANE_W24_Z3/`,
sealed, 4 scripts.

### **THE STRUCTURAL FINDING, AND IT EXPLAINS A CLASS OF FAILURES AT ONCE**

**The boundary links of every carrier in this program lie on NO CYCLE.** Measured directly: links
0, 1, 2 on a cycle `True`; links 3 and 4 — the boundary — `False`.

> **No cycle means no flux can thread it, no persistent current can circulate on it, and nothing can
> flow there. By construction, at any group, at any coupling.**

**This is the common cause of results the register has been carrying as separate findings.** The
boundary charge was frozen (`[H,Q] = 0.000e+00`) because nothing in the bulk could drive it. The
chirality could not reach the surface because there was no channel. Every boundary measurement came
back **forced** because a fringe has no degrees of freedom of its own to be otherwise.

**THE BOUNDARY OF A REGION IS A CLOSED CURVE. OURS WAS A SET OF LOOSE ENDS.** That is what makes a
surface a surface, and every carrier here — K1, the theta graph, tri_chain12, Heawood, the open
graphs of the erratum against W-22 — has had dangling ends where it needed a cycle. **This is one
error, not five, and it has one fix.**

### ORIENTATION: WHAT Z_2 STRUCTURALLY CANNOT EXPRESS

**Z_2 forces every charge to be its own inverse**, so the two senses of a circuit are the same
object and the Wilson loop is its own adjoint. **`Im<W> = 0.000000` at every `theta`, structurally,
not numerically.** At **Z_3** the identical term with the identical coefficient gives
`-0.059329`, `-0.151602`, `-0.222120` as `theta` runs `0.4 -> 2.0`. **The ground state acquires a
sense of rotation.** `theta` is dimensionless, periodic, and sets no scale — the first slot of that
kind to appear in this program, and it appeared without being installed.
**AND IT IS CONFINED TO THE BULK:** every boundary observable reads `0.000000` at every `theta`,
for the reason above — the dangling links touch no plaquette.

### AND THE PAIRS ARE OURS, WITH THE SCOPE STATED

Tripling the group changed **nothing** in the measured structure: same three-way operator split
(4 invariant, 3 junk, 1 sector-changing), same Gauss identity at machine zero, same
`I(S:F)/H(S) = 1.000000`. **So the two-valuedness is the registrar's, not the carrier's** (W-15).
**Scope, stated:** every quantity compared there is one a constraint already fixes, so the
comparison could not have separated the groups. **No evidence the pairs come from `Z_2`; no strong
evidence either way.** What Z_3 *does* add is orientation, and that is a different axis.

### CORRECTIONS CARRIED

`-dE_0/dtheta` is **not** a chirality and is **larger at Z_2** (`-0.229`) than at Z_3 (`-0.065`):
`theta` reweights the plaquette in both groups. **`Im<W>` is the orientation-sensitive quantity;
the energy derivative is not.** Recorded because the registrar conflated them.
Five measurements in this stretch were of quantities a constraint already fixes — including
`<psi|i[H,Q]|psi>`, which vanishes in any eigenstate identically. **Standing rule adopted here:
before measuring any quantity, compute whether the constraint determines it. One commutator, first,
every time.**

### **NEXT STEP — ONE FIX FOR THE WHOLE CLASS**

**A carrier whose boundary is a CLOSED CURVE**: an interior region separated from an exterior by a
cycle of links, so that the boundary has flux through it, can carry a persistent current, and is a
surface in the sense this program has been assuming it had. **On that carrier, and only there, the
open questions become askable** — does the boundary carry a current; does the chirality reach it;
does the charge stay correlated with what wrote it; and is there a regime that is *allowed but held*.

**REOPENS IF:** a fringe is shown to carry a flux or a current after all · or the closed-curve
carrier reproduces the same forced results, which would move the cause elsewhere.

---

## W-25 — THE CLOSED-CURVE CARRIER: **THE FIX WORKS. THE BOUNDARY CARRIES CHIRALITY AND IS MORE THAN ITS PIECES. AND THE WHEEL CANNOT DISTINGUISH AREA FROM PERIMETER, WHICH IS THE NEXT CARRIER'S REQUIREMENT.**

Built at the principal's direction after W-24 identified that every prior carrier put its boundary on
dangling links. `LANE_W25_CLOSED_CURVE/`, sealed, 4 scripts.

### THE CARRIER, AND THE TEST EVERY PREVIOUS ONE FAILED

**Wheel `W_n`:** hub = interior, **rim = a closed curve = the boundary**, `n` spokes, `n` interior
triangles. **All rim links lie on cycles** — so a flux can thread the boundary and a current can
circulate on it, which no earlier carrier permitted. Physical dims `16 / 81 / 256` at `Z_2/Z_3/Z_4`.
**DISCRETE STOKES VERIFIED EXACTLY AT EVERY GROUP:** `W_rim = T_0 T_1 T_2 T_3`, residual `0.00e+00`.
**The boundary flux equals the enclosed curvature** — the first relation in this program between a
boundary quantity and a bulk one.

### THE CHIRALITY REACHES THE BOUNDARY

`Im<W_rim>` is `-0.006997 / +0.002180 / -0.011434` at `Z_3` and **exactly `0.000000` at `Z_2` at
every theta, structurally** — the Wilson loop is its own adjoint there, so the two senses of the
circuit are one object. **On every fringe carrier this read zero because the boundary touched no
plaquette. The rim touches `n`.**

### THE BOUNDARY IS MORE THAN ITS PIECES — AND IT GROWS

The four `<T_k>` are **identical** (`+0.23272-0.06043j`), forced by the wheel's `n`-fold symmetry,
**so they cannot cancel by sign.** And `|<W_rim>|` **exceeds** the product of the individual
expectations in every case. Enhancement `|<W_rim>| / |prod_k <T_k>|`, at `theta = 1.0`:

```
   rim n :      3       4       5       6
   Z_2   :   2.17    2.66    3.34    4.23
   Z_3   :   2.13    2.60    3.24    4.09
```

**Monotone growth, ~`1.25^n` — a constant connected correlation per boundary link.**
**AND IT IS THE SAME AT BOTH GROUPS**, while the chirality is zero at one of them. **So the
enhancement is NOT a chirality effect; it is generic plaquette correlation and would appear with no
orientation at all.** Recorded so it is not later attributed to orientation.

### REGISTRAR ERROR CORRECTED — A NEW CLASS

The previous turn claimed `Im<W_rim>` was "20-30x smaller than a single triangle despite Stokes, so
they must nearly cancel." **False, and the reasoning was invalid.** Stokes is an **operator**
identity and says nothing about expectations; comparing a four-operator product to a single operator
is comparing `<x^4>` to `<x>`. **Not a forced measurement this time — INCOMMENSURABLE OBJECTS
COMPARED.** Both hypotheses it generated (screening; resolution floor) are dead: symmetry forbids
cancellation, and exact diagonalisation at dim 81 has no resolution floor.

### **THE NEXT CARRIER'S REQUIREMENT, AND IT IS STRUCTURAL**

**On a wheel the enclosed area and the perimeter are BOTH `n`.** They scale identically, so this
carrier **cannot distinguish an area law from a perimeter law** — and that distinction is the
boundary-versus-bulk question, the one gravitational signature phrased natively in information
language. **Same failure mode as the fringe: the question is sound, the carrier cannot carry it.**

> **NEXT: a carrier where AREA AND PERIMETER SCALE DIFFERENTLY** — a 2-D lattice patch, where a
> region of side `k` has area `~k^2` and perimeter `~k`. Only there can `<W>` be tested for an area
> law against a perimeter law, and only there is the boundary-versus-bulk scaling askable.

**REOPENS IF:** the enhancement is shown to be a finite-size artefact rather than a correlation per
link · or a wheel is exhibited on which area and perimeter can be separated.

---

## W-26 — **THE BOUNDARY TOTAL IS A MIRROR, NOT A MEMORY. THE DISTRIBUTION IS NOT.** — and a relabelling that inverts three earlier readings

Directed by the principal, 2026-08-17: *"It's the boundary."*

### THE RELABELLING, AND IT CHANGES WHAT EARLIER RESULTS SAID

**On a graph, a surface is a SEPARATOR.** Measured on wheels `n = 4,5,6`: removing the **rim** leaves
`1` component and removing the **spokes** gives `2`. **The rim separates nothing. The spokes are the
edge boundary of the region `{hub}`.** The registrar called the rim "the boundary" because it is
drawn as a circle — a picture, not a property of the graph.

**Three consequences.** (i) W-25's "the chirality reaches the boundary" was measured on the **rim**,
which is not the boundary — that result is about the exterior's internal structure and is
**relabelled, not withdrawn**. (ii) The residual order surviving thermalisation (one trit at
`T = 50`) is **on the boundary**, not the interior. (iii) **W-24's requirement was half right**: a
surface must be a separator **and** a cycle. The fringe was neither; the wheel's rim is a cycle that
separates nothing; the wheel's spokes separate but are not a cycle. **No carrier in this program has
had both.**

### "FORCED BY GAUSS" WAS BEING MISREAD, INCLUDING BY ME

The hub's Gauss law fixes the **sum** of the spokes — which is *flux through the boundary equals
charge enclosed*. **That is the physical law, not bookkeeping.** The registrar has been voiding
constraint-determined results as vacuous; this one is the constraint **being** the physics, and
W-19's `I(S:F) = H(S)` is the same statement. **Max `S(spokes) = 4.0000` exactly at `n=5, Z_3`
against `5` for the rim: the trit is real and it is what Gauss says.**

**BUT IT MAKES THE BOUNDARY A MIRROR, AND THIS NEEDS NO MEASUREMENT.** Gauss holds at **every
instant**, so the boundary total can never disagree with the present interior charge and can never
retain one the interior has left. **Perfect fidelity, zero retention. A memory must be able to
disagree with the present; the constraint forbids it.**

### **THE DISTRIBUTION IS UNCONSTRAINED, AND IT RETAINS**

Gauss fixes **one** quantity — the sum of `n` boundary links. **`n-1` boundary degrees of freedom
are unconstrained by the interior.** Every boundary measurement in this program has been of the
total. Measured (`wheel n=4`, `Z_3`, dim 81): two physical states with the **same** boundary total
`0` and **different** distributions `(0,0,0,0)` vs `(0,0,1,2)`:

```
   t        0.0     0.5     1.0     2.0     5.0    20.0   100.0
   boundary TD    1.000   0.895   0.833   0.908   0.725   0.891   0.908
   rim TD (ctrl)  1.000   0.949   0.890   0.975   0.778   0.954   0.970
```

**The boundary distribution stays distinguishable while the total is pinned. The first quantity in
this program that is neither forced nor absent.**

**CAVEAT, LOAD-BEARING:** the rim control behaves the same way, so **nothing decays on this
carrier** — it may simply not thermalise at dimension 81. *"Stays distinguishable"* may be a
statement about a small closed system. **Not established as a boundary property until run on a
carrier that thermalises, or an open one.**

### NEXT

1. **A carrier whose boundary is BOTH a separator AND a cycle** — a 2-D lattice patch, where the
   edge boundary of a block of plaquettes is a closed curve *and* removing it disconnects.
2. **Rerun the distribution test where things decay** — larger, or open (coupled to a bath), so the
   control can fail.
3. **An open system is now required for a further reason** (W-25/AA): the global entropy is exactly
   `0` and stays `0` under unitary evolution. **Nothing can be paid in a closed pure system**, and
   every physical account of record formation involves paying.

**REOPENS IF:** the distribution retention is shown to be non-thermalisation · or a carrier is
exhibited whose boundary is both separator and cycle and behaves differently.

---

## W-27 — THE PLANAR PATCH: **A BOUNDARY THAT IS BOTH A CYCLE AND A SEPARATOR, AND WHERE EM'S TWO HALVES FINALLY MEET**

3x3 planar patch, `Z_2`, 12 links, 4 plaquettes, physical dim 16 (cycle rank 4).
**PERIMETER = 8 links: all degrees `2` (a CYCLE) and removal gives `5` components (a SEPARATOR).**
The registrar had claimed duality forbids one set being both. **It does not, and the patch has the
object** — the wheel's rim was a cycle that separated nothing, the fringe was neither.
**Discrete Stokes exact:** product of the 4 plaquettes = the perimeter loop, `0.00e+00`.

### WHERE EM WAS WHILE WE WERE MISSING THIS — **WE HAD HALF OF IT AT A TIME**

Electric flux is a **cut** quantity (Gauss); magnetic flux is a **cycle** quantity (Wilson). Every
prior carrier could express only one at its boundary: the fringe's dangling links lie on no cycle,
so **no magnetic structure at the boundary at all**; the wheel's rim is a cycle with no cut and its
spokes a cut with no cycle. **So the mirror result (Gauss, electric) and the chirality result
(Wilson, magnetic) were measured on different sets and could never interact.**

### AND THE CONJUGATE PAIR IS NOT WHAT THE REGISTRAR EXPECTED

Electric flux on the **whole perimeter is the IDENTITY** on the physical sector — `nontrivial: False`
— because a closed `Z`-loop is fixed by the Gauss laws it encloses. And `[E, M]` on the **same** set
is `0.000e+00` always: any subset whose `X`-product stays physical is a closed loop, and closed loops
here touch an **even** number of links, so anticommutations pair off.
**THE CONJUGATE PAIR IS A LOOP AND A CUT THAT PIERCES IT AN ODD NUMBER OF TIMES:** `Z` on an odd
subset of the perimeter against the perimeter loop gives `8.000e+00`; even subsets give `0.000e+00`.
**A boundary cannot hold a definite loop flux and a definite odd-piercing cut at once. First
structural limit in this program on what a boundary can record that is not a Gauss identity.**

---

## W-28 — **THE PROCESS, ASSEMBLED AND WORKING: A RECORD SURVIVES THE ENTROPY BEING PAID, IFF IT COMMUTES WITH THE DYNAMICS**

`LANE_W28_PROCESS/`, sealed. Carrier W-27's patch; open system via Lindblad, `gamma = 0.5`, `T = 10`.
**RECORD = the perimeter flux sector.** Write `rho_A` on `W = +1`, `rho_B` on `W = -1`, evolve, read.

```
  phase      bath        ||[H,W]||   S(global)      TD    <W>_A    <W>_B
  magnetic   magnetic     0.00e+00      3.0000  1.0000   1.0000  -1.0000
  magnetic   electric     0.00e+00      3.0000  1.0000   1.0000  -1.0000
  mixed      magnetic     3.20e+01      4.0000  0.0000   0.0000  -0.0000
  mixed      electric     3.20e+01      4.0000  0.0016   0.0008  -0.0008
```

**IN THE MAGNETIC PHASE THE RECORD SURVIVES PERFECTLY WHILE THREE BITS OF ENTROPY ARE PAID.** The
written values hold at `+1.0000` and `-1.0000`; the states stay perfectly distinguishable. **In the
mixed phase, where the record no longer commutes with `H`, it is destroyed completely.** Same
carrier, same bath, same coupling — **one variable: whether the recorded quantity is conserved.**

**THE POINTER CONDITION, MEASURED NOT IMPORTED:** *a record survives dissipation exactly when the
recorded quantity commutes with the dynamics.*

**WHAT DOES NOT HOLD, AND IT WEAKENS THE CLAIM.** In the magnetic phase **the bath choice does not
matter** — magnetic and electric baths both give `TD = 1.0000`. If einselection were doing the work
the electric bath should have damaged it. **So the conservation is doing everything and the
monitoring nothing.** This is a conservation result, not a decoherence result, and the einselection
half of the pointer story is **untested here**.

**AND WHY W-27b FAILED, NOW DIAGNOSED.** That run dephased in `Z` on the cut while `H`'s plaquette
terms rotate states out of the electric basis: Hamiltonian moves the record, bath measures it, and
the pair random-walks the distinction to zero (`TD 0.019` at `gamma = 1`). **Not a construction
failure — the recorded quantity was not conserved.**

**NEXT.** (1) Separate conservation from monitoring: find a bath that damages a conserved record, or
show none can. (2) The record here is a **global** flux sector, not a local one — test whether a
record localised to part of the boundary also survives. (3) The `gamma = 0` control in W-27b drifted
`S` to `0.0035` where it must be exactly `0`; the Heun integrator does not conserve purity and the
small-`gamma` numbers there carry unbounded error.

---

## W-29 — THE WRITE STEP: **WRITABLE AND DURABLE ARE CONJUGATE, AND THEY ARE CONJUGATE BECAUSE THE GEOMETRY IS FROZEN**

`LANE_W29_WRITE/`, sealed. Carrier: W-27's patch. **W-28 measured PERSISTENCE, not formation** — the
value was put in by hand and stayed. **Formation means the record ACQUIRES a value it did not have.**

**W-29, the null and its cause.** Writer = `Z` on a cut link at the interior vertex.
`||[Z_interior, W]|| = 0.000e+00` — **they commute, so nothing can be written.** Separation `0` at
every setting. **The registrar used an operator W-27 had already shown cannot conjugate the loop.**

**W-29b, the right operator, still null.** The conjugate of a loop is a cut that **pierces it an odd
number of times** — `||[Z_S, W]|| = 8.000e+00`. Best separation `0.058`, and only in the **closed**
case; with dissipation `~1e-5`.

> **THE OBSTRUCTION, STATED SHARPLY. W-28: the record survives dissipation IFF it commutes with `H`.
> W-29b: it can be written ONLY by something that does NOT commute with it. So the property that
> makes a record durable is exactly the property that makes it unwritable.** Preparing the writer
> definite leaves the record maximally uncertain, and conservation keeps it there forever.

**W-29c, the dynamics on its own.** *"If the dynamics of the parts force the record creation that's
a process as well"* (the principal). Start **maximally mixed**, nothing written, and let it run:

```
   bath          g2   gamma    T        <W>     |<W>|   S(final)
   plaquette    1.0    0.50  20.0   0.000000  0.000000    4.0000
   electric     1.0    0.50  20.0   0.000000  0.000000    4.0000
   proj +1      0.0    0.50  20.0   0.999909  0.999909    3.0007
   none         1.0    0.00  20.0   0.000000  0.000000    4.0000
```

**`|<W>| = 0.999909` — the dynamics drives the record to a definite sector from nothing.** BUT
**only for the `proj +1` bath, which is a projector onto the answer.** Every physically motivated
bath leaves `<W>` at exactly `0.000000` with entropy pinned at the maximum. **This is dissipative
state preparation: the selection is in the registrar's Lindblad operator, not in the physics.**
**What is missing is specific — a bath with TWO absorbing sectors and a broken tie**, so the system
selects rather than being told. That is spontaneous symmetry breaking, and it is untested.

### THE DESCRIPTION OF THE PROCESS AND ITS TERMS — `PROCESS_DESCRIPTION_V001.md`

**EM: PRESENT** — both Hamiltonian terms **and** the constraint; electric, magnetic and Gauss are
Maxwell's equations. **ALPHA: PRESENT, AND THERE ARE TWO SLOTS** — `g^2`, and `theta` which is
dimensionless, periodic and does nothing at `Z_2`. **GRAVITY: ABSENT, AND THE TEMPTING FILL IS
WRONG** — assigning it to the Gauss constraint counts EM twice, because Gauss is Maxwell's first
equation.

> **WHAT IS MISSING IS THAT THE CARRIER IS STATIC.** The graph is fixed and nothing the field does
> changes it. Every gravitational functional signature reduces to this. **A gravitational term would
> be the carrier responding to the field.**

**AND THAT IS WHY THE WRITE FAILED.** With a frozen geometry the only conserved quantities are the
ones the geometry hands you, and those are precisely the ones nothing can write to.

**NEXT.** (1) A bath with **two** absorbing sectors — does the system select, and does the selection
correlate with anything in the initial state? (2) **A carrier that responds to the field** — the
missing term, and the first construction in this program with a dynamical geometry.

---

## W-30 — **THE W-29 OBSTRUCTION IS NOT A FEATURE OF OUR PATCH. IT IS GENERAL, AND NO VALUE OF THE COUPLING ESCAPES IT.**

`LANE_W30_DYNAMIC_CARRIER/`, sealed. W-29 found the conjugacy on one carrier. **A finding on one
carrier is a property of that carrier until tested elsewhere.**

### W-30a — FOUR CARRIERS, TWO GROUPS, ONE ANSWER

Carriers: the 3×3 patch, a single 2×2 plaquette, a triangle (3-cycle), and the theta graph (two
vertices, three parallel links, **two** independent cycles). Groups `Z_2` and `Z_3`. `R` = rim loop
built as the **product of plaquettes**, so orientation is correct by construction. At `g^2 = 0`:

```
  carrier         N  dim   ||R||  unitarity  ||[H,R]||  ||[L,R]||  max|d<R>/dt|
  3x3 patch       2   16   4.000  0.000e+00  0.000e+00  0.000e+00     1.110e-16
  3x3 patch       3   81   9.000  0.000e+00  0.000e+00  0.000e+00     2.252e-17
  2x2 plaquette   2    2   1.414  0.000e+00  0.000e+00  0.000e+00     0.000e+00
  triangle        3    3   1.732  0.000e+00  0.000e+00  0.000e+00     1.144e-16
  theta graph     3    9   3.000  0.000e+00  0.000e+00  0.000e+00     1.241e-16      (8 rows, all)
```

`R` is **unitary** (defect `0.000e+00`) and **nontrivial** — exactly `N` distinct eigenvalues on
every carrier, so it can label `N` sectors. It is a perfectly good record-bearing observable.
**`d<R>/dt` is zero to machine precision on every carrier, at every group, for every random state.**

> **THE NO-GO.** A perfectly conserved quantity cannot be written, **because that is what conserved
> means.** Its value does not change, so it was never acquired — it always was. And a quantity that
> is not conserved is one the bath can reach. **On a static carrier there is no third option.**

### W-30b — THE COUPLING IS THE DIAL BETWEEN THE HORNS, AND IT DOES NOT REACH A WAY OUT

Exact Lindblad exponential (RK4 agrees to `~1e-14`, CONTROL D — this also **retires the `w27b_open`
Heun defect** for this lane). `LASTS` = time-average of `<R>` over `t ∈ [10,20]`; `OSC` = the same at
`gamma = 0`, separating precession from survival.

```
     g^2   MOVES=||[H,R]||/dim   LASTS(gam=.5)   OSC(gam=0)   LASTS-OSC
    0.00                0.0000        1.000000     1.000000    0.000000
    0.01                0.0400        0.952075     0.998332   -0.046257
    0.03                0.1200        0.646844     0.985428   -0.338584
    0.10                0.4000        0.016774     0.856770   -0.839995
    1.00                4.0000        0.000076     0.211098   -0.211023
   20.00               80.0000        0.240293     0.319244   -0.078951
```

**At `g^2 = 0` the record is untouchable: `<R>` = `1.000000000000` with the bath ON** (CONTROL A).
**The bath cannot damage a conserved record at all.** And `MOVES` is exactly `0` — nothing can write
it. **Switch the coupling on and both change together.** By `g^2 = 0.1` survival is `1.7%`.

**The apparent revival at large `g^2` is not survival.** There `LASTS ≈ OSC` — the same value with
the bath off — so it is the **dephasing residue** (the projection onto the commutant of a
now-dominant electric term), and it is flat at `~0.3`, not rising. The `MOVES×LASTS` product peaks
at the **edge** of the scan only because `MOVES = 4g^2` grows without bound; that is the metric
inflating, not a record persisting.

> **NO VALUE OF THE COUPLING SOLVES IT. THE OBSTRUCTION IS TOTAL ON A STATIC CARRIER.**

### WHAT THIS ESTABLISHES — GRAVITY'S ROLE, DERIVED FROM ITS ABSENCE

The charter asks for the roles of EM, gravity and alpha in record formation. **EM** supplies the
field, the boundary that is both cycle and separator, and the loop that is a candidate record.
**Alpha (`g^2`)** is the dial along the axis of the obstruction — it is what trades writability
against durability, and it lands on that axis without being put there. **And the missing ingredient
is now specified by what it must do:**

> **A record can only form if a quantity can BECOME conserved. On a fixed `H` the conserved
> quantities are fixed for all time, so nothing can become one. The carrier must respond.**

**This is not "gravity would be nice to have." It is: without a responsive carrier there are no
records at all.** Record formation REQUIRES the ingredient we have been calling gravity, and the
requirement is derived from a null, not assumed.

**CAVEAT, STANDING.** The escape not yet tested is **approximate** conservation — einselection uses
pointer states stable on observation timescales, not forever. That is a hierarchy of timescales, and
a timescale ratio is a pure number. **Whether alpha sets that ratio is untested and is the obvious
next question if the dynamical carrier fails.**

---

## W-32 — **THE WINDOW EXISTS. W-30's "NO COUPLING ESCAPES" IS WRONG, AND THE COUPLING SETS THE WIDTH.**

`LANE_W32_TIMESCALES/`, sealed. Opened at the principal's direction: real records are never exactly
conserved — domains flip, crystals anneal, DNA mutates — so the question is not whether `<R>`
survives to a fixed time but **how many times the record can be written before it decays.** That
ratio is dimensionless.

### ERRATUM TO W-30

**W-30b measured survival at a FIXED ABSOLUTE TIME (`T=10`) while the WRITE timescale was itself
changing with `g^2`.** Two moving quantities, one held still. Corrected here by taking both rates
from the Lindbladian spectrum, where no observation window enters.

> **W-30's sentence "no value of the coupling solves it; the obstruction is total on a static
> carrier" IS WITHDRAWN.** The no-go on *exact* conservation (W-30a) stands untouched and is
> unaffected. What fails is the extension of it to approximate records.

### THE MEASUREMENT

`Gamma_R` = slowest decay rate among left-eigenmodes of the Lindbladian **whose overlap with `R`
exceeds `1e-6`** — so it is the record's lifetime, not an unrelated mode's. Overlap is `1.000` at
weak coupling; **at `g^2=0` NO mode carries `R` at all** (`1.3e-16`), confirming exact conservation
there and exposing that W-32's first pass had reported an unrelated mode's rate.

```
  BATH ON THE CUT -- disjoint from the rim, [L,R] = 0
       g^2        omega        Gamma_R        RATIO   d ln(RATIO)/d ln(g2)
   0.00050     0.002000   7.999981e-06     250.0006
   0.00100     0.004000   3.199970e-05     125.0012                 -1.000
   0.00200     0.008000   1.279952e-04      62.5024                 -1.000
   0.00500     0.020000   7.998119e-04      25.0059                 -1.000
   0.01000     0.040000   3.196988e-03      12.5118                 -0.999

  BATH ON THE RIM -- the environment sees the record directly, [L,R] = 8
   0.00050     0.002000   5.999266e+00       0.0003
   0.00500     0.020000   4.992902e+00       0.0040                 +1.001
   0.05000     0.200000   1.022078e+00       0.1957                 +2.166
```

**`omega = 4g^2` (first order). `Gamma_R = 32g^4` (second order). `RATIO = 1/(8g^2)`, slope `-1.000`
across two decades.** The record can be written `250` times over at `g^2 = 0.0005`, and the count
**doubles every time the coupling halves.**

### THE EXPONENT IS PHYSICS, NOT GENERIC COUNTING

The obvious deflation is that `omega ~ g^2` and `Gamma ~ g^4` is just golden-rule order counting.
**It is not, and the control settles it:**

```
  bath on CUT :  d ln(Gamma)/d ln(g2) = +1.996     RATIO ~ 1/g^2   protection DIVERGES as g2 -> 0
  bath on RIM :  d ln(Gamma)/d ln(g2) = -0.310     RATIO ~ g^2     protection VANISHES as g2 -> 0
```

**Same Hamiltonian, same record, same coupling strength — opposite exponents.** The difference is
only **where the environment touches.** When the bath commutes with the record it can reach it only
*through* the `H`-mixing, so decay is one order higher than drive and the record is shielded. When
the bath sees the record directly, `Gamma_R` is flat at `~5` regardless of coupling and there is no
protection at any `g^2`.

> **THE CONDITION FOR A RECORD IS GEOMETRIC. The record lives on the rim; the environment couples to
> the cut. `[L,R] = 0` is the shielding condition, and the W-27 patch — whose perimeter is
> simultaneously a CYCLE and a SEPARATOR — is exactly the structure that supplies it.** That
> construction was built at W-27 and only now does work.

### WHAT THIS DOES TO THE THREE TERMS

**ALPHA GAINS A SECOND, SHARPER ROLE.** Beyond being the dial along the obstruction axis (W-30b),
`1/g^2` **is the number of times a record can be written and re-read before it decays.** A pure
number, measured, not fitted.

**AND THE W-30 CLAIM ABOUT GRAVITY IS WEAKENED AND MUST BE RESTATED.** W-30 concluded "without a
responsive carrier there are no records at all." **That is now false as stated** — approximate
records form and persist on a static carrier, and approximate is what real records are.

**WHAT REMAINS OPEN IS SELECTION, NOT DURABILITY.** W-32 shows a record, once written, persists
`1/(8g^2)` write-times. **W-29c showed nothing writes it spontaneously** — every physically
motivated bath left `<W>` at exactly `0.000000` with entropy pinned at maximum, and only a projector
that already knew the answer broke the tie. **So the live question is no longer "can a record last."
It is "what decides WHICH value it takes."** That is the tie-breaking / symmetry-breaking problem,
and it is where the responsive carrier must now be aimed.

---

## W-33 — **THE CARRIER IS A MEMORY, NOT A RECORDER. IT PRESERVES EXACTLY; IT NEVER AMPLIFIES.**

`LANE_W33_SELECTION/`, sealed.

**STEADY-STATE COUNT — decisive before any dynamics is run, because a Lindbladian with a UNIQUE
steady state cannot select anything from any initial condition.**

```
   g^2      bath   #zero modes
  0.0000     cut         2      DEGENERATE -> a record manifold exists
  0.0010     cut         1      unique                 (next rate 3.200e-05)
  0.0100     cut         1      unique                 (next rate 3.197e-03)
  any        rim         1      unique                 (next rate ~1.0)
  any    all links       1      unique
```

**The record manifold is EXACT only at `g^2 = 0` with a shielded bath.** For `g^2 > 0` it is lifted
— but at rate `32g^4`, so for small coupling the manifold is approximate and long-lived, which is
W-32's window seen from the spectral side.

**WHAT THE STEADY STATE REMEMBERS.** Start from `rho(b) = (I + b R)/dim` and take `t -> infinity`
exactly by projecting onto the zero modes:

```
  bath=cut g2=0.0000   slope d<R>(inf)/d<R>(0) = +1.000000000   PRESERVES EXACTLY
      -0.90->-0.900  -0.50->-0.500  -0.20->-0.200  0.00->0.000  +0.20->+0.200  +0.50->+0.500  +0.90->+0.900
  bath=cut g2=0.0010   slope = -0.000000000                     erases (at t=inf; 1/(8g^2) cycles first)
  bath=rim any g2      slope = +0.000000000                     erases
```

> **THE SLOPE IS EXACTLY `1` OR EXACTLY `0`. IT IS NEVER GREATER THAN `1`, AT ANY SETTING TESTED.
> THE DYNAMICS NEVER BREAKS A TIE. IT KEEPS WHAT IT IS GIVEN AND MANUFACTURES NOTHING.**

**AND THIS REFRAMES W-29c'S NULL.** W-29c started from a maximally mixed state — which carries no
information — and found `<R> = 0.000000`. **That is the correct answer, not a failure: a record of
nothing is nothing.** Reproduced here exactly (`+0.000000000000` at every setting). The construction
was never failing to hold a record; it was being asked to invent one.

## W-34 — **THE PREDICTABILITY SIEVE SELECTS THE BOUNDARY, AND WE NOMINATED NOTHING**

Every prior lane NOMINATED the rim loop and then asked how it fared — which installs the answer.
Here nothing is nominated. The left-eigenmodes of the Liouvillian **are** the observables with
definite decay rates, so the slowest ones **are** the pointer observables. Read them off and ask
what they are.

```
  g2=0.001, bath on CUT          rate        best operator match
        (identity)          2.674e-15        identity
        SLOWEST OBSERVABLE  3.200e-05        RIM LOOP (all 4)      overlap 1.000
        next                5.803e-01        no clean match
        next                9.841e-01        no clean match
```

**CROSS-CHECK — all 16 magnetic operators ranked by their own decay rate, `g^2=0.01`:**

```
    7.157165e-02   RIM LOOP (all 4)          <-- 28x slower than anything else
    1.983944e+00   plaquette (0,)
    1.985041e+00   plaquette (1,)
    2.010041e+00   3-plaquette loop (0,1,3)
    2.029031e+00   2-plaquette loop (0,2)
    3.860385e+00   2-plaquette loop (0,3)
```

> **THE DYNAMICS PICKS THE FULL BOUNDARY, AND ONLY THE FULL BOUNDARY.** Not "bigger is slower" —
> single, double and triple loops are all clustered at `~2.0`. The rim wins by a factor of 28.

**AND THE MECHANISM IS GEOMETRIC, NOT ENERGETIC. The rim loop is the unique gauge-invariant loop
whose support (`PERIM`) is DISJOINT from where the environment couples (`CUT`).** Every other loop
runs through a cut link. **The pointer observable is the one the environment cannot touch.**

**CONFIRMED BY MOVING THE BATH.** With the bath on the RIM instead, no magnetic operator survives at
all — the slow modes become **electric-type (diagonal)**, all at rate `~1.0`. **The pointer basis
follows the environment's coupling.** That is einselection behaving exactly as advertised, and it
means the boundary is not privileged a priori: it is privileged *given where the environment is*.

**THIS RETIRES THE CIRCULARITY.** The rim loop was chosen by the registrar at W-27 and used as the
record ever since. It is now independently re-derived as what this construction actually protects.

---

## W-35 — **THE RECORD FORMS ON EVERY RUN. THE INSTRUMENT WAS BLIND TO IT FROM W-28 TO W-34.**

`LANE_W35_TRAJECTORIES/`, sealed. Raised by a design agent in the W-31 panel and tested here.

**Every lane from W-28 onward measured `<R> = Tr(R rho)` — an ENSEMBLE AVERAGE. A fair coin has mean
zero and every single flip is definite.** So `d<R>/dt = 0` is consistent with every individual run
acquiring a perfectly definite record. Unravel the same Lindbladian and measure the second moment.

```
  HOMODYNE, 600 trajectories, T=15, every start has <R> = 0 EXACTLY
     g^2  gamma      E[<R>]     +-SE    E[<R>^2]   |<R>|>0.9  |<R>|>0.99
   0.000   0.00    0.000000   0.0000    0.000000       0.000       0.000   <- GATE 0
   0.000   0.10    0.003189   0.0324    0.631305       0.432       0.113
   0.000   0.30   -0.014558   0.0397    0.945506       0.922       0.785
   0.000   1.00   -0.019572   0.0408    0.999078       0.998       0.995
   0.100   0.30   -0.001021   0.0258    0.398010       0.082       0.002
   1.000   0.30    0.006214   0.0122    0.090022       0.003       0.000

  GROWTH IN TIME (g2=0, gamma=0.3):  t=0: 0.000  t=1: 0.196  t=3: 0.445
                                     t=6: 0.700  t=10: 0.867  t=15: 0.946  t=25: 0.991
```

> **`E[<R>] = -0.0146 +- 0.0397` (zero, as every previous lane found) WHILE `E[<R>^2] = 0.9455` AND
> 78.5% OF RUNS SIT AT `|<R>| > 0.99`. At `gamma=1.0`, `E[<R>^2] = 0.999078` and 99.5% of runs are
> definite. THE TIE IS BROKEN ON EVERY SINGLE TRAJECTORY. THE VALUE IS A FAIR COIN.**

**GATES.** `gamma=0` gives `E[<R>^2] = 2.7e-29` (machine zero — the integrator invents nothing);
stable across `dt` from `0.008` to `0.001`; and the trajectory ensemble reproduces the exact master
equation (`<R>_exact = -0.00000000`, `<R>_traj = -0.0146`, within one standard error).

### THE SAME LINDBLADIAN, UNRAVELLED THE OTHER WAY, HAS NO RECORD AT ALL

```
  JUMP unravelling, identical generator:  g2=0, gamma=0.3 -> E[<R>^2] = 0.000000
                                          g2=0, gamma=1.0 -> E[<R>^2] = 0.000000
```

**Both unravellings give the same `rho(t)`.** With unitary `L`, jump rates are state-independent
(`gamma||L psi||^2 = gamma` always), so a jump carries **no information** and every trajectory keeps
`<R>` fixed. Homodyne extracts information continuously. **Record formation is therefore not a
property of the master equation. It is a fact about what the environment does with what it learns.**

### ERRATUM — W-29c AND W-33

**W-29c reported that no physically motivated bath selects a sector (`<W> = 0.000000` exactly) and
W-33 reported that the dynamics "preserves but never amplifies" (slope exactly `1` or `0`, never
above `1`). Both numbers are correct and both readings are WITHDRAWN.** They measured the *fairness*
of the coin. The selection they were looking for was happening on every run.

**W-32 and W-34 are NOT affected** — decay rates and pointer selection are properties of the
generator and stand as measured.

### AND IT RESOLVES THE MECHANISM

`[L,R] = 0` was read at W-32 as *shielding* — the environment cannot reach the record. **It is more
than that. `[L,R] = 0` means the environment cannot DISTURB `R`, while `L` and `R` remain CORRELATED
in the state, so the environment LEARNS `R` anyway. That is a quantum non-demolition measurement,
and QND is exactly the condition for a record.** The coupling `g^2` spoils it by making `R`
non-conserved: `E[<R>^2]` falls `0.9455 -> 0.398 -> 0.090` as `g^2` goes `0 -> 0.1 -> 1.0`.

**THE PROCESS, ASSEMBLED FROM WHAT IS MEASURED.** `[H,R] = 0` makes the record durable (W-30a).
`[L,R] = 0` makes it undisturbable while still legible (W-35). Correlation between `L` and `R` makes
the environment learn it. Continuous monitoring collapses each run to a definite value (W-35).
`1/(8g^2)` sets how long it lasts (W-32). And the dynamics selects which observable this happens to
— the full boundary, by a factor of 28 (W-34).

---

## W-36 — **THE RECORD IS NOT OBJECTIVE, AND THE OBSTRUCTION IS LOCALITY**

`LANE_W36_REDUNDANCY/`, sealed. W-35 left record formation depending on how the environment is
read, which is not good enough. Zurek's criterion settles it without choosing an unravelling: a
record is **objective** when many disjoint fragments EACH carry a full copy. That is a property of
the joint state, so the environment is built explicitly and never traced out.

```
  I(R:F) in bits, 6 explicit environment qubits, kappa=8.0, T=16.0, ceiling 1 bit
  |F|                                    0      1      2      3      4      5      6
  A  env on Z at the CUT   (LOCAL)    0.000  0.044  0.100  0.187  0.351  0.658  0.929
  B  env on R itself    (NONLOCAL)    0.000  0.999  1.000  1.000  1.000  1.000  1.000
  C  env on one plaquette  (LOCAL)    0.000  0.037  0.037  0.037  0.037  0.037  0.037
  CONTROL kappa = 0                  -0.000 -0.000 -0.000 -0.000 -0.000 -0.000 -0.000
```

**A — the configuration every earlier lane used — has NO redundancy.** `I(|F|=1)/I(all) = 0.047`:
information rises only as the fragment approaches the whole environment. **One delocalised copy. An
observer holding any part of the environment knows essentially nothing.**

**B is the textbook quantum-Darwinism plateau — `0.999` bits from a SINGLE qubit, flat thereafter.
Every fragment independently holds the whole record. That is objectivity.**

**C reads `ratio 0.999` and is NOT redundancy — the total is `0.037` bits, so every fragment knows
equally nothing.** The ratio test is meaningless without a floor on the total, and it is gated in
`PUBLISHED_CONVENTIONS`.

### THE OBSTRUCTION, STATED

> **`R` IS A PRODUCT OVER 8 RIM LINKS. IT IS NONLOCAL. Redundancy requires the environment to couple
> to the record observable itself, and a LOCAL environment cannot couple to a nonlocal one.**

**And this is the same shape as W-30's obstruction, one level up.** W-30: writable vs durable.
**W-36: legible vs protected.** `[L,R] = 0` is what protects the record (W-32) and it is also what
stops any local fragment from carrying a copy. **The observable that survives is the one nobody can
read.**

**Note that B is not smuggling.** `L = R` is symmetric between the `+1` and `-1` sectors, so it
cannot decide which wins — W-35's coin stays fair — and `R` was selected independently by the W-34
sieve, not nominated here.

### THE NEXT STEP, AND IT IS PHYSICAL

**In gauge theory a nonlocal observable is measured by transporting a probe charge around the loop —
Aharonov–Bohm.** The interaction is local at every instant; the nonlocality is in the PATH, paid for
in TIME proportional to the loop's length. **So the environment CAN read the rim loop, by a local
process extended in time, and case B is the effective description of exactly that.**

**BUILD: a probe charge transported around the rim, and measure whether the record it acquires is
redundantly copied.** If it is, the process closes: the record is written by monitoring, protected
by `[L,R]=0`, selected by the sieve, persists `1/(8g^2)` write-times, and is made objective by
transport. **If the transport time itself destroys the record, that is a new obstruction and the
timescale competition becomes the next quantity to measure.**

---

## W-37 — **A LOCAL PROBE READS THE NONLOCAL RECORD, AND THE READING IS PURELY TOPOLOGICAL**

`LANE_W37_TRANSPORT/`, sealed. W-36's obstruction was that `R` is nonlocal and a local environment
cannot copy it. Gauge theory's own answer is transport.

**EXPLICIT MATTER.** A probe charge hopping the 8 rim vertices, a static anti-charge at the centre.
**The constraint bites and was checked rather than assumed: with the probe alone the total `Z_2`
charge is odd and the physical sector contains `0` states.** Gauss becomes matter-dependent —
`div(s)_v = [probe site] + [centre]` — so **the physical sector MOVES with the probe.** Physical
dim `128 = 8 × 16`. The only interaction is the hop; it touches **one link at a time** and never
mentions `R`.

```
      T   I(R:probe) bits    P(probe back at 0)
   0.00        0.000000                1.0000
   1.00        0.002606                0.0501
   2.00        0.270865                0.1578
   3.00        0.669357                0.0355
  12.00        0.829907                0.4985

  CONTROL A  tau = 0 (no transport)          I = 2.403e-16
  CONTROL B  RING CUT (one hop removed)      I = 0.000000 at T = 4, 8, 12, 20
```

> **`I(R:probe)` reaches `0.83` bits through hops that each touch a single link. CUT THE RING —
> keeping every local interaction and deleting only the closed path — AND IT IS EXACTLY `0.000000`
> AT EVERY TIME. The information is topological. It exists only because a closed path exists.**

**`||[R, H_hop]|| = 0.0e+00`.** The transport **cannot disturb the record.** So the reading is QND,
and **many probes may each read it independently** — which is redundancy, and therefore objectivity.
**W-36's condition 5 is met, by transport.**

`I` oscillates rather than saturating because the dynamics here is coherent and nothing decoheres the
probe; the information sloshes between probe and field. **Locking it in requires decohering the
probe's position, which is untested.**

### THE CLOSING RELATION, AND WHAT IS STILL UNMEASURED IN IT

Readout takes a traversal: `T_read ~ perimeter / tau`. The record decays at `Gamma = 32 g^4` (W-32).
**A region can be read before it forgets only if `32 g^4 · perimeter / tau << 1`** — a relation
between the coupling, the size of the region, and the transport rate.

**THE PERIMETER SCALING IS NOT MEASURED.** One patch, one ring size, so `T_read ∝ perimeter` is an
expectation and not a result. **Rings of several sizes are the next build, and if the scaling holds
the relation above is the program's next pure number.**

## W-37b/c — **THE PERIMETER LAW HOLDS, AND READABILITY DEPENDS ON THE PARITY OF THE LOOP**

```
  bare cycle graphs, tau = 1.0
     n   dim  ||[R,H]||   T_read (I=0.5 bit)   max I     T_read/n
     4     8    0.0e+00          1.25          0.9998      0.312
     5    10    0.0e+00         never          0.0000        --
     6    12    0.0e+00          1.75          0.9827      0.292
     7    14    0.0e+00         never          0.0000        --
     8    16    0.0e+00          2.50          0.9901      0.312
    10    20    0.0e+00          3.00          0.9329      0.300
    12    24    0.0e+00          3.50          0.9338      0.292

  d ln(T_read) / d ln(n) = +0.960        CUT RING, every n: max I <= 4.4e-16, never reaches 0.5
```

> **`T_read/n` IS CONSTANT AT `0.30` ACROSS A FACTOR OF THREE IN RING SIZE. THE PERIMETER LAW IS
> MEASURED, NOT ASSUMED** (W-37 had flagged it as an expectation). And the cut-ring control returns
> `<= 4.4e-16` at every size — **no closed path, no information, at any scale.**

**SO THE CLOSING RELATION IS NOW A MEASURED ONE.** `T_read = 0.30 · perimeter / tau` and
`Gamma = 32 g^4`, so a region can be read before it forgets iff

> **`9.6 · g^4 · perimeter / tau  <<  1`** — a pure number joining the coupling, the SIZE of the
> region, and the transport rate. **A region has a maximum legible size set by the coupling.**

### AND AN UNEXPECTED ONE — ODD LOOPS CARRY NOTHING

`n = 5` and `n = 7` give **exactly** `0.0000`, and TEST 1 confirms it across **every start site and
every time** (`4.4e-16`, `2.2e-16`). An exact zero is a symmetry.

**TEST 3 finds it.** For a ring, `R=+1` is periodic hopping and `R=-1` is antiperiodic. At `n=5` the
periodic spectrum is `[-2, -0.618, -0.618, 1.618, 1.618]` and the antiperiodic is
`[-1.618, -1.618, 0.618, 0.618, 2]` — **the exact negative.** For odd `n` the antiperiodic momenta
ARE the periodic momenta shifted by `pi`, and `cos(k+pi) = -cos(k)`.

> **ON AN ODD LOOP THE FLUX IS ABSORBED INTO A SHIFT OF THE PROBE'S MOMENTUM. THIS PROBE CANNOT SEE
> IT — not because the record is absent, but because the carrier's ARITHMETIC hides it.**

**Stated carefully: this is a fact about a single particle hopping to nearest neighbours. Whether
some other probe can read an odd loop is untested,** and it is a sharp, cheap question.

**A DEFECT, LOGGED.** `w37c` TEST 2 projects a localised start into each sector and that projection
is degenerate at every `n`, so it printed "(a sector is empty)" throughout and measured nothing. It
is kept as written. The result rests on TEST 1 and TEST 3.

---

## ERRATUM (VECTORISATION) — **A REAL BUG IN W-32/33/34, AND IT MOVES NOTHING**

Found by a W-31 adversary. `numpy`'s `reshape(-1)` is **row-major**, so `vec(AXB) = (A ⊗ Bᵀ) vec X`
and the Lindblad generator is `-i(H⊗I - I⊗Hᵀ) + γ Σ (L⊗L* - I⊗I)`. Those lanes used the
**column-major** form `-i(I⊗H - Hᵀ⊗I)` with `kron(L*,L)`.

**Verified against RK4 with a deliberately non-symmetric complex state:** column-major is wrong by
`2.8e-02`; row-major matches to `5.1e-15`.

**All four scripts were patched and re-run. Every registered number is unchanged.** The diffs are
machine-precision only (`1.32e-16 → 2.41e-16`), signs of exact zeros, and reordering among
degenerate eigenvalues. **The reason is structural: the two conventions differ by the swap
permutation, so they are SIMILAR and every eigenvalue is identical** — and W-32's rates, W-34's
sieve ranking and W-33's steady-state counts are all spectral. **W-32's `1/(8g²)`, W-34's 28× margin
and W-33's slope `1.000000000` stand as registered.** W-35's main result used trajectory evolution
and never touched the generator; W-36 and W-37 evolve state vectors.

## W-31 — **ROUTE B: THE RESPONSIVE CARRIER, REFUTED 3/3. THE TOPOLOGY DID NO WORK.**

`LANE_W31_RESPONSIVE_CARRIER/`, sealed. Commissioned as 4 independent designs under distinct lenses,
3 judges, 1 build, 3 adversaries. **The registrar's own hypothesis was withheld from the designers.
The gr-faithful lens found it independently and sharpened it** — from "H depends on the state" to a
topology change: `||[R_E, H_capped]|| = 16/32/64 > 0` (writable while the cap exists) and
`||[R_E, H_punctured]|| = 0.000e+00` (exactly conserved once the hole opens). All three judges
picked it unanimously (`breaks_nogo = 9,9,9`).

**ALL THREE ADVERSARIES REFUTED IT, INDEPENDENTLY, AFTER REPRODUCING THE ARITHMETIC EXACTLY.**

> **ADVERSARY 2 BUILT A TOPOLOGY-FREE TWIN** — no hole, no cap, `b₁ = 0` everywhere, no carrier
> response: two copies of the same capped disk, the second with the writer set to `κ=0` by hand,
> joined by a leak. **It reproduces EVERY headline number to `1e-16`**, spread `0.107789` included.
> Because `ιᵀ H_disk(κ=0) ι = H_annulus + c·I` to `8.9e-16`: **the "punctured annulus" IS the
> zero-spoke-flux sector of the capped disk with the writer switched off.**

**ADVERSARY 3** independently built a 32-dimensional control — no hole, no second Gauss law, no Betti
number — and got the same numbers. **ADVERSARY 1** showed the headline metric is decoupled from what
it claims to measure: the `g²=0` null is forced by a weak symmetry `S` that anticommutes with `R_E`
and leaves all five "unbiased" starts exactly invariant; a record **does** form at `g²=0`
(`±0.145299`), contradicting the build's own reading; and the frozen value decomposes as
`+0.034219` population `+ 0.081116` **coherence** — a continuous readout of the initial relative
phase, so **nothing discrete was recorded at all.**

### WHAT THIS SETTLES

> **THE RESPONSIVE CARRIER CONTRIBUTED EXACTLY NOTHING. The load-bearing content was "the
> destination has no writer, and the transfer commutes with the record" — both hand-set.**

**And it converges with V003 from the opposite direction.** V003 concluded that nothing measured yet
REQUIRES the carrier to respond. Route B was the attempt to build a case that does, and the strongest
design three judges could pick reduces to a hand-set constant. **Two independent routes, same
answer: the gravity line is still empty, and the emptiness is now a result rather than a gap in
effort.**

## W-37d — **THE ODD-LOOP NULL IS ABOUT THE PROBE, NOT THE RECORD. OBJECTIVITY IS PROBE-RELATIVE.**

W-37c registered odd loops as carrying "exactly nothing" and the registrar billed it as a fact about
the carrier. **Tested; it is a fact about the probe.**

**MECHANISM (P1, confirmed).** A `Z_2` flux of `pi` and a sign flip of every hopping are the SAME
operation when `n` is odd, because negating all `n` amplitudes shifts the total flux by `n·pi`.
`spectrum(-H(flux 0)) == spectrum(H(flux pi))` at `n = 5,7,9,11` and NOT at `n = 4,6,8,10`. With
real hopping the position distribution is time-symmetric, so the two sectors become identical **by
position** — which is all the probe reported.

**DECISIVE (P3).** Same odd ring, same record, probe given a complex hopping phase:

```
   n      phase    max I(R:probe) bits
   5     0.0000               0.000000
   5     0.3000               0.922109
   5     1.5708               0.973503
   7     0.0000               0.000000
   7     1.5708               0.874649
```

> **NOTHING ABOUT THE RECORD CHANGED. THE NULL WAS THE PROBE'S.**

**CONSEQUENCE FOR THE PROCESS DESCRIPTION.** W-37 claimed condition 5 (redundant copying →
objectivity) was met. **It is met RELATIVE TO A PROBE.** A record can be present, protected, and
invisible to one probe while fully legible to another. **`V003` §1 condition 5 must be read with
that qualifier**, and "objective" in this program means "objective to probes that can couple to it."

**A FAILED PREDICTION, LOGGED AND NOT BUILT ON.** P2 predicted `Z_3` odd rings would read fine.
`Z_3` returned **exactly `0.000000` at every ring size including EVEN ones**, where `Z_2` gives
`0.999792`. That contradicts the `Z_2` result rather than extending it, so it indicates a defect in
the `Z_3` construction — most likely the partial trace, which matches gauge configurations across
probe positions that satisfy DIFFERENT Gauss laws. **Unresolved. Nothing is inferred from it.**

**ALSO CORRECTED.** `w37c`'s docstring asserts the odd-`n` spectra "coincide as SETS"; its own TEST 3
prints `different` at every `n`, and the true relation is that they are exact NEGATIVES. The register
entry for W-37c stated the negative relation correctly; the code comment is wrong and is superseded here.

---

## W-38 — **THE REGION IS NOT GIVEN. THE ENVIRONMENT DECIDES WHAT COUNTS AS A REGION.**

`LANE_W38_REGION/`, sealed. Opened on the principal's reading: **a region is another way of thinking
about a boundary.** On a graph that is an identity — a region IS a cycle that separates — so "which
region carries the record" and "which boundary is selected" are one question, and W-34 already
answered a special case of it while the registrar noted the general one in passing and moved on.

**15 regions** (subsets of the 4 plaquettes), each with its boundary loop. Sieve run seven times,
with the environment in a different place each time.

```
              bath       winner        rate    runner-up   margin   winner's boundary
                                                                    touching the bath
      CUT (centre) (0, 1, 2, 3)  7.1572e-02         (1,)    27.7x        0 of 8
 one corner plaq 0         (3,)  6.0883e-02         (2,)    16.4x        0 of 4
 one corner plaq 3         (0,)  6.0883e-02         (2,)    16.4x        0 of 4
     single link 0         (3,)  2.7401e-03         (1,)    14.1x        0 of 4
       left column         (3,)  4.5939e-02         (1,)     1.0x        0 of 4   (tie: two regions avoid it)
   RIM (perimeter)         (3,)  1.9812e+00         (0,)     1.0x        2 of 4   (no selection)
         all links         (0,)  3.9455e+00         (3,)     1.0x        4 of 4   (CONTROL: spread only 2.03x)
```

> **THE SURVIVING REGION IS ALWAYS THE ONE WHOSE BOUNDARY IS DISJOINT FROM THE ENVIRONMENT, AND IT
> MOVES WHEN THE ENVIRONMENT MOVES.** Bath at the centre → the WHOLE patch survives. Bath on corner
> plaquette `0` → the OPPOSITE corner `(3)` survives. Bath on corner `3` → corner `(0)` survives.
> Bath on a single link → a region avoiding it. **Same carrier, same Hamiltonian, same field. Only
> the environment's location changed.**

**THE CONTROLS FIRE CORRECTLY AND THEY ARE WHAT MAKE THIS READABLE.** With the bath on **all** links
no region is spared and the spread collapses to `2.03x` — no winner. With the bath on the **rim**,
which touches nearly every boundary, the margin is `1.0x` — no selection. And **"left column"
produces an honest TIE**, because regions `(1)` and `(3)` both have boundaries disjoint from
`[0,4,6,9]`. **When several regions avoid the environment, several regions tie.** The effect is
not "some region always wins"; it is specifically avoidance of the environment.

### WHAT THIS ESTABLISHES

**This is the first result in the program where a GEOMETRIC notion is an OUTPUT rather than an
input.** Every region in this program until now was chosen by the registrar — the 3×3 patch, its
perimeter, the rim loop. Here the carrier offers 15 regions and says nothing about which is real;
**the environment's placement picks one, with a margin of 14–28×.**

**Stated carefully, because the temptation to call this gravity is strong and premature:**

- **WHAT IS MEASURED:** what counts as a record-bearing region is determined by the matter/environment
  content, not fixed in advance. **Content decides geometry** — which is general relativity's central
  move, appearing here for the first time.
- **WHAT IS NOT MEASURED:** the graph is still static. Nothing deforms. What is selected is which
  SUBSET of a fixed complex is a region, not the complex itself. **This is not dynamical geometry
  and must not be reported as such.**

**And it is consistent with the whole night: it happened on a STATIC carrier.** W-31 (refuted 3/3),
V003, and now W-38 all point the same way — **the responsive carrier has not been needed for
anything yet, including for geometry to become an output.**

---

## W-39 — **RECORDS CROWD EACH OTHER OUT. THE CARRIER HOLDS 3, NOT 4, AND PROTECTION DEGRADES WITH EACH ONE.**

`LANE_W39_CAPACITY/`, sealed. Opened on the principal's objection: **dynamical geometry is unlikely
to appear at the level of a SINGLE record — a single record is a test particle, and test particles
never source geometry.** It would emerge where MANY records must share a carrier.

**CAPACITY.** The algebraic bound is `4` and is **FORCED** — it is the cycle rank, i.e. the dimension
of the physical space, and is reported only so it cannot be mistaken for a finding. The geometric
bound — can a set be simultaneously protected, i.e. does a nonempty bath avoid every boundary at
once — is **3**.

> **THE TWO DISAGREE. THE LIMIT IS PACKING, NOT DIMENSION.** `k=4`: **zero** independent sets are
> simultaneously protectable. Example of a pair that provably cannot coexist: regions `(0,3)` and
> `(1,2)`, whose boundaries together cover all 12 links, so **no bath avoids both.**

**AND THE CAPACITY NUMBER UNDERSTATES IT — PROTECTION DEGRADES AS RECORDS ARE ADDED:**

```
   k  |bath|   gamma  gamma*|bath|  worst member  best non-member   margin
   1       8   0.500          4.00    6.7285e-02       2.9460e+00    43.8x
   2       5   0.800          4.00    5.3364e-02       1.0612e-01     2.0x
   3       2   2.000          4.00    5.8497e-02       6.2405e-02     1.1x
```

**THE CONFOUND WAS REAL AND IS CONTROLLED.** The protecting bath necessarily SHRINKS as records are
added (8 → 5 → 2 links), so a weaker environment could explain the collapse by itself. **RUN 2 holds
the total dissipation fixed at `gamma·|bath| = 4.00` and the margin collapses identically: `43.8x →
2.0x → 1.1x`.** The crowding is not the bath weakening. (RUN 3, holding `|bath|` literally constant
at 2, is confounded in the other direction — with so small a bath many regions are protected at once
— and is retained but not read.)

> **AT `k=3` THE MARGIN IS `1.1x`: TWO OF THE THREE RECORDS ARE BARELY DISTINGUISHED FROM
> UNPROTECTED OPERATORS. THE USABLE CAPACITY IS NEARER 2 THAN 3.**

### WHAT THIS DOES TO THE GRAVITY CONCLUSION

**V003 and W-31 concluded that nothing measured REQUIRES a responsive carrier, and the registrar
offered the principal two readings: either records genuinely do not need gravity, or we had not
found the forcing measurement. The principal supplied a third, and it is better than both:
WE HAD NOT BUILT THE STAGE.**

**Every result up to W-38 used exactly ONE record on a passive background.** Nothing there could
have detected backreaction, in the same way a test particle cannot detect that it curves spacetime.

**W-39 is the first measurement in which records constrain EACH OTHER** — they compete for
environment placement, the competition is measured and controlled, and it has a hard limit at 3 and
a practical limit near 2. **That mutual constraint is the precondition for backreaction, and it is
where a responsive carrier could first do work.** Whether it does is untested.

**NEXT.** The carrier's capacity is a property of a FIXED graph. **Ask what happens when the packing
constraint has no solution** — `k=4` is already such a case here. On a static carrier the answer is
that the fourth record simply fails. **A responsive carrier is exactly the thing that could answer
differently, and that is the first question in this program where it could matter.**

---

## W-40 — **REDUNDANCY MEASURED. READING IS FREE; STORING IS NOT.**

`LANE_W40_REDUNDANCY/`, sealed. W-37 claimed condition 5 and **argued** it from `||[R,H_hop]|| = 0`.
W-37e's two-probe attempt returned zero, which was exclusion — hard-core probes never complete a
circuit, and a return trip flips every link twice giving `R² = I`. **Let them pass and measure it.**

```
    n   k   dim        I per probe (bits)          I joint
    4   1     8      0.9333                        0.9333
    4   2    32      0.9333  0.9333                0.9924
    4   3   128      0.9333  0.9333  0.9333        0.9989
    5   3   250      0.7342  0.7342  0.7342        0.9737
    6   3   432      0.8079  0.8079  0.8079        0.9869
```

> **EACH PROBE ALONE HOLDS THE BIT, AND ADDING PROBES DOES NOT COST THE OTHERS ANYTHING — the
> per-probe value is IDENTICAL at `k = 1, 2, 3`, to four decimals.** The joint value climbs toward
> `1` (`0.9333 → 0.9924 → 0.9989`), so the copies are independent samples rather than one shared
> copy. **Condition 5 is measured.**

**WHY THE PER-PROBE VALUE IS EXACTLY UNCHANGED, since an exact identity always deserves suspicion:**
the probes couple to each other only through the gauge field, and what they read — the global flux —
**is conserved (`[R,H] = 0`)**. So each reads the same undisturbed quantity independently. That is
the QND property doing its work, and it is what we set out to verify, not an artifact.

**CONTROLS.** `tau = 0`: `||H|| = 0.000e+00`, nothing moves, nothing is learned. **Real hopping:
`n=4` gives `0.9998`, `n=6` gives `0.9827`, and `n=5` gives EXACTLY `0.0` for every probe** — W-37d's
probe-relativity reproduced here, and it applies per-probe.

### THE ASYMMETRY, AND IT IS THE POINT

> **MANY OBSERVERS CAN READ ONE RECORD WITHOUT INTERFERING (W-40: no per-probe cost at all).
> MANY RECORDS CANNOT SHARE ONE CARRIER WITHOUT DEGRADING EACH OTHER (W-39: margin `43.8x → 2.0x →
> 1.1x`, controlled for bath strength).**
>
> **READING IS FREE. STORING IS NOT.**

The two measurements use the same machinery and differ only in what is multiplied. **Storage is what
has a capacity, and capacity is what could force a carrier to respond.** Observation cannot, and
W-40 is what rules it out rather than leaving it assumed.

**`PROCESS_DESCRIPTION_V003` §1 condition 5 is upgraded from ARGUED to MEASURED**, and keeps W-37d's
qualifier: **objective to probes that can couple to it.**

---

## W-41 — **k=4: THE CARRIER EVICTS. IT DOES NOT SHARE THE DAMAGE.**

`LANE_W41_K4/`, sealed. W-39 established the packing constraint has no solution at `k=4`. Demand it
and watch the failure mode, because the failure mode is the physics.

**A DERIVATION THAT ALSO MAKES THE SCAN CHEAP.** For the adjoint Lindbladian the coherent term
contributes **nothing** to the first-order rate — `Tr(O† i[H,O])` is purely imaginary for Hermitian
`H` — and at `Z_2` a Wilson loop and a link operator either commute or anticommute. So

> **`rate(O) = 2·gamma·|bath ∩ boundary(O)| + O(g⁴)`. The decay rate is literally a COUNT of how many
> bath links sit on the boundary.** W-32, W-38 and W-39 in one line.

**VALIDATED AGAINST THE FULL SPECTRUM before being used:** `|overlap| = 1` → `0.92–0.99`;
`|overlap| = 0` → `0.003–0.06`, the `g⁴` residue.

**THE RESULT.** Of 793 baths, **zero** protect all four — the constraint genuinely has no solution.
And the failure is not shared:

```
  bath [0], confirmed on the full 256x256 spectrum
     region (0,)   rate 9.2242e-01   EVICTED
     region (1,)   rate 3.8571e-02   protected
     region (2,)   rate 3.8571e-02   protected
     region (3,)   rate 2.7401e-03   protected
```

> **THREE RECORDS SURVIVE AT THE `g⁴` RESIDUE AND ONE IS THROWN OUT, A FACTOR OF 24–340x WORSE.
> THE CARRIER HAS A HARD CAPACITY OF 3 AND ENFORCES IT BY EVICTION, NOT BY DEGRADING EVERYONE.**

**AND WHICH ONE DIES IS SET BY GEOMETRY.** Each link lies on 1 or 2 of the four boundaries. A bath
on a **multiplicity-1** link evicts exactly **one** record; on a **multiplicity-2** link it evicts
**two**. **Where the environment sits decides how many records die**, and the count is exact.

### WHAT THIS OPENS

**On a static carrier the fourth record simply fails.** There is no negotiation and no deformation —
the capacity is enforced instantly and cleanly. **This is the first question in the program where a
responsive carrier could give a DIFFERENT answer that would mean something:** instead of evicting,
a carrier that could grow — one more plaquette, one more link — could accommodate the fourth.

**That makes the gravity-shaped mechanism concrete and testable for the first time: CONTENT
EXCEEDING CAPACITY FORCES THE CARRIER TO EXPAND.** The prerequisite measurement is whether capacity
actually grows with carrier size — `4` plaquettes give capacity `3`; what do `6` or `9` give? If
capacity tracks size, "demand exceeds capacity → the carrier grows" is a well-posed dynamics rather
than a wish. **Untested.**

---

## W-42 — **CAPACITY = AREA − 1, EXACTLY. AND IT CORRECTS W-39.**

`LANE_W42_CAPACITY_LAW/`, sealed. Prerequisite for "content exceeding capacity forces the carrier to
expand": capacity must actually depend on size.

W-41's counting formula removes the Hilbert space entirely. The boundary map is **linear over GF(2)**;
a set of records is simultaneously protectable iff **some link lies on none of their boundaries**;
and that is **one linear condition**. So capacity `= max_L dim ker(f_L)`.

```
     patch  verts  links  plaq m  cycle rank  capacity   = m-1
       3x3      9     12       4           4         3    yes
       4x3     12     17       6           6         5    yes
       4x4     16     24       9           9         8    yes
       5x4     20     31      12          12        11    yes
       5x5     25     40      16          16        15    yes
```

**Cross-checked against W-39's explicit brute force at `3x3`: both give `3`.**

> **CAPACITY = m − 1 = AREA − 1**, where `m` is the number of plaquettes, which for a planar patch is
> the area in lattice cells and equals the cycle rank.

### ERRATUM TO W-39

**W-39 reported capacity `3` as a GEOMETRIC bound, contrasted with the FORCED algebraic bound `4`,
and concluded "the limit is packing, not dimension." THAT CONCLUSION IS WITHDRAWN.** `m − 1` is the
kernel dimension of a single linear functional — **it is forced too**, just one less than the space's
dimension rather than equal to it. The registrar drew the forced/not-forced line in the wrong place.

**W-39's other results are unaffected**: the *measured* crowding (margin `43.8x → 2.0x → 1.1x`, with
total dissipation held fixed) is dynamics, not counting, and stands. So does W-41's eviction.

### AND THE HONEST READING, WHICH CUTS AGAINST AN EASY STORY

```
    3x3: area  4  perimeter  8  capacity  3     cap/area 0.750   cap/perimeter 0.375
    4x4: area  9  perimeter 12  capacity  8     cap/area 0.889   cap/perimeter 0.667
    5x5: area 16  perimeter 16  capacity 15     cap/area 0.938   cap/perimeter 0.938
```

> **CAPACITY TRACKS AREA, NOT PERIMETER.** `cap/area → 1` while `cap/perimeter` grows without bound.
> **This is NOT a holographic or boundary-law count and must not be reported as one** — the registrar
> has already conflated two different area laws once in this program (confinement vs Bekenstein) and
> will not do it again.

**WHAT IT BUYS.** The dynamics is now well-posed: **to hold `k` records a carrier needs area `≥ k+1`.**
Below that it evicts (W-41). **"Demand exceeds capacity → the carrier must grow" is a statement with
a number in it**, and the number is one plaquette per record.

---

## CORRESPONDENCE AUDIT — **THE EMPTY GRAVITY LINE IS A PROPERTY OF OUR GRAPH, NOT OF RECORDS**

`GRAVACLE_CORRESPONDENCE_V001.md`, sealed. The principal directed that the program not re-cover
ground already in `10.5281/zenodo.21238968` (Gravacle v337), the guard against fitting being to test
for it afterwards rather than to avoid the source.

**FITTING LINE, VERIFIED BY SEARCH:** before today the project contained **zero** substantive
references to the paper — the only hits are the GitHub org name. **W-27 to W-42 are uncontaminated.
W-43 onward is not, and must be flagged.**

**THE ONE THING THAT CHANGES A CONCLUSION.** The paper's Field-Registration Principle: *EM enters
records through gauge/action phase; gravity enters through metric/proper-time action.* **Our carrier
has gauge phase and a Gauss law, and its links have no length, no metric, no proper time, no clock.**

> **"GRAVITY IS ABSENT FROM THE PROCESS" IS NOT A RESULT ABOUT RECORDS. IT IS A RESULT ABOUT OUR
> GRAPH. We never installed the thing gravity enters through.** V003 and the W-31 refutation both
> reported the empty line as though it were a finding about record formation. **Withdrawn.**

**AND IT DIAGNOSES ROUTE B.** W-31 made the carrier's **topology** dynamical. Gravity is said to
enter through **metric**. **Route B varied the wrong degree of freedom** — which is why three
adversaries found the topology contributed exactly nothing.

**INDEPENDENT CORRESPONDENCES, all produced before reading:** the paper's **Action-Holonomy Sign
Theorem** (the interference sign IS the holonomy of a closed oriented action loop) against W-37's
cut-ring control at exactly `0.000000`; and its **Single-Run Boundary Closure Theorem** (exact
selection is single-run closure, Born capacity is the public ensemble law) against W-35's
`E[<R>] = -0.015 ± 0.040` beside `E[<R>²] = 0.999`. **The second match is precise and it overturned
four of this program's own lanes before the paper was consulted.**

**WHERE THIS PROGRAM SITS IN THE PAPER'S OPEN PROBLEM.** Appendix C names the **accumulation proof**
as unsolved: (1) a local record-to-geometry map, (2) compatibility across overlapping record
patches, (3) the large-record limit. **W-39/41/42 are step 2 with numbers in it.** But our records
carry no geometry, so **capacity constrains how many records fit and says nothing about a metric.
Step 2 has the right shape and no content because step 1 is missing.**

**THE MISSING TERM, CONCRETELY.** Links need a **length / proper-time weight**. Then a transported
probe accumulates a gauge **product round a closed loop** (homotopy-invariant — we have it) and a
**sum of lengths along its path** (not homotopy-invariant — we do not). **Structurally different,
and the difference is measurable: one is blind to path length, the other is not.**

**FITTING RISK FLAGGED IN ADVANCE.** W-37d's complex hopping phase is formally what a proper-time
weight contributes, and it is tempting to claim we already saw gravity. **We did not.** That phase
was **uniform and installed by hand** to break a time-reversal symmetry; a metric weight must **vary
across links and be sourced**. The separating test is uniform versus varying phase, and whether
anything sources the variation.

---

## W-43 — **LOCAL RECORDS ARE OBJECTIVE FOR FREE; THE ONE THE DYNAMICS SELECTS IS NOT**

`LANE_W43_METRIC/`, sealed. **CONTAMINATED LANE — designed after reading Gravacle v337.**

**W-43's chirality test failed twice and both failures are logged rather than hidden.** The measure
was broken (the antipodal site weighted `-n/2` when it is equally `+n/2`, so a symmetric
distribution scored `-0.419` — exactly the "nothing added" baseline). And the concept was wrong: **a
`Z_2` flux cannot break time reversal, because `0` and `pi` are both TR-invariant**, so chirality
could never separate a `Z_2` gauge record from a metric one. **Nothing is inferred from it.**

**W-43b, the measurement that works.** Same environment, same coupling, same machinery as W-36:

```
  GAUGE record  (rim Wilson loop, a product over 8 links)
    |F| : 0      1      2      3      4      5      6
    I   : 0.000  0.044  0.100  0.187  0.351  0.658  0.929      ratio 0.047

  LOCAL record  (a two-level variable at one site)
    I   : 0.000  0.999  1.000  1.000  1.000  1.000  1.000      ratio 0.999

  CONTROL kappa = 0, both:  <= 9.5e-15 at every fragment size
```

### THE HONEST LIMITATION, BEFORE THE RESULT

> **THIS IS NOT A GRAVITY MEASUREMENT.** The "metric record" is a two-level mass at one site with
> `H_sys = 0` and the environment coupled to it **directly**. **The potential it would source never
> enters the dynamics.** So W-43b compares a LOCAL observable against a NONLOCAL one — and that a
> local thing is locally copyable is close to a tautology.

### WHAT IS NOT TAUTOLOGICAL

> **THE RECORD THE DYNAMICS SELECTS IS THE NONLOCAL ONE.** W-34's sieve, nominating nothing, picked
> the full boundary by a factor of **28**. So the carrier protects precisely the record that **cannot
> be made objective**, while the record that is objective for free is **not** the one it protects.

**PROTECTION AND OBJECTIVITY PULL APART, AND THE TWO SIT ON OPPOSITE SIDES.** W-36 named this as
legible-vs-protected; W-43b puts both numbers on one page under identical conditions: `0.047`
against `0.999`. **And the only structural difference is whether the observable lives on a closed
loop or at a point.**

**NEXT, AND IT IS THE REAL TEST.** Make the potential do work: let it source a phase that a
transported probe accumulates, so the metric enters the dynamics rather than being a label the
environment reads off. **Then ask whether THAT is redundantly readable, and whether the sieve still
prefers the loop.** Until then the gravity line has a candidate term and no measurement.

---

## W-45 — **CAPACITY PASSES ALL FOUR OF GRAVITY'S FUNCTIONAL MARKERS, WITH NO METRIC ANYWHERE**

`LANE_W45_FUNCTIONAL/`, sealed. Opened on the principal's standing correction, given twice and
ignored twice by the registrar: **do not expect the full footprint of CLASSICAL gravity at the
record level; it will look different while performing the same function.** *"An early embryo doesn't
look like a human."* The W-44 workflow specifies gravity by metric, proper time and clocks — its
classical form — and therefore imports precisely the assumption we were told to stop importing.

**So test by FUNCTION.** Gravity's signature, stripped of its classical shape, against **capacity**
— a quantity that arose from this program's own measurements (W-39, W-41, W-42) and not from gravity.

```
  F1  UNIVERSAL SOURCING   every record consumes; there is no neutral record
        15 of 15 records consume exactly 1; consumption values observed: [1]         PASS
  F2  NO SCREENING         a record cannot hide behind others and avoid paying
        0 screened cases of 1152 with capacity genuinely available                   PASS
  F3  ONE SIGN             nothing ever frees capacity; there is no negative record
        no record of any size increases what is available to others                  PASS
  F4  ARENA, NOT FORCE     it never appears in H and exerts nothing; it sets what
        is POSSIBLE -- k=1,2,3 possible, k=4 IMPOSSIBLE, no arrangement exists       PASS
```

**F2 FAILED AS FIRST WRITTEN AND THE FAILURE WAS THE TEST.** `108 of 1260` additions appeared to pay
nothing — but in every one of them **capacity was already exhausted**, so there was nothing left to
charge. That is **saturation, not screening.** Gating on capacity `> 0` beforehand gives **zero**
screened cases. Both versions are retained in the lane.

**AND THE CONTRAST IS SHARP.** The `Z_2` charge on the same carrier fails **all four**: it has two
signs, charges annihilate in pairs (`1+1 = 0 mod 2`), and a neutral sector exists (`8` of `16`
configurations are globally neutral). **EM screens. Capacity cannot.**

### WHAT IT IS

> **CAPACITY IS ANOTHER WORD FOR ROOM. And `capacity = area - 1` exactly (W-42), so it is not a
> metaphor for space — IT IS THE AREA, MINUS ONE.**

Everything takes some, nothing makes more, nothing can avoid paying, and it never pushes anything
around — it only decides what can coexist. **That is what space does with respect to matter, derived
here from record dynamics on a pure gauge theory with no metric, no proper time and no clock
installed anywhere.**

**STATED AT ITS PROPER STRENGTH.** This is a functional match, not a metric. Capacity is presently a
**static bookkeeping quantity of a fixed graph** — it has no dynamics of its own, nothing sources a
change in it, and nothing expands. **It is the embryonic form the principal named, and it should be
reported as that and not as gravity.**

**NOT HOLOGRAPHIC, AND THE REGISTRAR WILL NOT CLAIM OTHERWISE.** Capacity tracks the extensive
measure (plaquettes), not the boundary (W-42: `cap/area -> 1` while `cap/perimeter` grows without
bound).

### THE NEXT STEP, AND IT IS NOW MOTIVATED FROM INSIDE

**MAKE CAPACITY DYNAMICAL.** Not "let the graph's topology respond" (Route B, refuted 3/3) and not
"install a metric" (W-44, an imported classical shape). **Let the carrier's ROOM change in response
to demand.** W-41 already measured the alternative: on a fixed carrier, demand beyond capacity is
resolved by **eviction**. So the dichotomy is sharp and it is ours, not borrowed:

> **EITHER THE CARRIER GROWS, OR RECORDS ARE EVICTED. Both are measurable, and only one of them
> looks like an expanding arena.**

---

## W-46 — **STORAGE IS A VOLUME LAW. LEGIBILITY IS A BOUNDARY LAW. BOTH READINGS WERE RIGHT.**

`LANE_W46_HOLOGRAPHY/`, sealed. The principal proposed that gravity and capacity converge exactly at
the Bekenstein/holographic bound. W-42 had measured `capacity = area - 1`, a **volume** law, and the
registrar recorded that as a disagreement. **It was not a disagreement. Two different quantities had
been run together.**

- **STORAGE** — how many independent records a region can HOLD.
- **LEGIBILITY** — how many an observer confined OUTSIDE the region can DISTINGUISH.

**Holography is a claim about what escapes, not about what fits.** Measured, for an `n x n` block of
plaquettes:

```
   |A|  perimeter  STORAGE  LEGIBLE   storage/|A|   legible/perimeter
     1        4        0        1         0.000            0.250
     4        8        3        3         0.750            0.375
     9       12        8        5         0.889            0.417
    16       16       15        7         0.938            0.438
```

> **STORAGE = |A| - 1 = n² - 1 — the BULK.**
> **LEGIBILITY = 2n - 1 = perimeter/2 - 1 — the BOUNDARY.**

`storage/|A| -> 1` while `legible/perimeter -> 1/2`. **Storage scales with area, legibility with
perimeter, on the same regions in the same measurement.**

### WHAT THIS SETTLES

> **A REGION CAN HOLD MORE THAN IT CAN EVER TELL YOU, AND THE GAP GROWS WITH ITS SIZE.** At `n=4`,
> `15` records are storable and `7` are externally distinguishable. The ratio diverges as `n`.

**AND IT IS THE SAME FACT THIS PROGRAM KEPT MEASURING FROM OTHER DIRECTIONS.** W-36: a gauge record
has no redundancy in a local environment (`0.047`). W-43b: `0.047` against `0.999` for a local
observable, identical conditions. W-37: reading requires transport around a CLOSED path, and cutting
the ring gives exactly `0.000000`. **All four are the boundary law seen from different sides.**

**THE CORRECTION IS THE REGISTRAR'S.** W-42 stated "capacity tracks area, not perimeter — this is
NOT a holographic count" and warned against conflating two area laws. **That warning was right about
STORAGE and wrong to close the question**, because the holographic claim was never about storage.

**FORCED-OR-NOT, STATED PLAINLY.** This is a rank computation, so it IS a counting fact: only the
boundary layer of `A` pairs with anything outside, so the rank is the perimeter's. **That is not a
defect — it is the mechanism.** Holography here is the statement that *the outside touches only the
edge*, and the arithmetic says so exactly.

---

## W-47 — **THE COMPOSED LAW: FORM SURVIVES, CONSTANT IS OFF BY ~4x**

`LANE_W47_COMBINED/`, sealed. **The first lane in this program that composes measured relations into
a single prediction and tries to break it.** Every prior lane measured one relation.

**THE COMPOSITION.** EM gives a holonomy on a closed boundary of perimeter `P` (W-27, W-34). Alpha
gives its decay rate `Gamma = 32 g^4` (W-32) and its readout time `T_read = 0.30 P/tau` (W-37b).
Capacity gives the boundary law `P/2 - 1` for what escapes (W-46). A record is exportable only if it
is read before it decays:

> **`X = Gamma · T_read = 9.6 · g⁴ · P / tau  <  1`**
>
> **The falsifiable content is that the surviving fraction depends on `X` ALONE** — three parameters
> must collapse onto one curve.

```
   X band                        n    mean I/I0     sd
   X < 0.5   predicted LEGIBLE   16     0.4826    0.2130
   0.5 - 1.5 transition          15     0.1599    0.1051
   1.5 - 4.0                     16     0.0305    0.0259
   X > 4     predicted ILLEGIBLE 17     0.0036    0.0032

   spread at fixed X:  X=0.37 -> 0.394   X=1.0 -> 0.169   X=2.7 -> 0.022   X=20 -> 0.001
```

**WHAT PASSES.** The surviving fraction falls **130x monotonically** across the `X` bands, and points
at the same `X` from different `(P, tau)` agree ever more closely as `X` grows — absolute spread
`0.394 -> 0.001`. **The composition has the right form: readability is gated by the single product
`Gamma·T_read`, not by `P`, `g` and `tau` separately.**

**WHAT FAILS, AND IT WAS DECLARED IN ADVANCE AS THE WEAKER FAILURE MODE.** Half the readability is
lost at **`X ~ 0.25`**, not at the predicted `X = 1` — **the constant is off by about a factor of 4.**
And the RELATIVE spread does not shrink (`sd/mean` runs `0.44, 0.66, 0.85, 0.89`), so the collapse is
**suggestive, not clean.**

> **VERDICT: the composed law is CONFIRMED IN FORM AND NOT IN CONSTANT.** It is one relation binding
> EM, alpha and the capacity/boundary law, it makes a prediction, and the prediction is about a
> quarter right in magnitude. **That is a partial pass and is registered as one.**

**TWO DEFECTS, BOTH CAUGHT BEFORE READING ANYTHING.** v1's decay channel was `P+ - P-`, which **is**
`R` — W-30a proves such a channel cannot move the record, and indeed `I` was identical to six
decimals across every `g^2`. v2 compared raw `I` across `n` when bare readability differs by `n`
independently of decay. **Both are in the lane conventions; the reported numbers are v3.**

---

## W-48 — **THE DECOHERENCE HALF OF OUR LAW IS TEXTBOOK. THE POWER LAW WAS OURS, NOT NATURE'S.**

`LANE_W48_VS_STANDARD/`, sealed. Run **before** touching any experimental data, on the principal's
own guard: if the composed law is `exp(-Gamma T)` in lattice clothing, fitting it to C70 or Panda
data would be curve-matching a textbook result and calling it confirmation.

```
  I_max   (max over t -- what W-47 reported)      n=48
     EXPONENTIAL  rms log-resid 1.2783
     POWER LAW    rms log-resid 0.6425      <- power law wins by 0.64

  I_fixed (at T_read, NO optimisation)            n=44
     EXPONENTIAL  rms log-resid 1.0229      <- exponential wins by 0.08
     POWER LAW    rms log-resid 1.0991
```

> **THE POWER LAW EXISTS ONLY IN THE OPTIMISED SIGNAL. Taking a max over time of a decaying-but-
> still-rising quantity manufactures it. The un-optimised signal is exponential in `X`, which is
> exactly what standard decoherence predicts.**

**SO ON THIS AXIS OUR COMPOSED LAW ADDS NOTHING.** `X = Gamma·T_read` gating legibility is
`exp(-Gamma T)` with `T` fixed by the loop's perimeter — and that a holonomy is unreadable once the
ring exceeds the coherence length is textbook mesoscopic physics. **Taking W-47 to C70 or Panda data
would have fitted a known result and called it confirmation. It is not done.**

**HONEST ABOUT THE STRENGTH.** Exponential beats power law by only `0.08` in rms log-residual, and
**both fits are poor** (`~1.0` in log-residual is a factor of `e` of scatter). The correct statement
is *consistent with exponential, marginally better than power law, neither fitting well* — **not
"confirmed exponential."**

### WHAT THIS LEAVES, AND IT IS A SHARPENING RATHER THAN A LOSS

**The composed law has two halves and they are not equally novel.**

- **THE ALPHA / DECOHERENCE HALF IS STANDARD.** Decay rate times readout time gating survival. W-48
  says so plainly. **Do not take it to data as though it were new.**
- **THE CAPACITY HALF IS NOT A DECOHERENCE STATEMENT AT ALL.** `storage = area - 1`,
  `legibility = perimeter/2 - 1`, capacity enforced by **eviction**, and capacity satisfying all
  four of gravity's functional markers while EM's charge fails all four. **Standard decoherence says
  nothing whatever about how much a region can hold versus export.**

> **THE NOVEL CONTENT, IF THERE IS ANY, IS ON THE GRAVITY/CAPACITY SIDE — AND IT HAS NEVER BEEN
> TESTED AGAINST ANYTHING.**

**AND THAT CHANGES WHICH LITERATURE IS RELEVANT.** C70 and Panda measure coherence against transit —
the half we now know is standard. **The storage-versus-legibility claim belongs to topological
codes**, where the count of logical operators, their nonlocality, and the readout of nonlocal
observables are the actual physics. That is a different corpus from the one in the paper's registry.

---

## W-44 — **A METRIC THAT GENUINELY DOES WORK BUYS NOTHING. AND IT KILLS W-43b.**

`LANE_W44_METRIC_DYNAMICS/`. **CONTAMINATED LANE, and specified wrongly on purpose-in-hindsight:**
it defines gravity by its CLASSICAL form — metric, proper time, clocks — which is the import the
principal had already flagged twice. **Read under that rule: a classical-shaped null reads two ways.**

Three independent builds (link length, site potential, clock), registrar's expectation withheld,
then three adversaries. **The adversaries verified independently that the metric really does work
this time:** it enters `H` only multiplied by the coupling knob, no jump operator touches it
(`max ||[N_M⊗I, L_k]|| = 0.000e+00`), and coupling-zero reproduces separately built no-metric
carriers to `~1e-14`. **That is precisely where W-43b failed, and it is not repeated.**

### REDUNDANCY — UNANIMOUS NULL, AND IT OVERTURNS W-43b

```
  single-fragment / whole-environment Holevo ratio, same state, same fragments
  build                                 metric    gauge
  link-length, length on a rim link     0.0758   0.0478
  link-length, length on a cut link     0.0385   0.0509
  site-potential, env on cut links      0.0090   0.0046
  site-potential, env at the mass       0.0013   0.0835
  clock                                 0.0135   0.0105
```

> **EVERY VALUE IS FAR BELOW `0.5`, AND THE ORDERING FLIPS WITH PLACEMENT. A metric degree of freedom
> that actually participates in the dynamics acquires NO redundancy, and NO advantage over the gauge
> loop.** Three independent constructions agree.

**W-43b's HEADLINE IS WITHDRAWN.** Its metric record scored `0.999` — a full bit from a single
fragment — **because the environment was coupled directly to the mass operator while `H_sys = 0`.**
The registrar flagged that limitation in the lane conventions and in the register at the time. **The
flag was right and the number was an artifact.** With the metric doing real work it scores `0.0013`
to `0.0758`.

**WHAT SURVIVES OF W-43b:** nothing of the metric side. **W-46 is unaffected** — storage-vs-legibility
is a GF(2) rank computation, not a redundancy measurement.

### SELECTION — THE RIM STILL WINS, BUT THE CONTEST IS TWO OPERATORS

At every reported operating point the slowest operator is still the full-boundary Wilson loop
tensored trivially on the metric (`3.16e-03` clock, `2.94e-03` site, `1.66e-02` link). **But an
adversary counted the bath commutant directly: only `8` of `64` dictionary operators commute with
all four jump operators — `56` are excluded by COUNTING before any dynamics runs.**

**And the winner is not robust.** Sweeping the metric's tunnelling rate `eps`, which the builds
declared unmotivated and never swept:

```
  lam=0.6, g2=0.01:  eps=0.01  RIM 3.164e-03  I(x)sz 4.954e-04   METRIC WINS
                     eps=0.02  RIM 3.164e-03  I(x)sz 1.978e-03   METRIC WINS
                     eps=0.05  RIM 3.164e-03  I(x)sz 1.217e-02   rim wins  3.85x
                     eps=0.15  RIM 3.163e-03  I(x)sz 9.192e-02   rim wins 29.06x  <- reported
```

At `eps=0.02` the slowest non-trivial left-eigenmode of the full `1024x1024` Lindbladian matches
`I⊗sz_M` at overlap `0.975`, so **this is not a dictionary artifact.** And the margin itself is
estimator-dependent: `2.56` at `T=1`, `30.9` at `T=10`, `3.95` at `T=1000`.

> **"THE SIEVE SELECTS THE GAUGE BOUNDARY EVEN WITH A METRIC PRESENT" IS A REPORT ON THE RATIO
> `eps/g²`, NOT A ROBUST FACT.** The winner is stable in `T` but not in `eps`.

### DEFECTS THE ADVERSARIES FOUND

**Site-potential transport stage: REFUTED** — a controlled-phase between the readout register and
the projector onto the very quantity whose Holevo is reported. **Clock section 5: REFUTED twice** —
`np.diag(np.arange(DC)) - (DC-1)/2.0` subtracts from **every matrix element**, not the diagonal
(shipped `<C> = -4.500` where the design requires `0.000`), and the fringe-contrast statistic is
blind: with the bug fixed, contrast is `1.000000` at every `omega` while `S(clock)` rises to `1.61`
bits. **Link-length transport: weakened** — the length register is frozen there, so the reading is a
probe seeing a static classical defect. **All three failures are in the third, UNREQUESTED stage.**

### THE READING, AT ITS PROPER STRENGTH

**Installing a metric — three ways, with the metric verified to participate — changed nothing about
redundancy and does not robustly change what the carrier protects.** Under the principal's standing
rule this reads two ways: either a metric is not what performs gravity's function at the record
level, **or** our carrier is not one where a classical-shaped metric can do anything. **It is not
evidence that gravity is absent, and it is not a baseline.**

**AND IT LEAVES W-45 STANDING AS THE BETTER ROUTE.** Capacity — which required no metric at all —
passes all four functional markers. **The classical-shaped install bought nothing; the
functional identification bought four for four.**

---

## W-49 — **CAPACITY ACQUIRES DYNAMICS. THE CARRIER'S SIZE IS SET BY WHAT IS RECORDED IN IT.**

`LANE_W49_CAPACITY_DYNAMICS/`, sealed. W-45 left capacity passing all four of gravity's functional
markers as **static bookkeeping with no equation of motion**. This gives it one.

### A METHODOLOGICAL FINDING THAT REACHES BACK THROUGH THE PROGRAM

**w49 failed and its controls diagnosed why: `<Nhat> = 0.50000` at every `mu` from `0` to `8`.** The
energy cost did nothing.

> **EVERY BATH IN THIS PROGRAM, FROM W-28 ONWARD, USED UNITARY JUMP OPERATORS. A unitary jump
> operator is unital: it drives the state toward maximally mixed. THAT IS AN INFINITE-TEMPERATURE
> BATH.** It is adequate for the dephasing and pointer questions the program has been asking, **but
> an energy cost cannot register in it, so no energetic question could ever have been asked.**

(A second failure: extracting a steady state as a null eigenvector returns an arbitrary element of a
degenerate manifold, which need not be a state — the `Delta=0` control returned `<Nhat> = -0.208`,
impossible for a projector. Steady states are now obtained by **evolving a physical state**.)

### WITH A THERMAL BATH, AND THE RECORD ACTUALLY CONSERVED

A Davies bath with detailed balance at temperature `T`, coupled to the **structure alone**
(`A = I_gauge ⊗ sigma_x`), so every plaquette flux is exactly conserved and records persist.
`<Nhat>` then tracks Gibbs to `1e-9`, and `T -> infinity` recovers the old `0.5`.

**The record must be conserved or the sectors mix and the test is void.** `g²·ELEC` does **not**
commute with the plaquettes (`||[R2,H]|| = 1.265` at `g2=0.05`), so the first pass leaked and showed
a spurious `4e-4`. At `g2 = 0` the commutator is `0.000e+00`:

```
     mu     T    g2  gate   <N | R2=-1>   <N | R2=+1>    difference
   1.00   0.5  0.05  1.00      0.760505      0.760886     3.811e-04   <- leaking, void
   1.00   0.5  0.00  1.00      0.011863      0.852879     8.410e-01
   2.00   0.5  0.00  1.00      0.005835      0.500000     4.942e-01
   1.00   0.5  0.00  0.50      0.035575      0.500000     4.644e-01
   1.00   0.5  0.00  0.00      0.147121      0.147121     1.943e-16   <- CONTROL
```

> **THE CARRIER IS ESSENTIALLY ABSENT (`0.012`) OR ESSENTIALLY PRESENT (`0.853`) DEPENDING ON A
> SINGLE BIT OF RECORD CONTENT. Remove the structure-content coupling and the difference is exactly
> zero.** Content determines geometry, dynamically, measured.

### AT ITS PROPER STRENGTH

**Nothing in `H` references a record, a boundary or a capacity.** The only structure-content coupling
is the gating `-Nhat ⊗ (W2 + h.c.)` — **which is what it MEANS for the carrier to have that
plaquette**, not an extra assumption. No rule says "grow when full". The correlation is a consequence.

**But this is thermodynamics, not new physics.** A degree of freedom whose energy depends on a
conserved quantity will have a conserved-quantity-dependent equilibrium; that is Boltzmann.
**What it changes is the program's ledger, not the physics literature:** item 14 moves from *static
bookkeeping* to *has an equation of motion*, and this is the first construction here in which the
carrier's structure responds to its content at all.

**AND IT ANSWERS W-41's DICHOTOMY.** *Either the carrier grows, or records are evicted.* **On a
carrier that can pay for growth, it grows** — by `0.84` of a plaquette for one bit of content.
Eviction was the fixed-carrier answer, not the only one.

---

## W-50 — **THE FOUR MARKERS ARE FORCED BY COUNTING. W-45 IS WEAKER THAN IT WAS REPORTED.**

`LANE_W50_NECESSITY/`, sealed. Proof obligation **C**: showing records *can* form this way is not
showing this is *how* they form. The gravity identification rests on capacity satisfying four
functional markers, so strip the structure that produced them.

**Two carriers, same records, same cycle rank, same physical dimension `16`, differing only in
whether record boundaries share links:**

```
                        CONNECTED (3x3 patch)      DISCONNECTED (4 separate squares)
  link multiplicity     [1, 2]  SHARED             [1]  DISJOINT
  capacity              3                          3
  F1 universal          PASS  consumption [1]      PASS  consumption [1, 3]
  F2 no screening       PASS  0 of 1152            PASS  0 of 792
  F3 one sign           PASS                       PASS
  F4 arena              PASS  k=4 impossible       PASS  k=4 impossible
  W-41 eviction         [0,0,0,1] -> 3 kept, 1 out [0,0,0,1] -> 3 kept, 1 out
```

> **ALL FOUR MARKERS SURVIVE INTACT WHEN THE SHARING IS REMOVED, AND SO DOES THE EVICTION PROFILE,
> EXACTLY. The gravity-like character does not depend on the carrier's structure at all.**

### THE CORRECTION, AND IT IS TO THE REGISTRAR'S OWN HEADLINE

**W-45 reported "capacity passes all four of gravity's functional markers" as though it were a
substantive structural finding. It is not.** All four follow from one fact: *independent records
live in a finite-dimensional space and each uses one dimension.*

- **F1** every record consumes — because every independent record uses a dimension
- **F2** no screening — forced by linearity
- **F3** one sign — adding constraints can only shrink a kernel
- **F4** arena — capacity is finite

> **THAT IS LINEAR ALGEBRA, NOT PHYSICS. ANY FINITE RESOURCE PASSES ALL FOUR.** Disk space passes
> all four. **The markers are necessary and nearly vacuous.**

**WHAT REAL CONTENT SURVIVES.** The contrast with EM's charge, which fails all four, is genuine — a
**resource** and a **two-signed conserved charge** are different kinds of object, and gravity is
resource-like where charge is not. **But that is far more modest than "capacity has gravity's
functional signature," and the registrar presented the stronger claim.**

**AND IT NAMES WHAT THE MARKERS CANNOT SEE.** Attraction. The equivalence principle. A field
equation. Any relation between content and a metric. **We have none of these, and the four markers
were never capable of detecting their absence.**

**ONE STRUCTURAL DIFFERENCE DID SHOW.** Consumption values are `[1]` connected and `[1, 3]`
disconnected — a record spanning several components blocks more links. So the carrier's structure is
visible in *how much* is consumed, though not in *whether* the markers pass.

### WHAT THIS DOES TO PROOF OBLIGATION C

**Partially answered, in the deflating direction.** The gravity-like structure is **necessary in the
sense of unavoidable** — no carrier with `m` independent records can lack it. **But its necessity is
trivial: it follows from counting rather than from anything about the carrier, the field, or the
process.** So it does not establish that anything gravity-like is *required for records* in a way
that carries physical content. **Item 14 must be read at this reduced strength.**

---

## MAP ERRATUM — **STEP 3 IMPORTED THE CLASSICAL FOOTPRINT AGAIN**

`MAP_TO_THE_RELATION_V001.md` step 3, as first written, read *"our capacity story is repulsive,
gravity is attractive"* and proposed that finding no attraction would be a **standing objection** to
the gravity identification.

> **THAT MAKES CLASSICAL ATTRACTION THE STANDARD THE RECORD SURFACE MUST MEET.** It is the same
> import the principal has now flagged **three times** — classical-shaped nulls read two ways
> (before W-44), *"an early embryo doesn't look like a human"* (before W-45), and again now.
> **The registrar has made this error at W-44, at W-45's framing, and here.**

**CORRECTED.** The functional content of attraction is **positive feedback**: content makes more
content more likely in the same place. The step now asks that, with all three outcomes read in
advance and **none of them treated as a test gravity must pass**:

- **positive feedback** — content concentrates
- **negative feedback** — already measured (W-39 crowding, W-41 eviction); **NOT a refutation**,
  since crowding may be the record-level form of something that only becomes attraction after
  accumulation
- **neither** — the strongest of the three, since it would mean the resource has no dynamics at all

**Steps 1, 2 and 4 were checked for the same defect.** Step 1 (reconstruct adjacency from records)
and Step 2 (alpha as resolution limit) carry no classical gravitational shape. **Step 4 does head
toward classical geometry, and that is correct** — convergence to classical geometry at large scale
is the expected emergent endpoint, which is a different thing from demanding classical form at the
record level, and it is last in the order for that reason.

---

## W-52 — **SPATIAL STRUCTURE IN RECORD FORMATION IS REAL, CREATED BY ALPHA, AND ITS SIGN OSCILLATES**

`LANE_W52_FEEDBACK/`, sealed. Map step 3 in its corrected form: not *do records attract* (the
classical import) but **does a record make another more likely nearby**.

Records `p0` and `p1` share one link; `p0` and `p3` share none. Bath `[4,9]` is the set of links
avoided by all three, so every comparison runs under the same bath.

```
   g^2    d(neighbour)         d(far)            difference        sigma
  0.00  0.03629 +-0.00305  0.03584 +-0.00313   +0.00045 +-0.00437    0.1
  0.02  0.03072 +-0.00382  0.04565 +-0.00330   -0.01493 +-0.00504    3.0
  0.05  0.09756 +-0.00429  0.05370 +-0.00492   +0.04386 +-0.00653    6.7
  0.10  0.07516 +-0.00522  0.16637 +-0.00591   -0.09121 +-0.00788   11.6
  0.20  0.01712 +-0.00309  0.06358 +-0.00365   -0.04645 +-0.00478    9.7
```

**TWO THINGS ARE ESTABLISHED AND ONE IS REFUSED.**

**ESTABLISHED — the spatial structure is real and it is ALPHA that creates it.** At `g² = 0` the
difference is `0.1 sigma`: **exactly indifferent, as the mechanism requires.** The mechanism was
stated before the run and holds: plaquette operators commute, `||[H(g²=0), W_p1]|| = 0.000e+00`, so
at zero coupling the records are independent degrees of freedom and **no spatial structure is
possible**. Switch the electric term on and `||[H(0.05), W_p1]|| = 0.980` — **it is the only channel
that couples plaquettes, and neighbouring plaquettes couple through the link they share.**

**ESTABLISHED — the effect is large.** Up to `11.6 sigma`, six independent seeds.

**REFUSED — that it is FEEDBACK in either direction. THE SIGN OSCILLATES WITH THE COUPLING:**
negative at `0.02`, positive at `0.05`, negative again at `0.10` and `0.20`.

> **THE REGISTRAR RAN `g² = 0.05` ALONE FIRST AND WOULD HAVE REGISTERED "POSITIVE, SPATIALLY
> STRUCTURED FEEDBACK". It is one point of a sign-oscillating quantity.** The likely cause is the
> same class as two earlier defects in this program: a **fixed sampling time** (`T = 15`) reading a
> coherently precessing quantity — cf. W-32's fixed-`T` survival and W-47's max-over-time.

**AND THE PRE-REGISTERED OUTCOME LIST WAS INCOMPLETE.** The map named three readings in advance —
positive, negative, indifferent. **The measurement returned a fourth: oscillating.** Pre-registering
outcomes is the right discipline and it did not save the reading; **only the sweep did.**

**WHAT IT LEAVES.** Records are **not** indifferent to one another — that possibility is dead at
`11.6 sigma`. **What couples them is alpha, through shared links, and nothing else.** Whether the
coupling concentrates or disperses content is **not answered**, and cannot be until the measurement
is made at something other than a fixed time.

---

## PHASE A — **THREE OF THE FOUR NECESSITY RESULTS ARE NOW THEOREMS, NOT MEASUREMENTS**

`THEOREMS_V001.md` + `LANE_PA_THEOREMS/`, sealed. The program had **measured** four results that are
provable. Measuring a theorem on one lattice is strictly weaker than proving it.

**T1 — durability ⟹ unwritable. PROVED, AND STRONGER THAN WHAT WAS MEASURED.** For any GKSL
generator with **arbitrary** jump operators and `R` **normal**: `[H,R]=0` and `[L_k,R]=0` give
`d⟨R⟩/dt = 0` for every state. **No lattice, no gauge group, no dimension, and no unitarity of the
jumps enters** — the program had only ever measured the unitary-jump case. Verified at the edges
with random operators and exact simultaneous diagonalisation: `≤ 1.2e-15` across `D = 4…12`, with
each hypothesis breaking it when dropped. **And normality is nearly free: if `H` is Hermitian with
non-degenerate spectrum and `[H,R]=0`, then `R` is automatically normal.**

**T2 — reading requires a closed path. PROVED.** An open path's holonomy transforms as
`W → g_a W g_b⁻¹` and admits no gauge-invariant function of the link variables alone; a closed one
transforms by conjugation. **W-37's cut-ring `0.000000` is an instance of the theorem, not evidence
for it.**

**T3 — `capacity = m − 1`. PROVED FOR ARBITRARY COMPLEXES.** `bd` is linear over `GF(2)`, so
protectability via a link `ℓ` is the kernel of ONE linear functional. Verified on the `3×3` patch,
disconnected squares, a **tetrahedron**, a **cube surface**, and six random complexes — **10 of 10**,
with the orphan-link hypothesis confirmed as the one doing the work. **Planarity, dimension, lattice
structure and gauge group play no role.**

**T4 — `legibility = perimeter/2 − 1`. NOT PROVED, and it is the one Phase C depends on.** The
mechanism is clear (only a region's boundary layer pairs with anything outside) but the coefficient
is underived and the class of regions unestablished.

> **WHAT PHASE A BUYS: three of the program's four necessity results are no longer statements about
> a `Z_2` lattice.** T1 in particular is a general fact about open quantum systems. **What it does
> not buy: the carrier is still imported, and T4 — the one with any prospect of a distinguishing
> prediction — is still a measurement on rectangular planar blocks.**

---

## PHASE A / T4 — **W-46's `2n−1` IS WITHDRAWN. THE BULK/BOUNDARY SPLIT SURVIVES AND IS STRONGER.**

Reducing the pairing shows legibility is the **GF(2) rank of the outside–inside plaquette
adjacency**, whose kernel is exactly the region-combinations whose boundary is interior to `A`.

```
    n   |A|=n²   legibility   invisible   (n−2)²
    2       4            4           0        0
    3       9            8           1        1
    4      16           12           4        4
    5      25           16           9        9
```

**`invisible = (n−2)²` exactly — the plaquettes touching no outside plaquette. `legibility = 4(n−1)`.**

**W-46 REPORTED `2n−1 = perimeter/2 − 1`. WITHDRAWN. Its region sat at the lattice CORNER**, so part
of `A`'s boundary faced the lattice edge instead of outside plaquettes. Re-running that geometry
reproduces `3, 5, 7`; the same regions placed in the interior give `4, 8, 12`.

> **THE CONCLUSION IS UNCHANGED AND BETTER FOUNDED: storage scales with AREA, legibility with
> PERIMETER — and the mechanism is now exact. What a region cannot export is precisely its interior.**

**T4 REMAINS OPEN AS A CLOSED FORM.** `legibility = boundary-touching count` is exact for rectangles
and **false in general**: the `plus` pentomino has all five cells touching the outside and legibility
`4`; `3×3`-minus-centre has all eight touching and legibility `8`. The kernel characterisation is
general; the formula is not.

**PHASE A CLOSES AT THREE THEOREMS AND ONE CHARACTERISATION.** T1 (any GKSL generator, arbitrary
jumps, `R` normal), T2 (gauge theory), T3 (arbitrary complexes) are proved. T4 has an exact
mechanism and a shape-dependent coefficient.

---

## PHASE B — **UNIVERSALITY: T1 AND T2 ARE GENERAL. T3 IS ABELIAN-SPECIFIC BY PROOF.**

`LANE_PB_UNIVERSALITY/`, sealed.

**B1 — T1 needs no Markov approximation at all.** Modelling system **and** environment explicitly and
unitarily: if `R` acts on the system, `[H_sys,R]=0` and `[H_int, R⊗I]=0`, then `R⊗I` commutes with the
full `H_tot`, since `H_env` acts on the other factor. **`⟨R⟩` is then exactly conserved with no weak
coupling, no memorylessness and no bath spectral assumption.** Verified with an explicit **finite**
environment (so it recurs and has memory) at strong coupling `3.0`: drift `≤ 6.6e-13` over `t = 0…40`
at `D_S×D_E` up to `5×8`; each hypothesis dropped breaks it. **GKSL was never required — T1 is a
conservation statement.**

**B2/B3 — T3 is field- and dimension-independent.** `capacity = m−1` verified over **GF(2), GF(3),
GF(5)** on plaquette and triangle complexes, and in **three dimensions** with regions as cubes and
boundaries as faces (`2×1×1`, `2×2×1`, `2×2×2`). The proof used only linearity of the boundary map
over a field, and neither the characteristic nor the dimension appears.

**WHERE T3 STOPS, ASSESSED NOT TESTED.** It needs an **abelian** group: a non-abelian holonomy depends
on ordering and base point, so "the product of the plaquettes" is not a linear function of the region
and the kernel argument fails. **Untested and open.** **T2 is group-independent** — its proof
(`W → g_a W g_b⁻¹`) holds for any `G`.

---

## PHASE C — **NO DISTINGUISHING PREDICTION. THE SURVIVING CANDIDATE IS LOCALITY.**

`LANE_PC_PREDICTION/`, sealed. W-48 had already killed the decoherence half of the composed law by
showing the un-optimised signal is `exp(−ΓT)`. **The one surviving candidate was the
storage/legibility split.** The plan's obligation was to establish it predicts something the existing
account does not, **before** touching data.

**The deflation to rule out: an operator supported strictly inside a region commutes with everything
outside it. That is microcausality.** So: is legibility exactly the number of independent records
whose support touches the region's boundary?

```
   shape              |A|   legibility   locality bound   differ?
   square 3x3           9            8                8
   square 4x4          16           12               12
   rect 3x4            12           10               10
   plus pentomino       5            4                4
   3x3 minus centre     8            8                8
   ... 13 shapes tested, differences observed: [0]
```

> **LEGIBILITY IS THE LOCALITY BOUND IN EVERY CASE. THE STORAGE/LEGIBILITY SPLIT IS MICROCAUSALITY
> EXPRESSED IN GAUGE-THEORETIC LANGUAGE.** It predicts nothing the existing account lacks, and it is
> **not** taken to data.

### WHAT THIS MEANS FOR THE PROGRAM, STATED WITHOUT SOFTENING

**Phase D is not reachable.** It was contingent on C producing a prediction that could fail, and C
produced none.

**Every component of the account now reduces to known physics.** The decoherence half is `exp(−ΓT)`
(W-48). The capacity half is locality (C1). Single-run versus ensemble is unravelling-dependence,
standard. W-49's carrier response is Boltzmann. W-52's spatial structure is the electric term coupling
plaquettes that share a link.

> **THE PROGRAM HAS PRODUCED A CORRECT, INTERNALLY VERIFIED, THEOREM-BACKED ACCOUNT OF RECORD
> FORMATION IN LATTICE GAUGE THEORY THAT IS ENTIRELY COMPOSED OF KNOWN PHYSICS. There is no novel
> prediction in it.** That is the ending `PLAN_TO_PROOF_V001.md` named in advance as possible, and it
> is the one that occurred.

**AND THE REASON IS STRUCTURAL, NOT ACCIDENTAL.** The only slot that was ever empty is **gravity**,
and three independent routes to it failed — dynamical topology (refuted 3/3), a classical metric
(three builds, bought nothing), and the functional identification (forced by counting). **An account
built from EM and a coupling alone can only reproduce what EM and a coupling already explain.**

---

## PHASE C ERRATUM — **C1 TESTED LEGIBILITY, NOT CAPACITY. THE VERDICT WAS OVERSTATED.**

The Phase C entry wrote *"the capacity half is locality (C1)"*. **C1 tested no such thing.** It tested
**LEGIBILITY** — how many of a region's records are distinguishable from outside — and found it equals
the locality bound in 13 of 13 shapes. **That is a result about what ESCAPES a region, not about how
many records it can HOLD.** The two are separate claims and the registrar collapsed them.

**WHAT IS ACTUALLY ESTABLISHED, ITEM BY ITEM:**

| claim | status | by |
|---|---|---|
| legibility (what escapes) reduces to locality | **established** | C1, 13 shapes |
| capacity's four gravity-like markers are forced by counting | **established** | W-50 |
| the carrier's response to content is Boltzmann | **established** | W-49 |
| **capacity + boundary formation are not gravity's record-level form** | **NOT ESTABLISHED** | — |

**The last line was never tested and the Phase C verdict implied it.**

### THE PRINCIPAL'S FRAMING, AND WHAT WOULD ACTUALLY TEST IT

> *"we don't expect gravity to appear at the record level looking the same as classical gravity. It
> will likely have other features such as capacity and boundary formation."*

Everything this program has ruled out was ruled out **at a single scale**: one region, one carrier,
one bath placement. **The four markers being forced by counting says they carry no information ABOUT
THE CARRIER. It does not say a large collection of capacity constraints carries no information.**

**THE UNTESTED ROUTE IS ACCUMULATION** — and it is the paper's own Appendix C step 3, the large-record
limit, which this program has never attempted. W-39/41/42 are step 2 (compatibility across
overlapping regions) with numbers in them. **Step 3 asks whether a sufficiently rich collection of
capacity and boundary constraints determines something that is not present in any one of them.**

**Why this is not special pleading.** Every deflation so far — `exp(−ΓT)`, locality, Boltzmann,
counting — is a statement about ONE constraint. **Emergence claims are claims about MANY, and none of
the tests run so far could have detected one.** That is a real gap in the coverage, not a rescue.

**IT IS ALSO FALSIFIABLE.** If a large collection of capacity/boundary constraints determines nothing
beyond what each determines separately, accumulation adds nothing and the negative verdict stands at
every scale. **That is a decidable question and it has not been asked.**

---

## W-51 — **GEOMETRY IS NOT IN THE RECORDS. IT IS IN THE ENVIRONMENT COUPLING.**

`LANE_W51_RECORD_TO_GEOMETRY/`. Map step 1: can the carrier's geometry be reconstructed from records
alone. Three independent methods, expectation withheld, three adversaries. **All three reported
fidelity `1.000` on the true carrier against `~0.5` shuffled. The headline is nonetheless withdrawn.**

### THE RECONSTRUCTION IS AN EXACT IDENTITY, SO NOTHING COULD HAVE FAILED

At `g² = 0` the Hamiltonian is a sum of records and commutes with all of them, so
`Γ(S) = 2γ·|bdy(S)|` exactly. Boundaries are GF(2)-additive, hence
`|bdy(p⊕q)| = |bdy(p)| + |bdy(q)| − 2|bdy(p)∩bdy(q)|`, so **the pair excess IS the shared-link count
by polarisation** — and "shares a link" is by definition that intersection being non-empty. **Two
lines of set algebra, no dynamics.**

**VERIFIED BY EXECUTION, NOT ARGUMENT.** A reconstructor with **no Hamiltonian, no Lindbladian and no
time evolution — just set XOR** — reproduces the lane exactly: block and ring exact, and across 200
count-preserving shuffles it matched **each shuffle's own truth 200/200**, scoring mean `0.5057`
against the original. **Those are the lane's published numbers.**

> **THE MANDATED SHUFFLE CONTROL HAS NO DISCRIMINATING POWER.** A physics-free reconstructor passes it
> with an identical signature. One adversary scored 60 shuffles with and without dynamics:
> `max|difference| = 0.0e+00`.

**Build 1 (spectral) is circular outright** — its input contains the full incidence matrix
(`rank(B) = L`, so `A` follows from one least-squares solve; recovered at `2e-15` by two adversaries
working independently, recovering *which* links, not just how many).

### THE RESULT THAT SURVIVES IS A NULL, AND IT IS THE STRONGEST CLAIM IN THE LANE

> **ADJACENCY IS NOT IN THE RECORD ALGEBRA.** Every record is a product of link operators, all
> commutators are `0.00e+00`, and the multiplication table is `(Z_2)^m` — **identical on block, chain
> and ring, which have three different dual graphs.**
>
> **Geometry enters only through `Z_k` — the electric term and the bath coupling.**

**TWO NON-VACUOUS SURVIVORS.** The **isolation test**: swap the stored incidence for BLOCK while the
operators are RING and fidelity is `1.000` vs RING, `0.667` vs BLOCK, and the mirror image on
reversal — **the answer follows the dynamics and ignores the stored structure entirely.** And the
**`γ` sweep**: fidelity `1.0000` at `γ = 0.5, 0.1, 0.01, 1e-4`, collapsing to `0.5333` at `γ = 0`.

> **REMOVE THE ENVIRONMENT AND NOTHING COMES BACK. THE RECORD SECTOR ON ITS OWN IS GEOMETRY-BLIND.**

### MAP STEP 2 — ALPHA DOES NOT SET THE RESOLUTION

`g²` sets the **validity window of the rate estimator**, not a resolution. It is **not monotone in
`g²` and not monotone in `g²·t*`**: `(g²=1.0, t*=0.2)` fails completely while `(g²=0.5, t*=1.0)` — a
larger product — is exact. Breakdown coincides with the measured coherent leakage rate rising to the
size of the dissipative rates. **The predicted coincidence with the legibility bound does not occur.**

### WHAT THIS DOES TO THE MAP

**Steps 1 and 2 are answered, both negatively.** Geometry is not recoverable from records in any
inferential sense, and the coupling is not a resolution limit. **The three-way relation the map was
built to obtain does not exist by this route.**

**AND IT SHARPENS THE PRINCIPAL'S FRAMING RATHER THAN CONTRADICTING IT.** If boundary formation is
part of gravity's record-level character, **W-51 locates it on the ENVIRONMENT side, not the record
side.** The records carry no geometry at all; what carries it is where the environment couples —
which is also what W-38 measured when the environment's placement selected which region was a region.
**Two independent lanes now say the same thing: geometry lives in the coupling, not in the content.**

---

## W-53 — **THE MIRROR TEST IS VOID, AND THE ASYMMETRY IT EXPOSES IS THE REAL POINT**

`LANE_W53_MIRROR/`, sealed. **Forced-or-not check run BEFORE any dynamics**, on W-51's lesson that a
reconstruction can succeed perfectly and mean nothing.

**VOID.** For a single-link bath `{ℓ}` the protected set is `{S : ℓ ∉ bd(S)} = ker f_ℓ`, and the
family `{ker f_ℓ}` determines every row of the boundary map. Confirmed in code on four carriers:
**the incidence matrix is recovered exactly, by pure set algebra, with no dynamics and no link
labels.** The mirror test is W-51's identity in different clothing. **No dynamics was built.**

### THE ASYMMETRY, WHICH IS NOT VOID

```
  RECORD ALGEBRA (W-51)      identical on block, chain and ring -> GEOMETRY-BLIND
  PROTECTED-SET FAMILY       distinct on all 4 carriers, 6 of 6 pairs -> GEOMETRY-FULL
```

**The same carriers that the record algebra cannot tell apart, the environment-coupling structure
separates completely.** That is W-51's finding stated at its sharpest: **the records carry no
geometry; the coupling carries all of it.**

### WHAT THIS SETTLES ABOUT THE PRINCIPAL'S FRAMING

> *"capacity + boundary formation may well be the signature of gravity at the record level"*

**The first half is where geometry actually lives, and that is now established twice.** But in this
carrier the coupling's geometry is **installed by the registrar**, so what W-38 measured was
**selection among given regions, not formation of them.** Boundary *formation* — the set of possible
boundaries being an output rather than an input — **cannot be tested on a carrier whose incidence is
an input to every quantity we can compute.**

> **EVERY RECONSTRUCTION QUESTION IS CIRCULAR ON THIS CARRIER. That is a structural fact about the
> construction, not a result about records or gravity.**

**So the framing is not refuted and is not testable here.** Testing it needs a carrier in which
adjacency is not given — which is the emergent-geometry problem proper, and is beyond this
construction. **That is the honest boundary of what this program can reach.**

---

## W-54 — **ACCUMULATION ADDS NOTHING. CAPACITY IS ENTIRELY LOCAL.**

`LANE_W54_ACCUMULATION/`, sealed. Road item 20 — the last route this program could decide on its own.
Every prior deflation was obtained at a single scale, so none could have detected a collective effect.

**MADE WELL-POSED FIRST.** On a fixed finite carrier every quantity is an exact function of the
incidence matrix, so emergence in the strong sense cannot arise — that needs a limit. The sharpest
version that **can** fail: **global capacity is `m−1`; how much of it is reachable using only records
supported inside a single part of a partition, and how much requires records that span parts and
belong to none?**

**PREDICTION, STATED IN ADVANCE: deficit `= k−1`, independent of partition shape.**

```
   partition                m    k   global   local-only   deficit   predicted
   9 singletons             9    9        8            8         0          8
   3 rows                   9    3        8            8         0          2
   3 columns                9    3        8            8         0          2
   uneven 5+3+1             9    3        8            8         0          2
   6 singletons             6    6        5            5         0          5
   3 vertical pairs         6    3        5            5         0          2
```

> **THE PREDICTION IS WRONG AND THE RESULT IS A CLEAN NEGATIVE. DEFICIT IS ZERO FOR EVERY PARTITION,
> INCLUDING SINGLETONS.** The full capacity is achievable with **single-plaquette records** — the most
> local objects the carrier has. **No part of a carrier's capacity is collective.**

**Why the prediction was wrong:** the parts are not independent carriers but subsets of one, and
single-plaquette records are already local to singletons. Their independent rank is `m`, and those
avoiding a link on exactly one plaquette give `m−1` directly.

### WHAT THIS CLOSES

**Item 20 is closed, negatively. Accumulation adds nothing on this carrier, and it was the last
question the program could decide by itself.**

Combined with what is already sealed: the decoherence half is `exp(−ΓT)` (W-48); legibility is
locality (C1); the capacity markers are counting (W-50); the carrier's response is Boltzmann (W-49);
geometry is not in the records and every reconstruction question here is circular (W-51, W-53); and
now **capacity has no collective content either.**

> **THE PROGRAM'S HONEST END STATE: a complete, internally verified, three-theorem account of record
> formation in lattice gauge theory, every component of which reduces to known physics, containing no
> distinguishing prediction and no collective content.**

**WHAT REMAINS IS NOT DECIDABLE HERE.** Item 24 — boundary *formation*, the set of possible
boundaries being an output rather than an input — requires a carrier whose adjacency is not installed.
**That is the emergent-geometry problem proper, and no experiment on this construction can reach it.**

---

## W-55 — **ITEM 24: BOUNDARIES ARE NOT DERIVABLE FROM GENERIC DYNAMICS. THEY MUST BE AN INPUT.**

`LANE_W55_FACTORIZATION/`, sealed. Item 24 in its sharpest form: in a bare Hilbert space with no
installed structure, a "boundary" is a **tensor factorisation**. So — can one be an output of the
dynamics rather than an input?

**W-55's optimiser search is VOID and is retained as such.** It plateaued near `0.94` whatever it was
given: at `eps=0` the planted Hamiltonian is exactly local (`1.0000`) and the search reached only
`0.9461`, while at `eps=0.40` it reported `0.936` — **above** the planted truth of `0.875`. The
optimiser, not the physics, set the number.

### THE EXACT QUESTION NEEDS NO OPTIMISER

`H` is exactly local under **some** factorisation iff its spectrum is a **SUMSET** `{a_i + b_j}` —
conjugation preserves the spectrum, so this is a property of the eigenvalues alone.

> **AND IT IS SETTLED BY DIMENSION COUNTING, RIGOROUSLY.** A `d_A × d_B` sumset has `d_A + d_B − 1`
> free parameters while the spectrum has `d_A·d_B` entries. For `4×4`: **7 parameters describing 16
> numbers.** The sumset spectra form a 7-dimensional subvariety of `R^16` — **measure zero. A generic
> Hamiltonian admits NO factorisation making it local.**

```
   planted sumsets           2.3e-16, 7.0e-16, 2.3e-02*     (*one optimiser failure, logged)
   generic random Hermitian  4.4e-02, 4.9e-02, 5.2e-02, 4.6e-02   never near zero, 60 restarts each
   OUR 3x3 patch, H_magnetic 7.7e-17                        EXACTLY a sumset
```

### THE NUMBER THAT ANSWERS ITEM 24

> **OUR OWN CARRIER'S HAMILTONIAN IS EXACTLY A SUMSET — `7.7e-17`. A factorisation exists for it.**
>
> **And it exists because we BUILT it that way:** `H_mag = −Σ_p (W_p + h.c.)` is a sum of **commuting**
> plaquette terms, so its spectrum is a sum of independent contributions by construction.

**THAT IS WHERE THE BOUNDARIES CAME FROM.** Not from the dynamics — from the registrar's choice of a
Hamiltonian assembled out of commuting local pieces. **The moment that choice was made, the set of
boundaries was fixed, and every lane afterwards was reading back an input.**

### ITEM 24, CLOSED

**Boundaries are an output only for Hamiltonians that already carry the structure. For generic
dynamics they are not derivable at all, and must be supplied.** This is a statement about dynamics,
not about our lattice — it is why W-53 found every reconstruction question circular, and it would
have been circular on any carrier built the same way.

**HONEST LIMIT.** One planted control failed at `2.3e-02`, the same order as the generic values, so
the numerics alone separate planted from generic by only ~2×. **The rigorous content is the dimension
count, not the optimisation;** the numbers illustrate it and do not carry it.

---

## SURFACE ASSUMPTION AUDIT — **TWO IMPORTED ASSUMPTIONS. ONE MEASURED TO FAIL, ONE NEVER EXAMINED.**

`SURFACE_ASSUMPTION_AUDIT_V001.md`, sealed. The principal asked whether an imported assumption is the
root cause of not finding a process. **The audit had never been done.** Eight assumptions in the
five-step process, classified.

**A1 — THE HILBERT SPACE FACTORISES. IMPORTED AFTER BEING MEASURED TO FAIL.** The register already
carries it as **IMP-2**: a gauge theory's gauge-invariant algebra does not factorise across a region
boundary. **W-37d broke on this directly.** Yet every Holevo quantity, every partial trace and every
"reduced environment state" in this program presupposes it. **And Gap 2 is the same assumption seen
from the other side** — W-55 asked where factorisations come from and answered: from the modeller.

**A2 — A BACKGROUND TIME PARAMETER. IMPORTED, LOAD-BEARING EVERYWHERE, NEVER EXAMINED ONCE.**
**125 of 125 lanes evolve in `t`. Zero derive it.** A project-wide search for relational or emergent
time returns nothing.

> **GRAVITY IS THE DYNAMICS OF SPACE-TIME STRUCTURE. WITH `t` A FIXED EXTERNAL PARAMETER THERE IS
> NOWHERE FOR GRAVITY TO ACT — AND NONE OF THE THREE FAILED ROUTES TOUCHED IT.** W-31 varied the
> graph; W-44 added lengths, potentials and a clock **read against `t`**; W-45/W-50 varied nothing
> dynamical. **All three left `t` alone.**

**And it matches the paper's own statement of where gravity enters — metric/PROPER-TIME action. We
have no proper time.** One global parameter shared by every part of the carrier, which is precisely
the structure general relativity does not have. **W-44's clock is the sharpest evidence: a clock that
reads the background parameter cannot register a proper-time difference, because there is only one.**

**A8 (Markovian bath) is CLEARED** — Phase B/B1 showed T1 needs no Markov assumption.

**STATUS.** That A2 is imported and unexamined is **established**. That removing it would produce a
process is a **hypothesis, not a result** — and it is falsifiable in the useful direction: build a
carrier with no background time, where evolution is relational, and see whether anything
gravity-shaped has somewhere to act. **If it still does not, A2 was not the obstruction.**

---

## SURFACE AUDIT / A1 — **THE TENSOR PRODUCT WAS THE WRONG FRAMEWORK, AND GAP 2 MAY BE AN ARTIFACT OF IT**

The principal: *"Why would you assume this at the record level?"* — of the step *"to say anything about
a system being monitored by its environment, you first have to cut that space into parts."*

**The tensor product is a POSTULATE ABOUT COMPOSITE SYSTEMS.** It says: given two things, each with a
state space, the pair gets `H_A ⊗ H_B`. **It presupposes that the parts already exist.** Importing it
where parts are the thing to be derived assumes the answer.

### THE FRAMEWORK THAT DOES NOT PRESUPPOSE PARTS

**Algebraic quantum theory.** A subsystem is a **SUBALGEBRA** of observables; its environment is the
**COMMUTANT**. No factorisation is required, and relative entropy is defined for these objects
without one.

> **AND IN THAT FRAMEWORK A REGION'S ALGEBRA IS NOT A FREE CHOICE. It is determined by the field
> content and the constraint: the gauge-invariant operators supported there. THE GAUSS LAW FIXES
> THEM.**

**THIS PROGRAM ALREADY BUILT THAT AND THEN WALKED AWAY FROM IT.** The physical sector was constructed
**from the constraint** — that is algebraic. Then the registrar laid tensor products on top to compute
partial traces, Holevo quantities and redundancy, **while the register carried IMP-2 stating the
algebra does not factorise.**

### WHAT THIS DOES TO GAP 2

> **"Gap 2: supplied, not derived" MAY BE AN ARTIFACT OF THE CHOSEN FRAMEWORK RATHER THAN A FACT
> ABOUT RECORDS.**

**W-55 asked which TENSOR FACTORISATION makes `H` local and answered: none, generically.** It never
asked **which SUBALGEBRAS the constraint singles out.** Those are different questions, and only the
second is the one a gauge theory actually poses.

### WHAT IT COSTS

**Every redundancy result obtained through a partial trace — W-36, W-40, W-43b — was computed in a
framework the system does not have.** The numbers stand as computed; **their interpretation assumed a
split that is not there.** W-37d/W-37e already hit this from the other side and the fix was to abandon
partial traces for classical mutual information — **that fix was applied locally and never propagated.**

**RE-POSED GAP 2:** *given the constraint algebra, is the decomposition into region subalgebras
determined?* **Untested. It is a different question from W-55's and it is the one this carrier can
actually be asked.**

---

## W-56 — **GIVEN A REGION, ITS BOUNDARY DATA IS DERIVED. THE ALGEBRAIC FRAMEWORK NEEDS NO FACTORISATION.**

`LANE_W56_SUBALGEBRA/`, sealed. Road item 25, in the framework that does not presuppose parts.
**NO PARTIAL TRACE IS TAKEN ANYWHERE IN THIS LANE.**

The physical space is `2^4` and every operator in play is a Pauli string on the four effective
qubits, so the question is exact `GF(2)` linear algebra: a plaquette shift is `X_p`; an electric `Z`
on link `ℓ` is the product of `Z_i` over the cycles containing `ℓ`. For a link set `S`, `A(S)` is
generated by `{Z_ℓ : ℓ ∈ S}` and `{W_p : p entirely inside S}`, and its **CENTRE** is
`V ∩ V^⊥` under the symplectic form — **determined by `A(S)` alone, chosen by no one.**

```
   region              inside      dim A(S)   centre generators   |centre|
   one plaquette p0    [0]              4     Z1, Z2                    4
   two adjacent p0,p1  [0,1]            6     Z3, Z2                    4
   two diagonal p0,p3  [0,3]            6     Z1, Z2                    4
   three p0,p1,p2      [0,1,2]          7     Z3                        2
   all four            [0,1,2,3]        8     (identity only)           1
```

> **IN EVERY CASE THE CENTRE IS GENERATED BY THE ELECTRIC CHARGES OF THE CYCLES OUTSIDE THE REGION,
> AND NOTHING ELSE.** Checked explicitly: the centre lives only on the outside cycles for every
> region tested. **It shrinks monotonically as the region grows — `2, 2, 2, 1, 0` generators — and the
> whole algebra is a FACTOR, with trivial centre.**

**This is the Casini–Huerta–Rosabal / Donnelly–Wall structure computed on this program's own carrier
rather than cited.** The registrar's earlier guesses were wrong in both directions: the **magnetic
loop is never central**, and the naive "electric flux summed over boundary links" is central only for
some regions. **The centre was reported rather than guessed, and it is the complement's data.**

### WHAT THIS DOES TO GAP 2 — IT NARROWS IT, IT DOES NOT CLOSE IT

**DERIVED, and this is genuine:** the framework needs **no tensor factorisation at all**; a subsystem
is a subalgebra and its environment the commutant. **Given a region, its boundary data — the centre —
follows from the algebra with nothing chosen.** W-55's negative result was about *tensor
factorisations*, and it does not apply here.

**STILL SUPPLIED:** which links constitute the region. Gap 2 moves from *"which of a continuous family
of factorisations"* to *"which set of links"* — **a far weaker choice, and in a field theory naming a
region of space is a natural act rather than an arbitrary one. But it is still an input.**

**AND THE COST STANDS.** W-36, W-40 and W-43b computed redundancy through partial traces, which this
lane shows are unnecessary and which IMP-2 says are unavailable. **Their numbers were obtained in a
structure the system does not have; the algebraic versions have not been computed.**

---

## W-57 — **REDUNDANCY WITHOUT PARTIAL TRACES: THE RECORD HAS NO PARTIAL LEGIBILITY AT ALL**

`LANE_W57_ALGEBRAIC_REDUNDANCY/`, sealed. Paying the debt from W-56: W-36, W-40 and W-43b computed
"what a fragment knows" through **partial traces**, which presuppose a factorisation IMP-2 says the
system does not have.

**THE ALGEBRAIC QUANTITY NEEDS NO FACTORISATION.** For a subalgebra `A`,
`d_A = max_{a ∈ A, ||a||≤1} |Tr(a ρ_+) − Tr(a ρ_−)| = ||P_A(ρ_+ − ρ_−)||_1`, with `P_A` the
Hilbert-Schmidt conditional expectation. **Only the algebra enters — no factor, no trace-out, no
environment Hilbert space.**

```
   region S                    |S|   dim A      d_A
   one cut link                  1       2   0.000000
   all four cut links            4       8   0.000000
   one perimeter link            1       2   0.000000
   half the perimeter            4      16   0.000000
   PERIMETER MINUS ONE LINK      7      16   0.000000
   the whole perimeter           8      32   1.000000
   everything                   12     256   1.000000
```

> **AN EXACT STEP FUNCTION. Every proper sub-region resolves the record's sectors at `0.000000` —
> INCLUDING seven of the eight perimeter links — and the whole loop resolves them at `1.000000`.
> THERE IS NO PARTIAL LEGIBILITY.**

**This differs from W-36's graded curve** (`0.044, 0.100, 0.187, 0.351, 0.658, 0.929`), and the two
are asking different questions: W-36 measured how much an explicit environment had *learned* through
a coupling; W-57 measures what a region's own algebra *can resolve at all*. **The algebraic answer is
exact where the other was approximate.**

**HONEST LIMIT ON THE STRENGTH.** The sectors are defined by `R` and the states are maximally mixed
within each, so **only `R` itself distinguishes them** — and `R ∈ A(S)` iff its support lies in `S`.
**For this state the result is close to tautological**, and it is the algebraic form of T2 rather than
a new fact. **Its real content is the EXACTNESS: a topological record admits no partial legibility,
not even at 7 links of 8.** With a state carrying correlations, other operators could contribute; that
is untested.

**THE FIRST VERSION FAILED ITS OWN CONTROL** and the definition it exposed is worth keeping: `A(S)`
must be generated by the electric operators on `S` **and by every gauge-invariant Wilson loop whose
SUPPORT lies inside `S`** — not merely by the whole plaquettes inside `S`. Under the narrower
definition `A(PERIM)` does not contain the rim loop at all.

---

## W-58 — **RELATIONAL TIME: THE CONSTRUCTION IS BUILDABLE, THE TEST OF A2 IS NOT SETTLED**

`LANE_W58_RELATIONAL_TIME/`, sealed. Surface audit **A2**: 125 of 125 lanes evolve in an external
parameter and zero derive it. **This is the first lane with no background time anywhere** — a
constraint `H_tot|Ψ⟩ = 0`, time as a correlation between a clock and the rest.

**WHAT IS ESTABLISHED.** The construction works: the constraint is satisfied to `2.7e-16`, and
conditioning on the clock reproduces unitary evolution to `2.6e-15` **from a stationary state with no
external time**. So relational time is realisable on this program's machinery.

**WHAT IS NOT ESTABLISHED, AND THE LANE SAYS SO.** The gravity-relevant case is a **clock whose rate
depends on what is there** — an interaction inside the constraint, making time non-universal. Every
`λ > 0` run is **VOID**: the kernel collapses to dimension 1 and the conditional state is **static**
(`largest variation 0.000`), so there is no history to condition on.

**THREE GATES WERE NEEDED AND THE FIRST TWO WERE WRONG.** The sign (`H_C` generates `+t`, so the
history satisfies `H_C − H_S`); the spectrum (negative `H_S` against non-negative `H_C` leaves only
`n=0`, a static eigenstate). **And the third gate still does not discriminate: participation over
clock readings is `12 of 12` for a single energy eigenstate too, because `|⟨t_k|n⟩|² = 1/d` uniformly
by Fourier duality.** The correct gate is whether the conditional **state** varies with `k`.

> **HONEST STATUS OF A2: UNRESOLVED, WITH A NAMED CAUSE.** A history state needs the clock and system
> spectra to match across many levels, and a generic interaction destroys that matching in finite
> dimensions — the kernel drops from `4` to `1`. **Whether that is a real obstruction or an artifact
> of a `12 × 4` toy is exactly what this construction is too small to say.**

**A2 THEREFORE REMAINS THE ONE UNEXAMINED SURFACE ASSUMPTION.** It has now been *approached* rather
than tested: the timeless framework is buildable, and the state-dependent-rate case — the only one
where gravity could enter — needs a carrier large enough to keep a history alive under interaction.

---

## W-59 — **THE OBSTRUCTIONS COHERE. SIX REDUCE TO FOUR, AND THREE OF THE FOUR ARE NOT DYNAMICAL.**

`LANE_W59_OBSTRUCTION_STRUCTURE/`, sealed. The principal: *"maybe these failures are a feature."*
Every obstruction has the form **"X and Y cannot both be had"** — the shape physics takes when the
content is real. If they are content rather than six accidents they should **cohere**. Decidable, so
decided.

**REDUCTION (a) — W-57 FOLLOWS FROM T2.** Across every subset size `|S| = 0…12`, resolution is an
**exact indicator** of "contains the loop's support", with **zero partial values anywhere**. A
gauge-invariant record is a closed loop (T2); the sectors are distinguished only by it; so a region
resolves them iff it holds the whole loop. **W-57 is T2 plus the state choice, not an independent
obstruction.**

**REDUCTION (b) — W-55 AND W-58 ARE ONE FACT.**

```
   lambda   constraint kernel dim   sumset residual of H_tot
     0.00                       4                  3.022e-16
     0.05                       1                  5.204e-02
     0.15                       1                  4.568e-02
     0.40                       1                  4.297e-02
```

**The kernel collapse and the sumset failure occur together and discontinuously, at the same
coupling.** W-55's "no factorisation" and W-58's "no global history" are both
**SPECTRAL COINCIDENCE IS NON-GENERIC.** One fact, met twice.

### THE STRUCTURE THAT REMAINS

```
   1  T1   a conserved quantity is conserved              DYNAMICAL -- and near-tautological
   2  T2   gauge invariance forces closed loops           KINEMATIC      => implies W-57
   3  T3   capacity = m - 1                               COMBINATORIAL  => implies eviction
   4       spectral coincidence is non-generic            MEASURE-THEORETIC => W-55 and W-58
```

> **SIX OBSTRUCTIONS REDUCE TO FOUR, AND THREE OF THE FOUR ARE NOT DYNAMICAL AT ALL.** Only T1 is,
> and it is close to a tautology about superselection.

**THE HONEST READING, WHICH IS NEITHER THE DEFLATION NOR THE ELEVATION.** The obstructions **do**
cohere — the principal's suspicion is confirmed — but **what they cohere into is structure, not
dynamics.** Record formation in this account is fixed almost entirely by gauge kinematics, counting,
and genericity. **And that is why gravity never entered: gravity would have to be dynamical, and
there is almost no dynamical content in the picture at all.**

**THE ONE SUBSTANTIVE MEMBER IS THE FOURTH.** *The structures records need — a factorisation, a global
history — are NON-GENERIC.* They exist only for special Hamiltonians: sums of commuting local terms,
spectra that match across many levels. **So "a world that has records" is a statement that its
Hamiltonian is special, and the content of the whole account is in WHICH way it is special.** That is
a sharper question than any this program has yet asked, and nothing here answers it.

---

## W-60 — **WHICH HAMILTONIANS ADMIT RECORDS: EXACTLY THE DEGENERATE ONES. AND THAT MEANS EXACT SYMMETRY.**

`LANE_W60_WHICH_HAMILTONIANS/`, sealed. W-59 left the one substantive obstruction — the structures
records need are non-generic — and the content is in **which** way a Hamiltonian must be special.

**THE CRITERION, DERIVED FROM T1 RATHER THAN GUESSED.** A record needs `[H,R] = 0` and must be
non-trivial — not merely a function of the energy, or it carries nothing beyond which level the
system is in. If every eigenvalue of `H` is simple, everything commuting with `H` **is** a function
of `H`. Therefore:

> **A RECORD EXISTS ⟺ THE COMMUTANT OF `H` IS NON-ABELIAN ⟺ `H` HAS A DEGENERATE EIGENVALUE.**

```
                                   multiplicities        dim comm  (predicted)   records beyond f(H)
  GENERIC random Hermitian      [1,1,1,1,1,1,1,1]           8      (8)                    0
  planted [2,2,2,2]                     [2,2,2,2]          16     (16)                   12
  planted [4,4]                             [4,4]          32     (32)                   30
  our 3x3 patch, H_magnetic           [1,4,6,4,1]          70     (70)                   65
  our patch + 0.05 * H_electric  [1,1,2,1,1,3,1,...]        26     (26)                   14
```

**`dim commutant = Σ (multiplicity)²` exactly in every case.**

**RECORDS ARE NON-GENERIC, AND THE REASON IS SPECTRAL: `0 of 400` random Hermitian matrices had ANY
degenerate eigenvalue.**

### THE SYMMETRY LINK, AND IT IS BRUTAL

By Wigner, spectral degeneracy in a physical Hamiltonian comes from **symmetry**. Break it:

```
     eps        multiplicities        dim commutant   records beyond f(H)
   0e+00              [2,2,2,2]                  16                   12
   1e-06      [1,1,1,1,1,1,1,1]                   8                    0
   1e-03      [1,1,1,1,1,1,1,1]                   8                    0
```

> **A SYMMETRY-BREAKING PERTURBATION OF `1e-06` DESTROYS EVERY RECORD. Not gradually — immediately.**

**SO RECORDS REQUIRE EXACT SYMMETRY, NOT APPROXIMATE SYMMETRY.** An approximate symmetry gives
*none*. This is why real records exist at all: the symmetries that carry them — gauge invariance,
particle identity — are **exact**, and exact symmetries are not spoiled by small perturbations. **It
is also why this program's carrier has records: its gauge symmetry is exact.**

### AND A NEW, SHARP ROLE FOR ALPHA

**The electric term breaks the degeneracy.** `H_magnetic` alone carries **65** records; adding *any*
electric term drops it to **14** — and `g² = 0.05` and `g² = 0.5` give the identical `26`-dimensional
commutant. **The reduction is immediate, not gradual.**

> **ALPHA'S ROLE, RESTATED: it is what breaks the degeneracy that records live in. At `g² = 0` the
> carrier holds 65; at any `g² > 0` it holds 14.** That is consistent with every earlier lane — the
> exact records survive only at `g² = 0` — and it says why in one line.

**WHAT THIS ANSWERS AND WHAT IT DOES NOT.** It answers W-59's question: the special Hamiltonians are
the degenerate ones, degeneracy means symmetry, and it must be exact. **It does not say why the
world's Hamiltonian has the exact symmetries it has** — which is now the only question left standing
in this program, and is not one the carrier can be asked.

---

## W-61 — **THE SPACE EMITS THE DEGENERACY. THE PRINCIPAL'S HYPOTHESIS IS CONFIRMED.**

`LANE_W61_TOPOLOGICAL/`, sealed. W-60 left a fatal fragility: records need exact degeneracy, and
SYMMETRY-induced degeneracy dies at `1e-06`. The principal: *"maybe gravity at the record level emits
symmetry."* There is a kind of degeneracy that behaves exactly that way, and W-60 supplies the
contrast that makes the test decisive.

**Same theory, same perturbation, two topologies.**

```
  TORUS 2x2 (genus 1)   unperturbed degeneracy 4
        eps       ground splitting     splitting / gap
    0.0e+00              2.665e-15           6.661e-16
    1.0e-06              4.867e-13           1.217e-13
    1.0e-03              1.282e-06           3.208e-07
    1.0e-01              5.925e-03           1.593e-03

  DISK 3x3 (genus 0)    unperturbed degeneracy 1 -- no splitting possible (control)

  SYMMETRY degeneracy (W-60), same perturbation sizes
    eps = 1e-06   ground splitting  2.003e-06     <- LINEAR in the perturbation
    eps = 1e-03   ground splitting  2.002e-03     <- LINEAR
```

> **AT `eps = 1e-06`: SYMMETRY DEGENERACY SPLITS BY `2.0e-06` — LINEARLY, THE SAME ORDER AS THE
> PERTURBATION. TOPOLOGICAL DEGENERACY SPLITS BY `4.9e-13` — FOUR MILLION TIMES SMALLER.**

**The degeneracy records need does NOT have to be supplied by a symmetry. It can be EMITTED BY THE
TOPOLOGY OF THE SPACE**, and when it is, it is protected in a way symmetry-induced degeneracy is not:
the splitting is suppressed to high order because only a perturbation that WRAPS the space can lift
it. **The disk shows the same theory with no such degeneracy at all — the difference is genus, not
dynamics.**

### WHAT THIS DOES TO THE PROGRAM

**W-60's fragility is resolved and its conclusion is corrected.** *"Records require exact symmetry"*
was too narrow. **Records require exact degeneracy, and TOPOLOGY supplies it more robustly than
symmetry does.**

**AND THIS IS THE FIRST PLACE GRAVITY ENTERS THIS PROGRAM NON-TRIVIALLY AND IN ITS OWN SHAPE.** Not a
metric (W-44: bought nothing), not dynamical topology as a variable (W-31: refuted 3/3), not a
functional resource analogy (W-50: forced by counting). **Here the GENUS OF THE SPACE — a purely
geometric fact, with no symmetry imposed and no metric anywhere — determines whether records can
exist at all, and how well they survive.**

**The principal predicted this shape three times before it appeared:** gravity at the record level
will look different while performing the same function. **The function is supplying the exact
degeneracy records live in. The form is topology, not force, not metric, not curvature.**

**AND IT MEETS HARLOW–OOGURI FROM THE OTHER SIDE.** Quantum gravity forbids exact *global* symmetries,
so a gravitating world cannot use symmetry to protect its records. **Topological degeneracy needs no
symmetry — which is precisely what such a world has left.**

---

## W-62 — **RECORD CREATION. ALL THREE TERMS IN ONE OPERATION, AND THE OLDEST OBSTRUCTION IS RESOLVED.**

`LANE_W62_CREATION/`, sealed. **The logical operators are COMPUTED from the structure, not nominated**
— magnetic ones are cycles that are not boundaries; electric ones are link sets with even overlap
with every plaquette that are not vertex stars. Both by `GF(2)` linear algebra. (The first attempt
guessed them and got `||[W,H]|| = 22.6`.)

```
  TORUS 2x2, physical dim 32, ground-space degeneracy 4, gap 4.0000

  1  the computed logicals commute with H       ||[M,H]|| = 0.00e+00   ||[Z,H]|| = 0.00e+00
  2  and ANTICOMMUTE on the ground space        ||{M,Z}|| = 1.26e-15   ||[M,Z]|| = 4.00
     the record label Z has eigenvalues         [-1, -1, +1, +1]
  3  CONTROL: a CONTRACTIBLE loop does not write            ||[C,Z]|| = 3.01e-17
  4  CONTROL: no LOCAL operator writes          worst single-link Z commutator = 1.59e-16
  5  DISK 3x3: NON-CONTRACTIBLE cycles = 0      same theory, NO record
```

### THE OLDEST OBSTRUCTION IS RESOLVED

W-29/W-30 established that **writable and durable are conjugate** — and it stood for the whole
program. **It was a search for a LOCAL writer for a record only a NON-LOCAL operation can set.**

> **THE RECORD `Z` IS UNTOUCHED BY EVERY LOCAL OPERATOR (`1.6e-16`) AND BY EVERY CONTRACTIBLE LOOP
> (`3.0e-17`), AND IS WRITTEN BY THE NON-CONTRACTIBLE TRANSPORT `M`, WHICH ANTICOMMUTES WITH IT
> (`4.00`). WRITABLE AND DURABLE AT THE SAME TIME.**

T1's no-go applies to operations **commuting** with the record. Transport around a non-contractible
cycle does not commute with it — **and is not local, so noise cannot mimic it.**

### THE THREE TERMS, IN ONE OPERATION

**EM** supplies the gauge field and the holonomy that **is** the record — `Z` on a computed
even-overlap set, with clean `±1` labels on the ground space.

**GRAVITY** supplies the possibility of the operation at all. **On the torus there are non-contractible
cycles; on the disk there are ZERO.** Same theory, same Hamiltonian, same everything — **genus 0 has
no record and no way to write one.** The genus is the difference.

**ALPHA** is what costs the record:

```
     g^2    ground degeneracy    record labels     |{M,Z}|
  0.0000                    4    [-1,-1,+1,+1]    1.26e-15
  0.0001                    3    [-1,+1,+1]       6.02e-08
  0.0100                    1    [+1]             3.94e-23
  0.5000                    1    [+1]             9.97e-32
```

**At `g² = 0` the record has a 4-fold home, clean labels, and a working writer. By `g² = 0.01` the
degeneracy is gone, the label collapses to a single value, and the writer no longer acts.**

> **EM CARRIES THE RECORD. THE GENUS MAKES IT WRITABLE AND PROTECTS IT. ALPHA DESTROYS IT.**
> That is the three-way combination the charter asked for, in one operation, with each role measured
> and each control firing.

---

## W-63 / PROOF_V001 — **THE RESULTS ARE DERIVED, AND THE DERIVATION PREDICTS NUMBERS IT WAS NOT FITTED TO**

`LANE_W63_DERIVATION/` + `PROOF_V001.md`, sealed. Four theorems with proofs, then checked against
what the lanes already measured.

**A — record space dimension `= |H_1| = 2^{2g}`.** Ground states are constant on `B_1`-cosets in
`Z_1`. Derived `4` for the `2×2` and `2×3` tori; **measured `4`.** Disk: `H_1 = 0` ⟹ dimension `1`,
no record; measured **zero** non-contractible cycles.

**B — logicals are `H_1` and `H^1`, and a writer always exists.** `Z(c)M(z) = (−1)^{⟨c,z⟩}M(z)Z(c)`,
and by **Poincaré duality the intersection pairing is non-degenerate**, so every record has an
anticommuting partner. Measured `‖{M,Z}‖ = 1.26e-15`, `‖[M,Z]‖ = 4.00`.

**C — nothing supported on a contractible region acts.** Measured `1.59e-16` and `3.01e-17`.

**D — splitting is `O(ε^d)`, `d` the minimal non-contractible cycle. THE PREDICTION THAT COULD HAVE
FAILED.** Derived `d = 2`; **measured slope `2.000, 2.000, 2.000, 2.000` across two decades**, and
`ε² = 1e-12` against W-61's `4.9e-13`. Symmetry degeneracy has `d = 1` and measured linear.

> **THE THREE ROLES, PROVED: EM supplies the field and the holonomy that IS the record. GENUS supplies
> the record space (A), the writer (B iii), and the protection (C, D). ALPHA is a sum of local terms
> and therefore splits the record space at order `d` and destroys it (D).**

**AND IT RESOLVES W-29/W-30 WITH A PROOF RATHER THAN A WORKAROUND.** Writable and durable are
conjugate **for local operations** — Theorems C and D. B(iii) supplies a **non-local** writer, and C
shows noise cannot mimic it.

**ATTRIBUTION, STATED IN THE PROOF ITSELF.** Theorems A–D **are Kitaev's toric code** and are not
this program's discovery. **What is this program's is the chain making them necessary:** W-60 (a
record exists iff `H` is degenerate) and W-61 (symmetry-sourced degeneracy dies at `1e-06`,
topology-sourced survives by `4×10⁶`), which together force the degeneracy to be topological.

---

## P1 — **THE DEFINITION, AND IT IS NOW THE PROGRAM'S ANCHOR**

`LANE_P1_DEFINITION/` + `CORE_FRAMEWORK_V001.md`, sealed. The program had used "record" in at least
three senses. **One definition, five clauses, and the content is carried by the definition so the
theorems come out clean.**

`R` is a **record** for `(H,{L_k})` if it is **(i)** a bit, **(ii)** durable (`[H,R]=0`, `[L_k,R]=0`),
**(iii)** non-trivial — *not constant on some eigenspace of `H`*, i.e. it distinguishes states of the
**same energy** — **(iv)** writable by some admissible operation, and **(v)** protected: **no**
contractible operation writes it.

**IT DISCRIMINATES.** Toric-code record: **all five pass** (`4.899` non-triviality, `11.314`
protection). Break the degeneracy with a local perturbation: **(iii) collapses to `0.000`.** Take `R`
a function of `H`: **(iii) and (iv) both fail.**

**THREE PROPOSITIONS FALL OUT IMMEDIATELY.** (iii) ⟹ `H` is degenerate — **P2's forward direction is
now one line.** (ii)+(iv) ⟹ the writer is built from neither `H` nor the jumps. **(iv)+(v) ⟹ THE
WRITER IS NON-LOCAL — so W-29/W-30's obstruction is a consequence of the definition rather than an
accident of our carrier.**

**AND THE FRAMEWORK NOW CARRIES THE OBSTRUCTION LIST.** Five solvable (O1 the converse — the one that
decides theorem vs conjecture; O2 the symmetry half; O3 exhaustiveness, to be side-stepped rather
than solved; O4 defining "admissible"; O5 an approximate version of (ii)) and four out of reach and
named as such (lattice-not-spacetime topology; no response to matter; the outcome problem; zero
empirical contact). **Order of attack: O4 → O3 → O2 → O1 → O5.**

---

## ERRATUM — **"TOPOLOGY, NOT SYMMETRY" WAS THE WRONG DICHOTOMY. IT IS A HIGHER-FORM SYMMETRY.**

The principal asked two questions: *why do we assume symmetry?* and *with the number of twins at the
record level, do we have a form of symmetry that goes by another name?* **Both land, and the second
renames the result.**

**ON THE ASSUMPTION.** Statement 2 presupposed symmetry as the candidate source of degeneracy —
standard Wigner lore. **The lore holds one way only.** *Symmetry ⟹ degeneracy* is a theorem;
*degeneracy ⟹ symmetry* is not. **W-61 is a counterexample this program produced itself**, and the
registrar imported the converse anyway.

**ON THE NAME — CONFIRMED IN THE LITERATURE.** Topological order **is** the spontaneous breaking of a
**1-FORM SYMMETRY**. In the toric code the Wilson and 't Hooft loops are **exact 1-form symmetries**,
both commute with `H`, and **the ground states transform non-trivially under the operators on
non-contractible loops — which IS the robust degeneracy.**

> **THE RECORD IS THE ORDER PARAMETER OF A SPONTANEOUSLY BROKEN 1-FORM SYMMETRY.**

### WHAT THIS CORRECTS AND WHAT IT IMPROVES

| was | now |
|---|---|
| statement 2: *symmetry cannot supply the degeneracy* | **ORDINARY (0-form) symmetry cannot** |
| statement 3: *topology can* | **a 1-form symmetry can — and topology is what supplies one** |
| statement 4: *records require genus* | **records require a HIGHER-FORM symmetry; on a manifold that means non-trivial topology** |

**AND IT EXPLAINS THE ROBUSTNESS INSTEAD OF ONLY MEASURING IT.** A 0-form symmetry is generated by
operators acting at points, so a local perturbation couples to it directly and splits the multiplet at
first order. **A 1-form symmetry is generated by EXTENDED (loop) operators, which no local operator
can reach — a perturbation must act along an entire loop.** That is precisely Theorems C and D, and
it is now a *symmetry* statement rather than a topological accident.

**THIS IS A BETTER RESULT, NOT A RETREAT.** "Records require topology" becomes **"records require a
symmetry whose generators are extended"** — which says *why* the protection exists.

### AND IT CREATES A TENSION THAT MUST BE FACED

Harlow–Ooguri forbids **exact global symmetries** in quantum gravity, and the modern statement
**includes higher-form global symmetries**. If records require a 1-form symmetry, and gravity forbids
exact global ones, records would be impossible — which is absurd.

**RESOLUTION, AND IT SHARPENS THE CLAIM.** In the toric code the 1-form symmetry is **EMERGENT** — a
symmetry of the low-energy theory, not an exact microscopic global symmetry. **Emergent symmetries
are not what Harlow–Ooguri forbids.** So:

> **RECORDS REQUIRE AN EMERGENT HIGHER-FORM SYMMETRY.** Exact microscopic symmetry is forbidden by
> gravity and is fatally fragile in any case (W-61); an emergent 1-form symmetry is neither.

**STATUS: this is a re-reading of results already sealed, not a new measurement.** Every number stands
— `slope 2.000`, `4.9e-13` against `2.0e-06`, `2^{2g}`. **What changes is what they are called and why
they hold.** The emergence claim (that the 1-form symmetry here is emergent rather than exact
microscopic) is **asserted from the literature and not yet verified on this carrier.**

---

## O6 — **THE 1-FORM SYMMETRY IS EXACT IN OUR CARRIER, NOT EMERGENT. THE HARLOW–OOGURI TENSION IS LIVE.**

`LANE_O6_EMERGENCE/`, sealed. The emergence claim was asserted from the literature and carries the
whole gravity connection, so it was tested. **It failed.**

```
  1-FORM (toric code, generic local gauge-invariant perturbation)
      eps        ||[W, H+epsV]||   slope      ground splitting   slope
    1.0e-04            0.000e+00                    1.613e-08
    1.0e-03            0.000e+00     --             1.613e-06   2.000
    1.0e-02            0.000e+00     --             1.613e-04   2.000

  0-FORM CONTROL
    1.0e-04            8.517e-04                    1.598e-04
    1.0e-03            8.517e-03   1.000            1.598e-03   1.000
    1.0e-02            8.517e-02   1.000            1.593e-02   0.998
```

**THE CONTROL FIRES CORRECTLY:** a 0-form symmetry is broken at slope `1.000` and its degeneracy
splits at slope `1.000` — **no gap**, exactly as an exact-but-fragile symmetry should behave.

> **BUT THE 1-FORM SYMMETRY IS NOT BROKEN AT ALL: `0.000e+00` AT EVERY PERTURBATION STRENGTH.** It is
> **EXACT**, not emergent — while its degeneracy still splits at order `2`.

**WHY, AND IT IS NOT AN ACCIDENT.** In the gauge-invariant sector the only available local operators
are the electric `Z_k` and the plaquettes. **Both commute with the Wilson-loop generator.** A
single-link shift would not, but it is not gauge-invariant and does not act on the physical space.
**Gauge invariance itself forbids every operator that could break the 1-form symmetry.**

### WHAT THIS DOES TO THE CLAIM

**The proposed resolution of the Harlow–Ooguri tension is withdrawn.** It rested on the 1-form
symmetry being emergent. **In this carrier it is exact — and it is exact precisely because we imposed
gauge invariance.** That is the same "we installed it" defect as A1 and A2.

> **THE TENSION IS LIVE, NOT RESOLVED.** Records require a 1-form symmetry (statements 1–4). Quantum
> gravity forbids exact global symmetries **of every form**. Our carrier supplies an **exact** one.
> **Something must give, and this construction cannot say what.**

**THE THREE LIVE POSSIBILITIES, none decidable here.** (a) The 1-form symmetry of a *gauge* theory is
not "global" in the Harlow–Ooguri sense — the likeliest resolution, and a question about their result
rather than ours. (b) It must be **emergent** in any real realization, which needs the gauge structure
itself to emerge from a non-gauge microscopic model — **not our carrier, where gauge invariance is an
input.** (c) The requirement is weaker than statements 1–4 assert.

**O6 IS NOT CLOSED. It is upgraded from "asserted" to "tested and failed", and the gravity connection
now rests on an open question rather than a citation.**

---

## G1 — **THE RECORD COUNT IS SET BY GENUS, NOT AREA. `capacity = area − 1` NEVER COUNTED RECORDS.**

`LANE_G1_RECORD_COUNT/`, sealed. Two capacity numbers were in the register with different scalings —
`T3: area − 1` (linear in area) and `Thm A: 2^{2g}` (exponential in genus). **Only one can be counting
records under the five-clause definition.** Counted properly: on the **ground space**, and **modulo
vertex stars**, since two link sets differing by a star act identically there.

```
   carrier              area   ground degeneracy   INDEPENDENT records (dim H¹)   T3's area−1
   torus 2x2 (g=1)         4                   4                            2              3
   torus 2x3 (g=1)         6                   4                            2              5
```

> **THE RECORD COUNT IS `2` ON BOTH — IT DOES NOT MOVE WHEN THE AREA GOES FROM 4 TO 6. It equals
> `2g`, exactly as Theorem A predicts, giving a record space of `2^{2g} = 4`.**

**AND T3's OBJECTS DO NOT SURVIVE AS RECORDS.** Of `7` and `31` plaquette-boundary loops, exactly
**1** is independent modulo stars in each case — against T3's counts of `3` and `5`.

> **`capacity = area − 1` IS STILL TRUE AS A STATEMENT ABOUT WHICH PLAQUETTE SUBSETS ARE
> SIMULTANEOUSLY PROTECTABLE FROM A BATH. IT IS NOT A RECORD COUNT.** Those objects fail the
> definition: they act trivially on the ground space, so they are not records under clause (iii).

### WHAT THIS SETTLES ABOUT GRAVITY'S ROLE

**The number of records a world can hold is a TOPOLOGICAL invariant.** It is `2g`. **It does not
depend on the area, the volume, the lattice spacing, the number of degrees of freedom, or any local
geometric quantity.** Double the area and nothing changes; change the genus and everything does.

> **AND THIS IS EXACTLY THE PRINCIPAL'S POINT THAT GRAVITY AT THE RECORD LEVEL LOOKS NOTHING LIKE
> CLASSICAL GRAVITY.** Classical gravity is local, metric, and curvature-driven. **What sets the
> record count here is a GLOBAL, NON-METRIC, TOPOLOGICAL property, and no local measurement can see
> it.** Same function — determining what geometry permits — in a form with no classical resemblance
> whatever.

**TWO MEASUREMENT DEFECTS, BOTH THE SAME CLASS AND BOTH LOGGED.** The first pass checked protection on
the **full Hilbert space**, where a short string maps out of the ground space and so cannot flip
anything — making clause (v) vacuously true for every candidate. The second counted rank in the
**link space** rather than in the quotient `Z¹/B¹`, giving `5` and `7` instead of `2`. **Both are
"counted in the wrong space"; the record lives on the ground space, modulo stars, and nowhere else.**

---

## G2 — **GRAVITY DEFINED, AND AN EXACT BRIDGE TO CURVATURE**

`LANE_G2_GRAVITY_DEF/` + `GRAVITY_DEFINED_V001.md`, sealed.

**DEFINITION: gravity at the record level is `( H₁(Σ), ⟨·,·⟩ )`** — the first homology of the carrier
and its intersection form. The count, the space, the writer and the protection all come from that
pair and nothing else. **No metric, no curvature tensor, no force, no local content.**

**AND THE BRIDGE, WHICH HAD NOT BEEN STATED.** For a closed orientable surface `dim H₁ = 2 − χ`, and
Gauss–Bonnet gives `∫K dA = 2πχ`, so

> **number of independent records `= 2 − (1/2π) ∫ K dA`**

```
   torus (chi=0)         2 − chi = 2      measured degeneracy 4 = 2^2
   sphere (chi=2)        2 − chi = 0      measured degeneracy 1 = 2^0
   disk (boundary)       formula N/A      measured degeneracy 1  (H_1 = 0)
```

> **THE RECORD COUNT IS TOPOLOGICAL, NON-LOCAL AND NON-METRIC — AND IT EQUALS AN INTEGRAL OF THE
> GAUSSIAN CURVATURE.** That is how a count no local measurement can see and the classical
> gravitational quantity turn out to be the same fact.

**FIVE FEATURES OF THE EMERGENCE PICTURE:** the count is an integrated curvature (exact); the
protection order `d` is a **length**; it is global and non-metric; it is a spontaneously broken 1-form
symmetry; and **a simply-connected world has no records at all** — the sphere gives degeneracy `1`.

> **THE SHARPEST STATEMENT THE PROGRAM HAS PRODUCED: a world whose space is a sphere cannot remember
> anything. Not remembers poorly — cannot hold one record.**

**AND THE HONEST LIMIT.** **Content does not change `χ`.** In classical gravity matter curves space;
here nothing records do alters the topology. **We have the half of gravity that says what geometry
permits, and not the half that says how content shapes geometry. Until X2 closes, this is a
correspondence and not an emergence.**

---

## X2 — **χ IS A PARAMETER, NOT AN OBSERVABLE. AND UNDER ALLOW/REQUIRE THAT IS NOT A DEFECT.**

`LANE_X2_BACKREACTION/`, sealed.

**THE STRUCTURAL FACT.** `χ = V − E + F` is a property of the **cell complex**. The Hilbert space is
**built from** that complex. Every operator — unitary, dissipative, any — is a map on that space.
**So `χ` has no eigenvalues, no expectation value and no equation of motion.** Verified: the torus
`2×2` and the cube surface have the **same Hilbert dimension (32)** and different `χ` (`0` vs `2`);
the torus `2×2` and `2×3` differ in **size** but share `χ` — **and share the degeneracy (4)**. The
record count follows `χ` and nothing else, and **nothing inside the theory varies `χ`.**

> **X2 IS NOT A MISSING MEASUREMENT. In any framework with a fixed complex, content CANNOT shape
> geometry, because geometry is not the sort of thing the theory contains.**

### THE PRINCIPAL'S CORRECTION — ALLOW / REQUIRE

The registrar wrote: *"we have the half of gravity that says what geometry permits, and not the half
that says how content shapes geometry — so this is a correspondence, not an emergence."*

**The principal: isn't this the allow/require distinction?** — and it is, and the registrar had it
backwards.

> **"MATTER CURVES SPACE" IS THE CLASSICAL FORM OF GRAVITY. Demanding it AT THE RECORD LEVEL is the
> same import flagged four times before.** At the record level there is no metric to curve. **What
> gravity does there is ALLOW: it says what can exist.** And that is not half a job — **it is the
> whole job at that level.**

**AND X2's RESULT IS EVIDENCE FOR THIS, NOT AGAINST IT.** `χ` is not an observable *at this level*.
Requiring it to respond to content is requiring a level-inappropriate quantity to have dynamics.
**ALLOW is the record-level face of gravity; REQUIRE is its macroscopic face.** Asking for both at
once is a category error — the embryo again.

**THE STATEMENT, CORRECTED:**

> **At the record level gravity ALLOWS. `records = 2 − (1/2π)∫K dA` says exactly what may exist, and
> a sphere permits nothing. That is complete as a record-level statement.** Whether allow becomes
> require at larger scales is a **separate question at a separate level**, and it is the paper's own
> Appendix C step 3.

**AND IT KEEPS A FALSIFIER, so this is not a way of never owing backreaction.** If accumulated
allow-statements never yield anything require-like at any scale, then allow is all there is and the
identification with gravity weakens accordingly. **That is decidable, later, and not here.**

**STATUS CHANGE: X2 moves from "out of reach — the gravity claim is only a correspondence" to
"category error — allow is the record-level face; require belongs to a level this construction does
not reach."**

---

## X2 ERRATUM — **ALLOW IS POTENTIAL AT EVERY LEVEL. THERE IS NO ALLOW→REQUIRE TRANSITION.**

The principal: *"allow is always potential. Why would we expect that to be different at the collective
records level?"*

**The registrar had just written that ALLOW is gravity's record-level face and REQUIRE its macroscopic
face. That is wrong, and the correction is structural.**

> **EINSTEIN'S FIELD EQUATIONS ARE CONSTRAINTS.** In the Hamiltonian formulation they are literally
> that — the Hamiltonian and momentum constraints, `H ≈ 0`, conditions on admissible initial data.
> They do not say what happens; they say which `(geometry, matter)` pairs are **permitted**.
> **"Matter curves space" was never a REQUIRE. It is an ALLOW written in differential form.**

**SO THE CORRECTION IS DOUBLE.** There is no level at which allow becomes require, because require
does not exist at any level. **And the falsifier the registrar attached to X2 — "if accumulated allow
never yields require, the identification weakens" — CAN NEVER FIRE, and is therefore worthless.** It
is withdrawn.

### THE QUESTION THAT REPLACES IT, AND IT IS SHARPER

If gravity is **allow** at every level, the live question is no longer *when does allow become
require* but:

> **IS OUR RECORD-LEVEL ALLOW THE SAME CONSTRAINT AS THE GRAVITATIONAL ONE, SEEN AT A DIFFERENT
> SCALE — OR TWO UNRELATED CONSTRAINTS THAT HAPPEN BOTH TO BE PERMISSIVE?**

**And there is a real structural parallel to test, not merely assert.** Both are **Gauss laws** in the
Dirac sense: conditions on data generated by the theory's own gauge symmetry, with the physical
content appearing as a **boundary term**.

| | our record level | classical gravity |
|---|---|---|
| the constraint | Gauss law at each vertex, `∂s = 0` | Hamiltonian + momentum constraints, `H ≈ 0` |
| generated by | the gauge symmetry | diffeomorphisms |
| physical content is a | **boundary term** — the record is a boundary holonomy | **boundary term** — ADM energy is a surface integral, *because of* the constraint |
| permits | `records = 2 − (1/2π)∫K dA` | which `(g, T)` pairs may exist |

**THE NEW FALSIFIER, and unlike the old one it can fire.** If the record-level constraint and the
gravitational constraint have **different structure** — different generators, or physical content that
is *not* a boundary term on one side — then they are two unrelated permissive statements and the
identification is a coincidence of form. **That is checkable.**

**STATUS: X2's reclassification stands (χ is a parameter, not an observable, and demanding backreaction
imports the classical form). What is withdrawn is the allow→require story layered on top of it, and
the unfalsifiable test that came with it.**

---

## THE IMPORT RULE — **A RESCUE IS ALWAYS DOWNSTREAM OF AN IMPORT**

The principal: *"We won't have to rescue measures if we avoid importing classical gravity measures
into the record level."*

The registrar had responded to the level-migration failure by building a **rescue-detector** into
lane X-6. That treats the symptom. **The rescue was only ever needed because an import was there to
be rescued.** Remove the import and the guard is unnecessary.

### THE IMPORT WAS IN X-6'S OWN FALSIFIER, AND THE REGISTRAR WROTE IT

> *"Different generators, or physical content that is not a boundary term on one side, kills the
> identification."*

This makes **classical gravity's constraint algebra the criterion for whether our object is real.**
A mismatch would read as OUR defect. That is a classical measure applied to the record level as a
test. **It is the same error as demanding backreaction (X-2) and as expecting require (X2 erratum),
in its third costume.**

### THE VECTOR IS THE NAME

**Calling `H₁(Σ)` "gravity" is what lets classical properties in.** Every time the name is used, the
classical field's attributes arrive with it and have to be individually refused — backreaction, an
equation of motion, a require-half, a matching constraint algebra. Each refusal has cost a lane.

**The charter asks for the ROLE of gravity in record formation. It does not ask that classical
gravity be present.** A role-assignment is not a claim of identity, and no property of the classical
field transfers through the name.

### THE RULE, STANDING

> **NO CLASSICAL MEASURE MAY BE APPLIED TO THE RECORD LEVEL AS A CRITERION.**
> A classical quantity may be **computed** and **reported as a relationship**. It may never be the
> **test our object must pass.**
> **TELL:** any falsifier of the form *"our object fails to exhibit ⟨classical feature⟩"* is an
> import and is **void**. A falsifier must be stated in the record level's own terms.

**AND THE OBJECTION THIS RAISES, ANSWERED.** If nothing can falsify "our constraint is gravity's",
that identification is unfalsifiable — the sin corrected hours ago. **The resolution is that the
identification is not claimed.** A claim not made needs no falsifier. What is claimed is the role,
and the role has its own falsifiers already registered (a record on a genus-0 carrier; a local
writer).

### EFFECT ON LANE X-6

- **The COMPUTE half stands.** Measuring our own constraint algebra — commutators, closure,
  structure constants vs functions, the boundary operator, `Z₁/B₁` — imports nothing. It is our
  object described in its own terms, and it is the actual work.
- **The IDENTIFICATION half is WITHDRAWN as an import.** "Is it the same constraint as gravity's"
  is not a charter question.
- What gravity's algebra is remains worth **knowing** and will be recorded as a **relationship**,
  carrying no verdict on our object.
- **The rescue-detector is retired as a permanent fixture.** It was scaffolding for an import.

**CLARIFICATION, same day.** The principal: *"I don't have a problem with using the term 'gravity' as
long as we don't use it as an excuse to import a variety of classical metrics that aren't
applicable."* **The word is restored as a permitted label.** The registrar had over-corrected by
purging it; the purge was never the point. **What survives, and is the whole of the gain: Γ is
defined by R1–R3 and by nothing else, no attribute transfers through the name (D-1, D-2), the
identification with `H₁` is a CANDIDATE and not a definition (G-1 withdrawn → G-5), and uniqueness
is now askable (G-6).** Say "gravity"; never let it do work.

---

## G-7 — **THE CAPACITY LAW IS AN INDEX.** `2 − χ` WAS THE WRONG FORMULA, AND THE DISK SHOWED IT

`LANE_G7_INDEX/g7_index.py` · ranks computed over GF(2) by explicit Gaussian elimination; `∂₁∂₂ = 0`
asserted on every carrier.

> ### `N = dim H₀ + dim H₂ − χ`   —   **holds on all five carriers**
> ### `N = 2 − χ`   —   **FAILS on the disk**

| carrier | V | E | F | χ | H₀ | **H₁** | H₂ | `H₀+H₂−χ` | `2−χ` |
|---|---|---|---|---|---|---|---|---|---|
| torus 2×2 | 4 | 8 | 4 | 0 | 1 | **2** | 1 | **2** ✔ | 2 ✔ |
| torus 2×3 | 6 | 12 | 6 | 0 | 1 | **2** | 1 | **2** ✔ | 2 ✔ |
| torus 3×3 | 9 | 18 | 9 | 0 | 1 | **2** | 1 | **2** ✔ | 2 ✔ |
| tetrahedron (sphere) | 4 | 6 | 4 | 2 | 1 | **0** | 1 | **0** ✔ | 0 ✔ |
| **disk 3×3** | 9 | 12 | 4 | 1 | 1 | **0** | **0** | **0** ✔ | **1 ✘** |

### WHAT THIS CORRECTS

**The disk was registered as a SCOPE CAVEAT — "the relation holds for closed orientable surfaces."
It was never a scope limit.** `2 − χ` is the wrong formula; the index form gets the disk right by the
same mechanism it gets the torus right. **`GRAVITY_DEFINED_V001.md` stated the caveat as if the
carrier were the problem. The formula was the problem.**

### WHERE THE `2` CAME FROM

**`2 = dim H₀ + dim H₂`** — one for **connected**, one for **closed and orientable**. It was never a
constant. On the disk `H₂ = 0`, so the `2` was never there, and the surface law fails by exactly
that `1`.

### WHY THIS MATTERS BEYOND THE ARITHMETIC

**`N` is an INDEX** — an alternating sum of dimensions. `N = 2 − (1/2π)∫K dA` is **Gauss–Bonnet**,
and Gauss–Bonnet is the two-dimensional case of **Euler–Poincaré**. Consequences, in order of weight:

1. **The capacity law needs no geometry.** It is stated for an arbitrary chain complex. Genus,
   curvature and area are special-case vocabulary.
2. **Area-independence stops being an observation and becomes structural.** `χ` is unchanged by local
   refinement of the complex, so `N` cannot move under it. **This is what G1 measured** (2 records at
   area 4 and at area 6) — now with a reason, and confirmed here at a third area.
3. **Curvature enters only through its integral**, which is the same statement seen through the
   surface case.

### THIS DOES NOT DEPEND ON G-6

G-7 is established on **our own carriers**, independently of whether `H₁` is unique. **G-6 decides
whether homology is FORCED by R1–R3; G-7 says that where homology is what carries the record, the
count is an index and not a surface formula.** Both branches of G-6 leave G-7 standing.

**AND THE REGISTRAR'S FRAMING OF G-6 WAS WRONG AND IS CORRECTED.** It was written as if
NOT-UNIQUE were "the bigger result". **UNIQUE is at least as large: it makes R1–R3 force homology,
turns statement 4 from "records ADMIT" into "records REQUIRE", and closes O-3.** The two branches
are large in different directions. Ranking them was the two-valued-object tic, not a finding.

---

## G-4 — **THE RECORD-LEVEL CONSTRAINT, IN ITS OWN TERMS.** `LANE_X6_CONSTRAINT_ALGEBRA`

Stated with no external referent, as **D-1** requires.

> ### The Gauss constraints close into **`Z₂^(NV−1)`** — elementary abelian, exponent 2, with **STRUCTURE CONSTANTS.**

| measured | 2×2 torus | 2×3 torus | meaning |
|---|---|---|---|
| `max ‖[G_v,G_w]‖` over **all** pairs | **0** exact (6 pairs) | **0** exact (15 pairs) | **abelian.** Integer zero in both representations, not a tolerance |
| positive control `‖[σˣ,σᶻ]‖_F` | 32.0000000000 vs closed form 32.0000000000 | — | **a reported zero is a measurement, not a broken norm** |
| `max‖G_v² − I‖` | 0 | 0 | every generator an **involution**; group exponent 2 |
| `∏_all v G_v` | **I** | **I** | one global relation — the generators are **not independent** |
| relation subgroup, **exhaustive** over all `2^NV` words | **2 elements** | **2 elements** | that relation is the **only** one; corank exactly 1 |
| GF(2) incidence rank | 3 (`NV`=4) | 5 (`NV`=6) | `NV−1` — the statement that **the graph is connected**. Scale-free |
| group order, explicit enumeration | 8 = 2³ | 32 = 2⁵ | `Z₂^(NV−1)`, by two independent methods |
| **state-to-state variation in the structure coefficients** | **0** | **0** | **STRUCTURE CONSTANTS.** 48 generic + 48 physical + 48 superposition states, plus the operator identity on every ordered pair |

**The composition law was read off FROM each state, not asserted: on every sampled state the matching
word set was exactly `{e_v + e_w + K}`, the predicted coset, with no state-dependent member.**

### THE CONSTRAINT MAKES A REGION'S TOTAL A SURFACE OPERATOR — AND HIDES THE INTERIOR

| | result |
|---|---|
| `‖P_R − ∏_{l∈∂R} σ_l‖` | **0.000000e+00**, every region, both bases — the bulk constraints cancel **exactly** |
| `support(P_R)`, `\|R\| = 1,4,9,16,25` | 4, 8, 12, 16, 20 — **= 4√\|R\| exactly**, residuals `[0,0,0,0,0]`, fitted exponent **0.500000**: a **PERIMETER law** |
| uncancelled product | `2k² + 2k` — area **plus** perimeter. **The cancellation removes an area's worth of operator content** |
| placement dependence | **min = max in all 14 rows** — the W-46 corner defect is guarded, not merely avoided |
| **hidden interior information** | **`(k−1)²` bits — an AREA — and NONE of it is visible in the boundary operator** |

**So the constraint puts the TOTAL at the boundary while the CONTENT stays an area inside it, invisible
from outside. Boundary support `4k` against hidden interior `(k−1)²`: perimeter against area, diverging
as `k/4`.**

### THE RECORD SPACE AS A QUOTIENT — EXACTLY

| lattice | `dim Z₁` | `dim B₁` | `dim H₁` | |
|---|---|---|---|---|
| torus 2×2 | 5 | 3 | **2** | |
| torus 2×3 | 7 | 5 | **2** | |
| torus 3×3 | 10 | 8 | **2** | |
| | **moves with area** (`F+1`) | **moves with area** (`F−1`) | **CONSTANT** | |

Gauss law ⟺ `ker ∂₁`, exhaustive **match rate 1.0**. Cosets `= 2^{dim H₁}`; every coset has exactly
`2^{dim B₁}` elements. **The record space IS `Z₁/B₁` — the Gauss-law sector modulo plaquettes.**
**Both the numerator and the denominator grow with area; only the quotient does not.**

---

## AND THE WITHDRAWN X-6 QUESTION — **IT WOULD HAVE FAILED. RECORDED SO THE WITHDRAWAL IS NOT AN ESCAPE.**

X-6 was withdrawn as an import (**D-1**) *before* these results returned. Both adversarial roles then
returned the same verdict, and the register states it rather than burying it:

- **prosecution: FIRES.** Ours is a finite abelian *group* with structure constants; gravity's is a Lie
  **algebroid** with structure **functions** (`h^{ab}` in `{H⊥,H⊥}`), and by **Hojman–Kuchař–Teitelboim**
  the hypersurface-deformation algebra plus locality plus `h_ij` **uniquely determines** the ADM
  constraints. Our side has **no analogue of `H⊥` at all**.
- **rescue-detector: FIRES**, and it pre-named the three escapes — the Chern–Simons presentation swap,
  level migration, and the continuum deferral.

**AND THE SPECIFICITY AUDIT KILLED HALF MY ORIGINAL PARALLEL INDEPENDENTLY.** "Physical content appears
as a boundary term" is **generic**: by **Noether's second theorem** *any* local gauge symmetry gives
`Q = ∮k`, off-shell and identically. **Likelihood ratio ≈ 1.** So even had X-6 run as stated, one of its
two legs carried **no evidential weight whatsoever** — I had built a test half of which could not
discriminate anything.

### RECORDED AS A RELATIONSHIP, CARRYING NO VERDICT (permitted by D-1)

| | 3+1 gravity | 2+1 gravity | ours |
|---|---|---|---|
| structure | **functions** (`h^{ab}`) | **constants** in Chern–Simons form | **constants** (measured 0 variation) |
| Lie algebra? | **no** — algebroid | yes, in CS form | yes — abelian group |
| local propagating dof | 2 | **0** | 0 |
| physical space | `C̄/gauge` | `Hom(π₁(Σ),G)/G` | `Hom(π₁(Σ),Z₂) = Z₁/B₁` |

**The last row is the same functor at a different group, and it is logged as an observation with zero
evidential weight for identity. It is not a claim, and nothing in the program rests on it.**

---

## G-10 / G-11 — **THE INDEX SURVIVES OFF MANIFOLDS. THE PROTECTION DOES NOT.** `LANE_G10_NONMANIFOLD`

Eleven carriers: manifolds, a non-orientable surface, a surface with boundary, a pure graph, three
pinch-point complexes, two tori wedged at a point, and **an abstract chain complex with no geometric
realisation at all.** All six self-checks PASS.

| | result |
|---|---|
| **(A) index identity** `dim H₁ = dim H₀ + dim H₂ − χ` | **HOLDS on all 11**, including the abstract complex |
| **(B) surface law** `dim H₁ = 2 − χ` | **FAILS on 6 of 11** — disk, theta graph, all three bouquets, two wedged tori, abstract complex |
| **(C) capacity** — Z₂ gauge degeneracy vs `2^{dim H₁}` | **MATCH on all four tested, INCLUDING the non-manifolds.** Bouquet: 8 = 8. Theta: 4 = 4. Two tori wedged: 16 = 16 |
| **(D) protection** — does `d` grow? | **NO. It does not track the index.** |

### (D) IS THE RESULT, AND IT SEPARATES R1 FROM R3

| family | size | `dim H₁` | **distance `d`** |
|---|---|---|---|
| torus L×L (manifold) | 2, 3, 4 | 2, 2, 2 | **2, 3, 4 — grows** |
| **bouquet of k triangles (pinch)** | 2, 3, 4 | **2, 3, 4 — grows** | **3, 3, 3 — CONSTANT** |

> **The bouquet has UNBOUNDED CAPACITY WITH FIXED PROTECTION.** As many records as you like, every one
> of them destroyable by a **weight-3** operator. The torus is the mirror image: capacity fixed at 2,
> protection growing without limit.

### WHAT THIS DOES TO G-7 — **IT SCOPES IT, AND THE SCOPE IS THE POINT**

**"The capacity law needs no geometry" is TRUE.** Verified on a complex with no geometric realisation.

**"The record needs no geometry" is FALSE.** A record requires **R1 AND R3**. The index delivers R1 on
any chain complex whatsoever. **It says NOTHING about R3.**

> **Γ IS NOT AN INDEX. Γ is an index PLUS a condition that the nontrivial classes have MINIMUM WEIGHT
> GROWING WITH THE CARRIER.** The first half is pure algebra and geometry-free. **The second half is
> where something metric re-enters — and it is the half that does the protecting.**

**The registrar's own next-step framing was therefore half right and half wrong.** "If the law survives
on a non-manifold, needs-no-geometry is a claim" — it survives, and the claim is real **for capacity
only**. The test as posed would have returned a clean PASS and concealed the separation. **It was the
distance column, which the framing did not ask for, that carried the finding.**

### AND IT BEARS ON G-6

**R1 is cheap** — satisfied by graphs, pinch points, and abstract complexes alike. **If uniqueness is
to hold, R3 must be what forces the structure.** G-6 should be read against R3, not R1.

---

## G-6 — **UNIQUENESS IS BROKEN. `H₁(Σ)` IS ONE STRUCTURE SATISFYING R1–R3, NOT THE ONLY ONE.**

Five lenses, each adversarially adjudicated against the five clauses. **One lens (fracton) died on an
API error and did not return — that region of the space is UNSWEPT and is recorded as a coverage gap,
not as a null.**

| survivor | strength | why it is not a surface |
|---|---|---|
| **`H_k` of an ARBITRARY `F₂` chain complex** | **DECISIVE** | the general CSS template; contains members no surface realises |
| **asymptotically good qLDPC** — Panteleev–Kalachev, Leverrier–Zémor, Dinur–Hsieh–Lin–Vidick | **STRONG** | `[[n, Θ(n), Θ(n)]]`. **Not a manifold ANYWHERE**: every vertex link is a `Δ×Δ` grid of `Δ²` squares, and a 2-manifold vertex link is a single cycle. Independently, **Delfosse's bound** `kd² ≤ C(log k)²n` holds for *any* tiling of *any* compact 2-manifold; here `kd² = Θ(n³)` — **exceeded by `Θ(n²/log²n)`** |
| **3D toric code**, `H₁(T³)` | **STRONG** | `[[3L³, 3, L]]`, two distances `L` and `L²`; geometric locality unchanged |
| hypergraph product (Tillich–Zémor) | weak | the toric code *is* the hypergraph product of two repetition codes — closest to the anchor in disguise |

### WHAT DIES AND WHAT SURVIVES

- **DEAD: "records require a SURFACE."** Statement 4's manifold reading is refuted.
- **SURVIVES: "records require non-trivial homology."** Every survivor is `H_k` of an `F₂` complex.
- **BUT NOT ESTABLISHED EITHER:** R1–R3 nowhere require a **CSS** structure. Non-CSS stabiliser codes
  and **non-abelian topological order (Levin–Wen, Fibonacci — writers are Wilson loops, not `F₂`
  chain-complex objects)** are untested. **The general condition is AT MOST this narrow and possibly
  wider.**

### THE CONDITION, AS FAR AS IT IS ESTABLISHED

> **Γ = a length-2 `F₂` chain complex whose `k`-SYSTOLE *and* `k`-COSYSTOLE both exceed the local scale.**

**R3 NEEDS TWO NUMBERS, NOT ONE.** `d_Z` and `d_X` coincide **only in self-dual cases such as the
square-lattice torus** — which is exactly why this program measured a single `d` and derived `d = 2`.
**The single distance was an artefact of the carrier we chose.**

### ERRATUM TO THEOREM B'S ATTRIBUTION — **THE WRITER IS NOT POINCARÉ DUALITY**

`GRAVITY_DEFINED_V001` credited the writer to *"non-degeneracy of the intersection form (Poincaré
duality)"*. **Over-credited.** The pairing `H_k(C) × H^k(C) → F₂` is **perfect for ANY chain complex
over a field** (universal coefficients; `Ext` vanishes over a field). **R2 needs no manifold, no
orientability and no Poincaré duality — it is free linear algebra.** Theorem B's *result* stands;
its stated *source* was too strong.

### AND IT COMPOSES WITH G-11

**G-11 showed R1 is cheap — graphs and pinch points supply capacity freely. G-6 now shows the
homology class is also not unique.** What discriminates is **R3**, and R3 is the systole/cosystole
condition. **The program's centre of gravity moves from "which topology" to "what makes minimum
weight grow."**

---

## D-4 — **ALLOW/REQUIRE IS THE PROGRAM'S TYPE DISCIPLINE, NOT A FACT ABOUT GRAVITY.** F-2 AND F-3 WERE THE SAME ERROR ONE LEVEL DOWN

The principal, on the registrar's *"nothing in the program says how a system comes to have growing
minimum weight"*: **"Isn't this the whole allow/require framework?"**

**It is, and the registrar had walked straight back into it — an hour after registering D-1.**

> **MINIMUM WEIGHT (systole / cosystole) IS A PROPERTY OF THE CARRIER, SO IT IS AN ALLOW.** It says a
> record **may** exist. **Asking what PRODUCES growing minimum weight is asking what produces a
> PERMISSION — the identical category error as demanding backreaction at X-2.**

### THE DISCIPLINE, GENERALISED FROM GRAVITY TO THE WHOLE PROGRAM

| type | what it is | the legitimate question | the category error |
|---|---|---|---|
| **ALLOW** | a permission — what the world's structure **admits**. `records = index`; `systole > local scale`; the Gauss law; Einstein's constraints | **what does it permit, and what is the permitted set?** | *"what produces it?"* |
| **OCCUPANCY** | which of the permitted possibilities **is actually the case** | **what dynamics puts the system HERE rather than elsewhere in the allowed set?** — legitimate, and it is where the arrow lives | — |

**This is not a claim that nothing in physics is a production rule.** Maxwell has evolution equations;
Lindblad dynamics is a production rule. **The discipline is about TYPE-MATCHING: dynamics produces
OCCUPANCIES, never PERMISSIONS.**

### WHAT THIS DOES TO THE RFP, REGISTERED THE SAME DAY IT WAS WRITTEN

| | as written | ruling |
|---|---|---|
| **F-1** arrow | a before with no record, an after with one | **STANDS** — about states, not permissions |
| **F-2** *"the protection TURNS ON"* | min weight goes from `O(1)` to growing | **WITHDRAWN AS WRITTEN.** It asks what produces the systole. **CORRECTED: given a carrier whose systole permits a record, what makes the system occupy a DEFINITE record instead of a mixture?** |
| **F-3** *"the carrier ARISES"* | the complex is not assumed | **RECLASSIFIED, same class as X-2.** It asks what produces the allow. `T-II.5` inherits the reclassification |
| **F-4** dynamics not an act | follows from `(H, {L_k})` | **STANDS** |
| **F-5** measurement | is a measurement the creation of a record? | **COLLAPSES INTO CORRECTED F-2** |

### AND THE COLLAPSE IS A RESULT, NOT A TIDY-UP

> **ONCE FORMATION IS AN OCCUPANCY QUESTION, THE FORMATION QUESTION AND THE MEASUREMENT QUESTION ARE
> THE SAME QUESTION:** *what selects one element of an allowed, exactly degenerate set?*
> **F-5 was never a separate item.** `T-V.3` — "the natural bridge to the measurement problem" — is
> not a bridge to a neighbouring problem. **It is the problem, under the program's own name for it.**

### THE PREDICTION SURVIVES, AND IT WAS ALREADY THE RIGHT SHAPE

The registered advance prediction — *a Davies bath reaches an equal mixture of all `2^{2g}` record
states, which is no record* — **is an OCCUPANCY statement.** It asks what the dynamics puts in the
allowed space, not what makes the space. **It was correctly typed before the requirements around it
were.** It stands unchanged as the RFP's first test.

---

## F-7 FIRST TEST — **THE PREDICTION HOLDS, AND THE REAL RESULT IS STRONGER THAN THE PREDICTION**

`LANE_F7_OCCUPANCY/f7_davies.py`. Toric code, 2×2 torus, `L=8`, dim 256, ground space **4**.
**Logicals COMPUTED from `Z₁/B₁` and `Z¹/B¹`, not nominated** — `‖[Z̄,H₀]‖ = ‖[X̄,H₀]‖ = 0.000e+00`,
`‖{Z̄,X̄}‖ = 0.000e+00`, `Z̄` splits the ground space 2–2.
**DEFECT CAUGHT AND FIXED BEFORE ANY NUMBER WAS READ:** the first build *nominated* the logicals and
self-check reported `‖[Z̄,H₀]‖ = 64` — the identical W-62 defect. Every number below is from the
computed operators.
**BATH POSITIVE CONTROL (the W-49 trap):** the Davies generator reproduces Boltzmann exactly at
`β = 0.5` and `β = 2.0` — `[0.622459, 0.377541]` and `[0.880797, 0.119203]`.

### 1. THE REGISTERED PREDICTION IS CONFIRMED EXACTLY

| β | `‖L(ρ_Gibbs)‖` | `⟨Z̄⟩` | `⟨Z̄₂⟩` | `⟨X̄⟩` | purity on code | `S/ln4` |
|---|---|---|---|---|---|---|
| 0.5 | 7.557e-16 | **0.000000** | 0.000000 | 0.000000 | **0.250000** | **1.000000** |
| 1.0 | 5.913e-16 | **0.000000** | 0.000000 | 0.000000 | **0.250000** | **1.000000** |
| 2.0 | 2.642e-17 | **0.000000** | 0.000000 | 0.000000 | **0.250000** | **1.000000** |
| 5.0 | 4.086e-22 | **0.000000** | 0.000000 | 0.000000 | **0.250000** | **1.000000** |

**Purity `= 1/4` and entropy `= ln 4` exactly: the maximally mixed state on the code space.
THERMAL RELAXATION CANNOT FORM A RECORD.** And a **definite** record decays —
`d⟨Z̄⟩/dt = −1.341e-03` from `⟨Z̄⟩ = +1`.

### 2. AND THE PREDICTION'S *REASON* WAS INCOMPLETE — LIFTING THE DEGENERACY IS NOT ENOUGH

The prediction blamed degeneracy. **Lift it and the record still does not form.**

| perturbation | `‖[pert, X̄]‖` | splitting `ΔE` at `ε=0.2` | `⟨Z̄⟩` | selects? |
|---|---|---|---|---|
| `Σ X_l` | **0.00e+00** | 1.691e-01 | **0.000000** | **no** |
| `Σ Z_l` | 4.53e+01 | 1.691e-01 | **0.211000** | **YES** |

**`Σ X_l` commutes with the WRITER `X̄`, and `X̄` anticommutes with `Z̄` — so any `X̄`-symmetric state
has `⟨Z̄⟩ = 0` IDENTICALLY, at any splitting.** The `Σ Z_l` row is the positive control: the
measurement **can** see selection, so the zeros are a symmetry and not a broken instrument.

> ### **THE CONDITION FOR FORMATION IS NOT "LIFT THE DEGENERACY". IT IS "BREAK THE SYMMETRY GENERATED BY THE WRITER."**

### 3. THE EXCHANGE RATE — **SELECTION IS BOUGHT WITH PROTECTION, AT PARITY**

| `ε` | `ΔE` | `⟨Z̄⟩` | `⟨Z̄⟩/ΔE` |
|---|---|---|---|
| 0.010 | 4.0006e-04 | 0.000500 | **1.2500** |
| 0.020 | 1.6010e-03 | 0.002001 | **1.2500** |
| 0.050 | 1.0037e-02 | 0.012547 | **1.2500** |
| 0.100 | 4.0594e-02 | 0.050742 | **1.2500** |
| 0.200 | 1.6914e-01 | 0.211000 | 1.2475 |

`⟨Z̄⟩ = ½ tanh(βΔE/2)`, with `⟨Z̄⟩/ΔE = β/4` **constant to four figures over three decades.**

| log–log slope in `ε` | |
|---|---|
| splitting `ΔE` | **2.0262** — Theorem D's `ε^d` with derived `d = 2` |
| selection `⟨Z̄⟩` | **2.0210** |

> ### **THE BIAS IN THE RECORD AND THE LOSS OF ITS PROTECTION ARE THE SAME QUANTITY, WITH THE SAME EXPONENT `d`. THERE IS NO `ε` AT WHICH ONE APPEARS WITHOUT THE OTHER.**

### WHAT IS EXCLUDED, AND THE NEXT STEP

**Energy-based dynamics is excluded as the RFP** — not because degeneracy blocks it, but because
**anything that biases the record by an amount unprotects it by the same amount at the same order
`d`. Protection and formability are one quantity with opposite sign.**

**NEXT STEP, CONCRETE.** The escape must be dynamics that is **not** energy-based, i.e. not
detailed-balance. Three candidates, in order of cost:
1. **a driven / non-equilibrium bath** — no Gibbs steady state, so the parity argument does not apply
2. **measurement-induced selection** — which by **D-4** is the *same question*, not a different one
3. **a NON-LOCAL bath coupling** — and note this is exactly what **P-3** says the writer must be

**Test 3 first: it is the only one the program already has a theorem about.**

---

## F-10 — PREDICTION REGISTERED BEFORE THE RUN

**THE QUESTION, POSED SO IT CANNOT BE ANSWERED BY CONSTRUCTION.** Handing a bath a logical operator
and observing that a record forms proves nothing — it assumes the answer. The non-circular question is
a **necessary condition**:

> **With `H = H₀` exactly degenerate and fully protected, what is the MINIMUM WEIGHT of a bath
> coupling that produces `d⟨Z̄⟩/dt ≠ 0` from the maximally mixed code state?**

**PREDICTION.** By **Knill–Laflamme**, any operator `C` with `wt(C) < d` satisfies `P_g C P_g ∝ P_g`,
so it cannot distinguish record states. Therefore:

> ### **NO COUPLING OF WEIGHT `< d` CAN FORM A RECORD, AND THE FORMATION THRESHOLD EQUALS THE DESTRUCTION THRESHOLD. BOTH ARE `d`.**

At `d = 2`: **weight-1 couplings give exactly zero; weight-2 can be non-zero.**

**LIMITATION, STATED IN ADVANCE.** The 2×2 torus is the only carrier of feasible dimension (256), and
it has `d = 2`. **The threshold is therefore tested at ONE value of `d`. The scaling claim is NOT
tested by this run** — 3×3 would need dim 262144.

---

## F-10 RESULT — **THE FORMATION THRESHOLD IS THE CODE DISTANCE. BOTH ARE `d`.**

`LANE_F7_OCCUPANCY/f10b_threshold.py`.
**v1 WAS MIS-DESIGNED AND THE MEASUREMENT CAUGHT IT:** it measured `d⟨Z̄⟩/dt` under a **Davies** bath,
but **F-9 already excludes every detailed-balance bath at any weight** — so the answer was fixed by
the bath class before weight could matter, and it returned `threshold = None`. **A test whose result
is entailed by a prior result is not a test.** v2 separates the necessary condition (bath-independent)
from sufficiency (which must use a non-equilibrium bath by construction).

### A. NECESSARY — Knill–Laflamme, and it cannot inherit F-9's exclusion because no bath appears in it

| weight | # couplings | `max ‖P_g C P_g − (trC/4)P_g‖` | |
|---|---|---|---|
| **1** | 24 | **6.344e-16** | **acts as a SCALAR — cannot form** |
| **2** | 252 | **2.000e+00** | **distinguishes code states** |
| 3 | 1512 | 2.000e+00 | distinguishes code states |

> **MINIMUM WEIGHT WITH ANY CODE-SPACE ACTION = 2 = `d`. THE PRE-REGISTERED PREDICTION IS CONFIRMED.**

### B. SUFFICIENT — a non-equilibrium bath, and it costs NO protection

Single non-Hermitian jump operator `σ⁻ = P₊X̄P₋` (**support 3 qubits**; `‖σ⁻ − σ⁻†‖ = 16.0`, and **a
Hermitian jump operator can only dephase, never select**). `H` is untouched — still `H₀`, still exactly
degenerate.

| `t` | `⟨Z̄⟩` | code-space weight | purity on code |
|---|---|---|---|
| 0.00 | 0.000000 | **1.000000** | 0.250000 |
| 1.00 | 0.635830 | **1.000000** | 0.351070 |
| 4.00 | 0.982412 | **1.000000** | 0.491283 |
| 6.00 | **0.997667** | **1.000000** | 0.498835 |

**THE RECORD FORMS WITHOUT EVER LEAVING THE CODE SPACE AND WITHOUT LIFTING THE DEGENERACY.** Purity
tends to `1/2`, not 1, because selecting one logical qubit leaves the other mixed — as it must on a
4-dimensional code space.

### C. CONTROL — 576 operators `(A + iB)/2` with `A,B` of weight 1: `max |d⟨Z̄⟩/dt| = 5.551e-16`

**No weight-1 jump operator forms a record, Hermitian or not.**

---

### WHAT THIS ESTABLISHES, AND IT REVERSES F-9's VERDICT UNDER ONE CONDITION

**F-9:** under **energy-based** dynamics, selection is bought with protection at parity — same
exponent `d`, no `ε` where one appears without the other.
**F-10:** **drop detailed balance and selection costs NO protection at all** — `⟨Z̄⟩ → 0.9977` with
code-space weight `1.000000` throughout and `H` untouched. **But the coupling must weigh at least `d`.**

> ### **THE RFP HAS TWO REQUIREMENTS, BOTH NECESSARY:**
> ### **(1) NON-EQUILIBRIUM DYNAMICS — no detailed balance. (2) A COUPLING OF WEIGHT ≥ `d`.**
> ### **AND `d` IS THE SAME `d` THAT PROTECTS. A RECORD IS EXACTLY AS HARD TO FORM AS IT IS TO DESTROY.**

Requirement (2) extends **P-3 / Theorem B** from the writer to the *bath*: **whatever forms a record
must couple at least as non-locally as the record is protected.**

### THE OBSTRUCTION THIS CREATES — the program's problem has MOVED, not closed

**A bath coupling of weight ≥ `d` is itself a non-local process.** The question is no longer *how does
a record form* but **WHAT SUPPLIES A NON-LOCAL, NON-EQUILIBRIUM COUPLING?** That is the next
obstruction and it is sharper than what it replaced.

**LIMITATION, AS STATED IN ADVANCE:** `d = 2` only. **The threshold is tested at ONE value of `d`;
the scaling `threshold = d` is NOT tested by this run.**

---

## A-EM CONSOLIDATED — **EM'S ROLE IS LARGER THAN "SUPPLIES THE CARRIER", AND TODAY'S RESULTS SHOW HOW MUCH**

The principal: *"I think we saw in the past how EM has a powerful role in this process."*
**The register under-credited it.** `A-EM` read *"supplies the gauge field and the holonomy that IS
the record"*; `THE_CLAIM` reads *"EM provides the carrier"*. Four measured results say more.

### 1. EM SUPPLIES BOTH BOUNDARY MAPS — THE WHOLE CHAIN COMPLEX, NOT JUST THE FIELD

| the complex | `C₂ →^{∂₂} C₁ →^{∂₁} C₀` | what it is physically |
|---|---|---|
| `∂₁` | **IS the Gauss law** | EM's constraint |
| `∂₂` | **IS the plaquette / field-strength term** | EM's magnetic action |
| `C₁` | the links | where the gauge field lives |

**The carrier supplies the SETS (which cells exist and how they meet). EM supplies the MAPS.**
`H₁ = Z₁/B₁` needs both — measured at match rate **1.0** (G-9).

### 2. THE RECORD SPACE IS THE HOMOLOGY OF EM'S OWN COMPLEX — G-9, MEASURED

`dim Z₁ = 5,7,10` and `dim B₁ = 3,5,8` both move with area; `dim H₁ = 2,2,2` does not.
**Both the Gauss-law sector and the plaquette space are EM structures. The record is what survives
their quotient.**

### 3. **NEW — THE MINIMAL COUPLING THAT CAN FORM A RECORD IS AN EM HOLONOMY.** `f10c_whatformsit.py`

Of 252 weight-2 Paulis, **8 distinguish code states. All 8 commute with every stabiliser — all 8 are
LOGICAL.** And every one is an EM holonomy on a non-contractible cycle:

| | count |
|---|---|
| all-`Z` supported on a **CYCLE** — magnetic Wilson loop | **4** |
| all-`X` supported on a **COCYCLE** — electric Wilson loop | **4** |
| **neither** | **0** |

### 4. SO THREE DIFFERENT THINGS TURN OUT TO BE THE SAME KIND OF OBJECT

> ### **THE RECORD IS AN EM HOLONOMY ON A NON-CONTRACTIBLE CYCLE (Thm B ii).**
> ### **ITS WRITER IS AN EM HOLONOMY ON A NON-CONTRACTIBLE CYCLE (W-62).**
> ### **AND THE MINIMAL THING THAT CAN FORM IT IS AN EM HOLONOMY ON A NON-CONTRACTIBLE CYCLE — 8 of 8, ZERO EXCEPTIONS.**

**The program did not put these three in the same class. They arrived there separately, by three
different measurements, and the third is new today.**

### WHAT THIS DOES TO Γ

**Γ is not an independent field alongside EM.** It is **the statement that EM's own complex has
non-trivial homology.** The carrier decides which cells exist; EM decides the maps; **Γ is a property
of the pair, and of neither alone.**

**CAUTION, STATED:** this does NOT dissolve Γ into EM. The cell structure is not EM's — G-6 already
showed the homology can be that of complexes with no manifold behind them at all, and G-11 showed
protection needs something metric that homology alone does not supply. **Γ remains a distinct
requirement (R1–R3). What changes is that it is a requirement ON EM'S COMPLEX, not a separate
ingredient beside it.**

---

## D-5 / A-EM4 — **AT THE RECORD LEVEL EM HAS NO LOCAL CONTENT. THE FIELD-VALUE PICTURE DOES NOT TRANSFER.**

The principal: *"at the record level EM is a massive object and not just some field measured in other
contexts."* **`LANE_EM_EXTENT/em_extent.py`** — closed forms verified against explicit GF(2) ranks at
`L = 2,3,4,5` (all PASS).

### HOW BIG IS THE EM OBJECT THAT CARRIES ONE RECORD

| `L` | links | `dim Z₁` | `dim H₁` | **record fraction** | `d` (support) | `d`/links |
|---|---|---|---|---|---|---|
| 2 | 8 | 5 | 2 | 4.000000e-01 | 2 | 0.250000 |
| 4 | 32 | 17 | 2 | 1.176471e-01 | 4 | 0.125000 |
| 16 | 512 | 257 | 2 | 7.782101e-03 | 16 | 0.031250 |
| 256 | 131072 | 65537 | 2 | **3.051711e-05** | **256** | 0.001953 |

**THREE THINGS ARE TRUE AT ONCE, AND NO FAMILIAR OBJECT SATISFIES ALL THREE:**

| | |
|---|---|
| **absolute extent `d = L → ∞`** | the object carrying one bit grows **without bound**. It is not a local excitation |
| **record fraction `2/(L²+1) → 0`** | it is a **vanishing fraction** of EM's own gauge-invariant content. It is not "all of EM" either |
| **local content ZERO** | **no local operator reads it** — Thm C, `1.59e-16` |

### THE IMPORT GUARD

> **D-5. MEASURES THAT APPLY TO EM AS A LOCAL FIELD HAVE NO REFERENT AT THE RECORD LEVEL.**
> Field strength at a point, energy density, local flux, a value `E(x)` or `B(x)` — **none of these
> has a record-level meaning.** Using one is an import of the same class as classical gravity's, and
> **D-1 applies to it unchanged.**
> **TELL:** any sentence that assigns EM a value, density or strength **at a place**.

### AND IT CORRECTS A WEAK FRAMING THE PROGRAM HAS BEEN CARRYING

**`THE_CLAIM_V001` says "EM provides the carrier."** A *carrier* suggests a medium that something is
written **on**. **There is no such separation here: the record IS the holonomy.** The object and its
carrier are the same thing, that thing is extended, and **it has no local parts to be written on.**

**Found by the principal's prompt, not by the running import audit** — recorded so the audit's eventual
coverage is not overstated.

---

# IMPORT AUDIT — **FOUR LOAD-BEARING IMPORTS BESIDES CLASSICAL GRAVITY. THREE ARE HIGH AND ALL FOUR HIT ROWS REGISTERED TODAY.**

Six lenses, adversarially adjudicated, each testing not *is this imported* but **would a registered
PROVED row change if it were dropped.** The principal asked the question; the audit answered it against us.

## 1. **THE ENSEMBLE AVERAGE `⟨R⟩ = Tr(Rρ)` USED AS THE TEST FOR FORMATION** — HIGH, and the program had ALREADY CAUGHT IT

**`SURFACE_ASSUMPTION_AUDIT_V001.md:16`** — *"**A4** | "durable" = expectation value constant | imported;
**W-35 showed ensembles are blind to record formation**"* — **and A4 is NOT CLEARED.**
**`REGISTER:2611`** — *"Every lane from W-28 onward measured `⟨R⟩ = Tr(Rρ)` — an ENSEMBLE AVERAGE.
**A fair coin has mean zero and every single flip is definite.**"* W-35 unravelled the same Lindbladian
and found a **definite record on 78.5% of runs.** `REGISTER:2646` closes the escape: *"Both
unravellings give the same `ρ(t)`"* — so purity `1/4` and `S/ln4 = 1` do not rescue the reading either.

**AND THE ENTIRE F-SERIES WAS RUN THROUGH THAT SAME DOOR TODAY.** `grep` of `LANE_F7_OCCUPANCY/` for
`unravel|trajector|single-run|ensemble` returns **ZERO hits.** The audit naming A4 is timestamped
**11:28**; `f7_davies.py` is **14:37**; `f10c` is **14:56**. **Six PROVED rows registered on top of a
named, uncleared assumption within three hours of it being named.**

**WHY IT WAS INVISIBLE, and this is the structural lesson:**
1. **A4 lives in prose, not in the ledger.** It has no row and no D-entry. **The ledger is what gets checked.**
2. **A blind instrument returns exact zeros.** `⟨Z̄⟩ = 0.000000` and `S/ln4 = 1.000000` read as strong
   confirmation. **They are the signature of the blindness.**

## 2. **"MINIMUM WEIGHT MUST GROW WITH THE CARRIER"** — HIGH. **CLAUSE (v) IS BINARY AND SAYS NO SUCH THING**

Clause (v) reads *"NO contractible operation does."* **Binary.** Nothing in the corpus argues that a
record's protection must **grow**. That requirement is imported from **asymptotic good-code-family
theory**, where `d = Θ(n)` is the quality criterion.

**AND IT INVERTS G-11's CONCLUSION.** The bouquet's cycles ARE the generators of `H₁` — they are
**non-contractible**. So a weight-3 operator traversing one is **not a contractible operation**, and
**the bouquet SATISFIES clause (v) as written.** G-11's *numbers* stand (capacity 2,3,4 against `d`
3,3,3). Its *conclusion* — that protection fails there, so Γ needs a growing-weight condition — is the
import.

**G-12's refuter cannot fire:** `"the local scale"` occurs three times in the corpus and **none of them
is a definition.** That is precisely the defect withdrawn at `GRAVITY_DEFINED_V001`: *"it can never
fire… an unfalsifiable test is worthless."*

## 3. **GIBBS STANDING IN FOR "ENERGY-BASED DYNAMICS" AS A CLASS** — HIGH. **THE VERDICT SURVIVES; THE STATED CONTENT DOES NOT**

Re-run on our own carrier. The registered row reproduces exactly (ratio `1.2500`, slopes `2.0262/2.0210`)
— **and the parity is canonical-only.** The constant is `−f'(0)/f(0)/4`:

| occupation `f(ΔE)` | constant | selection slope | splitting slope |
|---|---|---|---|
| Gibbs `β=5` | 1.2500 | 2.0210 | 2.0262 |
| Gibbs `β=1` | 0.2408 | 2.0187 | 2.0262 |
| power law `1/(1+ΔE)³` | 0.6552 | 1.9795 | 2.0262 |
| **Gaussian `exp(−2ΔE²)`** | 0.0050 | **4.0465** | 2.0262 |

**The exponent equality FAILS for the Gaussian, which is equally a function of `H` alone.**
**And `β = 0` is canonical and gives `⟨Z̄⟩ = −0.000000` (spread `1.97e-13`) at splitting `1.691e-01`
— SPLITTING WITHOUT SELECTION, which the register asserts does not exist.**

**WHAT SURVIVES, with its correct warrant:** for **any** occupation that is a function of `H` alone,
exact degeneracy forces flatness on the code space — verified `⟨Z̄⟩ ≤ 1.1e-14` across Gibbs `β=5` and
`β=0.5`, power-law `1/(1+ΔE)⁷`, a non-monotone `0.1+|sin 3ΔE|`, and a **random** function of the
eigenvalue. **The exclusion holds. "The same quantity, same exponent, no `ε` where one appears without
the other" is WITHDRAWN.**

## 4. **DETAILED BALANCE AS A CLASS IN A NECESSITY CLAIM** — MEDIUM, **AND F-13 MAY CONTRADICT ITSELF**

The evidence base is **one Davies generator built from weight-1 site operators.** `grep` for
`weak coupling|Born-Markov|secular|mean force` across the whole corpus returns **nothing relevant.**

**THE SELF-CONTRADICTION, and it is internal:** a strong-coupling thermal environment relaxes to the
**mean-force Gibbs state** `Tr_B[e^{−βH_tot}]/Z`, which depends on the **coupling operators** and not
on `H_S` alone. **For a coupling of weight ≥ `d` — exactly what F-13 clause (2) REQUIRES — that state
is not degeneracy-blind.** So a bath in detailed balance with respect to `H_tot` **meets F-13's own
registered falsifier while satisfying clause (2).**

---

# THE REMAINING WORK — **SIX ROWS ATTACKED, ALL SIX ADVERSARIALLY VERIFIED. O-1 IS SETTLED.**

## O-1 — **THE CONVERSE IS FALSE.** AND THE REPAIRED CONVERSE IS A THEOREM

**REGISTRY HYGIENE, FLAGGED:** P-2 reads *"(ii)+(iv) ⟹ the writer is built from neither `H` nor
`{L_k}`."* **The O-1 row describes P-1's converse, not P-2's.** Recorded, not silently fixed.

### LEMMA 0 — CLAUSE (ii) IS FAR STRONGER THAN IT READS

For `R` Hermitian, `[L,R] = 0 ⟺ [L†,R] = 0`. **So clause (ii) forces `R` into the commutant of the
`*`-ALGEBRA `A = alg{I, H, L_k, L_k†}` — not of the set `{H, L_k}`.**

### THE CONVERSE FAILS ON BOTH GAPS, WITH CONTROLS

| | |
|---|---|
| `H = diag(0,0,1,1)`, **one generic jump** | Hermitian commutant = scalars only **400/400**; **record exists 0/400** |
| **positive control**, same `H`, **diagonal** jumps | **record found 400/400** |
| `H = 0` (maximal degeneracy, the most favourable case there is), `d = 2,3,4,6,8` | `dim A = d²` in **all 10** rows; commutant dim **1**; record **no** in all 10 |
| **codimension** of `{L : [L,R]=0}` | `dim = p²+q²`, **6/6 exact**. Codim `2pq` — **measure zero for every `d ≥ 2`** |
| instability | `‖[L(ε),R]‖/ε = 4.9456` at every `ε` from `1e-8` to `1e-1` — **first order. 0/200 records.** Clause (ii) cannot be rescued by "small enough noise" |

**THE TRAP THAT KILLS THE OBVIOUS REPAIR.** A case with `max‖P_λ L_k P_λ‖ = 0.000e+00` — the noise
identically zero **inside every eigenspace** — and **still no record**, because `P A P` is irreducible.
Control with the return legs deleted: record **YES**, `(ii) = 0.00e+00`. **Second-order excursions out
of the eigenspace and back destroy records, so inspecting `P L P` alone can NEVER prove the converse.**

### THEOREM O1-A
> **A record non-constant on `E_λ` exists ⟺ the compression `P_λ A P_λ` is a PROPER subalgebra of `B(E_λ)`.** 6/6 against brute force.

### THEOREM O1-B — **AND IT REFUTES THE CONJECTURE IN MY OWN BRIEF**

I wrote that flippability was probably free on any degenerate subspace. **False.**
> **`U` with `U†RU = −R` exists ⟺ `p = q` ⟺ `Tr R = 0`.**

Balanced `(2,1,1),(4,2,2),(6,3,3),(8,4,4)`: `‖U†RU+R‖ = 9.6e-16 … 2.5e-15`. Unbalanced: **no flipper**,
and the positive control — **best of 4000 random unitaries each** — never drops below `2.000116`.

> ### **COROLLARY O1-B1: IN AN ODD-DIMENSIONAL HILBERT SPACE, `Tr R` IS ODD, SO NO RECORD IS EVER WRITABLE.**

### THE REPAIRED CONVERSE — **A THEOREM, 210/210 BY EXHAUSTIVE ENUMERATION**
> **A record satisfying (i)–(iv) exists ⟺ the commutant of `alg{I,H,L_k,L_k†}` contains a projection
> that is (H1) non-trivial on some eigenspace of `H` and (H2) trace-balanced.**

---

## O-4 — **"ADMISSIBLE" DEFINED, AND THE UNDEFINED WORD WAS LOAD-BEARING IN CLAUSE (v), NOT (iv)**

> **ADMISSIBLE `U` ≝ a unitary with `[U,H] = 0`.**

Consequences: **clause (iv) becomes exactly `Tr(P_E R) = 0` on every eigenspace of `H`** — agreeing
with O1-B independently. **P-3 survives but its conclusion weakens** to *"no admissible writer fits
inside ONE contractible region."*

> ### **AND THE STING: under the trivial reading "any unitary", THE TORIC CODE FAILS THE PROGRAM'S OWN CLAUSE (v).** The word was never optional.

## O-2 / O-3

**O-2 PROVED under stated hypotheses.** First-order splitting **is exactly** the non-vanishing of
`Φ(V) = PVP − (tr PVP/n)P`; guaranteed at weight 1 for a **connected** 0-form symmetry because its
charge is a sum of local densities. **And the topological vanishing to order `d−1` IS the
Knill–Laflamme condition** — slope `= d` confirmed on **five codes**.
**O-3's registered side-step is REFUTED by two explicit counterexamples.**

## O-5 — the relaxation must be a **spectral width**, never a norm bound

> **Clause (ii) relaxes as a width bound on the cluster carrying the DRESSED record.** The tolerance is
> then **exactly an inverse lifetime**, `T(η) ≥ η/δ` with `δ = c·p^d/Δ^{d−1}`.
> **W-61's four-million-fold separation SURVIVES. Read the naive way — a norm bound on a fixed `R` — it collapses to 0.35.**

## O-7 / O-8 — both gaps swept

- **FRACTONS WIDEN NOTHING.** X-cube, checkerboard **and Haah's cubic code are all length-2 `F₂` chain complexes.**
- **NON-CSS stabiliser codes provably ESCAPE that class.**
- **NON-ABELIAN order FAILS CLAUSE (iv) OUTRIGHT**, on a parity obstruction an `F₂` complex can never trigger.

## O-10 — **CLOSED, AND UPGRADED FROM A MEASUREMENT TO A PROOF**

**Threshold `= d` on all six codes at `d = 1, 2, 3`, including a non-CSS and a degenerate code — and the
equality is now PROVED for every stabiliser code and every coupling, not only Paulis.**
**F-11/F-13 no longer rest on the single `d=2` carrier.**

---

## O-12 CLOSED — **CLAUSE (v) WAS FALSE OF OUR OWN CARRIER. THE WORD FIXES IT.**

**AS PREVIOUSLY WRITTEN, clause (v) carried no admissibility restriction and was FALSE OF THE TORIC
CODE.** A single-edge `X` on `(0,0,H)` satisfies `X†RX = −R` **exactly**.

| | `L=3` | `L=5` |
|---|---|---|
| single-region Paulis that flip `R` (**positive control**) | **1158** | **2.815e15** |
| **ADMISSIBLE** flippers, `[U,H]=0`, over certified-contractible regions | **0** (36 regions) | **0** (300 regions, to 24 edges) |

Exact by `F₂` linear algebra, not sampling — both sets are subspaces, so the count is `0` or exactly
half. **Independently reproduced by a second implementation** (brute-force symplectic enumeration over
all `4^{|T|}` Paulis, no shared code): `1158` and `0`, min admissible weight **3 = d**.

### P-3's STRONG READING IS NOT UNPROVEN. **IT IS FALSE.**

Clause (v) quantifies over **single** regions, never **products**. A (v′) forbidding products would
give *"nothing local can write it"* — **and (v′) is false by the program's own Theorem D**, which says
`d` local terms **do** reach the record. **The gap cannot be closed by strengthening (v).**

**The minimum-weight admissible writer at `L=3` is weight `3 = d`, and it is a PRODUCT of three
single-edge operators each individually INADMISSIBLE (`‖[X_e,H]‖ = 45.255`). The record is written by
an operation that is admissible only AS A WHOLE.** This is the same `d` as F-13's forming coupling.

### DEFECTS IN THE LANE THAT PRODUCED THIS, CAUGHT BY THE VERIFIER

1. **THE REGISTERED FALLBACK WAS A MISCOUNT.** DEF-A′'s published table is arithmetically impossible
   in 3 of 7 rows — it reports **more** admissible unitaries than *all* unitaries. Corrected, DEF-A′
   equals DEF-A except on one row, and it fails the very case its own docstring says it rescues.
   **DEF-A has NO working escape hatch.** The row I registered was wrong.
2. **5 OF THE 108 SELF-CHECKS WERE `check(..., True)` — LITERAL TAUTOLOGIES**, and one of them passed
   the demonstrably wrong table above. **103 real checks, 5 that cannot fail.**
3. My previous entry said the anchor was amended with the definition. **The patch matched nothing and
   silently wrote the file unchanged.** The clauses are in a table; I matched against a plain-text
   form that does not exist in the file. **Amended now, and verified by re-reading.**

> ### **D-8. A CHECK THAT CANNOT FAIL IS NOT A CHECK.** `check(..., True)`, an assertion with no
> ### computed left-hand side, or a control that returns the same value under every input, is
> ### decoration. Here five of them let an arithmetically impossible table through.
> ### **TELL: a self-check whose expected value is a literal.**

---

## O-11 — PREDICTION REGISTERED BEFORE THE RUN

**THE CHARGE.** F-13 clause (1) says the RFP needs dynamics **outside** the detailed-balance class. But
a **strong-coupling** thermal environment does not relax to `e^{−βH_S}`; it relaxes to the **mean-force
Gibbs state** `ρ_MF = Tr_B[e^{−βH_tot}]/Z`, which depends on the **coupling operators**. If `ρ_MF` is
non-flat on the code space for a coupling of weight ≥ `d` — exactly what clause (2) **requires** — then
a **genuinely thermal environment forms a record**, and F-13 contradicts itself.

**THE SYMMETRY THAT DECIDES IT.** `H_tot = H_S + H_B + λ A⊗B` commutes with `X̄` iff `[A,X̄] = 0`, and
`{Z̄,X̄} = 0` forces `⟨Z̄⟩ = 0` whenever it does. **So the coupling must break the writer's symmetry —
which is F-8, reappearing at the level of the total Hamiltonian.**

**THE CLOSED FORM, DERIVED IN ADVANCE.** Take `A = Z̄` (weight `2 = d`), with `H_B` and `B` diagonal.
`Z̄` commutes with `H_S`, so the total blocks by `Z̄ = s = ±1`, and the bath contributes
`Z_B(s) = Σ_j e^{−β(E_j + sλb_j)}`. On the degenerate code space `ρ_MF ∝ Z_B(s)`, hence

> ### `⟨Z̄⟩ = [Z_B(+1) − Z_B(−1)] / [Z_B(+1) + Z_B(−1)]`

**PREDICTIONS:**
1. **`λ → 0` gives `⟨Z̄⟩ → 0`** for every coupling — the weak-coupling limit recovers ordinary Gibbs and
   F-9's exclusion. **This is the control that makes any nonzero a strong-coupling effect.**
2. **A weight-`d` coupling gives `⟨Z̄⟩ ≠ 0` at finite `λ`, rising at order `λ¹`**, matching the closed
   form above. **If so, F-13 clause (1) IS CONTRADICTED: a detailed-balance environment forms a record.**
3. **A weight-1 coupling gives a bias SUPPRESSED as `λ^d`** — nonzero, but at high order, because the
   system must virtually leave the code space. **If so, clause (2) SURVIVES in a scaling form** —
   weight-`d` couplings act at first order, weight-1 only at order `d`.

**WHAT EACH OUTCOME MEANS.** If 2 and 3 both hold, **F-13's clause (1) is wrong as stated and must be
replaced by a statement about coupling weight and order, not about detailed balance at all.** If 2
fails, the charge is dismissed and clause (1) stands.

---

## O-11 RESULT — **F-13 CLAUSE (1) IS CONTRADICTED. A THERMAL ENVIRONMENT FORMS A RECORD.**

`LANE_O11_MEANFORCE`. Toric code 2×2 (dim 256) + a 4-level bath, `β=2`, mean-force state
`ρ_MF = Tr_B[e^{−βH_tot}]/Z`. **This is a genuinely thermal state of the TOTAL — detailed balance with
respect to `H_tot`.**

| `λ` | `⟨Z̄⟩`, coupling `A = Z̄` (weight `2 = d`) | closed form | `A = Z_e` (weight 1) |
|---|---|---|---|
| 0.000 | 0.00000000 | 0.00000000 | 7.197e-15 |
| 0.100 | **−0.15512514** | −0.15512514 | 4.053e-15 |
| 0.400 | **−0.54832339** | −0.54832339 | −2.489e-15 |
| 0.800 | **−0.82457285** | −0.82457285 | −4.504e-15 |

**CLOSED FORM AGREEMENT `7.197e-15`.** `λ→0` control **PASSES** — weak coupling recovers ordinary
Gibbs and F-9's exclusion, so the effect is genuinely strong-coupling. Slope in `λ` = **0.9947**, as
predicted at order 1.

> ### **PREDICTION 2 CONFIRMED. F-13'S CLAUSE (1) IS WRONG: THE RFP DOES NOT NEED NON-EQUILIBRIUM DYNAMICS. THE "NON-EQUILIBRIUM" REQUIREMENT WAS AN ARTEFACT OF THE WEAK-COUPLING (DAVIES) LIMIT, NOT A FACT ABOUT RECORDS.**

### PREDICTION 3 WAS FALSIFIED — **AND IN THE PROGRAM'S FAVOUR**

I predicted weight-1 couplings would bias at order `λ^d`. **They do not bias at all.**

| all 24 weight-1 couplings, `λ = 0.8` | count | `max\|⟨Z̄⟩\|` |
|---|---|---|
| **commute** with the writer `X̄` | 20 | 7.177e-15 |
| **do NOT commute** with `X̄` (`‖[A,X̄]‖ = 32`) | **4** | **5.317e-15** |

Noise floor (no coupling at all): `1.550e-15`.

**THE FOUR THAT BREAK THE WRITER'S SYMMETRY STILL GIVE ZERO.** So **the suppression is
KNILL–LAFLAMME, not symmetry** — and **clause (2) survives in its STRONGEST form: a weight-1 coupling
produces no bias at any coupling strength, not merely a suppressed one.**

### AND A THIRD CONDITION F-13 DID NOT HAVE

| weight-2 coupling | `‖[A,X̄]‖` | `\|⟨Z̄⟩\|` |
|---|---|---|
| `Z̄` | 32.000 | **0.824573** |
| `Z̄₂` | **0.000** | **0.000000** |

**Weight `≥ d` is NECESSARY BUT NOT SUFFICIENT. The coupling must also break the symmetry generated by
the writer** — F-8, reappearing at the level of `H_tot`.

> ### **THE RFP'S REQUIREMENTS, CORRECTED:**
> ### **(1) a coupling of weight ≥ `d`, and (2) that coupling must not commute with the writer.**
> ### **NOTHING ABOUT EQUILIBRIUM. AN ORDINARY THERMAL ENVIRONMENT SUFFICES.**

### THIS DISSOLVES O-9's PREMISE

O-9 asked *"what supplies a non-local, NON-EQUILIBRIUM coupling?"* **The non-equilibrium half is gone.**
What remains is the ordinary physical question **"what supplies a coupling extended over a region at
least as large as `d`?"** — a question about environments, not a demand for exotic dynamics.

### DEFECT CAUGHT IN THE FIRST RUN

The v1 script fitted a log-log slope to the weight-1 column and its **auto-verdict announced "clause (2)
does NOT survive."** That column was **noise** — magnitude `~1e-15` with the sign oscillating — so the
slope `0.2741` was a fit to zero. **No noise floor had been established.** This is **D-8 in a second
costume: a computed verdict is not a check unless the noise floor is measured alongside it.**

---

## F-1 / T-III.6 — **THE ARROW.** PREDICTION AND IMPORT-GUARD REGISTERED BEFORE THE RUN

**WHY THIS DECIDES THE PROGRAM.** Every clause of the anchor is **time-symmetric** — `R=R†`, `R²=I`,
`[H,R]=0`, `[L_k,R]=0`, writable, protected. **Not one mentions time.** Without an arrow we have a
theory of **protected bits**, not of **records**, and the charter's own word is unearned.

### THE DESIGN PROBLEM, AND WHY THE OBVIOUS TEST IS VOID

"Run the formation backwards" does not work as stated. With coupling `A = Z̄`, **`[Z̄, H_tot] = 0`**, so
`⟨Z̄⟩` is a **constant of motion** and unitary evolution never forms anything. The mean-force state is
reached by **relaxation**, and `ρ_MF` is itself **stationary** — an equilibrium state has no arrow.
**Testing reversibility of a stationary state would have measured nothing.**

### THE TEST THAT REPLACES IT

`ρ_MF` blocks by `Z̄ = s`, with bath state `e^{−β(H_B + sλB)}/Z_B(s)` in each block. **Those are
DIFFERENT bath states.** So the bath's state depends on the record's value, and the quantity is the
**Holevo information the bath holds about the record**:

> `χ(Z̄ : B) = S(ρ_B) − Σ_s p_s S(ρ_B^s)`

**AND THE ARROW STATEMENT:** local unitaries **cannot change** `χ`. If the bath holds a copy, **no
system-only operation can unmake it.** Forming needed a coupling; **unforming needs that coupling AND
control of the bath, which the system does not have.** That asymmetry is the candidate arrow.

### **THE IMPORT GUARD — AND IT IS THE POINT OF THE LANE (D-1)**

**If `χ(O:B) > 0` for EVERY observable, this is DECOHERENCE and we have imported einselection.** The
program's concept list already contains einselection and quantum Darwinism. **A general arrow is not
ours.**

> ### **DISCRIMINATOR: does the arrow carry the RECORD'S OWN THRESHOLD?**
> ### **`χ(O:B) > 0` only for observables of weight ≥ `d` ⟹ the arrow is RECORD-LEVEL.**
> ### **`χ(O:B) > 0` for weight-1 observables too ⟹ it is ambient decoherence, and F-1 is NOT ours.**

**PREDICTIONS:**
1. `χ(Z̄:B) > 0` for a weight-`d` coupling, matching a closed form from `Z_B(±1)`.
2. `χ(Z̄:B) = 0` for a weight-1 coupling — **the bath learns nothing about the record.**
3. **Local unitaries leave `χ` exactly invariant** — measured, not asserted.
4. **THE DECIDING ONE:** sweeping observables by weight, `χ` turns on at **weight `= d`** and is **zero
   below it**. If instead low-weight observables also carry `χ`, **the arrow is imported and this lane
   returns a negative result.**

---

## F-1 / T-III.6 RESULT — **THE ARROW IS REAL AND IT CARRIES THE RECORD'S OWN THRESHOLD. IT IS NOT AMBIENT DECOHERENCE.**

`LANE_F1_ARROW`. Toric code 2×2 + 4-level bath, mean-force state, `β=2`, `λ=0.8`.

### 1–2. THE BATH HOLDS INFORMATION ABOUT THE RECORD

| coupling | weight | `‖[A,X̄]‖` | **`χ(Z̄:B)` bits** | closed form |
|---|---|---|---|---|
| `Z̄` (logical) | 2 = `d` | 32.000 | **0.11448276** | **0.11448276** |
| `Z̄₂` (logical) | 2 = `d` | **0.000** | **0.00000000** | |
| `Z_e` (single site) | 1 | 32.000 | **0.00000000** | |
| identity | 0 | 0.000 | 0.00000000 | |

Closed form **exact**. And the same two conditions as F-13/F-15 reappear: **weight ≥ `d` AND
non-commuting with the writer.**

### 3. THE IRREVERSIBILITY, MEASURED

| | |
|---|---|
| covariance instrument check | `9.992e-16` **PASS** |
| **`I(S:B)` under 12 random system-only unitaries** | **`3.686e-14` — EXACTLY INVARIANT** |
| `χ` about the *fixed* label `Z̄` | `1.145e-01` — **moves, and does** |

> **A system-only operation can MOVE which observable the bath knows about. It cannot change `I(S:B)`
> at all. THE CORRELATION IS NOT ERASABLE FROM INSIDE.**

### 4. **THE DECIDING TEST — AND THIS IS THE RESULT**

Sweeping **all** observables by weight in the same state:

| observable weight | # swept | `max χ(O:B)` |
|---|---|---|
| **1** | 24 | **0.00000000** |
| **2 = `d`** | 252 | **0.11448276** |

**AND THE DISCRIMINATION THE IMPORT GUARD WAS BUILT FOR** — from part (c):

| coupling | weight | `I(S:B)` | `χ(Z̄:B)` |
|---|---|---|---|
| `Z_e` | **1** | **0.04549256** | **0.00000000** |
| `Z̄` | 2 = `d` | 0.11448276 | 0.11448276 |

> ### **A WEIGHT-1 COUPLING DOES ENTANGLE THE ENVIRONMENT WITH THE SYSTEM — `I(S:B) = 0.045` — AND TRANSFERS **ZERO BITS ABOUT THE RECORD**.**
> ### **AMBIENT DECOHERENCE IS PRESENT AT WEIGHT 1. THE RECORD'S ARROW IS NOT. IT TURNS ON AT `d`.**

**So the arrow is not einselection wearing our vocabulary. Decoherence and the record's arrow are
separated by a measured threshold, and that threshold is the record's own.**

### THE STATEMENT, IN THE RECORD LEVEL'S OWN TERMS

> **A record is of the past because the environment holds Holevo information about it, and no
> operation the system can perform removes that. Acquiring it requires the environment to couple over
> a region at least as large as `d` — the same `d` that protects it and the same `d` that forms it.**

### WHAT IS **NOT** ESTABLISHED — stated plainly, and it keeps this PARTIAL

1. **This is an equilibrium state, not a dynamical history.** F-1 asked for *"a before with no record
   and an after with one."* **One-wayness is measured; the temporal formation was not simulated.**
2. **The bath is 4 levels.** Redundancy — many environment fragments each holding a copy, which is
   what makes a record *objective* — is **untested**.
3. **One-wayness is necessary for "of the past". Whether it is sufficient is not settled here.**

### DEFECT IN v1, CAUGHT BY ITS OWN CONTROL

v1 reported *"prediction 3 falsified — local unitaries change `χ`."* **It conjugated the observable
the wrong way (`Q†Z̄Q` where covariance needs `QZ̄Q†`), so it measured a rotated observable.** The
corrected covariance check returns `9.992e-16`. **The reported falsification was mine, not nature's.**

---

# THE MODEL — **RECORDS CONSTRUCTED FROM `(H, {L_k})` AND NOTHING ELSE**

`model/record_model.py`. The principal: *"We need a model that can be used to construct records from
first principles without importing values."*

**INPUT: a Hamiltonian and a set of Lindblad operators. NO lattice, NO gauge group, NO temperature, NO
coupling constant, NO code, NO geometry.** Every step is a registered theorem: **C-9** (clause (ii)
puts `R` in the commutant of the `*`-algebra), **C-10** (proper-subalgebra criterion), **C-11 / O-4**
(admissible flipper ⟺ `Tr(P_E R) = 0`), **C-12** (the repaired converse).

### VALIDATION — **12 PASS, 0 FAIL**, each case a registered result re-derived from first principles

| | |
|---|---|
| odd dimension ⟹ no record is ever writable (C-11 corollary) | `C³, C⁵, C⁷` → 0 records |
| O-1 gap (a): degenerate `H` + **one generic jump** kills every record | 0 over 6 draws |
| **positive control**: same `H`, **diagonal** jumps | records in all 6 |
| O-1's **trap**: `max‖P_E L P_E‖ = 0` yet **no record** | 0 records |
| **control**: return legs deleted | records exist |
| gauge carrier handed in as a bare matrix | records found |

## **C-14 — THE RECORD-COUNT LAW, FROM THE CLAUSES ALONE**

> ### `k = min over eigenspaces E of v₂(m_E)` — the **2-adic valuation** of the multiplicities

**22 spectra, 22 PASS, 0 FAIL**, including odd, non-power-of-2 and mixed multiplicities.

**THE DERIVATION IS THE CLAUSES.** A record must be non-trivial (iii) and **trace-balanced on every
eigenspace** (iv = C-11); an independent **family** must split every joint block **evenly**. Each
independent record therefore **halves every eigenspace**, so the family size is bounded by how many
times every multiplicity can be halved.

**AND THE CONTROL SHOWS IT HAS CONTENT.** The naive `floor(log₂ min m_E)` predicts `[3,3] → 1`,
`[6,6] → 2`, `[5,5] → 2`. **Measured: 0, 1, 0.** The count is set by the **valuation**, not the size.

> ### **IT DERIVES THE GAUGE RESULT WITHOUT TOPOLOGY.** The toric code's ground multiplicity is `4`,
> ### so `k = v₂(4) = 2 = 2g`. **The index law G-7 is a CONSEQUENCE on that carrier, not the source
> ### of the count.**

## **C-15 — FOUR OF THE FIVE CLAUSES ARE CARRIER-FREE. THE FIFTH IS NOT.**

Clauses **(i)–(iv) are computable from `(H,{L_k})` alone.** **Clause (v) is not** — it requires a
**locality structure**, which is carrier data and is not derivable from the pair. The model **refuses
to guess it** and raises rather than supplying a default.

**Any claim resting on (v) therefore inherits a carrier.** That is the precise scope boundary of the
whole program, and it is now enforced in code rather than remembered.

## MULTI-RECORD — **PF-7 SUBSTANTIALLY ADVANCED**

| carrier | records | **independent family** | `log₂ dim` | **independently writable** |
|---|---|---|---|---|
| `C⁴` | 3 | **2** | 2 | **2 of 2** |
| `C⁸` | 35 | **3** | 3 | **3 of 3** |
| `C¹⁶` | 6435 | **4** | 4 | **4 of 4** |

**F₂-independence of sign vectors is the WRONG notion** — `k` independent bits need `2^k` joint
eigenspaces. The right criterion is that a new record **splits every existing joint block evenly**,
and under it **every member of the family is independently writable: one record can be flipped without
disturbing any other.**

**DEFECTS CAUGHT DURING THE BUILD, all by the validation:** the writer test was over-strict (asked one
unitary to fix 35 records at once); the XOR basis stored **unreduced** vectors so nothing ever reduced
to zero; and F₂-independence gave `n−2` where the answer is `log₂ n`.

---

## PF-3 — **DEF-A SURVIVES PART A UNIVERSALLY. PART B FOUND AN ERROR IN G-11 INSTEAD.**

### PART A — THE BALANCE LEMMA IS CARRIER-FREE

Tested **through the model**, so no carrier is involved at all: for every spectrum tried
(`[4] [8] [4,4] [2,2] [6,2] [4,2] [2,2,2,2] [8,4] [4,4,4]`), **every record the model builds is
balanced on every eigenspace, and an admissible writer is constructed and verified for every one.**

**POSITIVE CONTROL, and it fires:** an **unbalanced** involution (`Tr R = +2`) has **no** flipper —
best `‖U†RU + R‖` over **4000** random unitaries is **2.0000**, exactly the proved bound.

> **Clause (iv) under DEF-A ⟺ `Tr(P_E R) = 0` is UNIVERSAL, not a fact about the toric code.**

### PART B — **VOID AS A TEST OF DEF-A, AND IT INDICTS G-11 INSTEAD**

On a **bouquet of two triangles** (pinch-point non-manifold, dim 64), contractible regions = forests:

| | |
|---|---|
| any-unitary flippers inside a contractible region (**control**) | **1830** |
| **admissible** flippers inside a contractible region | **171** |

**That reads as DEF-A failing. It is not.** The carrier itself fails clause (v):

| carrier | links | **Z-systole** | **X-cosystole** | **d = min** |
|---|---|---|---|---|
| torus 2×2 | 8 | 2 | 2 | **2** |
| torus 3×3 | 18 | 3 | 3 | **3** |
| **bouquet of 2 triangles** | 6 | 3 | **1** | **1** |
| **bouquet of 3 triangles** | 9 | 3 | **1** | **1** |

## **ERRATUM TO G-11 — ITS DISTANCE COLUMN REPORTED THE SYSTOLE, NOT THE CODE DISTANCE**

**G-11 registered `d = 3,3,3` for the bouquet family. The true distance is `1`.** The bouquet has no
2-cells, so every single-edge `X` is a non-trivial `X`-logical — **cosystole 1**. G-12 had already
stated that **R3 needs two numbers**; G-11's own table was never corrected to use them.

**AND THE CORRECTION STRENGTHENS G-11's POINT.** The bouquet does not have *"unbounded capacity with
fixed protection."* It has **unbounded capacity with NO protection at all.** Capacity and protection
are independent — which is what G-11 was for, now stated correctly.

**PF-3 STATUS: part A closed, part B still owed.** The test needs a **non-manifold carrier with BOTH
distances ≥ 2**; the one chosen was invalid, and the run said so rather than reporting a verdict.

---

## PF-2 — PREDICTION REGISTERED BEFORE THE RUN

F-1 is PARTIAL for two stated reasons: it measured one-wayness on an **equilibrium** state, and
**redundancy** was untested (O-13). Both are addressed dynamically here.

### THE DESIGN POINT THAT MAKES A DYNAMICAL TEST POSSIBLE

With coupling `A = Z̄`, `[Z̄, H_tot] = 0`, so **`⟨Z̄⟩` is a constant of motion** — the record's *value*
never changes under unitary evolution. **That is not an obstacle, because forming a record is not the
record changing value. It is the ENVIRONMENT COMING TO HOLD INFORMATION ABOUT IT.** So the observable
is `χ(Z̄:B)(t)`, which starts at **exactly 0** for a product state and can grow while `⟨Z̄⟩` is fixed.

> **BEFORE = no record information anywhere. AFTER = the environment holds it. That is the history
> F-1 asked for, and it is measurable.**

**PREDICTIONS:**
1. `χ(Z̄:B)(t)` grows from **exactly 0** at `t=0` while `⟨Z̄⟩` stays constant — **the before/after.**
2. **CONTROL: full time reversal returns `χ → 0`.** Unitary evolution is exactly reversible, so if the
   instrument reported irreversibility here it would be manufacturing it. **The arrow must NOT show up
   in the closed dynamics.**
3. **REDUNDANCY (O-13): each bath fragment independently holds `χ > 0`** — many fragments each with a
   copy, which is what makes a record objective instead of merely correlated.
4. **Redundancy carries the record's threshold: a weight-1 coupling gives every fragment ZERO** —
   F-18 extended from the whole bath to its parts.

**WHAT THE ARROW THEN IS, STATED HONESTLY IN ADVANCE:** if 1–3 hold and 2 confirms the closed dynamics
is reversible, the arrow is **RELATIVE, not absolute** — one-way with respect to what the *system* can
do (F-19: `I(S:B)` is invariant under every system-only operation), while the joint evolution remains
reversible. **That is a weaker and more honest claim than "time has a direction", and it is the one
the measurements can support.**

---

## PF-2 RESULT — **THE HISTORY IS REAL, THE CLOSED DYNAMICS IS REVERSIBLE, AND THE RECORD IS REDUNDANT**

`LANE_PF2_DYNAMICAL`. Toric code 2×2 + a **3-qubit** bath, so the bath **has fragments**. Initial
state is a **product**, so `χ = 0` exactly at `t=0`.

### 1. THE BEFORE/AFTER F-1 ASKED FOR

| `t` | `⟨Z̄⟩` | **`χ(Z̄:B)` bits** | `I(S:B)` |
|---|---|---|---|
| 0.00 | −0.000000 | **0.00000000** | −0.00000000 |
| 0.50 | −0.000000 | 0.81447230 | 0.81447230 |
| 1.00 | −0.000000 | **0.97527192** | 0.97527192 |
| 4.00 | −0.000000 | 0.90811968 | 0.90811968 |

> **`⟨Z̄⟩` is EXACTLY CONSTANT throughout — the record's value never changes. What changes is that the
> environment comes to hold it.** Before: no record information anywhere. After: ~1 bit in the bath.
> **That is the history, and forming a record was never the record changing value.**

### 2. THE CONTROL — **AND IT MUST PASS, OR THE RESULT IS MANUFACTURED**

| | |
|---|---|
| `χ` at `t = +4.0` | 0.90811968 |
| `χ` at `t = 0` | **0.00000000** |
| `χ` at `t = −4.0` | **0.90811968** — equal to `t=+4` |

**The closed dynamics is exactly reversible and time-symmetric.** Had the instrument reported
irreversibility here it would have been inventing it.

### 3. REDUNDANCY — **O-13 CLOSED**

| coupling | weight | whole bath | frag {0} | frag {1} | frag {2} |
|---|---|---|---|---|---|
| `Z̄` (logical) | 2 = `d` | **0.90811968** | **0.789366** | 0.048377 | **0.678602** |
| `Z_e` (single site) | 1 | **0.00000000** | 0.000000 | 0.000000 | 0.000000 |

**Fragments independently hold the record — two of three carry most of a bit each.** Separate observers
can each learn it without communicating, which is what makes a record **objective** and not merely
correlated. **And redundancy carries the record's threshold: a weight-1 coupling gives EXACTLY ZERO to
the whole bath AND to every fragment.** F-18 extends from the bath to its parts.

### **THE ARROW, STATED AT THE STRENGTH THE MEASUREMENTS SUPPORT**

> **It is RELATIVE, not absolute.** The joint evolution is exactly reversible (control above). The
> one-wayness is **with respect to what the SYSTEM can do**: `I(S:B)` is invariant under every
> system-only operation (F-19), so nothing the system does removes the environment's copy — and by the
> redundancy above, there are several copies to remove.
> **This is weaker than "time has a direction" and it is what was actually measured.**

### TWO LIMITS, STATED

1. **`χ(t)` is NON-MONOTONIC** — 0.975 at `t=1`, 0.787 at `t=2`, 0.908 at `t=4`. **Partial recurrence,
   because a 3-qubit bath is small.** A macroscopic environment would not recur on any relevant scale,
   but that is asserted here, not measured.
2. **Redundancy is UNEVEN** — 0.789, 0.048, 0.679. The fragments have different energies and acquire
   information at different rates. **Redundancy is present; uniform redundancy is not claimed.**

---

## D-9 / F-20 ERRATUM — **DO NOT GRADE A FINDING AGAINST AN EXTERNAL STANDARD. STATE THE MECHANISM.**

The principal, on the registrar's *"that's weaker than 'time has a direction'"*: **"What does this mean
'weaker'???? We just need to know how it works."**

**"Weaker" compares our result to a yardstick nobody asked for.** The charter asks how records form.
*"Time has a direction"* is a cosmological question the program never posed, and ranking against it is
**D-1 applied to framing instead of to physics** — an imported measure used as a grade.

> **D-9. A RESULT IS DESCRIBED BY ITS MECHANISM, NEVER RANKED AGAINST AN EXTERNAL STANDARD.**
> **TELL:** the words *weaker*, *stronger*, *only*, *merely*, *falls short of*, or any comparison to a
> claim the program did not make.

### F-20 RESTATED — **THE MECHANISM, WITH NOTHING GRADED**

**1. The coupling makes the environment's evolution CONDITIONAL on the record.**
`H_tot = H_S + H_B + λ Z̄ ⊗ ΣX_j`. The bath evolves under `H_B + λΣX_j` when `Z̄ = +1` and under
`H_B − λΣX_j` when `Z̄ = −1`. **Two different bath Hamiltonians, selected by the record.**

**2. The two branches drive the bath to different states. That divergence IS the copy.**
`χ` measures exactly how far apart they are: `0.00000000` at `t=0`, `0.97527192` at `t=1`.

**3. The record's value never changes — it is being READ, not written.**
`[Z̄, H_tot] = 0`, so `⟨Z̄⟩` is exactly constant. **Every bath degree of freedom reads the same bit
and none of them can alter it.**

**4. Several copies form because the coupling touches several bath degrees of freedom.**
Each qubit's evolution is separately conditioned, so each separately diverges: `0.789 / 0.048 / 0.679`
bits. **Redundancy is not an extra property. It is what happens when more than one thing is coupled.**

**5. Only a weight-≥`d` operator can make the environment's evolution conditional at all.**
By Knill–Laflamme a weight-`<d` operator acts as a **scalar** on the code space, so it takes the same
value on both branches, so the bath's Hamiltonian is the **same** either way, so **nothing diverges and
no copy forms** — measured as exactly `0.00000000`, whole bath and every fragment.
**The threshold is `d` because `d` is the point at which an operator can tell the record's two values
apart.**

**6. The system cannot undo it because the system's operators do not act on the environment.**
`I(S:B)` is invariant under every system-only unitary (`3.686e-14`). Undoing the divergence means
running the bath's evolution backwards, which is not an operation the system has.

> ### **THAT IS HOW IT WORKS.** A record forms when an operator large enough to distinguish its values
> ### makes the environment's evolution depend on which value it has, so the environment falls into a
> ### different state for each — in every degree of freedom the coupling reaches, and irreversibly
> ### from the system's side, because the system has no operators there.

---

## PF-5 / O-9 — PREDICTION REGISTERED BEFORE THE RUN

**THE CHARGE.** Formation needs a coupling of weight ≥ `d` (F-13). **Physical environments couple
LOCALLY** — the standard model of an environment is `H_int = Σ_i λ_i A_i ⊗ B_i` with every `A_i` acting
on one site. **If only extended couplings form records, records do not form in nature**, and the whole
process is proved conditional on something the world does not supply.

**THE RESOLUTION I EXPECT, AND IT IS ALREADY IMPLIED BY C-13.** Clause (v) forbids a **single**
contractible operation. It says nothing about **products**, and C-13 recorded that *products of
contractible operations DO reach the record* — that is Theorem D. **A sum of local couplings is not one
weight-1 operator; at order `d` in `λ` its terms act together and reach the record.**

**PREDICTIONS — three couplings, one carrier:**
1. **A SINGLE local term**, `λ Z_e ⊗ X_0`: `χ = 0` at **every** `λ`. Powers of one weight-1 operator
   are `I` or `Z_e`; they never reach the logical. **This is the control.**
2. **A SUM of local terms**, `λ Σ_l Z_l ⊗ X_{l mod 3}`: **`χ > 0`**, suppressed — log-log slope in `λ`
   of order **`2d`** (the amplitude to reach the record is `~λ^d`, and `χ` is second order in it).
3. **A weight-`d` coupling**, `λ Z̄ ⊗ ΣX_j`: `χ > 0` at **first** order, slope ~2.

**IF 1–3 HOLD, PF-5 IS ANSWERED AND O-9's PREMISE DISSOLVES.** The world does not need to supply an
exotic extended coupling. **An ordinary local environment forms records, at order `λ^d`** — so the same
`d` that protects the record, gates its formation, and sets the environment's threshold also **sets how
slowly it forms.** That is one quantity doing four jobs, and it is measurable.

**IF 2 RETURNS ZERO,** local environments genuinely cannot form records, O-9 stands as a real
obstruction, and the process is conditional on an input nothing supplies.

---

## PF-5 / O-9 RESULT — **AN ORDINARY LOCAL ENVIRONMENT DOES FORM RECORDS. THE INPUT IS SUPPLIED.**

`LANE_PF5_LOCALBATH`. Toric code 2×2 + 3-qubit bath, unitary from a product state, `t = 4.0`.
Noise floor `8.882e-16`.

| `λ` | **1. SINGLE local term** | **2. SUM of local terms** | **3. weight-`d` coupling** |
|---|---|---|---|
| 0.020 | 1.221e-15 | 0.00000031 | 0.00224457 |
| 0.100 | 0.000e+00 | 0.00016142 | 0.05127919 |
| 0.400 | 0.000e+00 | 0.01906127 | 0.47151593 |
| **log-log slope in `λ`** | **exactly zero** | **3.8988** *(predicted `2d = 4`)* | **1.9457** *(predicted 2)* |

### THE DISCRIMINATOR IS THE WHOLE POINT

**One local term gives EXACTLY ZERO at every coupling strength.** Powers of a single weight-1 operator
are `I` or that operator — **they never reach the logical, at any order.**
**A SUM of local terms gives `λ^{2d}`.** You need **`d` different local terms acting together**.

> ### **THAT IS C-13 MEASURED DYNAMICALLY.** Clause (v) forbids a *single* contractible operation;
> ### **products of them DO reach the record**, and here is the rate.

### O-9's PREMISE IS DISSOLVED

O-9 asked *"what supplies a coupling extended over a region of size ≥ `d`?"* **Nothing has to.** The
standard model of an environment — `Σ_i λ_i A_i ⊗ B_i`, every term on one site — **forms records
already.** What it costs is **suppression by `λ^{2d}`.**

> ### **AND THE SAME `d` NOW DOES FOUR JOBS, ALL MEASURED:**
> | | |
> |---|---|
> | **protects** the record | splitting `~ ε^d` (Thm D, slope 2.000) |
> | **gates formation** | threshold `= d` (F-11, six codes, `d = 1,2,3`) |
> | **gates what the environment can learn** | `χ` turns on at weight `d` (F-17, F-22) |
> | **sets how slowly a local environment forms it** | `χ ~ λ^{2d}` (here, slope 3.8988) |
>
> **One quantity, four roles. None of them was put in by hand.**

### **AND IT IS A CANDIDATE FOR PF-6 — THE FIRST ONE THE PROGRAM HAS HAD**

`X-4`, `T-VI.3` and `T-VI.4` have stood at **BLOCKED, "none identified"** since they were written.
**`χ ~ λ^{2d}` is quantitative, is about ordinary local environments, and involves only codes and baths
— objects that exist on hardware.** A record of distance `d` is **exponentially slow in `d`** to form.
Registered as the first candidate; **whether it distinguishes this account from any other is not yet
established and PF-6 stays open.**

---

## D-8, THIRD INSTANCE — **IN THE REPRODUCTION HARNESS ITSELF**

The first full-corpus run from a snapshot copy failed, and **both failures are D-8**:

1. **`ROOT` was derived from `$BASH_SOURCE`.** A copy in `/tmp` resolved `ROOT` to `/`, so the run
   `cd`'d to the filesystem root, found no lanes, and reported *"no sealed output — SKIP"* for
   everything. **It looked like a clean run over an empty corpus.** Fixed: `WAC_ROOT` override, and the
   script now **refuses to start** unless it can see `CORE_FRAMEWORK_V001.md`.
2. **THE GRID CHECK ANNOUNCED SUCCESS WITH ITS INPUTS MISSING.** `shasum` failed on both sides, both
   variables were the empty string, `"" = ""` compared true, and the harness printed *"regenerating the
   grid reproduces it byte-for-byte."* **A check that passes when the thing it checks is absent.**
   Fixed: an empty hash is now a failure.

**The verification harness had the exact defect it was built to catch.** Recorded because D-8 has now
fired three times — `check(..., True)` in the O-4 lane, a log-log slope fitted to a `1e-15` column in
O-11, and here.

---

# PF-6 — **ANSWERED NO. NOTHING IN THIS ACCOUNT IS UNPREDICTED BY A RIVAL, AND TWO REGISTERED CLAIMS ARE FALSE.**

Six rival accounts, every claimed difference adversarially adjudicated. **The deflationary lens won.**

## WHAT IS NOT OURS — the attribution, in full

| our result | whose it already is |
|---|---|
| formation threshold `= d` | **Knill–Laflamme** (PRA 55, 900, 1997) + Bravyi–Terhal cleaning lemma. **It is the DEFINITION of `d`** |
| protection `ε^d` | **Bravyi–Hastings–Michalakis**, J. Math. Phys. 51, 093512 (2010) |
| `χ` exactly 0 below `d`, non-zero at `d` | **Cleve–Gottesman–Lo** (1999); and now inside quantum Darwinism itself — **Girard–Cheng–Cao, arXiv:2606.06588, Thm 4** |
| weight-1 entangles but transfers 0 record bits | **information–disturbance**, Schumacher–Nielsen (1996) |
| redundancy across fragments | **quantum Darwinism** (Zurek, Zwolak–Zurek) + GCC Thm 4 |
| `χ ~ λ^{2d}` | `d`-th order degenerate perturbation theory, squared |
| a single local term gives exactly 0 at all orders | **group theory** — powers of one weight-1 Pauli are `I` or itself. **Algebra, not measurement** |
| `I(S:B)` invariant under system-only unitaries | mutual information is invariant under local unitaries. **Arithmetic** |
| clause (v)'s topological restriction | **TQO-1**, Bravyi–Hastings–Michalakis; stated in Kitaev 2003 |

> ### **"ONE QUANTITY, FOUR ROLES" IS ONE THEOREM READ FOUR TIMES. It is not four coincidences, and reporting the convergence is a check that cannot fail. F-24 IS WITHDRAWN.**

## **D-10 — A VERIFICATION WITH NO FAILING BRANCH IS AN INSTRUMENT CHECK, NOT A RESULT**

**`threshold = d` was verified on six codes at three distances — and `threshold = d` IS the definition
of `d`.** The verification had no branch on which it could have come out otherwise. **It confirmed the
instrument, not the account.**
**TELL:** a "verification" whose negation would be a contradiction in terms.

## THREE REGISTERED CLAIMS ARE FALSE OR MISSCOPED

### **F-25 WITHDRAWN — the eleventh withdrawn claim.** The exponent is `2n*`, not `2d`

`n*` = the minimum number of the environment's **available** local terms whose product acts
non-scalar on the code space. On `[[5,1,3]]` with `Z`-only couplings the measured slope is
**9.4619 → 9.8435 → 9.9612 → 9.9966 → 10 = 2n*`**, against `2d = 6`. `[[4,2,2]]` control: 3.9854 → 4.

**THE CAN-FAIL CHECK THAT DECIDES IT:** on the *same* code and estimator, widening the single-site set
from `{Z}` to `{X,Y,Z}` moves the slope to **5.9078 → 5.6933 → 6**, while `{Z}` stays at 10.
**The exponent tracks the COUPLING ALGEBRA, not the code.**

**And our own file closes the escape:** `LANE_PF5_LOCALBATH/pf5_localbath.py:38` uses a `Z`-type,
record-commuting sum and its own docstring calls it *"an ORDINARY LOCAL environment."*

### **C-14 IS FALSE AS STATED.** Verified counterexample

`model/count_law.py` only ever tested the law with an **empty noise set**, so *"from the spectrum
alone"* was never exercised against a non-trivial commutant.

| collective decoherence, `H = J²`, `L = {Jx,Jy,Jz}` | multiplicities | `min v₂(m_E)` | **model `k`** |
|---|---|---|---|
| **n = 3** | `[4,4]` | **2** | **0** |

**C-14 predicts 2. The truth is 0.** And these are **noiseless subsystems built and operated
experimentally** — Kwiat et al. (2000), Kielpinski et al. (2001), Viola et al. (2001).
**C-14 must be restricted to empty `L`, or restated in terms of the minimal projections of `A′`** —
which is what the model actually computes.

### THE PROTECTION LAW IS RESCOPED

`ε^d` holds for a **generic** single-site perturbation — which is what `LANE_O23_SPLITTING` actually
measured (3.0000 on `[[5,1,3]]` and Steane). Under a **`Z`-only** perturbation the slope is
**5.0001 = n\***. **The lane is right; the unqualified verbal form "a local perturbation splits as
`ε^d`" is not.**

### PF-5's EXPONENT NEEDS ITS BATH PREMISE NAMED

`slope(χ)` is **2.0064** at bath purity `q=0` and **3.9665** at `q=0.02`, while the amplitude-level
`slope(D)` holds at **1.9987–1.9992** throughout. **The factor of two is set by how the environment is
prepared. The record-intrinsic exponent is the amplitude exponent.**

## X-4, T-VI.3, T-VI.4 — **THEY DO NOT MOVE**

All three stay **BLOCKED**. **PF-6 is answered NO, and what it gains is an adjudicated reason in place
of silence.** Every arm run so far is closed unitary evolution with a small bath; **the one untested
place where amplitudes and rates need not stand in the squared relation is a genuinely dissipative
bath.**

---

## PF-6 ERRATUM — **THE REGISTRAR APPLIED A TEST THE REGISTER ITSELF REJECTS**

The principal: *"it's expected that we will use conventional physics… What physics has not accomplished
is to describe how those discoveries and observations combine to create records."*

**`STATUS_V001.md` section C already says this, and I ignored it:**

> **"THE TEST APPLIED: *does the assembly answer an open question?* — NOT *can the parts be named?*
> Naming parts does not make an assembly a restatement."**

**PF-6's finding stands: no measurement distinguishes this account from a rival, so X-4, T-VI.3 and
T-VI.4 remain BLOCKED.** That is a real limit on **empirical contact** and nothing here softens it.

**What does NOT follow is the conclusion I drew from it.** I wrote that the result *"materially reduces
what the program can claim as its own."* **It does not.** The charter never asked for a new physical
law. It asks **how EM, gravity and alpha COMBINE to form records** — and every component of a correct
assembly *should* have a named owner. **Finding that each part is known is what a right answer looks
like, not evidence against one.**

> ### **THE CLAIM WAS NEVER "A NEW LAW". IT IS "A PROCESS NOBODY HAS DESCRIBED".**
> ### **PF-6 tested the first and found nothing. It never tested the second.**

**WHAT SURVIVES PF-6 UNTOUCHED:** the assembly itself — that a record is a trace-balanced projection in
the commutant (C-12), written by an operation admissible only as a whole (C-13), carried by an EM
holonomy that is also its writer and its minimal forming coupling (A-EM3), formed when a coupling makes
the environment's evolution conditional on it (F-20), with the environment's copy unremovable from
inside (F-19). **No cited rival assembles these into a formation process; each owns one part.**

**AND THE THREE FALSE CLAIMS STAY FALSE.** F-25, C-14 and the protection scope were genuine errors,
caught by genuine counterexamples, and this erratum does not touch them.

**D-9 fired here and I missed it on myself:** *"materially reduces what the program can claim"* is a
**grade**, not a mechanism.

---

## O-16 — PREDICTION REGISTERED BEFORE THE RUN

**THE ONE UNTESTED ARM.** Every arm so far is **closed unitary evolution with a small bath**, where
information measures are the square of an amplitude by construction. **A genuinely dissipative
(Lindblad) bath is the one place amplitudes and rates need not stand in that relation.**

**THE MECHANISM THAT MAKES ME EXPECT A DIFFERENCE.** For a Pauli jump operator `L = P_i` with
`{P_i, Z̄} = 0`, the Lindblad term gives `Tr(Z̄ · (P_i ρ P_i − ρ)) = −2 Tr(Z̄ρ)` — **first order in the
rate `γ`, with no `d` in it at all.** A single-site jump decays the record immediately, where in the
unitary arm a single local term gave **exactly zero at every order**.

**PREDICTIONS:**
1. **The record's decay exponent under single-site Lindblad jumps is 2 in `λ` (first order in
   `γ = λ²`), INDEPENDENT of `d`.**
2. **The static splitting exponent under the same operator set is `n*`** — as PF-6 established.
3. **Therefore the two exponents DIFFER, and the rival's squared rule does not carry over.**

**IF 1–3 HOLD, THE FINDING IS PHYSICAL AND SHARP:** **distance protects against STATIC perturbation and
does NOT protect against DISSIPATION.** Uncorrected, a distance-`d` code loses its record at the
physical rate — which is why quantum error correction needs an **agent actively correcting**, and this
account has none.

**IF THE EXPONENTS MATCH,** the rival's rule survives the dissipative arm too, O-16 closes negative,
and the account has no arm left in which it differs from perturbation theory.
