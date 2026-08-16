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
