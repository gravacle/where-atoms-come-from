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
