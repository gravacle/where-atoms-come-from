# W-19 RULING — THE CARRIER THRESHOLD

**Register head W-17. Program: WHERE ATOMS COME FROM, attempt two, first computation.**
Question put to the lanes: *what is the smallest gauge carrier on which a redundancy plateau can
actually be exhibited?*

**LINEAGE, DECLARED.** I am Opus 5. Lane A, lane B, the refuter of lane A and the refuter of lane B
are Opus 5. This is one block, continuous from W-07 — same model, same lineage, no independent
adversary anywhere in it. Every "independent" check below, including my own, is a same-lineage
check, and the standing correction of `feedback_read_sealed_corrections_first` applies to this
ruling as much as to the lanes: chains of this shape fail by being under-read, not by being
under-adversarial, and same-lineage adversaries share failure modes. Read the agreements in this
document as weaker evidence than the disagreements.

**MY OWN LANE.** `/Users/bgm/MB Work/where-atoms-come-from/W19_RULING/` — 6 scripts, 5 captured
outputs, sealed in `SEALS.sha256` (verifies 11/11). `rule_verify.py` is written from the mathematics
and imports nothing from any lane: no `zn_gauge.py`, no `lib_b.py`, no `rlib.py`, no `carriers.py`.
It reproduces the sealed LANE_T1_NEW_PROGRAM null by a third route —
**I(0:1) = 0.690763, I(0:2) = 0.384496, I(0:{1,2}) = 1.075259** — before any claim below is made
(`OUT_rule_main.txt` block 0). All four lane directories verify their own seals: **22/22, 18/18,
12/12, 14/14** (refuter A reports lane A as 21/21; lane A's `SEALS.sha256` has 22 entries and all 22
verify — a miscount, not a mismatch, and logged here rather than quietly fixed).

---

## 0. THE FRAME OBJECTIONS, REPORTED FIRST AND SEPARATELY

Standing gate, adopted 2026-08-17: frame objections are reported before the verdict and are never
counted as votes. Both refuters filed one. They are not the same objection and both are sustained,
one of them with a number this ruling produced.

**REFUTER A'S FRAME OBJECTION, verbatim in substance:** the question asks for a minimum without
naming the class it is a minimum over, so the class does all the work and the brief supplies none of
it — degree-2 vertices allowed gives 9, min degree ≥ 3 with multi-edges gives 9, simple gives 12,
girth ≥ 4 gives 15, girth ≥ 6 gives 21, nested-from-both-endpoints gives unreachable. Every one is a
correct answer to the question as written.

**REFUTER B'S FRAME OBJECTION, in substance:** read existentially over states, the answer is a
combinatorial triviality (its number: L ≥ girth + 4) saturable on demand with a stabiliser state on
any group for any pointer; read dynamically — the smallest carrier on which *the theory's own
Hamiltonian* produces a plateau — the question has content and neither lane answered it. And
"report a THRESHOLD, not an impression" is the instruction that did the damage: it rewards a number
and penalises the finding that the number is not the interesting object.

**RULING ON THE FRAME OBJECTIONS: BOTH SUSTAINED.** I add the measurement that settles the first one
quantitatively, which neither lane nor either refuter ran. Hold the carrier, the state, the coupling,
the system algebra, the number of fragments and δ all fixed, and move **only the partition** of the
environment into disjoint fragments:

| carrier | H_elec(S) | R_δ on lane A's BFS-cut partition | R_δ over 200 random equal-size disjoint partitions |
|---|---|---|---|
| dbl_chain9 (L=9) | 0.610000609 | **5 of 5** (sizes 2,1,2,1,2) | mean **2.185**, max 4 |
| tri_chain12 (L=12) | 0.663269043 | **5 of 5** (sizes 2,2,1,2,2) | mean **1.500**, max 4 |
| heawood (L=21) | 0.994827216 | **5 of 5** (sizes 2,4,8,4,2) | mean **0.400**, max 2, **zero in 127 of 200** |

`W19_RULING/OUT_rule_main2.txt` block 7. On lane A's own threshold carrier, at lane A's own coupling,
in the gauge-invariant channel, the reported R_δ = 5 becomes 0 in 63.5% of randomly chosen partitions
of the same environment into the same number of disjoint fragments. **R_δ as computed in this round
is a property of the partition, and the partition was chosen — correctly and knowingly — to be the
one the Gauss law guarantees.** That is refuter A's objection with a number on it.

I record one overreach in each objection. Refuter A writes that no value of the number, its own or
lane A's, is a step toward record formation; that is rhetorically true and evidentially too strong —
the round did produce one durable structural fact (the Gauss-surface identity, §5 below), which is a
real statement about how a Z_N gauge theory copies electric flux into its environment. It is a
theorem about the wrong quantity, not a null. Refuter B's objection still frames its own reading-1
"six" as an answer, when its own `what_falls` correctly says that six-link plateau is the redundancy
of a non-observable — a fact I confirm independently in §4.

---

## 1. VERDICT ON THE TWO LANES AND THEIR REFUTERS

| | verdict | why |
|---|---|---|
| **LANE A (the sweep)** | **INSTRUMENT SOUND, HEADLINE FALSE, SELF-FLAGS CORRECT** | Its 21-link threshold F1 is false as stated; its own §4 already declared the plateau weightless and its CF-1 already declared R_δ a graph distance. It was right about itself and wrong about the number. |
| **REFUTER A** | **SUSTAINED, and independently reproduced here** | Every load-bearing number of its minimality attack reproduces under my code. |
| **LANE B (the criterion)** | **CORRECT ON THE CHOICE-DEPENDENCE, WRONG ON BOTH THRESHOLDS** | Its B-7 (the verdict is not fixed by the link partition alone) is the finding of the round. Its reading-1 six-link threshold has no gauge content; its reading-2 twelve is over-derived. |
| **REFUTER B** | **SUSTAINED WITH ONE CORRECTION** | Its Perron-Frobenius theorem is right and I extend it to carriers it did not test. Its claim that lane A's plateau is "the redundancy of a non-observable" is true on theta and **false on lane A's own carriers**; see §4. |

### What I reproduced independently, and what I only seal-checked

**REPRODUCED under my own code, agreeing to the last printed digit:**

* Refuter A's two counterexamples. `tri_chain12` (V=8, L=12, cubic, girth 3, d=5): H(S) =
  **0.663269043** at g² = 0.50, four rule-A plateau points, R_δ = 5 on cuts [2,2,1,2,2].
  `dbl_chain9` = mg_chain(5) (V=6, L=9, all degrees exactly 3, d=5): H(S) = **0.610000609**, four
  points, R_δ = 5 on cuts [2,1,2,1,2]. `OUT_rule_main.txt` block 3.
* Refuter A's exhaustive search, re-run from scratch: connected simple graphs with min degree ≥ 3 —
  **V=5: 26 graphs, max d = 2; V=6: 1858, max d = 3; V=7: 236926, max d = 4**; and on V=8,
  **19355 labelled cubic graphs, 2520 of them reach d ≥ 5 at L = 12, every one of girth 3**.
  Counts identical to refuter A's. `OUT_rule_main4.txt`.
* Lane A's own live control, C1 arm B, on an independently constructed Heawood graph with different
  link labels and a differently built arm B: **I_ext/H(S) = 1.213848345** against arm A's
  1.000000000, H(S) = **0.307291435** at g² = 1.00. `OUT_rule_main3.txt` block 10b. Lane A's
  arithmetic is exact.
* Refuter B's Perron-Frobenius theorem, and **extended to lane A's min-degree-3 carriers, which
  refuter B did not test**: on theta_6, dbl_chain9 and tri_chain12 at g² = 0.30/1.00/3.00, the
  ground state of the *unconstrained* 2^L Hamiltonian has |⟨free|projected⟩| = **1.000000000** and
  min⟨G_v⟩ = **1.000000**, and E₀(full 2^L) = E₀(physical) to nine decimals. Control that can fail
  and does: flip the electric sign and min⟨G_v⟩ = **−1.000000000**. `OUT_rule_main.txt` block 2.
* Refuter A's latent instrument defect, confirmed by inspection: `zn_gauge.py:242` reads
  `return T @ T.T` — no conjugation. Latent on real states, fatal for any complex Hamiltonian.

**NOT INDEPENDENTLY VERIFIED — cited on the refuters' authority and their seals only:** refuter A's
above-ceiling orbit-RDM runs (Möbius-Kantor through Tutte-Coxeter, L = 24…45); refuter B's Galois
enumeration of gauge-invariant region algebras (2 / 16 / 67) and its "eight distinct verdicts across
the 16 algebras of a plaquette"; refuter B's theta_6 stabiliser broadcast and its L ≥ girth + 4.
These are the claims a genuinely independent adversary should be pointed at next, and this round has
none.

---

## 2. THE THRESHOLD

The brief demands the carrier, its link count, its cut, and the plateau — or the statement that no
carrier of this class works and the naming of a class that would. **This ruling gives both, because
the two are answers to two different questions and the round has conflated them.**

### 2a. THE THRESHOLD TO THE QUESTION AS ASKED — a number, delivered, and marked

Fixing the conventions the brief left open exactly as lane A fixed them (extended Hilbert space,
rule-A nested fragments, δ = 0.10, four points, min degree ≥ 3, simple graph):

> **CARRIER: `tri_chain12`. V = 8, L = 12 LINKS, cubic, girth 3, d = dist_{G−l}(u,v) = 5, C = 5,
> dim_phys = 2⁵ = 32.**
> Edges, link 0 first: (0,7) | (0,1) (0,2) (1,2) (1,3) (2,3) | (3,4) | (4,5) (4,6) (5,6) (5,7) (6,7).
> **CUT: S = link 0; E = the other 11 links.**
> **PLATEAU at g² = 0.50, Z₂ ground state:** H(S) = **0.663269043** bits;
> I(S:F)/H(S) = **1.000000, 1.000000, 1.000000, 1.000000** at |F| = 3, 5, 6, 9; then 2.000000 at
> |F| = 11 (purification). **Four plateau points. R_δ = 5** from five pairwise-disjoint genuine
> u–v cuts of sizes [2,2,1,2,2].
> **Multigraph floor: `dbl_chain9` (V=6, L=9, all degrees exactly 3, d=5), H(S) = 0.610000609,
> four points, R_δ = 5.** The general floor is L_min = ⌈3(P+2)/2⌉ for P plateau points.

**Lane A's 21 is refuted.** Its F1 proof (out_04_threshold.txt block 4b) counts BFS levels 1,2,4 on
each side of l, which requires no 4-cycle through u — global girth ≥ 5 — while the hypothesis it
states is girth *through l* ≥ 6, which constrains nothing about cycles avoiding l. Heawood is minimal
in the class of girth-6 cubic graphs, a hypothesis lane A's own criterion never uses. The exhaustive
search above closes it: no simple min-degree-3 graph on ≤ 7 vertices reaches d = 5, and 2520 of the
19355 labelled cubic graphs on 8 vertices do, at L = 12.

**AND THE NUMBER IS MARKED WEIGHTLESS, BY THE STANDING RULE.** Lane A said so of Heawood; it is true
of tri_chain12 and of every carrier in this class. Three Haar-random *physical* states on
tri_chain12 (complex Ginibre, no Hamiltonian anywhere in their construction) give the **identical**
curve to the ground state — 1.000000, 1.000000, 1.000000, 1.000000, 2.000000, four points, R_δ = 5 —
at H_elec(S) = 0.997064307 / 0.947928640 / 0.811703716. `OUT_rule_main.txt` block 5. The plateau
could not have failed. It is arithmetic about d, not a measurement of a carrier.

### 2b. THE THRESHOLD TO THE QUESTION WORTH ASKING — NO CARRIER OF THIS CLASS WORKS

Fix the criterion as it must be fixed for the answer to be about the theory rather than about the
bookkeeping. Five clauses, each forced by something measured in this round:

* **P1 — the system algebra is gauge-invariant.** On a single link the *only* non-trivial
  gauge-invariant subalgebra is alg{X_l}, the electric flux (lane B's T1; refuter B's exhaustion).
* **P2 — plateau fragments are pairwise disjoint**, which is the brief's own definition of R_δ.
* **P3 — |F| ≤ |E|/2**, the standard partial-information-plot convention.
* **P4 — at least four points.**
* **P5 — the plateau must be able to fail**: there must exist a state on the *same* carrier, cut,
  algebra and *pre-declared* partition that does not plateau.

**RULING: under P1–P5 there is no carrier in the class "Z_N pure lattice gauge theory on a finite
graph, system region a single link (or any forest)". Not at 9 links, not at 12, not at 21, not at
45, and not at any L.** The obstruction is not computational and does not sit beyond anybody's
ceiling; it is structural, and it has two halves that between them exhaust the class:

1. **Where the Gauss law carries the plateau, P5 fails identically.** For any physical state, any
   coupling, any N: if F contains a u–v cut in G−l and S ∪ F contains no cycle through l, then
   X_l = X(cut)⁻¹ exactly, so I(S:F) = H(S) exactly. Measured across four arms that should have
   separated and did not — tri_chain12 and heawood, CL channel, BFS-cut fragments: ground state at
   g² = 0.50, ground state at g² = 3.00 (where H_elec(S) = 0.001283123 and 0.008196771 — an *empty*
   record), and two Haar-random physical states, **all giving 1.000000 ×5, points = 5**.
   `OUT_rule_main2.txt` block 8a. Worse for lane A: its single live control dies here too. On arm B
   of its own C1 enclosure test the gauge-invariant channel reads **CHI/H = CL/H = 1.000000000, the
   same as arm A**, while only the gauge-variant extended-algebra channel moves to 1.213848345
   (`OUT_rule_main3.txt` block 10b). **Lane A has no live control in the gauge-invariant channel at
   all.**
2. **Where P5 can be satisfied, the plateau has no gauge content.** On d = 1 carriers the Gauss law
   constrains only the global parity, so a plateau on proper fragments must be carried by the state
   and *can* fail. It does fail, for everything the theory produces: on theta_6 with five disjoint
   single-link fragments, R_δ = **5** for the constructed electric GHZ and **0** for the ground state
   at g² = 0.10, 0.50 and 3.00 and for two Haar physical states; on theta_8 with seven disjoint
   fragments, R_δ = **7** for the electric GHZ and **0** for the magnetic GHZ, the ground state and
   Haar (`OUT_rule_main2.txt` block 8b/8c). So the criterion discriminates — and the thing it
   certifies is a classical repetition code in the X basis: CL(X₀:X_j) = 1.000000 for every j, on a
   state whose curves are identical on L bare qubits with no Gauss law and no plaquettes
   (`OUT_rule_main2.txt` block 9). **Zero variables moved with respect to gauge invariance.**

There is no third case. A forest region's gauge-invariant algebra is abelian and purely electric;
its only redundancy in the physical sector is the Gauss-surface identity, which is a theorem; and
any *additional* redundancy is a bare classical code that the gauge structure neither creates nor
forbids.

### 2c. THE CLASS THAT WOULD WORK, NAMED

Two changes, both required, both cheap, and they are one experiment:

* **A NON-TRIVIAL CHARGE SECTOR.** Move the physical sector from G_v = +1 everywhere to
  **G_v = −1 at two chosen vertices** (a static charge pair), then to dynamical Z_N matter. This is
  the one move that breaks the Perron-Frobenius inertness measured in §1: with a −1 sector the
  ground state is no longer the positive vector, so the Gauss law stops being a fact about which
  subspace the state happens to lie in and starts being a fact about which state is selected.
  Minimal arm-diff: the sign pattern in the Gauss constraints, and nothing else.
* **A SYSTEM REGION CARRYING A CYCLE**, |S| ≥ girth, so that its gauge-invariant algebra is
  non-abelian and has magnetic content. This is where the algebra choice is not a convention but a
  live disagreement of one full bit (§4), and therefore the only place the criterion has anything to
  decide.

**The carrier to run it on is `tri_chain12`** — V = 8, L = 12, cubic, girth 3, d = 5, so it carries a
three-link magnetic system region (the triangle on links 1,2,3) *and* nine environment links, enough
for four disjoint fragments, *and* it is the simple floor. With Z₂ matter the physical dimension is
2^(L+1) = **8192**, against 2^22 = 4194304 on Heawood: a 512× saving that removes the instrumental
blocker lane A declared, before any orbit-basis RDM is applied.

---

## 3. RULING: DOES THE CRITERION DISCRIMINATE?

**IT DISCRIMINATES BETWEEN STATES. IT DOES NOT DISCRIMINATE BETWEEN CARRIERS. AND ON THE CARRIERS
THE SWEEP RANKED, IT DOES NOT DISCRIMINATE AT ALL.**

* **On lane A's sweep: NO, and therefore the sweep measured nothing about record formation.** Three
  Haar-random physical states score identically to the ground state on tri_chain12; two more do so
  on heawood; a coupling at which the record is empty (H_elec(S) = 0.0082 bits) scores identically
  to one where it is nearly a full bit. Every arm reads exactly 1.000000. `OUT_rule_main.txt` block
  5, `OUT_rule_main2.txt` block 8a. **I say so plainly, as the brief requires: the sweep measured
  d = dist_{G−l}(tail,head), a graph metric, and nothing else.** Lane A's own §4 and CF-1 said this
  about itself; the ruling is that they were right and that the threshold number rests on it.
* **On lane B's discrimination test: YES, but half of it is an artefact.** Lane B's B-1 reports a
  7-vs-0 separation on theta_8 between broadcast and scrambled arms, and refuter B re-scored it
  under three definitions and passed it. Neither noticed that **one of the two broadcast arms flips
  sides when the algebra is made gauge-invariant**: the magnetic GHZ scores R_δ = 7 under the
  extended algebra and **R_δ = 0 under alg{X_l}**, joining the Haar arm
  (`OUT_rule_main2.txt` block 8c). Only the electric GHZ is a gauge-invariant broadcast. B-1's
  separation survives; its arm list does not.
* **On the partition: NO.** §0's table. Same state, same algebra, same fragment count — R_δ moves
  from 5 to a mean of 0.400 on heawood. A criterion whose answer is set by a choice made after the
  state is known is not discriminating between states.

**A null reads two ways, and I score neither as confirmation.** The identical Haar/ground result
reads either as "the criterion is degenerate here" or as "the Gauss law really does broadcast the
electric flux and every physical state inherits it". The second reading is *true* — it is the
theorem of §5 — and it is exactly why the first reading is also true. Both hold at once; that is the
character of a constraint-forced observable, and it is why no link count derived from it is a
threshold.

---

## 4. RULING: THE BOUNDARY-ALGEBRA DEPENDENCE — THIS IS THE FINDING

The brief: *if I(S:F) depends on the choice, that IS the residue W-04 located, reached for the first
time with a real dynamics, and it is the finding.*

**IT DEPENDS. AND THE DEPENDENCE IS NOT A NUISANCE TO BE CONVENTIONALISED AWAY — IT IS ZERO ON ONE
KIND OF RECORD AND TOTAL ON THE OTHER, AND WHICH KIND YOU GET IS SET BY THE COUPLING.**

The measurement holds the carrier, the state, the fragments and δ fixed and moves exactly one
object: the algebra assigned to the system link. EXT = the full 2×2 matrix algebra in the extended
Hilbert space (what lane A and lane B reading 1 both compute); CHI = the Holevo quantity of
alg{X_l}, the only non-trivial gauge-invariant subalgebra of a single link; CL = electric flux on
both sides. `OUT_rule_main3.txt`.

| record | carrier / state | EXT / H(S) | gauge-invariant / H | **EXT − gauge-invariant** |
|---|---|---|---|---|
| **ELECTRIC**, non-enclosing F | dbl_chain9, tri_chain12, petersen, heawood; ground and Haar; g² = 0.10…3.00 | 1.000000000 | 1.000000000 | **0.000000000 bits, exactly, on every fragment** |
| **ELECTRIC**, enclosing F | same carriers, same states | 2.000000000 | 1.000000000 | H(S) — the purification jump, which is gauge-variant |
| lane A's own C1 arm B | heawood, g² = 1.00, \|F\| = 14 fixed | 1.213848345 | 1.000000000 | 0.2138 — **lane A's only live control lives entirely in the gauge-variant part** |
| **MAGNETIC** | theta_6 and theta_8, magnetic GHZ = the exact g² → 0 ground state | 1.000000000 at every proper fragment | **0.000000000** at every proper fragment | **1.000000000 bits. Full swing.** |

Three consequences, each of which changes what gets built.

1. **REFUTER B IS RIGHT ABOUT THETA AND WRONG ABOUT LANE A.** Its adjudication —
   I(alg{Z₀}:alg{Z₁}) = 1.000000 and I(alg{X₀}:alg{X₁}) = 0.000000, so the reading-1 plateau is the
   redundancy of an operator no measurement in this theory can perform — I confirm exactly on
   theta_6 and theta_8. But it does **not** transfer to lane A's carriers. There the plateau sits at
   1.000000 in *both* channels, because on a graph with d ≥ 2 the Gauss surfaces carry the electric
   flux, which is gauge-invariant. **Lane A's plateau is an observable. Lane B's reading-1 plateau is
   not.** The two lanes' "thresholds" are thresholds for different objects and were never comparable.
2. **THE ALGEBRA CHOICE MOVES THE PLATEAU COUNT, AND THEREFORE THE NUMBER, A FOURTH TIME.** Under
   alg{X_l} the purification jump disappears (a classical bit's mutual information with anything is
   bounded by its own entropy), so the last fragment stays at 1.000000 and the count rises from
   d − 1 to d: dbl_chain9 4 → 5, tri_chain12 4 → 5, heawood 4 → 5, and **petersen 3 → 4**, i.e.
   MARGINAL becomes EXHIBITED at L = 15. Pushing this to the floor gives four points at d = 4,
   V ≥ 5, **L = 8** — measured: `mg_chain(4)`, V = 5, L = 8, d = 4, H_elec(S) = **0.738303126**,
   EXT gives 3 points and alg{X_l} gives **4** (`OUT_rule_main5.txt`). **The answers now on the table
   are 8, 9, 12, 15, 21 and "unreachable", all correct, all under different unstated conventions.**
   This is the frame objection, now with the algebra axis added to the fragment axis.
3. **THE RESIDUE IS PHASE-INDEXED.** Small g² is where the ground state is magnetic and where the
   algebra disagreement is one full bit; large g² is where it is electric and the disagreement is
   exactly zero. So the question "which boundary algebra" is not a free convention — it is
   *undetermined precisely in the phase the dynamics selects at small coupling*, and determined
   (trivially, because the answer is Gauss-forced) in the other. **That is W-04's residue, located,
   and located with a Hamiltonian and a coupling underneath it for the first time.** It is the only
   thing in this round that required attempt two's genuine dynamics to say.

One correction to refuter B's phrasing that matters for the next step: Perron-Frobenius shows the
Gauss law is inert **on the state**. It is not inert on the algebra — that is where it bites, and it
bites for one full bit. Adding static charges will fix the state-side inertness; it will not create
the algebra-side residue, which is already here and already measured.

---

## 5. WHAT SURVIVES THIS ROUND

* **The instrument.** Three independent constructions of the physical sector — lane A's spanning-tree
  gauge fixing, refuter A's orbit averaging in the full N^L space, and my min-weight-cycle-basis
  orbit labelling — agree, and all three reproduce the sealed T1 null to six decimals. Gauss residual
  0.000e+00 on every lifted state I built, up to L = 21.
* **THE GAUSS-SURFACE THEOREM** (lane A's F3, refuter A's most durable result, independently
  re-derived and re-measured here). For any physical state of Z_N pure gauge theory, any coupling:
  if F contains a u–v cut in G−l and S ∪ F contains no cycle through l, then I(S:F) = H(S) exactly.
  Rule-A plateau points = d − 1 under the extended algebra and d under alg{X_l}; R_δ = d. Zero
  deviations on every carrier either lane or I have measured. **This is the round's one durable
  fact, and it is the reason the round has no measurement: it makes the criterion a theorem.**
* **Lane B's B-7**, that the redundancy verdict is not fixed by the link partition alone. Strengthened
  three ways: by refuter B's enumeration (8 verdicts across the 16 gauge-invariant algebras of a
  plaquette), by my EXT-vs-alg{X_l} table (§4), and by my partition randomisation (§0).
* **Lane B's T1**, magnetic content = cycle space, hence a single link has no localised
  gauge-invariant holonomy and a magnetic record needs |S| ≥ girth. This is what forces §2c's second
  clause.
* **Lane B's B-2/B-3**, the two trivial-satisfaction traps: H(S) = 0 scores maximally, and R_δ is
  blind to how much is recorded. My g² = 3.00 arms make the second one concrete on lane A's own
  carriers: H_elec(S) = 0.001283123 bits on tri_chain12, R_δ = 5, plateau perfect.
* **Refuter B's Perron-Frobenius theorem**, extended here to min-degree-3 carriers.
* **Refuter A's corrected floor arithmetic**, L_min = ⌈3(P+2)/2⌉, and its exhaustive search, both
  independently reproduced.

## 6. WHAT FALLS

* **Lane A's F1 (status PROVED): 21 links.** False as stated; the simple floor is 12 and the
  multigraph floor 9. Its secondary floor ("3 points ⇒ V ≥ 10, L ≥ 15, Petersen") is also false.
* **Lane A's F8 (status PROVED)**, that girth 4 fixes d = 3 and only girth lengthens the plateau.
  Refuted by refuter A's sq_chain15; and girth does not even order d (tri_chain12: girth 3, d = 5;
  Petersen: girth 5, d = 4).
* **Lane A's claim that rule B independently confirms the plateau.** Rule B is a band-crossing count,
  not a flatness test, and its count moves with L at fixed d.
* **Lane B's reading-1 threshold of six links.** Gauge-variant in its entirety: EXT = 1.000000 at
  every proper fragment where the gauge-invariant channel reads 0.000000. Confirmed independently.
* **Lane B's reading-2 threshold of twelve links**, on refuter B's argument (a purity discard that
  does not apply to a restricted pointer algebra, and the false step "each fragment must contain a
  cycle"). I did not independently rebuild its six-link replacement and do not certify that number
  either.
* **Refuter B's generalisation that "the gauge structure enters nowhere".** True of lane B's
  estimator, false as a statement about the theory: the Gauss law is inert on the state and decisive
  on the algebra, by one full bit.
* **The whole family of link-count answers, as thresholds.** 8, 9, 12, 15, 21 and "unreachable" are
  all correct answers to the question as written. A quantity that a convention can move by a factor
  of two and a partition can move from 5 to 0 is not a threshold.

## 7. MY OWN CONFOUNDS, RECORDED RATHER THAN FIXED

* **CF-R1.** My CL channel restricts *both* sides to electric flux, so it is a lower bound on the
  full gauge-invariant I(S:F) — the gauge-invariant algebra of a multi-link fragment also contains
  that fragment's Wilson loops. CHI (Holevo, environment unrestricted) upper-bounds it. On every
  carrier of §2a the two agree to 1e-9; on Haar states on theta they do **not** (theta_8 seed 11,
  |F| = 4: CHI = 0.656208 against CL = 0.074777). **I did not fix which of CHI and CL is "the"
  gauge-invariant criterion.** That is a sixth undeclared convention and it is mine.
* **CF-R2.** My random partitions are equal-size round-robin over a shuffle; they do not search for
  the best partition. The claim is that R_δ is partition-relative, not that no good partition exists
  — the Gauss cuts are a good partition and always exist. The null reads two ways and I score
  neither.
* **CF-R3.** Z₂ only, single-link system regions only, in everything I computed. All statements
  about magnetic *system regions* are lane B's and refuter B's; I verified only the theta magnetic
  GHZ.
* **CF-R4.** My plaquette set is a min-weight cycle basis by matroid greedy over the full cycle
  space; ties are broken by enumeration order, so it is canonical only up to ties. For the T1
  reproduction I had to override it with T1's specific pair (W₀₁, W₁₂), which is lane A's CF-6 in a
  new place.
* **CF-R5.** Ceiling: full-vector partial traces at N^L ≤ 2^21, always taken on the smaller side.
  I did not go above it and I did not verify refuter A's above-ceiling runs. **The threshold in §2b
  does not sit beyond any ceiling — it is a structural statement, not a search result — so the
  ceiling does not bear on it.**
* **CF-R6, the largest.** Same lineage throughout. My agreement with a refuter is weak evidence.

---

## 8. NEXT STEP — MANDATORY, AND SPECIFIED

**ONE RUN. Z₂ GAUGE + STATIC CHARGE ON `tri_chain12`, WITH A MAGNETIC SYSTEM REGION, AND THE
CRITERION PINNED IN WRITING BEFORE ANY STATE IS BUILT.**

1. **THE CARRIER AND THE CUT, FIXED NOW.** `tri_chain12`: V = 8, L = 12, cubic, girth 3, d = 5,
   edges as listed in §2a. **S = the triangle {links 1, 2, 3}** — the smallest region on a simple
   min-degree-3 carrier whose gauge-invariant algebra is non-abelian and carries a Wilson loop.
   E = the remaining 9 links, partitioned into **four pairwise-disjoint fragments declared before the
   state is computed**, sizes 2,2,2,3. Physical dimension 2⁵ = 32 pure gauge, 2¹³ = 8192 with matter.
   Everything is exactly diagonalisable; there is no ceiling in this run.
2. **THE ISOLATED VARIABLE IS THE CHARGE SECTOR AND NOTHING ELSE.** Arm 1: G_v = +1 at all eight
   vertices. Arm 2: G_v = −1 at vertices 0 and 4 (a static charge pair straddling the bridge). Arms
   byte-identical apart from the sign pattern. Sweep g² across the crossover. **Diff the arms and
   print the diff** — the commonest fatal defect of record in this program is a control whose arms
   are identical.
3. **THE CRITERION, PINNED IN ADVANCE, ALL SIX AXES.** (i) system algebra: named, from the
   gauge-invariant lattice on S, with the pointer subalgebra understood as an instance of the choice
   and not a separate axis; (ii) fragment algebra: named, and CHI-vs-CL settled (my CF-R1);
   (iii) fragments pairwise disjoint; (iv) |F| ≤ |E|/2; (v) join or union for composite fragments;
   (vi) H(S) > 0 reported beside every R_δ. **And the partition declared before the state.**
4. **THE QUESTION IS INVERTED. Do not ask where the plateau appears — ask where it can fail.** The
   pre-registered falsifier: on the pre-declared partition, does the charged-sector ground state
   plateau where the vacuum-sector ground state does not, and does a Haar physical state on the same
   carrier, cut, algebra and partition fail? **A plateau that a Haar state also produces is not
   reported as a result.** That single rule would have voided lane A's headline before it was
   written, and it is the one thing this round most needed.
5. **BANK, AND DROP.** Bank refuter A's orbit-basis RDM with its recorded correction (it lifts |A|,
   not L), and fix `zn_gauge.py:242` (`T @ T.T` → `T @ T.conj().T`) before any complex Hamiltonian
   goes through it. Drop the subdivision-invariant normalised redundancy: R_δ = d identically on
   every reachable carrier in both lanes and here, so R_δ/d ≡ 1 and carries no information. That
   sub-route is closed by arithmetic, not by a run.

**AND THE LOG ENTRY THE CHARTER ASKS FOR.** What did not work: the carrier-size question. It does not
work not because the answer was hard to find but because the criterion it asks about is, in pure
gauge theory on a forest region, a theorem — the Gauss-surface identity — wearing the clothes of a
measurement. The invented next step is the one above: put a charge in, so that the Gauss law selects
a state instead of merely containing one; put a cycle in the system, so that the algebra choice has
something to disagree about; and declare the partition before the state, so that the criterion can
fail. The residue this round located — that I(S:F) differs by one full bit between boundary algebras
on exactly the magnetic states the dynamics produces at small coupling, and by exactly zero on the
electric ones — is not a defect to be conventionalised away. It is the first thing in this program
that a Hamiltonian and a coupling were needed to say, and it is where record formation, if it is
anywhere here, is hiding.

---

*Files: `/Users/bgm/MB Work/where-atoms-come-from/W19_RULING/` — `rule_verify.py`, `rule_main.py`,
`rule_main2.py`, `rule_main3.py`, `rule_main4.py`, `rule_main5.py`, `OUT_rule_main.txt`,
`OUT_rule_main2.txt`, `OUT_rule_main3.txt`, `OUT_rule_main4.txt`, `OUT_rule_main5.txt`,
`SEALS.sha256` (verifies 11/11).*
