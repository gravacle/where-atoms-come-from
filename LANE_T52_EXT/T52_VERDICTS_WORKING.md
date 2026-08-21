# T-52 PER-ROW WORKING — carrier verdicts for the three gap rows
Date: 2026-08-21. Method FIXED per LANE_T9_AUDIT/METHOD.md — applied, not improved.
Read-only on the repo; THE REGISTRAR appends and reseals. Companion: GAP_ENUMERATION.md,
gap_rows.tsv, REGISTRAR_DRAFTS.md, append_rows.tsv (ready-to-append lines).

The method's text I audit against (METHOD.md, quoted):
  * "STRUCTURALLY DIFFERENT meaning different algebra/mechanism, never the same model at
    another size or parameter set" (lines 5-6).
  * "Probes defaulted to SINGLE-CARRIER; every TWO-CARRIER claim was attacked by an
    adversarial refuter checking the cited evidence exists, concerns the row's actual
    result, and the carriers genuinely differ" (lines 6-8).
  * NOT-CARRIER-SHAPED = "rows whose content is not a computation on a carrier —
    meta-classifications, literature statements, definitional consequences; used
    sparingly" (lines 13-14).
Calibration rows from the sealed TSV: C-88 (SINGLE — the 2-D lattice a spot-check
"dimension discriminator", demoted because the full gate set ran ONLY on the 3-D venue);
C-87 (TWO-CARRIER toric-2x2; spin-chain — cross-dimension venues ACCEPTED because the
row's full result ran on each, three instruments, disjoint windows); C-78 (four carriers
listed); A-EM/A-AL (role rows with their own measured instances audited carrier-shaped).

---

## A-PR — VERDICT: NOT-CARRIER-SHAPED

KIND: T-15 restatement row. No computation of its own; every clause's carrier support is
inherited through cited rows (per-clause table, LANE_T15_ROLES/A_PR_RESTATED.md), each
citation already carrying its own sealed T-9 verdict.

WHY NOT-CARRIER-SHAPED AND NOT AN INHERITED TWO-CARRIER:
1. Method fit. The row's content "is not a computation on a carrier" — METHOD's own
   definition, verbatim. Nothing of A-PR's was ever run on any carrier; what carriers
   carry are the CITED rows' results, already audited under their own IDs. This is the
   same shape as C-85 ("a reproduction/validation claim ... not a carrier computation")
   and G-13 in the sealed table — a 15th meta row, within "sparingly".
2. The refuter standard blocks inheritance. TWO-CARRIER requires cited evidence that
   THIS ROW'S actual result was established on two carriers. A-PR's actual result is the
   assembled three-roles statement; the assembly was established by a citation table and
   a row-fidelity verifier (A_PR_RESTATED.md, VERIFICATION), not on carriers. Only its
   parts are carrier-borne, and they are audited elsewhere.
3. An inherited mark rots silently when a citation moves — and one HAS moved: C-72 is
   now PARTIAL. A verdict that re-states the citations' verdicts would double-count
   sealed evidence and go stale without any audit action. NOT-CARRIER-SHAPED with the
   inheritance recorded in the evidence cell keeps T-16's grep honest (the support is
   visible, not claimed as the row's own).
4. The SINGLE-CARRIER default does not fit either: the default is for carrier-shaped
   rows, and no single carrier could be named without arbitrariness (the citation union
   spans 8: toric-2x2, bouquet, DD4, steane-713, spin-chain, rank2-abelian,
   macrospin-CoCrPt, NAND-floating-gate).
Note the precedent boundary: A-EM/A-AL are role rows WITH their own measured instances
(f10c Paulis, fitted exponents) — carrier-shaped. A-PR has no instance of its own.

SIDE FINDING (evidence for the registrar, not part of the verdict): A-PR's ledger
bracket and A_PR_RESTATED.md assert "no load-bearing citation WITHDRAWN/FAILED/PARTIAL"
— clean at landing (VERIFICATION section, round 1) — but C-72, load-bearing in the
EM-world clause, is now PARTIAL. C-71 (still PROVED, TWO-CARRIER) alone carries that
clause's two-carrier support, but the row's own gate sentence is no longer true of the
current ledger. Registrar follow-up needed; separate from the carrier mark.

POINTERS: ledger/status_ledger.tsv row A-PR; REGISTER_V001.md:10026-10047 (T-15);
LANE_T15_ROLES/A_PR_RESTATED.md (per-clause table + VERIFICATION);
LANE_T15_ROLES/T15_draft_and_citations.json.

---

## C-90 — VERDICT: SINGLE-CARRIER (macrospin-CoCrPt)

KIND: computed result (T-44 world lane): D=3 critical row, mu_c = 1/6 exact, power-law
class with exponent bracket containing 1.

THE ONE VENUE: LANE_T44_B_WORLD/PUBLISHED_CONVENTIONS.txt (MODEL paragraph): the walk
sum is computed "on the venue's own lattice -- the grain-adjacency graph of the census
access geometry (GR1 grains, adjacency = shares a face)". That is the T-29/T-43-B HDD
access model — macrospin-CoCrPt census lineage, the carrier name the sealed audit
already uses for exactly this venue (C-88's row). Every one of the row's claims —
mu_c = 1/6 three ways, the power-law class in the INV window, G finite with Watson
inside, the coefficient bracket onto 3/(2pi) — lives on this lattice and nowhere else.

CANDIDATE SECOND CARRIERS, EACH REFUSED UNDER THE METHOD:
1. The cross-dimension discriminator (FINDINGS item 4, labeled D-15: D=1 linear, D=2
   log, D=3 power, disjoint windows). The D=1/D=2 venues carry DIFFERENT results (1/2
   linear, 1/4 log — those are C-87's results, already audited TWO-CARRIER under C-87's
   own ID); they discriminate dimensions, they do not re-establish the D=3 result. This
   is the C-88 demotion precedent verbatim — "a weakened test, not a second carrier."
   The contrast with C-87 is exact: there each venue carried THAT row's full result;
   here the other venues carry a different row's result.
2. The adversarial verifier's independent Z^3 DP, kernel deepening 1400->2800, m=1000
   divergence witness, forbidden continuum Fourier instrument — independent MACHINERY
   on the SAME venue. The audit counts carriers, not instruments.
3. Owner constants (Watson G(0), Spitzer 3/(2pi), OZ rates) — comparison-only AFTER the
   class, by the lane's own declaration ("Owner anchors AFTER the class").
No other candidate exists in the lane or the register entry. SINGLE-CARRIER stands —
which is what the T-52 plan row itself anticipated ("the chain's head stands on one
venue"). The mark makes the debt printable; it does not discharge it.

POINTERS: ledger/status_ledger.tsv row C-90; REGISTER_V001.md:9961-9999 (T-44 LANDS);
LANE_T44_B_WORLD/PUBLISHED_CONVENTIONS.txt (MODEL paragraph; FINDINGS 2 and 4; NEXT
STEP paragraph "1/6 world, 1/4 corner, 1/2 chain"); LANE_T44_B_WORLD/VERIFY.

---

## C-91 — VERDICT: TWO-CARRIER (toric-2x2; macrospin-CoCrPt)

KIND: computed result (T-48, three lanes + three adversarial verifiers + judge):
masslessness is measure conservation; mass gap in closed form.

THE TIERS QUESTION, ANSWERED FROM THE METHOD'S TEXT: tiers of one construction are NOT
carriers — a tier is instrument stratification, and independent machinery on one venue
is machinery (the C-90 refusal above). The TWO-CARRIER claim here does NOT rest on
"kernel/corner/world are three tiers"; it rests on the fact that two of the tiers ran
the row's actual result, full-suite, on carrier builds that differ in ALGEBRA/MECHANISM
— METHOD's own definition of structurally different — not in size or parameter set.

CARRIER 1 — toric-2x2 family, the stabiliser writer algebra (corner tier,
LANE_T48_B_CORNER, 112/112, own adversarial verifier). Venue: plaquette-adjacency
graphs of the toric carriers (4,6) and (3,7), built from plaquette supports alone,
writers the carrier's own single-edge X operators gated against the O-54-C admissible
coset. The row's result established here: link-uniformity EARNED from the elementary
writer's invariant tuple (1,2,2,2,1) identical on every link (PUBLISHED_CONVENTIONS
"WRITER-ALGEBRA ROUTE (earned)"); conservation deg*t = 1 has unique solution
t* = 1/deg == mu_c EXACTLY, re-located in-lane; conserving anisotropic members exactly
critical; beta-biased ensembles land computably off conservation AND off criticality
(open gap, nonsingular resolvent) — both directions of the identity.

CARRIER 2 — macrospin-CoCrPt lineage, written media (world tier, LANE_T48_C_WORLD,
42 gates, own adversarial verifier). Venue: the census access geometry's deg-6 grain
lattice; writers are the (iv') energy-conserving dilation writers on written media with
the activation convention of model/project_model.py (u = exp(-E_b/kT), v = u*b) —
t48c_world.py lines 8-31. The row's result established here: E1 (iv' literal transport)
conserving and critical at EVERY dE and barrier, m = 1/6 = mu_c, criticality structural;
E3 (the model's own erase channel) NEVER critical, mass gap in closed form
ln(mu_c/mu) = ln(1 + e^{dE/kT}/l) with f0 and E_b dropping out — both directions again.

WHY THIS SURVIVES THE REFUTER STANDARD:
1. Evidence exists: two sealed lanes, each with its own full gate suite and its own
   independent-machinery adversarial verifier, judge-checked (REGISTER 10049-10075).
2. It concerns the row's ACTUAL result: both carriers compute the SAME registered
   identity — the induced per-link amplitude of the measure-conserving writer ensemble
   equals mu_c = 1/deg identically (each venue its own number: 1/4, 1/6), and
   non-conserving ensembles are computably massive. Not different sub-claims on
   different builds: the identity and its converse land on each carrier separately.
3. The carriers genuinely differ — different algebra/mechanism: a stabiliser syndrome
   algebra with Gamma-priced X writers versus thermally-activated Boltzmann dilation
   writers on written media. This is the same carrier-class distance the sealed audit
   already recognizes (C-78 lists toric-2x2 and macrospin-CoCrPt as distinct carriers
   of one row). It is NOT a size/parameter change of one model.
4. The C-88 demotion precedent does NOT apply: neither venue is a spot-check — each ran
   the complete suite (112/112; 42 gates), each was separately adversarially verified.
   The acceptance precedent is C-87 (TWO-CARRIER, toric-2x2; spin-chain — venues of
   different dimension count when the row's full result runs on each).

NOT COUNTED (self-refuted before the refuter gets there):
  * Toric (4,6) vs (3,7): same model at another size — never counts (METHOD line 6).
  * The world tier's b-sweep (b = 1 CoCrPt orientation / b < 1 NAND occupancy): the
    encodings enter as PARAMETER SETS of one two-state activation kernel on one lattice
    build — so NAND-floating-gate is NOT claimed as a carrier of this row.
  * Kernel tier venues (C8 Ising-ring bond graph, T3/T4 toric, Z27 grain; 126/126) and
    the corner lane's chain C_24: supporting evidence for the abstract identity, but C8
    is chain-lineage in a lane where the chain also serves as the D=1 discriminator, so
    spin-chain is not claimed — the two carriers above suffice and are unweakened.
  * The square (5,5) venue: symmetry control only, by the lane's own label.

POINTERS: ledger/status_ledger.tsv row C-91; REGISTER_V001.md:10049-10075 (T-48 LANDS);
LANE_T48_A_DERIVATION/PUBLISHED_CONVENTIONS.txt (VENUES block: deg computed never
declared); LANE_T48_B_CORNER/PUBLISHED_CONVENTIONS.txt (VENUE block; GUARD; K1; D-15
controls); LANE_T48_C_WORLD/t48c_world.py lines 1-64 (venue, writers, E1/E2/E3,
declared samples); each lane's VERIFY subdir; LANE_T48_D_JUDGE.

---

## CONSEQUENCE LINE FOR THE CHAIN'S HEAD (for the register entry)
C-88 SINGLE-CARRIER, C-90 SINGLE-CARRIER, C-91 TWO-CARRIER: the measure-conservation
identity now stands on two structurally different carriers, but the D=3 critical row
itself (C-90) still stands on the one census venue. The single-venue debt at the head
is now PRINTED, not discharged; a second structurally different D=3 world venue for the
C-90 result is the named open debt.
