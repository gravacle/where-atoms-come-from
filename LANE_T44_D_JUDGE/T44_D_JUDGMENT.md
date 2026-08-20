# T-44-D — THE JUDGE. Reachable classes of a Γ-constrained coupling; composition; the one Newton comparison.

Date: 2026-08-20. Judge lane for T-44 (opened by O-54: "which falloff exponents a Γ-constrained
coupling can REACH"). Inputs: three sealed, adversarially verified lanes —

- LANE_T44_A_CORNER  (D=2 exponent classes)   — 70/70 gates, verifier not-refuted (high), fresh-code recheck 24/24
- LANE_T44_B_WORLD   (D=3 exponent classes)   — 58/58 gates, verifier not-refuted (high), exact deepening M=1400→2800
- LANE_T44_C_COMPOSE (composition classes)    — 36/36 gates, verifier not-refuted (high), independent linear-solve 22/22

Judge's own verification (this lane, logged below): SEALS.sha256 re-checked OK in all three lanes;
gate summaries re-counted from the sealed OUTs (70 PASS 0 FAIL / 58 PASS 0 FAIL / 36/36 PASS).
Every number quoted below traces to a sealed lane artifact; the judge adds no new measurement.

Supersession note: this lane already holds the pre-rerun judgment (JUDGMENT.txt + D24_AUDIT.txt,
2026-08-20 13:53, retained untouched), issued when only lanes A and C existed — the D=3 probe had
overflowed (O-58 N1). It ruled the Newton comparison NEITHER on exactly one earned-vs-attributed
ground, and PRE-REGISTERED its own upgrade rule verbatim: *"A computed 1/d bracket earns MATCHES at
member level"* (JUDGMENT.txt, open piece N1), prescribing the exact method (Perron/sector route,
landmark 1/6, INVERSE-DISTANCE window beside the log and linear windows, Watson/Spitzer anchors as
comparison-only). LANE_T44_B_WORLD then executed that prescription to the letter and survived
adversarial verification with one caveat sharpened. This judgment applies the pre-registered rule to
the completed record; the register's T-44 entry (C-87/C-88) is amended accordingly. Everything else
in the prior judgment stands, including its open pieces N2 and N3.

---

## 1. THE REACHABLE-CLASS STATEMENT (C-87, completed)

**For a coupling between records at earned separation d — implemented as a weighted sum over
Γ-admissible strings, per-link amplitude μ (the coupling tier's one declared parameter), a weight-w
string contributing μ^w, computed EXACTLY on the venue's own lattice — the reachable classes are
exactly three, and each is located by computation:**

**(a) EXPONENTIAL decay for μ < μ_c.** Leading term N_min·μ^d exactly — the w_min = d confinement
stratum (C-80/O-54; Wegner 1971 / Wilson 1974), gated equal in weight AND count to the true
Γ-admissible-writer coset (2^22 enumeration) before any sum. Rates approach the Ornstein–Zernike
owner values; the exponential family degenerates monotonically toward μ_c (gated).

**(b) At μ_c exactly — the earned dimension's OWN lattice Green's-function class, and nothing else.**
μ_c = 1/deg is the venue's own computed number (1/4 on the D=2 corner venue, 1/6 on the D=3 world
venue, 1/2 on the D=1 chain), located three independent ways (Perron row-sum sandwich; exact sector
sandwich; exact-rational resolvent pole — singular AT μ_c, solvable beside it). The critical member
by earned dimension, one instrument, three disjoint declared windows:
- D=1: LINEAR kernel (doubling ratio in [9/5, 11/5]) — control venue;
- D=2: MARGINAL/LOG — G itself diverges (Pólya marginality); regularized kernel logarithmic, computed
  intervals containing the exact owner anchors 4/π, 4−8/π, increment (2/π)ln2. **No power law is
  reachable in D=2 and none was forced.**
- D=3: **POWER LAW, exponent 1 — the inverse-separation law.** G finite at criticality (the venue's
  own transience). Class call (INV window [2/5,3/5], all three rays) depth-stable at M=1400, M=2800,
  and under an independent continuum instrument; d·G(d) = 0.5147, 0.4869, 0.4794, 0.4779 at
  d = 2, 4, 8, 16, converging onto a constant — the computed statement "exponent 1" — with the
  coefficient bracket tightening onto the attributed owner value 3/(2π) = 0.477465 (M=2800 bracket
  [0.476369, 0.487321], owner inside).

**(c) DIVERGENT for μ > μ_c** — term-by-term growth, no mediated coupling in the venue limit.

**No further class appeared in the swept family, in any venue.** The critical class is the venue's
own: μ = 1/6 is critical on D=3 and subcritical on D=2; the D=1/D=2/D=3 critical kernels occupy
disjoint computed windows.

### Hypotheses the statement needs (each named, none silent)
- **H1 (string model, declared):** walks with backtracking are the admissible-string ensemble. The
  leading strata are gated equal to the true Γ-coset (weight and count, every pair, both tiers); the
  small-μ behavior agrees on Γ's actual admissible set; beyond that the ensemble is model-declared.
- **H2 (venue limit):** critical and supercritical labels are venue-limit statements of the walk
  model — a finite venue's coset sum is a polynomial, always finite.
- **H3 (declared windows and scope):** class labels are computed booleans against windows declared
  before measurement; power exclusion below critical carries the declared scope p ≤ 8.
- **H4 (μ swept, never fitted):** nothing in Γ yet places μ at μ_c — see Section 3.

### Owners for the standard parts (comparison-only, applied after each class was computed)
Lattice/random-walk Green's functions and radius of convergence: Spitzer 1964, Lawler. Recurrence /
transience by dimension: Pólya 1921. Potential-kernel anchors (4/π, 4−8/π, (2/π)ln2): Stöhr 1950,
Spitzer 1964, McCrea–Whipple 1940. G(0) anchor D=3: Watson 1939. Asymptote 3/(2π): Spitzer P26.1.
Subcritical rates: Ornstein–Zernike. Walk bijections: Feller I. Perron–Frobenius/Gershgorin;
Collatz–Wielandt certificate (Varga). Confinement character of w_min = d: Wegner/Wilson (standing
C-80/O-54 attribution). OURS: the Γ-priced coupling (μ^d leading term from the earned w_min = d),
dimension entering only as the venue, and the reachable-class statement under those constraints.

---

## 2. THE COMPOSITION FACTS (C-88) and what they mean for a source term under C-72's split

C-72's standing encoding split: **occupancy** (one-signed presence) and **orientation** (two-signed
values) are two structurally different record encodings. The composition lane ran two declared
admissibility rules beside each other: **R1** (transparent endpoints — sources are only endpoints)
and **R2** (opaque sources — walks through the other source are lost).

- **F1 — Superposition is a property of the ADMISSIBILITY RULE, not of the walk sum.** Under R1 the
  defect is identically zero at every order — a per-order integer identity. Under R2 the coupling is
  strictly subadditive (self-shadowing defect < 0, shrinking with separation).
- **F2 — The mediator kernel is strictly positive (gated), so C-72's sign structure survives
  mediation UNCHANGED — the mediator contributes magnitude only.** Orientation-encoded sources
  screen (exact mirror-pair zero; dipole partially screened, depth μ-dependent). Occupancy-encoded
  sources accumulate one-signed at every probe with nothing to cancel. Opacity reduces counts and
  never flips a sign: screening-by-sign and shadowing-by-opacity are computationally distinct.
- **F3 — Shadowing exists only under R2/λ and equals a LENGTHENING of earned separation** (shadow
  onset order = punctured-venue BFS distance; C-80's w_min = d survives obstacles). Under R1 there
  is no shadowing at all.

**For a source term: the profile a Newtonian source needs — exact superposition, one-signed
accumulation, no screening, no shadowing — is reached at EXACTLY ONE design point: R1 + occupancy.**
Orientation gives the electromagnetic-like two-signed screenable profile; R2 the absorptive one.
All facts reproduce in D=2 (venue-structural, magnitudes dimension-dependent) and hold at μ rows
above and below the landmarks (six rows total including adversarial fresh rows) — classes of the
certified-convergent regime as a whole, not artifacts of any μ choice.

**Which rule and encoding the record surface itself enforces is NOT decided by these lanes** — it is
named in Section 3 as part of the next piece.

---

## 3. FINAL COMPARISON — the only place Newton is named (D-1)

**Question:** is Newton's form the unique critical member in earned D = 3?

**VERDICT: MATCHES — at member/form level, at one computed design point** (the pre-rerun judge's
own pre-registered rule applied: the computed 1/d bracket now exists).

The grounds, all computed:
1. In earned D=3 the critical member of the reachable family IS the inverse-separation law: the
   class call is depth-stable on all rays, d·G(d) converges onto a constant, and the coefficient
   bracket tightens onto the attributed 3/(2π). It is the UNIQUE critical member — the family holds
   exactly one class at μ_c, and it is this one.
2. The composed source term reaches the full Newtonian profile (exact superposition, one-signed
   accumulation, unscreenable, unshadowed) at exactly one design point: R1 + occupancy.
3. Newton's form was assumed nowhere: no gravitational form in any construction lane; μ_c located
   by computation; dimension enters only as the venue; owner values confined to comparison-after-class.

So: **at (R1, occupancy, μ = μ_c, earned D = 3) the Γ-constrained mediated coupling is Newton's
form — 1/d, superposing, accumulating, unscreenable — and at no other point in the swept design
space.** The earlier NEITHER is superseded on its own stated ground: the D=3 member it named as
owner-attributed is now computed (and the D-24 objection — substituting an owner for a computation —
no longer applies). **O-58 piece N1 CLOSES.**

**What MATCHES does not say (scope, kept):**
- (i) Exponent 1 is an asymptotic statement pinned by coefficient convergence; the finite-pair
  exponent bracket at depth excludes exactly 1 with a finite-d drift ~ +0.04 (the sharpened T-44-B
  caveat, carried verbatim).
- (ii) The 1/d law is exact in an EMERGENT Euclidean norm over the earned L1 separation (√2/√3
  coefficient factors computed); which norm the surface's own physics selects is not settled.
- (iii) The verdict is conditional on OCCUPYING the design point: nothing in Γ yet forces μ to sit
  AT μ_c, and nothing yet forces the surface to enforce R1 + occupancy sourcing. MATCHES is a
  statement about the critical member's form and uniqueness — never that gravity is derived.

### What would EARN criticality — the masslessness condition (the named next piece)

μ_c = 1/deg is exactly the per-link normalization at which the writer measure is conserved — the
venue's row sums equal deg, so μ·deg = 1 is the Perron normalization and the resolvent's pole.
**Labeled hypothesis (carried from the register's T-44 entry, still unproved): if Γ's writer kernel
conserves measure — per-link amplitude redistributed, never created or lost — then μ sits at μ_c
structurally, and masslessness is measure conservation, not tuning.** The falsifiable next
increment: compute, from the record surface's OWN writer ensemble (not a declared μ), whether the
induced per-link amplitude is the measure-conserving one. Beside it, the two companion pieces:
whether the surface enforces R1 + occupancy (the design point), and which norm the surface selects
(emergent Euclidean over earned L1). These are O-58's remaining pieces N2 and N3, plus the norm
question this judgment adds.

---

## 4. INCREMENT VERDICT per C-77: **SATISFIED**

C-77's increment standard: derive, from record-surface boundary structure alone — no classical form
assumed, every concept earning its place per D-24 — a scale-level relation that matches a known
gravitational one; the claim fails if every such derivation requires assuming the classical form or
yields relations incompatible with the known ones.

**Delivered by T-44:** from the earned w_min = d (C-78/C-80, boundary-crossing cost) and one
declared mediator parameter, the reachable-class family was computed exactly; its unique critical
member in the earned dimension (C-79's computed 3) is the 1/d potential — Newton's form — reached
without assuming it, with the Newtonian composition profile located at exactly one computed design
point. The derivation did not require the classical form, and the relation it yields is compatible
with the known one. **The third increment is SATISFIED — at family-and-form level, with the
conditional (occupancy of the critical point) carried in the verdict sentence itself.**

Not closed by this verdict (per the standing discipline, the route stays open with its next step
named): what earns criticality — masslessness as measure conservation (O-58 N2), surface enforcement
of the design point (O-58 N3), and norm selection. Disposition of the register amendment
(upgrading C-87's D=3 member from owner-attributed to computed, and re-issuing the comparison
verdict MATCHES over NEITHER) awaits the principal.

---

## 5. D-24 AUDIT (judge lane)

- **The judge adds no new measurement.** Own verification limited to: SEALS.sha256 re-check in all
  three lanes (OK), gate re-counts from sealed OUTs (70/0, 58/0, 36/36), verifier artifacts
  inspected. Logged in JUDGE_CHECKS.txt.
- **Every number quoted traces to a sealed artifact**; no number re-derived or adjusted here.
- **MATCHES is issued at form level only**, against computed brackets and depth-stable class calls;
  owner values (3/(2π), Watson, Spitzer anchors) appear as comparison anchors only, never on any
  measurement path.
- **The conditional structure is load-bearing and kept in the verdict sentence**: unique critical
  member (computed) + occupancy of criticality (NOT earned) + design-point enforcement (NOT earned).
- **D-1 kept program-wide:** Newton absent from all three construction lanes; named only in
  Section 3 of this judgment.
- **Supersession is grounded, not editorial:** the earlier NEITHER named its own missing piece (the
  D=3 computation); that piece now exists and survived adversarial verification; no other ground of
  the earlier entry is touched.
- **INSERTED and labeled (inherited from lanes):** μ rows, the λ opacity dial, the R1/R2 rule pair,
  the walk ensemble (H1). **EARNED (inherited):** d (C-78), dimension (C-79), w_min = d (C-80),
  μ_c (located three ways), every class label (computed booleans, pre-declared windows).
- **Relevance test:** the borrowed Green's-function machinery was applied to the named variable
  G_μ(d) — the Γ-priced admissible-writer sum — and to nothing else.

## Files
- /Users/bgm/MB Work/where-atoms-come-from/LANE_T44_D_JUDGE/T44_D_JUDGMENT.md (this file — the completed-record judgment)
- /Users/bgm/MB Work/where-atoms-come-from/LANE_T44_D_JUDGE/JUDGE_CHECKS.txt (this run's verification log)
- /Users/bgm/MB Work/where-atoms-come-from/LANE_T44_D_JUDGE/JUDGMENT.txt (pre-rerun judgment, retained — carries the pre-registered upgrade rule)
- /Users/bgm/MB Work/where-atoms-come-from/LANE_T44_D_JUDGE/D24_AUDIT.txt (pre-rerun run's audit, retained)
- /Users/bgm/MB Work/where-atoms-come-from/LANE_T44_D_JUDGE/SEALS.sha256 (covers all of the above)
