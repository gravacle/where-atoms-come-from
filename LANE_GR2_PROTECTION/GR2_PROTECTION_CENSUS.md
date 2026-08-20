# LANE_GR2_PROTECTION — IS CLAUSE (v) FALSE OF THE WORLD'S RECORDS?

Lane opened 2026-08-20. Sibling lane LANE_GR1_CENSUS was empty at open time (running in
parallel), so this lane fields its own ten-record census; GR1's list should be reconciled
against this one at merge.

Companion arithmetic: `gr2_arrhenius_numbers.py` / `.OUT.txt` (all derived numbers below
recomputed there from CODATA constants; literature inputs cited inline).

Clause (v) under test, as written:
  (v) no admissible operation on a SINGLE CONTRACTIBLE REGION has U†RU = −R,
      ADMISSIBLE U := [U,H] = 0 (O-4, PROVISIONAL).
Program context that binds this lane: O-4/C-60 (admissible = does no work), O-5 (clause (ii)
relaxes to a width bound, tolerance = inverse lifetime), O-12 (admissibility was what saved
(v) on the toric code), T-11-resolved/D-4-extended (clauses are permissions; carriers satisfy
them or not), O-16 (ε^d protects against static perturbation only; dissipation kills at γ,
no d anywhere in the exponent), O-49/O-50 (toric code canonical for (v)).

Framed per D-4: the question is not "is clause (v) true" but WHICH CARRIERS IN THE WORLD
SATISFY IT, and whether the mechanism that actually keeps the world's records alive is the
mechanism clause (v) names.

---

## 0. THE STRUCTURAL FACT THAT DECIDES PART 1 BEFORE ANY NUMBERS

Clause (v) is satisfiable only by records whose two states are locally indistinguishable —
this is a small theorem, not an intuition:

  THEOREM (support-compression). Let the two record states |0⟩, |1⟩ have identical reduced
  density matrices outside some contractible region A (pure global states). Then there exists
  a unitary supported ON A mapping |0⟩ to |1⟩. (Uhlmann/purification: equal complements ⟹
  the two states differ by a unitary on A. Standard.)

  COROLLARY. If a record's two states differ only inside a contractible region A — which is
  what "the bit is stored in the grain / the gate / the base pair" means — then a flipper
  supported on the single contractible region A EXISTS. If additionally the two states are
  degenerate under H (the precondition for ANY admissible flipper to exist, i.e. for clause
  (iv) to hold at all, by C-60's dE = 0), the flipper on A can be taken admissible on the
  record subspace. Hence for contractible-support records: (iv) ⟹ ¬(v).
  The pair {(iv), (v)} with one admissible class is UNSATISFIABLE off non-contractible
  supports. The toric code escapes precisely because its logical states have equal reduced
  states on EVERY contractible region (TQO-1) — the record's support is a non-contractible
  cycle no contractible region contains.

Every record in the census below stores its bit in a simply-connected lump and is READ by a
local probe (read head, sense amp, polymerase, eye). Local readability = the two states have
different reduced states on a small region = exactly what TQO forbids. The world's records
are locally readable BY DESIGN; clause-(v) records are locally unreadable BY CONSTRUCTION
(the register already knows the model's readout takes a traversal, W-37). That is the
mismatch in one sentence, before any energy scale is quoted.

The T-11 ruling does not rescue this. T-11's discriminator was region/system FRACTION (a
region ≥ half the carrier is not local; 22% with d=3 admitted no flipper). For a hard-drive
grain the flipping region is the grain: 4×10⁻²⁵ m³ against a ~10⁻⁸ m³ platter — fraction
~10⁻¹⁷ — and it flips the bit with certainty. Real records fail the fraction criterion by
seventeen orders of magnitude, not marginally.

---

## 1. THE CENSUS — TEN REAL RECORDS, THE LOCAL OPERATION THAT KILLS EACH, AND THE SUPPORT COMPARISON

Constants: kT(300 K) = 4.14×10⁻²¹ J = 25.9 meV. Landauer floor kT ln2 = 2.87×10⁻²¹ J.
Arrhenius throughout: τ = τ₀ exp(E_b/kT), attempt time τ₀ ~ 10⁻⁹–10⁻¹³ s (stated per row).

### R1. Hard-disk bit (CoCrPt-oxide perpendicular media grain)
- Support: grain ~8 nm dia × 8 nm, V ≈ 4×10⁻²⁵ m³; a user bit = majority over ~10–20 grains.
- Barrier: K_u ≈ 3×10⁵ J/m³ (CoCrPt alloy; literature range 2–6×10⁵) ⟹ K_uV ≈ 0.75 eV ≈
  29 kT per grain; industry design rule Δ = K_uV/kT ≥ 40–60 (Weller & Moser, IEEE Trans.
  Magn. 35, 4423 (1999)). τ₀ ≈ 10⁻⁹ s: Δ=40 → τ ≈ 7.5 yr; Δ=60 → τ ≈ 3.6×10⁹ yr.
- Local killer: write-head field ~1 T over ~50 nm (Zeeman 2μB ≈ 2.5 eV > barrier, computed
  from Ms = 5×10⁵ A/m); or HAMR laser spot ~30 nm at ~700 K (kills Δ thermally). Extent =
  the bit's own support. Write dissipation ~10 pJ ≈ 3×10⁹ × Landauer.
- Support comparison: the operation acts on the grain; the grain is contractible; region =
  record support; fraction of carrier ~10⁻¹⁷. Clause (v)'s comparison FAILS decisively.
- Mechanism: ENERGETIC (Neel-Arrhenius). Note the barrier is EXTENSIVE — K_u·V grows with
  grain volume. Nature buys protection with volume, i.e. with redundancy of spins.

### R2. NAND flash cell (floating gate / charge trap)
- Support: cell ~20 nm; stored charge ~10–100 electrons behind 6–8 nm tunnel oxide.
- Barrier: Si/SiO₂ conduction-band offset 3.1 eV = 120 kT (standard). Intrinsic thermionic
  escape τ ~ 10⁻¹³·e¹²⁰ s ~ 10³⁹ s — irrelevant; REAL retention (~10 yr at 55 °C, JEDEC) is
  set by trap-assisted TUNNELING through oxide defects, modeled with Ea ≈ 1.1 eV. Kinetic-
  plus-tunneling leakage, still barrier physics; part of the protection is barrier WIDTH
  (tunneling suppression), the quantum cousin of Arrhenius.
- Local killer: Fowler-Nordheim tunneling at ~10 MV/cm (~18 V on the control gate), extent =
  the cell; or heat (retention bake). Erase energy ~fJ–pJ ≈ 10⁶–10⁹ × Landauer.
- Support comparison: operation on the cell itself. (v) FAILS.
- Mechanism: ENERGETIC/KINETIC (barrier height + width).

### R3. DNA base pair (the genetic record)
- Support: ~1 nm, one base; 2 bits per position.
- Barrier: covalent backbone ~3.6 eV; the WEAK channel is chemistry: cytosine deamination
  Ea ≈ 1.1–1.25 eV, t½ ≈ 200 yr per base in ssDNA, ~140× slower in duplex (~3×10⁴ yr)
  (Lindahl, Nature 362, 709 (1993)); depurination ~10⁴ events per human cell per day.
- Local killer: one 4.9 eV UV photon (254 nm) forms a pyrimidine dimer on a ~1 nm site;
  a hydroxyl radical; a polymerase misincorporation. Extent = the base itself.
- Support comparison: (v) FAILS at the single-base scale.
- Mechanism: ENERGETIC — plus, uniquely in this census, an ACTIVE AGENT: repair machinery
  (proofreading + mismatch repair) takes the per-generation error rate from ~10⁻⁵ raw
  polymerase to ~10⁻⁹–10⁻¹⁰. Where the barrier alone is too low for the required lifetime,
  nature installs an error-correcting agent — precisely O-16's missing agent, existing in
  biology and absent in the model. Redundancy also: the record is replicated across ~10¹³
  cells; no single contractible region reaches the ORGANISM's genome record.

### R4. Ink on paper (carbon ink manuscript)
- Support: mark ~100 μm; pigment = graphitic carbon black.
- Barrier: carbon oxidation is negligible at 300 K (graphitic C–C ~4.9 eV; effective
  combustion activation ~1.5–2 eV needs ~500 K) — why 2000-year carbon-ink texts survive.
  The SUBSTRATE is the weak link: cellulose depolymerization Ea ≈ 1.0–1.25 eV (paper-aging
  kinetics literature), centuries at 300 K, dry.
- Local killer: mechanical abrasion (~mJ over the mark — razor scrape; palimpsests are the
  historical record of exactly this local operation); locally applied bleach/oxidant.
- Support comparison: the scrape covers the mark; (v) FAILS.
- Mechanism: ENERGETIC/KINETIC (chemical barriers + absence of agents: dryness, dark).

### R5. Photographic silver grain (developed film / daguerreotype)
- Support: filamentary Ag grain ~0.2–2 μm.
- Barrier: Ag is a noble metal; oxidation/sulfidation Ea ~ 1 eV scale (tarnish is the slow
  channel); archival film rated centuries.
- Local killer: ferricyanide bleach on the grain (standard photographic bleach chemistry);
  a ~μJ focused pulse melts it (Ag melting point 1235 K).
- Support comparison: (v) FAILS on the grain.
- Mechanism: ENERGETIC.

### R6. Optical disc mark (pressed pit / CD-RW phase-change)
- Support: pit ~800×400×100 nm in polycarbonate, or amorphous GST mark ~500 nm.
- Barrier: GST (Ge₂Sb₂Te₅) crystallization Ea ≈ 2.3 eV = 89 kT ⟹ τ(300 K, τ₀=10⁻¹³ s) ~
  10¹⁸ yr intrinsic (real archival life ~30–100 yr, set by dye/reflector chemistry at
  ~0.8–1.2 eV, not by GST); polycarbonate Tg ≈ 420 K.
- Local killer: the drive's own write laser, ~10 mW × 50 ns = 0.5 nJ on the mark (melts GST
  at ~900 K); or a scratch (~mJ, mechanical).
- Support comparison: (v) FAILS on the mark.
- Mechanism: ENERGETIC (deep well; the record IS a frozen non-equilibrium state).

### R7. STT-MRAM cell (CoFeB magnetic tunnel junction)
- Support: free layer ~40 nm dia × 1.5 nm.
- Barrier: design rule Δ ≥ 60 for 10-yr retention at ppm error rates (industry standard;
  same Neel-Arrhenius as R1) ⟹ E_b ≈ 1.55 eV.
- Local killer: spin-transfer torque current ~100 μA × 10 ns through the junction itself,
  ~0.5 pJ ≈ 1.7×10⁸ × Landauer.
- Support comparison: the write current flows through the record's own pillar. (v) FAILS.
- Mechanism: ENERGETIC.

### R8. Tree ring (dendrochronological archive)
- Support: annual layer ~0.1–5 mm through the trunk cross-section.
- Barrier: cellulose/lignin covalent bonds 3–4 eV; effective decay channels are catalyzed
  (fungal enzymes drop effective barriers to ~0.5–1 eV; combustion self-sustains above
  ~500 K, releasing 17 MJ/kg).
- Local killer: rot at the ring (requires water + microbes = supplying the catalyst), or
  fire. Extent = the ring's own material.
- Support comparison: (v) FAILS.
- Mechanism: KINETIC in the purest sense — the barrier is fixed; protection is managed by
  excluding the AGENTS (dryness ⟹ bristlecone ring records ~5×10³ yr; waterlogged anoxic
  ⟹ ~10⁴ yr).

### R9. Lunar crater (geological record)
- Support: bowl ~1 km, depth ~200 m.
- Barrier: to erase it locally you must resurface it: a comparable impact needs a ~50 m
  impactor at 17 km/s, KE ≈ 3×10¹⁶ J (π-group crater scaling, Melosh, order of magnitude);
  or ~km³ of lava. The "barrier" is astronomically high in absolute joules and ZERO in
  topology — any sufficiently energetic local event on the crater's own footprint erases it.
- Lifetime: ~10⁹ yr on the Moon (no atmosphere, no water — no agents); the SAME record shape
  on Earth lasts ~10⁷ yr (erosion agents present). The record is identical; the lifetime is
  set entirely by the agent census of the environment. Pure KINETIC protection.
- Support comparison: (v) FAILS.

### R10. Superconducting flux quantum (persistent current in a closed loop)
  — included as the world's BEST candidate for clause (v); see Part 4. Even here:
- Support: the winding number n of the order-parameter phase around the loop — genuinely
  non-contractible support. The closest thing nature offers to the program's Z̄.
- BUT the killer is local anyway: a PHASE SLIP is an event at ONE cross-section of the wire
  (a coherence-length segment goes normal, the phase unwinds by 2π there). Barrier =
  condensation energy of that segment, ΔF ≈ (B_c²/2μ₀)·A·ξ: for a 100×100 nm Al wire,
  ΔF ≈ 0.27 eV ⟹ Δ ≈ 760 at 4.2 K, Δ ≈ 3200 at 1 K. Measured bound: persistent current in
  a superconducting solenoid, decay time > ~10⁵ yr (File & Mills, PRL 10, 93 (1963)).
- Local killer in practice: a heat pulse of ~ΔF (fractions of an eV!) on one wire segment,
  or a persistent-switch heater — extent ~ξ³ ~10⁻²¹ m³, a vanishing fraction of the loop.
- Support comparison: the RECORD's support is non-contractible (the whole loop), but the
  FLIPPER's support is a single contractible segment. Topology reroutes the flip through a
  barrier; it does not remove the local flipper. Under clause (v) as written the phase-slip
  unitary is inadmissible ([U,H] ≠ 0 — it does work against the condensate), so (v) "holds"
  — but only in the same vacuous way it holds for a hard-drive grain with split wells: by
  the admissibility filter, not by topology, and clause (iv) fails by the same filter.
- Mechanism: ENERGETIC, topology-ASSISTED (topology localizes and raises the barrier; the
  suppression is still exp(−ΔF/kT), finite, Arrhenius).

### Census table

| # | record | support | flip operation | flip extent vs support | flip energy | barrier E_b | Δ = E_b/kT | lifetime | mechanism |
|---|---|---|---|---|---|---|---|---|---|
| R1 | HDD grain | 8 nm | head field 1 T / HAMR | = support | 2.5 eV Zeeman, 10 pJ dissipated | 1.0–1.6 eV | 40–60 | 7.5 yr – 3.6 Gyr | energetic |
| R2 | flash gate | 20 nm | FN tunneling 18 V | = support | ~fJ–pJ | 3.1 eV (height), 1.1 eV eff. | 120 / ~43 | ~10 yr @ 55 °C | energetic+tunneling |
| R3 | DNA base | 1 nm | UV photon / hydrolysis | = support | 4.9 eV (one photon) | 1.1–1.25 eV | ~46 | 200 yr (ss), 30 kyr (ds) | energetic + AGENT |
| R4 | ink mark | 100 μm | scrape / bleach | = support | ~mJ | 1.0–1.25 eV (substrate) | ~45 | centuries | kinetic |
| R5 | Ag grain | 1 μm | bleach / melt | = support | ~μJ | ~1 eV | ~40 | centuries | energetic |
| R6 | GST mark | 500 nm | own write laser | = support | 0.5 nJ | 2.3 eV | 89 | 10¹⁸ yr intrinsic | energetic |
| R7 | MRAM MTJ | 40 nm | STT current | = support | 0.5 pJ | 1.55 eV | 60 | 10 yr by design | energetic |
| R8 | tree ring | mm | rot / fire | = support | self-sustaining | 3–4 eV, cat. ~0.5–1 | agent-set | 10³–10⁴ yr | kinetic |
| R9 | crater | km | later impact | = support | 3×10¹⁶ J | — (absolute J) | agent-set | 10⁹ yr (Moon) | kinetic |
| R10 | flux quantum | loop (non-contractible!) | phase slip at ONE segment | ≪ support (ξ segment) | ~0.3 eV heat | ΔF ≈ 0.27 eV | 760 @ 4 K | >10⁵ yr measured | energetic, topology-assisted |

In every row the destroying operation is supported on a single contractible region, and in
nine of ten rows that region is the record's own (contractible) support. In the tenth (flux
quantum) the region is far SMALLER than the record's non-contractible support — the strongest
possible violation of the intuition behind (v): even nature's one non-contractible record is
flipped by a local operation, paying only a finite barrier.

---

## 2. THE TWO KINDS OF PROTECTION, AND WHICH ONE THE WORLD USES

- TOPOLOGICAL (clause (v) as written / TQO-1, Bravyi–Hastings–Michalakis): NO local operator
  connects the sectors AT ALL. The connecting matrix element is identically zero for every
  operator on a contractible region; splitting and mixing are suppressed as λ^d in the
  region size; the two states are locally indistinguishable everywhere.
- ENERGETIC/KINETIC (Neel–Arrhenius/Kramers): a local connector EXISTS; its dynamical weight
  is exp(−E_b/kT) (or exp(−S/ħ) for tunneling). The two states are locally distinguishable
  — that is what makes the record readable — and the flip is rare, not impossible.

Which does the world use? All ten rows above: energetic/kinetic. Δ ranges 40–90 for
engineered media (the 10-year × ppm design band), ~46 for DNA chemistry, effectively
~700–3000 for persistent currents, "agent-limited" for the geological rows where the barrier
is so high that lifetime is set by the census of energetic events in the environment.

Two structural points, both load-bearing:

(a) THE BATH IS THE PROTECTOR, NOT THE THREAT. Energetic protection WORKS BECAUSE of
    dissipation: Kramers escape is rare precisely because friction re-thermalizes the state
    into the well after every sub-barrier excursion; and einselection (the bath continuously
    monitoring the pointer observable m_z) is what makes the record classical and objective.
    Topological protection is the exact opposite: O-16 measured it — under a dissipative
    bath the model's record dies at rate γ (exponent 2, no d), and the register itself says
    the ε^d law is protection against STATIC perturbation only. So the model's mechanism is
    destroyed by the very environment that STABILIZES the world's mechanism.

(b) NATURE'S BARRIER IS EXTENSIVE; THE MODEL'S DISTANCE IS NOT A BARRIER. K_uV grows with
    grain volume; ΔF_slip grows with wire cross-section; a book's protection grows with the
    number of printed copies; a bit grows safer with the number of grains in it. Meanwhile
    the 2D toric code at finite temperature has an O(1) barrier (create one anyon pair,
    ~2×gap, then diffuse the pair apart for free): its lifetime is INDEPENDENT of L — no
    self-correcting quantum memory in <4D (Alicki–Fannes–Horodecki 2009; Bravyi–Terhal
    2009; Landon-Cardinal–Poulin 2013). Physically instantiated at finite T, the program's
    CANONICAL carrier (O-49/O-50) is a WORSE record than one CoCrPt grain, and gets no
    better as it grows. A hard disk is a better record than a toric code by the model's own
    O-16 criterion.

---

## 3. THE MISMATCH, SAID PLAINLY, AND WHAT HANGS ON IT

THE WORLD USES ENERGETIC AND KINETIC PROTECTION. THE PROGRAM ASSUMES TOPOLOGICAL PROTECTION.
THAT IS A FUNDAMENTAL MISMATCH.

It is not a tolerance gap (O-5's kind, where an exact clause has a natural approximate form
the world satisfies). Clause (v) names the wrong MECHANISM: it demands the non-existence of
a local flipper, where the world's records all possess a local flipper and survive by its
dynamical suppression. No tolerance parameter connects "no operator exists" to "the operator
exists with weight e^(−Δ)" — except a rate bound, which is exactly what clause (v), a
statement in operator algebra with no temperature and no time in it (T-22: nothing in the
corpus is dimensionful), cannot express.

Rows that depend on clause (v) as written (searched REGISTER_V001.md + STATUS_LEDGER_V001.md):

| row | claim | fate under the mismatch |
|---|---|---|
| A-GR (PROVED) | Γ supplies record space, writer, AND the protection (R3) | NARROWS to the model class. Γ protects nothing at finite T in a dissipative world (O-16 is the program's own proof). "Supplies the protection" is false of every record in this census. |
| A-GR2 (PROVED) | Γ supplies the formation CHANNEL via intersection pairing | Formation claim, not protection — logically survives, but inherits model-only scope until a real carrier exhibits the pairing. |
| G-5 (PROVED) | (H₁, pairing) satisfies R1–R3, candidate third term | Survives as mathematics; falls as a candidate mechanism for the WORLD's records — R3-as-homology is instantiated by none of the ten. |
| G-14 (PROVED) | Γ = non-triviality of the homology of EM's complex | Same fate as G-5. |
| G-15 (PROVED) | the condition in clause (v)'s own binary terms: no contractible operation reaches a non-trivial class | FALLS as a world-claim. The world's protection is not a binary and not homological; it is a graded exponent Δ. This is the row most directly contradicted by the census. |
| Theorem D / ε^d protection law | distance suppresses splitting as ε^d | Already rescoped by O-16 to static perturbations; under the census it is a T=0 statement with NO natural referent — none of the ten records' threats are static admissible perturbations. |
| O-12 (PROVED) | admissibility rescues (v) on the toric code | Survives internally; but note the rescue mechanism (energy filter) is exactly the mechanism the WORLD uses — see the amendment. |
| T-11/O-49/O-50 | contractible-region conventions; toric canonical | Survive as model conventions. O-50's "canonical benchmark" carrier is, at finite T, protected by an O(1) barrier — canonical for (v), anti-canonical for real protection. |
| F-9/F-10 (exchange rate) | selection is bought with protection, same exponent d; parity | The model-internal statement survives. NOTE A GENUINE REAL-WORLD ECHO: by detailed balance, the barrier a write must climb IS the barrier that protects (every row of the census: write energy ≳ E_b). "As hard to form as to destroy" is TRUE of the world — in Arrhenius form, not in exponent-d form. Worth registering as the one place the program's structure and the world's agree in shape. |
| G-11 / W-39 / capacity rows | capacity vs protection trade-offs in d | Model mathematics; survive; silent about the world. |
| O-5 (PROVED) | clause (ii) relaxes to width bound, tolerance = inverse lifetime | SURVIVES AND IS PROMOTED: it is the template the amendment generalizes. Checked against the census: yes, it suffices for real records — every row's [H,R] ≠ 0 and [L_k,R] ≠ 0 is exactly a finite-lifetime statement, and the O-5 tolerance 1/T(η) is the Arrhenius rate. Clause (ii) in O-5 form FITS the world. |
| O-16 (PROVED) | distance does not protect against dissipation | VINDICATED — it is the row that discovered this mismatch from inside the model, one lane early. |

The writable/protected tension (months of program history) is resolved by the census in the
world's favor and at the definition's expense: nature never satisfies (iv) and (v) with one
operation class. Every real write is INADMISSIBLE — it does work, 10⁶–10¹⁰ × Landauer in the
engineered rows — and every real protection is a barrier against the THERMAL operation class.
The tension was an artifact of demanding both properties of the same admissible set.

---

## 4. IS ANYTHING IN NATURE TOPOLOGICALLY PROTECTED IN THE PROGRAM'S SENSE?

| candidate | topological? | writable (iv)? | verdict |
|---|---|---|---|
| SC flux quantum / persistent current | winding number: yes; but a LOCAL phase slip flips it over a finite barrier (0.27 eV example above) | yes — flux pump, heater (both inadmissible, do work) | energetic, topology-assisted. FAILS (v)-as-TQO: local connector exists. The closest nature comes, and it is still Arrhenius. |
| Quantized Hall conductance (Chern number) | yes — genuinely no local operator moves it within the gapped phase | NO — it is a property of the Hamiltonian/phase, set by B and density; changing it crosses a phase transition; nothing admissible writes it; arguably fails (iii) too (it is not an operator on one H's eigenspaces but a label of H itself) | protected constant, NOT a record. |
| Total baryon number | nearest to exact: T=0 sphaleron suppression e^(−4π/α_W) ≈ 10⁻¹⁶⁴ | NO for the total (nothing writes it); the WRITABLE version — B in a region — is flipped by ordinary local transport (carry a proton out), i.e. UNPROTECTED | the writable/protected tension realized in nature: the protected quantity is unwritable, the writable one unprotected. |
| Skyrmions (chiral magnets) | continuum topology only; on the lattice, annihilation barriers ~0.1–1 eV, routinely written/deleted with local currents | yes | energetic. Marketing says "topologically protected"; the measured lifetimes are Arrhenius. |
| Vortices / winding in superfluids, domain walls | same structure as flux quantum | yes | energetic (phase slips, finite barriers). |
| FQH / engineered topological order | yes in the ground space | in principle (braiding) | DOES NOT OCCUR AS A NATURAL RECORD; the entire field of topological quantum memory exists because these must be built and actively corrected; finite-T lifetime O(1) in system size (<4D). |

Pattern, and it answers the writable/protected question empirically: in nature, exactly the
quantities that approach clause-(v) protection stop being writable (Chern number, total B),
and everything writable is protected energetically. Where topology appears in a real record
(flux quantum), its role is to LOCALIZE AND RAISE AN ENERGY BARRIER — topology in the service
of Arrhenius, not instead of it.

One more world-mechanism the clause set does not name: REDUNDANCY. Define R as the majority
observable over N copies (the N grains of one HDD bit, the 10¹³ cells carrying one genome,
the surviving copies of one book). Then no operation on a single small region flips R —
a (v)-LIKE clause holds with "distance" ~ N/2. But this is the classical repetition code:
scattered support (excluded by T-11's "single region" reading), locally READABLE (violating
TQO), zero homology. The world's only working (v)-analogue is repetition/extensivity — the
same mechanism as K_uV, one level up. Zurek-style redundancy (the program's W-36/W-40
territory) is the natural-record protection story, and it is not homological.

---

## 5. VERDICT

Clause (v) as written — topological protection, no admissible single-contractible-region
flipper, homology as the mechanism — DESCRIBES NONE OF THE WORLD'S RECORDS.

- Not the engineered ones (R1–R7): their flippers act on their own contractible supports.
- Not the natural archives (R3, R8, R9): same, with agent-limited lifetimes.
- Not the topology-flavored ones (R10, skyrmions, vortices): a local phase slip always
  exists; protection is a finite barrier; topology only shapes the barrier.
- The quantities in nature that DO satisfy (v)'s spirit (Chern number, total baryon number)
  fail clause (iv) and are constants, not records.

Under the admissibility filter [U,H]=0 the failure re-expresses itself but does not go away:
for every census record either the two states are degenerate (then a local admissible
flipper exists on the record's support and (v) fails exactly, by the support-compression
theorem) or they are split (then NO admissible operation writes the record anywhere and (iv)
fails, C-60). The pair {(iv),(v)} is unsatisfiable on contractible supports, and the world
stores its bits on contractible supports because that is what makes them locally readable.

What IS true, and the program should bank rather than mourn:
1. Clause (ii) in O-5's width form fits the world perfectly — tolerance = inverse Arrhenius
   lifetime. O-5 generalizes; it is the template.
2. The write/protect parity (F-9) has a true world analogue via detailed balance: the write
   pays the barrier that protects.
3. O-16 found this mismatch from inside; its scoping of ε^d to static threats is exactly
   right and is the hinge for the amendment.
4. Redundancy/extensivity is a real (v)-like protection the clause set could name honestly.

## 6. MINIMAL AMENDMENT

Replace the binary, homological clause (v) with a graded, dynamical one, on O-5's pattern:

  (v-E) PROTECTED, energetic form: for the record's own open dynamics (H, {L_k}) and any
  bounded perturbation supported on a single contractible region small relative to the
  record's support, the induced flip rate obeys
        Γ_flip ≤ τ₀⁻¹ · exp(−Δ),
  with stability exponent Δ declared alongside the record. Equivalently, in the model's
  operator terms: every product of single-region admissible-or-not operations connecting R
  to −R passes through intermediate states costing energy ≥ E_b = Δ·kT (a barrier
  statement), OR the flipping region must be a non-vanishing fraction of the record's
  support (the repetition/extensivity statement — T-11's own fraction criterion, promoted
  from carrier diagnostic to clause).

  Protection becomes a NUMBER Δ (dimensionless, so T-22's discipline survives), records form
  a one-parameter family, and the world's records enter at Δ ≈ 40–90 (engineered), ~46
  (DNA), ~10²–10³ (persistent currents). Topological carriers are the Δ = ∞ limit — a limit
  the census shows nature does not instantiate and the no-self-correction theorems say
  cannot be instantiated passively below 4D at finite T.

Cost of (v-E) to the corpus: G-15 falls as a world-claim; A-GR's "supplies the protection"
narrows to the model class; Theorem D becomes the Δ=∞/T=0 corner case; G-5/G-14 survive as
mathematics with their world-relevance an OPEN claim needing a carrier that realizes the
pairing; O-5, O-16, C-60/O-4, F-9's parity, and all capacity mathematics survive; the
writable/protected tension dissolves (writes are inadmissible barrier-crossings; protection
is the same barrier faced thermally). The 162 PROVED rows keep their proofs; what changes is
which of them are ABOUT anything.

Next step (per program rule that no route closes without one): commission the counter-lane —
build the model's OWN Arrhenius record (a double-well / repetition carrier with a declared Δ
inside the existing (H,{L_k}) machinery) and check which of Theorems A–D survive with Δ in
place of d. That is the test of whether the program's structure needs homology or only
needed A BARRIER, and it is checkable with the corpus's existing tools.

## NUMBERS I AM LEAST SURE OF (stated per the lane's honesty rule)
- Flash retention effective Ea (1.1 eV) is a JEDEC modeling convention, not a single
  physical barrier; true leakage is defect-dominated and device-specific.
- CoCrPt K_u spans 2–6×10⁵ J/m³ across generations; per-grain Δ at the low end is ~29 kT
  and the bit relies on multi-grain redundancy.
- Al coherence length/cross-section in R10 are illustrative; ΔF varies by orders of
  magnitude across materials — the CONCLUSION (finite, local, Arrhenius) does not.
- Skyrmion annihilation barriers (0.1–1 eV) are from simulation literature, temperature-
  and material-dependent.
- File & Mills persistent-current bound is remembered as >~10⁵ yr; the exact figure should
  be checked against the 1963 paper before any register entry cites it.
- Crater impactor KE uses order-of-magnitude π-scaling only.
No web access was used; all values are standard-literature figures from training knowledge,
recomputed for consistency in gr2_arrhenius_numbers.py.
