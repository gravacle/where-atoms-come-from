# LANE_GR1_CENSUS — THE FIVE CLAUSES AGAINST THIRTEEN REAL RECORDS

Date: 2026-08-20. Inputs honored: O-5 (approximate (ii) = spectral width, tolerance T(η) ≥ η/δ,
δ = inverse lifetime), O-4/DEF-A (ADMISSIBLE := [U,H]=0, provisional), C-60 ([U,H]=0 ⇒ dE=0),
C-12 (converse: projection in commutant, non-trivial on an eigenspace, trace-balanced), O-12
(clause (v) admissibility repair), T-22 (nothing in corpus dimensionful).

Arithmetic in `numbers.txt` (script inline there). kT(300 K) = 4.14e-21 J = 25.9 meV;
kT ln 2 = 2.87e-21 J = 17.9 meV.

## THE CENSUS SET (13 records, 24 orders of magnitude in lifetime-energy product)

1. **HDD grain** — CoCrPt-SiO2 perpendicular media, cylinder d≈8 nm h≈10 nm, V=5.0e-25 m³,
   K_u≈5e5 J/m³ (industry range 2–7e5), M_s≈6e5 A/m. K_uV = 2.5e-19 J = **61 kT** (design target
   ~60 kT for 10-yr retention; Sharrock/Néel-Arrhenius with τ₀~1e-9 s gives τ ~ 2e17 s for the
   median grain; the 10-yr spec is set by the distribution tail and demag fields, not the median).
   Switching field μ₀H_K = 2K_u/M_s ≈ 1.7 T; head fields ~1–1.8 T over ~30 nm.
2. **Silver-halide grain** — AgBr microcrystal 0.5–2 μm; latent image = Ag_n cluster, n_c ≈ 3–6
   atoms at a sensitivity site; write ≈ 4–10 absorbed photons × ~2.6 eV ≈ 10–26 eV; development
   gain ~1e9. Latent-image regression time ~months–years (τ ~ 3e7 s, poorly controlled number).
3. **Zircon U-Pb** — ²⁰⁶Pb/²³⁸U ratio over ~1e15 atoms per crystal; Pb-diffusion Ea ≈ 550 kJ/mol
   = 5.7 eV/atom (Cherniak & Watson); Boltzmann factor at 300 K = e^-220; closure ~900 °C;
   demonstrated retention 4.4 Gyr = 1.4e17 s.
4. **CMB photon polarisation** — single free photon; polarisation doublet exactly degenerate at
   fixed k (E = ħω both states); free-streaming 13.8 Gyr = 4.4e17 s; width ħ/τ ~ 1.5e-33 eV.
5. **Apparatus pointer** — latched CMOS bistable (SRAM cell / readout latch): two circuit states
   symmetric by design; barrier ~40–60 kT while powered (thermal flips unobservable; SER is
   radiation-dominated); write ≈ CV²/2 ~ 5e-16 J = **1.7e5 × kT ln 2**; VOLATILE — the record
   exists only in the driven steady state.
6. **Cosmic-ray / fission track in mica** — damage trail ~10 μm × 5 nm, ~1e4–1e5 displaced atoms,
   stored strain/defect energy ~ 5e5 eV = 2e7 kT (few % of the ~100 MeV deposited); muscovite
   retention ~1e9 yr = 3e16 s (apatite anneals at ~110 °C, zircon ~240 °C, muscovite higher).
7. **DNA base pair** — 4-valued (2 bits); bare durability: depurination t½ ≈ 730 yr = 2.3e10 s at
   37 °C (Lindahl); ancient-DNA survival 1e4–1e6 yr cold/dry; in vivo maintained indefinitely by
   ACTIVE repair. Write: polymerase + proofreading ~2 ATP + pyrophosphate ~ 20–100 kT/base.
8. **Potentiated synapse** — weight is graded (EPSP amplitude, receptor count 0–200); substrate
   turns over in hours–days; persistence is a self-reinforcing biochemical attractor (CaMKII
   bistability — CONTESTED in the literature — plus receptor-recycling loops); potentiation costs
   ~1e4–1e5 ATP ~ 3e3–3e4 eV.
9. **Tree ring** — ring width/density: continuous analog; preserved wood ~1e4 yr (bristlecone
   chronology ~9e3 yr) = 3e11 s; written by a season of photosynthesis: ~0.5 kg biomass ≈ 8e6 J
   = **3e27 × kT ln 2** — the least admissible write in the census.
10. **Lava-flow TRM** — SD magnetite grains ~60 nm, M_s = 4.8e5 A/m; shape-anisotropy KV/kT ≈ 780
    (τ beyond astronomical; Gyr paleorecords exist); grain-level bit = ±M along easy axis;
    ensemble record (field direction) is continuous.
11. **Flash floating gate** — ~100 e⁻ (modern 3D NAND: tens) behind 6–8 nm SiO2, barrier 3.1 eV;
    retention spec 10 yr = 3e8 s; stored electrostatic energy ≈ QΔV/2 ≈ 2.4e-17 J = 150 eV =
    5.8e3 kT; write ≥ QΔV ~ 4e4 × kT ln 2 at the gate (~1e8 × with peripherals).
12. **Ice-core gas bubble** — CO2 ppm in trapped air: continuous; EPICA record 800 kyr = 2.5e13 s;
    degradation by clathrate formation and diffusion; written by firn compaction (dissipative).
13. **Footprint in rock** — analog shape; ~1e8 yr = 3e15 s for lithified trace fossils; grains
    rearranged: per-grain ΔE ~ mgh ~ 1.7e-9 J = 4e11 kT — configurationally "cheap" but nowhere
    near degenerate on any microscopic tolerance; write = a footstep, 1–100 J dissipated.

## THE TABLE — record × clause

P = pass, AP(x) = approximate pass with stated tolerance, F = fail. The number in each cell is the
one that decides it.

| # | Record | (i) BIT | (ii) DURABLE (O-5 form) | (iii) NON-TRIVIAL (same-energy) | (iv) WRITABLE + ADMISSIBLE (DEF-A) | (v) PROTECTED |
|---|---|---|---|---|---|---|
| 1 | HDD grain | **P** (Ising ± below blocking) | **AP** τ=2e17 s, δ=ħ/τ=4.8e-52 J | **AP** split 0.007 kT (Earth) to ~7 kT (50 mT stray) ≪ 61 kT barrier | **F** write dissipates ~2K_uV = 175 × kT ln2; dE≠0 | **F** 30-nm head region, 1.7 T, flips one grain; ~5e-19 J |
| 2 | AgBr grain | **P** (developable threshold n≥n_c) | **AP** τ~3e7 s, δ~3e-42 J | **F** exposed vs unexposed differ ~10–26 eV chemical | **F** write = +16 eV absorbed ≈ 870 × kT ln2 | **F** local heat (~100 °C, minutes) disperses Ag_n; ~eV/atom |
| 3 | Zircon U-Pb | **F** continuous ratio (analog) | **AP** τ=1.4e17 s, δ~7.5e-52 J | **F** U-config vs Pb-config differ by MeV (nuclide masses) | **F** "write" = α-decay, 4–5 MeV each ≈ 1e8 × kT ln2 | **F** local heating >900 °C, μm³ ≈ 3e-9 J = 1e12 kT |
| 4 | CMB photon pol. | **P** exact qubit | **P exact** [H_free,R]=0, δ→0 (ħ/τ=1.5e-33 eV en route) | **P exact** polarisation doublet degenerate at fixed k | **P** waveplate unitary on degenerate doublet: dE=0, ADMISSIBLE | **F** the same local waveplate flips it at ~zero energy |
| 5 | Latched pointer | **P** (by design) | **AP** driven: τ~e^60×ns~1e17 s while powered; **0 s unpowered** | **F/AP** symmetric by design but mismatch ΔE ~ 1e3 kT; barrier-metastable, not degenerate | **F** CV²/2 = 1.7e5 × kT ln2 | **F** the word line: local, ~fJ |
| 6 | Track in mica | **P** (present/absent per site) | **AP** τ=3e16 s, δ~3.5e-51 J | **F** track stores ~0.5 MeV = 2e7 kT of strain | **F** write = 100 MeV fragment ≈ 1e10 × kT ln2 | **F** local anneal: μm³ to 400 °C ~ 1e-9 J |
| 7 | DNA base pair | **P** (4-valued → 2 bits, discrete) | **AP** bare τ=2.3e10 s; in vivo held BY jump operators (repair) | **F** A·T vs G·C ≠ same molecule; tautomer/mutation ΔE ~ 0.2–1 eV | **F** ~20–100 kT dissipated per base written | **F** one UV photon (4.9 eV) or one mutagen molecule at the site |
| 8 | Synapse (LTP) | **F** graded weight (analog attractor) | **AP** only as dissipative attractor; substrate τ ~ days, record τ ~ years | **F** potentiated state = more protein, ΔE ~ 3e4 eV | **F** ~1e5 ATP ≈ 1e6 × kT ln2 | **F** local low-freq stimulation / ZIP depotentiates |
| 9 | Tree ring | **F** continuous width/density | **AP** τ=3e11 s, δ~3.3e-46 J | **F** different mass distribution, ΔE ~ MJ | **F** 3e27 × kT ln2 | **F** local carving/burning |
| 10 | Lava TRM (grain) | **P** grain ±M; **F** ensemble direction (continuous) | **AP** τ~e^780×1e-9 s, δ≈0 | **AP** grain: Earth-split 2.5 kT ≪ 780 kT barrier | **F** written by cooling through T_B ≈ 855 K: heat flow ~ J/cm³ | **F** local reheat >580 °C or local ~0.1 T field |
| 11 | Flash gate | **P** SLC (MLC: 4–16 discrete levels) | **AP** τ=3e8 s, δ=3.3e-43 J | **F** charged vs empty differ 150 eV = 5.8e3 kT electrostatic | **F** ≥ 4e4 × kT ln2 at the gate | **F** 254-nm UV (4.9 eV > 3.1 eV barrier) — that IS the EPROM eraser |
| 12 | Ice-core bubble | **F** continuous ppm | **AP** τ=2.5e13 s, δ~4e-48 J | **F** composition = particle-number sectors, not degenerate | **F** compaction work, ≫kT per molecule sorted | **F** local melt, ~J/cm³ |
| 13 | Footprint | **F** analog shape | **AP** τ=3e15 s, δ~3.5e-50 J | **F** ΔE ~ 4e11 kT per displaced grain (smallest fractional offender, still fails any width tolerance) | **F** footstep, 1–100 J ≈ 1e20 × kT ln2 | **F** local erosion/smoothing, ~mJ |

## TALLIES

- **(i) BIT: 8/13 pass** discretely (grain-level readings), **5/13 fail** — the analog records
  (zircon ratio, synapse weight, tree ring, ice core, footprint; TRM ensemble direction also).
  For the discrete 8, multi-valued cases (DNA 4-ary, MLC flash 16-ary) are harmless coarse-
  graining: replace the involution by a commuting spectral family of projections — C-12 is
  ALREADY stated at projection level, so the corpus's own converse survives this generalization.
  For the analog 5 it is not harmless: no involution exists to test; the definition is silent on
  half the world's records rather than false of them.
- **(ii) DURABLE: NEVER FAILS in the O-5 approximate form — 12 AP + 1 exact pass.** Lifetimes
  run 3e7 s (latent image) to 2e17 s (HDD median, zircon, CMB); the O-5 width δ = ħ/τ runs
  1e-41 to 1e-52 J — comfortably below every relevant energy scale by ≥ 20 orders. O-5's
  relaxation is VINDICATED as necessary and sufficient in form. One caveat with teeth: for the
  three MAINTAINED records (pointer, in-vivo DNA, synapse) durability is delivered BY the jump
  operators, not despite them — [L_k, R] = 0 is the wrong sign; the L_k actively restore R.
  These need the dissipative-attractor reading (O-16/F-13 territory), not O-5's passive form.
- **(iii) NON-TRIVIAL: FAILS MOST OFTEN — 9 clean fails, 3 approximate passes, 1 exact pass.**
  Passes exactly only for the CMB photon (true degenerate doublet); approximately only for the
  magnetisation pairs (HDD grain, TRM grain: time-reversal-paired ±M, splitting 0.007–7 kT ≪
  barrier 61–780 kT) and arguably the symmetric latch. Every chemical, structural, charge, and
  compositional record stores its value in states of grossly different energy: 0.2 eV (DNA
  tautomer) → 150 eV (flash) → 0.5 MeV (track) → MJ (tree ring).
- **(iv) ADMISSIBLE-WRITABLE under DEF-A: 12/13 FAIL; 1 passes** (waveplate on the photon's
  degenerate doublet — the single genuinely workless write in the census). Real write costs run
  60 × kT ln 2 (DNA) to 3e27 × kT ln 2 (tree ring). Under DEF-A + C-60 (dE = 0 identically), NO
  real record except the photon is writable in the program's sense. All thirteen are obviously
  writable physically — the world wrote them.
- **(v) PROTECTED: 13/13 FAIL as physical operations.** Every census member has an explicit
  single-region eraser with a finite energy price: write head (5e-19 J), UV photon (4.9 eV),
  local anneal (1e-9 J), mutagen (eV), waveplate (~0). Real protection is a COST FLOOR
  (E_barrier ≫ kT), never a kinematic impossibility. Under DEF-A's admissible reading (v)
  passes VACUOUSLY for every barrier record — any local flipper must do work ≥ E_b and is
  therefore inadmissible — but that same vacuity is exactly what makes (iv) fail. As written,
  (iv)+(v) jointly demand "flippable free globally, unflippable free locally"; real records are
  "flippable at cost everywhere, with a locality-independent floor E_b."
- **ALL FIVE AS WRITTEN: 0/13.** Best performer: the CMB photon, exact on (i)–(iv), killed by
  (v) — and it is a degenerate case in the other direction: unprotected, unread, no redundancy.
  Best material record: the HDD grain — P, AP, AP, F, F.

## THREE STRUCTURAL FINDINGS (the census's actual content)

**FINDING 1 — exact (ii) + exact (iii) are JOINTLY unsatisfiable in generic real matter.** If
[R,H] = 0 exactly and the spectrum of H is non-degenerate, R is a function of H, hence constant
on every eigenspace — (iii) fails identically. Generic macroscopic Hamiltonians have no exact
degeneracies (level repulsion); exact degeneracy survives only via symmetry (time reversal →
the ±M pairs) or engineered topology (the corpus's codes). This is WHY the census's only
(iii)-passers are the magnetisation pairs and the photon: they are the census's only
symmetry-protected doublets. The corpus's carriers pass (ii)+(iii) exactly because they were
BUILT degenerate; the world's records are built METASTABLE instead. Two different protection
mechanisms: degeneracy-protection (corpus) vs barrier-protection (world). The five clauses
formalize the first; nature overwhelmingly uses the second.

**FINDING 2 — the failed clauses (iii) and (iv) are BOTH cured by one move: draw the system
boundary around carrier + local bath.** Microcanonically, "track present + hotter phonons" and
"track absent + cooler phonons" have the SAME total energy: the configuration energy is exported
to the bath, and R (dressed) is non-constant on the closed system's energy shell. Quantitatively
the enlargement always works: even a 1 cm³ phonon bath has energy-shell width √(C k_B T²) ≈
1.6e-9 J ≈ 10 GeV — 2e4 × the worst configuration offset in the census (track, 0.5 MeV) and
1e8 × flash's 150 eV. Simultaneously, the physical write IS an energy-conserving unitary on
carrier + writer + bath (this is O-4's untested "or a physical channel" reading): admissibility
holds at the dilation level, C-60 survives (total dE = 0), and Landauer's kT ln 2 reappears
correctly as the FREE-energy cost paid into the bath rather than as an energy non-conservation.
The clauses were written for the carrier alone; real records satisfy their natural forms only on
carrier ⊗ bath. This also answers H-6 from an unexpected side: the bath must be large enough to
thermalize the write energy — records are macroscopic because the CLAUSES (amended) demand a
macroscopic bath, not because the carrier must be big.

**FINDING 3 — clause (v) must become an energy statement or remain vacuous.** No real record
has kinematic protection; every one has energetic protection, measured by one number, E_b/kT:
61 (HDD), ~120 (3.1 eV flash barrier), ~50 (DNA glycosidic 1.3 eV), 220 (zircon Pb diffusion),
780 (SD magnetite), ~40–60 (latch). The empirical form of (v) is: every single-region operation
that flips R must transfer energy ≥ E_b, with E_b/kT ≫ 1. The corpus's toric-code result (O-12:
zero admissible single-region flippers, minimum admissible writer weight d) is the E_b → ∞
kinematic limit of this. Retention times then follow from (v') + Arrhenius, which is exactly
O-5's tolerance — clauses (ii) and (v), amended, become two faces of ONE number, the barrier.

## MINIMAL AMENDMENTS AND THEIR COST TO THE CORPUS

- **(i') SPECTRAL FAMILY**: replace the involution R by a commuting family of projections {P_a}
  (analog records: a POVM/pointer decomposition). COST: LOW. C-12 already lives at projection
  level. O1-B/O1-B1 (trace parity, odd-dimension corollary) apply per projection pair; the F₂
  chain-complex corpus (O-7/O-8, O-10) survives as the discrete subclass. Analog records are
  new territory, not counterexamples.
- **(ii') = O-5 as registered**, extended with the dissipative-attractor case ([L_k,R] restoring
  rather than commuting) for maintained records. COST: NONE to existing rows; O-11/O-16 move
  from side gallery to load-bearing.
- **(iii') ON-SHELL, ENLARGED**: R (dressed) non-constant on an energy shell of carrier ⊗ local
  bath; shell width √(C k_B T²) ≥ ΔE_config, verified above for all 13. COST: MEDIUM. Any row
  that reads (iii) as "records carry zero configuration energy" is FALSIFIED by the census —
  records generically carry configuration energy (up to MeV here) and EXPORT it. Rows using
  carrier-level degeneracy as an input (the W-4x record-counting/metric chain, W-42 capacity,
  W-43/W-44 metric, W-51 record-to-geometry) survive only if their record COUNT is insensitive
  to well-depth asymmetry — needs an explicit re-check, not a presumption. GLOSSARY's "we know
  of no standard name" for (iii) now has an answer: the standard object is a metastable well
  pair, and its name in nature is not degeneracy.
- **(iv') DILATION ADMISSIBILITY** (O-4's untested second reading): a write is admissible iff it
  is implemented by an energy-conserving unitary on carrier ⊗ writer ⊗ bath; work accounting via
  free energy, floor kT ln 2 + ΔE_config (Landauer). COST: MEDIUM-LOW. C-60 survives verbatim at
  the total level. O-12's counts, O1-B's trace-balance criterion, and C-12's (H2) were derived
  for system-level unitaries; whether trace balance survives dilation-writability is OPEN — flag
  for a dedicated lane. DEF-A rows survive as exact statements about a definition real writes
  never satisfied.
- **(v') COST FLOOR**: every single-contractible-region channel flipping R requires energy
  transfer ≥ E_b, E_b ≫ kT. COST: LOW. O-12/T-11/P-3-weak survive as the kinematic limit;
  nothing in the corpus asserts real-world kinematic protection (H-13 already registers the
  negative as clause-relative).

**Net cost:** none of the 162 PROVED rows is falsified AS A STATEMENT ABOUT ITS OWN OBJECTS.
What the census kills is any reading of them as statements about the world's records: under the
clauses AS WRITTEN, the world contains approximately ZERO records, and the corpus's carriers are
the only known members of the defined class. The amended clauses (i')–(v') are satisfied by all
13 census members with sane tolerances — and they are not a retreat to vagueness: each has one
number (δ = ħ/τ; shell width √(C k_B T²); E_b/kT) that decides it.

## NUMBERS I AM UNSURE OF (bounds, not citations)

- CoCrPt K_u: used 5e5 J/m³ (lit. range 2–7e5); conclusion insensitive within range.
- Latent-image regression τ: ~3e7 s is a rough figure; anywhere in 1e6–1e9 s changes nothing.
- Fission-track stored energy: "few % of 100 MeV" is an estimate; 1e4–1e6 eV bounds it; (iii)
  fails at any point in that range.
- CaMKII bistability as THE synaptic tag is contested; the attractor picture (and every cell
  entry) survives whichever molecular substrate wins.
- SRAM barrier 40–60 kT and mismatch asymmetry ~1e3 kT: device-dependent order of magnitude.
- Flash electron count: 3D NAND is now tens of electrons; used 100 (older planar); scales cells
  by <10×, no verdict changes.
- Muscovite track retention 1e9 yr: extrapolated from annealing kinetics, order of magnitude.

## NEXT STEP (no route closes without one)

Commission LANE_GR2-adjacent work: (a) re-derive C-12's trace-balance condition (H2) under
dilation admissibility (iv') — it is the one corpus theorem whose survival is genuinely in
doubt; (b) re-run one W-4x record-count result with asymmetric wells (ΔE up to the shell width)
to test whether the gravity chain needs (iii) or only (iii'); (c) put the HDD grain — the best
real performer — through the model as an (H, {L_k}) instance: two-state macrospin, Arrhenius
L_k, and check which clauses the MODEL itself certifies. That is the first PROOF grounded in a
world record.
