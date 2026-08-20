# LANE_T39_A_ENUM — ALPHA-THROUGH-E_B (probe 1) AND THE FULL CANDIDATE ENUMERATION

Date: 2026-08-20. Inputs honored: GR1 census (E_b/kT column), C-69/C-70 (instrument, demoted),
C-71/C-72, C-76 (corrected), C-79/C-80/C-81/C-82/C-83, C-41, C-12, PF-6 (the prior FAILED run
of this same audit). Arithmetic in `t39_alpha_amplification.py` / `.txt`.

## JOB 1 — THE DERIVATION

Step 1, mechanism-correct scalings of E_b with the fine structure constant (m_e c^2 held fixed;
per-atom energies so lattice-constant factors cancel; owners named):

- s = d ln E_b/d ln alpha = 2 — chemical bond, diffusion activation, band offset, electrostatic:
  all Hartree/Rydberg scale, E_h = alpha^2 m_e c^2. Owner: standard atomic/molecular physics;
  operative throughout the varying-constants literature (Uzan 2003 review; Flambaum's enhanced-
  sensitivity program; King et al. 2012 many-multiplet).
- s = 4 — shape (magnetostatic/dipolar) anisotropy: mu_B^2/a0^3 = alpha^4 m_e c^2 / 4, exact
  identity. Owner: textbook dipolar scaling; Neel for the anisotropy application.
- s = 6 — uniaxial magnetocrystalline anisotropy: 2nd order in spin-orbit, K ~ xi^2/W with
  xi ~ Z^2 alpha^4 m_e c^2 and W ~ Hartree ~ alpha^2 m_e c^2. Owner: van Vleck 1937; Bruno 1989.
- s = 10 — cubic magnetocrystalline anisotropy (4th order in SOC, xi^4/W^3). Listed for
  completeness; the census magnetite barrier is SHAPE anisotropy (s = 4), not cubic MAE.
- PVED: weak-scale (G_F, ~Z^5 enhancement; Letokhov 1975, Quack). Not an Arrhenius barrier in
  the census; excluded from the amplification table, recorded here as instructed.

Step 2, the record-law consequence. tau = f0^-1 exp(E_b/kT) gives
  d ln tau / d ln alpha = (E_b/kT) * s  +  O(1)  [the O(1) is d ln f0^-1/d ln alpha; attempt
  frequencies are phonon/gyromagnetic, alpha-dependence O(1-10) in the log, negligible against
  the 86-3120 exponent terms].

Step 3, numbers (from the census's own E_b/kT column; full table in t39_alpha_amplification.txt):

| record | E_b/kT | s | A = d ln tau/d ln alpha | 1% alpha shift multiplies tau by |
|---|---|---|---|---|
| HDD CoCrPt grain | 61 | 6 | 366 | 39 |
| SD magnetite TRM | 780 | 4 | 3120 | 3.6e13 |
| Zircon Pb closure | 220 | 2 | 440 | 81 |
| Flash (3.1 eV barrier) | 120 | 2 | 240 | 11 |
| DNA depurination | 50 | 2 | 100 | 2.7 |
| AgBr latent image | 45 | 2 | 90 | 2.5 |
| CMOS latch (powered) | 50 | 2 | 100 | 2.7 |
| generic 60 kT chemical | 60 | 2 | 120 | 3.3 |
| generic 60 kT uniaxial-SOC | 60 | 6 | 360 | 37 |

The brief's example: a 1% alpha shift moves a 60 kT chemical barrier's lifetime by e^1.2 = 3.3x,
and the same 60 kT barrier's lifetime by e^3.6 = 37x if the barrier is uniaxial-SOC magnetic.

Derived side results:
- RATIO SIGNATURE. For coexisting records of one epoch the log-lifetime shifts stand in fixed
  ratios set only by (E_b/kT)*s: magnetite-vs-zircon = 3120/440 = 7.1; the mechanism exponents
  alone stand as 2:4:6.
- TUNNELING FLATNESS. Flash's tunneling channel has exponent kappa*d with kappa ~ alpha m_e c/hbar
  and d = N a0 ~ N/(alpha): kappa*d ~ N, alpha-INDEPENDENT at fixed oxide atom count. A
  tunneling-limited retention channel is alpha-flat while every Arrhenius channel is
  alpha-amplified.
- C-83 COMPOSITION. W = tau/t_epoch inherits the full amplification, so the certifiability
  crossover n* ~ 6W moves as exp[(E_b/kT) s (dalpha/alpha)] — the one alpha sentence that runs
  through program-owned machinery (no rival owns CERT_W).

HONEST OWNERSHIP OF THE LOGIC: "a natural record's survival amplifies sensitivity to varying
constants" is OWNED in the decay-constant domain — Peebles & Dicke 1962 (Re-187), Shlyakhter
1976 / Damour & Dyson 1996 (Oklo, ~10^6-10^7 amplification via a 0.1 eV resonance), Olive et al.
2002. The Arrhenius-barrier version for material record media, and the 2:4:6 mechanism-ratio
signature, we did not find stated by any rival — but every ingredient is owned, and any bound
extractable from it is uncompetitive against Oklo/quasar (|dalpha/alpha| < 1e-7 at Gyr). The
advocates decide whether "equivalent statement" covers it.

## JOB 2 — the candidate list is returned as structured data to the audit (T-39); this file is
the lane's derivation record. Candidates CAND-1..CAND-12 with exposure marks are in the return.

## ERRATA NOTED AGAINST THE BRIEF'S INVENTORY (verified in the ledger before use)
- "w_min = d (C-80)" — the ledger's C-80 is the O-54 falloff judgment. The weight-d material
  (threshold = d, splitting ~ eps^d, chi = 0.1145 bits at weight d and 0 below) lives in
  F-13/F-24/PF-6/A-PR, and PF-6 attributed its substance to Bravyi-Hastings-Michalakis (H-2).
- C-71/C-72 carry PROVED status with the registrar note "VERIFIED ON ONE MECHANISM ONLY" per
  the strengthened 2026-08-20 bar (each law's second mechanism is the other's first; the
  cross-mechanism kernel awaits the written-vs-erased NAND run). Used with that caveat.

## NEXT STEP (no route closes without one)
Whichever alpha candidate survives the advocates: the sharpest falsifier is the flash
Arrhenius-vs-tunneling channel discrimination — one retention-bake dataset (Arrhenius plot
curvature vs oxide thickness) separates the alpha-amplified channel from the alpha-flat one,
turning ALPHA-3 from a scaling remark into a measured channel split. Commission it as a
literature-extraction lane (no lab needed; JEDEC retention data exists).
