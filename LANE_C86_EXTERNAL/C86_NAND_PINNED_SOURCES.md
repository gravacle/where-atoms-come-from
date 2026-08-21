# C-86 EXTERNAL-DATA RUN — NAND/MLC ARM — SOURCE-PINNED (2026-08-21)

Target: falsifier (i) of the C-86 row — JEDEC-class MLC bake data, retirement-staircase
positions, margin-free. House style follows LANE_T41_EXTERNAL/CITATIONS.md: every pin carries
SOURCE, NUMBER, UNCERTAINTY, and SEMANTICS. A number with no traceable primary source is
pinned as CLASS, never as datum. The azobenzene lesson governs throughout: the killers here
are semantic (what "failure" means, what units the axis is in, what temperature the clock ran
at), not arithmetic. Arithmetic on these pins: `c86_nand_arm.py` → `c86_nand_arm.txt`.

## A. THE JEDEC FRAME — what "retirement" means in the class the falsifier names

1. **JESD218B Table 1 (SSD classes and requirements)** — read first-hand from the standard
   (JEDEC JESD218B, scanned public copy, p. 6): Client — retention **1 year at 30 °C**
   power-off, UBER **≤ 1e-15**, FFR **≤ 3%**; Enterprise — retention **3 months at 40 °C**
   power-off, UBER **≤ 1e-16**, FFR **≤ 3%**. UNCERTAINTY: none — normative spec values.
   SEMANTICS: a "data error" is a drive returning wrong/unrecoverable data to the host
   (§3.2) — i.e., AFTER the controller's ECC and at the vendor's read-reference margins.
   **JEDEC-class retirement is UBER/margin-defined, drive-level, ECC-downstream. It is not
   a record census and is margin-laden by construction.**

2. **JESD218B §6.1.3 + Table 3 + Table 4 (the bake machinery)** — same document, pp. 11–13:
   high-temperature retention stress legs **96 h at ≥ 66 °C** or **500 h at ≥ 52 °C** (both
   classes); target stress-temperature table computed "assuming an activation energy of
   **1.1 eV**" (§6.1.3 text, Annex A normative), with the standard's own stated motive of
   adding "margin against possible inaccuracies in the 1.1 eV acceleration model".
   UNCERTAINTY: the standard treats 1.1 eV as an assumption, not a measurement.
   SEMANTICS: the 1.1 eV maps stress time to use time; it is the clock conversion the whole
   JEDEC class runs on.

3. **JESD218B §3.19 (the standard's own two-mechanism warning)** — p. 4: retention failure
   mechanisms are "accelerated in different ways by temperature"; charge leak through the
   transfer dielectric "can be weakly accelerated or even decelerated by high temperature",
   detrapping "can be highly temperature accelerated (see JEP122)"; "it is not known a
   priori whether high temperature is worse than low temperature." SEMANTICS: **the
   standard's own text denies that one Ea covers retention** — the 1.1 eV is the detrapping
   channel only. Any bake-positioned staircase inherits this caveat at full strength.

## B. ACTIVATION ENERGIES — the JEP122-class numbers

4. **Detrapping Ea ≈ 1.1 eV, measured on a 3D FG array (open-access primary)** — Malavena
   et al. (Politecnico di Milano + Micron Vimercate), *Analysis of High-Temperature Data
   Retention in 3D Floating-Gate NAND Flash Memory Arrays*, IEEE J. Electron Devices Soc.
   11, 524–530 (2023), DOI 10.1109/JEDS.2023.3320722, open PDF via re.public.polimi.it.
   Bakes 75–150 °C, RT reads; post-cycling ΔVT transients merge onto one trend under time
   shifts whose Arrhenius plot gives "a single activation energy EA nearly equal to
   **1.1 eV**", called "another specific signature of charge detrapping observed over all
   Flash technology generations". UNCERTAINTY: "nearly equal"; no error bar quoted.
   SEMANTICS: ΔVT of monitored cells between RT reads bracketing bakes; **voltages
   "normalized to the same arbitrary constant"** — even Micron's own published array data
   withholds absolute volts.

5. **Detrapping Ea ≈ 1.0 eV constant across generations and cycling (CLASS pin)** —
   *Activation Energies (Ea) of Failure Mechanisms in Advanced NAND Flash Cells for
   Different Generations and Cycling*, IEEE TED (2013); paywalled (ResearchGate/IEEE),
   abstract-level content only: detrapping Ea "almost the same (approximately **1.0 eV**)
   regardless of the generation or cycling times". PINNED AS CLASS. This is the row that
   T-41's citation #3 (≈ 1.0 eV, constant) traces to at the same access depth.

6. **The industry bake-equivalence table inverts to exactly 1.1 eV (computed control)** —
   Li, Ye, Kuo, Xue, *How the Common Retention Acceleration Method of 3D NAND Flash Memory
   Goes Wrong?*, HotStorage 2021 (Best Paper; open slides at hotstorage.org), slide 5:
   1 year at 25 °C declared equivalent to 97.65 h @ 60 °C, 11.16 h @ 80 °C, 1.61 h @ 100 °C,
   0.28 h @ 120 °C. Inversion (instrument §1): implied Ea = **1.100, 1.100, 1.100, 1.101 eV**
   — one number to the table's own quantization, agreeing with pin 2. UNCERTAINTY: ±0.001 eV
   from table rounding. SEMANTICS: this is the conversion clock, pinned from the practicing
   side rather than the standards side.

## C. LEVEL-RESOLVED RETENTION STRUCTURE — the data the falsifier actually needs

7. **Per-state mean/σ vs retention time, numeric table (the best public census surrogate)** —
   Cai, Ghose, Haratsch, Luo, Mutlu, *Error Characterization, Mitigation, and Recovery in
   Flash-Memory-Based Solid-State Drives*, Proc. IEEE 105(9) (2017); open PDF arXiv:1706.08642,
   Appendix Table 5. Real TLC chips, FPGA read-retry platform, **2,000 P/E** (1-day row is
   byte-identical to Table 4's 2,000-cycle row), retention **1 day → 1 year**, room
   temperature (the companion MLC study, HPCA 2015, states 20 °C explicitly), **normalized
   Vt: nominal max = 512, 0 = GND — "absolute threshold voltage values are proprietary
   information to flash vendors"**. Transcribed in full into the instrument. Computed
   margin-free (instrument §2, floors stated there):
   - per-state 1-day→1-year mean drift: ER +23.3, P1 +10.0, P2 +6.1, P3 +3.3, P4 +0.4,
     P5 −2.3, P6 −4.7, P7 −7.3 (units of the normalized scale);
   - every programmed adjacent gap SHRINKS; d(gap)/yr from −2.4 (P5-P6) to −3.9 (P1-P2);
   - resolvability: P1-P2 falls 3.2 → 2.3 combined-σ over the year; **no pair merges**.
   UNCERTAINTY: table quantization 0.1 units; per-state σ 8.5–12.8 units is published.
   SEMANTICS: distribution means from full read-retry sweeps — margin-free in exactly the
   falsifier's sense; but ROOM temperature, one anonymized chip, normalized axis.

8. **MLC level-resolved shape (the falsifier's own device class)** — Cai, Luo, Haratsch,
   Mai, Mutlu, *Data Retention in MLC NAND Flash Memory: Characterization, Optimization,
   and Recovery*, HPCA 2015 (open PDF, ETH/CMU): 2y-nm (20–24 nm) MLC, **20 °C**, retention
   1–40 days, up to 8k P/E, 5 s dwell; Finding 3: **ΔP3 > ΔP2 > ΔP1**, P1 mean "remains
   almost constant"; P2, P3 shift down, all distributions widen. Numbers live in figures
   (Figs. 2–4), normalized axis. UNCERTAINTY: figure-resolution. SEMANTICS: same margin-free
   distribution semantics as pin 7; the paper REJECTS Arrhenius bake acceleration for
   characterization on record: it "may exaggerate some causes of retention loss over others".

9. **Level-resolved failure structure under REAL 1-year retention vs bake (3D)** — Li et al.
   (pin 6's paper), slides 7–15: four real 3D chips (3 TLC + 1 QLC; 2 FG + 2 CT), YEESTOR
   9083 platform, **one full year at 25 °C** side-by-side with 60/80/100/120 °C bakes.
   Quantified findings: bakes **underestimate** real long-retention RBER (CT worse than FG;
   RBER axes 1e-5…1e-1 per wordline in slides); per-read-voltage (V1…V7 / V1…V15) RBER
   ordering **changes with temperature** — "different voltage states suffer from retention
   errors in different degrees"; real-time per-voltage RBER time series over 366 days exist
   (slide 15). UNCERTAINTY: figure-resolution (numbers in plots, not tables). SEMANTICS:
   RBER at fixed read voltages = margin-crossing counts, level-resolved but margin-DEFINED;
   the temperature distortion is measured directly against a real one-year clock.

10. **3D early retention + per-state kinetics** — Luo, Ghose, Cai, Haratsch, Mutlu,
    *Improving 3D NAND Flash Memory Lifetime by Tolerating Early Retention Loss and Process
    Variation*, SIGMETRICS 2018 (open PDF arXiv:1807.05140): 3D MLC charge-trap, 20 °C,
    nine retention points 7 min → 24 days, 0–10K P/E: RBER rises **an order of magnitude in
    ~3 h and another by ~11 days** at 10K P/E; optimal Vc moves 5 voltage steps in the first
    3 h and 5 more by 11 days; Va approximately constant; per-state means follow
    V = A·log(t) + B (Fig. 22, fitted lines in figures). SEMANTICS: layer counts and
    voltage axes deliberately anonymized "to protect" vendor information; charge-trap 3D
    kinetics are log-t from minutes out — a different small-t shape than a single-barrier
    escape, pinned here so the model-side lane meets it with eyes open.

11. **Level-dependence of the mechanisms on 3D (sign structure)** — Malavena et al. (pin 4),
    Figs. 2–4: on FRESH 3D FG arrays, ΔVT after 125 °C bake is **nonmonotonic in level**
    (E-BP), because polysilicon trap **depassivation adds a POSITIVE ΔVT with "a more marked
    dependence on the cell VT level"** while detrapping is negative and grows with level;
    depassivation is "almost negligible" for L1/L2. All ΔVT negative post-cycling.
    SEMANTICS: on 3D, the bake moves levels by TWO mechanisms of opposite sign with
    different level laws — a single-escape retirement order read off bake data would be
    confounded at the mechanism level.

12. **Mielke et al., the canonical JEDEC-class MLC bake paper (CLASS pin)** — Mielke,
    Marquart, et al., *Bit Error Rate in NAND Flash Memories*, IRPS 2008, IEEE doc 4558857.
    Paywalled everywhere located (IEEE, academia.edu, RG); abstract-level content only:
    MLC parts from four manufacturers, raw error data, UBER "a strong function" of P/E
    cycling and retention time. PINNED AS CLASS. **Its per-level failure structure, if any,
    is not publicly readable** — noted as a gap, not quoted.

## D. WHAT DOES NOT EXIST PUBLICLY (findings, per the commissioning note)

- **No located public source gives level-resolved BAKE retirement times** (t at which a
  level pair merges or a level is retired, per level, margin-free). Drive/component-level
  JEDEC qual results are UBER/margin/ECC-defined (pins 1, 12); academic level-resolved data
  is room-temperature (pins 7, 8, 10) or bake-distorted-by-mechanism-mixing (pins 9, 11).
- **No located public source gives the absolute calibration** the parameter-free drop-time
  formula t*_i = f0⁻¹ e^{(B_i−dE_i)/kT}/(1+e^{−dE_i/kT}) needs: normalized-Vt→eV conversion
  (dE_i), per-level barriers B_i, and attempt frequency f0 are all vendor-proprietary, by
  the sources' own statements (pins 4, 7, 10).
- **The first staircase drop is outside every public window**: at RT/2,000 P/E no adjacent
  pair merges within 1 year (instrument §2c) — public data shows the approach, not a step.

## E. EXECUTABILITY VERDICT — falsifier (i) on public data

**PARTIALLY EXECUTABLE.** Split by what the falsifier asks:

- **Retirement ORDER and approach SHAPE: executable now.** Pins 7–9 give margin-free,
  level-resolved drift (numeric at pin 7) and one-year real-time per-voltage error series
  (pin 9). A model-side lane can test ordinal/shape predictions (which pair closes first,
  drift-vs-level monotonicity, the two-sided small-t form) against these without any
  proprietary number. Caveats that must ride along: planar TLC/MLC shows clean
  TAT-dominated high-level ordering (pins 7, 8); 3D adds an opposite-sign second mechanism
  that breaks level-monotonicity (pin 11) and log-t early kinetics (pin 10).
- **Absolute staircase POSITIONS t*_i: not executable on public data as found.** The
  needed calibration (B_i, dE_i in eV, f0, absolute volts) is withheld by every source, and
  level-resolved bake retirement times are unpublished. Executing this half needs either
  vendor data under NDA or an in-house read-retry characterization (a purchasable FPGA
  platform + retail chips — the pin-7/pin-9 methodology is fully described and replicable).
- **The bake clock itself is contested at the source level**: the JEDEC class runs on one
  1.1 eV detrapping Ea (pins 2, 6) while JESD218B's own §3.19 (pin 3), HPCA 2015 (pin 8),
  and the HotStorage measurements (pin 9) each independently deny that this single clock
  reproduces room-temperature level structure. Any comparison of C-86's staircase against
  bake-positioned times inherits this as a first-class systematic, azobenzene-style: the
  convention (which Ea, which temperature law) moves positions by orders of magnitude.

The accumulation directive applies: nothing above was selected for agreement; the
same-function/different-shape reading and the no-gravity-present reading both remain open
on this arm — the data located does not yet reach the regime where the staircase either
appears or fails to.

## NEXT STEP (named, per discipline)

Commission the model-side comparison lane on the ORDER/SHAPE half only: map C-86's
record-pair structure onto adjacent-level pairs, derive its ordinal predictions (first-closing
pair, drift-sign structure, small-t form), and score them against pin 7's transcribed table
and pin 9's one-year series — absolute-position claims explicitly out of scope until a
calibrated dataset exists. In parallel, decide whether an in-house read-retry census
(pin 7/9 methodology, retail chips) is worth the bench cost, since it is the only located
route to margin-free absolute staircase positions.
