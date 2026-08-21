# T-50 DESIGN ONE V002 — THE SAME-READ CONTRAST EXPONENT, EIGHT CLOSURES — 2026-08-21

LANE: LANE_T50/V002, commissioned by the T-50 judgment (LANE_T50/JUDGMENT.txt §2, the
eight-point repair list). READ-FIRST MATERIAL CONSUMED IN FULL before construction: the
judgment; refuter A's verification run and its three escalations (VERIFY_A/); refuter B's
Design One verdict and attack run (VERIFY_B/d1_*); the sealed V001 (t50_design.md,
t50_contrast.py); ERRATUM_REFUTED.json (the 32 prior defects); C72_BRIEF_V001.md; the
C-72 GROUNDED cell as it now reads.

BINDING AND HONOURED: the principal's ruling (baseline-free differential; the exponent OF
the written-minus-unwritten contrast); the three computed constraints in the T-50 plan
row; the accumulation directive (every outcome below is a measured state with a name —
nothing is framed as failure against an imported standard); D-1, D-8, D-15. The bar is
unchanged: a physicist anywhere runs this against their own data; every input is on a
datasheet or readable from the part — in this design, NOTHING is taken from a datasheet
at all; every falsifier is triggerable only by the claim being false.

Verification: `v2_pipeline.py` -> sealed `v2_run.txt`; `v2_mutations.py` -> sealed
`v2_mutations.txt` (sha256 sidecars in this directory). Nothing outside LANE_T50/V002 and
the scratchpad was written. No git, no reproduce.sh, no r3.sh, no ledger or register edit.

---

## 1. THE OBSERVABLE — THE CORE, KEPT

Both refuters kept the core, and V002 does not move it.

**Preparation.** On ONE part, the reader prepares an interleaved partition by the part's
own write/erase operations: sectors of equal length **4096 cells** (pinned), alternating
WRITTEN and UNWRITTEN on the **registered address-balanced role map**: sector s takes the
written role iff (s + s//16) is even, so every address class a = s mod 16 appears in both
roles across the part (this map is what makes the rung-1 fixed-pattern estimator exist;
closure C5). At least 256 sectors. The erase map is the reader's own preparation record.

**Read.** ONE read pass, one instrument, one configuration. v_i = the per-cell RAW read
value in the instrument's own units (NAND: per-cell threshold voltage via read-retry /
margin-read; magnetic: the polar-Kerr analyzer-difference signal per resolved cell —
closure C3). The two registered treatments of closure C5 (address-class de-trend from
the same read's unwritten sectors; optionally the all-erased calibration read) are the
only permitted operations on v, and both are built from the part's own unwritten state.

**The contrast.** For block size N on the PINNED grid {16, 32, 64, 128, 256, 512, 1024,
2048, 4096} (closure C7), over ALL disjoint adjacent written–unwritten sector pairs
(≥ 16 pairs on the pinned geometry), FOUR uniform-random in-sector placements per pair
per grid point, blocks wholly inside one sector:

    D_k(N) = sum_{i in W-block k} v_i − sum_{i in U-block k} v_i

Per-cell terms common to the pair cancel by construction (cN − cN). No V_t,neutral, no
C_fg, no datasheet quantity exists anywhere in the definition.

**The statistic (constraint 2) — THE POOL, PINNED.** A_WU(N) = median |D(N)| over the
whole pool (every disjoint pair × placement). The guard, the control fit, the fitted
exponents, and the point certificate all read pools of exactly this construction — one
object, one moment, one population, the C1 discipline applied to the statistic itself.
(A K-subsample statistic beside a pool-level certificate would be two different
objects, and the seam leaks between them: measured during construction — a subset read
fitted 0.8843 while its pool certified — and closed by this pinning.) The centred
variance is EXCLUDED and displayed as the non-discriminator it is (sealed run R9: ~N
for both encodings; the discriminating moments are E|Q|, exponents 1 vs 1/2, and the
uncentred second moment, 2 vs 1).

**Same-read controls (D-15, D-8).** A_UU(N) = median |D| over the pool of ALL disjoint
unwritten–unwritten pairs of the same read (the noise floor, the guard's denominator,
and the fixed-pattern detector), with A_WW displayed beside it.

**The observable** is the fitted log-log slope pair (beta_WU ± SE, beta_UU ± SE) over
the admissible window, plus xi = beta_WU − beta_UU.

**Density held fixed in the definition (constraint 3)** — blocks subsample sectors, so
the written density at every N is the sector's own; the measurable condition is closure
C1 below.

---

## 2. THE EIGHT CLOSURES

Each closure names the defect verbatim from the record, the registered structure that
closes it, and the sealed measurement.

### C1 — The density condition reads the statistic's own moment (kills A-D1)

**The defect.** V001's density condition pinned the MEAN per-block f̂ while the
registered statistic median_k|D_k| reads the MEDIAN block. A legal declared skew mask
(multiplicative cascade; two-level) passed every V001 condition while registered
falsifier (b) fired 30/30 OK reads on correct physics (VERIFY_A, strict grid escalation)
— the constraint-3 confound in a new costume, and the automatic-REFUTED trigger.

**The closure — structural, not a patch.** For a one-carrier write,
median_k|D_k(N)| = N_E · N · median_k f̂_k(N) + O(noise) — a monotone transform. So the
registered condition now pins EXACTLY the moment the statistic reads, over EXACTLY the
population the statistic samples: **the MEDIAN per-block programmed fraction over ALL
in-sector block placements of every written sector** (computed from the read's own 0/1
decisions; deterministic in the declared pattern) **must be trend-free:
|OLS slope on log10 N, N ≥ 64| ≤ 0.02**, else the read is VOID_DENSITY_MEDIAN and no
clause applies. A one-carrier pattern cannot trend the statistic's typical block without
trending the condition's own object — the channel is closed identically, not at a
measured power.

**Co-tuning (the judge's tolerance note).** A sub-tolerance median trend can move the
fitted exponent by at most the tolerance itself: 0.02, five times inside the 0.1 band
half-width. Measured gap (sealed R2): honest patterns worst |slope| 0.0070 over ten
pattern families × 20 seeds (tolerance 2.9x above); smallest skew-mask slope 0.0361
(tolerance 1.8x below); every skew mask VOID at rate 1.000 over 50 seeds; clause fires
on the masks: 0. The window starts at N = 64 because the discrete binomial-median bias
of honest sparse patterns lives at N = 16–32 and is not a density trend (measured).

**The mean condition is REPLACED, not kept — and the reason is registered.** Over the
placement population the mean per-block f̂ equals the sector density identically at
every N: a mean-trend condition has NO FAILING BRANCH (D-8/INST-14 — the defect class
this program removes). The M3 family (fixed record in growing blocks) is caught instead
by the guard and the median condition: READS_UNWRITTEN / VOID_DENSITY_MEDIAN at rate
1.000 (sealed R2, suite MU3).

**The quantifier, re-derived from the new condition (the judge's instruction):** every
occupancy prediction and clause below quantifies over *guard-passing, density-median-
trend-free, control-in-band* reads — nothing else.

### C2 — The DC-free-coded in-read crosstalk voider (kills A-D2, and retires V001's biased-erase scope clause)

**The defect.** A class-correlated read offset (written sectors read with a small
additive non-magnetic signal — topography/charging crosstalk) fired clause (c) 13/40 at
0.10/grain and 34/40 at 0.20/grain on a correctly-screening medium.

**The closure.** The orientation interleave carries FOUR classes: DATA, DC-FREE-CODED,
DC-SATURATED, AC-ERASED — cycle [DATA, U, DCF, U, DC, U]. Physics: the DCF sectors
screen exactly (that is what the code is). Any additive in-read artifact — crosstalk on
written sectors, a biased erase on unwritten ones — enters every written-minus-U
contrast identically, so it drives beta_DCF-U out of the control band. **Clause (c) is
VOID unless the same read's beta_DCF-U sits in the control band [0.35, 0.65]** (state
INCONCLUSIVE_CROSSTALK otherwise). Refuter A verified the detector arms at 0.498 on
honest reads and predicted it voids both measured fire channels; the sealed run measures
it: fire_c 0/40 at every crosstalk level 0.02–0.20 (states INCONCLUSIVE_CROSSTALK),
arming on honest reads with beta_DCF-U median +0.5043 in band, and the SAME detector
catches the biased erase at 1.000 (suite MU12) — V001's scope clause ("AC-erased BY
PROCEDURE") is replaced by a measurement the read makes itself.

### C3 — The orientation half carried by a NAMED M-reading instrument (kills B-K1)

**The defect.** Every V001-named magnetic reader (MFM "is standard", scanning Hall, NV)
maps STRAY FIELD; every stray-field map has transfer exp(−|k|d)(1−exp(−|k|t)) = 0 at
k = 0, so a uniformly magnetized sector's interior is invisible, the DC positive control
can never pass (bDC-U = nan, 0/20 at every standoff), and clause (c) was a falsifier no
named instrument can trigger.

**The closure — the instrument named, and what it measures verified before naming.**
**Polar Kerr microscopy (wide-field magneto-optic Kerr imaging).** The polar Kerr
rotation of reflected polarized light is proportional to the LOCAL out-of-plane
magnetization M_z within the optical spot (surface/penetration-depth weighted) — NOT to
the stray field. Its response at k = 0 is nonzero: the uniform-film polar-Kerr
hysteresis loop — the standard laboratory measurement of M(H) on uniformly magnetized
perpendicular films — IS the demonstration, since a saturated uniform film (pure k = 0)
returns the full Kerr rotation. Perpendicular CoCrPt media are the polar geometry.
v_i = the raw analyzer-difference signal per resolved cell, in the instrument's own
units; the magneto-optic proportionality is a per-read prefactor and cancels from the
exponent.

**The resolution requirement is a preparation choice, not a device constant:** the
bit-cell length must be at least the microscope's declared optical resolution, and the
READER'S OWN WRITER sets the bit length (spin-stand or drive write clock). Writing long
cells (~1 µm and up) is the reader's own act; no datasheet enters.

**Model-side verification (sealed R5):** the Kerr read modeled as grain-averaged local
M through an optical PSF (k=0-preserving) + read noise. DC-sector interior map value
+1.00015 vs AC-erased -0.00007 under Kerr; +0.00007 under the stray transfer (K1's
kill reproduced in the same table). All four exponents land: beta_DATA-U +0.4973,
beta_UU +0.5049, beta_DCF-U +0.5043 (armed), beta_DC-U +1.0000 (medians, 50 reps) —
B3 true 50/50. And the
requirement is SELF-MEASURED: the stray-field transfer applied to the same surface
drives the read to INCONCLUSIVE_DC_CONTROL at rate 1.000 (suite MU9) — a reader on the
wrong instrument class is told so by the read's own positive control, never silently.

**Clause (c) is triggerable (the K1 complaint answered by measurement):** genuinely
accumulating data (75% one-way) through the same Kerr model fires clause (c) at rate
1.000 with all guards green (suite MU11).

### C4 — Per-clause guard scope (kills B-K2)

**The defect.** V001's "every clause carries the SAME guard, density check and scope
conditions" put the orientation prediction's own subject (a screening DATA sector,
guard ratio 1.000) inside the void region, while the sealed pipeline silently exempted
orientation — two honest readers computed different reachable sets, and the magnetic
"density check" named an operation that does not exist.

**The closure.** The guard-scope table (§5) states, per clause, exactly which
conditions guard it. The DATA sector carries NO void guard — its predicted state is
screening. Orientation admissibility is carried by the same read's DC-saturated
positive control (guard margin + accumulation band), the U-U control band, and the C2
voider. The density condition is registered as OCCUPANCY-ONLY; the orientation
analogue is the declared write pattern itself, verified from the reader's own write
record (there being no read-side programmed fraction on a magnetic surface). Code
equals text: the pipeline implements exactly the table.

### C5 — The registered fixed-pattern treatment (kills B-K3 and A-D3)

**The defect.** Sector-scale fixed pattern — the real instrument's dominant systematic
(NAND wordline/layer mean-V_t structure; scan-line offsets) — at 0.1% of the programmed
charge threw 12/30 reads INCONCLUSIVE on a perfectly accumulating part; real parts sit
orders above; de-trending was forbidden by V001's own text; and the sealed fire_a fired
falsifier (a) on the correctly-accumulating part 28/30 (see C6). A-D3: iid per-sector
offsets at 0.5 e/cell made the protocol "INCONCLUSIVE forever."

**The closure — two exactly-specified rungs, both built from the part's own unwritten
state, plus the detector and the named outcome.**

- **RUNG 1 (always applied, same read):** subtract, from every sector of address class
  a, the mean per-cell value of the SAME READ'S unwritten sectors of class a. The
  address-balanced role map (§1) guarantees every class has unwritten observations.
  This removes ADDRESS-STRUCTURED pattern — exactly the named dominant systematics
  (wordline/layer position; scan lines) — with nothing but the read itself. Measured:
  a deck-profile pattern at amplitudes 0.5–50 e/cell (up to 500x the K3 kill level) is
  fully restored to OK at rate 1.000 (sealed R4).
- **RUNG 2 (registered escalation):** erase the whole part, read it (the all-erased
  CALIBRATION READ — "de-trending from the part's own unwritten sectors", all of
  them), then prepare the interleave and read; subtract each sector's calibration
  mean. This removes ANY static per-sector pattern. Measured: iid per-sector offsets
  0.05–1.0 e/cell (10x above K3's fatal level) — rung 2 restores B1 at rate 1.000 with
  worst |beta_WU − 1| ≤ 0.0033 (sealed R4), and applied to honest reads it does not
  distort (max beta_WU shift 0.0019, sealed R4).
- **THE DETECTOR AND THE NAMED OUTCOME:** the same-read U-U control band IS the
  fixed-pattern detector (a per-sector artifact drives beta_UU from 1/2 toward 1 —
  refuter B's own measurement). Its trip is the named state INCONCLUSIVE_CONTROL: an
  insufficient-instrument outcome with a registered remedy (apply rung 2; if the
  control still leaves the band, the read cannot resolve the claim on this instrument
  — recorded, never fired).
- **SCOPE, STATED HONESTLY (the judge's sector-response-uniformity instruction):** a
  per-sector-random offset that is NOT static across the two preparations is
  observationally equivalent, inside one read, to a physical sector-mean difference —
  no within-read treatment can separate them; the control band detects it and the read
  says INCONCLUSIVE_CONTROL. Response GAIN nonuniformity scales pair prefactors; the
  same control catches it when it matters. Both are measured outcomes, not clause fires.

The suite carries the judge's fixed-pattern member: MU7 (iid 0.5 e/cell) flips the
detector at rate 1.000 with fire_a = 0, and rung 2 recovers B1 at rate 1.000.

### C6 — INCONCLUSIVE precedence wired into the fire booleans (kills B-K4's fire half, A-D6)

**The defect.** V001's sealed fire_a never consulted the control band: code and text
disagreed on the load-bearing boolean, and falsifier (a) fired on a correctly
accumulating part 28/30 under fixed pattern.

**The closure.** There is ONE implementation. The read STATE is computed first, in the
registered order — INCONCLUSIVE_RAILED, READS_UNWRITTEN, VOID_DENSITY_MEDIAN,
INCONCLUSIVE_CONTROL, then OK/SEAM (occupancy); INCONCLUSIVE_RAILED,
INCONCLUSIVE_DC_CONTROL, INCONCLUSIVE_CONTROL, INCONCLUSIVE_CROSSTALK, then OK
(orientation) — and no falsifier boolean can be True unless the state is OK or SEAM
(occupancy) / OK (orientation), by control flow, not by convention. The mutation suite
ASSERTS the invariant on every read of every member and aborts nonzero on violation.
Measured: fire_a 0/30 at every fixed-pattern level where V001 fired 26–29/30.

### C7 — Grid, pairing and placements pinned; the seam claim re-derived per reader (kills A-D4, B-K4, B-repair-7)

**The defect.** The registered-minimum reader freedom (K = 8, 8 points, any geometric
grid ≥ 1.5 decades) moved the seam to where the worst guard-passing beta_WU was
+0.8022 — the EVERY-quantifier was false for registered-compliant readers. And kappa=8's
seam property was a property of the lane's own grid AND noise law (Laplace on the
lane's grid: 0.8921) — no fixed constant transfers to the reader's instrument.

**The closure, in two registered moves.**

1. **The grid, pairing and placements are PINNED**: N ∈ {16·2^j, j = 0..8}; ALL
   disjoint adjacent pairs (≥ 16 available), four placements each — the statistic IS
   the pool; sector 4096 cells, ≥ 256 sectors, guard window ≥ 6 surviving points. A
   part that cannot support this geometry is out of scope — stated. The
   reader-freedom attack surface is gone.
2. **Above the void guard (kappa_void = 8, unchanged) there is NO registered constant.**
   The point-band sentence is asserted only for reads whose OWN CERTIFICATE is clean:
   B = 400 surrogate ladders, each a full-size bootstrap replica of the read's own
   written–unwritten pool — THE VERY OBJECT THE STATISTIC REPORTS, carrying the
   read's own noise law AND its own programmed-count shot noise (the seam's dominant
   fluctuation) — each surrogate run through EXACTLY the pipeline's guard and fit
   against the read's own A_UU. One surrogate below the band edge and the read is a
   SEAM read: measured, recorded, the 2-SE falsifier still armed, the point-band
   sentence NOT asserted. Because statistic and certificate read the same object, a
   read whose own ladder tilts below the band cannot certify — its surrogates inherit
   the tilt. Nothing but the reader's own draws enters — this is the judgment's
   "kappa re-derived from the reader's own A_UU and measured noise law", made
   operational.

**Measured (sealed R3):** seam ensembles at f = 0.015–0.10 under uniform, Laplace, and
Student-t3 noise at matched sigma (refuter B's attack-5 laws): the worst CERTIFIED (OK)
beta_WU anywhere is +0.9865 — and NO ensemble read needed the SEAM state: on the pinned
pool statistic a read either reads as unwritten or certifies. Refuter B's +0.8022 and
refuter A's 0.8815–0.8954 guard-passing families were creatures of the K = 8–16 subset
statistic, and the pinning removed them. The SEAM state remains registered and
EXERCISED: a genuinely screening surface cannot certify, lands SEAM, and the 2-SE
falsifier fires there at rate 1.000 (the suite's two-signed member). Falsifier fires
across every honest ensemble: 0 (the 2-SE rule held, as both refuters also measured).

### C8 — The railed-population branch (kills B-K5, A-D5)

**The defect.** A railed erased population (consumer read-retry floors at GND —
custody, cai_procieee ~line 1780) gave A_UU = 0: vacuous guard pass, log10(0) fit,
bUU = nan, no registered branch.

**The closure.** A_UU(16) = 0 is the FIRST state checked: INCONCLUSIVE_RAILED —
insufficient access, nothing fitted, no guard evaluated, no nan reaches any decision.
The registered text tells the reader what suffices: a per-cell analog access mode that
resolves the erased distribution (characterization tester or margin-read with
negative-range support); plain consumer 0/1 access supplies only the density condition
and BER control, not the observable. Measured: rate 1.000, fires 0, no nan escapes
(sealed R6, suite MU8).

---

## 3. PREDICTION (the C-72 clause, re-derived quantifiers)

On an OCCUPANCY-ENCODED surface (one-carrier write; NAND page against the same part's
erased sectors), for EVERY read that is guard-passing, density-median-trend-free, and
control-in-band (the C1 quantifier):

- **the falsifier-grade claim (all such reads, OK and SEAM):** beta_WU does not fall
  below the accumulation band by more than 2 fitted SE, and xi does not fall below
  0.25 by more than 2 SE — the exponent is 1 and the excess over the same-read control
  is real;
- **the point-band claim (certified reads — state OK):** beta_WU = 1 within
  [0.9, 1.1], beta_UU = 1/2 within [0.35, 0.65], xi ≥ 0.25, for every certified
  guard-passing data pattern — the exponent is pattern-INDEPENDENT (prefactor free);
- SEAM reads (certificate not clean) are recorded with their measured beta_WU; the
  point sentence is not asserted there. Measured: no honest seam-ensemble read needed
  SEAM (reads void or certify; worst certified +0.9865) — the refuters' 0.80–0.89
  guard-passing families no longer exist on the pool statistic — while a genuinely
  screening surface lands SEAM with the 2-SE falsifier armed and firing (suite MU1).

On an ORIENTATION-ENCODED surface (two-signed write; magnetic medium read by polar
Kerr microscopy), for reads whose DC positive control, U-U control, and DCF voider are
green: data-bearing sectors SCREEN — beta_DATA-U in the control band with xi < 0.25,
for random AND DC-free-coded data — while the SAME read's DC-saturated sectors give
beta_DC-U = 1 within [0.9, 1.1] (the within-read positive control: the instrument
demonstrably sees k = 0 accumulation).

SCOPE: NAND — the read returns the written pattern (BER < 1%, the read's own control).
Magnetic — the erase artifact and any in-read crosstalk are DETECTED by the C2 voider
(INCONCLUSIVE_CROSSTALK), not policed by a scope clause.

## 4. FALSIFIER

All clauses: the pinned grid and K, the registered conditions of §5 for that clause,
and the offending inequality by MORE THAN 2 FITTED SE. States other than OK/SEAM
(occupancy) or OK (orientation) RECORD and never fire (C6).

- **(a)** An occupancy read in state OK or SEAM whose beta_WU + 2·SE < 0.9, or whose
  xi + 2·SE_xi < 0.25. (Accumulation absent where the write mechanism is one-carrier.)
- **(b)** Two declared data patterns on the SAME part, both reads state OK, whose
  beta_WU differ by more than 0.2 beyond 2·SE of the difference. (Pattern-dependence
  of the exponent.)
- **(c)** An orientation read in state OK (DC control ≥ 0.9 with its guard passed, U-U
  control in band, DCF voider armed) whose xi − 2·SE_xi ≥ 0.25 while
  beta_DC-U − 2·SE ≥ 0.9. (Accumulation where screening is predicted, with the
  instrument proven able to see accumulation in the same read.)

Measured false-fire rate over 600 honest reads through the full V002 machinery: 0
(sealed R7). Measured trigger rate where the named physics IS present: clause (a)
1.000 under the two-signed write; clause (c) 1.000 under one-way data (suite MU1,
MU11) — every clause has a measured firing branch and a measured non-firing honest
ensemble.

## 5. GUARD-SCOPE TABLE (C4 — normative)

| registered condition                       | clause (a) | clause (b) | clause (c) | point claim |
|--------------------------------------------|-----------|-----------|-----------|-------------|
| occ void guard (kappa_void = 8)            | yes       | yes       | —         | yes         |
| occ density-median condition (C1)          | yes       | yes       | —         | yes         |
| occ control band / FP detector (C5)        | yes       | yes       | —         | yes         |
| occ point certificate (C7)                 | —         | OK-only   | —         | yes         |
| occ railed branch (C8)                     | yes       | yes       | —         | yes         |
| ori DC positive control (C3/C4)            | —         | —         | yes       | —           |
| ori U-U control band                       | —         | —         | yes       | —           |
| ori DC-free voider (C2)                    | —         | —         | yes       | —           |

Clause (a) fires in states OK and SEAM (the 2-SE rule is its own protection — measured
0 false fires in every seam ensemble of both refuters and this run); clause (b)
requires both classes OK; clause (c) requires state OK. The DATA sector carries no
void guard: its predicted state is screening (the K2 closure). The density condition
is occupancy-only; the orientation analogue is the declared write pattern, verified
from the reader's own write record.

## 6. INSTRUMENT PROTOCOL

**Occupancy (NAND).**
1. Choose a part with per-cell analog read access (characterization tester,
   open-channel controller, or vendor read-retry/margin-read). The access mode must
   resolve the erased distribution: if A_UU(16) = 0 the read is INCONCLUSIVE_RAILED
   and the text says so (C8).
2. Erase; program the interleave on the registered address-balanced role map with a
   declared data pattern; record the map (the reader's own preparation record).
3. One read pass: per-cell V_t (any units). BER control: the read returns the pattern
   (< 1%).
4. Apply rung 1 (address-class de-trend from the read's own unwritten sectors).
   Compute the pools, the guard, the density-median ladder from the read's own 0/1
   decisions, the control fit, the certificate. If INCONCLUSIVE_CONTROL: apply rung 2
   (all-erased calibration read) and re-run (C5).
5. Report the state; in OK/SEAM report (beta_WU ± SE, beta_UU ± SE, xi ± SE), the
   margin, and the certificate outcome. Clauses per §4/§5.

**Orientation (magnetic, polar Kerr).**
1. Any perpendicular recording medium plus a writer (spin-stand or drive) and a
   polar/wide-field Kerr microscope with declared optical resolution. Write the
   four-class interleave [DATA, U, DCF, U, DC, U] at bit-cell length ≥ the declared
   resolution (the writer's clock — the reader's own choice). AC-erase the U sectors.
2. One Kerr image pass of the track; v_i = raw analyzer-difference per resolved cell.
   Blocks keep ≥ 4 cells from sector boundaries (PSF margin).
3. Compute the DC positive control (guard + band), the U-U control, the DCF voider,
   then DATA-U. States and clause (c) per §5. Scan-line structure: rung 1/rung 2
   apply unchanged (the calibration read is the AC-erased-everywhere image).
4. VSM remains out (saturation destroys the state; whole-sample integral); MFM /
   scanning Hall / NV remain out and are SELF-DETECTED by the DC control (C3) — the
   protocol tells a reader on the wrong instrument that the read is inconclusive
   rather than letting a void pass silently.

## 7. WHAT A READER NEEDS — COMPLETE LIST, EVERY INPUT SOURCED

1. An occupancy part with per-cell analog read access — sourced: the part itself plus
   a published, vendor-independent access technique; the railed branch names the
   failure mode of consumer access (C8).
2. An orientation part: any perpendicular magnetic medium + writer + polar Kerr
   microscope; bit length set by the reader's writer to ≥ the declared resolution —
   sourced: the reader's own preparation.
3. The declared analysis constants, ALL in this text (C7): the pinned grid
   {16..4096, ×2}; the pool statistic (all disjoint adjacent pairs, ≥ 16, four
   placements each); sector 4096, ≥ 256 sectors, kappa_void = 8, ≥ 6 surviving
   points, density-median tolerance 0.02 on N ≥ 64, bands [0.9, 1.1] / [0.35, 0.65],
   xi ≥ 0.25, the 2-SE rule, certificate B = 400, PSF edge margin 4 cells, BER < 1%.
4. The two registered treatments (C5) and the certificate procedure (C7) — both built
   from the reader's own read(s), specified operation-by-operation in this text and
   implemented in the sealed pipeline.
5. NOTHING from any datasheet: no V_t,neutral, no C_fg, no Q_p, no Delta, no epsilon,
   no squareness class, no noise constant, no kappa above the void guard. The
   cross-encoding comparison is of dimensionless exponents; no cross-instrument
   calibration, no matched N.

---

## 8. THE DEFECTS TABLE

Every named defect from BOTH T-50 refuter verdicts, then the 32 prior defects
(ERRATUM_REFUTED.json), each CLOSED (how, with the sealed measurement) or CONCEDED
(why, with the mitigation named). Sealed evidence: v2_run.txt (R-blocks),
v2_mutations.txt (MU-members).

### Refuter A (computation), Design One

| # | defect | disposition |
|---|--------|-------------|
| A-D1 | median-vs-mean density confound: skew masks fire clause (b) 30/30 on correct physics | **CLOSED (C1).** The condition now reads the statistic's own moment over its own population; every named mask (cascade 1.5/0.5 … 1.7/0.3, two-level) is VOID_DENSITY_MEDIAN at 1.000 with 0 clause reaches (R2, MU4–6); honest/mask separation 2.9x / 1.8x around the 0.02 tolerance; sub-tolerance trends bounded at 0.02 of exponent, inside the band. |
| A-D2 | class-correlated read offset fires clause (c) 13/40 @0.10, 34/40 @0.20 | **CLOSED (C2).** The DCF voider: fire_c 0/40 at every level; states INCONCLUSIVE_CROSSTALK; armed on honest reads (R5, MU10). |
| A-D3 | 0.5 e/cell iid sector offsets: 100% INCONCLUSIVE forever | **CLOSED (C5).** Rung 2 restores B1 at 1.000 up to 1.0 e/cell; the un-escalated state is named and never fires; the one-read identifiability limit is stated as scope (R4). |
| A-D4 | kappa=8 seam underpowered: honest guard-passing beta 0.8815–0.8954 | **CLOSED (C7).** The family no longer exists on the pinned pool statistic (reads void or certify); worst CERTIFIED beta anywhere +0.9865; falsifier fires: 0 (R3). |
| A-D5 | A(N)=0 undefined branch | **CLOSED (C8).** INCONCLUSIVE_RAILED at 1.000, no nan, no vacuous guard (R6, MU8). |
| A-D6 | fire_a omits the registered INCONCLUSIVE precedence | **CLOSED (C6).** States before fires by control flow; one implementation; invariant asserted on every suite read (MU-all). |

### Refuter B (instrument), Design One

| # | defect | disposition |
|---|--------|-------------|
| B-K1 | no named stray-field instrument can trigger clause (c); orientation half has no runnable side | **CLOSED (C3).** Polar Kerr microscopy named, k=0 response verified and stated; DC interior +1.00015 on the Kerr model vs +0.00007 on the stray transfer (both displayed); clause (c) fires at 1.000 on genuinely accumulating data (MU11); the stray-field instrument is self-detected at 1.000 (MU9). |
| B-K2 | shared guard voids the orientation prediction's own subject; magnetic density check undefined | **CLOSED (C4).** Per-clause guard-scope table; DATA carries no void guard; density condition registered occupancy-only with the write-record analogue named. |
| B-K3 | fixed pattern at 0.1% of programmed charge: INCONCLUSIVE forever + fire_a on correct physics 28/30 | **CLOSED (C5 + C6).** Two registered rungs from the part's own unwritten state; detector + named state; fire_a 0/30 at every level; rung-2 recovery 1.000; suite member at 1.000 (R4, MU7). |
| B-K4 | worst guard-passing beta +0.8022 at registered minima; EVERY-quantifier false | **CLOSED (C7).** Grid/K pinned (the minima no longer exist); the point sentence covers exactly the certified set; the falsifier's 2-SE form covers the rest and measured 0 false fires (R3). |
| B-K5 | railed erased population: undefined protocol on the named access mode | **CLOSED (C8).** (R6, MU8.) |
| B-r7 | kappa's seam property is noise-law-dependent; no constant transfers | **CLOSED (C7).** No constant above the void guard remains; the certificate is built from the reader's own draws; verified under uniform, Laplace, t3 (R3). |

(B's repairs 1–6 are K1–K5 restated plus the precedence repair; each is covered above.)

### The 32 prior defects (ERRATUM_REFUTED.json: computation C-0..13 = COMP-1..14, instrument I-0..17 = INST-1..18)

V001 closed these and both refuters re-verified the closures that matter (the three
binding constraints, the anchor, the moment table). V002 carries every V001 closure
forward — the observable, the contrast, the statistic and its excluded moments, the
unit-freedom, the anchor — and repairs the four the refuters showed V001 had re-opened.
Dispositions, one line each:

| defect | V002 disposition |
|---|---|
| COMP-1 (N-decay clause unguarded in f) | CLOSED — density fixed in the definition + C1 median condition + guard; M3 family voided at 1.000 (R2, MU3). |
| COMP-2 (void condition vs falsifier set) | CLOSED — one state machine, per-clause scope table (C4, C6); the orientation-side recurrence (K2) closed. |
| COMP-3 (clauses disagree on reachable sets) | CLOSED — the scope table IS the reachable-set statement; code equals text. |
| COMP-4 (rho ≠ M_r/M_s under tilt) | CLOSED — no squareness anywhere; 30° tilt moves nothing (R5 tilt row: prefactor). |
| COMP-5 (baseline not available) | CLOSED — structural: cN − cN; R8 sweep incl. ±0.5 e and 50 e with the raw shadow displayed. |
| COMP-6 (datasheet supplies nothing) | CLOSED — zero datasheet inputs (§7.5). |
| COMP-7 (nominal Q_p) | CLOSED — no Q_p; prefactors free (R10: exponent flat f=0.10–0.90). |
| COMP-8 (floor slack artifact) | CLOSED — no floor exists; mu=2/20 invariance measured with shadow (MU2). |
| COMP-9 (closed-form noise constants) | CLOSED — noise scale is the read's own A_UU pool; no constant registered (and C7 removed the last one above the void guard). |
| COMP-10 (checks are theorems about the generator) | CLOSED — every decision has a designated killer at 1.000 (MU matrix); the one unfailable candidate (a mean-density condition) was identified and REMOVED for exactly this reason (C1). |
| COMP-11 (orientation literals in model/) | FLAGGED, UNCHANGED — outside this lane's write scope; the registrar's T-50(b) landing repointed the occupancy gate; the five orientation literals remain the registrar's (named, not hidden). |
| COMP-12 (mixed statistics) | CLOSED — matched ensembles, medians with IQR throughout (R1, R5). |
| COMP-13 (one shared residual field) | CLOSED — every replicate its own field (all R-blocks). |
| COMP-14 (0.00096 vs 0.0027) | CLOSED — neither number is used; all numbers computed fresh under declared seeds. |
| INST-1 (V_t,neutral circular) | CLOSED — no neutral reference exists in the definition. |
| INST-2 (Delta spread/offset conflation) | CLOSED — no Delta; sector-mean physics joins the signal; instrument sector offsets handled by C5. |
| INST-3 (C_fg hunt) | CLOSED — unit freedom measured to 1e-10 (R9). |
| INST-4 (epsilon unfalsifiable) | CLOSED — no epsilon. |
| INST-5 (rho = 1−2eps not physics) | CLOSED — no such identification; tilt row (R5). |
| INST-6 (borrowed squareness class) | CLOSED — no borrowed figure anywhere. |
| INST-7 (datasheet claim false) | CLOSED — §7.5; N_E is model-side only; scale-free exponent. |
| INST-8 (VSM destroys the state) | CLOSED — Kerr read is non-destructive, no saturation step, no whole-sample integral (C3). |
| INST-9 ("one instrument" spans two) | CLOSED — one instrument per surface, stated; cross-encoding comparison is of dimensionless exponents. |
| INST-10 (clauses are corollaries) | CLOSED — three clauses, three prepared surfaces, scope table; no floor so no no-N corollary. |
| INST-11 (no ensemble/confidence) | CLOSED — pinned K, grid, 2-SE rule; false fire 0/600 measured (R7). |
| INST-12 (DC-free coded wrong value) | CLOSED — DCF is now a REGISTERED CLASS with its own role (the C2 voider) and lands in the control band (R5). |
| INST-13 (sqrt(2/piN) constant) | CLOSED — as COMP-9. |
| INST-14 (floor check cannot fail) | CLOSED — no floor; and the same discipline removed the mean-density condition (C1). |
| INST-15 (control A equality half) | CLOSED — no equality-to-literal control; all controls are populated branches with variance; the certificate replaces the unfailable-constant pattern. |
| INST-16 (18% mutation power) | CLOSED — every designated cell at 1.000 (MU matrix). |
| INST-17 (two documents disagree) | CLOSED — one pipeline computes every number once; the design quotes the sealed txt. |
| INST-18 (orphaned gate cells) | CLOSED BY THE REGISTRAR at the T-50(b) landing (R8 repointing map in the sealed GATE); V002 adds the two judge-carried members to ITS OWN suite and flags them as the gate's appendable members (judgment c1) — the gate append remains the registrar's. |

---

## 9. HONEST RISKS (named, not narrated away)

1. **The Kerr instrument model is a model.** Its k=0 claim is standard magneto-optics
   (uniform-film polar loops), but the sealed run verifies a MODEL of local-M imaging
   (grain averaging + PSF + additive artifacts), not a Kerr microscope. Refuter B said
   V002 "would face me on the Kerr-instrument raw-value law": the raw-value law
   registered here is v_i ∝ local mean M_z + additive instrument terms, with the
   additive terms detected by C2/C5 and the proportionality a prefactor. A real Kerr
   channel with a MULTIPLICATIVE state-correlated artifact (reflectivity change on
   written tracks entering the Kerr channel) would evade the additive voider —
   partially mitigated by the DCF class (any state-correlated artifact rides DCF too,
   whatever its algebra, if it is class-level); a cell-level multiplicative artifact
   correlated with the DATA pattern itself remains the named residual exposure.
2. **The certificate is a bootstrap.** Because it resamples the very object the
   statistic reports, a read whose own ladder tilts below the band cannot certify (the
   surrogates inherit the tilt) — the leak this design measured and closed during its
   own construction, when a K = 16 subset statistic beside a pool certificate let a
   0.8843 read certify. What remains is finite-B resolution: a read whose true
   below-band resampling probability is below ~1/B can certify and land below the
   band at that residual rate. The falsifier is protected regardless (2-SE, 0 fires
   measured everywhere); the exposure is confined to the point-band sentence and
   shrinks as the reader raises B (400 is a floor, not a ceiling).
3. **Rung 2 assumes the fixed pattern is static across two preparations.** Drift
   between the calibration and analysis reads at sector scale re-enters as an
   un-removed offset; the control band catches it (INCONCLUSIVE_CONTROL), so the
   failure mode is a lost verdict, never a false one — but a part whose FP drifts is a
   part this protocol cannot conclude on, and that is stated.
4. **The orientation half's crosstalk voider is conservative.** It voids at 0.02/grain
   crosstalk (measured), so clause (c) is reachable only on reasonably clean Kerr
   reads. The cost of closing A-D2 structurally is that dirty instruments yield
   INCONCLUSIVE_CROSSTALK rather than verdicts — the accumulation directive prefers a
   named non-verdict to a false fire, and the alternative (a tolerance) is exactly the
   V001 defect class.
5. **Same-lineage exposure.** This lane constructed both the closures and their
   verification; the commissioning rule stands: V002 faces BOTH refuters again, with
   this design, both verdicts, ERRATUM_REFUTED.json, and the judgment as their
   read-first material, before anything is registered into C-72.

**NEXT STEP NAMED (no route closes without one):** stage LANE_T50 VERIFY_A2 / VERIFY_B2
— the two refuters re-commissioned against this V002 (computation lens; instrument
lens), read-first material as above; on survival, the registrar lands the C-72
PREDICTION/FALSIFIER rescope from §3–§4 and appends the two judge-carried members to
the sealed gate's suite (judgment caveat c1).
