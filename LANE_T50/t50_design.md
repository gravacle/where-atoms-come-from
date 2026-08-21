# T-50 DESIGN ONE — THE SAME-READ CONTRAST EXPONENT — V001 — 2026-08-21

LANE: LANE_T50. Directed repair of C-72's PREDICTION/FALSIFIER per the principal's ruling of
2026-08-20 ("stop trying to rescue an absolute single-page ratio") and the registrar's three
computed constraints in the T-50 plan row. Verification: `t50_contrast.py` -> sealed
`t50_contrast_run.txt` (sha256 sidecars in this directory). Nothing outside this directory was
written. No git, no reproduce.sh, no r3.sh.

---

## 1. THE OBSERVABLE

**Name: the SAME-READ CONTRAST EXPONENT.** Dimensionless.

**Preparation (the reader's own, no vendor data).** On ONE part, the reader prepares an
INTERLEAVED partition by the part's own write/erase operations: sectors of equal length L_s
cells, alternating WRITTEN (programmed with a declared data pattern; magnetic: data-written)
and UNWRITTEN (NAND: erased; magnetic: AC-erased by procedure). On the magnetic surface the
interleave carries THREE classes: DATA-written, DC-SATURATED (one-way written), AC-ERASED —
so the surface carries its own positive and negative branches in one read (D-15). The erase
map is the reader's own preparation record.

**Read.** ONE read pass, one instrument, one configuration. v_i = the per-cell RAW read value
in the instrument's own units:
- NAND: per-cell threshold voltage from read-retry / margin-read (volts or DAC counts);
- magnetic: scanning-magnetometer map value at bit-cell resolution (MFM / scanning Hall / NV).
No conversion, no calibration, no reference subtraction.

**The contrast.** For block size N on a geometric grid spanning >= 1.5 decades (every block
wholly inside ONE sector; K >= 8 disjoint pairs per N, paired by sector ADJACENCY):

    D_k(N) = sum_{i in W-block k} v_i  −  sum_{i in U-block k} v_i

Any per-cell term common to the two sectors in the read — reference choice, offset, drift slow
on the pair scale — enters both sums as the same +cN and CANCELS BY CONSTRUCTION. There is no
V_t,neutral, no C_fg, no datasheet quantity anywhere in the definition.

**The statistic (constraint 2).** A_WU(N) = median_k |D_k(N)| — the first-absolute-moment
statistic; equivalently the UNCENTRED median_k D_k(N)^2, whose exponent is 2β. The CENTRED
variance is EXCLUDED: it scales as N for both encodings and discriminates nothing (measured
in the sealed run, same table).

**Same-read controls (D-15, D-8).** A_UU(N) from disjoint unwritten-unwritten pairs and
A_WW(N) from written-written pairs of the SAME read — the noise floor carried beside every
fit, and the demonstration that the signal is the sector-mean difference, not within-sector
fluctuation.

**The observable is the fitted pair of log-log slopes** (β_WU ± SE, β_UU ± SE) over the
declared grid, plus ξ = β_WU − β_UU.

**Density fixed IN THE DEFINITION (constraint 3).** Blocks subsample sectors, so the written
density at every N is the sector's own — a fixed-length record in a growing block is not a
configuration of this observable. Enforced measurably: the read's own 0/1 decisions give the
per-block programmed fraction f̂; the observable is DEFINED only while mean f̂(N) shows no
trend in N (|OLS slope of f̂ on log10 N| <= 0.05 — 5x above honest fit noise, 4x below the
counterexample's smallest trend; measured in the sealed run). A pattern that populates fewer than half
the blocks at some N reads as unwritten there (measured outcome).

**Guard (self-measured void condition).** A written sector READS AS UNWRITTEN — no clause
applies — when median|D_WU(N_min)| <= 8 × A_UU(N_min) at every admissible N_min leaving
>= 1.5 decades (>= 6 grid points). The guard is computed from the same read; no datasheet
threshold exists anywhere. κ = 8 was MEASURED into place, not decreed: the sealed run's seam
scan (R5) shows κ = 5 admits a guard-passing sector below the accumulation band (worst β
0.8893) while κ = 8 is the smallest scanned value that does not (worst β 0.9240).

## 2. WHY THIS HONOURS THE THREE CONSTRAINTS

1. **Contrast AND exponent, not alternatives.** The registered quantity is the exponent OF the
   written-minus-unwritten contrast. The common-mode immunity is exact and structural (cN − cN),
   verified in the sealed run by an offset sweep (including the registrar's 0.5 e/cell case and a
   read-long drift) with the corrupted RAW single-sector exponent shown in the same table.
2. **E|Q| / uncentred Q².** Primary statistic median|D| (exponents 1 vs 1/2); equivalence table
   for uncentred D² (2 vs 1); centred variance shown non-discriminating (≈1 for both) beside it.
3. **Density held fixed in the definition** — blocks subsample sectors; the density-trend check
   has a real failing branch and catches refuter A's counterexample (mutation M3, flip rate in
   the sealed run).

## 3. PREDICTION (the C-72 GROUNDED cell clause)

On an OCCUPANCY-ENCODED surface (one-carrier write; NAND page against the same part's erased
sectors): β_WU = 1 within [0.9, 1.1] for EVERY guard-passing data pattern — the exponent is
pattern-INDEPENDENT; the prefactor (device, density, units) is free — with the same-read
control β_UU = 1/2 within [0.35, 0.65] and excess ξ >= 0.25.

On an ORIENTATION-ENCODED surface (two-signed write; magnetic medium): data-bearing sectors
SCREEN — β_WU sits in the control band with ξ < 0.25, for random AND DC-free-coded data —
while the SAME part's DC-saturated sectors in the SAME read give β_DC-U = 1 within [0.9, 1.1]
(the within-read positive control: the instrument demonstrably sees accumulation).

SCOPE: NAND — the read must return the written pattern (BER < 1%; the read's own instrument
control). Magnetic — the unwritten sectors must be AC-erased BY PROCEDURE; a DC or biased
erase is one-way writing of that sector, and the claim then predicts accumulation (exponent 1)
for its contrast as well — the outcome is measured either way (sealed run R3-M2: a +0.25/grain
erase bias fires clause (c) on 34/50 seeds if this scope clause is ignored — the clause is
load-bearing and measured, not decorative).

## 4. FALSIFIER

All clauses: K >= 8 pairs per point, >= 8 points, >= 1.5 decades, and the offending inequality
must hold by more than 2 fitted SE. Every clause carries the SAME guard, density check and
scope conditions as the prediction — no clause reaches into the guarded-void region.

(a) An occupancy-encoded part — guard-passing, density-trend-free, BER < 1% — whose β_WU
    falls below 0.9, or whose ξ falls below 0.25. (Accumulation absent where the write
    mechanism is one-carrier.)
(b) Two guard-passing data patterns at declared densities on the SAME occupancy part whose
    β_WU differ by more than 0.2. (Pattern-dependence of the exponent.)
(c) An orientation-encoded medium, AC-erased by procedure, whose data-vs-unwritten ξ >= 0.25
    WHILE the same read's DC-saturated control shows β_DC-U >= 0.9. (Accumulation where
    screening is predicted, with the instrument proven able to see accumulation.)

An unwritten-unwritten control leaving [0.35, 0.65], or a failed DC positive control on the
magnetic side, marks the READ INCONCLUSIVE — recorded, not fired (accumulation directive:
outcomes are measured).

## 5. WHAT A READER NEEDS — COMPLETE LIST

1. An occupancy part: any NAND device with per-cell analog read access (characterization
   tester, open-channel controller, or vendor read-retry/margin-read commands — published,
   vendor-independent technique). FLAG: plain consumer 0/1 read access does NOT suffice for
   the analog observable; it supplies only the density check and BER control.
2. An orientation part: any magnetic recording medium plus a writer (spin-stand or drive) for
   the DATA / DC-SATURATED / AC-ERASED interleave, and a bit-cell-resolving scanning
   magnetometer (MFM is standard). VSM is NOT used: no saturation step, no whole-sample
   integral (the two defects that killed the previous magnetic protocol).
3. The declared analysis constants, all in the registered text: the N-grid, K >= 8, guard
   κ = 8, density tolerance 0.05, the bands [0.9,1.1] / [0.35,0.65], ξ threshold 0.25, the
   2-SE rule.
4. NOTHING from any datasheet: no V_t,neutral, no C_fg, no Q_p, no Δ, no ε, no squareness
   class. The cross-encoding comparison is of dimensionless exponents, so no cross-instrument
   calibration and no matched N between surfaces is required.

## 6. THE 32 RESIDUAL DEFECTS

Tested one by one in the sealed output (`t50_contrast_run.txt`, section DEFECT TABLE), each
ADDRESSED with the run that shows it or CONCEDED with the mitigation named. Summary of the
three disqualifying groups: (a) baseline — eliminated structurally, verified by offset/drift
sweep; (b) clauses firing on correct physics — one shared self-measured guard on every clause,
density fixed in the definition, scope conditions explicit, false-fire rate measured; (c) the
orientation identification — squareness, ε, VSM all removed; the magnetic observable is the
same map contrast with its own within-read positive control.

## 7. FOR THE REGISTRAR (T-50 part b) — NOT THIS LANE'S WRITE SCOPE

This lane's R1/R3 suite (baseline booleans + mutation flip rates measured per seed) is a
candidate repointing target for the gate cells of P-FORM-9 and P-ROLES-2 (INST-18): every
decision boolean has a measured failing branch (M1 two-signed write, M3 density-falling), the
non-zero-mean residual mutation is shown STRUCTURALLY NEUTRALIZED by the contrast (invariance
measured, not asserted), and no check compares against a literal. The five orientation literal
checks in model/validate_geometry.py (COMP-11) remain the registrar's to replace; this design
does not depend on them.

## 8. NEXT STEP NAMED (no route closes without one)

Commission the instrument half: run the protocol of §5 on one real NAND part via read-retry
per-cell V_t — the model-side pipeline in this lane consumes a per-cell value array and needs
no change. The census flag (N_E = 100 e is planar-era) does not touch the design: the
exponent is scale-free, so "tens of electrons" moves the prefactor only.
