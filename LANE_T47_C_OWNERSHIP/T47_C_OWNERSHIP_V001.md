# LANE_T47_C_OWNERSHIP — THE OWNERSHIP AND FALSIFIER CHECK FOR k(t_m)

Date: 2026-08-20. Lane C of T-47 (the C-84 audit's creator step). Assignment: merciless
ownership check of the candidate count law `k(t_m) = #{i : dE_i <= B_i - kT ln(f0 t_m)}`
against the three nearest rivals, with concrete falsifiers per accessible surface.
Fairness rule as in C-84/T-39: **a rival owns a statement only if it already made it** —
adjacency, machinery, or one-derivation-step-away is not ownership, and is reported as such.

Convention map (project_model.RecordSurface, solidity-review-corrected): the model's `E_b`
is the activation energy from the METASTABLE (upper) well; the brief's `B` (barrier above
the LOWER well) = `E_b + dE`, so the brief's `B - dE` **is** the model's `E_b`. The exact
survival crossing is `E_b = kT [ln(f0 t_m) + ln(1 + e^{-dE/kT})]`; the brief's
`delta(t_m) = B - kT ln(f0 t_m)` is the leading form, low by at most `kT ln 2` (attained
exactly at the symmetric corner).

## (1) SHARROCK / STREET-WOOLLEY / NEEL — magnetic viscosity

**Their count objects, precisely.** Street-Woolley (Proc. Phys. Soc. A 62, 562 (1949)):
time-dependent remanent magnetization `M(t) = M_0 - S ln t` — an integral over a barrier
distribution of the survival probability of each moment **in its as-written state**.
Sharrock (J. Appl. Phys. 76, 6413 (1994)): time-dependent coercivity/switching field from
the field-reduced barrier `E_B(H) = K_u V (1 - H/H_K)^n` — the time it takes an applied
field plus thermal activation to REVERSE a moment. Neel (Ann. Geophys. 5, 99 (1949)):
the rate itself, `tau = f0^{-1} exp(E_B/kT)`. Every observable in this family is
**magnetization-weighted survival of the occupied well.**

**Where k(t_m) departs.** The count law's unit is the two-valued record: a pair counts only
if BOTH values are durable to t_m, so the record's dropout time is the MINORITY-well time
`t* = f0^{-1} exp((B - dE)/kT)` regardless of which value is written. The viscosity
observables never form this object: on a biased carrier the favored-branch remanence decays
at the same rate `g_u + g_l` but to the asymptote `tanh(dE/2kT)` — it reads "stable"
forever while the record count is 0 past `t*`. Computed on a declared exchange-bias-class
carrier (B = 60 kT, dE = 10 kT, T = 300 K, f0 = 1e9 Hz, t_m = 10 yr): favored-branch
retained signal 0.999909 as t → ∞; record dead past `t* = 5.18e12 s`. Departure term in
step position: `Delta ln t* = -dE/kT` (here -10.0; time factor 4.5e-5).

**Is the departure already in their formalism under another name? Checked honestly: the
BARRIER TERM yes, the COUNT OBJECT no.** Sharrock's field-reduced barrier and the
Preisach-Arrhenius thermal-aftereffect models (Korman-Mayergoyz lineage; hysterons biased
by an interaction field — IEEE Trans. Magn. 38, 2422-ff class) carry exactly the
asymmetry-reduced barrier `B - dE`, under the names "field-reduced energy barrier" /
"interaction-field-shifted hysteron". What no publication in this family states is the
count: #{units bistable-to-spec} as a function of the retention spec. Their integrals are
all magnetization-weighted as-written survival. The `B - dE` RATE is theirs; the
worst-case-over-values COUNT LAW is not.

**The symmetric corner is theirs outright.** Charap-Lu-He (IEEE Trans. Magn. 33, 978
(1997)) and Weller-Moser (IEEE Trans. Magn. 35, 4423 (1999), already source-pinned in
T-41) state the bit-stability criterion `K_u V / kT >~ ln(f0 t_m)` (35 onset, ~60 for
10-year) and derive density limits from a 10-year retention spec. That IS
`B >= kT ln(f0 t_m)` at `dE = 0` — the dE = 0 corner of k(t_m) was published in 1997-99
and the register must concede it by name.

## (2) MLC FLASH — retention level-merging

**Their objects.** JEDEC-class retention bake (JESD22-A117 / JESD47 family; Arrhenius
acceleration is standard practice — Kioxia/Macronix retention briefs) measures Vt
distributions vs bake; Cai et al. (HPCA 2015, "Data Retention in MLC NAND Flash Memory")
characterize per-level decay and establish that **the mean-Vt slope with retention age is
steeper for higher-voltage states** — the asymmetry ORDERING (which level dies first) is
empirically theirs. Worst-pattern retention testing means engineering PRACTICE already
probes worst-case-over-values durability empirically. Information-theoretic treatments
(capacity/mutual-information vs retention and wear, e.g. dynamic-voltage-allocation
arXiv:1403.4333 and NAND channel-capacity estimation) state a continuous
capacity-vs-retention TRADEOFF.

**Does any published treatment DERIVE the usable-level count as a staircase in the
retention spec from an activation model?** Not found. The search covered retention
relaxation schemes (FAST-2012 class and IEEE 8607232), channel-capacity treatments, and
reliability characterization. What exists: empirical level retirement, continuous
capacity curves, per-level decay slopes. What does not exist as a made statement: integer
level count `k(t_m)` with derived step positions `t*_i = f0^{-1} exp((B_i - dE_i)/kT)` in
`ln t_m` from the carrier's activation parameters. **The staircase LAW is unclaimed on
this surface; the ordering and the practice are theirs.** (Caveat: search-level check,
not an exhaustive full-text sweep — same provisional grade as T-39's S-1 concession.)

## (3) QEC MEMORY-TIME (Alicki et al.)

Their objects: STORAGE-TIME bounds — Alicki-Fannes-Horodecki (2D Kitaev: relaxation time
bounded by a size-independent constant ∝ exp(2Δ/kT)); Alicki-Horodecki^3 (4D Kitaev
thermally stable); Bravyi-Terhal (general relaxation-rate upper bound; no self-correction
in 2D stabilizer memories) — and the COUNT-vs-DISTANCE tradeoff `k d^2 <= c n`
(Bravyi-Poulin-Terhal), which contains no time. The composition (count vs retention time
through d(t_m)) is one derivation step away and **was never made**. Under the fairness
rule: UNCLAIMED. No count-law-grade data class exists on this surface today (experiments
measure the lifetime of a few logical qubits, not a retirement staircase) — the QEC
falsifier is honestly N/A for now.

## (4) THE FALSIFIERS, per accessible surface

**F1 — MLC/TLC flash (data class exists: JEDEC retention-bake Vt distributions).**
Step 1: from SHORT bakes only, extract per-level effective activation (B_i, dE_i-analogue
= level-dependent barrier lowering from field-enhanced leakage) and f0 — extraction, not
tuning. Step 2: predict the retirement staircase — positions `ln t*_i = (B_i - dE_i)/kT -
ln f0` and the order (highest level first). Step 3: compare against LONG-bake empirical
retirement (adjacent-level overlap crossing the ECC margin). **Falsified if** the observed
retirement times/order deviate from the predicted ln-t staircase beyond the stated
extraction precision. Published datasets of Cai-et-al. class already suffice for a first
pass.

**F2 — asymmetric / exchange-biased magnetic media (data class exists: both-branch
loop + viscosity measurements).** Write BOTH polarities; measure (a) remanence decay per
branch, (b) the writable-bit count (worst-polarity read-back) vs time. Prediction: the bit
count collapses at `t* = f0^{-1} exp((B - dE)/kT)` while the favored branch holds signal
`tanh(dE/2kT)` to times `exp(dE/kT)` longer — for the declared carrier a factor 2.2e4 in
time with asymptote 0.99991. **Falsified if** the writable-bit count tracks the
favored-branch (Sharrock-inferred) lifetime instead of the computed minority-well
lifetime.

**F3 — QEC memories:** no existing data class at count-law grade; stated, not claimed.

## (5) VERDICT: **PARTIAL** — with the owned fragment stated both ways

**Conceded to rivals (they made the statements):** the dE = 0 corner criterion
`B >= kT ln(f0 t_m)` and its density-limit count consequence (Neel 1949 rate;
Charap-Lu-He 1997; Weller-Moser 1999); the asymmetry-reduced barrier `B - dE` as a RATE
(Sharrock 1994 field-reduced form; Preisach-Arrhenius biased hysterons); the flash
retirement ORDERING (higher levels first — Cai et al. 2015) and the worst-pattern testing
PRACTICE (JEDEC); the storage-time bounds and `k d^2 <= c n` (Alicki et al.;
Bravyi-Terhal; BPT).

**Wholly program-owned fragment (no rival made it):** the two-valued-record COUNT OBJECT
(a pair counts only while BOTH values are durable — worst case over values, not survival
of the as-written value) stated as a LAW: `k(t_m)` a decreasing staircase in `ln t_m`
with per-record derived dropout `t*_i = f0^{-1} exp((B_i - dE_i)/kT)`; PLUS the dual-width
unification — the one clause (ii') on the record's own Liouvillian mode yields the
population width `B - kT ln(f0 t_m)` on the diagonal record and the coherence width
`hbar/t_m` on the off-diagonal record (verified numerically below). The unification is
what makes it a derivation rather than a patch, and no rival account contains either half
as a count statement.

## CONTROLS (D-15) and the numeric record — `control_departure.py`, ALL PASS

- **C1** record-mode rate through the Liouvillian = `g_u + g_l` (rel. err 1.3e-16); no
  closed form smuggled.
- **C2** swept-dE durability crossing lands on the derived exact width to 5.8e-15 kT;
  brief's leading form verified within its `kT ln(1+e^{-dE/kT}) <= kT ln 2` correction
  (numerically 0.0000 kT at this crossing, dE* = 19.71 kT).
- **C3** both polarities relax at the SAME rate to DIFFERENT asymptotes
  (+-tanh(dE/2kT)); departure term computed, max err 3.3e-16.
- **C4** symmetric control: flat-staircase corner, `tau = exp(B/kT)/(2 f0)` exact; the
  ownership-detection control FIRES (corner = the owned Charap/Weller-Moser criterion,
  `B >= kT ln(2 f0 t_m)` — note the ln 2 the leading form drops).
- **C5** coherence-type control: off-diagonal record's own mode (sigma+) has
  `|lambda| = sqrt((dE/hbar)^2 + ((g_u+g_l)/2)^2)`; durable inside `dE = 0.5 hbar/t_m`,
  not durable at `20 hbar/t_m` — **hbar/t_m reappears as the off-diagonal width.**
  Error caught during the build, logged per D-16 practice: the first draft evaluated
  sigma_x, which is NOT an eigenmode — its Rayleigh quotient drops the rotation and read
  8.8e-18 for both cases; `spectrum()`'s own docstring caution ("ask for the record's own
  mode") is the fix.

## D-24 AUDIT of this lane

Every number declared or carrier-typical, none tuned: T = 300 K, f0 = 1e9 Hz (Sharrock's
canonical attempt frequency), B = 60 kT (Weller-Moser 10-year class, source-pinned in
T-41), dE = 10 kT (declared example), t_m = 10 yr (retention spec). No chosen tolerance
enters any conclusion: the crossing is bisected to machine precision against the exact
derived form; clause_ii's 1e-9 is grounded.py's documented NUMERICAL floor, not a
physical width. The lane introduces no new width anywhere.

## NEXT STEP (no route closes without one)

Commission the F1 first pass on published Cai-et-al.-class Vt-vs-bake data: extraction
from short bakes, prediction of the retirement staircase, comparison at long bakes — the
cheapest existing-data test of the wholly-owned fragment. F2 is the discriminating
laboratory test if F1 is confounded by trap-spectrum dispersion.

Files: control_departure.py (controls, exit 0 = all pass), this record.
