# DWELL SEARCH — ARM 1: SPINTRONICS — SOURCE-PINNED (2026-08-21)

Target: OPEN external data for the two record laws
- **C-69 LIFETIME**: tau = exp(E_b/kT) / (2 f0 cosh(dE/2kT)) — in dwell language 1/tau* = 1/tau_up + 1/tau_dn
- **C-70 STEADY**: <R>_ss = tanh(dE/2kT), equivalently tau_up/tau_dn = exp(-dE/kT) by detailed balance

The object sought: thermally activated TWO-STATE telegraph with PUBLISHED asymmetric dwell times
(tau_P AND tau_AP or their distributions), known temperature, best of all swept in T or in a bias
playing dE's role. Benchmark of what we want: the RESTRICTED Funatsu Zenodo deposit (6767828).
House style per LANE_T41_EXTERNAL/CITATIONS.md and LANE_C86_EXTERNAL/*: source, value,
uncertainty, SEMANTICS, ACCESS STATUS verified by actually loading the resource. Every access
claim below states the act that verified it.

Local extracted texts (for re-verification): scratchpad/bielefeld_1908.02139.txt,
kaneko_apex2024.txt, koenders_prapplied2025.txt, safranski_2010.14393.txt,
vodenicarevic_1706.05262.txt, rehm_2209.01480.txt, coupling_2307.15165.txt.

---

## PINNED — OPEN RAW DATA (the one open deposit found)

1. **Raw SMTJ telegraph time series, single condition** — DATUM, raw data OPEN (verified by
   download). F. Kammerbauer & L. Schnitzspan (JGU Mainz), *Time series sMTJ switching data*,
   **Zenodo 15222667**, doi:10.5281/zenodo.15222667 (2025), **CC-BY-4.0**. 28 parquet chunks of
   1e8 rows (~100 MB each, ~2.8 GB total) + 206-byte header, columns (timestep, voltage), series
   name `#140_r1c3p14_timeseries_0.1V`. Header verbatim content: sampling time 5.0e-6 s
   (200 kHz), 1.00e3 acquisitions, measured 2024/02/05, applied voltage 1.000000e-3 with
   "10k Ohm pre-resistance" (filename says 0.1 V — source-vs-junction voltage bookkeeping left
   ambiguous by the header), Keithley 2400 + Tektronix DPO 7354.
   ACCESS: record JSON loaded via Zenodo API (access: open); header file downloaded and read;
   chunk 0 verified downloadable by HTTP HEAD → **200 OK, content-length 100,213,314**.
   COMPANION PAPER (defines the device): K. Koenders, L. Schnitzspan, F. Kammerbauer, S. Shu,
   G. Jakob, M. Kläui, J. H. Mentink, N. Ahmad, M. van Gerven, *Noise-based local learning using
   stochastic magnetic tunnel junctions*, **Phys. Rev. Applied 23, 054035 (2025)** (Editors'
   Suggestion); dataset cited as its Ref. [41] — verified by loading the publisher-version PDF
   from the Mainz OA repository (openscience.ub.uni-mainz.de) and reading the reference verbatim;
   also arXiv:2412.12783. Device: **in-plane** low-uniaxial-anisotropy circular SMTJ, **60 nm**
   diameter, stack Ta(10)/Ru(10)/Ta(10)/PtMn(20)/CoFe(2.2)/Ru(0.8)/CoFeB(2.4)/MgO(1.1)/
   CoFeB(3.0)/Ta(10)/Ru(30), TMR > 100 %, RA ≈ 15 Ω·µm²; paper's own noise traces taken for
   2500 s at 40 kHz under in-plane fields ≈ −3 to +3 mT (three bias conditions).
   SEMANTICS: this is the only OPEN raw equivalent of the restricted Funatsu deposit found
   anywhere — but it is ONE condition (one bias, no field sweep in the deposit; deposit sampling
   200 kHz differs from the paper's 40 kHz figure data, so the deposit is a separate longer/faster
   run on the same device family). **Measurement temperature is stated nowhere** (ambient
   implied); any kT-bearing number extracted from it inherits an assumed T ≈ 295 K. On a single
   condition the occupancy-ratio = dwell-ratio identity is definitional, not a physics test.
   EXECUTABILITY (C-69/C-70): input-hygiene grade, and genuinely valuable at that grade —
   two-state verification, exponentiality of BOTH dwell distributions, stationarity over ~1.4e4 s,
   and the both-values factor 1/tau* = 1/tau_P + 1/tau_AP measured directly from raw data with
   ~1e6+ transitions. Knowns: sampling grid, bias, device class. Unknowns: T (assumed ambient),
   field during the deposit run, per-state calibration (from voltage histogram). No sweep → no
   cosh/tanh curve from this deposit alone.

## PINNED — OPEN FULL-TEXT SOURCES WITH TWO-SIDED DWELL DATA (figure-grade, text-anchored)

2. **All three handles (T, H, V) on one device family, both dwell times** — DATUM at
   table/figure grade, full text OPEN. G. Reiss, J. Ludwig, K. Rott (Bielefeld),
   *Tuning superparamagnetism in perpendicular magnetic tunnel junctions*, **arXiv:1908.02139v3**
   (2019; journal publication not located — treat as preprint). Perpendicular CoFeB(1.1 nm)
   free layer, exchange-biased reference, circular 140 nm pillars, MgO 1.2/1.4/1.6 nm.
   Measures tau_P AND tau_AP separately: vs perpendicular field at MULTIPLE temperatures
   (Fig. 2; example condition 48 °C, −190 Oe, 100 mV, ms-range switching), vs bias voltage at
   FIVE fields at 65 °C (Fig. 3a); exponential dwell histograms shown (Fig. 1b). Text/table
   anchors for digitization: **Table I** (VE/A·MS = 9.28/8.33/6.29 A·nm², ΔE = 2.5/1.3/2.4 eV,
   K* = 26/13.5/24.9 kJ/m³ for tMgO = 1.2/1.4/1.6 nm); spin-torque slope
   Δ− = −(20 ± 2) pJ/V × exp(−(13 ± 6)/nm · tMgO), anisotropy term β = (30 ± 15) fJ/V·m;
   Gaussian tuning-curve anchors H0 = 157.8 Oe (σ = 6.29 Oe) at +80 mV and H0 = 180.6 Oe
   (σ = 6.95 Oe) at −70 mV; entropic intercept ln w ≈ 35.
   ACCESS: full text extracted from the arXiv PDF (local copy scratchpad/bielefeld_1908.02139.txt).
   No dataset deposit; raw numbers live in figures, anchored by the above.
   SEMANTICS (three traps, all load-bearing): (a) the measured prefactor is **tau_0/w with
   ln w ≈ 35** (entropy of ~10^15 switching pathways) — C-69's f0 is NOT a bare FMR attempt
   frequency here; (b) the 140 nm electrode is NOT a macrospin — activation volume radius
   ≈ (17 ± 4) nm; two-valuedness at readout level still holds (exponential dwells, two resistance
   states); (c) the field axis must be shifted by H_comp (stray-field coupling to the reference)
   before H' plays dE's role. Their Eq. 3 decomposition is literally the record laws' split:
   d/dU [ln tau_P + ln tau_AP] = barrier part (C-69's symmetric factor),
   d/dU [ln tau_P − ln tau_AP] = asymmetry part (C-70's ratio) — the FIELD version of the ratio
   is clean Zeeman/Boltzmann; the VOLTAGE version is spin-torque, a nonequilibrium drive, NOT
   detailed-balance dE (see semantics control #9).
   EXECUTABILITY: strongest open figure-grade object for BOTH laws. C-70: digitize Fig. 2,
   slope of ln(tau_P/tau_AP) vs H' at each T → test exp(−dE/kT) with dE = 2·mu0·(VA·MS)·H'
   against Table I's VE/A·MS, T known per curve. C-69: symmetric part
   sqrt(tau_P·tau_AP) vs 1/T at H' = 0 → activation E_b with the entropic prefactor as the
   registered deviation from a naive f0. Knowns: T per curve, Table I anchors. Unknowns: raw
   dwell values (digitization), w's T-dependence.

3. **The benchmark's OWN open-access paper (data restricted, figures open)** — DATUM at figure
   grade. T. Funatsu, S. Kanai, J. Ieda, S. Fukami, H. Ohno, *Local bifurcation with
   spin-transfer torque in superparamagnetic tunnel junctions*, **Nat. Commun. 13, 4079 (2022)**,
   doi:10.1038/s41467-022-31788-1 — paper OPEN ACCESS, verified loaded at **PMC9283488**
   (Europe PMC flags it open). Data Availability VERBATIM: "The data that support the plots
   within this paper have been deposited in Zenodo at https://zenodo.org/record/6767828".
   **That deposit is RESTRICTED** — verified TODAY via Zenodo API record 6767828: access
   status "restricted", no files visible (title: "Dataset: Local bifurcation with spin-transfer
   torque in superparamagnetic tunnel junctions", Kanai/Funatsu/Ieda/Fukami/Ohno, 2022-06-28).
   This confirms the C-86 refuter's catch; the words "deposited in Zenodo" in the paper do NOT
   mean open.
   Content: perpendicular CoFeB/MgO SMTJs, device A **34 nm** (TMR ~73 %, RA 5.5 Ω·µm²),
   device B 28 nm (TMR ~74 %, RA 8.1 Ω·µm²); **tau_P and tau_AP separately vs perpendicular
   field Hz and vs bias** (expected switching times vs Hz in Fig. 4e; thermal stability factors
   Δ_P, Δ_AP vs (Hz, V) in Fig. 5a); dwell times ~0.3 ms to seconds; **tau_0 = 1 ns ASSUMED,
   not measured** — every Δ inherits ln f0.
   EXECUTABILITY: C-70 clean via the Hz sweep (Zeeman dE, detailed balance intact): tanh
   occupancy and exp ratio testable by digitizing Fig. 4e, anchored by the stated device
   parameters and dwell range; C-69 via 1/tau* vs Hz cosh shape at fixed T. Temperature not
   stated as a number in main text (room ambient; the twin APEX paper below carries the
   T-semantics). Raw-grade rerun requires the restricted deposit → that route stays with the
   benchmark.

4. **Temperature-SWEPT lifetime law on perpendicular SMTJs** — DATUM at figure grade, OPEN
   (CC-BY). H. Kaneko, R. Ota, K. Kobayashi, S. Kanai, M. Elyasi, G. E. W. Bauer, H. Ohno,
   S. Fukami, *Temperature dependence of the properties of stochastic magnetic tunnel junction
   with perpendicular magnetization*, **Appl. Phys. Express 17, 053001 (2024)**,
   doi:10.35848/1882-0786/ad43b0. ACCESS: full text extracted from the IOP OA PDF (local copy
   scratchpad/kaneko_apex2024.txt). No data deposit, no availability statement.
   Content: T = **20–130 °C sweep** of relaxation time for devices with electrical diameter
   14.6 ± 0.7, 16.6 ± 0.7, 18.6 ± 0.8, 28.4 ± 1.3 nm; RA = 11.3 ± 1.0 Ω·µm², TMR 70–85 %,
   tau_ave ~100 ms at RT. tau_P and tau_AP defined per state; **the T-sweep is taken at the
   sigmoid shift field where tau_P = tau_AP** — i.e. the dE = 0 axis of C-69's law, where
   tau = exp(E_b/kT)/(2 f0) exactly. Arrhenius fit result stated in text: **tau_0 comes out
   1e-18 to 1e-14 s** (vs 100 ps–10 ns expected), data deviate from the fit line at high T, and
   the resolved (via Fokker–Planck Eq. 3 with alpha = 0.005) ΔE and tau_0 are BOTH strongly
   T-dependent (spin-wave/magnon mechanism proposed). Time-averaged side: <r>–Hz sigmoid
   measured (Fig. 1d); its center slope s_H obeys their Eq. 2 — which is exactly the
   tanh(dE/2kT) slope law with dE = 2 mu0 MS V Hz re-expressed — and s_H(T) decreasing gives
   MS(T) falling ~10 % by 130 °C.
   SEMANTICS: the cleanest open statement that in real s-MTJs **(E_b, f0) are themselves
   T-dependent** — a naive fixed-parameter Arrhenius/C-69 comparison fails by construction, and
   the failure is the registered content, not a defect of the instrument (the record-mode
   extraction must carry per-T parameters). tau measured at the tau_P = tau_AP point means this
   source tests C-69's symmetric factor only; the asymmetry handle lives in the Hz sigmoid.
   EXECUTABILITY: C-69 — digitize Fig. 2b (Arrhenius plots, 4 devices), anchored by ~100 ms RT
   and the stated tau_0 range; run the two-state instrument with (E_b, f0) free per T and
   register the drift. C-70 — digitize Fig. 2c (s_H vs T), test the 1/T slope law of the tanh
   with the MS(T) anchor. Knowns: T (hotplate + thermocouple), diameters, RA, TMR. Unknowns:
   per-point tau values (figures), exact shift fields.

5. **Both-state rates vs current, fitted parameters stated in text** — DATUM (fitted parameters),
   OPEN. A. Mizrahi, T. Hirtzlin, A. Fukushima, H. Kubota, S. Yuasa, J. Grollier, D. Querlioz,
   *Neural-like computing with populations of superparamagnetic basis functions*,
   **Nat. Commun. 9, 1533 (2018)** — OA, loaded at PMC5906599. In-plane elliptic
   60 × 120 nm² CoFeB(1.7 nm)/MgO junctions; escape rates of both states modified by
   spin-transfer torque; fitted per-junction values stated: **ΔE/kT ≈ 8.87–18.68,
   I_c = 85–550 nA, attempt frequency 1 GHz ASSUMED**; natural fluctuation rates kHz–70 kHz.
   Data Availability VERBATIM: "The datasets generated and analyzed during this study are
   available from the corresponding author on reasonable request" — NOT deposited.
   Measurement T not stated as a number (room). SEMANTICS: current as dE-analog is STT
   (nonequilibrium), same caveat as #2/#9; assumed 1 GHz prefactor puts ln f0 inside every ΔE.
   EXECUTABILITY: C-70 ratio-vs-current at figure grade with stated fit anchors; weaker than
   #2/#3 (one handle, no T sweep, request-only data).

6. **Zero-field asymmetric telegraph, single condition, open** — DATUM (single point).
   D. Vodenicarevic et al., *Low-energy truly random number generation with superparamagnetic
   tunnel junctions for unconventional computing*, **Phys. Rev. Applied 8, 054045 (2017)**,
   arXiv:1706.05262v3 — full text extracted (scratchpad/vodenicarevic_1706.05262.txt).
   In-plane 50 × 150 nm² elliptic MRAM-type junctions, room temperature, NO applied field;
   exponential dwell histograms in BOTH states (their Fig. 2a), mean dwell asymmetry attributed
   to the pinned-layer stray field (a built-in dE); stated: F_MTJ = 1.66 kHz,
   **tau_1 + tau_0 ≈ 604 µs**; rate equations r_{0→1} = f0 exp(−ΔE_{0→1}/kT) etc. given.
   No data availability statement (none existed for PRA 2017). EXECUTABILITY: hygiene-grade +
   one fixed-dE point; the built-in stray-field dE is not independently measured — pin as the
   open in-plane twin of entry 1's data class.

7. **Prefactor semantics under bias/field (Meyer–Neldel), open** — CLASS (semantics control
   with numbers). L. Soumah, ..., P. Talatchian (SPINTEC + NIST), *Entropy-assisted nanosecond
   stochastic operation in perpendicular superparamagnetic tunnel junctions*,
   **arXiv:2402.03452v2**, published Phys. Rev. Applied (2024/2025; NIST mirror
   tsapps.nist.gov pub_id 959930); HTML full text loaded. 50 ± 20 nm perpendicular SMTJs at
   T_RT = 294 K; mean dwell times from ~10 ms down to **2.7 ns** under in-plane field
   25–70 mT and bias; Langer-theory analysis yields **Arrhenius prefactors 0.1–10 fs** and a
   confirmed **Meyer–Neldel compensation** tau_0 ∝ exp(−ΔE·const); barriers ~5 kT (ns regime)
   to ~25 kT (µs regime). No data availability statement found.
   SEMANTICS: with entries 2 and 4, this closes the case that across the SMTJ literature the
   measured prefactor spans **~1e-18 to 1e-9 s and co-varies with the barrier** — any C-69
   comparison that fixes f0 a priori manufactures its own refutation; f0 must be an envelope
   parameter, exactly as T-41 treated it.

8. **NIST/SPINTEC single-device dwell statistics, open preprints, no deposits** — CLASS.
   (a) L. A. Pocher et al., *Measurement-driven Langevin modeling of superparamagnetic tunnel
   junctions*, arXiv:2403.11988 — ~20 nm perpendicular device, T = 300 K stated, 505 mV bias,
   124.5 mT field at 18°; dwell-time distributions of BOTH states with ~3e5 transitions; no
   data availability statement (HTML full text checked). (b) T. N. Adeyeye et al., *Sampling
   from exponential distributions in the time domain with superparamagnetic tunnel junctions*,
   arXiv:2412.10317 (Phys. Rev. Applied 23, 044047 (2025)) — P→AP first-passage times only
   (AP dwell relegated to a drift appendix), 918–930 µA current steps, room T assumed; no data
   statement in the arXiv version (checked). SEMANTICS: the working NIST practice confirms the
   arm-wide finding — even measurement-driven modeling papers ship no raw deposits.

## PINNED — SEMANTICS CONTROLS AND DISQUALIFIED CLASSES

9. **Voltage/current asymmetry is NOT Boltzmann dE** — CLASS control. L. Schnitzspan, M. Kläui,
   G. Jakob, *Electrical coupling of superparamagnetic tunnel junctions mediated by
   spin-transfer-torques*, arXiv:2307.15165 (full text extracted). In-plane SMTJ pairs; dwell
   times of both states vs applied voltage; the dwell-time ratio is driven by the STT factor
   (1 ∓ V/V_c) in the exponent — a nonequilibrium torque, not a detailed-balance level
   splitting. Data "available from the corresponding author" (verbatim), not deposited.
   RULE FOR THE LANE: C-70's exp(−dE/kT) is cleanly tested only on FIELD sweeps (Zeeman dE);
   voltage sweeps test an effective-dE extension and must be registered as such.

10. **Joule heating breaks "known temperature" at operating bias** — CLASS control.
    L. Schnitzspan, M. Kläui, G. Jakob, *Nanosecond true random number generation with
    superparamagnetic tunnel junctions — identification of Joule heating and
    spin-transfer-torque effects*, arXiv:2301.05694. 50 nm in-plane SMTJs, RA 15 Ω·µm²,
    dwell times < 10 ns at ~1 MA/cm²; the paper's own finding: switching-rate increase is
    dominated by **elevated temperature at the tunneling site**, i.e. junction T ≠ ambient T
    at high bias. RULE: any kT extracted from a high-bias sweep carries an unmeasured ΔT —
    the "unknown temperature" disqualifier applies WITHIN otherwise-good datasets.

11. **No barrier → no law to test (scope boundary)** — CLASS, disqualified. C. Safranski,
    J. Kaiser, P. Trouilloud, P. Hashemi, G. Hu, J. Z. Sun (IBM), *Demonstration of nanosecond
    operation in stochastic magnetic tunnel junctions*, Nano Lett. 21, 2040 (2021),
    arXiv:2010.14393 — full text extracted. Easy-plane MTJ: fluctuation quantified ONLY by
    magnetoresistance autocorrelation width (~2–5 ns, bias-dependent), ambient ~20 °C; no
    two-state dwell distributions, and by design no activation barrier (diffusive in-plane
    rotation). This is the commissioned disqualifier "two states not thermally activated over
    a barrier" made concrete: C-69/C-70 have no object here. Same exclusion covers the
    easy-plane-dominant PRB 108, 064418 class.

12. **Pulse-actuated (driven) telegraph** — CLASS, disqualified for the free laws. Kent-group
    SMART devices: L. Rehm et al., *Stochastic magnetic actuated random transducer devices
    based on perpendicular magnetic tunnel junctions*, Phys. Rev. Applied 19, 024035 (2023),
    arXiv:2209.01480 (full text extracted — switching PROBABILITY vs ns pulse amplitude at
    T_bath = 295 K and 4 K, Δ = 39 at RT, no free dwell times); A. Sidi El Valli, M. Tsao,
    D. Chen, A. D. Kent, *Tunable random telegraph noise in stable perpendicular MTJs*,
    arXiv:2509.13458 (40 nm devices, room T, zero field; tau_0-state tuned 2.3 µs → 28.8 ns by
    alternating-polarity ns pulses; stated barriers Δ_AP→P ≈ 26, Δ_P→AB ≈ 51; no data
    statement). The two-state statistics here is engineered by the drive protocol; equilibrium
    detailed balance is not present to test. Usable only under an explicit drive-model
    extension — register any use as such.

13. **Paywalled-only, no open mirror (access-blocked)** — CLASS. (a) K. Hayakawa, S. Kanai,
    T. Funatsu, J. Igarashi, B. Jinnai, W. A. Borders, H. Ohno, S. Fukami, *Nanosecond random
    telegraph noise in in-plane magnetic tunnel junctions*, **Phys. Rev. Lett. 126, 117202
    (2021)** — no arXiv posting exists (arXiv API searched); PRL paywalled; APS serves a free
    accepted-manuscript PDF at link.aps.org/accepted/10.1103/PhysRevLett.126.117202 — verified
    to exist (a real-browser load offers the PDF as a file download) but machine fetch is
    Cloudflare-blocked (HTTP 403; curl gets the challenge page), so no numbers were extracted
    here beyond the abstract: in-plane easy-axis SMTJs with relaxation times down to **8 ns**
    at negligible bias, ~5 orders faster than typical perpendicular SMTJs. Its theory twin
    (S. Kanai et al., PRB 103, 094423 (2021)) is likewise arXiv-absent. (b) W. A. Borders,
    A. Z. Pervaiz, S. Fukami, K. Y. Camsari, H. Ohno, S. Datta, *Integer factorization using
    stochastic magnetic tunnel junctions*, **Nature 573, 390 (2019)** — nature.com
    robot-blocked (303 to auth), no arXiv or OA mirror found; data availability unverified.
    Neither is needed for execution: the same group's open members (entries 3, 4) carry the
    two-sided dwell data.

## REPOSITORY SWEEP (negative results, each attempted directly)

- **Zenodo**: API query "superparamagnetic tunnel junction" → 0 open hits; targeted queries
  found exactly ONE open deposit (entry 1) and the restricted benchmark 6767828 (verified
  restricted via API). No other SMTJ telegraph deposits.
- **NIST data portal (data.nist.gov)**: RMM API query → 10 records, all Josephson/metrology;
  zero SMTJ records despite the NIST papers in entry 8.
- **OSF**: API title filter "magnetic tunnel junction" → 0 nodes.
- **IEEE DataPort**: two targeted web searches → no SMTJ/telegraph datasets surfaced (site is
  JS-gated; not directly enumerable — weakest of the four negatives).
- **Materials Cloud**: search/API endpoints returned 404 to direct fetch; no MTJ hits via web
  search. Negative, low confidence, low prior (repository is DFT/simulation-centric).

---

## EXECUTABILITY VERDICT — ARM 1

**The commissioned object exists in the open literature, at two grades:**

- **Raw-data grade (ONE deposit)**: Zenodo 15222667 (entry 1) — open, verified downloadable,
  ~2.8e9 samples of genuine SMTJ telegraph. Single condition, temperature unstated: it grounds
  C-69/C-70's INPUTS (two-valuedness, exponential dwells, the both-values factor
  1/tau* = 1/tau_P + 1/tau_AP at one point) but yields no cosh/tanh sweep by itself.
- **Figure-grade with text anchors (three strong sources)**: Reiss arXiv:1908.02139 (entry 2:
  tau_P AND tau_AP vs H at multiple T and vs V at five fields, Table I anchors — the only open
  source with all three handles); Funatsu Nat Commun 2022 (entry 3: the benchmark's own OA
  figures, field+voltage sweeps, tau_0 = 1 ns assumed); Kaneko APEX 2024 (entry 4: the
  T-swept Arrhenius axis at dE = 0 plus the tanh-slope law vs T). Digitization is permitted for
  all three under the lane's anchoring rule — each states key values in text/tables.

**No open dataset matches the restricted benchmark at raw-sweep grade.** The open record's
shape: raw data without sweeps (entry 1), sweeps without raw data (entries 2–4).

**Semantics that must ride with any comparison** (entries 7, 9, 10): measured prefactors are
entropic/Meyer–Neldel (1e-18–1e-9 s, co-varying with barrier) so f0 is an envelope parameter,
never a constant; voltage/current asymmetry is STT-nonequilibrium, so clean C-70 ratio tests
use FIELD sweeps only; junction temperature at operating bias exceeds ambient (Joule heating),
so kT from high-bias data carries unmeasured ΔT. Cross-references already pinned in
LANE_C86_EXTERNAL/MAGNETIC_ARM_CITATIONS_V001.md: Wernsdorfer N=1 telegraph (entry 1 there),
Krause per-island f0 spread 1e13–1e16 Hz (entry 2 there), and the Funatsu access correction
this file re-verifies.

**Named next steps (route stays open):**
(a) run the record-mode instrument on Zenodo 15222667 chunk 0 — threshold, dwell histograms
    both states, exponentiality, 1/tau* both-values factor; ambient-T assumption declared;
    mind the sampling-grid artifact class (finite 5 µs grid vs sub-grid dwells) before quoting
    any rate;
(b) digitize Reiss Fig. 2 anchored by Table I: C-70 ratio-slope vs H' per temperature and
    C-69 symmetric part vs 1/T — the one open dataset where both laws run on the same device
    with T known per curve;
(c) digitize Funatsu Fig. 4e / Kaneko Fig. 2b as the perpendicular-device replicates (field
    sweep and T sweep respectively), registering the assumed-tau_0 and T-dependent-(E_b, f0)
    semantics as findings, not nuisances.
