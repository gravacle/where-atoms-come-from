# LANE_DWELL_SEARCH — RANKED VERDICT (ACCESS-AND-SEMANTICS REFUTER PASS, 2026-08-21)

Refuter mandate: default REFUTED; every claim below survived an independent re-fetch THIS
session, or is marked otherwise. Targets: the top candidates of ARM1_SPINTRONICS.md,
ARM2_CHARGE.md, ARM3_KRAMERS.md, against C-69/C-70's object — two thermal states over one
barrier, BOTH dwell times published, temperature known, dE swept or independently sourced.

Fresh evidence copies: scratchpad `refuter/` (arxiv_1206.7049|1703.07699|1908.02139|
1912.09314|2511.17125 .pdf/.txt, koenders_2412.12783, zenodo_15222667.json, smtj.hdr,
figshare_24794955.json, lyons_head.bin). Retrieved Woodside SI Table 1:
LANE_DWELL_SEARCH/WOODSIDE_TABLE1_RETRIEVED.md.

## A. VERIFICATION RESULTS PER CANDIDATE

### Arm 1.1 — Zenodo 15222667 (Kammerbauer & Schnitzspan SMTJ raw telegraph)
- ACCESS: CONFIRMED OPEN. API record re-loaded (access open, cc-by-4.0); header file
  re-downloaded (206 B) and read byte-identical to the pin (sampling 5.0E-6 s, 1.00E+3
  acquisitions, 2024/02/05, applied voltage 1.000000E-3, "10k Ohm pre-resistance",
  Keithley[sic "Keihley"] 2400, Tektronix DPO 7354, 200k); chunk 0 ranged GET -> HTTP 206,
  bytes 0-15/100,213,314 (size exactly as pinned).
- Companion linkage independently re-verified via arXiv:2412.12783: "available on Zenodo [41]",
  ref [41] = 10.5281/zenodo.15222667.
- DISCREPANCIES (2, citation-level): (i) the record's own title is "Superparamagnetic MTJ time
  series for NP noise" — "Time series sMTJ switching data" is the companion paper's citation
  title, not the deposit's; (ii) the deposit holds 29 parquet chunks (chunk_0..chunk_28), not
  28 — ~2.9e9 samples, ~2.9 GB.
- SEMANTICS: temperature stated NOWHERE (metadata + header re-checked) — the brief's
  unknown-temperature disqualifier applies to LAW contact; single bias, no sweep. The arm's
  own "input-hygiene grade" framing is correct and is the ceiling.

### Arm 1.2 — Reiss/Ludwig/Rott arXiv:1908.02139
- NUMBERS: ALL CONFIRMED verbatim in a fresh arXiv PDF: Fig. 2 caption = mean dwell times
  tau_P/AP vs external field for VARYING TEMPERATURE (tMgO = 1.4 nm); Fig. 3a caption = five
  fields at 65 C; example condition 48 C / -190 Oe / 100 mV, exponential dwell distributions;
  Table I = VE/A*MS {9.28, 8.33, 6.29} A nm^2, dE {2.5, 1.3, 2.4} eV, K* {26, 13.5, 24.9}
  kJ/m^3; Delta- = -(20±2) pJ/V x exp(-(13±6)/nm * tMgO); beta = beta'*tMgO = (30±15) fJ/V m;
  H0 = 157.8 Oe (sigma 6.29) at +80 mV, H0 = 180.6 Oe (sigma 6.95) at -70 mV; ln w ~ 35;
  activation radius (17±4) nm; H_comp field-shift procedure.
- NEW HAZARD FOUND: the source contradicts itself on Fig. 3a's device — caption says
  tMgO = 1.4 nm, text says 1.2 nm ("for tMgO = 1.2nm at 65 C"). Any digitization of Fig. 3a
  must resolve which device it is before a number enters the register.
- ACCESS: arXiv open, re-downloaded. Preprint-only status unchanged (not re-searched).

### Arm 2.1 — Saira et al., PRL 109, 180601 (2012), arXiv:1206.7049
- NUMBERS: ALL CONFIRMED verbatim in a fresh arXiv PDF (v3, main + supplement): Delta =
  218 ± 3 ueV, Ec/kB = 1.94 ± 0.05 K, RT = 100 ± 13 MOhm; 70 Hz at 214 mK; JE route
  Ec/kB = 1.91 ± 0.03 K; Fig. S1 caption verbatim "bath temperatures 197, 218, and 238 mK.
  Open and filled symbols correspond to 0 -> 1 and 1 -> 0 transitions"; S1(b) detailed-balance
  Ec estimates; S1(c) degeneracy rate vs T with Eq. (1) fit; detailed balance
  Gamma(-E) = e^(-E/kBT) Gamma(E) stated; Q-data 66-214 mK with bracketed electron-T
  corrections 66(77)/82(90)/126(130) mK.
- ACCESS: arXiv open, re-downloaded. SEMANTICS: exactly the commissioned object — both
  directional rates, gate-swept dE = 2Ec(ng-1/2), three known bath temperatures, activation
  in Delta at degeneracy; dE independent TWICE (detailed balance Ec AND the Jarzynski route).
  Prefactor is sqrt(T)-bearing, not constant f0 — carries as a registered semantic, not a
  defect. NO discrepancy found.

### Arm 2.2 — Mukherjee et al., arXiv:2511.17125 (ReS2-hBN single defect)
- NUMBERS: CONFIRMED in a fresh arXiv PDF: Ea ~ 208 meV (twice in text); two-level window
  "~ 88-150 K"; gate sweep at Vds = 50 mV, T = 110 K, Vbg 30-70 V; ratio law
  tau_c/tau_e ∝ exp((ET-EF)/kBT); trap-depth Eq. (1.6) and d ~ 3 nm inside hBN; Lorentzian
  corner = combined rate. NO data-availability statement (confirmed absent).
- REFUTER QUALIFICATION: the pinned "ET-EF = -4...-9 meV" appears NOWHERE in text — Fig. 3f is
  its only carrier. The arm's own uncertainty line said figure-level; the executable-brief
  wording should not call it an anchor. Anchors that exist in text: 208 meV, 110 K, 50 mV,
  30-70 V, d ~ 3 nm.

### Arm 3.1 — Rondin et al., arXiv:1703.07699 (levitated nanoparticle)
- NUMBERS: ALL CONFIRMED verbatim: "T = 300 K is the temperature of the gas"; "we measure the
  energy barriers UA ≈ 4kBT and UC ≈ 5kBT"; |omega_S^B|/2pi ≈ 51 kHz; Gamma/Pgas ≈ 51 Hz/Pa;
  pressure swept 200 Pa -> 2e4 Pa; "R = RAC + RCA. Note that R is the relaxation rate of a
  non-equilibrium population in the wells towards equilibrium" — literally C-69's corrected
  object. No data-availability statement (confirmed absent). ETH handle status remains
  undetermined (not re-attempted).

### Arm 3.2 — Zijlstra et al., arXiv:1912.09314 (bistable optical trap)
- NUMBERS: ALL CONFIRMED verbatim: "T = 295 K is the" (Methods); "dwell times for transitions
  in both directions were analyzed separately"; relaxation 4.4 ± 0.1 s vs inverse sum of rate
  coefficients 3.9 ± 0.3 s (Fig. 1 caption region). No deposit statement found — confirmed.

### Arm 3 Tier-2 — Lyons/Woodside PRX 14, 011017 + figshare 24794955
- DEPOSIT ACCESS: RE-CONFIRMED OPEN — figshare API (CC BY 4.0, single ZIP, 163,149,357 B),
  ranged GET -> HTTP 206 + ZIP magic PK\x03\x04.
- PAPER: the APS 403 was BEATEN this session — full PRX text loaded in the browser pane, and
  the SI PDF (Lyons_etal_SI.pdf, 2,149,370 B) fetched in-page and text-extracted.
- REFUTATION OF THE ARM'S UNBLOCK PLAN: the temperature is NOT IN the paper. Zero occurrences
  of "temperature"/degree-C/Kelvin values in the full main text (53,953 chars rendered) or in
  the readable SI caption text. Arm 3's named next step "load the PRX text to pin T" CANNOT
  succeed — T must be imported from the group's protocol refs (their Ref. [10]), a mapping cost.
- NEW SEMANTICS FOUND (Methods, verbatim): constant-force measurements were done "adjusting the
  power of the zero-stiffness trap ... until the hairpin spent approximately equal amounts of
  time in the folded and unfolded states" — the CF data are TUNED TO dE ≈ 0; there is no CF
  bias sweep. Forces stated in the SI: 13.9 pN (hairpin-free REFERENCE construct), 11.1 and
  16.1 pN (non-constant-force handle/bead calibrations). Per-molecule CF forces: not stated.
- NET: demoted from "best candidate for a sweep" to raw-hygiene + one dE≈0 point with imported T.

### Arm 3 §F — Woodside et al., PNAS 103, 6190 (2006): UNBLOCKED, CLASS -> DATUM
- The PMC proof-of-work gate was cleared by real-browser navigation (the page's own script ran
  as designed; no challenge computed out-of-band, no CAPTCHA solved). SI Table 1 turned out to
  be an HTML anchor (#T1) INSIDE pnas_0511048103_index.html — not a PDF — which is why every
  PDF route missed it.
- VERIFIED: T = 23 ± 0.5 C verbatim in Supporting Text, INCLUDING the 2.0 ± 0.5 C laser-heating
  correction; full 20-hairpin Table 1 retrieved (Dx, F_1/2, DG, ln k_u0, ln t_1/2, Dx‡_f,
  Dx‡_u, experiment AND model, with uncertainties); SI Fig. 4 text-grade measured pair:
  unfolded-state tau = 41 ms @ 14.3 pN and 11 ms @ 13.5 pN (hairpin 20R55/T4).
- Full table + transcription cautions: WOODSIDE_TABLE1_RETRIEVED.md.

### Benchmark control — Funatsu Zenodo 6767828
- RE-CONFIRMED RESTRICTED (third independent check): API access_right "restricted", no files.

## B. RANKED TABLE (survivors, by first-target priority; READINESS per the commissioned scale)

| # | Source (arm) | Readiness | Both dwell branches | dE independent? | T known | Access (verified how) | What still costs |
|---|---|---|---|---|---|---|---|
| 1 | **Saira PRL 109, 180601 / arXiv:1206.7049 (A2)** | **NEEDS-DIGITIZATION** (anchored: Ec, Delta, RT, 70 Hz@214 mK, JE cross-check) | YES — Gamma(0->1), Gamma(1->0), Fig. S1a | YES, twice (detailed-balance Ec + JE Ec) | 197/218/238 mK bath | arXiv PDF re-downloaded | log-axis digitization of S1a, S1c |
| 2 | Woodside PNAS 2006 (A3) | **READY** on the parameter table (retrieved, numeric, uncertainties); NEEDS-DIGITIZATION for law-vs-DATA curves (Fig. 3D-F) | YES — t_1/2 + Dx‡_f/Dx‡_u parameterize both branches; SI Fig. 4 two measured points | YES — F_1/2, Dx, DG per hairpin, passive clamp (earned CF) | 23 ± 0.5 C incl. laser heating | PMC via real-browser PoW; SI table extracted in-page | fits-vs-fits circularity if Table 1 alone is used; Fig. 3 digitization for measured curves |
| 3 | Reiss arXiv:1908.02139 (A1) | NEEDS-DIGITIZATION | YES — tau_P and tau_AP vs H at multiple T | Field version yes (Zeeman, after H_comp shift); voltage version is STT (excluded for clean C-70) | Per curve (48-65 C range shown) | arXiv PDF re-downloaded | Fig. 2 digitization; resolve the Fig. 3a 1.2-vs-1.4 nm source conflict; preprint-only |
| 4 | Mukherjee arXiv:2511.17125 (A2) | NEEDS-DIGITIZATION + NEEDS-MAPPING | YES — tau_c, tau_e (gate-swept at 110 K); combined rate T-swept | Gate plays dE at fixed T; Ea is NMP capture activation, mapping to E_b must be earned | 88-150 K window | arXiv PDF re-downloaded | Fig. 2d/3e-f digitization with only 208 meV/110 K as text anchors; ET-EF values are figure-only |
| 5 | Zijlstra arXiv:1912.09314 (A3) | READY — one text-anchored point (4.4±0.1 vs 3.9±0.3 s, 1.58 sigma; T = 295 K; f0 independent) | YES (stated; per-particle values unpublished) | YES (Boltzmann inversion) | 295 K | arXiv PDF re-downloaded | not expandable without author data — single point is the ceiling |
| 6 | Rondin arXiv:1703.07699 (A3) | NEEDS-DIGITIZATION — single point | NO — sum R only (which IS C-69's tau) | YES (shared saddle, 4 vs 5 kT) | 300 K | arXiv PDF re-downloaded | no C-70 contact; f0-not-constant finding carries |
| 7 | Zenodo 15222667 (A1) | NOT-USABLE for law contact (T unknown = brief's disqualifier; no sweep); READY as input-hygiene instrument feed under a DECLARED ambient-T assumption | Extractable from raw | n/a (single condition — ratio is definitional) | NOT STATED anywhere | API + header + ranged GET re-verified | record-mode run is hygiene-grade only |
| 8 | Lyons figshare 24794955 (A3) | NEEDS-MAPPING (T import from protocol refs; Dx handle correction; CF is dE≈0 by construction) | YES (fold/unfold folders; segments only) | CF earned in principle, but tuned to balance | ABSENT from paper+SI (refuter finding) | figshare API + ranged GET; PRX text loaded | demoted from sweep candidate |
| 9 | Bercy NAR 2015 (A3) | NOT-USABLE as a law test (the §A circle — identity to machine zero); consistency floor only | YES, numeric in caption | NO | 29 C measured | not re-verified this session (arm's verification stands unchallenged) | none — floor role only |

## C. THE NAMED FIRST TARGET

**Fire the Saira row first** (ranked #1). Reasons, all verified this session: it is the only
object in the three arms where BOTH record laws are touchable in ONE instrument-grade system
with dE independently sourced TWICE (detailed-balance Ec = 1.94 ± 0.05 K vs Jarzynski
Ec = 1.91 ± 0.03 K — the circle that voids naive C-70 tests is broken inside the source
itself), both directional rates published at three known temperatures, and the activation form
carried at degeneracy with stated Delta and a declared non-constant (sqrt-T) prefactor. Its
only cost is anchored digitization, bounded by four text values. Woodside 2006 (#2) fires
second as the 20-system swept replicate — its independent-dE ingredients are now numeric and
in hand (WOODSIDE_TABLE1_RETRIEVED.md), but its dwell content is fit-parameter-grade until
Fig. 3D-F is digitized, and the exponential force law is assumed by those fits.

**Nothing located requires falling back to the Funatsu access request as the first move** —
but file it in parallel anyway: 6767828 stays the only raw-sweep-grade SMTJ benchmark
(re-confirmed restricted), and no open equivalent at that grade exists in any arm.

## D. NEXT STEPS (no route closes without one)

1. Grounding lane on Saira: digitize Fig. S1a (2 directions x 3 temperatures) and S1c against
   the four text anchors; carry the Zander/Arm-2-pin-12 sampling-artifact test quantity.
2. Second lane on Woodside 2006: re-read #T1 in a browser (one load) to clear the transcription
   cautions, then digitize main Fig. 3D-F for law-vs-data contact; Table 1 supplies dE and T.
3. Registrar note for Arm 1: correct the deposit title and chunk count in ARM1_SPINTRONICS.md
   entry 1 (two citation-level errata, access verdict unchanged).
4. Registrar note for Arm 3: mark §C.1's unblock item (1) as REFUTED-as-formulated (T is not in
   the PRX text); Lyons T must be imported from the protocol lineage if the deposit is ever used.
