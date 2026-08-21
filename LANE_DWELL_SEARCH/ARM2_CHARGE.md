# LANE_DWELL_SEARCH — ARM 2: CHARGE & DEVICE TELEGRAPHY — SOURCE-PINNED (2026-08-21)

Target: OPEN external contact for the two record laws (rows C-69, C-70, post-solidity-review
convention):

    LIFETIME:  tau = 1/(g_up+g_dn) = exp(E_b/kT) / (2 f0 cosh(dE/2kT))
    STEADY:    <R>_ss = tanh(dE/2kT)  <=>  tau_up/tau_dn = exp(-dE/kT)   (detailed balance)

QUALIFIES: thermally activated two-state system, PUBLISHED asymmetric dwell/residence times
(tau_up AND tau_dn or full distributions), temperature known — best if T or a dE-playing bias is
SWEPT. DISQUALIFIES: single-lifetime-only; unknown T; figure-only with no text/table anchors;
paywalled with no mirror; two states not thermally activated over a barrier.

House style per LANE_T41_EXTERNAL/CITATIONS.md and LANE_C86_EXTERNAL/C86_NAND_PINNED_SOURCES.md:
every pin carries SOURCE, NUMBER, UNCERTAINTY, SEMANTICS, ACCESS. Every ACCESS status below was
verified THIS SESSION by loading the resource (method stated per pin); nothing is assumed
(the C-86 "raw data public" lesson governs).

---

## A. SINGLE-ELECTRON BOX — thermally activated 1e transitions (STRONGEST ARM-2 CLASS)

1. **Saira, Yoon, Tanttu, Möttönen, Averin, Pekola — PRL 109, 180601 (2012); arXiv:1206.7049v3.**
   *THE ROW. The closest open equivalent of the restricted Funatsu benchmark found in this arm.*
   - NUMBERS (main + supplement): metallic single-electron box (Cu/AlOx/Al NIS junction), SET
     charge readout, two charge states n = 0,1. **Fig. S1(a): both directional rates Γ(0→1) and
     Γ(1→0), 10–1000 Hz, vs gate offset n_g−1/2 over ±0.10, at bath T = 197, 218, 238 mK.**
     Detailed balance applied explicitly: k_BT·ln(Γ01/Γ10) = 2E_c(n_g−1/2) → **E_c/k_B = 1.94 ±
     0.05 K** (independent JE route gives 1.91 ± 0.03 K — a consistency check inside the paper).
     Rate AT degeneracy vs T (Fig. S1c, ~180–240 mK) fits Γ(0) = sqrt(πk_BTΔ/2)·e^(−Δ/k_BT)/(e²R_T)
     with **Δ = 218 ± 3 µeV, R_T = 100 ± 13 MΩ**; text anchor: **70 Hz at 214 mK**. Q-distribution
     data at bath T = 66–214 mK.
   - UNCERTAINTY: stated errors above; the rate-vs-n_g data are figure-level (log axis) but
     anchored by E_c, Δ, R_T, and the 70 Hz@214 mK text value.
   - SEMANTICS: STEADY law exact — dE = 2E_c(n_g−1/2) swept by gate, both dwell branches
     published at three known temperatures. LIFETIME law in activation form — the degeneracy-point
     rate is activated in the superconducting gap Δ (with a sqrt(T) prefactor, not constant f0:
     pin the prefactor difference). Authors' own caveat: below ~150 mK electron T decouples from
     bath (their bracketed 77/90/130 mK corrections) — the detailed-balance test lives in the
     197–238 mK window where thermal activation is verified.
   - ACCESS: arXiv open; **full PDF incl. supplement loaded and read this session**.
   - **VERDICT: EXECUTABLE.** Both laws touchable in one instrument-grade object. Cost: anchored
     digitization of Fig. S1(a),(c).

2. **Koski, Maisi, Pekola, Averin — PNAS 111, 13786 (2014); arXiv:1402.5907 (Szilard engine).**
   - Same SEB class, ~0.1 K, dwell statistics under feedback. ACCESS: arXiv abs loaded; PNAS
     direct fetch returned 403 (do not cite publisher access as open without the arXiv mirror).
   - **VERDICT: SECONDARY** — redundant with pin 1 for the two laws; feedback protocol adds
     semantics the laws don't need. Use only if pin 1 fails registrar scrutiny.

---

## B. QUANTUM-DOT COUNTING — dE swept by gate at fixed T (STEADY-law testbed only)

These are tunneling rates with NO thermal barrier: they exercise C-70's detailed-balance identity
(occupancy = Fermi function ⇔ <R>_ss = tanh(dE/2kT)) but CANNOT test C-69's activation prefactor.

3. **Gustavsson et al., "Electron counting in quantum dots", Surf. Sci. Rep. 64, 191 (2009);
   arXiv:0905.4675.**
   - NUMBERS: Fig. 10 — measured Γ_in, Γ_out vs dot-level detuning fitted by Γ_in = Γ·f(E),
     Γ_out = Γ·(1−f(E)) with **Γ = 9.2 kHz, T = 230 mK** (T extracted from the Fermi-function
     width; lever arm from Coulomb diamonds). Occupancy identity stated in text:
     f(ΔE/k_BT) = ⟨τ_out⟩/(⟨τ_in⟩+⟨τ_out⟩).
   - UNCERTAINTY: figure-level rates; Γ and T anchored in caption/text.
   - SEMANTICS: exactly the steady law after the occupancy↔(±1) mapping; QPC detector present
     (see pin 6 before trusting exponents).
   - ACCESS: ar5iv full-text HTML loaded this session.
   - **VERDICT: EXECUTABLE (STEADY only).**

4. **Gustavsson et al., PRB 74, 195305 (2006); arXiv:cond-mat/0605365.** Rates through both
   barriers individually, time-resolved QPC counting. ACCESS: arXiv abs loaded; numbers
   figure-level in full text (not pulled this session). **VERDICT: SECONDARY** behind pin 3.

5. **Hofmann et al., PRL 117, 206803 (2016); arXiv:1610.00928.**
   - NUMBERS: waiting times τ_in, τ_out vs detuning at **T ≈ 50 mK** (fixed), feedback-enhanced;
     ratio used as W_in/W_out = p_N/p_(N−1) = m/n to measure LEVEL DEGENERACY (orbital and spin,
     tuned by E and B fields).
   - SEMANTICS — refuter-facing pin: with degeneracy the ratio is g·exp(−dE/kT), NOT exp(−dE/kT).
     The bare two-state tanh law carries no multiplicity factor; any lane meeting degenerate
     carriers must carry this generalization or record a false null. That is this pin's value.
   - UNCERTAINTY: figure-level; no data-availability statement; single T.
   - ACCESS: ar5iv full text loaded this session.
   - **VERDICT: MID** — not the cleanest data row; pinned for the degeneracy trap.

6. **Utsumi, Golubev, Marthaler, Saito, Fujisawa, Schön — PRB 81, 125331 (2010); arXiv:0908.0229.**
   - WARNING PIN (not a data row): bidirectional single-electron counting through a double dot;
     the fluctuation-theorem exponent FAILS with the bath temperature and is repaired only when
     QPC detector back-action (shot noise, level fluctuations) is modeled — i.e. the effective T
     in exp(−dE/kT) is NOT the thermometer T for QPC-monitored dots.
   - ACCESS: arXiv abs loaded this session.
   - **VERDICT: SEMANTICS GUARD for all of section B** (and for pin 1's SET readout, where the
     authors bounded back-action by design).

---

## C. RTN — single-oxide/dielectric-defect telegraphy (LIFETIME-law testbed via T sweep)

7. **Mukherjee et al., "Single-Defect Spectroscopy via RTN in Graphene-Contacted ReS2-hBN
   Heterostructures", arXiv:2511.17125; PRB (2026).** *The strongest open RTN row found.*
   - NUMBERS: clean two-level RTN in one device (RG1), **88–150 K window** (V_ds = 50 mV,
     V_bg = 65 V). Lorentzian corner f_c = 1/τ shifts **mHz → Hz** with T; Arrhenius ln f_c vs 1/T
     → **E_a ≈ 208 meV** (Fig. 2d; text value). Gate sweep at T = 110 K: τ_c and τ_e PDFs and
     means at V_bg = 30–70 V (Fig. 3c–e; τ_c ~0.1–1 s falls with V_bg, τ_e ~1 s near-constant);
     **E_T − E_F from τ_c/τ_e ∝ exp((E_T−E_F)/k_BT): −4 to −9 meV across V_bg = 30–70 V**
     (Fig. 3f); trap located **d ≈ 3 nm inside hBN** via d ≈ −t_sub·(k_BT/q)·δln(τ_c/τ_e)/δV_bg;
     D_it ~ 4e11 cm⁻²eV⁻¹. Temperature-dependent lifetimes Arrhenius in SI (Fig. S11).
   - UNCERTAINTY: "~208 meV", no error bar; τ and E_T−E_F values figure-level (meV scale legible).
   - SEMANTICS: the Arrhenius object is the COMBINED rate 1/τ = 1/τ_c + 1/τ_e — structurally the
     SAME object as C-69's corrected τ = 1/(g_up+g_dn) with its cosh(dE/2kT) factor; the ratio
     τ_c/τ_e is the detailed-balance/dE handle, gate-swept at known T. E_a is capture-cross-section
     activation (NMP two-state double-well, their Fig. 4d) — map it to E_b with eyes open, not by
     name-matching. NO data-availability statement in the preprint (checked to the references).
   - ACCESS: arXiv PDF (31 pp) loaded; pp. 1–14 read this session.
   - **VERDICT: EXECUTABLE** — T-swept combined rate AND bias-swept ratio in one device; both
     law shapes touchable. Cost: digitization of Figs. 2d, 3e–f anchored on the text values.

8. **Grasser, Reisinger, Wagner, Schanovsky, Goes, Kaczer — TDDS, IRPS 2010, pp. 16–25; open PDF
   at iue.tuwien.ac.at (CP2010_Grasser_1.pdf).**
   - NUMBERS: production 2.2 nm plasma-nitrided pMOSFET; 13 individual defects tracked over
     3 months. **Fig. 8 per-defect Arrhenius with E_A printed on the plot** — τ_e (at −0.55 V):
     0.6, 0.79, 0.62, 0.71, 1.14, 1.05, 1.39, 0.89, 1.09, 1.03 eV; τ_c (at −1.7 V stress): 0.56,
     0.46, 0.46, 0.56, 0.51, 0.73, 0.4, 0.99 eV. T window 1000/T = 2.1–2.9 K⁻¹ (~345–476 K);
     τ spans 1e−6…1e4 s; τ_c(V), τ_e(V) at 125/175 °C over nearly full operating range (Figs. 6,
     7, 9).
   - SEMANTICS — TWO DISQUALIFIERS FOR THE STEADY LAW, pinned so nobody trips: (i) τ_c is measured
     at STRESS bias, τ_e at RECOVERY bias — the pair never shares an operating point, so τ_c/τ_e
     is NOT a detailed-balance ratio at any single dE (the paper says this explicitly, contrasting
     with RTN); (ii) the paper itself demonstrates metastable states (tRTN, switching traps,
     disappearing defects — a 4-state Markov model): the defect is NOT a two-state GKSL system.
   - UNCERTAINTY: "approximate activation energies"; no error bars.
   - ACCESS: open PDF loaded and **read in full** this session.
   - **VERDICT: CONDITIONALLY EXECUTABLE** — activation-FORM contact for individual escape
     channels (a spread of real E_b values with real τ ranges); **DISQUALIFIED for the tanh/steady
     law** by (i)+(ii).

9. **Campbell, Qin, Cheung, Yu, Suehle, Oates, Sheng — IIRW 2008, pp. 105–109; open PDF at
   tsapps.nist.gov (pub_id 33210).**
   - NUMBERS: 0.085 × 0.055 µm SiON nMOSFET (1.4 nm oxide), ROOM temperature (stated), source
     −50 mV. τ_capture and τ_emission vs gate overdrive −150…+100 mV at three bandwidths;
     text anchors: **τ_c ≈ 9 ms at V_G−V_TH = −150 mV**; ranges **τ_c = 1e−3…1e−2 s,
     τ_e = 1e−1…1e+1 s**; ratio τ_c/τ_e exponential in overdrive spanning 1e−5…1e0 (Fig. 7).
   - SEMANTICS: both dwell branches at one operating point (true RTN semantics, unlike pin 8),
     ratio exponential in a dE-playing bias at known T — steady-law-shaped. But single
     temperature: no E_b. Class warning carried by the paper itself: the standard
     bulk-oxide-tunneling RTN model FAILS on these devices (τ values orders too slow); the
     microscopic identity of the two states (interface state vs bulk trap) is contested even
     where two-state statistics are clean. Bandwidth systematics bounded in-paper (factor ~2).
   - ACCESS: open PDF loaded and **read in full** this session. Companion (extended) version:
     *Random Telegraph Noise in Highly Scaled nMOSFETs*, IRPS 2009, p. 382 — also open at
     tsapps.nist.gov (pub_id 901584; loaded, p. 1 read: same group, same 1.4 nm devices, same
     elastic-tunneling refutation; no temperature sweep in it either).
   - **VERDICT: MEDIUM** — bias-swept ratio at known T; no T sweep; digitization beyond the two
     text anchors.

10. **Ma, Bi, et al., Nanomaterials 12, 4344 (2022); PMC9741056 (CC-BY).**
    - REGIME-BOUNDARY PIN: 22-nm FDSOI MOSFET, τ_c and τ_e vs V_g at **10–100 K**;
      tunneling-dominated — τ0 not constant, thermal-activation correction coefficient χ falls
      0.34 → 0 as T drops 100 K → 10 K; trap depth 0.13 nm.
    - SEMANTICS: cryogenic RTN LEAVES the thermally-activated class — the boundary the search
      must respect (same reason the 14 K Coulomb-blockade RTN row, arXiv:2206.09086 / IEEE EDL 43,
      5 (2022), is not pinned as a law contact).
    - ACCESS: PMC full text loaded this session. **VERDICT: DISQUALIFIED as law contact; pinned
      to mark the regime boundary.**

11. **Kirton & Uren, Adv. Phys. 38, 367–468 (1989).**
    - The canonical T-swept single-defect τ_c/τ_e lineage (thermally activated capture AND
      emission, E_B ~ 0.6 eV scale) that pins 7–9 descend from.
    - ACCESS: **paywalled (Taylor & Francis); no open mirror found — searched this session
      (re-searched on the second pass, same outcome).**
    - **VERDICT: CLASS PIN ONLY.** The numbers are not open; must not be cited as executable.

12. **Zander, "The Detector Before the System — Spurious Rate Anomalies in Discretely Sampled
    Threshold Detection", Zenodo preprint, DOI 10.5281/zenodo.21935750 (v2, 2026-08-14, CC-BY-4.0).**
    - EXTRACTION-SYSTEMATICS GUARD (not a data row): closed-form framework showing that
      threshold detection with hysteresis on a DISCRETELY SAMPLED noisy two-state signal
      produces apparent rate changes (~7% in its worked case) purely from sampling-grid
      effects — no process change required; supplies a test quantity for flagging the artifact.
      Aimed at superparamagnetic-tunnel-junction readout but generic to Schmitt-trigger dwell
      extraction.
    - SEMANTICS: this is the same systematic family Saira's supplement bounds by Monte Carlo
      (pin 1: missed back-and-forth transitions, threshold-crossing delay) and Campbell bounds
      empirically (pin 9: factor ~2 across 300 Hz–30 kHz). Any digitized-dwell fit in this arm
      should carry its test quantity or an equivalent bound before an exponent is trusted.
    - ACCESS: Zenodo record page loaded this session; single open PDF, CC-BY.
    - **VERDICT: METHODS GUARD for pins 1, 7, 9** — pin before firing any extraction.

---

## D. DATASET REPOSITORIES — searched this session, outcome

- **Zenodo**, query "random telegraph noise" (API, 25 records reviewed): best dataset hit is the
  RTNinja deposit (10.5281/zenodo.16750729) — **SYNTHETIC** Monte-Carlo RTN plus framework
  outputs (companion arXiv:2507.08424). Not measured; disqualified as data (possibly useful as
  dwell-extraction tooling). Remaining hits: presentations, simulations, TRNG circuit data.
  **No open measured dwell-time dataset found under this query.**
- **IEEE DataPort**: no open-access measured RTN dataset identified at search level; full catalog
  needs a (free) login — logged as residual, not as absence.
- **Zenodo second-pass query** "random telegraph" AND (dwell OR capture OR emission) (API): ONE
  hit — the Zander threshold-detection preprint, now pin 12 (methods, not data). A targeted web
  pass for data-availability statements in RTN papers (Zenodo/figshare/DataPort) surfaced no open
  measured dwell-time dataset; the 13,000-trace array-chip RTN study (Microelectron. Eng. 2025,
  ScienceDirect S0167931725001261) is paywalled with no visible open deposit — logged as residual
  to check for a repository link if a mirror appears.

## E. DISQUALIFIED CLASSES (logged so the next lane does not re-walk them)

- NV-center charge-state telegraphy (NV⁻/NV⁰): photo-driven ionization/recombination, not thermal
  activation over a barrier.
- Cryogenic Coulomb-blockade / deep-cryo MOSFET RTN: tunneling regime (pins 10 and the 14 K EDL
  row).
- Synthetic RTN datasets (Zenodo RTNinja class).

## F. NEXT STEPS (no route closes without one)

1. **Fire the Saira row first**: digitize Fig. S1(a) (both directions × three temperatures) and
   S1(c), anchored on E_c/k_B = 1.94 K, Δ = 218 µeV, R_T = 100 MΩ, 70 Hz@214 mK → direct
   tanh/steady fit AND activation-form fit in one object. Highest value per hour in this arm.
2. **ReS2 row second**: digitize Figs. 2d + 3e–f anchored on 208 meV / mHz–Hz / −4…−9 meV; when
   the PRB published version appears, check its SI for data files (preprint has no statement).
3. Residuals: IEEE DataPort logged-in catalog pass; Song et al., Nat. Commun. 8, 2121 (2017)
   (open, MoS2 noise nanospectroscopy) as a second 2D row if wanted; Ensslin-group ETH Research
   Collection for any deposited dot-counting raw data.
4. Cross-arm note: pin 1 (both dwell branches, multi-T, gate-swept dE, open) is this arm's answer
   to the Funatsu benchmark; if Arm 1/3 find nothing stronger, it should lead the combined brief.

---

## G. RE-VERIFICATION PASS (2026-08-21, second independent session)

Every ACCESS claim and every load-bearing number above was independently re-verified in a second
session on 2026-08-21, by loading each resource again (method per pin). Nothing above rests on the
first pass's word alone.

- **Pin 1 (Saira arXiv:1206.7049)**: abs + ar5iv full text + arXiv v3 PDF pp. 3–4 and supplement
  pp. 1–4 read. CONFIRMED first-hand: Fig. S1(a) shows BOTH directional rates (open = 0→1,
  filled = 1→0), 10–1000 Hz, n_g−1/2 ∈ ±0.10, at bath T = 197, 218, 238 mK (caption quoted);
  S1(b) E_c = 1.94 ± 0.05 K; S1(c) degeneracy rate vs T ~180–240 mK with Eq. (1) fit; main text
  Δ = 218 ± 3 µeV, R_T = 100 ± 13 MΩ, 70 Hz at 214 mK, JE route E_c/k_B = 1.91 ± 0.03 K,
  Q-data at 66–214 mK with bracketed electron-T corrections (77/90/130 mK). All match the pin.
- **Pin 2 (Koski arXiv:1402.5907)**: abs loaded — title/authors/SEB confirmed.
- **Pin 3 (Gustavsson review arXiv:0905.4675)**: ar5iv loaded — Fig. 10(c) caption quoted
  verbatim: "Γ=9.2 kHz and T=230 mK"; Eq. (13) occupancy identity f = ⟨τ_out⟩/(⟨τ_in⟩+⟨τ_out⟩)
  confirmed.
- **Pin 4 (cond-mat/0605365)**: abs loaded — rates-through-both-barriers claim confirmed.
- **Pin 5 (Hofmann arXiv:1610.00928)**: abs loaded — detailed-balance degeneracy method confirmed.
- **Pin 6 (Utsumi arXiv:0908.0229)**: abs loaded — QPC back-action "strongly modifies the
  tunneling statistics" quoted from abstract.
- **Pin 7 (Mukherjee arXiv:2511.17125)**: abs + PDF pp. 3–8 read. CONFIRMED: Fig. 2a traces at
  88/96/103/112/123/135/150/230 K (the 88–150 K window stands; the abstract rounds to 90–150 K);
  Eq. (1.5) f_c = f0·exp(−E_a/k_BT) with E_a ~ 208 meV from Fig. 2d (text value); Fig. 3c–e τ_c,
  τ_e PDFs and means at V_bg = 30–70 V, T = 110 K, τ_e ≈ 1 s near-constant, τ_c falling; Fig. 3f
  E_T−E_F axis −4…−9 meV; Eq. (1.6) trap-depth formula and d ~ 3 nm inside hBN; S11 Arrhenius
  lifetimes cited in text. Note for the semantics line: f_c = 1/τ is the Lorentzian corner —
  the COMBINED rate, as pinned.
- **Pin 8 (Grasser TDDS IRPS 2010)**: open PDF re-downloaded from iue.tuwien.ac.at; pp. 1, 4–5
  read. CONFIRMED: title/authors; Fig. 8 left (τ_e, −0.55 V recovery) printed E_A values
  {0.6, 0.79, 0.62, 0.71, 1.03, 1.14, 1.05, 1.39, 1.39, 1.09, 0.89} eV and right (τ_c, −1.7 V
  stress) {0.4, 0.73, 0.99, 0.51, 0.56, 0.46, 0.46, 0.56} eV — same sets as pinned; 1000/T axis
  2.1–2.9 K⁻¹; caption says "approximate activation energies"; stress-vs-recovery bias split,
  tRTN, disappearing/switching traps and metastable-state language all on the pages read —
  both steady-law disqualifiers confirmed in the source's own words.
- **Pin 9 (Campbell IIRW 2008)**: open PDF re-downloaded from tsapps.nist.gov; pp. 1–3 read.
  CONFIRMED: room temperature, source −50 mV, τ_capture ~9 ms printed in Fig. 2(b) at
  V_G−V_TH = −150 mV; text ranges τ_c = 1e−3…1e−2 s, τ_e = 1e−1…1e1 s; three bandwidths
  300 Hz/3 kHz/30 kHz with error "of order of only a factor of two"; elastic-tunneling model
  refuted in text. IRPS 2009 companion access verified (pub_id 901584 loads).
- **Pin 10 (PMC9741056)**: PMC full text loaded — 10–100 K, χ falls 0.34 → 0, trap depth
  0.13 nm, τ_c/τ_e gate dependence: all confirmed.
- **Pin 11 (Kirton & Uren)**: paywall re-confirmed; targeted search again found no open mirror.
- **Section D**: Zenodo RTNinja deposit (16750729) re-loaded — description states synthetic
  Monte-Carlo generation (~1.6 TB, CC-BY-NC); "not measured" confirmed.

Discrepancies found by the pass: NONE that change any verdict. One refinement recorded at pin 7
(88 K vs the abstract's rounded 90 K — resolved in the pin's favor by Fig. 2a). Pin 12 and the
IRPS 2009 companion note were ADDED by this pass.
