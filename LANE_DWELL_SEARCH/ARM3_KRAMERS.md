# ARM 3 — MECHANICAL AND SOFT-MATTER KRAMERS — SOURCE-PINNED (2026-08-21)

Target: OPEN external data for the two record laws of C-69 / C-70 —

    LIFETIME:  tau      = exp(E_b/kT) / (2 f0 cosh(dE/2kT))     [= 1/(g_up + g_dn)]
    STEADY:    <R>_ss   = tanh(dE/2kT)   <=>   tau_up/tau_dn = exp(-dE/kT)

The object sought is a thermally activated TWO-STATE system with published
asymmetric dwell times at known temperature, ideally swept.

House style follows `LANE_T41_EXTERNAL/CITATIONS.md` and
`LANE_C86_EXTERNAL/C86_NAND_PINNED_SOURCES.md`: every pin carries SOURCE, NUMBER,
UNCERTAINTY, SEMANTICS, and **ACCESS STATUS verified by actually loading the
resource**. The C-86 refuter's catch — "raw data public" asserted but false — is
the governing lesson; §E below is a literal transcript of what was loaded, with
HTTP evidence, so that any claim of openness here can be re-checked in one
command. Arithmetic on these pins: `arm3_numbers.py` → `arm3_numbers.txt`.

---

## A. THE FINDING THAT REORDERS THE ARM — the circle in C-70

Before any source: a pair of dwell times **alone cannot test C-70**. Given only
`tau_up` and `tau_dn`, the only available definition of the asymmetry is
`dE/kT := ln(tau_dn/tau_up)`, and then `tanh(dE/2kT) == (tau_dn - tau_up)/(tau_dn +
tau_up)` is an algebraic identity for any two positive numbers. `arm3_numbers.py`
§1 demonstrates this to machine zero on real pinned dwell times.

This is the same failure mode the solidity review already found inside C-70
(detailed balance forces `tanh` before any physics is consulted). It reappears on
the data side: **importing measured dwell times does not by itself convert the
identity into a measurement.** What breaks the circle is a `dE` that is *not*
derived from the dwell times — from a Boltzmann-inverted potential, or from a
control parameter whose mapping to `dE` has been earned.

The arm is therefore ranked below by whether `dE` is independent, not by how
pretty the dwell-time statistics are. This inverts the naive ranking: the largest
open dataset in the arm (§C.1) is the *least* immediately executable, and a
figure caption in a 2017 paper (§B.1) is the most.

---

## B. TIER 1 — INDEPENDENT `dE`, THE CIRCLE BROKEN

### B.1 Rondin, Gieseler, Ricci, Quidant, Dellago, Novotny — levitated nanoparticle, Kramers turnover

**SOURCE.** *Direct Measurement of Kramers' Turnover with a Levitated
Nanoparticle*, Nature Nanotechnology **12**, 1130–1133 (2017);
open preprint **arXiv:1703.07699v2** (11 Aug 2017), loaded in full including
Supplementary Sections A–B.

**NUMBERS.**
- `T = 300 K` — stated verbatim in Methods: *"where k_B is the Boltzmann constant
  and T = 300 K is the temperature of the gas"*. UNCERTAINTY: none quoted; the
  paper's own closing sentence names improved **temperature control** as a
  needed upgrade, so treat as ±few K, uncontrolled.
- Barriers `U_A ≈ 4 k_B T`, `U_C ≈ 5 k_B T` — *"Experimentally, we measure the
  energy barriers U_A ≈ 4k_BT and U_C ≈ 5k_BT"*. UNCERTAINTY: "≈", no error bar;
  the paper attributes its ≤20% theory/experiment mismatch to exactly these.
- Saddle curvature `|ω_S^B|/2π ≈ 51 kHz`; `Γ/P_gas ≈ 51 Hz/Pa`.
- Pressure swept `200 Pa → 2×10^4 Pa` stepwise (Γ from ≈10 kHz to ≈1 MHz).
- Rate `R` from ~1.2×10³ to ~2.4×10³ s⁻¹ across the sweep — **DIGITIZED from
  Fig. 3**, permitted here only because the barriers, T, Γ/P and the turnover
  condition are all stated numerically in text to anchor it.

**SEMANTICS — the reason this is Tier 1.** `R` is obtained by fitting an
exponential to the autocorrelation of the binary population operator `h_A`, and
the paper states plainly: **`R = R_AC + R_CA`**, *"the sum of the rate constants
for the forward and backward processes … the relaxation rate of a non-equilibrium
population in the wells towards equilibrium."* That is **exactly** the corrected
C-69 object `tau = 1/(g_up + g_dn)` — the same relaxation eigenvalue of a
two-state Liouvillian, not the naive `exp(E_b/kT)/f0` the row was demoted for.
And because `U_A` and `U_C` are measured **to a shared saddle B** by Boltzmann
inversion of the position histogram, `dE = U_C − U_A ≈ 1 k_BT` is **independent
of any dwell time**. The circle of §A is broken here and nowhere else so cleanly.

**HONEST NOTE — what the sweep is, and the warning it carries.** The swept
parameter is gas pressure, i.e. the **damping Γ — the prefactor**, not `dE` and
not `T`. With barriers held fixed, rate ratio = prefactor ratio, so their own
Fig. 3 shows `f0` moving by a factor of ~2 across the measured decades and
falling away without bound on the plotted asymptotes. **C-69 assumes a single
`f0`.** Rondin is thus simultaneously the best anchor in the arm and a direct
warning: any C-69 anchor treating `f0` as a material constant carries an unstated
damping-regime assumption. This is the arm's contribution to the register, and it
argues *against* the row's convention rather than for it.

**ACCESS STATUS.** Paper: **OPEN, VERIFIED** — arXiv PDF loaded and read
end-to-end (main text, Methods, Supplementary A/B).
Raw data: **NOT OPEN / NOT VERIFIABLE.** The arXiv version carries no data
availability statement; Acknowledgements are followed directly by Methods. The
ETH Research Collection record (handle 20.500.11850/221251) returned a
**hard block** — an "Access Restricted" page stating access *"is temporarily
restricted from your location / your provider"* due to scraping volume. That is
an **IP block, not evidence about the record**: this lane could not determine
whether that handle is a dataset or a publication record, and must not claim
either. A DataCite lookup on the guessed DOI `10.3929/ethz-b-000221251` returned
empty attributes — also not evidence.

**EXECUTABILITY VERDICT: EXECUTABLE NOW, single-point, from text + one digitized
figure.** `dE ≈ 1 kT` independent; `T = 300 K` stated; measured quantity is
literally C-69's `tau`. Not executable as a swept test of `dE` or `T`.
Blocking item if a *sweep* is wanted: the ETH record's true nature, which needs a
non-blocked route or an author request.

### B.2 Zijlstra, Nettels, Satija, Makarov, Schuler — dielectric bead in a bistable optical trap

**SOURCE.** *Transition Path Dynamics of a Dielectric Particle in a Bistable
Optical Trap*, Phys. Rev. Lett. **125**, 146001 (2020); open preprint
**arXiv:1912.09314**, loaded in full including Methods.

**NUMBERS.**
- `T = 295 K` — stated in Methods, Eq. (1): *"where … T = 295 K is the
  temperature"*. Sample: fused-silica microspheres, 0.54 µm diameter, in water,
  *"at room temperature"*. UNCERTAINTY: no error bar on T.
- Barrier heights **tuned systematically ~2 → 8 k_BT** by laser power,
  *"corresponding to a hundredfold change in transition rate coefficients"*.
- One text-anchored pair (Fig. 1c caption): autocorrelation relaxation time
  **4.4 ± 0.1 s**; *"close to the inverse sum of the transition rate
  coefficients, **3.9 ± 0.3 s**"*. Diffusive component 2.2 ± 0.1 ms.
- 24 trajectories, **7637 transitions**; 33 µs resolution; runs up to ~1 hour.
- `D = 0.73–0.89 µm²/s` measured from MSD, ±0.01 precision.

**SEMANTICS — why this is the best-designed system in the arm.** Three properties
line up with the record laws in a way nothing else pinned here does:
1. **Separate directional dwell times.** Methods, verbatim: *"Most potentials were
   slightly asymmetric, resulting in different barrier heights for transitions
   starting from the left and right wells… For determining rate coefficients as a
   function of barrier height, **dwell times for transitions in both directions
   were analyzed separately**."* This is `tau_up` and `tau_dn`, on the same
   particle, by construction.
2. **`dE` independent.** The full 3D potential is obtained by Boltzmann inversion
   `G_3D(r) = −k_BT ln P_3D(r)` from 10⁶–10⁸ sampled positions — the asymmetry is
   read from the *histogram*, not from the rates.
3. **`f0` independent.** `k_0 = (D/2πk_BT)·√(κ_b κ_w)` is built from the measured
   diffusion coefficient and the measured curvatures at barrier and well. **The
   prefactor is not a fit parameter.** This is precisely the defect the solidity
   review found in C-69 ("the rates ARE the inserted jump strengths"): here the
   rate is *predicted* from independently measured landscape quantities and then
   compared to the observed rate (Fig. 2c).

The text-anchored pair above is, on its own, a **measured** instance of the
identity C-69 rests on — population relaxation eigenvalue vs. inverse sum of the
two rate constants — agreeing at **1.58 σ** (`arm3_numbers.py` §3).

**HONEST NOTE.** The reported barrier asymmetry is incidental ("slightly
asymmetric", uncontrolled), not a swept `dE`; the swept parameter is barrier
*height* via laser power. The per-particle rate values live in **Fig. 2b/2c
scatter only** — no table. Digitization of Fig. 2c is defensible under the brief's
rule (text states T, barrier range, D, and the one anchored rate pair) but yields
`(k_exp, k_Kramers)` pairs *without* their individual `dE`, which is what the test
needs. That is the blocking gap.

**ACCESS STATUS.** Preprint: **OPEN, VERIFIED** (arXiv PDF read in full).
Published PRL page and the group's hosted PDF: **NOT LOADABLE by this lane** —
`journals.aps.org` returned 403; the group-hosted PDF fetched but was not
text-extractable by the tool used, so the *published* version's data-availability
statement is **UNVERIFIED**, and no claim is made about it.
Raw data: **NO DEPOSIT FOUND.** The arXiv version has no data availability
statement (Acknowledgements → Author contributions → Competing interests →
References, checked directly). Searches of Zenodo returned nothing for this work.
**Do not record this as "data available".**

**EXECUTABILITY VERDICT: EXECUTABLE NOW for the C-69 relaxation-eigenvalue
identity at one text-anchored point (1.58 σ).** NOT executable as a swept `dE`
test without the underlying per-particle landscape parameters. Blocking item:
those parameters exist (they are in the authors' analysis) but are not published
in numeric form — an author request is the honest route.

---

## C. TIER 2 — OPEN RAW DATA, BUT THE ASYMMETRY HANDLE NEEDS EARNING

### C.1 Lyons, Devi, Hoffer, Woodside — DNA hairpin 30R50/T4 + bead, PRX 2024

**SOURCE.** Dataset: *Trajectories of non-productive attempts at thermally
activated energy-barrier crossing for beads and DNA hairpin 30R50/T4*,
figshare **DOI 10.6084/m9.figshare.24794955.v1**, published 2024-02-14,
**CC BY 4.0**, single file `Lyons_etal_PRX_2023_Data.zip`, **163,149,357 bytes**.
Paper: Phys. Rev. X **14**, 011017 (2024).

**CONTENTS — verified by reading the ZIP central directory directly, not from the
description.** Three measurement families, all space-delimited text, **1 MHz**
sampling, values in nm:
- `Bead Hopping Data/` — `BeadHopping_SampleTrajectory.txt`, plus
  `Forwards`/`Reverse` `Fluctuations` and `TransitionPaths` files.
  **The forward/reverse split is the directional handle.**
- `Hairpin Constant Force Data/Molecule 1..4/` — `HP_CF_Mol*_SampleTrajectory.txt`
  + `Folding`/`Unfolding` `Fluctuations` and `TransitionPaths`.
- `Hairpin High Stiffness Data/Molecule 1..5/` — `HP_HS_Mol*` equivalents.
- Reference-construct folders `RC_HS_11pN_*`, `RC_HS_16pN_*` (hairpin-free
  control, two forces).

**SEMANTICS.** The `*SampleTrajectory.txt` files are, per the deposit's own
description, *"a segment of the full trajectory"* — continuous hopping traces from
which dwell times **are** extractable, but segments, not the full records
(*"Additional data are available from the authors if needed"*). The
`*Fluctuations.txt` and `*TransitionPaths.txt` files are **within-well excursions
and barrier crossings** — they are not dwells and must not be counted as such.
The folder names `CF` vs `HS` are load-bearing: see the modality note below.

**HONEST NOTE ON THE dE-ANALOGUE — this is where the arm needs care.** For the
**constant-force (CF)** files the tilt is `dG(F) = dG_0 − F·Δx` with `F` externally
fixed by a passive force clamp; the mapping `dE := dG(F)` is **EARNED**, subject
to two named corrections: handle/linker compliance in `Δx`, and the two-state
assumption. For the **high-stiffness (HS)** files the trap separation is fixed, so
the trap spring stores energy when the hairpin opens; the tilt is **not** `F·Δx`
and the landscape is molecule-plus-spring. **For HS the mapping is NOT earned and
must not be asserted.** `arm3_numbers.py` §4 shows the scale that makes this
matter: at 5–16 pN and Δx = 10–20 nm the force term is *tens* of `k_BT`, so `dE`
sweeps through zero over well under 1 pN — a force calibration without an
uncertainty makes the implied `dE` worthless.

**ACCESS STATUS.** Dataset: **OPEN, VERIFIED** — figshare API returned CC BY 4.0
and the file manifest; a ranged GET returned HTTP 206 with a valid ZIP local-file
header (`50 4b 03 04`), and a second ranged GET of the final 120 kB returned the
central directory, from which the file listing above was read. No login, no
embargo. Paper: **NOT LOADABLE by this lane** — `journals.aps.org` 403 and the
APS accepted-manuscript route 403. **PRX is a fully open-access journal, so the
text is open in principle, but this lane did not load it.** Consequently
**temperature is UNPINNED for this source** — not "room temperature", *unpinned*.

**EXECUTABILITY VERDICT: PARTIALLY EXECUTABLE, blocked on two named items.**
The raw two-state hopping traces are genuinely open and directionally resolved.
Blocking: (1) the bath temperature and the constant-force values must be pinned
from the PRX text via a route that actually loads (the 403 is this lane's limit,
not the journal's policy); (2) `Δx` handle-corrected, for the CF→`dE` mapping.
Restrict any first pass to the **CF** folders. This is the arm's best candidate
for a *sweep* if (1) resolves, because CF forces differ across molecules.

### C.2 Bercy & Bockelmann — RNA vs DNA hairpins under tension

**SOURCE.** *Hairpins under tension: RNA versus DNA*, Nucleic Acids Research
**43**(20), 9928 (2015), Laboratoire de Nanobiophysique, ESPCI ParisTech.
**Open Access, CC BY 4.0** (license statement read in the article itself).

**NUMBERS — the only pinned pair of true asymmetric dwell times in text.**
- `T = 29 °C` — and pinned with unusual care by the authors: *"The temperature
  inside the sample during an experiment has been measured. It is 29°C."* This is
  a **measured in-sample** temperature, not an assumed room temperature — rare,
  and exactly what the brief asks for.
- RNA10: **τ_unfolded = 1.77 s, τ_folded = 3.1 s**.
- DNA10: **τ_unfolded = 0.11 s, τ_folded = 0.25 s**.
  Both quoted at extensions *"where both states have a similar probability of
  occupation"*. UNCERTAINTY: none quoted on the dwell times.
- Zero-force rates (Table 3): DNA10 `k_0 = (1.0–3.0)×10⁻⁴ s⁻¹`;
  RNA10 `k_0 = (2.9–5.5)×10⁻⁵ s⁻¹`. Loading rates 1–25 pN/s.

**SEMANTICS.** Constant-**extension** hopping: *"the molecular construct is brought
to a defined extension, which is then maintained constant during 10 s to 2 min"*,
with **stepwise increases** in extension — so the asymmetry *is* swept, but by
extension, in the **HS-like modality of §C.1**, where the force→`dE` map is **NOT
EARNED**. The two quoted dwell-time pairs sit at the near-balance point.

**WHAT THE NUMBERS ACTUALLY DO.** Fed to C-70 they reproduce the identity to
machine zero (`arm3_numbers.py` §1): DNA10 → `dE = 0.821 kT` (21.38 meV at 29 °C),
RNA10 → `dE = 0.560 kT` (14.59 meV). **This is the circle of §A, not a test.**
Their real value is as a **consistency floor**: any pipeline that reads these and
fails to return machine zero is broken. Recorded as such, not as an anchor.

**ACCESS STATUS.** **OPEN, VERIFIED** — article and Figure 7 caption loaded
directly from `academic.oup.com`; CC BY 4.0 statement read in the article. Raw
data: none deposited; dwell-time numbers exist **only in the figure caption**, and
the per-extension-step values are figure-only.

**EXECUTABILITY VERDICT: EXECUTABLE ONLY AS A CONSISTENCY FLOOR.** Two clean
dwell-time pairs at a *measured* temperature — genuinely valuable, but the
asymmetry is not independent and the modality's `dE` map is unearned. Blocking
item for promotion: an independent `dE` at each extension step, which the paper
does not provide.

---

## D. DISQUALIFIED — with the reason, verified rather than assumed

Each of these was pursued far enough to disqualify it **on evidence**, not on the
abstract. Recording them is the point: three of the four are the obvious
first-guess hits for this arm, and all three fail.

| Source | Verified reason for disqualification |
|---|---|
| **Hoffer, Neupane, Woodside — figshare 14897091**, *Transition-path trajectories in DNA hairpin folding*, CC BY 4.0, 780,782,293 B | **NO DWELL TIMES.** The obvious candidate — and it fails. Ranged read of the ZIP central directory shows the archive contains *only* `Raw Folding Transtions <hairpin>_N.dat` / `Raw Unfolding Transtions <hairpin>_N.dat` — isolated barrier-crossing segments — plus `Data Info` files, under four folders (`Hairpin 20R100T4`, `Hairpin 30R50T4`, `Handles only`, `Simulated`). Transition paths are the *excluded middle* of a dwell: no residence times are recoverable. Verified by reading the archive, not by trusting the description. |
| **Ricci, Rica, Spasenović, Gieseler, Rondin, Novotny, Quidant**, *Optically levitated nanoparticle as a model system for stochastic bistable dynamics*, Nat. Commun. **8**, 15141 (2017) | **DRIVEN, NOT THERMALLY ACTIVATED**, on three counts: bistability is created by parametric modulation of a Duffing oscillator; the noise is *injected*, with an **effective temperature `T_N` swept 10–30 K** that is not a bath temperature; and only an overall switching rate `Γ` is reported (`Γ_0 = 1814 ± 96 s⁻¹`), never `τ_A` and `τ_C` separately. Data availability, verbatim: *"available from the corresponding author on reasonable request."* Fails the brief's driven-oscillator exclusion and the unknown-temperature exclusion simultaneously. |
| **Badzey, Zolfagharkhani, Gaidarzhy, Mohanty**, *Temperature dependence of a nanomechanical switch*, Appl. Phys. Lett. **86**, 023106 (2005) | **DRIVEN DUFFING.** The two states are the stable oscillation points of a magnetomotively driven beam at 23.4973 MHz under external modulation — the paper's own Eq. (1) is the driven Duffing equation. Temperature (275–825 mK) is swept, but it degrades *switching fidelity*; the reported "residence time fraction" is a duty-cycle-like occupancy in **Fig. 3c only**, never a dwell-time distribution. Representative of the whole bistable-nanobeam / stochastic-resonance literature, which is driven by construction. |
| **Boneß, Margiani, Belzig, Eichler, Zilberberg — Zenodo 19697506**, *Nonequilibrium Kramers Turnover in a Kerr Parametric Oscillator*, CC BY 4.0, 2026-04-22 | **DRIVEN; NO BATH TEMPERATURE; NO DIRECTIONAL SPLIT.** Genuinely open — `exp_KPO_rates.csv` (26,346 B) was downloaded and its header read. Columns: `drives, detunings, resonances, frequencies, lambdas, muValues, effDissipation, effTemperature, scaledRate` — 147 rows. The temperature is an `effTemperature` in scaled units (≈0.012–0.024), the rate is a single `scaledRate` (no up/down pair), and the system is a parametrically driven Kerr oscillator. The one open "Kramers turnover" dataset in existence, and it is the wrong physics for this row. |
| **LUMICKS Pylake tutorial dataset — Zenodo 7729812**, `hairpin.h5`, 86,987,907 B, CC BY 4.0 | **UNKNOWN TEMPERATURE, NO PROVENANCE.** Verified genuinely open: a ranged GET returned HTTP 206 and the HDF5 magic number `89 48 44 46 0d 0a 1a 0a`. It is real two-state DNA-hairpin hopping force data at constant trap distance, and the tutorial fits a 2-state HMM whose transition matrix (`0.0366` vs `0.0822` off-diagonals) is an asymmetry handle. **But**: no temperature stated anywhere in the dataset or tutorial, no associated publication, no force-calibration provenance, and constant-distance (HS) modality. Assuming room temperature is exactly the move the C-86 refuter punished. **Rejected on unknown temperature.** Could be promoted only if LUMICKS documents the sample temperature. |

---

## E. ACCESS-VERIFICATION LOG — what was actually loaded

Recorded so that no openness claim above rests on an abstract or a description.

| Resource | Method | Result |
|---|---|---|
| arXiv:1703.07699 (Rondin) | PDF fetched, pages 1–8 read | **OPEN** — text, Methods, Suppl. A/B read |
| arXiv:1912.09314 (Zijlstra) | PDF fetched, pages 1–5 and 9–14 read | **OPEN** — Methods incl. `T = 295 K` read |
| ETH Research Collection 20.500.11850/221251 | HTTP GET, browser UA | **BLOCKED** — 1337-byte "Access Restricted" page: scraping-mitigation IP block. Nature of record **undetermined** |
| DataCite `10.3929/ethz-b-000221251` | REST API | Empty attributes — **not evidence either way** |
| figshare 14897091 (Hoffer/Woodside) | API metadata; ranged GET 0–3 (HTTP 206, `PK\x03\x04`); ranged GET final 64 kB → central directory | **OPEN**, contents read → **disqualified on contents** |
| figshare 24794955 (Lyons/Woodside) | API metadata; ranged GET final 120 kB → central directory | **OPEN**, CC BY 4.0, full file listing read |
| Zenodo 7729812 (`hairpin.h5`) | ranged GET 0–63 | **HTTP 206**, HDF5 magic confirmed — open, no login |
| Zenodo 19697506 (`exp_KPO_rates.csv`) | full GET, 26,346 B | **OPEN** — header + 147 rows inspected |
| PMC1458853 (Woodside PNAS 2006) | article fetched | **OPEN** — main text read |
| PMC1458853 SI (Table 1) | 5 supplementary PDFs; EuropePMC `supplementaryFiles`; `fullTextXML` | **NOT OBTAINED** — PMC serves a JS proof-of-work gate (returns 1.8 kB HTML, not PDF); EuropePMC returns empty/404. **SI numbers UNVERIFIED** |
| PNAS `doi/suppl/10.1073/pnas.0511048103` | HTTP HEAD | **403** |
| `academic.oup.com` NAR 43:9928 | article fetched | **OPEN** — Fig. 7 caption + license read |
| `journals.aps.org` PRX 14.011017 / PRL 125.146001 | fetch + accepted-manuscript route | **403 both** — published text **not loaded** |
| PMC4968418 (Woodside, instrumental artifacts) | fetched twice | **reCAPTCHA gate both times** — not read |

---

## F. THE NEAR-MISS WORTH NAMING — Woodside et al., PNAS 2006

Woodside, Behnke-Parks, Larizadeh, Travers, Herschlag & **Block**,
*Nanomechanical measurements of the sequence-dependent folding landscapes of
single nucleic acid hairpins*, PNAS **103**(16), 6190–6195 (2006),
DOI `10.1073/pnas.0511048103`, PMC1458853 (volume/page/authors confirmed against
Crossref) — **main text OPEN and read**. This is the Block-lineage source the
brief names. It is, on paper, the ideal Arm-3 source, and it is recorded here as
a **CLASS pin, not a datum**:

- **20 hairpins**, 6–30 bp stems, 0–100% GC, 3–30 nt loops.
- Folded and unfolded lifetimes measured **separately vs force**, with the
  explicit law `τ_f(F) = τ_f,0·exp(F·Δx_f‡/k_BT)` — force as an exponential tilt.
- `T = 23 ± 0.5 °C` stated, **with an uncertainty**.
- **Passive force clamp** — *"one of the traps was made weaker … pulled to the
  edge of the trap, where the local stiffness became zero"* — i.e. the **CF
  modality where the `dE` map is earned** (§C.1).
- Force-tilt used explicitly: `P_u(F) = {1 + exp[(F_1/2 − F)·Δx]}⁻¹`.
- Rates spanning **six orders of magnitude**; `k_u,0 = 10^(−15±2) s⁻¹`.

Everything the row needs — swept asymmetry, both dwell times, stated `T` with
error, earned modality — and **the numbers live in SI Table 1**, which the main
text refers to and which this lane **could not load** (PMC proof-of-work gate;
PNAS 403). Per the C-86 discipline this is pinned as **CLASS: a table of
per-hairpin `F_1/2`, `Δx‡`, and rate parameters is known to exist and to be
publicly posted; its contents are UNVERIFIED by this lane.** No number from it is
quoted here.

**This is the single highest-value unblocking action in Arm 3.** One successful
load of that SI would likely convert the arm from "one anchored point" to "a
20-system swept table with earned `dE`".

---

## G. VERDICT SUMMARY

| Source | `tau_up`/`tau_dn` split | `dE` independent? | `T` stated | Data open | Verdict |
|---|---|---|---|---|---|
| Rondin 2017 (levitated NP) | no (sum `R` only — **which is C-69's `tau`**) | **yes** (4 vs 5 kT, shared saddle) | 300 K | paper only | **EXECUTABLE, 1 point** |
| Zijlstra 2020 (bead double well) | **yes**, explicit | **yes** (Boltzmann inversion); `f0` also independent | 295 K | paper only | **EXECUTABLE, 1 anchored point (1.58 σ)** |
| Lyons 2024 (hairpin + bead) | **yes** (fold/unfold folders) | via `F·Δx`, **earned for CF only** | **unpinned** | **yes, CC BY** | **PARTIAL — blocked on `T`, `Δx`** |
| Bercy 2015 (RNA/DNA hairpins) | **yes**, numeric, measured `T` | **no** — circular | **29 °C measured** | no | **CONSISTENCY FLOOR ONLY** |
| Woodside 2006 (20 hairpins) | yes (per main text) | via `F·Δx`, earned (passive clamp) | 23 ± 0.5 °C | **SI unverified** | **CLASS PIN — top unblock target** |
| Ricci 2017 · Badzey 2005 · KPO 2026 · Pylake | — | — | effective / none | mixed | **DISQUALIFIED** (driven; unknown `T`) |

**Net for the register.** Arm 3 yields **two independently-`dE` anchored points**
(Rondin, Zijlstra) at stated temperatures, **one genuinely open raw-trajectory
dataset** (Lyons, CC BY, blocked on `T`), and **one consistency floor** (Bercy).
It also yields two findings that cut *against* the rows rather than for them:
the §A circle, which shows imported dwell times do not by themselves rescue C-70
from being an identity; and Rondin's damping sweep, which shows `f0` is not a
constant of the material — an assumption C-69 makes silently.

**Named next steps, in priority order.** (1) Load PNAS 2006 SI Table 1 by a route
that clears the PMC proof-of-work gate. (2) Load the PRX 14.011017 open text to
pin `T` and the CF forces for the already-open Lyons deposit. (3) Determine what
ETH handle 20.500.11850/221251 actually is, from a non-blocked network.
