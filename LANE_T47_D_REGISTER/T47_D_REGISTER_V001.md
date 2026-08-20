# T47_D_REGISTER_V001 — Registrar's draft: the count law k(t_m) registers

Lane: LANE_T47_D_REGISTER (the audit's creator step, T-47; commissioned by C-84).
Inputs: three sealed lanes, each independently verified adversarially, all surviving —
- LANE_T47_A_WIDTH (t47_a_width.py/.txt; VERIFY/verify_t47a.py/.txt, rerun_lane.txt — bit-identical rerun, exit 0)
- LANE_T47_B_STAIRCASE (t47b_staircase.py/.txt, 15 PASS 0 FAIL; VERIFY/verify_t47b.py/.txt — 47/47 independent checks, exit 0)
- LANE_T47_C_OWNERSHIP (T47_C_OWNERSHIP_V001.md, control_departure.py; VERIFY/verify_t47c.py/.txt — 27 checks, exit 0, plus control rerun)

---

## 1. VERDICT: REGISTER, PARTIAL OWNERSHIP

**The width is DERIVED.** No smuggled tolerance anywhere: three independent tokenizer/source
censuses of every numeric literal on the conclusion paths (lanes A, B, C VERIFY) found only carrier
parameters, SI constants, declared instrument floors (the sealed 1e-9 grounded.py floor, unmodified),
machine-precision bisection, and bracket constants proven non-load-bearing (crossing spread 0.0
under perturbation). The 1e-2/1e-3 strings survive only in prose describing the C-76 failure.
**The C-76 failure mode has no successor in any of the three lanes.**

**k(t_m) holds on every carrier tested:** the two-well grid (t_m 1 s..100 y × f0 1e9..1e13 × T
200..400 K; worst exact-formula residual 6.0e-15 relative, independently re-bisected on a full-
eigendecomposition instrument and on a decay-fit instrument that uses no eigen-mode selection at
all); uniform N=12 and lognormal-like N=10 ensembles (instrument count vs formula count discrepancy
0 at all 1201 grid points, each ensemble); the symmetric D-15 control (staircase flat at k=8, single
threshold exact to 1e-12); the T-31 generic-asymmetry [[4,2,2]] carrier rebuilt from scratch
(staircase 2→1→0 while exact-multiplicity v_2 reads 0 at every eps>0); and the coherence-face D-15
control (hbar/t_m governs, crossing exact to 2.1e-16, full-|λ| rule confirmed load-bearing).
60-digit mpmath root agrees with the closed form to 8e-62.

Ownership is PARTIAL (fragments conceded below), so the row registers the wholly-owned LAW with
the concessions written into it.

---

## 2. THE REGISTERED LAW (draft C-row text; number assigned by the registrar at append — C-85 expected)

> ### **THE SURVIVING-RECORD COUNT LAW — the width is the carrier's own, derived, two-faced by
> ### necessity; k(t_m) is the URM's wholly-owned falsifiable count law.**
>
> ONE criterion — clause (ii') applied to the record's OWN Liouvillian mode, |λ_record| ≤ 1/t_m —
> yields both widths with zero adjustable content:
>
> **Population face** (diagonal record, [H,R] = 0): with barrier B above the lower well, splitting
> dE, shared Arrhenius prefactor f0, the record's mode decays at g_u + g_l =
> f0 e^(−B/kT)(e^(dE/kT) + 1), and the EXACT registered width is
>
>     δ_pop(t_m) = kT ln( e^(B/kT)/(f0 t_m) − 1 )   [compute as kT·ln(expm1(B/kT − ln f0 t_m))]
>              = [B − kT ln(f0 t_m)] − kT ln(1 + e^(−dE*/kT)),
>
> i.e. the candidate B − kT ln(f0 t_m) minus a DERIVED two-sided-escape correction bounded in
> (0, kT ln 2], attained exactly at the symmetric corner. The exact form, not the naive form, is
> the law: near threshold the naive width errs by up to 0.27 kT, all of it the derived correction.
>
> **Coherence face** (off-diagonal record, [H,R] = −dE·R): the SAME criterion on the same carrier
> gives δ_coh(t_m) = ħ/t_m. With dissipation the coherence eigenvalue modulus is ONE expression,
> |λ| = sqrt( ((g_u+g_l)/2)² + (dE/ħ)² ), whose two corners are the two widths — 22 orders of
> magnitude apart on the same surface at the same t_m. This is C-76's error quantified: the
> coherence width ħ/t_m was the wrong width for population-type records by that factor.
>
> **THE COUNT LAW:**
>
>     k(t_m) = #{ i : dE_i ≤ kT ln( e^(B_i/kT)/(f0 t_m) − 1 ) }  ≡  #{ i : g_u(i) + g_l(i) ≤ 1/t_m }
>
> — a DECREASING STAIRCASE in ln t_m, one step per record, each record dropping out at the
> parameter-free time
>
>     t*_i = f0^(−1) e^((B_i − dE_i)/kT) / (1 + e^(−dE_i/kT)).
>
> A record counts only while BOTH values are durable — the worst case over values, the SHALLOWER
> well's escape — not while the as-written value survives.
>
> **Corner limits, all verified:** (i) symmetric carriers (dE = 0) — staircase flat, single
> threshold at t_m = e^(B/kT)/(2 f0) (the commission's f0^(−1) e^(B/kT) carries a derived factor-2
> correction: both wells escape); below it C-14's count is recovered exactly, min_E v_2(m_E) = N =
> the flat-staircase count, computed for N = 1..8. (ii) t_m → ∞ — δ_pop → −∞ (no thermally
> activated two-valued record survives forever at T > 0; DEF-A's E_b → ∞ demand recovered) and
> δ_coh → 0 (the exact commutant count, C-75). (iii) The no-crossing condition of the exact form IS
> the symmetric bound — the same inequality, proven.
>
> **Scope, declared:** thermally activated Markovian (GKSL/Davies) two-state carriers with shared
> Kramers prefactor f0; the model declines non-thermal surfaces (thermal=False). On multi-pathway
> carriers the one-f0 closed form leaves ln-residuals (0.62 ln-units on the T-31 carrier at
> eps = 0.16, exceeding the local step gap) — there the staircase and its exponent are carried by
> the record-mode instrument, not the one-prefactor closed form. The joint word of k records decays
> at the SUM of their rates: the family's joint-correlation horizon is shorter than each member's —
> a named open edge for any family-level claim.

**PROVED-bar statement (any physicist, their own data):** take your carrier's (B_i, dE_i, f0, T)
— measured, not fitted to this law — and a retention spec t_m; the law predicts the integer count
of functioning two-valued records and the ln-t positions of every step, with no adjustable width
anywhere. RECORDS VERIFIED candidates, named: (a) exchange-biased / asymmetric magnetic media
grains (MFM grain-by-grain census vs ln t; Weller–Moser/Sharrock-class media with set exchange
bias); (b) MLC NAND flash Vt levels under JEDEC-class retention bake (Cai-et-al. HPCA 2015-class
published Vt-vs-bake data; devices with established two-state Arrhenius detrapping only).

---

## 3. C-76 RECONCILIATION (the corner reconnecting to v_2)

**C-76's objection is answered in its own terms: the width that recovers a count is the carrier's
own, not a choice** — at the symmetric corner the derived width makes the staircase flat and
min_E v_2(m_E) = N equals the count exactly (computed, N = 1..8, and in the separated regime
[4,56,4] → v_2 = 2 = k); under generic asymmetry v_2 correctly reads 0 forever while k(t_m) is the
dated census that replaces T-31's binary kill (2→1→0 on the rebuilt T-31 carrier); clustered-v_2 is
hereby DEMOTED to a corner proxy — exact at dE = 0 and in the separated regime, broken by chain
merging (exhibited: dE = (0.90, 1.05)δ gives proxy 2 vs true k = 1) — and **the law is registered
on record modes, not on width-clustered spectra.** T-31's CLUSTERED runs are re-read as the
coherence-face count, whose ħ/t_m width is now derived, not borrowed.

---

## 4. OWNERSHIP VERDICT: PARTIAL — the LAW is owned; fragments conceded by name

Fairness rule (C-84): a rival owns a statement only if its account already made it.

**Conceded, by name:**
- Neel 1949; Street–Woolley 1949; Sharrock 1994 — the activation window kT ln(f0 t) as a
  remanence-DECAY device (viscosity S, fluctuation field, field-reduced barriers E_B(H)).
- Charap–Lu–He, IEEE Trans. Magn. 33, 978 (1997); Weller–Moser, IEEE Trans. Magn. 35, 4423 (1999)
  — the dE = 0 corner criterion B ≥ kT ln(f0 t_m) and its density-count consequence (the
  ownership-detection control fires OWNED at the corner, as it must; T-41 source-pinned).
- Neel/Brown field-tilted two-well relaxation and the Korman–Mayergoyz Preisach–Arrhenius /
  Preisach–Neel aftereffect lineage — (B − dE)-type effective barriers and two-sided rates AT THE
  RATE LEVEL; their aftereffect observable is a moment-WEIGHTED survivor integral with a fitted
  density, never an integer census.
- Standard two-state relaxation kinetics — the (1 + e^(−dE/kT)) both-rates factor itself
  (1/τ = forward + backward is chemical-relaxation textbook material).
- Exchange-bias literature (York group) — the QUALITATIVE recognition that one-directionally
  stable magnetization cannot store a bit; "stated as a LAW" is load-bearing in the owned fragment.
- MLC flash engineering — level-decay ordering (Cai et al., HPCA 2015), worst-pattern retention
  practice (JEDEC JESD22-A117/JESD47), and physics-based double-well retention models (ABDWT
  class) which remain parameter-fitted.
- Alicki–Fannes–Horodecki / Bravyi–Terhal / BPT — coherence-face storage-TIME bounds and kd² ≤ cn;
  the composition through a derived width was never made; adjacency is not ownership.

**Wholly program-owned (the registered fragment):** (a) the width DERIVED from the record
definition — clause (ii') on the record's own mode — rather than fitted or spec'd; (b) the
diagonal/coherence unification: BOTH widths, activation and ħ/t_m, from the ONE criterion, two
corners of one modulus — no rival states it; (c) the two-valued-record COUNT OBJECT as a LAW: the
discrete, margin-free integer census k(t_m) with parameter-free drop times t*_i, dying at the
shallower value's escape; (d) the computed departure term (below). C-84's shared gap — no rival
account owns a surviving-record COUNT law — SURVIVES all stretch directions tested (Sharrock/
Street–Woolley, Preisach thermal aftereffect, field-tilted Neel–Brown, MLC incl. ABDWT, Alicki/QEC).
Ownership calls on flash/Preisach full text remain search-grade (lane C caveat 1) — hardening pass
named below.

**The departure term, computed:** Sharrock/viscosity counting integrates surviving POPULATION;
k(t_m) counts FUNCTIONING TWO-VALUED RECORDS. A grain written into the deep well holds
m ≈ tanh(dE/2kT) indefinitely while the record is dead; departure = Σ_over-dead tanh(dE_i/2kT) —
zero on symmetric media, → N as dE/kT grows (4.84/6.70/8.09/9.05 of 10 across asymmetry scales
0.5/1/2/4; at the declared carrier B = 60 kT, dE = 10 kT: remanence 0.999909 forever while k = 0
past t* = 5.18e12 s). Equivalently: the record dies a factor e^(dE/kT) EARLIER than the
favored-branch (Sharrock-inferred) barrier implies — Δln t* = −dE/kT.

---

## 5. FALSIFIERS (existing data classes; no lab needed)

**F1 — MLC flash (JEDEC-class bake data):** from SHORT bakes extract per-level (B_i, dE_i-analogue,
f0); the law predicts the level-retirement staircase positions ln t*_i = (B_i − dE_i)/kT − ln f0
(the (1+e^(−dE/kT)) factor is untestable at dE ≫ kT — the dE-shift and the step SPACINGS
(dE_i − dE_j)/kT are what bite) and the order (highest level first). A measured dropout order or
spacing inconsistent with the recorded level asymmetries, beyond extraction precision, falsifies
k(t_m). Scope: devices with established two-state Arrhenius detrapping (the model declines
distributed-trap/stretched-exponential surfaces). First pass on published Cai-et-al.-class
Vt-vs-bake data — commissioned, named next step.

**F2 — exchange-biased / asymmetric media:** write both polarities, read both. The law predicts
record death at 1/(g_u + g_l) — a factor e^(dE/kT) before the along-bias remanence lifetime — while
remanence ≈ tanh(dE/2kT) persists. Falsified if against-bias retention matches along-bias retention
at measured dE ≫ kT. **Against a Preisach–Neel adversary, only the grain-by-grain integer census
vs ln t (MFM) plus the derived-vs-fitted width discriminates — the falsifier must be run as a
survivor COUNT, not as drop-time asymmetry alone.** The both-values factor (up to 2× in t*) is
measurable only where dE ≲ a few kT; the count-vs-magnitude and asymmetry-term departures are the
robust falsifiers.

**F3 — QEC:** honestly N/A today; no count-grade data class exists on that surface. Stated, not
claimed.

---

## 6. REGISTER NOTES — errors and instrument scope (per the errors-belong-in-the-register rule)

1. **Two derived corrections to the commission's own text stand:** the symmetric corner bound is
   e^(B/kT)/(2 f0), not f0^(−1) e^(B/kT); and t*_i carries the (1 + e^(−dE_i/kT)) factor. Neither
   alters the headline width.
2. **H-units convention (load-bearing for coherence):** `grounded.liouvillian` takes H in generator
   (angular-frequency) units; `project_model` passes joules — harmless for every sealed
   population-mode number ([H, σ_z] = 0) but coherence-face callers MUST pass H/ħ. Confirmed
   against the model's own code by two independent VERIFY passes.
3. **`grounded.clause_ii`'s Rayleigh quotient misreads multi-shell carriers** (0.497 vs true 3.6e3
   on the 16-dim carrier — uniform trace weighting swamped by excited-shell flips); on such
   carriers the record's mode must be taken in the basin-lumped slow sector. Same failure class as
   the spectrum() 2× lesson, one level up. (On the two-well carrier the quotient is provably exact:
   I-component and rotation are HS-orthogonal.)
4. **Errors caught and logged in the lanes, none surviving into conclusions:** lane B — a
   scale-relative mode-exclusion floor let the coherence frequency swallow the slow population mode
   (C-75's error class), and an instantaneous-outflow basin rate overcounted by the recrossing
   factor ~2; both replaced. Lane C — first coherence control used σ_x (not an eigenmode) and read
   durable for both cases; fixed to the record's own mode σ_+, failure reproduced by VERIFY. Lane
   B-VERIFY's own first decay-fit instrument had a bug (its own, found and logged; the lane's
   numbers never implicated).
5. **One check-slack flag (not a D-24 kill):** lane B's F-check "residual ≪ step separation" passes
   only via its +2.0 ln-unit slack; at eps = 0.16 the residual (0.62) exceeds the step gap (0.138)
   — the staircase order there is carried by the instrument, exactly as scope-clause registered.
6. **σ_x/probe hygiene:** C5's 0.5/20 probe points are probes, not widths; the coherence boundary
   independently computes to 1.0000000000000000 ħ/t_m (probed 0.9/1.1).

---

## 7. NAMED NEXT STEPS (no route closes without one)

- **N1 (commissioned by lane C):** F1 first pass on published Cai-et-al.-class Vt-vs-bake data —
  extract (B_i, dE_i, f0), predict the retirement staircase, compare.
- **N2:** unequal attempt frequencies f0_u ≠ f0_l shift δ_pop by kT ln(f0_u/f0_l) inside the same
  derivation — run it if data demands.
- **N3:** the joint-word decay edge (summed rates) for any family-level claim.
- **N4:** full-text ownership hardening pass on interacting/biased Preisach thermal-relaxation
  papers and flash channel-capacity papers (lane C caveat 1; same provisional grade as T-39's S-1
  concession).

## 8. DISPOSITION

Register the C-row above (head-number at append), mark T-47 DONE, carry notes 1–3 as register
notes on the model files, and log N1–N4. The three lane directories plus their VERIFY
subdirectories are the sealed evidence; this lane contains only this draft.
