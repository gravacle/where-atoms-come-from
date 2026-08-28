# HUST-2018 dual-method ordinary-gravity forward theorem

**Theorem ID:** `GRA-HUST18-DMGF-V001`

**Date:** 2026-08-27

**Claim class:** exact reproduction of the published processed-coefficient
forward equations; deterministic extraction from official figure-level
raw-like period and acceleration series; exact ownership and identifiability
ceiling against `GRA-GC-CFAGC-V001`.

**Status:**
`OFFICIAL_FIGURE_LEVEL_TOS_RESPONSE_RECOVERED__THREE_PROCESSED_AAF_SOURCE_RESPONSE_FORWARDS_REPRODUCED__FULL_GC16_NOT_EXECUTABLE__NO_NEW_G`

## 1. The closure gained

The official HUST release supplies a real, dual-channel finite-apparatus
instantiation that is stronger than the previous PDF-only readiness result.
The time-of-swing (ToS) channel measures a gravitational **operator/stiffness**
contrast, while angular-acceleration feedback (AAF) measures a gravitational
**source/forcing** response. In the notation of the calibrated finite-apparatus
model, they separately exercise the two ordinary-gravity ownership positions

\[
 d_\theta-Gk_g \quad\hbox{(ToS stiffness)},
 \qquad
 Ga \quad\hbox{(AAF source response)}.                       \tag{H01}
\]

This is a real-data check of the non-double-counting architecture. It is not
a derivation of either coefficient from record lineage, and it is not a full
execution of `GC16`.

## 2. Exact published response equations

For ToS the source masses change the gravitational torsion constant
\(K_{1g}=G C_g\), so the exact supplementary equation is

\[
 G={\Delta\omega^2\over \Delta C_g/I}
 \left[
 1-{\Delta K\over I\Delta\omega^2}
 +{I_mK^2\over I K_m^2}
 \right].                                                    \tag{H02}
\]

Here \(\Delta C_g/I\) is the processed finite-mass source coefficient,
\(\Delta\omega^2\) is the near/far stiffness response after background
handling, \(\Delta K\) owns the fibre correction, and the last term is the
independently listed two-stage magnetic-damper correction. The stiffness
belongs in the physical operator; it may not also be entered as a source
torque.

For AAF, feedback reduces fibre twist and the inertial angular acceleration
balances the gravitational torque. The exact supplementary equation is

\[
 G={\alpha_t(2\omega_d)\over
       \left|\sum_{l=2}^{\infty}P_{g,l,2}\right|}
 \left(1+{K\over K_m}{I_m\over I}\right).                    \tag{H03}
\]

\(P_{g,l,2}\) is a processed multipole source-response coefficient; it is not
a raw mass-coordinate file. The AAF response belongs to the inhomogeneous
source/forcing channel after the displayed mechanical correction.

## 3. ToS figure-level response extraction

Official Figure-2 source workbook cells `a!B3:D22` contain ten near and ten
far three-day period summaries with source masses present. Cells
`b!B3:D22` contain the corresponding ten-plus-ten background summaries with
the source masses absent. These are raw-like intermediate observations, not
the original 0.5-second angle stream.

Convert every period to angular-frequency squared,
\(\omega^2=(2\pi/T)^2\). For every internal near-far-near and
far-near-far triple, linearly interpolate the two equal-configuration
endpoints to the middle time, form the contrast, average the 18 overlapping
contrasts, and subtract the source-absent background. This implements the
published A-B-A drift cancellation at the released resolution and gives

\[
 \Delta\omega^2_{\mathrm{ABA}}
 =1.6626945111323172\times10^{-6}\ \mathrm{s}^{-2}.           \tag{H04}
\]

The printed repeated-fibre-1 value is
\(1.662699(18)\times10^{-6}\ \mathrm{s}^{-2}\); (H04) differs by
\(-0.2494\) of that printed standard uncertainty. Because the 18 triples
share endpoints, this packet does not pretend that they are independent
observations or attach a coverage claim.

Independent panel-wise regressions, each imposing one common quadratic
ageing drift on its 20 near/far summaries, followed by source-absent
background subtraction, give

\[
 \Delta\omega^2_{\mathrm{quad}}
 =1.6626989120180067\times10^{-6}\ \mathrm{s}^{-2},           \tag{H05}
\]

only \(-0.00489\) printed standard uncertainties from the table value. This
is a deliberately labelled diagnostic rather than a reconstruction of the
authors' unavailable event-level fit or covariance. Dividing (H04) and (H05)
by the printed processed source coefficient
\(\Delta C_g/I=24912.12\ \mathrm{kg\,m^{-3}}\) gives
\(6.6742393306\times10^{-11}\) and
\(6.6742569963\times10^{-11}\) SI, respectively. Those are
response/source quotients, not new estimates of \(G\): the full \(\Delta K\)
and row-level correction custody is absent.

Across all seven Supplementary-Table-2 rows, the packet evaluates
\(\Delta\omega^2/(\Delta C_g/I)\) and the magnetic term directly. The
remaining closure residual against the authors' already-derived row summaries
is \(-5.46\) to \(-8.37\) ppm. It is explicitly an **unowned correction
bracket**, not noise set to zero and not evidence for new physics.

## 4. AAF processed-coefficient closure

Supplementary Tables 1 and 3 provide all displayed inputs in (H03) for three
campaigns. Direct substitution, without any accepted value of \(G\), gives:

| Campaign | Recomputed \(G\) (SI) | Difference from authors' row summary |
|---|---:|---:|
| AAF-I | \(6.674532777558953\times10^{-11}\) | \(-0.18315\) ppm |
| AAF-II | \(6.674375348038621\times10^{-11}\) | \(+0.05215\) ppm |
| AAF-III | \(6.674535082204988\times10^{-11}\) | \(+0.01232\) ppm |

The differences are below the rounding resolution of the displayed source,
response, and mechanical inputs. Thus the three processed-coefficient AAF
forwards close exactly at published precision.

Workbook cells `f!B3:C6`, `f!D9:E18`, and `f!F21:G35` contain 4, 10,
and 15 already-derived \(G\) outcomes, respectively. Their unweighted
campaign means are \(6.67453375\), \(6.67437500\), and \(6.67453480\) in
units of \(10^{-11}\), reproducing the three printed central summaries.
These 29 rows are not raw accelerations and are never relabelled as such.

Cells `e!B3:C10001` supply 9,999 one-second acceleration samples. On the
caption-scored first two hours (`e!B3:C7202`), a fixed \(1/600\) Hz source
harmonic plus a separately fitted lab-fixed background harmonic yields

\[
 A_{\mathrm{source}}=461.993465\ \mathrm{nrad\,s^{-2}},
 \qquad
 A_{\mathrm{background}}=77.850059\ \mathrm{nrad\,s^{-2}},   \tag{H06}
\]

at background frequency \(0.7397\) mHz. This independently demonstrates the
frequency separation described in the paper (about \(462\) versus \(77\)
\(\mathrm{nrad\,s^{-2}}\)). It is a representative figure segment, not a
campaign-average \(\alpha_t\) extraction and not a fourth AAF result.

## 5. Cross-method source-model stress

The workbook's already-derived combined outcomes are
\(6.674184(78)\times10^{-11}\) for ToS (`c!E10:F10`) and
\(6.674484(78)\times10^{-11}\) for AAF (`f!H39:I39`). Their central
separation is

\[
 3.0000\times10^{-15}\ \mathrm{SI}
 =44.9483\ \mathrm{ppm}                                      \tag{H07}
\]

relative to their midpoint. If and only if the two quoted standard
uncertainties are treated as uncorrelated, the descriptive ratio is
\(2.7196\). No cross-method covariance matrix is released, so (H07) is
retained only as a source-model/systematics stress test. It is not a
new-physics claim, a new measurement, or permission to average the two
outcomes.

## 6. Exact identifiability ceiling

The release does **not** supply:

1. the 0.5-second ToS angle/environment event stream or full 20-kHz AAF
   encoder/controller records;
2. the finite mass-coordinate and density files used to calculate
   \(\Delta C_g\) and \(P_{g,l,m}\);
3. row-level trajectory and independent global source-scale calibration with
   joint covariance;
4. complete frequency-dependent torsion, support, auxiliary-mode, and readout
   transfer;
5. signed row-level corrections and physical remainders, including complete
   ownership of \(\Delta K\);
6. observation, calibration, campaign, and cross-method covariance;
7. a conserved complete apparatus stress ledger including drives and support
   reactions; or
8. the predeclared nuisance, null, holdout, and likelihood packet required by
   the prospective GC protocol.

Therefore the public release identifies processed response/source quotients
and validates their dual operator/source placement. It does not independently
rederive the finite-source kernels, isolate a global source scale, execute the
full dressed `GC16` likelihood, confirm RGRL/GFT, infer a lineage charge, or
measure \(G\) anew.

## 7. Strongest lawful successor

Obtain or independently reconstruct the authors' finite mass-coordinate files,
complete mechanical/readout transfer, row covariance, and correction/remainder
ledger. Then recompute \(\Delta C_g\) and \(P_{g,l,m}\) from geometry and run
the full `GC16` model prospectively. Until those fields exist, further
fitting of the 7 and 29 already-derived \(G\) rows cannot close the missing
physics.
