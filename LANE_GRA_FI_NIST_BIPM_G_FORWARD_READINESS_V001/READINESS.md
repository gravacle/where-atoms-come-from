# NIST/BIPM 2026 public-apparatus G-forward readiness audit

**Lane ID:** `GRA-FI-NIST-BIPM-GFR-V001`

**Date:** 2026-08-27

**Claim class:** exact primary-source custody; exact public-summary reduction
to the committed finite-apparatus torque forward; exact Jacobian and nuisance
rank audit; reproducible uncertainty/correlation diagnostics; exact missing-
field ledger

**Disposition:**
`PUBLIC_SUMMARY_REDUCED_TORQUE_FORWARD_EXACT__FINITE_CONTRAST_SOURCE_COLUMN_AND_NOMINAL_ZERO_SOURCE_STIFFNESS_OWNED__TABLE15_RATIO_AND_TABLE18_COVARIANCE_DIAGNOSTICS_REPRODUCED__FULL_GC16_REAL_APPARATUS_FIT_NOT_READY__NO_ACCEPTED_G_INPUT_NO_LINEAGE_OR_GRAVITY_EMERGENCE_CLAIM`

## 1. Result

The committed NIST/BIPM 2026 paper supports a real, apparatus-specific
**reduced torque forward map**:

\[
 \boxed{
 \Delta N_j=G A_j+r_j,
 \qquad
 A_j={16\Gamma_jm_{s,j}m_t\over R_{s,j}}.}           \tag{FI01}
\]

Here `Delta N_j` is the paper's already calibrated two-position torque
difference, and `A_j` is the finite two-state contrast analogue of the source
column `a` in the committed GC model.  It is not silently called the
infinitesimal trajectory derivative in GC06.  Table 15's four configuration
rows supply eight torque-mode observations; represented as eight analysis
rows, they give `A_j=209.4792946 kg^2/m` for sapphire and
`467.2522461--467.4580783 kg^2/m` for copper.  Thus

\[
 {\partial\Delta N_j\over\partial G}=A_j\ne0.         \tag{FI02}
\]

At the ideal source-mass torque extrema, the paper states that the source-mass
contribution has zero second derivative.  Therefore the nominal
source-induced GC torsion-gradient entry is

\[
 k_{g,j}=0                                             \tag{FI03}
\]

at those evaluation points.  Finite dither, miscentering, and support/local-
gravity corrections to (FI03) are not supplied by the public packet.

This is substantive progress: the public apparatus can now be attached to the
**source-column side** of the finite-apparatus G model without a synthetic
mass geometry or an accepted value of `G`.  It is not yet possible to execute
the full dressed GC16 likelihood or claim an independent real-`G` estimate,
because the paper publishes calibrated summaries rather than the complete raw
geometry, transfer, covariance, and physical-remainder packet required by the
prospective protocol.

## 2. Primary-source custody

The sole empirical source is the committed official open-access PDF:

- `LANE_GRA_SPAG_PUBLIC_DATA_SUBSTITUTE_V001/SOURCE/nist_bipm_2026.pdf`
- SHA-256
  `c79552d62f4d4f4e85cfbbb00f135c1d985b596d9cdcde9bee57cfe4618f33dc`
- 31 PDF pages
- Schlamminger et al., *Metrologia* **63** (2026) 025012,
  “Redetermination of the gravitational constant with the BIPM torsion
  balance at NIST.”

Relevant pages were text-extracted and visually rendered from the pinned
bytes.  PDF-page numbers are one greater than the article's printed page in
this file.

| Public object | Exact custody |
|---|---|
| four-source/four-test geometry and torque law | PDF p. 6 / printed p. 5, Figure 3 and equations (1)--(5) |
| free-deflection transfer and zero source curvature at torque extrema | PDF p. 7 / printed p. 6, equations (6)--(10) and text below equation (7) |
| method-dependent relative sensitivities | PDF p. 10 / printed p. 9, Table 1 and equation (21) |
| nominal masses, dimensions, radii, period, capacitances, and capacitance gradients | PDF p. 11 / printed p. 10, Table 2 |
| mass-integration geometry sensitivity and uncertainty | PDF p. 15 / printed p. 14, Table 6 and equations (26)--(33) |
| full/partial mass-integration comparisons and density orientation | PDF p. 18 / printed p. 17, Table 10 and sections 6.7--6.9 |
| background-torque and moment-of-inertia checks | PDF p. 23 / printed p. 22, Tables 12--13 and sections 8.2--8.4 |
| autocollimator non-linearity summary | PDF p. 25 / printed p. 24, Table 14 and equation (62) |
| eight torque summaries and four derived G summaries | PDF p. 26 / printed p. 25, Tables 15--16 |
| category uncertainties and four-result correlations | PDF p. 27 / printed p. 26, Tables 17--18 and equation (63) |

No paper figure, abstract, or prose summary is promoted to event-level raw
measurement.  Table 15 is explicitly treated as eight published calibrated
torque summaries with Type-A uncertainties only.

## 3. Exact public parameter extraction

From Tables 2 and 15:

\[
 m_t=1.150156\ {\rm kg},\qquad
 m_s^{\rm Cu}=11.19168\ {\rm kg},\qquad
 m_s^{\rm Sa}=5.01560\ {\rm kg}.                    \tag{FI04}
\]

The row-specific `R_s`, `Gamma`, torque, and Type-A uncertainty give:

| Configuration | Mode | `A_j` (`kg^2/m`) | `Delta N` (`nN m`) | summary ratio `Delta N/A_j` (SI) |
|---|---:|---:|---:|---:|
| Sapphire | free | 209.479294635 | 13.9799(2) | `6.673642865e-11` |
| Sapphire | servo | 209.479294635 | 13.9778(3) | `6.672640379e-11` |
| Copper 0 deg | free | 467.458078302 | 31.1979(3) | `6.673946060e-11` |
| Copper 0 deg | servo | 467.458078302 | 31.1962(4) | `6.673582391e-11` |
| Copper 120 deg | free | 467.256225338 | 31.1842(3) | `6.673897170e-11` |
| Copper 120 deg | servo | 467.256225338 | 31.1828(6) | `6.673597549e-11` |
| Copper 240 deg | free | 467.252246069 | 31.1856(2) | `6.674253631e-11` |
| Copper 240 deg | servo | 467.252246069 | 31.1836(3) | `6.673825597e-11` |

These ratios are a deterministic reconstruction of the paper's summary
algebra, not a new estimate.  The denominator is an already-integrated
published `Gamma`, not an independently reconstructed finite mass measure.

Table 1 additionally supplies the relative sensitivities of inferred `G`:

| input | free | servo |
|---|---:|---:|
| `m_s` | -1 | -1 |
| `m_t` | -0.1 | -1 |
| `R_s` | -5.4 | -5.4 |
| `R_t` | +2.4 | +4.4 |
| `I_disk` | +0.11 | 0 |
| `phi_t` | +1 | -1 |

The mode contrast is physically important: calibration errors cannot be
owned by one generic scalar nuisance without checking which method they enter.

## 4. Exact GC ownership ledger

| GC object | NIST/BIPM public mapping | Sole ownership in the reduced forward | Public status |
|---|---|---|---|
| source measure `dmu_s` | four source cylinders; average masses/dimensions in Table 2; row `R_s,Gamma` in Table 15 | inside `A_j` once | nominal aggregate only; no source point cloud or density file |
| detector measure `dmu_d` | four test cylinders plus the 299-shape disk model summarized in Tables 7--10 | inside `A_j` once | test-cylinder averages public; complete coordinate/shape file absent |
| source trajectory | two azimuthal extrema separated by 37.67 deg, located by dither | finite contrast defining `A_j` | exact run-level positions, dither path, and timestamps absent |
| metric-mediated source torque | `G A_j` | inhomogeneous physical source column only | summary-level exact |
| source gravitational stiffness `G k_g` | source second derivative vanishes at ideal extrema | dressed operator only | nominal `k_g=0`; finite-offset correction absent |
| bare torsion transfer `d_theta` | period, inertia, flexure, anelasticity, and free fit | physical response operator | only partial/calibrated summaries; no raw row transfer |
| auxiliary/support transfer `d_x,lambda` | gimbal, support, damping, controller, and apparatus modes | Schur term only | no GC-normalized auxiliary calibration supplied |
| readout `C` | autocollimator for free mode; capacitance/voltage servo calibration for servo mode | after physical solution once | Table summaries only; raw transfer/calibration covariance absent |
| physical torque remainder `r_theta` | gas thermal torque, empty-carousel/disk residual, local gravity, drive/support reaction, electromagnetic and servo heating effects | before response inverse | categories and some bounds public, signed row columns absent |
| auxiliary remainder `r_x` | unmodeled support/controller forcing | auxiliary physical column | absent |
| homogeneous data `d_h` | free oscillation/initial state and finite-gain servo residual | homogeneous solution column | fitted into published torque; raw data absent |
| readout nuisance/noise `B eta + epsilon` | Type-A scatter and analysis/readout nuisances | after physical solution | row variances only; nuisance templates absent |
| observation covariance `Sigma_y` | eight Table-15 torque summaries | likelihood only | **not published**; Table 18 is for four already-derived G values |
| calibration covariance `Sigma_nu` | 23 inputs grouped into Table-17 categories | propagated once or explicit nuisance | full covariance and calibration records absent |
| global source scale `s` | nongravitational mass calibration | `p=Gs`, never duplicated in `A_j` | nominal masses and category uncertainties public; independent scale set/covariance absent |
| complete conserved apparatus source | masses, vacuum can, carousel, drive, supports, fields, gas, controller | physical source/remainder ledger | not released as a complete stress/source model |

The published `Delta N` values already contain method-specific transfer and
calibration.  Setting `C=1` in (FI01) means “take published torque as the
observation”; it does **not** assert that the raw autocollimator or voltage
readout has unit gain.

## 5. Jacobian and identifiability results

With calibration and remainders fixed, the eight-by-one source Jacobian has
rank one and any nonzero row identifies the scalar product entering (FI01).
The exact ceilings appear when physical ownership is restored:

1. With a free global source scale, the two parameter columns for `(G,s)` are
   collinear.  Their Jacobian rank is one, so only `p=Gs` is identified.
2. With one arbitrary physical torque remainder per row, the design
   `[A | I_8]` has rank 8 rather than 9; `A` lies in the nuisance span and `G`
   is not identifiable.
3. Even one arbitrary remainder per **configuration** is enough.  The free
   and servo rows share the same `A_j` within each configuration, so `A` is
   exactly a linear combination of the four configuration indicators.  The
   five-column design has rank four.
4. If one makes the much stronger assumption of only two common additive
   method offsets, `[A | B_free | B_servo]` has rank three and `G` is
   algebraically identifiable.  The paper does not justify that restriction:
   it reports configuration-dependent uncertainty and an unexplained
   free/servo discrepancy.

The last point is why the missing remainder ledger is causal, not clerical.
The observed free-minus-servo differences are `1.4--2.1 pN m`; assigning them
to the wrong layer can manufacture precision.

## 6. Covariance and sensitivity diagnostics actually supported

The analyzer reproduces Table 17's four combined relative uncertainties by
root-sum-square to within `0.064 ppm`, the expected printed rounding.  The
four central values are the Table-16 values.  Its displayed uncertainty
column is rounded to whole ppm `(23,30,38,94)`, whereas the covariance below
uses the tenth-ppm standard uncertainties `(23.2,30.3,37.5,93.9)` printed on
the Table-17 combined row and Table-18 diagonal.  This source distinction is
retained explicitly in the result schema.

The Table-18 four-result correlation matrix is positive definite with eigenvalues

\[
 (0.5000073,\ 0.7657793,\ 0.9020189,\ 1.8321945),    \tag{FI05}
\]

rank four, and condition number `3.66434`.  It can therefore be used for a
diagnostic on the four **already-derived** Table-16 G summaries.  A formal GLS
calculation gives weights

\[
 (0.617817,\ 0.254323,\ 0.126685,\ 0.001175)         \tag{FI06}
\]

and a formal standard uncertainty `1.4154e-15 SI` (`21.21 ppm`).  This is not
the covariance of the eight torque rows and excludes the paper's subsequent
dark-uncertainty layer.

For comparison only, pretending that the eight Table-15 Type-A errors are
independent, all denominators exact, and all remainders zero gives the formal
diagonal diagnostic

\[
 G_{\rm diag}=6.673920939\times10^{-11},\qquad
 u_{\rm diag}=2.4773\times10^{-16}\quad(3.712\ {\rm ppm}).     \tag{FI07}
\]

Equation (FI07) is deliberately labelled
`ALGEBRA_DIAGNOSTIC_ONLY__NOT_A_G_ESTIMATE`.  The gap between 3.7 ppm and the
paper's 23--94 ppm per-result combined uncertainties quantitatively shows why
Table-15 Type A values cannot stand in for the missing covariance and physical
remainders.

The paper's Table-19/final-consensus hierarchy is not used anywhere in this
packet because its stated prior is centered on CODATA.  No accepted `G`, paper
consensus, or synthetic `G` enters extraction, calibration, weights, ranks, or
verification.

## 7. Exact fields preventing a full real-G execution

The committed GC protocol cannot be executed on this public PDF until the
following are released or independently reconstructed:

1. event-level free-deflection angle/time/frequency data with run and
   source-position labels;
2. event-level servo voltage, capacitance-gradient, residual-twist,
   controller, and timing data;
3. complete source/test/disk finite-element mass-coordinate and density files
   used to compute `Gamma` and its trajectory derivatives;
4. row-level source/test mass and geometry calibration records, their joint
   covariance, and an independent global source-scale interval;
5. full mechanical transfer: total inertia, frequency-dependent torsion
   response, damping, gimbal/support modes, and auxiliary-mode couplings;
6. the raw autocollimator calibration curve/ensemble and row-level readout
   transfer covariance;
7. the eight-row torque observation covariance--Table 18 is only a four-by-
   four covariance summary of already-derived G values;
8. signed row-level physical remainder columns or bounds for gas thermal,
   background, disk, local-gravity, drive/support, electromagnetic, and
   controller effects;
9. the complete conserved apparatus source/stress ledger and finite-dither
   correction to nominal `k_g=0`; and
10. the predeclared nuisance design, held-out rows, likelihood domain, and
    covariance treatment required by the prospective GC protocol.

## 8. Scientific status

This packet closes a real no-lab increment: one public apparatus now has a
source-custodied, unit-correct, parameter-owned reduced forward and an exact
readiness obstruction.  It proves that the missing step is no longer the
ordinary source coefficient; it is the raw transfer/covariance/remainder
completion needed to test that coefficient without aliasing it.

It does not provide an independent numerical `G` estimate, test a lineage
source, confirm RGRL or Gravity Formation Theory, or derive gravity.  It also
does not weaken the finite-apparatus model: the public paper passes the source-
column interface and fails only where the prospective protocol intentionally
requires information that a paper-level summary cannot substitute for.

## 9. Reproduction

Run:

```text
python3 LANE_GRA_FI_NIST_BIPM_G_FORWARD_READINESS_V001/analyze_nist_bipm_g_readiness.py
```

The zero-argument run verifies all dependency hashes and page-level custody,
recomputes the source columns, ranks, Table-17 RSS values, and Table-18
covariance diagnostics, and returns `24/24`.  `--json` reproduces
`RESULT.json` byte-for-byte.
