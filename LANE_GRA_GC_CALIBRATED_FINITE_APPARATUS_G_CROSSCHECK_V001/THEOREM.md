# Calibrated finite-apparatus source model for a (G) cross-check

**Theorem ID:** `GRA-GC-CFAGC-V001`

**Date:** 2026-08-27

**Claim class:** exact finite-source Newtonian kernel inside the declared
weak-field Einstein endpoint; exact two-mode dressed torsion transfer; exact
derivative ownership and same-data zero theorem; exact product-identifiability
and global source-scale ceiling; conditional covariance/profile identified set;
synthetic executable validation.

**Status:**
`FINITE_CALIBRATED_MASS_AND_TRAJECTORY_TO_DRESSED_TORSION_OBSERVABLE_FORWARD_MODEL_CLOSED__POINT_AND_EXTENDED_SOURCE_KERNELS_EXPLICIT__BALANCED_SOURCE_HAS_ZERO_MONOPOLE_AND_DIPOLE_TANGENT_WITH_NONZERO_STF_QUADRUPOLE_TANGENT__SOURCE_REMAINDER_DATA_AND_READOUT_OWNED_ONCE__G_TIMES_UNCALIBRATED_GLOBAL_SOURCE_SCALE_ONLY_IDENTIFIABLE__INDEPENDENT_SOURCE_CALIBRATION_YIELDS_A_G_IDENTIFIED_SET__SYNTHETIC_15_OF_15_PASS__NO_G_MEASUREMENT_OR_LINEAGE_CHARGE`

## 1. Result and boundary

This lane supplies the finite-apparatus calculation requested by
`GRAVITY_RGRL_LOCAL_RESPONSE_G_IDENTIFIABILITY_CEILING_V002.md`.  It takes
nongravitationally calibrated source and detector masses, their finite spatial
profiles, a source trajectory, a calibrated torsion mode, one calibrated
coupled support mode, readout transfer, covariance, and prospectively declared
nuisances to a predicted observable

\[
 \boxed{\mathbf y_{\rm pred}=\mathbf F(G,\nu).}                 \tag{GC01}
\]

It is a physical source-calibration and endpoint-
coefficient cross-check model.  It is not a record-lineage source model.  No
lineage charge is inferred, scalar gamma is not treated as mass or stress, and
the accepted numerical value of (G) is not used to calibrate any model
coefficient.  The executable result is synthetic and is not a measurement of
(G).

The model preserves the adopted clarification:

\[
 S_{\rm dr}u=R_{\rm dr}u=\mathbb C_{dL}u=0
 \quad\Longrightarrow\quad \delta y_{\rm phys}=0.              \tag{GC02}
\]

Thus this lane does not create a force from the off-shell RGRL-C tangent.  It
constructs the independently calibrated ordinary-source lane against which a
future, separately derived lineage stress column could be tested.

## 2. Complete physical source and declared approximation

Let (T_{\rm complete}^{\mu\nu}) contain the source bodies, detector bodies,
supports, drive, and every other physical apparatus stress required by

\[
 \partial_\mu T_{\rm complete}^{\mu\nu}=0.                    \tag{GC03}
\]

An accelerated point mass by itself is not a conserved source.  In a real
application the drive and support stresses must therefore be included in the
complete retained system or bounded in the physical remainder.  This lane
scores settled, weak-field, nonrelativistic response after the source
trajectory and support state have been independently measured.  The exact
linearized retarded metric kernel is

\[
 \bar h_{\mu\nu}(t,\mathbf x)
 ={4G\over c^4}\int
 {T^{\rm complete}_{\mu\nu}
  (t-|\mathbf x-\mathbf X|/c,\mathbf X)
  \over |\mathbf x-\mathbf X|}\,d^3X .                         \tag{GC04}
\]

The executable torsion model takes its quasistatic Newtonian component.  Terms
from source velocity, stress active mass, retardation, higher weak-field
orders, support motion, and unmodeled modes are not silently discarded: they
belong to the declared remainder and receive a bound in an actual protocol.

## 3. Exact finite-source kernels

Let (d\mu_s(q)) and (d\mu_d(r)) be independently measured, compact,
positive Lagrangian source and detector mass measures with fixed weights and
disjoint spatial supports.  Their calibrated embedding maps are
(\mathbf X(q;u)) and (\mathbf x(r)).  The detector is rigidly rotated by
(R_z(\theta)), and (u) is a calibrated small source-trajectory coordinate.
This fixed-weight transport assumption is what makes the trajectory derivative
below complete.  At the leading Einstein/Newtonian order,

\[
 U_N(\theta,u)
 =-G\iint
 {d\mu_s(q)d\mu_d(r)
  \over |\mathbf X(q;u)-R_z(\theta)\mathbf x(r)|}.              \tag{GC05}
\]

Define the exact geometry function and the two load-bearing derivatives

\[
 K(\theta,u):=-{1\over G}{\partial U_N\over\partial\theta},
 \qquad
 a:={\partial K\over\partial u}\bigg|_{0,0},
 \qquad
 k_g:={\partial K\over\partial\theta}\bigg|_{0,0}.             \tag{GC06}
\]

For discrete calibrated mass elements (M_i) at (mathbf R_i(u)) and (m_j)
at (mathbf r_j), write

\[
 \mathbf d_{ij}=\mathbf R_i-\mathbf r_j,
 \quad n_{ij}=\hat{\mathbf z}\!\cdot
 (\mathbf r_j\times\mathbf R_i),
 \quad \mathbf V_i={d\mathbf R_i\over du}.
\]

Then the point-element kernel and source derivative are exactly

\[
 K=\sum_{ij}M_i m_j{n_{ij}\over |\mathbf d_{ij}|^3},           \tag{GC07}
\]

\[
 a=\sum_{ij}M_i m_j\left[
 {\hat{\mathbf z}\cdot(\mathbf r_j\times\mathbf V_i)
  \over|\mathbf d_{ij}|^3}
 -{3n_{ij}(\mathbf d_{ij}\cdot\mathbf V_i)
  \over|\mathbf d_{ij}|^5}
 \right].                                                       \tag{GC08}
\]

With (mathbf r'_j=\hat{\mathbf z}\times\mathbf r_j), the gravitational
torsion stiffness is

\[
 k_g=\sum_{ij}M_i m_j\left[
 {\hat{\mathbf z}\cdot(\mathbf r'_j\times\mathbf R_i)
  \over|\mathbf d_{ij}|^3}
 +{3n_{ij}(\mathbf d_{ij}\cdot\mathbf r'_j)
  \over|\mathbf d_{ij}|^5}
 \right].                                                       \tag{GC09}
\]

Equations (GC05)--(GC09) are also the exact extended-source kernel: replace
the sums by the calibrated Lagrangian mass integrals.  A finite element cloud is exact
for the declared discrete mass measure.  For mutually nonoverlapping uniform
spheres, Newton's shell theorem reduces every pair exactly to (GC07)--(GC09)
at its centers.  No point-source approximation is then being made.

The executable uses two equal antipodal source masses with fixed weights.  For
their trajectory tangent (\mathbf V_i), the source contrast obeys exactly

\[
 \delta M=0,
 \qquad
 \delta D_a=\sum_i M_iV_{ia}=0,                               \tag{GC09a}
\]

while its symmetric trace-free quadrupole tangent

\[
 \delta Q_{ab}=\sum_iM_i\left[
 3(V_{ia}R_{ib}+R_{ia}V_{ib})
 -2(\mathbf R_i\!\cdot\!\mathbf V_i)\delta_{ab}\right]       \tag{GC09b}
\]

is trace-free and nonzero.  Consequently the Fourier transform of the
Newtonian mass-density tangent has no (k^0) or (k^1) term and begins at
(O(|\mathbf k|^2)).  This finite profile avoids treating the calculation as a
literal density-(k=0) inversion; it does not establish the full conserved
(T^{\mu\nu}) zero-mode/range conditions (GI25), which remain an explicit
real-apparatus requirement.  This is not a claim to have derived general
finite-(k) gravity.  The symmetric torsion read is likewise differential:
constant potential and common translation do not enter (K).  A relativistic
application must express the same read in a detector-frame tidal observable.

## 4. Full coupled and dressed apparatus operator

Use Fourier amplitudes proportional to (e^{+i\omega t}).  Let the calibrated
bare torsion and auxiliary-mode denominators be

\[
 d_\theta(\omega)=\kappa-I\omega^2+i b\omega,
 \qquad
 d_x(\omega)=\kappa_x-I_x\omega^2+i b_x\omega,                 \tag{GC10}
\]

and let (lambda) be their independently calibrated reciprocal coupling.
The auxiliary coordinate (x) is a calibrated dimensionless generalized mode
normalized so that (d_x) and the reciprocal off-diagonal (lambda) have the
displayed stiffness units.  With another coordinate convention, the two
off-diagonal conversions must be retained separately and (\lambda^2/d_x)
replaced by their calibrated product over (d_x).
If the global source-mass calibration is (s>0), define

\[
 p:=Gs.                                                         \tag{GC11}
\]

Because both (a) and (k_g) are linear in the source measure, the exact
linearized two-mode apparatus problem is

\[
 \begin{pmatrix}
 d_\theta-pk_g&-\lambda\\
 -\lambda&d_x
 \end{pmatrix}
 \binom{\theta}{x}
 =\binom{pa+r_\theta}{r_x}.                                    \tag{GC12}
\]

Here the metric has already been eliminated through the leading Einstein
kernel.  Consequently (pa) is the metric-mediated torque, not a second copy
of the material source.  The (pk_g) term is the implicit gravitational
stiffness in the operator.  It may not also be added to the explicit source
column.

Eliminating the calibrated auxiliary mode gives the exact Schur complement

\[
 D_{\rm dr}(\omega;p)
 =d_\theta(\omega)-pk_g
 -{\lambda^2\over d_x(\omega)},                                 \tag{GC13}
\]

\[
 r_{\rm dr}=r_\theta+{\lambda r_x\over d_x},
 \qquad
 \theta_{\rm part}={pa+r_{\rm dr}\over D_{\rm dr}}.           \tag{GC14}
\]

This is the finite-apparatus realization of V002's (A_{\rm dr}^R),
(S_{\rm dr}), (R_{\rm dr}), and retarded inverse.  With readout transfer

\[
 C(\omega)=g_{\rm ro}e^{-i\omega\tau_{\rm ro}},                \tag{GC15}
\]

homogeneous initial/boundary/branch response (d_h), linear readout nuisance
templates (B\eta), and detector noise (epsilon), the complete observable
is

\[
 \boxed{
 y(\omega)
 =C(\omega)\left[
 {pa+r_\theta+\lambda r_x/d_x\over
  d_\theta-pk_g-\lambda^2/d_x}
 +d_h\right]
 +B(\omega)\eta+\epsilon(\omega).}                             \tag{GC16}
\]

The static quotient is stable and has no zero mode in the declared two-mode
sector if

\[
 \kappa-pk_g>0,qquad
 (\kappa-pk_g)\kappa_x-\lambda^2>0.                            \tag{GC17}
\]

An actual scan must check (GC17) over its entire admitted (p) and geometry
domain rather than only at a fitted value.

## 5. Non-double-counted ownership

| Physical object | Sole location in (GC12)--(GC16) |
|---|---|
| calibrated source mass, geometry, and trajectory | (a), multiplied once by (p=Gs) |
| gravitational torque gradient | (-pk_g) inside (D_{\rm dr}) |
| implicit support/detector feedback | (-\lambda^2/d_x) inside (D_{\rm dr}) |
| explicit physical unmodeled torques | (r_\theta,r_x), before the retarded inverse |
| initial, incoming, boundary, or solution-branch change | (d_h), homogeneous solution column |
| detector/readout/analysis nuisance | (B\eta+\epsilon), after the physical solution |

This table implements V002's derivative ownership exactly once.  In
particular, an electromagnetic force, support displacement, heat load, gas
effect, or controller reaction is physical and cannot be relabeled as
post-solution readout noise.

## 6. Theorem GC-1: same-data zero

On the stable quotient (GC17), if

\[
 a=0,qquad r_\theta=r_x=0,qquad d_h=0,                         \tag{GC18}
\]

then the physical response in (GC16) is exactly zero for every (G), even
though an independently adopted off-shell metric tangent may have full rank.

**Proof.**  The inhomogeneous column and homogeneous solution data vanish.
The invertible matrix in (GC12) therefore has the unique zero solution.
Readout noise can create an observed residual but not a physical torque or
metric column.  QED.

## 7. Theorem GC-2: what the apparatus identifies

First suppose (s), (a_n), (k_{g,n}), (d_{\theta,n}), (d_{x,n}),
(lambda), and (C_n) are fixed independently, and all remainder, data, and
nuisance columns are either zero or known.  Define

\[
 d_{0n}:=d_{\theta,n}-{\lambda^2\over d_{x,n}},
 \qquad
 F_n(p)=C_n{p a_n\over d_{0n}-pk_{g,n}}.                        \tag{GC19}
\]

If at least one scored row has (C_na_nd_{0n}\ne0) and its denominator does
not vanish, that row is globally injective in positive (p).  Indeed,

\[
 {p_1a_n\over d_{0n}-p_1k_{g,n}}
 ={p_2a_n\over d_{0n}-p_2k_{g,n}}
 \Longrightarrow (p_1-p_2)d_{0n}=0
 \Longrightarrow p_1=p_2.                                     \tag{GC20}
\]

Thus an exact noiseless calibrated observation identifies (p).  It
identifies (G) only when (s) is independently fixed.  Without that source
normalization,

\[
 \boxed{F(G,s)=F(Gq,s/q)\quad(q>0),}                            \tag{GC21}
\]

so no amount of same-apparatus gravity data separates (G) from a free global
source-mass scale.  This is an exact, not statistical, nonidentifiability.

If a positive product interval ([p_-,p_+]) is identified and independent
mass metrology supplies (s\in[s_-,s_+]), then the exact quotient set is

\[
 \boxed{
 G\in\left[{p_-\over s_+},{p_+\over s_-}\right].}             \tag{GC22}
\]

This is the finite-apparatus specialization of V002 (GI29)--(GI32).

## 8. Nuisance and covariance treatment

For complex lock-in data, stack real and imaginary parts into one real vector.
Let (Sigma_y) be the prospectively frozen observation covariance.  For
small independent calibration errors (delta\nu) with covariance
(Sigma_\nu), first-order propagation is

\[
 \Sigma_{\rm eff}(p)
 =\Sigma_y+J_\nu(p)\Sigma_\nu J_\nu(p)^T.                       \tag{GC23}
\]

The same calibration uncertainty may appear only once.  If the linearization
is inadequate, the nuisance must instead be retained in the full implicit
forward model.  If (Sigma_{\rm eff}) depends materially on (p), the
Gaussian likelihood must include both its quadratic form and log determinant;
holding it fixed after seeing the gravity result is not allowed.

In the executable, the gain/delay Jacobian is frozen at the predeclared
absolute reference (p_{\rm cov}=7.0\times10^{-11}), not at the synthetic truth
and not at a fitted value.  This makes (GC23) a prospective first-order noise
model rather than a hidden truth-dependent weight.

For frozen covariance and linear nuisance templates, define

\[
 \chi^2(p,\eta)
 =[\mathbf y-\mathbf F(p)-B\eta]^T
 \Sigma_{\rm eff}^{-1}
 [\mathbf y-\mathbf F(p)-B\eta]
 +\eta^T\Lambda\eta.                                            \tag{GC24}
\]

The code profiles (eta) analytically at every (p).  Global
identifiability with an unpenalized nuisance requires

\[
 \mathbf F(p_1)-\mathbf F(p_2)\notin\operatorname{col}B
 \quad\hbox{for every }p_1\ne p_2                               \tag{GC25}
\]

in the admitted scan.  A likelihood interval derived from a chosen
(Delta\chi^2) is a statistical interval, not an exact theorem.  For bounded
remainders or non-Gaussian calibration sets, the authoritative object remains
the V002 identified set

\[
 \mathcal I_G=\left\{G>0:\exists\nu\in\mathcal N,
 \|\mathbf y-\mathbf F(G,\nu)\|_W\le\epsilon\right\}.           \tag{GC26}
\]

## 9. Executable synthetic validation

`verify_calibrated_finite_apparatus_g.py` uses an arbitrary generator-only
truth (G_{\rm syn}=7.314159265358979\times10^{-11}) SI.  This deliberately
is not an accepted reference value.  The estimator receives only the
synthetic observations, calibrated apparatus, covariance, and scan domain.
The absolute scan (p\in[10^{-12},1.5\times10^{-10}]) with 6001 points is
declared independently of the generator, and (GC17) is checked at both
endpoints for every geometry row.  Because its determinant is affine in (p),
the endpoint check covers the entire scan.  Truth is used only after fitting
for coverage and error diagnostics.

The deterministic run passes `15/15` checks, including:

1. exact zero monopole/dipole and nonzero STF quadrupole source tangent;
2. analytic point torque versus the potential derivative;
3. analytic source and torsion-gradient derivatives;
4. the exact nonoverlapping uniform-sphere reduction;
5. stable quotient over the full absolute scan;
6. equality of the full coupled solve and Schur complement;
7. the same-source/remainder/data zero theorem;
8. exact (G\)-source-scale degeneracy;
9. covariance propagation without duplicate ownership or truth-dependent
   weights;
10. blinded product recovery and nominal profile-grid interval coverage;
11. propagation of the independent source-scale interval to (G);
12. held-out geometry/frequency prediction without refitting; and
13. a detectable `0.5429%` synthetic calibration shift when both the
    auxiliary Schur term and gravitational stiffness are incorrectly omitted.

For the frozen synthetic realization,

\[
 p_{\rm true}=7.318547760918193\times10^{-11},
\]

\[
 p_{95\%,\rm grid}\in
 [7.311600000000000,7.324016666666666]\times10^{-11},           \tag{GC27}
\]

and (s\in[0.999,1.001]) gives

\[
 G_{95\%,\rm grid}\in
 [7.304295704295705,7.331348014681347]\times10^{-11}.           \tag{GC28}
\]

These are nominal profile-grid, not exact-coverage, intervals.  The held-out
Mahalanobis value is a deterministic implementation diagnostic, not a
calibrated independent-study p-value.  The numbers validate the estimator
implementation only and carry no empirical weight.

## 10. What this closes and what remains

This lane closes the bounded no-lab calculation from a finite calibrated
ordinary source through an implicit, dressed apparatus operator to a (G)
identified set.  It proves that an uncalibrated source normalization cannot be
mistaken for a derivation of (G), and it provides a concrete forward model
for an actual torsion dataset when one is supplied.

It does **not**:

- measure or derive the realized numerical (G);
- derive a lineage stress-energy functional;
- prove GI21's lineage-to-physical-source type join;
- turn the off-shell RGRL-C tangent into an on-shell force;
- empirically confirm RGRL, EIR, GFT, or universal coupling; or
- replace a complete conserved-source and systematics model for a real
  apparatus.

A future lineage-(G) cross-check may reuse (GC05)--(GC26) only after a
lineage-dependent complete stress column has been derived and calibrated
nongravitationally.  Until then, the correct matched-lineage physical
prediction remains (GC02).

## 11. Custody

This theorem uses, without modifying:

- `GRAVITY_RGRL_LOCAL_RESPONSE_G_IDENTIFIABILITY_CEILING_V002.md`, SHA-256
  `140eb379f22c369b0442b5585a4828122ac4f5858b54ae3f2bf6391866ac84a3`;
- `GRAVITY_RGRL_ONSHELL_OFFSHELL_AND_SPAG_CLARIFICATION_V001.md`, SHA-256
  `50e2be7f06e79d943eddb6c37c63094dd758bb8b6798987e73b8eacf8ef448df`;
- `GRAVITY_RGRL_ONSHELL_OFFSHELL_CLARIFICATION_ADOPTION_V001.md`, SHA-256
  `4959f99898b216edc7da3e212ce2e26422287899fcf8f3b41cd34ef5d8bb3ff8`;
- `GRAVITY_SPAG_EIR_CONSISTENCY_AND_EMPIRICAL_PATH_V001.md`, SHA-256
  `1bab3e5901a1566d287a2a7b04563718fbfab321b75ba41e48aa6cd152862ae1`;
- `GRAVITY_RGRL_IR_ENDPOINT_CLOSURE_THEOREM_V001.md`, SHA-256
  `c883c4c9f3816e453766846a1691ef27cb50d6ea7e5676bc52ed1928617f82bf`.

## Disposition

`FINITE_CALIBRATED_SOURCE_GEOMETRY_AND_TRAJECTORY_DEFINE_EXACT_NEWTONIAN_TORQUE_SOURCE_AND_STIFFNESS_KERNELS__THE_COUPLED_APPARATUS_REDUCES_TO_ONE_EXACT_RETARDED_SCHUR_COMPLEMENT__SOURCE_REMAINDER_HOMOGENEOUS_DATA_AND_READOUT_ARE_OWNED_ONCE__THE_FORWARD_MAP_IS_GLOBALLY_INJECTIVE_IN_P_EQUALS_G_TIMES_S_ON_ITS_STABLE_DOMAIN__A_FREE_GLOBAL_SOURCE_SCALE_IS_EXACTLY_DEGENERATE_WITH_G__INDEPENDENT_SOURCE_CALIBRATION_CONVERTS_THE_PRODUCT_SET_TO_A_G_SET__SYNTHETIC_VALIDATION_PASSES__NO_EMPIRICAL_OR_LINEAGE_CLAIM`
