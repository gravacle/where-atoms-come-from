# Point-local lineage-response / Newton-coefficient identifiability theorem

**Theorem ID:** `RGRL-PLGI-V001`

**Date:** 2026-08-27

**Claim class:** exact linear identifiability theorem at the already earned
point-local / zero-spatial-momentum scope; exact separation of constitutive
lineage response from Einstein source response; exact minimal-input statement
for inferring the Einstein--Hilbert coefficient or Newton's constant.

**Status:**
`POINT_LOCAL_PHYSICAL_LINEAGE_RESPONSE_EQUALS_ELL_F_SQUARED_TIMES_THE_S4_FORM_FACTOR__COMPATIBLE_O3_GIVES_ONE_TRACE_AND_ONE_SHEAR_RESPONSE__EINSTEIN_HILBERT_FIXES_THE_METRIC_OPERATOR_BUT_NOT_THE_LINEAGE_SOURCE_OR_BRANCH_MAP__G_IS_IDENTIFIABLE_FROM_H_ONLY_AFTER_INDEPENDENT_COMPLETE_STRESS_ZERO_HOMOGENEOUS_DATA_REMAINDER_AND_GAUGE_CALIBRATION__MATCHED_KEEP_BREAK_TESTS_RGRL_ANCESTRY_BUT_CANNOT_CALIBRATE_G`

**Not claimed:** a numerical derivation of `G`; a new finite-momentum
classification; a universal ratio between trace and shear; that a lineage
intervention is automatically a stress source; that the matched SPAG contrast
measures `G`; or that a gauge/constraint zero mode can be inverted without an
independent boundary and observable prescription.

## 1. Frozen scope and typed objects

Work at the q4 tetrahedral symmetric point and at the derivative-zero or
`k=0` scope of `RGRL-LMRK-V001`.  Let

\[
 X:=D\,\delta L\in\operatorname{Sym}^2(V)          \tag{PL01}
\]

be the tensorized lineage-intervention direction.  This is only a convenient
transport of the source coordinate through the exact invertible kinematic map
`D`; it does not make `X` a stress tensor.  The physical metric response is

\[
 \delta s_{\rm sp}
 =\ell_F^2D H^R\delta L.                           \tag{PL02}
\]

For `r in {A,E,T}`, Schur reduction gives

\[
 P_r\delta s_{\rm sp}
 =\underbrace{\ell_F^2 h_r^R}_{\displaystyle\kappa_r^R}
   P_rX,
 \qquad
 \boxed{\kappa_r^R=\ell_F^2h_r^R}.                \tag{PL03}
\]

Thus the singular weights of `D` do not appear in the physical response
eigenvalue when both input and output are expressed in the same transported
tensor coordinate `X`.  They appear only if one compares an uncalibrated
Euclidean norm on edge coordinates with the Frobenius norm on metric tensors.

If the compatible physical `O(3)` action required by `LM26a` is present on
both source and output, then exactly

\[
 \boxed{h_E^R=h_T^R=:h_{\rm sh}^R,
 \qquad \kappa_E^R=\kappa_T^R=:\kappa_{\rm sh}^R,} \tag{PL04}
\]

while `h_A=h_tr` remains an independent trace response.  This is the strongest
relation among `(h_A,h_E,h_T)` supplied by the present symmetry and
kinematics alone.

## 2. The endpoint equation and the three distinct lineage channels

Linearize the closed endpoint equation on one admitted background:

\[
 \delta\!\left(G_{\mu\nu}+\Lambda_{\rm eff}g_{\mu\nu}\right)
 ={8\pi G_{\rm eff}\over c^4}\,
   \delta T_{\mu\nu}^{\rm complete}
 +\delta\Delta_{\mu\nu}^{\rm rem}.                \tag{PL05}
\]

Fix the same gauge/boundary quotient, zero-mode prescription, and retarded
solution class used to define the physical response.  On every projected
sector for which that problem has a unique retarded inverse, denote the
unit-normalized linearized Einstein solution map by `mathscr G_r^R`.  It
inverts the geometric left side of (PL05); the factor `8 pi G/c^4` is not
included in its definition.

A lineage intervention can enter the solution in three logically different
ways:

1. through a complete, conserved physical stress variation
   `delta T_complete = mathscr S_r delta X_r`, where `mathscr S_r` is an
   independently calibrated lineage-to-stress map;
2. through homogeneous, initial, boundary, or solution-branch data, denoted
   by the metric-response column `b_r^R`; and
3. through the uniquely owned remainder response, denoted by `rho_r^R`.

The retarded linear solution therefore has the exact projected form

\[
 \boxed{
 \kappa_r^R
 =b_r^R
 +{8\pi G_{\rm eff}\over c^4}\,\Phi_r^R
 +\rho_r^R,}
 \qquad
 \Phi_r^R:=\mathcal P_r\mathscr G_r^R\mathscr S_r. \tag{PL06}
\]

Here `mathcal P_r` includes the declared output projection and normalization
relative to `P_rX`.  Equation (PL06) is a decomposition, not an assumption
that every lineage response is source-driven.  If retained matter reacts to
the metric, that feedback must either be included in the independently solved
complete source map or in the retarded operator; it may not be counted in
both.

Using the SI action coefficient

\[
 C_{R,\rm SI}^{\rm eff}={c^3\over16\pi G_{\rm eff}},
 \qquad
 {8\pi G_{\rm eff}\over c^4}={1\over2cC_{R,\rm SI}^{\rm eff}}, \tag{PL07}
\]

the same relation is

\[
 \kappa_r^R-b_r^R-\rho_r^R
 ={\Phi_r^R\over2cC_{R,\rm SI}^{\rm eff}}.         \tag{PL08}
\]

## 3. Theorem PLGI-1 -- what is and is not identifiable now

The sealed Gravity Formation Theory fixes the following two facts separately:

1. the point-local lineage response has the form (PL03)--(PL04); and
2. the endpoint metric has the Einstein--Hilbert coefficient
   `C_R,SI^eff=c^3/(16 pi G_eff)>0`, calibrated by the common gravitational
   endpoint.

Those facts do **not** determine a further numerical relation between
`h_tr`, `h_sh`, and `G_eff`.  The missing join is not another representation
calculation.  It is the physical decomposition `(mathscr S_r,b_r,rho_r)` in
(PL06).

### Proof

First, at fixed `G_eff` and fixed endpoint operator, changing the map from a
lineage intervention into homogeneous/branch data changes `b_r` and hence
`h_r` without changing the Einstein--Hilbert coefficient.  Therefore `G`
does not determine `h`.

Second, even in a source-only model with `b_r=rho_r=0`, (PL06) is invariant
under

\[
 G_{\rm eff}\mapsto\lambda G_{\rm eff},
 \qquad
 \mathscr S_r\mapsto\lambda^{-1}\mathscr S_r,
 \qquad \lambda>0.                                \tag{PL09}
\]

Thus `h_r` determines only the product `G_eff mathscr S_r` until the complete
lineage-to-stress normalization is fixed independently.  An independent
rescaling of the lineage coordinate or of `ell_F` supplies the analogous
calibration degeneracy unless those quantities are physically frozen.

Finally, Einstein's diffeomorphism Ward identity constrains the complete
four-dimensional conserved source and propagates the constraints.  It does
not turn the spatial lineage coordinate `X` into such a source and does not
eliminate (PL09).  QED.

## 4. Corollary PLGI-2 -- the exact source-calibrated `G` relation

Suppose, for at least one sector and frequency in the declared point-local
response band, all of the following are established independently:

1. `ell_F` and the physical amplitude of `X=D delta L` are absolutely
   calibrated;
2. `mathscr S_r` gives the complete conserved four-stress response, including
   writer, controller, support, work, reservoir, EM, interaction, and boundary
   ownership exactly once;
3. the intervention supplies no unmodelled homogeneous, incoming,
   initial-data, boundary, or branch response, so `b_r` is known (zero in the
   source-only design);
4. `rho_r` is independently calculated or bounded below the inference scale;
   and
5. the gauge/constraint quotient and retarded Einstein solution map are fixed,
   and `Phi_r` is finite and nonzero.

Then `G_eff` is identifiable and obeys

\[
 \boxed{
 G_{\rm eff}
 ={c^4\over8\pi}
 {\ell_F^2h_r^R-b_r^R-\rho_r^R\over\Phi_r^R}.}    \tag{PL10}
\]

The right side must be the same positive, real, frequency-independent constant
for every independently calibrated sector and frequency after the declared
retarded phases and contact terms are handled.  Equivalently,

\[
 {\ell_F^2h_r^R-b_r^R-\rho_r^R\over\Phi_r^R}
 ={8\pi G_{\rm eff}\over c^4}.                    \tag{PL11}
\]

Equations (PL10)--(PL11) are a genuine calculation and a strong consistency
test, but not a parameter-free prediction: `Phi_r` contains the independently
measured lineage-to-complete-stress calibration.  Ratios between two sectors
cancel `G` and test the source/response construction,

\[
 {\kappa_r-b_r-\rho_r\over\kappa_s-b_s-\rho_s}
 ={\Phi_r\over\Phi_s},                            \tag{PL12}
\]

but cannot recover the common absolute scale if all `Phi` values share one
unknown normalization.

No universal `h_tr/h_sh` ratio follows from the Einstein--Hilbert action
alone.  The trace channel depends on the conserved four-source embedding,
constraints, gauge/observable choice, and background.  Once those data are
fixed, the Einstein retarded map calculates the ratio through (PL06); without
them, importing the spin-two pole numerator as a ratio of the two spatial
lineage form factors is a type error.

## 5. Corollary PLGI-3 -- the matched ancestry experiment cannot measure `G`

The RGRL-C / SPAG KEEP-versus-BREAK contrast is deliberately matched in the
complete ordinary stress and collateral channels.  For the ideal matched
contrast,

\[
 \delta T_{\mu\nu}^{\rm complete}=0,
 \qquad \mathscr S_r=0.                            \tag{PL13}
\]

Therefore (PL06) reduces to

\[
 \boxed{\kappa_r^R=b_r^R+\rho_r^R,}               \tag{PL14}
\]

and the ratio in (PL10) is undefined rather than a measurement of `G`.  A
nonzero, qualified matched response tests whether lineage supplies a
constitutive geometry/branch column after ordinary channels are removed.  It
does not say that the extra response is proportional to record count or gamma,
and it does not calibrate Newton's constant.

If complete stress, remainder, homogeneous/initial/boundary data, and every
physical state variable were all literally identical in a well-posed
single-solution endpoint, uniqueness would force the response to vanish.  The
scientific content of the ancestry test is precisely whether authenticated
lineage is a missing constitutive state/branch variable.  Calling it a stress
source in advance would assume away that test.

## 6. Exact remaining point-local calculation

The open calculation is now minimal and physical:

1. choose either the **matched ancestry lane**, which measures `b_r+rho_r` and
   tests RGRL but cannot determine `G`; or a separate **source-calibrated
   lane**, which fixes `mathscr S_r` and can apply (PL10);
2. construct the full conserved `mathscr S_r` for that particular apparatus or
   physical model, rather than equating it with gamma, record number, or an
   uncalibrated lineage label;
3. calculate the point-local retarded Einstein projection `Phi_r` on the
   declared gauge/boundary quotient; and
4. compare all available sector/frequency estimates using the constant-ratio
   condition (PL11).

No further group-theory machinery is needed at this scope.  A microscopic
derivation of the RGRL branch map would predict `b_r`; a calibrated
lineage-to-stress model would predict `Phi_r`.  At least one of those physical
maps must be supplied before the numerical `h` functions can be predicted.

## 7. Dependency and custody ledger

This theorem uses, without changing:

- `GRAVITY_FORMATION_THEORY_CLOSURE_V001.md`, SHA-256
  `63c37c85442fa96739591a1380a41ee29a9f6f66a6ca0afd5b3470d22fdce028`;
- `GRAVITY_RGRL_S4_LINEAGE_METRIC_RESPONSE_KERNEL_THEOREM_V001.md`, SHA-256
  `49e97e9cd3c9d8c75c65f3717156071bfcc0d88b3be3118aa442f74fb711f50d`;
- `GRAVITY_RGRL_IR_ENDPOINT_CLOSURE_THEOREM_V001.md`, SHA-256
  `c883c4c9f3816e453766846a1691ef27cb50d6ea7e5676bc52ed1928617f82bf`;
- `FINAL_GRAVITY_REAL_WORLD_THEOREM_V001.md`, SHA-256
  `1caabded24b861932b319ed715556a5d4123b2cff5ea3004676e12c4c76de155`;
- `GRAVITY_RGRL_SPAG_PROSPECTIVE_PROTOCOL_V001.md`, SHA-256
  `9495ca2b9edf3ebf1133e077d746e77b78ebda6e0fc061c178c80109506386b9`.

## 8. Disposition

`AT_POINT_LOCAL_SCOPE_O3_FIXES_H_E_EQUALS_H_T_BUT_NOT_TRACE_TO_SHEAR__THE_PHYSICAL_RESPONSE_EIGENVALUE_IS_ELL_F_SQUARED_H_R__EINSTEIN_HILBERT_SUPPLIES_THE_COMMON_METRIC_OPERATOR_AND_G_BUT_NOT_THE_LINEAGE_TO_COMPLETE_STRESS_OR_BRANCH_MAP__ONLY_THE_SOURCE_CALIBRATED_ZERO_HOMOGENEOUS_CONTROLLED_REMAINDER_JOIN_IDENTIFIES_G__THE_MATCHED_KEEP_BREAK_COLUMN_TESTS_CONSTITUTIVE_LINEAGE_ANCESTRY_AND_DOES_NOT_MEASURE_G`
