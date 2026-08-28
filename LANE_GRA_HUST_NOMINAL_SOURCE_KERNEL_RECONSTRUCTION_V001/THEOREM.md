# HUST nominal source-kernel reconstruction and public identifiability ceiling

**Theorem ID:** `GRA-HUST-NSKR-V001`

**Date:** 2026-08-27

**Claim class:** exact conditional homogeneous finite-source functional with
converged numerical evaluation; public calibration sensitivity reconstruction;
post-calculation processed-coefficient comparison; finite public-field
identifiability ceiling.

**Status:**
`AAF_CONDITIONAL_HOMOGENEOUS_SOURCE_KERNEL_RECONSTRUCTED__TOS_CONDITIONAL_HOMOGENEOUS_STIFFNESS_KERNEL_RECONSTRUCTED__PUBLIC_GEOMETRY_SENSITIVITIES_RECOVERED__FULL_CENTRAL_REMAINDER_NOT_IDENTIFIED_FROM_UNPROCESSED_PUBLIC_FIELDS__GC16_NOT_CLOSED__NO_NEW_G`

## 1. Result

The official HUST-2018 public dimensions, masses, sphere-centre separations,
campaign temperatures, and an explicit nominal pairwise-centred placement
premise are sufficient to calculate two useful conditional finite-apparatus
objects without using a measured value of \(G\), a gravity response, or a
published processed source coefficient:

1. the AAF \(m=2\) nominal source-forcing coefficient; and
2. the seven ToS nominal near-minus-far stiffness coefficients
   \(\Delta C_g/I\).

The Newtonian functionals are exact for the declared **conditional homogeneous
geometry**: a uniform centred cuboid and mutually disjoint homogeneous source
spheres at the assumed nominal centres. Newton's shell theorem makes each
sphere exactly a point mass for every cuboid element in this domain; an
analytic minimization over source azimuth gives positive surface clearances of
69.87 mm for AAF and 4.17 mm for ToS. The displayed coefficients, however, are
converged Gauss--Legendre/DFT evaluations of those exact integrals, not
certified exact numerical values. The real spheres' reported nonsphericity,
density inhomogeneity, coating, clamps, ferrules, and unpublished
three-dimensional CMM offsets remain outside the conditional domain.

Extended Data Figure 3 and Table 4 publish four AAF pair separations, not the
individual three-dimensional CMM coordinates relative to the rotation axis.
Even after fixing the overall centroid, those four differences leave two shear
degrees of freedom. Equation (H01) therefore adds the explicit premise that
each horizontal pair and each vertical column is separately centred. It is the
natural nominal symmetric realization, but it is not uniquely entailed by the
public pair distances.

This is a source-kernel reconstruction, not a new measurement of \(G\).

## 2. AAF source coefficient reconstructed first

Let the uniform pendulum core occupy the box \(B\), with

\[
 I_0={M(L^2+W^2)\over 12}.
\]

The nominal source centres are realized directly from the four reported
horizontal and vertical coordinate separations:

\[
\begin{aligned}
R_7&=(-S_{7,9}/2,0,+S_{7,10}/2),&
R_9&=(+S_{7,9}/2,0,+S_{9,12}/2),\\
R_{10}&=(-S_{10,12}/2,0,-S_{7,10}/2),&
R_{12}&=(+S_{10,12}/2,0,-S_{9,12}/2).
\end{aligned}                                                     \tag{H01}
\]

This preserves the small upper/lower and left/right mismatches instead of
silently replacing them by a perfect rectangle. For relative source azimuth
\(\phi\), the angular acceleration per \(G\) is

\[
 a_0(\phi)={1\over I_0}\sum_s M_s
 \int_B\rho_0\,
 {xR_{s,y}(\phi)-yR_{s,x}(\phi)\over
 |R_s(\phi)-r|^3}\,d^3r,                                      \tag{H02}
\]

and the independently calculated source coefficient is

\[
 K^{(0)}_{\rm AAF}
 =2\left|{1\over2\pi}\int_0^{2\pi}
 a_0(\phi)e^{-2i\phi}\,d\phi\right|.                          \tag{H03}
\]

Tensor Gauss--Legendre orders 12 and 16 and 256 versus 512 azimuth samples
agree below \(4\times10^{-12}\ {\rm kg\,m^{-3}}\). This is a strong empirical
convergence check, not a rigorous quadrature-error certificate. The calculated
values are:

| campaign | temperature | \(K^{(0)}_{\rm AAF}\) (kg m\(^{-3}\)) | public-input RSS \(u\) | processed comparator, opened afterward | processed minus nominal |
|---|---:|---:|---:|---:|---:|
| AAF-I | 22.8 °C | 6926.660438859 | 0.074057725 | 6926.352(74) | -0.308438859 (-44.531 ppm) |
| AAF-II | 23.7 °C | 6926.700007763 | 0.074029014 | 6926.334(75) | -0.366007763 (-52.843 ppm) |
| AAF-III | 23.7 °C | 6926.700007763 | 0.074029014 | 6926.415(74) | -0.285007763 (-41.148 ppm) |

AAF-I uses only the published upper-pair temperature coefficient,
\(-1.9(1)\ \mu{\rm m}\,{^\circ\rm C}^{-1}\), to transport \(S_{7,9}\) from
23.7 °C to the campaign's published 22.8 °C mean. Other AAF-I geometry fields
remain at their published 23.7 °C reference values because the release does
not provide their full campaign transport. Thus AAF-I is a partial public
temperature transport, not a complete campaign-temperature mass map. The
processed values in the fifth column come from Supplementary Table 3. The
executable quarantines them in a separate file and reads them only after
equations (H01)--(H03), convergence, and input sensitivities are computed;
this proves code-level input independence, not unknowable historical cognitive
independence.

The independently propagated AAF relative sensitivities are 0.157 ppm from
pendulum dimensions, 0.309 ppm from source masses, 8.975--8.980 ppm from the
two horizontal separations, and 5.792 ppm from the two vertical separations.
They reproduce the corresponding official main-Table-1 classes of 0.16,
0.31--0.32, 8.98, and 5.79 ppm. Thus the 41--53 ppm
processed-minus-nominal scalar discrepancy is not quadrature error and is not
removed by the local public-input sensitivity scale. It locates the net
difference between this declared nominal core and the authors' processed
coefficient; it does not identify a unique physical mass-map remainder.

## 3. ToS stiffness coefficients reconstructed second

For each public source separation \(S\), take the two nominal source centres
at \(\pm S/2\) and rotate their line by \(\phi\). Define

\[
 N=xR_y-yR_x,\qquad P=xR_x+yR_y,\qquad d=|R-r|.
\]

Here \(\phi\) is the source-line azimuth relative to the pendulum. A positive
pendulum deflection decreases that relative azimuth, so the conventional
positive gravitational torsion constant satisfies
\(C_g/G=-\partial(N_g/G)/\partial\theta_{\rm pend}
=\partial(N_g/G)/\partial\phi\). This fixes the stiffness sign before the
near-minus-far subtraction.

The exact azimuth derivative of the torque kernel is

\[
 {d\over d\phi}{N\over d^3}={P\over d^3}-{3N^2\over d^5}.     \tag{H04}
\]

Consequently the nominal ToS coefficient is evaluated without a finite-angle
derivative:

\[
 K^{(0)}_{\rm ToS}={\Delta C_g\over I_0}
 =a_0'(0)-a_0'(\pi/2).                                        \tag{H05}
\]

Gauss--Legendre orders 16 and 20 agree below
\(2\times10^{-10}\ {\rm kg\,m^{-3}}\). Again, this establishes numerical
convergence at the shown scale but is not a certified error bound. The seven
results are:

| run | \(K^{(0)}_{\rm ToS}\) (kg m\(^{-3}\)) | public-input RSS \(u\) | processed comparator, opened afterward | processed minus nominal |
|---|---:|---:|---:|---:|
| TOS-I-F1 first | 24915.183963140 | 0.222910472 | 24912.86(23) | -2.323963140 (-93.284 ppm) |
| TOS-I-F1 repeat | 24914.429047781 | 0.222900386 | 24912.12(23) | -2.309047781 (-92.688 ppm) |
| TOS-I-F2 | 24914.429047781 | 0.222900386 | 24912.15(23) | -2.279047781 (-91.483 ppm) |
| TOS-I-F3 first | 24914.243522929 | 0.216659663 | 24911.70(22) | -2.543522929 (-102.102 ppm) |
| TOS-I-F3 repeat | 24914.294701994 | 0.216660328 | 24911.72(22) | -2.574701994 (-103.353 ppm) |
| TOS-II-F4 first | 25005.342718853 | 0.248129925 | 25003.05(25) | -2.292718853 (-91.698 ppm) |
| TOS-II-F4 repeat | 25005.259029544 | 0.241955574 | 25002.95(25) | -2.309029544 (-92.350 ppm) |

The direct sensitivity calculation recovers the official dimension classes
(1.815 versus 1.82 ppm for ToS-I; 2.726 versus 2.73 ppm for ToS-II), the
horizontal-distance classes (8.730, 8.474, 9.526 and 9.268 ppm), and the mass
class for ToS-II (0.545 versus 0.55 ppm). The direct independent-scalar
ToS-I mass propagation is 0.727 ppm; the official campaign table assigns
0.73 ppm to fibres 1/2 and 0.55 ppm to fibre 3, so the latter difference is
retained as campaign/covariance ownership rather than overwritten.

The TOS-I pendulum and source dimensions are published at the 20.2 °C
reference while its run-specific separations are tabulated at 20.1 or 20.3 °C;
TOS-II is at 21.5 °C throughout. No unpublished thermal coefficients are used
to transport the remaining TOS-I dimensions. The seven values are therefore
conditional evaluations of the exact listed public fields, not complete
thermal reconstructions of each campaign.

## 4. Calibration envelope and identified set

For public inputs \(x_i\) with quoted standard uncertainties \(u_i\), the
executable evaluates the symmetric sensitivities

\[
 s_i={|K(x_i+u_i)-K(x_i-u_i)|\over2},\qquad
 u_{K,\rm lin}=\sqrt{\sum_i s_i^2}.                            \tag{H06}
\]

It also reports \(\sum_i s_i\) as a local linearized axis-box diagnostic.
Neither number is promoted to an exact coverage theorem: a quoted standard
uncertainty does not define compact support, and the public release does not
supply the full covariance.

For any explicitly chosen input domain \({\cal X}\), the formal conditional
image is

\[
 {\cal I}_{\rm hom}({\cal X})=\{K_{\rm hom}(X;P):X\in{\cal X}\}, \tag{H07}
\]

where \(P\) includes the nominal placement premise. The public standard
uncertainties do not themselves define a compact-support domain
\({\cal X}\), so this lane does not promote (H07) to an exact public coverage
set. What is numerically supplied is the nominal point plus local symmetric
sensitivities.

The full-apparatus object has the form

\[
 K_{\rm full}=K_{\rm hom}(X)+r_{\rm clamp}+r_{\rm coat}
 +r_{\rho}+r_{\rm CMM}+r_{\rm shape}+\cdots .                 \tag{H08}
\]

The pinned unprocessed public fields do not locate or numerically bound the
central remainder in (H08). The official one-standard-uncertainty classes
quantify uncertainty around the experiment's processed model; they do not
publish the missing central spatial maps.

This non-identifiability is structural. Extended Data Table 1 specifies the
bare cuboid, whereas Supplementary Table 1 gives the inertia of the full
pendulum assembly. Two missing clamp/coating distributions can have the same
total \(I_z\) yet different azimuthal multipoles: rotating a non-axisymmetric
missing distribution about the fibre leaves \(I_z\) invariant while changing
its \(q_{lm}\) phase and therefore its source torque. Hence a core numerator
divided by the full-assembly inertia is not the full coefficient. The
executable displays that forbidden mixed normalization only as a collision
diagnostic, never as a result.

Supplementary Tables 2 and 3 do publish processed scalar coefficients. Once
read by the executable after the nominal calculation, they locate a scalar
processed-minus-nominal discrepancy and permit the comparison tables above.
They do not turn that scalar back into an independently reconstructed mass
map. The reported remainder bands in
`RESULT.json` conservatively add the two quoted one-standard-uncertainty half
widths because their covariance is unavailable; those bands are neither
confidence intervals nor independent source reconstructions.

## 5. What this supplies to GC16—and what it does not

This lane supplies the strongest public-data partial source kernel now
available for the calibrated finite-apparatus model:

- AAF supplies a nominal finite-source \(m=2\) forcing coefficient in the
  homogeneous domain;
- ToS supplies nominal finite-source near/far stiffness entries; and
- both supply a checked geometry sensitivity ledger and a quantified
  post-calculation central-remainder diagnostic.

It does **not** close the real-data GC16 map. Still absent are the full
detector/source/support/drive mass and stress measures, row-level source
trajectory and phase, complete damped and coupled mechanical transfer,
readout transfer, covariance, signed physical remainders, and prospective
null/held-out rows. In particular, a prescribed moving source mass is not the
complete conserved \(T^{\mu\nu}\) of its drive and supports.

The sharp lawful successor is therefore `HUST-FULL-MASS-MAP-GC16`: obtain or
release the experiment's finite-element/CMM mass-coordinate-density files,
matched full-assembly inertia and numerator, transfer calibration, row-level
phase-referenced observations, covariance, and signed remainder ledger; then
evaluate the already frozen GC16 architecture without fitting a geometry
coefficient to gravity data. Until those fields exist publicly, the nominal
kernel plus the proved ceiling is the maximal no-lab result.

## 6. Claim ceiling

This theorem does not claim a new or re-estimated \(G\), confirmation of RGRL
or Gravity Formation Theory, complete conserved-stress ownership, equality of
the homogeneous model and the HUST apparatus, or a thermodynamic/general
gravity result. It does not claim that the AAF placement is uniquely fixed by
the public pair separations, a complete campaign thermal transport, or a
certified exact numerical quadrature value. It does not use the published
processed coefficients as calculation inputs.
