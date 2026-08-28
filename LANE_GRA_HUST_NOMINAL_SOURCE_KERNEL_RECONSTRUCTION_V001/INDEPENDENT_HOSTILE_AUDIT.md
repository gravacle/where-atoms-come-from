# Independent hostile audit

**Lane:** `LANE_GRA_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_V001`

**Disposition:** `PASS_AFTER_MATERIAL_NARROWING_AND_REPAIR`

**Independent executable:** `hostile_audit_hust_nominal_source_kernels.py`

**Transcript:** `HOSTILE_AUDIT_TRANSCRIPT.txt`

## 1. Audit question

I treated the production reconstruction, its coordinate realization, and its
use of the word "exact" as untrusted. I independently rebuilt both Newtonian
functionals from `SOURCE_FIELDS.json`, checked the official figures and tables
visually, differentiated the ToS torque by an independent finite difference,
used different cubature orders, and constructed a second AAF coordinate set
with identical published pair distances and identical overall centroid.

The production numbers and physical signs survive. Three material claim
defects did not: numerical cubature had been called exact, the AAF coordinate
realization had been presented as if uniquely public, and the AAF-I coefficient
had been described too nearly as a complete campaign-temperature
reconstruction. Those defects are repaired below.

## 2. Official source semantics

The pinned official Supplementary Information formula pages and Tables 1--3
were rendered and read visually. The pinned official Extended Data Figure 3
and Tables 1, 2, and 4 were also inspected at their native resolution.

The source establishes the following:

- the ToS near position has the pendulum in line with the source pair and the
  far position is a 90-degree rotation;
- the AAF signal is the amplitude at twice the relative turntable frequency;
- Extended Data Figure 3 and Table 4 publish two ToS pair separations and four
  AAF horizontal/vertical pair separations;
- they do not publish the individual three-dimensional CMM coordinates of all
  sphere centres relative to the rotation axis;
- the official Extended Data Table 4 caption gives the upper AAF horizontal
  temperature coefficient as -1.9(1) micrometre/C, says the lower horizontal
  and vertical distances are constant within 2 micrometres over 4 C, and gives
  the ToS coefficient ceiling; and
- Supplementary Table 3 gives the AAF campaign temperatures and processed
  coefficients.

The official caption at
`https://www.nature.com/articles/s41586-018-0431-5/tables/5` was checked live
on 2026-08-27. Its numerical transcription is source-located in
`SOURCE_FIELDS.json`, but the caption page itself is not a separately
hash-pinned raw object. That custody ceiling is now explicit.

## 3. Independent reconstruction

The hostile executable imports no production module. At Gauss--Legendre order
20 with 768 independent azimuth samples it obtains:

- AAF-I: 6926.660438859112 kg m^-3;
- AAF-II: 6926.700007763448 kg m^-3; and
- AAF-III: 6926.700007763448 kg m^-3.

At Gauss--Legendre order 24 it independently obtains all seven ToS values,
from 24914.243522929188 to 25005.342718853015 kg m^-3, agreeing with the stored
order-20 values within 2.3e-10 kg m^-3. For the first ToS row the independently
resolved components are

\[
 C_{g,n}/I=18939.796309048641,
 \qquad C_{g,f}/I=-5975.387654091569,
\]

so near minus far is positive and has the stored factor and sign. A central
finite difference of the torque at both 0 and pi/2 independently reproduces
the analytic derivative

\[
 {d\over d\phi}{N\over d^3}={P\over d^3}-{3N^2\over d^5}.
\]

The sign is not inferred from the processed coefficient. With \(\phi\) the
source-line azimuth relative to the pendulum, a positive pendulum displacement
changes \(\phi\) by the opposite amount, so
\(C_g/G=-\partial(N_g/G)/\partial\theta_{\rm pend}
=\partial(N_g/G)/\partial\phi\). The positive near derivative, negative far
derivative, and near-minus-far convention therefore reproduce the physical
stiffness sign directly from the source geometry.

The AAF factor of two was separately checked against a unit sine at harmonic
two. The core inertia in every row equals
\(M(L^2+W^2)/12\). The displayed core-numerator/full-inertia values reproduce
that deliberately mixed arithmetic, but remain labelled collision diagnostics
and are not promoted to apparatus coefficients.

## 4. Shell-theorem domain

For a source centre on a circular orbit of horizontal radius \(R\), the exact
minimum horizontal distance to the rotating rectangular cuboid is obtained by
aligning the orbit with the cuboid's farthest horizontal corner. Combining that
with the independent vertical gap gives an analytic minimum over all azimuths,
not a sampled one. The resulting minimum surface clearances are

- AAF: 0.069871243962 m; and
- ToS: 0.004170645924 m.

Both are positive. Therefore Newton's shell theorem lawfully replaces each
*declared homogeneous spherical* source by its centre point mass throughout
the cuboid. This exact reduction does not apply to the real nonsphericity,
density multipoles, mounts, or unreported offsets; those remain in the missing
full-apparatus remainder.

## 5. Numerical exactness defect and repair

The order comparisons are extraordinarily stable: the independent AAF
order-20 minus order-16 differences are about 1.55e-11 kg m^-3, and the ToS
order-24 minus order-20 differences are about -2.3e-10 kg m^-3. They establish
convergence at a scale vastly below the public input uncertainties.

They are not, however, a certified quadrature error bound. No interval
arithmetic, analytic remainder theorem, or other rigorous cubature certificate
was supplied. The theorem and result now attach exactness only to the declared
Newtonian functional, the ToS derivative identity, and the conditional
shell-theorem reduction. The reported coefficients are correctly described as
converged numerical evaluations.

## 6. AAF coordinate non-uniqueness defect and repair

Four pair differences plus an overall zero centroid do not fix eight horizontal
and vertical coordinates. Two shear degrees of freedom remain. The hostile
audit constructs a second coordinate set by shifting the upper and lower rows
oppositely in the horizontal direction and the left and right columns
oppositely in the vertical direction. With 30-micrometre shears it preserves
all four published pair distances and the overall centroid exactly, yet changes
the nominal coefficient from 6926.700007763433 to 6926.702624350222 kg m^-3.

This is a direct identifiability collision. It does not refute the natural
pairwise-centred nominal realization; it proves that realization is an explicit
conditional premise rather than a unique consequence of the public pair data.
The theorem, source fields, result, and claim ceiling now say so.

## 7. Temperature transport and sensitivities

The AAF-I calculation correctly transports only \(S_{7,9}\):

\[
0.3422874+[-1.9\times10^{-6}](22.8-23.7)
=0.34228911\ {\rm m}.
\]

All other AAF geometry remains at its published 23.7 C reference. ToS-I uses
pendulum/source dimensions at the published 20.2 C reference and run-specific
separations at 20.1 or 20.3 C; ToS-II is at 21.5 C. The lane now calls these
partial public transports and does not imply a complete campaign-temperature
mass map.

The hostile executable independently redoes the symmetric one-standard-
uncertainty perturbations. For AAF-II it obtains 0.156575 ppm from pendulum
dimensions, 0.308787 ppm from source masses, 8.974932 ppm from horizontal
distances, and 5.792504 ppm from vertical distances. It also reproduces
representative published ToS-I fibres-1/2 and ToS-II dimension, mass, and
distance classes. As recorded in the theorem, independent scalar propagation
gives 0.727 ppm for the ToS-I source masses in every fibre, whereas the campaign
table assigns 0.55 ppm to fibre 3; that unresolved campaign/covariance
ownership is retained rather than called a reconstruction. These checks
support the local sensitivity ledger, not an exact coverage set: standard
uncertainties do not define compact support and the covariance is incomplete.

## 8. Comparator quarantine and remainder semantics

No accepted value of \(G\), measured response, or processed coefficient occurs
in `SOURCE_FIELDS.json` or in the kernel code. The production executable builds
AAF, then ToS, and only then reads `PUBLISHED_COMPARATORS.json`; mutating that
file cannot change either independent functional. This proves code-level input
quarantine. It cannot prove a historical mental state, so the theorem no longer
overstates what executable ordering establishes.

All processed-minus-nominal discrepancies have the reproduced negative sign.
They locate scalar differences between the conditional nominal calculation and
the authors' processed coefficients. They do not identify a unique spatial
clamp/coating/CMM/density remainder. The summed one-standard-uncertainty bands
remain explicitly non-coverage diagnostics. Equation (H07) is now a formal
image for an explicitly chosen input domain, not an exact public identified set
that the quoted standard uncertainties fail to define.

## 9. Scientific ceiling

After repair, the lane lawfully supplies a strong no-lab partial source model:
three conditional AAF harmonic coefficients, seven conditional ToS stiffness
coefficients, exact shell-theorem domain clearance, and independently checked
local sensitivities. It still does not provide the full mass/stress map,
individual CMM coordinates, transfer functions, covariance, row-level phase
data, signed remainders, or prospective nulls required by GC16. It does not
derive or re-estimate \(G\), confirm RGRL or Gravity Formation Theory, or close
conserved-stress ownership.

With those boundaries explicit, the hostile audit passes 94/94 checks.
