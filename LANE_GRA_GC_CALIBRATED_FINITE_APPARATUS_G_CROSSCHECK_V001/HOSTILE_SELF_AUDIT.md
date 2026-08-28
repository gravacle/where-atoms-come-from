# Hostile self-audit: calibrated finite-apparatus (G) cross-check V001

**Audit target:** `GRA-GC-CFAGC-V001` theorem, protocol, implementation, and
synthetic result

**Date:** 2026-08-27

**Disposition:** `PASS_WITH_EXPLICIT_SCIENTIFIC_CEILINGS`

## 1. Exact mathematics checked

1. From
   (U/G=-\sum_{ij}M_im_j/|\mathbf R_i-R_z(\theta)\mathbf r_j|),
   direct differentiation gives (K=-G^{-1}\partial_\theta U) with the sign in
   (GC07).  Independent centered finite differences agree.
2. Differentiating (K) with respect to the source tangent and detector torsion
   angle gives (GC08) and (GC09), including the opposite signs in their radial
   derivatives.  Independent finite differences agree.
3. Newton's shell theorem supplies the center-mass reduction only for mutually
   nonoverlapping spherical bodies; the implementation enforces this domain.
4. With the declared (e^{+i\omega t}) convention, damping contributes
   (+ib\omega) and a positive readout delay contributes
   (e^{-i\omega\tau}).
5. Direct inversion of the full calibrated two-mode operator equals the Schur
   complement to numerical roundoff.  The implicit gravitational stiffness
   appears only in the operator; it is not duplicated in the source column.
6. The same-source/remainder/data zero follows from invertibility, independent
   of any adopted off-shell metric tangent.
7. Cross multiplication proves global one-row injectivity in (p=Gs) on the
   nonsingular domain.  The transformation ((G,s)\mapsto(Gq,s/q)) leaves both
   source and implicit stiffness unchanged, proving exact global-scale
   nonidentifiability.  Positive interval division gives (GC22).
8. The static determinant is affine in (p).  Endpoint positivity, positive
   (\kappa_x), and the declared positive damping/inertia conditions establish
   the two-mode stable scan domain used by the executable.

## 2. Source and finite-spatial-profile scope

The equal antipodal source pair has fixed total mass, zero mass-dipole
trajectory tangent, and a nonzero symmetric trace-free mass-quadrupole tangent.
Therefore only the Fourier transform of its Newtonian mass-density tangent is
proved to start at (O(|\mathbf k|^2)).  No claim is made that every component
of the complete conserved (\delta T^{\mu\nu}) starts at that order.  Momentum,
drive, and support stresses remain necessary for the relativistic Ward
identity, and V002 GI25 remains a real-apparatus condition.

The extended-source derivative is exact for a fixed-weight Lagrangian mass
measure transported by the declared embedding.  Density creation, loss, or
unmodeled deformation is outside that exact formula and must enter a larger
source model or the bounded remainder.

## 3. Units and operator type

(K), (a), and (k_g) have units of mass squared per length when (u) and
(\theta) are dimensionless.  Multiplication by (G) gives torque or generalized
stiffness as appropriate.  The auxiliary coordinate is explicitly normalized
as a dimensionless generalized mode so the reciprocal coupling in (GC12) is
well typed.  A differently normalized physical coordinate requires two
calibrated conversion factors and their product in the Schur complement.

## 4. Derivative and uncertainty ownership

The audit found no duplicate derivative:

- source mass/trajectory enters once through (pa);
- gravitational feedback enters once through (-pk_g);
- the auxiliary mode enters once through (-\lambda^2/d_x);
- physical torques enter before the inverse;
- homogeneous data enter as solution data; and
- readout bias and detector noise enter after the physical response.

The observation covariance and independent gain/delay calibration covariance
are combined once.  Their Jacobian is frozen at a declared absolute reference,
not at the generator truth or fitted value.  This is only a first-order
calibration model; material nonlinearity requires explicit nuisance inference
and the likelihood log determinant when covariance depends on (p).

## 5. Synthetic noncircularity

The first implementation incorrectly set the product scan to fixed percentages
of the hidden synthetic truth and checked stability at another truth-derived
point.  Hostile review rejected that construction.  The repaired implementation
declares an absolute scan ([10^{-12},1.5\times10^{-10}]) and an independent
source-scale calibration interval before the arbitrary generator constants,
then checks stability at both scan endpoints.  The generator truth is now used
only to create observations and, after inference, to score coverage and error.
The accepted empirical value of (G) is nowhere used.

The first covariance implementation also evaluated its calibration Jacobian at
the generator truth.  It was replaced by a prospectively fixed absolute
reference (p_{\rm cov}=7.0\times10^{-11}).

## 6. Statistical ceiling

The nominal profile-grid interval uses a one-parameter (\Delta\chi^2)
threshold.  It is an implementation diagnostic, not a finite-sample coverage
proof.  The held-out score is likewise a deterministic diagnostic because the
single synthetic realization and its correlated covariance do not constitute
an independent prospective experiment.  For a real dataset, bounded
remainders and calibration sets must be propagated into (GC26).

## 7. Independent hostile finding and repair

An independent hostile reviewer found and rechecked three seal-relevant issues:

1. truth-dependent scan/stability bounds -- repaired with an absolute frozen
   domain;
2. overbroad inference from density multipoles to complete conserved
   (T^{\mu\nu}) -- narrowed to the Newtonian density tangent with GI25 open;
3. stale theorem counts and result intervals after repair -- synchronized to
   the final 15/15 execution.

The reviewer then returned `FINAL PASS`: kernel signs,
conservation/approximation ceiling, finite-density spatial-momentum scope,
Schur complement, identifiability, covariance ownership, fixed scan, and
synthetic generator/estimator separation are clean.

## 8. Scientific ceilings retained

This package does not:

- measure (G) from real data;
- derive the numerical value or microscopic origin of (G);
- derive a lineage-dependent complete stress-energy column;
- prove GI21's lineage-to-source compatibility join;
- turn the off-shell RGRL-C ancestry tangent into an on-shell response;
- establish the complete real-apparatus Ward identity or GI25; or
- empirically confirm RGRL, EIR, or GFT.

Within those ceilings, it closes the requested no-lab step: a calibrated finite
ordinary source, exact declared Newtonian geometry, full retained two-mode
dressed operator, transfer, nuisance/covariance model, and source-scale-aware
identified set form one non-double-counted executable forward calculation.
