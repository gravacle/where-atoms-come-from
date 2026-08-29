# Independent hostile audit: PMICS

**Date:** 2026-08-29  
**Target:** repaired PMICS bytes pinned by `TARGET_CUSTODY.sha256`  
**Verdict:** **PASS, at the exact flat-background/principal-symbol ceiling**

## Scope and independence

This audit reconstructed the q4 observable-memory map, Fourier Riemann,
Ricci, and Einstein symbols, gauge kernel, trace/tidal census, exact FZ
direction witness, and dependency custody without importing the builder
verifier or any of its helpers.  The independent verifier uses Walsh-parity
identities, a Leibniz determinant and rank-by-minors calculation, and a
directional four-linear curvature contraction rather than the builder's row
reduction and component-loop implementation.

Two material issues were found while the target was mutable and were repaired
before the bytes audited here were frozen:

1. The pure-gradient zero-kernel statement is exact for a flat-reference
   linearization, or for the two-derivative principal symbol in a locally
   frozen frame.  It is not the literal kernel of the full linearized
   curvature operator on a curved background, where a diffeomorphism changes
   background-curvature components by a Lie derivative.  The repaired target
   now states this ceiling explicitly.
2. EW's natural parameter `J_ab` is a preparation/control coordinate, not an
   observable retained record.  The repaired target uses the observable pair
   expectations `C_ab=<Y_ab>` throughout.  It proves the exact statewise
   identity

   `XX^T = I_V + sum_ab CARRIER[Y_ab](v_a odot v_b)`

   before expectation, hence

   `F_theta(C) = I_V + sum_ab C_ab(v_a odot v_b)`.

   The equality of the old and new tangent matrices is justified only because
   `D_J C=I_6` at the uniform point.  The target now explicitly says that
   RGRL's older word “controls” is instantiated here by observable `C`, not by
   preparation coordinate `J`.

## Independent reconstruction

- All five target files and all ten declared parent dependencies replayed at
  their pinned SHA-256 digests.  The dependency ledger contains exactly those
  ten entries, with no omitted or extra path.
- The tetrahedral contrast vectors independently give the six columns
  `v_a odot v_b`.  Their six-by-six determinant in coordinate order
  `(xx,yy,zz,xy,xz,yz)` is exactly `-1/2`, proving the EW map is an
  isomorphism.  Walsh parity reproduces every column as the derivative of the
  localization Fisher tensor.  Direct evaluation on all sixteen carrier
  states verifies the stronger statewise observable-memory identity above.
- The uniform pair-score covariance is exactly `I_6`.  Therefore
  `D_J C=I_6`, the `C` chart is locally valid, and the local complete-query
  expansion `-log gamma_Q=(1/4)||delta C||^2+O(||delta C||^3)` has the stated
  normalization.
- In a canonical orthonormal frame with `k=q e_3`, the independently rebuilt
  Riemann symbol is

  `R_(a n b n)=(q^2/2) h_ab`,

  and selects exactly `(h_11,h_12,h_22)`.  Its rank is three.  The sign agrees
  with the frozen derivative convention `partial -> i k`.
- The zero transverse block has dimension three and is exactly
  `h=k odot xi`.  Since the q4 map is invertible, its pullback exhausts the
  three-dimensional pair-memory kernel.  This is the flat/principal-symbol
  metric-curvature kernel; it does not say that distinct observable memory
  states are themselves physically gauge-identical in every other sector.
- The scalar symbol is
  `R^(3)=|k|^2 tr(P h P)` with the declared positive sign.  The transverse
  symmetric block splits into one trace and two trace-free directions.  Two
  explicit independent trace-free tensors have zero scalar curvature but
  nonzero Ricci/Einstein curvature.
- The independently contracted three-dimensional Einstein symbol is
  transverse, has rank three, and on the transverse plane equals

  `G_ab=(|k|^2/2)(H_ab-P_ab tr H)`.

  Two-dimensional trace reversal is invertible, so it neither adds nor loses
  a curvature direction.
- At `r=(7,15,-17)`, the orthogonal frame
  `u=(15,-7,0)`, `w=(-119,-255,-274)` reproduces the complete displayed
  three-by-six matrix.  Its `(12,13,23)` minor is exactly
  `-173782321152`; the sign and nonorthonormal-frame scaling are correct.
- At zero momentum all six first-order metric deformations have zero Fourier
  curvature, as required.

The independent verifier reports `112/112 PASS`.

## Hostile physical-claim review

No remaining material sign, rank, kernel, basis-normalization, dependency, or
operator-typing defect was found in the repaired target.  The PASS is bounded
as follows:

- PMICS proves a finite-dimensional kinematic curvature-capacity theorem.  It
  does not prove a curvature equation of motion or a gravity response.
- The physical interpretation uses adopted RGRL-B--C.  It is an axiomatic
  consequence inside that working theory, not an empirical derivation of
  RGRL and not a microscopic F3 realization.
- The full curved-background operator requires the omitted lower-derivative
  background-curvature terms and their Lie derivatives.
- `C` is observable pair memory; `J` remains a preparation coordinate.
- Scalar gamma certifies local distinguishability.  It is not curvature,
  stress, stiffness, `G`, or a force magnitude.
- The arbitrary conversion scale `ell_F`, physical gluing and refinement,
  the same-parent `C/J/j/Q` join, temporal/current/contact completion, the
  stress Ward identity, and positive effective-action stiffness remain open.
- The FZ rank-two TT source response and the PMICS metric perturbation remain
  different tensor types.  Their matching representation dimensions do not
  authorize substitution.
- No particle pole, graviton, Newtonian limit, Einstein equation, or gravity
  theorem follows from PMICS alone.

## Audited ceiling

The repaired packet proves that the six observable q4 pair-memory expectation
directions span every infinitesimal deformation of the candidate localization
metric, and that at each nonzero momentum its exact flat-reference or
principal two-derivative intrinsic-curvature symbol has rank three modulo the
three pure-gradient metric directions.  The quotient is one scalar-curvature
plus two trace-free tidal-curvature directions.  Its physical interpretation
is conditional on adopted RGRL; its F3 dynamics and gravity closure remain
open.

