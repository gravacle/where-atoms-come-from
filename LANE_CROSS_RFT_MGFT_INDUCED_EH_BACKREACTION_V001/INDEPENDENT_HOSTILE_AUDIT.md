# Independent hostile audit

**Lane:** `LANE_CROSS_RFT_MGFT_INDUCED_EH_BACKREACTION_V001`

**Date:** 2026-08-26

**Final verdict:** `ACCEPT_AFTER_REPAIR`

## Defects found in the first draft

1. The proper-time integral extended through the infrared while the retained
   slow functional was also present, so the split could double-count modes.
2. `C_R^pre=0` was treated as sufficient for a strict induced-origin claim
   even though finite local terms can move under a matching convention.
3. The composite-metric chain rule omitted the force from explicit dependence
   of the collective action on the underlying variables.
4. Nonsurjectivity was treated as equivalent to a nontrivial adjoint kernel.
   That implication is false for a dense, nonclosed range in infinite
   dimensions.

## Repairs verified

- The determinant now occupies only the finite fast shell
  `s in [kappa_R^-2, mu^-2]`; the complementary infrared/nonlocal contribution
  is assigned separately.  The shell moments and light-mass expansion are
  correct.
- Strict induced origin is conditioned on a physically distinguished
  microscopic regulator/matching prescription in which the absent preterm is
  fixed by an exact condition or symmetry.  Otherwise the claim is limited to
  an induced contribution/renormalization and the matched total.
- The exact underlying stationarity equation is
  `(Dg)^* E_g + F_explicit = 0`.  A full metric equation therefore also needs
  the independent explicit force to vanish or to be absorbed into a proved
  enlarged-field residual.
- The physical tangent criterion is stated correctly as
  `ker((Dg)^*)=0`, equivalently dense range of `Dg` in the Hilbert setting.
  Surjectivity is only a sufficient condition; finite-mode claims use full row
  rank.
- The universality claim is limited to one common leading coefficient
  multiplying the complete variational stress.  The packet does not promote
  that statement to WEP, SEP, or absence of additional forces.

## Verification

The repaired verifier passes `37/37`.  No remaining material defect was found.

This verdict accepts the conditional induced Einstein--Hilbert and
back-reaction mechanism on an already-earned common metric.  It does not
accept a derivation of that metric, the sign or observed value of `G`, the
smallness of `Lambda`, or the required composite tangent realization from the
current F3 parent.
