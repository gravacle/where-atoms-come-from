# Self-audit

## Adversarial checks performed

- Rebuilt the physical `D_a` and normalized root dyads with their distinct
  `1/3` and `1/8` normalizations.
- Reparsed all 720 paths and verified all 3600 proper-prefix gaps, the five
  gap classes, and `63/8`.
- Subtracted the source derivative of the initial ice energy from every
  forward virtual gap and repeated the calculation from the final endpoint.
- Used the Hermitian forward/reverse average; a one-sided Bloch entry was not
  promoted to an operator.
- Exhausted all four missing-label orientations and all `2^6` allowed local
  external-link assignments.
- Kept the unchanged FS one-edge numerator source and the FU DPAR denominator
  source in the same pre-Feshbach calculation.
- Exposed `S10 / FV-PURE` as an additional premise: FU `S1`--`S9` permit
  residual source kernels and do not alone make the ideal-Coulomb term the
  complete nonidentity derivative.
- Constructed four actual global `G_5` ice transitions and two actual
  direction-covering diagonal differences, so coefficient rank was not
  mistaken for operator rank.
- Removed identities only through functionals that annihilate them; no scalar
  shift was counted.
- Checked the exact cancellation slopes `0`, `2/5`, and `3/5`, specifically
  confirming that Coulomb `-1/2` is not one.
- Classified lower folds as diagonal by endpoint topology rather than
  silently dropping them.

## Remaining ceilings

The exact conclusion is formal source rank through H8.  A finite nonzero `h`
claim outside an unspecified small-coupling neighborhood would require a
convergence bound or exclusion of a tuned algebraic zero after higher orders.
No such numerical bound is claimed.

Off-shell operator rank is necessary but not sufficient for retarded or CTP
rank.  Ward constraints, state choice, spectral residues, tensor poles,
continuum gluing, common-metric coupling, gravity, and `G` remain outside this
lane.

The exact coefficients and rank are conditional on `FV-PURE`.  A residual
pair, cross-node, boundary, controller, or other source derivative with a
different normalized strain law is not covered and can change the H6 tensor
or its special slopes.
