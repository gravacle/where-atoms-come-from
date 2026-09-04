# GL6CO Author Self-Audit

## Exact checks

- Constructed all 24 tetrahedral port actions from the four coordinate
  vectors and verified that each is orthogonal and a signed permutation.
- Verified the exact cycle decomposition `A1+T2` under all 24 actions.
- Used a character calculation to prove that the quadratic invariant space
  has dimension five and the constant symmetric invariant space dimension
  two.
- Checked covariance of every displayed quadratic basis tensor under all
  24 group elements with rational polynomial arithmetic.
- Re-derived the four GL6CL tensor-writer rows from the centered incidence
  geometry and composed each invariant cycle block explicitly.
- Fixed the response normalization to the bare cycle susceptibility and
  applied the GL6CL writer scale exactly once; this prevents a hidden
  `lambda_T^2`/`mu^2` double count.
- Kept the `GL6CL` coordinate `j_+=(j_P+j_C)/2` distinct from the normalized
  common coordinate `(j_P+j_C)/sqrt(2)`.  The latter gives a `mu^2/2` cycle
  Hessian when it is compared with the normalized `GL6BV` contact.
- Verified the three coefficient formulas in the pullback and that cycle
  `A1-A1/A1-T2` quadratic terms first enter above order `k^2`.
- Derived the necessary-and-sufficient `SO(3)` symmetric-tensor restriction
  plane.  This prevents the earlier error of demanding ordinary-vector
  isotropy from `T2` alone.
- Supplied an explicit positive analytic coefficient witness satisfying the
  one matching relation.
- Reconstructed all four defect-frame contact projectors and proved the
  order-`k^2` contact formula without importing the GL6BV compact result.
- Reconstructed the full inherited pair-to-symmetric-tensor solder and a
  generic-momentum six-by-six linearized-Einstein reference matrix.  Its
  rank-three quotient and nonzero `A1/E2-T2` blocks make the `T2`-only claim
  ceiling an executable kill test.
- Checked the relative strong-lock scaling of the spectral cycle response
  and the separately typed contact before stating the finite-ratio/collective
  implication as power counting only.

## Scope challenges

- The stationary bulk symbol is an assumption being classified, not a result
  derived from the finite-component Perron--Frobenius theorem.
- The repaired GL6CL tensor rows are used; its unclassified arbitrary-profile
  `A1/E` direct-history components are not used.
- The stronger `D-O` comparison is labeled a held-out algebraic reference,
  not a Ricci identification.
- Contact and cycle response are not added as facts.  Their sum is written
  only as the equation a future same-state source-first calculation would
  have to satisfy.
- The response has not been inverted to a 1PI operator.  Full `E2-T2` and
  `E2-E2` blocks, diagonal order-six vertices, phase selection, continuum
  calibration, causal continuation, gravity, and `G` remain open.

This self-audit does not replace an independent hostile audit.
