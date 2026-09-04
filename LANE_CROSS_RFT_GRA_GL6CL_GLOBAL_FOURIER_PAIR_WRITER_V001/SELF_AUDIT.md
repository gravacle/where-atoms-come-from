# GL6CL Author Self-Audit

## Checks performed

- Derived all four inherited-coordinate hexagon centers and all twelve parent/child
  offsets from the infinite-parent node positions.
- Verified every parent offset is the negative of the same-pair child offset
  and every centered radius squared is `11/4`.
- Constructed `A1`, `E`, and `T2` projectors with exact rational arithmetic.
- Reconstructed the zero-mode complete tensor writer and locked read independently;
  checked their normal operators, ranks, determinant, and left inverses.
- Expanded the exact cosine symbol and independently evaluated the
  `(A,t1,t2,t3)` determinant through degree four.
- Evaluated the full `T2` normal through degree two.  This check caught an
  important possible overstatement: although the aggregate second offset
  moment and scalar determinant coefficient are isotropic, the tensor normal
  contains a cubic term already at degree two.
- A representation-scope review caught that `T2` alone is not closed under
  `SO(3)`: its cubic-looking `T-T` block may be the restriction of a
  rotationally covariant `E2+T2` operator.  The physical-anisotropy claim was
  removed and replaced by the exact unresolved completion diagnostic.
- A parent scope correction caught a second important overstatement before
  the repaired seal: the full canonical direct gradient cannot be promoted
  to a complete arbitrary-profile `A1/E` order-six writer.  All physical
  access and inverse claims were rebuilt from `D` on `A1+E` and `B P_T` on
  `T2`; the determinant changed from `3,670,016` to `524,288` while rank six
  and the smooth-field theorem survived.
- Derived all four leading relative-sector minors and their Cauchy--Binet sum
  of squares.  A Hadamard recombination proves the generic and face-diagonal
  rank statements.
- Verified exact sine-column dependencies on all six signed Cartesian
  face-diagonal families.
- Proved an exact finite-momentum common-sector rank loss without numerical
  trigonometric approximations.
- Checked the uniform pair/storage identity on all sixteen local bit words
  and matched the denominator derivative to six local writer vertices.

## Claim guards

This author audit does not substitute for independent hostile audit.  The
packet proves the complete `T2` Fourier source-access map and its obstructions.  It does not
claim that the source is record-generated, that a stationary susceptibility
or phase exists, that the cubic term disappears under refinement, or that
the symbol is spacetime, a metric, Ricci curvature, gravity, or `G`.

The canonical-direct unprojected Fourier row is not described as the full
effective writer.  Its `E` component, and nonuniform `A1` completion, remain
unclassified.  The common-field condition is a stated physical restriction, not something
derived from ring incidence.  The exact unsoldered rank ceiling is retained
to prevent that restriction from being hidden.
