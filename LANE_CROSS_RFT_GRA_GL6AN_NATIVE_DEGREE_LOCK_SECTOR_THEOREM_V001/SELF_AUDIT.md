# Self-audit

## Scope controls

- The theorem starts from the exact GL6AK Hamiltonian and authenticated
  parent/shared-child incidence.
- The lock line is tied back to the inherited coefficient:
  `Delta=4*U_d*(d_star-2)`.  The strict positive-`Delta`, positive-`U_d`
  subfamily requires `d_star>2` and has the explicit witness
  `d_star=3, Delta=4*U_d`.
- The mutable GL6AL N=4 data motivate the parameter ray but are not a premise.
- No conventional phase name, photon, graviton, ice rule, Ricci tensor, or
  continuum field equation is imported.
- The exact pair-sector statement is explicitly restricted to `q_v=0`.
- The perturbative statement is explicitly restricted to `U_d>0` and
  `h/U_d -> 0`; it is not promoted to the finite point `(-6,1)`.
- The linear-charge no-go is not overextended to every nonlinear or emergent
  conservation law.
- Complement symmetry is asserted for the infinite degree-four product
  automorphism and degree-four periodic quotients, not generic finite-open
  boundary products.
- Translation characters are not called physical momentum.
- `4-|s|` is typed as a constraint-Gram eigenvalue/squared singular value;
  the unsquared singular value vanishes linearly.
- Strong-lock projectors and extensive coefficients are defined on the
  explicit period-four, girth-at-least-six, degree-four quotient checked by
  the verifier.  They are not promoted to smaller quotients with wrapped
  four-cycles.  The hexagon coefficient is a finite linked matrix element
  embedded in the infinite incidence, not an infinite locked-sector
  projector.
- The six flipped active-link occupations are not promoted to six
  independently authenticated records.

## Exact checks

`verify_native_degree_lock.py` performs 1056 exact checks using only the Python
standard library.  They cover:

1. the completed-square identity for all local occupations `k=0,...,4`, five
   independent global occupation patterns, and exact inherited-domain
   admissibility;
2. full degree-four complement symmetry, an explicit finite-open partial-flip
   counterexample, and the exact Pauli commutator;
3. original-node degree, link endpoints, unique pair ownership, line degree,
   incidence rank, cycle-space dimension, and the linear-charge kernel;
4. the rank-two local pair kernel, its exact `E` eigenvalues, all six locked
   local states, and exact uniform covariance `(0,8/3,0)`;
5. exact Gaussian-rational samples of the incidence-character symbol and the
   quadratic centered-phase identity;
6. a chordless native hexagon, the absence of two/four-cycle port
   coincidences, and a separate modulo-four no-wrapped-four-cycle check;
7. every intermediate subset energy, all 720 sixth-order paths, an independent
   subset recursion, and the exact coefficient `-63/8`;
8. a deterministic period-four bipartite b-matching whose infinite lift is a
   degree-two background containing the target alternating hexagon; and
9. the exact finite-quotient adjacency occupation census and canonical scalar
   coefficients `H2/M=-h^2/(2U_d)` and
   `H4/M=-(7/24)h^4/U_d^3`.

Current result:

```text
GL6AN exact verification: PASS (1056/1056)
```

## Known open risks

1. A fresh independent hostile replay of the repaired bytes is pending.
2. The complete degenerate effective operator through sixth order is not
   derived.  The second- and fourth-order terms are scalar; sixth/higher
   diagonal and loop terms may gap or otherwise reorganize the locked-sector
   response.
3. The thermodynamic spectral threshold of the constrained effective theory
   is unknown.
4. The continuous incidence flat directions need not be quantum zero modes.
5. The finite `r_U=1` comparator is outside controlled strong-lock
   perturbation theory.
6. No physical momentum calibration, common cone, complete-stress Ward law,
   Ricci comparison, gravity identification, or `G` calculation is present.

## Freeze decision

The repaired author bytes are frozen by the fail-closed manifest/seal.  That
seal establishes custody only.  Do not promote the result until a fresh
independent hostile audit checks the repaired algebra, perturbative
coefficients, finite-projector typing, locked-background construction, and
stated ceilings.
