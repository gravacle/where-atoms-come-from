# GA encoded-current self-audit

## Algebra and equivalence

- The encoder is an isometry onto the complete zero-total-charge subspace of
  one link plus one reservoir qubit.
- Equivalence is claimed only on that invariant encoded subspace; dark states
  in the four-dimensional pair Hilbert space are not discarded silently.
- Tensor-product composition is checked on all four q4 links, including the
  degree square, flips, pair source, and flip source.
- Feshbach inheritance is conditional on an encoded invariant projector and
  the same source-dependent intertwiner.

## Current and placement

- The derived operator is a scalar link-to-reservoir U(1) exchange current.
- It is not called a spatial bond current, vertex divergence, stress current,
  `T^{0j}`, or metric Ward current.
- Co-located midpoint support is one explicit witness, not a selected
  physical reservoir placement.
- A displaced reservoir is checked to require a nonzero connector current at
  `m=1`, while global `m=0` charge remains conserved.
- A zero outer current applies only during the closed hold.  The explicit
  active outer-port model retains its boundary current and is shown to leak
  from the code.

## Source discipline

- `FV-PURE/FY` preservation is stated only under
  `GA-CLOSED-FULL-CODE-SCALAR-HOLD`: every added term is one common scalar on
  the full encoded `P+Q` Hilbert space at source off and has at most one common
  identity first spatial derivative.
- Scalar terms are retained modulo reference identity, not called literally
  zero.  The Feshbach reference shifts with `c[j]`, virtual gaps cancel it, and
  only nonidentity coefficients/ranks/commutators/responses are unchanged.
- Source independence is explicitly rejected as a sufficient condition; the
  constant `mu Q_R -> -mu q_* Z` counterexample changes source-off gaps.
- Local quadratic reservoir charging, linear bias, shared charging, active
  outer ports, and source-dependent transfer/support terms are distinguished.
- The exact negative examples prevent promotion from an existence witness to
  an arbitrary physical grounded-reservoir theorem.

## Claim ceiling

No hardware implementation, autonomous reset, full FU `S3--S4`, visible QED,
spatial gauge Ward identity, stress tensor, continuum limit, gravity, or `G`
is claimed.
