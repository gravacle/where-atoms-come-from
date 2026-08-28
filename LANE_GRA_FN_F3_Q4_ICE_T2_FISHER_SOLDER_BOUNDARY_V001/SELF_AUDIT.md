# Self-audit -- q4/F3 ice `T2` Fisher-solder boundary

**Lane:** `GRA-FN-F3-Q4-ITFSB-V001`

**Date:** 2026-08-27

## Verdict

`ACCEPT_AS_EXACT_LOCAL_SECOND_JET_AND_DERIVABILITY_BOUNDARY__DO_NOT_PROMOTE_TO_PHYSICAL_METRIC_SOLDER_OR_TENSOR_MODE`

## 1. What was independently recomputed inside this packet

The verifier reconstructs the six two-in/two-out states and proves that they
are the three antipodal doubled axes in `V=1^perp`.  It then computes the
uniform second, third, and fourth moments with rational arithmetic.  All 81
components of the Fisher Hessian are checked against the exact cumulant
formula.  Its `A1`, `E`, and `T2` eigenvalues are checked on explicit rational
bases, and its six-dimensional polarized rank is checked without numerical
tolerances.

For an arbitrary rational interior distribution, the verifier checks

\[
 F=4\operatorname{diag}(w)-mm^{\mathsf T},
\]

the complemented-state covariance, and the exact total-covariance identity
for an unlabelled complement mixture.  It checks rank two at a
complement-preserving normalized point, rank three after adding a scalar,
rank five at a generic broken normalized point, and rank six only after the
scalar is added.  The determinant `2m_1m_2m_3` of the off-diagonal response
block is reproduced as an exact symbolic monomial calculation.

## 2. Most likely overclaim and how it is blocked

The tempting overclaim is: "the second derivative spans `Sym^2(V)`, so the
physical metric solder is complete."  That is false.  The Hessian has a
six-dimensional **polarized second-jet span**, while a state with one source
vector still lies on a three-dimensional quadratic image.  More decisively,
the exact all-orders formula shows that every off-diagonal `T2` entry is
`-m_i m_j`.  At a complement-symmetric state `m=0`, the first-order `T2`
tangent is absent for every probability family, not just the exponential
curve at the uniform point.

At a generic broken background a full augmented rank is real, but it is
obtained by linearizing a vector dyad.  It requires nonzero mean on all three
ice axes and an independently supplied scalar.  It is therefore not an
independent tensor degree of freedom, and it does not evade the Maxwell
spin-one pole classification when that infrared input applies.

## 3. Scope and dependency risks

The finite theorem depends only on the final `FK` theorem and its independent
audit.  The conditional pole interpretation additionally uses the final
`FL` theorem and independent audit, whose bytes are now pinned separately in
the theorem and verifier.  That interpretation is explicitly separated from
the exact local theorem and remains conditional on `MAXWELL-IR`.

For an arbitrary interior state `p`, `F(p)` is the covariance/Fisher matrix
of the fixed three-parameter one-link exponential tilt based at `p`.  The
generic rank-five/rank-six result concerns the map from background-state
controls (and an independently supplied scalar) into the six covariance
components.  It is not a six-parameter Fisher metric and does not count
propagating modes.

The result is local and diagonal.  It does not exclude a derived nonlocal or
noncommuting information query.  It also does not establish state
preparation, scalar ownership, thermodynamic stability, `O(3)` restoration,
or a continuum tensor Ward identity.

## 4. Promotion ceiling

Safe canonical language is:

> The odd one-link source first changes the local ice Fisher metric at second
> order.  Its polarized Hessian has full `A1+E+T2` tetrahedral rank, but the
> `T2` component is exactly the dyad of a complement-breaking vector mean.
> Every complement-preserving local diagonal family remains `T2`-silent at
> first order.  A generic nonzero-flux background plus an independent scalar
> gives algebraic rank six, but only as vector-background response, not as an
> independently derived metric or graviton.

Unsafe language includes "metric derived," "emergent graviton," "helicity
two," "rotationally invariant nonlinear response," or "gravity closed."
