# Builder self-audit

## Claim discipline

- The lane is conditional on `FV-PURE`; it does not silently weaken or
  discharge FV `S1`--`S10`.
- “Complete” is restricted to `Q_diag^(2,4,6)` plus the frozen irreducible
  ring, at homogeneous momentum, on the selected 180-state FO component,
  through H6, modulo `H_id`.
- Static coefficient rank is never called retarded rank.  Operator,
  source-to-commutator, ground retarded, and `M_1` ranks are separated.
- The exact dynamic proof is the Fraction-exact pair-plus-identity reduction
  composed with audited FW algebra.  Float eigensystem work is an independent
  checksum only.

## Completeness gates

- Every even multiplicity partition through length six is enumerated.
- `(4)` and `(6)` are explicitly checked to have zero irreducible orders.
- `(4,2)` and `(2,4)` are separately included.
- The `(2,2,2)` family includes all ninety unique multiset orders before
  intermediate-`P` filtering.
- Exact endpoint-referenced gap derivatives and all numerator derivatives
  are included.
- The complete BW/Feshbach derivative folds through H6 are implemented from
  the self-consistency equation, including `2 a2 da2 k22`.
- Translation is proven label- and incidence-preserving before one
  representative per equal-size orbit is used.
- The normalized zero-momentum block is checked with `1/sqrt(30)` weights.

## Scale gates

- The enumerator works at `U_d=1`; the theorem restores `x=h/U_d` and
  `J_6=(63/8)U_dx^6` order by order.
- `Q=-2 dH/dj` is applied before comparing with the pair source.
- The exact Hilbert identities are retained, displayed, and removed only
  from nonidentity response ranks.

## Root discipline

`f_E(x)` is first treated as a formal power-series unit.  The unique finite
through-H6 cancellation is separately proved by exact monotonicity and an
exact rational bracket; NumPy supplies only its decimal checksum.  It is not
a physical threshold, a convergence statement, or a gravity onset.

## Explicit ceilings

This packet does not prove a local/nonzero-momentum source, full-ice-sector
universality, H8 completion, CTP response, a Ward identity, RGRL-B, a
thermodynamic massless tensor, gravity, `G`, or Newton's constant.
