# Independent hostile audit: GD translation-owning recoil parent

**Date:** 2026-08-28  
**Target:** repaired and frozen GD/TORP core pinned by `TARGET_CUSTODY.sha256`  
**Verdict:** **PASS, at the bounded B1 algebraic-existence ceiling**

## Scope and repaired defects

This audit independently tested the reciprocal-lattice domain, Weyl signs,
half-kick encoder, complete inherited `P+Q` reduction, kinetic scalar,
factor-edge momentum and charge currents, Feshbach common-reference shift,
conditional controller/boundary class, and same-mechanical-state no-go.  It
did not modify the author packet.

Before freeze, hostile review required the following repairs:

1. all paired recoil factors were placed on one explicitly shared auxiliary
   torus, with `kappa_e in 2 Lambda^*` and `p_e/hbar in Lambda^*`;
2. the no-go was narrowed from all product states to reuse of one identical
   mechanical state by both logical codewords;
3. the exact current was typed factor-edge-local on an auxiliary torus, not
   as an already localized physical diamond-space current;
4. the controller Hamiltonian was distinguished from F3's Coulomb `H_C`;
5. controller/boundary closure was stated as a conditional admitted class,
   not an explicit active controller or energy-current construction;
6. the internal transfer in the bare Weyl-current formula was required to
   commute with both mechanical momenta, with extra commutators retained
   otherwise; and
7. projector commutation, Hermiticity, off-sector invariance, and the full
   analytic square-summability argument were made explicit.

The frozen bytes retain all seven repairs.

## Independent results

- Direct three-dimensional momentum-state action, using several unrelated
  kick and center vectors, reproduced the two Weyl commutator signs, fixed
  pair momentum, charge neutrality, exact `X` intertwining, and equal
  equal-mass kinetic energy.
- A separate two-codeword Heisenberg calculation reproduced
  `dot P_L+J=0`, `dot P_R-J=0`, and
  `J_P=-(hbar kappa/2q_*) I_Q` component by component.
- A four-link enumeration covered all sixteen inherited configurations,
  including the off-ice/virtual sector.  Unequal kick vectors still gave one
  full-code recoil scalar and preserved rows containing diagonal fields,
  pair interactions, and every single-link flip.
- An independent rational Schur-complement calculation verified that a
  common `P/Q` reference shift cancels from the Feshbach resolvent and adds
  only one identity to the effective Hamiltonian.  A deliberately unequal
  `P/Q` shift changed the virtual denominator, confirming the stated hold
  boundary.
- An explicit multi-edge factor graph closed total impulse.  The bilateral
  shift recurrence gave partial norm
  `(2N+1)|c|^2`, confirming that a nonzero normalizable state cannot be reused
  unchanged by both logical codewords.  This does not exclude the displayed
  link/reservoir product half-kick states.

## Hostile claim review

No material defect remains.  The source and response inheritance is exactly
conditional on the full-code scalar hold: source-dependent kicks, masses,
potentials, controllers, or ports require a new source/rank calculation.
Likewise, a closed auxiliary recoil ledger does not derive the kick scale,
physical placement, a localized current, or a spacetime stress tensor.

The construction is therefore a valid B1 algebraic existence witness.  Its
arbitrary kick and delocalized momentum-sector code are inputs to B2, not a
substitute for B2.

## Sealed ceiling

GD proves that the inherited GA charge exchange can be dressed by an exact
translation-owning relative Weyl recoil while preserving complete GA/FV/FY
nonidentity physics modulo one common full-code reference identity.  It does
**not** prove physical diamond-space placement, a complete spacetime source,
`T^{0j}`, stress-Ward closure, a tensor cone, gravity, or `G`.
