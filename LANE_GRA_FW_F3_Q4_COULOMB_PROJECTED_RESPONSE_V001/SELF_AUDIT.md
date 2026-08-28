# Self-audit

## Scope accepted

This packet makes one finite source-suboperator calculation only: the
homogeneous `FV-WITNESS` sum of `Q_pair^(0)` and the irreducible differentiated
`Q_ring^(6)` on FO's exact 180-state leading-`H6` component and unique sector
ground state.  It keeps the FU/FV `S1`--`S10` and `FV-PURE` conditionality
explicit.

## Main hostile checks

1. **The witness was not called the complete source.**  FV retains generated
   `Q_diag^(2,4,6)` derivatives and differentiated folds whose values were not
   needed for its rank theorem.  FW omits them explicitly.  Its rank two is
   not an upper bound on the complete fixed-order response.
2. **Static rank was not counted as response rank.**  The packet separately
   reports rank five modulo identity, source-to-commutator rank three, and
   ground-state retarded/`M1` rank two.
3. **Identities were removed.**  One uniform `E` combination is exactly a
   scalar on the selected connected component and is rejected from
   nonidentity rank.
4. **Conservation was not called propagation.**  The `A1` source is an energy
   rescaling plus a reference, and one `T2` source commutes with the finite
   Hamiltonian.  Neither contributes ground-state spectral weight.
5. **A dark state was not called a constraint.**  `Q_dark` fails to excite the
   selected ground state but does not commute with `H`; it is not a Ward
   identity.
6. **Degenerate eigenvalues were handled algebraically.**  The exact integer
   six-orbit Hamiltonian, closed-form ground vector, and rational polynomial
   energy projectors prove the two residue matrices and the absence of a
   third pole.  A complete numerical degeneracy sum is only an independent
   replay; no arbitrary eigenvector inside a degenerate subspace is scored.
7. **The source scale was not tuned.**  The analytic answer retains
   `rho=Ud/J6`; executable replays at `rho=1` and `rho=2` verify the separate
   direct/ring factors.
8. **FV's theorem was not contradicted.**  FV proves rank six using family and
   local matrix-element witnesses.  This packet restricts to one component
   and one uniform source, where one family direction becomes scalar.
9. **Finite poles were not called particles.**  The two algebraic gaps are
   discrete finite-sector poles, not a thermodynamic dispersion or massless
   tensor.
10. **No local-source structure was invented.**  FV owns a homogeneous affine
   derivative.  A nonzero-momentum response requires a separately frozen
   microscopic block ownership and is left to the next calculation.
11. **No H8 or Ward promotion was made.**  The present parent is leading H6.
    Full through-H8 CTP contacts, ports, constraints, limit order, and
    thermodynamic scaling remain open.

## Numerical boundary

The Hamiltonian, source support, exact commutators, and finite-component
identities are replayed directly.  The zero-momentum proof uses the exact
integer `6 x 6` Hamiltonian, a closed algebraic ground vector, exact rational
polynomial energy projectors, and exact denominator-twelve source blocks.
Floating evaluation of those identities uses conservative `5e-11` or tighter
tolerances.  Complete double-precision diagonalization of the 180-state
matrix, with inherited eigensystem residual below `2e-13`, independently
replays the proof.

## Withheld conclusions

The packet does not claim the complete fixed-order source response, a
complete-sector or thermodynamic no-go, a local or nonzero-momentum response,
a six-channel CTP phase, a Ward algebra, RGRL-B, a graviton, gravity, or `G`.
