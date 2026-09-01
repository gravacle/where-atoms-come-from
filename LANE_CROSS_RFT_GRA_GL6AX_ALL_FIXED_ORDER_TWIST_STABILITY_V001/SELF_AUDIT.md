# GL6AX self-audit

## Exactness attacks

1. **Disconnected symmetric differences:** the affine-moment proof uses the
   complete endpoint difference, not connectedness; disjoint multi-loop
   changes conserve every port total.
2. **Winding counterexample:** a periodic `j/3` alternating row has length
   `2L_j` and changes `N_j,N_3` by `-/+L_j`.  Port conservation is therefore
   not falsely extended to all periodic updates.
3. **Sharp threshold:** every parent-to-parent two-edge step moves each
   lifted coordinate by at most one, so nonzero winding needs at least
   `2L_min` edges; the row construction attains it.
4. **Folded and Feshbach terms:** fixed-order conservation uses only the final
   locked endpoints.  Locality is stated separately and is not inferred for
   a raw reducible Feshbach representation.  The Schrieffer--Wolff/Kato claim
   uses the canonical word-based gauge; an arbitrary extra locked-subspace
   unitary can obscure the displayed charge blocks.
5. **Diagonal terms:** arbitrary configuration-dependent diagonal terms
   commute with the twist.  They can change sector selection but cannot add
   twist energy within the declared centered sector.
6. **Multi-loop terms:** bounded-support simultaneous loop changes are
   included through the local double-commutator density.  Arbitrarily
   separated products require a linked or quasi-local interaction bound.
7. **Complex amplitudes:** `+/-` twist averaging cancels all odd-in-twist
   contributions.  No reality, zero-current, or time-reversal premise is
   used.
8. **Interaction tails:** the quasi-local statement includes the explicit
   wrapping-tail norm `T_L`; divergent second moments and system-spanning
   terms are not silently included.

## Promotion attacks

1. Fixed-order stability is not called a convergent all-orders theorem.
2. A formal order-by-order conserved port total is not called an exact
   finite-coupling symmetry of the microscopic link-flip Hamiltonian.
3. Centered-sector closure is not called selection of the actual ground.
4. An anisotropic finite-size dichotomy is not called selected-GNS
   gaplessness or an isotropic low-character mode.
5. No finite twisted trial is called a record-authenticated physical mode.
6. No photon, graviton, Ricci tensor, Einstein equation, gravity law, or `G`
   is inferred.

## Order-of-limits attack

The theorem fixes the perturbative truncation `K` before sending the box size
to infinity.  It does not exchange the thermodynamic limit with
`K -> infinity`.  Doing so requires a uniform quasi-local block
diagonalization, a dressed locked-sector map, and explicit wrapping-tail
control.
