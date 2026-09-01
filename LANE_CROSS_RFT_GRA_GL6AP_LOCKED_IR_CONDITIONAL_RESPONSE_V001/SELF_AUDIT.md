# Self-audit

## Exact scope controls

- Only sealed GL6AN and its distinct hostile audit are logical inputs.
- The strict inherited domain remains `U_d>0`, `d_star>2`,
  `Delta=4*U_d*(d_star-2)`.
- The link kernel and local pair sector are kept in their inequivalent
  representations; a shared dimension two at generic character is not used
  as an identification.
- `4-|s|` is a Gram eigenvalue/squared singular value.  The singular value is
  linear.
- The loop nonconservation result uses a complete degree-two `Q4` background,
  not an arbitrary local assignment.
- The upgrade from one nonzero loop witness to no conserved uniform `E`
  component uses the `S4`-invariance of the commutator kernel and the
  irreducibility of `E`.
- The general quadratic kernel is conditional on a selected invariant state
  and character analyticity.  The even reciprocal expansion is separately
  conditional on a reciprocal kernel.
- Fixed `Q4` is not assigned a continuum or infrared limit; all
  long-wavelength statements require a growing-quotient or thermodynamic
  completion.
- Every locked projector is finite-volume `P_Q`; infinite-volume constraint
  statements use local null constraints and no global projector.
- Ground-state spectral criteria are not silently applied to arbitrary KMS
  or merely stationary states.
- A positive atom at each character is distinguished from a residue bounded
  away from zero in the infrared limit.
- The single-mode quotient excludes zero-frequency elastic weight from both
  moments before bounding the least positive support point.
- No complete sixth-order effective Hamiltonian is claimed.
- The nonconservation equation is written for `H_eff^(6)` and is not promoted
  to an uncomputed all-orders operator.
- No character is called physical momentum.
- No gauge variable or named gauge/photon/graviton phase is imported.

## Exact checks

`verify_locked_ir_conditions.py` independently checks:

1. the `S4` characters of the port `T2` and opposite-pair `E`, including
   `Hom_S4(T2,E)=0`;
2. `Sym^2(E)=A1+E`, `Sym^2(T2)=A1+E+T2`, the absence of a linear-gradient
   invariant, and the allowed `E` cubic invariant;
3. the exact rank-three trivial and rank-two generic incidence kernels;
4. the centered quadratic character identities and Gram/singular scaling;
5. the three locked opposite-pair types and exact `E` plane;
6. a complete degree-two period-four background with an alternating target
   hexagon whose loop toggle changes type counts by `(2,0,-2)`; and
7. theorem custody and strict claim ceilings.

## Open physics

1. The full connected sixth-order diagonal and loop operator is neither
   imported nor rederived within this GL6AN-only lane.
2. No thermodynamic ground/KMS/other physical state is selected.
3. The `E` mass, stiffnesses, damping, threshold, and spectral residue are
   unknown.
4. No physical embedding/refinement calibrates character as momentum.
5. No common cone, stress law, Ricci comparison, gravity, or `G` follows.

## Freeze decision

The author packet is frozen and sealed after its exact replay.  Do not promote
it before an independent hostile audit checks the representation mismatch,
loop nonconservation witness, quadratic invariant count, spectral criteria,
finite-volume gates, and ceilings.
