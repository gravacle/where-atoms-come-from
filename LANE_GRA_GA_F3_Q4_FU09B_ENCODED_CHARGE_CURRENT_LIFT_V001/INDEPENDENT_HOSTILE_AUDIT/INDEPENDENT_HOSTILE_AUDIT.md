# Independent hostile audit: GA encoded charge/current lift

**Date:** 2026-08-28  
**Target:** repaired and re-sealed GA core pinned by `TARGET_CUSTODY.sha256`  
**Verdict:** **PASS, at the bounded GA claim ceiling**

## Scope and repaired defects

This audit independently tested the FU09b encoder, internal and outer current
signs, reservoir-placement boundary, conditional source inheritance, and all
listed failure modes.  It did not modify the author packet.

Two source-preservation defects were found before the audited target was
re-frozen:

1. source independence alone was incorrectly sufficient under the first hold;
   a constant `mu Q_R` is source independent but restricts to `-mu q_* Z` and
   changes the source-off Hamiltonian and virtual gaps; and
2. the first repair allowed nonzero encoded identity terms while retaining a
   literal Hamiltonian/Feshbach equality.  A common `c[j]I` instead requires a
   common scalar across the full encoded `P+Q` Hilbert space, a consistent
   Feshbach reference shift `z -> z+c[j]`, and claims only modulo the resulting
   identity/reference term.

The repaired theorem, result, self-audit, machine-readable result, verifier,
manifest, and seal now retain both corrections.  This PASS applies only to
those repaired bytes.

## Independent results

- The one-link isometry spans exactly the complete two-dimensional zero-total-
  charge eigenspace.  It intertwines both inherited generators,
  `Z_tilde V=VZ` and `X_tilde V=VX`, while the two off-code charge sectors are
  dark.  One reservoir charge state cannot cancel both link charges, so a
  two-level paired reservoir is dimension-minimal for this fixed-total-charge
  witness.
- For `H_flip=-h X_tilde`, independent commutator algebra reproduced
  `qdot_L=2ihq_*(A-B)/hbar`, `Qdot_R=-qdot_L`, and the orientation
  `I_(L->R)=-qdot_L`.  The current operator is nonzero and Hermitian, and its
  encoded action intertwines the inherited link-charge motion.
- Tensoring four independent pairs reproduced all q4 generator
  intertwiners.  With the frozen conventions, all six coordinates of
  `16 Q_pair/U_d` and `3 Q_flip/h` intertwine at the correct signs and
  normalizations.
- Adding a nonzero six-component identity source changed only the encoded
  identity/reference part.  Exact centering restored the original source in
  every coordinate and preserved its nonidentity Gram matrix and rank.
- A separate rational one-dimensional `P/Q` Feshbach calculation proved
  `H_eff_tilde(z+c)=H_eff(z)+cI` for one common full-code scalar.  Giving `P`
  and `Q` different scalars changed the virtual denominator, confirming that
  the repaired full-code condition is load bearing.
- The co-located four-support construction closes the Fourier-weighted
  link-to-reservoir exchange equation and conserves the full `m=1` density.
  Assigning all reservoirs one common displaced phase leaves a nonzero
  mismatch, while global `m=0` charge remains conserved.  This demonstrates
  allocation dependence; it does not choose a connector path.
- The explicit link/reservoir/exterior model reproduced all three oriented
  boundary equations and total U(1) conservation.  Its reservoir-only outer
  exchange has a nonzero off-code block, so an active port does not preserve
  the closed encoded hold.
- The failure classification is exact: `Q_R^2` is scalar in the minimal qubit;
  `Q_R` restricts to `-Z`; `Q_(R,1)Q_(R,2)` restricts to `Z_1Z_2`; and a
  strain-dependent Hermitian transfer derivative supplies a new nonidentity
  encoded source.

## Hostile claim review

No remaining material sign, normalization, leakage, source, or operator-
typing defect was found.  The final packet correctly states that:

- the link-to-reservoir exchange is an internal scalar U(1) current, not a
  diamond-bond current, vertex divergence, stress current, or `T^{0j}`;
- FU09b does not select reservoir placement, connector ownership, or a
  code-preserving active port;
- exact FV--FY inheritance concerns nonidentity physics modulo one declared
  full-code identity/reference shift, not literal equality unless `c[j]=0`;
- independent paired qubits are an existence witness, not a theorem about an
  arbitrary shared reservoir, grounded conductor, reset mechanism, or support;
  and
- visible QED, a Ward packet, continuum locality, gravity, and `G` remain
  outside this result.

## Sealed ceiling

GA proves a finite fixed-total-charge realization of the FU09b dressed flip,
its internal and outer U(1) balance, placement nonuniqueness, and conditional
preservation of the established nonidentity source physics.  It does **not**
derive a spatial bond current, vertex divergence, `T^{0j}`, a gauge or metric
Ward identity, continuum propagation, gravity, or `G`.

