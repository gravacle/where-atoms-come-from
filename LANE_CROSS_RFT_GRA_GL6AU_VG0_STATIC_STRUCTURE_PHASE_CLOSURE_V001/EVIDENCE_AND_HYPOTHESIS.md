# GL6AU evidence and hypothesis ledger

## A. Proved in this lane

1. The exact first-character cycle norm is
   `||C^*u||^2=12sin^2(pi/L)`.
2. PF positivity bounds every orientation expectation by `0<=t_d<=1`.
3. The same-component single-mode estimate is
   `Delta_C<=6Jsin^2(pi/L)/S_T(2pi/L)`.
4. Any static exponent below two closes the selected finite-component gap.
5. For `S_T>=s/L`, a real quadrature has variance at least `(s/2)L^2` and
   the gap bound is `O(J/L)`.
6. The PF ground-state transform makes the missing variance a property of
   the exact reversible stationary law `pi=psi^2`.
7. RK equality, the elementary cubic half-filled twist, sector-tower closure,
   and continuous-link rigorous gauge theorems do not transfer the missing
   lower bound without additional premises.

## B. Controlled scientific evidence

The exact GL6AT crosswalk places the order-six model at the same `v/g=0`
fully-packed-loop quantum-ice Hamiltonian studied by Shannon et al. and
Benton et al.

- Shannon et al.: zero-temperature GFMC plus finite ED place zero well inside
  the numerically inferred liquid region and find Coulomb-form flux scaling.
- Benton et al.: equal-time QMC and finite-size energy fits calibrate a
  Gaussian lattice theory with two transverse linear modes and link weight
  proportional to frequency, giving the expected `S_T(q)~|q|` law.

These are direct phase evidence and a quantitatively calibrated effective
description.  They are not a rigorous all-`L` lower bound, a connectivity
proof, or direct real-time microscopic spectroscopy.

## C. Named hypotheses available for adoption

`ICE0-STATIC`:

```text
selected translation-stable global ground components exist, and
S_T(2pi/L)>=s/L for all sufficiently large L with s>0.
```

Consequence already proved: `Delta_C(L)<=6pi^2J/(sL)`.

`ICE0-GNS-BRIDGE`:

```text
the selected PF states converge compatibly to a ground state, and their
finite-character trials yield normalized form-domain excitations outside
the full zero-energy GNS subspace with energy tending to zero.
```

Consequence with `ICE0-STATIC`: zero GNS spectral gap in the selected
pure-model representation.

Neither named hypothesis is silently counted as a theorem.  The evidence
supports them at different strengths: strong numerical/effective support for
`ICE0-STATIC`, and field-theory rather than microscopic proof for the GNS
bridge.

## D. Forbidden promotions

The comparison-model phase label does not supply the all-orders F3 phase,
calibrated physical momentum or speed, a common physical cone, stress
coupling, Ricci/Einstein response, gravity, or `G`.
