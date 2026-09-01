# GL6AO Result

Starting only from sealed GL6AN, the complete canonical locked-sector
Hamiltonian on the declared `Q_4` quotient is

```text
H_eff = C P
  -(M/2)(h^2/U_d)P
  -(7M/24)(h^4/U_d^3)P
  -(h^6/U_d^5)[
      (893M/1080)P
      +(63/8) sum_{c in Hex(Q_4)} P(product_{e in c}X_e)P
    ]
  +O(h^8/U_d^7).
```

Here `M=256`.  A distinct locked configuration is reachable in six flips if
and only if the symmetric difference is one alternating six-cycle.  Exact
enumeration of all `720` flip orders gives its universal amplitude `-63/8`.

The diagonal calculation is also complete.  All repeated-pair and
three-distinct-link return words were classified, and the canonical folds
cancel every `M^3` and `M^2` term.  The residual is the common scalar
`-(893/1080)M`; it is independent of the locked configuration.  Thus order
six contains no diagonal flippability potential.

The cycle toggle has a finite-support formulation using local degree-two
projectors, so it defines a formal uniformly finite-range linked interaction
on the infinite incidence.  This is a microscopic collective-dynamics
theorem, not yet a theorem of an infrared phase, pole, physical cone,
gravity, or `G`.
