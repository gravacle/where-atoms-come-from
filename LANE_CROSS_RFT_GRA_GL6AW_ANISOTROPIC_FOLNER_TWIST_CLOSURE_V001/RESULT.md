# GL6AW result

The exact pure hexagon model has a new finite-size closure route that does
not assume the GL6AU static-structure lower bound.  On a rectangular torus
with even twist length, odd transverse area, and centered port occupation
`N_0=V/2`, a one-port large twist changes translation character by `-1`.

For a translation-stable PF ground component, the twisted vector stays in
the same component and is orthogonal to its ground state.  Exactly two
hexagon orientations per cell acquire the twist phase, giving

```text
Delta_C <=2JV[1-cos(2pi/L_parallel)]
        <=4pi^2J A_perp/L_parallel.
```

With periods `(m,2m^3,m)` for odd `m`, all dimensions tend to infinity and
`Delta_C<=2pi^2J/m`.  If the controlling component is not translation
stable, its translate is an exactly degenerate controlling component
instead.  Hence the centered sector has an exact anisotropic Følner
dichotomy: ground-component degeneracy or same-component gap closure.

This bypasses `ICE0-STATIC` only for finite-size anisotropic closure.  It
does not prove isotropic scaling, a low-character pole, or infinite-volume
selected-GNS gaplessness; the separate GNS bridge remains open.  No physical
cone, stress/Ricci law, gravity, or `G` is derived.
