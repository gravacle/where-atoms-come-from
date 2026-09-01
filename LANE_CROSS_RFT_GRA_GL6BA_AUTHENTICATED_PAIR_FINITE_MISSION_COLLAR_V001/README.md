# GL6BA — authenticated-pair finite-mission collar

This packet executes the direct moderate-`R` branch identified by `GL6AZ`.
It compares the exact authenticated pair read under the full F3 Hamiltonian
with the same full Hamiltonian restricted to a finite physical collar.  It
does not use the out-of-domain `GL6AY` prethermal approximation.

For every finite coupling ratio `R=U_d/h`, finite dimensionless mission time
`sigma_obs`, and desired binary-probability tolerance, the theorem supplies a
finite collar radius with an explicit state-independent error.  The admitted
members `R=2` and `R=5/2` are therefore directly treatable; neither is chosen
as the physical member.

The primary theorem compares every larger complete finite all-formed/`MATCH`
FPSS exterior
with its induced finite collar using the same postformation state and its
exact reduction.  The collar cut is a proof device, not a separately
performed physical switch.  The infinite-volume statement is only the
quasi-local mathematical corollary supplied by `GL6AK`, not an infinite
authenticated record.

The main new refinement is a boundary-only Duhamel split.  The exact `A_3`
cell-collar boundary has `12(3L^2+3L+1)` inherited crossing terms, so the
resulting pair-marginal certificate is

```text
D_TV(p^Omega,p^(L)) <= min(1,
  (3 L^2+3 L+1) T_(2L+1)(48 R |sigma_obs|)).
```

No graviton, Ricci target, Einstein equation, gravity identification, or
Newton constant enters.
