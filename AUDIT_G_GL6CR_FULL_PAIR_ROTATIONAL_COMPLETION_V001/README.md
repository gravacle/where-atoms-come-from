# GL6CR independent hostile audit

This packet independently reconstructs the complete algebraic content of
`GL6CR` without importing or executing the author derivation.

Disposition: `PASS`.

The replay enumerates all 24 tetrahedral port permutations and all 126 raw
symmetric-kernel/quadratic-momentum seeds.  It independently obtains the
nine-dimensional `S4` response space, the four-dimensional rotational
subspace, the five exact completion tests, and the rank-eight direct Ward
system with the unique Einstein/Fierz--Pauli null ray.

Run:

```text
python3 independent_gl6cr_audit.py
python3 verify_packet.py
```

This is an algebraic classification result.  It does not derive the Ward
identity from F3, calculate the complete same-state response, construct its
physical 1PI/quotient kernel, prove gravity, or calculate `G`.
