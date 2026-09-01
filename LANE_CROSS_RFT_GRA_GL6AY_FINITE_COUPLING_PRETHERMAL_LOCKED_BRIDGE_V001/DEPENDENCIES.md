# GL6AY local dependencies and custody

GL6AY uses three frozen, sealed, independently audited local inputs:

- `GL6AN` for the exact native reduction
  `H=U_d sum_v(k_v-2)^2-h sum_e X_e+C` and the finite local constraint
  geometry;
- `GL6AO` for the canonical locked perturbative coefficients through order
  six, including the first non-scalar coefficient
  `-(63/8)(h^6/U_d^5) sum_c T_c`; and
- `GL6AX` for exact contractible locked-to-locked conservation of all four
  port totals, the sharp `2L_min` wrapping threshold, and the controlled
  quasi-local twist theorem.

For every local input, `DEPENDENCIES.sha256` pins the author theorem,
manifest, and seal plus the distinct hostile-audit report, manifest, and
seal.  The exact external theorem versions are separately pinned in
`PRIMARY_SOURCES.md`.

No mutable later lane enters the proof.

