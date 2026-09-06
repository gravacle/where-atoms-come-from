# GL6CY development — CU V002 raw source composition

**Status:** exact development calculation; author-complete only after the
fast verifier below passes.  The imported `GL6CU V002` bytes are
author-frozen but remain development evidence here until their distinct
hostile audit passes.

This directory is intentionally separate from
`LANE_CROSS_RFT_GRA_GL6CY_COMPLETE_PHYSICAL_SOURCE_REEMBEDDING_RESIDUAL_TARGET_V001`.
No byte in that target was changed.

The executable reconstructs the complete selected-Q4 raw Hamiltonian source
Jets from the pinned CU V002 engine and evaluates

\[
 {\cal R}^{\rm raw}_{sr}
 =H_{,AB}\eta_{A,s}\eta_{B,r}+H_{,A}\eta_{A,sr}
\]

for all six tensor test rows, all three native center directions, and all 64
characters of the inherited Q4 regulator.  It separates the bare term, the
`h2`, `h4`, and `h6` diagonal terms, and each of the 64 active alternating
six-cycle `h6` owners.  The nonlinear second source `eta_sr` is constructed
and contracted; it is never silently set to zero.

The same run compiles each arbitrary-nonuniform source Hessian into the exact
finite-range translation stencil

\[
 {\cal R}^{(n)}_{sr}(m)
 =\sum_{\Delta\in\mathbb Z_4^3}
 K^{(n)}_{sr}(\Delta)i^{m\cdot\Delta}.
\]

That stencil is the reusable local input for an overlap-family calculation
on `G_L`.  It is not a cell-level Ward inference and does not presume a
microscopic continuous position law.

## Scope ceiling

The computed object is a **raw Hamiltonian source-family composition**.  It
is not a connected response, a CTP kernel, a 1PI kernel, a Schur/Legendre
quotient, a physical Ward residual, an Einstein/Fierz--Pauli kernel, gravity,
`C_R`, or `G`.

Run the fast custody/formula check with:

```text
python3 verify_development.py
```

Run the full exact owner and 64-character replay with:

```text
python3 derive_raw_source_composition.py
```

Use `--json` to emit the exact character-function classes, the exact
translation-difference stencils, and hashes for all 64 individual cycle
functions.  The full replay completed on 2026-09-04 with `21/21` checks in
`439.229` seconds.

