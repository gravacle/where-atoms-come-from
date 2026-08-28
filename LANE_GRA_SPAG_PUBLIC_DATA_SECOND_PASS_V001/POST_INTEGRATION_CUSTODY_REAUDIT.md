# Post-integration custody re-audit

**Lane:** `GRA-SPAG-PUBLIC-DATA-SECOND-PASS-V001`  
**Date:** 2026-08-27  
**Verdict:** `ACCEPT_NONCIRCULAR_CONTEXT_BASELINE_REPAIR`

## Defect and repair

The first frozen verifier treated both repository experiment registers as live
hash-pinned dependencies. That is circular: the registers were legitimate
search-time inputs to the duplicate screen, but acceptance of this lane then
adds the lane's result to those same mutable registers. The integration changed
the register bytes without changing the search execution or result, so the old
verifier rejected a lawful downstream integration.

The repair preserves the exact two search-time register seals as
`contextual_register_baselines` and removes only their live-byte comparison.
The historical seals are:

- `GRAVITY_EMERGENCE_EXPERIMENT_REGISTER_V001.md`:
  `ae7d2e672b3ba59f9c93d160c8562c95541b8e4c00d68fb21bf27d5315c2b58c`;
- `PROGRAM_EXPERIMENT_REGISTER_V001.md`:
  `084ed543ddaa5c55fd90e12a76c43688016aedd9218fe8850c34a5b30514c0d5`.

Both hashes independently reproduce from the immutable pre-integration commit
`a7e54f1cc5295bcfe885415f006521356f683627`. The verifier now requires that
commit identifier, those exact historical values, and both live paths, while
continuing to recompute every stable logical dependency. It does not replace
the baselines with moving post-integration hashes.

## Scientific-content invariance

The frozen public report, result, and first hostile audit remain byte-identical:

- `PUBLIC_DATA_SECOND_PASS.md`:
  `3d4300b9c2998aab4a485771f097f860e570a3931b8a948be9e1b034925931a8`;
- `RESULT.json`:
  `d643d38c60822be6bfb348082ce280ce482f075b01c93e9fe0a7a2053b4daa2f`;
- `INDEPENDENT_HOSTILE_AUDIT.md`:
  `668abf36aea57e3465562465a6b77dff1ce65f0a4a41ae548fbea2f20a90c5e5`.

All 28 query strings, both retained component roots, the zero retained
same-parent lineage-root count, all inclusion and exclusion gates, and the
protected Panda inventory are unchanged. In short, the queries, retained roots, and Panda holdout gates are unchanged.

## Ceiling re-audit

The repair does not add returned-hit-list custody, exhaustive-repository
custody, a physical null, a lineage estimand, an independent `G` estimate, or
permission to open the Panda response holdout. It changes dependency typing,
not evidence or scientific claims. Future edits to stable dependencies still
fail the verifier; future register integrations do not invalidate the exact
historical context against which this search was performed.
