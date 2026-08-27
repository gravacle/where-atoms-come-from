# AQ4DL terminal custody seal

**Final custody status:**
`FINAL_SEALED_EXACT_CONDITIONAL_SELECTOR__NOT_AN_EM_OR_GRAVITY_DERIVATION`

**Date:** 2026-08-27

## Two-layer custody

The eight-file builder packet was frozen and independently hostile-audited
without changing its theorem, verifier, verification record, or builder
manifest.  Their embedded `MUTABLE` / `NOT_SEALED` labels record the stage at
which those bytes were frozen.  This additive custody layer is the
authoritative post-audit status record.

Frozen anchors:

- theorem:
  `ae0ed7bdd758f14e830612a5f6f7dc0207efbe4b4ad0b0a8152fe793c1a99a0d`
- deterministic verifier:
  `6f94e3e3a88bd188b8b94d0ab10facfc37897d80367fe7a5eb4fdc4168e8a747`
- builder verification record:
  `c63856edca54475cc446dbcaed6d0791e18bbb67fde7c581db8c5a53e4887bcc`
- frozen eight-file builder manifest:
  `6c862b17a3e4409b99687292dc5edc3dc19ad83c09d54947b40cca2b987dc7c8`
- independent hostile audit:
  `55d17928fbea3e2681b9da6fc580013699111352a3d97354864ae6a273429911`

The frozen verifier reported `PASS 238/238` on the builder packet and
`PASS 241/241` after the independent audit was added.  With this custody
record, terminal checker, and terminal manifest also present, its unchanged
dynamic byte-hygiene census reports `PASS 250/250`; the original eight-entry
builder manifest remains frozen and valid.

`TERMINAL_MANIFEST.sha256` is an additive full-packet ledger.  It covers all
eight frozen builder sources, the builder manifest, independent audit, this
custody record, and the terminal checker.  The terminal checker authenticates
both custody layers.

The audit's nonblocking convention note is retained: in the displayed
alternative gauge normalization, `g_D=e_D` after fixing the same physical
charge convention.  This is a normalization identification, not a new
scientific premise or a post-freeze theorem edit.

## Claim ceiling

The sealed theorem is an exact conditional REQUIRE-side selector:
`QFRONT-DIM` plus a same-front canonical Maxwell sector whose inherited charge
is classically marginal without an engineering-scale compensator implies
`D=q=4`.  It does not prove those actual-world premises, produce four reusable
operations, derive electromagnetism, calculate alpha's numerical value, or
derive gravity.
