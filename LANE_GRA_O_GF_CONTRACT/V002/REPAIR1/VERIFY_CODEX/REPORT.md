# Independent scientific-contract QA — V002 Repair1

**Audit:** `GRA-O-GF-CONTRACT-V002-REPAIR1-VERIFY-CODEX`  
**Date:** 2026-08-22  
**Disposition:** **REFUTED AT FIRST SCIENTIFIC-STATE INCONSISTENCY**  
**Scientific result:** `NO_PROOF_OUTPUT`

## Controls reproduced

- Sealed V002 manifest: **PASS, 12/12**.
- Repair1 manifest: **PASS, 6/6**.
- Exact empty-platform regression: sealed V002 accepted it, derived
  `actual_platform_present=true`, and promoted GF0–GF3 plus GE1–GE2; Repair1
  refused it with `PLATFORM_INSTANCE_SCHEMA_INVALID`, set
  `actual_platform_present=false`, and emitted seven `NO_PROOF_OUTPUT` values.
- Fixed Repair1 suite: **PASS, 47/47 twice**, zero failures, identical payload
  SHA-256 `6126960565ddb29f9754811dca16dacbac2062499b0e08a04f1dacbe849ef61b`.
- The four published platform regressions (exact empty packet, empty fields,
  unknown surface, and cross-package platform binding) all refuse as documented.

## First scientific-state inconsistency

Repair1 does not require every formation observation to be an observation object.
Starting from its own closed `joined_platform_candidate()`, the audit changed only
`REAL.PLATFORM.FORMATION.observations` to `[0]`, then recomputed that inline
artifact's canonical hash and parent closure with the suite's own helper.

Observed result:

```text
accepted                    = true
custody.disposition         = QUALIFIED
platform_repair1.valid      = true
actual_platform_present     = true
FORM.ALLOW0                 = PASS
GF0, GF1, GF2, GF3          = PASSES_DECLARED_DOMAIN
GE1, GE2                    = PASSES_DECLARED_DOMAIN
UGE                         = NO_PROOF_OUTPUT
```

Expected result under Repair1 clauses 8 and 56–58 is refusal,
`actual_platform_present=false`, `PLATFORM_FORMATION_EVIDENCE_INVALID`, and seven
`NO_PROOF_OUTPUT` values.

The cause is the predicate in `validator_v002_repair1.py:197`: it checks bad values
only for entries satisfying `isinstance(item, dict)`. A nonempty list containing no
dictionaries therefore passes the `any(...)` test vacuously. The base artifact
validator constrains the top-level `GATE_EVIDENCE` keys but does not close the nested
observation shape, so no earlier check rejects `[0]`.

This is an incomplete platform package that Repair1 classifies as actual and uses to
emit authoritative GF/GE proof outputs. Consequently the full actual-platform
conjunction and the synthetic/incomplete-package authority barrier are not closed.
The audit stops at this first inconsistency; later conjunction claims cannot receive
an overall PASS from this version.

## Scope and custody

No V002 or Repair1 subject file was modified. All audit bytes are confined to
`V002/REPAIR1/VERIFY_CODEX/`; the adversarial candidate was generated and evaluated
in memory.
