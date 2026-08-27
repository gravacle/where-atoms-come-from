# Independent scientific-contract QA — V002 Repair2

**Audit:** `GRA-O-GF-CONTRACT-V002-REPAIR2-VERIFY-CODEX`  
**Date:** 2026-08-22  
**Disposition:** **REFUTED AT FIRST REMAINING SCIENTIFIC-STATE INCONSISTENCY**  
**Scientific result:** `NO_PROOF_OUTPUT`

## Controls reproduced

- V002 manifest: **PASS, 12/12**.
- Repair1 manifest: **PASS, 6/6**.
- Repair2 manifest: **PASS, 6/6**.
- Fixed Repair2 suite: **PASS, 54/54 twice**, zero failures, identical result
  payload SHA-256
  `564429f111bd975572a9bee443c334fe4630c7d3abc2c1d2cd28898700dfb947`.
- The exact Repair1 scalar escape was reproduced with
  `REAL.PLATFORM.FORMATION.observations=[0]`: Repair1 accepted the candidate,
  derived `actual_platform_present=true`, and promoted GF0–GF3 plus GE1–GE2.
  Repair2 refused the same candidate with
  `FORMATION_OBSERVATION_MEMBER_NOT_OBJECT`, derived
  `actual_platform_present=false`, and emitted seven `NO_PROOF_OUTPUT` values.
- The typed three-stage positive control passed. Fifteen independent mutations
  confirmed refusal for event and surface instability, unknown surface, invalid
  stage/role, naive or unordered time, uncertified or invalid value, unstable unit,
  unjoined or unstable source identity, incomplete or over-complete lifecycle,
  duplicate predicate identity, and non-reproducibility.
- The sealed synthetic fixture, a platform under `SYNTHETIC_TEST`, a platform with
  `synthetic_only=true`, and an incomplete lifecycle all remained non-authoritative.
- GF0 and GF1 passed in the typed compatibility control without a gravitational
  geometry, spacetime metric, or stress-energy field. The transport `W` calibration
  is an information/transport metric, not a gravitational metric requirement.

## First remaining scientific-state inconsistency

Repair2 can silently repair an invalid submitted formation-to-source custody join
and then promote it as an actual platform.

Starting from Repair2's own typed platform control, the audit changed only this
submitted parent reference:

```text
artifact:       REAL.PLATFORM.FORMATION
parent:         REAL.PLATFORM.SOURCE
field:          parents[].sha256
correct value:  252769ceef4757cb00b12f94a561b118fe4aba37dbaf5a8f611acd112bfa107b
submitted value:0000000000000000000000000000000000000000000000000000000000000000
```

No payload byte, artifact digest, binding, observation, or scientific field was
changed. Repair1 correctly evaluated the submitted candidate as:

```text
accepted                    = false
custody.disposition         = REFUSE
diagnostic                  = ARTIFACT_PARENT_HASH_MISMATCH:
                              REAL.PLATFORM.FORMATION:REAL.PLATFORM.SOURCE
actual_platform_present     = false
GF0 through UGE             = NO_PROOF_OUTPUT
```

Repair2 evaluated that same submitted candidate as:

```text
accepted                    = true
custody.disposition         = QUALIFIED
custody.errors              = []
formation_observation_valid = true
actual_platform_present     = true
GF0, GF1, GF2, GF3          = PASSES_DECLARED_DOMAIN
GE1, GE2                    = PASSES_DECLARED_DOMAIN
UGE                         = NO_PROOF_OUTPUT
```

The cause is the interaction among
`validator_v002_repair2.py:72`, `:195-202`, and `:226-228`:

1. the original artifact validator reports the parent-hash mismatch in
   `artifact_errors`, but observation validation does not fail on it because the
   evidence payload still decodes;
2. `_project_for_repair1` then rewrites **every** resolved `parents[].sha256` in the
   copied package to the registry digest, including the malformed submitted
   formation/source edge; and
3. Repair1 receives only the normalized copy, so its custody check can no longer see
   the original mismatch.

This violates the source-custody join, the incomplete-package non-promotion barrier,
and the stated actual-platform conjunction. Repair2 therefore cannot receive an
overall PASS even though its direct observation-member checks behave as documented.

## Scope and custody

The audit stops at this first authority-bearing inconsistency. No V002, Repair1, or
Repair2 subject byte was modified. All audit bytes are confined to
`V002/REPAIR2/VERIFY_CODEX/`; the adversarial candidate was generated and evaluated
in memory.
