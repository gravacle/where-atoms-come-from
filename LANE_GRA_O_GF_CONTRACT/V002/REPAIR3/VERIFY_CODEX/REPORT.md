# Independent terminal QA — V002 Repair3

**Audit ID:** `GRA-O-GF-CONTRACT-V002-REPAIR3-VERIFY-CODEX`  
**Disposition:** `PASS`  
**Scientific result:** `NO_PROOF_OUTPUT`  
**First decisive defect:** none

## Acceptance result

The bounded independent audit found no surviving acceptance bypass in Repair3.
Repair3 evaluates the original candidate through Repair1 exactly once and first,
then gives the same original object to the observation validator exactly once. The
observation check preserved input bytes, the whole evaluation preserved input bytes,
and traps placed on Repair2's compatibility projection and projecting evaluator were
not triggered. The added layer therefore behaved monotonically: it removed
eligibility on observation failure and never repaired upstream custody.

## Exact Repair2 regression

The audit independently changed only the recorded digest on the
`REAL.PLATFORM.FORMATION -> REAL.PLATFORM.SOURCE` link. Repair1 refused the original
candidate with `ARTIFACT_PARENT_HASH_MISMATCH`. Repair2 reproduced the laundering
escape, accepted the projected copy, asserted an actual platform, and promoted GF0
through GE2. Repair3 evaluated the unmodified candidate, preserved the Repair1
refusal, set `actual_platform_present=false`, and emitted `NO_PROOF_OUTPUT` for every
GF0–UGE claim. Its observation layer was valid in this case but could not promote the
upstream refusal.

## Fixed and independent checks

- All six subject manifests verified: V001, the V001 independent verifier, V002,
  Repair1, Repair2, and Repair3.
- The sealed Repair3 suite ran twice: 59/59 passed on both runs, with byte-identical
  result payloads and payload digest
  `043761bf06a90025dff7940e91bd1af4ed72db515555eeece004de8e8f613223`.
- Independent wrong-child-digest, unknown-parent, self-cycle, orphan, scalar `[0]`,
  and mixed object/scalar cases all refused, preserved input bytes, and emitted only
  `NO_PROOF_OUTPUT`.
- The complete typed lifecycle passed as the required constructed positive control.
  It is a test-control fixture, not an admitted empirical platform and carries no
  scientific weight.
- With the platform binding absent, `actual_platform_present=false` and every
  authoritative GF0–UGE output is `NO_PROOF_OUTPUT`.

## Nongeometric boundary

GF0 requires custody, a platform candidate, Gamma, seed, calibrated transport, and
ancestry. GF1 adds `FORM.ALLOW0`, `GAMMA.PROCESS`, `SEED.PHI`,
`SCALE.CALIBRATED_TRANSPORT`, and `ANCESTRY.FULL_PATH`. Neither milestone requires a
gravity-characteristic gate or candidate field for spacetime geometry, a
gravitational metric, stress-energy, or `T_lab`.

The validator does contain an internal covariance/Gram compatibility calculation
named `metric_residual` inside calibrated transport. This is linear-algebraic
transport compatibility, not a spacetime or source-geometry premise, and no
candidate `metric` field is admitted or required at GF0/GF1.

## Terminal conclusion

Repair3 closes the concrete Repair2 custody-laundering defect under the fixed audit
scope. No further repair request is warranted. This is a mechanics result only:
without an admitted platform-specific `JOINT_SEED` packet, the authoritative
scientific disposition remains GF0 through UGE `NO_PROOF_OUTPUT`.
