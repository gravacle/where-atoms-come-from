# V002 Repair3 — immutable monotone custody evaluation

**Repair:** `GRA-O-GF-CONTRACT-V002-REPAIR3`  
**Date:** 2026-08-22  
**Trigger:** independently reproduced parent-hash laundering through Repair2's compatibility projection  
**Scope:** evaluation order and artifact-custody preservation only  
**Scientific result:** `NO_PROOF_OUTPUT`

## Refutation reproduced

In an otherwise qualified typed platform candidate, the recorded digest on the
`REAL.PLATFORM.FORMATION -> REAL.PLATFORM.SOURCE` parent link was replaced by an
incorrect 64-hex digest. Repair2 validated the rich observations, constructed a
compatibility copy, and rewrote every parent digest from the copied registry before
calling Repair1. That rewrite erased the custody contradiction. Repair2 accepted the
candidate, set `actual_platform_present=true`, and promoted GF0 through GE2.

Repair3 preserves that exact candidate as a regression. It first reproduces the
Repair2 escape, then evaluates the unmodified input through Repair1. The original
parent mismatch remains visible as `ARTIFACT_PARENT_HASH_MISMATCH`; the candidate is
refused, actual-platform authority is false, and every GF0–UGE proof output is
`NO_PROOF_OUTPUT`.

## Immutable evaluation order

Repair3 applies one monotone conjunction:

1. Call the sealed Repair1 evaluator exactly once with the original instance.
2. Run Repair2's typed all-member formation-observation predicate directly against
   that same original instance and its original content-addressed registry.
3. If Repair1 refuses, preserve that refusal and force actual-platform authority
   false. Observation success cannot repair it.
4. If Repair1 qualifies but observation typing fails, add those diagnostics, force
   actual-platform authority false, make every gate and milestone unscoreable, and
   emit no proof output.
5. If both antecedents pass, preserve Repair1's derived scientific state and attach
   the read-only observation result.

No candidate is reconstructed. No payload, artifact digest, byte length, parent
digest, parent identity, binding, or registry entry is rewritten. No projected
candidate is evaluated. Repair2's public observation diagnostic key is retained only
as a compatibility alias for the same single read-only observation check; Repair2's
projecting evaluator is never called by Repair3.

## Bounded custody guards

The fixed regressions require refusal for:

- the exact corrupted formation-to-source parent digest;
- a wrong digest on the formation child artifact;
- an unknown formation parent;
- a formation self-parent/cycle; and
- a correctly encoded but unreachable orphan artifact.

Each Repair3 case also asserts byte identity of the input before and after evaluation.
The inherited scalar-list and mixed-list attacks still refuse, while the complete
typed three-stage lifecycle remains a positive control. All prior 54 cases run under
Repair3 unchanged.

## Scientific and theory boundary

Repair3 changes no gate predicate, proof conjunction, JOINT_SEED selection, Gamma
process representation, calibrated/refinement transport rule, gravity characteristic,
data split, platform admission requirement, or product-reproduction requirement. It
adds no metric, geometry, tensor, stress-energy, mass, or microscopic gravitational-
source requirement at GF0 or GF1. The sealed nongeometric clarification remains:

```text
LANE_GRA_S_JOINT_SEED_DECISION/PRINCIPAL_CLARIFICATION_001.md
sha256 2268762250f69c1ee8297ecd22cc3e67d490ccfb44062bf94ec18db9f768b8cb
```

## Version custody

V001, its verifier, V002, Repair1, and Repair2 remain byte-identical to their sealed
manifests. Repair3 is confined to `V002/REPAIR3`. There are no shared URM, model,
ledger, register, proof, or plan edits.

No actual admitted platform packet is supplied by this repair. The authoritative
scientific disposition remains GF0 through UGE `NO_PROOF_OUTPUT`.
