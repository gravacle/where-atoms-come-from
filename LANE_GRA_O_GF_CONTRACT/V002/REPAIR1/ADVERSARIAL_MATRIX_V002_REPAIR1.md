# Repair1 adversarial matrix

| Case | Legacy behavior / attack | Repair1 requirement |
|---|---|---|
| exact accepted counterexample | sealed V002 accepts and promotes GF0–GE2 | reproduce legacy promotion, then refuse closed-schema failure |
| empty platform variant | empty ID, surface list, and map list | refuse nonempty-identity conjunction |
| unknown surface variant | platform lists a surface without measured process/map identity | refuse exact process/surface/map join |
| cross-package variant | otherwise joined platform claims another package | refuse package-scope mismatch |

All four cases require `actual_platform_present=false` and seven authoritative
`NO_PROOF_OUTPUT` results. The sealed 43-case V002 suite also runs unchanged under
the overlay, making 47 fixed cases total.
