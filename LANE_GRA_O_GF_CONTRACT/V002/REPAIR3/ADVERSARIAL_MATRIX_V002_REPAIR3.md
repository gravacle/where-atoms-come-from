# Repair3 adversarial matrix

| Case | Mutation | Required result |
|---|---|---|
| exact parent laundering | corrupt the digest on `REAL.PLATFORM.FORMATION -> REAL.PLATFORM.SOURCE` | reproduce Repair2 promotion; Repair3 exposes parent-hash mismatch and refuses |
| wrong child digest | replace `REAL.PLATFORM.FORMATION.sha256` without changing its bytes | refuse artifact hash mismatch |
| unknown parent | replace the formation source parent identity with `UNKNOWN.PARENT` | refuse unresolved parent |
| parent cycle | add the formation artifact as its own parent | refuse self-loop/cycle |
| orphan | append a valid, content-addressed, unbound source artifact | refuse unreachable artifact |

Every new case asserts `accepted=false`, `actual_platform_present=false`, GF0–UGE
`NO_PROOF_OUTPUT`, the exact expected diagnostic class, and byte identity of the
candidate across evaluation.

The inherited 54-case inventory is rerun through Repair3. In particular, `[0]` and a
mixed typed/scalar observation list still refuse. A fully typed
FORMATION/CLOSURE/PERSISTENCE lifecycle is an explicit positive control and remains
qualified; this prevents a reject-all implementation from passing. Total counted
cases: 59.
