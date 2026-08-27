# V002 Repair2 — typed all-member formation observations

**Repair:** `GRA-O-GF-CONTRACT-V002-REPAIR2`  
**Date:** 2026-08-22  
**Trigger:** independently reproduced vacuous observation-list promotion  
**Scope:** actual-platform formation-observation members only  
**Scientific result:** `NO_PROOF_OUTPUT`

## Refutation reproduced

Repair1 accepted an otherwise joined platform after replacing
`REAL.PLATFORM.FORMATION.observations` with:

```json
[0]
```

Its generator examined only members already known to be dictionaries. The scalar was
filtered out, the empty generator returned no failure, the platform remained valid,
and GF0 through GE2 promoted. Repair2 preserves that exact candidate as a regression:
Repair1's escape is first reproduced, then Repair2 refuses it with
`FORMATION_OBSERVATION_MEMBER_NOT_OBJECT`.

## Closed observation member

For each actual-platform `FORM.ALLOW0` evidence artifact, `observations` must be a
nonempty list. Before reading any member field, Repair2 requires **every** member to
be an object. It then requires every member to have exactly:

```text
predicate_id  event_id  surface_id  stage  role  time
value         unit      source_artifact_id  reproducible
```

There is no filtering, coercion, defaulting, or truthiness fallback. Unknown and
missing keys are refused. Each identity is a nonempty string; the source resolves to
the evidence source list; the surface resolves to the platform; the timestamp is
timezone-aware and precedes the platform freeze; the value is a certified boolean or
finite real and must be exactly the frozen passing value; and `reproducible` is true.
Nonfinite JSON is refused before scientific evaluation.

## Lifecycle join

One scored formation event contains exactly one observation at each required stage:

```text
FORMATION   role FORMATION
CLOSURE     role CLOSURE
PERSISTENCE role PERSISTENCE
```

Across those three members, event ID, surface ID, unit identity, and source artifact
identity are stable; predicate IDs are unique; and timestamps strictly increase in
the listed order while remaining pre-freeze. Mixing lifecycle stages from different
events or surfaces cannot manufacture a formation pass.

Repair2 projects a successfully validated rich observation into Repair1's sealed
three-field compatibility view only after all-member validation. This preserves prior
V002/Repair1 scoring without letting their narrow view become the observation
authority.

## Nongeometric boundary recorded, not expanded

Repair2 adds no metric, geometry, tensor, stress-energy, mass, or microscopic
gravitational-source requirement to GF0 or GF1. Microscopic `JOINT_SEED` blocks remain
typed formation/process/current candidates; source-responsive geometry is a
conditional IR identification only.

The following sealed principal clarification is recorded as a normative parent for
the later gamma-flow schema, not as a new Repair2 predicate:

```text
LANE_GRA_S_JOINT_SEED_DECISION/PRINCIPAL_CLARIFICATION_001.md
sha256 2268762250f69c1ee8297ecd22cc3e67d490ccfb44062bf94ec18db9f768b8cb
```

## Version custody

Every V001, V001-verifier, V002, and Repair1 byte remains unchanged. Repair2 is an
isolated overlay with its own validator, suite, result, audit, and manifest. It
changes no seed selection, platform admission rule beyond observation typing, gate or
proof conjunction, physics endpoint, dataset, shared URM code, register, ledger, or
plan.

No actual admitted platform has been instantiated. GF0 through UGE therefore remain
`NO_PROOF_OUTPUT`.
