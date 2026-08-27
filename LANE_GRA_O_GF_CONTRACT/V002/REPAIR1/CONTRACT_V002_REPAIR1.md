# V002 Repair1 — actual-platform conjunction

**Repair:** `GRA-O-GF-CONTRACT-V002-REPAIR1`  
**Date:** 2026-08-22  
**Trigger:** independently reproduced synthetic-to-scientific platform promotion  
**Scope:** `PLATFORM_INSTANTIATION` semantics only  
**Scientific result:** `NO_PROOF_OUTPUT`

## Refutation reproduced

Sealed V002 accepted this correctly hashed artifact after changing only
`package_mode` to `SCIENTIFIC` and binding it:

```json
{
  "platform_id": "",
  "surface_ids": [],
  "synthetic_only": false,
  "platform_map_artifact_ids": [],
  "freeze_time": ""
}
```

The legacy evaluator returned `accepted=true`,
`actual_platform_present=true`, and promoted GF0 through GE2. Repair1 preserves that
exact candidate as a regression and refuses it at the platform schema/conjunction.

## Derived actual-platform predicate

`actual_platform_present` is true if and only if every following antecedent passes:

1. package mode is `SCIENTIFIC`, package custody is qualified, and the platform is
   explicitly non-synthetic;
2. platform ID, package ID, surface list, map list, formation-evidence list, process
   list, custody list, and timezone-aware freeze time are nonempty and closed;
3. the platform package ID equals the candidate package ID;
4. process IDs equal the bound independently measured `Gamma_rec` processes and the
   platform surfaces equal their unique measured surface IDs;
5. every surface has exactly one resolved `PLATFORM_MAP` of the correct kind;
6. every map joins that surface's process, the sealed `JOINT_SEED` definition, full
   D dilation, calibrated V004 transport, formation evidence, and source-custody
   object to the same platform and package;
7. seed, dilation, and transport IDs equal the top-level bindings; the dilation is
   measured on a platform surface and has an exchange ledger; the transport's
   physical surface is in the platform;
8. formation evidence is a reproducible `FORM.ALLOW0` observation whose source bytes
   belong to the platform custody object;
9. every `SOURCE_DATA` leaf in the process/seed/dilation/transport/formation ancestry
   is covered by a non-synthetic, content-addressed custody object with acquisition,
   specimen, independent-unit, provenance, and license identities;
10. the platform freeze is valid, no later than the package freeze, and before both
    response access and validation numeric access; and
11. FORM, GAMMA, SEED, SCALE, and ANCESTRY already derive `PASS` from the sealed V002
    evaluator.

Failure of any antecedent sets `actual_platform_present=false`, refuses the candidate,
and derives every milestone as `UNSCOREABLE` with authoritative
`NO_PROOF_OUTPUT`. The platform boolean and labels carry no authority by themselves.

## Version custody

Every sealed V002 file and `MANIFEST_V002.sha256` remains byte-identical. Repair1 is
an isolated overlay in this directory with its own validator, suite, result, audit,
and manifest. It changes no theory, seed choice, gate definition, proof conjunction,
shared URM code, dataset, register, ledger, or plan.

The actual admitted platform remains absent. `SELECT_JOINT_SEED` is still sealed, and
GF0 through UGE remain `NO_PROOF_OUTPUT`.
