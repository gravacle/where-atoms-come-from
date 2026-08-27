# Final T-53A scientific data-integrity recheck

Run from any directory:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -B "/Users/bgm/MB Work/where-atoms-come-from/LANE_T53_A_WORLD_WRITE/INGEST_VERIFY_FINAL/verify_final.py"
```

This final, default-refuted QA pass freezes the repaired input implementation, reconstructs
the 760-row table independently from the ten raw Lake Shore files, compares it byte for byte
with the registered adapter and stored normalized table, repeats the official checks/CLI,
and evaluates the complete retained set of malformed and substituted input fixtures.

It emits deterministic `RESULT.json`, `RESULT.txt`, and `D24.md`. It does not score a theory
or authorize a record-formation proof, universal claim, independent reproduction, or public
URM registration.
