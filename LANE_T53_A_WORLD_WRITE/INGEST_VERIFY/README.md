# T-53A world-observation ingest verifier

Run from any directory:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -B "/Users/bgm/MB Work/where-atoms-come-from/LANE_T53_A_WORLD_WRITE/INGEST_VERIFY/verify_ingest.py"
```

The script freezes the reviewed source/artifact hashes, runs the shipped checks and actual
CLI three times, reconstructs all 760 normalized rows directly from the ten raw instrument
files without importing the project VSM adapter, and executes 36 mutation attacks in owned
temporary directories.

Outputs:

- `RESULT.json`: complete machine-readable evidence, including every refusal/acceptance.
- `RESULT.txt`: short default-refuted verdict.
- `D24.md`: initial failures, corrections, current refuting findings, and declared boundaries.

The verifier neither scores record formation nor authorizes a proof, universal claim, or
public URM registration.
