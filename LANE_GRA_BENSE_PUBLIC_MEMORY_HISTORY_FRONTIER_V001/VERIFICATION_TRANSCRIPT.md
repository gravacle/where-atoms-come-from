# Verification transcript

**Date:** 2026-08-27  
**Protocol:** `GRA-BENSE-PMHF-V001`

## Deterministic physics replay

Command:

```text
python3 analyze_public_history.py --check
```

Output:

```text
PASS: RESULT.json matches deterministic recomputation
```

## Independent verification

Command:

```text
python3 verify_lane.py
```

Output:

```text
PASS: 432/432 independent checks passed
```

The independent verifier does not import the analyzer. It reparses the four
raw CSVs, independently reconstructs branch directions and boundaries,
recomputes all 97 opposite-branch comparisons and their pointwise digests,
recomputes force-event indices and return extrema, checks causal typing and
claim ceilings, verifies source/core manifests and the lane seal, and then
runs the analyzer's deterministic replay as a separate check.

## Syntax, JSON, and control-byte scan

Command class:

```text
PYTHONPYCACHEPREFIX=/private/tmp/bense_pmfh_pycache python3 -m py_compile analyze_public_history.py verify_lane.py
parse every lane JSON
scan all text/code/data/manifest files for forbidden control bytes
```

Output:

```text
PASS: Python compile, JSON parse, and text control-byte scan
```

The binary primary-source PDF is hash-checked rather than interpreted as text.

## Sealed load-bearing hashes

```text
754b0d4d1f6b6dc64751cf5c359bda1bef498377f2bb72de65061d5596a8fd70  SOURCE_MANIFEST.sha256
87a5d67c58508ccad984e7b39c3e7ebb4bbca5c48a8164f01cffc3515769274b  CORE_MANIFEST.sha256
4063539f62a6698beb4c82b15e2ce13ed0978f30353a49c9cac3836265bbbb6a  RESULT.json
9eef510038cf541f70b22082613231024e42d81e17818fb0cf8b866df540f75c  THEOREM.md
48dceb6cb61dc1d1fad60805b0fe26734223b03ca11758bb04ef2e999f95941b  HOSTILE_AUDIT.md
53bdf3dd2a1417e1c267c58e65afda50cb343a19dfea42f6bbfc7b85015072e4  analyze_public_history.py
fefdca4e29fabf76a9841bbb6e89e7323fd66c7ab317a71381a72c1b34b1e21c  verify_lane.py
```

`LANE_SEAL.json` reproduces and binds these hashes. No canonical model,
program register, or git file was edited by this isolated lane.
