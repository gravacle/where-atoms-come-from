#!/usr/bin/env python3
"""Run the fixed 82-check public-data proof-frontier validation chain."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
CHECKS = (
    (
        "MISSING-DATA PROOF FRONTIER CORE",
        "checks_proof_frontier.py",
        "PROOF_FRONTIER_CORE_CHECKS: 76/76 PASS",
    ),
    (
        "PUBLIC URM PROOF FRONTIER DELEGATES",
        "checks_proof_frontier_urm.py",
        "PROOF_FRONTIER_URM_CHECKS: 6/6 PASS",
    ),
)


def main() -> int:
    passed = 0
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment.setdefault("PYTHONPYCACHEPREFIX", "/private/tmp/wac-proof-frontier-pycache")
    print("VALIDATE PUBLIC-DATA PROOF FRONTIER — PROOF-FIRST, ZERO PROOF PROMOTION")
    print("=" * 78)
    for label, script, expected in CHECKS:
        result = subprocess.run(
            [sys.executable, "-B", str(HERE / script)],
            cwd=HERE,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        output = result.stdout.rstrip("\n")
        if output:
            print(output)
        exact = result.returncode == 0 and output.splitlines()[-1:] == [expected]
        print(f"  {'PASS' if exact else 'FAIL'}  {label}")
        passed += int(exact)
    overall = passed == len(CHECKS)
    print("=" * 78)
    print(
        f"PROOF FRONTIER CONTRACT: {'PASS' if overall else 'FAIL'} "
        f"({passed}/{len(CHECKS)} blocks; 82 checks; authoritative proof outputs ZERO)"
    )
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
