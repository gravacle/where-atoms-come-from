#!/usr/bin/env python3
"""Validate the URM's measurement-input door without scoring any science.

The 44 checks cover the registered instrument adapter (11), the generic closed
observation contract (28), and public URM delegation (5).  Passing this executable
input gate is necessary infrastructure; it is not evidence that either theory is true.
"""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
CHECKS = (
    ("LAKESHORE VSM ADAPTER", "checks_lakeshore_vsm.py", "LAKESHORE_VSM_CHECKS: 11/11 PASS"),
    ("WORLD OBSERVATION CONTRACT", "checks_world_observation.py", "WORLD_OBSERVATION_CHECKS: 28/28 PASS"),
    ("PUBLIC URM OBSERVATION DOOR", "checks_world_observation_urm.py", "URM_WORLD_OBSERVATION_CHECKS: 5/5 PASS"),
)


def main() -> int:
    passed = 0
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    print("VALIDATE URM WORLD-OBSERVATION INPUT — NO SCIENTIFIC VERDICT")
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
        if exact:
            passed += 1
    overall = passed == len(CHECKS)
    print("=" * 78)
    print(
        f"WORLD INPUT CONTRACT: {'PASS' if overall else 'FAIL'} "
        f"({passed}/{len(CHECKS)} blocks; 44 checks; scientific weight ZERO)"
    )
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
