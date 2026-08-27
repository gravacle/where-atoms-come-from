#!/usr/bin/env python3
"""Run the fixed 48-check origin-neutral gamma-flow validation chain."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
CHECKS = (
    (
        "GENERIC GAMMA FLOW CORE",
        "checks_gamma_flow.py",
        "GAMMA_FLOW_CORE_CHECKS: 42/42 PASS",
    ),
    (
        "PUBLIC URM GAMMA FLOW DELEGATES",
        "checks_gamma_flow_urm.py",
        "GAMMA_FLOW_URM_CHECKS: 6/6 PASS",
    ),
)


def main() -> int:
    passed = 0
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    print("VALIDATE ORIGIN-NEUTRAL URM GAMMA FLOW — REPAIR3 BOUNDED")
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
        f"GAMMA FLOW CONTRACT: {'PASS' if overall else 'FAIL'} "
        f"({passed}/{len(CHECKS)} blocks; 48 checks; synthetic proof weight ZERO)"
    )
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
