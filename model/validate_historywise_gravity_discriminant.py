#!/usr/bin/env python3
"""Run the fixed 84-check historywise-gravity formal-discriminant chain."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
CHECKS = (
    (
        "PINNED FORMAL DISCRIMINANT CORE",
        "checks_historywise_gravity_discriminant.py",
        "HISTORYWISE_GRAVITY_DISCRIMINANT_CORE_CHECKS: 60/60 PASS",
    ),
    (
        "ZERO-INPUT URM DELEGATES AND INTEGRATED CONJUNCTION",
        "checks_historywise_gravity_discriminant_urm.py",
        "HISTORYWISE_GRAVITY_DISCRIMINANT_URM_CHECKS: 24/24 PASS",
    ),
)


def main() -> int:
    passed = 0
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment["PYTHONPYCACHEPREFIX"] = "/private/tmp/wac-historywise-discriminant-pycache"
    print("VALIDATE HISTORYWISE-GRAVITY FORMAL DISCRIMINANT — ZERO PHYSICAL PROOF PROMOTION")
    print("=" * 86)
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
    print("=" * 86)
    print(
        f"HISTORYWISE-GRAVITY FORMAL DISCRIMINANT: {'PASS' if overall else 'FAIL'} "
        f"({passed}/{len(CHECKS)} blocks; 84 checks; physical/empirical proof weight ZERO)"
    )
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())

