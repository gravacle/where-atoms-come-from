#!/usr/bin/env python3
"""Run the fixed 38-check generic formation-input validation chain."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
CHECKS = (
    (
        "GENERIC FORMATION INPUT CORE",
        "checks_formation_input.py",
        "FORMATION_INPUT_CORE_CHECKS: 32/32 PASS",
    ),
    (
        "PUBLIC URM FORMATION DELEGATES",
        "checks_formation_input_urm.py",
        "FORMATION_INPUT_URM_CHECKS: 6/6 PASS",
    ),
)


def main() -> int:
    passed = 0
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    print("VALIDATE GENERIC URM FORMATION INPUT — NO SCIENTIFIC VERDICT")
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
        f"FORMATION INPUT CONTRACT: {'PASS' if overall else 'FAIL'} "
        f"({passed}/{len(CHECKS)} blocks; 38 checks; scientific weight ZERO)"
    )
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
