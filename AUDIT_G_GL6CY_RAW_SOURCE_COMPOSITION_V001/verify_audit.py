#!/usr/bin/env python3
"""Fail-closed verifier for the bounded independent GL6CY source audit."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "GL6CY_DEVELOPMENT_CU_V002_RAW_SOURCE_COMPOSITION"
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


required = {"README.md", "independent_reconstruction.py", "verify_audit.py",
            "TARGET.sha256", "MANIFEST.sha256", "SEAL.sha256"}
for name in sorted(required):
    check((HERE / name).is_file(), f"audit file {name}")

target_lines = [line for line in (HERE / "TARGET.sha256").read_text().splitlines()
                if line.strip()]
check(len(target_lines) == 7, "seven target custody entries")
for line in target_lines:
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    check(path.is_file() and digest(path) == expected, f"target byte {relative}")

run = subprocess.run([sys.executable, "-B", str(HERE / "independent_reconstruction.py")],
                     cwd=ROOT, capture_output=True, text=True, check=False)
check(run.returncode == 0, "independent reconstruction exits zero")
for token in (
    "PASS__INDEPENDENT_GL6CY_RAW_SOURCE_COMPOSITION__",
    "REPLAY=FRESH_FULL_TARGET_64_CHARACTER_COMPARISON",
    "SOURCE=H_AB_ETA_S_ETA_R_PLUS_H_A_ETA_SR__SECOND_SOURCE_NONZERO_AND_CONTRACTED",
    "K0=H0_H2_ZERO__H4_32128_OVER_27_S__H6_DIAG_PLUS_CYCLE_7212448_OVER_6075_S",
    "CEILING=RAW_HAMILTONIAN_ONLY__NO_CONNECTED_CTP_1PI_WARD_GL6CR_GRAVITY_CR_OR_G",
):
    check(token in run.stdout, f"reconstruction output {token}")

print(f"PASS__AUDIT_G_GL6CY_RAW_SOURCE_COMPOSITION__{checks}/{checks}")
print("SCOPE=EXACT_DECLARED_RAW_SOURCE_COMPOSITION_ONLY")
