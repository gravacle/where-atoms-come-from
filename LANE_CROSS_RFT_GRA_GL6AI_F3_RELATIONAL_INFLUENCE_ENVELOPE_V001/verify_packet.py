#!/usr/bin/env python3
"""Custody and scope verifier for a frozen GL6AI packet."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path


root = Path(__file__).resolve().parent
repo = root.parent
checks = 0

required = [
    "THEOREM.md",
    "RESULT.md",
    "README.md",
    "SELF_AUDIT.md",
    "PRESCREEN_REQUEST.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "ENVELOPE_LEDGER.json",
    "verify_f3_relational_envelope.py",
    "VERIFICATION.txt",
    "MANIFEST.sha256",
    "SEAL.sha256",
]
for name in required:
    assert (root / name).is_file(), name
    checks += 1

ledger = json.loads((root / "ENVELOPE_LEDGER.json").read_text())
assert ledger["certified_constants"]["lambda_F3"] == (
    "4 J Delta_L/hbar = 48 |Ud|/hbar"
)
assert "not exact finite-speed support" in ledger["ceiling"]
checks += 2

theorem = (root / "THEOREM.md").read_text()
for token in [
    "lambda_{\\rm F3}:={4J\\Delta_L\\over\\hbar}",
    "={48|U_d|\\over\\hbar}",
    "d_L(e,f)\\ge d_{\\mathcal G_N}(m,n)",
    "d_{\\rm cell}(s,Y)",
    "exponentially quasi-local",
    "exact finite-speed support or a stationary bulk mode",
    "No sharper constant is claimed",
]:
    assert token in theorem, token
    checks += 1

run = subprocess.run(
    ["python3", str(root / "verify_f3_relational_envelope.py")],
    check=True,
    text=True,
    capture_output=True,
)
assert "PASS 123777/123777" in run.stdout
checks += 1

for line in (root / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, rel = line.split("  ", 1)
    data = (repo / rel).read_bytes()
    assert hashlib.sha256(data).hexdigest() == expected, rel
    checks += 1

manifest_lines = [
    line.split("  ", 1) for line in (root / "MANIFEST.sha256").read_text().splitlines()
]
assert len(manifest_lines) == 11
checks += 1
for expected, rel in manifest_lines:
    data = (repo / rel).read_bytes()
    assert hashlib.sha256(data).hexdigest() == expected, rel
    checks += 1

seal_expected, seal_rel = (root / "SEAL.sha256").read_text().strip().split("  ", 1)
assert seal_rel == f"{root.name}/MANIFEST.sha256"
assert seal_expected == hashlib.sha256((root / "MANIFEST.sha256").read_bytes()).hexdigest()
checks += 1

print(f"PASS__GL6AI_PACKET__{checks}/{checks}")
