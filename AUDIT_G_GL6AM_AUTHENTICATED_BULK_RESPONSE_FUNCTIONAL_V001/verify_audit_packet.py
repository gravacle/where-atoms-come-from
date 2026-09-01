#!/usr/bin/env python3
"""Fail-closed custody verifier for the sealed GL6AM hostile audit."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


required = (
    "README.md",
    "AUDIT.md",
    "AUDITED_TARGETS.sha256",
    "independent_gl6am_replay.py",
    "VERIFICATION.txt",
    "verify_audit_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
)
for name in required:
    check((HERE / name).is_file(), f"required audit file: {name}")

for ledger in ("AUDITED_TARGETS.sha256", "MANIFEST.sha256"):
    members = set()
    for line in (HERE / ledger).read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        target = ROOT / relative
        check(target.is_file(), f"{ledger} target exists: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"{ledger} hash: {relative}")
        members.add(relative)
    if ledger == "MANIFEST.sha256":
        for name in required:
            if name not in ("MANIFEST.sha256", "SEAL.sha256"):
                check(f"{HERE.name}/{name}" in members, f"manifest coverage: {name}")

seal_lines = [line for line in (HERE / "SEAL.sha256").read_text().splitlines() if line.strip()]
check(len(seal_lines) == 1, "one seal entry")
expected, relative = seal_lines[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets audit manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest() == expected,
      "seal hash")

audit = (HERE / "AUDIT.md").read_text()
normalized_audit = " ".join(audit.split())
check("Hostile verdict: PASS" in normalized_audit, "audited disposition")
check("no physical momentum calibration" in normalized_audit, "physics ceiling")
check("generally not invariant under a defect generator" in normalized_audit,
      "nonequilibrium ceiling")

print(f"PASS__GL6AM_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
