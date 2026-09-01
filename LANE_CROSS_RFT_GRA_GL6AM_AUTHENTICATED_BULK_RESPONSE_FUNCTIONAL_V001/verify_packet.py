#!/usr/bin/env python3
"""Fail-closed structural and custody check for the mutable GL6AM packet."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


required_files = (
    "README.md",
    "THEOREM.md",
    "RESULT.md",
    "SELF_AUDIT.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "RESPONSE_LEDGER.json",
    "verify_authenticated_bulk_response.py",
    "verify_packet.py",
    "VERIFICATION.txt",
    "MANIFEST.sha256",
)
for name in required_files:
    check((HERE / name).is_file(), f"required file: {name}")

manifest_paths = set()
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"manifest target exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"manifest hash: {relative}")
    manifest_paths.add(relative)

for name in required_files:
    if name != "MANIFEST.sha256":
        relative = f"{HERE.name}/{name}"
        check(relative in manifest_paths, f"manifest coverage: {name}")

for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"dependency exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"dependency hash: {relative}")

ledger = json.loads((HERE / "RESPONSE_LEDGER.json").read_text())
check(ledger["lane"] == "GL6AM V001", "ledger lane")
check("INDEPENDENT_REVIEW_REQUIRED" in ledger["status"], "mutable status")

theorem = (HERE / "THEOREM.md").read_text()
self_audit = (HERE / "SELF_AUDIT.md").read_text()
result = (HERE / "RESULT.md").read_text()
check("independent hostile review\nrequired before freeze" in theorem,
      "theorem must remain unsealed")
check("Remaining blocker" in self_audit, "hostile blocker disclosure")
check("generally nonequilibrium" in result, "defect scope in result")
check(not (HERE / "SEAL.sha256").exists(), "self-reviewed packet must not self-seal")

print(f"PASS__GL6AM_PACKET__{checks}/{checks}")
