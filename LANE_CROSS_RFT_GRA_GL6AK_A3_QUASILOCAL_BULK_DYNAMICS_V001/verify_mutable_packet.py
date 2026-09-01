#!/usr/bin/env python3
"""Frozen author-packet custody and scope verifier for GL6AK V001."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
CHECKS = 0


def check(condition: bool, message: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(message)
    CHECKS += 1


required = {
    "THEOREM.md",
    "RESULT.md",
    "README.md",
    "SELF_AUDIT.md",
    "PRESCREEN_REQUEST.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "PRESCREEN_AUDIT.sha256",
    "BULK_DYNAMICS_LEDGER.json",
    "verify_a3_bulk_dynamics.py",
    "verify_mutable_packet.py",
    "VERIFICATION.txt",
    "MANIFEST.sha256",
    "SEAL.sha256",
}
for name in sorted(required):
    path = HERE / name
    check(path.is_file() and path.stat().st_size > 0, f"missing/empty {name}")

ledger = json.loads((HERE / "BULK_DYNAMICS_LEDGER.json").read_text())
check(ledger["status"].startswith("AUTHOR_FROZEN"), "ledger freeze status")
check(ledger["edges"]["site_degree"] == 6, "ledger degree")
check("48 |Ud|/hbar" in ledger["interaction"]["lambda_F3"], "ledger lambda")
check("no pole" in ledger["physics_ceiling"], "ledger ceiling")

theorem = (HERE / "THEOREM.md").read_text()
for token in (
    "author frozen after independent hostile pre-freeze review",
    "post-freeze custody/replay audit required",
    "3\\|A\\||X|\\sum_{r=R}^{\\infty}",
    "\\lambda_{\\rm F3}=24J/\\hbar",
    "\\bar\\omega\\tau_t=\\bar\\omega",
    "\\mu=\\mu_{A_1}P_{A_1}+\\mu_EP_E+\\mu_{T_2}P_{T_2}",
    "E=\\hbar\\nu",
    "CHARACTER_NOT_PHYSICAL_MOMENTUM",
    "AUTHOR_FROZEN_POSTFREEZE_AUDIT_REQUIRED",
):
    check(token in theorem, f"theorem token: {token}")

physics = subprocess.run(
    ["python3", str(HERE / "verify_a3_bulk_dynamics.py")],
    check=True,
    text=True,
    capture_output=True,
)
check("PASS 6304/6304" in physics.stdout, "author physics replay")

for ledger_name in ("DEPENDENCIES.sha256", "PRESCREEN_AUDIT.sha256"):
    for line in (HERE / ledger_name).read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = REPO / relative.removeprefix("../")
        check(path.is_file(), f"pinned member exists: {relative}")
        check(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
              f"pinned member hash: {relative}")

audit_text = (
    REPO / "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001" / "AUDIT.md"
).read_text()
audit_verification = (
    REPO / "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001" / "VERIFICATION.txt"
).read_text()
check("38cb58ef9fc52e1252e0b0d3415c54488c0471c0ddf25a35b7adf5aba41bccc9" in audit_text,
      "reviewed mutable theorem hash")
check("Hostile pre-freeze verdict: CLEAN" in audit_text, "clean prescreen verdict")
check("PASS 79644/79644" in audit_verification, "independent prescreen count")

manifest = [
    line.split("  ", 1)
    for line in (HERE / "MANIFEST.sha256").read_text().splitlines()
]
check(len(manifest) == 12, "manifest member count")
for expected, relative in manifest:
    path = REPO / relative
    check(path.is_file(), f"manifest member exists: {relative}")
    check(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
          f"manifest member hash: {relative}")

seal_hash, seal_relative = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal_relative == f"{HERE.name}/MANIFEST.sha256", "seal target")
check(seal_hash == hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest(),
      "seal hash")

print(f"PASS {CHECKS}/{CHECKS}")
