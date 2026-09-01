#!/usr/bin/env python3
"""Fail-closed custody verifier for the sealed GL6AN hostile audit."""

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
    "independent_gl6an_replay.py",
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
        check(target.is_file() and not target.is_symlink(),
              f"{ledger} regular target: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"{ledger} hash: {relative}")
        members.add(relative)
    if ledger == "AUDITED_TARGETS.sha256":
        check(len(members) == 11, "exact audited author target count")
    else:
        for name in required:
            if name not in ("MANIFEST.sha256", "SEAL.sha256"):
                check(f"{HERE.name}/{name}" in members,
                      f"manifest coverage: {name}")

seal_lines = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
              if line.strip()]
check(len(seal_lines) == 1, "one seal entry")
expected, relative = seal_lines[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets audit manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest() == expected,
      "seal hash")

normalized = " ".join((HERE / "AUDIT.md").read_text().split())
check("Hostile verdict: PASS" in normalized, "audited disposition")
check("actual singular value is the square root" not in normalized,
      "no stale unsourced phrasing")
check("corresponding singular value is the square root" in normalized,
      "Gram/singular distinction")
check("not a six-record operation" in normalized, "record ceiling")
check("not an authenticated finite open mission" in
      " ".join((ROOT / "LANE_CROSS_RFT_GRA_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/THEOREM.md").read_text().split()),
      "finite-quotient ceiling")
check("gravity, or `G`" in (HERE / "AUDIT.md").read_text(), "physics ceiling")

print(f"PASS__GL6AN_HOSTILE_AUDIT_PACKET__{checks}/{checks}")

