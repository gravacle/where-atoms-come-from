#!/usr/bin/env python3
"""Fail-closed custody verification for the distinct GL6AO audit packet."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001"
AUDIT_FILES = {
    "README.md", "AUDIT.md", "AUDITED_TARGETS.sha256",
    "independent_gl6ao_replay.py", "VERIFICATION.txt", "verify_audit_packet.py",
}
AUTHOR_FILES = {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md", "DEPENDENCIES.md",
    "DEPENDENCIES.sha256", "VERIFICATION.txt", "verify_sixth_order_hamiltonian.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def rows(path):
    parsed = []
    for raw in path.read_text().splitlines():
        if not raw.strip():
            continue
        expected, relative = raw.split(maxsplit=1)
        check(len(expected) == 64 and all(c in "0123456789abcdef" for c in expected),
              f"valid digest syntax: {relative}")
        check(relative not in {item[1] for item in parsed}, f"unique row: {relative}")
        target = ROOT / relative
        check(target.is_file(), f"target exists: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        parsed.append((expected, relative))
    return parsed


targets = rows(HERE / "AUDITED_TARGETS.sha256")
check(len(targets) == 11, "exact audited target row count")
check({Path(relative).name for _, relative in targets} == AUTHOR_FILES,
      "audited targets cover exact author packet")
check({path.name for path in AUTHOR.iterdir() if path.is_file()} == AUTHOR_FILES,
      "author packet has no undeclared files")

author_manifest = rows(AUTHOR / "MANIFEST.sha256")
check(len(author_manifest) == 9, "author manifest exact row count")
check({Path(relative).name for _, relative in author_manifest} == AUTHOR_FILES - {"MANIFEST.sha256", "SEAL.sha256"},
      "author manifest exact coverage")
seal_rows = rows(AUTHOR / "SEAL.sha256")
check(len(seal_rows) == 1 and Path(seal_rows[0][1]).name == "MANIFEST.sha256",
      "author seal targets author manifest")

dependency_rows = rows(AUTHOR / "DEPENDENCIES.sha256")
check(len(dependency_rows) == 6, "exact transitive dependency count")
check(all("GL6AN" in relative and "GL6AL" not in relative for _, relative in dependency_rows),
      "transitive custody confined to sealed GL6AN")

audit_manifest = rows(HERE / "MANIFEST.sha256")
check(len(audit_manifest) == 6, "audit manifest exact row count")
check({Path(relative).name for _, relative in audit_manifest} == AUDIT_FILES,
      "audit manifest exact coverage")
check({path.name for path in HERE.iterdir() if path.is_file()} == AUDIT_FILES | {"MANIFEST.sha256", "SEAL.sha256"},
      "audit directory has no undeclared files")
audit_seal = rows(HERE / "SEAL.sha256")
check(len(audit_seal) == 1 and Path(audit_seal[0][1]).name == "MANIFEST.sha256",
      "audit seal targets audit manifest")

replay = (HERE / "independent_gl6ao_replay.py").read_text()
check("import verify_sixth_order_hamiltonian" not in replay,
      "independent replay does not import author verifier")
audit = (HERE / "AUDIT.md").read_text()
for marker in (
    "K_6=T_6-bX_4+b^2A_3-dA_2",
    "-{893\\over1080}M",
    "-{63\\over8}",
    "18-link collar",
    "No order-six flippable-hexagon diagonal potential",
    "**Hostile verdict: PASS.**",
):
    check(marker in audit, f"audit claim/scope marker: {marker}")
check("gravity, or `G`" in audit, "strict gravity/G ceiling")

print(f"PASS__GL6AO_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
