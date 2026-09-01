#!/usr/bin/env python3
"""Fail-closed custody verification for the distinct GL6AP audit packet."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR_NAME = "LANE_CROSS_RFT_GRA_GL6AP_LOCKED_IR_CONDITIONAL_RESPONSE_V001"
AUTHOR = ROOT / AUTHOR_NAME
AUDIT_FILES = {
    "README.md", "AUDIT.md", "AUDITED_TARGETS.sha256",
    "independent_gl6ap_replay.py", "VERIFICATION.txt", "verify_audit_packet.py",
}
AUTHOR_FILES = {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md", "DEPENDENCIES.md",
    "DEPENDENCIES.sha256", "VERIFICATION.txt", "verify_locked_ir_conditions.py",
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
    seen = set()
    for raw in path.read_text().splitlines():
        if not raw.strip():
            continue
        expected, relative = raw.split(maxsplit=1)
        relative_path = Path(relative)
        check(len(expected) == 64 and all(c in "0123456789abcdef" for c in expected),
              f"valid digest syntax: {relative}")
        check(not relative_path.is_absolute() and ".." not in relative_path.parts,
              f"safe relative path: {relative}")
        check(relative not in seen, f"unique row: {relative}")
        target = ROOT / relative_path
        check(target.is_file(), f"target exists: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        seen.add(relative)
        parsed.append((expected, relative))
    return parsed


targets = rows(HERE / "AUDITED_TARGETS.sha256")
check(len(targets) == 11, "exact audited target row count")
check({Path(relative).name for _, relative in targets} == AUTHOR_FILES,
      "audited targets cover exact author packet")
check(all(Path(relative).parent == Path(AUTHOR_NAME) for _, relative in targets),
      "audited targets confined to author packet")
check({path.name for path in AUTHOR.iterdir() if path.is_file()} == AUTHOR_FILES,
      "author packet has no undeclared files")

author_manifest = rows(AUTHOR / "MANIFEST.sha256")
check(len(author_manifest) == 9, "author manifest exact row count")
check({Path(relative).name for _, relative in author_manifest}
      == AUTHOR_FILES - {"MANIFEST.sha256", "SEAL.sha256"},
      "author manifest exact coverage")
check(all(Path(relative).parent == Path(AUTHOR_NAME)
          for _, relative in author_manifest), "author manifest path confinement")
author_seal = rows(AUTHOR / "SEAL.sha256")
check(len(author_seal) == 1
      and author_seal[0][1] == f"{AUTHOR_NAME}/MANIFEST.sha256",
      "author seal targets author manifest")

dependencies = rows(AUTHOR / "DEPENDENCIES.sha256")
check(len(dependencies) == 7, "exact author dependency count")
check(all("GL6AN" in relative for _, relative in dependencies),
      "dependencies confined to GL6AN author/audit custody")
check(all(marker not in relative for _, relative in dependencies
          for marker in ("GL6AL", "GL6AO", "GL6AQ")),
      "mutable and parallel lanes excluded")
check(any(relative.endswith("GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/THEOREM.md")
          for _, relative in dependencies), "GL6AN theorem pinned")
check(any(relative.endswith("AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/AUDIT.md")
          for _, relative in dependencies), "GL6AN hostile audit pinned")

audit_manifest = rows(HERE / "MANIFEST.sha256")
check(len(audit_manifest) == 6, "audit manifest exact row count")
check({Path(relative).name for _, relative in audit_manifest} == AUDIT_FILES,
      "audit manifest exact coverage")
check(all(Path(relative).parent == Path(HERE.name) for _, relative in audit_manifest),
      "audit manifest path confinement")
check({path.name for path in HERE.iterdir() if path.is_file()}
      == AUDIT_FILES | {"MANIFEST.sha256", "SEAL.sha256"},
      "audit directory has no undeclared files")
audit_seal = rows(HERE / "SEAL.sha256")
check(len(audit_seal) == 1
      and audit_seal[0][1] == f"{HERE.name}/MANIFEST.sha256",
      "audit seal targets audit manifest")

replay = (HERE / "independent_gl6ap_replay.py").read_text()
for forbidden_import in (
    "import verify_locked_ir_conditions",
    "import verify_native_degree_lock",
):
    check(forbidden_import not in replay, "independent replay imports no author verifier")

audit = (HERE / "AUDIT.md").read_text()
for marker in (
    "\\operatorname{Hom}_{S_4}(T_2,E)=0",
    "twelve `+1` and twelve `-1`",
    "(N_1,N_2,N_3):(1,2,3)\\longmapsto(3,2,1)",
    "-{63\\over8}",
    "exactly two quadratic spatial covariants",
    "a mass is symmetry allowed",
    "restriction to `(0,infinity)` is therefore necessary",
    "Physical momentum additionally requires",
    "**Hostile verdict: PASS.**",
):
    check(marker in audit, f"audit claim/scope marker: {marker}")
check("gravity, or `G` is derived" in audit, "strict gravity/G ceiling")

verification = (HERE / "VERIFICATION.txt").read_text()
check("PASS__GL6AP_INDEPENDENT_REPLAY__" in verification,
      "independent replay recorded")
check("PASS__GL6AP_HOSTILE_AUDIT_PACKET__" in verification,
      "audit packet replay recorded")
check("NORMAL_AND_OPTIMIZED" in verification, "optimized-mode replay recorded")

print(f"PASS__GL6AP_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
