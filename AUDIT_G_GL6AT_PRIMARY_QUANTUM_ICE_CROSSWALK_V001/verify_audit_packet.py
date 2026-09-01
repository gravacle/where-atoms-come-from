#!/usr/bin/env python3
"""Fail-closed custody verifier for the distinct GL6AT hostile audit."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR_NAME = "LANE_CROSS_RFT_GRA_GL6AT_PRIMARY_QUANTUM_ICE_CROSSWALK_V001"
AUTHOR = ROOT / AUTHOR_NAME
AUTHOR_FILES = {
    "README.md", "RESULT.md", "PRIMARY_SOURCES.md", "EVIDENCE_LADDER.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "SELF_AUDIT.md",
    "VERIFICATION.txt", "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
AUDIT_FILES = {
    "README.md", "AUDIT.md", "AUDITED_TARGETS.sha256",
    "independent_gl6at_replay.py", "VERIFICATION.txt", "verify_audit_packet.py",
}
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def rows(path: Path):
    parsed = []
    seen = set()
    for raw in path.read_text().splitlines():
        if not raw.strip():
            continue
        expected, relative = raw.split(maxsplit=1)
        relative_path = Path(relative)
        check(len(expected) == 64 and all(c in "0123456789abcdef" for c in expected),
              f"valid digest: {relative}")
        check(not relative_path.is_absolute() and ".." not in relative_path.parts,
              f"safe path: {relative}")
        check(relative not in seen, f"unique row: {relative}")
        target = ROOT / relative_path
        check(target.is_file(), f"target exists: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        seen.add(relative)
        parsed.append((expected, relative))
    return parsed


targets = rows(HERE / "AUDITED_TARGETS.sha256")
check(len(targets) == 11, "exact audited target count")
check({Path(relative).name for _, relative in targets} == AUTHOR_FILES,
      "all and only author files pinned")
check(all(Path(relative).parent == Path(AUTHOR_NAME) for _, relative in targets),
      "audited targets confined to author packet")
check({path.name for path in AUTHOR.iterdir() if path.is_file()} == AUTHOR_FILES,
      "author packet exact inventory")

author_manifest = rows(AUTHOR / "MANIFEST.sha256")
check(len(author_manifest) == 9, "author manifest row count")
check({Path(relative).name for _, relative in author_manifest}
      == AUTHOR_FILES - {"MANIFEST.sha256", "SEAL.sha256"},
      "author manifest exact coverage")
author_seal = rows(AUTHOR / "SEAL.sha256")
check(len(author_seal) == 1 and
      author_seal[0][1] == f"{AUTHOR_NAME}/MANIFEST.sha256",
      "author seal targets author manifest")

dependencies = rows(AUTHOR / "DEPENDENCIES.sha256")
check(len(dependencies) == 12, "exact dependency count")
check(all(("GL6AO" in relative or "GL6AP" in relative)
          for _, relative in dependencies), "dependencies confined to GL6AO/AP")
check(sum("GL6AO" in relative for _, relative in dependencies) == 6,
      "six GL6AO dependency rows")
check(sum("GL6AP" in relative for _, relative in dependencies) == 6,
      "six GL6AP dependency rows")

audit_manifest = rows(HERE / "MANIFEST.sha256")
check(len(audit_manifest) == 6, "audit manifest row count")
check({Path(relative).name for _, relative in audit_manifest} == AUDIT_FILES,
      "audit manifest exact coverage")
check(all(Path(relative).parent == Path(HERE.name) for _, relative in audit_manifest),
      "audit manifest path confinement")
check({path.name for path in HERE.iterdir() if path.is_file()}
      == AUDIT_FILES | {"MANIFEST.sha256", "SEAL.sha256"},
      "audit exact inventory")
audit_seal = rows(HERE / "SEAL.sha256")
check(len(audit_seal) == 1 and
      audit_seal[0][1] == f"{HERE.name}/MANIFEST.sha256",
      "audit seal targets audit manifest")

replay = (HERE / "independent_gl6at_replay.py").read_text()
for forbidden in ("import verify_packet", "import verify_complete", "import verify_locked"):
    check(forbidden not in replay, f"independent replay excludes {forbidden}")

audit = " ".join((HERE / "AUDIT.md").read_text().split())
for marker in (
    "Z_e=2n_e-1=2S_e^z", "with no factor of two", "mu=v/g", "`v/g=0`",
    "published Eq. (55c)/arXiv-v1 Eq. (57c)",
    "published Eq. (70c)/arXiv-v1 Eq. (72c)", "`Omega^4`",
    "cannot be relabeled as the sealed local pair `E`", "unknown higher F3 orders",
    "gravity, or Newton's constant `G` is derived", "**Hostile verdict: PASS.**",
):
    check(marker in audit, f"audit scope marker: {marker}")

verification = (HERE / "VERIFICATION.txt").read_text()
check("PASS__GL6AT_INDEPENDENT_REPLAY__" in verification,
      "independent replay recorded")
check("PASS__GL6AT_HOSTILE_AUDIT_PACKET__" in verification,
      "audit packet replay recorded")
check("NORMAL_AND_OPTIMIZED" in verification, "optimized modes recorded")

print(f"PASS__GL6AT_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
