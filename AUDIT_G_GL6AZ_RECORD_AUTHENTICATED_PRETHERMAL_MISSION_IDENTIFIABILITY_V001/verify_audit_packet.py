#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6AZ hostile audit."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(f"FAIL [{checks + 1}] {label}")
    checks += 1


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> dict[str, str]:
    output: dict[str, str] = {}
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        parts = raw.split("  ", 1)
        check(len(parts) == 2, f"two-column row {path.name}:{lineno}")
        digest, relative = parts
        check(bool(re.fullmatch(r"[0-9a-f]{64}", digest)),
              f"digest syntax {path.name}:{lineno}")
        check(relative not in output, f"unique path {path.name}:{lineno}")
        output[relative] = digest
    return output


required = {
    "README.md", "AUDIT.md", "AUDITED_TARGETS.sha256", "VERIFICATION.txt",
    "independent_gl6az_replay.py", "verify_audit_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}
check({path.name for path in HERE.iterdir() if path.is_file()} == required,
      "audit file set is exact")

targets = rows(HERE / "AUDITED_TARGETS.sha256")
check(len(targets) == 13, "thirteen frozen author targets")
for relative, expected in sorted(targets.items()):
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(), f"regular author target: {relative}")
    check(path.parent == AUTHOR, f"author target confined: {relative}")
    check(sha256(path) == expected, f"author target hash: {relative}")

expected_author_names = {
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "IDENTIFIABILITY_LEDGER.json",
    "MANIFEST.sha256", "README.md", "RESULT.md", "SEAL.sha256",
    "SEARCH_LEDGER.md", "SELF_AUDIT.md", "THEOREM.md", "VERIFICATION.txt",
    "verify_packet.py", "verify_record_authenticated_prethermal_mission.py",
}
check({Path(relative).name for relative in targets} == expected_author_names,
      "author target names exact")

author_manifest = rows(AUTHOR / "MANIFEST.sha256")
check(set(author_manifest) == expected_author_names - {"MANIFEST.sha256", "SEAL.sha256"},
      "author manifest covers exact pre-manifest packet")
for relative, expected in sorted(author_manifest.items()):
    check(sha256(AUTHOR / relative) == expected, f"author manifest resolves: {relative}")
author_seal = rows(AUTHOR / "SEAL.sha256")
check(author_seal == {"MANIFEST.sha256": sha256(AUTHOR / "MANIFEST.sha256")},
      "author seal resolves")

dependencies = rows(AUTHOR / "DEPENDENCIES.sha256")
check(len(dependencies) == 33, "thirty-three frozen dependency targets")
for relative, expected in sorted(dependencies.items()):
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(), f"regular dependency: {relative}")
    check(sha256(path) == expected, f"dependency hash: {relative}")
for marker in ("GL6AY", "GL6AM", "GL6V", "GL6AN", "GL5ZZF", "UDCL"):
    check(any(marker in relative for relative in dependencies),
          f"dependency family present: {marker}")

manifest = rows(HERE / "MANIFEST.sha256")
check(set(manifest) == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "audit manifest covers exact pre-manifest packet")
for relative, expected in sorted(manifest.items()):
    check(sha256(HERE / relative) == expected, f"audit manifest resolves: {relative}")
seal = rows(HERE / "SEAL.sha256")
check(seal == {"MANIFEST.sha256": sha256(HERE / "MANIFEST.sha256")},
      "audit seal resolves")

audit = " ".join((HERE / "AUDIT.md").read_text(encoding="utf-8").split()).lower()
verification = " ".join((HERE / "VERIFICATION.txt").read_text(encoding="utf-8").split())
for token in (
    "spurious low-ratio branch",
    "complete source-proof domain",
    "1861.32559690908",
    "18,665,728.0078",
    "grouped strong support cannot cancel",
    "selected ready/`match` factor",
    "binary pair marginal",
    "immediately before the terminal read dilation",
    "every positive `r`",
    "same-parent clock bind",
    "inside-domain branch",
    "outside-domain branch",
    "hostile verdict: pass",
):
    check(token in audit, f"audit scope token: {token}")
for token in (
    "PASS: 519660/519660 GL6AZ constructive checks",
    "PASS: 240/240 GL6AZ packet checks",
    "PASS__INDEPENDENT_GL6AZ_POST_REPAIR_REPLAY__",
    "PASS__GL6AZ_HOSTILE_AUDIT_PACKET__",
):
    check(token in verification, f"verification token: {token}")
for forbidden in (
    "proves gravity",
    "is gravity",
    "derives newton's constant",
    "proves a necessary prethermal threshold",
    "proves physical phase failure",
    "selects one numerical r",
    "bounds the full retained output",
):
    check(forbidden not in audit, f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AZ_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
