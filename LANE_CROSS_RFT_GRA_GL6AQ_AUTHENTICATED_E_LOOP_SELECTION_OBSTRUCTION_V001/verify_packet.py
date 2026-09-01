#!/usr/bin/env python3
"""Fail-closed packet and dependency verification for GL6AQ V001."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
EXPECTED = {
    "README.md",
    "THEOREM.md",
    "RESULT.md",
    "SELF_AUDIT.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "verify_authenticated_e_loop_selection.py",
    "verify_packet.py",
    "VERIFICATION.txt",
}
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def verify_hash_rows(path: Path, allowed_root: Path | None = None):
    rows = []
    for raw in path.read_text().splitlines():
        if not raw.strip():
            continue
        expected, relative = raw.split(maxsplit=1)
        check(len(expected) == 64 and all(c in "0123456789abcdef" for c in expected),
              f"valid SHA-256 syntax: {relative}")
        check(relative not in rows, f"unique hash row: {relative}")
        target = ROOT / relative
        check(target.is_file(), f"hash target exists: {relative}")
        if allowed_root is not None:
            check(target.parent == allowed_root, f"manifest target remains inside packet: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        rows.append(relative)
    return rows


manifest_rows = verify_hash_rows(HERE / "MANIFEST.sha256", HERE)
check({Path(row).name for row in manifest_rows} == EXPECTED, "manifest file set is exact")
check(len(manifest_rows) == len(EXPECTED), "manifest row count is exact")
check({path.name for path in HERE.iterdir() if path.is_file()} == EXPECTED | {"MANIFEST.sha256"},
      "packet directory has no undeclared files")

dependency_rows = verify_hash_rows(HERE / "DEPENDENCIES.sha256")
check(len(dependency_rows) == 11, "dependency row count is exact")
check(all("GL6AM" in row or "GL6AN" in row for row in dependency_rows),
      "dependencies are confined to sealed GL6AM/GL6AN custody")

check(not (HERE / "SEAL.sha256").exists(), "author packet does not self-seal")
print(f"PASS__GL6AQ_PACKET__{checks}/{checks}")
