#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6AQ hostile audit."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


required = {
    "README.md",
    "AUDIT.md",
    "AUDITED_TARGETS.sha256",
    "VERIFICATION.txt",
    "independent_gl6aq_replay.py",
    "verify_audit_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
}
check({path.name for path in HERE.iterdir() if path.is_file()} == required,
      "audit directory file set exact")


def verify_rows(path, allowed_parent=None):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        check(len(expected) == 64 and all(char in "0123456789abcdef" for char in expected),
              f"hash syntax: {relative}")
        check(relative not in rows, f"unique row: {relative}")
        target = ROOT / relative
        check(target.is_file(), f"target exists: {relative}")
        if allowed_parent is not None:
            check(target.parent == allowed_parent, f"target confined: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        rows.append(relative)
    return rows


target_rows = verify_rows(HERE / "AUDITED_TARGETS.sha256", AUTHOR)
check(len(target_rows) == 10, "ten frozen author targets")
check({Path(row).name for row in target_rows} == {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256",
    "verify_authenticated_e_loop_selection.py", "verify_packet.py",
    "VERIFICATION.txt", "MANIFEST.sha256",
}, "author target set exact")

author_manifest = verify_rows(AUTHOR / "MANIFEST.sha256", AUTHOR)
check(len(author_manifest) == 9, "author manifest has nine rows")
check(set(author_manifest) == set(target_rows) - {f"{AUTHOR.name}/MANIFEST.sha256"},
      "author manifest reconciles with audited targets")

dependencies = verify_rows(AUTHOR / "DEPENDENCIES.sha256")
check(len(dependencies) == 11, "eleven frozen dependency rows")
check(all("GL6AM" in row or "GL6AN" in row for row in dependencies),
      "dependencies confined to GL6AM/GL6AN")
check(any(row.endswith("AUDIT_G_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001/SEAL.sha256")
          for row in dependencies), "GL6AM hostile seal pinned")
check(any(row.endswith("AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/SEAL.sha256")
          for row in dependencies), "GL6AN hostile seal pinned")

manifest_rows = verify_rows(HERE / "MANIFEST.sha256", HERE)
check(len(manifest_rows) == 6, "audit manifest has six rows")
check({Path(row).name for row in manifest_rows} == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "audit manifest file set exact")

seal_rows = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
             if line.strip()]
check(len(seal_rows) == 1, "one audit seal row")
expected, relative = seal_rows[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets audit manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest() == expected,
      "audit seal hash")

audit = " ".join((HERE / "AUDIT.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
for token in (
    "C|_E=(8/3) I_E",
    "P_Q X_e P_Q = 0",
    "P_E K(t)=0",
    "||delta M||^2=16",
    "-(63/8) h^6/U_d^5",
    "product_{e in C} kappa_e",
    "Query/source confusion",
    "Finite/bulk confusion",
    "Universal/existential confusion",
    "not a locked-state witness",
    "universal nonzero claim is refuted",
    "Hostile verdict: PASS",
):
    check(token in audit, f"audit scope token: {token}")
for token in (
    "PASS__INDEPENDENT_GL6AQ_HOSTILE_REPLAY__307976/307976",
    "PASS__GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION__1010/1010",
    "PASS__GL6AQ_PACKET__95/95",
    "PASS__INDEPENDENT_GL6AM_HOSTILE_REPLAY__1476/1476",
    "PASS__GL6AM_HOSTILE_AUDIT_PACKET__54/54",
    "PASS__INDEPENDENT_GL6AN_HOSTILE_REPLAY__3939/3939",
    "PASS__GL6AN_HOSTILE_AUDIT_PACKET__58/58",
    "PASS__GL6AQ_AUDIT_PACKET__",
):
    check(token in verification, f"verification token: {token}")

for forbidden in (
    "the pair query is the physical k source",
    "the product trace is a locked state",
    "no stationary locked state exists",
    "the loop is a graviton",
    "proves gravity",
    "derives g",
):
    check(forbidden not in audit.lower(), f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AQ_AUDIT_PACKET__{checks}/{checks}")
