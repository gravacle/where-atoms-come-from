#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6AX hostile audit."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AX_ALL_FIXED_ORDER_TWIST_STABILITY_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(label)
    checks += 1


required = {
    "README.md", "AUDIT.md", "AUDITED_TARGETS.sha256", "VERIFICATION.txt",
    "independent_gl6ax_replay.py", "verify_audit_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}
check({path.name for path in HERE.iterdir() if path.is_file()} == required,
      "audit directory file set exact")


def verify_rows(path: Path, allowed_parent: Path | None = None):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        check(len(expected) == 64 and
              all(char in "0123456789abcdef" for char in expected),
              f"hash syntax: {relative}")
        check(relative not in rows, f"unique row: {relative}")
        target = ROOT / relative
        check(target.is_file() and not target.is_symlink(),
              f"regular target: {relative}")
        if allowed_parent is not None:
            check(target.parent == allowed_parent, f"target confined: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        rows.append(relative)
    return rows


targets = verify_rows(HERE / "AUDITED_TARGETS.sha256", AUTHOR)
check(len(targets) == 11, "eleven frozen author targets")
check({Path(row).name for row in targets} == {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "VERIFICATION.txt",
    "verify_all_fixed_order_twist_stability.py", "verify_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}, "author target set exact")

author_manifest = verify_rows(AUTHOR / "MANIFEST.sha256", AUTHOR)
check(len(author_manifest) == 9, "author manifest has nine rows")
check(set(author_manifest) == set(targets) - {
    f"{AUTHOR.name}/MANIFEST.sha256", f"{AUTHOR.name}/SEAL.sha256"
}, "author manifest reconciles with frozen targets")
author_seal = [line for line in (AUTHOR / "SEAL.sha256").read_text().splitlines()
               if line.strip()]
check(len(author_seal) == 1, "one author seal row")
expected, relative = author_seal[0].split(maxsplit=1)
check(relative == f"{AUTHOR.name}/MANIFEST.sha256", "author seal target")
check(hashlib.sha256((AUTHOR / "MANIFEST.sha256").read_bytes()).hexdigest()
      == expected, "author seal resolves")

dependencies = verify_rows(AUTHOR / "DEPENDENCIES.sha256")
check(len(dependencies) == 18, "eighteen frozen dependencies")
for marker in ("GL6AN", "GL6AO", "GL6AW"):
    check(sum(marker in row for row in dependencies) == 6,
          f"six {marker} author/audit custody objects")
    check(sum(f"AUDIT_G_{marker}" in row for row in dependencies) == 3,
          f"three distinct {marker} audit objects")
for seal_relative in (row for row in dependencies if row.endswith("/SEAL.sha256")):
    seal_path = ROOT / seal_relative
    rows = [line for line in seal_path.read_text().splitlines() if line.strip()]
    check(len(rows) == 1, f"one dependency seal row: {seal_relative}")
    expected, manifest_relative = rows[0].split(maxsplit=1)
    check(manifest_relative.endswith("/MANIFEST.sha256"),
          f"dependency seal targets a manifest: {seal_relative}")
    manifest_path = ROOT / manifest_relative
    check(manifest_path.is_file() and not manifest_path.is_symlink(),
          f"dependency manifest is regular: {manifest_relative}")
    check(hashlib.sha256(manifest_path.read_bytes()).hexdigest() == expected,
          f"dependency seal resolves: {seal_relative}")

manifest = verify_rows(HERE / "MANIFEST.sha256", HERE)
check(len(manifest) == 6, "audit manifest has six rows")
check({Path(row).name for row in manifest} == required - {
    "MANIFEST.sha256", "SEAL.sha256"
}, "audit manifest file set exact")
seal = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
        if line.strip()]
check(len(seal) == 1, "one audit seal row")
expected, relative = seal[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "audit seal target")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest()
      == expected, "audit seal resolves")

audit = " ".join((HERE / "AUDIT.md").read_text().split()).lower()
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
for token in (
    "unimodular affine simplex",
    "degree-four vertices, self-touching even subgraphs",
    "sharp minimum is thus exactly `2l_min`",
    "canonical physical-basis word representation",
    "termwise `u(1)` symmetrization",
    "no seam or origin term was omitted",
    "all odd terms, including a nonzero equilibrium-current expectation",
    "the coefficient and factor of two are correct",
    "not a missing factor four",
    "the order is fixed before the volume limit",
    "compatible selected-state/gns bridge remains open",
    "hostile verdict: pass",
):
    check(token in audit, f"audit scope token: {token}")

for token in (
    "PASS__GL6AX_ALL_FIXED_ORDER_TWIST_STABILITY__542791/542791",
    "PASS__GL6AX_PACKET__155/155",
    "PASS__INDEPENDENT_GL6AX_HOSTILE_REPLAY__35867/35867",
    "PASS__GL6AX_HOSTILE_AUDIT_PACKET__270/270",
):
    check(token in verification, f"verification token: {token}")

for forbidden in (
    "proves exact finite-coupling port conservation",
    "proves centered-sector selection",
    "proves selected-gns gaplessness",
    "proves an isotropic cone",
    "is a physical photon",
    "is a graviton",
    "is gravity",
    "derives newton's constant",
):
    check(forbidden not in audit, f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AX_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
