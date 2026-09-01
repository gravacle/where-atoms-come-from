#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6AU hostile audit."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AU_VG0_STATIC_STRUCTURE_PHASE_CLOSURE_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(label)
    checks += 1


required = {
    "README.md", "AUDIT.md", "AUDITED_TARGETS.sha256", "VERIFICATION.txt",
    "independent_gl6au_replay.py", "verify_audit_packet.py",
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
check(len(targets) == 12, "twelve frozen author targets")
check({Path(row).name for row in targets} == {
    "README.md", "THEOREM.md", "RESULT.md", "EVIDENCE_AND_HYPOTHESIS.md",
    "SELF_AUDIT.md", "DEPENDENCIES.md", "DEPENDENCIES.sha256",
    "VERIFICATION.txt", "verify_vg0_static_structure_closure.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}, "author target set exact")

author_manifest = verify_rows(AUTHOR / "MANIFEST.sha256", AUTHOR)
check(len(author_manifest) == 10, "author manifest has ten rows")
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
check(len(dependencies) == 32, "thirty-two dependencies")
for marker in ("GL6AO", "GL6AR", "GL6AS"):
    check(sum(marker in row for row in dependencies) == 6,
          f"six {marker} custody objects")
check(sum("GL6AT" in row for row in dependencies) == 14,
      "fourteen GL6AT author/audit objects")
check(sum("AUDIT_G_GL6AT" in row for row in dependencies) == 3,
      "three distinct GL6AT audit objects")

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
    "exact transverse plane",
    "same connected component",
    "threshold is strictly `alpha<2`",
    "PF transform",
    "low-character variance of the nonuniform PF law",
    "does not prove isotropic finite-size scaling",
    "anisotropic Følner",
    "separate GNS bridge",
    "named hypotheses",
    "No pole, physical cone, stress/Ricci law, gravity, or `G` is derived",
    "Hostile verdict: PASS",
):
    check(token.lower() in audit, f"audit scope token: {token}")

for token in (
    "PASS__GL6AU_VG0_STATIC_STRUCTURE__8111/8111",
    "PASS__GL6AU_PACKET__163/163",
    "PASS__INDEPENDENT_GL6AU_HOSTILE_REPLAY__41046/41046",
    "Post-freeze extension",
):
    check(token in verification, f"verification token: {token}")

for forbidden in (
    "ice0-static is proved", "ice0-gns-bridge is proved",
    "the v/g=0 phase is rigorously proved", "is a physical photon",
    "is a graviton", "is gravity", "derives the einstein equation",
    "derives newton's constant",
):
    check(forbidden not in audit, f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AU_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
