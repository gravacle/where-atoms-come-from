#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6AO author packet."""

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
    "THEOREM.md",
    "RESULT.md",
    "SELF_AUDIT.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "VERIFICATION.txt",
    "verify_sixth_order_hamiltonian.py",
    "verify_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
)
for name in required:
    check((HERE / name).is_file(), f"required file: {name}")

dependency_paths = set()
for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"dependency exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"dependency hash: {relative}")
    dependency_paths.add(relative)
check(len(dependency_paths) == 6, "exact dependency count")
check(all("GL6AL" not in path for path in dependency_paths),
      "mutable GL6AL is excluded")
check(all("GL6AN" in path for path in dependency_paths),
      "only GL6AN author/audit custody is imported")
check("LANE_CROSS_RFT_GRA_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/THEOREM.md"
      in dependency_paths, "GL6AN theorem pinned")
check("AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/AUDIT.md"
      in dependency_paths, "GL6AN hostile audit pinned")

manifest_paths = set()
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"manifest target exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"manifest hash: {relative}")
    manifest_paths.add(relative)
for name in required:
    if name not in ("MANIFEST.sha256", "SEAL.sha256"):
        check(f"{HERE.name}/{name}" in manifest_paths,
              f"manifest coverage: {name}")
check(len(manifest_paths) == 9, "exact author manifest count")

seal_lines = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
              if line.strip()]
check(len(seal_lines) == 1, "one seal row")
expected, relative = seal_lines[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest()
      == expected, "seal hash")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
required_scope = (
    "K_6 = T_6 - b X_4 + b^2 A_3 - d A_2",
    "matrix element through order six is an alternating",
    "<s'|K_6|s> = -63/8",
    "<s|K_6|s> = -(893/1080) M",
    "no order-six diagonal potential proportional to the number of flippable six-cycles",
    "M^3` and `M^2` terms cancel exactly",
    "No global locked projector is inserted into the infinite quasi-local algebra",
    "Phi_c^(6) = -(63/8)(h^6/U_d^5) tau_c",
    "does **not** prove",
    "gravity, or `G`",
)
for token in required_scope:
    check(token in theorem, f"theorem scope token: {token}")
check("ONLY_ALTERNATING_HEXAGON" in verification,
      "verification records offdiagonal classification")
check("H6_DIAG=-893M/1080" in verification,
      "verification records diagonal coefficient")
check("PASS__GL6AO_COMPLETE_SIXTH_ORDER__120304/120304" in verification,
      "physics replay count recorded")
check("PASS__GL6AO_PACKET__82/82" in verification,
      "packet replay count recorded")
check("no flippable-cycle-count potential" in readme,
      "README states the absent diagonal potential")
check("formal uniformly finite-range linked interaction" in result,
      "result keeps thermodynamic claim formal")
check("not called a photon, graviton, or gravity" in self_audit,
      "self-audit blocks conventional promotion")

for forbidden in (
    "is a graviton",
    "is gravity",
    "derives Newton's constant",
    "proves a gapless pole",
    "proves the Einstein equation",
    "configuration-dependent diagonal term at order six exists",
):
    check(forbidden not in theorem + " " + result + " " + readme,
          f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AO_PACKET__{checks}/{checks}")
