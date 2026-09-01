#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the frozen GL6AN author packet."""

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
    "verify_native_degree_lock.py",
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
check(len(dependency_paths) == 7, "exact authoritative dependency count")
check(all("GL6AL" not in path for path in dependency_paths),
      "mutable GL6AL is not a proof dependency")
check("LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/THEOREM.md"
      in dependency_paths, "AK theorem pinned")
check("AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/POSTFREEZE_AUDIT.md"
      in dependency_paths, "AK postfreeze audit pinned")
check("AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/SEAL.sha256"
      in dependency_paths, "AK audit seal pinned")

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
        check(f"{HERE.name}/{name}" in manifest_paths, f"manifest coverage: {name}")

seal_lines = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
              if line.strip()]
check(len(seal_lines) == 1, "one seal entry")
expected, relative = seal_lines[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest() == expected,
      "seal hash")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
required_scope = (
    "\\Delta=4U_d(d_\\star-2)",
    "d_\\star=3",
    "generic finite open product",
    "constraint-Gram eigenvalue",
    "singular value itself vanishes linearly",
    "explicit period-four quotient",
    "No global projector in the infinite quasi-local algebra",
    "bi-infinite path",
    "H_{{\\rm eff},\\mathcal Q}^{(4)}",
    "six independently authenticated records",
    "sixth and higher diagonal/loop terms remain unclassified",
)
for token in required_scope:
    check(token in theorem, f"repaired theorem scope: {token}")
for forbidden in (
    "|\\mathbb L|h^2",
    "multi-record collective operation",
    "formation-lineage-to-collective-dynamics",
    "state-dependent diagonal fourth",
    "quadratically soft singular value",
):
    check(forbidden not in theorem + " " + result + " " + readme,
          f"forbidden stale promotion: {forbidden}")
check("all-formed lineage branch" in result, "result retains branch typing")
check("squared singular value" in readme, "README spectral typing")
check("PASS (1056/1056)" in readme, "README verifier count")

print(f"PASS__GL6AN_PACKET__{checks}/{checks}")
