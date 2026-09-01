#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6AP author packet."""

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
    "verify_locked_ir_conditions.py",
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
check(len(dependency_paths) == 7, "exact dependency count")
check(all("GL6AN" in path for path in dependency_paths),
      "only GL6AN author/audit custody is imported")
check(all(marker not in path for path in dependency_paths
          for marker in ("GL6AL", "GL6AO", "GL6AQ")),
      "mutable or parallel lanes excluded")
check("LANE_CROSS_RFT_GRA_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/THEOREM.md"
      in dependency_paths, "GL6AN theorem pinned")
check("AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/AUDIT.md"
      in dependency_paths, "GL6AN hostile audit pinned")

manifest_paths = set()
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    relative_path = Path(relative)
    check(relative_path.parent == Path(HERE.name),
          f"manifest target stays in packet: {relative}")
    target = ROOT / relative_path
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
readme = " ".join((HERE / "README.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())

for document in (theorem, readme):
    check("author frozen and sealed" in document, "frozen author status")
check("independent hostile audit" in theorem.lower(),
      "theorem requires independent hostile audit")
check("PASS__GL6AP_LOCKED_IR_CONDITIONS__434/434" in verification,
      "exact replay count recorded")

required_scope = (
    "Hom}_{S_4}(T_2,E)=0",
    "(N_1,N_2,N_3)_i=(1,2,3)",
    "(N_1,N_2,N_3)_f=(3,2,1)",
    "checks all 24 graph automorphisms",
    "no nonzero uniform linear combination of pair-`E` is conserved",
    "scalar mass `r_E` is symmetry allowed",
    "exactly two invariant spatial quadratic contractions",
    "has a discrete character set and by itself has no literal",
    "does not insert a global locked projector",
    "a pole sequence can mathematically have residues tending to zero",
    "removes any elastic zero-frequency atom",
    "no all-orders effective Hamiltonian is being asserted here",
    "physical embedding `X:A3 -> R^3`",
    "Nothing in GL6AP assumes or derives a gauge phase",
)
for token in required_scope:
    check(token in theorem, f"theorem scope token: {token}")
check("fixed `Q4` alone has no infrared sequence" in readme,
      "README finite-volume gate")
check("neither imported nor rederived within this GL6AN-only lane"
      in self_audit, "self-audit parallel-lane independence")
check("positive-frequency weight" in result, "result spectral-weight scope")

aggregate = " ".join((theorem, readme, result, self_audit)).lower()
for forbidden in (
    "the e-sector pole is proved",
    "character is physical momentum",
    "is a gauge theory",
    "is a photon",
    "is a graviton",
    "is gravity",
    "derives newton's constant",
    "proves a common cone",
    "[h_{\\rm eff},\\mathbf n_e]\\ne0",
):
    check(forbidden not in aggregate, f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AP_PACKET__{checks}/{checks}")
