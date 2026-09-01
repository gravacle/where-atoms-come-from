#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6AS author packet."""

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
    "verify_native_hexagon_collective.py",
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
check(len(dependency_paths) == 17, "exact dependency count")
check(sum("GL6AO" in path for path in dependency_paths) == 6,
      "six AO custody objects")
check(sum("GL6AP" in path for path in dependency_paths) == 6,
      "six AP custody objects")
check(sum("GL6AQ" in path for path in dependency_paths) == 5,
      "five AQ custody objects")
check(not any(
    path == "LANE_CROSS_RFT_GRA_GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION_V001/SEAL.sha256"
    for path in dependency_paths
), "missing AQ author seal is not invented")
check(all(marker not in path for path in dependency_paths
          for marker in ("GL6AL", "GL6AR")),
      "mutable and unrelated lanes excluded")

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

seal_lines = [
    line for line in (HERE / "SEAL.sha256").read_text().splitlines()
    if line.strip()
]
check(len(seal_lines) == 1, "one seal row")
expected, relative = seal_lines[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest()
      == expected, "seal hash")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
dependencies = " ".join((HERE / "DEPENDENCIES.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())

for document in (theorem, readme):
    check("author frozen and sealed" in document, "frozen author status")
check("independent hostile audit" in theorem.lower(),
      "theorem requires hostile audit")
check("PASS__GL6AS_NATIVE_HEXAGON_COLLECTIVE__1704/1704" in verification,
      "physics replay count recorded")
check("deliberately has no author `SEAL.sha256`" in dependencies,
      "AQ seal absence documented")

required_scope = (
    "[H_{\\rm hex},N_a]=0",
    "\\operatorname{im}C(\\chi)=\\ker B(\\chi)",
    "C_1(\\theta)C_1(\\theta)^T",
    "conditional diagnostic, not a dispersion theorem",
    "\\omega_1^2=\\omega_2^2",
    "cannot be set equal to the bare",
    "exact quadratic composite",
    "\\langle {\\rm one}\\text{-}T_2|O(c)|0\\rangle=0",
    "not calibrated physical momenta",
    "support thresholds only if the form factor",
    "K_{E\\leftarrow\\mathrm{loop}}(\\omega,1)=0",
    "\\sum_c\\tau_c=-H_{\\rm hex}/J",
    "leading Fock overlap is therefore a two-mode channel",
    "\\operatorname{Sym}^2_0(T_2)=E\\oplus T_2",
    "normalized product trace remains an exact universal counterexample",
    "Nothing here assumes a conventional gauge phase",
)
for token in required_scope:
    check(token in theorem, f"theorem scope token: {token}")
check("No native cone or `E` pole follows" in result,
      "result begins with no-go")
check("Fixed `Q4` is used for exact microscopic checks only" in readme,
      "README finite-volume gate")
check("The always-conserved `A1` energy density" in self_audit,
      "self-audit limits only-candidate quantifier")

aggregate = " ".join((theorem, readme, result, self_audit)).lower()
for forbidden in (
    "the t2 mode is proved gapless",
    "the e pole is proved",
    "is a gauge theory",
    "is a photon",
    "is a graviton",
    "is gravity",
    "derives a physical cone",
    "derives newton's constant",
    "fixed q4 has an infrared limit",
    "g=j",
):
    check(forbidden not in aggregate, f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AS_PACKET__{checks}/{checks}")
