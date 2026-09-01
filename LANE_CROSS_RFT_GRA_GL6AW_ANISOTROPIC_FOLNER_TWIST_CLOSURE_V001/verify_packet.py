#!/usr/bin/env python3
"""Fail-closed custody, inventory, and scope verifier for GL6AW."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(label)
    checks += 1


required = {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "VERIFICATION.txt",
    "verify_anisotropic_folner_twist.py", "verify_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}
check({path.name for path in HERE.iterdir() if path.is_file()} == required,
      "author file set exact")


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


dependencies = verify_rows(HERE / "DEPENDENCIES.sha256")
check(len(dependencies) == 18, "eighteen exact dependencies")
for marker in ("GL6AR", "GL6AS", "GL6AU"):
    check(sum(marker in row for row in dependencies) == 6,
          f"six {marker} author/audit objects")
check(not any("GL6AV" in row or "GL6AW" in row for row in dependencies),
      "no mutable or circular lane dependency")

manifest = verify_rows(HERE / "MANIFEST.sha256", HERE)
check(len(manifest) == 9, "author manifest has nine rows")
check({Path(row).name for row in manifest} == required - {
    "MANIFEST.sha256", "SEAL.sha256"
}, "author manifest file set exact")
seal = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
        if line.strip()]
check(len(seal) == 1, "one author seal row")
expected, relative = seal[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest()
      == expected, "seal resolves")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())

for document in (theorem, readme):
    check("author frozen and sealed" in document, "frozen author status")
check("independent hostile audit" in theorem.lower() or
      "independent hostile audit" in readme.lower(),
      "distinct hostile audit required")
for token in (
    "The centered sector is nonempty",
    "Y U_0 Y^-1=exp(-2pi i N_0/L1) U_0",
    "exp(-pi i L0 L2)=-1",
    "The twisted trial is therefore normalized, lies in the same connected component",
    "exactly two orientations per cell",
    "Delta_C <=2JV[1-cos(2pi/L1)]",
    "Ground-component dichotomy",
    "(L0,L1,L2)=(m,2m^3,m)",
    "Delta_C(m)<=2pi^2 J/m",
    "not a derived small-character collective dispersion",
    "full zero-energy projection",
    "not a selected-GNS mode",
):
    check(token in theorem, f"theorem scope token: {token}")
check("bypass" in result.lower() and "finite-size anisotropic closure" in result,
      "result states exact bypass scope")
check("Term count" in self_audit and "Følner scope" in self_audit,
      "self-audit attacks counting and geometry")
check("PASS__GL6AW_ANISOTROPIC_FOLNER_TWIST__" in verification,
      "exact replay recorded")

aggregate = " ".join((theorem, readme, result, self_audit)).lower()
for forbidden in (
    "proves isotropic gaplessness", "proves gns gaplessness",
    "is a physical photon", "is a graviton", "is gravity",
    "derives a physical cone", "derives the einstein equation",
    "derives newton's constant",
):
    check(forbidden not in aggregate, f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AW_PACKET__{checks}/{checks}")
