#!/usr/bin/env python3
"""Fail-closed custody, inventory, and scope verifier for GL6AX."""

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
    "verify_all_fixed_order_twist_stability.py", "verify_packet.py",
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
for marker in ("GL6AN", "GL6AO", "GL6AW"):
    check(sum(marker in row for row in dependencies) == 6,
          f"six {marker} author/audit objects")
check(not any("GL6AX" in row for row in dependencies),
      "no circular GL6AX dependency")

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
check("distinct independent hostile audit" in theorem.lower() or
      "distinct independent hostile audit" in readme.lower(),
      "distinct hostile audit required")
for token in (
    "The four augmented columns `(1,d_a)` are affinely independent",
    "Delta N_j=L_j W_j",
    "minimum wrapping Hamming distance is consequently exactly `2L_min`",
    "Fixed-order conservation theorem",
    "raw energy-dependent Feshbach expression",
    "It is no restriction to take each term to commute",
    "`+/-` twist average without a current assumption",
    "Delta_L <= 2pi^2 D_2(L) L0 L2/L1",
    "An exponentially decaying interaction",
    "no finite fixed order of controlled local corrections",
    "The order of limits is fixed",
    "full zero-energy projection",
):
    check(token in theorem, f"theorem scope token: {token}")
check("complex Hermitian amplitudes." in theorem,
      "complex-amplitude scope explicit")
check("minimum winding change has exactly `2L_min` links" in result,
      "result states sharp wrapping threshold")
check("Interaction tails" in self_audit and "Order-of-limits attack" in self_audit,
      "self-audit attacks tails and order of limits")
check("PASS__GL6AX_ALL_FIXED_ORDER_TWIST_STABILITY__" in verification,
      "exact replay recorded")

aggregate = " ".join((theorem, readme, result, self_audit)).lower()
for forbidden in (
    "proves exact finite-coupling port conservation",
    "proves selected-gns gaplessness",
    "proves an isotropic mode",
    "is a physical photon",
    "is a graviton",
    "is gravity",
    "derives the einstein equation",
    "derives newton's constant",
):
    check(forbidden not in aggregate, f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AX_PACKET__{checks}/{checks}")
