#!/usr/bin/env python3
"""Custody, replay, and scope checks for GL6AF."""

import hashlib
import json
import pathlib
import subprocess


root = pathlib.Path(__file__).resolve().parent
repo = root.parent
checks = 0


def require(condition, message):
    global checks
    if not condition:
        raise AssertionError(message)
    checks += 1


for line in (root / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = repo / relative
    require(path.is_file(), f"missing dependency {relative}")
    require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
            f"dependency drift {relative}")

required = {
    "THEOREM.md", "RESULT.md", "README.md", "SELF_AUDIT.md",
    "FORMATION_SOURCE_LEDGER.json", "DEPENDENCIES.sha256",
    "verify_formation_pattern_e2_source.py", "verify_packet.py",
}
require(required.issubset({path.name for path in root.iterdir()}),
        "packet file census")

result = subprocess.run(
    ["python3", str(root / "verify_formation_pattern_e2_source.py")],
    check=True, capture_output=True, text=True,
)
require("PASS GL6AF exact formation-pattern checks" in result.stdout,
        "exact formation-pattern replay")

ledger = json.loads((root / "FORMATION_SOURCE_LEDGER.json").read_text())
require(ledger["rank_by_formed_count"] ==
        {"0": 0, "1": 0, "2": 1, "3": 2, "4": 2},
        "rank ledger")
require(len(ledger["distinct_two_record_covectors"]) == 3,
        "covector ledger")
require("two_formed_E_restriction" in ledger,
        "fixed-E restriction ledger")

theorem = (root / "THEOREM.md").read_text()
for token in (
    "One common physical parent, all branches retained",
    "No event label is inserted",
    "all sixteen retained patterns",
    "0\\hbox{ or }1",
    "two links `a,b` be formed",
    "w_{ab}^{T}w_{ab}",
    "B_E^{(\\boldsymbol\\kappa)}",
    "need not preserve the `E` plane",
    "[\\Pi_{\\boldsymbol\\kappa},U_{\\rm src}(j)]=0",
    "probability of a retained `kappa` branch is independent",
    "not the semantic",
    "authentication-sensitive force",
    "preserving the number of formed records",
    "-{E_\\star^2\\over2\\hbar^2}D",
    "pair observable is itself a record",
    "not collective stiffness",
    "BROKEN_S4_NO_INVARIANT_BLOCK_CLAIM",
    "PHYSICAL_K_NOT_SEMANTIC_REC_DYNAMICS",
    "NO_COLLECTIVE_STIFFNESS_BULK_STRESS_RICCI_GRAVITY_OR_G_CLAIM",
):
    require(token in theorem, f"missing theorem token {token}")

manifest = root / "MANIFEST.sha256"
if manifest.exists():
    for line in manifest.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = repo / relative
        require(path.is_file(), f"missing manifest member {relative}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
                f"manifest drift {relative}")

print(f"PASS GL6AF packet checks {checks}/{checks}")
