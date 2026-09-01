#!/usr/bin/env python3
"""Custody and scope checks for the GL6AA packet."""

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
    "THEOREM.md",
    "RESULT.md",
    "README.md",
    "SELF_AUDIT.md",
    "ATLAS_LEDGER.json",
    "DEPENDENCIES.sha256",
    "verify_record_authenticated_atlas.py",
    "verify_packet.py",
}
require(required.issubset({p.name for p in root.iterdir()}),
        "packet file census")

result = subprocess.run(
    ["python3", str(root / "verify_record_authenticated_atlas.py")],
    check=True,
    capture_output=True,
    text=True,
)
require("PASS GL6AA exact checks 1686208/1686208" in result.stdout,
        "exact atlas replay")
json.loads((root / "ATLAS_LEDGER.json").read_text())
require(True, "ledger JSON")

manifest = root / "MANIFEST.sha256"
if manifest.exists():
    for line in manifest.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = repo / relative
        require(path.is_file(), f"missing manifest member {relative}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
                f"manifest drift {relative}")

seal = root / "SEAL.sha256"
if seal.exists():
    expected, relative = seal.read_text().strip().split("  ", 1)
    path = repo / relative
    require(path == manifest, "seal must pin this packet manifest")
    require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
            "manifest seal drift")

theorem = (root / "THEOREM.md").read_text()
for token in (
    "response bytes alone",
    "source/controller path",
    "conditional finite derivative",
    "flags need not be mutually",
    "`K_e=1` **iff**",
    "queried injectivity of the independent physical-site source `sigma`",
    "postselected on `MATCH`",
    "record-authenticated relational incidence atlas",
    "wavevector remain calibration data",
    "Ricci ansatz",
):
    require(token in theorem, f"missing scope token {token}")

print(f"PASS GL6AA packet checks {checks}/{checks}")
