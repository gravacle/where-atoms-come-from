#!/usr/bin/env python3
"""Frozen-packet custody, integrity, and scope checks for GL6AG."""

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
checks = 0


def require(condition, message):
    global checks
    if not condition:
        raise AssertionError(message)
    checks += 1


required = (
    "THEOREM.md",
    "README.md",
    "RESULT.md",
    "SELF_AUDIT.md",
    "PRESCREEN_REQUEST.md",
    "DEPENDENCIES.sha256",
    "EXACT_MATCHED_LEDGER.json",
    "verify_structure_and_ledger.py",
    "verify_n1_matched_formation_propagation.cpp",
    "VERIFICATION.txt",
    "verify_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
)
for name in required:
    path = HERE / name
    require(path.is_file() and path.stat().st_size > 0, f"required file: {name}")

require(not (HERE / "AUDIT.md").exists(), "post-freeze audit remains separate")

theorem = (HERE / "THEOREM.md").read_text()
result = (HERE / "RESULT.md").read_text()
request = (HERE / "PRESCREEN_REQUEST.md").read_text()
verification = (HERE / "VERIFICATION.txt").read_text()
source = (HERE / "verify_n1_matched_formation_propagation.cpp").read_text()
ledger = json.loads((HERE / "EXACT_MATCHED_LEDGER.json").read_text())

for token in (
    "author frozen after independent hostile pre-freeze review",
    "matched contrast",
    "w_{ab}:=E_{(ab),:}^{T}\\in\\mathbb R^2",
    "nonzero exactly when",
    "term-ablation",
    "not presented\nas an authenticated physical intervention",
    "not the semantic terminal predicate",
    "No statement about the absolute receiver mean",
):
    require(token in theorem, f"theorem token: {token}")

require("INDEPENDENT_HOSTILE_PRESCREEN_CLEAN__SOURCE_FREEZE_AUTHORIZED_AND_COMPLETED"
        in request, "clean prescreen and freeze status")
require("Do not edit the mutable builder files" in request,
        "independent-review boundary")
require("does not claim an authenticated controller" in result,
        "term-ablation scope in result")
require(ledger["matched_intervention"]["absolute_neighbor_formula_used"] is False,
        "ledger excludes absolute-neighbor formula")
require(ledger["bridge_off"]["authenticated_physical_switch_claimed"] is False,
        "ledger excludes physical switch")

for token in (
    "constexpr int kLinks = 16",
    "constexpr int kDimension = 1 << kLinks",
    "Integer(63371264)",
    "Integer(123422773248)",
    "for (int mask = 0; mask < 16; ++mask)",
    "shared_bridges ? 30 : 24",
):
    require(token in source, f"exact replay token: {token}")

require("PASS GL6AG structure/ledger checks 144/144" in verification,
        "fast verification recorded")
require("PASS GL6AG exact matched-formation checks 3955/3955" in verification,
        "full verification recorded")
require("PASS GL6AG frozen packet checks 62/62" in verification,
        "packet verification recorded")
require("post-\nfreeze audit remains required" in verification,
        "post-freeze audit boundary")

manifest_lines = (HERE / "MANIFEST.sha256").read_text().splitlines()
require(len(manifest_lines) == 11, "manifest member census")
for line in manifest_lines:
    expected, relative = line.split("  ", 1)
    path = HERE.parent / relative
    require(path.is_file(), f"manifest member exists: {relative}")
    require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
            f"manifest member hash: {relative}")

seal_hash, seal_relative = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
require(seal_relative ==
        f"{HERE.name}/MANIFEST.sha256", "seal target")
require(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest() ==
        seal_hash, "seal hash")

print(f"PASS GL6AG frozen packet checks {checks}/{checks}")
