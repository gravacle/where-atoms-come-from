#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6AS hostile audit."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AS_NATIVE_HEXAGON_COLLECTIVE_RESPONSE_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


required = {
    "README.md", "AUDIT.md", "AUDITED_TARGETS.sha256", "VERIFICATION.txt",
    "independent_gl6as_replay.py", "verify_audit_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}
check({path.name for path in HERE.iterdir() if path.is_file()} == required,
      "audit directory file set exact")


def verify_rows(path, allowed_parent=None):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        check(len(expected) == 64 and all(char in "0123456789abcdef" for char in expected),
              f"hash syntax: {relative}")
        check(relative not in rows, f"unique row: {relative}")
        target = ROOT / relative
        check(target.is_file() and not target.is_symlink(), f"regular target: {relative}")
        if allowed_parent is not None:
            check(target.parent == allowed_parent, f"target confined: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        rows.append(relative)
    return rows


targets = verify_rows(HERE / "AUDITED_TARGETS.sha256", AUTHOR)
check(len(targets) == 11, "eleven frozen author targets")
check({Path(row).name for row in targets} == {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "VERIFICATION.txt",
    "verify_native_hexagon_collective.py", "verify_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}, "author target file set exact")

author_manifest = verify_rows(AUTHOR / "MANIFEST.sha256", AUTHOR)
check(len(author_manifest) == 9, "author manifest has nine rows")
check(set(author_manifest) == set(targets) - {
    f"{AUTHOR.name}/MANIFEST.sha256", f"{AUTHOR.name}/SEAL.sha256"
}, "author manifest reconciles with audited targets")
author_seal = [line for line in (AUTHOR / "SEAL.sha256").read_text().splitlines()
               if line.strip()]
check(len(author_seal) == 1, "one author seal row")
expected, relative = author_seal[0].split(maxsplit=1)
check(relative == f"{AUTHOR.name}/MANIFEST.sha256", "author seal targets manifest")
check(hashlib.sha256((AUTHOR / "MANIFEST.sha256").read_bytes()).hexdigest()
      == expected, "author seal resolves")

dependencies = verify_rows(AUTHOR / "DEPENDENCIES.sha256")
check(len(dependencies) == 17, "seventeen dependency rows")
check(sum("GL6AO" in row for row in dependencies) == 6, "six AO objects")
check(sum("GL6AP" in row for row in dependencies) == 6, "six AP objects")
check(sum("GL6AQ" in row for row in dependencies) == 5, "five AQ objects")
check(not any(row.endswith(
    "LANE_CROSS_RFT_GRA_GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION_V001/SEAL.sha256"
) for row in dependencies), "missing AQ author seal not invented")

manifest = verify_rows(HERE / "MANIFEST.sha256", HERE)
check(len(manifest) == 6, "audit manifest has six rows")
check({Path(row).name for row in manifest} == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "audit manifest file set exact")
seal = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
        if line.strip()]
check(len(seal) == 1, "one audit seal row")
expected, relative = seal[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "audit seal targets manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest() == expected,
      "audit seal resolves")

audit = " ".join((HERE / "AUDIT.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
for token in (
    "[H_hex,N_a]=0",
    "im C(z)=ker B(z)",
    "C_1 C_1^T = 4",
    "structure factor is indispensable",
    "do not prove a gapless mode or pole",
    "conditional harmonic result",
    "not call it physical momentum",
    "exact locked mean-square operator overlap",
    "two-mode form factor",
    "K_(E<-loop)(omega,1)=0",
    "Sym^2_0(T2)=E+T2",
    "finite-group composite channel",
    "not a selected locked ground-state witness",
    "Hostile verdict: PASS",
):
    check(token.lower() in audit.lower(), f"audit scope token: {token}")

for token in (
    "PASS__INDEPENDENT_GL6AS_HOSTILE_REPLAY__11059/11059",
    "PASS__GL6AS_NATIVE_HEXAGON_COLLECTIVE__1704/1704",
    "PASS__GL6AS_PACKET__125/125",
    "PASS__GL6AO_COMPLETE_SIXTH_ORDER__120304/120304",
    "PASS__GL6AO_HOSTILE_AUDIT_PACKET__156/156",
    "GL6AP exact verification: PASS (434/434)",
    "PASS__GL6AP_HOSTILE_AUDIT_PACKET__208/208",
    "PASS__GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION__1010/1010",
    "PASS__GL6AQ_AUDIT_PACKET__209/209",
    "PASS__GL6AS_HOSTILE_AUDIT_PACKET__",
):
    check(token in verification, f"verification token: {token}")

for forbidden in (
    "t2 is a photon", "t2 is a graviton", "the character is physical momentum",
    "the composite is gravity", "derives a physical cone", "derives stress",
    "derives ricci", "derives newton's constant", "therefore g=j",
):
    check(forbidden not in audit.lower(), f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6AS_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
