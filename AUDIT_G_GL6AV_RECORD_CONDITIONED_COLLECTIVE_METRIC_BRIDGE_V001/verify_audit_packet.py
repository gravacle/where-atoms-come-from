#!/usr/bin/env python3
"""Fail-closed custody verifier for the distinct GL6AV hostile audit."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR_NAME = "LANE_CROSS_RFT_GRA_GL6AV_RECORD_CONDITIONED_COLLECTIVE_METRIC_BRIDGE_V001"
AUTHOR = ROOT / AUTHOR_NAME
AUTHOR_FILES = {
    "README.md", "RESULT.md", "THEOREM.md", "DEPENDENCIES.md",
    "SELF_AUDIT.md", "VERIFICATION.txt",
    "verify_record_conditioned_metric_bridge.py", "verify_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}
AUDIT_FILES = {
    "README.md", "AUDIT.md", "AUDITED_TARGETS.sha256",
    "independent_gl6av_replay.py", "VERIFICATION.txt", "verify_audit_packet.py",
}
DEPENDENCIES = {
    "LANE_CROSS_RFT_GRA_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/THEOREM.md":
        "f75edcb115c3f7c86c6598f4597366b36e363df2d03ad919cc607b57dfb6b20c",
    "LANE_CROSS_RFT_GRA_GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION_V001/THEOREM.md":
        "1d1b01380ec8fd7ce83c69d45b68d9bde36bbe1dacdd32e3a5909ee6723a5ace",
    "LANE_CROSS_RFT_GRA_GL6AS_NATIVE_HEXAGON_COLLECTIVE_RESPONSE_V001/THEOREM.md":
        "bfe36071a24ccc7d6d7a16afeeea1b5554a95562ae91ac59c709db478000db9f",
    "LANE_CROSS_RFT_GRA_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001/THEOREM.md":
        "8407cee5196bfa4240f02159a5f59f941903dcf7a10e2baa18cf52a01ac8f743",
    "LANE_CROSS_RFT_GRA_GL6AT_PRIMARY_QUANTUM_ICE_CROSSWALK_V001/RESULT.md":
        "1c2b48c40a88000a88f8446fa2aea5116acb0896d94848274c87a22003bbcad9",
}
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def rows(path: Path):
    parsed = []
    seen = set()
    for raw in path.read_text().splitlines():
        if not raw.strip():
            continue
        expected, relative = raw.split(maxsplit=1)
        relative_path = Path(relative)
        check(len(expected) == 64 and all(c in "0123456789abcdef" for c in expected),
              f"valid digest: {relative}")
        check(not relative_path.is_absolute() and ".." not in relative_path.parts,
              f"safe path: {relative}")
        check(relative not in seen, f"unique row: {relative}")
        target = ROOT / relative_path
        check(target.is_file(), f"target exists: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        seen.add(relative)
        parsed.append((expected, relative))
    return parsed


targets = rows(HERE / "AUDITED_TARGETS.sha256")
check(len(targets) == 10, "exact audited target count")
check({Path(relative).name for _, relative in targets} == AUTHOR_FILES,
      "all and only author files pinned")
check(all(Path(relative).parent == Path(AUTHOR_NAME) for _, relative in targets),
      "audited targets confined to author packet")
check({path.name for path in AUTHOR.iterdir() if path.is_file()} == AUTHOR_FILES,
      "author packet exact inventory")
check(not any(path.is_dir() for path in AUTHOR.iterdir()), "author packet has no subdirectories")

author_manifest = rows(AUTHOR / "MANIFEST.sha256")
check(len(author_manifest) == 8, "author manifest row count")
check({Path(relative).name for _, relative in author_manifest}
      == AUTHOR_FILES - {"MANIFEST.sha256", "SEAL.sha256"},
      "author manifest exact coverage")
check(hashlib.sha256((AUTHOR / "MANIFEST.sha256").read_bytes()).hexdigest()
      == "283e1010be399b5e31cf93da34dee4c075dc29f9a98d9b4804ba6ea2e411073b",
      "frozen author manifest identity")
author_seal = rows(AUTHOR / "SEAL.sha256")
check(len(author_seal) == 1 and
      author_seal[0][1] == f"{AUTHOR_NAME}/MANIFEST.sha256",
      "author seal targets author manifest")
check(hashlib.sha256((AUTHOR / "SEAL.sha256").read_bytes()).hexdigest()
      == "c2ea7f69d86542583f58dffa5b8f6207c1b28b77c941ee016e1ac5f69c76b191",
      "frozen author seal-file identity")

for relative, expected in DEPENDENCIES.items():
    target = ROOT / relative
    check(target.is_file(), f"dependency exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"dependency hash: {relative}")

audit_manifest = rows(HERE / "MANIFEST.sha256")
check(len(audit_manifest) == 6, "audit manifest row count")
check({Path(relative).name for _, relative in audit_manifest} == AUDIT_FILES,
      "audit manifest exact coverage")
check(all(Path(relative).parent == Path(HERE.name) for _, relative in audit_manifest),
      "audit manifest path confinement")
check({path.name for path in HERE.iterdir() if path.is_file()}
      == AUDIT_FILES | {"MANIFEST.sha256", "SEAL.sha256"},
      "audit exact inventory")
check(not any(path.is_dir() for path in HERE.iterdir()), "audit has no subdirectories")
audit_seal = rows(HERE / "SEAL.sha256")
check(len(audit_seal) == 1 and
      audit_seal[0][1] == f"{HERE.name}/MANIFEST.sha256",
      "audit seal targets audit manifest")

replay = (HERE / "independent_gl6av_replay.py").read_text()
for forbidden in ("import verify_packet", "import verify_record_conditioned", "runpy"):
    check(forbidden not in replay, f"independent replay excludes {forbidden}")

env = os.environ.copy()
env["PYTHONPYCACHEPREFIX"] = "/private/tmp/gl6av_distinct_audit_pycache"
commands = (
    ("independent", HERE / "independent_gl6av_replay.py",
     "PASS__GL6AV_INDEPENDENT_REPLAY__"),
    ("author_exact", AUTHOR / "verify_record_conditioned_metric_bridge.py",
     "PASS__GL6AV_EXACT_REPLAY_COMPLETE__67/67"),
    ("author_packet", AUTHOR / "verify_packet.py",
     "PASS__GL6AV_AUTHOR_PACKET__62/62"),
)
for name, script, marker in commands:
    for optimized in (False, True):
        command = [sys.executable]
        if optimized:
            command.append("-O")
        command.append(str(script))
        completed = subprocess.run(command, cwd=ROOT, env=env,
                                   capture_output=True, text=True, check=False)
        check(completed.returncode == 0,
              f"{name} exit {'optimized' if optimized else 'normal'}")
        check(marker in completed.stdout,
              f"{name} marker {'optimized' if optimized else 'normal'}")

audit = " ".join((HERE / "AUDIT.md").read_text().split())
for marker in (
    "not the complete nonuniform effective Hamiltonian",
    "time-dependent unitary to finite volume",
    "unique multilinear extension",
    "does not authenticate one infinite homogeneous query",
    "lower-order diagonal shifts can be configuration dependent",
    "L_q=q^6L_1",
    "Independently selected `q`-dependent stationary states",
    "not for a purported global Hamiltonian",
    "determinant `-48`",
    "one half of its stated signed four-term combination",
    "not four authenticated retained-source directions",
    "pulse/read coordinate",
    "two modes with distinct speeds",
    "ds_q^2=-q^{12}v_1^2dt^2+h_{ij}dx^idx^j",
    "not physical two-way back-reaction",
    "No authenticated global homogeneous `q`",
    "Newton's constant `G`",
    "**Hostile verdict: PASS.**",
):
    check(marker in audit, f"audit scope marker: {marker}")

verification = (HERE / "VERIFICATION.txt").read_text()
check("PASS__GL6AV_INDEPENDENT_REPLAY__" in verification,
      "independent replay recorded")
check("PASS__GL6AV_AUTHOR_PACKET__62/62" in verification,
      "author packet replay recorded")
check("PASS__GL6AV_HOSTILE_AUDIT_PACKET__" in verification,
      "audit packet replay recorded")
check("NORMAL_AND_OPTIMIZED" in verification, "optimized modes recorded")

print(f"PASS__GL6AV_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
