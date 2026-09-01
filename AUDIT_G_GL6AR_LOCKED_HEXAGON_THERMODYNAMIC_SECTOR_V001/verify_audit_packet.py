#!/usr/bin/env python3
"""Fail-closed custody verifier for the distinct GL6AR hostile audit."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AR_LOCKED_HEXAGON_THERMODYNAMIC_SECTOR_V001"
AUTHOR_FILES = {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "VERIFICATION.txt",
    "verify_locked_hexagon_thermodynamics.py", "verify_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}
AUDIT_CORE = {
    "README.md", "AUDIT.md", "AUDITED_TARGETS.sha256",
    "independent_gl6ar_replay.py", "VERIFICATION.txt", "verify_audit_packet.py",
}
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise RuntimeError(label)
    checks += 1


def rows(path):
    parsed = []
    seen = set()
    for raw in path.read_text().splitlines():
        if not raw.strip():
            continue
        expected, relative = raw.split(maxsplit=1)
        check(len(expected) == 64 and all(c in "0123456789abcdef" for c in expected),
              f"valid digest syntax: {relative}")
        check(relative not in seen, f"unique row: {relative}")
        seen.add(relative)
        target = ROOT / relative
        check(target.is_file(), f"target exists: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        parsed.append((expected, relative))
    return parsed


targets = rows(HERE / "AUDITED_TARGETS.sha256")
check(len(targets) == 11, "exact audited target count")
check({Path(relative).name for _, relative in targets} == AUTHOR_FILES,
      "all and only author files pinned")
check({path.name for path in AUTHOR.iterdir() if path.is_file()} == AUTHOR_FILES,
      "author directory inventory exact")

author_manifest = rows(AUTHOR / "MANIFEST.sha256")
check(len(author_manifest) == 9, "author manifest count")
check({Path(relative).name for _, relative in author_manifest}
      == AUTHOR_FILES - {"MANIFEST.sha256", "SEAL.sha256"},
      "author manifest coverage exact")
author_seal = rows(AUTHOR / "SEAL.sha256")
check(len(author_seal) == 1 and Path(author_seal[0][1]).name == "MANIFEST.sha256",
      "author seal targets manifest")

dependencies = rows(AUTHOR / "DEPENDENCIES.sha256")
check(len(dependencies) == 12, "exact dependency count")
check(all("GL6AO" in relative or "GL6AN" in relative
          for _, relative in dependencies), "dependencies confined to AO/AN")
check(sum("GL6AO" in relative for _, relative in dependencies) == 6,
      "six AO custody rows")
check(sum("GL6AN" in relative for _, relative in dependencies) == 6,
      "six AN custody rows")
check(all(not any(later in relative for later in ("GL6AP", "GL6AQ", "GL6AR", "GL6AS", "GL6AT"))
          for _, relative in dependencies), "no later mutable lane in custody")

audit_manifest = rows(HERE / "MANIFEST.sha256")
check(len(audit_manifest) == 6, "audit manifest count")
check({Path(relative).name for _, relative in audit_manifest} == AUDIT_CORE,
      "audit manifest coverage exact")
check({path.name for path in HERE.iterdir() if path.is_file()}
      == AUDIT_CORE | {"MANIFEST.sha256", "SEAL.sha256"},
      "audit directory inventory exact")
audit_seal = rows(HERE / "SEAL.sha256")
check(len(audit_seal) == 1 and Path(audit_seal[0][1]).name == "MANIFEST.sha256",
      "audit seal targets audit manifest")

replay = (HERE / "independent_gl6ar_replay.py").read_text()
check("import verify_locked_hexagon_thermodynamics" not in replay,
      "independent replay imports no author verifier")
audit = (HERE / "AUDIT.md").read_text()
for marker in (
    "an 18-link support for every `Q_4` term",
    "three plane-independent integer coordinate cut fluxes",
    "Delta_C=t(rho_C-lambda_2(C))",
    "= 18t||w||_infinity^2 L/V_L",
    "finite periodic subsequence bound",
    "P_0=1_{\\{0\\}}(H_omega)",
    "**Hostile verdict: PASS.**",
):
    check(marker in audit, f"audit theorem marker: {marker}")

theorem = " ".join((AUTHOR / "THEOREM.md").read_text().split())
result = " ".join((AUTHOR / "RESULT.md").read_text().split())
for marker in (
    "does not identify that closure with infinite-volume GNS spectral gaplessness",
    "does not prove retention of a nonlocal asymptotic sector",
    "If `V_L>0`",
    "projection onto the full zero-energy subspace",
    "convergence of the energy density or existence of a thermodynamic energy-density limit is not proved",
):
    check(marker in theorem, f"repaired theorem ceiling: {marker}")
check("not a fourth derived geometric cut" in result,
      "result distinguishes dependent port count from geometric cuts")

combined = (audit + " " + theorem + " " + result).lower()
for forbidden in (
    "proves gns gaplessness",
    "the mode is a photon",
    "the mode is a graviton",
    "derives physical momentum",
    "derives gravity",
    "derives newton's constant",
):
    check(forbidden not in combined, f"forbidden promotion absent: {forbidden}")

verification = (HERE / "VERIFICATION.txt").read_text()
for marker in (
    "PASS__GL6AR_INDEPENDENT_REPLAY__",
    "PASS__GL6AR_LOCKED_HEXAGON_THERMODYNAMICS__72161/72161",
    "PASS__GL6AR_PACKET__145/145",
    "PASS__GL6AR_HOSTILE_AUDIT_PACKET__",
):
    check(marker in verification, f"verification marker: {marker}")

print(f"PASS__GL6AR_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
