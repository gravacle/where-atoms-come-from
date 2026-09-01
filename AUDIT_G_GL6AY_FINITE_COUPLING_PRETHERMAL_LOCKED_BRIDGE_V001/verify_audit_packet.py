#!/usr/bin/env python3
"""Fail-closed verifier for the independent GL6AY hostile audit packet."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(label)
    checks += 1


required = {
    "README.md", "AUDIT.md", "PRIMARY_SOURCE_REPLAY.md",
    "AUDITED_TARGETS.sha256", "VERIFICATION.txt",
    "independent_gl6ay_replay.py", "verify_audit_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}
check({path.name for path in HERE.iterdir() if path.is_file()} == required,
      "audit inventory exact")


def rows(path: Path, parent: Path | None = None):
    answer = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        check(len(expected) == 64 and
              all(char in "0123456789abcdef" for char in expected),
              f"hash syntax {relative}")
        check(relative not in answer, f"unique row {relative}")
        target = ROOT / relative
        check(target.is_file() and not target.is_symlink(),
              f"regular target {relative}")
        if parent is not None:
            check(target.parent == parent, f"target parent {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash resolves {relative}")
        answer.append(relative)
    return answer


targets = rows(HERE / "AUDITED_TARGETS.sha256", AUTHOR)
check(len(targets) == 12, "all twelve frozen author files pinned")
check({Path(item).name for item in targets} ==
      {path.name for path in AUTHOR.iterdir() if path.is_file()},
      "target set equals author set")

author_manifest_hash = hashlib.sha256((AUTHOR / "MANIFEST.sha256").read_bytes()).hexdigest()
check(author_manifest_hash ==
      "d51c6aea006c1b5cdc7a75023dfd59e0cdb363549565c92aab0a5e1ee5083710",
      "frozen author manifest")
author_seal = (AUTHOR / "SEAL.sha256").read_text().strip().split(maxsplit=1)
check(len(author_seal) == 2 and author_seal[0] == author_manifest_hash,
      "author seal resolves")

manifest = rows(HERE / "MANIFEST.sha256", HERE)
check(len(manifest) == 7, "seven audit content rows")
check({Path(item).name for item in manifest} == required - {
    "MANIFEST.sha256", "SEAL.sha256"
}, "audit manifest inventory exact")
seal = (HERE / "SEAL.sha256").read_text().strip().split(maxsplit=1)
check(len(seal) == 2, "one audit seal row")
check(seal[1] == f"{HERE.name}/MANIFEST.sha256", "audit seal target")
check(seal[0] == hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest(),
      "audit seal resolves")

audit = " ".join((HERE / "AUDIT.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
source = " ".join((HERE / "PRIMARY_SOURCE_REPLAY.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
author_theorem = " ".join((AUTHOR / "THEOREM.md").read_text().split())

disposition = "FAIL__REPAIR_REQUIRED__GLOBAL_LOCK_PROJECTOR_AND_DRESSED_SUBSPACE_SCOPE"
check(disposition in audit and disposition in readme, "failure disposition explicit")
for token in (
    "Material defect one",
    "no global spectral projector `P`",
    "P_S^0=chi(N_S=0)",
    "[D_hat(S),N_S]=0",
    "preserves the twist estimate",
    "Material defect two",
    "supplies no volume-uniform estimate",
    "small local rotations",
    "First winding coefficient — pass",
    "Whole-band boundary — pass",
    "Required author repair",
):
    check(token in audit, f"audit finding token {token}")
for token in (
    "1509.05386v3", "Theorem 3.1", "Theorem 3.3",
    "1704.08703v2", "Appendix A",
    "1105.0675v1", "Section 4",
):
    check(token in source, f"source replay token {token}")
check("P D_hat P has exact termwise port U(1)^4" in author_theorem,
      "frozen target contains first disputed promotion")
check("physical dressed space `Y^*P`" in author_theorem,
      "frozen target contains second disputed promotion")
check("PASS__INDEPENDENT_GL6AY_HOSTILE_REPLAY__" in verification,
      "independent replay recorded")
check("PASS__GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE__" in verification,
      "author replay recorded")
check("PASS__GL6AY_PACKET__" in verification,
      "author packet replay recorded")

print(f"PASS__GL6AY_HOSTILE_AUDIT_PACKET__{checks}/{checks}")
