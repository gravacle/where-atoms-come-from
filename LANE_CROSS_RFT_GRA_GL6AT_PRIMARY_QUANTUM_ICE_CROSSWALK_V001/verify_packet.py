#!/usr/bin/env python3
"""Fail-closed local custody and claim-structure verifier for GL6AT.

This verifier checks frozen internal bytes, packet inventory, exact elementary
crosswalk arithmetic, primary-source metadata, and advertised claim ceilings.
It deliberately does not claim to reproduce or authenticate external papers.
"""

from __future__ import annotations

from fractions import Fraction
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
    "RESULT.md",
    "PRIMARY_SOURCES.md",
    "EVIDENCE_LADDER.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "SELF_AUDIT.md",
    "VERIFICATION.txt",
    "verify_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
)

for name in required:
    check((HERE / name).is_file(), f"required file: {name}")
check(not any(path.is_dir() for path in HERE.iterdir()), "no packet subdirectories")
check(sorted(path.name for path in HERE.iterdir() if path.is_file())
      == sorted(required), "exact packet inventory")


def verify_hash_rows(filename: str) -> set[str]:
    paths: set[str] = set()
    for line in (HERE / filename).read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        check(len(expected) == 64 and all(ch in "0123456789abcdef" for ch in expected),
              f"well-formed hash: {relative}")
        check(relative not in paths, f"unique hash target: {relative}")
        target = ROOT / relative
        check(target.is_file(), f"hash target exists: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"hash matches: {relative}")
        paths.add(relative)
    return paths


dependency_paths = verify_hash_rows("DEPENDENCIES.sha256")
check(len(dependency_paths) == 12, "exact dependency count")
check(all(("GL6AO" in path or "GL6AP" in path) for path in dependency_paths),
      "only sealed GL6AO/GL6AP targets imported")
check(sum("GL6AO" in path for path in dependency_paths) == 6,
      "six GL6AO author/audit targets")
check(sum("GL6AP" in path for path in dependency_paths) == 6,
      "six GL6AP author/audit targets")
for token in (
    "LANE_CROSS_RFT_GRA_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/THEOREM.md",
    "AUDIT_G_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/AUDIT.md",
    "LANE_CROSS_RFT_GRA_GL6AP_LOCKED_IR_CONDITIONAL_RESPONSE_V001/THEOREM.md",
    "AUDIT_G_GL6AP_LOCKED_IR_CONDITIONAL_RESPONSE_V001/AUDIT.md",
):
    check(token in dependency_paths, f"pinned dependency: {token}")

manifest_paths = verify_hash_rows("MANIFEST.sha256")
expected_manifest = {
    f"{HERE.name}/{name}" for name in required
    if name not in ("MANIFEST.sha256", "SEAL.sha256")
}
check(manifest_paths == expected_manifest, "exact manifest coverage")
check(len(manifest_paths) == 9, "exact manifest count")

seal_lines = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
              if line.strip()]
check(len(seal_lines) == 1, "one seal row")
seal_hash, seal_target = seal_lines[0].split(maxsplit=1)
check(seal_target == f"{HERE.name}/MANIFEST.sha256", "seal targets manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest()
      == seal_hash, "seal hash matches")

# Exact elementary crosswalk arithmetic, independent of prose.
check(Fraction(63, 8) > 0, "positive ring coefficient")
tetra = []
for a in range(4):
    tetra.append(tuple(Fraction(1 if a == b else 0) - Fraction(1, 4)
                       for b in range(4)))
for a in range(4):
    check(sum(tetra[a]) == 0, f"tetra vector in A3 plane: {a}")
    for b in range(4):
        dot = sum(tetra[a][j] * tetra[b][j] for j in range(4))
        check(dot == (Fraction(3, 4) if a == b else Fraction(-1, 4)),
              f"tetra Gram entry: {a},{b}")

readme = " ".join((HERE / "README.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
sources = " ".join((HERE / "PRIMARY_SOURCES.md").read_text().split())
ladder = " ".join((HERE / "EVIDENCE_LADDER.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
all_claims = " ".join((readme, result, sources, ladder, self_audit))

for token in (
    "g={63\\over8}{h^6\\over U_d^5}>0",
    "v/g=0",
    "standard diamond net",
    "two dimers touching every",
    "fully-packed loop",
    "not the standard degree-one diamond quantum dimer model",
    "The exactly soluble RK point is `v/g=1`",
    "uncomputed `O(h^8/U_d^7)` completion",
    "fixed `Q_4` theorem by itself has no literal infrared limit",
    "No primary result located in this screen proves",
):
    check(token in readme + " " + result, f"exact scope token: {token}")

for doi in (
    "10.1103/PhysRevB.69.064404",
    "10.1103/PhysRevLett.108.067204",
    "10.1103/PhysRevB.86.075154",
    "10.1103/PhysRevB.96.035136",
):
    check(doi in sources, f"primary DOI recorded: {doi}")
for version in (
    "cond-mat/0305401v3",
    "1105.4196v3",
    "1204.1325v2",
    "1703.03836v1",
):
    check(version in sources, f"versioned arXiv source: {version}")

for token in (
    "finite numerical evidence",
    "not a rigorous thermodynamic theorem",
    "effective-theory outputs",
    "approximate composite-operator mechanism",
    "not the sealed two-dimensional local pair `E` channel",
    "absence is a screen result, not a proof",
):
    check(token in all_claims, f"evidence ceiling token: {token}")

for token in (
    "`mu=v/g`",
    "R_{01}-R_{23}",
    "M_{01}+M_{23}",
    "Z_e=2n_e-1=2S_e^z",
    "three-dimensional `T2g`",
    "two-dimensional centered complementary-pair sums",
    "Eq. (57c)",
    "Eq. (55c)",
    "Omega^4",
    "weight tending to zero",
):
    check(token in result + " " + sources + " " + ladder,
          f"operator-overlap token: {token}")

for forbidden in (
    "the `v/g=0` phase is rigorously proved",
    "GL6AO proves a microscopic pole",
    "GL6AO derives a physical cone",
    "GL6AO derives gravity",
    "GL6AO derives Newton's constant",
    "local pair `E` and `T2g` are identical",
    "A3 character is calibrated physical momentum",
):
    check(forbidden not in all_claims, f"forbidden promotion absent: {forbidden}")

check("PASS__GL6AT_PACKET__" in verification,
      "verification records packet pass")
check("external primary-paper content not locally rederived" in verification,
      "verification states honest external-source limit")

print(f"PASS__GL6AT_PACKET__{checks}/{checks}")
