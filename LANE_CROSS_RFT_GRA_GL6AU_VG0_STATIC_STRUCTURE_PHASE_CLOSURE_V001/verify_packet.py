#!/usr/bin/env python3
"""Fail-closed custody, inventory, and claim-scope verifier for GL6AU."""

from __future__ import annotations

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
    "THEOREM.md",
    "RESULT.md",
    "EVIDENCE_AND_HYPOTHESIS.md",
    "SELF_AUDIT.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "VERIFICATION.txt",
    "verify_vg0_static_structure_closure.py",
    "verify_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
)
for name in required:
    check((HERE / name).is_file(), f"required file: {name}")


# Exact upstream custody.
dependency_paths = set()
for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"dependency exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"dependency hash: {relative}")
    dependency_paths.add(relative)
check(len(dependency_paths) == 32, "exact dependency count")
for marker in ("GL6AO", "GL6AR", "GL6AS"):
    check(sum(marker in path for path in dependency_paths) == 6,
          f"six {marker} author/audit custody objects")
check(sum("GL6AT" in path for path in dependency_paths) == 14,
      "eleven GL6AT author objects plus three audit custody objects")
check(sum("AUDIT_G_GL6AT" in path for path in dependency_paths) == 3,
      "three GL6AT hostile-audit custody objects")
check(all(marker not in path for path in dependency_paths
          for marker in ("GL6AL", "GL6AV", "GL6AW")),
      "mutable and later lanes excluded")

at_expected = {
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "EVIDENCE_LADDER.md",
    "MANIFEST.sha256", "PRIMARY_SOURCES.md", "README.md", "RESULT.md",
    "SEAL.sha256", "SELF_AUDIT.md", "VERIFICATION.txt", "verify_packet.py",
}
at_seen = {
    Path(path).name for path in dependency_paths
    if path.startswith("LANE_CROSS_RFT_GRA_GL6AT_PRIMARY_QUANTUM_ICE_CROSSWALK")
}
check(at_seen == at_expected, "complete exact GL6AT author snapshot")


# Exact local inventory and seal.
manifest_paths = set()
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    relative_path = Path(relative)
    check(relative_path.parent == Path(HERE.name),
          f"manifest target stays in packet: {relative}")
    target = ROOT / relative_path
    check(target.is_file(), f"manifest target exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"manifest hash: {relative}")
    manifest_paths.add(relative)
for name in required:
    if name not in ("MANIFEST.sha256", "SEAL.sha256"):
        check(f"{HERE.name}/{name}" in manifest_paths,
              f"manifest coverage: {name}")
check(len(manifest_paths) == 10, "exact author manifest count")

seal_lines = [
    line for line in (HERE / "SEAL.sha256").read_text().splitlines()
    if line.strip()
]
check(len(seal_lines) == 1, "one seal row")
expected, relative = seal_lines[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest()
      == expected, "seal hash")


# Claim-class and theorem-scope replay.
theorem = " ".join((HERE / "THEOREM.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
ledger = " ".join((HERE / "EVIDENCE_AND_HYPOTHESIS.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
dependencies = " ".join((HERE / "DEPENDENCIES.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())

for document in (theorem, readme):
    check("author frozen and sealed" in document, "frozen author status")
check("independent hostile audit" in theorem.lower(),
      "theorem requires distinct hostile audit")
check("PASS__GL6AU_VG0_STATIC_STRUCTURE__8111/8111" in verification,
      "physics replay count recorded")
check("both author-frozen and independently audited" in dependencies,
      "AT author/audit custody typed explicitly")

required_theorem_tokens = (
    "preserved by all three unit cell translations",
    "||C(z)^*u||^2=3|1-exp(iq_L)|^2",
    "0<=t_d<=1",
    "Delta_C(L) <= 6J sin^2(pi/L)/S_(u,L)(q_L)",
    "Static-exponent closure theorem.",
    "any exponent below two suffices",
    "Var(F_c)+Var(F_s)=N S_(u,L)(q_L)",
    "Extensive variance was sufficient, not necessary",
    "f(q)=O(q^2) does not imply Delta(q)->0",
    "H_0=-JA=H_RK-JD",
    "Flux-sector towers are not a GNS proof",
    "No exact transfer-matrix, measure, or coupling map",
    "GL6AT evidence classification has passed its distinct hostile audit",
    "ICE0-STATIC",
    "ICE0-GNS-BRIDGE",
    "Neither hypothesis is proved here",
)
for token in required_theorem_tokens:
    check(token in theorem, f"theorem scope token: {token}")

check("The phase problem has been reduced to one measured static exponent"
      in result, "result gives exact residual gate")
check("strong numerical/effective support" in ledger,
      "evidence strength is typed")
check("finite numerical evidence" in self_audit,
      "self-audit preserves numerical ceiling")

aggregate = " ".join((theorem, readme, result, ledger, self_audit)).lower()
for forbidden in (
    "ice0-static is proved",
    "ice0-gns-bridge is proved",
    "the v/g=0 phase is rigorously proved",
    "the pole is proved",
    "the character is physical momentum",
    "is a physical photon",
    "is a graviton",
    "is gravity",
    "therefore derives a physical cone",
    "we derive the einstein equation",
    "we derive newton's constant",
):
    check(forbidden not in aggregate, f"forbidden promotion absent: {forbidden}")


print(f"PASS__GL6AU_PACKET__{checks}/{checks}")
