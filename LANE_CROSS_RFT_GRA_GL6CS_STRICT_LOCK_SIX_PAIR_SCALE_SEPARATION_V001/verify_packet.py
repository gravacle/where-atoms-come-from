#!/usr/bin/env python3
"""Fail-closed custody, replay, and scope verifier for GL6CS."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def unique_object(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


required = {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "EXACT_LEDGER.json",
    "VERIFICATION.txt", "derive_strict_lock_scale_separation.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required packet file {name}")

dependency_lines = [
    line for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
    if line.strip()
]
check(len(dependency_lines) == 36, "thirty-six exact dependency bytes pinned")
dependency_names = set()
for line in dependency_lines:
    expected, relative = line.split("  ", 1)
    check(relative not in dependency_names, f"unique dependency {relative}")
    dependency_names.add(relative)
    candidate = ROOT / relative
    check(candidate.is_file(), f"dependency exists {relative}")
    check(digest(candidate) == expected, f"dependency hash {relative}")

ledger = json.loads((HERE / "EXACT_LEDGER.json").read_text(),
                    object_pairs_hook=unique_object)
ring = ledger["local_ring_flip_census"]
check(ring["count"] == 24 and len(ring["rows"]) == 24,
      "all twenty-four locked ring incidences")
check(all(sum(int(value) * int(value) for value in row["difference"]) == 16
          for row in ring["rows"]), "every ring difference norm squared sixteen")
scales = ledger["scales"]
check(scales["coupling_exponents_at_fixed_component"] == {
          "EE": -6, "ET": 0, "TT_contact": 2, "TT_writer": 6,
      }, "four exact fixed-component powers")
check(scales["no_orientation_mixing_required_collective_enhancement"] == {
          "contact_relative_to_E": "O(r^-8)",
          "writer_susceptibility_relative_to_E": "O(r^-12)",
      }, "required collective enhancement powers")
check("fixed-frame repetition" in ledger["accumulation_boundary"],
      "fixed-frame repetition boundary")
check("no thermodynamic phase" in ledger["scope"] and "gravity" in ledger["scope"],
      "promotion ceiling")

replay = subprocess.run(
    [sys.executable, "-B", str(HERE / "derive_strict_lock_scale_separation.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(replay.returncode == 0, "exact replay exits zero")
for token in (
    "PASS__GL6CS_STRICT_LOCK_SCALE_SEPARATION__110/110",
    "LOCKED_RING_CHANGE=PURE_E2_NONZERO;24/24",
    "FIXED_COMPONENT_EXPONENTS=EE_-6;ET_0;TT_CONTACT_2;TT_WRITER_6",
    "FIXED_FRAME_STRICT_LOCK_ROTATIONAL_EQUALITY=IMPOSSIBLE_FOR_BOUNDED_COEFFICIENTS",
):
    check(token in replay.stdout, f"exact replay token {token}")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
for token in (
    "K_{EE}=O(r^{-6}/U_d)", "K_{ET}=O(1/U_d)",
    "K_{TT}^{\\rm ct}=O(r^2/U_d)", "K_{TT}^{\\rm wr}=O(r^6/U_d)",
    "`O(r^-8)`", "`O(r^-12)`",
    "does not rule out", "does not select any of those possibilities",
):
    check(token in theorem, f"theorem claim/scope token {token}")
check("cannot satisfy" in result and "does not choose" in result,
      "result preserves accumulation ceiling")
check("does not replace an independent hostile audit" in self_audit,
      "self-audit independence ceiling")

manifest_lines = [
    line for line in (HERE / "MANIFEST.sha256").read_text().splitlines()
    if line.strip()
]
manifest_names = set()
for line in manifest_lines:
    expected, name = line.split("  ", 1)
    check(Path(name).parent == Path("."), f"manifest local path {name}")
    check(name not in manifest_names, f"manifest unique path {name}")
    manifest_names.add(name)
    check(digest(HERE / name) == expected, f"manifest hash {name}")
check(manifest_names == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "manifest covers every non-custody packet byte")
seal = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal == [digest(HERE / "MANIFEST.sha256"), "MANIFEST.sha256"],
      "seal names and hashes packet manifest")

print(f"PASS__GL6CS_PACKET__{checks}/{checks}")
