#!/usr/bin/env python3
"""Fail-closed custody, replay, and scope verifier for GL6CR."""

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
    "VERIFICATION.txt", "derive_full_pair_rotational_completion.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required packet file {name}")

dependency_lines = [
    line for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
    if line.strip()
]
check(len(dependency_lines) == 21, "twenty-one exact dependency bytes pinned")
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
check(ledger["metric_solder"]["rank"] == 6, "solder rank six")
check(ledger["quadratic_symbols"]["S4_dimension_reynolds"] == 9,
      "cubic quadratic dimension nine")
check(ledger["quadratic_symbols"]["SO3_self_adjoint_dimension"] == 4,
      "rotational quadratic dimension four")
check(ledger["quadratic_symbols"]["matching_codimension"] == 5,
      "rotational matching codimension five")
direct = ledger["direct_cubic_gauge_shortcut"]
check(direct["constraint_rank_on_S4_nine_space"] == 8,
      "direct Ward constraint rank eight")
check(direct["null_dimension"] == 1 and direct["equals_Einstein_ray"] is True,
      "unique direct Ward null is Einstein ray")
check(direct["unique_null_ray_in_invariant_basis"] ==
      ["0", "2", "-2", "0", "-2", "-2", "4", "-1", "1"],
      "frozen direct Ward null coordinates")
check("not yet derived" in direct["guard"], "physical Ward derivation ceiling")

replay = subprocess.run(
    [sys.executable, "-B", str(HERE / "derive_full_pair_rotational_completion.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(replay.returncode == 0, "exact replay exits zero")
for token in (
    "PASS__GL6CR_FULL_PAIR_ROTATIONAL_COMPLETION__599/599",
    "QUADRATIC_DIMENSIONS=S4_9;SO3_4;FIVE_MATCHING_CONDITIONS",
    "T2_RESTRICTION_DIMENSION=2;FULL_E2_T2_COMPLETION_MANDATORY",
    "DIRECT_CUBIC_GAUGE_NULL=RANK8;UNIQUE_EINSTEIN_RAY;SO3_FOLLOWS",
):
    check(token in replay.stdout, f"exact replay token {token}")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
for token in (
    "\\operatorname{rank}A_{\\rm Ward}=8", "\\dim\\ker A_{\\rm Ward}=1",
    "rotational completion follows", "from the parent relational redundancy",
    "Imposing (CR22) by hand", "gravity, or `G`",
):
    check(token in theorem, f"theorem claim/scope token {token}")
check("rank eight" in result and "must be derived" in result,
      "result preserves direct-Ward ceiling")
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

print(f"PASS__GL6CR_PACKET__{checks}/{checks}")
