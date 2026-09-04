#!/usr/bin/env python3
"""Fail-closed verifier for the independent GL6CR hostile audit."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CR_FULL_PAIR_ROTATIONAL_COMPLETION_V001"
TARGET_MANIFEST_DIGEST = "825cccc8232bbc292b55f9bab72410a60e209243ec3822c116658c2a2a55ab4c"
TARGET_SEAL_FILE_DIGEST = "ed06c8c30d04d04e067de34dadff2bdcc9d93542faef630ac37bfd7ea9d127a5"
checks = 0


def check(condition, label):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict_json(path):
    def no_duplicates(pairs):
        out = {}
        for key, value in pairs:
            if key in out:
                raise ValueError(f"duplicate JSON key: {key}")
            out[key] = value
        return out

    return json.loads(path.read_text(), object_pairs_hook=no_duplicates)


required = {
    "README.md",
    "AUDIT_REPORT.md",
    "TARGET.sha256",
    "INDEPENDENT_RESULT.json",
    "VERIFICATION.txt",
    "independent_gl6cr_audit.py",
    "verify_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required audit file: {name}")


# Pin every byte of the frozen author packet.
target_paths = set()
for line in (HERE / "TARGET.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    candidate = ROOT / relative
    check(candidate.is_file(), f"target exists: {relative}")
    check(digest(candidate) == expected, f"target hash: {relative}")
    check(relative not in target_paths, f"target row unique: {relative}")
    target_paths.add(relative)
check(len(target_paths) == 12, "exact target byte count")
check(digest(TARGET / "MANIFEST.sha256") == TARGET_MANIFEST_DIGEST,
      "frozen target manifest digest")
check(digest(TARGET / "SEAL.sha256") == TARGET_SEAL_FILE_DIGEST,
      "frozen target seal-file digest")

target_manifest_names = set()
for line in (TARGET / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, name = line.split(maxsplit=1)
    candidate = TARGET / name
    check(candidate.is_file(), f"target manifest file exists: {name}")
    check(digest(candidate) == expected, f"target manifest hash: {name}")
    check(name not in target_manifest_names, f"target manifest row unique: {name}")
    check(Path(name).parent == Path("."), f"target manifest local path: {name}")
    target_manifest_names.add(name)
check(len(target_manifest_names) == 10, "exact target manifest count")
target_seal = [line for line in (TARGET / "SEAL.sha256").read_text().splitlines()
               if line.strip()]
check(len(target_seal) == 1, "one target seal row")
expected, name = target_seal[0].split(maxsplit=1)
check(name == "MANIFEST.sha256", "target seal names local manifest")
check(expected == digest(TARGET / "MANIFEST.sha256"), "target seal closes")

# Recheck every transitive dependency selected by the author packet.
dependency_paths = set()
for line in (TARGET / "DEPENDENCIES.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    candidate = ROOT / relative
    check(candidate.is_file(), f"target dependency exists: {relative}")
    check(digest(candidate) == expected, f"target dependency hash: {relative}")
    check(relative not in dependency_paths, f"target dependency unique: {relative}")
    dependency_paths.add(relative)
check(len(dependency_paths) == 21, "exact target dependency count")


# Confirm both author replays on the frozen bytes.  They are custody checks;
# the independent mathematical verdict below does not use their results.
author = subprocess.run(
    [sys.executable, "-B", str(TARGET / "derive_full_pair_rotational_completion.py")],
    cwd=ROOT, text=True, capture_output=True, timeout=90, check=False)
check(author.returncode == 0, "author science replay exits zero")
for token in (
    "PASS__GL6CR_FULL_PAIR_ROTATIONAL_COMPLETION__599/599",
    "QUADRATIC_DIMENSIONS=S4_9;SO3_4;FIVE_MATCHING_CONDITIONS",
    "DIRECT_CUBIC_GAUGE_NULL=RANK8;UNIQUE_EINSTEIN_RAY;SO3_FOLLOWS",
):
    check(token in author.stdout, f"author science replay token: {token}")

author_packet = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_packet.py")],
    cwd=ROOT, text=True, capture_output=True, timeout=90, check=False)
check(author_packet.returncode == 0, "author packet replay exits zero")
check("PASS__GL6CR_PACKET__129/129" in author_packet.stdout,
      "author packet replay token")


# Independent exact result and fresh independent replay.
result = strict_json(HERE / "INDEPENDENT_RESULT.json")
check(result["audit"] == "GL6CR independent hostile audit",
      "independent result identity")
check(result["method"] ==
      "independent exact rational reconstruction; author module neither imported nor executed",
      "independent method identity")
check(result["disposition"] == "PASS", "independent disposition")
check(result["checks"] == 668, "independent mathematical check count")
check(result["reynolds_raw_seed_count"] == 126, "all raw Reynolds seeds")
check(result["s4_quadratic_dimension_character"] == 9,
      "character dimension nine")
check(result["s4_quadratic_dimension_reynolds"] == 9,
      "Reynolds dimension nine")
check(result["s4_constant_dimension"] == 3, "constant S4 dimension three")
check(result["so3_quadratic_dimension"] == 4, "SO3 dimension four")
check(result["so3_constant_dimension"] == 2, "constant SO3 dimension two")
check(result["rotational_codimension"] == 5, "rotational codimension five")
check(result["direct_ward_constraint_rows"] == 180,
      "complete Ward row census")
check(result["direct_ward_rank"] == 8, "direct Ward rank eight")
check(result["direct_ward_nullity"] == 1, "direct Ward nullity one")
check(result["direct_ward_ray"] ==
      ["0", "2", "-2", "0", "-2", "-2", "4", "-1", "1"],
      "unique Ward ray")
check(result["einstein_coordinates"] ==
      ["0", "1/8", "-1/8", "0", "-1/8", "-1/8", "1/4", "-1/16", "1/16"],
      "Einstein Reynolds coordinates")
check(result["einstein_generic_rank"] == 3, "Einstein quotient rank three")
check(result["target_manifest_sha256"] == TARGET_MANIFEST_DIGEST,
      "independent result pins frozen target")
check("not derived" in result["ceiling"], "independent result physical ceiling")

independent = subprocess.run(
    [sys.executable, "-B", str(HERE / "independent_gl6cr_audit.py"),
     "--frozen-target-hash", TARGET_MANIFEST_DIGEST],
    cwd=ROOT, text=True, capture_output=True, timeout=90, check=False)
check(independent.returncode == 0, "fresh independent replay exits zero")
for token in (
    "PASS__GL6CR_INDEPENDENT_HOSTILE_REPLAY__668/668",
    "SPACE=S4_9;SO3_4;ROTATIONAL_CODIMENSION_5;T2_PROJECTION_RANK_1",
    "DIRECT_WARD=180x9;RANK_8;NULLITY_1;UNIQUE_EINSTEIN_RAY",
    "DISPOSITION=PASS;ALGEBRAIC_CLASSIFIER_ONLY;PHYSICAL_F3_WARD_1PI_GRAVITY_G_OPEN",
):
    check(token in independent.stdout, f"independent replay token: {token}")


# Claim typing: the exact algebraic implication is allowed; a physical Ward,
# Ricci, gravity, or G derivation is not.
theorem = " ".join((TARGET / "THEOREM.md").read_text().split())
target_result = " ".join((TARGET / "RESULT.md").read_text().split())
for token in (
    "The Ward identity is a target condition in this packet, not yet a derived property",
    "Equation (CR20) is a target classifier, not a derivation of the target from the parent",
    "complete same-state F3 1PI/quotient kernel must still be shown",
    "Imposing (CR22) by hand would merely insert the desired answer",
    "prove Ricci/Einstein dynamics, gravity, or `G`",
):
    check(token in theorem, f"target theorem claim ceiling: {token}")
check("This does not yet prove that the microscopic parent lands on that ray" in
      target_result, "target result physical ceiling")
check("must be derived from the parent's relational redundancy rather than imposed" in
      target_result, "target result Ward derivation ceiling")

audit = " ".join((HERE / "AUDIT_REPORT.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
for token in (
    "**PASS**",
    "21 x 6 = 126` raw symmetric-pair/quadratic-momentum seeds",
    "complete `180 x 9` exact constraint matrix",
    "exactly sixteen times the Reynolds-basis coordinates",
    "The target does not claim that the physical parent already satisfies the Ward identity",
    "does **not** establish the second bracket from F3",
    "No material defect remains on the declared algebraic surface",
):
    check(token in audit, f"audit report token: {token}")
check("Disposition: `PASS`" in readme, "README disposition")
check("does not derive the Ward identity from F3" in readme,
      "README physical ceiling")
check("PASS__GL6CR_INDEPENDENT_HOSTILE_REPLAY__668/668" in verification,
      "verification records independent replay")
check("PASS__AUDIT_G_GL6CR_PACKET__" in verification,
      "verification records audit packet replay")

lower_docs = (audit + " " + readme).lower()
for forbidden in (
    "we have proved gravity",
    "f3 satisfies the ward identity",
    "derives the physical ward identity",
    "calculates the value of g",
    "is a graviton",
):
    check(forbidden not in lower_docs, f"forbidden promotion absent: {forbidden}")


# Audit manifest and one-file seal.
manifest_paths = set()
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    candidate = ROOT / relative
    check(candidate.is_file(), f"audit manifest file exists: {relative}")
    check(digest(candidate) == expected, f"audit manifest hash: {relative}")
    check(relative not in manifest_paths, f"audit manifest row unique: {relative}")
    manifest_paths.add(relative)
for name in required - {"MANIFEST.sha256", "SEAL.sha256"}:
    check(f"{HERE.name}/{name}" in manifest_paths,
          f"audit manifest covers: {name}")
check(len(manifest_paths) == 7, "exact audit manifest count")

seal_lines = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
              if line.strip()]
check(len(seal_lines) == 1, "one audit seal row")
expected, relative = seal_lines[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "audit seal target")
check(expected == digest(HERE / "MANIFEST.sha256"), "audit seal closes")

print(f"PASS__AUDIT_G_GL6CR_PACKET__{checks}/{checks}")
