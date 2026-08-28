#!/usr/bin/env python3
"""Final FY audit custody plus independent and frozen-packet replays."""

from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent
ROOT = LANE.parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def parse_manifest(path):
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line:
            expected, relative = line.split("  ", 1)
            rows.append((expected, relative))
    return rows


core = parse_manifest(AUDIT / "CORE_CUSTODY.sha256")
check(len(core) == 24 and len({relative for _, relative in core}) == 24,
      "core custody lists thirteen FY files and eleven dependencies exactly")
for expected, relative in core:
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(),
          f"core custody regular file: {relative}")
    check(digest(path) == expected, f"core custody hash: {relative}")

result = json.loads((AUDIT / "AUDIT_RESULT.json").read_text(encoding="utf-8"))
check(result["disposition"] == "PASS", "machine-readable disposition is PASS")
check(result["independent_checks"] == {"passed": 184, "total": 184},
      "machine-readable independent count is 184/184")
check(result["frozen_builder_checks"] == {
          "full_replay": "72/72", "packet_custody": "58/58",
          "exact_h2": "4/4"},
      "machine-readable frozen builder counts are complete")
check(result["remaining_material_defects"] == [],
      "machine-readable audit has no remaining material defect")
check(result["source_off_coefficients"] ==
      {"a2": "-60", "a4": "-35", "a6": "-893/9"},
      "machine-readable source-off coefficients are exact")
check(result["exact_m1_diagonal_lift"] ==
      {"H2": "-1", "H4": "-37/12", "H6": "-16247/900"},
      "machine-readable H2/H4/H6 lift is exact")
check(all(sample["ranks_operator_adH_spectral_static_M1_TT"] ==
          [6, 6, 6, 6, 6, 2] for sample in result["samples"]),
      "both machine-readable sampled rank packets are exact")
check(all(sample["pole_ranks"] == [1, 3, 1, 1]
          for sample in result["samples"]),
      "both machine-readable pole-rank packets are exact")
check(all(word in result["claim_ceiling"] for word in
          ("finite", "sampled", "Ward", "continuum", "massless",
           "RGRL-B", "gravity", "G", "Newton")),
      "machine-readable ceiling retains every nonpromotion boundary")

report = " ".join((AUDIT / "INDEPENDENT_HOSTILE_AUDIT.md").read_text(
    encoding="utf-8").split())
for phrase in ("**Disposition:** `PASS`", "184/184", "72/72", "58/58",
               "count-state path dynamic program",
               "formal native dual-series fixed point", "Phi_240",
               "All eighteen", "420", "fourteen", "m=29",
               "6 -> 6 -> 6 -> 6 -> 6", "sampled finite ranks", "Ward",
               "continuum locality", "massless graviton", "RGRL-B",
               "gravity emergence", "Newton's law"):
    check(phrase in report, f"audit report retains: {phrase}")

environment = dict(os.environ)
environment["PYTHONDONTWRITEBYTECODE"] = "1"
independent = subprocess.run(
    [sys.executable, str(AUDIT / "audit_native_support_m1.py")],
    cwd=LANE, env=environment, text=True, capture_output=True, check=False)
check(independent.returncode == 0, "independent hostile replay exits successfully")
check("SUMMARY 184/184 independent FY hostile-audit checks passed" in
      independent.stdout, "independent hostile replay reproduces 184/184")
check("orbit 5 exact H6 lift over Q(zeta_240)" in independent.stdout,
      "independent replay reaches the eighteenth exact lift gate")
check("VERDICT PASS; exact finite-graph native-support theorem only" in
      independent.stdout, "independent replay retains its finite ceiling")

packet = subprocess.run(
    [sys.executable, str(LANE / "verify_packet.py")], cwd=LANE,
    env=environment, text=True, capture_output=True, check=False)
check(packet.returncode == 0, "frozen FY packet replay exits successfully")
check("SUMMARY 58/58 FY packet checks passed" in packet.stdout,
      "frozen FY packet replay reproduces 58/58")
check("no Ward/locality/gravity claim" in packet.stdout,
      "frozen FY packet replay retains its scientific ceiling")

h2 = subprocess.run(
    [sys.executable, str(LANE / "verify_exact_m1_h2_lift.py")], cwd=LANE,
    env=environment, text=True, capture_output=True, check=False)
check(h2.returncode == 0, "frozen exact H2 replay exits successfully")
check("SUMMARY 4/4 exact H2 m=1 lift checks passed" in h2.stdout,
      "frozen exact H2 replay reproduces 4/4")
check("Qdiag2(m=1)=-Qpair(m=1) exactly" in h2.stdout,
      "frozen H2 replay retains the exact cyclotomic claim")

manifest = parse_manifest(AUDIT / "AUDIT_MANIFEST.sha256")
expected_files = {
    "AUDIT_RESULT.json", "CORE_CUSTODY.sha256",
    "INDEPENDENT_HOSTILE_AUDIT.md", "VERIFICATION.txt",
    "audit_native_support_m1.py", "verify_independent_audit.py",
}
check({relative for _, relative in manifest} == expected_files,
      "audit manifest lists exactly six frozen audit payloads")
for expected, relative in manifest:
    path = AUDIT / relative
    check(path.is_file() and not path.is_symlink(),
          f"audit manifest regular file: {relative}")
    check(digest(path) == expected, f"audit manifest hash: {relative}")

seal = parse_manifest(AUDIT / "AUDIT_SEAL.sha256")
check(seal == [(digest(AUDIT / "AUDIT_MANIFEST.sha256"),
                "AUDIT_MANIFEST.sha256")],
      "audit seal binds exactly the audit manifest")

allowed = {9, 10, 13}
for relative in expected_files:
    payload = (AUDIT / relative).read_bytes()
    check(not any(byte < 32 and byte not in allowed for byte in payload),
          f"audit payload has no forbidden control byte: {relative}")

print(f"SUMMARY {checks}/{checks} independent-audit custody checks passed")
print("FINAL PASS")
