#!/usr/bin/env python3
"""Fail-closed verifier for the frozen GL6CU V002 PASS audit."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V002"
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
    "README.md", "AUDIT_REPORT.md", "INDEPENDENT_RESULT.json",
    "TARGET.sha256", "REFERENCES.sha256", "VERIFICATION.txt", "independent_reconstruction.py",
    "verify_audit.py", "verify_target_custody_from_completed_replay.py",
    "MANIFEST.sha256", "SEAL.sha256",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required audit file {name}")


target_lines = [line for line in (HERE / "TARGET.sha256").read_text().splitlines()
                if line.strip()]
check(len(target_lines) == 12, "twelve final target bytes")
target_names = set()
for line in target_lines:
    expected, relative = line.split("  ", 1)
    check(len(expected) == 64, f"target SHA width {relative}")
    check(relative not in target_names, f"unique target byte {relative}")
    target_names.add(relative)
    candidate = ROOT / relative
    check(candidate.is_file(), f"target exists {relative}")
    check(digest(candidate) == expected, f"target hash {relative}")
check(target_names == {
    f"{TARGET.name}/{name}" for name in (
        "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
        "DEPENDENCIES.md", "DEPENDENCIES.sha256", "EXACT_LEDGER.json",
        "VERIFICATION.txt", "derive_complete_six_pair_h6_source_jet.py",
        "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256")
}, "complete target byte set")


result = json.loads((HERE / "INDEPENDENT_RESULT.json").read_text(),
                    object_pairs_hook=unique_object)
check(result["disposition"] == "PASS", "audit PASS disposition")
check(result["independent_checks"] == 2769, "independent check count")
check(set(result["v001_repairs"]) == {
    "packet_verifier", "support_scope", "compact_ledger", "ctp_custody"
}, "four exact V001 repairs")
check(result["independent_reconstruction"]["literal_stencil_records_tested"] == 8,
      "eight literal stencils independently tested")
check(result["independent_reconstruction"]["mixed_hessian_polarization_per_literal_stencil"] is True,
      "mixed Hessian independently tested")
check(result["ctp_and_half_source"]["physical_retarded_direct_contact"] ==
      "-(E_star^2/2)K'' delta(t-s)", "negative half-source pullback")
check(result["compact_ledger_alone_is_full_coefficient_export"] is False and
      result["sealed_executable_reconstructs_all_committed_classes"] is True,
      "compact ledger and executable boundary")
for key in ("stationary_connected_response_earned", "one_PI_quotient_earned",
            "ward_null_earned", "gl6cr_invocation_allowed",
            "einstein_shape_earned", "gravity_earned", "C_R_earned",
            "G_earned"):
    check(result[key] is False, f"ceiling {key}")


independent = subprocess.run(
    [sys.executable, "-B", str(HERE / "independent_reconstruction.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(independent.returncode == 0, "independent reconstruction exits zero")
for token in (
    "PASS__INDEPENDENT_GL6CU_V002_HOSTILE_RECONSTRUCTION__2769/2769",
    "REPAIRS=V001_R1_R2_R3_R4_CLOSED",
    "CENSUS=4096_ENUMERATED__3904_COMPUTED__192_INACTIVE",
    "ARITHMETIC=EIGHT_LITERAL_STENCILS_INDEPENDENTLY_SCALAR_PULLBACK_TESTED",
    "HESSIAN=ACTUAL_SECOND_DERIVATIVE_AND_MIXED_POLARIZATION_PASS",
    "CTP=-SIGMA_KPP__HALF_SOURCE_RETARDED=-E2_OVER_2_KPP",
    "CEILING=OPERATOR_SOURCE_JET_ONLY__NO_CONNECTED_1PI_WARD_GRAVITY_OR_G",
):
    check(token in independent.stdout, f"independent output {token}")


# The audit-required exact 21849-check replay is frozen in VERIFICATION.txt.
# Reuse that one result while executing every other assertion in the unmodified
# target packet verifier; this avoids repeating a 407-second sealed ledger.
target_verify = subprocess.run(
    [sys.executable, "-B", str(HERE / "verify_target_custody_from_completed_replay.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(target_verify.returncode == 0, "target custody harness exits zero")
check("PASS__GL6CU_V002_PACKET__188/188" in target_verify.stdout,
      "target packet custody passes 188/188")
check("PASS__GL6CU_V002_CUSTODY_WITH_COMPLETED_REPLAY__ONE_REPLAY_REUSED" in
      target_verify.stdout, "exact completed replay reused once")


ledger = json.loads((TARGET / "EXACT_LEDGER.json").read_text(),
                    object_pairs_hook=unique_object)
check(ledger["schema"] == "GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V002" and
      ledger["status"] == "PASS" and ledger["checks"] == 21849,
      "target ledger schema status and exact check count")
census = ledger["owner_census"]
check(census["graph_wide_geometric_supports_enumerated"] == 4096 and
      census["selected_row_nonbare_operator_jets_computed"] == 3904 and
      census["inactive_six_cycle_supports_enumerated_not_differentiated"] == 192,
      "target exact support scope")
storage = ledger["stencil_storage_contract"]
check(storage["all_reconstructed_local_stencil_classes_committed"] == 1654 and
      storage["family_order_records"] == 8 and
      storage["literal_first_exact_stencils_stored"] == 8,
      "target compact storage census")
check(ledger["source_convention"]["ctp_branch_contact"].startswith(
      "direct same-branch connected contact is -sigma K''_AB delta(t-s)"),
      "target branchwise CTP sign")


# Repair provenance and the defining BV convention are exact dependencies.
dependency_lines = [line for line in
                    (TARGET / "DEPENDENCIES.sha256").read_text().splitlines()
                    if line.strip()]
check(len(dependency_lines) == 16, "sixteen target dependency bytes")
dependency_names = set()
for line in dependency_lines:
    expected, relative = line.split("  ", 1)
    check(relative not in dependency_names, f"unique dependency {relative}")
    dependency_names.add(relative)
    candidate = ROOT / relative
    check(candidate.is_file() and digest(candidate) == expected,
          f"dependency hash {relative}")
for token in ("GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001/THEOREM.md",
              "GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V001/MANIFEST.sha256",
              "AUDIT_G_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V001/AUDIT_REPORT.md"):
    check(any(token in name for name in dependency_names),
          f"repair/custody dependency {token}")

# The downstream negative-half-source statement is tied to the frozen GL6W
# theorem and its independent hostile audit, not introduced by this audit.
reference_lines = [line for line in
                   (HERE / "REFERENCES.sha256").read_text().splitlines()
                   if line.strip()]
check(len(reference_lines) == 5, "five GL6W reference custody bytes")
reference_names = set()
for line in reference_lines:
    expected, relative = line.split("  ", 1)
    check(relative not in reference_names, f"unique reference {relative}")
    reference_names.add(relative)
    candidate = ROOT / relative
    check(candidate.is_file() and digest(candidate) == expected,
          f"reference hash {relative}")
check(any(name.endswith("GL6W_Q4_NORMALIZED_OPERATOR_NEUTRAL_CTP_KERNEL_V001/THEOREM.md")
          for name in reference_names) and
      any(name.endswith("AUDIT_G_GL6W_Q4_NORMALIZED_OPERATOR_NEUTRAL_CTP_KERNEL_V001/AUDIT.md")
          for name in reference_names), "GL6W theorem and hostile audit pinned")


report = " ".join((HERE / "AUDIT_REPORT.md").read_text().split())
for token in (
    "PASS__REPAIRED_COMPLETE_SELECTED_BRANCH",
    "2769/2769", "21849/21849", "188/188",
    "4096_{\\rm geometric\\ supports}",
    "3904_{\\rm computed\\ selected-row\\ Jets}",
    "does **not** literally export all coefficient tables",
    "K^{{\\rm direct},\\sigma\\sigma}_{AB} =-\\sigma K''_{AB}\\delta(t-s)",
    "j_A=-{E_\\star\\over2}J_A",
    "{\\cal G}^{R,\\rm direct}=2W_{ar}^{\\rm direct}",
    "No stationary connected response, `1PI` kernel, physical Ward null",
):
    check(token in report, f"audit report token {token}")


verification = (HERE / "VERIFICATION.txt").read_text()
for token in (
    "PASS__INDEPENDENT_GL6CU_V002_HOSTILE_RECONSTRUCTION__2769/2769",
    "PASS__GL6CU_V002_COMPLETE_SIX_PAIR_H6_SOURCE_JET__21849/21849",
    "PASS__GL6CU_V002_PACKET__188/188",
    "FINAL_DISPOSITION=PASS_AT_OPERATOR_SOURCE_JET_SCOPE",
    "NO_PROMOTION=CONNECTED_1PI_WARD_GL6CR_EINSTEIN_GRAVITY_CR_G",
):
    check(token in verification, f"verification transcript token {token}")


manifest_lines = [line for line in (HERE / "MANIFEST.sha256").read_text().splitlines()
                  if line.strip()]
manifest_names = set()
for line in manifest_lines:
    expected, name = line.split("  ", 1)
    check(Path(name).parent == Path("."), f"manifest local path {name}")
    check(name not in manifest_names, f"manifest unique {name}")
    manifest_names.add(name)
    check((HERE / name).is_file(), f"manifest file exists {name}")
    check(digest(HERE / name) == expected, f"manifest hash {name}")
check(manifest_names == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "manifest covers all non-custody audit bytes")
seal = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal == [digest(HERE / "MANIFEST.sha256"), "MANIFEST.sha256"],
      "seal pins manifest")


print(f"PASS__AUDIT_GL6CU_V002_PACKET__{checks}/{checks}")
