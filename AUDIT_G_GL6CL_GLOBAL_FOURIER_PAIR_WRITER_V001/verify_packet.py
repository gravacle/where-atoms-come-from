#!/usr/bin/env python3
"""Fail-closed custody, independent replay, and scope audit for GL6CL."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CL_GLOBAL_FOURIER_PAIR_WRITER_V001"
checks = 0


def check(condition, label):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_hash_list(name):
    count = 0
    seen = set()
    for line in (HERE / name).read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        path = (HERE / relative.strip()).resolve()
        check(path.is_file(), f"{name} target exists: {relative}")
        check(digest(path) == expected, f"{name} hash: {relative}")
        check(relative not in seen, f"{name} has no duplicate: {relative}")
        seen.add(relative)
        count += 1
    return count


target_count = verify_hash_list("TARGET.sha256")
context_count = verify_hash_list("CONTEXT.sha256")
check(target_count == 12, "twelve target bytes pinned")
check(context_count == 8, "eight GL6CH target/audit context bytes pinned")

target_science = subprocess.run(
    [sys.executable, "-B", str(TARGET / "derive_global_fourier_pair_writer.py")],
    cwd=TARGET, check=True, text=True, stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
)
check("PASS__GL6CL_GLOBAL_FOURIER_PAIR_WRITER__122/122" in target_science.stdout,
      "fresh target science replay passes 122/122")
target_packet = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_packet.py")],
    cwd=TARGET, check=True, text=True, stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
)
check("PASS__GL6CL_PACKET__129/129" in target_packet.stdout,
      "fresh target packet verifier passes 129/129")

independent = subprocess.run(
    [sys.executable, "-B", str(HERE / "verify_gl6cl_independent.py")],
    cwd=HERE, check=True, text=True, stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
)
for token in (
    "PASS GL6CL independent geometry/offset moments: 4 orientations, 12 offsets",
    "PASS zero-mode scopes/ranks: D=A1+E rank3; BPT=T2 rank3; combined rank6 det524288",
    "PASS common expansion and smooth ball: |k|^4<32/363",
    "PASS T2 rotational guard: full E2+T2 completion required",
    "PASS relative sector: zero0/generic3/face-leading2; six exact dependencies",
    "PASS finite-momentum loss: common tensor rank1 at q=(pi,0,0)",
    "PASS uniform A1 storage/writer identity: 315/4",
    "AUDIT DISPOSITION: PASS",
):
    check(token in independent.stdout, f"independent output: {token}")

result = json.loads((HERE / "INDEPENDENT_RESULT.json").read_text())
check(result["schema"] == "AUDIT_G_GL6CL_INDEPENDENT_V001" and
      result["disposition"] == "PASS", "frozen independent schema/disposition")
geometry = result["geometry"]
check(geometry["orientations"] == 4 and geometry["centered_offsets"] == 12 and
      geometry["all_radius_squared"] == "11/4", "orientation/offset census")
check(geometry["second_moment"] ==
      [["11", "0", "0"], ["0", "11", "0"], ["0", "0", "11"]],
      "isotropic aggregate second moment")
check(geometry["fourth_diagonal"] == ["83/4"] * 3 and
      geometry["fourth_mixed"] == ["19/4"] * 3,
      "exact fourth offset moments")

zero = result["zero_mode"]
check(zero["canonical_direct_rank"] == 4 and
      zero["canonical_direct_normal"] == "24 P_A + 8 P_T",
      "canonical-direct rank-four bookkeeping row")
check(zero["complete_tensor_writer_rank"] == 3 and
      zero["complete_tensor_writer_normal"] == "8 P_T",
      "complete projected tensor writer rank three")
check(zero["locked_read_rank"] == 3 and
      zero["locked_read_normal"] == "4 P_A + 16 P_E",
      "locked read alone supplies A1+E")
check(zero["combined_rank"] == 6 and zero["combined_determinant"] == "524288" and
      zero["combined_normal"] == "4 P_A + 16 P_E + 8 P_T",
      "typed combined map has rank six and exact determinant")
check(zero["unsoldered_rank"] == 9 and
      zero["unsoldered_all_momentum_rank_ceiling"] == 10,
      "unsoldered zero rank and all-momentum ceiling")

scope = result["writer_scope"]
check(scope["complete_arbitrary_profile"] == "B_plus P_T and B_minus P_T only" and
      scope["A1_plus_E_access"] == "locked diagonal read D" and
      scope["unprojected_canonical_direct_row"] == "bookkeeping only",
      "writer/read/bookkeeping scope is exact")
check("autonomous source" in scope["not_established"] and
      "gravity" in scope["not_established"], "downstream claims remain open")

common = result["common_expansion"]
check(common["canonical_direct_AT_determinant"].startswith("768-1408|k|^2+1072|k|^4"),
      "canonical-direct determinant expansion")
check(common["orthonormal_T_normal"] ==
      "8I-2|k|^2I+12kk^T-28diag(k_i^2)+O(|k|^4)",
      "exact CL23 T2 block")
check("neither physical isotropy nor anisotropy" in common["scope_conclusion"] and
      "not SO(3)-closed" in common["T2_SO3_scope"],
      "CL23 is not promoted to an SO3 diagnosis")

smooth = result["smooth_ball"]
check(smooth["frobenius_bound_squared_coefficient"] == "363/4" and
      smooth["minimum_zero_T_singular_value_squared"] == "8" and
      smooth["full_rank_condition"] == "|k|^4<32/363",
      "rigorous smooth-ball constants")
check(smooth["analytic_left_inverse"] == "[C(k)^*C(k)]^-1 C(k)^*",
      "analytic left inverse on smooth ball")

relative = result["relative_sector"]
check(relative["zero_rank"] == 0 and relative["generic_leading_rank"] == 3 and
      relative["face_diagonal_leading_rank"] == 2,
      "relative leading ranks")
check(relative["exact_face_diagonal_column_dependencies"] == 6 and
      "58 kx^2 ky^2 kz^2" in relative["minor_sum"],
      "six exact relative dependencies and minor sum")

finite = result["finite_momentum"]
check(finite["q"] == "(pi,0,0)" and
      finite["complete_common_T_rank"] == 1 and
      finite["all_three_columns_identical"] is True,
      "exact finite-momentum rank loss")

scalar = result["uniform_A1"]
check(scalar["identity"] == "sum_{a<b} Z_a Z_b=2(n-2)^2-2" and
      scalar["occupancy_values"] == {"0": 6, "1": 0, "2": -2, "3": 0, "4": 6},
      "uniform A1 identity on all occupancy sectors")
check(scalar["denominator_derivative"] == scalar["six_vertex_sum"] ==
      "+315/4 h^6/U_d^6", "uniform denominator/six-vertex coefficient match")

ledger = json.loads((TARGET / "EXACT_LEDGER.json").read_text())
check(ledger["exact_symbol"]["complete_tensor_writer"] ==
      "delta a_d^T=mu[B_plus P_T j_plus+B_minus P_T j_minus]",
      "target ledger limits complete writer to projected tensor source")
check("arbitrary-profile A1/E" in ledger["exact_symbol"]["scope_guard"] and
      "unclassified" in ledger["exact_symbol"]["scope_guard"],
      "target ledger guards arbitrary-profile A1/E completion")
check(ledger["locked_read_and_common_writer"]["C_definition"] ==
      "C(k)=(D;B_+(k)P_T)", "target composition uses D plus projected writer")
check("cannot distinguish physical rotational anisotropy" in
      ledger["common_sector"]["rotation_ceiling"], "target ledger guards SO3 claim")

theorem = " ".join((TARGET / "THEOREM.md").read_text().split())
readme = " ".join((TARGET / "README.md").read_text().split())
result_text = " ".join((TARGET / "RESULT.md").read_text().split())
for token in (
    "without `P_T` are canonical-direct bookkeeping",
    "The arbitrary-profile `A1/E` off-diagonal completion is not classified",
    "`A1+E` through the locked read and `T2` through the complete writer",
    "does **not**, by itself, prove physical rotational anisotropy or isotropy",
    "rank exactly nine",
    "rank one there",
    "same-parent storage-energy/future-writer linkage",
):
    check(token in theorem, f"target theorem boundary: {token}")
check("The locked read supplies `A1+E`; the complete writer supplies `T2`" in readme,
      "README carries typed access split")
check("complete arbitrary-profile writer claim is rank three on `T2`" in result_text,
      "result carries projected-writer scope")

report = " ".join((HERE / "AUDIT_REPORT.md").read_text().split())
for token in (
    "**Disposition:** **PASS**",
    "complete arbitrary-profile order-six writer is only",
    "D^*D=4P_A+16P_E",
    "determinant `524288`",
    "|k|^4<{32\\over363}",
    "cannot by itself diagnose physical rotational anisotropy",
    "rank one at this point",
    "315\\over4",
):
    check(token in report, f"audit report token: {token}")

verification = (HERE / "VERIFICATION.txt").read_text()
for token in ("DISPOSITION: PASS", "524288", "32/363", "SO(3)", "rank loss to 1", "315/4"):
    check(token in verification, f"verification summary token: {token}")

manifest_count = verify_hash_list("MANIFEST.sha256")
check(manifest_count == 8, "eight audit payloads in manifest")
seal_lines = [line for line in (HERE / "SEAL.sha256").read_text().splitlines() if line.strip()]
check(len(seal_lines) == 1, "one audit seal row")
expected, name = seal_lines[0].split(maxsplit=1)
check(name == "MANIFEST.sha256" and expected == digest(HERE / name),
      "audit seal authenticates manifest")

print(f"PASS target custody: {target_count}/{target_count}")
print(f"PASS GL6CH context custody: {context_count}/{context_count}")
print("PASS fresh target replays: science 122/122; packet 129/129")
print("PASS independent GL6CL physics and hostile scope checks")
print(f"PASS audit manifest: {manifest_count}/{manifest_count}")
print(f"PASS total checks: {checks}/{checks}")
print("AUDIT DISPOSITION: PASS")
