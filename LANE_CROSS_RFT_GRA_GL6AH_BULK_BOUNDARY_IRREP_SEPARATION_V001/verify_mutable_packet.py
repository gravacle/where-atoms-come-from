#!/usr/bin/env python3
"""Run GL6AH exact replays and enforce frozen custody/scope boundaries."""

import hashlib
import json
import pathlib
import shutil
import subprocess
import tempfile


HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent
checks = 0


def require(condition, message):
    global checks
    if not condition:
        raise AssertionError(message)
    checks += 1


required = {
    "THEOREM.md",
    "RESULT.md",
    "README.md",
    "SELF_AUDIT.md",
    "PRESCREEN_REQUEST.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "EXACT_IRREP_LEDGER.json",
    "verify_local_bulk_boundary.py",
    "verify_n1_connector_supports.cpp",
    "verify_mutable_packet.py",
    "VERIFICATION.txt",
    "MANIFEST.sha256",
    "SEAL.sha256",
}
for name in required:
    path = HERE / name
    require(path.is_file() and path.stat().st_size > 0,
            f"required mutable packet file: {name}")

require(not (HERE / "AUDIT.md").exists(), "post-freeze audit remains separate")

for target, audit in (
    ("LANE_CROSS_RFT_GRA_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001",
     "AUDIT_G_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001"),
    ("LANE_CROSS_RFT_GRA_GL6AB_E2_MULTI_CONNECTOR_TRIANGLE_V001",
     "AUDIT_G_GL6AB_E2_MULTI_CONNECTOR_TRIANGLE_V001"),
    ("LANE_CROSS_RFT_GRA_GL6AF_FORMATION_PATTERN_E2_SOURCE_V001",
     "AUDIT_G_GL6AF_FORMATION_PATTERN_E2_SOURCE_V001"),
    ("LANE_CROSS_RFT_GRA_GL6AG_N1_FORMATION_NEIGHBOR_PROPAGATION_V001",
     "AUDIT_G_GL6AG_N1_FORMATION_NEIGHBOR_PROPAGATION_V001"),
    ("LANE_CROSS_RFT_GRA_GL6Z_MULTICELL_CTP_K2_JET_V001",
     "AUDIT_G_GL6Z_MULTICELL_CTP_K2_JET_V001"),
):
    require((REPO / target / "THEOREM.md").is_file(), f"dependency {target}")
    audit_text = (REPO / audit / "AUDIT.md").read_text()
    require("**Audit verdict: PASS.**" in audit_text, f"audited dependency {audit}")

local = subprocess.run(
    ["python3", "-B", str(HERE / "verify_local_bulk_boundary.py")],
    check=True, capture_output=True, text=True,
)
require("PASS__GL6AH_LOCAL_BULK_BOUNDARY__177/177" in local.stdout,
        "local exact replay count")
for token in (
    "DIRECT_Q6=-128_H4_UD2_U_B_IN;EXPECTATION=+8/45",
    "TAG_POLYNOMIAL=-128_G2",
    "DIRECT_IRREPS=EQUAL_A1_T2;E_EXACT_NULL",
    "CHAIN_Q12=+63371264_W_BD;CHAIN_Q16=-123422773248_W_BD",
    "ALTERNATE_NO_DIRECT=ZERO_BELOW_Q14;Q14=-14721024_U_1_IN;E_NULL",
    "HOMOGENEOUS_ETA_ZERO",
):
    require(token in local.stdout, f"local replay token {token}")

cxx = shutil.which("c++")
require(cxx is not None, "C++ compiler available")
with tempfile.TemporaryDirectory(prefix="gl6ah-") as temporary:
    binary = pathlib.Path(temporary) / "verify_gl6ah_supports"
    command = [cxx, "-O3", "-std=c++17"]
    if pathlib.Path("/opt/homebrew/include/gmpxx.h").is_file():
        command.extend(["-I/opt/homebrew/include", "-L/opt/homebrew/lib"])
    command.extend([
        str(HERE / "verify_n1_connector_supports.cpp"),
        "-lgmpxx", "-lgmp", "-o", str(binary),
    ])
    subprocess.run(command, check=True, capture_output=True, text=True)
    support = subprocess.run(
        [str(binary)], check=True, capture_output=True, text=True,
    )
require("PASS__GL6AH_N1_CONNECTOR_SUPPORTS__1293/1293" in support.stdout,
        "full connector census count")
for token in (
    "Q6_SUPPORTS=1:DIRECT_01:-128_U0_IN",
    "Q12_SUPPORTS=2:CHAIN_01_12:+63371264_W02;CHAIN_01_13:+63371264_W03",
    "Q16_SUPPORTS=2:CHAIN_01_12:-123422773248_W02;CHAIN_01_13:-123422773248_W03",
    "FULL_Q12=-63371264_W01;FULL_Q16=+123422773248_W01",
):
    require(token in support.stdout, f"support replay token {token}")

ledger = json.loads((HERE / "EXACT_IRREP_LEDGER.json").read_text())
require(ledger["status"].startswith("author frozen"), "frozen ledger")
require(ledger["direct_edge"]["raw_order_6"] ==
        "-128 h^4 Ud^2 u_b^in", "ledger q6")
require(len(ledger["N1_connector_support_census"]["q6_nonzero_supports"]) == 1,
        "ledger one q6 support")
require(len(ledger["N1_connector_support_census"]["q12_nonzero_supports"]) == 2,
        "ledger two q12 supports")
require(len(ledger["N1_connector_support_census"]["q16_nonzero_supports"]) == 2,
        "ledger two q16 supports")
require(ledger["scope"]["physical_K_not_semantic_REC"], "K/REC scope")
require(not ledger["scope"]["bulk_shear_claimed"], "no bulk shear")
require(not ledger["scope"]["stationary_mode_claimed"], "no stationary mode")
require(not ledger["scope"]["common_cone_claimed"], "no common cone")
require(not ledger["scope"]["Ricci_or_gravity_or_G_claimed"],
        "no Ricci/gravity/G")

theorem = (HERE / "THEOREM.md").read_text()
for token in (
    "author frozen after independent hostile pre-freeze review",
    "not inferred from the `N=1` cell graph `K4`",
    "-128h^4U_d^2u_b^{\\rm in}",
    "{8\\over45}",
    "-128h^4U_d^2g^2u_b^{\\rm in}",
    "u_b^{\\rm in}={1\\over2}\\mathbf1+{1\\over2}s_b",
    "E^Tx_q=0\\quad\\hbox{for every order }q",
    "+63371264\\,w_{bd}",
    "-123422773248\\,w_{bd}",
    "r_{n,d}:=\\#\\{c\\ne d:n_c>0\\}",
    "\\sum_{d\\ne b}w_{bd}=0",
    "term\nablations and tags are algebraic diagnostics",
    "not the semantic fact `REC`",
    "No stationary mode",
    "open detuning neighborhood",
    "AUTHOR_FROZEN_POSTFREEZE_AUDIT_REQUIRED",
):
    require(token in theorem, f"theorem scope token {token}")

for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = REPO / relative
    require(path.is_file(), f"dependency member exists: {relative}")
    require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
            f"dependency member hash: {relative}")

manifest_lines = (HERE / "MANIFEST.sha256").read_text().splitlines()
require(len(manifest_lines) == 12, "manifest member census")
for line in manifest_lines:
    expected, relative = line.split("  ", 1)
    path = REPO / relative
    require(path.is_file(), f"manifest member exists: {relative}")
    require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
            f"manifest member hash: {relative}")

seal_hash, seal_relative = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
require(seal_relative == f"{HERE.name}/MANIFEST.sha256", "seal target")
require(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest() ==
        seal_hash, "seal hash")

print(f"PASS__GL6AH_FROZEN_PACKET__{checks}/{checks}")
print(local.stdout, end="")
print(support.stdout, end="")
