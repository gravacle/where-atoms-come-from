#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for GL6CM."""

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
    "VERIFICATION.txt", "verify_global_h6_stationary_writer_response.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required packet file {name}")

dependency_lines = [line for line in
                    (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
                    if line.strip()]
check(len(dependency_lines) == 21, "twenty-one exact dependency bytes pinned")
for line in dependency_lines:
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    check(path.is_file(), f"dependency exists {relative}")
    check(digest(path) == expected, f"dependency hash {relative}")

for packet in (
    "LANE_CROSS_RFT_GRA_GL6AR_LOCKED_HEXAGON_THERMODYNAMIC_SECTOR_V001",
    "AUDIT_G_GL6AR_LOCKED_HEXAGON_THERMODYNAMIC_SECTOR_V001",
    "LANE_CROSS_RFT_GRA_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001",
    "AUDIT_G_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001",
    "LANE_CROSS_RFT_GRA_GL6CJ_SAME_PARENT_SIX_DIRECTION_PAIR_OPERATOR_COMPOSITION_V001",
    "LANE_CROSS_RFT_GRA_GL6CK_TWO_OVERLAP_H6_TENSOR_WRITER_RESPONSE_V001",
    "AUDIT_G_GL6CK_TWO_OVERLAP_H6_TENSOR_WRITER_RESPONSE_V001",
):
    manifest = ROOT / packet / "MANIFEST.sha256"
    seal = (ROOT / packet / "SEAL.sha256").read_text().strip().split("  ", 1)
    check(seal[0] == digest(manifest) and
          seal[1] in {"MANIFEST.sha256", f"{packet}/MANIFEST.sha256"},
          f"upstream seal closes exact manifest {packet}")

ledger = json.loads((HERE / "EXACT_LEDGER.json").read_text(),
                    object_pairs_hook=unique_object)
check(ledger["schema"] == "GL6CM_GLOBAL_H6_STATIONARY_WRITER_RESPONSE_V001",
      "ledger schema")
check(ledger["verdict"] ==
      "PASS_WRITER_ONLY_SPECTRAL_RESPONSE_WITH_STRICT_TWO_OVERLAP_WITNESS",
      "ledger verdict")
check(ledger["spectral_response"] == {
          "formula": "K_spec(j,k)=2 lambda_T^2 Re <0|B_j Q(H0-E0)^-1 Q B_k|0>",
          "sign": "positive semidefinite",
          "reciprocity": True,
          "kernel": "K_spec(j,j)=0 iff Q B_j|0>=0",
          "common_cycle_rescaling": "null",
      }, "exact spectral response ledger")
check(ledger["two_overlap"]["star_curvature"] ==
      "-sqrt(2)(w0-w1)^2/(4J)" and
      ledger["two_overlap"]["literal_Q4_negative_energy_hessian"] ==
      "-(175sqrt(2)/32)h^6/U_d^7",
      "two-overlap response ledger")
check("possible diagonal order-h6 first-source vertex" in ledger["excluded"] and
      "record authentication or autonomous record generation of the source" in
      ledger["excluded"] and
      "metric, Ricci, Einstein, gravity, or Newton G" in ledger["excluded"],
      "scope exclusions ledger")

replay = subprocess.run(
    [sys.executable, "-B",
     str(HERE / "verify_global_h6_stationary_writer_response.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(replay.returncode == 0, "physics replay exits zero")
for token in (
    "PASS__GL6CM_EXACT_ALGEBRA__16/16",
    "SPECTRAL_RESPONSE=RECIPROCAL_PSD;KERNEL=QB|0>_ZERO",
    "COMMON_RING_RESCALING=NULL;ISOLATED_K2=NULL",
    "TWO_OVERLAP=175SQRT2_OVER32_H6_UD7_STRICT",
    "CONTACT_BULK_REALTIME_RICCI_GRAVITY_G=OPEN",
):
    check(token in replay.stdout, f"physics replay token {token}")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
for token in (
    "complete off-diagonal source-linear order-six writer block",
    "explicitly writer-only linear family",
    "possible diagonal order-six first-source vertex",
    "two-writer-vertex spectral contribution",
    "record authentication",
    "Gamma(phi)=s phi-W(s)=E(s)+s phi",
    "neither a Ricci coefficient nor `G`",
):
    check(token in theorem, f"theorem scope guard {token}")
check("stationary writer-sector response" in result and
      "record authentication" in result,
      "result states advance and open authentication")
check("not mislabeled as the full" in self_audit and
      "does not sum two-overlap answers" in self_audit,
      "self-audit guards omitted terms and accumulation")

manifest_lines = [line for line in
                  (HERE / "MANIFEST.sha256").read_text().splitlines()
                  if line.strip()]
manifest_names = set()
for line in manifest_lines:
    expected, name = line.split("  ", 1)
    manifest_names.add(name)
    check(Path(name).parent == Path("."), f"manifest local path {name}")
    check(digest(HERE / name) == expected, f"manifest hash {name}")
check(manifest_names == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "manifest covers every non-custody packet byte")
seal = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal == [digest(HERE / "MANIFEST.sha256"), "MANIFEST.sha256"],
      "seal names and hashes packet manifest")

print(f"PASS GL6CM packet {checks}/{checks}")
