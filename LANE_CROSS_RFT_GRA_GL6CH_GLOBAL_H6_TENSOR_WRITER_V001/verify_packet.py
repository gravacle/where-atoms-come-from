#!/usr/bin/env python3
"""Fail-closed custody, exact-ledger, and scope verifier for GL6CH."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition, label):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(label)


required = (
    "README.md",
    "THEOREM.md",
    "RESULT.md",
    "SELF_AUDIT.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "EXACT_LEDGER.json",
    "VERIFICATION.txt",
    "derive_global_h6_tensor_writer.py",
    "verify_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
)
for name in required:
    check((HERE / name).is_file(), f"required file: {name}")

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
check(len(dependency_paths) == 12, "exact dependency count")
check(all(("GL6AO" in path or "GL6BX" in path)
          for path in dependency_paths),
      "only GL6AO/GL6BX author and audit custody imported")
check(any(path.endswith("GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/THEOREM.md")
          for path in dependency_paths), "GL6AO theorem pinned")
check(any(path.endswith("GL6BX_EIGHTEEN_LINK_ACTIVE_COLLAR_FULL_PAIR_SOURCE_V001/THEOREM.md")
          for path in dependency_paths), "GL6BX theorem pinned")
check(any(path.startswith("AUDIT_G_GL6AO") and path.endswith("AUDIT.md")
          for path in dependency_paths), "GL6AO hostile audit pinned")
check(any(path.startswith("AUDIT_G_GL6BX") and
          path.endswith("INDEPENDENT_HOSTILE_AUDIT.md")
          for path in dependency_paths), "GL6BX hostile audit pinned")
check(all("GL6CG" not in path for path in dependency_paths),
      "unsealed GL6CG is not imported")

manifest_paths = set()
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"manifest target exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"manifest hash: {relative}")
    manifest_paths.add(relative)
for name in required:
    if name not in ("MANIFEST.sha256", "SEAL.sha256"):
        check(f"{HERE.name}/{name}" in manifest_paths,
              f"manifest coverage: {name}")
check(len(manifest_paths) == 10, "exact manifest count")

seal_lines = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
              if line.strip()]
check(len(seal_lines) == 1, "one seal row")
expected, relative = seal_lines[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest()
      == expected, "seal hash")

ledger = json.loads((HERE / "EXACT_LEDGER.json").read_text())
direct = ledger["direct_hexagon"]
geometry = ledger["global_geometry"]
lower = ledger["lower_order_exclusion"]
operator = ledger["operator"]
check(ledger["lane"] == "GL6CH", "ledger lane")
check(direct["amplitude"] == "-63/8", "source-free amplitude")
check(direct["canonical_pair_gradient_coefficient"] == "105/8",
      "canonical full pair gradient")
check(direct["tensor_vector_coefficient"] == "105/16",
      "T2 vector coefficient")
check(direct["theta_derivative"] == "105/8", "Theta directional derivative")
check(direct["context_count"] == 288, "complete local context count")
check(sum(row[1] for row in direct["energy_profile_histogram"]) == 720,
      "energy histogram exhausts 720 orders")
check(len(direct["energy_profile_histogram"]) == 9,
      "nine exact energy profiles")
check(sum(row[1] for row in direct["canonical_gradient_histogram"]) == 288,
      "canonical-gradient histogram exhausts contexts")
check(len(direct["canonical_gradient_histogram"]) == 6,
      "six cycle-pair canonical gradients")
check(lower["locked_local_words"] == 6, "six local locked words")
check(lower["radius_one_neighborhoods"] == 486,
      "complete radius-one census")
check(lower["h2_identity"] == "V2_v=-M_v", "h2 identity")
check(lower["h4_identity"] == "V4_v=-(4/9)1_6-(37/12)M_v",
      "h4 identity")
check(geometry["q4_hexagons"] == 256, "Q4 hexagon count")
check(geometry["q4_four_cycles"] == 0, "declared Q4 has no four-cycle")
check(geometry["hexagons_per_link"] == 6, "Q4 cycles per link")
check(geometry["orientation_rank"] == 3, "orientation writer rank")
check(sorted(row[1] for row in geometry["orientation_counts"]) == [64] * 4,
      "balanced orientation counts")
check(operator["source_free"] == "-(63/8)(h^6/U_d^5) sum_c T_c",
      "operator source-free term")
check(operator["first_T_source"] ==
      "+(105/16)(h^6/U_d^6) sum_c T_c sum_{v in c} j_v.Theta_{v,c}",
      "operator first-T term")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
for token in (
    "Theta_{v,c}=e_{ab}-e_{cd}",
    "g_v={105\\over8}e_{ab}",
    "P_Tg_v={105\\over16}\\Theta_{v,c}",
    "V_v^{(2)}=-M_v",
    "V_v^{(4)}=-{4\\over9}{\\bf1}_6-{37\\over12}M_v",
    "[H_{\\rm eff}(j_T)]_{\\rm off}^{(6)}",
    "\\sum_du_du_d^T=32P_T",
    "candidate-field-dependent future writer",
    "does **not** establish",
    "No graviton",
):
    check(token in theorem, f"theorem token: {token}")
check("five non-trace directions `E2+T2`" in result,
      "result retains five-shear access statement")
check("does not promote tensor source access" in readme,
      "README guards promotion")
check("tempting but false full-vector shortcut" in theorem,
      "hostile full-gradient correction recorded")
check("Hostile correction made during authorship" in self_audit,
      "self-audit records correction")
check("PASS__GL6CH_GLOBAL_H6_TENSOR_WRITER__60307/60307" in verification,
      "physics replay recorded")
check("PASS__GL6CH_PACKET__122/122" in verification,
      "packet replay recorded")

for forbidden in (
    "is a graviton",
    "is gravity",
    "proves gravity",
    "derives the Ricci tensor",
    "derives Newton's constant",
    "calculates G",
    "proves a stationary phase",
    "proves a causal cone",
):
    check(forbidden not in theorem + " " + result + " " + readme,
          f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6CH_PACKET__{checks}/{checks}")
