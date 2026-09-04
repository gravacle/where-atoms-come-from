#!/usr/bin/env python3
"""Fail-closed custody, ledger, and scope verifier for GL6CJ."""

from __future__ import annotations

import hashlib
import json
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
    "derive_same_parent_pair_composition.py",
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
check(len(dependency_paths) == 9, "exact dependency count")
check(all(("GL6CH" in path or "GL6AV" in path)
          for path in dependency_paths),
      "only GL6CH and GL6AV custody imported")
check(any(path.endswith("GL6CH_GLOBAL_H6_TENSOR_WRITER_V001/THEOREM.md")
          for path in dependency_paths), "GL6CH theorem pinned")
check(any(path.endswith("GL6AV_RECORD_CONDITIONED_COLLECTIVE_METRIC_BRIDGE_V001/THEOREM.md")
          for path in dependency_paths), "GL6AV theorem pinned")
check(any(path.startswith("AUDIT_G_GL6AV") and path.endswith("AUDIT.md")
          for path in dependency_paths), "GL6AV hostile audit pinned")
check(all("EW_Q4" not in path for path in dependency_paths),
      "mutable EW is not imported")

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
read = ledger["diagonal_locked_read"]
writer = ledger["h6_tensor_writer"]
combined = ledger["combined"]
check(ledger["lane"] == "GL6CJ", "ledger lane")
check(ledger["single_source"] ==
      "H(j)=H(0)+sum_v j_v.M_v before Feshbach elimination",
      "single pre-Feshbach source")
check(read["rank"] == 3, "diagonal read rank")
check(len(read["locked_words"]) == 6, "six locked words")
check(len(read["unique_pair_words"]) == 3, "three complement classes")
check(sorted(row[1] for row in read["unique_pair_words"]) == [2, 2, 2],
      "two locked words per pair class")
check(writer["rank"] == 3, "writer rank")
check(writer["q4_cells"] == 64, "Q4 cells")
check(writer["q4_nodes"] == 128, "Q4 nodes")
check(writer["q4_links"] == 256, "Q4 links")
check(writer["q4_hexagons"] == 256, "Q4 hexagons")
check(writer["total_node_cycle_incidences"] == 1536,
      "Q4 node-cycle incidences")
check(writer["cycles_per_node"] == 12, "twelve cycles per node")
check(writer["cycles_per_local_pair"] == 2, "two cycles per local pair")
check(sorted(row[1] for row in writer["pair_counts"]) == [2] * 6,
      "representative pair census")
check(combined["rank"] == 6, "combined rank six")
check(combined["dressed_writer_scale"] ==
      "lambda_T=(105/16)h^6/U_d^6", "dressed writer scale")
check(combined["physical_amplitude_inverse"] ==
      "j_T=(2/105)(U_d^6/h^6) sum_c delta_a_c Theta_{v,c}",
      "physical coefficient inverse")
check(combined["operator_support"] ==
      "locked read diagonal; ring writer off-diagonal",
      "operator support separation")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
dependencies = " ".join((HERE / "DEPENDENCIES.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
for token in (
    "one microscopic source",
    "{\\cal D}^*{\\cal D}=\\sum_{s:\\,|s|=2}M(s)M(s)^T =4P_A+16P_E",
    "{\\cal W}_v^*{\\cal W}_v =\\sum_{c\\ni v}\\Theta_{v,c}\\Theta_{v,c}^T=8P_T",
    "exactly twelve elementary hexagons meet the node",
    "each of the six unordered local port pairs occurs exactly twice",
    "{\\cal C}_v^*{\\cal C}_v=4P_A+16P_E+8P_T",
    "R_C C_v=I_6",
    "algebraic source/read type split is closed at the operator-jet level",
    "AV-CONSTITUTIVE` and `AV-UPDATE` gates remain open",
    "does **not** prove",
):
    check(token in theorem, f"theorem token: {token}")
check("former **algebraic source/read type split" in theorem,
      "typed split qualification")
check("does not turn an external query field into an autonomous physical field"
      in readme, "README autonomy guard")
check("does not classify any order-six diagonal source term" in self_audit and
      "typed projections rather than a formal equality" in self_audit,
      "order-six completeness and typed-projection ceiling")
check("truncation identity, not an asymptotic equality" in theorem and
      "No equality for the complete first-source vertex at order six is asserted"
      in theorem,
      "CJ09 and CJ19 are dimensionally and operationally typed")
check("mutable `EW` packet is not a dependency" in dependencies,
      "EW dependency ceiling")
check("PASS__GL6CJ_SAME_PARENT_PAIR_COMPOSITION__4253/4253" in verification,
      "science replay recorded")
check("PASS__GL6CJ_PACKET__116/116" in verification,
      "packet replay recorded")
check("Autonomous field generation" in result and
      "record authentication" in result and
      "unclassified order-six vertex pieces" in result,
      "result retains physical and completeness ceilings")

for forbidden in (
    "is a metric",
    "is gravity",
    "proves Ricci",
    "proves the Einstein equation",
    "derives Newton's constant",
    "calculates G",
    "closes AV-CONSTITUTIVE",
    "closes AV-UPDATE",
    "autonomously generates j",
):
    check(forbidden not in theorem + " " + result + " " + readme,
          f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6CJ_PACKET__{checks}/{checks}")
