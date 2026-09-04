#!/usr/bin/env python3
"""Custody and scope verifier for the independent GL6CJ hostile audit."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CJ_SAME_PARENT_SIX_DIRECTION_PAIR_OPERATOR_COMPOSITION_V001"
checks = 0


def require(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON member {key!r}")
        result[key] = value
    return result


def strict_json(path):
    return json.loads(path.read_text(), object_pairs_hook=unique_object)


required = {
    "AUDIT_REPORT.md", "INDEPENDENT_RESULT.json", "MANIFEST.sha256",
    "README.md", "SEAL.sha256", "TARGET.sha256", "UPSTREAM.sha256",
    "VERIFICATION.txt", "verify_gl6cj_independent.py", "verify_packet.py",
}
for name in sorted(required):
    require((HERE / name).is_file(), f"required audit file {name}")

target_pins = [line for line in (HERE / "TARGET.sha256").read_text().splitlines()
               if line.strip()]
require(len(target_pins) == 12, "all twelve repaired-target bytes are pinned")
pin_names = set()
for line in target_pins:
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    require(path.is_file(), f"target byte exists {relative}")
    require(digest(path) == expected, f"target byte hash {relative}")
    pin_names.add(path.name)
require(pin_names == {path.name for path in TARGET.iterdir() if path.is_file()},
        "target pin set is complete and unique")

target_manifest_names = set()
for line in (TARGET / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    require(path.parent == TARGET, f"target manifest path stays in target {relative}")
    require(digest(path) == expected, f"target manifest hash {path.name}")
    target_manifest_names.add(path.name)
require(target_manifest_names == pin_names - {"MANIFEST.sha256", "SEAL.sha256"},
        "target manifest covers every non-custody byte")
target_seal = (TARGET / "SEAL.sha256").read_text().strip().split("  ", 1)
require(target_seal == [digest(TARGET / "MANIFEST.sha256"),
                        f"{TARGET.name}/MANIFEST.sha256"],
        "target seal closes and names repaired manifest")
require(target_seal[0] ==
        "20b0d3bd4b9fd89860dc2754e1a4a05ac6fb42017eec589a6400809123e4e993" and
        digest(TARGET / "SEAL.sha256") ==
        "5712bfd7745a4ffe660b0dc69dfae49a1a24ad157e11055222c47b06765a7466",
        "requested repaired target manifest and seal-file pins")

target_dependencies = {}
for line in (TARGET / "DEPENDENCIES.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    target_dependencies[relative] = expected
    require(path.is_file(), f"target dependency exists {relative}")
    require(digest(path) == expected, f"target dependency hash {relative}")
require(len(target_dependencies) == 9,
        "target has nine unique CH/AV author-and-audit dependencies")

upstream_lines = [line for line in (HERE / "UPSTREAM.sha256").read_text().splitlines()
                  if line.strip()]
require(len(upstream_lines) == 3, "three GL6CH hostile-audit custody bytes pinned")
for line in upstream_lines:
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    require(path.is_file(), f"upstream audit byte exists {relative}")
    require(digest(path) == expected, f"upstream audit byte hash {relative}")
upstream_seal = (ROOT / "AUDIT_G_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001/SEAL.sha256").read_text().strip().split("  ", 1)
require(upstream_seal == [
            digest(ROOT / "AUDIT_G_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001/MANIFEST.sha256"),
            "AUDIT_G_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001/MANIFEST.sha256",
        ], "GL6CH hostile-audit seal closes its manifest")

target_ledger = strict_json(TARGET / "EXACT_LEDGER.json")
independent = strict_json(HERE / "INDEPENDENT_RESULT.json")
require(target_ledger["lane"] == "GL6CJ" and
        independent["schema"] ==
        "AUDIT_G_GL6CJ_SAME_PARENT_SIX_DIRECTION_PAIR_OPERATOR_COMPOSITION_V001",
        "target and audit schemas")
require(independent["mathematical_verdict"] == "PASS" and
        independent["checks"] == 16456,
        "frozen independent mathematics PASS and count")
require(independent["projector_ranks"] == {"A1": 1, "E": 2, "T2": 3},
        "independent projector ranks 1+2+3")
require(independent["diagonal"]["rank"] == 3 and
        len(independent["diagonal"]["kernel"]) == 3 and
        independent["diagonal"]["normal"] ==
        target_ledger["diagonal_locked_read"]["normal"],
        "independent diagonal rank, kernel, and normal")
require(independent["writer"]["all_nodes_checked"] == 128 and
        independent["writer"]["canonical_hexagons"] == 256 and
        independent["writer"]["all_simple_hexagons"] == 256 and
        independent["writer"]["cycles_per_node"] == 12 and
        independent["writer"]["cycles_per_pair"] == 2,
        "independent complete Q4 writer incidence")
require(independent["writer"]["representative"]["rank"] == 3 and
        len(independent["writer"]["representative"]["kernel"]) == 3 and
        independent["writer"]["representative"]["normal"] ==
        target_ledger["h6_tensor_writer"]["normal"],
        "independent writer rank, kernel, and normal")
require(independent["writer"]["representative"]["combined_rank"] == 6 and
        independent["writer"]["representative"]["combined_normal"] ==
        target_ledger["combined"]["normal"],
        "independent combined rank-six normal")
require(independent["writer"]["amplitude_inverse"] ==
        "(2/105)(U_d^6/h^6)" and
        target_ledger["combined"]["physical_amplitude_inverse"] ==
        "j_T=(2/105)(U_d^6/h^6) sum_c delta_a_c Theta_{v,c}",
        "independent dressed-amplitude inverse normalization")
require("selected operator derivatives" in
        independent["operator_typing"]["combined_status"] and
        "endogenous or autonomous source field" in
        independent["operator_typing"]["not_established"],
        "independent operator-jet and autonomy boundary")

replay = subprocess.run(
    [sys.executable, "-B", str(HERE / "verify_gl6cj_independent.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
require(replay.returncode == 0, "independent mathematics replay exits zero")
for token in (
    "PASS GL6CJ independent mathematics 16458/16458",
    "DIAGONAL=A1_PLUS_E;RANK3;NORMAL_4PA_PLUS16PE;KERNEL_T2",
    "WRITER=T2;ALL_128_NODES;12_CYCLES;GRAM_8PT;KERNEL_A1_PLUS_E",
    "COMBINED=RANK6;EXACT_RECONSTRUCTION;SAME_PREFESHBACH_SOURCE",
    "STATUS=SELECTED_OPERATOR_JET_ONLY;NO_AUTONOMOUS_FIELD_RESPONSE_METRIC_GRAVITY",
):
    require(token in replay.stdout, f"independent replay token {token}")

theorem = " ".join((TARGET / "THEOREM.md").read_text().split())
result = " ".join((TARGET / "RESULT.md").read_text().split())
self_audit = " ".join((TARGET / "SELF_AUDIT.md").read_text().split())
for token in (
    r"V_{v,\rm diag}^{[0,2,4]}(j;s)",
    "Equation (CJ09) is a truncation identity, not an asymptotic equality",
    r"\tag{CJ19a}",
    r"\tag{CJ19b}",
    "No equality for the complete first-source vertex at order six is asserted",
    "diagonal order-six vertex and the `A1/E` off-diagonal order-six pieces remain unclassified",
    "operator-jet level",
    "not a reciprocal stationary response",
):
    require(token in theorem, f"repaired theorem scope token {token}")
require("unclassified order-six vertex pieces" in result and
        "gravity, and `G` remain outside this theorem" in result,
        "result retains incomplete-vertex and gravity ceilings")
require("not a dimensionally incomplete `O(r^6)` asymptotic formula" in self_audit and
        "typed projections rather than a formal equality with undefined complements" in self_audit,
        "self-audit records both hostile repairs")
for forbidden in (
    "+O(r^6)",
    "unclassified displayed-order complements",
    "is an autonomous bulk field",
    "is a reciprocal stationary response",
    "is a metric",
    "proves gravity",
):
    require(forbidden not in theorem, f"repaired theorem excludes {forbidden}")

report = " ".join((HERE / "AUDIT_REPORT.md").read_text().split())
for token in (
    "**Final verdict:** `PASS` after repair",
    "independent mathematical replay passes `16458/16458` exact checks",
    "bare dimensionless `O(r^6)`",
    "repaired `CJ19a` and `CJ19b`",
    r"\ker {\cal D}=T_2",
    r"\ker {\cal W}_v=A_1\oplus E",
    "all `128` constraint nodes",
    "six diagonal rows are basis-state operator entries",
    "only an operator-jet/source-access closure",
    "NO_COMPLETE_H6_VERTEX_AUTONOMOUS_FIELD_RESPONSE_PHASE_METRIC_RGRLB_RICCI_GRAVITY_OR_G",
):
    require(token in report, f"audit report proof/scope guard {token}")

manifest_lines = [line for line in (HERE / "MANIFEST.sha256").read_text().splitlines()
                  if line.strip()]
manifest_names = set()
for line in manifest_lines:
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    require(Path(relative).parent == Path(HERE.name),
            f"audit manifest path remains within packet {relative}")
    manifest_names.add(path.name)
    require(digest(path) == expected, f"audit manifest hash {path.name}")
require(manifest_names == required - {"MANIFEST.sha256", "SEAL.sha256"},
        "audit manifest covers every non-custody byte")
seal = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
require(seal == [digest(HERE / "MANIFEST.sha256"),
                 f"{HERE.name}/MANIFEST.sha256"],
        "audit seal closes and names audit manifest")

print(f"PASS GL6CJ independent hostile audit packet {checks}/{checks}")
