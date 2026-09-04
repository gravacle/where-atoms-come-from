#!/usr/bin/env python3
"""Custody and scope verifier for the independent GL6CH hostile audit."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001"
checks = 0


def require(condition, label):
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
            raise ValueError(f"duplicate JSON member {key!r}")
        answer[key] = value
    return answer


def strict_json(path):
    return json.loads(path.read_text(), object_pairs_hook=unique_object)


required = {
    "AUDIT_REPORT.md", "INDEPENDENT_RESULT.json", "MANIFEST.sha256",
    "README.md", "SEAL.sha256", "TARGET.sha256", "VERIFICATION.txt",
    "verify_gl6ch_independent.py", "verify_packet.py",
}
for name in sorted(required):
    require((HERE / name).is_file(), f"required audit file {name}")

target_pins = [line for line in (HERE / "TARGET.sha256").read_text().splitlines()
               if line.strip()]
require(len(target_pins) == 12, "all twelve target bytes are pinned")
pin_names = set()
for line in target_pins:
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    require(path.is_file(), f"target byte exists {relative}")
    require(digest(path) == expected, f"target byte hash {relative}")
    pin_names.add(path.name)
require(pin_names == {path.name for path in TARGET.iterdir() if path.is_file()},
        "target pin set is complete and unique")

target_manifest_rows = {}
for line in (TARGET / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    target_manifest_rows[path.name] = expected
    require(path.parent == TARGET, f"target manifest path stays in packet {relative}")
    require(digest(path) == expected, f"target manifest hash {path.name}")
require(set(target_manifest_rows) == pin_names - {"MANIFEST.sha256", "SEAL.sha256"},
        "target manifest covers every non-custody byte")
target_seal = (TARGET / "SEAL.sha256").read_text().strip().split("  ", 1)
require(target_seal == [digest(TARGET / "MANIFEST.sha256"),
                        f"{TARGET.name}/MANIFEST.sha256"],
        "target seal closes and names target manifest")
require(target_seal[0] ==
        "a895ecdb1ab5340634808c0d6d379e96a0b161f8756d44ba2460ac5c404a34e5" and
        digest(TARGET / "SEAL.sha256") ==
        "e61abdab6ea225bb1c52c3d2f4e2050dc3dae4f70c0be1ba7712e369cc7ff61a",
        "requested stable target manifest and seal-file pins")

dependency_rows = {}
for line in (TARGET / "DEPENDENCIES.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    dependency_rows[relative] = expected
    require(path.is_file(), f"target dependency exists {relative}")
    require(digest(path) == expected, f"target dependency hash {relative}")
require(len(dependency_rows) == 12,
        "target has twelve unique AO/BX author-and-audit dependencies")

target_ledger = strict_json(TARGET / "EXACT_LEDGER.json")
require(target_ledger["lane"] == "GL6CH",
        "target ledger lane")
require(target_ledger["direct_hexagon"]["amplitude"] == "-63/8" and
        target_ledger["direct_hexagon"]["canonical_pair_gradient_coefficient"] ==
        "105/8" and
        target_ledger["direct_hexagon"]["tensor_vector_coefficient"] == "105/16",
        "target frozen direct coefficients")
require(target_ledger["lower_order_exclusion"]["h2_identity"] == "V2_v=-M_v" and
        target_ledger["lower_order_exclusion"]["h4_identity"] ==
        "V4_v=-(4/9)1_6-(37/12)M_v",
        "target frozen lower identities")
require(target_ledger["global_geometry"]["q4_four_cycles"] == 0 and
        target_ledger["global_geometry"]["q4_hexagons"] == 256 and
        target_ledger["global_geometry"]["orientation_rank"] == 3,
        "target frozen girth, cycle, and rank results")

replay = subprocess.run(
    [sys.executable, "-B", str(HERE / "verify_gl6ch_independent.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
require(replay.returncode == 0, "independent replay exits zero")
for token in (
    "PASS GL6CH independent hostile audit 41466/41466",
    "DIRECT=720x2;-63/8;FULL_GRADIENT=105/8_E_AB;T2=105/16_THETA",
    "LOWER=H0_H2_H4_T2_ZERO;DIFFERENTIATED_FINITE_STAR_486_CASES",
    "GEOMETRY=Q4_GIRTH6;256_ALL_SIMPLE_HEXAGONS;6_OWNERS_PER_EDGE",
    "ORIENTATION=SUM0;NORM24;INNER-8;GRAM32PT;RANK3",
    "DIMENSIONS_REMAINDERS_FOLDS=PASS;GENERIC_GRAPH_SCOPE=NARROW",
    "WORDING=CANDIDATE_FUTURE_WRITER_ONLY;NO_PHASE_RICCI_GRAVITY_G_GRAVITON",
):
    require(token in replay.stdout, f"independent replay token {token}")

independent = strict_json(HERE / "INDEPENDENT_RESULT.json")
require(independent["schema"] == "AUDIT_G_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001" and
        independent["verdict"] == "PASS" and
        independent["checks"] == 41464 and
        independent["material_defects"] == [],
        "frozen independent PASS verdict and count")
require(independent["independent_direct_history"] == {
            "canonical_gradient": "(105/8)e_ab",
            "energy_profiles": {
                "(2, 2, 2, 2, 2)": 96, "(2, 2, 2, 4, 2)": 48,
                "(2, 2, 4, 2, 2)": 48, "(2, 2, 4, 4, 2)": 96,
                "(2, 4, 2, 2, 2)": 48, "(2, 4, 2, 4, 2)": 24,
                "(2, 4, 4, 2, 2)": 96, "(2, 4, 4, 4, 2)": 192,
                "(2, 4, 6, 4, 2)": 72,
            },
            "literal_theta_derivative": "105/8",
            "local_contexts": 288,
            "orders_per_phase": 720,
            "phases": 2,
            "source_free": "-63/8",
            "t2_gradient": "(105/16)Theta_ab",
        }, "frozen independent direct-history result")
require(independent["independent_lower_orders"]["radius_one_neighborhoods"] == 486 and
        independent["independent_lower_orders"]["q_only_path_count_histogram"] ==
        {"480": 486} and
        independent["independent_lower_orders"]["h4"] ==
        "-(4/9)1_6-(37/12)M_v",
        "frozen differentiated finite-star lower-order result")
require(independent["independent_geometry"]["four_cycles"] == 0 and
        independent["independent_geometry"]["simple_six_cycles"] == 256 and
        independent["independent_geometry"]["canonical_cycles"] == 256 and
        independent["independent_geometry"]["cycles_per_edge"] == 6 and
        independent["independent_geometry"]["orientation_rank"] == 3,
        "frozen independent global geometry result")
require(set(independent["folding_dimensions_scope"]
            ["dimension_of_every_displayed_term"].values()) == {1},
        "frozen independent dimensions all equal energy")
require("self-generation and record authentication not proved" in
        independent["folding_dimensions_scope"]["writer_scope"] and
        "disciplined; no phase, metric, Ricci, gravity, G, or graviton promotion" ==
        independent["folding_dimensions_scope"]["interpretive_verdict"],
        "frozen candidate-writer interpretation ceiling")

report = " ".join((HERE / "AUDIT_REPORT.md").read_text().split())
for token in (
    "independent replay passes `41466/41466` exact checks",
    "all `480` four-flip identity words",
    "separate graph DFS obtains exactly `256` simple six-cycles",
    "all twelve ordered nonzero differences",
    "have energy dimension one",
    "`Arbitrary locked state` means every locked basis state of that parent",
    "not a new record-authentication theorem in GL6CH",
    "NO_PHASE_RECORD_AUTH_METRIC_RICCI_GRAVITY_G_OR_GRAVITON",
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

print(f"PASS GL6CH independent hostile audit packet {checks}/{checks}")
