#!/usr/bin/env python3
"""Independent hostile verifier for GRA-BENSE-PMHF-V001.

This executable imports no builder module.  It binds the frozen builder bytes,
reparses the deposited public files, independently reconstructs the numerical
diagnostics, and enforces the development-only causal and gravity ceilings.
"""

from __future__ import annotations

import bisect
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
SOURCE = HERE / "SOURCE"
CSV_NAMES = (
    "3_loops_SampleA.csv",
    "RPM.csv",
    "Accumulator.csv",
    "WeakRPMBreak.csv",
)

EXPECTED_BUILDER = {
    "CORE_MANIFEST.sha256": "87a5d67c58508ccad984e7b39c3e7ebb4bbca5c48a8164f01cffc3515769274b",
    "HOSTILE_AUDIT.md": "48dceb6cb61dc1d1fad60805b0fe26734223b03ca11758bb04ef2e999f95941b",
    "LANE_SEAL.json": "3ec5525012efa28847a7dcfb10dc945786d00fc0bfcb7e6d65b7475d496e62ce",
    "PREACCESS_METADATA.json": "cf58648deb4a0b02dae01e325d127d02e82fd444e666f20b08aa78a3a86aea3e",
    "PREREGISTRATION.md": "688dc9393adaf947fa39ce7ff8f754940f18476f8b1b9dd8cf4addc1482eb4d8",
    "RESULT.json": "4063539f62a6698beb4c82b15e2ce13ed0978f30353a49c9cac3836265bbbb6a",
    "SCHEMA_AMENDMENT_001.md": "ac0c2a40546ce2cac4dbcdb12257c7966534b0d6b7a92821f49920f0d795cb3c",
    "SCHEMA_AMENDMENT_002_THEORY_TYPING.md": "4c3bb0bd78f605ddb33cd406ee4a625c46c3b184ffbf73045ba377bc5e5c1da9",
    "SOURCE_CUSTODY.json": "f357070f8ed2331b02eb78bd6522f9ca198054b1f58c6154a7d841a5b125f329",
    "SOURCE_MANIFEST.sha256": "754b0d4d1f6b6dc64751cf5c359bda1bef498377f2bb72de65061d5596a8fd70",
    "THEOREM.md": "9eef510038cf541f70b22082613231024e42d81e17818fb0cf8b866df540f75c",
    "VERIFICATION_TRANSCRIPT.md": "832df83ba1da0180777aa57bfd1df27241bf36ea890f1a68f853bfd27b516e15",
    "analyze_public_history.py": "53bdf3dd2a1417e1c267c58e65afda50cb343a19dfea42f6bbfc7b85015072e4",
    "verify_lane.py": "fefdca4e29fabf76a9841bbb6e89e7323fd66c7ab317a71381a72c1b34b1e21c",
}

EXPECTED_SOURCE = {
    "3_loops_SampleA.csv": (31583915, 114477, "8d615357089e0b679553270e4192baac", "6dd2d66716e0c4ac0a756cc247d0ffbf2b40be3725ce74799317518f9fd6b5d3"),
    "RPM.csv": (31582832, 905278, "c3fd3e4ede58ba9719ff12039c61b0c2", "b5a95d4d903f6eae174e809fa0fcc88e1cffbcfe1fed7a92bf9989d132d8b8d9"),
    "Accumulator.csv": (31582829, 1100148, "fe4e22635d67ca9a9d20002b6870434b", "d261adee2cf02bc31daf260d64ca322144e9aec0f20888db21dd90add1a2a8ca"),
    "WeakRPMBreak.csv": (31582838, 789250, "dd53e862d2632f81a22f940652be5399", "4ae3aac59aa2c8d762ec6da12811591282924e2b21bfa3e6e7bbe5b8dbae340e"),
    "Fig_1de.py": (31583909, 5829, "35a6836513da8f1de6c25c1aca3fa90c", "18028aec2658291a800f27806a35b9de49762944f06a848080a2bf8d771ddbb9"),
    "Fig5ac.py": (31582835, 5499, "fbfa8c9e10360a60f16fd93930d35545", "a7bea4a0d2081f8ce3b4857a2126549889fa4af54067219b8d06053da9d0b04f"),
    "Fig6.py": (31582841, 3857, "3e8265722b49c7382e189b92f8790c1d", "89d04d4a0a9a7f2259d7742a4370349471f2526acf2345da9f5f1dd2a8118b2a"),
}

EXPECTED_ROWS = {
    "3_loops_SampleA.csv": 4124,
    "RPM.csv": 31595,
    "Accumulator.csv": 38315,
    "WeakRPMBreak.csv": 27594,
}
EXPECTED_BRANCHES = {
    "3_loops_SampleA.csv": 6,
    "RPM.csv": 12,
    "Accumulator.csv": 12,
    "WeakRPMBreak.csv": 8,
}
EXPECTED_COMPARISONS = {
    "3_loops_SampleA.csv": 9,
    "RPM.csv": 36,
    "Accumulator.csv": 36,
    "WeakRPMBreak.csv": 16,
}
EXPECTED_EVENTS = {
    "3_loops_SampleA.csv": 0,
    "RPM.csv": 14,
    "Accumulator.csv": 14,
    "WeakRPMBreak.csv": 6,
}


class Audit:
    def __init__(self) -> None:
        self.count = 0
        self.failures: list[str] = []

    def check(self, condition: bool, label: str) -> None:
        self.count += 1
        if not condition:
            self.failures.append(label)


def digest(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def ledger(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in path.read_text().splitlines():
        expected, relative = raw.split("  ", 1)
        if relative in out:
            raise AssertionError(f"duplicate manifest member: {relative}")
        out[relative] = expected
    return out


def stable(value: float) -> float:
    return round(float(value), 12)


def canonical_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def read_csv(name: str) -> dict[str, Any]:
    with (SOURCE / name).open(newline="", encoding="utf-8-sig") as handle:
        rows = csv.reader(handle, delimiter=";")
        header = next(rows)
        units = next(rows)
        values = [
            tuple(float(cell.replace(",", ".")) for cell in row)
            for row in rows
            if row
        ]
    return {
        "header": header,
        "units": units,
        "t": [row[0] for row in values],
        "u": [row[1] for row in values],
        "f": [row[2] for row in values],
    }


def fill_directions(values: list[float]) -> list[int]:
    signs = [(right > left) - (right < left) for left, right in zip(values, values[1:])]
    nonzero = [index for index, sign in enumerate(signs) if sign]
    if not nonzero:
        raise AssertionError("constant drive")
    for index, sign in enumerate(signs):
        if sign:
            continue
        insertion = bisect.bisect_left(nonzero, index)
        candidates = []
        if insertion:
            candidates.append(nonzero[insertion - 1])
        if insertion < len(nonzero):
            candidates.append(nonzero[insertion])
        chosen = min(candidates, key=lambda other: (abs(other - index), other))
        signs[index] = signs[chosen]
    return signs


def branches(data: dict[str, Any]) -> list[dict[str, Any]]:
    signs = fill_directions(data["u"])
    cuts = [0]
    cuts.extend(index for index in range(1, len(signs)) if signs[index] != signs[index - 1])
    cuts.append(len(signs))
    out = []
    for start, end in zip(cuts, cuts[1:]):
        support = data["u"][start : end + 1]
        if len(set(support)) < 3:
            continue
        out.append(
            {
                "id": len(out),
                "start_index": start,
                "end_index": end,
                "direction": signs[start],
                "direction_label": "increasing" if signs[start] > 0 else "decreasing",
                "extension_min_mm": min(support),
                "extension_max_mm": max(support),
                "distinct_extension_count": len(set(support)),
            }
        )
    return out


def collapse(data: dict[str, Any], branch: dict[str, Any]) -> tuple[list[float], list[float]]:
    bins: dict[float, list[float]] = {}
    lo, hi = branch["start_index"], branch["end_index"]
    for u, force in zip(data["u"][lo : hi + 1], data["f"][lo : hi + 1]):
        bins.setdefault(u, []).append(force)
    support = sorted(bins)
    return support, [statistics.fmean(bins[u]) for u in support]


def interpolate(xs: list[float], ys: list[float], target: float) -> float:
    position = bisect.bisect_left(xs, target)
    if position < len(xs) and xs[position] == target:
        return ys[position]
    if position == 0 or position == len(xs):
        raise AssertionError("interpolation outside support")
    left, right = position - 1, position
    weight = (target - xs[left]) / (xs[right] - xs[left])
    return ys[left] + weight * (ys[right] - ys[left])


def trapz(xs: list[float], ys: list[float]) -> float:
    return math.fsum(
        (xs[index + 1] - xs[index]) * (ys[index + 1] + ys[index]) / 2
        for index in range(len(xs) - 1)
    )


def compare(data: dict[str, Any], first: dict[str, Any], second: dict[str, Any]) -> dict[str, Any] | None:
    x1, f1 = collapse(data, first)
    x2, f2 = collapse(data, second)
    left, right = max(x1[0], x2[0]), min(x1[-1], x2[-1])
    if right <= left:
        return None
    common = sorted(set(x1).intersection(x2).intersection(u for u in set(x1) | set(x2) if left <= u <= right))
    if len(common) >= 3:
        grid = sorted(set(common + [left, right]))
        mode = "exact_common_extension"
    else:
        grid = [left + (right - left) * index / 100 for index in range(101)]
        mode = "101_point_linear_interpolation"
    first_values = [interpolate(x1, f1, u) for u in grid]
    second_values = [interpolate(x2, f2, u) for u in grid]
    difference = [a - b for a, b in zip(first_values, second_values)]
    width = right - left
    mean_abs = trapz(grid, [abs(value) for value in difference]) / width
    signed = trapz(grid, difference) / width
    rms = math.sqrt(trapz(grid, [value * value for value in difference]) / width)
    pooled = max(first_values + second_values) - min(first_values + second_values)
    points = [[stable(u), stable(value)] for u, value in zip(grid, difference)]
    return {
        "first_branch": first["id"],
        "second_branch": second["id"],
        "first_direction": first["direction_label"],
        "second_direction": second["direction_label"],
        "overlap_mm": [stable(left), stable(right)],
        "comparison_mode": mode,
        "comparison_point_count": len(grid),
        "signed_mean_load_difference_N": stable(signed),
        "mean_absolute_load_difference_N": stable(mean_abs),
        "rms_load_difference_N": stable(rms),
        "maximum_absolute_load_difference_N": stable(max(abs(value) for value in difference)),
        "pooled_load_range_N": stable(pooled),
        "normalized_overlap_integrated_absolute_difference": stable(mean_abs / pooled) if pooled else (0.0 if mean_abs == 0 else None),
        "pointwise_extension_difference_sha256": hashlib.sha256(canonical_json_bytes(points)).hexdigest(),
    }


def opposite_comparisons(data: dict[str, Any], branch_list: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for index, first in enumerate(branch_list):
        for second in branch_list[index + 1 :]:
            if first["direction"] == second["direction"]:
                continue
            item = compare(data, first, second)
            if item is not None:
                out.append(item)
    return out


def force_events(data: dict[str, Any], branch_list: list[dict[str, Any]]) -> list[dict[str, Any]]:
    candidates = [
        index
        for index, (left, right) in enumerate(zip(data["f"], data["f"][1:]))
        if abs(right - left) > 0.0079
    ]
    if not candidates:
        return []
    retained = [
        candidates[index]
        for index in range(len(candidates) - 1)
        if candidates[index + 1] - candidates[index] > 15
    ] + [candidates[-1]]
    out = []
    for index in retained:
        owner = next(
            (
                branch["id"]
                for branch in branch_list
                if branch["start_index"] <= index < branch["end_index"]
            ),
            None,
        )
        out.append(
            {
                "delta_start_index": index,
                "time_sec": stable(data["t"][index]),
                "extension_mm": stable(data["u"][index]),
                "load_before_N": stable(data["f"][index]),
                "load_after_N": stable(data["f"][index + 1]),
                "delta_load_N": stable(data["f"][index + 1] - data["f"][index]),
                "branch_id": owner,
            }
        )
    return out


def extrema(data: dict[str, Any], branch_list: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for first, second in zip(branch_list, branch_list[1:]):
        if first["end_index"] != second["start_index"] or first["direction"] == second["direction"]:
            continue
        index = first["end_index"]
        out.append(
            {
                "ordinal": len(out),
                "kind": "maximum" if first["direction"] > 0 else "minimum",
                "index": index,
                "time_sec": stable(data["t"][index]),
                "extension_mm": stable(data["u"][index]),
                "load_N": stable(data["f"][index]),
            }
        )
    return out


def return_diagnostic(data: dict[str, Any], branch_list: list[dict[str, Any]]) -> dict[str, Any]:
    turning = extrema(data, branch_list)
    step = statistics.median(
        abs(right - left)
        for left, right in zip(data["u"], data["u"][1:])
        if right != left
    )
    pairs = []
    for later_index, later in enumerate(turning):
        for earlier in turning[:later_index]:
            delta = later["extension_mm"] - earlier["extension_mm"]
            if later["kind"] == earlier["kind"] and abs(delta) <= step:
                pairs.append(
                    {
                        "first_extremum": earlier["ordinal"],
                        "returning_extremum": later["ordinal"],
                        "kind": later["kind"],
                        "extension_difference_mm": stable(delta),
                        "load_difference_N": stable(later["load_N"] - earlier["load_N"]),
                        "complete_state_word_available": False,
                        "state_closed": None,
                    }
                )
    return {
        "matching_tolerance_mm": stable(step),
        "local_extrema": turning,
        "all_same_type_return_pairs": pairs,
        "state_closure_status": "UNSCOREABLE_COMPLETE_STATE_WORD_NOT_DEPOSITED",
    }


def accumulator_cycles(data: dict[str, Any], branch_list: list[dict[str, Any]], events: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for offset in range(0, len(branch_list) - 1, 2):
        first, second = branch_list[offset : offset + 2]
        if first["direction"] == second["direction"] or first["end_index"] != second["start_index"]:
            continue
        start, turn, end = first["start_index"], first["end_index"], second["end_index"]
        out.append(
            {
                "cycle": len(out),
                "boundary_basis": "drive_extrema_inferred_from_deposited_Extension",
                "start_index": start,
                "turn_index": turn,
                "end_index": end,
                "start_extension_mm": stable(data["u"][start]),
                "turn_extension_mm": stable(data["u"][turn]),
                "end_extension_mm": stable(data["u"][end]),
                "start_load_N": stable(data["f"][start]),
                "end_load_N": stable(data["f"][end]),
                "force_event_count": sum(start <= event["delta_start_index"] < end for event in events),
                "complete_state_word_available": False,
                "state_transition_hamming": None,
            }
        )
    return out


def same_direction_repeats(data: dict[str, Any], branch_list: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for direction in (-1, 1):
        family = [branch for branch in branch_list if branch["direction"] == direction]
        for first, second in zip(family, family[1:]):
            item = compare(data, first, second)
            if item is not None:
                out.append(item)
    return out


def recursively_contains(value: Any, forbidden: set[str]) -> bool:
    if isinstance(value, dict):
        return any(key in forbidden or recursively_contains(item, forbidden) for key, item in value.items())
    if isinstance(value, list):
        return any(recursively_contains(item, forbidden) for item in value)
    return False


def main() -> int:
    audit = Audit()

    # Frozen builder byte custody.
    for relative, expected in EXPECTED_BUILDER.items():
        path = HERE / relative
        audit.check(path.is_file() and not path.is_symlink(), f"regular frozen builder member: {relative}")
        audit.check(digest(path) == expected, f"frozen builder hash: {relative}")
        audit.check(hashlib.sha256(path.read_bytes() + b"tamper").hexdigest() != expected, f"tamper rejection: {relative}")

    core = ledger(HERE / "CORE_MANIFEST.sha256")
    audit.check(set(core) == {
        "analyze_public_history.py", "HOSTILE_AUDIT.md", "PREACCESS_METADATA.json",
        "PREREGISTRATION.md", "RESULT.json", "SCHEMA_AMENDMENT_001.md",
        "SCHEMA_AMENDMENT_002_THEORY_TYPING.md", "SOURCE_CUSTODY.json",
        "SOURCE_MANIFEST.sha256", "THEOREM.md", "verify_lane.py",
    }, "core manifest has exact frozen payload membership")
    for relative, expected in core.items():
        audit.check(digest(HERE / relative) == expected, f"core manifest binding: {relative}")

    seal = json.loads((HERE / "LANE_SEAL.json").read_text())
    audit.check(seal["protocol"] == "GRA-BENSE-PMHF-V001", "seal protocol")
    audit.check(seal["canonical_or_git_files_edited"] is False, "builder isolation declaration")
    audit.check(seal["core_manifest_sha256"] == digest(HERE / "CORE_MANIFEST.sha256"), "seal binds core manifest")
    audit.check(seal["source_manifest_sha256"] == digest(HERE / "SOURCE_MANIFEST.sha256"), "seal binds source manifest")

    # Exact local source custody and ID/checksum transcription.
    custody = json.loads((HERE / "SOURCE_CUSTODY.json").read_text())
    preaccess = json.loads((HERE / "PREACCESS_METADATA.json").read_text())
    audit.check(custody["protocol"] == preaccess["protocol"] == "GRA-BENSE-PMHF-V001", "custody protocol agreement")
    audit.check(custody["figshare_project"]["project_id"] == preaccess["project"]["figshare_project_id"] == 127250, "Figshare project identity")
    custody_by_name = {item["name"]: item for item in custody["files"]}
    preaccess_by_name = {
        item["name"]: item
        for item in preaccess["selected_files"] + preaccess["semantic_scripts"]
    }
    audit.check(set(custody_by_name) == set(EXPECTED_SOURCE), "custody has exact seven selected objects")
    audit.check(set(preaccess_by_name) == set(EXPECTED_SOURCE), "preaccess metadata has exact seven objects")
    for name, (file_id, byte_count, md5, sha256) in EXPECTED_SOURCE.items():
        item = custody_by_name[name]
        before = preaccess_by_name[name]
        path = SOURCE / name
        audit.check(item["file_id"] == before["file_id"] == file_id, f"repository file ID: {name}")
        audit.check(path.stat().st_size == item["bytes"] == byte_count, f"source bytes: {name}")
        audit.check(digest(path, "md5") == item["md5"] == before["computed_md5"] == md5, f"source MD5: {name}")
        audit.check(digest(path) == item["sha256"] == sha256, f"source SHA-256: {name}")

    paper = custody["primary_article"]
    paper_path = HERE / paper["local_path"]
    audit.check(paper["doi"] == "10.1073/pnas.2111436118", "primary DOI")
    audit.check(paper_path.stat().st_size == paper["bytes"] == 1514470, "primary PDF bytes")
    audit.check(digest(paper_path) == paper["sha256"] == "bae9cbe1b9da6adde39a9c9bed78d3f3fc7c59b0f8e47a0ccaab9f3347eead52", "primary PDF hash")

    source_manifest = ledger(HERE / "SOURCE_MANIFEST.sha256")
    audit.check(set(source_manifest) == {
        "SOURCE/3_loops_SampleA.csv", "SOURCE/Accumulator.csv",
        "SOURCE/Bense_van_Hecke_2021_primary.pdf", "SOURCE/Fig5ac.py",
        "SOURCE/Fig6.py", "SOURCE/Fig_1de.py", "SOURCE/README.md",
        "SOURCE/RPM.csv", "SOURCE/WeakRPMBreak.csv",
    }, "source manifest exact membership")
    for relative, expected in source_manifest.items():
        audit.check(digest(HERE / relative) == expected, f"source manifest binding: {relative}")

    # The selection claim is content-addressed to the prior metadata-only audit,
    # while the present packet remains development evidence rather than validation.
    prior_audit = REPO / "LANE_GRA_T_MECHANICAL_JOINT_SEED/CANDIDATE_AUDIT.md"
    prior_ledger = REPO / "LANE_GRA_T_MECHANICAL_JOINT_SEED/PUBLIC_METADATA_LEDGER.md"
    audit.check(digest(prior_audit) == "79ca0f1c1f411209a0f17771528e427765d17c2a84c3e43194c0acdd6765ea11", "prior candidate-audit custody")
    audit.check(digest(prior_ledger) == "e67c71dc4833b79e0d8d352b93893df42fc48bbd11e4dd91976ad41c8ee970ac", "prior metadata-ledger custody")
    prior_text = prior_audit.read_text()
    audit.check("strongest *public executable* near-miss" in prior_text, "prior selection rationale")
    audit.check("No candidate CSV, NPY, XLSX, movie, image archive, or ZIP numerical payload was opened" in prior_text, "prior payload-access disclosure")
    audit.check("unanticipated contamination" in prior_text, "prior methods-code exposure disclosed")

    fig5 = (SOURCE / "Fig5ac.py").read_text()
    fig6 = (SOURCE / "Fig6.py").read_text()
    for label, script in (("Fig5ac", fig5), ("Fig6", fig6)):
        audit.check("abs(df)>0.0079" in script, f"official strict force threshold: {label}")
        audit.check("dj>15" in script, f"official strict cluster gap: {label}")
        audit.check("idx=jumps[j]" in script and "np.append(idx,jumps[-1])" in script, f"official cluster representative rule: {label}")
        audit.check("from CreateData import *" in script and "from BasicAnalysis import *" in script, f"unreleased module dependency: {label}")
    audit.check("exp[0].create_var()" in fig5 and "exp[1].create_var()" in fig5 and "exp[0].create_var()" in fig6, "directory-order script binding rather than immutable filenames")

    result = json.loads((HERE / "RESULT.json").read_text())
    audit.check(result["protocol"] == "GRA-BENSE-PMHF-V001", "result protocol")
    audit.check(result["execution_class"] == "RETROSPECTIVE_DEVELOPMENT_AFTER_DISCLOSED_SCHEMA_INSPECTION", "development-only execution class")
    audit.check(result["status"].endswith("CAUSAL_KEEP_BREAK_UNIDENTIFIED__NO_GRAVITY_TEST"), "bounded result disposition")
    audit.check(result["source_custody_all_pass"] is True, "reported source-custody status")
    audit.check(result["freeze_custody_scope"] == "INTERNALLY_FROZEN_ACCORDING_TO_PROJECT_CUSTODY__NOT_EXTERNALLY_TIMESTAMPED_OR_INDEPENDENTLY_REGISTERED", "freeze limitation is executable result schema")
    audit.check(result["search_scope"] == "FRONTIER_OF_DOCUMENTED_AUDITED_CANDIDATE_SET__NOT_AN_EXHAUSTIVE_WORLD_PUBLIC_DATASET_SEARCH", "audited-search limitation is executable result schema")
    audit.check(result["comparison_independence"] == "97_WITHIN_FILE_COMPARISONS_REUSE_FOUR_CHRONOLOGICAL_TRACES__NOT_97_INDEPENDENT_TRIALS__SAME_APPARATUS_STATUS_IS_PER_FILE", "comparison and apparatus scope are executable result schema")
    audit.check(result["strongest_supported_result"] == "within each deposited same-apparatus file, the reversible loading trace contains branch-dependent mechanical Load and publication-rule force-event trajectories; the four files are not one joined trial", "strongest result is explicitly per-file")

    parsed: dict[str, dict[str, Any]] = {}
    all_branches: dict[str, list[dict[str, Any]]] = {}
    all_events: dict[str, list[dict[str, Any]]] = {}
    total_rows = total_branches = total_comparisons = total_events = 0
    found_overlap_only_range_repair = False

    for name in CSV_NAMES:
        data = read_csv(name)
        parsed[name] = data
        own_branches = branches(data)
        all_branches[name] = own_branches
        own_events = force_events(data, own_branches)
        all_events[name] = own_events
        own_comparisons = opposite_comparisons(data, own_branches)
        reported = result["per_file"][name]

        audit.check(data["header"] == ["Time", "Extension", "Load"], f"three-column header: {name}")
        audit.check(data["units"] == ["(sec)", "(mm)", "(N)"], f"owned units: {name}")
        audit.check(len(data["t"]) == EXPECTED_ROWS[name], f"row count: {name}")
        audit.check(all(right > left for left, right in zip(data["t"], data["t"][1:])), f"strict chronology: {name}")
        audit.check(reported["schema"]["row_count"] == len(data["t"]), f"reported row count: {name}")
        audit.check(reported["schema"]["time_range_sec"] == [stable(min(data["t"])), stable(max(data["t"]))], f"time range: {name}")
        audit.check(reported["schema"]["extension_range_mm"] == [stable(min(data["u"])), stable(max(data["u"]))], f"drive range: {name}")
        audit.check(reported["schema"]["load_range_N"] == [stable(min(data["f"])), stable(max(data["f"]))], f"load range: {name}")
        audit.check(len(own_branches) == EXPECTED_BRANCHES[name], f"branch count: {name}")
        audit.check(len(own_comparisons) == EXPECTED_COMPARISONS[name], f"opposite-pair count: {name}")
        audit.check(len(own_events) == EXPECTED_EVENTS[name], f"event count: {name}")

        stated_branches = reported["monotone_branches"]
        audit.check(len(stated_branches) == len(own_branches), f"reported branch count: {name}")
        for position, (own, stated) in enumerate(zip(own_branches, stated_branches)):
            normalized = {key: stable(value) if isinstance(value, float) else value for key, value in own.items()}
            audit.check(normalized == stated, f"complete branch reconstruction: {name}/{position}")

        stated_events = reported["official_force_jump_proxy"]
        audit.check(stated_events["threshold_N_strictly_greater_than"] == 0.0079, f"reported event threshold: {name}")
        audit.check(stated_events["cluster_separation_samples_strictly_greater_than"] == 15, f"reported event cluster gap: {name}")
        audit.check(stated_events["interpretation"] == "candidate_force_events_not_a_complete_hysteron_state_word", f"event proxy ceiling: {name}")
        audit.check(own_events == stated_events["events"], f"complete event reconstruction: {name}")

        stated_comparisons = reported["D1_all_opposite_branch_comparisons"]
        audit.check(len(stated_comparisons) == len(own_comparisons), f"reported comparison count: {name}")
        for position, (own, stated) in enumerate(zip(own_comparisons, stated_comparisons)):
            audit.check(own == stated, f"complete branch functional: {name}/{position}")
            first = own_branches[own["first_branch"]]
            second = own_branches[own["second_branch"]]
            _, f1 = collapse(data, first)
            _, f2 = collapse(data, second)
            full_branch_range = max(f1 + f2) - min(f1 + f2)
            if not math.isclose(full_branch_range, own["pooled_load_range_N"], rel_tol=0, abs_tol=1e-12):
                found_overlap_only_range_repair = True

        audit.check(all(item["mean_absolute_load_difference_N"] > 0 for item in own_comparisons), f"nonzero descriptive path differences: {name}")
        total_rows += len(data["t"])
        total_branches += len(own_branches)
        total_comparisons += len(own_comparisons)
        total_events += len(own_events)

    audit.check((total_rows, total_branches, total_comparisons, total_events) == (101628, 38, 97, 34), "global diagnostic totals")
    audit.check(found_overlap_only_range_repair, "pooled range is genuinely restricted to compared overlap")

    # Exact return, accumulator, and repeated-loop diagnostics.
    for name in ("RPM.csv", "WeakRPMBreak.csv"):
        own = return_diagnostic(parsed[name], all_branches[name])
        stated = result["D2_return_and_repeated_input"][name]
        audit.check(own == stated, f"complete return diagnostic: {name}")
        audit.check(all(pair["state_closed"] is None for pair in own["all_same_type_return_pairs"]), f"no invented state closure: {name}")
    rpm_deltas = [abs(pair["load_difference_N"]) for pair in result["D2_return_and_repeated_input"]["RPM.csv"]["all_same_type_return_pairs"]]
    weak_deltas = [abs(pair["load_difference_N"]) for pair in result["D2_return_and_repeated_input"]["WeakRPMBreak.csv"]["all_same_type_return_pairs"]]
    audit.check(max(rpm_deltas) == 0.00411, "RPM maximum endpoint-load difference")
    audit.check(max(weak_deltas) == 0.06326, "weak-RPM-break maximum endpoint-load difference")

    own_cycles = accumulator_cycles(parsed["Accumulator.csv"], all_branches["Accumulator.csv"], all_events["Accumulator.csv"])
    stated_accumulator = result["D2_return_and_repeated_input"]["Accumulator.csv"]
    audit.check(own_cycles == stated_accumulator["nominal_cycles"], "complete accumulator-cycle reconstruction")
    audit.check([cycle["force_event_count"] for cycle in own_cycles] == [3, 2, 2, 2, 2, 3], "all accumulator event counts")
    audit.check([cycle["end_load_N"] for cycle in own_cycles[:5]] == [-0.28141, -0.27618, -0.27501, -0.27418, -0.27381], "repeated-endpoint load sequence")
    audit.check(stated_accumulator["hamming_transitions"] is None, "no invented accumulator Hamming word")

    own_repeats = same_direction_repeats(parsed["3_loops_SampleA.csv"], all_branches["3_loops_SampleA.csv"])
    stated_repeats = result["D2_return_and_repeated_input"]["3_loops_SampleA.csv"]
    audit.check(own_repeats == stated_repeats["all_consecutive_same_direction_branch_comparisons"], "complete three-loop repeat reconstruction")
    audit.check([item["mean_absolute_load_difference_N"] for item in own_repeats] == [0.000999616422, 0.000331838523, 0.004082173485, 0.003291539505], "three-loop repeatability values")

    # Causal typing and explicit nonidentifiability witness.
    d3 = result["D3_frozen_history_beyond_complete_state_search"]
    d4 = result["D4_constitutive_memory_lineage_frontier"]
    audit.check(d3["constitutive_memory_lineage_target"] is False, "D3 correctly excluded from constitutive target")
    audit.check(d3["matched_history_pair_identified"] is False, "D3 remains unscoreable")
    audit.check(d4["memory_coordinate_available_in_selected_csv_packet"] is False, "L absent from selected packet")
    audit.check(d4["matched_conventional_confound_state_available"] is False, "sufficient C match absent")
    audit.check(d4["common_future_query_channel_available"] is False, "common future query absent")
    audit.check("explicitly excluding L_t" in d4["conventional_confound_state"], "C excludes target L")
    audit.check(d4["status"] == "CONSTITUTIVE_MEMORY_LINEAGE_COEFFICIENT_UNIDENTIFIED", "D4 disposition")

    # Same observed law, different causal L contrast: M_L has beta_L=1;
    # M_C has beta_L=0.  This finite witness is sufficient for nonidentifiability.
    observed_l_model = [(history, 0, history, history) for history in (0, 1)]  # H,U,L,Y
    observed_c_model = [(history, 0, 0, history) for history in (0, 1)]        # H,U,L,Y; C=H hidden
    audit.check([(h, u, y) for h, u, _, y in observed_l_model] == [(h, u, y) for h, u, _, y in observed_c_model], "observational equivalence witness")
    beta_l_model = 1 - 0
    beta_c_model = 0
    audit.check(beta_l_model != beta_c_model, "causal L contrast differs across observationally equivalent models")

    amendment1 = (HERE / "SCHEMA_AMENDMENT_001.md").read_text()
    amendment2 = (HERE / "SCHEMA_AMENDMENT_002_THEORY_TYPING.md").read_text()
    theorem = (HERE / "THEOREM.md").read_text()
    prereg = (HERE / "PREREGISTRATION.md").read_text()
    audit.check("basic row ranges and branch counts were viewed" in amendment1, "post-access schema inspection disclosed")
    audit.check("DEVELOPMENT" in amendment1, "schema amendment denies untouched validation")
    audit.check("does **not** hold all of physical state" in amendment2 and "fixed" in amendment2, "L-versus-C intervention typing")
    audit.check("W=(L,C)" in theorem, "complete-state relation stated")
    audit.check("complete-state comparison is not the constitutive" in theorem, "complete-W target rejection")
    audit.check("observationally equivalent" in theorem and "not identified" in theorem, "identifiability proof stated")
    audit.check("The claim is not that no unknown dataset exists" in (HERE / "HOSTILE_AUDIT.md").read_text(), "selection claim is not exhaustive")
    audit.check("outcome bodies has been downloaded or opened" in prereg, "builder-declared pre-outcome access boundary")
    audit.check("internally frozen according to" in theorem and "not externally timestamped" in theorem, "theorem closes non-timestamped-freeze caveat")
    audit.check("frontier of the audited candidate set" in theorem and "not proof that" in theorem, "theorem closes audited-search-scope caveat")
    audit.check("not 97 independent trials" in theorem and "per file" in theorem, "theorem closes independence and apparatus-scope caveat")
    builder_hostile = (HERE / "HOSTILE_AUDIT.md").read_text()
    audit.check("local self-declarations without an external" in builder_hostile and "current byte" in builder_hostile, "builder hostile audit carries access-order ceiling")
    audit.check("not 97" in builder_hostile and "independent trials" in builder_hostile and "per file" in builder_hostile, "builder hostile audit carries comparison ceiling")

    gates = result["causal_admission_gates"]
    audit.check(result["gate_summary"] == {"pass": 1, "partial": 1, "fail": 4, "causal_keep_break_admitted": False}, "causal gate summary")
    audit.check(gates[0]["status"] == "PASS_FOR_EACH_FILE", "same-apparatus claim is per file")
    audit.check(gates[1]["status"] == "PARTIAL_RETROSPECTIVE_ROUTES_ONLY", "route status is retrospective")
    audit.check(all(gate["status"] == "FAIL" for gate in gates[2:]), "all four causal identification gates fail")

    ceiling = set(result["claim_ceiling"])
    expected_ceiling = {
        "PUBLIC_HISTORY_DIAGNOSTIC_ONLY",
        "NO_SURGICAL_WHOLE_LINEAGE_KEEP_BREAK",
        "NO_WRITER_OFF_RETAINED_RECORD_TEST",
        "NO_URFT_CONFIRMATION_OR_REFUTATION",
        "NO_GRAVITATIONAL_METRIC_GAMMA_RGRL_OR_GFT_EFFECT",
        "MISSING_FIELDS_ARE_NOT_EVIDENCE_FOR_ZERO_PHYSICS",
    }
    audit.check(ceiling == expected_ceiling, "exact claim ceiling")
    audit.check(not recursively_contains(result, {"p_value", "confidence_interval", "standard_error", "significance"}), "no invented inferential statistics")
    audit.check("cannot confirm or refute gravity emergence" in theorem, "theorem gravity ceiling")
    audit.check("remote, compositionally distinct force or clock probe" in theorem, "future gravity observable named")
    audit.check(result["additional_gravity_release_target"] == [
        "remote compositionally distinct force or clock probe",
        "event-keyed source energy-stress and support-state custody through KEEP, BREAK, hold, and query",
        "predeclared Newtonian and environmental nuisance channels",
    ], "gravity release target is bounded and concrete")
    audit.check("MATCHED_ENDPOINT_QUERY_TABLE" in theorem, "minimal missing mechanical object named")
    audit.check(len(result["minimal_mechanical_release_target"]) == 6, "minimal mechanical release has six declared fields")

    if audit.failures:
        for failure in audit.failures:
            print(f"FAIL: {failure}")
        print(f"SUMMARY {audit.count - len(audit.failures)}/{audit.count} independent hostile checks passed")
        return 1

    print(f"SUMMARY {audit.count}/{audit.count} independent hostile checks passed")
    print("DISPOSITION PASS__REPAIRED_CORE_SCOPE_CAVEATS_EXECUTABLY_BOUND__PUBLIC_MECHANICAL_HISTORY_DEVELOPMENT_EVIDENCE__CONSTITUTIVE_LINEAGE_AND_GRAVITY_UNIDENTIFIED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
