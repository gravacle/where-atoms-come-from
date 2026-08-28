#!/usr/bin/env python3
"""Independent verifier for the public mechanical-memory frontier lane.

This verifier does not import the production analyzer.  It reparses the raw
files, independently reconstructs branches, force events, comparison metrics,
return extrema, custody, and the claim ceiling.  It additionally asks the
production executable to replay its frozen RESULT.json.
"""

from __future__ import annotations

import hashlib
import json
import math
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "SOURCE"
NAMES = ["3_loops_SampleA.csv", "RPM.csv", "Accumulator.csv", "WeakRPMBreak.csv"]
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
EXPECTED_EVENTS = {
    "3_loops_SampleA.csv": 0,
    "RPM.csv": 14,
    "Accumulator.csv": 14,
    "WeakRPMBreak.csv": 6,
}


class Checks:
    def __init__(self) -> None:
        self.total = 0
        self.failures: list[str] = []

    def test(self, condition: bool, label: str) -> None:
        self.total += 1
        if not condition:
            self.failures.append(label)


def file_hash(path: Path, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def q(x: float) -> float:
    return round(float(x), 12)


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def parse_raw(name: str) -> dict[str, Any]:
    lines = (SRC / name).read_text(encoding="utf-8-sig").splitlines()
    header = lines[0].split(";")
    units = lines[1].split(";")
    numeric = [[float(cell.replace(",", ".")) for cell in line.split(";")] for line in lines[2:] if line]
    return {
        "header": header,
        "units": units,
        "t": [row[0] for row in numeric],
        "x": [row[1] for row in numeric],
        "f": [row[2] for row in numeric],
    }


def directions(x: list[float]) -> list[int]:
    raw = [(b > a) - (b < a) for a, b in zip(x, x[1:])]
    nonzero = [i for i, value in enumerate(raw) if value]
    for i, value in enumerate(raw):
        if not value:
            nearest = min(nonzero, key=lambda j: (abs(j - i), j))
            raw[i] = raw[nearest]
    return raw


def make_branches(data: dict[str, Any]) -> list[dict[str, Any]]:
    d = directions(data["x"])
    starts = [0] + [i for i in range(1, len(d)) if d[i] != d[i - 1]]
    stops = starts[1:] + [len(d)]
    branches = []
    for start, stop in zip(starts, stops):
        points = data["x"][start : stop + 1]
        if len(set(points)) < 3:
            continue
        branches.append(
            {
                "id": len(branches),
                "start": start,
                "end": stop,
                "dir": d[start],
                "xmin": min(points),
                "xmax": max(points),
            }
        )
    return branches


def collapse(data: dict[str, Any], branch: dict[str, Any]) -> tuple[list[float], list[float]]:
    grouped: dict[float, list[float]] = {}
    for x, force in zip(
        data["x"][branch["start"] : branch["end"] + 1],
        data["f"][branch["start"] : branch["end"] + 1],
    ):
        grouped.setdefault(x, []).append(force)
    xs = sorted(grouped)
    return xs, [sum(grouped[x]) / len(grouped[x]) for x in xs]


def interp(xs: list[float], ys: list[float], target: float) -> float:
    if target == xs[0]:
        return ys[0]
    if target == xs[-1]:
        return ys[-1]
    lo = 0
    hi = len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= target:
            lo = mid
        else:
            hi = mid
    if xs[lo] == target:
        return ys[lo]
    weight = (target - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] * (1 - weight) + ys[hi] * weight


def integrate(xs: list[float], ys: list[float]) -> float:
    answer = 0.0
    for i in range(len(xs) - 1):
        answer += (xs[i + 1] - xs[i]) * (ys[i + 1] + ys[i]) / 2
    return answer


def comparison(data: dict[str, Any], a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any] | None:
    ax, af = collapse(data, a)
    bx, bf = collapse(data, b)
    left = max(ax[0], bx[0])
    right = min(ax[-1], bx[-1])
    if right <= left:
        return None
    shared = sorted(x for x in set(ax).intersection(bx) if left <= x <= right)
    if len(shared) >= 3:
        grid = sorted(set(shared + [left, right]))
        mode = "exact_common_extension"
    else:
        grid = [left + i * (right - left) / 100 for i in range(101)]
        mode = "101_point_linear_interpolation"
    one = [interp(ax, af, x) for x in grid]
    two = [interp(bx, bf, x) for x in grid]
    delta = [u - v for u, v in zip(one, two)]
    width = right - left
    mean_abs = integrate(grid, [abs(v) for v in delta]) / width
    signed = integrate(grid, delta) / width
    rms = math.sqrt(integrate(grid, [v * v for v in delta]) / width)
    span = max(one + two) - min(one + two)
    point_object = [[q(x), q(v)] for x, v in zip(grid, delta)]
    return {
        "first_branch": a["id"],
        "second_branch": b["id"],
        "first_direction": "increasing" if a["dir"] > 0 else "decreasing",
        "second_direction": "increasing" if b["dir"] > 0 else "decreasing",
        "overlap_mm": [q(left), q(right)],
        "comparison_mode": mode,
        "comparison_point_count": len(grid),
        "signed_mean_load_difference_N": q(signed),
        "mean_absolute_load_difference_N": q(mean_abs),
        "rms_load_difference_N": q(rms),
        "maximum_absolute_load_difference_N": q(max(abs(v) for v in delta)),
        "pooled_load_range_N": q(span),
        "normalized_overlap_integrated_absolute_difference": q(mean_abs / span) if span else 0.0,
        "pointwise_extension_difference_sha256": hashlib.sha256(canonical_bytes(point_object)).hexdigest(),
    }


def opposite_comparisons(data: dict[str, Any], branches: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    for i, a in enumerate(branches):
        for b in branches[i + 1 :]:
            if a["dir"] == b["dir"]:
                continue
            item = comparison(data, a, b)
            if item is not None:
                output.append(item)
    return output


def events(data: dict[str, Any]) -> list[int]:
    candidates = [i for i, (a, b) in enumerate(zip(data["f"], data["f"][1:])) if abs(b - a) > 0.0079]
    if not candidates:
        return []
    return [candidates[i] for i in range(len(candidates) - 1) if candidates[i + 1] - candidates[i] > 15] + [candidates[-1]]


def extrema_and_returns(data: dict[str, Any], branches: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    extrema = []
    for a, b in zip(branches, branches[1:]):
        if a["end"] == b["start"] and a["dir"] != b["dir"]:
            i = a["end"]
            extrema.append(
                {
                    "ordinal": len(extrema),
                    "kind": "maximum" if a["dir"] > 0 else "minimum",
                    "index": i,
                    "time_sec": q(data["t"][i]),
                    "extension_mm": q(data["x"][i]),
                    "load_N": q(data["f"][i]),
                }
            )
    step = statistics.median(abs(b - a) for a, b in zip(data["x"], data["x"][1:]) if b != a)
    pairs = []
    for j, current in enumerate(extrema):
        for prior in extrema[:j]:
            dx = current["extension_mm"] - prior["extension_mm"]
            if current["kind"] == prior["kind"] and abs(dx) <= step:
                pairs.append(
                    {
                        "first_extremum": prior["ordinal"],
                        "returning_extremum": current["ordinal"],
                        "kind": current["kind"],
                        "extension_difference_mm": q(dx),
                        "load_difference_N": q(current["load_N"] - prior["load_N"]),
                        "complete_state_word_available": False,
                        "state_closed": None,
                    }
                )
    return extrema, pairs


def verify_manifests(checks: Checks) -> None:
    for manifest_name in ("SOURCE_MANIFEST.sha256", "CORE_MANIFEST.sha256"):
        path = ROOT / manifest_name
        checks.test(path.exists(), f"{manifest_name} exists")
        if not path.exists():
            continue
        for line in path.read_text().splitlines():
            digest, relative = line.split("  ", 1)
            target = ROOT / relative
            checks.test(target.is_file(), f"manifest target exists: {relative}")
            if target.is_file():
                checks.test(file_hash(target) == digest, f"manifest hash: {relative}")


def verify_seal(checks: Checks) -> None:
    path = ROOT / "LANE_SEAL.json"
    checks.test(path.exists(), "lane seal exists")
    if not path.exists():
        return
    seal = json.loads(path.read_text())
    checks.test(seal["protocol"] == "GRA-BENSE-PMHF-V001", "seal protocol")
    checks.test(seal["canonical_or_git_files_edited"] is False, "seal isolation declaration")
    expected = {
        "source_manifest_sha256": "SOURCE_MANIFEST.sha256",
        "core_manifest_sha256": "CORE_MANIFEST.sha256",
        "result_sha256": "RESULT.json",
        "theorem_sha256": "THEOREM.md",
        "hostile_audit_sha256": "HOSTILE_AUDIT.md",
        "analyzer_sha256": "analyze_public_history.py",
        "verifier_sha256": "verify_lane.py",
    }
    for field, relative in expected.items():
        checks.test(seal[field] == file_hash(ROOT / relative), f"seal hash {relative}")


def main() -> int:
    checks = Checks()
    result = json.loads((ROOT / "RESULT.json").read_text())
    custody = json.loads((ROOT / "SOURCE_CUSTODY.json").read_text())

    checks.test(result["protocol"] == "GRA-BENSE-PMHF-V001", "protocol")
    checks.test(result["execution_class"].startswith("RETROSPECTIVE_DEVELOPMENT"), "development class")
    checks.test(result["status"].endswith("CAUSAL_KEEP_BREAK_UNIDENTIFIED__NO_GRAVITY_TEST"), "bounded disposition")
    checks.test(result["source_custody_all_pass"] is True, "reported custody pass")

    expected = {item["name"]: item for item in custody["files"]}
    for name, item in expected.items():
        path = SRC / name
        checks.test(path.stat().st_size == item["bytes"], f"bytes {name}")
        checks.test(file_hash(path, "md5") == item["md5"], f"md5 {name}")
        checks.test(file_hash(path) == item["sha256"], f"sha256 {name}")
    paper = custody["primary_article"]
    paper_path = ROOT / paper["local_path"]
    checks.test(paper_path.stat().st_size == paper["bytes"], "paper bytes")
    checks.test(file_hash(paper_path) == paper["sha256"], "paper sha256")

    all_data: dict[str, dict[str, Any]] = {}
    all_branches: dict[str, list[dict[str, Any]]] = {}
    for name in NAMES:
        data = parse_raw(name)
        all_data[name] = data
        branches = make_branches(data)
        all_branches[name] = branches
        reported = result["per_file"][name]
        checks.test(data["header"] == ["Time", "Extension", "Load"], f"header {name}")
        checks.test(data["units"] == ["(sec)", "(mm)", "(N)"], f"units {name}")
        checks.test(len(data["t"]) == EXPECTED_ROWS[name], f"row count independent {name}")
        checks.test(all(b > a for a, b in zip(data["t"], data["t"][1:])), f"time order {name}")
        checks.test(reported["schema"]["row_count"] == len(data["t"]), f"row result {name}")
        checks.test(reported["schema"]["time_range_sec"] == [q(min(data["t"])), q(max(data["t"]))], f"time range {name}")
        checks.test(reported["schema"]["extension_range_mm"] == [q(min(data["x"])), q(max(data["x"]))], f"extension range {name}")
        checks.test(reported["schema"]["load_range_N"] == [q(min(data["f"])), q(max(data["f"]))], f"load range {name}")
        checks.test(len(branches) == EXPECTED_BRANCHES[name], f"branch count independent {name}")
        checks.test(len(reported["monotone_branches"]) == len(branches), f"branch count result {name}")
        for own, stated in zip(branches, reported["monotone_branches"]):
            checks.test(stated["id"] == own["id"], f"branch id {name}/{own['id']}")
            checks.test(stated["start_index"] == own["start"] and stated["end_index"] == own["end"], f"branch custody {name}/{own['id']}")
            checks.test(stated["direction"] == own["dir"], f"branch direction {name}/{own['id']}")
        independent_events = events(data)
        stated_events = reported["official_force_jump_proxy"]["events"]
        checks.test(len(independent_events) == EXPECTED_EVENTS[name], f"event count independent {name}")
        checks.test(len(stated_events) == len(independent_events), f"event count result {name}")
        checks.test([x["delta_start_index"] for x in stated_events] == independent_events, f"event indices {name}")
        for index, stated in zip(independent_events, stated_events):
            checks.test(stated["delta_load_N"] == q(data["f"][index + 1] - data["f"][index]), f"event delta {name}/{index}")
        independent_pairs = opposite_comparisons(data, branches)
        stated_pairs = reported["D1_all_opposite_branch_comparisons"]
        checks.test(len(independent_pairs) == len(stated_pairs), f"D1 pair count {name}")
        for position, (own, stated) in enumerate(zip(independent_pairs, stated_pairs)):
            checks.test(own == stated, f"D1 complete comparison {name}/{position}")

    for name in ("RPM.csv", "WeakRPMBreak.csv"):
        own_extrema, own_pairs = extrema_and_returns(all_data[name], all_branches[name])
        stated = result["D2_return_and_repeated_input"][name]
        checks.test(stated["local_extrema"] == own_extrema, f"extrema {name}")
        checks.test(stated["all_same_type_return_pairs"] == own_pairs, f"return pairs {name}")
        checks.test(stated["state_closure_status"].startswith("UNSCOREABLE"), f"state closure ceiling {name}")

    accumulator = result["D2_return_and_repeated_input"]["Accumulator.csv"]
    checks.test(len(accumulator["nominal_cycles"]) == 6, "accumulator all inferred cycles")
    checks.test(accumulator["complete_state_word_available"] is False, "accumulator no state word")
    checks.test(accumulator["hamming_transitions"] is None, "accumulator no invented Hamming")
    checks.test(
        [cycle["force_event_count"] for cycle in accumulator["nominal_cycles"]] == [3, 2, 2, 2, 2, 3],
        "accumulator event accounting",
    )

    d3 = result["D3_frozen_history_beyond_complete_state_search"]
    checks.test(d3["constitutive_memory_lineage_target"] is False, "D3 is not constitutive target")
    checks.test(d3["complete_state_word_available_in_selected_csv_packet"] is False, "D3 W absent in selected packet")
    checks.test(d3["internal_stress_field_available"] is False, "D3 stress absent")
    checks.test(d3["common_future_query_channel_available"] is False, "D3 query absent")
    checks.test(d3["matched_history_pair_identified"] is False, "D3 no pair")
    checks.test(d3["status"].endswith("AND_NOT_THE_PROGRAM_TARGET"), "D3 disposition")
    d4 = result["D4_constitutive_memory_lineage_frontier"]
    checks.test(d4["memory_coordinate_available_in_selected_csv_packet"] is False, "D4 L absent")
    checks.test(d4["matched_conventional_confound_state_available"] is False, "D4 C match absent")
    checks.test(d4["common_future_query_channel_available"] is False, "D4 query absent")
    checks.test(d4["status"] == "CONSTITUTIVE_MEMORY_LINEAGE_COEFFICIENT_UNIDENTIFIED", "D4 disposition")
    checks.test(result["gate_summary"] == {"pass": 1, "partial": 1, "fail": 4, "causal_keep_break_admitted": False}, "gate summary")
    checks.test(result["freeze_custody_scope"].startswith("INTERNALLY_FROZEN"), "internal freeze scope")
    checks.test("NOT_EXTERNALLY_TIMESTAMPED" in result["freeze_custody_scope"], "no external timestamp promotion")
    checks.test(result["search_scope"].startswith("FRONTIER_OF_DOCUMENTED_AUDITED_CANDIDATE_SET"), "audited candidate-set scope")
    checks.test("NOT_97_INDEPENDENT_TRIALS" in result["comparison_independence"], "comparison independence ceiling")
    checks.test("PER_FILE" in result["comparison_independence"], "same-apparatus per-file scope")

    ceiling = set(result["claim_ceiling"])
    for token in (
        "NO_SURGICAL_WHOLE_LINEAGE_KEEP_BREAK",
        "NO_WRITER_OFF_RETAINED_RECORD_TEST",
        "NO_URFT_CONFIRMATION_OR_REFUTATION",
        "NO_GRAVITATIONAL_METRIC_GAMMA_RGRL_OR_GFT_EFFECT",
        "MISSING_FIELDS_ARE_NOT_EVIDENCE_FOR_ZERO_PHYSICS",
    ):
        checks.test(token in ceiling, f"claim ceiling {token}")

    theorem = (ROOT / "THEOREM.md").read_text()
    for phrase in (
        "causal memory-lineage or gravity-history contrast",
        "observationally equivalent",
        "not identified",
        "does **not** say that history or memory is unreal",
        "complete-state comparison is not the constitutive",
        "MATCHED_ENDPOINT_QUERY_TABLE",
        "remote, compositionally distinct force or clock probe",
        "cannot confirm or refute gravity emergence",
        "internally frozen according to",
        "project custody",
        "not externally timestamped",
        "frontier of the audited candidate set",
        "not 97 independent trials",
    ):
        checks.test(phrase in theorem, f"theorem ceiling phrase: {phrase}")
    checks.test("0.0079" in (SRC / "Fig5ac.py").read_text(), "official event threshold script 1")
    checks.test("dj>15" in (SRC / "Fig5ac.py").read_text(), "official event gap script 1")
    checks.test("0.0079" in (SRC / "Fig6.py").read_text(), "official event threshold script 2")
    checks.test("dj>15" in (SRC / "Fig6.py").read_text(), "official event gap script 2")
    amendment = (ROOT / "SCHEMA_AMENDMENT_002_THEORY_TYPING.md").read_text()
    checks.test("constitutive memory-lineage estimand" in amendment, "theory-typing amendment present")
    checks.test("does **not** demand that all future-active physical state" in theorem, "constitutive target does not freeze all state")

    replay = subprocess.run(
        [sys.executable, str(ROOT / "analyze_public_history.py"), "--check"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    checks.test(replay.returncode == 0, "production deterministic replay")
    checks.test("PASS" in replay.stdout, "production replay emits PASS")

    verify_manifests(checks)
    verify_seal(checks)

    if checks.failures:
        for failure in checks.failures:
            print(f"FAIL: {failure}")
        print(f"FAIL: {checks.total - len(checks.failures)}/{checks.total} checks passed")
        return 1
    print(f"PASS: {checks.total}/{checks.total} independent checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
