#!/usr/bin/env python3
"""Internally frozen public mechanical-memory diagnostic for GRA-BENSE-PMHF-V001.

The executable deliberately stops short of a causal KEEP/BREAK claim.  It
checks source custody, executes the internally frozen reversible-history
diagnostics on every eligible branch, and records why the deposited fields do
not identify a lineage-only coefficient or a gravitational response. The
local freeze is not an externally timestamped registration.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "SOURCE"
CSV_NAMES = (
    "3_loops_SampleA.csv",
    "RPM.csv",
    "Accumulator.csv",
    "WeakRPMBreak.csv",
)
EXPECTED_HEADERS = ["Time", "Extension", "Load"]
EXPECTED_UNITS = ["(sec)", "(mm)", "(N)"]
EVENT_THRESHOLD_N = 0.0079
EVENT_CLUSTER_GAP_SAMPLES = 15
GRID_SIZE = 101


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def digest_file(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def q(value: float) -> float:
    """Stable report precision; far below the deposited five-decimal resolution."""
    return round(float(value), 12)


def verify_custody() -> list[dict[str, Any]]:
    custody = json.loads((ROOT / "SOURCE_CUSTODY.json").read_text())
    rows: list[dict[str, Any]] = []
    for item in custody["files"]:
        path = SOURCE / item["name"]
        actual = {
            "name": item["name"],
            "bytes": path.stat().st_size,
            "md5": digest_file(path, "md5"),
            "sha256": digest_file(path, "sha256"),
        }
        actual["pass"] = (
            actual["bytes"] == item["bytes"]
            and actual["md5"] == item["md5"]
            and actual["sha256"] == item["sha256"]
        )
        rows.append(actual)
    paper = custody["primary_article"]
    paper_path = ROOT / paper["local_path"]
    paper_actual = {
        "name": paper_path.name,
        "bytes": paper_path.stat().st_size,
        "sha256": digest_file(paper_path, "sha256"),
    }
    paper_actual["pass"] = (
        paper_actual["bytes"] == paper["bytes"]
        and paper_actual["sha256"] == paper["sha256"]
    )
    rows.append(paper_actual)
    if not all(row["pass"] for row in rows):
        raise RuntimeError("source-custody mismatch")
    return rows


def parse_decimal_comma(token: str) -> float:
    return float(token.replace(",", "."))


def read_series(name: str) -> dict[str, Any]:
    with (SOURCE / name).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader)
        units = next(reader)
        if header != EXPECTED_HEADERS or units != EXPECTED_UNITS:
            raise RuntimeError(f"unexpected schema in {name}: {header!r} / {units!r}")
        rows = [tuple(parse_decimal_comma(cell) for cell in row) for row in reader if row]
    if not rows:
        raise RuntimeError(f"empty series: {name}")
    time = [row[0] for row in rows]
    extension = [row[1] for row in rows]
    load = [row[2] for row in rows]
    if any(b <= a for a, b in zip(time, time[1:])):
        raise RuntimeError(f"time is not strictly increasing: {name}")
    positive_steps = [abs(b - a) for a, b in zip(extension, extension[1:]) if b != a]
    if not positive_steps:
        raise RuntimeError(f"drive coordinate is constant: {name}")
    return {
        "name": name,
        "time": time,
        "extension": extension,
        "load": load,
        "median_positive_extension_step_mm": statistics.median(positive_steps),
    }


def filled_interval_directions(extension: list[float]) -> list[int]:
    raw = [1 if b > a else -1 if b < a else 0 for a, b in zip(extension, extension[1:])]
    if not any(raw):
        raise RuntimeError("no nonzero drive increment")
    left: list[int | None] = [None] * len(raw)
    last: int | None = None
    for i, direction in enumerate(raw):
        if direction:
            last = i
        left[i] = last
    right: list[int | None] = [None] * len(raw)
    nxt: int | None = None
    for i in range(len(raw) - 1, -1, -1):
        if raw[i]:
            nxt = i
        right[i] = nxt
    filled = raw[:]
    for i, direction in enumerate(filled):
        if direction:
            continue
        li, ri = left[i], right[i]
        if li is None:
            chosen = ri
        elif ri is None:
            chosen = li
        else:
            chosen = li if i - li <= ri - i else ri
        assert chosen is not None
        filled[i] = raw[chosen]
    return filled


def monotone_branches(series: dict[str, Any]) -> list[dict[str, Any]]:
    directions = filled_interval_directions(series["extension"])
    groups: list[tuple[int, int, int]] = []
    start_interval = 0
    for i in range(1, len(directions)):
        if directions[i] != directions[i - 1]:
            groups.append((start_interval, i, directions[i - 1]))
            start_interval = i
    groups.append((start_interval, len(directions), directions[-1]))
    branches: list[dict[str, Any]] = []
    for interval_start, interval_stop, direction in groups:
        point_start, point_end = interval_start, interval_stop
        values = series["extension"][point_start : point_end + 1]
        if len(set(values)) < 3:
            continue
        branches.append(
            {
                "id": len(branches),
                "start_index": point_start,
                "end_index": point_end,
                "direction": direction,
                "direction_label": "increasing" if direction > 0 else "decreasing",
                "extension_min_mm": min(values),
                "extension_max_mm": max(values),
                "distinct_extension_count": len(set(values)),
            }
        )
    return branches


def collapsed_branch(series: dict[str, Any], branch: dict[str, Any]) -> tuple[list[float], list[float]]:
    accum: dict[float, list[float]] = {}
    lo, hi = branch["start_index"], branch["end_index"]
    for x, force in zip(series["extension"][lo : hi + 1], series["load"][lo : hi + 1]):
        accum.setdefault(x, []).append(force)
    xs = sorted(accum)
    forces = [statistics.fmean(accum[x]) for x in xs]
    return xs, forces


def interpolate(xs: list[float], ys: list[float], x: float) -> float:
    if x < xs[0] or x > xs[-1]:
        raise ValueError("interpolation outside branch")
    if x == xs[0]:
        return ys[0]
    if x == xs[-1]:
        return ys[-1]
    low, high = 0, len(xs) - 1
    while high - low > 1:
        mid = (low + high) // 2
        if xs[mid] <= x:
            low = mid
        else:
            high = mid
    if xs[low] == x:
        return ys[low]
    fraction = (x - xs[low]) / (xs[high] - xs[low])
    return ys[low] + fraction * (ys[high] - ys[low])


def trapezoid(xs: list[float], ys: list[float]) -> float:
    return sum((b - a) * (u + v) / 2 for a, b, u, v in zip(xs, xs[1:], ys, ys[1:]))


def compare_branches(
    series: dict[str, Any], first: dict[str, Any], second: dict[str, Any]
) -> dict[str, Any] | None:
    x1, y1 = collapsed_branch(series, first)
    x2, y2 = collapsed_branch(series, second)
    low, high = max(x1[0], x2[0]), min(x1[-1], x2[-1])
    if high <= low:
        return None
    exact = sorted(x for x in set(x1) & set(x2) if low <= x <= high)
    if len(exact) >= 3:
        # The integral is over the *closed* overlap.  If one support endpoint
        # was not deposited by both branches, interpolate only that endpoint;
        # every interior value remains an exact shared Extension value.
        grid = sorted(set(exact + [low, high]))
        mode = "exact_common_extension"
    else:
        grid = [low + (high - low) * i / (GRID_SIZE - 1) for i in range(GRID_SIZE)]
        mode = "101_point_linear_interpolation"
    first_load = [interpolate(x1, y1, x) for x in grid]
    second_load = [interpolate(x2, y2, x) for x in grid]
    differences = [a - b for a, b in zip(first_load, second_load)]
    abs_difference = [abs(value) for value in differences]
    width = high - low
    pooled_range = max(first_load + second_load) - min(first_load + second_load)
    mean_abs = trapezoid(grid, abs_difference) / width
    signed_mean = trapezoid(grid, differences) / width
    rms = math.sqrt(trapezoid(grid, [value * value for value in differences]) / width)
    if pooled_range == 0:
        normalized = 0.0 if mean_abs == 0 else None
    else:
        normalized = mean_abs / pooled_range
    point_digest_input = [[q(x), q(value)] for x, value in zip(grid, differences)]
    return {
        "first_branch": first["id"],
        "second_branch": second["id"],
        "first_direction": first["direction_label"],
        "second_direction": second["direction_label"],
        "overlap_mm": [q(low), q(high)],
        "comparison_mode": mode,
        "comparison_point_count": len(grid),
        "signed_mean_load_difference_N": q(signed_mean),
        "mean_absolute_load_difference_N": q(mean_abs),
        "rms_load_difference_N": q(rms),
        "maximum_absolute_load_difference_N": q(max(abs_difference)),
        "pooled_load_range_N": q(pooled_range),
        "normalized_overlap_integrated_absolute_difference": None if normalized is None else q(normalized),
        "pointwise_extension_difference_sha256": hashlib.sha256(canonical_bytes(point_digest_input)).hexdigest(),
    }


def all_opposite_branch_comparisons(
    series: dict[str, Any], branches: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    comparisons: list[dict[str, Any]] = []
    for i, first in enumerate(branches):
        for second in branches[i + 1 :]:
            if first["direction"] == second["direction"]:
                continue
            comparison = compare_branches(series, first, second)
            if comparison is not None:
                comparisons.append(comparison)
    return comparisons


def official_force_events(series: dict[str, Any], branches: list[dict[str, Any]]) -> list[dict[str, Any]]:
    load = series["load"]
    candidates = [i for i, (a, b) in enumerate(zip(load, load[1:])) if abs(b - a) > EVENT_THRESHOLD_N]
    if not candidates:
        return []
    retained = [candidates[i] for i in range(len(candidates) - 1) if candidates[i + 1] - candidates[i] > EVENT_CLUSTER_GAP_SAMPLES]
    retained.append(candidates[-1])
    events: list[dict[str, Any]] = []
    for index in retained:
        branch_id = None
        for branch in branches:
            if branch["start_index"] <= index < branch["end_index"]:
                branch_id = branch["id"]
                break
        events.append(
            {
                "delta_start_index": index,
                "time_sec": q(series["time"][index]),
                "extension_mm": q(series["extension"][index]),
                "load_before_N": q(load[index]),
                "load_after_N": q(load[index + 1]),
                "delta_load_N": q(load[index + 1] - load[index]),
                "branch_id": branch_id,
            }
        )
    return events


def extrema(branches: list[dict[str, Any]], series: dict[str, Any]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for previous, following in zip(branches, branches[1:]):
        if previous["end_index"] != following["start_index"]:
            continue
        if previous["direction"] == following["direction"]:
            continue
        index = previous["end_index"]
        kind = "maximum" if previous["direction"] > 0 else "minimum"
        out.append(
            {
                "ordinal": len(out),
                "kind": kind,
                "index": index,
                "time_sec": q(series["time"][index]),
                "extension_mm": q(series["extension"][index]),
                "load_N": q(series["load"][index]),
            }
        )
    return out


def return_revisits(series: dict[str, Any], branches: list[dict[str, Any]]) -> dict[str, Any]:
    local_extrema = extrema(branches, series)
    tolerance = series["median_positive_extension_step_mm"]
    pairs: list[dict[str, Any]] = []
    for j, current in enumerate(local_extrema):
        for prior in local_extrema[:j]:
            delta_x = current["extension_mm"] - prior["extension_mm"]
            if current["kind"] == prior["kind"] and abs(delta_x) <= tolerance:
                pairs.append(
                    {
                        "first_extremum": prior["ordinal"],
                        "returning_extremum": current["ordinal"],
                        "kind": current["kind"],
                        "extension_difference_mm": q(delta_x),
                        "load_difference_N": q(current["load_N"] - prior["load_N"]),
                        "complete_state_word_available": False,
                        "state_closed": None,
                    }
                )
    return {
        "matching_tolerance_mm": q(tolerance),
        "local_extrema": local_extrema,
        "all_same_type_return_pairs": pairs,
        "state_closure_status": "UNSCOREABLE_COMPLETE_STATE_WORD_NOT_DEPOSITED",
    }


def event_count(events: list[dict[str, Any]], start: int, end: int) -> int:
    return sum(start <= event["delta_start_index"] < end for event in events)


def inferred_cycles(
    series: dict[str, Any], branches: list[dict[str, Any]], events: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    cycles: list[dict[str, Any]] = []
    for offset in range(0, len(branches) - 1, 2):
        first, second = branches[offset], branches[offset + 1]
        if first["direction"] == second["direction"] or first["end_index"] != second["start_index"]:
            continue
        start, turn, end = first["start_index"], first["end_index"], second["end_index"]
        cycles.append(
            {
                "cycle": len(cycles),
                "boundary_basis": "drive_extrema_inferred_from_deposited_Extension",
                "start_index": start,
                "turn_index": turn,
                "end_index": end,
                "start_extension_mm": q(series["extension"][start]),
                "turn_extension_mm": q(series["extension"][turn]),
                "end_extension_mm": q(series["extension"][end]),
                "start_load_N": q(series["load"][start]),
                "end_load_N": q(series["load"][end]),
                "force_event_count": event_count(events, start, end),
                "complete_state_word_available": False,
                "state_transition_hamming": None,
            }
        )
    return cycles


def repeatability_comparisons(
    series: dict[str, Any], branches: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    by_direction = {
        direction: [branch for branch in branches if branch["direction"] == direction]
        for direction in (-1, 1)
    }
    for direction in (-1, 1):
        family = by_direction[direction]
        for first, second in zip(family, family[1:]):
            comparison = compare_branches(series, first, second)
            if comparison is not None:
                output.append(comparison)
    return output


def series_result(series: dict[str, Any]) -> dict[str, Any]:
    branches = monotone_branches(series)
    events = official_force_events(series, branches)
    return {
        "schema": {
            "headers": EXPECTED_HEADERS,
            "units": EXPECTED_UNITS,
            "row_count": len(series["time"]),
            "time_range_sec": [q(min(series["time"])), q(max(series["time"]))],
            "extension_range_mm": [q(min(series["extension"])), q(max(series["extension"]))],
            "load_range_N": [q(min(series["load"])), q(max(series["load"]))],
            "time_strictly_increasing": True,
            "median_positive_extension_step_mm": q(series["median_positive_extension_step_mm"]),
        },
        "monotone_branches": [
            {key: q(value) if isinstance(value, float) else value for key, value in branch.items()}
            for branch in branches
        ],
        "D1_all_opposite_branch_comparisons": all_opposite_branch_comparisons(series, branches),
        "official_force_jump_proxy": {
            "threshold_N_strictly_greater_than": EVENT_THRESHOLD_N,
            "cluster_separation_samples_strictly_greater_than": EVENT_CLUSTER_GAP_SAMPLES,
            "event_count": len(events),
            "events": events,
            "interpretation": "candidate_force_events_not_a_complete_hysteron_state_word",
        },
        "_branches_internal": branches,
        "_events_internal": events,
    }


def build_result() -> dict[str, Any]:
    custody = verify_custody()
    parsed = {name: read_series(name) for name in CSV_NAMES}
    per_file = {name: series_result(series) for name, series in parsed.items()}

    rpm_branches = per_file["RPM.csv"].pop("_branches_internal")
    per_file["RPM.csv"].pop("_events_internal")
    weak_branches = per_file["WeakRPMBreak.csv"].pop("_branches_internal")
    per_file["WeakRPMBreak.csv"].pop("_events_internal")
    accumulator_branches = per_file["Accumulator.csv"].pop("_branches_internal")
    accumulator_events = per_file["Accumulator.csv"].pop("_events_internal")
    repeat_branches = per_file["3_loops_SampleA.csv"].pop("_branches_internal")
    per_file["3_loops_SampleA.csv"].pop("_events_internal")

    d2 = {
        "RPM.csv": return_revisits(parsed["RPM.csv"], rpm_branches),
        "WeakRPMBreak.csv": return_revisits(parsed["WeakRPMBreak.csv"], weak_branches),
        "Accumulator.csv": {
            "nominal_cycles": inferred_cycles(
                parsed["Accumulator.csv"], accumulator_branches, accumulator_events
            ),
            "cycle_boundary_custody": "DRIVE_EXTREMA_INFERRED_SCRIPT_STATE_BINDING_UNAVAILABLE",
            "complete_state_word_available": False,
            "hamming_transitions": None,
        },
        "3_loops_SampleA.csv": {
            "cycle_count_from_branch_pairs": len(repeat_branches) // 2,
            "all_consecutive_same_direction_branch_comparisons": repeatability_comparisons(
                parsed["3_loops_SampleA.csv"], repeat_branches
            ),
            "role": "apparatus_repeatability_development_check_not_independent_validation",
        },
    }

    gates = [
        {
            "gate": 1,
            "requirement": "one physical apparatus and chronological source schedule",
            "status": "PASS_FOR_EACH_FILE",
            "basis": "publication-owned Time, Extension, and Load acquisition",
        },
        {
            "gate": 2,
            "requirement": "prospectively identifiable alternative routes",
            "status": "PARTIAL_RETROSPECTIVE_ROUTES_ONLY",
            "basis": "reversible histories are deposited but this execution is DEVELOPMENT",
        },
        {
            "gate": 3,
            "requirement": "authenticated memory-lineage coordinate changed by KEEP/BREAK while conventional confound state is matched separately",
            "status": "FAIL",
            "basis": "the selected CSV packet has no event-bound groove-state/transition-lineage coordinate and no sufficient conventional-state ledger",
        },
        {
            "gate": 4,
            "requirement": "common later query schedule",
            "status": "FAIL",
            "basis": "no event-keyed paired route endpoints followed by one bound query channel",
        },
        {
            "gate": 5,
            "requirement": "matched work, heat, stress, support, actuator, geometry, and read disturbance",
            "status": "FAIL",
            "basis": "only Time, Extension, and Load are deposited",
        },
        {
            "gate": 6,
            "requirement": "independent trials or joint covariance/noise model",
            "status": "FAIL",
            "basis": "no trial IDs or joint covariance/noise object are deposited",
        },
    ]
    return {
        "protocol": "GRA-BENSE-PMHF-V001",
        "execution_class": "RETROSPECTIVE_DEVELOPMENT_AFTER_DISCLOSED_SCHEMA_INSPECTION",
        "status": "PUBLIC_SAME_APPARATUS_REVERSIBLE_HISTORY_DIAGNOSTIC_EXECUTED__CAUSAL_KEEP_BREAK_UNIDENTIFIED__NO_GRAVITY_TEST",
        "freeze_custody_scope": "INTERNALLY_FROZEN_ACCORDING_TO_PROJECT_CUSTODY__NOT_EXTERNALLY_TIMESTAMPED_OR_INDEPENDENTLY_REGISTERED",
        "search_scope": "FRONTIER_OF_DOCUMENTED_AUDITED_CANDIDATE_SET__NOT_AN_EXHAUSTIVE_WORLD_PUBLIC_DATASET_SEARCH",
        "comparison_independence": "97_WITHIN_FILE_COMPARISONS_REUSE_FOUR_CHRONOLOGICAL_TRACES__NOT_97_INDEPENDENT_TRIALS__SAME_APPARATUS_STATUS_IS_PER_FILE",
        "source_custody_all_pass": all(row["pass"] for row in custody),
        "source_custody": custody,
        "per_file": per_file,
        "D2_return_and_repeated_input": d2,
        "D3_frozen_history_beyond_complete_state_search": {
            "estimand_class": "DIRECT_RESIDUAL_HISTORY_EFFECT_CONDITIONAL_ON_COMPLETE_PRESENT_STATE",
            "constitutive_memory_lineage_target": False,
            "complete_state_word_available_in_selected_csv_packet": False,
            "internal_stress_field_available": False,
            "common_future_query_channel_available": False,
            "matched_history_pair_identified": False,
            "status": "DIRECT_RESIDUAL_HISTORY_COEFFICIENT_UNIDENTIFIED_AND_NOT_THE_PROGRAM_TARGET",
            "reason": "the frozen D3 endpoint and query fields are absent; moreover, if complete future-active physical state were truly identical, state sufficiency makes a same-channel direct history residual trivial rather than a constitutive-lineage test",
        },
        "D4_constitutive_memory_lineage_frontier": {
            "target_memory_coordinate": "L_t = authenticated event-keyed groove-state/transition-lineage coordinate deliberately changed by KEEP versus BREAK",
            "conventional_confound_state": "C_t = Extension, baseline Load, internal stress-strain, geometry/damage, energy/work/heat/acoustic loss, actuator/support/read/EM state, explicitly excluding L_t",
            "memory_coordinate_available_in_selected_csv_packet": False,
            "matched_conventional_confound_state_available": False,
            "common_future_query_channel_available": False,
            "status": "CONSTITUTIVE_MEMORY_LINEAGE_COEFFICIENT_UNIDENTIFIED",
            "reason": "the observed path dependence can be assigned either to an authenticated memory coordinate or to unmeasured conventional present state with the same deposited Time/Extension/Load law",
        },
        "causal_admission_gates": gates,
        "gate_summary": {
            "pass": 1,
            "partial": 1,
            "fail": 4,
            "causal_keep_break_admitted": False,
        },
        "strongest_supported_result": "within each deposited same-apparatus file, the reversible loading trace contains branch-dependent mechanical Load and publication-rule force-event trajectories; the four files are not one joined trial",
        "claim_ceiling": [
            "PUBLIC_HISTORY_DIAGNOSTIC_ONLY",
            "NO_SURGICAL_WHOLE_LINEAGE_KEEP_BREAK",
            "NO_WRITER_OFF_RETAINED_RECORD_TEST",
            "NO_URFT_CONFIRMATION_OR_REFUTATION",
            "NO_GRAVITATIONAL_METRIC_GAMMA_RGRL_OR_GFT_EFFECT",
            "MISSING_FIELDS_ARE_NOT_EVIDENCE_FOR_ZERO_PHYSICS",
        ],
        "minimal_mechanical_release_target": [
            "event-keyed randomized paired KEEP and surgical BREAK routes in the same specimen that deliberately differ in authenticated lineage L_t",
            "independent event-bound measurement of L_t, such as the declared per-groove state/transition ancestry vector",
            "matched conventional confound state C_t excluding L_t: Extension, baseline Load, internal stress-strain, geometry/damage, support, actuator state, energy/work/heat, acoustic loss, read disturbance, and EM conditions",
            "matched declared marginal occupation if the lineage intervention changes assignment, pairing, correlation, or ancestry at fixed count",
            "writer off plus registered settling interval followed by one identical subthreshold query schedule",
            "independent block IDs and a joint covariance/noise object",
        ],
        "additional_gravity_release_target": [
            "remote compositionally distinct force or clock probe",
            "event-keyed source energy-stress and support-state custody through KEEP, BREAK, hold, and query",
            "predeclared Newtonian and environmental nuisance channels",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, help="write canonical JSON result")
    parser.add_argument("--check", action="store_true", help="compare recomputation to RESULT.json")
    args = parser.parse_args()
    result = build_result()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        expected_path = ROOT / "RESULT.json"
        if not expected_path.exists():
            raise SystemExit("RESULT.json missing")
        if expected_path.read_text() != rendered:
            raise SystemExit("RESULT.json does not match deterministic recomputation")
        print("PASS: RESULT.json matches deterministic recomputation")
        return 0
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
