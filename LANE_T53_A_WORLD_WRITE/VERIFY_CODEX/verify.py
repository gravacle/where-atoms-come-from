#!/usr/bin/env python3
"""Independent, stdlib-only verifier for the T53-A raw VSM measurements.

This verifier intentionally does not import any project model or lane analyzer.
Its verdict is REFUTED unless every predeclared integrity, protocol, numerical,
and scope predicate executes and succeeds.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import traceback
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
SOURCE = RAW / "SOURCE.json"

STEP1_SUFFIX = "-hys-dcd-forc1 Step 1 Hysteresis Measurement.csv"
STEP2_SUFFIX = "-hys-dcd-forc1 Step 2 Remanence Curves.csv"

REQUIRED_PREDICATES = (
    "source_manifest_parseable",
    "source_identity_fixed",
    "manifest_has_ten_unique_measurements",
    "raw_directory_membership_exact",
    "manifest_sizes_all_match",
    "manifest_md5_all_match",
    "manifest_sha256_all_match",
    "five_physical_sample_pairs_exact",
    "hysteresis_headers_define_applied_field_protocol",
    "hysteresis_rows_all_finite_and_good",
    "hysteresis_paths_traverse_both_directions",
    "hysteresis_zero_field_interpolants_unique",
    "hysteresis_writer_off_remanence_opposes_by_history",
    "dcd_headers_define_writer_off_read_protocol",
    "dcd_rows_all_finite_and_good",
    "dcd_reverse_pulse_paths_complete",
    "dcd_low_field_retention_computable",
    "dcd_sign_change_threshold_unique",
    "dcd_independent_threshold_agrees_with_embedded_extraction",
    "dcd_final_reverse_pulse_reverses_sign",
    "scientific_scope_explicitly_bounded",
)


def digest_bytes(data: bytes, algorithm: str) -> str:
    return hashlib.new(algorithm, data).hexdigest()


def digest_file(path: Path, algorithm: str = "sha256") -> str:
    return digest_bytes(path.read_bytes(), algorithm)


def read_measurement(path: Path) -> dict[str, Any]:
    """Parse one Lake Shore export without relying on project machinery."""

    raw = path.read_bytes()
    text = raw.decode("latin-1")
    lines = text.splitlines()
    table_index = next(i for i, line in enumerate(lines) if line.startswith("Step,"))
    header_lines = lines[:table_index]
    rows: list[dict[str, Any]] = []
    for cells in csv.reader(lines[table_index + 1 :]):
        if not cells or not cells[0].strip():
            continue
        if len(cells) < 8:
            raise ValueError(f"short data row in {path.name}: {cells!r}")
        rows.append(
            {
                "step": int(cells[0]),
                "iteration": int(cells[1]),
                "segment": int(cells[2]),
                "field_T": float(cells[3]),
                "moment_Am2": float(cells[4]),
                "time_s": float(cells[5]),
                "field_status": cells[6].strip(),
                "moment_status": cells[7].strip(),
            }
        )

    key_values: dict[str, str] = {}
    for line in header_lines:
        if ": " in line:
            key, value = line.split(": ", 1)
            key_values[key.strip()] = value.strip()

    return {
        "name": path.name,
        "raw_sha256": digest_bytes(raw, "sha256"),
        "header_lines": header_lines,
        "key_values": key_values,
        "rows": rows,
    }


def has_line(measurement: dict[str, Any], exact_line: str) -> bool:
    return exact_line in measurement["header_lines"]


def interpolate_y_at_x_zero(points: list[tuple[float, float]]) -> tuple[float, tuple[float, float]]:
    candidates: list[tuple[float, tuple[float, float]]] = []
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 == 0.0:
            candidates.append((y0, (x0, x0)))
        elif x1 == 0.0:
            candidates.append((y1, (x1, x1)))
        elif x0 * x1 < 0.0:
            y_at_zero = y0 + (0.0 - x0) * (y1 - y0) / (x1 - x0)
            candidates.append((y_at_zero, (x0, x1)))
    if len(candidates) != 1:
        raise ValueError(f"expected one x=0 bracket, found {len(candidates)}")
    return candidates[0]


def interpolate_x_at_y_zero(points: list[tuple[float, float]]) -> tuple[float, tuple[float, float]]:
    candidates: list[tuple[float, tuple[float, float]]] = []
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if y0 == 0.0:
            candidates.append((x0, (y0, y0)))
        elif y1 == 0.0:
            candidates.append((x1, (y1, y1)))
        elif y0 * y1 < 0.0:
            x_at_zero = x0 + (0.0 - y0) * (x1 - x0) / (y1 - y0)
            candidates.append((x_at_zero, (y0, y1)))
    if len(candidates) != 1:
        raise ValueError(f"expected one y=0 bracket, found {len(candidates)}")
    return candidates[0]


def render_text(result: dict[str, Any]) -> str:
    lines = [
        "INDEPENDENT VSM RAW-MEASUREMENT VERIFICATION",
        f"Verdict: {result['verdict']}",
        f"Scope: {result['verdict_scope']}",
        "",
        "Default rule: REFUTED unless every declared predicate executes and succeeds.",
        f"Predicates executed: {result['predicate_summary']['executed']}/{result['predicate_summary']['required']}",
        f"Predicates succeeded: {result['predicate_summary']['succeeded']}/{result['predicate_summary']['required']}",
        "",
        "Independent numerical results:",
    ]
    for sample in result.get("samples", []):
        h = sample["hysteresis"]
        d = sample["dcd"]
        lines.extend(
            [
                f"- {sample['sample_id']}",
                f"  hysteresis H=0 remanence after +field history: {h['descending_zero_field_remanence_Am2']:.16e} A m^2",
                f"  hysteresis H=0 remanence after -field history: {h['ascending_zero_field_remanence_Am2']:.16e} A m^2",
                f"  branch-magnitude asymmetry: {h['branch_magnitude_asymmetry_fraction']:.16e}",
                f"  DCD retained fraction at largest reverse pulse <=1 mT: {d['low_field_retained_fraction']:.16e}",
                f"  DCD sign-change pulse field: {d['computed_sign_change_field_T']:.16e} T",
                f"  DCD final reverse pulse: {d['final_reverse_pulse_field_T']:.16e} T; retained moment {d['final_writer_off_moment_Am2']:.16e} A m^2",
            ]
        )
    lines.extend(["", "What these files support:"])
    lines.extend(f"- {item}" for item in result["adjudication"]["supports"])
    lines.extend(["", "What these files do not support:"])
    lines.extend(f"- {item}" for item in result["adjudication"]["does_not_support"])
    lines.extend(["", "Predicate ledger:"])
    for predicate in result["predicates"]:
        state = "SUCCEEDED" if predicate["succeeded"] else "FAILED"
        lines.append(f"- {predicate['name']}: {state}")
    if result["errors"]:
        lines.extend(["", "Verifier errors:"])
        lines.extend(f"- {error['type']}: {error['message']}" for error in result["errors"])
    return "\n".join(lines) + "\n"


def render_d24(result: dict[str, Any]) -> str:
    failures = [p for p in result["predicates"] if not p["succeeded"]]
    lines = [
        "# D24 audit — independent T53-A VSM verification",
        "",
        f"- Verdict: `{result['verdict']}`",
        f"- Verdict scope: `{result['verdict_scope']}`",
        "- Default: `REFUTED` until every declared predicate executes and succeeds.",
        f"- Predicates executed: {result['predicate_summary']['executed']} of {result['predicate_summary']['required']}",
        f"- Predicates succeeded: {result['predicate_summary']['succeeded']} of {result['predicate_summary']['required']}",
        "",
        "## Errors and predicate failures",
        "",
    ]
    if not result["errors"] and not failures:
        lines.append("None observed in this bounded raw-integrity and measurement recomputation.")
    else:
        for error in result["errors"]:
            lines.append(f"- `{error['type']}`: {error['message']}")
        for failure in failures:
            lines.append(f"- Predicate `{failure['name']}` failed: `{json.dumps(failure['evidence'], sort_keys=True)}`")
    lines.extend(
        [
            "",
            "## Scope audit",
            "",
            "The non-refutation verdict, if reached, applies only to exact raw-file integrity, protocol semantics exposed by the raw headers, and the independently recomputed numerical measurements. It is not a record-formation law, a universal result, a gravity result, or a program-status decision.",
            "",
            "The dataset is retrospective and restricted to five labeled magnetic sediment specimens. It contains no randomized intended-message assignment, no blind prediction, no no-write control cohort, no long common-hold survival study, no independent laboratory reproduction, no cross-surface coverage, and no gravity observable.",
            "",
            "## Reproducibility",
            "",
            "Run `PYTHONDONTWRITEBYTECODE=1 python3 -B verify.py` from this directory. The verifier uses only the Python standard library and writes deterministic result artifacts with no clock or machine-specific path fields.",
        ]
    )
    return "\n".join(lines) + "\n"


def verify() -> dict[str, Any]:
    result: dict[str, Any] = {
        "schema": "wacf.independent-vsm-raw-verification.v1",
        "verdict": "REFUTED",
        "verdict_scope": "raw-integrity-and-measurement-only",
        "default_rule": "REFUTED unless every declared predicate executes and succeeds",
        "verifier_sha256": digest_file(Path(__file__)),
        "source_manifest_sha256": digest_file(SOURCE),
        "predicates": [],
        "samples": [],
        "errors": [],
        "adjudication": {
            "supports": [
                "The supplied files are byte-identical to the ten measurements pinned by SOURCE.json.",
                "Five distinctly labeled physical sediment specimens underwent externally applied magnetic-field histories.",
                "Both hysteresis branches retain nonzero, history-dependent moments at interpolated zero applied field.",
                "The DCD protocol explicitly records moments away from the applied reverse-pulse field, and retained response remains measurable after low-field pulses.",
                "Each DCD series has an independently recomputable sign-change threshold and an opposite-sign retained response after the final reverse pulse.",
            ],
            "does_not_support": [
                "A general or universal law of record formation, or universal formation terms.",
                "Necessity or sufficiency of any proposed record-formation criterion.",
                "Prospective or blinded prediction, randomized message writing, a no-write control, or long common-hold survival.",
                "Coverage of any surface class beyond these five magnetic sediment specimens, much less any and every bona-fide record surface.",
                "Redundant-carrier formation, independent-physicist reproduction, or a released-URM reproduction.",
                "Gravity emergence or any gravity observable.",
                "Completion of T53, a ledger status change, registration, or a program-level proof claim.",
            ],
        },
    }

    def predicate(name: str, succeeded: bool, evidence: Any) -> None:
        if name not in REQUIRED_PREDICATES:
            raise KeyError(f"undeclared predicate: {name}")
        if any(item["name"] == name for item in result["predicates"]):
            raise KeyError(f"duplicate predicate: {name}")
        result["predicates"].append(
            {"name": name, "executed": True, "succeeded": bool(succeeded), "evidence": evidence}
        )

    try:
        source = json.loads(SOURCE.read_text(encoding="utf-8"))
        predicate("source_manifest_parseable", isinstance(source, dict), {"type": type(source).__name__})

        expected_identity = {
            "record_id": "14564186",
            "doi": "10.5281/zenodo.14564186",
            "license": "cc-by-4.0",
            "related_publication_doi": "10.1029/2025PA005360",
            "archive_sha256": "202a932bcf30af20087c21277c9d19f2b7217a69d849a198db9016f170f330ab",
        }
        observed_identity = {key: source.get(key) for key in expected_identity}
        predicate("source_identity_fixed", observed_identity == expected_identity, observed_identity)

        selected = source.get("selected_files", [])
        names = [item.get("file") for item in selected if isinstance(item, dict)]
        predicate(
            "manifest_has_ten_unique_measurements",
            len(selected) == 10 and len(names) == 10 and len(set(names)) == 10 and all(isinstance(n, str) for n in names),
            {"entries": len(selected), "unique_names": len(set(names))},
        )

        observed_raw_names = sorted(path.name for path in RAW.iterdir())
        expected_raw_names = sorted(["SOURCE.json", *names])
        predicate(
            "raw_directory_membership_exact",
            observed_raw_names == expected_raw_names and all((RAW / name).is_file() and not (RAW / name).is_symlink() for name in expected_raw_names),
            {"observed": observed_raw_names, "expected": expected_raw_names},
        )

        size_results: dict[str, dict[str, Any]] = {}
        md5_results: dict[str, dict[str, Any]] = {}
        sha_results: dict[str, dict[str, Any]] = {}
        for entry in selected:
            path = RAW / entry["file"]
            data = path.read_bytes()
            size_results[path.name] = {"expected": entry["bytes"], "observed": len(data)}
            md5_results[path.name] = {
                "expected": entry["md5_repository"],
                "observed": digest_bytes(data, "md5"),
            }
            sha_results[path.name] = {
                "expected": entry["sha256"],
                "observed": digest_bytes(data, "sha256"),
            }
        predicate(
            "manifest_sizes_all_match",
            all(v["expected"] == v["observed"] for v in size_results.values()),
            size_results,
        )
        predicate(
            "manifest_md5_all_match",
            all(v["expected"] == v["observed"] for v in md5_results.values()),
            md5_results,
        )
        predicate(
            "manifest_sha256_all_match",
            all(v["expected"] == v["observed"] for v in sha_results.values()),
            sha_results,
        )

        step1_names = sorted(name for name in names if name.endswith(STEP1_SUFFIX))
        step2_names = sorted(name for name in names if name.endswith(STEP2_SUFFIX))
        step1_samples = {name[: -len(STEP1_SUFFIX)] for name in step1_names}
        step2_samples = {name[: -len(STEP2_SUFFIX)] for name in step2_names}
        sample_ids = sorted(step1_samples)
        predicate(
            "five_physical_sample_pairs_exact",
            len(step1_names) == 5
            and len(step2_names) == 5
            and step1_samples == step2_samples
            and len(sample_ids) == 5
            and all(re.fullmatch(r"U1537A-[A-Za-z0-9-]+cm", sample) for sample in sample_ids),
            {"step1_samples": sorted(step1_samples), "step2_samples": sorted(step2_samples)},
        )

        hys_by_sample = {
            sample: read_measurement(RAW / f"{sample}{STEP1_SUFFIX}") for sample in sample_ids
        }
        dcd_by_sample = {
            sample: read_measurement(RAW / f"{sample}{STEP2_SUFFIX}") for sample in sample_ids
        }

        hys_header_evidence: dict[str, Any] = {}
        hys_headers_ok = True
        for sample, measurement in hys_by_sample.items():
            evidence = {
                "measurement_marker": has_line(measurement, "#HYSTERESIS MEASUREMENT"),
                "max_field": measurement["key_values"].get("Max field"),
                "initial_curve": measurement["key_values"].get("Include initial curve"),
                "acquisition_mode": measurement["key_values"].get("Acquisition mode"),
                "table": next(
                    (line for line in measurement["header_lines"] if line.startswith("##DATA TABLE")), None
                ),
            }
            ok = (
                evidence["measurement_marker"]
                and evidence["max_field"] == "1.5 T"
                and evidence["initial_curve"] == "False"
                and evidence["acquisition_mode"] == "Continuous"
                and evidence["table"] is not None
                and "Moment" in evidence["table"]
                and "Field" in evidence["table"]
            )
            hys_headers_ok = hys_headers_ok and ok
            hys_header_evidence[sample] = evidence
        predicate("hysteresis_headers_define_applied_field_protocol", hys_headers_ok, hys_header_evidence)

        hys_rows_evidence: dict[str, Any] = {}
        hys_rows_ok = True
        for sample, measurement in hys_by_sample.items():
            rows = measurement["rows"]
            finite = all(
                math.isfinite(row[key])
                for row in rows
                for key in ("field_T", "moment_Am2", "time_s")
            )
            statuses_good = all(
                row["field_status"] == "GOOD" and row["moment_status"] == "GOOD" for row in rows
            )
            structure = all(row["step"] == 1 and row["iteration"] == 0 for row in rows)
            timestamps = all(b["time_s"] >= a["time_s"] for a, b in zip(rows, rows[1:]))
            ok = len(rows) > 1000 and finite and statuses_good and structure and timestamps
            hys_rows_ok = hys_rows_ok and ok
            hys_rows_evidence[sample] = {
                "rows": len(rows),
                "finite": finite,
                "statuses_good": statuses_good,
                "step_iteration_expected": structure,
                "timestamps_nondecreasing": timestamps,
            }
        predicate("hysteresis_rows_all_finite_and_good", hys_rows_ok, hys_rows_evidence)

        hys_path_evidence: dict[str, Any] = {}
        hys_paths_ok = True
        for sample, measurement in hys_by_sample.items():
            segments: dict[int, list[dict[str, Any]]] = {}
            for row in measurement["rows"]:
                segments.setdefault(row["segment"], []).append(row)
            segment_ids = sorted(segments)
            descending = segments.get(0, [])
            ascending = segments.get(1, [])
            endpoints = {
                "segment_ids": segment_ids,
                "descending_start_T": descending[0]["field_T"] if descending else None,
                "descending_end_T": descending[-1]["field_T"] if descending else None,
                "ascending_start_T": ascending[0]["field_T"] if ascending else None,
                "ascending_end_T": ascending[-1]["field_T"] if ascending else None,
            }
            ok = (
                segment_ids == [0, 1]
                and endpoints["descending_start_T"] > 1.49
                and endpoints["descending_end_T"] < -1.49
                and endpoints["ascending_start_T"] < -1.49
                and endpoints["ascending_end_T"] > 1.49
            )
            hys_paths_ok = hys_paths_ok and ok
            hys_path_evidence[sample] = endpoints
        predicate("hysteresis_paths_traverse_both_directions", hys_paths_ok, hys_path_evidence)

        hys_numerical: dict[str, Any] = {}
        hys_interpolants_ok = True
        hys_opposition_ok = True
        for sample, measurement in hys_by_sample.items():
            segments: dict[int, list[tuple[float, float]]] = {0: [], 1: []}
            for row in measurement["rows"]:
                segments[row["segment"]].append((row["field_T"], row["moment_Am2"]))
            descending_m, descending_bracket = interpolate_y_at_x_zero(segments[0])
            ascending_m, ascending_bracket = interpolate_y_at_x_zero(segments[1])
            denominator = (abs(descending_m) + abs(ascending_m)) / 2.0
            asymmetry = abs(abs(descending_m) - abs(ascending_m)) / denominator
            hys_numerical[sample] = {
                "descending_zero_field_remanence_Am2": descending_m,
                "descending_zero_field_bracket_T": list(descending_bracket),
                "ascending_zero_field_remanence_Am2": ascending_m,
                "ascending_zero_field_bracket_T": list(ascending_bracket),
                "branch_magnitude_asymmetry_fraction": asymmetry,
            }
            hys_interpolants_ok = hys_interpolants_ok and all(
                math.isfinite(value) for value in (descending_m, ascending_m, asymmetry)
            )
            hys_opposition_ok = hys_opposition_ok and descending_m > 0.0 and ascending_m < 0.0
        predicate("hysteresis_zero_field_interpolants_unique", hys_interpolants_ok, hys_numerical)
        predicate(
            "hysteresis_writer_off_remanence_opposes_by_history",
            hys_opposition_ok,
            {
                sample: {
                    "after_positive_field_history_Am2": values["descending_zero_field_remanence_Am2"],
                    "after_negative_field_history_Am2": values["ascending_zero_field_remanence_Am2"],
                }
                for sample, values in hys_numerical.items()
            },
        )

        dcd_header_evidence: dict[str, Any] = {}
        dcd_headers_ok = True
        for sample, measurement in dcd_by_sample.items():
            kv = measurement["key_values"]
            evidence = {
                "measurement_marker": has_line(measurement, "#REMANENCE CURVES MEASUREMENT"),
                "irm": kv.get("IRM"),
                "dcd": kv.get("DCD"),
                "saturation_field": kv.get("Saturation field"),
                "points": kv.get("Number of points"),
                "final_field": kv.get("Final field"),
                "measure_at_applied_fields": kv.get("Measure moment at applied fields"),
                "pause_at_zero_field": kv.get("Pause at zero field"),
            }
            ok = (
                evidence["measurement_marker"]
                and evidence["irm"] == "False"
                and evidence["dcd"] == "True"
                and evidence["saturation_field"] == "2.5 T"
                and evidence["points"] == "150"
                and evidence["final_field"] == "2.5 T"
                and evidence["measure_at_applied_fields"] == "False"
                and evidence["pause_at_zero_field"] == "0 s"
            )
            dcd_headers_ok = dcd_headers_ok and ok
            dcd_header_evidence[sample] = evidence
        predicate("dcd_headers_define_writer_off_read_protocol", dcd_headers_ok, dcd_header_evidence)

        dcd_rows_evidence: dict[str, Any] = {}
        dcd_rows_ok = True
        for sample, measurement in dcd_by_sample.items():
            rows = measurement["rows"]
            finite = all(
                math.isfinite(row[key])
                for row in rows
                for key in ("field_T", "moment_Am2", "time_s")
            )
            statuses_good = all(
                row["field_status"] == "GOOD" and row["moment_status"] == "GOOD" for row in rows
            )
            structure = all(
                row["step"] == 2 and row["iteration"] == 0 and row["segment"] == 0 for row in rows
            )
            timestamps = all(b["time_s"] >= a["time_s"] for a, b in zip(rows, rows[1:]))
            ok = len(rows) == 150 and finite and statuses_good and structure and timestamps
            dcd_rows_ok = dcd_rows_ok and ok
            dcd_rows_evidence[sample] = {
                "rows": len(rows),
                "finite": finite,
                "statuses_good": statuses_good,
                "step_iteration_segment_expected": structure,
                "timestamps_nondecreasing": timestamps,
            }
        predicate("dcd_rows_all_finite_and_good", dcd_rows_ok, dcd_rows_evidence)

        dcd_path_evidence: dict[str, Any] = {}
        dcd_paths_ok = True
        for sample, measurement in dcd_by_sample.items():
            rows = measurement["rows"]
            magnitudes = [abs(row["field_T"]) for row in rows]
            nondecreasing = all(b >= a for a, b in zip(magnitudes, magnitudes[1:]))
            first_is_nearest_zero = magnitudes[0] == min(magnitudes)
            final_field = rows[-1]["field_T"]
            ok = nondecreasing and first_is_nearest_zero and final_field < 0 and abs(abs(final_field) - 2.5) < 1e-4
            dcd_paths_ok = dcd_paths_ok and ok
            dcd_path_evidence[sample] = {
                "first_pulse_field_T": rows[0]["field_T"],
                "final_pulse_field_T": final_field,
                "absolute_pulse_magnitude_nondecreasing": nondecreasing,
                "first_point_nearest_zero": first_is_nearest_zero,
            }
        predicate("dcd_reverse_pulse_paths_complete", dcd_paths_ok, dcd_path_evidence)

        dcd_numerical: dict[str, Any] = {}
        low_retention_ok = True
        threshold_ok = True
        embedded_agreement_ok = True
        final_reversal_ok = True
        for sample, measurement in dcd_by_sample.items():
            rows = measurement["rows"]
            baseline = rows[0]["moment_Am2"]
            low_rows = [row for row in rows if abs(row["field_T"]) <= 0.001]
            low_endpoint = max(low_rows, key=lambda row: abs(row["field_T"]))
            retained_fraction = low_endpoint["moment_Am2"] / baseline
            points = [(row["field_T"], row["moment_Am2"]) for row in rows]
            sign_change_field, moment_bracket = interpolate_x_at_y_zero(points)
            embedded_hr = float(measurement["key_values"]["Hr [T]"])
            threshold_delta = sign_change_field - embedded_hr
            final = rows[-1]
            dcd_numerical[sample] = {
                "initial_near_zero_pulse_field_T": rows[0]["field_T"],
                "initial_writer_off_moment_Am2": baseline,
                "low_field_reverse_pulse_T": low_endpoint["field_T"],
                "low_field_writer_off_moment_Am2": low_endpoint["moment_Am2"],
                "low_field_retained_fraction": retained_fraction,
                "computed_sign_change_field_T": sign_change_field,
                "sign_change_moment_bracket_Am2": list(moment_bracket),
                "instrument_embedded_Hr_T": embedded_hr,
                "computed_minus_embedded_Hr_T": threshold_delta,
                "final_reverse_pulse_field_T": final["field_T"],
                "final_writer_off_moment_Am2": final["moment_Am2"],
            }
            low_retention_ok = (
                low_retention_ok
                and len(low_rows) >= 2
                and baseline > 0.0
                and low_endpoint["moment_Am2"] > 0.0
                and math.isfinite(retained_fraction)
                and retained_fraction > 0.0
            )
            threshold_ok = threshold_ok and sign_change_field < 0.0 and abs(sign_change_field) < 2.5
            embedded_agreement_ok = embedded_agreement_ok and abs(threshold_delta) <= 5e-15
            final_reversal_ok = final_reversal_ok and baseline * final["moment_Am2"] < 0.0
        predicate("dcd_low_field_retention_computable", low_retention_ok, dcd_numerical)
        predicate(
            "dcd_sign_change_threshold_unique",
            threshold_ok,
            {sample: values["computed_sign_change_field_T"] for sample, values in dcd_numerical.items()},
        )
        predicate(
            "dcd_independent_threshold_agrees_with_embedded_extraction",
            embedded_agreement_ok,
            {
                sample: values["computed_minus_embedded_Hr_T"] for sample, values in dcd_numerical.items()
            },
        )
        predicate(
            "dcd_final_reverse_pulse_reverses_sign",
            final_reversal_ok,
            {
                sample: {
                    "initial_Am2": values["initial_writer_off_moment_Am2"],
                    "final_Am2": values["final_writer_off_moment_Am2"],
                }
                for sample, values in dcd_numerical.items()
            },
        )

        scope = result["adjudication"]
        scope_ok = (
            len(scope["supports"]) >= 5
            and len(scope["does_not_support"]) >= 7
            and any("universal" in item for item in scope["does_not_support"])
            and any("Gravity" in item for item in scope["does_not_support"])
            and any("T53" in item for item in scope["does_not_support"])
        )
        predicate("scientific_scope_explicitly_bounded", scope_ok, scope)

        result["samples"] = [
            {
                "sample_id": sample,
                "files": {
                    "hysteresis": hys_by_sample[sample]["name"],
                    "dcd": dcd_by_sample[sample]["name"],
                },
                "hysteresis": hys_numerical[sample],
                "dcd": dcd_numerical[sample],
            }
            for sample in sample_ids
        ]
    except Exception as exc:  # Preserve a deterministic audit rather than losing the failed run.
        result["errors"].append(
            {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            }
        )

    executed_names = [item["name"] for item in result["predicates"]]
    executed_exactly_once = len(executed_names) == len(set(executed_names))
    all_executed = set(executed_names) == set(REQUIRED_PREDICATES) and len(executed_names) == len(REQUIRED_PREDICATES)
    succeeded = sum(1 for item in result["predicates"] if item["succeeded"])
    result["predicate_summary"] = {
        "required": len(REQUIRED_PREDICATES),
        "executed": len(executed_names),
        "succeeded": succeeded,
        "all_executed_exactly_once": executed_exactly_once and all_executed,
    }
    if all_executed and executed_exactly_once and succeeded == len(REQUIRED_PREDICATES) and not result["errors"]:
        result["verdict"] = "NOT_REFUTED"
    return result


def write_outputs(result: dict[str, Any]) -> None:
    json_path = HERE / "VERDICT.json"
    text_path = HERE / "VERDICT.txt"
    audit_path = HERE / "D24_AUDIT.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    text_path.write_text(render_text(result), encoding="utf-8")
    audit_path.write_text(render_d24(result), encoding="utf-8")

    output_hashes = {
        "schema": "wacf.independent-vsm-output-hashes.v1",
        "verifier": {"file": "verify.py", "sha256": digest_file(Path(__file__))},
        "outputs": {
            path.name: {"bytes": path.stat().st_size, "sha256": digest_file(path)}
            for path in (audit_path, json_path, text_path)
        },
    }
    (HERE / "OUTPUT_HASHES.json").write_text(
        json.dumps(output_hashes, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    verification = verify()
    write_outputs(verification)
    print(f"{verification['verdict']} {verification['predicate_summary']['succeeded']}/{verification['predicate_summary']['required']}")

