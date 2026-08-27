#!/usr/bin/env python3
"""Retrospective dimensionless DCD response analysis on five physical specimens."""

from __future__ import annotations

from bisect import bisect_right
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MODEL = ROOT / "model"
if str(MODEL) not in sys.path:
    sys.path.insert(0, str(MODEL))

from lakeshore_vsm import PAIR_TOKEN, dcd_writer_off, read_trace


PROTOCOL = HERE / "EXPLORATORY_NORMALIZED_RESPONSE_PROTOCOL_V002.md"
RAW = HERE / "raw"
GRID = tuple(round(index * 0.05, 2) for index in range(1, 31))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interpolate(xs: list[float], ys: list[float], point: float) -> float:
    if point < xs[0] or point > xs[-1]:
        raise ValueError(f"point {point} is outside [{xs[0]}, {xs[-1]}]")
    right = bisect_right(xs, point)
    if right == 0:
        return ys[0]
    if right == len(xs):
        return ys[-1]
    left = right - 1
    if xs[left] == point:
        return ys[left]
    fraction = (point - xs[left]) / (xs[right] - xs[left])
    return ys[left] + fraction * (ys[right] - ys[left])


def main() -> int:
    paths = sorted(RAW.glob("* Step 2 Remanence Curves.csv"))
    if len(paths) != 5:
        raise ValueError(f"expected five DCD curves, found {len(paths)}")

    curves: dict[str, dict] = {}
    for path in paths:
        sample = path.name.split(PAIR_TOKEN, 1)[0]
        trace = read_trace(path, sample_id=sample)
        strict_crossings = [
            index for index in range(len(trace.moment_Am2) - 1)
            if trace.moment_Am2[index] * trace.moment_Am2[index + 1] < 0
        ]
        if len(strict_crossings) != 1:
            raise ValueError(
                f"expected one strict retained-moment sign change for {sample}, "
                f"found {len(strict_crossings)}"
            )
        result = dcd_writer_off(trace)
        initial = trace.moment_Am2[0]
        zero_field = result["remanent_sign_change_field_T"]
        if initial == 0 or zero_field == 0:
            raise ValueError(f"zero normalization for {sample}")
        xs = [abs(field / zero_field) for field in trace.field_T]
        rs = [moment / initial for moment in trace.moment_Am2]
        if any(b < a for a, b in zip(xs, xs[1:])):
            raise ValueError(f"normalized writer field is not monotone for {sample}")
        values = [interpolate(xs, rs, point) for point in GRID]
        curves[sample] = {
            "source_path": str(path.relative_to(ROOT)),
            "source_sha256": sha256(path),
            "row_count": len(xs),
            "initial_remanent_moment_Am2": initial,
            "sign_change_field_T": zero_field,
            "grid_values": values,
        }

    grid_rows = []
    residuals: list[float] = []
    for index, point in enumerate(GRID):
        values = [curves[sample]["grid_values"][index] for sample in sorted(curves)]
        mean = statistics.fmean(values)
        spread = statistics.pstdev(values)
        residuals.extend(value - mean for value in values)
        grid_rows.append({
            "normalized_writer_field": point,
            "values": values,
            "mean": mean,
            "population_std": spread,
            "range": max(values) - min(values),
        })

    payload = {
        "analysis": "RETROSPECTIVE_NORMALIZED_DCD_RESPONSE_V001",
        "protocol_sha256": sha256(PROTOCOL),
        "sample_count": len(curves),
        "samples": sorted(curves),
        "curves": curves,
        "grid": grid_rows,
        "summary": {
            "rms_deviation_from_gridwise_mean": math.sqrt(
                statistics.fmean(value * value for value in residuals)
            ),
            "largest_population_std": max(row["population_std"] for row in grid_rows),
            "largest_range": max(row["range"] for row in grid_rows),
        },
        "scientific_scope": "ONE_MECHANISM_RETROSPECTIVE_DESCRIPTION",
        "formation_verdict": "NONE_NOT_SCORED",
        "universal_claim_authorized": False,
        "gravity_claim_authorized": False,
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    (HERE / "NORMALIZED_RESPONSE.json").write_text(rendered, encoding="utf-8")

    lines = [
        "EXPLORATORY NORMALIZED WRITER-OFF RESPONSE",
        "=" * 58,
        f"samples: {payload['sample_count']}",
        f"grid: {GRID[0]:.2f}..{GRID[-1]:.2f} in steps of 0.05",
        ("RMS deviation from gridwise mean: "
         f"{payload['summary']['rms_deviation_from_gridwise_mean']:.12g}"),
        f"largest population std: {payload['summary']['largest_population_std']:.12g}",
        f"largest range: {payload['summary']['largest_range']:.12g}",
        "",
        "sample                              H_zero [T]       m_initial [A m^2] rows",
    ]
    for sample in sorted(curves):
        curve = curves[sample]
        lines.append(
            f"{sample:<35} {curve['sign_change_field_T']: .12g}  "
            f"{curve['initial_remanent_moment_Am2']: .12g} {curve['row_count']:4d}"
        )
    lines.extend([
        "",
        "VERDICT: NONE_NOT_SCORED — retrospective one-mechanism physics only.",
    ])
    (HERE / "NORMALIZED_RESPONSE.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
