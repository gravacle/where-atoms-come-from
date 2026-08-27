#!/usr/bin/env python3
"""Measure writer-off remanence in the five raw U1537 VSM sample pairs.

This is a retrospective measurement reader.  It intentionally emits no theory verdict and
does not score T-53.  The only fixed diagnostic scale is the declared 1 mT low-field probe;
its retention fraction is printed descriptively, not tested against a learned tolerance.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
MODEL = HERE.parent / "model"
if str(MODEL) not in sys.path:
    sys.path.insert(0, str(MODEL))

from lakeshore_vsm import dcd_writer_off, hysteresis_writer_off, read_trace


CONTROL_FIELD_T = 1.0e-3
PAIR_TOKEN = "-hys-dcd-forc1 "


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, default=HERE / "raw")
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()

    source = json.loads((args.raw / "SOURCE.json").read_text(encoding="utf-8"))
    frozen = {row["file"]: row["sha256"] for row in source["selected_files"]}
    files = sorted(path for path in args.raw.glob("*.csv"))
    if {path.name for path in files} != set(frozen):
        raise SystemExit("REFUSED: raw CSV membership differs from SOURCE.json")
    for path in files:
        if sha256(path) != frozen[path.name]:
            raise SystemExit(f"REFUSED: raw SHA-256 mismatch for {path.name}")

    grouped: dict[str, dict[str, Path]] = {}
    for path in files:
        if PAIR_TOKEN not in path.name:
            raise SystemExit(f"REFUSED: unrecognized source filename {path.name}")
        sample, kind = path.name.split(PAIR_TOKEN, 1)
        grouped.setdefault(sample, {})[kind] = path

    sample_rows = []
    for sample in sorted(grouped):
        pair = grouped[sample]
        if set(pair) != {
            "Step 1 Hysteresis Measurement.csv",
            "Step 2 Remanence Curves.csv",
        }:
            raise SystemExit(f"REFUSED: incomplete measurement pair for {sample}")
        hysteresis = read_trace(pair["Step 1 Hysteresis Measurement.csv"], sample_id=sample)
        remanence = read_trace(pair["Step 2 Remanence Curves.csv"], sample_id=sample)
        if hysteresis.instrument != remanence.instrument:
            raise SystemExit(f"REFUSED: cross-instrument pair for {sample}")
        h = hysteresis_writer_off(hysteresis)
        d = dcd_writer_off(remanence, control_field_T=CONTROL_FIELD_T)
        sample_rows.append(
            {
                "sample_id": sample,
                "instrument": hysteresis.instrument,
                "hysteresis": h,
                "reverse_pulse_curve": d,
            }
        )

    sign_observations = [
        sign
        for row in sample_rows
        for sign in row["hysteresis"]["writer_off_signs"]
    ]
    final_reversals = sum(
        row["reverse_pulse_curve"]["baseline_remanent_moment"]
        * row["reverse_pulse_curve"]["final_writer_off_moment"]
        < 0
        for row in sample_rows
    )
    result = {
        "schema": "T53A_WORLD_WRITE_RESULT_V001",
        "scope": "RETROSPECTIVE_REAL_SURFACE_PROCESS_OBSERVATION",
        "scientific_disposition": "EVIDENCE_ONLY_NO_THEORY_VERDICT",
        "t53_status_change_authorized": False,
        "dataset": {
            "doi": source["doi"],
            "related_publication_doi": source["related_publication_doi"],
            "license": source["license"],
            "archive_sha256": source["archive_sha256"],
            "source_manifest_sha256": sha256(args.raw / "SOURCE.json"),
        },
        "declared_diagnostic": {
            "low_field_control_T": CONTROL_FIELD_T,
            "semantics": "descriptive retention after reverse pulses no larger than 1 mT",
        },
        "aggregate_measurements": {
            "physical_samples": len(sample_rows),
            "writer_off_branch_signs": sign_observations,
            "opposite_writer_sign_pairs": sum(
                row["hysteresis"]["writer_off_signs"] == [1, -1]
                for row in sample_rows
            ),
            "maximum_branch_magnitude_asymmetry": max(
                row["hysteresis"]["magnitude_asymmetry"] for row in sample_rows
            ),
            "minimum_low_field_retention_fraction": min(
                row["reverse_pulse_curve"]["minimum_low_field_retention_fraction"]
                for row in sample_rows
            ),
            "remanent_sign_change_field_T_range": [
                min(abs(row["reverse_pulse_curve"]["remanent_sign_change_field_T"])
                    for row in sample_rows),
                max(abs(row["reverse_pulse_curve"]["remanent_sign_change_field_T"])
                    for row in sample_rows),
            ],
            "final_reverse_pulse_sign_reversals": final_reversals,
        },
        "samples": sample_rows,
        "limitations": [
            "outcomes were public before this exploratory analysis",
            "not randomized or blinded",
            "no independent no-write cohort",
            "no long-time retention series",
            "no observer-fragment redundancy measurement",
            "one magnetic mechanism only",
            "cannot establish clause necessity, a general formation law, or universality",
        ],
    }

    args.output.mkdir(parents=True, exist_ok=True)
    json_path = args.output / "RESULT.json"
    text_path = args.output / "RESULT.txt"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    a = result["aggregate_measurements"]
    lines = [
        "T-53A — RAW WORLD-SURFACE WRITE OBSERVATION",
        "=" * 72,
        f"scope: {result['scope']}",
        f"disposition: {result['scientific_disposition']}",
        f"dataset: {source['doi']} ({source['license']})",
        f"physical samples: {a['physical_samples']}",
        f"opposite writer-off sign pairs: {a['opposite_writer_sign_pairs']}/{a['physical_samples']}",
        f"writer-off branch signs: {a['writer_off_branch_signs']}",
        "maximum |+Mr| vs |-Mr| asymmetry: "
        f"{a['maximum_branch_magnitude_asymmetry']:.9f}",
        f"declared low-field control: {CONTROL_FIELD_T:.6g} T",
        "minimum remanence fraction after <= control field: "
        f"{a['minimum_low_field_retention_fraction']:.9f}",
        "remanent sign-change field magnitude range: "
        f"[{a['remanent_sign_change_field_T_range'][0]:.9g}, "
        f"{a['remanent_sign_change_field_T_range'][1]:.9g}] T",
        f"final reverse-pulse sign reversals: {a['final_reverse_pulse_sign_reversals']}/"
        f"{a['physical_samples']}",
        "",
        "READ: actual VSM data exhibit writer-direction-dependent remanence after the",
        "writer field is removed, retention under the declared small reverse probe, and",
        "a finite sign-change field on all five samples. This is one real-surface",
        "observation of candidate formation terms. It is retrospective, single-mechanism,",
        "and supplies no general-theory or T-53 verdict.",
    ]
    text_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
