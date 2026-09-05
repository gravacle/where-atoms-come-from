#!/usr/bin/env python3
"""Independent result-level audit for the closed GL6FJ L=4 m001 packet.

This is deliberately a reconstruction audit, not a response promotion.  It
does not import the author's reconstructor, invert a response, perform a
Schur complement, or name a Ward null.  It authenticates the closed result,
parses all 300 raw scalar projections, reconstructs each of the six physical
response owners by polarization, and compares that independently reconstructed
finite response to the promoted packet.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import stat
import sys

import numpy as np


sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.resolve()
RESULT = ROOT / "RESULT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
TARGET = ROOT / "DEVELOPMENT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
PRELAUNCH_AUDIT = ROOT / "AUDIT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
CHANNELS = 24
UPPER = CHANNELS * (CHANNELS + 1) // 2
OWNERS = (
    "stay_contact", "writer_contact", "stay_spectral", "writer_spectral",
    "stay_writer_interference", "disconnected",
)
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def unique(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key " + key)
        result[key] = value
    return result


def strict_json(path: Path):
    value = json.loads(
        path.read_bytes().decode("utf-8"), object_pairs_hook=unique,
        parse_constant=lambda token: (_ for _ in ()).throw(
            ValueError("nonfinite JSON " + token)))

    def visit(node):
        if type(node) is dict:
            for child in node.values():
                visit(child)
        elif type(node) is list:
            for child in node:
                visit(child)
        elif type(node) is int:
            check(-(2 ** 63) <= node <= 2 ** 64 - 1, "JSON integer range")
        elif type(node) is float:
            check(math.isfinite(node), "JSON finite number")
        else:
            check(type(node) in (str, bool, type(None)), "JSON scalar type")
    visit(value)
    return value


def manifest_rows(path: Path):
    rows = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        check(line.count("  ") == 1, "manifest delimiter " + str(number))
        expected, name = line.split("  ", 1)
        pure = PurePosixPath(name)
        check(len(expected) == 64 and set(expected) <= set("0123456789abcdef"),
              "manifest digest " + str(number))
        check(name not in rows and name and not pure.is_absolute() and
              all(part not in {"", ".", ".."} for part in pure.parts),
              "manifest safe unique name " + str(number))
        rows[name] = expected
    return rows


def verify_closed(directory: Path, *, expected_payload: set[str] | None = None):
    mode = os.lstat(directory).st_mode
    check(stat.S_ISDIR(mode) and not stat.S_ISLNK(mode), "ordinary directory")
    files = {}
    for root, dirs, names in os.walk(directory, followlinks=False):
        for name in dirs:
            p = Path(root) / name
            check(not p.is_symlink(), "no symlink directory " + str(p))
        for name in names:
            p = Path(root) / name
            check(not p.is_symlink() and stat.S_ISREG(os.lstat(p).st_mode),
                  "ordinary regular file " + str(p))
            files[p.relative_to(directory).as_posix()] = p
    manifest = manifest_rows(directory / "MANIFEST.sha256")
    visible = set(files) - {"MANIFEST.sha256", "SEAL.sha256"}
    check(set(manifest) == visible, "manifest exact payload census")
    if expected_payload is not None:
        check(set(manifest) == expected_payload, "expected payload census")
    for name, expected in manifest.items():
        check(digest(files[name]) == expected, "manifest content pin " + name)
    seal = manifest_rows(directory / "SEAL.sha256")
    check(seal == {"MANIFEST.sha256": digest(directory / "MANIFEST.sha256")},
          "seal binds manifest")
    return digest(directory / "MANIFEST.sha256"), digest(directory / "SEAL.sha256")


def pair(index: int):
    cursor = 0
    for i in range(CHANNELS):
        for j in range(i, CHANNELS):
            if cursor == index:
                return i, j
            cursor += 1
    raise AssertionError("unreachable ray index")


def reconstruct(values: list[float]):
    check(len(values) == UPPER and all(math.isfinite(x) for x in values),
          "finite 300 scalar values")
    out = np.zeros((CHANNELS, CHANNELS), dtype=np.float64)
    for ray, value in enumerate(values):
        i, j = pair(ray)
        if i == j:
            out[i, i] = value
    for ray, value in enumerate(values):
        i, j = pair(ray)
        if i != j:
            out[i, j] = out[j, i] = value - 0.5 * (out[i, i] + out[j, j])
    return out


def packed(matrix: np.ndarray):
    return [float(matrix[i, j]) for i in range(CHANNELS)
            for j in range(i, CHANNELS)]


def ao_transform():
    result = np.zeros((24, 24), dtype=np.float64)
    scale = 1.0 / math.sqrt(2.0)
    for sector in range(2):
        sign = 1.0 if sector == 0 else -1.0
        for pair_index in range(6):
            for quad in range(2):
                row = sector * 12 + 2 * pair_index + quad
                p = 2 * pair_index + quad
                c = 2 * (6 + pair_index) + quad
                result[row, p] = scale
                result[row, c] = sign * scale
    check(np.linalg.norm(result @ result.T - np.eye(24)) < 3e-15,
          "independent acoustic-optical transform orthogonal")
    return result


def as_matrix(rows):
    result = np.asarray(rows, dtype=np.float64)
    check(result.shape == (24, 24) and np.isfinite(result).all(),
          "finite 24 by 24 result matrix")
    return result


def close(left: np.ndarray, right: np.ndarray, label: str, tol: float = 6e-11):
    scale = max(float(np.linalg.norm(left)), float(np.linalg.norm(right)), 1.0)
    check(float(np.linalg.norm(left - right)) <= tol * scale, label)


def main() -> int:
    expected_result = {f"RAW/m001_ray_{i:03d}.json" for i in range(UPPER)} | {
        "FULL_RESPONSE_RECONSTRUCTION.json", "RAW_MANIFEST.sha256",
        "RAW_SEAL.sha256", "LAUNCH_CUSTODY.json",
        "reconstruct_m001_full_response_v002.py", "EXPECTED_HISTORY.json",
        "FROZEN_DEPENDENCIES.sha256", "FROZEN_DEPENDENCY_CONTRACT.json",
        "admission_contract_v003.py", "admission_contract_v002.py",
        "gl6fj_m001_runner_frozen.bin", "verify_result_packet.py",
    }
    result_manifest, result_seal = verify_closed(RESULT, expected_payload=expected_result)
    target_manifest, target_seal = verify_closed(TARGET)
    prelaunch_manifest, prelaunch_seal = verify_closed(PRELAUNCH_AUDIT)
    check(result_manifest == "09420cb679195fd2c68f976def1e4a861a98bc56408a3d2a0fafb270587719f9",
          "result manifest root fingerprint")
    check(result_seal == "6ab9053b9090d53d66e23356015d13009412359078c0c7d23df0634c7db37ca7",
          "result seal root fingerprint")
    check(target_manifest == "685844e9dbee9691fa6c576f2ffb6ec37a99827b29cc1a39e118aa9685b2c0b5",
          "target manifest root fingerprint")
    check(target_seal == "2204fe01c92bfcec176e183102addf65aef81f721cafe1aa92e9d0c9840b296b",
          "target seal root fingerprint")
    check(prelaunch_manifest == "4ed08f75432659d7b224c4f4b2fac405cc4b47cabb2226964b1a524523b6ff62",
          "prelaunch audit manifest root fingerprint")
    check(prelaunch_seal == "2f5dc997c8d099cc9603ea7544d3c5824c9535e751d8d1bc53a7c35b8830f0d4",
          "prelaunch audit seal root fingerprint")

    packet = strict_json(RESULT / "FULL_RESPONSE_RECONSTRUCTION.json")
    check(packet["status"] == "RECONSTRUCTED_SELECTED_COMPONENT_RESPONSE__NOT_PHYSICAL_WARD",
          "result claim ceiling status")
    check(packet["quotient"] == {
        "available_routine": "lawful_range_inverse",
        "optical_elimination_routine": "lawful_optical_schur",
        "performed": False,
        "reason": "NO_EXTERNALLY_CERTIFIED_SOURCE_VISIBLE_RANGE_OR_ERROR_BOUND",
    }, "no unauthorized quotient")
    check(packet["certified_ward_rank"] is None and
          packet["physical_ward_null_claimed"] is False,
          "no physical Ward promotion")
    check(packet["ceiling"] == "SELECTED_COMPONENT_L4_H6_M001_RESPONSE_ONLY__NO_M123_COMPOSITION_ORBIT_COMPLETION_REFINEMENT_PHYSICAL_WARD_EINSTEIN_GRAVITY_C_R_OR_G",
          "no gravity promotion")
    checkpoints = packet["checkpoints"]
    check(type(checkpoints) is list and len(checkpoints) == 18,
          "result has exact checkpoint count")

    raw = [strict_json(RESULT / "RAW" / f"m001_ray_{i:03d}.json") for i in range(UPPER)]
    reference_history = None
    for ray, record in enumerate(raw):
        i, j = pair(ray)
        check(record["ray_index"] == ray and record["ray_left"] == i and
              record["ray_right"] == j and record["character"] == [0, 0, 1],
              "ray catalog and m001 character " + str(ray))
        check(record["real_source_chart_unaliased"] is True and
              record["finite_run_anomalous_zero_imposed"] is False and
              record["physical_ward_rank_claimed"] is False,
              "raw claim boundary " + str(ray))
        history = (record["burn_digest"], tuple(
            (c["window"], c["path_steps"], c["decision_digest"],
             c["surviving_ancestors"], c["minimum_active"], c["maximum_active"])
            for c in record["checkpoints"]))
        if reference_history is None:
            reference_history = history
        else:
            check(history == reference_history, "same-state history " + str(ray))

    transform = ao_transform()
    worst_owner = 0.0
    worst_packet = 0.0
    final_eigenvalues = None
    final_ratio = None
    last_lineages = []
    for position, promoted in enumerate(checkpoints):
        source_rows = [record["checkpoints"][position] for record in raw]
        check((promoted["window"], promoted["path_steps"]) ==
              (source_rows[0]["window"], source_rows[0]["path_steps"]),
              "checkpoint position " + str(position))
        total = reconstruct([float(row["density_response"]) for row in source_rows])
        owner_sum = np.zeros((24, 24), dtype=np.float64)
        for owner in OWNERS:
            owner_sum += reconstruct([
                float(row["lambda_geometric"]) * float(row[owner]) /
                (float(row["path_steps"]) * 8.0 * 64.0 * (63.0 / 8.0))
                for row in source_rows])
        owner_error = float(np.linalg.norm(total - owner_sum) / max(np.linalg.norm(total), 1.0))
        worst_owner = max(worst_owner, owner_error)
        close(total, owner_sum, "six owner response recombination " + str(position))
        reported = np.asarray(promoted["density_total"]["packed_real_upper"], dtype=np.float64)
        expected = np.asarray(packed(total), dtype=np.float64)
        packet_error = float(np.linalg.norm(reported - expected) / max(np.linalg.norm(expected), 1.0))
        worst_packet = max(worst_packet, packet_error)
        check(reported.shape == (UPPER,) and np.isfinite(reported).all(),
              "reported packed finite response " + str(position))
        check(packet_error <= 6e-11, "independent packed response match " + str(position))
        reported_ao = as_matrix(promoted["density_total"]["real_ao_24x24"])
        close(reported_ao, transform @ total @ transform.T,
              "independent AO response match " + str(position))
        check(abs(float(promoted["raw_owner_total_recombination_residual"])) <= 3e-8,
              "reported scalar owner residual bounded " + str(position))
        check(abs(float(promoted["linear_rate_reconstruction_residual"])) == 0.0,
              "exact linear-rate reconstruction " + str(position))
        last_lineages.append(int(promoted["surviving_ancestors"])) if promoted["path_steps"] == 256 else None
        if position == len(checkpoints) - 1:
            final_eigenvalues = np.linalg.eigvalsh(total)
            final_ratio = float(promoted["density_total"]["anomalous_over_normal_frobenius"])

    check(final_eigenvalues is not None and final_ratio is not None,
          "final metrics available")
    check(min(last_lineages) >= 1 and len(last_lineages) == 6,
          "six final genealogy windows survive")
    output = {
        "schema": "GL6FJ_M001_INDEPENDENT_RESULT_AUDIT_V001",
        "result_packet": RESULT.name,
        "result_manifest_sha256": result_manifest,
        "result_seal_file_sha256": result_seal,
        "target_manifest_sha256": target_manifest,
        "target_seal_file_sha256": target_seal,
        "prelaunch_audit_manifest_sha256": prelaunch_manifest,
        "prelaunch_audit_seal_file_sha256": prelaunch_seal,
        "independent_checks": CHECKS,
        "disposition": "PASS__AUTHENTICATED_L4_M001_FINITE_SELECTED_COMPONENT_RESPONSE__NO_QUOTIENT_OR_PHYSICAL_WARD",
        "final_window": {"window": 5, "path_steps": 256},
        "final_response_eigenvalue_min": float(final_eigenvalues[0]),
        "final_response_eigenvalue_max": float(final_eigenvalues[-1]),
        "final_anomalous_over_normal_frobenius": final_ratio,
        "T256_surviving_ancestors_by_window": last_lineages,
        "worst_relative_six_owner_recombination_residual": worst_owner,
        "worst_relative_promoted_reconstruction_residual": worst_packet,
        "physical_ward_null_claimed": False,
        "gravity_claimed": False,
        "claim_ceiling": "FINITE_L4_M001_SELECTED_COMPONENT_RESPONSE_ONLY__NO_CERTIFIED_RANGE_1PI_SCHUR_PHYSICAL_WARD_EINSTEIN_IR_LIMIT_EXCHANGE_GRAVITY_C_R_OR_G",
    }
    print(json.dumps(output, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError, json.JSONDecodeError) as error:
        print("FAIL=" + str(error), file=sys.stderr)
        raise SystemExit(1)
