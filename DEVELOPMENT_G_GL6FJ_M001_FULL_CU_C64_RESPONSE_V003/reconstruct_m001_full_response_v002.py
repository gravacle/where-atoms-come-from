#!/usr/bin/env python3
"""GL6FJ V002 fixed-m001 reconstruction and quotient-ready transforms.

This program never launches the long calculation.  It consumes the 300
subordinate scalar JSON payloads only after they exist, authenticates their
same-state histories, reconstructs each raw checkpoint and physical owner,
and preserves the complete 24-real response.  A range inverse is available
only when an external certified range basis and error bound are supplied.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import stat
import sys
from typing import Any, Iterable

import numpy as np
import admission_contract_v002 as admission

sys.dont_write_bytecode = True

CHANNELS = 24
UPPER = CHANNELS * (CHANNELS + 1) // 2
MOMENTUM = (0, 0, 1)
MOMENTUM_LABEL = "m001"
QUADRATIC_OWNERS = (
    "stay_contact",
    "writer_contact",
    "stay_spectral",
    "writer_spectral",
    "stay_writer_interference",
    "disconnected",
)
LINEAR_OWNERS = ("stay_score_mean", "writer_score_mean", "first_mean")
DENOMINATOR = 8.0 * 64.0 * (63.0 / 8.0)
RAW_CEILING = ("SELECTED_COMPONENT_L4_H6_M001__FULL_REAL_RESPONSE_"
               "RECONSTRUCTIBLE_ONLY_AFTER_ALL_300_SAME_STATE_RAYS__NO_"
               "RANGE_QUOTIENT_OR_WARD_CLAIM")
EXPECTED_HISTORY_PATH = Path(__file__).resolve().parent / "EXPECTED_HISTORY.json"


Failure = admission.Failure


def check(condition: bool, message: str) -> None:
    if not condition:
        raise Failure(message)


def finite_tree(value: Any) -> None:
    admission.validate_finite_tree(value)


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Reject JSON objects whose duplicate names would otherwise be hidden."""
    result: dict[str, Any] = {}
    for key, value in pairs:
        check(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def strict_json_loads(payload: str) -> Any:
    return admission.strict_json_bytes(payload.encode("utf-8"))


def strict_json_file(path: Path) -> Any:
    return admission.strict_json(path)


HISTORY_CHECKPOINT_KEYS = (
    "window", "path_steps", "lambda_geometric", "stay_transitions",
    "writer_transitions", "surviving_ancestors", "minimum_active",
    "maximum_active", "decision_digest",
)


def source_zero_history(raw: dict[str, Any]) -> dict[str, Any]:
    return {
        "burn_stays": raw["burn_stays"],
        "burn_writers": raw["burn_writers"],
        "burn_digest": raw["burn_digest"],
        "checkpoints": [
            {key: checkpoint[key] for key in HISTORY_CHECKPOINT_KEYS}
            for checkpoint in raw["checkpoints"]
        ],
    }


def validated_expected_history() -> dict[str, Any]:
    payload = strict_json_file(EXPECTED_HISTORY_PATH)
    check(type(payload) is dict and set(payload) == {
        "schema", "source", "burn_digest", "burn_stays", "burn_writers",
        "checkpoints",
    }, "expected-history exact schema")
    check(payload["schema"] == "GL6FE_FROZEN_SOURCE_ZERO_HISTORY_V002" and
          payload["source"] ==
          "independently audited GL6FA and GL6FB c64 pilot results" and
          all(type(payload[key]) is int for key in (
              "burn_digest", "burn_stays", "burn_writers")) and
          type(payload["checkpoints"]) is list and
          len(payload["checkpoints"]) == 18,
          "expected-history exact outer types")
    for key in ("burn_digest", "burn_stays", "burn_writers"):
        admission.exact_u64(payload[key], "expected-history uint64 " + key)
    expected_order = [(window, steps) for window in range(6)
                      for steps in (64, 128, 256)]
    for position, row in enumerate(payload["checkpoints"]):
        check(type(row) is dict and set(row) == set(HISTORY_CHECKPOINT_KEYS),
              "expected-history checkpoint exact schema")
        integers = set(HISTORY_CHECKPOINT_KEYS) - {"lambda_geometric"}
        check(all(type(row[key]) is int for key in integers) and
              type(row["lambda_geometric"]) is float,
              "expected-history checkpoint exact types")
        for key in integers:
            admission.exact_u64(
                row[key], "expected-history checkpoint uint64 " + key)
        admission.exact_number(
            row["lambda_geometric"],
            "expected-history finite positive lambda", lower=0.0)
        check((row["window"], row["path_steps"]) == expected_order[position] and
              row["stay_transitions"] >= 0 and
              row["writer_transitions"] >= 0 and
              row["stay_transitions"] + row["writer_transitions"] ==
              1024 * row["path_steps"] and
              1 <= row["surviving_ancestors"] <= 1024 and
              0 < row["minimum_active"] <= row["maximum_active"] <= 64 and
              row["decision_digest"] >= 0 and
              row["lambda_geometric"] > 0.0,
              "expected-history checkpoint domains")
    return {
        "burn_stays": payload["burn_stays"],
        "burn_writers": payload["burn_writers"],
        "burn_digest": payload["burn_digest"],
        "checkpoints": payload["checkpoints"],
    }


def upper_index(left: int, right: int) -> int:
    check(0 <= left <= right < CHANNELS, "upper-index domain")
    return left * CHANNELS - left * (left - 1) // 2 + right - left


def ray_pair(index: int) -> tuple[int, int]:
    check(0 <= index < UPPER, "ray index domain")
    cursor = 0
    for left in range(CHANNELS):
        for right in range(left, CHANNELS):
            if cursor == index:
                return left, right
            cursor += 1
    raise AssertionError("unreachable")


def canonical_direction(index: int) -> np.ndarray:
    left, right = ray_pair(index)
    result = np.zeros(CHANNELS)
    if left == right:
        result[left] = 1.0
    else:
        result[[left, right]] = 1.0 / math.sqrt(2.0)
    return result


def reconstruct_quadratic(values: Iterable[float]) -> np.ndarray:
    """Reconstruct H from q(e_i) and q((e_i+e_j)/sqrt(2))."""
    rows = np.asarray(tuple(values), dtype=np.float64)
    check(rows.shape == (UPPER,) and np.isfinite(rows).all(),
          "quadratic projection shape/finite")
    result = np.zeros((CHANNELS, CHANNELS), dtype=np.float64)
    for index in range(UPPER):
        left, right = ray_pair(index)
        if left == right:
            result[left, left] = rows[index]
    for index in range(UPPER):
        left, right = ray_pair(index)
        if left != right:
            value = rows[index] - 0.5 * (
                result[left, left] + result[right, right])
            result[left, right] = result[right, left] = value
    check(np.isfinite(result).all(), "finite reconstructed Hessian")
    return result


def pack_upper(matrix: np.ndarray) -> list[float]:
    array = np.asarray(matrix, dtype=np.float64)
    check(array.shape == (CHANNELS, CHANNELS), "pack shape")
    check(np.linalg.norm(array - array.T) <=
          2e-12 * max(np.linalg.norm(array), 1.0), "pack symmetric")
    return [float(array[left, right]) for left in range(CHANNELS)
            for right in range(left, CHANNELS)]


def reconstruct_linear(values: Iterable[float], tolerance: float = 2e-8
                       ) -> tuple[np.ndarray, float]:
    rows = np.asarray(tuple(values), dtype=np.float64)
    check(rows.shape == (UPPER,) and np.isfinite(rows).all(),
          "linear projection shape/finite")
    vector = np.asarray(
        [rows[upper_index(index, index)] for index in range(CHANNELS)])
    residual = 0.0
    for index in range(UPPER):
        left, right = ray_pair(index)
        expected = (vector[left] if left == right else
                    (vector[left] + vector[right]) / math.sqrt(2.0))
        residual = max(residual, abs(float(rows[index]) - expected))
    scale = max(float(np.max(np.abs(rows))), 1.0)
    check(residual <= tolerance * scale, "same-state first-score polarization")
    return vector, float(residual)


def complex_payload(matrix: np.ndarray) -> list[list[list[float]]]:
    array = np.asarray(matrix, dtype=np.complex128)
    check(np.isfinite(array.real).all() and np.isfinite(array.imag).all(),
          "finite complex matrix")
    return [[[float(value.real), float(value.imag)] for value in row]
            for row in array]


def real_payload(matrix: np.ndarray) -> list[list[float]]:
    array = np.asarray(matrix, dtype=np.float64)
    check(np.isfinite(array).all(), "finite real matrix")
    return [[float(value) for value in row] for row in array]


def normal_anomalous(real_hessian: np.ndarray
                     ) -> tuple[np.ndarray, np.ndarray]:
    raw = np.asarray(real_hessian, dtype=np.float64)
    check(raw.shape == (24, 24), "real Hessian dimension")
    normal = np.empty((12, 12), dtype=np.complex128)
    anomalous = np.empty((12, 12), dtype=np.complex128)
    for left in range(12):
        ql, il = 2 * left, 2 * left + 1
        for right in range(12):
            qr, ir = 2 * right, 2 * right + 1
            normal[left, right] = (
                raw[ql, qr] + raw[il, ir]
                + 1j * (raw[il, qr] - raw[ql, ir]))
            anomalous[left, right] = (
                raw[ql, qr] - raw[il, ir]
                + 1j * (raw[il, qr] + raw[ql, ir]))
    return normal, anomalous


def real_from_normal_anomalous(normal: np.ndarray,
                               anomalous: np.ndarray) -> np.ndarray:
    """Inverse of ``normal_anomalous`` in the fixed q/iq convention."""
    n = np.asarray(normal, dtype=np.complex128)
    a = np.asarray(anomalous, dtype=np.complex128)
    check(n.shape == (12, 12) and a.shape == (12, 12),
          "normal/anomalous inverse dimensions")
    check(np.isfinite(n.real).all() and np.isfinite(n.imag).all() and
          np.isfinite(a.real).all() and np.isfinite(a.imag).all(),
          "normal/anomalous inverse finite")
    raw = np.empty((24, 24), dtype=np.float64)
    for left in range(12):
        ql, il = 2 * left, 2 * left + 1
        for right in range(12):
            qr, ir = 2 * right, 2 * right + 1
            raw[ql, qr] = (n[left, right].real + a[left, right].real) / 2.0
            raw[il, ir] = (n[left, right].real - a[left, right].real) / 2.0
            raw[il, qr] = (n[left, right].imag + a[left, right].imag) / 2.0
            raw[ql, ir] = (a[left, right].imag - n[left, right].imag) / 2.0
    return raw


def complex_ao_transform() -> np.ndarray:
    identity = np.eye(6)
    return np.block([[identity, identity],
                     [identity, -identity]]) / math.sqrt(2.0)


def real_ao_transform() -> np.ndarray:
    """Rows: AC pair0 Q/I...pair5 Q/I, then OP in same order."""
    result = np.zeros((24, 24), dtype=np.float64)
    scale = 1.0 / math.sqrt(2.0)
    for sector in range(2):
        sign = 1.0 if sector == 0 else -1.0
        for pair in range(6):
            for quadrature in range(2):
                row = sector * 12 + 2 * pair + quadrature
                p = 2 * pair + quadrature
                c = 2 * (6 + pair) + quadrature
                result[row, p] = scale
                result[row, c] = sign * scale
    check(np.linalg.norm(result @ result.T - np.eye(24)) < 3e-15,
          "real AO transform orthogonal")
    return result


def transformed_payload(matrix: np.ndarray) -> dict[str, Any]:
    raw = np.asarray(matrix, dtype=np.float64)
    normal, anomalous = normal_anomalous(raw)
    cao = complex_ao_transform()
    rao = real_ao_transform()
    normal_ao = cao @ normal @ cao.T
    anomalous_ao = cao @ anomalous @ cao.T
    real_ao = rao @ raw @ rao.T
    return {
        "packed_real_upper": pack_upper(raw),
        "normal_12x12": complex_payload(normal),
        "anomalous_12x12": complex_payload(anomalous),
        "normal_ao_12x12": complex_payload(normal_ao),
        "anomalous_ao_12x12": complex_payload(anomalous_ao),
        "real_ao_24x24": real_payload(real_ao),
        "real_ao_blocks": {
            "acoustic_acoustic_12x12": real_payload(real_ao[:12, :12]),
            "acoustic_optical_12x12": real_payload(real_ao[:12, 12:]),
            "optical_acoustic_12x12": real_payload(real_ao[12:, :12]),
            "optical_optical_12x12": real_payload(real_ao[12:, 12:]),
        },
        "anomalous_symmetry_residual": float(
            np.linalg.norm(anomalous - anomalous.T)),
        "anomalous_over_normal_frobenius": float(
            np.linalg.norm(anomalous) / max(np.linalg.norm(normal), 1e-300)),
    }


def lawful_range_inverse(response: np.ndarray, range_basis: np.ndarray,
                         certified_operator_error: float
                         ) -> tuple[np.ndarray, dict[str, float]]:
    """Invert only an externally certified source-visible response range.

    This routine does not discover a range from the same noisy point estimate.
    The caller must supply the audited basis and an operator-norm error bound.
    """
    x = np.asarray(response, dtype=np.float64)
    u = np.asarray(range_basis, dtype=np.float64)
    check(x.shape == (24, 24) and u.ndim == 2 and u.shape[0] == 24,
          "range inverse dimensions")
    check(np.isfinite(x).all() and np.isfinite(u).all(),
          "range inverse finite inputs")
    check(certified_operator_error >= 0 and
          math.isfinite(certified_operator_error), "certified error bound")
    symmetry_residual = float(np.linalg.norm(x - x.T, 2))
    check(symmetry_residual <= certified_operator_error,
          "response symmetry within certified operator error")
    orthogonal_residual = float(np.linalg.norm(u.T @ u - np.eye(u.shape[1])))
    check(orthogonal_residual <= 2e-10, "certified basis orthonormality")
    projector = u @ u.T
    complement = np.eye(24) - projector
    left_complement = float(np.linalg.norm(complement @ x, 2))
    right_complement = float(np.linalg.norm(x @ complement, 2))
    check(left_complement <= certified_operator_error and
          right_complement <= certified_operator_error,
          "full response complement exceeds certified range error")
    reduced = u.T @ x @ u
    singular = np.linalg.svd(reduced, compute_uv=False)
    minimum = float(singular[-1]) if singular.size else math.inf
    check(minimum > certified_operator_error,
          "range block not separated from estimator error")
    gamma = u @ np.linalg.inv(reduced) @ u.T
    check(np.isfinite(gamma).all(), "finite certified range inverse")
    inverse_error = (certified_operator_error /
                     (minimum * (minimum - certified_operator_error))
                     if certified_operator_error else 0.0)
    return gamma, {
        "response_symmetry_operator_norm": symmetry_residual,
        "orthogonal_residual": orthogonal_residual,
        "left_complement_operator_norm": left_complement,
        "right_complement_operator_norm": right_complement,
        "reduced_minimum_singular": minimum,
        "certified_operator_error": certified_operator_error,
        "certified_range_inverse_error_bound": inverse_error,
        "range_inverse_residual": float(
            np.linalg.norm(projector @ (x @ gamma) - projector)),
    }


def lawful_optical_schur(kernel_real_ao: np.ndarray,
                          certified_full_kernel_operator_error: float
                          ) -> tuple[np.ndarray, dict[str, float]]:
    """Eliminate optical modes under a full-kernel operator error bound."""
    kernel = np.asarray(kernel_real_ao, dtype=np.float64)
    check(kernel.shape == (24, 24), "AO kernel dimension")
    check(np.isfinite(kernel).all() and
          math.isfinite(certified_full_kernel_operator_error),
          "finite Schur inputs")
    epsilon = certified_full_kernel_operator_error
    check(epsilon >= 0.0, "nonnegative full-kernel operator error")
    symmetry_residual = float(np.linalg.norm(kernel - kernel.T, 2))
    check(symmetry_residual <= epsilon,
          "AO kernel symmetry within certified error")
    aa, ao = kernel[:12, :12], kernel[:12, 12:]
    oa, oo = kernel[12:, :12], kernel[12:, 12:]
    block_transpose_residual = float(np.linalg.norm(ao - oa.T, 2))
    check(block_transpose_residual <= epsilon,
          "AO/OA transpose relation within certified error")
    singular = np.linalg.svd(oo, compute_uv=False)
    minimum = float(singular[-1])
    check(minimum > epsilon,
          "optical block not certified invertible")
    inverse = np.linalg.inv(oo)
    schur = aa - ao @ inverse @ oa
    check(np.isfinite(inverse).all() and np.isfinite(schur).all(),
          "finite optical inverse and Schur complement")
    norm_ao = float(np.linalg.norm(ao, 2))
    norm_oa = float(np.linalg.norm(oa, 2))
    inverse_error = (epsilon / (minimum * (minimum - epsilon))
                     if epsilon else 0.0)
    schur_error = (
        epsilon
        + epsilon / (minimum - epsilon) * (norm_oa + epsilon)
        + norm_ao * inverse_error * (norm_oa + epsilon)
        + norm_ao / minimum * epsilon
        if epsilon else 0.0)
    schur_symmetry = float(np.linalg.norm(schur - schur.T, 2))
    roundoff = (64.0 * np.finfo(np.float64).eps
                * max(float(np.linalg.norm(schur, 2)), 1.0))
    check(schur_symmetry <= 2.0 * schur_error + roundoff,
          "Schur symmetry within propagated full-kernel error")
    return schur, {
        "kernel_symmetry_operator_norm": symmetry_residual,
        "ao_oa_transpose_operator_norm": block_transpose_residual,
        "optical_minimum_singular": minimum,
        "certified_full_kernel_operator_error": epsilon,
        "certified_optical_inverse_error_bound": inverse_error,
        "certified_schur_error_bound": schur_error,
        "schur_symmetry_operator_norm": schur_symmetry,
        "floating_point_roundoff_allowance": roundoff,
        "optical_inverse_residual": float(
            np.linalg.norm(oo @ inverse - np.eye(12))),
    }


def load_raw_directory(path: Path) -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}
    mode = os.lstat(path).st_mode
    check(stat.S_ISDIR(mode) and not stat.S_ISLNK(mode),
          "raw input is an ordinary unlinked directory")
    entries = list(os.scandir(path))
    check(all(stat.S_ISREG(entry.stat(follow_symlinks=False).st_mode)
              and not entry.is_symlink() for entry in entries),
          "closed raw directory regular unlinked files only")
    sources = sorted((Path(entry.path) for entry in entries), key=lambda row: row.name)
    expected_names = {f"m001_ray_{index:03d}.json" for index in range(UPPER)}
    check({source.name for source in sources} == expected_names,
          "closed raw directory exact 300-file census")
    for source in sources:
        try:
            raw = strict_json_file(source)
        except (json.JSONDecodeError, UnicodeDecodeError) as error:
            raise Failure(f"invalid raw JSON {source}: {error}") from error
        check(type(raw) is dict and
              raw.get("schema") == "GL6FJ_SCALAR_PROJECTION_RAW_V001",
              "exact GL6FJ canonical raw schema")
        index = raw.get("ray_index")
        check(type(index) is int and 0 <= index < UPPER,
              "raw ray index domain")
        check(source.name == f"m001_ray_{index:03d}.json",
              "raw filename/index identity")
        check(index not in rows, "duplicate canonical projection")
        rows[index] = raw
    check(set(rows) == set(range(UPPER)), "exact 300 canonical census")
    return rows


def validate_common_raw(rows: dict[int, dict[str, Any]]) -> None:
    check(type(rows) is dict and set(rows) == set(range(UPPER)),
          "exact 300 canonical raw census")
    outer_keys = {
        "schema", "mode", "L", "volume", "component", "action_order",
        "momentum_label", "character", "channels", "symmetric_entries",
        "design_kind", "ray_index", "ray_left", "ray_right", "ray_kind",
        "population", "burn", "window_length", "windows", "stay_units",
        "writer_units", "stay_weight", "random_seed",
        "support_cache_dense_entries", "burn_stays", "burn_writers",
        "burn_digest", "checkpoints", "support_cache_hits",
        "support_cache_misses", "writer_cache_hits", "writer_cache_misses",
        "sampled_cold_residual", "production_execution",
        "physical_ward_rank_claimed",
        "translation_forced_anomalous_expectation_zero",
        "finite_run_anomalous_zero_imposed", "real_source_chart_dimension",
        "real_source_chart_unaliased", "ceiling",
    }
    checkpoint_keys = {
        "window", "path_steps", "lambda_geometric", "stay_transitions",
        "writer_transitions", "surviving_ancestors", "minimum_active",
        "maximum_active", "decision_digest", "stay_score_mean",
        "writer_score_mean", "first_mean", "first_score_rate",
        "stay_contact", "writer_contact", "stay_spectral",
        "writer_spectral", "stay_writer_interference", "disconnected",
        "owner_total", "log_hessian_rate", "raw_connected_perron_hessian",
        "density_response",
    }
    exact_history = validated_expected_history()
    reference: dict[str, Any] | None = None
    reference_cache: list[int] | None = None
    expected_checkpoint_pairs = [(window, steps) for window in range(6)
                                 for steps in (64, 128, 256)]
    for ray, raw in rows.items():
        finite_tree(raw)
        check(type(raw) is dict and set(raw) == outer_keys,
              "exact raw schema keys")
        left, right = ray_pair(ray)
        exact = {
            "schema": "GL6FJ_SCALAR_PROJECTION_RAW_V001",
            "mode": "subordinate_raw", "L": 4, "volume": 64,
            "component": "GL6CC_SELECTED_DENSE_H6_COMPONENT",
            "action_order": "through_h6", "momentum_label": "m001",
            "character": [0, 0, 1], "channels": 24,
            "symmetric_entries": 300,
            "design_kind": "canonical_full_response",
            "ray_index": ray, "ray_left": left, "ray_right": right,
            "ray_kind": ("unit" if left == right else "normalized_pair_sum"),
            "population": 1024, "burn": 2560, "window_length": 256,
            "windows": 6, "stay_units": 4096, "writer_units": 64,
            "stay_weight": 64, "random_seed": 96511001,
            "support_cache_dense_entries": 4008960,
            "production_execution": False,
            "physical_ward_rank_claimed": False,
            "translation_forced_anomalous_expectation_zero": True,
            "finite_run_anomalous_zero_imposed": False,
            "real_source_chart_dimension": 24,
            "real_source_chart_unaliased": True,
            "ceiling": RAW_CEILING,
        }
        for key, value in exact.items():
            admission.deep_exact(raw.get(key), value,
                                 f"raw fixed field {key} at ray {ray}")
        integer_fields = {
            "L", "volume", "channels", "symmetric_entries", "ray_index",
            "ray_left", "ray_right", "population", "burn", "window_length",
            "windows", "stay_units", "writer_units", "stay_weight",
            "random_seed", "support_cache_dense_entries", "burn_stays",
            "burn_writers", "burn_digest", "support_cache_hits",
            "support_cache_misses", "writer_cache_hits", "writer_cache_misses",
            "real_source_chart_dimension",
        }
        check(all(type(raw[key]) is int for key in integer_fields),
              "raw exact integer fields excluding Boolean")
        string_fields = {"schema", "mode", "component", "action_order",
                         "momentum_label", "design_kind", "ray_kind", "ceiling"}
        check(all(type(raw[key]) is str for key in string_fields),
              "raw exact string fields")
        boolean_fields = {"production_execution", "physical_ward_rank_claimed",
                          "translation_forced_anomalous_expectation_zero",
                          "finite_run_anomalous_zero_imposed",
                          "real_source_chart_unaliased"}
        check(all(type(raw[key]) is bool for key in boolean_fields),
              "raw exact Boolean fields")
        check(type(raw["character"]) is list and
              all(type(value) is int for value in raw["character"]),
              "raw character exact integer vector")
        admission.exact_number(raw["sampled_cold_residual"],
                               "bounded finite cold residual",
                               lower=0.0, upper=3e-7)
        for key in ("burn_stays", "burn_writers", "burn_digest",
                    "support_cache_hits", "support_cache_misses",
                    "writer_cache_hits", "writer_cache_misses"):
            admission.exact_u64(raw[key], "uint64 counter/digest " + key)
        cache = [raw[key] for key in (
            "support_cache_dense_entries", "support_cache_hits",
            "support_cache_misses", "writer_cache_hits",
            "writer_cache_misses")]
        if reference_cache is None:
            reference_cache = cache
        else:
            admission.deep_exact(cache, reference_cache,
                                 "all 300 same-history cache/access counters")
        check(raw["burn_stays"] + raw["burn_writers"] == 1024 * 2560,
              "burn transition census")
        checkpoints = raw.get("checkpoints")
        check(type(checkpoints) is list and len(checkpoints) == 18,
              "exact 18-checkpoint census")
        observed_pairs = []
        for checkpoint in checkpoints:
            check(type(checkpoint) is dict and set(checkpoint) == checkpoint_keys,
                  "exact checkpoint schema keys")
            integer_checkpoint = {
                "window", "path_steps", "stay_transitions",
                "writer_transitions", "surviving_ancestors", "minimum_active",
                "maximum_active", "decision_digest",
            }
            check(all(type(checkpoint[key]) is int for key in integer_checkpoint),
                  "checkpoint exact integer fields")
            for key in checkpoint_keys - integer_checkpoint:
                admission.exact_number(checkpoint[key],
                                       "checkpoint finite number " + key)
            for key in integer_checkpoint:
                admission.exact_u64(checkpoint[key],
                                    "checkpoint uint64 field " + key)
            observed_pairs.append((checkpoint["window"], checkpoint["path_steps"]))
            steps = checkpoint["path_steps"]
            check(checkpoint["stay_transitions"] >= 0 and
                  checkpoint["writer_transitions"] >= 0 and
                  checkpoint["stay_transitions"] + checkpoint["writer_transitions"]
                  == 1024 * steps, "checkpoint transition census")
            check(1 <= checkpoint["surviving_ancestors"] <= 1024 and
                  0 < checkpoint["minimum_active"] <=
                  checkpoint["maximum_active"] <= 64 and
                  checkpoint["decision_digest"] >= 0 and
                  checkpoint["lambda_geometric"] > 0.0,
                  "checkpoint genealogy/activity/growth domains")
            owner_scale = max(abs(float(checkpoint["owner_total"])), 1.0)
            owner_sum = sum(float(checkpoint[key]) for key in QUADRATIC_OWNERS)
            check(abs(owner_sum - float(checkpoint["owner_total"])) <=
                  3e-11 * owner_scale, "scalar owner-once recombination")
            first_scale = max(abs(float(checkpoint["first_mean"])), 1.0)
            check(abs(float(checkpoint["first_mean"]) -
                      float(checkpoint["stay_score_mean"]) -
                      float(checkpoint["writer_score_mean"])) <=
                  3e-11 * first_scale, "scalar first-score owner identity")
            check(abs(float(checkpoint["first_score_rate"]) -
                      float(checkpoint["first_mean"]) / steps) <=
                  3e-11 * max(first_scale / steps, 1.0),
                  "first-score rate identity")
            check(abs(float(checkpoint["disconnected"]) +
                      float(checkpoint["first_mean"]) ** 2) <=
                  3e-11 * max(float(checkpoint["first_mean"]) ** 2, 1.0),
                  "disconnected owner identity")
            check(abs(float(checkpoint["log_hessian_rate"]) -
                      float(checkpoint["owner_total"]) / steps) <=
                  3e-11 * max(owner_scale / steps, 1.0),
                  "log-Hessian rate identity")
            expected_raw = (float(checkpoint["lambda_geometric"]) *
                            float(checkpoint["owner_total"]) / steps)
            check(abs(float(checkpoint["raw_connected_perron_hessian"]) -
                      expected_raw) <= 3e-11 * max(abs(expected_raw), owner_scale),
                  "Perron Hessian identity")
            expected_density = expected_raw / DENOMINATOR
            check(abs(float(checkpoint["density_response"]) - expected_density) <=
                  3e-11 * max(abs(expected_density), 1.0),
                  "density normalization identity")
        check(observed_pairs == expected_checkpoint_pairs,
              "exact window/checkpoint ordering")
        history = source_zero_history(raw)
        admission.deep_exact(history, exact_history,
                             "exact frozen source-zero history")
        if reference is None:
            reference = history
        else:
            admission.deep_exact(history, reference,
                                 "all 300 projections share exact histories")


def reconstruct(rows: dict[int, dict[str, Any]]) -> dict[str, Any]:
    validate_common_raw(rows)
    result: dict[str, Any] = {
        "schema": "GL6FJ_M001_FULL_RESPONSE_RECONSTRUCTION_V002",
        "status": "RECONSTRUCTED_SELECTED_COMPONENT_RESPONSE__NOT_PHYSICAL_WARD",
        "design": {
            "channels": 24, "symmetric_entries": 300,
            "momentum_label": "m001", "character": [0, 0, 1],
            "self_conjugate": False, "same_state_history_required": True,
            "real_source_chart_dimension": 24,
            "real_source_chart_unaliased": True,
            "translation_forced_anomalous_expectation_zero": True,
            "finite_run_anomalous_zero_imposed": False,
            "reconstruction_stage":
                "EACH_RAW_WINDOW_AND_CHECKPOINT_BEFORE_AGGREGATION",
        },
        "checkpoints": [],
        "quotient": {
            "performed": False,
            "reason": "NO_EXTERNALLY_CERTIFIED_SOURCE_VISIBLE_RANGE_OR_ERROR_BOUND",
            "available_routine": "lawful_range_inverse",
            "optical_elimination_routine": "lawful_optical_schur",
        },
        "certified_ward_rank": None,
        "physical_ward_null_claimed": False,
        "ceiling": (
            "SELECTED_COMPONENT_L4_H6_M001_RESPONSE_ONLY__NO_M123_COMPOSITION_"
            "ORBIT_COMPLETION_REFINEMENT_PHYSICAL_WARD_EINSTEIN_GRAVITY_C_R_OR_G"
        ),
    }
    first_raw = rows[0]
    for checkpoint_index, reference in enumerate(first_raw["checkpoints"]):
        quadratic: dict[str, np.ndarray] = {}
        for owner in QUADRATIC_OWNERS + ("owner_total", "density_response"):
            quadratic[owner] = reconstruct_quadratic(
                rows[ray]["checkpoints"][checkpoint_index][owner]
                for ray in range(UPPER))
        owner_sum = sum((quadratic[key] for key in QUADRATIC_OWNERS),
                        np.zeros((24, 24)))
        owner_residual = float(np.linalg.norm(owner_sum - quadratic["owner_total"]))
        owner_scale = max(float(np.linalg.norm(quadratic["owner_total"])), 1.0)
        check(owner_residual <= 2e-8 * owner_scale,
              "owner-once matrix recombination")
        scale = (float(reference["lambda_geometric"]) /
                 (float(reference["path_steps"]) * DENOMINATOR))
        density_residual = float(np.linalg.norm(
            scale * quadratic["owner_total"] - quadratic["density_response"]))
        density_scale = max(float(np.linalg.norm(quadratic["density_response"])), 1.0)
        check(density_residual <= 2e-8 * density_scale,
              "density normalization matrix reconstruction")

        linear_vectors: dict[str, np.ndarray] = {}
        linear: dict[str, Any] = {}
        for owner in LINEAR_OWNERS + ("first_score_rate",):
            vector, residual = reconstruct_linear(
                rows[ray]["checkpoints"][checkpoint_index][owner]
                for ray in range(UPPER))
            linear_vectors[owner] = vector
            linear[owner] = {
                "vector_24": [float(value) for value in vector],
                "maximum_pair_relation_residual": residual,
            }
        linear_owner_residual = float(np.linalg.norm(
            linear_vectors["first_mean"] - linear_vectors["stay_score_mean"] -
            linear_vectors["writer_score_mean"]))
        linear_rate_residual = float(np.linalg.norm(
            linear_vectors["first_score_rate"] -
            linear_vectors["first_mean"] / float(reference["path_steps"])))
        linear_scale = max(float(np.linalg.norm(linear_vectors["first_mean"])), 1.0)
        check(linear_owner_residual <= 2e-8 * linear_scale,
              "linear owner matrix identity")
        check(linear_rate_residual <= 2e-8 * linear_scale,
              "linear rate matrix identity")

        result["checkpoints"].append({
            "window": reference["window"],
            "path_steps": reference["path_steps"],
            "lambda_geometric": reference["lambda_geometric"],
            "surviving_ancestors": reference["surviving_ancestors"],
            "minimum_active": reference["minimum_active"],
            "maximum_active": reference["maximum_active"],
            "decision_digest": reference["decision_digest"],
            "linear_scores": linear,
            "linear_owner_recombination_residual": linear_owner_residual,
            "linear_rate_reconstruction_residual": linear_rate_residual,
            "raw_owner_total_recombination_residual": owner_residual,
            "density_normalization_reconstruction_residual": density_residual,
            "density_total": transformed_payload(quadratic["density_response"]),
            "density_owners": {
                owner: transformed_payload(scale * quadratic[owner])
                for owner in QUADRATIC_OWNERS
            },
            "owner_typing": (
                "SAME_CHECKPOINT_OWNER_ONCE_LINEAR_RECONSTRUCTION__"
                "NO_COMPONENTWISE_MEDIAN_RECOMBINATION"
            ),
            "anomalous_typing": (
                "TRANSLATION_FORCES_ZERO_IN_THE_EXACT_EXPECTATION__FINITE_RUN_"
                "ESTIMATE_RETAINED_AS_MEASURED_CONTROL_NOT_PROJECTED_TO_ZERO"
            ),
        })
    finite_tree(result)
    return result


def selftest() -> dict[str, Any]:
    rng = np.random.default_rng(96512001)
    matrix = rng.normal(size=(24, 24))
    matrix = (matrix + matrix.T) / 2.0
    first = rng.normal(size=24)
    q = np.zeros(UPPER)
    linear = np.zeros(UPPER)
    for index in range(UPPER):
        left, right = ray_pair(index)
        direction = np.zeros(24)
        if left == right:
            direction[left] = 1.0
        else:
            direction[[left, right]] = 1.0 / math.sqrt(2.0)
        q[index] = direction @ matrix @ direction
        linear[index] = direction @ first
    recovered = reconstruct_quadratic(q)
    recovered_first, first_residual = reconstruct_linear(linear, 2e-13)
    reconstruction_residual = float(np.linalg.norm(recovered - matrix))
    check(reconstruction_residual < 2e-13, "synthetic matrix reconstruction")
    check(np.linalg.norm(recovered_first - first) < 2e-13,
          "synthetic first reconstruction")

    normal, anomalous = normal_anomalous(matrix)
    normal_hermiticity = float(np.linalg.norm(normal - normal.conj().T))
    anomalous_symmetry = float(np.linalg.norm(anomalous - anomalous.T))
    transform_inverse = float(np.linalg.norm(
        real_from_normal_anomalous(normal, anomalous) - matrix))
    check(normal_hermiticity < 2e-13, "normal block Hermitian")
    check(anomalous_symmetry < 2e-13, "anomalous block symmetric")
    check(transform_inverse < 2e-13, "normal/anomalous transform invertible")
    rao = real_ao_transform()
    real_ao_inverse = float(np.linalg.norm(rao.T @ (rao @ matrix @ rao.T) @ rao
                                           - matrix))
    cao = complex_ao_transform()
    complex_ao_inverse = float(np.linalg.norm(
        cao.T @ (cao @ normal @ cao.T) @ cao - normal))
    check(real_ao_inverse < 2e-13, "real AO transform invertible")
    check(complex_ao_inverse < 2e-13, "complex AO transform invertible")

    # Test the lawful range inverse and AO Schur on a known rank-7 response.
    raw_basis, _ = np.linalg.qr(rng.normal(size=(24, 7)))
    response = raw_basis @ np.diag(np.arange(1.0, 8.0)) @ raw_basis.T
    gamma, diagnostics = lawful_range_inverse(response, raw_basis, 1e-10)
    projector = raw_basis @ raw_basis.T
    range_residual = float(np.linalg.norm(response @ gamma - projector))
    check(range_residual < 2e-12, "synthetic lawful range inverse")
    kernel_ao = np.diag(np.arange(1.0, 25.0))
    schur, schur_diagnostics = lawful_optical_schur(kernel_ao, 1e-10)
    check(np.linalg.norm(schur - np.diag(np.arange(1.0, 13.0))) < 1e-13,
          "synthetic optical Schur")
    transformed = transformed_payload(matrix)
    return {
        "schema": "GL6FJ_RECONSTRUCTION_SELFTEST_V002",
        "ray_count": UPPER,
        "reconstruction_residual": reconstruction_residual,
        "redundant_diagnostic_rows": 0,
        "first_pair_residual": first_residual,
        "normal_hermiticity_residual": normal_hermiticity,
        "anomalous_symmetry_residual": anomalous_symmetry,
        "normal_anomalous_inverse_residual": transform_inverse,
        "real_ao_inverse_residual": real_ao_inverse,
        "complex_ao_inverse_residual": complex_ao_inverse,
        "range_inverse_residual": range_residual,
        "range_diagnostics": diagnostics,
        "schur_diagnostics": schur_diagnostics,
        "normal_shape": [len(transformed["normal_12x12"]),
                         len(transformed["normal_12x12"][0])],
        "result": "PASS",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("selftest")
    reconstruction = subparsers.add_parser("reconstruct")
    reconstruction.add_argument("input_directory", type=Path)
    reconstruction.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.command == "selftest":
        print(json.dumps(selftest(), sort_keys=True, indent=2))
        return 0
    rows = load_raw_directory(args.input_directory)
    payload = reconstruct(rows)
    check(not args.output.exists(), "write-once reconstruction output")
    args.output.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
