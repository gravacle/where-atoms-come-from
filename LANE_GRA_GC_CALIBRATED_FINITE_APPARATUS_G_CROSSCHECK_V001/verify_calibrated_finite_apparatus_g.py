#!/usr/bin/env python3
"""Deterministic checks for the calibrated finite-apparatus G cross-check."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from finite_apparatus_g_model import (  # noqa: E402
    MechanicalCalibration,
    connected_grid_interval,
    effective_covariance,
    forward_complex,
    full_coupled_complex,
    geometry_coefficients,
    point_potential_per_G,
    profile_product_grid,
    profiled_quadratic,
    source_scale_identified_interval,
    stack_complex,
    stack_complex_columns,
    static_stability_margin,
    uniform_sphere_coefficients,
)


# These absolute search and covariance-linearization domains are frozen without
# reference to the synthetic generator below.
PRODUCT_SCAN_SI = (1.0e-12, 1.5e-10)
PRODUCT_SCAN_POINTS = 6001
SOURCE_SCALE_CALIBRATION_INTERVAL = (0.9990, 1.0010)
# Frozen before the synthetic truth is generated.  It is used only to
# linearize independently measured gain/delay uncertainty into covariance.
COVARIANCE_LINEARIZATION_PRODUCT = 7.0e-11

# This is an arbitrary hidden-generator value, not CODATA, an accepted-G input,
# or a value used by the estimator.  It is disclosed only after estimation.
SYNTHETIC_TRUTH_G = 7.314159265358979e-11
SYNTHETIC_SOURCE_SCALE = 1.0006


checks: list[tuple[str, bool, str]] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    checks.append((name, bool(condition), detail))


def source_cloud(phi: float):
    radius = 0.22
    angles = np.array((phi, phi + np.pi))
    positions = np.column_stack(
        (
            radius * np.cos(angles),
            radius * np.sin(angles),
            np.zeros(angles.size),
        )
    )
    return np.array((11.0, 11.0)), positions


detector_masses = np.array((0.25, 0.25))
detector_positions = np.array(((0.07, 0.0, 0.0), (-0.07, 0.0, 0.0)))


# ---------------------------------------------------------------------------
# Exact finite-source kernels and their two required derivatives.
phi_test = 0.47
source_masses, source_positions = source_cloud(phi_test)
dphi = 0.04
z_cross_source = np.column_stack(
    (-source_positions[:, 1], source_positions[:, 0], np.zeros(2))
)
source_tangents = dphi * z_cross_source

# The balanced two-source modulation is a finite-support shear source: fixed
# masses give zero monopole derivative, the equal antipodal pair gives zero
# dipole derivative, and its symmetric trace-free quadrupole derivative is
# nonzero.  This is an exact statement about the declared discrete measure.
dipole_tangent = np.sum(source_masses[:, None] * source_tangents, axis=0)
quadrupole_tangent = np.zeros((3, 3))
for mass, position, tangent in zip(
    source_masses, source_positions, source_tangents
):
    quadrupole_tangent += mass * (
        3.0
        * (
            np.outer(tangent, position)
            + np.outer(position, tangent)
        )
        - 2.0 * np.dot(position, tangent) * np.eye(3)
    )
check(
    "balanced source has zero monopole/dipole derivative and nonzero STF quadrupole derivative",
    np.linalg.norm(dipole_tangent) < 1.0e-14
    and abs(np.trace(quadrupole_tangent)) < 1.0e-14
    and np.linalg.norm(quadrupole_tangent) > 1.0e-6,
    (
        f"dipole_norm={np.linalg.norm(dipole_tangent):.6e}, "
        f"quadrupole_trace={np.trace(quadrupole_tangent):.6e}, "
        f"quadrupole_norm={np.linalg.norm(quadrupole_tangent):.6e}"
    ),
)
kernel, source_column_test, gradient_test = geometry_coefficients(
    source_masses,
    source_positions,
    source_tangents,
    detector_masses,
    detector_positions,
)

epsilon = 1.0e-6
potential_plus = point_potential_per_G(
    source_masses,
    source_positions,
    detector_masses,
    detector_positions,
    theta=epsilon,
)
potential_minus = point_potential_per_G(
    source_masses,
    source_positions,
    detector_masses,
    detector_positions,
    theta=-epsilon,
)
torque_from_potential = -(potential_plus - potential_minus) / (2.0 * epsilon)
check(
    "point torque is minus the angle derivative of the exact potential",
    np.isclose(kernel, torque_from_potential, rtol=2e-9, atol=2e-10),
    f"analytic={kernel:.16e}, finite_difference={torque_from_potential:.16e}",
)

_, positions_plus = source_cloud(phi_test + epsilon * dphi)
_, positions_minus = source_cloud(phi_test - epsilon * dphi)
kernel_plus = geometry_coefficients(
    source_masses,
    positions_plus,
    dphi
    * np.column_stack((-positions_plus[:, 1], positions_plus[:, 0], np.zeros(2))),
    detector_masses,
    detector_positions,
)[0]
kernel_minus = geometry_coefficients(
    source_masses,
    positions_minus,
    dphi
    * np.column_stack((-positions_minus[:, 1], positions_minus[:, 0], np.zeros(2))),
    detector_masses,
    detector_positions,
)[0]
source_finite_difference = (kernel_plus - kernel_minus) / (2.0 * epsilon)
check(
    "source trajectory column is the analytic geometry derivative",
    np.isclose(source_column_test, source_finite_difference, rtol=2e-8, atol=2e-10),
    f"analytic={source_column_test:.16e}, finite_difference={source_finite_difference:.16e}",
)

gradient_plus = geometry_coefficients(
    source_masses,
    source_positions,
    source_tangents,
    detector_masses,
    detector_positions,
    theta=epsilon,
)[0]
gradient_minus = geometry_coefficients(
    source_masses,
    source_positions,
    source_tangents,
    detector_masses,
    detector_positions,
    theta=-epsilon,
)[0]
gradient_finite_difference = (gradient_plus - gradient_minus) / (2.0 * epsilon)
check(
    "gravitational stiffness is the analytic torsion-angle derivative",
    np.isclose(gradient_test, gradient_finite_difference, rtol=2e-8, atol=2e-9),
    f"analytic={gradient_test:.16e}, finite_difference={gradient_finite_difference:.16e}",
)

sphere_coefficients = uniform_sphere_coefficients(
    source_masses,
    source_positions,
    source_tangents,
    (0.018, 0.018),
    detector_masses,
    detector_positions,
    (0.009, 0.009),
)
check(
    "nonoverlapping uniform-sphere kernel equals the center-mass kernel",
    np.allclose(sphere_coefficients, (kernel, source_column_test, gradient_test), rtol=0, atol=0),
)


# ---------------------------------------------------------------------------
# Prospectively fixed synthetic finite apparatus.
frequencies_hz = np.array((0.0035, 0.0050, 0.0067, 0.0084, 0.0113, 0.0138, 0.0162, 0.0190))
omega = 2.0 * np.pi * frequencies_hz
phis = np.array((0.28, 0.39, 0.51, 0.64, 0.78, 0.91, 1.03, 1.16))
trajectory_amplitudes = np.array((0.035, -0.040, 0.032, -0.038, 0.041, -0.034, 0.037, -0.039))
source_columns = []
torsion_gradients = []
for phi, amplitude in zip(phis, trajectory_amplitudes):
    masses, positions = source_cloud(float(phi))
    tangents = amplitude * np.column_stack(
        (-positions[:, 1], positions[:, 0], np.zeros(2))
    )
    _, column, gradient = geometry_coefficients(
        masses,
        positions,
        tangents,
        detector_masses,
        detector_positions,
    )
    source_columns.append(column)
    torsion_gradients.append(gradient)
source_columns = np.asarray(source_columns, dtype=complex)
torsion_gradients = np.asarray(torsion_gradients)

calibration = MechanicalCalibration(
    torsion_inertia=5.0e-4,
    torsion_stiffness=2.0e-6,
    torsion_damping=6.5e-11,
    auxiliary_inertia=1.5e-4,
    auxiliary_stiffness=8.0e-6,
    auxiliary_damping=1.2e-7,
    coupling=4.0e-7,
    readout_gain=1.0003,
    readout_delay=0.23,
)

product_truth = SYNTHETIC_TRUTH_G * SYNTHETIC_SOURCE_SCALE
stability_at_bounds = np.vstack(
    tuple(
        static_stability_margin(
            product_bound, torsion_gradients, calibration
        )
        for product_bound in PRODUCT_SCAN_SI
    )
)
check(
    "declared positive-G scan remains inside the static stable quotient",
    bool(np.all(stability_at_bounds > 0.0)),
    f"minimum endpoint determinant={np.min(stability_at_bounds):.6e}",
)

remainder_theta = 2.5e-14 * np.exp(1j * np.linspace(0.2, 1.1, omega.size))
remainder_aux = 8.0e-15 * np.exp(-1j * np.linspace(0.1, 0.8, omega.size))
dressed = forward_complex(
    product_truth,
    omega,
    source_columns,
    torsion_gradients,
    calibration,
    remainder_torque=remainder_theta,
    auxiliary_remainder=remainder_aux,
)
full = full_coupled_complex(
    product_truth,
    omega,
    source_columns,
    torsion_gradients,
    calibration,
    remainder_torque=remainder_theta,
    auxiliary_remainder=remainder_aux,
)
check(
    "Schur-complement response equals the full coupled two-mode solve",
    np.allclose(dressed, full, rtol=2e-13, atol=2e-18),
    f"max_abs_difference={np.max(np.abs(dressed-full)):.6e}",
)

zero_response = forward_complex(
    product_truth,
    omega,
    np.zeros_like(source_columns),
    torsion_gradients,
    calibration,
)
check(
    "same-source same-remainder same-data response is exactly zero",
    np.array_equal(zero_response, np.zeros_like(zero_response)),
)

rescaling = 1.071
scale_left = forward_complex(
    SYNTHETIC_TRUTH_G * SYNTHETIC_SOURCE_SCALE,
    omega,
    source_columns,
    torsion_gradients,
    calibration,
)
scale_right = forward_complex(
    (SYNTHETIC_TRUTH_G * rescaling)
    * (SYNTHETIC_SOURCE_SCALE / rescaling),
    omega,
    source_columns,
    torsion_gradients,
    calibration,
)
check(
    "uncalibrated global source scale is exactly degenerate with G",
    np.array_equal(scale_left, scale_right),
)


# ---------------------------------------------------------------------------
# Correlated synthetic observation; truth is not supplied to the estimator.
signal = forward_complex(
    product_truth,
    omega,
    source_columns,
    torsion_gradients,
    calibration,
)
nuisance_complex = np.column_stack(
    (
        1.0e-8 * np.exp(0.31j) * np.ones(omega.size),
        1.0e-8j * omega / np.max(omega),
    )
)
nuisance_design = stack_complex_columns(nuisance_complex)
nuisance_truth = np.array((0.72, -0.43))

base_sigma = 2.0e-8
diagonal = (base_sigma * (1.0 + 0.2 * np.arange(2 * omega.size) / (2 * omega.size - 1))) ** 2
common_vector = np.tile((1.0, -0.35), omega.size)
observation_covariance = np.diag(diagonal) + (4.0e-9**2) * np.outer(
    common_vector, common_vector
)

# Independent gain and delay calibration uncertainty is propagated once.  Its
# Jacobian is frozen at the declared reference product above, not evaluated at
# the synthetic truth and not updated during the G scan.
signal_real = stack_complex(signal)
covariance_reference_signal = forward_complex(
    COVARIANCE_LINEARIZATION_PRODUCT,
    omega,
    source_columns,
    torsion_gradients,
    calibration,
)
gain_jacobian = (
    stack_complex(covariance_reference_signal) / calibration.readout_gain
)
delay_jacobian = stack_complex(-1j * omega * covariance_reference_signal)
calibration_jacobian = np.column_stack((gain_jacobian, delay_jacobian))
calibration_covariance = np.diag((2.0e-4**2, 1.5e-3**2))
total_covariance = effective_covariance(
    observation_covariance, calibration_jacobian, calibration_covariance
)
check(
    "observation and independent calibration covariance are combined once",
    np.min(np.linalg.eigvalsh(total_covariance)) > 0.0
    and np.all(np.diag(total_covariance) >= np.diag(observation_covariance)),
)

rng = np.random.default_rng(20260827)
noise = rng.multivariate_normal(np.zeros(2 * omega.size), total_covariance)
observed = signal_real + nuisance_design @ nuisance_truth + noise

product_grid = np.linspace(
    PRODUCT_SCAN_SI[0], PRODUCT_SCAN_SI[1], PRODUCT_SCAN_POINTS
)


def prediction_real(product: float) -> np.ndarray:
    return stack_complex(
        forward_complex(
            product,
            omega,
            source_columns,
            torsion_gradients,
            calibration,
        )
    )


grid, scores, nuisance_profile = profile_product_grid(
    observed,
    total_covariance,
    product_grid,
    prediction_real,
    nuisance_design=nuisance_design,
)
p_low, p_high, p_best, best_index = connected_grid_interval(
    grid, scores, delta=3.841458820694124
)
check(
    "blinded synthetic product lies in the profiled covariance interval",
    p_low <= product_truth <= p_high,
    f"interval=[{p_low:.16e},{p_high:.16e}], truth={product_truth:.16e}",
)

source_scale_interval = SOURCE_SCALE_CALIBRATION_INTERVAL
g_low, g_high = source_scale_identified_interval(
    (p_low, p_high), source_scale_interval
)
check(
    "independently calibrated source-scale interval propagates to G",
    g_low <= SYNTHETIC_TRUTH_G <= g_high,
    f"G_interval=[{g_low:.16e},{g_high:.16e}]",
)

relative_product_error = abs(p_best - product_truth) / product_truth
check(
    "blinded synthetic product point estimate is accurate",
    relative_product_error < 3.0e-3,
    f"relative_error={relative_product_error:.6e}",
)


# Held-out frequency/configuration prediction with no refit of p.
train_rows = np.arange(0, 6)
hold_rows = np.arange(6, 8)


def real_indices(rows):
    return np.ravel(np.column_stack((2 * rows, 2 * rows + 1)))


train_index = real_indices(train_rows)
hold_index = real_indices(hold_rows)


def train_prediction(product: float) -> np.ndarray:
    return prediction_real(product)[train_index]


train_grid, train_scores, train_nuisance = profile_product_grid(
    observed[train_index],
    total_covariance[np.ix_(train_index, train_index)],
    product_grid,
    train_prediction,
    nuisance_design=nuisance_design[train_index],
)
train_best_index = int(np.argmin(train_scores))
train_best_product = float(train_grid[train_best_index])
train_eta = train_nuisance[train_best_index]
hold_prediction = (
    prediction_real(train_best_product)[hold_index]
    + nuisance_design[hold_index] @ train_eta
)
hold_score, _ = profiled_quadratic(
    observed[hold_index],
    hold_prediction,
    total_covariance[np.ix_(hold_index, hold_index)],
)
check(
    "held-out configurations agree without refitting G or nuisance",
    hold_score < 13.276704135987622,
    f"four-real-dof Mahalanobis={hold_score:.6f}",
)


# A deliberately incomplete bare transfer must not reproduce the dressed fit.
bare_calibration = MechanicalCalibration(
    torsion_inertia=calibration.torsion_inertia,
    torsion_stiffness=calibration.torsion_stiffness,
    torsion_damping=calibration.torsion_damping,
    auxiliary_inertia=calibration.auxiliary_inertia,
    auxiliary_stiffness=calibration.auxiliary_stiffness,
    auxiliary_damping=calibration.auxiliary_damping,
    coupling=1.0e-30,
    readout_gain=calibration.readout_gain,
    readout_delay=calibration.readout_delay,
)


def bare_prediction(product: float) -> np.ndarray:
    return stack_complex(
        forward_complex(
            product,
            omega,
            source_columns,
            np.zeros_like(torsion_gradients),
            bare_calibration,
        )
    )


_, bare_scores, _ = profile_product_grid(
    observed,
    total_covariance,
    product_grid,
    bare_prediction,
    nuisance_design=nuisance_design,
)
bare_best = float(product_grid[int(np.argmin(bare_scores))])
bare_shift = abs(bare_best - p_best) / p_best
check(
    "omitting the dressed operator produces a detectable calibration shift",
    bare_shift > 2.0e-3,
    f"relative_shift={bare_shift:.6e}",
)


passed = sum(item[1] for item in checks)
failed = len(checks) - passed
result = {
    "schema": "WAC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001",
    "claim_class": "SYNTHETIC_FORWARD_MODEL_VALIDATION_ONLY",
    "accepted_G_used_as_input": False,
    "lineage_charge_inferred": False,
    "synthetic_truth_G_SI": SYNTHETIC_TRUTH_G,
    "synthetic_source_scale": SYNTHETIC_SOURCE_SCALE,
    "covariance_linearization_product": COVARIANCE_LINEARIZATION_PRODUCT,
    "prospective_product_scan_SI": list(PRODUCT_SCAN_SI),
    "prospective_product_scan_points": PRODUCT_SCAN_POINTS,
    "synthetic_product_truth": product_truth,
    "profiled_product_best": p_best,
    "profiled_product_interval_95_grid": [p_low, p_high],
    "source_scale_interval": list(source_scale_interval),
    "propagated_G_interval_95_grid_SI": [g_low, g_high],
    "profiled_nuisance_at_best": nuisance_profile[best_index].tolist(),
    "held_out_mahalanobis_four_real_dof": hold_score,
    "bare_operator_relative_product_shift": bare_shift,
    "checks_passed": passed,
    "checks_total": len(checks),
    "checks": [
        {"name": name, "pass": condition, "detail": detail}
        for name, condition, detail in checks
    ],
    "scientific_ceiling": (
        "NO_MEASUREMENT_OF_G__NO_LINEAGE_SOURCE__NO_RGRL_EMPIRICAL_PASS__"
        "FINITE_APPARATUS_FORWARD_AND_IDENTIFIABILITY_MODEL_ONLY"
    ),
}
print(json.dumps(result, indent=2, sort_keys=True))
raise SystemExit(0 if failed == 0 else 1)
