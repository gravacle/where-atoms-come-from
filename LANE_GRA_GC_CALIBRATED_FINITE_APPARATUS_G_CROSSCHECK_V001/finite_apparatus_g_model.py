"""Finite-apparatus source-to-torsion forward model for a G cross-check.

The model is deliberately independent of any record/lineage charge.  Its
source is a nongravitationally calibrated compact mass measure and trajectory.
It retains the gravitational torque gradient and one calibrated auxiliary
mechanical mode in the dressed operator, and keeps physical remainders,
homogeneous data, and readout terms in distinct columns.

Fourier convention: real parts of amplitudes times exp(+i omega t).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np


class ForwardModelError(ValueError):
    """The declared finite-apparatus calibration is ill typed or singular."""


def _vector(name: str, value, *, complex_ok: bool = False) -> np.ndarray:
    dtype = complex if complex_ok else float
    result = np.asarray(value, dtype=dtype)
    if result.ndim != 1 or result.size == 0:
        raise ForwardModelError(f"{name} must be a nonempty vector")
    if not np.all(np.isfinite(result)):
        raise ForwardModelError(f"{name} contains a nonfinite value")
    return result


def _points(name: str, value) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if result.ndim != 2 or result.shape[1] != 3 or result.shape[0] == 0:
        raise ForwardModelError(f"{name} must have shape (n,3)")
    if not np.all(np.isfinite(result)):
        raise ForwardModelError(f"{name} contains a nonfinite value")
    return result


def rotation_z(angle: float) -> np.ndarray:
    """Right-handed rotation about the torsion axis."""
    c = float(np.cos(angle))
    s = float(np.sin(angle))
    return np.array(((c, -s, 0.0), (s, c, 0.0), (0.0, 0.0, 1.0)))


def rotate_detector(points, theta: float) -> np.ndarray:
    return _points("detector positions", points) @ rotation_z(theta).T


def point_potential_per_G(
    source_masses,
    source_positions,
    detector_masses,
    detector_positions,
    *,
    theta: float = 0.0,
) -> float:
    """Exact Newtonian U/G for two finite discrete positive mass measures."""
    sm = _vector("source masses", source_masses)
    sp = _points("source positions", source_positions)
    dm = _vector("detector masses", detector_masses)
    dp = rotate_detector(detector_positions, theta)
    if sm.size != sp.shape[0] or dm.size != dp.shape[0]:
        raise ForwardModelError("mass and position counts disagree")
    if np.any(sm <= 0.0) or np.any(dm <= 0.0):
        raise ForwardModelError("all calibrated mass elements must be positive")

    total = 0.0
    for mass_s, pos_s in zip(sm, sp):
        for mass_d, pos_d in zip(dm, dp):
            distance = float(np.linalg.norm(pos_s - pos_d))
            if distance <= 0.0:
                raise ForwardModelError("source and detector supports overlap")
            total -= mass_s * mass_d / distance
    return total


def geometry_coefficients(
    source_masses,
    source_positions,
    source_tangents,
    detector_masses,
    detector_positions,
    *,
    theta: float = 0.0,
) -> tuple[float, float, float]:
    """Return (K, a, k_g) for a finite calibrated mass measure.

    K is torque/G, a=d(K)/du for the declared source-trajectory tangent,
    and k_g=d(K)/d(theta).  These are analytic derivatives of the exact
    discrete Newton kernel, not fitted gravity amplitudes.
    """
    sm = _vector("source masses", source_masses)
    sp = _points("source positions", source_positions)
    sv = _points("source tangents", source_tangents)
    dm = _vector("detector masses", detector_masses)
    dp = rotate_detector(detector_positions, theta)
    if sm.size != sp.shape[0] or sp.shape != sv.shape or dm.size != dp.shape[0]:
        raise ForwardModelError("mass, position, or trajectory counts disagree")
    if np.any(sm <= 0.0) or np.any(dm <= 0.0):
        raise ForwardModelError("all calibrated mass elements must be positive")

    z_axis = np.array((0.0, 0.0, 1.0))
    torque = 0.0
    source_column = 0.0
    torsion_gradient = 0.0
    for mass_s, pos_s, vel_s in zip(sm, sp, sv):
        for mass_d, pos_d in zip(dm, dp):
            displacement = pos_s - pos_d
            radius = float(np.linalg.norm(displacement))
            if radius <= 0.0:
                raise ForwardModelError("source and detector supports overlap")
            weight = mass_s * mass_d
            cross_sr = float(np.cross(pos_d, pos_s)[2])
            inv3 = radius ** -3
            inv5 = radius ** -5

            torque += weight * cross_sr * inv3
            source_column += weight * (
                float(np.cross(pos_d, vel_s)[2]) * inv3
                - 3.0 * cross_sr * float(displacement @ vel_s) * inv5
            )

            detector_tangent = np.cross(z_axis, pos_d)
            torsion_gradient += weight * (
                float(np.cross(detector_tangent, pos_s)[2]) * inv3
                + 3.0
                * cross_sr
                * float(displacement @ detector_tangent)
                * inv5
            )
    return float(torque), float(source_column), float(torsion_gradient)


def uniform_sphere_coefficients(
    source_masses,
    source_centers,
    source_tangents,
    source_radii,
    detector_masses,
    detector_centers,
    detector_radii,
    *,
    theta: float = 0.0,
) -> tuple[float, float, float]:
    """Exact external coefficients for nonoverlapping uniform spheres.

    Newton's shell theorem reduces every nonoverlapping spherical pair to its
    center mass.  The explicit separation check is part of the model domain.
    """
    sp = _points("source centers", source_centers)
    dp = rotate_detector(detector_centers, theta)
    sr = _vector("source radii", source_radii)
    dr = _vector("detector radii", detector_radii)
    if sr.size != sp.shape[0] or dr.size != dp.shape[0]:
        raise ForwardModelError("sphere radius and center counts disagree")
    if np.any(sr < 0.0) or np.any(dr < 0.0):
        raise ForwardModelError("sphere radii must be nonnegative")
    for i, pos_s in enumerate(sp):
        for j, pos_d in enumerate(dp):
            if np.linalg.norm(pos_s - pos_d) <= sr[i] + dr[j]:
                raise ForwardModelError("uniform spheres are not disjoint")
    return geometry_coefficients(
        source_masses,
        sp,
        source_tangents,
        detector_masses,
        detector_centers,
        theta=theta,
    )


@dataclass(frozen=True)
class MechanicalCalibration:
    """Nongravitational calibration of two normalized dimensionless modes.

    The auxiliary coordinate is normalized so the reciprocal off-diagonal
    coupling has the same generalized-stiffness convention in both rows.
    """

    torsion_inertia: float
    torsion_stiffness: float
    torsion_damping: float
    auxiliary_inertia: float
    auxiliary_stiffness: float
    auxiliary_damping: float
    coupling: float
    readout_gain: float = 1.0
    readout_delay: float = 0.0

    def __post_init__(self) -> None:
        positive = (
            self.torsion_inertia,
            self.torsion_stiffness,
            self.torsion_damping,
            self.auxiliary_inertia,
            self.auxiliary_stiffness,
            self.auxiliary_damping,
            self.readout_gain,
        )
        if not all(np.isfinite(item) and item > 0.0 for item in positive):
            raise ForwardModelError("mechanical calibration must be finite and positive")
        if not np.isfinite(self.coupling) or not np.isfinite(self.readout_delay):
            raise ForwardModelError("coupling and delay must be finite")


def bare_denominators(omega, calibration: MechanicalCalibration):
    omega = _vector("angular frequencies", omega)
    d_theta = (
        calibration.torsion_stiffness
        - calibration.torsion_inertia * omega**2
        + 1j * calibration.torsion_damping * omega
    )
    d_aux = (
        calibration.auxiliary_stiffness
        - calibration.auxiliary_inertia * omega**2
        + 1j * calibration.auxiliary_damping * omega
    )
    if np.any(np.abs(d_aux) == 0.0):
        raise ForwardModelError("auxiliary retarded denominator is singular")
    return d_theta, d_aux


def dressed_denominator(
    product: float,
    omega,
    torsion_gradient,
    calibration: MechanicalCalibration,
) -> np.ndarray:
    """Exact two-mode Schur complement after the metric is integrated out.

    product is p=G*s, where s is a global source-mass scale calibration.
    """
    gradient = _vector("torsion gradient", torsion_gradient)
    d_theta, d_aux = bare_denominators(omega, calibration)
    if gradient.size != d_theta.size:
        raise ForwardModelError("gradient and frequency counts disagree")
    return (
        d_theta
        - float(product) * gradient
        - calibration.coupling**2 / d_aux
    )


def readout_transfer(omega, calibration: MechanicalCalibration) -> np.ndarray:
    omega = _vector("angular frequencies", omega)
    return calibration.readout_gain * np.exp(-1j * omega * calibration.readout_delay)


def forward_complex(
    product: float,
    omega,
    source_column,
    torsion_gradient,
    calibration: MechanicalCalibration,
    *,
    remainder_torque=None,
    auxiliary_remainder=None,
    homogeneous_angle=None,
    readout_bias=None,
) -> np.ndarray:
    """Complete reduced forward map with every derivative owned once."""
    if not np.isfinite(product) or product < 0.0:
        raise ForwardModelError("G times source scale must be finite and nonnegative")
    omega = _vector("angular frequencies", omega)
    source = _vector("source column", source_column, complex_ok=True)
    gradient = _vector("torsion gradient", torsion_gradient)
    if source.size != omega.size or gradient.size != omega.size:
        raise ForwardModelError("source, gradient, and frequency counts disagree")

    zeros = np.zeros(omega.size, dtype=complex)
    r_theta = zeros if remainder_torque is None else _vector(
        "remainder torque", remainder_torque, complex_ok=True
    )
    r_aux = zeros if auxiliary_remainder is None else _vector(
        "auxiliary remainder", auxiliary_remainder, complex_ok=True
    )
    data = zeros if homogeneous_angle is None else _vector(
        "homogeneous angle", homogeneous_angle, complex_ok=True
    )
    bias = zeros if readout_bias is None else _vector(
        "readout bias", readout_bias, complex_ok=True
    )
    if any(item.size != omega.size for item in (r_theta, r_aux, data, bias)):
        raise ForwardModelError("all response columns must match the frequency count")

    denominator = dressed_denominator(
        product, omega, gradient, calibration
    )
    if np.any(np.abs(denominator) == 0.0):
        raise ForwardModelError("dressed retarded denominator is singular")
    _, d_aux = bare_denominators(omega, calibration)
    dressed_remainder = r_theta + calibration.coupling * r_aux / d_aux
    angle = (product * source + dressed_remainder) / denominator + data
    return readout_transfer(omega, calibration) * angle + bias


def full_coupled_complex(
    product: float,
    omega,
    source_column,
    torsion_gradient,
    calibration: MechanicalCalibration,
    *,
    remainder_torque=None,
    auxiliary_remainder=None,
) -> np.ndarray:
    """Solve the full two-mode system before the Schur complement."""
    omega = _vector("angular frequencies", omega)
    source = _vector("source column", source_column, complex_ok=True)
    gradient = _vector("torsion gradient", torsion_gradient)
    zeros = np.zeros(omega.size, dtype=complex)
    r_theta = zeros if remainder_torque is None else _vector(
        "remainder torque", remainder_torque, complex_ok=True
    )
    r_aux = zeros if auxiliary_remainder is None else _vector(
        "auxiliary remainder", auxiliary_remainder, complex_ok=True
    )
    d_theta, d_aux = bare_denominators(omega, calibration)
    result = np.empty(omega.size, dtype=complex)
    for index in range(omega.size):
        operator = np.array(
            (
                (d_theta[index] - product * gradient[index], -calibration.coupling),
                (-calibration.coupling, d_aux[index]),
            ),
            dtype=complex,
        )
        column = np.array(
            (product * source[index] + r_theta[index], r_aux[index]),
            dtype=complex,
        )
        result[index] = np.linalg.solve(operator, column)[0]
    return readout_transfer(omega, calibration) * result


def static_stability_margin(
    product: float,
    torsion_gradient,
    calibration: MechanicalCalibration,
) -> np.ndarray:
    """Static positive-definiteness determinant for each geometry row."""
    gradient = _vector("torsion gradient", torsion_gradient)
    return (
        (calibration.torsion_stiffness - product * gradient)
        * calibration.auxiliary_stiffness
        - calibration.coupling**2
    )


def stack_complex(values) -> np.ndarray:
    values = _vector("complex response", values, complex_ok=True)
    return np.column_stack((values.real, values.imag)).reshape(-1)


def stack_complex_columns(columns) -> np.ndarray:
    columns = np.asarray(columns, dtype=complex)
    if columns.ndim != 2 or columns.shape[0] == 0:
        raise ForwardModelError("complex nuisance columns must have shape (n,k)")
    return np.stack((columns.real, columns.imag), axis=1).reshape(
        2 * columns.shape[0], columns.shape[1]
    )


def effective_covariance(observation_covariance, jacobian, calibration_covariance):
    """First-order independent calibration-covariance propagation."""
    sigma_y = np.asarray(observation_covariance, dtype=float)
    jacobian = np.asarray(jacobian, dtype=float)
    sigma_nu = np.asarray(calibration_covariance, dtype=float)
    if sigma_y.ndim != 2 or sigma_y.shape[0] != sigma_y.shape[1]:
        raise ForwardModelError("observation covariance must be square")
    if jacobian.shape[0] != sigma_y.shape[0]:
        raise ForwardModelError("calibration Jacobian has the wrong row count")
    if sigma_nu.shape != (jacobian.shape[1], jacobian.shape[1]):
        raise ForwardModelError("calibration covariance has the wrong shape")
    result = sigma_y + jacobian @ sigma_nu @ jacobian.T
    if np.min(np.linalg.eigvalsh(result)) <= 0.0:
        raise ForwardModelError("effective covariance is not positive definite")
    return result


def profiled_quadratic(
    observed,
    predicted,
    covariance,
    nuisance_design=None,
    nuisance_precision=None,
) -> tuple[float, np.ndarray]:
    """Profile linear readout nuisances with optional Gaussian precision."""
    observed = _vector("observed real response", observed)
    predicted = _vector("predicted real response", predicted)
    covariance = np.asarray(covariance, dtype=float)
    if observed.shape != predicted.shape or covariance.shape != (
        observed.size,
        observed.size,
    ):
        raise ForwardModelError("observation, prediction, and covariance disagree")
    weight = np.linalg.inv(covariance)
    residual = observed - predicted
    if nuisance_design is None:
        return float(residual @ weight @ residual), np.zeros(0)

    design = np.asarray(nuisance_design, dtype=float)
    if design.ndim != 2 or design.shape[0] != observed.size:
        raise ForwardModelError("nuisance design has the wrong shape")
    precision = (
        np.zeros((design.shape[1], design.shape[1]))
        if nuisance_precision is None
        else np.asarray(nuisance_precision, dtype=float)
    )
    if precision.shape != (design.shape[1], design.shape[1]):
        raise ForwardModelError("nuisance precision has the wrong shape")
    normal = design.T @ weight @ design + precision
    rhs = design.T @ weight @ residual
    nuisance = np.linalg.pinv(normal, rcond=1e-14) @ rhs
    final = residual - design @ nuisance
    score = float(final @ weight @ final + nuisance @ precision @ nuisance)
    return score, nuisance


def profile_product_grid(
    observed,
    covariance,
    product_grid,
    prediction: Callable[[float], np.ndarray],
    *,
    nuisance_design=None,
    nuisance_precision=None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Deterministically profile a positive p=G*s grid."""
    grid = _vector("product grid", product_grid)
    if np.any(np.diff(grid) <= 0.0) or np.any(grid <= 0.0):
        raise ForwardModelError("product grid must be strictly increasing and positive")
    scores = np.empty(grid.size)
    nuisance_rows = []
    for index, product in enumerate(grid):
        score, nuisance = profiled_quadratic(
            observed,
            prediction(float(product)),
            covariance,
            nuisance_design=nuisance_design,
            nuisance_precision=nuisance_precision,
        )
        scores[index] = score
        nuisance_rows.append(nuisance)
    return grid, scores, np.asarray(nuisance_rows)


def connected_grid_interval(
    grid,
    scores,
    *,
    delta: float,
) -> tuple[float, float, float, int]:
    """Return the connected grid identified interval around the best fit."""
    grid = _vector("grid", grid)
    scores = _vector("scores", scores)
    if grid.size != scores.size or delta < 0.0:
        raise ForwardModelError("grid, scores, or threshold are invalid")
    best = int(np.argmin(scores))
    accepted = scores <= scores[best] + delta
    left = best
    right = best
    while left > 0 and accepted[left - 1]:
        left -= 1
    while right + 1 < grid.size and accepted[right + 1]:
        right += 1
    if np.any(accepted[:left]) or np.any(accepted[right + 1 :]):
        raise ForwardModelError("identified grid set is disconnected")
    if left == 0 or right == grid.size - 1:
        raise ForwardModelError("identified interval reaches the scan boundary")
    return float(grid[left]), float(grid[right]), float(grid[best]), best


def source_scale_identified_interval(
    product_interval: tuple[float, float],
    source_scale_interval: tuple[float, float],
) -> tuple[float, float]:
    """Exact positive quotient interval {p/s}."""
    p_low, p_high = map(float, product_interval)
    s_low, s_high = map(float, source_scale_interval)
    if not (0.0 < p_low <= p_high and 0.0 < s_low <= s_high):
        raise ForwardModelError("product and source-scale intervals must be positive")
    return p_low / s_high, p_high / s_low
