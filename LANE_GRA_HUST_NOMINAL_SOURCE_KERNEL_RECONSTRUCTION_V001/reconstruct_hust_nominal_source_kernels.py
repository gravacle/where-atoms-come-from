#!/usr/bin/env python3
"""Reconstruct HUST-2018 nominal source kernels from public geometry.

The calculation consumes only transcribed public dimensions, masses, centre
separations, campaign temperatures, and independent mechanical inertias.  It
first constructs the conditional homogeneous AAF forcing coefficient and the
conditional homogeneous ToS stiffness coefficients.  Only after those values
and their public-input sensitivities are frozen does it open the separate file
of published processed coefficients as a post-calculation comparator.

No value of G, measured response, processed source coefficient, or fitted
gravity result enters either kernel.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIELDS = HERE / "SOURCE_FIELDS.json"
COMPARATORS = HERE / "PUBLISHED_COMPARATORS.json"
RESULT = HERE / "RESULT.json"

DEPENDENCIES = {
    "SOURCE/41586_2018_431_MOESM1_ESM.pdf": (
        HERE / "SOURCE/41586_2018_431_MOESM1_ESM.pdf",
        "5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb",
    ),
    "SOURCE/41586_2018_431_Tab1_ESM.jpg": (
        HERE / "SOURCE/41586_2018_431_Tab1_ESM.jpg",
        "c84bc71bcd0115ccbbbdddd70dd1d755ffdabd37c91772f30cb51a264a200195",
    ),
    "SOURCE/41586_2018_431_Tab2_ESM.jpg": (
        HERE / "SOURCE/41586_2018_431_Tab2_ESM.jpg",
        "96c40827af03f4de0715ea77bc69c0612d3ef94dc15eff035213e8a4dc0649c1",
    ),
    "SOURCE/41586_2018_431_Tab4_ESM.jpg": (
        HERE / "SOURCE/41586_2018_431_Tab4_ESM.jpg",
        "567e5c8b953cba86e642bfd01b3880b262873606e2c69af44fa13a9ff4f629ce",
    ),
    "SOURCE/41586_2018_431_Fig6_ESM.jpg": (
        HERE / "SOURCE/41586_2018_431_Fig6_ESM.jpg",
        "f29cdc1909149fc3a03264299889ef7f8346eb0fae6d7415cebbcf3d14163312",
    ),
    "SOURCE/nature_main_table1_error_budget.html": (
        HERE / "SOURCE/nature_main_table1_error_budget.html",
        "23436d4be7600a7a9dffa02cc4167a20b6eea032a181e77899bb57bb90aa02e9",
    ),
    "GC16 theorem": (
        ROOT / "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/THEOREM.md",
        "cbf0733633ba93756b08dded7486a9be76beb572807693455c761bd36a8f0f5b",
    ),
    "GC16 protocol": (
        ROOT / "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/PROTOCOL.md",
        "6ec7d8f0ce9a184d25612107dbfc294dd22d124ebe859831401e9cc0c8e8b819",
    ),
    "GC16 executable": (
        ROOT / "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/finite_apparatus_g_model.py",
        "6c17498d2d65f6420498ac559a97a2c3bbf49e110dd971da34b4c9c9bea2e4e4",
    ),
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def verify_dependencies() -> dict[str, str]:
    observed = {name: digest(path) for name, (path, _) in DEPENDENCIES.items()}
    expected = {name: expected for name, (_, expected) in DEPENDENCIES.items()}
    if observed != expected:
        raise RuntimeError(f"source/dependency custody failure: {observed}")
    return observed


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def cuboid_rule(params: dict[str, float], order: int) -> tuple[np.ndarray, ...]:
    node, weight = leggauss(order)
    x_grid, y_grid, z_grid = np.meshgrid(
        node * params["L"] / 2.0,
        node * params["W"] / 2.0,
        node * params["H"] / 2.0,
        indexing="ij",
    )
    # The weights integrate a uniform body average: sum(weights)=1.
    weights = np.einsum("i,j,k->ijk", weight, weight, weight).ravel() / 8.0
    return x_grid.ravel(), y_grid.ravel(), z_grid.ravel(), weights


def core_inertia(params: dict[str, float]) -> float:
    return params["M"] * (params["L"] ** 2 + params["W"] ** 2) / 12.0


def aaf_centres(params: dict[str, float]) -> np.ndarray:
    """Conditional pairwise-centred realization of the four reported distances.

    The Extended Data table calls S7,9 and S10,12 horizontal coordinate
    distances and S7,10 and S9,12 vertical coordinate distances.  The small
    left/right mismatches are retained rather than projected onto a perfect
    rectangle.  Pairwise centring is an explicit premise: the four separations
    and one overall centring condition do not publish the two remaining shear
    coordinates of the individual sphere centres.
    """

    return np.asarray(
        (
            (-params["S79"] / 2.0, 0.0, +params["S710"] / 2.0),
            (+params["S79"] / 2.0, 0.0, +params["S912"] / 2.0),
            (-params["S1012"] / 2.0, 0.0, -params["S710"] / 2.0),
            (+params["S1012"] / 2.0, 0.0, -params["S912"] / 2.0),
        ),
        dtype=float,
    )


def aaf_harmonic(params: dict[str, float], *, order: int = 16, angles: int = 512) -> float:
    """Return the m=2 angular-acceleration coefficient per G, in kg/m^3."""

    x, y, z, weights = cuboid_rule(params, order)
    inertia = core_inertia(params)
    centres = aaf_centres(params)
    source_masses = np.asarray((params["m7"], params["m9"], params["m10"], params["m12"]))
    phase = np.arange(angles, dtype=float) * (2.0 * math.pi / angles)
    acceleration_per_g = np.empty(angles, dtype=float)
    for index, phi in enumerate(phase):
        c, s = math.cos(phi), math.sin(phi)
        rotated = centres.copy()
        rotated[:, 0] = c * centres[:, 0] - s * centres[:, 1]
        rotated[:, 1] = s * centres[:, 0] + c * centres[:, 1]
        torque_per_g = 0.0
        for mass, (rx, ry, rz) in zip(source_masses, rotated):
            distance2 = (rx - x) ** 2 + (ry - y) ** 2 + (rz - z) ** 2
            numerator = x * ry - y * rx
            torque_per_g += mass * params["M"] * float(
                np.sum(weights * numerator / distance2**1.5)
            )
        acceleration_per_g[index] = torque_per_g / inertia
    coefficient = 2.0 * abs(np.mean(acceleration_per_g * np.exp(-2j * phase)))
    return float(coefficient)


def tos_stiffness(params: dict[str, float], *, order: int = 20) -> tuple[float, float, float]:
    """Return (near-minus-far stiffness, near, far), each divided by G I.

    For source azimuth phi, N=x R_y-y R_x and P=x R_x+y R_y.  The derivative
    of the torque kernel is P/d^3-3N^2/d^5.  Its value at phi=0 is C_g,n/I;
    its value at phi=pi/2 is C_g,f/I in the paper's sign convention.
    """

    x, y, z, weights = cuboid_rule(params, order)
    inertia = core_inertia(params)
    values = []
    for phi in (0.0, math.pi / 2.0):
        c, s = math.cos(phi), math.sin(phi)
        derivative = 0.0
        for mass, sign in ((params["m1"], -1.0), (params["m2"], +1.0)):
            rx = sign * params["S"] * c / 2.0
            ry = sign * params["S"] * s / 2.0
            distance2 = (rx - x) ** 2 + (ry - y) ** 2 + z**2
            n_cross = x * ry - y * rx
            p_dot = x * rx + y * ry
            integrand = p_dot / distance2**1.5 - 3.0 * n_cross**2 / distance2**2.5
            derivative += mass * params["M"] * float(np.sum(weights * integrand)) / inertia
        values.append(derivative)
    return float(values[0] - values[1]), float(values[0]), float(values[1])


def exact_minimum_clearance_over_azimuth(
    params: dict[str, float], centres: np.ndarray, radii: np.ndarray
) -> float:
    """Exact minimum sphere-surface/cuboid clearance over source azimuth.

    Every declared centre starts in the x-z plane and rotates about z.  For a
    fixed orbital radius, the closest point of the rotating centre to the
    rectangular x-y cross-section is obtained by aligning it with the farthest
    cuboid corner.  The vertical separation is independent.  This removes the
    earlier sampled-clearance qualifier from the shell-theorem domain check.
    """

    if not np.allclose(centres[:, 1], 0.0, atol=0.0, rtol=0.0):
        raise ValueError("analytic clearance requires declared centres in the x-z plane")
    half_xy_radius = math.hypot(params["L"] / 2.0, params["W"] / 2.0)
    half_height = params["H"] / 2.0
    clearances = []
    for centre, radius in zip(centres, radii):
        orbit_radius = math.hypot(float(centre[0]), float(centre[1]))
        horizontal_gap = max(orbit_radius - half_xy_radius, 0.0)
        vertical_gap = max(abs(float(centre[2])) - half_height, 0.0)
        clearances.append(math.hypot(horizontal_gap, vertical_gap) - float(radius))
    return min(clearances)


def central_sensitivities(function, params: dict[str, float], steps: dict[str, float]) -> dict:
    components = {}
    nominal = float(function(params))
    for name, step in steps.items():
        plus = dict(params)
        minus = dict(params)
        plus[name] += step
        minus[name] -= step
        half_change = abs(float(function(plus)) - float(function(minus))) / 2.0
        components[name] = {
            "input_standard_uncertainty": step,
            "coefficient_standard_sensitivity_kg_m-3": half_change,
            "relative_ppm": half_change / nominal * 1e6,
        }
    rss = math.sqrt(sum(item["coefficient_standard_sensitivity_kg_m-3"] ** 2 for item in components.values()))
    l1 = sum(item["coefficient_standard_sensitivity_kg_m-3"] for item in components.values())
    return {
        "method": "symmetric one-standard-uncertainty finite perturbations",
        "components": components,
        "rss_standard_uncertainty_kg_m-3": rss,
        "rss_relative_ppm": rss / nominal * 1e6,
        "linearized_axis_box_half_width_kg_m-3": l1,
        "linearized_axis_box_relative_ppm": l1 / nominal * 1e6,
        "coverage_ceiling": (
            "RSS is first-order standard-uncertainty propagation. The L1 number is a local "
            "linearized axis-box diagnostic. Neither is an exact confidence or coverage theorem."
        ),
    }


def group_rss(sensitivities: dict, names: tuple[str, ...]) -> float:
    components = sensitivities["components"]
    return math.sqrt(sum(components[name]["relative_ppm"] ** 2 for name in names))


def aaf_parameters(fields: dict, campaign: dict) -> tuple[dict[str, float], dict[str, float]]:
    pendulum = fields["pendulums"]["AAF"]
    spheres = {item["id"]: item for item in fields["source_spheres"]["AAF"]}
    distance = fields["aaf_distances_23p7C"]
    temperature_offset = campaign["average_temperature_C"] - 23.7
    s79 = (
        distance["S7_9_m"]
        + distance["upper_horizontal_temperature_coefficient_m_per_C"] * temperature_offset
    )
    u_s79 = math.hypot(
        distance["u_S7_9_m"],
        temperature_offset * distance["u_upper_horizontal_temperature_coefficient_m_per_C"],
    )
    params = {
        "L": pendulum["length_m"],
        "W": pendulum["width_m"],
        "H": pendulum["height_m"],
        "M": pendulum["mass_kg"],
        "m7": spheres["7"]["mass_kg"],
        "m9": spheres["9"]["mass_kg"],
        "m10": spheres["10"]["mass_kg"],
        "m12": spheres["12"]["mass_kg"],
        "S79": s79,
        "S1012": distance["S10_12_m"],
        "S710": distance["S7_10_m"],
        "S912": distance["S9_12_m"],
    }
    steps = {
        "L": pendulum["u_length_m"],
        "W": pendulum["u_width_m"],
        "H": pendulum["u_height_m"],
        "M": pendulum["u_mass_kg"],
        "m7": spheres["7"]["u_mass_kg"],
        "m9": spheres["9"]["u_mass_kg"],
        "m10": spheres["10"]["u_mass_kg"],
        "m12": spheres["12"]["u_mass_kg"],
        "S79": u_s79,
        "S1012": distance["u_S10_12_m"],
        "S710": distance["u_S7_10_m"],
        "S912": distance["u_S9_12_m"],
    }
    return params, steps


def reconstruct_aaf(fields: dict) -> list[dict]:
    output = []
    pendulum = fields["pendulums"]["AAF"]
    sphere_radii = np.asarray(
        [item["diameter_m"] / 2.0 for item in fields["source_spheres"]["AAF"]]
    )
    for campaign in fields["aaf_campaigns"]:
        params, steps = aaf_parameters(fields, campaign)
        coefficient = aaf_harmonic(params, order=16, angles=512)
        coarse = aaf_harmonic(params, order=12, angles=256)
        sensitivity = central_sensitivities(
            lambda trial: aaf_harmonic(trial, order=12, angles=256), params, steps
        )
        groups = {
            "pendulum_dimensions_ppm": group_rss(sensitivity, ("L", "W", "H")),
            "source_masses_ppm": group_rss(sensitivity, ("m7", "m9", "m10", "m12")),
            "horizontal_distances_ppm": group_rss(sensitivity, ("S79", "S1012")),
            "vertical_distances_ppm": group_rss(sensitivity, ("S710", "S912")),
        }
        source_centres = aaf_centres(params)
        core_i = core_inertia(params)
        full_i = pendulum["supplement_full_apparatus_I_kg_m2"]
        output.append(
            {
                "id": campaign["id"],
                "temperature_C": campaign["average_temperature_C"],
                "temperature_corrected_S7_9_m": params["S79"],
                "nominal_homogeneous_coefficient_kg_m-3": coefficient,
                "definition": "2*abs(mean_phi[(torque/(G*I0))*exp(-2i*phi)])",
                "uniform_cuboid_I0_kg_m2": core_i,
                "supplement_full_apparatus_I_kg_m2": full_i,
                "core_torque_divided_by_full_I_forbidden_mix_kg_m-3": coefficient * core_i / full_i,
                "normalization_rule": (
                    "The numerator and inertia must describe the same mass distribution. "
                    "The displayed mixed number is an identifiability diagnostic, not an alternative result."
                ),
                "exact_minimum_homogeneous_sphere_to_cuboid_clearance_over_azimuth_m": exact_minimum_clearance_over_azimuth(
                    params, source_centres, sphere_radii
                ),
                "quadrature": {
                    "primary_tensor_Gauss_order": 16,
                    "primary_azimuth_samples": 512,
                    "coarse_tensor_Gauss_order": 12,
                    "coarse_azimuth_samples": 256,
                    "coarse_coefficient_kg_m-3": coarse,
                    "primary_minus_coarse_kg_m-3": coefficient - coarse,
                },
                "public_input_sensitivity": sensitivity,
                "sensitivity_groups": groups,
                "main_table1_group_comparator_ppm": fields[
                    "main_table1_geometry_sensitivity_ppm"
                ][campaign["id"]],
                "identified_set": {
                    "conditional_homogeneous_functional": (
                        "K_AAF(X | P_pairwise-centred), where X is a chosen input domain and "
                        "P_pairwise-centred is the explicit nominal placement premise"
                    ),
                    "full_apparatus_from_unprocessed_public_fields": (
                        "K_AAF(X)+r_missing; the pinned release does not numerically bound "
                        "the central clamp/coating/density/3D-coordinate remainder"
                    ),
                },
                "temperature_transport": (
                    "For AAF-I, only S7,9 is transported from the published 23.7 C reference "
                    "using the published upper-pair coefficient; all other geometry fields remain "
                    "at their 23.7 C table values. For AAF-II/III the campaign and reference "
                    "temperatures coincide. This is a partial public temperature transport; no "
                    "complete campaign-temperature mass map is claimed."
                ),
            }
        )
    return output


def tos_parameters(fields: dict, run: dict) -> tuple[dict[str, float], dict[str, float]]:
    pendulum = fields["pendulums"][run["apparatus"]]
    spheres = fields["source_spheres"][run["apparatus"]]
    params = {
        "L": pendulum["length_m"],
        "W": pendulum["width_m"],
        "H": pendulum["height_m"],
        "M": pendulum["mass_kg"],
        "m1": spheres[0]["mass_kg"],
        "m2": spheres[1]["mass_kg"],
        "S": run["sphere_center_distance_m"],
    }
    steps = {
        "L": pendulum["u_length_m"],
        "W": pendulum["u_width_m"],
        "H": pendulum["u_height_m"],
        "M": pendulum["u_mass_kg"],
        "m1": spheres[0]["u_mass_kg"],
        "m2": spheres[1]["u_mass_kg"],
        "S": run["u_sphere_center_distance_m"],
    }
    return params, steps


def reconstruct_tos(fields: dict) -> list[dict]:
    output = []
    for run in fields["tos_runs"]:
        params, steps = tos_parameters(fields, run)
        coefficient, near, far = tos_stiffness(params, order=20)
        coarse, _, _ = tos_stiffness(params, order=16)
        sensitivity = central_sensitivities(
            lambda trial: tos_stiffness(trial, order=16)[0], params, steps
        )
        groups = {
            "pendulum_dimensions_ppm": group_rss(sensitivity, ("L", "W", "H")),
            "source_masses_ppm": group_rss(sensitivity, ("m1", "m2")),
            "horizontal_distance_ppm": group_rss(sensitivity, ("S",)),
        }
        pendulum = fields["pendulums"][run["apparatus"]]
        spheres = fields["source_spheres"][run["apparatus"]]
        base_centres = np.asarray(((-params["S"] / 2.0, 0.0, 0.0), (+params["S"] / 2.0, 0.0, 0.0)))
        radii = np.asarray([item["diameter_m"] / 2.0 for item in spheres])
        core_i = core_inertia(params)
        full_i = pendulum["supplement_full_apparatus_I_kg_m2"]
        output.append(
            {
                "id": run["id"],
                "apparatus": run["apparatus"],
                "temperature_C": run["temperature_C"],
                "sphere_center_distance_m": params["S"],
                "nominal_homogeneous_Delta_Cg_over_I_kg_m-3": coefficient,
                "near_Cg_over_I_kg_m-3": near,
                "far_Cg_over_I_kg_m-3": far,
                "definition": "d(torque/(G*I0))/dphi at phi=0 minus its value at phi=pi/2",
                "uniform_cuboid_I0_kg_m2": core_i,
                "supplement_full_apparatus_I_kg_m2": full_i,
                "core_curvature_divided_by_full_I_forbidden_mix_kg_m-3": coefficient
                * core_i
                / full_i,
                "normalization_rule": (
                    "The curvature numerator and inertia must describe the same mass distribution. "
                    "The displayed mixed number is an identifiability diagnostic, not an alternative result."
                ),
                "exact_minimum_homogeneous_sphere_to_cuboid_clearance_over_azimuth_m": exact_minimum_clearance_over_azimuth(
                    params, base_centres, radii
                ),
                "quadrature": {
                    "primary_tensor_Gauss_order": 20,
                    "coarse_tensor_Gauss_order": 16,
                    "coarse_coefficient_kg_m-3": coarse,
                    "primary_minus_coarse_kg_m-3": coefficient - coarse,
                },
                "public_input_sensitivity": sensitivity,
                "sensitivity_groups": groups,
                "main_table1_group_comparator_ppm": fields[
                    "main_table1_geometry_sensitivity_ppm"
                ][run["id"]],
                "identified_set": {
                    "conditional_homogeneous_functional": (
                        "K_TOS(X | P_axis-centred), where X is a chosen input domain and "
                        "P_axis-centred is the explicit nominal placement premise"
                    ),
                    "full_apparatus_from_unprocessed_public_fields": (
                        "K_TOS(X)+r_missing; the pinned release does not numerically bound "
                        "the central clamp/coating/density/3D-coordinate remainder"
                    ),
                },
                "temperature_transport": (
                    "The run-specific source separation is used at its table temperature. "
                    "TOS-I pendulum/source dimensions remain at their published 20.2 C reference "
                    "while runs are at 20.1 or 20.3 C; TOS-II fields and runs are at 21.5 C. "
                    "No unpublished full thermal transport is inserted."
                ),
            }
        )
    return output


def attach_post_calculation_comparators(results: list[dict], published: list[dict], key: str) -> None:
    """Open the processed table only after the geometry-first result exists."""

    by_id = {item["id"]: item for item in published}
    for result in results:
        comparator = by_id[result["id"]]
        nominal = result[key]
        processed = comparator["processed_coefficient_kg_m-3"]
        u_public = comparator["u_processed_coefficient_kg_m-3"]
        u_ideal = result["public_input_sensitivity"]["rss_standard_uncertainty_kg_m-3"]
        remainder = processed - nominal
        conservative_half = u_public + u_ideal
        result["post_calculation_processed_comparator"] = {
            "role": "POST_CALCULATION_ONLY__NOT_A_KERNEL_INPUT",
            **comparator,
            "processed_minus_nominal_remainder_kg_m-3": remainder,
            "processed_minus_nominal_relative_ppm": remainder / processed * 1e6,
            "unknown_covariance_conservative_sum_of_1sigma_half_widths_kg_m-3": conservative_half,
            "remainder_band_from_summed_1sigma_half_widths_kg_m-3": [
                remainder - conservative_half,
                remainder + conservative_half,
            ],
            "band_ceiling": (
                "This is interval addition of two quoted one-standard-uncertainty half widths "
                "under unknown covariance. It is not an independent reconstruction or coverage interval."
            ),
        }


def build_result() -> dict:
    source_hashes = verify_dependencies()
    fields = read_json(FIELDS)

    # These two calls are the entire geometry-first reconstruction.  The
    # processed-comparator file has deliberately not been opened yet.
    aaf = reconstruct_aaf(fields)
    tos = reconstruct_tos(fields)

    comparators = read_json(COMPARATORS)
    attach_post_calculation_comparators(
        aaf, comparators["AAF"], "nominal_homogeneous_coefficient_kg_m-3"
    )
    attach_post_calculation_comparators(
        tos, comparators["TOS"], "nominal_homogeneous_Delta_Cg_over_I_kg_m-3"
    )

    return {
        "schema": "WAC_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_RESULT_V001",
        "status": (
            "AAF_CONDITIONAL_HOMOGENEOUS_SOURCE_KERNEL_RECONSTRUCTED__"
            "TOS_CONDITIONAL_HOMOGENEOUS_STIFFNESS_KERNEL_RECONSTRUCTED__"
            "PUBLIC_GEOMETRY_SENSITIVITIES_RECOVERED__FULL_CENTRAL_REMAINDER_"
            "NOT_IDENTIFIED_FROM_UNPROCESSED_PUBLIC_FIELDS__GC16_NOT_CLOSED__NO_NEW_G"
        ),
        "source_and_dependency_hashes": source_hashes,
        "calculation_order": [
            "AAF ideal homogeneous coefficient from Extended Data Tables 1, 2 and 4",
            "ToS ideal homogeneous stiffness coefficients from Extended Data Tables 1, 2 and 4",
            "public-input sensitivity propagation and main Table 1 comparison",
            "post-calculation comparison to processed Supplementary Tables 2 and 3",
        ],
        "conditional_model_domain": {
            "detector": "uniform centred rectangular cuboid with Extended Data Table 1 dimensions and mass",
            "sources": (
                "mutually disjoint homogeneous spheres; shell theorem makes each source exactly a "
                "point mass at its stated centre for this domain"
            ),
            "layout": (
                "explicit nominal pairwise-centred realization consistent with Extended Data Figure 3; "
                "reported horizontal and vertical separations are retained independently, but the "
                "unpublished individual 3D CMM coordinates are not thereby reconstructed"
            ),
            "physics": "quasistatic Newtonian component of the declared weak-field Einstein endpoint",
            "numerics": (
                "tensor Gauss-Legendre cubature plus the exact analytic ToS azimuth derivative; "
                "the reported values are converged numerical evaluations without a certified "
                "quadrature error bound"
            ),
        },
        "AAF": aaf,
        "TOS": tos,
        "identifiability_ceiling": {
            "proved_collision": (
                "Extended Data Table 1 gives only the core cuboid while Supplementary Table 1 gives "
                "the inertia of the full pendulum assembly. Many unreported clamp/coating/ferrule mass "
                "maps share those scalars but produce different torque/curvature numerators. Mixing the "
                "core numerator with full-assembly inertia is therefore not a lawful reconstruction."
            ),
            "missing_public_geometry_fields": [
                "three-dimensional mass-coordinate or density map for clamp and ferrule",
                "coating mass/thickness map and density-inhomogeneity field",
                "individual three-dimensional CMM sphere centres relative to the rotation axis",
                "the two AAF shear coordinates left free by four pair separations plus overall centring",
                "source nonsphericity and density multipoles beyond the scalar summaries",
                "campaign covariance tying dimensions, masses, distances, alignment, and corrections",
            ],
            "full_kernel_identified_set": (
                "K_full=K_homogeneous(X)+r_missing is non-singleton and is not numerically bounded "
                "by the pinned unprocessed public fields. Supplementary Tables 2 and 3 publish a "
                "processed scalar image, but using it locates rather than independently reconstructs r_missing."
            ),
        },
        "residual_gap_to_GC16": {
            "now_supplied": [
                "AAF nominal m=2 finite-source forcing coefficient in the conditional homogeneous domain",
                "seven ToS nominal near-minus-far finite-source stiffness coefficients in that domain",
                "public geometry sensitivity ledger and a post-calculation central-remainder diagnostic",
            ],
            "still_missing": [
                "full detector/source/support/drive mass and stress measures with conserved-stress ownership",
                "row-level source trajectory and phase referenced event stream",
                "complete torsion and auxiliary-mode transfer including damping, coupling, and calibration covariance",
                "readout gain, delay, filters, nuisance design, and row-level observation covariance",
                "signed physical remainder ledger and bounds before the retarded inverse",
                "prospective held-out/null rows and independent source-scale covariance",
            ],
            "conclusion": (
                "The public release now supports a geometry-first nominal kernel and a sharp remainder "
                "ceiling, but not a real-data evaluation of GC16."
            ),
        },
        "claim_ceiling": [
            "no new or re-estimated value of G",
            "no RGRL or Gravity Formation Theory confirmation",
            "no claim of complete conserved-stress ownership",
            "no claim that a uniform nominal kernel equals the full HUST apparatus coefficient",
            "no use of a published processed coefficient as a calculation input",
            "no certified exact numerical value or rigorous quadrature error bound",
            "no claim that the conditional pairwise-centred AAF placement is uniquely fixed by the public pair distances",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=RESULT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build_result()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != encoded:
            raise SystemExit("stored RESULT.json does not match a clean reconstruction")
        print("RESULT.json matches clean public-geometry reconstruction")
    else:
        args.output.write_text(encoded, encoding="utf-8")
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
