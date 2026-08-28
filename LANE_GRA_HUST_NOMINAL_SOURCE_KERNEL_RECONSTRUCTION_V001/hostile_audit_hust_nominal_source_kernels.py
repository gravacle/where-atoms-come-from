#!/usr/bin/env python3
"""Hostile independent audit of the HUST nominal source-kernel lane.

This executable intentionally does not import the production reconstruction.
It rebuilds the conditional Newtonian functionals from SOURCE_FIELDS.json,
uses different quadrature orders, checks the torque derivative with a finite
difference, proves the analytic shell-theorem clearance in the declared
geometry, and constructs a pair-distance-preserving AAF coordinate collision.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent

SOURCE_HASHES = {
    "SOURCE/41586_2018_431_MOESM1_ESM.pdf": "5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb",
    "SOURCE/41586_2018_431_Tab1_ESM.jpg": "c84bc71bcd0115ccbbbdddd70dd1d755ffdabd37c91772f30cb51a264a200195",
    "SOURCE/41586_2018_431_Tab2_ESM.jpg": "96c40827af03f4de0715ea77bc69c0612d3ef94dc15eff035213e8a4dc0649c1",
    "SOURCE/41586_2018_431_Tab4_ESM.jpg": "567e5c8b953cba86e642bfd01b3880b262873606e2c69af44fa13a9ff4f629ce",
    "SOURCE/41586_2018_431_Fig6_ESM.jpg": "f29cdc1909149fc3a03264299889ef7f8346eb0fae6d7415cebbcf3d14163312",
    "SOURCE/nature_main_table1_error_budget.html": "23436d4be7600a7a9dffa02cc4167a20b6eea032a181e77899bb57bb90aa02e9",
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def box_rule(p: dict[str, float], order: int) -> tuple[np.ndarray, ...]:
    nodes, weights_1d = leggauss(order)
    x, y, z = np.meshgrid(
        nodes * p["L"] / 2,
        nodes * p["W"] / 2,
        nodes * p["H"] / 2,
        indexing="ij",
    )
    weights = np.einsum("i,j,k->ijk", weights_1d, weights_1d, weights_1d).ravel() / 8
    return x.ravel(), y.ravel(), z.ravel(), weights


def inertia(p: dict[str, float]) -> float:
    return p["M"] * (p["L"] ** 2 + p["W"] ** 2) / 12


def aaf_parameters(fields: dict, campaign: dict) -> dict[str, float]:
    pendulum = fields["pendulums"]["AAF"]
    spheres = {row["id"]: row for row in fields["source_spheres"]["AAF"]}
    distance = fields["aaf_distances_23p7C"]
    delta_t = campaign["average_temperature_C"] - 23.7
    return {
        "L": pendulum["length_m"],
        "W": pendulum["width_m"],
        "H": pendulum["height_m"],
        "M": pendulum["mass_kg"],
        "m7": spheres["7"]["mass_kg"],
        "m9": spheres["9"]["mass_kg"],
        "m10": spheres["10"]["mass_kg"],
        "m12": spheres["12"]["mass_kg"],
        "S79": distance["S7_9_m"]
        + distance["upper_horizontal_temperature_coefficient_m_per_C"] * delta_t,
        "S1012": distance["S10_12_m"],
        "S710": distance["S7_10_m"],
        "S912": distance["S9_12_m"],
    }


def aaf_steps(fields: dict, campaign: dict) -> dict[str, float]:
    pendulum = fields["pendulums"]["AAF"]
    spheres = {row["id"]: row for row in fields["source_spheres"]["AAF"]}
    distance = fields["aaf_distances_23p7C"]
    delta_t = campaign["average_temperature_C"] - 23.7
    return {
        "L": pendulum["u_length_m"],
        "W": pendulum["u_width_m"],
        "H": pendulum["u_height_m"],
        "M": pendulum["u_mass_kg"],
        "m7": spheres["7"]["u_mass_kg"],
        "m9": spheres["9"]["u_mass_kg"],
        "m10": spheres["10"]["u_mass_kg"],
        "m12": spheres["12"]["u_mass_kg"],
        "S79": math.hypot(
            distance["u_S7_9_m"],
            delta_t * distance["u_upper_horizontal_temperature_coefficient_m_per_C"],
        ),
        "S1012": distance["u_S10_12_m"],
        "S710": distance["u_S7_10_m"],
        "S912": distance["u_S9_12_m"],
    }


def aaf_centres(p: dict[str, float], shear_x: float = 0.0, shear_z: float = 0.0) -> np.ndarray:
    """Family preserving all four pair separations and the overall centroid."""

    return np.asarray(
        [
            [-p["S79"] / 2 + shear_x, 0, +p["S710"] / 2 + shear_z],
            [+p["S79"] / 2 + shear_x, 0, +p["S912"] / 2 - shear_z],
            [-p["S1012"] / 2 - shear_x, 0, -p["S710"] / 2 + shear_z],
            [+p["S1012"] / 2 - shear_x, 0, -p["S912"] / 2 - shear_z],
        ],
        dtype=float,
    )


def aaf_coefficient(
    p: dict[str, float], *, order: int, angles: int, shear_x: float = 0.0, shear_z: float = 0.0
) -> float:
    x, y, z, weights = box_rule(p, order)
    centres = aaf_centres(p, shear_x, shear_z)
    masses = np.asarray([p["m7"], p["m9"], p["m10"], p["m12"]])
    phases = np.arange(angles) * (2 * math.pi / angles)
    values = np.empty(angles)
    for index, phi in enumerate(phases):
        c, s = math.cos(float(phi)), math.sin(float(phi))
        rotated = centres.copy()
        rotated[:, 0] = c * centres[:, 0]
        rotated[:, 1] = s * centres[:, 0]
        torque_per_g = 0.0
        for mass, (rx, ry, rz) in zip(masses, rotated):
            d2 = (rx - x) ** 2 + (ry - y) ** 2 + (rz - z) ** 2
            torque_per_g += mass * p["M"] * np.sum(weights * (x * ry - y * rx) / d2**1.5)
        values[index] = torque_per_g / inertia(p)
    return float(2 * abs(np.mean(values * np.exp(-2j * phases))))


def tos_parameters(fields: dict, run: dict) -> dict[str, float]:
    pendulum = fields["pendulums"][run["apparatus"]]
    spheres = fields["source_spheres"][run["apparatus"]]
    return {
        "L": pendulum["length_m"],
        "W": pendulum["width_m"],
        "H": pendulum["height_m"],
        "M": pendulum["mass_kg"],
        "m1": spheres[0]["mass_kg"],
        "m2": spheres[1]["mass_kg"],
        "S": run["sphere_center_distance_m"],
    }


def tos_steps(fields: dict, run: dict) -> dict[str, float]:
    pendulum = fields["pendulums"][run["apparatus"]]
    spheres = fields["source_spheres"][run["apparatus"]]
    return {
        "L": pendulum["u_length_m"],
        "W": pendulum["u_width_m"],
        "H": pendulum["u_height_m"],
        "M": pendulum["u_mass_kg"],
        "m1": spheres[0]["u_mass_kg"],
        "m2": spheres[1]["u_mass_kg"],
        "S": run["u_sphere_center_distance_m"],
    }


def symmetric_components(function, p: dict[str, float], steps: dict[str, float]) -> dict[str, float]:
    output = {}
    for name, step in steps.items():
        plus, minus = dict(p), dict(p)
        plus[name] += step
        minus[name] -= step
        output[name] = abs(function(plus) - function(minus)) / 2
    return output


def relative_group(components: dict[str, float], names: tuple[str, ...], nominal: float) -> float:
    return math.sqrt(sum(components[name] ** 2 for name in names)) / nominal * 1e6


def tos_torque(p: dict[str, float], phi: float, *, order: int) -> float:
    x, y, z, weights = box_rule(p, order)
    total = 0.0
    for mass, sign in ((p["m1"], -1.0), (p["m2"], +1.0)):
        rx = sign * p["S"] * math.cos(phi) / 2
        ry = sign * p["S"] * math.sin(phi) / 2
        d2 = (rx - x) ** 2 + (ry - y) ** 2 + z**2
        total += mass * p["M"] * np.sum(weights * (x * ry - y * rx) / d2**1.5)
    return float(total / inertia(p))


def tos_derivative(p: dict[str, float], phi: float, *, order: int) -> float:
    x, y, z, weights = box_rule(p, order)
    total = 0.0
    for mass, sign in ((p["m1"], -1.0), (p["m2"], +1.0)):
        rx = sign * p["S"] * math.cos(phi) / 2
        ry = sign * p["S"] * math.sin(phi) / 2
        d2 = (rx - x) ** 2 + (ry - y) ** 2 + z**2
        cross = x * ry - y * rx
        dot = x * rx + y * ry
        total += mass * p["M"] * np.sum(weights * (dot / d2**1.5 - 3 * cross**2 / d2**2.5)
        )
    return float(total / inertia(p))


def tos_coefficient(p: dict[str, float], *, order: int) -> tuple[float, float, float]:
    near = tos_derivative(p, 0.0, order=order)
    far = tos_derivative(p, math.pi / 2, order=order)
    return near - far, near, far


def exact_clearance(p: dict[str, float], centres: np.ndarray, radii: np.ndarray) -> float:
    corner_radius = math.hypot(p["L"] / 2, p["W"] / 2)
    values = []
    for centre, radius in zip(centres, radii):
        horizontal = max(math.hypot(float(centre[0]), float(centre[1])) - corner_radius, 0.0)
        vertical = max(abs(float(centre[2])) - p["H"] / 2, 0.0)
        values.append(math.hypot(horizontal, vertical) - float(radius))
    return min(values)


def pair_data(centres: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            centres[1, 0] - centres[0, 0],
            centres[3, 0] - centres[2, 0],
            centres[0, 2] - centres[2, 2],
            centres[1, 2] - centres[3, 2],
        ]
    )


def main() -> None:
    checks: list[str] = []

    def ck(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    for relative, expected in SOURCE_HASHES.items():
        ck(f"source custody {relative}", digest(HERE / relative) == expected)

    fields_text = (HERE / "SOURCE_FIELDS.json").read_text(encoding="utf-8")
    fields = json.loads(fields_text)
    result = json.loads((HERE / "RESULT.json").read_text(encoding="utf-8"))
    comparators_text = (HERE / "PUBLISHED_COMPARATORS.json").read_text(encoding="utf-8")
    comparators = json.loads(comparators_text)
    production = (HERE / "reconstruct_hust_nominal_source_kernels.py").read_text(encoding="utf-8")
    theorem = (HERE / "THEOREM.md").read_text(encoding="utf-8")

    ck("no processed coefficient in source fields", "processed_coefficient" not in fields_text)
    ck("no accepted G in source fields", "6.674" not in fields_text)
    ck("temperature coefficient transcription", fields["aaf_distances_23p7C"]["upper_horizontal_temperature_coefficient_m_per_C"] == -1.9e-6)
    ck("official temperature caption located", fields["source_locations"]["distance_temperature_caption"].startswith("https://www.nature.com/"))
    ck("AAF placement premise explicit", "two shear degrees of freedom" in fields["conditional_layout_premise"]["AAF"])

    stored_aaf = {row["id"]: row for row in result["AAF"]}
    aaf_metrics = []
    for campaign in fields["aaf_campaigns"]:
        p = aaf_parameters(fields, campaign)
        high = aaf_coefficient(p, order=20, angles=768)
        medium = aaf_coefficient(p, order=16, angles=512)
        stored = stored_aaf[campaign["id"]]
        ck(f"AAF rebuilt {campaign['id']}", math.isclose(high, stored["nominal_homogeneous_coefficient_kg_m-3"], abs_tol=6e-9))
        ck(f"AAF independent convergence {campaign['id']}", abs(high - medium) < 6e-9)
        ck(f"AAF core inertia {campaign['id']}", math.isclose(stored["uniform_cuboid_I0_kg_m2"], inertia(p), abs_tol=1e-18))
        aaf_metrics.append((campaign["id"], high, high - medium))

    p_aaf = aaf_parameters(fields, fields["aaf_campaigns"][1])
    base_centres = aaf_centres(p_aaf)
    sheared_centres = aaf_centres(p_aaf, shear_x=30e-6, shear_z=-30e-6)
    ck("AAF pair data invariant under shear collision", np.allclose(pair_data(base_centres), pair_data(sheared_centres), atol=2e-17, rtol=0))
    ck("AAF centroid invariant under shear collision", np.allclose(np.mean(sheared_centres, axis=0), 0.0, atol=2e-17, rtol=0))
    base_kernel = aaf_coefficient(p_aaf, order=16, angles=512)
    sheared_kernel = aaf_coefficient(p_aaf, order=16, angles=512, shear_x=30e-6, shear_z=-30e-6)
    ck("AAF pair distances do not identify kernel", abs(sheared_kernel - base_kernel) > 5e-4)
    ck("AAF-I S79 temperature transport", math.isclose(aaf_parameters(fields, fields["aaf_campaigns"][0])["S79"], 0.34228911, abs_tol=1e-15))
    ck("AAF partial temperature ceiling stored", all("partial" in row["temperature_transport"] for row in result["AAF"]))

    campaign = fields["aaf_campaigns"][1]
    p_sens = aaf_parameters(fields, campaign)
    aaf_components = symmetric_components(
        lambda trial: aaf_coefficient(trial, order=12, angles=256),
        p_sens,
        aaf_steps(fields, campaign),
    )
    aaf_nominal_sens = aaf_coefficient(p_sens, order=12, angles=256)
    aaf_groups = {
        "pendulum_dimensions": relative_group(aaf_components, ("L", "W", "H"), aaf_nominal_sens),
        "source_masses": relative_group(aaf_components, ("m7", "m9", "m10", "m12"), aaf_nominal_sens),
        "horizontal_distance": relative_group(aaf_components, ("S79", "S1012"), aaf_nominal_sens),
        "vertical_distance": relative_group(aaf_components, ("S710", "S912"), aaf_nominal_sens),
    }
    official_aaf = fields["main_table1_geometry_sensitivity_ppm"][campaign["id"]]
    ck("AAF independent dimension sensitivity", abs(aaf_groups["pendulum_dimensions"] - official_aaf["pendulum_dimensions"]) < 0.004)
    ck("AAF independent mass sensitivity", abs(aaf_groups["source_masses"] - official_aaf["source_masses"]) < 0.006)
    ck("AAF independent horizontal sensitivity", abs(aaf_groups["horizontal_distance"] - official_aaf["horizontal_distance"]) < 0.006)
    ck("AAF independent vertical sensitivity", abs(aaf_groups["vertical_distance"] - official_aaf["vertical_distance"]) < 0.004)

    phase = np.arange(512) * (2 * math.pi / 512)
    ck("AAF Fourier amplitude factor", math.isclose(2 * abs(np.mean(np.sin(2 * phase) * np.exp(-2j * phase))), 1.0, abs_tol=2e-15))

    stored_tos = {row["id"]: row for row in result["TOS"]}
    tos_metrics = []
    for run in fields["tos_runs"]:
        p = tos_parameters(fields, run)
        high, near, far = tos_coefficient(p, order=24)
        medium, _, _ = tos_coefficient(p, order=20)
        stored = stored_tos[run["id"]]
        ck(f"ToS rebuilt {run['id']}", math.isclose(high, stored["nominal_homogeneous_Delta_Cg_over_I_kg_m-3"], abs_tol=6e-8))
        ck(f"ToS independent convergence {run['id']}", abs(high - medium) < 6e-8)
        ck(f"ToS near/far sign {run['id']}", near > 0 and far < 0)
        ck(f"ToS core inertia {run['id']}", math.isclose(stored["uniform_cuboid_I0_kg_m2"], inertia(p), abs_tol=1e-18))
        tos_metrics.append((run["id"], high, near, far, high - medium))

    for run in (fields["tos_runs"][0], fields["tos_runs"][5]):
        p_sens = tos_parameters(fields, run)
        tos_components = symmetric_components(
            lambda trial: tos_coefficient(trial, order=16)[0],
            p_sens,
            tos_steps(fields, run),
        )
        tos_nominal_sens = tos_coefficient(p_sens, order=16)[0]
        groups = {
            "pendulum_dimensions": relative_group(tos_components, ("L", "W", "H"), tos_nominal_sens),
            "source_masses": relative_group(tos_components, ("m1", "m2"), tos_nominal_sens),
            "horizontal_distance": relative_group(tos_components, ("S",), tos_nominal_sens),
        }
        official = fields["main_table1_geometry_sensitivity_ppm"][run["id"]]
        ck(f"ToS independent dimension sensitivity {run['id']}", abs(groups["pendulum_dimensions"] - official["pendulum_dimensions"]) < 0.01)
        ck(f"ToS independent mass sensitivity {run['id']}", abs(groups["source_masses"] - official["source_masses"]) < 0.01)
        ck(f"ToS independent distance sensitivity {run['id']}", abs(groups["horizontal_distance"] - official["horizontal_distance"]) < 0.01)

    p_tos = tos_parameters(fields, fields["tos_runs"][0])
    step = 1e-5
    for label, phi in (("near", 0.0), ("far", math.pi / 2)):
        finite = (tos_torque(p_tos, phi + step, order=20) - tos_torque(p_tos, phi - step, order=20)) / (2 * step)
        analytic = tos_derivative(p_tos, phi, order=20)
        ck(f"ToS torque derivative sign/factor {label}", abs(finite - analytic) < 1e-5)

    aaf_radii = np.asarray([row["diameter_m"] / 2 for row in fields["source_spheres"]["AAF"]])
    aaf_clearance = exact_clearance(p_aaf, base_centres, aaf_radii)
    ck("AAF analytic shell clearance positive", aaf_clearance > 0.069)
    ck("AAF analytic clearance stored", math.isclose(aaf_clearance, stored_aaf["AAF-II"]["exact_minimum_homogeneous_sphere_to_cuboid_clearance_over_azimuth_m"], abs_tol=2e-15))
    aaf_stored = stored_aaf["AAF-II"]
    ck("AAF forbidden inertia mix arithmetic", math.isclose(aaf_stored["core_torque_divided_by_full_I_forbidden_mix_kg_m-3"], aaf_stored["nominal_homogeneous_coefficient_kg_m-3"] * aaf_stored["uniform_cuboid_I0_kg_m2"] / aaf_stored["supplement_full_apparatus_I_kg_m2"], abs_tol=2e-12))
    ck("AAF forbidden inertia mix not promoted", "not an alternative result" in aaf_stored["normalization_rule"])
    tos_clearances = []
    for run in fields["tos_runs"]:
        p = tos_parameters(fields, run)
        centres = np.asarray([[-p["S"] / 2, 0, 0], [p["S"] / 2, 0, 0]])
        radii = np.asarray([row["diameter_m"] / 2 for row in fields["source_spheres"][run["apparatus"]]])
        clearance = exact_clearance(p, centres, radii)
        tos_clearances.append(clearance)
        ck(f"ToS analytic clearance stored {run['id']}", math.isclose(clearance, stored_tos[run["id"]]["exact_minimum_homogeneous_sphere_to_cuboid_clearance_over_azimuth_m"], abs_tol=2e-15))
    ck("ToS analytic shell clearance positive", min(tos_clearances) > 0.004)
    tos_stored = stored_tos[fields["tos_runs"][0]["id"]]
    ck("ToS forbidden inertia mix arithmetic", math.isclose(tos_stored["core_curvature_divided_by_full_I_forbidden_mix_kg_m-3"], tos_stored["nominal_homogeneous_Delta_Cg_over_I_kg_m-3"] * tos_stored["uniform_cuboid_I0_kg_m2"] / tos_stored["supplement_full_apparatus_I_kg_m2"], abs_tol=2e-12))
    ck("ToS forbidden inertia mix not promoted", "not an alternative result" in tos_stored["normalization_rule"])

    index_aaf = production.index("aaf = reconstruct_aaf(fields)")
    index_tos = production.index("tos = reconstruct_tos(fields)")
    index_comparator = production.index("comparators = read_json(COMPARATORS)")
    ck("production code-level comparator quarantine", index_aaf < index_tos < index_comparator)
    ck("no comparator literal in production kernel", "6926.352" not in production and "24912.86" not in production)
    ck("comparator file explicitly post-calculation", "POST_CALCULATION_ONLY" in comparators["hard_rule"])
    ck("processed-minus-nominal signs AAF", all(row["post_calculation_processed_comparator"]["processed_minus_nominal_remainder_kg_m-3"] < 0 for row in result["AAF"]))
    ck("processed-minus-nominal signs ToS", all(row["post_calculation_processed_comparator"]["processed_minus_nominal_remainder_kg_m-3"] < 0 for row in result["TOS"]))
    ck("remainder bands disclaim coverage", all("not an independent reconstruction or coverage interval" in row["post_calculation_processed_comparator"]["band_ceiling"] for group in (result["AAF"], result["TOS"]) for row in group))

    ck("result numerical ceiling", "without a certified" in result["conditional_model_domain"]["numerics"])
    ck("theorem numerical ceiling", "not a rigorous quadrature-error certificate" in theorem)
    ck("theorem layout ceiling", "not uniquely entailed" in theorem and "public pair distances" in theorem)
    ck("full apparatus not identified", "non-singleton" in result["identifiability_ceiling"]["full_kernel_identified_set"])
    ck("GC16 remains open", "not a real-data evaluation of GC16" in result["residual_gap_to_GC16"]["conclusion"])
    ck("conserved stress remains open", any("conserved-stress" in item for item in result["claim_ceiling"]))
    ck("no new G remains explicit", result["claim_ceiling"][0] == "no new or re-estimated value of G")
    ck("no accepted G literal in production", "6.674" not in production)

    print("HUST_NOMINAL_SOURCE_KERNEL_HOSTILE_AUDIT: PASS")
    print(f"Checks: {len(checks)}/{len(checks)}")
    for name, value, delta in aaf_metrics:
        print(f"AAF {name}: independent={value:.12f} kg m^-3; high-minus-medium={delta:+.3e}")
    for name, value, near, far, delta in tos_metrics:
        print(f"ToS {name}: independent={value:.12f}; near={near:.12f}; far={far:.12f}; high-minus-medium={delta:+.3e}")
    print(f"AAF exact conditional clearance={aaf_clearance:.12f} m")
    print(f"ToS minimum exact conditional clearance={min(tos_clearances):.12f} m")
    print(f"AAF pair-distance collision: baseline={base_kernel:.12f}; sheared={sheared_kernel:.12f}; delta={sheared_kernel-base_kernel:+.12e} kg m^-3")
    print("AAF independent sensitivity classes (ppm): " + ", ".join(f"{key}={value:.6f}" for key, value in aaf_groups.items()))
    print("Numerical status: converged conditional evaluation; no certified quadrature error bound")
    print("Placement status: pairwise-centred AAF realization is an explicit premise, not uniquely public")
    print("Temperature status: public partial transport only; no complete campaign-temperature mass map")
    print("GC16/new G/RGRL-GFT/conserved-stress closure: not claimed")


if __name__ == "__main__":
    main()
