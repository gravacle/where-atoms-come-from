#!/usr/bin/env python3
"""Independent verifier for the HUST nominal source-kernel lane.

This file does not import the reconstruction executable.  It reimplements the
finite homogeneous AAF harmonic and analytic ToS stiffness from the pinned
public fields, checks the post-calculation quarantine, and verifies the lane
manifest and seal.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

EXPECTED_SOURCE_HASHES = {
    "SOURCE/41586_2018_431_MOESM1_ESM.pdf": "5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb",
    "SOURCE/41586_2018_431_Tab1_ESM.jpg": "c84bc71bcd0115ccbbbdddd70dd1d755ffdabd37c91772f30cb51a264a200195",
    "SOURCE/41586_2018_431_Tab2_ESM.jpg": "96c40827af03f4de0715ea77bc69c0612d3ef94dc15eff035213e8a4dc0649c1",
    "SOURCE/41586_2018_431_Tab4_ESM.jpg": "567e5c8b953cba86e642bfd01b3880b262873606e2c69af44fa13a9ff4f629ce",
    "SOURCE/41586_2018_431_Fig6_ESM.jpg": "f29cdc1909149fc3a03264299889ef7f8346eb0fae6d7415cebbcf3d14163312",
    "SOURCE/nature_main_table1_error_budget.html": "23436d4be7600a7a9dffa02cc4167a20b6eea032a181e77899bb57bb90aa02e9",
}

EXPECTED_GC_HASHES = {
    "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/THEOREM.md": "cbf0733633ba93756b08dded7486a9be76beb572807693455c761bd36a8f0f5b",
    "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/PROTOCOL.md": "6ec7d8f0ce9a184d25612107dbfc294dd22d124ebe859831401e9cc0c8e8b819",
    "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/finite_apparatus_g_model.py": "6c17498d2d65f6420498ac559a97a2c3bbf49e110dd971da34b4c9c9bea2e4e4",
}

MANIFEST_FILES = {
    "HOSTILE_AUDIT_TRANSCRIPT.txt",
    "INDEPENDENT_HOSTILE_AUDIT.md",
    "PDF_SOURCE_READING.md",
    "PUBLISHED_COMPARATORS.json",
    "RESULT.json",
    "SELF_AUDIT.md",
    "SOURCE_CUSTODY.json",
    "SOURCE_FIELDS.json",
    "THEOREM.md",
    "hostile_audit_hust_nominal_source_kernels.py",
    "reconstruct_hust_nominal_source_kernels.py",
    "verify_hust_nominal_source_kernels.py",
    *EXPECTED_SOURCE_HASHES.keys(),
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def rule(p: dict, order: int) -> tuple[np.ndarray, ...]:
    node, weight = leggauss(order)
    x, y, z = np.meshgrid(
        node * p["L"] / 2.0,
        node * p["W"] / 2.0,
        node * p["H"] / 2.0,
        indexing="ij",
    )
    weights = np.einsum("i,j,k->ijk", weight, weight, weight).ravel() / 8.0
    return x.ravel(), y.ravel(), z.ravel(), weights


def inertia(p: dict) -> float:
    return p["M"] * (p["L"] ** 2 + p["W"] ** 2) / 12.0


def independent_aaf(p: dict, order: int = 12, angles: int = 256) -> float:
    x, y, z, w = rule(p, order)
    centres = np.array(
        [
            [-p["S79"] / 2, 0, p["S710"] / 2],
            [p["S79"] / 2, 0, p["S912"] / 2],
            [-p["S1012"] / 2, 0, -p["S710"] / 2],
            [p["S1012"] / 2, 0, -p["S912"] / 2],
        ],
        dtype=float,
    )
    masses = [p[name] for name in ("m7", "m9", "m10", "m12")]
    phases = np.arange(angles) * (2 * np.pi / angles)
    values = []
    for phi in phases:
        c, s = np.cos(phi), np.sin(phi)
        rotated = centres.copy()
        rotated[:, 0] = c * centres[:, 0]
        rotated[:, 1] = s * centres[:, 0]
        torque = 0.0
        for mass, (rx, ry, rz) in zip(masses, rotated):
            d2 = (rx - x) ** 2 + (ry - y) ** 2 + (rz - z) ** 2
            torque += mass * p["M"] * np.sum(w * (x * ry - y * rx) / d2**1.5)
        values.append(torque / inertia(p))
    return float(2 * abs(np.mean(np.asarray(values) * np.exp(-2j * phases))))


def independent_tos(p: dict, order: int = 16) -> tuple[float, float, float]:
    x, y, z, w = rule(p, order)
    values = []
    for phi in (0.0, np.pi / 2):
        c, s = np.cos(phi), np.sin(phi)
        total = 0.0
        for mass, sign in ((p["m1"], -1.0), (p["m2"], 1.0)):
            rx, ry = sign * p["S"] * c / 2, sign * p["S"] * s / 2
            d2 = (rx - x) ** 2 + (ry - y) ** 2 + z**2
            cross = x * ry - y * rx
            dot = x * rx + y * ry
            total += (
                mass
                * p["M"]
                * np.sum(w * (dot / d2**1.5 - 3 * cross**2 / d2**2.5))
                / inertia(p)
            )
        values.append(float(total))
    return values[0] - values[1], values[0], values[1]


def aaf_params(fields: dict, campaign: dict) -> tuple[dict, dict]:
    pend = fields["pendulums"]["AAF"]
    spheres = {row["id"]: row for row in fields["source_spheres"]["AAF"]}
    dist = fields["aaf_distances_23p7C"]
    dt = campaign["average_temperature_C"] - 23.7
    p = {
        "L": pend["length_m"], "W": pend["width_m"], "H": pend["height_m"], "M": pend["mass_kg"],
        "m7": spheres["7"]["mass_kg"], "m9": spheres["9"]["mass_kg"],
        "m10": spheres["10"]["mass_kg"], "m12": spheres["12"]["mass_kg"],
        "S79": dist["S7_9_m"] + dist["upper_horizontal_temperature_coefficient_m_per_C"] * dt,
        "S1012": dist["S10_12_m"], "S710": dist["S7_10_m"], "S912": dist["S9_12_m"],
    }
    u = {
        "L": pend["u_length_m"], "W": pend["u_width_m"], "H": pend["u_height_m"], "M": pend["u_mass_kg"],
        "m7": spheres["7"]["u_mass_kg"], "m9": spheres["9"]["u_mass_kg"],
        "m10": spheres["10"]["u_mass_kg"], "m12": spheres["12"]["u_mass_kg"],
        "S79": math.hypot(dist["u_S7_9_m"], dt * dist["u_upper_horizontal_temperature_coefficient_m_per_C"]),
        "S1012": dist["u_S10_12_m"], "S710": dist["u_S7_10_m"], "S912": dist["u_S9_12_m"],
    }
    return p, u


def tos_params(fields: dict, run: dict) -> tuple[dict, dict]:
    pend = fields["pendulums"][run["apparatus"]]
    spheres = fields["source_spheres"][run["apparatus"]]
    p = {
        "L": pend["length_m"], "W": pend["width_m"], "H": pend["height_m"], "M": pend["mass_kg"],
        "m1": spheres[0]["mass_kg"], "m2": spheres[1]["mass_kg"], "S": run["sphere_center_distance_m"],
    }
    u = {
        "L": pend["u_length_m"], "W": pend["u_width_m"], "H": pend["u_height_m"], "M": pend["u_mass_kg"],
        "m1": spheres[0]["u_mass_kg"], "m2": spheres[1]["u_mass_kg"], "S": run["u_sphere_center_distance_m"],
    }
    return p, u


def independent_sensitivities(function, p: dict, u: dict) -> tuple[dict[str, float], float]:
    nominal = function(p)
    components = {}
    for name, step in u.items():
        plus, minus = dict(p), dict(p)
        plus[name] += step
        minus[name] -= step
        components[name] = abs(function(plus) - function(minus)) / 2
    return components, math.sqrt(sum(value**2 for value in components.values())) / nominal * 1e6


def main() -> None:
    checks: list[str] = []

    def ck(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    for relative, expected in EXPECTED_SOURCE_HASHES.items():
        ck(f"source hash {relative}", digest(HERE / relative) == expected)
    for relative, expected in EXPECTED_GC_HASHES.items():
        ck(f"GC dependency hash {relative}", digest(ROOT / relative) == expected)

    fields_text = (HERE / "SOURCE_FIELDS.json").read_text(encoding="utf-8")
    fields = json.loads(fields_text)
    comparator_text = (HERE / "PUBLISHED_COMPARATORS.json").read_text(encoding="utf-8")
    comparators = json.loads(comparator_text)
    result = json.loads((HERE / "RESULT.json").read_text(encoding="utf-8"))
    ck("fields exclude processed coefficient", "processed_coefficient" not in fields_text)
    ck("comparator quarantine label", "POST_CALCULATION_ONLY" in comparators["hard_rule"])
    ck("pendulum transcription", fields["pendulums"]["AAF"]["length_m"] == 0.09105243)
    ck("sphere transcription", fields["source_spheres"]["AAF"][0]["mass_kg"] == 8.5435826)
    ck("distance transcription", fields["aaf_distances_23p7C"]["S7_9_m"] == 0.3422874)
    ck("ToS row count", len(fields["tos_runs"]) == 7)
    ck("AAF campaign count", len(fields["aaf_campaigns"]) == 3)
    ck("AAF conditional placement declared", "two shear degrees of freedom" in fields["conditional_layout_premise"]["AAF"])
    ck("temperature caption source located", "/tables/5" in fields["source_locations"]["distance_temperature_caption"])

    expected_aaf = {
        "AAF-I": 6926.660438859097,
        "AAF-II": 6926.700007763433,
        "AAF-III": 6926.700007763433,
    }
    expected_aaf_u_ppm = {
        "AAF-I": 10.691692680507465,
        "AAF-II": 10.687486717039393,
        "AAF-III": 10.687486717039393,
    }
    aaf_by_id = {row["id"]: row for row in result["AAF"]}
    for campaign in fields["aaf_campaigns"]:
        name = campaign["id"]
        p, u = aaf_params(fields, campaign)
        calculated = independent_aaf(p)
        ck(f"AAF independent coefficient {name}", math.isclose(calculated, expected_aaf[name], abs_tol=2e-9))
        ck(
            f"AAF stored coefficient {name}",
            math.isclose(aaf_by_id[name]["nominal_homogeneous_coefficient_kg_m-3"], calculated, abs_tol=2e-9),
        )
        components, u_ppm = independent_sensitivities(independent_aaf, p, u)
        ck(f"AAF independent sensitivity {name}", math.isclose(u_ppm, expected_aaf_u_ppm[name], abs_tol=2e-8))
        ck(
            f"AAF stored sensitivity {name}",
            math.isclose(aaf_by_id[name]["public_input_sensitivity"]["rss_relative_ppm"], u_ppm, abs_tol=2e-8),
        )
        horizontal = math.hypot(components["S79"], components["S1012"]) / calculated * 1e6
        vertical = math.hypot(components["S710"], components["S912"]) / calculated * 1e6
        ck(f"AAF horizontal class {name}", abs(horizontal - 8.98) < 0.006)
        ck(f"AAF vertical class {name}", abs(vertical - 5.79) < 0.004)

    expected_tos = [
        24915.18396314044, 24914.429047780784, 24914.429047780784,
        24914.243522929413, 24914.29470199404, 25005.34271885324, 25005.25902954411,
    ]
    tos_by_id = {row["id"]: row for row in result["TOS"]}
    for run, expected in zip(fields["tos_runs"], expected_tos):
        p, _ = tos_params(fields, run)
        calculated, near, far = independent_tos(p)
        stored = tos_by_id[run["id"]]
        ck(f"ToS independent coefficient {run['id']}", math.isclose(calculated, expected, abs_tol=3e-8))
        ck(
            f"ToS stored coefficient {run['id']}",
            math.isclose(stored["nominal_homogeneous_Delta_Cg_over_I_kg_m-3"], calculated, abs_tol=3e-8),
        )
        ck(f"ToS near/far signs {run['id']}", near > 0 and far < 0)
        ck(f"ToS analytic decomposition {run['id']}", math.isclose(near - far, calculated, abs_tol=1e-10))

    ck("AAF processed differences bounded", all(41 < abs(row["post_calculation_processed_comparator"]["processed_minus_nominal_relative_ppm"]) < 54 for row in result["AAF"]))
    ck("ToS processed differences bounded", all(91 < abs(row["post_calculation_processed_comparator"]["processed_minus_nominal_relative_ppm"]) < 104 for row in result["TOS"]))
    ck("all exact azimuth clearances positive", all(row["exact_minimum_homogeneous_sphere_to_cuboid_clearance_over_azimuth_m"] > 0 for group in (result["AAF"], result["TOS"]) for row in group))
    ck("AAF convergence", max(abs(row["quadrature"]["primary_minus_coarse_kg_m-3"]) for row in result["AAF"]) < 1e-9)
    ck("ToS convergence", max(abs(row["quadrature"]["primary_minus_coarse_kg_m-3"]) for row in result["TOS"]) < 1e-8)
    ck("full kernel not identified", "not numerically bounded" in result["identifiability_ceiling"]["full_kernel_identified_set"])
    ck("GC16 not closed", "not a real-data evaluation" in result["residual_gap_to_GC16"]["conclusion"])
    ck("six GC16 gaps", len(result["residual_gap_to_GC16"]["still_missing"]) == 6)
    ck("seven claim ceilings", len(result["claim_ceiling"]) == 7)
    ck("numerical exactness ceiling", "without a certified" in result["conditional_model_domain"]["numerics"])
    ck("AAF partial temperature transport", all("partial" in row["temperature_transport"] for row in result["AAF"]))

    analyzer = (HERE / "reconstruct_hust_nominal_source_kernels.py").read_text(encoding="utf-8")
    index_aaf = analyzer.index("aaf = reconstruct_aaf(fields)")
    index_tos = analyzer.index("tos = reconstruct_tos(fields)")
    index_comparator = analyzer.index("comparators = read_json(COMPARATORS)")
    ck("calculation before comparator", index_aaf < index_tos < index_comparator)
    ck("no accepted G literal in analyzer", "6.674" not in analyzer)
    ck("no comparator values in kernel source", "6926.352" not in analyzer and "24912.86" not in analyzer)
    ck("tampered source hash rejected", hashlib.sha256((b"x" + (HERE / "SOURCE/41586_2018_431_Tab1_ESM.jpg").read_bytes())).hexdigest() != EXPECTED_SOURCE_HASHES["SOURCE/41586_2018_431_Tab1_ESM.jpg"])

    # Comparator mutation cannot alter an independently calculated kernel.
    p, _ = aaf_params(fields, fields["aaf_campaigns"][0])
    before = independent_aaf(p)
    mutated = json.loads(comparator_text)
    mutated["AAF"][0]["processed_coefficient_kg_m-3"] *= 2
    after = independent_aaf(p)
    ck("comparator tamper cannot alter kernel", before == after)

    theorem = (HERE / "THEOREM.md").read_text(encoding="utf-8")
    audit = (HERE / "SELF_AUDIT.md").read_text(encoding="utf-8")
    pdf_reading = (HERE / "PDF_SOURCE_READING.md").read_text(encoding="utf-8")
    ck("theorem AAF first", theorem.index("AAF source coefficient reconstructed first") < theorem.index("ToS stiffness coefficients reconstructed second"))
    ck("theorem no new G", "not a new measurement of \\(G\\)" in theorem)
    ck("theorem GC ceiling", "does **not** close the real-data GC16 map" in theorem)
    ck("theorem conserved stress ceiling", "complete conserved" in theorem)
    ck("theorem placement premise", "two shear" in theorem and "not uniquely entailed" in theorem)
    ck("theorem numerical ceiling", "not a rigorous quadrature-error certificate" in theorem)
    ck("self audit hostile disposition", "independent hostile audit" in audit)
    ck("PDF visual record", "Visual inspection: completed" in pdf_reading)

    text_files = [path for path in HERE.rglob("*") if path.is_file() and path.suffix.lower() in {".md", ".json", ".py", ".txt", ".sha256", ".html"}]
    ck("no carriage returns", all(b"\r" not in path.read_bytes() for path in text_files))
    ck("no backspace bytes", all(b"\b" not in path.read_bytes() for path in text_files))
    ck("no form-feed bytes", all(b"\f" not in path.read_bytes() for path in text_files))

    run = subprocess.run(
        [sys.executable, str(HERE / "reconstruct_hust_nominal_source_kernels.py"), "--check"],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    ck("clean reconstruction exits zero", run.returncode == 0)
    ck("clean reconstruction text", run.stdout.strip() == "RESULT.json matches clean public-geometry reconstruction")

    hostile = subprocess.run(
        [sys.executable, str(HERE / "hostile_audit_hust_nominal_source_kernels.py")],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    hostile_transcript = (HERE / "HOSTILE_AUDIT_TRANSCRIPT.txt").read_text(encoding="utf-8")
    ck("independent hostile audit exits zero", hostile.returncode == 0)
    ck("independent hostile transcript exact", hostile.stdout == hostile_transcript)

    manifest = {}
    for line in (HERE / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        manifest[relative] = expected
    ck("manifest set", set(manifest) == MANIFEST_FILES)
    ck("manifest hashes", all(digest(HERE / name) == expected for name, expected in manifest.items()))
    ck("manifest no symlinks", all(not (HERE / name).is_symlink() for name in manifest))
    seal_hash, seal_name = (HERE / "LANE_SEAL.sha256").read_text(encoding="utf-8").strip().split("  ", 1)
    ck("seal target", seal_name == "MANIFEST.sha256")
    ck("seal hash", seal_hash == digest(HERE / "MANIFEST.sha256"))
    lane_files = {path.relative_to(HERE).as_posix() for path in HERE.rglob("*") if path.is_file()}
    ck("lane file set", lane_files == MANIFEST_FILES | {"MANIFEST.sha256", "LANE_SEAL.sha256", "VERIFICATION.txt"})

    text = (
        "HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_CHECK: PASS\n"
        f"Checks: {len(checks)}/{len(checks)}\n"
        "Official primary source objects pinned: 6\n"
        "AAF: 3 conditional homogeneous m=2 source coefficients independently recomputed\n"
        "ToS: 7 conditional homogeneous stiffness coefficients independently recomputed\n"
        "Public geometry sensitivities and post-calculation remainder signs reproduced\n"
        "Full apparatus kernel, GC16, new G, RGRL/GFT confirmation, and conserved-stress closure: not claimed\n"
    )
    (HERE / "VERIFICATION.txt").write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
