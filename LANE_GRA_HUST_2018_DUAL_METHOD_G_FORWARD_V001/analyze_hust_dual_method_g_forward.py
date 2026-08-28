#!/usr/bin/env python3
"""Reproduce the bounded HUST-2018 dual-method G-forward audit.

Only the official Nature supplementary PDF and Figure-2 source workbook are
used.  No accepted or CODATA value of G enters a fit, calibration, scan, or
comparison.  The script distinguishes figure-level raw-like series from
already-derived G summaries and never promotes the released processed
coefficients into a full GC16 apparatus fit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import posixpath
import re
import zipfile
import xml.etree.ElementTree as ET

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PDF = HERE / "SOURCE/41586_2018_431_MOESM1_ESM.pdf"
XLSX = HERE / "SOURCE/41586_2018_431_MOESM3_ESM.xlsx"
RESULT = HERE / "RESULT.json"

EXPECTED_DEPENDENCIES = {
    "supplementary_pdf": "5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb",
    "figure2_xlsx": "331ea6ee8a6a2558d4f6a8ffa4cd52e4c5d52d47f2e2a697889d0d3c8a3ad27a",
    "gc_theorem": "cbf0733633ba93756b08dded7486a9be76beb572807693455c761bd36a8f0f5b",
    "gc_protocol": "6ec7d8f0ce9a184d25612107dbfc294dd22d124ebe859831401e9cc0c8e8b819",
    "gc_model": "6c17498d2d65f6420498ac559a97a2c3bbf49e110dd971da34b4c9c9bea2e4e4",
}

DEPENDENCIES = {
    "supplementary_pdf": PDF,
    "figure2_xlsx": XLSX,
    "gc_theorem": ROOT / "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/THEOREM.md",
    "gc_protocol": ROOT / "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/PROTOCOL.md",
    "gc_model": ROOT / "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/finite_apparatus_g_model.py",
}

NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS_PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_dependencies() -> dict[str, str]:
    observed = {name: sha256(path) for name, path in DEPENDENCIES.items()}
    if observed != EXPECTED_DEPENDENCIES:
        raise RuntimeError(f"dependency hash mismatch: {observed}")
    return observed


def column_number(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference)
    if letters is None:
        raise ValueError(f"invalid cell reference: {reference}")
    number = 0
    for character in letters.group(0):
        number = number * 26 + ord(character) - ord("A") + 1
    return number


class XlsxSource:
    """Small read-only OOXML reader with explicit sheet/cell custody."""

    def __init__(self, path: Path):
        self.archive = zipfile.ZipFile(path)
        self.shared_strings = self._shared_strings()
        self.sheets = self._sheet_targets()

    def _shared_strings(self) -> list[str]:
        try:
            payload = self.archive.read("xl/sharedStrings.xml")
        except KeyError:
            return []
        root = ET.fromstring(payload)
        values = []
        for item in root.findall(f"{{{NS_MAIN}}}si"):
            values.append("".join(node.text or "" for node in item.iter(f"{{{NS_MAIN}}}t")))
        return values

    def _sheet_targets(self) -> dict[str, str]:
        workbook = ET.fromstring(self.archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(self.archive.read("xl/_rels/workbook.xml.rels"))
        target_by_id = {
            node.attrib["Id"]: node.attrib["Target"]
            for node in relationships.findall(f"{{{NS_PKG_REL}}}Relationship")
        }
        result = {}
        for sheet in workbook.findall(f".//{{{NS_MAIN}}}sheet"):
            relationship_id = sheet.attrib[f"{{{NS_REL}}}id"]
            target = target_by_id[relationship_id]
            if target.startswith("/"):
                target = target.lstrip("/")
            else:
                target = posixpath.normpath(posixpath.join("xl", target))
            result[sheet.attrib["name"]] = target
        return result

    def dimensions(self, sheet_name: str) -> str:
        with self.archive.open(self.sheets[sheet_name]) as stream:
            for _, element in ET.iterparse(stream, events=("end",)):
                if element.tag == f"{{{NS_MAIN}}}dimension":
                    return element.attrib["ref"]
                if element.tag == f"{{{NS_MAIN}}}sheetData":
                    break
        raise RuntimeError(f"dimension absent: {sheet_name}")

    def cells(self, sheet_name: str) -> dict[str, object]:
        values: dict[str, object] = {}
        with self.archive.open(self.sheets[sheet_name]) as stream:
            for _, cell in ET.iterparse(stream, events=("end",)):
                if cell.tag != f"{{{NS_MAIN}}}c":
                    continue
                reference = cell.attrib["r"]
                kind = cell.attrib.get("t")
                value_node = cell.find(f"{{{NS_MAIN}}}v")
                if kind == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iter(f"{{{NS_MAIN}}}t"))
                elif value_node is None:
                    value = None
                elif kind == "s":
                    value = self.shared_strings[int(value_node.text)]
                elif kind in {"str", "e"}:
                    value = value_node.text or ""
                elif kind == "b":
                    value = value_node.text == "1"
                else:
                    raw = value_node.text or ""
                    number = float(raw)
                    value = int(number) if number.is_integer() else number
                values[reference] = value
                cell.clear()
        return values

    def close(self) -> None:
        self.archive.close()


def row_values(cells: dict[str, object], start: int, stop: int, columns: str) -> list[list[object]]:
    result = []
    for row in range(start, stop + 1):
        result.append([cells.get(f"{column}{row}") for column in columns])
    return result


def period_series(rows: list[list[object]]) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    near: list[tuple[float, float]] = []
    far: list[tuple[float, float]] = []
    for time_day, period_near, period_far in rows:
        time = float(time_day)
        if period_near is not None:
            near.append((time, (2.0 * math.pi / float(period_near)) ** 2))
        if period_far is not None:
            far.append((time, (2.0 * math.pi / float(period_far)) ** 2))
    if len(near) != 10 or len(far) != 10:
        raise RuntimeError(f"unexpected ToS series sizes: {len(near)}, {len(far)}")
    return near, far


def interpolate(left: tuple[float, float], right: tuple[float, float], time: float) -> float:
    t0, y0 = left
    t1, y1 = right
    return y0 + (y1 - y0) * (time - t0) / (t1 - t0)


def aba_contrast(near: list[tuple[float, float]], far: list[tuple[float, float]]) -> dict:
    nfn = [
        interpolate(near[index], near[index + 1], far[index][0]) - far[index][1]
        for index in range(9)
    ]
    fnf = [
        near[index][1] - interpolate(far[index - 1], far[index], near[index][0])
        for index in range(1, 10)
    ]
    values = np.asarray(nfn + fnf, dtype=float)
    return {
        "n_f_n_count": len(nfn),
        "f_n_f_count": len(fnf),
        "contrast_count": values.size,
        "mean_delta_omega2_s-2": float(np.mean(values)),
        "minimum_delta_omega2_s-2": float(np.min(values)),
        "maximum_delta_omega2_s-2": float(np.max(values)),
        "uncertainty_status": "OVERLAPPING_TRIPLETS__NO_INDEPENDENCE_OR_COVERAGE_CLAIM",
    }


def common_quadratic_contrast(
    near: list[tuple[float, float]], far: list[tuple[float, float]]
) -> dict:
    rows = [(time, value, 1.0) for time, value in near] + [
        (time, value, 0.0) for time, value in far
    ]
    rows.sort()
    time = np.asarray([item[0] for item in rows], dtype=float)
    response = np.asarray([item[1] for item in rows], dtype=float)
    indicator = np.asarray([item[2] for item in rows], dtype=float)
    scaled = (time - np.mean(time)) / np.std(time)
    design = np.column_stack((np.ones(time.size), scaled, scaled**2, indicator))
    beta, _, rank, _ = np.linalg.lstsq(design, response, rcond=None)
    residual = response - design @ beta
    return {
        "model": "omega2=b0+b1*t_scaled+b2*t_scaled^2+delta*near_indicator",
        "row_count": time.size,
        "design_rank": int(rank),
        "near_minus_far_delta_omega2_s-2": float(beta[-1]),
        "residual_rms_s-2": float(np.sqrt(np.mean(residual**2))),
        "claim": "MINIMAL_COMMON_CURVATURE_DRIFT_DIAGNOSTIC__NOT_AUTHORS_EVENT_LEVEL_FIT",
    }


TOS_MECHANICS = {
    "fiber_1": {"I_m": 2.18e-5, "K_m": 1.2e-5, "I": 4.7705e-5, "K": 12.2e-9},
    "fiber_2": {"I_m": 2.18e-5, "K_m": 1.2e-5, "I": 4.7706e-5, "K": 47.4e-9},
    "fiber_3": {"I_m": 2.18e-5, "K_m": 1.2e-5, "I": 4.7705e-5, "K": 10.1e-9},
    "fiber_4": {"I_m": 2.21e-5, "K_m": 1.4e-5, "I": 4.6477e-5, "K": 10.6e-9},
}

TOS_ROWS = [
    ("fiber_1_first", "fiber_1", 24912.86, 0.23, 1.662732, 0.000017, 6.674154, 0.000095),
    ("fiber_1_repeated", "fiber_1", 24912.12, 0.23, 1.662699, 0.000018, 6.674222, 0.000098),
    ("fiber_2", "fiber_2", 24912.15, 0.23, 1.662698, 0.000051, 6.674237, 0.000219),
    ("fiber_3_first", "fiber_3", 24911.70, 0.22, 1.662684, 0.000020, 6.674274, 0.000103),
    ("fiber_3_repeated", "fiber_3", 24911.72, 0.22, 1.662683, 0.000017, 6.674266, 0.000094),
    ("fiber_4_first", "fiber_4", 25003.05, 0.25, 1.668719, 0.000023, 6.674017, 0.000117),
    ("fiber_4_repeated", "fiber_4", 25002.95, 0.25, 1.668734, 0.000023, 6.674105, 0.000116),
]


def tos_summary_forward() -> list[dict]:
    output = []
    for name, fiber, coefficient, coefficient_u, delta_w2, delta_w2_u, reported, reported_u in TOS_ROWS:
        mechanics = TOS_MECHANICS[fiber]
        magnetic = (
            mechanics["I_m"] * mechanics["K"] ** 2
            / (mechanics["I"] * mechanics["K_m"] ** 2)
        )
        raw = delta_w2 * 1e-6 / coefficient
        magnetic_only = raw * (1.0 + magnetic)
        reported_si = reported * 1e-11
        output.append({
            "id": name,
            "fiber": fiber,
            "delta_Cg_over_I_kg_m-3": coefficient,
            "delta_Cg_over_I_standard_u_kg_m-3": coefficient_u,
            "mean_delta_omega2_s-2": delta_w2 * 1e-6,
            "mean_delta_omega2_standard_u_s-2": delta_w2_u * 1e-6,
            "raw_response_over_source_ratio_SI": raw,
            "magnetic_damper_factor_minus_one": magnetic,
            "magnetic_damper_ppm": magnetic * 1e6,
            "magnetic_only_ratio_SI": magnetic_only,
            "published_G_summary_SI_comparison_only": reported_si,
            "published_G_summary_standard_u_SI": reported_u * 1e-11,
            "magnetic_only_minus_published_relative_ppm": (magnetic_only / reported_si - 1.0) * 1e6,
            "unowned_full_bracket_residual_ppm": (
                reported_si / raw - (1.0 + magnetic)
            ) * 1e6,
        })
    return output


AAF_ROWS = [
    {
        "id": "AAF-I", "P_abs_kg_m-3": 6926.352, "P_u_kg_m-3": 0.074,
        "alpha_nrad_s-2": 462.0912, "alpha_u_nrad_s-2": 0.0016,
        "I_m": 2.401e-5, "K_m": 1.199e-5, "I": 2.776e-5, "K": 6.313e-9,
        "published_G_e-11": 6.674534, "published_G_u_e-11": 0.000083,
    },
    {
        "id": "AAF-II", "P_abs_kg_m-3": 6926.334, "P_u_kg_m-3": 0.075,
        "alpha_nrad_s-2": 462.0791, "alpha_u_nrad_s-2": 0.0012,
        "I_m": 2.401e-5, "K_m": 1.199e-5, "I": 2.776e-5, "K": 6.313e-9,
        "published_G_e-11": 6.674375, "published_G_u_e-11": 0.000082,
    },
    {
        "id": "AAF-III", "P_abs_kg_m-3": 6926.415, "P_u_kg_m-3": 0.074,
        "alpha_nrad_s-2": 462.2941, "alpha_u_nrad_s-2": 0.0006,
        "I_m": 2.404e-5, "K_m": 21.24e-5, "I": 2.776e-5, "K": 6.313e-9,
        "published_G_e-11": 6.674535, "published_G_u_e-11": 0.000075,
    },
]


def aaf_summary_forward() -> list[dict]:
    output = []
    for source in AAF_ROWS:
        correction = 1.0 + (source["K"] / source["K_m"]) * (source["I_m"] / source["I"])
        raw = source["alpha_nrad_s-2"] * 1e-9 / source["P_abs_kg_m-3"]
        calculated = raw * correction
        published = source["published_G_e-11"] * 1e-11
        output.append({
            **source,
            "raw_alpha_over_P_SI": raw,
            "magnetic_damper_factor": correction,
            "magnetic_damper_ppm": (correction - 1.0) * 1e6,
            "recomputed_G_SI": calculated,
            "published_G_summary_SI_comparison_only": published,
            "recomputed_minus_published_relative_ppm": (calculated / published - 1.0) * 1e6,
        })
    return output


def aaf_time_stream_fit(rows: list[list[object]]) -> dict:
    samples = [(float(time), float(value)) for time, value in rows if time is not None and value is not None]
    if len(samples) != 9999:
        raise RuntimeError(f"unexpected acceleration sample count: {len(samples)}")
    time_all = np.asarray([item[0] for item in samples], dtype=float)
    response_all = np.asarray([item[1] for item in samples], dtype=float)
    time = time_all[:7200]
    response = response_all[:7200]
    centered = (time - np.mean(time)) / (time[-1] - time[0])
    source_frequency = 1.0 / 600.0

    base = np.column_stack((
        np.ones(time.size), centered,
        np.sin(2.0 * math.pi * source_frequency * time),
        np.cos(2.0 * math.pi * source_frequency * time),
    ))
    base_beta, _, _, _ = np.linalg.lstsq(base, response, rcond=None)
    source_only_amplitude = float(math.hypot(base_beta[-2], base_beta[-1]))

    aaf_i = np.column_stack((
        np.ones(time.size), centered,
        np.sin(2.0 * math.pi * 0.0025 * time),
        np.cos(2.0 * math.pi * 0.0025 * time),
    ))
    aaf_i_beta, _, _, _ = np.linalg.lstsq(aaf_i, response, rcond=None)
    aaf_i_amplitude = float(math.hypot(aaf_i_beta[-2], aaf_i_beta[-1]))

    best: tuple[float, float, np.ndarray] | None = None
    for index in range(6001):
        background_frequency = 0.0005 + index * 1e-7
        design = np.column_stack((
            base,
            np.sin(2.0 * math.pi * background_frequency * time),
            np.cos(2.0 * math.pi * background_frequency * time),
        ))
        beta, _, _, _ = np.linalg.lstsq(design, response, rcond=None)
        residual = response - design @ beta
        rss = float(residual @ residual)
        if best is None or rss < best[0]:
            best = (rss, background_frequency, beta)
    assert best is not None
    rss, background_frequency, beta = best
    return {
        "workbook_available_sample_count": len(samples),
        "workbook_time_range_s": [float(time_all[0]), float(time_all[-1])],
        "caption_scored_two_hour_sample_count": time.size,
        "caption_scored_time_range_s": [float(time[0]), float(time[-1])],
        "source_frequency_Hz_fixed_from_AAF_II_III": source_frequency,
        "source_only_fit_amplitude_nrad_s-2": source_only_amplitude,
        "AAF_I_2p5mHz_control_amplitude_nrad_s-2": aaf_i_amplitude,
        "background_grid_Hz": {"minimum": 0.0005, "maximum": 0.0011, "step": 1e-7},
        "two_tone_best_background_frequency_Hz": background_frequency,
        "two_tone_source_amplitude_nrad_s-2": float(math.hypot(beta[2], beta[3])),
        "two_tone_background_amplitude_nrad_s-2": float(math.hypot(beta[4], beta[5])),
        "two_tone_residual_rms_nrad_s-2": math.sqrt(rss / time.size),
        "claim": "FIGURE_LEVEL_SPECTRAL_SEPARABILITY_DIAGNOSTIC__NOT_CAMPAIGN_ALPHA_OR_G_EXTRACTION",
    }


def workbook_outcome_summaries(cells_c: dict[str, object], cells_f: dict[str, object]) -> dict:
    tos = [
        {"number": int(cells_c[f"B{row}"]), "G_e-11": float(cells_c[f"C{row}"]), "u_e-11": float(cells_c[f"D{row}"])}
        for row in range(3, 10)
    ]
    groups = {
        "AAF-I": ("B", "C"),
        "AAF-II": ("D", "E"),
        "AAF-III": ("F", "G"),
    }
    aaf = {}
    for name, (value_column, uncertainty_column) in groups.items():
        rows = []
        for row in range(3, 40):
            value = cells_f.get(f"{value_column}{row}")
            if value is not None:
                rows.append({
                    "workbook_row": row,
                    "G_e-11": float(value),
                    "u_e-11": float(cells_f[f"{uncertainty_column}{row}"]),
                })
        values = np.asarray([item["G_e-11"] for item in rows], dtype=float)
        aaf[name] = {
            "count": len(rows),
            "arithmetic_mean_G_e-11": float(np.mean(values)),
            "minimum_G_e-11": float(np.min(values)),
            "maximum_G_e-11": float(np.max(values)),
            "rows": rows,
        }
    return {
        "ToS_seven_already_derived_G_rows": tos,
        "ToS_combined_G_e-11": float(cells_c["E10"]),
        "ToS_combined_u_e-11": float(cells_c["F10"]),
        "AAF_29_already_derived_G_rows_by_campaign": aaf,
        "AAF_combined_G_e-11": float(cells_f["H39"]),
        "AAF_combined_u_e-11": float(cells_f["I39"]),
    }


def analyze() -> dict:
    hashes = check_dependencies()
    source = XlsxSource(XLSX)
    try:
        dimensions = {name: source.dimensions(name) for name in ("a", "b", "c", "d", "e", "f")}
        cells_a = source.cells("a")
        cells_b = source.cells("b")
        cells_c = source.cells("c")
        cells_e = source.cells("e")
        cells_f = source.cells("f")
    finally:
        source.close()

    expected_dimensions = {
        "a": "B2:D22", "b": "B2:D22", "c": "B1:F10",
        "d": "B2:G129604", "e": "B2:C10001", "f": "A1:I39",
    }
    if dimensions != expected_dimensions:
        raise RuntimeError(f"workbook dimension mismatch: {dimensions}")
    expected_headers = {
        "a": [cells_a.get("B2"), cells_a.get("C2"), cells_a.get("D2")],
        "b": [cells_b.get("B2"), cells_b.get("C2"), cells_b.get("D2")],
        "e": [cells_e.get("B2"), cells_e.get("C2")],
    }
    if expected_headers != {
        "a": ["Time (day)", "Period at near position (s)", "Period at far position (s)"],
        "b": ["Time (day)", "Period at near position (s)", "Period at far position (s)"],
        "e": ["Time (s)", "Acceleration (nrad/s2)"],
    }:
        raise RuntimeError(f"workbook header mismatch: {expected_headers}")

    near_present, far_present = period_series(row_values(cells_a, 3, 22, "BCD"))
    near_background, far_background = period_series(row_values(cells_b, 3, 22, "BCD"))
    aba_present = aba_contrast(near_present, far_present)
    aba_background = aba_contrast(near_background, far_background)
    quadratic_present = common_quadratic_contrast(near_present, far_present)
    quadratic_background = common_quadratic_contrast(near_background, far_background)
    aba_net = (
        aba_present["mean_delta_omega2_s-2"]
        - aba_background["mean_delta_omega2_s-2"]
    )
    quadratic_net = (
        quadratic_present["near_minus_far_delta_omega2_s-2"]
        - quadratic_background["near_minus_far_delta_omega2_s-2"]
    )
    printed_delta = 1.662699e-6
    printed_delta_u = 18e-12
    repeated_coefficient = 24912.12

    summaries = workbook_outcome_summaries(cells_c, cells_f)
    tos_combined = summaries["ToS_combined_G_e-11"] * 1e-11
    aaf_combined = summaries["AAF_combined_G_e-11"] * 1e-11
    tos_u = summaries["ToS_combined_u_e-11"] * 1e-11
    aaf_u = summaries["AAF_combined_u_e-11"] * 1e-11
    cross_difference = aaf_combined - tos_combined

    return {
        "schema": "WAC_HUST_2018_DUAL_METHOD_G_FORWARD_V001",
        "status": (
            "OFFICIAL_FIGURE_LEVEL_TOS_RESPONSE_RECOVERED__THREE_PROCESSED_AAF_"
            "SOURCE_RESPONSE_FORWARDS_REPRODUCED__FULL_GC16_NOT_EXECUTABLE__NO_NEW_G"
        ),
        "source": {
            "paper": "Li_et_al_Nature_560_582_588_2018",
            "doi": "10.1038/s41586-018-0431-5",
            "article_url": "https://www.nature.com/articles/s41586-018-0431-5",
            "dependency_sha256": hashes,
            "pinned_official_files": [
                {
                    "file": "SOURCE/41586_2018_431_MOESM1_ESM.pdf",
                    "url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-018-0431-5/MediaObjects/41586_2018_431_MOESM1_ESM.pdf",
                    "role": "equations_and_supplementary_tables_1_to_4",
                },
                {
                    "file": "SOURCE/41586_2018_431_MOESM3_ESM.xlsx",
                    "url": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-018-0431-5/MediaObjects/41586_2018_431_MOESM3_ESM.xlsx",
                    "role": "official_source_data_for_main_Figure_2",
                },
            ],
            "pdf_page_custody": {
                "ToS_and_AAF_equations": "PDF pages 4-6; printed SI pages 3-5",
                "two_stage_mechanics": "PDF page 19; Supplementary Table 1",
                "ToS_coefficients": "PDF page 20; Supplementary Table 2",
                "AAF_coefficients": "PDF page 21; Supplementary Table 3",
            },
            "workbook_dimensions": dimensions,
        },
        "exact_equations": {
            "ToS": (
                "G=(Delta_omega2/(Delta_Cg/I))*"
                "[1-Delta_K/(I*Delta_omega2)+I_m*K^2/(I*K_m^2)]"
            ),
            "AAF": (
                "G=(alpha_t(2*omega_d)/abs(sum_l>=2 P_g,l,2))*"
                "[1+(K/K_m)*(I_m/I)]"
            ),
            "GC16_target": (
                "y=C*[(p*a+r_theta+lambda*r_x/d_x)/(d_theta-p*k_g-lambda^2/d_x)+d_h]"
                "+B*eta+epsilon"
            ),
        },
        "ToS_figure_level_response": {
            "source_present_panel_a": aba_present,
            "source_absent_background_panel_b": aba_background,
            "A_B_A_background_subtracted_delta_omega2_s-2": aba_net,
            "A_B_A_minus_printed_delta_omega2_s-2": aba_net - printed_delta,
            "A_B_A_minus_printed_in_printed_sigma": (aba_net - printed_delta) / printed_delta_u,
            "common_quadratic_present": quadratic_present,
            "common_quadratic_background": quadratic_background,
            "common_quadratic_background_subtracted_delta_omega2_s-2": quadratic_net,
            "common_quadratic_minus_printed_delta_omega2_s-2": quadratic_net - printed_delta,
            "common_quadratic_minus_printed_in_printed_sigma": (
                quadratic_net - printed_delta
            ) / printed_delta_u,
            "printed_repeated_fiber_1_delta_omega2_s-2": printed_delta,
            "printed_repeated_fiber_1_delta_omega2_standard_u_s-2": printed_delta_u,
            "printed_repeated_fiber_1_delta_Cg_over_I_kg_m-3": repeated_coefficient,
            "A_B_A_response_over_processed_source_coefficient_SI": aba_net / repeated_coefficient,
            "common_quadratic_response_over_processed_source_coefficient_SI": quadratic_net / repeated_coefficient,
            "ceiling": (
                "THE_40_PERIOD_VALUES_ARE_THREE_DAY_PERIOD_SUMMARIES_NOT_0p5_SECOND_ANGLE_DATA;_"
                "DELTA_Cg_OVER_I_IS_A_PROCESSED_GEOMETRY_COEFFICIENT"
            ),
        },
        "ToS_seven_processed_coefficient_forwards": tos_summary_forward(),
        "AAF_three_processed_coefficient_forwards": aaf_summary_forward(),
        "AAF_figure_level_acceleration_stream": aaf_time_stream_fit(
            row_values(cells_e, 3, 10001, "BC")
        ),
        "already_derived_outcome_rows": summaries,
        "cross_method_source_model_stress": {
            "ToS_combined_SI": tos_combined,
            "AAF_combined_SI": aaf_combined,
            "AAF_minus_ToS_SI": cross_difference,
            "difference_relative_to_midpoint_ppm": (
                cross_difference / ((aaf_combined + tos_combined) / 2.0) * 1e6
            ),
            "quoted_standard_u_SI": {"ToS": tos_u, "AAF": aaf_u},
            "z_if_and_only_if_cross_method_covariance_is_zero": (
                cross_difference / math.sqrt(tos_u**2 + aaf_u**2)
            ),
            "cross_method_covariance_released": False,
            "interpretation": (
                "SOURCE_MODEL_AND_SYSTEMATICS_STRESS_ONLY__NOT_NEW_PHYSICS_AND_NOT_A_NEW_G_ESTIMATE"
            ),
        },
        "public_field_inventory": {
            "raw_like_or_intermediate": {
                "ToS_source_present_period_summaries": 20,
                "ToS_source_absent_background_period_summaries": 20,
                "AAF_acceleration_1_second_samples": 9999,
                "AAF_residual_and_free_twist_PSD_rows": "sheet_d_large_frequency_domain_arrays",
            },
            "already_derived": {
                "ToS_G_rows": 7,
                "AAF_G_rows": 29,
                "ToS_combined_G_rows": 1,
                "AAF_combined_G_rows": 1,
            },
            "processed_source_response_coefficients_in_pdf": {
                "ToS_Delta_Cg_over_I_and_mean_Delta_omega2_rows": 7,
                "AAF_sum_P_and_mean_alpha_rows": 3,
            },
        },
        "GC16_execution": {
            "full_real_apparatus_execution": False,
            "closure_gained": [
                "one Figure-2 ToS source-present/background period contrast recovered by A-B-A",
                "same ToS contrast independently recovered by panel-wise common-quadratic drift fits and background subtraction",
                "seven ToS processed response/source quotients and the known magnetic term evaluated",
                "three AAF processed response/source/mechanical equations reproduced within printed rounding",
                "29 AAF derived outcomes regroup to the three published campaign central means",
                "two-tone fit separates the approximately 462 nrad/s2 source response from the approximately 77 nrad/s2 lab-fixed background in the released two-hour stream",
                "cross-method 44.95 ppm source-model stress quantified without accepted G",
            ],
            "missing_fields": [
                "event-level 0.5 s ToS angle and environmental channels with run/configuration labels",
                "full 20 kHz AAF encoder/controller records and campaign-level alpha time series",
                "finite source/detector mass-coordinate-density files used to compute Delta_Cg and P_g,l,m",
                "row-level trajectory, source-scale calibration, and their joint covariance",
                "complete frequency-dependent torsion/support/auxiliary-mode transfer and readout calibration",
                "signed row-level correction and physical-remainder ledger, including ToS Delta_K ownership",
                "observation, calibration, campaign, and cross-method covariance matrices",
                "complete conserved apparatus stress/source ledger including drive and support reactions",
                "predeclared nuisance design, likelihood domain, null rows, and held-out rows required by GC protocol",
            ],
            "identifiability": (
                "PROCESSED_COEFFICIENT_RATIOS_IDENTIFIED__RAW_GEOMETRY_TRANSFER_REMAINDER_AND_GLOBAL_SOURCE_SCALE_NOT_SEPARATELY_REIDENTIFIED"
            ),
        },
        "claim_ceiling": {
            "new_G_measurement": False,
            "accepted_G_used_as_input": False,
            "full_GC16_execution": False,
            "RGRL_or_GFT_confirmed": False,
            "lineage_charge_inferred": False,
            "raw_event_level_reanalysis": False,
            "strongest_claim": (
                "OFFICIAL_HUST_RELEASE_CLOSES_A_REAL_PROCESSED_SOURCE_RESPONSE_FORWARD_AND_"
                "LIMITED_FIGURE_LEVEL_RESPONSE_EXTRACTION__NOT_THE_CALIBRATED_FULL_APPARATUS_G_CROSSCHECK"
            ),
        },
        "next_step": (
            "OBTAIN_OR_RECONSTRUCT_THE_AUTHORS_FINITE_MASS_GEOMETRY_FILES_FULL_TRANSFER_AND_"
            "ROW_COVARIANCE;_THEN_RECOMPUTE_DELTA_Cg_AND_P_FROM_GEOMETRY_AND_RUN_GC16_WITH_"
            "THE_COMPLETE_REMAINDER_LEDGER"
        ),
    }


def serialized(result: dict) -> str:
    return json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="compare with the sealed RESULT.json")
    args = parser.parse_args()
    result = analyze()
    payload = serialized(result)
    if args.check:
        if not RESULT.is_file():
            raise RuntimeError("RESULT.json absent")
        if RESULT.read_text(encoding="utf-8") != payload:
            raise RuntimeError("RESULT.json is stale")
        print("HUST_DUAL_METHOD_ANALYZER_CHECK: PASS")
    else:
        RESULT.write_text(payload, encoding="utf-8")
        print(payload, end="")


if __name__ == "__main__":
    main()
