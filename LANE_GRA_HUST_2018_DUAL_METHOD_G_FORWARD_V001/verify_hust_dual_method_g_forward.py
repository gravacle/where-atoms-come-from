#!/usr/bin/env python3
"""Independent verifier for the bounded HUST-2018 G-forward lane."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import posixpath
import subprocess
import sys
import zipfile
import xml.etree.ElementTree as ET

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PDF = HERE / "SOURCE/41586_2018_431_MOESM1_ESM.pdf"
XLSX = HERE / "SOURCE/41586_2018_431_MOESM3_ESM.xlsx"
NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG = "http://schemas.openxmlformats.org/package/2006/relationships"

EXPECTED_HASHES = {
    "SOURCE/41586_2018_431_MOESM1_ESM.pdf": "5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb",
    "SOURCE/41586_2018_431_MOESM3_ESM.xlsx": "331ea6ee8a6a2558d4f6a8ffa4cd52e4c5d52d47f2e2a697889d0d3c8a3ad27a",
    "../LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/THEOREM.md": "cbf0733633ba93756b08dded7486a9be76beb572807693455c761bd36a8f0f5b",
    "../LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/PROTOCOL.md": "6ec7d8f0ce9a184d25612107dbfc294dd22d124ebe859831401e9cc0c8e8b819",
    "../LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/finite_apparatus_g_model.py": "6c17498d2d65f6420498ac559a97a2c3bbf49e110dd971da34b4c9c9bea2e4e4",
}

MANIFEST_FILES = {
    "INDEPENDENT_HOSTILE_AUDIT.md",
    "READINESS.md",
    "RESULT.json",
    "SELF_AUDIT.md",
    "SOURCE_CUSTODY.json",
    "THEOREM.md",
    "analyze_hust_dual_method_g_forward.py",
    "verify_hust_dual_method_g_forward.py",
    "verify_hust_dual_method_hostile_audit.py",
    "SOURCE/41586_2018_431_MOESM1_ESM.pdf",
    "SOURCE/41586_2018_431_MOESM3_ESM.xlsx",
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def workbook_cells(path: Path) -> tuple[dict[str, str], dict[str, dict[str, object]]]:
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in shared_root.findall(f"{{{NS}}}si"):
                shared.append("".join(node.text or "" for node in item.iter(f"{{{NS}}}t")))
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {
            item.attrib["Id"]: item.attrib["Target"]
            for item in relations.findall(f"{{{PKG}}}Relationship")
        }
        sheet_paths = {}
        for sheet in workbook.findall(f".//{{{NS}}}sheet"):
            target = targets[sheet.attrib[f"{{{REL}}}id"]]
            if target.startswith("/"):
                target = target.lstrip("/")
            else:
                target = posixpath.normpath(posixpath.join("xl", target))
            sheet_paths[sheet.attrib["name"]] = target

        dimensions = {}
        outputs = {}
        for name in ("a", "b", "c", "e", "f"):
            root = ET.fromstring(archive.read(sheet_paths[name]))
            dimensions[name] = root.find(f"{{{NS}}}dimension").attrib["ref"]
            cells = {}
            for cell in root.findall(f".//{{{NS}}}c"):
                reference = cell.attrib["r"]
                kind = cell.attrib.get("t")
                value_node = cell.find(f"{{{NS}}}v")
                if kind == "inlineStr":
                    value = "".join(node.text or "" for node in cell.iter(f"{{{NS}}}t"))
                elif value_node is None:
                    value = None
                elif kind == "s":
                    value = shared[int(value_node.text)]
                elif kind in {"str", "e"}:
                    value = value_node.text or ""
                else:
                    number = float(value_node.text)
                    value = int(number) if number.is_integer() else number
                cells[reference] = value
            outputs[name] = cells
        return dimensions, outputs


def period_arrays(cells: dict[str, object]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    tn, wn, tf, wf = [], [], [], []
    for row in range(3, 23):
        time = float(cells[f"B{row}"])
        if f"C{row}" in cells:
            tn.append(time)
            wn.append((2.0 * math.pi / float(cells[f"C{row}"])) ** 2)
        if f"D{row}" in cells:
            tf.append(time)
            wf.append((2.0 * math.pi / float(cells[f"D{row}"])) ** 2)
    return tuple(np.asarray(item, dtype=float) for item in (tn, wn, tf, wf))


def lerp(t0: float, y0: float, t1: float, y1: float, time: float) -> float:
    return y0 + (y1 - y0) * (time - t0) / (t1 - t0)


def aba(cells: dict[str, object]) -> float:
    tn, wn, tf, wf = period_arrays(cells)
    values = []
    for index in range(9):
        values.append(lerp(tn[index], wn[index], tn[index + 1], wn[index + 1], tf[index]) - wf[index])
    for index in range(1, 10):
        values.append(wn[index] - lerp(tf[index - 1], wf[index - 1], tf[index], wf[index], tn[index]))
    return float(np.mean(values))


def quadratic(cells: dict[str, object]) -> float:
    tn, wn, tf, wf = period_arrays(cells)
    time = np.concatenate((tn, tf))
    response = np.concatenate((wn, wf))
    flag = np.concatenate((np.ones(10), np.zeros(10)))
    order = np.argsort(time)
    time, response, flag = time[order], response[order], flag[order]
    scaled = (time - np.mean(time)) / np.std(time)
    design = np.column_stack((np.ones(20), scaled, scaled**2, flag))
    beta = np.linalg.lstsq(design, response, rcond=None)[0]
    return float(beta[-1])


def two_tone(cells: dict[str, object]) -> tuple[int, float, float, float, float]:
    time = np.asarray([float(cells[f"B{row}"]) for row in range(3, 7203)], dtype=float)
    response = np.asarray([float(cells[f"C{row}"]) for row in range(3, 7203)], dtype=float)
    scaled = (time - np.mean(time)) / (time[-1] - time[0])
    source_frequency = 1.0 / 600.0
    best = None
    for index in range(6001):
        background_frequency = 0.0005 + index * 1e-7
        design = np.column_stack((
            np.ones(7200), scaled,
            np.sin(2.0 * math.pi * source_frequency * time),
            np.cos(2.0 * math.pi * source_frequency * time),
            np.sin(2.0 * math.pi * background_frequency * time),
            np.cos(2.0 * math.pi * background_frequency * time),
        ))
        beta = np.linalg.lstsq(design, response, rcond=None)[0]
        residual = response - design @ beta
        rss = float(residual @ residual)
        if best is None or rss < best[0]:
            best = rss, background_frequency, beta
    assert best is not None
    rss, frequency, beta = best
    return (
        7200,
        frequency,
        float(math.hypot(beta[2], beta[3])),
        float(math.hypot(beta[4], beta[5])),
        math.sqrt(rss / 7200),
    )


def main() -> None:
    checks = []

    def ck(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    result = json.loads((HERE / "RESULT.json").read_text(encoding="utf-8"))
    custody = json.loads((HERE / "SOURCE_CUSTODY.json").read_text(encoding="utf-8"))

    ck("schema", result["schema"] == "WAC_HUST_2018_DUAL_METHOD_G_FORWARD_V001")
    ck("status", result["status"].endswith("FULL_GC16_NOT_EXECUTABLE__NO_NEW_G"))
    ck("doi", result["source"]["doi"] == custody["article"]["doi"] == "10.1038/s41586-018-0431-5")
    ck("article url", custody["article"]["article_url"] == "https://www.nature.com/articles/s41586-018-0431-5")
    ck("release count", len(custody["official_release"]) == 7)
    ck("release ids", [item["id"] for item in custody["official_release"]] == [f"MOESM{i}" for i in range(1, 8)])
    ck("pin counts", custody["pinned_file_count"] == 2 and custody["unpinned_file_count"] == 5)
    ck("pin flags", [item["id"] for item in custody["official_release"] if item["pinned_in_lane"]] == ["MOESM1", "MOESM3"])
    ck("official urls", all(item["url"].startswith("https://media.springernature.com/original/") for item in custody["official_release"]))
    ck("release hashes unique", len({item["sha256"] for item in custody["official_release"]}) == 7)
    ck("release sizes positive", all(item["bytes"] > 10000 for item in custody["official_release"]))

    for relative, expected in EXPECTED_HASHES.items():
        path = (HERE / relative).resolve()
        ck(f"dependency exists {relative}", path.is_file() and not path.is_symlink())
        ck(f"dependency hash {relative}", digest(path) == expected)
    ck("source sizes", PDF.stat().st_size == 2711453 and XLSX.stat().st_size == 4169004)

    dimensions, sheets = workbook_cells(XLSX)
    ck("dimensions a", dimensions["a"] == "B2:D22")
    ck("dimensions b", dimensions["b"] == "B2:D22")
    ck("dimensions c", dimensions["c"] == "B1:F10")
    ck("dimensions e", dimensions["e"] == "B2:C10001")
    ck("dimensions f", dimensions["f"] == "A1:I39")
    ck("period headers", sheets["a"]["C2"] == sheets["b"]["C2"] == "Period at near position (s)")
    ck("acceleration header", sheets["e"]["C2"] == "Acceleration (nrad/s2)")

    aba_present = aba(sheets["a"])
    aba_background = aba(sheets["b"])
    aba_net = aba_present - aba_background
    ck("A-B-A present", math.isclose(aba_present, 1.6621980172599122e-6, rel_tol=0.0, abs_tol=3e-19))
    ck("A-B-A background", math.isclose(aba_background, -4.964938724048171e-10, rel_tol=0.0, abs_tol=3e-19))
    ck("A-B-A net", math.isclose(aba_net, result["ToS_figure_level_response"]["A_B_A_background_subtracted_delta_omega2_s-2"], rel_tol=0.0, abs_tol=3e-19))
    ck("A-B-A published sigma", math.isclose((aba_net - 1.662699e-6) / 18e-12, -0.24938153793, abs_tol=1e-10))

    quadratic_net = quadratic(sheets["a"]) - quadratic(sheets["b"])
    ck("quadratic net", math.isclose(quadratic_net, 1.6626989120180067e-6, rel_tol=0.0, abs_tol=3e-19))
    ck("quadratic result custody", math.isclose(quadratic_net, result["ToS_figure_level_response"]["common_quadratic_background_subtracted_delta_omega2_s-2"], rel_tol=0.0, abs_tol=3e-19))
    ck("ToS response quotient", math.isclose(aba_net / 24912.12, 6.67423933062428e-11, rel_tol=2e-15))
    ck("ToS seven rows", len(result["ToS_seven_processed_coefficient_forwards"]) == 7)
    residuals = [row["unowned_full_bracket_residual_ppm"] for row in result["ToS_seven_processed_coefficient_forwards"]]
    ck("ToS residual signs", all(value < 0.0 for value in residuals))
    ck("ToS residual range", min(residuals) < -8.37 and max(residuals) < -5.46 and max(residuals) > -5.47)

    ck("AAF row count", len(result["AAF_three_processed_coefficient_forwards"]) == 3)
    expected_aaf = {
        "AAF-I": (6.674532777558953e-11, -0.1831500218152371, 455.3960973721072),
        "AAF-II": (6.674375348038621e-11, 0.05214549991805484, 455.3960973721072),
        "AAF-III": (6.674535082204988e-11, 0.01231621205022293, 25.739273134872676),
    }
    for row in result["AAF_three_processed_coefficient_forwards"]:
        correction = 1.0 + row["K"] / row["K_m"] * row["I_m"] / row["I"]
        calculated = row["alpha_nrad_s-2"] * 1e-9 / row["P_abs_kg_m-3"] * correction
        expected_g, expected_ppm, expected_correction = expected_aaf[row["id"]]
        ck(f"AAF formula {row['id']}", math.isclose(calculated, expected_g, rel_tol=2e-15))
        ck(f"AAF result {row['id']}", math.isclose(row["recomputed_G_SI"], expected_g, rel_tol=2e-15))
        ck(f"AAF rounding {row['id']}", math.isclose(row["recomputed_minus_published_relative_ppm"], expected_ppm, abs_tol=1e-9))
        ck(f"AAF correction {row['id']}", math.isclose(row["magnetic_damper_ppm"], expected_correction, abs_tol=1e-9))
    ck("AAF rounding ceiling", max(abs(row["recomputed_minus_published_relative_ppm"]) for row in result["AAF_three_processed_coefficient_forwards"]) < 0.2)

    groups = {
        "AAF-I": ("B", range(3, 7), 4, 6.67453375),
        "AAF-II": ("D", range(9, 19), 10, 6.674375),
        "AAF-III": ("F", range(21, 36), 15, 6.6745348),
    }
    for name, (column, rows, count, expected_mean) in groups.items():
        values = [float(sheets["f"][f"{column}{row}"]) for row in rows]
        ck(f"campaign count {name}", len(values) == count)
        ck(f"campaign mean {name}", math.isclose(float(np.mean(values)), expected_mean, abs_tol=5e-15))
        stored = result["already_derived_outcome_rows"]["AAF_29_already_derived_G_rows_by_campaign"][name]
        ck(f"campaign stored {name}", stored["count"] == count and math.isclose(stored["arithmetic_mean_G_e-11"], expected_mean, abs_tol=5e-15))
    ck("29 derived rows", sum(item[2] for item in groups.values()) == 29)
    ck("seven ToS cells", [sheets["c"][f"C{row}"] for row in range(3, 10)] == [6.674154, 6.674222, 6.674237, 6.674274, 6.674266, 6.674017, 6.674105])
    ck("combined cells", sheets["c"]["E10"] == 6.674184 and sheets["c"]["F10"] == 0.000078 and sheets["f"]["H39"] == 6.674484 and sheets["f"]["I39"] == 0.000078)

    count, frequency, source_amplitude, background_amplitude, residual_rms = two_tone(sheets["e"])
    stream = result["AAF_figure_level_acceleration_stream"]
    ck("stream count", count == stream["caption_scored_two_hour_sample_count"] == 7200)
    ck("stream available", stream["workbook_available_sample_count"] == 9999)
    ck("stream frequency", frequency == stream["two_tone_best_background_frequency_Hz"] == 0.0007397)
    ck("stream source", math.isclose(source_amplitude, 461.99346479501804, abs_tol=1e-9))
    ck("stream background", math.isclose(background_amplitude, 77.85005910756898, abs_tol=1e-9))
    ck("stream residual", math.isclose(residual_rms, 15.710326418514784, abs_tol=1e-9))
    ck("stream custody", math.isclose(source_amplitude, stream["two_tone_source_amplitude_nrad_s-2"], abs_tol=1e-9) and math.isclose(background_amplitude, stream["two_tone_background_amplitude_nrad_s-2"], abs_tol=1e-9))
    ck("stream not campaign", "NOT_CAMPAIGN_ALPHA" in stream["claim"])

    cross = result["cross_method_source_model_stress"]
    ck("cross difference", math.isclose(cross["AAF_minus_ToS_SI"], 3e-15, rel_tol=0.0, abs_tol=2e-28))
    ck("cross ppm", math.isclose(cross["difference_relative_to_midpoint_ppm"], 44.94830495447016, abs_tol=1e-10))
    ck("cross z", math.isclose(cross["z_if_and_only_if_cross_method_covariance_is_zero"], 2.719641466102003, abs_tol=1e-12))
    ck("cross covariance ceiling", cross["cross_method_covariance_released"] is False)
    ck("cross interpretation", "NOT_NEW_PHYSICS" in cross["interpretation"] and "NOT_A_NEW_G_ESTIMATE" in cross["interpretation"])

    inventory = result["public_field_inventory"]
    ck("raw-like inventory", inventory["raw_like_or_intermediate"]["ToS_source_present_period_summaries"] == 20 and inventory["raw_like_or_intermediate"]["AAF_acceleration_1_second_samples"] == 9999)
    ck("derived inventory", inventory["already_derived"]["ToS_G_rows"] == 7 and inventory["already_derived"]["AAF_G_rows"] == 29)
    gc = result["GC16_execution"]
    ck("GC16 false", gc["full_real_apparatus_execution"] is False)
    ck("closure list", len(gc["closure_gained"]) == 7)
    ck("missing list", len(gc["missing_fields"]) == 9)
    ck("geometry missing", any("mass-coordinate-density" in item for item in gc["missing_fields"]))
    ck("covariance missing", any("covariance" in item for item in gc["missing_fields"]))
    ck("stress missing", any("stress/source ledger" in item for item in gc["missing_fields"]))
    ck("identifiability ceiling", gc["identifiability"].endswith("NOT_SEPARATELY_REIDENTIFIED"))
    ceilings = result["claim_ceiling"]
    ck("all claim booleans false", all(value is False for key, value in ceilings.items() if key != "strongest_claim"))
    ck("strongest claim", "NOT_THE_CALIBRATED_FULL_APPARATUS_G_CROSSCHECK" in ceilings["strongest_claim"])

    theorem = (HERE / "THEOREM.md").read_text(encoding="utf-8")
    readiness = (HERE / "READINESS.md").read_text(encoding="utf-8")
    audit = (HERE / "SELF_AUDIT.md").read_text(encoding="utf-8")
    analyzer = (HERE / "analyze_hust_dual_method_g_forward.py").read_text(encoding="utf-8")
    text_payloads = {
        path.relative_to(HERE).as_posix(): path.read_bytes()
        for path in HERE.rglob("*")
        if path.is_file() and path.suffix.lower() not in {".pdf", ".xlsx"}
    }
    ck("no carriage-return bytes", all(b"\r" not in payload for payload in text_payloads.values()))
    ck("no backspace bytes", all(b"\b" not in payload for payload in text_payloads.values()))
    ck("no form-feed bytes", all(b"\f" not in payload for payload in text_payloads.values()))
    ck("operator/source statement", "operator/stiffness" in theorem and "source/forcing" in theorem)
    ck("theorem A-B-A", "1.6626945111323172" in theorem)
    ck("theorem AAF", "6.674532777558953" in theorem)
    ck("theorem cross ceiling", "If and only if" in theorem and "source-model/systematics stress" in theorem)
    ck(
        "delta omega squared uses inverse-second-squared",
        theorem.count(r"\Delta\omega^2_{\mathrm{ABA}}") == 1
        and theorem.count(r"\Delta\omega^2_{\mathrm{quad}}") == 1
        and theorem.count(r"\mathrm{s}^{-2}") >= 3,
    )
    ck("delta omega squared not linear acceleration", r"\mathrm{m\,s^{-2}}" not in theorem)
    ck("readiness ranges", "`a!B3:D22`" in readiness and "`e!B3:C7202`" in readiness and "`f!H39:I39`" in readiness)
    ck("readiness verdict", "FULL_GC16_NOT_READY" in readiness)
    ck("audit verdict", "ACCEPT_WITH_PROCESSED_COEFFICIENT_AND_FIGURE_LEVEL_CEILINGS" in audit)
    ck("no MOESM4 analyzer", "MOESM4" not in analyzer)
    ck("no accepted comparison input", "accepted numerical value" not in analyzer.lower())

    run = subprocess.run(
        [sys.executable, str(HERE / "analyze_hust_dual_method_g_forward.py"), "--check"],
        cwd=ROOT,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    ck("analyzer check exit", run.returncode == 0)
    ck("analyzer check text", run.stdout.strip() == "HUST_DUAL_METHOD_ANALYZER_CHECK: PASS")

    manifest_lines = (HERE / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines()
    manifest = {}
    for line in manifest_lines:
        expected, relative = line.split("  ", 1)
        manifest[relative] = expected
    ck("manifest file set", set(manifest) == MANIFEST_FILES)
    ck("manifest hashes", all(digest(HERE / relative) == expected for relative, expected in manifest.items()))
    ck("manifest no symlinks", all(not (HERE / relative).is_symlink() for relative in manifest))
    seal_line = (HERE / "LANE_SEAL.sha256").read_text(encoding="utf-8").strip()
    seal_hash, seal_target = seal_line.split("  ", 1)
    ck("seal target", seal_target == "MANIFEST.sha256")
    ck("seal hash", seal_hash == digest(HERE / "MANIFEST.sha256"))
    lane_files = {
        path.relative_to(HERE).as_posix()
        for path in HERE.rglob("*")
        if path.is_file()
    }
    ck(
        "lane file set",
        lane_files
        == MANIFEST_FILES
        | {"MANIFEST.sha256", "LANE_SEAL.sha256", "VERIFICATION.txt", "HOSTILE_VERIFICATION.txt"},
    )

    text = (
        "HUST_DUAL_METHOD_G_FORWARD_CHECK: PASS\n"
        f"Checks: {len(checks)}/{len(checks)}\n"
        "Pinned official sources: 2; release inventory: 7\n"
        "ToS: one A-B-A background-subtracted figure-level response recovered\n"
        "AAF: three processed source-response forwards reproduced within 0.2 ppm\n"
        "AAF stream: 461.993 nrad/s2 source and 77.850 nrad/s2 background separated\n"
        "Cross-method stress: 44.948 ppm; 2.720 only under zero cross covariance\n"
        "Full GC16, new G, RGRL/GFT confirmation, and lineage charge: not claimed\n"
    )
    (HERE / "VERIFICATION.txt").write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
