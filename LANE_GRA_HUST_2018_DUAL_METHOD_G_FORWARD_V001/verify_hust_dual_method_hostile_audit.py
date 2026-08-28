#!/usr/bin/env python3
"""Hostile replay of the HUST-2018 bounded public-data forward.

This verifier does not import the lane analyzer or its first verifier.  It
derives the scored quantities from the pinned primary files before comparing
them with RESULT.json.  The A-B-A replay uses a merged alternating sequence,
the quadratic replay uses QR, and the two-tone scan uses a residualized
two-column projection rather than the builder's repeated six-column solve.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import posixpath
import shutil
import tempfile
import zipfile
import xml.etree.ElementTree as ET

import numpy as np


HERE = Path(__file__).resolve().parent
PDF = HERE / "SOURCE/41586_2018_431_MOESM1_ESM.pdf"
XLSX = HERE / "SOURCE/41586_2018_431_MOESM3_ESM.xlsx"
NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG = "http://schemas.openxmlformats.org/package/2006/relationships"
BASE = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fs41586-018-0431-5/MediaObjects/"
)
RELEASE = {
    1: ("pdf", 2711453, "5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb"),
    2: ("xlsx", 1152664, "9e419d1150b6f7897a1352cee50a9721fed48254279e13fe2847200e677547ba"),
    3: ("xlsx", 4169004, "331ea6ee8a6a2558d4f6a8ffa4cd52e4c5d52d47f2e2a697889d0d3c8a3ad27a"),
    4: ("xlsx", 10330, "dda202fa672f2524c0cd33eb38c74c67f790ef24842b752afecd46311d7cf2aa"),
    5: ("xlsx", 36733, "08af72bca3f7861ac63c6e7695687cd71d3a3b9b631187ce9ce1417e09c66da3"),
    6: ("xlsx", 15013, "8d07c2cdfa92d5ff14a49fd8e0b876089a0df920996638a5c5e6ae5ac309aec7"),
    7: ("xlsx", 15635, "20081480399f3502894f1658236dd31841ebcbf49fb661aafea7655a66512903"),
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


def workbook(path: Path):
    with zipfile.ZipFile(path) as archive:
        if archive.testzip() is not None:
            raise AssertionError("xlsx CRC")
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = [
                "".join(node.text or "" for node in item.iter(f"{{{NS}}}t"))
                for item in root.findall(f"{{{NS}}}si")
            ]
        book = ET.fromstring(archive.read("xl/workbook.xml"))
        rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {
            node.attrib["Id"]: node.attrib["Target"]
            for node in rels.findall(f"{{{PKG}}}Relationship")
        }
        paths = {}
        for sheet in book.findall(f".//{{{NS}}}sheet"):
            target = targets[sheet.attrib[f"{{{REL}}}id"]]
            paths[sheet.attrib["name"]] = (
                target.lstrip("/") if target.startswith("/")
                else posixpath.normpath(posixpath.join("xl", target))
            )
        dimensions, outputs, formula_refs = {}, {}, []
        for name in ("a", "b", "c", "d", "e", "f"):
            root = ET.fromstring(archive.read(paths[name]))
            dimensions[name] = root.find(f"{{{NS}}}dimension").attrib["ref"]
            cells = {}
            for cell in root.findall(f".//{{{NS}}}c"):
                ref = cell.attrib["r"]
                if cell.find(f"{{{NS}}}f") is not None:
                    formula_refs.append(f"{name}!{ref}")
                kind = cell.attrib.get("t")
                value_node = cell.find(f"{{{NS}}}v")
                if kind == "inlineStr":
                    value = "".join(n.text or "" for n in cell.iter(f"{{{NS}}}t"))
                elif value_node is None:
                    value = None
                elif kind == "s":
                    value = shared[int(value_node.text)]
                elif kind in {"str", "e"}:
                    value = value_node.text or ""
                else:
                    number = float(value_node.text)
                    value = int(number) if number.is_integer() else number
                cells[ref] = value
            outputs[name] = cells
        return dimensions, outputs, formula_refs


def panel_observations(cells):
    rows = []
    for row in range(3, 23):
        time = float(cells[f"B{row}"])
        if cells.get(f"C{row}") is not None:
            rows.append((time, "near", (2.0 * math.pi / float(cells[f"C{row}"])) ** 2))
        if cells.get(f"D{row}") is not None:
            rows.append((time, "far", (2.0 * math.pi / float(cells[f"D{row}"])) ** 2))
    return sorted(rows)


def merged_aba(cells):
    rows = panel_observations(cells)
    values = []
    for index in range(1, len(rows) - 1):
        left, middle, right = rows[index - 1], rows[index], rows[index + 1]
        if not (left[1] == right[1] and left[1] != middle[1]):
            raise AssertionError("nonalternating A-B-A source sequence")
        predicted = left[2] + (right[2] - left[2]) * (
            (middle[0] - left[0]) / (right[0] - left[0])
        )
        values.append(middle[2] - predicted if middle[1] == "near" else predicted - middle[2])
    return np.asarray(values, dtype=float)


def qr_quadratic(cells):
    rows = panel_observations(cells)
    time = np.asarray([row[0] for row in rows], dtype=float)
    response = np.asarray([row[2] for row in rows], dtype=float)
    near = np.asarray([row[1] == "near" for row in rows], dtype=float)
    scaled = (time - time.mean()) / time.std()
    design = np.column_stack((np.ones(20), scaled, scaled * scaled, near))
    q, r = np.linalg.qr(design, mode="reduced")
    beta = np.linalg.solve(r, q.T @ response)
    return float(beta[-1]), float(np.sqrt(np.mean((response - design @ beta) ** 2)))


def qr_coefficients(design, response):
    q, r = np.linalg.qr(design, mode="reduced")
    return np.linalg.solve(r, q.T @ response)


def residualized_two_tone(cells):
    time = np.asarray([float(cells[f"B{row}"]) for row in range(3, 7203)])
    response = np.asarray([float(cells[f"C{row}"]) for row in range(3, 7203)])
    scaled = (time - time.mean()) / (time[-1] - time[0])
    source_frequency = 1.0 / 600.0
    base = np.column_stack((
        np.ones(7200), scaled,
        np.sin(2 * math.pi * source_frequency * time),
        np.cos(2 * math.pi * source_frequency * time),
    ))
    q, _ = np.linalg.qr(base, mode="reduced")
    y_residual = response - q @ (q.T @ response)
    base_rss = float(y_residual @ y_residual)
    best = None
    for index in range(6001):
        frequency = 0.0005 + index * 1e-7
        pair = np.column_stack((
            np.sin(2 * math.pi * frequency * time),
            np.cos(2 * math.pi * frequency * time),
        ))
        pair -= q @ (q.T @ pair)
        gram = pair.T @ pair
        rhs = pair.T @ y_residual
        rss = base_rss - float(rhs @ np.linalg.solve(gram, rhs))
        if best is None or rss < best[0]:
            best = (rss, frequency)
    assert best is not None
    frequency = best[1]
    design = np.column_stack((
        base,
        np.sin(2 * math.pi * frequency * time),
        np.cos(2 * math.pi * frequency * time),
    ))
    beta = qr_coefficients(design, response)
    residual = response - design @ beta
    return (
        frequency,
        math.hypot(beta[2], beta[3]),
        math.hypot(beta[4], beta[5]),
        math.sqrt(float(residual @ residual) / 7200),
    )


def read_manifest(root: Path):
    values = {}
    for line in (root / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        values[relative] = expected
    return values


def manifest_valid(root: Path) -> bool:
    values = read_manifest(root)
    return set(values) == MANIFEST_FILES and all(
        (root / relative).is_file()
        and not (root / relative).is_symlink()
        and digest(root / relative) == expected
        for relative, expected in values.items()
    )


def seal_valid(root: Path) -> bool:
    expected, relative = (root / "LANE_SEAL.sha256").read_text(encoding="utf-8").strip().split("  ", 1)
    return relative == "MANIFEST.sha256" and expected == digest(root / relative)


def main() -> None:
    checks = []

    def ck(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    result = json.loads((HERE / "RESULT.json").read_text(encoding="utf-8"))
    custody = json.loads((HERE / "SOURCE_CUSTODY.json").read_text(encoding="utf-8"))

    ck("release count", len(custody["official_release"]) == 7)
    for item, index in zip(custody["official_release"], range(1, 8)):
        extension, size, expected_hash = RELEASE[index]
        expected_url = BASE + f"41586_2018_431_MOESM{index}_ESM.{extension}"
        ck(f"release id {index}", item["id"] == f"MOESM{index}")
        ck(f"release url {index}", item["url"] == expected_url)
        ck(f"release size {index}", item["bytes"] == size)
        ck(f"release hash {index}", item["sha256"] == expected_hash)
    ck("pinned pdf bytes", PDF.stat().st_size == RELEASE[1][1] and digest(PDF) == RELEASE[1][2])
    ck("pinned xlsx bytes", XLSX.stat().st_size == RELEASE[3][1] and digest(XLSX) == RELEASE[3][2])

    dimensions, sheets, formula_refs = workbook(XLSX)
    ck("workbook dimensions", dimensions == {
        "a": "B2:D22", "b": "B2:D22", "c": "B1:F10",
        "d": "B2:G129604", "e": "B2:C10001", "f": "A1:I39",
    })
    ck("source ranges contain no formulas", formula_refs == [])
    ck("period units", sheets["a"]["B2"] == "Time (day)" and sheets["a"]["C2"] == "Period at near position (s)")
    ck("acceleration units", sheets["e"]["B2"] == "Time (s)" and sheets["e"]["C2"] == "Acceleration (nrad/s2)")
    ck("one-second stream", [sheets["e"][f"B{row}"] for row in range(3, 10002)] == list(range(1, 10000)))

    present = merged_aba(sheets["a"])
    background = merged_aba(sheets["b"])
    aba_net = float(present.mean() - background.mean())
    ck("18 present triples", present.size == 18)
    ck("18 background triples", background.size == 18)
    ck("A-B-A present", math.isclose(float(present.mean()), 1.6621980172599122e-6, abs_tol=3e-19))
    ck("A-B-A background", math.isclose(float(background.mean()), -4.964938724048171e-10, abs_tol=3e-19))
    ck("A-B-A net", math.isclose(aba_net, 1.6626945111323172e-6, abs_tol=3e-19))
    ck("A-B-A result custody", math.isclose(aba_net, result["ToS_figure_level_response"]["A_B_A_background_subtracted_delta_omega2_s-2"], abs_tol=3e-19))

    quad_present, rms_present = qr_quadratic(sheets["a"])
    quad_background, rms_background = qr_quadratic(sheets["b"])
    quad_net = quad_present - quad_background
    ck("quadratic present QR", math.isclose(quad_present, 1.6621944646749398e-6, abs_tol=4e-19))
    ck("quadratic background QR", math.isclose(quad_background, -5.044473430670186e-10, abs_tol=4e-19))
    ck("quadratic net QR", math.isclose(quad_net, 1.6626989120180067e-6, abs_tol=5e-19))
    ck("quadratic RMS", math.isclose(rms_present, 4.2671660643210314e-11, abs_tol=2e-22) and math.isclose(rms_background, 4.5740930903667815e-11, abs_tol=2e-22))

    mechanics = {
        "fiber_1": (2.18e-5, 1.2e-5, 4.7705e-5, 12.2e-9),
        "fiber_2": (2.18e-5, 1.2e-5, 4.7706e-5, 47.4e-9),
        "fiber_3": (2.18e-5, 1.2e-5, 4.7705e-5, 10.1e-9),
        "fiber_4": (2.21e-5, 1.4e-5, 4.6477e-5, 10.6e-9),
    }
    table2 = [
        ("fiber_1", 24912.86, 1.662732, 6.674154),
        ("fiber_1", 24912.12, 1.662699, 6.674222),
        ("fiber_2", 24912.15, 1.662698, 6.674237),
        ("fiber_3", 24911.70, 1.662684, 6.674274),
        ("fiber_3", 24911.72, 1.662683, 6.674266),
        ("fiber_4", 25003.05, 1.668719, 6.674017),
        ("fiber_4", 25002.95, 1.668734, 6.674105),
    ]
    inferred = []
    for stored, (fiber, coefficient, delta, published) in zip(result["ToS_seven_processed_coefficient_forwards"], table2):
        im, km, inertia, spring = mechanics[fiber]
        raw = delta * 1e-6 / coefficient
        magnetic = im * spring**2 / (inertia * km**2)
        residual = (published * 1e-11 / raw - (1 + magnetic)) * 1e6
        inferred.append(residual)
        ck(f"ToS raw {fiber} {delta}", math.isclose(stored["raw_response_over_source_ratio_SI"], raw, rel_tol=2e-15))
        ck(f"ToS magnetic {fiber} {delta}", math.isclose(stored["magnetic_damper_factor_minus_one"], magnetic, rel_tol=2e-15))
        ck(f"ToS bracket {fiber} {delta}", math.isclose(stored["unowned_full_bracket_residual_ppm"], residual, abs_tol=2e-9))
    ck("ToS bracket ceiling", min(inferred) < -8.37 and -5.47 < max(inferred) < -5.46)

    table3 = [
        ("AAF-I", 6926.352, 462.0912, 2.401e-5, 1.199e-5, 2.776e-5, 6.313e-9, 6.674534),
        ("AAF-II", 6926.334, 462.0791, 2.401e-5, 1.199e-5, 2.776e-5, 6.313e-9, 6.674375),
        ("AAF-III", 6926.415, 462.2941, 2.404e-5, 21.24e-5, 2.776e-5, 6.313e-9, 6.674535),
    ]
    for stored, row in zip(result["AAF_three_processed_coefficient_forwards"], table3):
        name, source, alpha, im, km, inertia, spring, published = row
        correction = 1 + spring / km * im / inertia
        calculated = alpha * 1e-9 / source * correction
        ppm = (calculated / (published * 1e-11) - 1) * 1e6
        ck(f"AAF identity {name}", stored["id"] == name)
        ck(f"AAF forward {name}", math.isclose(stored["recomputed_G_SI"], calculated, rel_tol=2e-15))
        ck(f"AAF rounding {name}", math.isclose(stored["recomputed_minus_published_relative_ppm"], ppm, abs_tol=2e-9) and abs(ppm) < 0.2)

    groups = {
        "AAF-I": ("B", range(3, 7), 6.67453375),
        "AAF-II": ("D", range(9, 19), 6.674375),
        "AAF-III": ("F", range(21, 36), 6.6745348),
    }
    for name, (column, rows, expected) in groups.items():
        values = [float(sheets["f"][f"{column}{row}"]) for row in rows]
        ck(f"derived campaign mean {name}", math.isclose(float(np.mean(values)), expected, abs_tol=5e-15))
    ck("derived row classification", result["public_field_inventory"]["already_derived"]["AAF_G_rows"] == 29)

    frequency, source_amplitude, background_amplitude, residual_rms = residualized_two_tone(sheets["e"])
    stream = result["AAF_figure_level_acceleration_stream"]
    ck("two-tone frequency", frequency == 0.0007397 == stream["two_tone_best_background_frequency_Hz"])
    ck("two-tone source", math.isclose(source_amplitude, 461.99346479501804, abs_tol=2e-9))
    ck("two-tone background", math.isclose(background_amplitude, 77.85005910756898, abs_tol=2e-9))
    ck("two-tone residual", math.isclose(residual_rms, 15.710326418514784, abs_tol=2e-9))
    ck("two-tone ceiling", "NOT_CAMPAIGN_ALPHA_OR_G_EXTRACTION" in stream["claim"])

    tos, aaf, uncertainty = 6.674184e-11, 6.674484e-11, 0.000078e-11
    difference = aaf - tos
    ppm = difference / ((aaf + tos) / 2) * 1e6
    z_zero_covariance = difference / math.sqrt(2 * uncertainty**2)
    cross = result["cross_method_source_model_stress"]
    ck("cross delta", math.isclose(difference, 3e-15, abs_tol=2e-28))
    ck("cross ppm", math.isclose(ppm, 44.94830495447016, abs_tol=1e-10) and math.isclose(cross["difference_relative_to_midpoint_ppm"], ppm, abs_tol=1e-10))
    ck("cross conditional z", math.isclose(z_zero_covariance, 2.719641466102003, abs_tol=1e-12) and cross["cross_method_covariance_released"] is False)

    gc = result["GC16_execution"]
    ck("GC16 not executed", gc["full_real_apparatus_execution"] is False)
    ck("GC16 missing geometry", any("mass-coordinate-density" in item for item in gc["missing_fields"]))
    ck("GC16 missing transfer", any("frequency-dependent" in item for item in gc["missing_fields"]))
    ck("GC16 missing covariance", any("covariance matrices" in item for item in gc["missing_fields"]))
    ck("GC16 missing stress", any("stress/source ledger" in item for item in gc["missing_fields"]))
    ck("no new-G claims", all(value is False for key, value in result["claim_ceiling"].items() if key != "strongest_claim"))

    theorem = (HERE / "THEOREM.md").read_text(encoding="utf-8")
    readiness = (HERE / "READINESS.md").read_text(encoding="utf-8")
    analyzer = (HERE / "analyze_hust_dual_method_g_forward.py").read_text(encoding="utf-8")
    ck("panel-wise wording", "Independent panel-wise regressions" in theorem and "its 20 near/far summaries" in theorem)
    ck("range custody", all(item in readiness for item in ("`a!B3:D22`", "`b!B3:D22`", "`e!B3:C7202`", "`f!H39:I39`")))
    ck("unit custody", r"\mathrm{s}^{-2}" in theorem and r"\mathrm{kg\,m^{-3}}" in theorem and "nrad" in theorem)
    ck("accepted G not loaded", "MOESM4" not in analyzer and "accepted_G_used_as_input\": False" in analyzer)
    text_files = [path for path in HERE.rglob("*") if path.is_file() and path.suffix.lower() not in {".pdf", ".xlsx"}]
    ck("control bytes", all(not any(byte in path.read_bytes() for byte in (0, 8, 12, 13)) for path in text_files))

    ck("manifest", manifest_valid(HERE))
    ck("seal", seal_valid(HERE))
    with tempfile.TemporaryDirectory(prefix="hust-hostile-") as temporary:
        copy = Path(temporary) / "lane"
        copy.mkdir()
        for relative in MANIFEST_FILES | {"MANIFEST.sha256", "LANE_SEAL.sha256"}:
            target = copy / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(HERE / relative, target)
        (copy / "THEOREM.md").write_bytes((copy / "THEOREM.md").read_bytes() + b"\nTAMPER\n")
        ck("content tamper rejected", not manifest_valid(copy))
        shutil.copy2(HERE / "THEOREM.md", copy / "THEOREM.md")
        (copy / "MANIFEST.sha256").write_bytes((copy / "MANIFEST.sha256").read_bytes() + b"\n")
        ck("manifest tamper rejected", not seal_valid(copy))

    text = (
        "HUST_DUAL_METHOD_HOSTILE_AUDIT: PASS\n"
        f"Checks: {len(checks)}/{len(checks)}\n"
        "Official release: 7/7 URLs, byte counts, and SHA-256 values independently checked\n"
        "ToS: merged-sequence A-B-A and independent QR panel fits reproduced\n"
        "AAF: three primary-table forwards and a residualized two-tone replay reproduced\n"
        "Cross-method stress: 44.948 ppm; 2.720 conditional on zero cross covariance\n"
        "Ceiling: processed/figure-level forward only; no new G and no GC16/RGRL/GFT test\n"
        "Tamper tests: content manifest and lane seal both reject mutation\n"
    )
    (HERE / "HOSTILE_VERIFICATION.txt").write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
