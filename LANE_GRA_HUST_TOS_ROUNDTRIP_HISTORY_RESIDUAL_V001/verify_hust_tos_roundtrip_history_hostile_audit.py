#!/usr/bin/env python3
"""Independent hostile audit for GRA-HUST18-TOS-RTHR-V001.

This executable does not import the lane analyzer or its first verifier.  It
re-opens the pinned OOXML workbook, rebuilds the endpoint histories, exact
integer weights, ranks, Gram matrices, telescoping identities, timing ledger,
and every published mean.  It then screens the claim ceiling.  ``--preseal``
is used only to create the deterministic transcript before the manifest is
rebuilt; normal execution also verifies the final manifest and outer seal.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import posixpath
import xml.etree.ElementTree as ET
import zipfile
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORKBOOK_REL = (
    "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/SOURCE/"
    "41586_2018_431_MOESM3_ESM.xlsx"
)
WORKBOOK = (HERE / WORKBOOK_REL).resolve()
WORKBOOK_SHA256 = "331ea6ee8a6a2558d4f6a8ffa4cd52e4c5d52d47f2e2a697889d0d3c8a3ad27a"

MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
OFFICE_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PACKAGE_REL = "http://schemas.openxmlformats.org/package/2006/relationships"

MANIFEST_MEMBERS = {
    "DEPENDENCIES.sha256",
    "HOSTILE_VERIFICATION.txt",
    "INDEPENDENT_HOSTILE_AUDIT.md",
    "README.md",
    "RESULT.json",
    "RESULT.md",
    "SELF_AUDIT.md",
    "SOURCE_CUSTODY.json",
    "THEOREM.md",
    "VERIFICATION.txt",
    "analyze_hust_tos_roundtrip_history.py",
    "verify_hust_tos_roundtrip_history.py",
    "verify_hust_tos_roundtrip_history_hostile_audit.py",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class Audit:
    def __init__(self) -> None:
        self.count = 0

    def require(self, condition: bool, label: str) -> None:
        if not condition:
            raise AssertionError(label)
        self.count += 1

    def equal(self, actual: object, expected: object, label: str) -> None:
        self.require(actual == expected, f"{label}: {actual!r} != {expected!r}")

    def close(self, actual: float, expected: float, label: str) -> None:
        self.require(
            math.isclose(float(actual), float(expected), rel_tol=2e-13, abs_tol=2e-24),
            f"{label}: {actual!r} != {expected!r}",
        )

    def close_list(self, actual: list[float], expected: list[float], label: str) -> None:
        self.equal(len(actual), len(expected), f"{label} length")
        for index, (left, right) in enumerate(zip(actual, expected)):
            self.close(left, right, f"{label}[{index}]")


def parse_hashes(path: Path) -> dict[str, str]:
    output: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, name = line.split(maxsplit=1)
        output[name.strip()] = digest
    return output


def sheet_targets(archive: zipfile.ZipFile) -> dict[str, str]:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    by_id = {
        node.attrib["Id"]: node.attrib["Target"]
        for node in relations.findall(f"{{{PACKAGE_REL}}}Relationship")
    }
    output: dict[str, str] = {}
    for sheet in workbook.findall(f".//{{{MAIN}}}sheet"):
        target = by_id[sheet.attrib[f"{{{OFFICE_REL}}}id"]]
        if target.startswith("/"):
            target = target[1:]
        else:
            target = posixpath.normpath(posixpath.join("xl", target))
        output[sheet.attrib["name"]] = target
    return output


def shared_strings(archive: zipfile.ZipFile) -> list[str]:
    try:
        root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    except KeyError:
        return []
    return [
        "".join(node.text or "" for node in item.iter(f"{{{MAIN}}}t"))
        for item in root.findall(f"{{{MAIN}}}si")
    ]


def read_cells(
    archive: zipfile.ZipFile, target: str, strings: list[str]
) -> dict[str, dict[str, object]]:
    root = ET.fromstring(archive.read(target))
    output: dict[str, dict[str, object]] = {}
    for cell in root.findall(f".//{{{MAIN}}}c"):
        reference = cell.attrib["r"]
        kind = cell.attrib.get("t")
        formula = cell.find(f"{{{MAIN}}}f")
        value_node = cell.find(f"{{{MAIN}}}v")
        raw = None if value_node is None else (value_node.text or "")
        if kind == "inlineStr":
            value: object = "".join(
                node.text or "" for node in cell.iter(f"{{{MAIN}}}t")
            )
        elif raw is None:
            value = None
        elif kind == "s":
            value = strings[int(raw)]
        elif kind == "b":
            value = raw == "1"
        elif kind in {"str", "e"}:
            value = raw
        else:
            value = float(raw)
        output[reference] = {
            "value": value,
            "raw": raw,
            "formula": None if formula is None else (formula.text or ""),
        }
    return output


def extract(cells: dict[str, dict[str, object]]) -> tuple[list[dict], list[dict]]:
    near: list[dict] = []
    far: list[dict] = []
    for row in range(3, 23):
        time = float(cells[f"B{row}"]["value"])
        for column, output in (("C", near), ("D", far)):
            entry = cells.get(f"{column}{row}")
            if entry is None or entry["value"] is None:
                continue
            period = float(entry["value"])
            output.append(
                {
                    "time": time,
                    "period": period,
                    "omega2": (2.0 * math.pi / period) ** 2,
                    "time_cell": f"B{row}",
                    "period_cell": f"{column}{row}",
                }
            )
    return near, far


def difference_matrix() -> list[list[int]]:
    output = [[0] * 10 for _ in range(9)]
    for index in range(9):
        output[index][index] = -1
        output[index][index + 1] = 1
    return output


def return_matrix(difference: list[list[int]]) -> list[list[int]]:
    output = [[0] * 40 for _ in range(36)]
    for block in range(4):
        for row in range(9):
            for column in range(10):
                output[block * 9 + row][block * 10 + column] = difference[row][column]
    return output


def differential_matrix(difference: list[list[int]]) -> list[list[int]]:
    output = [[0] * 40 for _ in range(18)]
    for orientation in range(2):
        for row in range(9):
            target = orientation * 9 + row
            present = orientation * 10
            background = 20 + orientation * 10
            for column in range(10):
                output[target][present + column] = difference[row][column]
                output[target][background + column] = -difference[row][column]
    return output


def matvec(matrix: list[list[int]], vector: list[float]) -> list[float]:
    return [sum(weight * value for weight, value in zip(row, vector)) for row in matrix]


def gram(matrix: list[list[int]]) -> list[list[int]]:
    return [
        [sum(a * b for a, b in zip(left, right)) for right in matrix]
        for left in matrix
    ]


def exact_rank(matrix: list[list[int]]) -> int:
    work = [[Fraction(value) for value in row] for row in matrix]
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [a - scale * b for a, b in zip(work[row], work[pivot_row])]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def summarize(values: list[float]) -> dict[str, float | int]:
    return {
        "count": len(values),
        "mean_s-2": sum(values) / len(values),
        "minimum_s-2": min(values),
        "maximum_s-2": max(values),
        "rms_s-2": math.sqrt(sum(value * value for value in values) / len(values)),
        "sum_s-2": sum(values),
    }


def check_summary(audit: Audit, actual: dict, values: list[float], label: str) -> None:
    expected = summarize(values)
    audit.equal(actual["count"], expected["count"], f"{label} count")
    for key in ("mean_s-2", "minimum_s-2", "maximum_s-2", "rms_s-2", "sum_s-2"):
        audit.close(actual[key], expected[key], f"{label} {key}")


def verify_manifest(audit: Audit) -> None:
    manifest = parse_hashes(HERE / "MANIFEST.sha256")
    audit.equal(set(manifest), MANIFEST_MEMBERS, "manifest member set")
    for name in sorted(MANIFEST_MEMBERS):
        audit.equal(sha256(HERE / name), manifest[name], f"manifest hash {name}")
    seal = parse_hashes(HERE / "LANE_SEAL.sha256")
    audit.equal(set(seal), {"MANIFEST.sha256"}, "outer seal member")
    audit.equal(sha256(HERE / "MANIFEST.sha256"), seal["MANIFEST.sha256"], "outer seal hash")

    result_bytes = bytearray((HERE / "RESULT.json").read_bytes())
    result_bytes[len(result_bytes) // 2] ^= 1
    audit.require(
        hashlib.sha256(result_bytes).hexdigest() != manifest["RESULT.json"],
        "mutated result rejected",
    )
    manifest_bytes = bytearray((HERE / "MANIFEST.sha256").read_bytes())
    manifest_bytes[0] ^= 1
    audit.require(
        hashlib.sha256(manifest_bytes).hexdigest() != seal["MANIFEST.sha256"],
        "mutated manifest rejected",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preseal", action="store_true")
    args = parser.parse_args()
    audit = Audit()

    audit.equal(sha256(WORKBOOK), WORKBOOK_SHA256, "workbook SHA-256")
    with zipfile.ZipFile(WORKBOOK) as archive:
        targets = sheet_targets(archive)
        strings = shared_strings(archive)
        audit.require({"a", "b"}.issubset(targets), "required sheets")
        present_cells = read_cells(archive, targets["a"], strings)
        background_cells = read_cells(archive, targets["b"], strings)

    for label, cells in (("present", present_cells), ("background", background_cells)):
        audit.equal(cells["B2"]["value"], "Time (day)", f"{label} time header")
        audit.equal(cells["C2"]["value"], "Period at near position (s)", f"{label} near header")
        audit.equal(cells["D2"]["value"], "Period at far position (s)", f"{label} far header")
        for row in range(3, 23):
            for column in ("B", "C", "D"):
                entry = cells.get(f"{column}{row}")
                if entry is not None:
                    audit.equal(entry["formula"], None, f"{label} {column}{row} no formula")
        audit.require(
            all(f"C{row}" in cells and f"D{row}" not in cells for row in range(3, 13)),
            f"{label} first block near-only",
        )
        audit.require(
            all(f"D{row}" in cells and f"C{row}" not in cells for row in range(13, 23)),
            f"{label} second block far-only",
        )

    pn, pf = extract(present_cells)
    bn, bf = extract(background_cells)
    groups = (pn, pf, bn, bf)
    for label, series in zip(("PN", "PF", "BN", "BF"), groups):
        audit.equal(len(series), 10, f"{label} endpoint count")
        audit.require(all(item["period"] > 0 for item in series), f"{label} positive periods")
        audit.require(
            all(series[i]["time"] < series[i + 1]["time"] for i in range(9)),
            f"{label} chronological order",
        )
        audit.require(
            all(abs(item["period"] * 1e5 - round(item["period"] * 1e5)) < 2e-6 for item in series),
            f"{label} nominal 1e-5 second period grid",
        )
    for label, near, far in (("present", pn, pf), ("background", bn, bf)):
        for index in range(9):
            audit.require(
                near[index]["time"] < far[index]["time"] < near[index + 1]["time"],
                f"{label} N-F-N chronology {index}",
            )
            audit.require(
                far[index]["time"] < near[index + 1]["time"] < far[index + 1]["time"],
                f"{label} F-N-F chronology {index}",
            )

    payload = json.loads((HERE / "RESULT.json").read_text(encoding="utf-8"))
    audit.equal(payload["schema"], "WAC_HUST_TOS_ROUNDTRIP_HISTORY_RESIDUAL_V001", "schema")
    audit.equal(payload["source"]["workbook_sha256"], WORKBOOK_SHA256, "payload source hash")
    audit.equal(payload["source"]["sheet_present"], "a", "present sheet")
    audit.equal(payload["source"]["sheet_background"], "b", "background sheet")
    audit.require("sequence ordinal" in payload["observable"]["cross_panel_pairing"], "ordinal pairing typed")
    audit.require("not a matched experimental trial" in payload["orientation_decomposition"]["pairing_ceiling"], "orientation pairing ceiling")

    payload_names = (
        "source_present_near",
        "source_present_far",
        "source_absent_background_near",
        "source_absent_background_far",
    )
    for name, independent in zip(payload_names, groups):
        actual = payload["inputs"][name]
        audit.equal(len(actual), 10, f"{name} length")
        for index, (left, right) in enumerate(zip(actual, independent)):
            audit.equal(left["cell_time"], right["time_cell"], f"{name}[{index}] time cell")
            audit.equal(left["cell_period"], right["period_cell"], f"{name}[{index}] period cell")
            audit.close(left["time_day"], right["time"], f"{name}[{index}] time")
            audit.close(left["period_s"], right["period"], f"{name}[{index}] period")
            audit.close(left["omega2_s-2"], right["omega2"], f"{name}[{index}] omega2")

    difference = difference_matrix()
    returns_weight = return_matrix(difference)
    differentials_weight = differential_matrix(difference)
    custody = payload["linear_custody"]
    audit.equal(custody["return_weight_matrix_36x40"], returns_weight, "return weights")
    audit.equal(custody["differential_weight_matrix_18x40"], differentials_weight, "differential weights")
    audit.equal(exact_rank(returns_weight), 36, "return exact rank")
    audit.equal(exact_rank(differentials_weight), 18, "differential exact rank")
    audit.equal(custody["return_weight_rank_exact"], 36, "reported return rank")
    audit.equal(custody["differential_weight_rank_exact"], 18, "reported differential rank")

    return_gram = gram(returns_weight)
    differential_gram = gram(differentials_weight)
    audit.equal(custody["return_unit_endpoint_overlap_gram_36x36"], return_gram, "return Gram")
    audit.equal(custody["differential_unit_endpoint_overlap_gram_18x18"], differential_gram, "differential Gram")
    for row in range(36):
        for column in range(36):
            same = row // 9 == column // 9
            expected = 2 if row == column else (-1 if same and abs(row - column) == 1 else 0)
            audit.equal(return_gram[row][column], expected, f"return Gram {row},{column}")
    for row in range(18):
        for column in range(18):
            same = row // 9 == column // 9
            expected = 4 if row == column else (-2 if same and abs(row - column) == 1 else 0)
            audit.equal(differential_gram[row][column], expected, f"differential Gram {row},{column}")

    endpoints = [item["omega2"] for series in groups for item in series]
    returns = matvec(returns_weight, endpoints)
    differentials = matvec(differentials_weight, endpoints)
    return_groups = {
        "source_present_N_F_N": returns[0:9],
        "source_present_F_N_F": returns[9:18],
        "source_absent_background_N_F_N": returns[18:27],
        "source_absent_background_F_N_F": returns[27:36],
    }
    differential_groups = {
        "present_minus_background_N_F_N": differentials[0:9],
        "present_minus_background_F_N_F": differentials[9:18],
    }
    for name, expected in return_groups.items():
        actual = payload["returns"][name]
        audit.close_list(actual["values_s-2"], expected, f"{name} values")
        check_summary(audit, actual["summary"], expected, f"{name} summary")
    for name, expected in differential_groups.items():
        actual = payload["panel_differentials"][name]
        audit.close_list(actual["values_s-2"], expected, f"{name} values")
        check_summary(audit, actual["summary"], expected, f"{name} summary")

    even = [(left + right) / 2.0 for left, right in zip(differentials[:9], differentials[9:])]
    odd = [(left - right) / 2.0 for left, right in zip(differentials[:9], differentials[9:])]
    actual_even = payload["orientation_decomposition"]["even_orientation_common_component"]
    actual_odd = payload["orientation_decomposition"]["odd_N_F_N_minus_F_N_F_half"]
    audit.close_list(actual_even["values_s-2"], even, "even values")
    audit.close_list(actual_odd["values_s-2"], odd, "odd values")
    check_summary(audit, actual_even["summary"], even, "even summary")
    check_summary(audit, actual_odd["summary"], odd, "odd summary")

    frozen_means = {
        "source_present_N_F_N": 3.647457758732247e-10,
        "source_present_F_N_F": 3.496735157677661e-10,
        "source_absent_background_N_F_N": 5.498090526271760e-10,
        "source_absent_background_F_N_F": 5.324953327068938e-10,
    }
    for name, expected in frozen_means.items():
        audit.close(sum(return_groups[name]) / 9.0, expected, f"frozen mean {name}")
    audit.close(sum(differentials[:9]) / 9.0, -1.8506327675395133e-10, "frozen dN mean")
    audit.close(sum(differentials[9:]) / 9.0, -1.8282181693912770e-10, "frozen dF mean")
    audit.close(sum(even) / 9.0, -1.8394254684653953e-10, "frozen even mean")
    audit.close(sum(odd) / 9.0, -1.1207299074118054e-12, "frozen odd mean")
    audit.close(math.sqrt(sum(value * value for value in odd) / 9.0), 6.524819132094064e-11, "frozen odd RMS")

    for label, values, series in (
        ("PN", returns[0:9], pn),
        ("PF", returns[9:18], pf),
        ("BN", returns[18:27], bn),
        ("BF", returns[27:36], bf),
    ):
        audit.close(sum(values), series[-1]["omega2"] - series[0]["omega2"], f"{label} telescope")
    audit.close(
        sum(differentials[:9]),
        (pn[-1]["omega2"] - pn[0]["omega2"]) - (bn[-1]["omega2"] - bn[0]["omega2"]),
        "dN telescope",
    )
    audit.close(
        sum(differentials[9:]),
        (pf[-1]["omega2"] - pf[0]["omega2"]) - (bf[-1]["omega2"] - bf[0]["omega2"]),
        "dF telescope",
    )

    durations: dict[str, list[float]] = {}
    for name, series in (
        ("source_present_N_F_N_days", pn),
        ("source_present_F_N_F_days", pf),
        ("source_absent_background_N_F_N_days", bn),
        ("source_absent_background_F_N_F_days", bf),
    ):
        durations[name] = [series[i + 1]["time"] - series[i]["time"] for i in range(9)]
    durations["present_minus_background_N_F_N_days"] = [
        a - b
        for a, b in zip(durations["source_present_N_F_N_days"], durations["source_absent_background_N_F_N_days"])
    ]
    durations["present_minus_background_F_N_F_days"] = [
        a - b
        for a, b in zip(durations["source_present_F_N_F_days"], durations["source_absent_background_F_N_F_days"])
    ]
    for name, expected in durations.items():
        audit.close_list(payload["timing"]["loop_duration_values"][name], expected, f"duration {name}")
    audit.close(max(abs(value) for value in durations["present_minus_background_N_F_N_days"]), 0.4445027, "N-F-N duration difference")
    audit.close(max(abs(value) for value in durations["present_minus_background_F_N_F_days"]), 0.0421588, "F-N-F duration difference")
    audit.require("ordinal correspondence is not an experimental match" in payload["timing"]["ceiling"], "timing matching ceiling")

    custody_doc = json.loads((HERE / "SOURCE_CUSTODY.json").read_text(encoding="utf-8"))
    audit.equal(custody_doc["source"]["sha256"], WORKBOOK_SHA256, "custody hash")
    audit.equal(custody_doc["source"]["copied_into_this_lane"], False, "source not copied")
    audit.require("sequence ordinal only" in custody_doc["panel_semantics"]["cross_panel_index"], "custody ordinal semantics")
    audit.require("do not encode a matched cross-panel trial" in custody_doc["panel_semantics"]["worksheet_headers_alone"], "headers do not own matching")

    theorem = (HERE / "THEOREM.md").read_text(encoding="utf-8")
    result_md = (HERE / "RESULT.md").read_text(encoding="utf-8")
    self_audit = (HERE / "SELF_AUDIT.md").read_text(encoding="utf-8")
    readme = (HERE / "README.md").read_text(encoding="utf-8")
    joined = " ".join("\n".join((theorem, result_md, self_audit, readme)).lower().split())
    for phrase in (
        "no matched no-excursion",
        "history-confound diagnostic",
        "sequence ordinal",
        "not measurement precision",
        "not an experimental match",
        "not empirical covariance",
        "telescopes",
        "no row covariance",
        "does **not** estimate `beta_tm`",
        "gravity emergence",
        "full `gc16`",
        "new value of \\(g\\)",
    ):
        audit.require(phrase in joined, f"required ceiling phrase: {phrase}")
    audit.require("no canonical model" in joined, "no canonical model edit")
    audit.require("cannot identify hysteresis or memory" in payload["interpretation"]["counterfactual_ceiling"], "JSON no-excursion ceiling")
    audit.require("computational digits" in payload["interpretation"]["printed_precision_ceiling"], "JSON precision ceiling")

    forbidden_keys = {
        "standard_error",
        "standard_errors",
        "p_value",
        "p_values",
        "confidence_interval",
        "confidence_intervals",
        "beta_tm_estimate",
        "lineage_charge_estimate",
        "coverage_probability",
    }
    keys: set[str] = set()

    def collect(value: object) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                keys.add(str(key).lower())
                collect(child)
        elif isinstance(value, list):
            for child in value:
                collect(child)

    collect(payload)
    audit.require(not (keys & forbidden_keys), "no unauthorized inferential fields")

    core_count = audit.count
    transcript = (
        "HUST_TOS_ROUNDTRIP_HISTORY_HOSTILE_AUDIT: PASS\n"
        f"core_checks_passed: {core_count}\n"
        "independent_official_xlsx_reparse: PASS\n"
        "chronology_and_orientation_labels: PASS\n"
        "exact_weight_ranks_grams_and_telescoping: PASS\n"
        "all_reported_means_and_duration_differences: PASS\n"
        "ordinal_cross_panel_pairing_ceiling: PASS\n"
        "printed_precision_ceiling: PASS\n"
        "matched_no_excursion_counterfactual: ABSENT\n"
        "manifest_and_lane_seal: PASS\n"
        "scientific_ceiling: HISTORY_CONFOUND_DIAGNOSTIC_ONLY__NO_MEMORY_LINEAGE_BETA_TM_COVERAGE_GRAVITY_EMERGENCE_GC16_OR_NEW_G_CLAIM\n"
    )
    if not args.preseal:
        verify_manifest(audit)
    (HERE / "HOSTILE_VERIFICATION.txt").write_text(transcript, encoding="utf-8")
    print(transcript, end="")


if __name__ == "__main__":
    main()
