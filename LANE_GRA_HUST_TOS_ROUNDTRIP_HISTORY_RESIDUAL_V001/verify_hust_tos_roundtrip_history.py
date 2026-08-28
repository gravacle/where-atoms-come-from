#!/usr/bin/env python3
"""Independent verifier for GRA-HUST18-TOS-RTHR-V001.

This verifier deliberately does not import the production analyzer.  It
re-opens the pinned OOXML workbook, reconstructs every scored value and exact
integer matrix, and checks the sealed payload and scientific ceiling.
"""

from __future__ import annotations

import hashlib
import json
import math
import posixpath
import xml.etree.ElementTree as ET
import zipfile
from fractions import Fraction
from pathlib import Path


LANE = Path(__file__).resolve().parent
ROOT = LANE.parent
WORKBOOK_REL = (
    "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/SOURCE/"
    "41586_2018_431_MOESM3_ESM.xlsx"
)
WORKBOOK = (LANE / WORKBOOK_REL).resolve()

MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
OD_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"

DEPENDENCIES = {
    WORKBOOK_REL: "331ea6ee8a6a2558d4f6a8ffa4cd52e4c5d52d47f2e2a697889d0d3c8a3ad27a",
    "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/SOURCE_CUSTODY.json":
        "76040b79f3ead7970298df5b784e3cab0ea637351b6dc881be1afe6796a6e474",
    "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/THEOREM.md":
        "44eb8c81a3d84dfa6829bcd6971d0261215877af0529318eab5cedbd3980c340",
    "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "e8625b4cbf67d73927db495e8111e3ffbd4e85f46b80183e55fbf8b4391d0b2e",
}

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


def parse_hash_file(path: Path) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, name = line.split(maxsplit=1)
        parsed[name.strip()] = digest
    return parsed


def workbook_sheet_targets(archive: zipfile.ZipFile) -> dict[str, str]:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {
        relation.attrib["Id"]: relation.attrib["Target"]
        for relation in relations.findall(f"{{{PKG_REL}}}Relationship")
    }
    output: dict[str, str] = {}
    for sheet in workbook.findall(f".//{{{MAIN}}}sheet"):
        target = targets[sheet.attrib[f"{{{OD_REL}}}id"]]
        if target.startswith("/"):
            target = target[1:]
        else:
            target = posixpath.normpath(posixpath.join("xl", target))
        output[sheet.attrib["name"]] = target
    return output


def workbook_shared_strings(archive: zipfile.ZipFile) -> list[str]:
    try:
        root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    except KeyError:
        return []
    return [
        "".join(node.text or "" for node in item.iter(f"{{{MAIN}}}t"))
        for item in root.findall(f"{{{MAIN}}}si")
    ]


def read_sheet(archive: zipfile.ZipFile, target: str, shared: list[str]) -> dict[str, object]:
    root = ET.fromstring(archive.read(target))
    output: dict[str, object] = {}
    for cell in root.findall(f".//{{{MAIN}}}c"):
        reference = cell.attrib["r"]
        kind = cell.attrib.get("t")
        value_node = cell.find(f"{{{MAIN}}}v")
        if kind == "inlineStr":
            value: object = "".join(node.text or "" for node in cell.iter(f"{{{MAIN}}}t"))
        elif value_node is None:
            value = None
        elif kind == "s":
            value = shared[int(value_node.text or "0")]
        elif kind == "b":
            value = value_node.text == "1"
        elif kind in {"str", "e"}:
            value = value_node.text or ""
        else:
            value = float(value_node.text or "nan")
        output[reference] = value
    return output


def extract_panel(cells: dict[str, object]) -> tuple[list[dict], list[dict]]:
    near: list[dict] = []
    far: list[dict] = []
    for row in range(3, 23):
        time = cells.get(f"B{row}")
        if f"C{row}" in cells:
            period = float(cells[f"C{row}"])
            near.append({
                "cell_time": f"B{row}",
                "cell_period": f"C{row}",
                "time_day": float(time),
                "period_s": period,
                "omega2_s-2": (2 * math.pi / period) ** 2,
            })
        if f"D{row}" in cells:
            period = float(cells[f"D{row}"])
            far.append({
                "cell_time": f"B{row}",
                "cell_period": f"D{row}",
                "time_day": float(time),
                "period_s": period,
                "omega2_s-2": (2 * math.pi / period) ** 2,
            })
    return near, far


def difference_matrix() -> list[list[int]]:
    output: list[list[int]] = []
    for row in range(9):
        values = [0] * 10
        values[row] = -1
        values[row + 1] = 1
        output.append(values)
    return output


def block_diagonal(block: list[list[int]], copies: int) -> list[list[int]]:
    rows = len(block)
    columns = len(block[0])
    output = [[0] * (columns * copies) for _ in range(rows * copies)]
    for copy in range(copies):
        for row in range(rows):
            for column in range(columns):
                output[copy * rows + row][copy * columns + column] = block[row][column]
    return output


def differential_weights(block: list[list[int]]) -> list[list[int]]:
    output = [[0] * 40 for _ in range(18)]
    for orientation in range(2):
        for row in range(9):
            output_row = orientation * 9 + row
            present_offset = orientation * 10
            background_offset = 20 + orientation * 10
            for column in range(10):
                output[output_row][present_offset + column] = block[row][column]
                output[output_row][background_offset + column] = -block[row][column]
    return output


def matvec(matrix: list[list[int]], vector: list[float]) -> list[float]:
    return [sum(weight * value for weight, value in zip(row, vector)) for row in matrix]


def gram(matrix: list[list[int]]) -> list[list[int]]:
    return [
        [sum(a * b for a, b in zip(left, right)) for right in matrix]
        for left in matrix
    ]


def rational_rank(matrix: list[list[int]]) -> int:
    work = [[Fraction(value) for value in row] for row in matrix]
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next((r for r in range(pivot_row, len(work)) if work[r][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(pivot_row + 1, len(work)):
            if not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [a - factor * b for a, b in zip(work[row], work[pivot_row])]
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


def check_summary(audit: Audit, actual: dict, expected_values: list[float], label: str) -> None:
    expected = summarize(expected_values)
    audit.equal(actual["count"], expected["count"], f"{label} count")
    for key in ("mean_s-2", "minimum_s-2", "maximum_s-2", "rms_s-2", "sum_s-2"):
        audit.close(actual[key], expected[key], f"{label} {key}")


def verify_manifest(audit: Audit) -> None:
    manifest = parse_hash_file(LANE / "MANIFEST.sha256")
    audit.equal(set(manifest), MANIFEST_MEMBERS, "manifest member set")
    for name in sorted(MANIFEST_MEMBERS):
        audit.equal(sha256(LANE / name), manifest[name], f"manifest hash {name}")
    seal = parse_hash_file(LANE / "LANE_SEAL.sha256")
    audit.equal(set(seal), {"MANIFEST.sha256"}, "outer seal target")
    audit.equal(sha256(LANE / "MANIFEST.sha256"), seal["MANIFEST.sha256"], "outer seal digest")

    # Negative sentinels prove the comparison would reject a one-byte mutation.
    result_bytes = (LANE / "RESULT.json").read_bytes()
    tampered_result = bytearray(result_bytes)
    tampered_result[len(tampered_result) // 2] ^= 1
    audit.require(
        hashlib.sha256(tampered_result).hexdigest() != manifest["RESULT.json"],
        "tampered RESULT sentinel rejected",
    )
    manifest_bytes = (LANE / "MANIFEST.sha256").read_bytes()
    tampered_manifest = bytearray(manifest_bytes)
    tampered_manifest[0] ^= 1
    audit.require(
        hashlib.sha256(tampered_manifest).hexdigest() != seal["MANIFEST.sha256"],
        "tampered manifest sentinel rejected",
    )


def main() -> None:
    audit = Audit()

    declared_dependencies = parse_hash_file(LANE / "DEPENDENCIES.sha256")
    audit.equal(declared_dependencies, DEPENDENCIES, "dependency declaration")
    for relative, expected in DEPENDENCIES.items():
        audit.equal(sha256((LANE / relative).resolve()), expected, f"dependency hash {relative}")
    workbook_bytes = WORKBOOK.read_bytes()
    mutated_workbook = bytearray(workbook_bytes)
    mutated_workbook[len(mutated_workbook) // 3] ^= 1
    audit.require(
        hashlib.sha256(mutated_workbook).hexdigest() != DEPENDENCIES[WORKBOOK_REL],
        "mutated workbook dependency rejected",
    )

    with zipfile.ZipFile(WORKBOOK) as archive:
        sheets = workbook_sheet_targets(archive)
        strings = workbook_shared_strings(archive)
        audit.require({"a", "b"}.issubset(sheets), "required workbook sheets")
        cells_present = read_sheet(archive, sheets["a"], strings)
        cells_background = read_sheet(archive, sheets["b"], strings)

    audit.equal(cells_present.get("B2"), "Time (day)", "sheet a time header")
    audit.equal(cells_present.get("C2"), "Period at near position (s)", "sheet a near header")
    audit.equal(cells_present.get("D2"), "Period at far position (s)", "sheet a far header")
    audit.equal(cells_background.get("B2"), "Time (day)", "sheet b time header")
    audit.equal(cells_background.get("C2"), "Period at near position (s)", "sheet b near header")
    audit.equal(cells_background.get("D2"), "Period at far position (s)", "sheet b far header")

    pn, pf = extract_panel(cells_present)
    bn, bf = extract_panel(cells_background)
    for label, series in (("PN", pn), ("PF", pf), ("BN", bn), ("BF", bf)):
        audit.equal(len(series), 10, f"{label} endpoint count")
        audit.require(all(item["period_s"] > 0 for item in series), f"{label} periods positive")
        audit.require(
            all(series[index]["time_day"] < series[index + 1]["time_day"] for index in range(9)),
            f"{label} chronological order",
        )
    for label, near, far in (("present", pn, pf), ("background", bn, bf)):
        for index in range(9):
            audit.require(
                near[index]["time_day"] < far[index]["time_day"] < near[index + 1]["time_day"],
                f"{label} N-F-N triple {index}",
            )
            audit.require(
                far[index]["time_day"] < near[index + 1]["time_day"] < far[index + 1]["time_day"],
                f"{label} F-N-F triple {index}",
            )

    payload = json.loads((LANE / "RESULT.json").read_text(encoding="utf-8"))
    audit.equal(payload["schema"], "WAC_HUST_TOS_ROUNDTRIP_HISTORY_RESIDUAL_V001", "result schema")
    audit.equal(payload["date"], "2026-08-27", "result date")
    audit.equal(payload["source"]["workbook_sha256"], DEPENDENCIES[WORKBOOK_REL], "result source hash")
    audit.equal(payload["source"]["sheet_present"], "a", "present sheet declaration")
    audit.equal(payload["source"]["sheet_background"], "b", "background sheet declaration")
    audit.equal(payload["source"]["range_each"], "B3:D22", "scored range declaration")
    audit.require(
        "sequence ordinal" in payload["observable"]["cross_panel_pairing"],
        "cross-panel ordinal pairing typed",
    )

    input_names = (
        "source_present_near",
        "source_present_far",
        "source_absent_background_near",
        "source_absent_background_far",
    )
    independent_series = (pn, pf, bn, bf)
    for name, expected_series in zip(input_names, independent_series):
        actual_series = payload["inputs"][name]
        audit.equal(len(actual_series), 10, f"result {name} length")
        for index, (actual, expected) in enumerate(zip(actual_series, expected_series)):
            audit.equal(actual["cell_time"], expected["cell_time"], f"{name}[{index}] time cell")
            audit.equal(actual["cell_period"], expected["cell_period"], f"{name}[{index}] period cell")
            audit.close(actual["time_day"], expected["time_day"], f"{name}[{index}] time")
            audit.close(actual["period_s"], expected["period_s"], f"{name}[{index}] period")
            audit.close(actual["omega2_s-2"], expected["omega2_s-2"], f"{name}[{index}] omega2")

    expected_labels = [
        *(f"P_N_{i}" for i in range(10)),
        *(f"P_F_{i}" for i in range(10)),
        *(f"B_N_{i}" for i in range(10)),
        *(f"B_F_{i}" for i in range(10)),
    ]
    audit.equal(payload["inputs"]["endpoint_order"], expected_labels, "endpoint order")
    audit.equal(payload["counts"]["endpoint_omega2_values"], 40, "endpoint count")
    audit.equal(payload["counts"]["all_panel_returns"], 36, "return count")
    audit.equal(payload["counts"]["present_minus_background_differentials"], 18, "differential count")

    d = difference_matrix()
    returns_matrix = block_diagonal(d, 4)
    diff_matrix = differential_weights(d)
    custody = payload["linear_custody"]
    audit.equal(custody["return_weight_matrix_36x40"], returns_matrix, "return matrix")
    audit.equal(custody["differential_weight_matrix_18x40"], diff_matrix, "differential matrix")
    audit.equal(len(returns_matrix), 36, "return matrix rows")
    audit.require(all(len(row) == 40 for row in returns_matrix), "return matrix columns")
    audit.equal(len(diff_matrix), 18, "differential matrix rows")
    audit.require(all(len(row) == 40 for row in diff_matrix), "differential matrix columns")
    audit.require(all(sum(abs(v) for v in row) == 2 for row in returns_matrix), "two endpoint return weights")
    audit.require(all(sum(row) == 0 for row in returns_matrix), "return weights annihilate constants")
    audit.require(all(sum(abs(v) for v in row) == 4 for row in diff_matrix), "four endpoint differential weights")
    audit.require(all(sum(row) == 0 for row in diff_matrix), "differential weights annihilate constants")
    audit.equal(rational_rank(returns_matrix), 36, "independent exact return rank")
    audit.equal(rational_rank(diff_matrix), 18, "independent exact differential rank")
    audit.equal(custody["return_weight_rank_exact"], 36, "reported return rank")
    audit.equal(custody["differential_weight_rank_exact"], 18, "reported differential rank")

    return_gram = gram(returns_matrix)
    diff_gram = gram(diff_matrix)
    audit.equal(custody["return_unit_endpoint_overlap_gram_36x36"], return_gram, "return Gram")
    audit.equal(custody["differential_unit_endpoint_overlap_gram_18x18"], diff_gram, "differential Gram")
    for row in range(36):
        for column in range(36):
            same_block = row // 9 == column // 9
            expected = 2 if row == column else (-1 if same_block and abs(row - column) == 1 else 0)
            audit.equal(return_gram[row][column], expected, f"return Gram[{row},{column}]")
    for row in range(18):
        for column in range(18):
            same_block = row // 9 == column // 9
            expected = 4 if row == column else (-2 if same_block and abs(row - column) == 1 else 0)
            audit.equal(diff_gram[row][column], expected, f"differential Gram[{row},{column}]")

    endpoint_values = [item["omega2_s-2"] for series in independent_series for item in series]
    return_values = matvec(returns_matrix, endpoint_values)
    differential_values = matvec(diff_matrix, endpoint_values)
    return_groups = {
        "source_present_N_F_N": return_values[0:9],
        "source_present_F_N_F": return_values[9:18],
        "source_absent_background_N_F_N": return_values[18:27],
        "source_absent_background_F_N_F": return_values[27:36],
    }
    differential_groups = {
        "present_minus_background_N_F_N": differential_values[0:9],
        "present_minus_background_F_N_F": differential_values[9:18],
    }
    for name, expected in return_groups.items():
        actual = payload["returns"][name]
        audit.close_list(actual["values_s-2"], expected, f"returns {name}")
        check_summary(audit, actual["summary"], expected, f"returns {name}")
    for name, expected in differential_groups.items():
        actual = payload["panel_differentials"][name]
        audit.close_list(actual["values_s-2"], expected, f"differentials {name}")
        check_summary(audit, actual["summary"], expected, f"differentials {name}")

    even = [(a + b) / 2 for a, b in zip(differential_values[:9], differential_values[9:])]
    odd = [(a - b) / 2 for a, b in zip(differential_values[:9], differential_values[9:])]
    even_actual = payload["orientation_decomposition"]["even_orientation_common_component"]
    odd_actual = payload["orientation_decomposition"]["odd_N_F_N_minus_F_N_F_half"]
    audit.close_list(even_actual["values_s-2"], even, "orientation even")
    audit.close_list(odd_actual["values_s-2"], odd, "orientation odd")
    check_summary(audit, even_actual["summary"], even, "orientation even")
    check_summary(audit, odd_actual["summary"], odd, "orientation odd")
    audit.close(
        payload["orientation_decomposition"]["mean_full_orientation_difference_s-2"],
        2 * sum(odd) / 9,
        "full orientation mean difference",
    )

    expected_means = {
        "source_present_N_F_N": 3.647457758732247e-10,
        "source_present_F_N_F": 3.496735157677661e-10,
        "source_absent_background_N_F_N": 5.498090526271760e-10,
        "source_absent_background_F_N_F": 5.324953327068938e-10,
    }
    for name, expected in expected_means.items():
        audit.close(sum(return_groups[name]) / 9, expected, f"frozen mean {name}")
    audit.close(sum(differential_values[:9]) / 9, -1.8506327675395133e-10, "frozen N-F-N differential mean")
    audit.close(sum(differential_values[9:]) / 9, -1.8282181693912770e-10, "frozen F-N-F differential mean")
    audit.close(sum(even) / 9, -1.8394254684653953e-10, "frozen even mean")
    audit.close(sum(odd) / 9, -1.1207299074118054e-12, "frozen odd mean")
    audit.close(math.sqrt(sum(value * value for value in odd) / 9), 6.524819132094064e-11, "frozen odd RMS")
    audit.require(
        "not a matched experimental trial"
        in payload["orientation_decomposition"]["pairing_ceiling"],
        "orientation pairing ceiling",
    )

    for name, values, endpoints in (
        ("PN", return_values[0:9], pn),
        ("PF", return_values[9:18], pf),
        ("BN", return_values[18:27], bn),
        ("BF", return_values[27:36], bf),
    ):
        audit.close(sum(values), endpoints[-1]["omega2_s-2"] - endpoints[0]["omega2_s-2"], f"{name} telescope")

    timing_groups: dict[str, list[float]] = {}
    for name, series in (
        ("source_present_N_F_N_days", pn),
        ("source_present_F_N_F_days", pf),
        ("source_absent_background_N_F_N_days", bn),
        ("source_absent_background_F_N_F_days", bf),
    ):
        timing_groups[name] = [series[i + 1]["time_day"] - series[i]["time_day"] for i in range(9)]
    timing_groups["present_minus_background_N_F_N_days"] = [
        a - b for a, b in zip(timing_groups["source_present_N_F_N_days"], timing_groups["source_absent_background_N_F_N_days"])
    ]
    timing_groups["present_minus_background_F_N_F_days"] = [
        a - b for a, b in zip(timing_groups["source_present_F_N_F_days"], timing_groups["source_absent_background_F_N_F_days"])
    ]
    for name, expected in timing_groups.items():
        actual = payload["timing"]["loop_duration_values"][name]
        audit.close_list(actual, expected, f"timing {name}")
        reported = payload["timing"]["loop_duration_summaries"][name]
        audit.equal(reported["count"], 9, f"timing summary {name} count")
        audit.close(reported["mean_day"], sum(expected) / 9, f"timing summary {name} mean")
        audit.close(reported["minimum_day"], min(expected), f"timing summary {name} min")
        audit.close(reported["maximum_day"], max(expected), f"timing summary {name} max")
        audit.close(reported["maximum_absolute_day"], max(abs(value) for value in expected), f"timing summary {name} abs max")
    audit.close(
        max(abs(value) for value in timing_groups["present_minus_background_N_F_N_days"]),
        0.4445027,
        "frozen maximum N-F-N duration mismatch",
    )
    audit.close(
        max(abs(value) for value in timing_groups["present_minus_background_F_N_F_days"]),
        0.0421588,
        "frozen maximum F-N-F duration mismatch",
    )
    audit.require(
        "ordinal correspondence is not an experimental match"
        in payload["timing"]["ceiling"],
        "timing ordinal-pairing ceiling",
    )

    custody_source = json.loads((LANE / "SOURCE_CUSTODY.json").read_text(encoding="utf-8"))
    audit.equal(custody_source["source"]["sha256"], DEPENDENCIES[WORKBOOK_REL], "custody workbook hash")
    audit.equal(custody_source["source"]["bytes"], WORKBOOK.stat().st_size, "custody workbook bytes")
    audit.equal(custody_source["source"]["copied_into_this_lane"], False, "source not copied")
    audit.equal(custody_source["read_only_spreadsheet_inspection"]["workbook_modified"], False, "spreadsheet read-only")
    audit.equal([(item["sheet"], item["range"]) for item in custody_source["scored_ranges"]], [("a", "B3:D22"), ("b", "B3:D22")], "custody scored ranges")
    audit.require(
        "sequence ordinal only" in custody_source["panel_semantics"]["cross_panel_index"],
        "custody cross-panel index ceiling",
    )

    theorem = (LANE / "THEOREM.md").read_text(encoding="utf-8")
    result_md = (LANE / "RESULT.md").read_text(encoding="utf-8")
    self_audit = (LANE / "SELF_AUDIT.md").read_text(encoding="utf-8")
    joined = " ".join((theorem + "\n" + result_md + "\n" + self_audit).split())
    for required in (
        "no matched no-excursion",
        "history-confound diagnostic",
        "sequence ordinal",
        "not measurement precision",
        "not an experimental match",
        "empirical covariance matrices",
        "telescopes",
        "no standard error",
        "no row covariance",
        "does **not** estimate `beta_TM`",
        "gravity emergence",
        "No canonical MODEL",
    ):
        audit.require(required.lower() in joined.lower(), f"required ceiling phrase: {required}")
    forbidden_keys = {
        "standard_error", "standard_errors", "p_value", "p_values",
        "confidence_interval", "confidence_intervals", "beta_tm_estimate",
        "lineage_charge_estimate", "coverage_probability",
    }
    all_keys: set[str] = set()

    def collect_keys(value: object) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                all_keys.add(str(key).lower())
                collect_keys(child)
        elif isinstance(value, list):
            for child in value:
                collect_keys(child)

    collect_keys(payload)
    audit.require(not (all_keys & forbidden_keys), "no unauthorized inferential result fields")
    audit.require("cannot identify hysteresis or memory" in payload["interpretation"]["counterfactual_ceiling"], "JSON counterfactual ceiling")
    audit.require(
        "computational digits" in payload["interpretation"]["printed_precision_ceiling"],
        "JSON printed-precision ceiling",
    )

    verify_manifest(audit)

    print("HUST_TOS_ROUNDTRIP_HISTORY_VERIFIER: PASS")
    print(f"checks_passed: {audit.count}")
    print("independent_source_reparse: PASS")
    print("all_18_present_minus_background_residuals: PASS")
    print("exact_return_matrix_rank: 36")
    print("exact_differential_matrix_rank: 18")
    print("exact_reused_endpoint_overlap_grams: PASS")
    print("telescoping_and_duration_ledger: PASS")
    print("tamper_sentinels: PASS")
    print("scientific_ceiling: HISTORY_CONFOUND_DIAGNOSTIC_ONLY__NO_COVARIANCE_LINEAGE_BETA_TM_COVERAGE_OR_GRAVITY_EMERGENCE_CLAIM")


if __name__ == "__main__":
    main()
