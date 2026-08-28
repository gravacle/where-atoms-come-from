#!/usr/bin/env python3
"""Deterministic HUST ToS equal-configuration round-trip history diagnostic.

This analyzer reads only the pinned official HUST Figure-2 workbook.  The
observable is linear in the transformed endpoint vector y=(2*pi/T)^2.  It is
not a lineage estimand and it assigns no uncertainty or statistical coverage.
"""

from __future__ import annotations

import hashlib
import json
import math
import posixpath
import re
import xml.etree.ElementTree as ET
import zipfile
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LANE = Path(__file__).resolve().parent
PARENT = ROOT / "LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001"
WORKBOOK = PARENT / "SOURCE" / "41586_2018_431_MOESM3_ESM.xlsx"

NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS_PKG_REL = "http://schemas.openxmlformats.org/package/2006/relationships"

EXPECTED_DEPENDENCIES = {
    "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/SOURCE/41586_2018_431_MOESM3_ESM.xlsx":
        "331ea6ee8a6a2558d4f6a8ffa4cd52e4c5d52d47f2e2a697889d0d3c8a3ad27a",
    "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/SOURCE_CUSTODY.json":
        "76040b79f3ead7970298df5b784e3cab0ea637351b6dc881be1afe6796a6e474",
    "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/THEOREM.md":
        "44eb8c81a3d84dfa6829bcd6971d0261215877af0529318eab5cedbd3980c340",
    "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "e8625b4cbf67d73927db495e8111e3ffbd4e85f46b80183e55fbf8b4391d0b2e",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check_dependencies() -> dict[str, str]:
    observed: dict[str, str] = {}
    for relative, expected in EXPECTED_DEPENDENCIES.items():
        path = (LANE / relative).resolve()
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(f"dependency hash mismatch: {relative}: {actual}")
        observed[relative] = actual

    declared: dict[str, str] = {}
    for line in (LANE / "DEPENDENCIES.sha256").read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, relative = line.split(maxsplit=1)
        declared[relative.strip()] = digest
    if declared != EXPECTED_DEPENDENCIES:
        raise RuntimeError("DEPENDENCIES.sha256 differs from analyzer dependency contract")
    return observed


def column_index(reference: str) -> int:
    letters = re.match(r"[A-Z]+", reference)
    if letters is None:
        raise ValueError(reference)
    number = 0
    for character in letters.group(0):
        number = number * 26 + ord(character) - ord("A") + 1
    return number


class XlsxSource:
    """Minimal read-only OOXML reader with explicit cell custody."""

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
        values: list[str] = []
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
        result: dict[str, str] = {}
        for sheet in workbook.findall(f".//{{{NS_MAIN}}}sheet"):
            relationship_id = sheet.attrib[f"{{{NS_REL}}}id"]
            target = target_by_id[relationship_id]
            if target.startswith("/"):
                target = target.lstrip("/")
            else:
                target = posixpath.normpath(posixpath.join("xl", target))
            result[sheet.attrib["name"]] = target
        return result

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
                    value: object = "".join(
                        node.text or "" for node in cell.iter(f"{{{NS_MAIN}}}t")
                    )
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


def period_series(cells: dict[str, object]) -> tuple[list[dict], list[dict]]:
    near: list[dict] = []
    far: list[dict] = []
    for row in range(3, 23):
        time = float(cells[f"B{row}"])
        if cells.get(f"C{row}") is not None:
            period = float(cells[f"C{row}"])
            near.append({
                "cell_time": f"B{row}", "cell_period": f"C{row}",
                "time_day": time, "period_s": period,
                "omega2_s-2": (2.0 * math.pi / period) ** 2,
            })
        if cells.get(f"D{row}") is not None:
            period = float(cells[f"D{row}"])
            far.append({
                "cell_time": f"B{row}", "cell_period": f"D{row}",
                "time_day": time, "period_s": period,
                "omega2_s-2": (2.0 * math.pi / period) ** 2,
            })
    if len(near) != 10 or len(far) != 10:
        raise RuntimeError("expected ten near and ten far summaries")
    if not all(near[i]["time_day"] < far[i]["time_day"] < near[i + 1]["time_day"] for i in range(9)):
        raise RuntimeError("near-far-near alternation failed")
    if not all(far[i]["time_day"] < near[i + 1]["time_day"] < far[i + 1]["time_day"] for i in range(9)):
        raise RuntimeError("far-near-far alternation failed")
    return near, far


def zeros(rows: int, columns: int) -> list[list[int]]:
    return [[0 for _ in range(columns)] for _ in range(rows)]


def first_difference_matrix() -> list[list[int]]:
    matrix = zeros(9, 10)
    for index in range(9):
        matrix[index][index] = -1
        matrix[index][index + 1] = 1
    return matrix


def return_matrix() -> list[list[int]]:
    difference = first_difference_matrix()
    matrix = zeros(36, 40)
    for block in range(4):
        for row in range(9):
            for column in range(10):
                matrix[block * 9 + row][block * 10 + column] = difference[row][column]
    return matrix


def differential_matrix() -> list[list[int]]:
    difference = first_difference_matrix()
    matrix = zeros(18, 40)
    for row in range(9):
        for column in range(10):
            matrix[row][column] = difference[row][column]
            matrix[row][20 + column] = -difference[row][column]
            matrix[9 + row][10 + column] = difference[row][column]
            matrix[9 + row][30 + column] = -difference[row][column]
    return matrix


def matvec(matrix: list[list[int]], vector: list[float]) -> list[float]:
    return [sum(coefficient * value for coefficient, value in zip(row, vector)) for row in matrix]


def gram(matrix: list[list[int]]) -> list[list[int]]:
    return [
        [sum(a * b for a, b in zip(row_a, row_b)) for row_b in matrix]
        for row_a in matrix
    ]


def exact_rank(matrix: list[list[int]]) -> int:
    work = [[Fraction(value) for value in row] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows) if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        factor = work[pivot_row][column]
        work[pivot_row] = [value / factor for value in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [a - factor * b for a, b in zip(work[row], work[pivot_row])]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def summary(values: list[float]) -> dict:
    mean = sum(values) / len(values)
    return {
        "count": len(values),
        "mean_s-2": mean,
        "minimum_s-2": min(values),
        "maximum_s-2": max(values),
        "rms_s-2": math.sqrt(sum(value * value for value in values) / len(values)),
        "sum_s-2": sum(values),
    }


def durations(series: list[dict]) -> list[float]:
    return [
        float(series[index + 1]["time_day"] - series[index]["time_day"])
        for index in range(9)
    ]


def timing_summary(values: list[float]) -> dict:
    return {
        "count": len(values),
        "mean_day": sum(values) / len(values),
        "minimum_day": min(values),
        "maximum_day": max(values),
        "maximum_absolute_day": max(abs(value) for value in values),
    }


def analyze() -> dict:
    dependencies = check_dependencies()
    workbook = XlsxSource(WORKBOOK)
    try:
        cells_present = workbook.cells("a")
        cells_background = workbook.cells("b")
    finally:
        workbook.close()

    present_near, present_far = period_series(cells_present)
    background_near, background_far = period_series(cells_background)
    series = [present_near, present_far, background_near, background_far]
    labels = [
        *(f"P_N_{index}" for index in range(10)),
        *(f"P_F_{index}" for index in range(10)),
        *(f"B_N_{index}" for index in range(10)),
        *(f"B_F_{index}" for index in range(10)),
    ]
    endpoint_vector = [item["omega2_s-2"] for group in series for item in group]

    matrix_return = return_matrix()
    matrix_differential = differential_matrix()
    returns = matvec(matrix_return, endpoint_vector)
    differentials = matvec(matrix_differential, endpoint_vector)

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
    even = [(a + b) / 2.0 for a, b in zip(differentials[0:9], differentials[9:18])]
    odd = [(a - b) / 2.0 for a, b in zip(differentials[0:9], differentials[9:18])]
    duration_groups = {
        "source_present_N_F_N_days": durations(present_near),
        "source_present_F_N_F_days": durations(present_far),
        "source_absent_background_N_F_N_days": durations(background_near),
        "source_absent_background_F_N_F_days": durations(background_far),
    }
    duration_groups["present_minus_background_N_F_N_days"] = [
        a - b
        for a, b in zip(
            duration_groups["source_present_N_F_N_days"],
            duration_groups["source_absent_background_N_F_N_days"],
        )
    ]
    duration_groups["present_minus_background_F_N_F_days"] = [
        a - b
        for a, b in zip(
            duration_groups["source_present_F_N_F_days"],
            duration_groups["source_absent_background_F_N_F_days"],
        )
    ]

    result = {
        "schema": "WAC_HUST_TOS_ROUNDTRIP_HISTORY_RESIDUAL_V001",
        "date": "2026-08-27",
        "status": (
            "EQUAL_CONFIGURATION_ENDPOINT_RETURNS_EXTRACTED__"
            "PRESENT_MINUS_BACKGROUND_ORIENTATIONS_SEPARATED__"
            "EXACT_REUSED_ENDPOINT_WEIGHT_AND_OVERLAP_MATRICES__"
            "DESCRIPTIVE_ONLY_NO_LINEAGE_OR_COVERAGE_CLAIM"
        ),
        "source": {
            "workbook_relative_path": (
                "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/SOURCE/"
                "41586_2018_431_MOESM3_ESM.xlsx"
            ),
            "workbook_sha256": dependencies[
                "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/SOURCE/"
                "41586_2018_431_MOESM3_ESM.xlsx"
            ],
            "sheet_present": "a",
            "sheet_background": "b",
            "range_each": "B3:D22",
            "classification": "THREE_DAY_PERIOD_SUMMARIES__RAW_LIKE_INTERMEDIATE__NOT_EVENT_LEVEL",
        },
        "observable": {
            "transformed_endpoint": "y=(2*pi/T)^2 in s^-2",
            "N_F_N_return": "y_N[i+1]-y_N[i], with F[i] chronologically between",
            "F_N_F_return": "y_F[i+1]-y_F[i], with N[i+1] chronologically between",
            "panel_differential": "return_source_present-return_source_absent_background",
            "cross_panel_pairing": (
                "The index i pairs the same sequence ordinal in separate source-present and "
                "source-absent panels; it is not a simultaneous, randomized, or authenticated "
                "matched-run pairing."
            ),
            "orientation_even": "(d_NFN+d_FNF)/2",
            "orientation_odd": "(d_NFN-d_FNF)/2",
            "middle_point_role": "certifies the opposite-configuration excursion; has zero endpoint weight",
        },
        "inputs": {
            "endpoint_order": labels,
            "source_present_near": present_near,
            "source_present_far": present_far,
            "source_absent_background_near": background_near,
            "source_absent_background_far": background_far,
        },
        "counts": {
            "endpoint_omega2_values": 40,
            "N_F_N_per_panel": 9,
            "F_N_F_per_panel": 9,
            "all_panel_returns": 36,
            "present_minus_background_differentials": 18,
        },
        "timing": {
            "loop_duration_values": duration_groups,
            "loop_duration_summaries": {
                name: timing_summary(values) for name, values in duration_groups.items()
            },
            "ceiling": (
                "The primary endpoint-return residual is not divided by loop duration. "
                "Same-ordinal present/background loop durations differ and remain a descriptive "
                "limitation; the ordinal correspondence is not an experimental match."
            ),
        },
        "returns": {
            name: {"values_s-2": values, "summary": summary(values)}
            for name, values in return_groups.items()
        },
        "panel_differentials": {
            name: {"values_s-2": values, "summary": summary(values)}
            for name, values in differential_groups.items()
        },
        "orientation_decomposition": {
            "even_orientation_common_component": {
                "values_s-2": even,
                "summary": summary(even),
            },
            "odd_N_F_N_minus_F_N_F_half": {
                "values_s-2": odd,
                "summary": summary(odd),
            },
            "mean_full_orientation_difference_s-2": 2.0 * summary(odd)["mean_s-2"],
            "pairing_ceiling": (
                "The component means depend only on orientation means, but individual c_i/h_i "
                "values and the reported h_i RMS use an ordinal cross-panel pairing that is not "
                "a matched experimental trial."
            ),
        },
        "linear_custody": {
            "domain": "40 transformed omega^2 endpoints in endpoint_order",
            "return_row_order": [
                *(f"P_NFN_{index}" for index in range(9)),
                *(f"P_FNF_{index}" for index in range(9)),
                *(f"B_NFN_{index}" for index in range(9)),
                *(f"B_FNF_{index}" for index in range(9)),
            ],
            "differential_row_order": [
                *(f"D_NFN_{index}" for index in range(9)),
                *(f"D_FNF_{index}" for index in range(9)),
            ],
            "return_weight_matrix_36x40": matrix_return,
            "return_weight_rank_exact": exact_rank(matrix_return),
            "return_unit_endpoint_overlap_gram_36x36": gram(matrix_return),
            "differential_weight_matrix_18x40": matrix_differential,
            "differential_weight_rank_exact": exact_rank(matrix_differential),
            "differential_unit_endpoint_overlap_gram_18x18": gram(matrix_differential),
            "overlap_semantics": (
                "Exact row-weight Gram matrices under identity endpoint covariance; they expose reused "
                "endpoints and are not empirical covariance matrices."
            ),
        },
        "interpretation": {
            "descriptive_facts": [
                "Both panel-differential orientation means are negative at the released summary resolution.",
                "Their mean difference is much smaller than the individual loop variation.",
                "Adjacent returns reuse one endpoint with opposite sign; the 18 values are not independent.",
                "Each orientation mean telescopes to its last-minus-first endpoint divided by nine.",
                "Cross-panel loop indices are sequence ordinals, not matched experimental trials.",
            ],
            "strongest_claim": (
                "A reproducible figure-level equal-configuration return diagnostic and exact overlap "
                "ledger are available from the official HUST release."
            ),
            "counterfactual_ceiling": (
                "Every return follows an opposite-configuration excursion, but there is no matched "
                "no-excursion trajectory. The observable cannot identify hysteresis or memory "
                "against ordinary time drift and is a history-confound diagnostic only."
            ),
            "printed_precision_ceiling": (
                "Reported decimal digits reproduce arithmetic on the printed workbook periods. "
                "Digits beyond the source precision are computational digits, not measurement "
                "precision or an uncertainty model."
            ),
            "not_claimed": [
                "record-lineage intervention or beta_TM",
                "lineage gravitational charge or gravity emergence",
                "statistical significance, confidence interval, uncertainty, or coverage",
                "separation of fibre ageing, controller, thermal, mechanical, or source-motion memory",
                "a new G or a full GC16 execution",
            ],
        },
        "dependency_sha256": dependencies,
    }
    return result


def main() -> None:
    result = analyze()
    output = LANE / "RESULT.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("HUST_TOS_ROUNDTRIP_HISTORY_ANALYZER: PASS")
    print("N-F-N differential mean s^-2:", result["panel_differentials"]["present_minus_background_N_F_N"]["summary"]["mean_s-2"])
    print("F-N-F differential mean s^-2:", result["panel_differentials"]["present_minus_background_F_N_F"]["summary"]["mean_s-2"])
    print("orientation-odd half mean s^-2:", result["orientation_decomposition"]["odd_N_F_N_minus_F_N_F_half"]["summary"]["mean_s-2"])
    print("Ceiling: descriptive only; no covariance, lineage, beta_TM, or gravity-emergence claim")


if __name__ == "__main__":
    main()
