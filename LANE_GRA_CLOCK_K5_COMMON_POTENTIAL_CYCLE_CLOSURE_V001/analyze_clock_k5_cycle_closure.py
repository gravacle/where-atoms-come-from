#!/usr/bin/env python3
"""Exact K5 common-node-potential diagnostic for the official Fig4 CSV.

The result is a downstream processed-data compatibility check only.  It does
not treat shared-clock outputs as statistically independent and does not fit a
gravity model, record lineage, or an emergent metric.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


LANE = Path(__file__).resolve().parent
ROOT = LANE.parent

SOURCE_HASHES = {
    "SOURCE/Fig1.csv": "5ff66b2229f9c3b57d0fd2fa27e38aa773e9bfb597441c22393f9a65a36fed61",
    "SOURCE/Fig2.csv": "2df463bc8840f1b366d6acd9e361bf58322b844b72f6309641fc66e60735bc30",
    "SOURCE/Fig3.csv": "a4384040d5678893b4df15ea0be2980661023305f5009de876179fab1dea632f",
    "SOURCE/Fig4.csv": "3c14df355b5cf1e6dcf138cf3b3de750f59ad270091e7874c70ad05204fa988d",
}

DEPENDENCY_HASHES = {
    "../LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/SEARCH_CUSTODY.json":
        "a9924636c082fa3177aed69e6e332c0e2b0b26464335bc4b49cc14318c205b37",
    "../LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/PUBLIC_DATA_SECOND_PASS.md":
        "3d4300b9c2998aab4a485771f097f860e570a3931b8a948be9e1b034925931a8",
    "../LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "668abf36aea57e3465562465a6b77dff1ce65f0a4a41ae548fbea2f20a90c5e5",
}

EXPECTED_MD5 = {
    "SOURCE/Fig1.csv": "a829e2044f0ef2dd450435d7b790e8c7",
    "SOURCE/Fig2.csv": "8edb5d957e81a07b07e9d350afe43e3a",
    "SOURCE/Fig3.csv": "257eab6f29cf3480c536cac73ed3998a",
    "SOURCE/Fig4.csv": "53987b91b94d75d844f865ac0a778e75",
}

EXPECTED_BYTES = {
    "SOURCE/Fig1.csv": 25550,
    "SOURCE/Fig2.csv": 2876,
    "SOURCE/Fig3.csv": 958,
    "SOURCE/Fig4.csv": 369,
}

EDGE_ORDER = [
    (1, 2), (2, 3), (3, 4), (4, 5), (1, 3),
    (2, 4), (3, 5), (1, 4), (2, 5), (1, 5),
]


def file_digest(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_hash_file(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, relative = line.split(maxsplit=1)
        result[relative.strip()] = digest
    return result


def verify_custody() -> dict:
    if parse_hash_file(LANE / "SOURCE.sha256") != SOURCE_HASHES:
        raise RuntimeError("SOURCE.sha256 differs from frozen source contract")
    if parse_hash_file(LANE / "DEPENDENCIES.sha256") != DEPENDENCY_HASHES:
        raise RuntimeError("DEPENDENCIES.sha256 differs from frozen dependency contract")
    for relative, expected in SOURCE_HASHES.items():
        path = LANE / relative
        if path.stat().st_size != EXPECTED_BYTES[relative]:
            raise RuntimeError(f"source byte-size mismatch: {relative}")
        if file_digest(path) != expected:
            raise RuntimeError(f"source SHA-256 mismatch: {relative}")
        if file_digest(path, "md5") != EXPECTED_MD5[relative]:
            raise RuntimeError(f"source MD5 mismatch: {relative}")
    for relative, expected in DEPENDENCY_HASHES.items():
        if file_digest((LANE / relative).resolve()) != expected:
            raise RuntimeError(f"dependency mismatch: {relative}")
    return {
        "source_sha256": SOURCE_HASHES,
        "source_md5": EXPECTED_MD5,
        "source_bytes": EXPECTED_BYTES,
        "dependency_sha256": DEPENDENCY_HASHES,
    }


def read_fig4() -> tuple[list[Fraction], list[Fraction], list[str]]:
    with (LANE / "SOURCE/Fig4.csv").open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.reader(stream))
    if len(rows) != 3 or any(len(row) != 20 for row in rows):
        raise RuntimeError("Fig4.csv must be exactly 3 x 20")
    labels: list[str] = []
    values: list[Fraction] = []
    uncertainties: list[Fraction] = []
    for index, edge in enumerate(EDGE_ORDER):
        label = rows[0][2 * index]
        expected_label = f"({edge[0]}, {edge[1]})"
        if label != expected_label or rows[0][2 * index + 1] != "Uncertainty":
            raise RuntimeError(f"unexpected Fig4 header at edge {edge}")
        if rows[1][2 * index:2 * index + 2] != ["cm", "cm"]:
            raise RuntimeError(f"unexpected Fig4 units at edge {edge}")
        labels.append(label)
        values.append(Fraction(rows[2][2 * index]))
        uncertainties.append(Fraction(rows[2][2 * index + 1]))
    return values, uncertainties, labels


def incidence() -> list[list[int]]:
    matrix = [[0] * 5 for _ in EDGE_ORDER]
    for row, (i, j) in enumerate(EDGE_ORDER):
        matrix[row][i - 1] = 1
        matrix[row][j - 1] = -1
    return matrix


def transpose(matrix: list[list[Fraction | int]]) -> list[list[Fraction | int]]:
    return [list(column) for column in zip(*matrix)]


def matmul(left: list[list[Fraction | int]], right: list[list[Fraction | int]]) -> list[list[Fraction]]:
    right_t = transpose(right)
    return [
        [sum(Fraction(a) * Fraction(b) for a, b in zip(row, column)) for column in right_t]
        for row in left
    ]


def matvec(matrix: list[list[Fraction | int]], vector: list[Fraction]) -> list[Fraction]:
    return [sum(Fraction(a) * b for a, b in zip(row, vector)) for row in matrix]


def exact_rank(matrix: list[list[Fraction | int]]) -> int:
    work = [[Fraction(value) for value in row] for row in matrix]
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next((row for row in range(pivot_row, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [a - factor * b for a, b in zip(work[row], work[pivot_row])]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def solve_square(matrix: list[list[Fraction]], vector: list[Fraction]) -> list[Fraction] | None:
    size = len(matrix)
    work = [row[:] + [vector[index]] for index, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return None
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [a - factor * b for a, b in zip(work[row], work[column])]
    return [work[row][-1] for row in range(size)]


def edge_coefficient(path_start: int, path_end: int) -> tuple[int, int]:
    edge = (min(path_start, path_end), max(path_start, path_end))
    return EDGE_ORDER.index(edge), (1 if path_start < path_end else -1)


def canonical_cycles() -> list[tuple[int, ...]]:
    cycles: list[tuple[int, ...]] = []
    for length in (3, 4, 5):
        for subset in itertools.combinations(range(1, 6), length):
            first = min(subset)
            remaining = tuple(node for node in subset if node != first)
            for ordering in itertools.permutations(remaining):
                cycle = (first,) + ordering
                reverse = (first,) + tuple(reversed(ordering))
                if cycle <= reverse:
                    cycles.append(cycle)
    return cycles


def cycle_row(cycle: tuple[int, ...]) -> list[int]:
    row = [0] * len(EDGE_ORDER)
    for start, end in zip(cycle, cycle[1:] + cycle[:1]):
        index, sign = edge_coefficient(start, end)
        row[index] += sign
    return row


def fraction_record(value: Fraction) -> dict[str, str | float]:
    return {"exact": str(value), "decimal": float(value)}


def matrix_exact(matrix: list[list[Fraction | int]]) -> list[list[str]]:
    return [[str(Fraction(value)) for value in row] for row in matrix]


def chebyshev_box_fit(
    values: list[Fraction], uncertainties: list[Fraction], b_matrix: list[list[int]]
) -> tuple[Fraction, list[Fraction], list[Fraction]]:
    """Exact min_z max_e |y_e-(Bz)_e|/sigma_e with z_5 fixed to zero."""
    best: tuple[Fraction, list[Fraction], list[Fraction]] | None = None
    for active_edges in itertools.combinations(range(10), 5):
        for signs in itertools.product((-1, 1), repeat=5):
            system: list[list[Fraction]] = []
            target: list[Fraction] = []
            for edge_index, sign in zip(active_edges, signs):
                system.append([
                    *(Fraction(value) for value in b_matrix[edge_index][:4]),
                    Fraction(sign) * uncertainties[edge_index],
                ])
                target.append(values[edge_index])
            solution = solve_square(system, target)
            if solution is None:
                continue
            rho = solution[-1]
            if rho < 0 or (best is not None and rho >= best[0]):
                continue
            nodes = solution[:4] + [Fraction(0)]
            residuals = [
                values[index] - sum(Fraction(weight) * node for weight, node in zip(row, nodes))
                for index, row in enumerate(b_matrix)
            ]
            if all(abs(residual) <= rho * sigma for residual, sigma in zip(residuals, uncertainties)):
                best = (rho, nodes, residuals)
    if best is None:
        raise RuntimeError("exact Chebyshev fit found no feasible vertex")
    return best


def analyze() -> dict:
    custody = verify_custody()
    values, uncertainties, labels = read_fig4()
    b_matrix = incidence()
    b_t = transpose(b_matrix)
    rank_b = exact_rank(b_matrix)
    laplacian = matmul(b_t, b_matrix)

    # Exact unweighted Euclidean cut/cycle decomposition.  For K5, the
    # Laplacian pseudoinverse on the zero-sum gauge subspace is (1/5)I.
    bt_y = matvec(b_t, values)
    nodes_ls = [value / 5 for value in bt_y]
    cut_component = matvec(b_matrix, nodes_ls)
    cycle_component = [value - cut for value, cut in zip(values, cut_component)]
    p_cut = [[Fraction(value) / 5 for value in row] for row in matmul(b_matrix, b_t)]
    identity = [[Fraction(int(row == column)) for column in range(10)] for row in range(10)]
    p_cycle = [[identity[row][column] - p_cut[row][column] for column in range(10)] for row in range(10)]

    # Six star-tree fundamental triangles form a cycle-space basis.
    fundamental_cycles = [(1, i, j) for i, j in itertools.combinations(range(2, 6), 2)]
    fundamental_matrix = [cycle_row(cycle) for cycle in fundamental_cycles]

    all_cycles = canonical_cycles()
    cycle_records: list[dict] = []
    for cycle in all_cycles:
        row = cycle_row(cycle)
        residual = sum(Fraction(weight) * value for weight, value in zip(row, values))
        sigma_l1 = sum(abs(weight) * sigma for weight, sigma in zip(row, uncertainties))
        cycle_records.append({
            "length": len(cycle),
            "cycle": "-".join(str(node) for node in cycle + cycle[:1]),
            "nodes": list(cycle),
            "edge_coefficients_in_source_order": row,
            "residual_cm": fraction_record(residual),
            "marginal_1sigma_l1_envelope_cm": fraction_record(sigma_l1),
            "absolute_residual_over_l1_envelope": float(abs(residual) / sigma_l1),
            "individually_closeable_inside_marginal_1sigma_box": abs(residual) <= sigma_l1,
        })

    rho, nodes_box, residuals_box = chebyshev_box_fit(values, uncertainties, b_matrix)
    corrected_edges_box = [value - residual for value, residual in zip(values, residuals_box)]

    exact_closed = all(not value for value in matvec(fundamental_matrix, values))
    rounding_witness_cycle = (2, 3, 4)
    rounding_row = cycle_row(rounding_witness_cycle)
    rounding_residual = sum(Fraction(weight) * value for weight, value in zip(rounding_row, values))
    # Under conventional nearest rounding, each of the three displayed
    # 0.01-cm values can differ from its undisplayed precursor by at most
    # 0.005 cm.  This is a rounding-model statement, not a claim about the
    # metrological uncertainty or the source's unpublished digits.
    rounding_half_unit_sum = Fraction(3, 200)
    if abs(rounding_residual) <= rounding_half_unit_sum:
        raise RuntimeError("rounding witness no longer excludes a rounded common-node vector")

    output = {
        "schema": "WAC_CLOCK_K5_COMMON_POTENTIAL_CYCLE_CLOSURE_V001",
        "date": "2026-08-27",
        "status": (
            "PAIRWISE_ESTIMATED_NOT_COMMON_FIT_DERIVED__STATISTICAL_INDEPENDENCE_NOT_ESTABLISHED__"
            "K5_CUT4_CYCLE6_EXACT__MARGINAL_INTERVAL_COMPATIBLE__NO_GRAVITY_OR_LINEAGE_PROMOTION"
        ),
        "source_custody": custody,
        "source_semantics": {
            "Fig4_quantity": "ten pairwise inferred relative clock-height differences in cm",
            "mapping": "Delta h=(delta f/f)c^2/g, assuming GR and one independently measured local g",
            "article_method": (
                "Each pairwise clock comparison was re-analyzed individually over the same 14-run raw data, "
                "with pair-specific systematic corrections and a weighted mean."
            ),
            "provenance_classification": (
                "SEPARATELY_ESTIMATED_PAIRWISE_OUTPUTS__NOT_ALGEBRAICALLY_GENERATED_FROM_ONE_COMMON_NODE_FIT__"
                "SHARED_CLOCK_AND_RUN_COVARIANCE_UNOWNED"
            ),
            "test_scope": (
                "Because the height mapping is one common scalar multiplier, zero-cycle compatibility is "
                "equivalent for the deposited inferred heights and their underlying pairwise redshifts."
            ),
        },
        "edges": [
            {
                "source_index": index,
                "label": label,
                "orientation": f"{edge[0]}->{edge[1]}",
                "value_cm": fraction_record(value),
                "reported_1sigma_cm": fraction_record(sigma),
            }
            for index, (edge, label, value, sigma) in enumerate(zip(EDGE_ORDER, labels, values, uncertainties))
        ],
        "preflight": {
            "exact_common_node_representation": exact_closed,
            "decision": "PROCEED_WITH_BOUNDED_K5_DIAGNOSTIC",
            "rounding_exclusion_witness": {
                "cycle": "2-3-4-2",
                "residual_cm": fraction_record(rounding_residual),
                "maximum_sum_of_half_last_digit_rounding_cm": fraction_record(rounding_half_unit_sum),
                "residual_over_rounding_bound": float(abs(rounding_residual) / rounding_half_unit_sum),
                "rounding_model": "independent conventional nearest rounding to 0.01 cm",
                "conclusion": "NOT_A_COMMON_NODE_VECTOR_UNDER_INDEPENDENT_NEAREST_0.01_CM_ROUNDING",
            },
            "statistical_independence": "NOT_ESTABLISHED_AND_NOT_ASSUMED",
        },
        "k5_exact_algebra": {
            "edge_order": [list(edge) for edge in EDGE_ORDER],
            "incidence_10x5": b_matrix,
            "incidence_rank": rank_b,
            "cut_space_dimension": rank_b,
            "cycle_space_dimension": 10 - rank_b,
            "laplacian_5x5": matrix_exact(laplacian),
            "cut_projector_10x10": matrix_exact(p_cut),
            "cycle_projector_10x10": matrix_exact(p_cycle),
            "cut_projector_rank": exact_rank(p_cut),
            "cycle_projector_rank": exact_rank(p_cycle),
            "fundamental_cycle_order": [list(cycle) for cycle in fundamental_cycles],
            "fundamental_cycle_matrix_6x10": fundamental_matrix,
            "fundamental_cycle_rank": exact_rank(fundamental_matrix),
            "fundamental_cycles_annihilate_incidence": matrix_exact(matmul(fundamental_matrix, b_matrix)),
        },
        "unweighted_descriptive_projection": {
            "gauge": "sum(node_potential_cm)=0",
            "node_potentials_cm": [fraction_record(value) for value in nodes_ls],
            "cut_component_cm": [fraction_record(value) for value in cut_component],
            "cycle_component_cm": [fraction_record(value) for value in cycle_component],
            "cycle_component_l2_cm": math.sqrt(sum(float(value * value) for value in cycle_component)),
            "maximum_absolute_cycle_component_cm": max(float(abs(value)) for value in cycle_component),
            "ceiling": "Euclidean projection only; not covariance-weighted and not a chi-square fit.",
        },
        "simple_cycles": {
            "counts": {
                "triangles": sum(record["length"] == 3 for record in cycle_records),
                "quadrilaterals": sum(record["length"] == 4 for record in cycle_records),
                "pentagons": sum(record["length"] == 5 for record in cycle_records),
                "all": len(cycle_records),
            },
            "all_unique_up_to_rotation_and_reversal": cycle_records,
            "maximum_absolute_residual_cm": max(float(abs(Fraction(record["residual_cm"]["exact"]))) for record in cycle_records),
            "maximum_absolute_residual_over_marginal_l1_envelope": max(record["absolute_residual_over_l1_envelope"] for record in cycle_records),
            "all_individually_closeable_inside_marginal_1sigma_box": all(record["individually_closeable_inside_marginal_1sigma_box"] for record in cycle_records),
        },
        "covariance_honest_compatibility": {
            "optimization": "rho*=min_z max_e |y_e-(Bz)_e|/sigma_e with node_5=0",
            "rho_star_exact": str(rho),
            "rho_star_decimal": float(rho),
            "compatible_with_all_reported_marginal_1sigma_intervals": rho <= 1,
            "one_exact_feasible_node_vector_cm_node5_zero": [fraction_record(value) for value in nodes_box],
            "edge_residuals_cm": [fraction_record(value) for value in residuals_box],
            "edge_standardized_absolute_residuals": [float(abs(value) / sigma) for value, sigma in zip(residuals_box, uncertainties)],
            "corrected_cut_edges_cm": [fraction_record(value) for value in corrected_edges_box],
            "interpretation": (
                "The marginal intervals have a simultaneous geometric intersection with the K5 cut space. "
                "This is not a joint-coverage probability, p-value, confidence region, or chi-square test."
            ),
        },
        "strongest_claim": (
            "The ten processed pairwise redshift-derived heights are not algebraically forced to close, and "
            "they are geometrically compatible with one five-node scalar height vector inside every reported "
            "marginal 1sigma interval. Under the source's already-assumed GR mapping and common g, that node "
            "scalar is proportional to a relative gravitational-potential vector."
        ),
        "not_claimed": [
            "statistical independence of the ten shared-clock pair estimates",
            "conventional chi-square, p-value, confidence region, or joint coverage",
            "an independent test of GR, gravity, a common metric, or local g",
            "record lineage, beta_TM, Gravity Formation Theory, or gravity emergence",
            "raw-data or run-wise closure; those vectors and covariance were not deposited",
        ],
    }
    return output


def main() -> None:
    result = analyze()
    (LANE / "RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    compatibility = result["covariance_honest_compatibility"]
    cycles = result["simple_cycles"]
    print("CLOCK_K5_COMMON_POTENTIAL_ANALYZER: PASS")
    print("preflight:", result["source_semantics"]["provenance_classification"])
    print("cut_rank: 4")
    print("cycle_dimension: 6")
    print("simple_cycles_checked:", cycles["counts"]["all"])
    print("rho_star:", compatibility["rho_star_decimal"])
    print("marginal_interval_compatible:", compatibility["compatible_with_all_reported_marginal_1sigma_intervals"])
    print("ceiling: PROCESSED_PAIRWISE_SCALAR_COMPATIBILITY_ONLY__NO_CHI_SQUARE_GRAVITY_METRIC_LINEAGE_OR_EMERGENCE_CLAIM")


if __name__ == "__main__":
    main()
