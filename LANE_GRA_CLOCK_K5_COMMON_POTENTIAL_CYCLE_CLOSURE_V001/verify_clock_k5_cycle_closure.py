#!/usr/bin/env python3
"""Independent hostile-style verifier for GRA-CLOCK-K5-CPC-V001.

The verifier imports no production analyzer code.  It reparses every official
CSV, independently rebuilds the exact K5 algebra and marginal-box optimum, and
checks the sealed result and its scientific ceiling.
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

SOURCE_CONTRACT = {
    "SOURCE/Fig1.csv": (25550, "a829e2044f0ef2dd450435d7b790e8c7", "5ff66b2229f9c3b57d0fd2fa27e38aa773e9bfb597441c22393f9a65a36fed61"),
    "SOURCE/Fig2.csv": (2876, "8edb5d957e81a07b07e9d350afe43e3a", "2df463bc8840f1b366d6acd9e361bf58322b844b72f6309641fc66e60735bc30"),
    "SOURCE/Fig3.csv": (958, "257eab6f29cf3480c536cac73ed3998a", "a4384040d5678893b4df15ea0be2980661023305f5009de876179fab1dea632f"),
    "SOURCE/Fig4.csv": (369, "53987b91b94d75d844f865ac0a778e75", "3c14df355b5cf1e6dcf138cf3b3de750f59ad270091e7874c70ad05204fa988d"),
}

DEPENDENCY_CONTRACT = {
    "../LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/SEARCH_CUSTODY.json": "a9924636c082fa3177aed69e6e332c0e2b0b26464335bc4b49cc14318c205b37",
    "../LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/PUBLIC_DATA_SECOND_PASS.md": "3d4300b9c2998aab4a485771f097f860e570a3931b8a948be9e1b034925931a8",
    "../LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/INDEPENDENT_HOSTILE_AUDIT.md": "668abf36aea57e3465562465a6b77dff1ce65f0a4a41ae548fbea2f20a90c5e5",
}

EDGES = [
    (1, 2), (2, 3), (3, 4), (4, 5), (1, 3),
    (2, 4), (3, 5), (1, 4), (2, 5), (1, 5),
]

MANIFEST_MEMBERS = {
    "DEPENDENCIES.sha256",
    "README.md",
    "RESULT.json",
    "RESULT.md",
    "SELF_AUDIT.md",
    "SOURCE.sha256",
    "SOURCE/Fig1.csv",
    "SOURCE/Fig2.csv",
    "SOURCE/Fig3.csv",
    "SOURCE/Fig4.csv",
    "SOURCE_CUSTODY.json",
    "THEOREM.md",
    "VERIFICATION.txt",
    "analyze_clock_k5_cycle_closure.py",
    "hostile_audit_clock_k5_cycle_closure.py",
    "HOSTILE_AUDIT_TRANSCRIPT.txt",
    "INDEPENDENT_HOSTILE_AUDIT.md",
    "verify_clock_k5_cycle_closure.py",
}


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
            math.isclose(float(actual), float(expected), rel_tol=2e-14, abs_tol=2e-14),
            f"{label}: {actual!r} != {expected!r}",
        )


def digest(path: Path, algorithm: str) -> str:
    state = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 18), b""):
            state.update(chunk)
    return state.hexdigest()


def parse_hashes(path: Path) -> dict[str, str]:
    output: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            value, name = line.split(maxsplit=1)
            output[name.strip()] = value
    return output


def read_csv(name: str) -> list[list[str]]:
    with (LANE / "SOURCE" / name).open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.reader(stream))


def rank_exact(matrix: list[list[Fraction | int]]) -> int:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for index in range(column, len(work[rank])):
            work[rank][index] /= pivot_value
        for row in range(len(work)):
            if row == rank:
                continue
            factor = work[row][column]
            if factor:
                for index in range(column, len(work[row])):
                    work[row][index] -= factor * work[rank][index]
        rank += 1
        if rank == len(work):
            break
    return rank


def dot(left: list[Fraction | int], right: list[Fraction | int]) -> Fraction:
    return sum(Fraction(a) * Fraction(b) for a, b in zip(left, right))


def product(left: list[list[Fraction | int]], right: list[list[Fraction | int]]) -> list[list[Fraction]]:
    columns = [list(column) for column in zip(*right)]
    return [[dot(row, column) for column in columns] for row in left]


def apply(matrix: list[list[Fraction | int]], vector: list[Fraction]) -> list[Fraction]:
    return [dot(row, vector) for row in matrix]


def solve(matrix: list[list[Fraction]], target: list[Fraction]) -> list[Fraction] | None:
    size = len(matrix)
    augmented = [matrix[row][:] + [target[row]] for row in range(size)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if augmented[row][column]), None)
        if pivot is None:
            return None
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [entry / scale for entry in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [a - factor * b for a, b in zip(augmented[row], augmented[column])]
    return [augmented[row][-1] for row in range(size)]


def build_incidence() -> list[list[int]]:
    output = [[0 for _ in range(5)] for _ in range(10)]
    for row, (start, end) in enumerate(EDGES):
        output[row][start - 1] = 1
        output[row][end - 1] = -1
    return output


def coefficient_for_step(start: int, end: int) -> tuple[int, int]:
    edge = (min(start, end), max(start, end))
    return EDGES.index(edge), 1 if start < end else -1


def weights_for_cycle(cycle: tuple[int, ...]) -> list[int]:
    output = [0] * 10
    closed = cycle + cycle[:1]
    for index in range(len(cycle)):
        edge_index, sign = coefficient_for_step(closed[index], closed[index + 1])
        output[edge_index] += sign
    return output


def enumerate_cycles() -> list[tuple[int, ...]]:
    found: set[tuple[int, ...]] = set()
    for length in (3, 4, 5):
        for path in itertools.permutations(range(1, 6), length):
            if len(set(path)) != length:
                continue
            rotations = [path[index:] + path[:index] for index in range(length)]
            reverse = tuple(reversed(path))
            rotations += [reverse[index:] + reverse[:index] for index in range(length)]
            found.add(min(rotations))
    return sorted(found, key=lambda cycle: (len(cycle), cycle))


def independent_box_optimum(
    values: list[Fraction], sigmas: list[Fraction], incidence: list[list[int]]
) -> tuple[Fraction, list[Fraction], list[Fraction]]:
    winner: tuple[Fraction, list[Fraction], list[Fraction]] | None = None
    # Five variables: four node coordinates (node 5 is gauge-fixed) and rho.
    for chosen in itertools.combinations(range(10), 5):
        for signs in itertools.product((1, -1), repeat=5):
            matrix = [
                [*(Fraction(value) for value in incidence[edge][:4]), Fraction(sign) * sigmas[edge]]
                for edge, sign in zip(chosen, signs)
            ]
            target = [values[edge] for edge in chosen]
            candidate = solve(matrix, target)
            if candidate is None or candidate[-1] < 0:
                continue
            rho = candidate[-1]
            if winner is not None and rho >= winner[0]:
                continue
            nodes = candidate[:4] + [Fraction(0)]
            residuals = [values[row] - dot(incidence[row], nodes) for row in range(10)]
            if all(abs(residual) <= rho * sigma for residual, sigma in zip(residuals, sigmas)):
                winner = rho, nodes, residuals
    if winner is None:
        raise AssertionError("box optimum not found")
    return winner


def check_fraction_record(audit: Audit, record: dict, expected: Fraction, label: str) -> None:
    audit.equal(record["exact"], str(expected), f"{label} exact")
    audit.close(record["decimal"], float(expected), f"{label} decimal")


def verify_manifest(audit: Audit) -> None:
    manifest = parse_hashes(LANE / "MANIFEST.sha256")
    audit.equal(set(manifest), MANIFEST_MEMBERS, "manifest members")
    for name in sorted(MANIFEST_MEMBERS):
        audit.equal(digest(LANE / name, "sha256"), manifest[name], f"manifest {name}")
    seal = parse_hashes(LANE / "LANE_SEAL.sha256")
    audit.equal(set(seal), {"MANIFEST.sha256"}, "seal member")
    audit.equal(digest(LANE / "MANIFEST.sha256", "sha256"), seal["MANIFEST.sha256"], "manifest seal")
    result_bytes = bytearray((LANE / "RESULT.json").read_bytes())
    result_bytes[len(result_bytes) // 2] ^= 1
    audit.require(hashlib.sha256(result_bytes).hexdigest() != manifest["RESULT.json"], "result tamper sentinel")
    source_bytes = bytearray((LANE / "SOURCE/Fig4.csv").read_bytes())
    source_bytes[-2] ^= 1
    audit.require(hashlib.sha256(source_bytes).hexdigest() != SOURCE_CONTRACT["SOURCE/Fig4.csv"][2], "source tamper sentinel")


def main() -> None:
    audit = Audit()

    declared_sources = parse_hashes(LANE / "SOURCE.sha256")
    audit.equal(set(declared_sources), set(SOURCE_CONTRACT), "source hash members")
    for name, (size, md5, sha) in SOURCE_CONTRACT.items():
        path = LANE / name
        audit.equal(path.stat().st_size, size, f"{name} bytes")
        audit.equal(digest(path, "md5"), md5, f"{name} md5")
        audit.equal(digest(path, "sha256"), sha, f"{name} sha256")
        audit.equal(declared_sources[name], sha, f"{name} declared sha256")

    declared_dependencies = parse_hashes(LANE / "DEPENDENCIES.sha256")
    audit.equal(declared_dependencies, DEPENDENCY_CONTRACT, "dependency declarations")
    for relative, expected in DEPENDENCY_CONTRACT.items():
        audit.equal(digest((LANE / relative).resolve(), "sha256"), expected, f"dependency {relative}")

    csv_rows = {name: read_csv(name) for name in ("Fig1.csv", "Fig2.csv", "Fig3.csv", "Fig4.csv")}
    audit.equal((len(csv_rows["Fig1.csv"]), len(csv_rows["Fig1.csv"][0])), (721, 5), "Fig1 shape")
    audit.equal((len(csv_rows["Fig2.csv"]), len(csv_rows["Fig2.csv"][0])), (12, 17), "Fig2 shape")
    audit.equal((len(csv_rows["Fig3.csv"]), len(csv_rows["Fig3.csv"][0])), (16, 6), "Fig3 shape")
    audit.equal((len(csv_rows["Fig4.csv"]), len(csv_rows["Fig4.csv"][0])), (3, 20), "Fig4 shape")
    numeric_counts = []
    for column in range(5):
        count = 0
        for row in csv_rows["Fig1.csv"]:
            try:
                float(row[column])
                count += 1
            except (ValueError, IndexError):
                pass
        numeric_counts.append(count)
    audit.equal(numeric_counts, [714, 719, 717, 716, 716], "Fig1 numeric counts")

    fig4 = csv_rows["Fig4.csv"]
    values: list[Fraction] = []
    sigmas: list[Fraction] = []
    for edge_index, edge in enumerate(EDGES):
        audit.equal(fig4[0][2 * edge_index], f"({edge[0]}, {edge[1]})", f"edge {edge} label")
        audit.equal(fig4[0][2 * edge_index + 1], "Uncertainty", f"edge {edge} uncertainty header")
        audit.equal(fig4[1][2 * edge_index:2 * edge_index + 2], ["cm", "cm"], f"edge {edge} units")
        values.append(Fraction(fig4[2][2 * edge_index]))
        sigmas.append(Fraction(fig4[2][2 * edge_index + 1]))
    audit.equal(values, [Fraction(text) for text in ("0.25", "0.43", "0.19", "0.23", "0.66", "0.53", "0.57", "0.73", "0.82", "1.1")], "Fig4 values")
    audit.equal(sigmas, [Fraction(text) for text in ("0.18", "0.18", "0.18", "0.19", "0.2", "0.2", "0.2", "0.23", "0.23", "0.27")], "Fig4 sigmas")

    payload = json.loads((LANE / "RESULT.json").read_text(encoding="utf-8"))
    audit.equal(payload["schema"], "WAC_CLOCK_K5_COMMON_POTENTIAL_CYCLE_CLOSURE_V001", "schema")
    audit.equal(payload["date"], "2026-08-27", "date")
    audit.equal(payload["preflight"]["exact_common_node_representation"], False, "not exact common-node vector")
    audit.equal(payload["preflight"]["statistical_independence"], "NOT_ESTABLISHED_AND_NOT_ASSUMED", "independence ceiling")
    audit.require("SHARED_CLOCK_AND_RUN_COVARIANCE_UNOWNED" in payload["source_semantics"]["provenance_classification"], "provenance covariance ceiling")
    for index, (edge, value, sigma) in enumerate(zip(EDGES, values, sigmas)):
        record = payload["edges"][index]
        audit.equal(record["source_index"], index, f"edge record {index} index")
        audit.equal(record["orientation"], f"{edge[0]}->{edge[1]}", f"edge record {index} orientation")
        check_fraction_record(audit, record["value_cm"], value, f"edge record {index} value")
        check_fraction_record(audit, record["reported_1sigma_cm"], sigma, f"edge record {index} sigma")

    b = build_incidence()
    audit.equal(rank_exact(b), 4, "independent incidence rank")
    b_t = [list(column) for column in zip(*b)]
    laplacian = product(b_t, b)
    audit.equal(laplacian, [[Fraction(4 if row == col else -1) for col in range(5)] for row in range(5)], "K5 Laplacian")
    p_cut = [[entry / 5 for entry in row] for row in product(b, b_t)]
    p_cycle = [[Fraction(int(row == col)) - p_cut[row][col] for col in range(10)] for row in range(10)]
    audit.equal(rank_exact(p_cut), 4, "independent cut projector rank")
    audit.equal(rank_exact(p_cycle), 6, "independent cycle projector rank")
    audit.equal(product(p_cut, p_cut), p_cut, "cut projector idempotence")
    audit.equal(product(p_cycle, p_cycle), p_cycle, "cycle projector idempotence")
    audit.equal(product(p_cut, p_cycle), [[Fraction(0)] * 10 for _ in range(10)], "projector orthogonality")

    algebra = payload["k5_exact_algebra"]
    audit.equal(algebra["edge_order"], [list(edge) for edge in EDGES], "result edge order")
    audit.equal(algebra["incidence_10x5"], b, "result incidence")
    for key, expected in (("incidence_rank", 4), ("cut_space_dimension", 4), ("cycle_space_dimension", 6), ("cut_projector_rank", 4), ("cycle_projector_rank", 6), ("fundamental_cycle_rank", 6)):
        audit.equal(algebra[key], expected, f"algebra {key}")
    audit.equal(algebra["laplacian_5x5"], [[str(value) for value in row] for row in laplacian], "result Laplacian")
    audit.equal(algebra["cut_projector_10x10"], [[str(value) for value in row] for row in p_cut], "result cut projector")
    audit.equal(algebra["cycle_projector_10x10"], [[str(value) for value in row] for row in p_cycle], "result cycle projector")

    fundamental = [weights_for_cycle((1, start, end)) for start, end in itertools.combinations(range(2, 6), 2)]
    audit.equal(rank_exact(fundamental), 6, "independent fundamental-cycle rank")
    audit.equal(product(fundamental, b), [[Fraction(0)] * 5 for _ in range(6)], "cycles annihilate incidence")
    audit.equal(algebra["fundamental_cycle_matrix_6x10"], fundamental, "result fundamental matrix")

    bt_y = apply(b_t, values)
    nodes_ls = [entry / 5 for entry in bt_y]
    cut = apply(b, nodes_ls)
    cycle_part = [value - fitted for value, fitted in zip(values, cut)]
    projection = payload["unweighted_descriptive_projection"]
    for index, expected in enumerate(nodes_ls):
        check_fraction_record(audit, projection["node_potentials_cm"][index], expected, f"LS node {index}")
    for index, expected in enumerate(cut):
        check_fraction_record(audit, projection["cut_component_cm"][index], expected, f"cut edge {index}")
    for index, expected in enumerate(cycle_part):
        check_fraction_record(audit, projection["cycle_component_cm"][index], expected, f"cycle edge {index}")
    audit.close(projection["cycle_component_l2_cm"], math.sqrt(sum(float(value * value) for value in cycle_part)), "cycle L2")
    audit.close(projection["maximum_absolute_cycle_component_cm"], max(float(abs(value)) for value in cycle_part), "cycle max component")

    cycles = enumerate_cycles()
    audit.equal(len(cycles), 37, "unique simple-cycle count")
    audit.equal(sum(len(cycle) == 3 for cycle in cycles), 10, "triangle count")
    audit.equal(sum(len(cycle) == 4 for cycle in cycles), 15, "quadrilateral count")
    audit.equal(sum(len(cycle) == 5 for cycle in cycles), 12, "pentagon count")
    stored_cycles = payload["simple_cycles"]["all_unique_up_to_rotation_and_reversal"]
    audit.equal(len(stored_cycles), 37, "stored cycle count")
    stored_by_nodes = {tuple(record["nodes"]): record for record in stored_cycles}
    audit.equal(set(stored_by_nodes), set(cycles), "stored cycle identity set")
    max_abs = Fraction(0)
    max_ratio = Fraction(0)
    for index, cycle in enumerate(cycles):
        row = weights_for_cycle(cycle)
        residual = dot(row, values)
        envelope = sum(abs(weight) * sigma for weight, sigma in zip(row, sigmas))
        max_abs = max(max_abs, abs(residual))
        max_ratio = max(max_ratio, abs(residual) / envelope)
        record = stored_by_nodes[cycle]
        audit.equal(record["nodes"], list(cycle), f"cycle {index} nodes")
        audit.equal(record["cycle"], "-".join(str(node) for node in cycle + cycle[:1]), f"cycle {index} name")
        audit.equal(record["edge_coefficients_in_source_order"], row, f"cycle {index} row")
        check_fraction_record(audit, record["residual_cm"], residual, f"cycle {index} residual")
        check_fraction_record(audit, record["marginal_1sigma_l1_envelope_cm"], envelope, f"cycle {index} envelope")
        audit.close(record["absolute_residual_over_l1_envelope"], float(abs(residual) / envelope), f"cycle {index} ratio")
        audit.equal(record["individually_closeable_inside_marginal_1sigma_box"], abs(residual) <= envelope, f"cycle {index} interval closeable")
    audit.equal(max_abs, Fraction(29, 100), "maximum exact cycle residual")
    audit.equal(max_ratio, Fraction(27, 82), "maximum exact cycle ratio")
    audit.close(payload["simple_cycles"]["maximum_absolute_residual_cm"], 0.29, "reported maximum residual")
    audit.close(payload["simple_cycles"]["maximum_absolute_residual_over_marginal_l1_envelope"], float(Fraction(27, 82)), "reported maximum ratio")
    audit.equal(payload["simple_cycles"]["all_individually_closeable_inside_marginal_1sigma_box"], True, "all cycles marginally closeable")

    witness_row = weights_for_cycle((2, 3, 4))
    witness = dot(witness_row, values)
    audit.equal(witness, Fraction(9, 100), "rounding witness residual")
    audit.require(witness > Fraction(3, 200), "rounding witness exceeds bound")
    witness_record = payload["preflight"]["rounding_exclusion_witness"]
    check_fraction_record(audit, witness_record["residual_cm"], Fraction(9, 100), "reported rounding witness")
    check_fraction_record(audit, witness_record["maximum_sum_of_half_last_digit_rounding_cm"], Fraction(3, 200), "reported rounding bound")
    audit.close(witness_record["residual_over_rounding_bound"], 6.0, "reported rounding ratio")
    audit.equal(
        witness_record["rounding_model"],
        "independent conventional nearest rounding to 0.01 cm",
        "reported rounding model",
    )

    rho, nodes_box, residuals_box = independent_box_optimum(values, sigmas, b)
    audit.equal(rho, Fraction(27, 82), "exact box optimum")
    audit.require(all(abs(residual) <= rho * sigma for residual, sigma in zip(residuals_box, sigmas)), "box solution feasible")
    compatibility = payload["covariance_honest_compatibility"]
    audit.equal(compatibility["rho_star_exact"], "27/82", "reported rho exact")
    audit.close(compatibility["rho_star_decimal"], float(rho), "reported rho decimal")
    audit.equal(compatibility["compatible_with_all_reported_marginal_1sigma_intervals"], True, "reported interval compatibility")
    for index, expected in enumerate(nodes_box):
        check_fraction_record(audit, compatibility["one_exact_feasible_node_vector_cm_node5_zero"][index], expected, f"box node {index}")
    for index, expected in enumerate(residuals_box):
        check_fraction_record(audit, compatibility["edge_residuals_cm"][index], expected, f"box residual {index}")
        audit.close(compatibility["edge_standardized_absolute_residuals"][index], float(abs(expected) / sigmas[index]), f"box standardized residual {index}")
        check_fraction_record(audit, compatibility["corrected_cut_edges_cm"][index], values[index] - expected, f"box corrected edge {index}")

    custody = json.loads((LANE / "SOURCE_CUSTODY.json").read_text(encoding="utf-8"))
    audit.equal(custody["repository"]["doi"], "10.5281/zenodo.8184043", "custody DOI")
    audit.equal(custody["repository"]["license_in_api_metadata"], "cc-by-4.0", "custody license")
    audit.equal(custody["spreadsheet_inspection"]["all_four_csvs_inspected_and_rendered"], True, "artifact-tool inspection")
    audit.equal(custody["spreadsheet_inspection"]["source_bytes_modified"], False, "source bytes unmodified")
    for entry in custody["files"]:
        contract = SOURCE_CONTRACT[f"SOURCE/{entry['name']}"]
        audit.equal(entry["bytes"], contract[0], f"custody {entry['name']} bytes")
        audit.equal(entry["published_md5"], contract[1], f"custody {entry['name']} md5")
        audit.equal(entry["sha256"], contract[2], f"custody {entry['name']} sha256")

    theorem = (LANE / "THEOREM.md").read_text(encoding="utf-8")
    result_md = (LANE / "RESULT.md").read_text(encoding="utf-8")
    self_audit = (LANE / "SELF_AUDIT.md").read_text(encoding="utf-8")
    corpus = " ".join((theorem + result_md + self_audit).replace("**", "").lower().split())
    for phrase in (
        "not statistically independent",
        "no joint covariance",
        "not a chi-square",
        "not an independent confirmation of gr",
        "necessary, never sufficient",
        "record lineage",
        "gravity emergence",
        "no canonical model",
    ):
        audit.require(phrase in corpus, f"required ceiling phrase: {phrase}")
    forbidden_keys = {
        "chi_square", "chi2", "p_value", "confidence_interval", "joint_coverage",
        "beta_tm_estimate", "metric_confirmation", "gravity_emergence_confirmation",
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
    audit.require(not (keys & forbidden_keys), "no unauthorized inferential result fields")

    verify_manifest(audit)

    print("CLOCK_K5_COMMON_POTENTIAL_VERIFIER: PASS")
    print(f"checks_passed: {audit.count}")
    print("official_four_csv_custody: PASS")
    print("preflight_not_common_fit_derived: PASS")
    print("statistical_independence_not_assumed: PASS")
    print("K5_cut_rank_cycle_dimension: 4 6")
    print("all_simple_cycles: 37")
    print("rho_star_exact: 27/82")
    print("marginal_interval_compatibility: PASS")
    print("tamper_sentinels: PASS")
    print("ceiling: NECESSARY_PROCESSED_COMMON_NODE_SCALAR_CHECK_ONLY__NO_CHI_SQUARE_GR_GRAVITY_METRIC_LINEAGE_OR_EMERGENCE_CLAIM")


if __name__ == "__main__":
    main()
