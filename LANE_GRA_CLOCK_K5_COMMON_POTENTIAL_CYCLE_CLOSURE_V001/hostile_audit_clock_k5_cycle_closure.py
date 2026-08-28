#!/usr/bin/env python3
"""Independent hostile audit of GRA-CLOCK-K5-CPC-V001.

This executable imports neither the production analyzer nor its verifier.  It
reparses the pinned CSV bytes, rebuilds the graph algebra, enumerates every K5
simple cycle, and certifies the Chebyshev optimum by a cycle-dual lower bound
plus an independently constructed exact feasible node vector.
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

EDGES = (
    (1, 2), (2, 3), (3, 4), (4, 5), (1, 3),
    (2, 4), (3, 5), (1, 4), (2, 5), (1, 5),
)

SOURCE_CONTRACT = {
    "SOURCE/Fig1.csv": (721, 5, 25550, "a829e2044f0ef2dd450435d7b790e8c7", "5ff66b2229f9c3b57d0fd2fa27e38aa773e9bfb597441c22393f9a65a36fed61"),
    "SOURCE/Fig2.csv": (12, 17, 2876, "8edb5d957e81a07b07e9d350afe43e3a", "2df463bc8840f1b366d6acd9e361bf58322b844b72f6309641fc66e60735bc30"),
    "SOURCE/Fig3.csv": (16, 6, 958, "257eab6f29cf3480c536cac73ed3998a", "a4384040d5678893b4df15ea0be2980661023305f5009de876179fab1dea632f"),
    "SOURCE/Fig4.csv": (3, 20, 369, "53987b91b94d75d844f865ac0a778e75", "3c14df355b5cf1e6dcf138cf3b3de750f59ad270091e7874c70ad05204fa988d"),
}

DEPENDENCY_CONTRACT = {
    "../LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/SEARCH_CUSTODY.json": "a9924636c082fa3177aed69e6e332c0e2b0b26464335bc4b49cc14318c205b37",
    "../LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/PUBLIC_DATA_SECOND_PASS.md": "3d4300b9c2998aab4a485771f097f860e570a3931b8a948be9e1b034925931a8",
    "../LANE_GRA_SPAG_PUBLIC_DATA_SECOND_PASS_V001/INDEPENDENT_HOSTILE_AUDIT.md": "668abf36aea57e3465562465a6b77dff1ce65f0a4a41ae548fbea2f20a90c5e5",
}

MANIFEST_MEMBERS = {
    "DEPENDENCIES.sha256",
    "HOSTILE_AUDIT_TRANSCRIPT.txt",
    "INDEPENDENT_HOSTILE_AUDIT.md",
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


def digest(path: Path, algorithm: str = "sha256") -> str:
    state = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 18), b""):
            state.update(block)
    return state.hexdigest()


def parse_hashes(path: Path) -> dict[str, str]:
    output: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            value, name = line.split(maxsplit=1)
            output[name.strip()] = value
    return output


def read_csv(relative: str) -> list[list[str]]:
    with (LANE / relative).open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.reader(stream))


def exact_rank(matrix: list[list[int | Fraction]]) -> int:
    work = [[Fraction(value) for value in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][column]
        work[rank] = [value / scale for value in work[rank]]
        for row in range(len(work)):
            if row == rank:
                continue
            scale = work[row][column]
            if scale:
                work[row] = [a - scale * b for a, b in zip(work[row], work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def incidence() -> list[list[int]]:
    matrix = [[0] * 5 for _ in EDGES]
    for row, (start, end) in enumerate(EDGES):
        matrix[row][start - 1] = 1
        matrix[row][end - 1] = -1
    return matrix


def dot(left: list[int | Fraction], right: list[int | Fraction]) -> Fraction:
    return sum(Fraction(a) * Fraction(b) for a, b in zip(left, right))


def matmul(
    left: list[list[int | Fraction]], right: list[list[int | Fraction]]
) -> list[list[Fraction]]:
    columns = [list(column) for column in zip(*right)]
    return [[dot(row, column) for column in columns] for row in left]


def cycle_canonical(path: tuple[int, ...]) -> tuple[int, ...]:
    rotations = [path[index:] + path[:index] for index in range(len(path))]
    reverse = tuple(reversed(path))
    rotations.extend(reverse[index:] + reverse[:index] for index in range(len(path)))
    return min(rotations)


def all_simple_cycles() -> list[tuple[int, ...]]:
    cycles: set[tuple[int, ...]] = set()
    for length in (3, 4, 5):
        for subset in itertools.combinations(range(1, 6), length):
            for ordering in itertools.permutations(subset):
                cycles.add(cycle_canonical(ordering))
    return sorted(cycles, key=lambda cycle: (len(cycle), cycle))


def cycle_row(cycle: tuple[int, ...]) -> list[int]:
    edge_lookup = {edge: index for index, edge in enumerate(EDGES)}
    output = [0] * len(EDGES)
    closed = cycle + cycle[:1]
    for start, end in zip(closed, closed[1:]):
        undirected = (min(start, end), max(start, end))
        output[edge_lookup[undirected]] += 1 if start < end else -1
    return output


def exact_feasible_nodes_by_difference_constraints(
    values: list[Fraction], sigmas: list[Fraction], rho: Fraction
) -> list[Fraction]:
    """Construct one exact feasible node vector without vertex enumeration."""
    constraints: list[tuple[int, int, Fraction]] = []
    for (start, end), value, sigma in zip(EDGES, values, sigmas):
        # z_start - z_end <= value + rho*sigma gives end -> start.
        constraints.append((end - 1, start - 1, value + rho * sigma))
        # z_start - z_end >= value - rho*sigma gives start -> end.
        constraints.append((start - 1, end - 1, -value + rho * sigma))

    # A zero-weight super-source is equivalent to initializing all distances
    # to zero.  Four complete relaxation passes suffice for five nodes when no
    # negative cycle exists.
    distance = [Fraction(0)] * 5
    for _ in range(4):
        for source, target, weight in constraints:
            distance[target] = min(distance[target], distance[source] + weight)
    if any(distance[target] > distance[source] + weight for source, target, weight in constraints):
        raise AssertionError("difference constraints have a negative cycle")
    gauge = distance[4]
    return [value - gauge for value in distance]


def fraction_record(audit: Audit, record: dict, expected: Fraction, label: str) -> None:
    audit.equal(record["exact"], str(expected), f"{label} exact")
    audit.close(record["decimal"], float(expected), f"{label} decimal")


def main() -> None:
    audit = Audit()

    declared_sources = parse_hashes(LANE / "SOURCE.sha256")
    audit.equal(set(declared_sources), set(SOURCE_CONTRACT), "source declaration members")
    parsed: dict[str, list[list[str]]] = {}
    for relative, (rows, columns, byte_size, md5, sha256) in SOURCE_CONTRACT.items():
        path = LANE / relative
        audit.equal(path.stat().st_size, byte_size, f"{relative} byte size")
        audit.equal(digest(path, "md5"), md5, f"{relative} MD5")
        audit.equal(digest(path), sha256, f"{relative} SHA-256")
        audit.equal(declared_sources[relative], sha256, f"{relative} declared SHA-256")
        table = read_csv(relative)
        parsed[relative] = table
        audit.equal(len(table), rows, f"{relative} row count")
        audit.require(all(len(row) == columns for row in table), f"{relative} rectangular {columns} columns")

    declared_dependencies = parse_hashes(LANE / "DEPENDENCIES.sha256")
    audit.equal(declared_dependencies, DEPENDENCY_CONTRACT, "dependency declaration")
    for relative, expected in DEPENDENCY_CONTRACT.items():
        audit.equal(digest((LANE / relative).resolve()), expected, f"dependency {relative}")

    fig4 = parsed["SOURCE/Fig4.csv"]
    values: list[Fraction] = []
    sigmas: list[Fraction] = []
    for index, edge in enumerate(EDGES):
        audit.equal(fig4[0][2 * index], f"({edge[0]}, {edge[1]})", f"Fig4 label {edge}")
        audit.equal(fig4[0][2 * index + 1], "Uncertainty", f"Fig4 uncertainty label {edge}")
        audit.equal(fig4[1][2 * index:2 * index + 2], ["cm", "cm"], f"Fig4 units {edge}")
        values.append(Fraction(fig4[2][2 * index]))
        sigmas.append(Fraction(fig4[2][2 * index + 1]))
    audit.equal(values, [Fraction(x) for x in ("0.25", "0.43", "0.19", "0.23", "0.66", "0.53", "0.57", "0.73", "0.82", "1.1")], "Fig4 values")
    audit.equal(sigmas, [Fraction(x) for x in ("0.18", "0.18", "0.18", "0.19", "0.2", "0.2", "0.2", "0.23", "0.23", "0.27")], "Fig4 uncertainties")
    audit.require(all(sigma > 0 for sigma in sigmas), "all marginal uncertainty widths positive")

    custody = json.loads((LANE / "SOURCE_CUSTODY.json").read_text(encoding="utf-8"))
    audit.equal(custody["repository"]["doi"], "10.5281/zenodo.8184043", "Zenodo DOI")
    audit.equal(custody["official_article"]["doi"], "10.1038/s41467-023-40629-8", "article DOI")
    audit.equal(custody["official_article"]["title"], "A lab-based test of the gravitational redshift with a miniature clock network", "article title")
    semantics = " ".join(custody["official_article"]["method_semantics"]).lower()
    for phrase in ("ten simultaneous pairwise comparisons", "same 14-run raw data", "assuming general relativity", "independently measured local g", "no ten-edge pair covariance matrix"):
        audit.require(phrase in semantics, f"source-custody semantic: {phrase}")

    b_matrix = incidence()
    audit.equal(exact_rank(b_matrix), 4, "K5 incidence rank")
    audit.equal(10 - exact_rank(b_matrix), 6, "K5 cycle dimension")
    for row, (start, end) in zip(b_matrix, EDGES):
        audit.equal(row[start - 1], 1, f"incidence + sign {start}->{end}")
        audit.equal(row[end - 1], -1, f"incidence - sign {start}->{end}")
        audit.equal(sum(row), 0, f"incidence row sum {start}->{end}")

    cycles = all_simple_cycles()
    audit.equal(len(cycles), 37, "all simple K5 cycles")
    audit.equal(sum(len(cycle) == 3 for cycle in cycles), 10, "triangle count")
    audit.equal(sum(len(cycle) == 4 for cycle in cycles), 15, "quadrilateral count")
    audit.equal(sum(len(cycle) == 5 for cycle in cycles), 12, "pentagon count")

    payload = json.loads((LANE / "RESULT.json").read_text(encoding="utf-8"))
    audit.equal(payload["schema"], "WAC_CLOCK_K5_COMMON_POTENTIAL_CYCLE_CLOSURE_V001", "result schema")
    audit.equal(payload["k5_exact_algebra"]["incidence_10x5"], b_matrix, "stored incidence")
    for key, expected in (("incidence_rank", 4), ("cut_space_dimension", 4), ("cycle_space_dimension", 6), ("cut_projector_rank", 4), ("cycle_projector_rank", 6), ("fundamental_cycle_rank", 6)):
        audit.equal(payload["k5_exact_algebra"][key], expected, f"stored {key}")

    stored_cycles = payload["simple_cycles"]["all_unique_up_to_rotation_and_reversal"]
    audit.equal(len(stored_cycles), 37, "stored cycle ledger length")
    stored = {tuple(record["nodes"]): record for record in stored_cycles}
    audit.equal(set(stored), set(cycles), "stored cycle identity set")

    cycle_results: dict[tuple[int, ...], tuple[list[int], Fraction, Fraction]] = {}
    maximum_residual = Fraction(0)
    maximum_ratio = Fraction(0)
    maximizers: list[tuple[int, ...]] = []
    for index, cycle in enumerate(cycles):
        row = cycle_row(cycle)
        reverse_row = cycle_row(tuple(reversed(cycle)))
        audit.equal(reverse_row, [-value for value in row], f"cycle orientation reversal {index}")
        audit.equal(matmul([row], b_matrix), [[Fraction(0)] * 5], f"cycle annihilates incidence {index}")
        residual = dot(row, values)
        envelope = sum(abs(weight) * sigma for weight, sigma in zip(row, sigmas))
        ratio = abs(residual) / envelope
        cycle_results[cycle] = row, residual, envelope
        maximum_residual = max(maximum_residual, abs(residual))
        if ratio > maximum_ratio:
            maximum_ratio = ratio
            maximizers = [cycle]
        elif ratio == maximum_ratio:
            maximizers.append(cycle)

        record = stored[cycle]
        audit.equal(record["cycle"], "-".join(str(node) for node in cycle + cycle[:1]), f"stored cycle name {index}")
        audit.equal(record["edge_coefficients_in_source_order"], row, f"stored cycle row {index}")
        fraction_record(audit, record["residual_cm"], residual, f"stored cycle residual {index}")
        fraction_record(audit, record["marginal_1sigma_l1_envelope_cm"], envelope, f"stored cycle envelope {index}")
        audit.close(record["absolute_residual_over_l1_envelope"], float(ratio), f"stored cycle ratio {index}")
        audit.equal(record["individually_closeable_inside_marginal_1sigma_box"], ratio <= 1, f"stored cycle interval flag {index}")

    audit.equal(maximum_residual, Fraction(29, 100), "maximum absolute cycle residual")
    audit.equal(maximum_ratio, Fraction(27, 82), "maximum simple-cycle ratio")
    audit.equal(maximizers, [(1, 3, 5, 4)], "unique dual maximizing cycle")
    audit.close(payload["simple_cycles"]["maximum_absolute_residual_cm"], 0.29, "stored maximum residual")
    audit.close(payload["simple_cycles"]["maximum_absolute_residual_over_marginal_l1_envelope"], float(Fraction(27, 82)), "stored maximum ratio")

    # Rounding-model preflight: all three source strings carry two displayed
    # decimal places, so conventional nearest rounding permits 0.005 cm each.
    witness_cycle = (2, 3, 4)
    _, witness_residual, _ = cycle_results[witness_cycle]
    half_ulp = [Fraction(1, 200) for _ in (fig4[2][2], fig4[2][4], fig4[2][10])]
    audit.equal([fig4[2][2], fig4[2][4], fig4[2][10]], ["0.43", "0.19", "0.53"], "rounding witness source strings")
    audit.equal(witness_residual, Fraction(9, 100), "rounding witness residual")
    audit.equal(sum(half_ulp), Fraction(3, 200), "rounding nearest half-ULP bound")
    audit.require(abs(witness_residual) > sum(half_ulp), "rounding witness excludes nearest-rounded closure")
    rounding = payload["preflight"]["rounding_exclusion_witness"]
    audit.equal(rounding["rounding_model"], "independent conventional nearest rounding to 0.01 cm", "stored rounding model")
    audit.equal(rounding["conclusion"], "NOT_A_COMMON_NODE_VECTOR_UNDER_INDEPENDENT_NEAREST_0.01_CM_ROUNDING", "stored rounding conclusion")

    # Independent rational primal/dual proof of rho*.  The maximizing cycle
    # gives the lower bound because its row annihilates every Bz.  Exact
    # difference constraints give a feasible upper-certificate vector.
    rho = Fraction(27, 82)
    dual_cycle = (1, 3, 5, 4)
    dual_row, dual_residual, dual_envelope = cycle_results[dual_cycle]
    audit.equal(dual_residual, Fraction(27, 100), "dual witness residual")
    audit.equal(dual_envelope, Fraction(41, 50), "dual witness envelope")
    audit.equal(abs(dual_residual) / dual_envelope, rho, "dual lower certificate")
    audit.equal(matmul([dual_row], b_matrix), [[Fraction(0)] * 5], "dual row annihilates cut space")

    nodes = exact_feasible_nodes_by_difference_constraints(values, sigmas, rho)
    audit.equal(nodes, [Fraction(4503, 4100), Fraction(1457, 1640), Fraction(2067, 4100), Fraction(2399, 8200), Fraction(0)], "independent feasible node vector")
    residuals = [value - dot(row, nodes) for value, row in zip(values, b_matrix)]
    standardized = [abs(residual) / sigma for residual, sigma in zip(residuals, sigmas)]
    audit.require(all(value <= rho for value in standardized), "independent primal vector feasible at rho")
    audit.equal(max(standardized), rho, "independent primal upper certificate")
    audit.require(abs(dual_residual) > (rho - Fraction(1, 10000)) * dual_envelope, "strict lower-rho dual sentinel")

    compatibility = payload["covariance_honest_compatibility"]
    audit.equal(compatibility["rho_star_exact"], "27/82", "stored rho exact")
    audit.close(compatibility["rho_star_decimal"], float(rho), "stored rho decimal")
    audit.equal(compatibility["compatible_with_all_reported_marginal_1sigma_intervals"], True, "stored marginal-box feasibility")
    payload_nodes = [Fraction(record["exact"]) for record in compatibility["one_exact_feasible_node_vector_cm_node5_zero"]]
    payload_residuals = [value - dot(row, payload_nodes) for value, row in zip(values, b_matrix)]
    audit.require(all(abs(residual) <= rho * sigma for residual, sigma in zip(payload_residuals, sigmas)), "stored primal certificate feasible")
    audit.equal(max(abs(residual) / sigma for residual, sigma in zip(payload_residuals, sigmas)), rho, "stored primal certificate reaches rho")

    # Exact projector checks are independent of the stored decimal diagnostics.
    b_transpose = [list(column) for column in zip(*b_matrix)]
    p_cut = [[value / 5 for value in row] for row in matmul(b_matrix, b_transpose)]
    identity = [[Fraction(int(row == column)) for column in range(10)] for row in range(10)]
    p_cycle = [[identity[row][column] - p_cut[row][column] for column in range(10)] for row in range(10)]
    audit.equal(exact_rank(p_cut), 4, "cut projector rank")
    audit.equal(exact_rank(p_cycle), 6, "cycle projector rank")
    audit.equal(matmul(p_cut, p_cut), p_cut, "cut projector idempotence")
    audit.equal(matmul(p_cycle, p_cycle), p_cycle, "cycle projector idempotence")
    audit.equal(matmul(p_cut, p_cycle), [[Fraction(0)] * 10 for _ in range(10)], "cut-cycle orthogonality")

    theorem = (LANE / "THEOREM.md").read_text(encoding="utf-8")
    result_md = (LANE / "RESULT.md").read_text(encoding="utf-8")
    self_audit = (LANE / "SELF_AUDIT.md").read_text(encoding="utf-8")
    hostile_md = (LANE / "INDEPENDENT_HOSTILE_AUDIT.md").read_text(encoding="utf-8")
    corpus = " ".join((theorem + result_md + self_audit + hostile_md).replace("**", "").lower().split())
    for phrase in (
        "conventional nearest rounding",
        "not statistically independent",
        "no joint covariance",
        "not a joint 68% coverage statement",
        "not an independent confirmation of gr",
        "necessary, never sufficient",
        "height vector is proportional",
        "record lineage",
        "gravity emergence",
        "no canonical model",
    ):
        audit.require(phrase in corpus, f"required scientific ceiling: {phrase}")

    manifest = parse_hashes(LANE / "MANIFEST.sha256")
    audit.equal(set(manifest), MANIFEST_MEMBERS, "manifest members")
    for member in sorted(MANIFEST_MEMBERS):
        audit.equal(digest(LANE / member), manifest[member], f"manifest hash {member}")
    seal = parse_hashes(LANE / "LANE_SEAL.sha256")
    audit.equal(set(seal), {"MANIFEST.sha256"}, "lane-seal member")
    audit.equal(seal["MANIFEST.sha256"], digest(LANE / "MANIFEST.sha256"), "lane-seal hash")

    tampered = bytearray((LANE / "SOURCE/Fig4.csv").read_bytes())
    tampered[-2] ^= 1
    audit.require(hashlib.sha256(tampered).hexdigest() != SOURCE_CONTRACT["SOURCE/Fig4.csv"][4], "source tamper sentinel")
    tampered_result = bytearray((LANE / "RESULT.json").read_bytes())
    tampered_result[len(tampered_result) // 2] ^= 1
    audit.require(hashlib.sha256(tampered_result).hexdigest() != manifest["RESULT.json"], "result tamper sentinel")

    print("CLOCK_K5_INDEPENDENT_HOSTILE_AUDIT: PASS")
    print(f"checks_passed: {audit.count}")
    print("official_pinned_csvs_reparsed: 4/4")
    print("source_semantics: PAIRWISE_REANALYZED__GR_AND_COMMON_G_ASSUMED__COVARIANCE_NOT_DEPOSITED")
    print("K5_cut_rank_cycle_dimension: 4 6")
    print("all_simple_cycles_recomputed: 37")
    print("rounding_witness_nearest_model_cm: 9/100 > 3/200")
    print("dual_cycle: 1-3-5-4-1")
    print("rho_star_exact_primal_dual: 27/82")
    print("marginal_box_intersection: PASS")
    print("ceiling: PROCESSED_HEIGHT_NODE_SCALAR_COMPATIBILITY_ONLY__NO_JOINT_STATISTICS_GR_GRAVITY_METRIC_LINEAGE_OR_EMERGENCE_PROMOTION")


if __name__ == "__main__":
    main()
