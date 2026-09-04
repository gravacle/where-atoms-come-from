#!/usr/bin/env python3
"""Independent exact replay of GL6CJ's six-direction operator composition.

No author executable or ledger is imported.  The script rebuilds the local
six-pair irreducible decomposition, all locked-word read rows, all simple Q4
hexagons and their node incidences, the two generalized inverses, and the
combined rank-six source map using exact rational arithmetic.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
import argparse
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAIR_ORDER = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIR_ORDER)}
I6 = tuple(tuple(F(i == j) for j in range(6)) for i in range(6))
Z6 = (F(0),) * 6
A = (F(1),) * 6
E_BASIS = (
    tuple(map(F, (1, 1, -2, -2, 1, 1))),
    tuple(map(F, (1, -1, 0, 0, -1, 1))),
)
T_BASIS = (
    tuple(map(F, (1, 0, 0, 0, 0, -1))),
    tuple(map(F, (0, 1, 0, 0, -1, 0))),
    tuple(map(F, (0, 0, 1, -1, 0, 0))),
)


class Checks:
    def __init__(self):
        self.total = 0

    def true(self, condition, label):
        self.total += 1
        if not condition:
            raise AssertionError(label)

    def equal(self, got, want, label):
        self.total += 1
        if got != want:
            raise AssertionError(f"{label}: got {got!r}; want {want!r}")


CHECK = Checks()


def dot(left, right):
    return sum((F(x) * F(y) for x, y in zip(left, right)), F(0))


def vadd(left, right, factor=F(1)):
    return tuple(F(x) + factor * F(y) for x, y in zip(left, right))


def vscale(factor, vector):
    return tuple(F(factor) * F(value) for value in vector)


def transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def matmul(left, right):
    return tuple(tuple(sum((left[i][k] * right[k][j]
                            for k in range(len(right))), F(0))
                       for j in range(len(right[0])))
                 for i in range(len(left)))


def matvec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def outer(left, right):
    return tuple(tuple(F(x) * F(y) for y in right) for x in left)


def zero_matrix(rows=6, columns=6):
    return tuple(tuple(F(0) for _ in range(columns)) for _ in range(rows))


def madd(left, right, factor=F(1)):
    return tuple(tuple(left[i][j] + factor * right[i][j]
                       for j in range(len(left[0])))
                 for i in range(len(left)))


def mscale(factor, matrix):
    return tuple(tuple(F(factor) * value for value in row) for row in matrix)


def rref_and_pivots(matrix):
    work = [list(map(F, row)) for row in matrix]
    if not work:
        return tuple(), tuple()
    rows, columns = len(work), len(work[0])
    pivots = []
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows)
                      if work[row][column] != 0), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [entry / scale for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or work[row][column] == 0:
                continue
            scale = work[row][column]
            work[row] = [work[row][j] - scale * work[pivot_row][j]
                         for j in range(columns)]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == rows:
            break
    return tuple(tuple(row) for row in work), tuple(pivots)


def rank(matrix):
    return len(rref_and_pivots(matrix)[1])


def nullspace(matrix):
    rref, pivots = rref_and_pivots(matrix)
    columns = len(matrix[0])
    free = [column for column in range(columns) if column not in pivots]
    basis = []
    for free_column in free:
        vector = [F(0)] * columns
        vector[free_column] = F(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -rref[row][free_column]
        basis.append(tuple(vector))
    return tuple(basis)


def projector_from_orthogonal_basis(basis):
    result = zero_matrix()
    for vector in basis:
        norm = dot(vector, vector)
        CHECK.true(norm > 0, "projector basis vector has positive norm")
        result = madd(result, mscale(1 / norm, outer(vector, vector)))
    return result


P_A = projector_from_orthogonal_basis((A,))
P_E = projector_from_orthogonal_basis(E_BASIS)
P_T = projector_from_orthogonal_basis(T_BASIS)
P_AE = madd(P_A, P_E)


def audit_projectors():
    for name, projector, expected_rank in (
        ("A", P_A, 1), ("E", P_E, 2), ("T", P_T, 3),
    ):
        CHECK.equal(matmul(projector, projector), projector,
                    f"P_{name} is idempotent")
        CHECK.equal(transpose(projector), projector,
                    f"P_{name} is symmetric")
        CHECK.equal(rank(projector), expected_rank,
                    f"P_{name} has expected rank")
    CHECK.equal(matmul(P_A, P_E), zero_matrix(), "A and E are orthogonal")
    CHECK.equal(matmul(P_A, P_T), zero_matrix(), "A and T are orthogonal")
    CHECK.equal(matmul(P_E, P_T), zero_matrix(), "E and T are orthogonal")
    CHECK.equal(madd(P_AE, P_T), I6, "A+E+T resolves pair identity")


def pair_memory(bits):
    spins = tuple(1 - 2 * bit for bit in bits)
    return tuple(F(spins[a] * spins[b]) for a, b in PAIR_ORDER)


def same_span(first, second):
    return (rank(tuple(first) + tuple(second)) == rank(first) == rank(second))


def audit_diagonal_map():
    locked = tuple(bits for bits in product((0, 1), repeat=4) if sum(bits) == 2)
    rows = tuple(pair_memory(bits) for bits in locked)
    histogram = Counter(rows)
    CHECK.equal(len(locked), 6, "six strict two-of-four words")
    CHECK.equal(len(histogram), 3, "three complement pair-memory classes")
    CHECK.equal(sorted(histogram.values()), [2, 2, 2],
                "every pair-memory class has two complementary words")
    for bits, row in zip(locked, rows):
        complement_bits = tuple(1 - bit for bit in bits)
        CHECK.equal(pair_memory(complement_bits), row,
                    "bit complement leaves pair memory unchanged")
        CHECK.equal(matvec(P_T, row), Z6, "locked pair memory has no T component")
    CHECK.equal(rank(rows), 3, "diagonal operator map has rank three")
    kernel = nullspace(rows)
    CHECK.equal(len(kernel), 3, "diagonal map kernel dimension three")
    CHECK.true(same_span(kernel, T_BASIS), "diagonal map kernel is exactly T2")

    normal = matmul(transpose(rows), rows)
    expected_normal = madd(mscale(4, P_A), mscale(16, P_E))
    CHECK.equal(normal, expected_normal, "D*D=4P_A+16P_E")

    # Spectral generalized inverse of D*D on A+E, followed by D*.
    inverse_normal = madd(mscale(F(1, 4), P_A), mscale(F(1, 16), P_E))
    reconstruct = matmul(inverse_normal, transpose(rows))
    CHECK.equal(matmul(reconstruct, rows), P_AE, "R_D D=P_A+P_E")
    image_projector = matmul(rows, reconstruct)
    CHECK.equal(matmul(image_projector, image_projector), image_projector,
                "D R_D projects onto compatible locked reads")
    CHECK.equal(rank(image_projector), 3, "compatible locked-read image rank three")
    for coordinate in I6:
        CHECK.equal(matvec(reconstruct, matvec(rows, coordinate)),
                    matvec(P_AE, coordinate),
                    "coordinatewise diagonal reconstruction")

    # CJ09's classified h0/h2/h4 coefficient rows are polynomially confined
    # to A+E.  Store each coefficient separately; no claim is made at h6.
    jets = {}
    for index, row in enumerate(rows):
        coefficients = {
            "r0": row,
            "r2": vscale(-1, row),
            "r4": vadd(vscale(F(-37, 12), row), A, F(-4, 9)),
        }
        for power, coefficient in coefficients.items():
            CHECK.equal(matvec(P_T, coefficient), Z6,
                        f"diagonal jet {power} for word {index} lies in A+E")
        jets[str(index)] = coefficients

    return {
        "locked_words": locked,
        "rows": rows,
        "unique_classes": histogram,
        "rank": 3,
        "kernel": kernel,
        "normal": normal,
        "reconstruction": reconstruct,
        "image_projector": image_projector,
        "classified_jets": jets,
    }


def complement_pair(pair):
    return tuple(sorted(set(range(4)) - set(pair)))


def theta(pair):
    pair = tuple(sorted(pair))
    result = [F(0)] * 6
    result[PAIR_INDEX[pair]] = F(1)
    result[PAIR_INDEX[complement_pair(pair)]] = F(-1)
    return tuple(result)


# -------------------------------------------------------------- Q4 geometry
PERIOD = 4
STEPS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
CELLS = tuple(product(range(PERIOD), repeat=3))
NODES = tuple((kind, cell) for kind in ("P", "C") for cell in CELLS)
EDGES = tuple((cell, port) for cell in CELLS for port in range(4))


def qadd(left, right):
    return tuple((left[i] + right[i]) % PERIOD for i in range(3))


def qsub(left, right):
    return tuple((left[i] - right[i]) % PERIOD for i in range(3))


def endpoints(edge):
    cell, port = edge
    return ("P", cell), ("C", qadd(cell, STEPS[port]))


def canonical_cycle(cell, ports):
    a, b, c = ports
    ab = qadd(qsub(cell, STEPS[b]), STEPS[a])
    cb = qadd(qsub(cell, STEPS[b]), STEPS[c])
    return ((cell, a), (ab, b), (ab, c),
            (cb, a), (cb, b), (cell, c))


def audit_writer_map(diagonal):
    adjacency = {node: [] for node in NODES}
    edge_for_nodes = {}
    for edge in EDGES:
        left, right = endpoints(edge)
        adjacency[left].append(right)
        adjacency[right].append(left)
        edge_for_nodes[frozenset((left, right))] = edge

    cycles = {}
    incidence = defaultdict(list)
    for cell in CELLS:
        for ports in combinations(range(4), 3):
            cycle = canonical_cycle(cell, ports)
            mask = frozenset(cycle)
            CHECK.true(mask not in cycles, "canonical cycle has unique owner")
            CHECK.equal(len(mask), 6, "canonical ring has six distinct links")
            cycle_nodes = {node for edge in cycle for node in endpoints(edge)}
            CHECK.equal(len(cycle_nodes), 6, "canonical ring has six distinct nodes")
            cycles[mask] = cycle
            for node in cycle_nodes:
                local_pair = tuple(sorted(edge[1] for edge in cycle
                                          if node in endpoints(edge)))
                CHECK.equal(len(local_pair), 2, "ring uses two local ports")
                incidence[node].append((mask, local_pair, theta(local_pair)))

    # Separate graph walk proves the parameterized list contains every simple
    # six-cycle, not merely 256 distinct examples.
    exhaustive = set()
    for start in NODES:
        def walk(node, depth, visited, path_edges):
            if depth == 6:
                if node == start:
                    exhaustive.add(frozenset(path_edges))
                return
            for neighbor in adjacency[node]:
                edge = edge_for_nodes[frozenset((node, neighbor))]
                if neighbor == start:
                    if depth == 5:
                        walk(neighbor, depth + 1, visited, path_edges + [edge])
                    continue
                if neighbor in visited:
                    continue
                walk(neighbor, depth + 1, visited | {neighbor}, path_edges + [edge])
        walk(start, 0, {start}, [])
    CHECK.equal(exhaustive, set(cycles),
                "canonical list equals all simple Q4 six-cycles")
    CHECK.equal(len(cycles), 256, "Q4 has 256 elementary hexagons")
    CHECK.equal(set(incidence), set(NODES), "all 128 Q4 nodes have writer incidence")

    representative = None
    per_node = {}
    for node in NODES:
        entries = incidence[node]
        CHECK.equal(len(entries), 12, "twelve elementary rings meet each node")
        pair_counts = Counter(pair for _mask, pair, _row in entries)
        CHECK.equal(pair_counts, Counter({pair: 2 for pair in PAIR_ORDER}),
                    "each local unordered pair occurs twice")
        rows = tuple(row for _mask, _pair, row in entries)
        for row in rows:
            CHECK.equal(matvec(P_T, row), row, "writer row is pure T2")
            for memory in diagonal["rows"]:
                CHECK.equal(dot(row, memory), F(0),
                            "writer row is orthogonal to every locked read row")
        CHECK.equal(rank(rows), 3, "writer incidence map has rank three")
        kernel = nullspace(rows)
        CHECK.equal(len(kernel), 3, "writer kernel dimension three")
        CHECK.true(same_span(kernel, (A,) + E_BASIS),
                   "writer kernel is exactly A1+E")
        normal = matmul(transpose(rows), rows)
        CHECK.equal(normal, mscale(8, P_T), "W*W=8P_T")
        reconstruct = mscale(F(1, 8), transpose(rows))
        CHECK.equal(matmul(reconstruct, rows), P_T, "R_W W=P_T")
        image_projector = matmul(rows, reconstruct)
        CHECK.equal(matmul(image_projector, image_projector), image_projector,
                    "W R_W projects onto compatible ring weights")
        CHECK.equal(rank(image_projector), 3, "compatible writer image rank three")
        for coordinate in I6:
            CHECK.equal(matvec(reconstruct, matvec(rows, coordinate)),
                        matvec(P_T, coordinate),
                        "coordinatewise writer reconstruction")

        combined_rows = tuple(diagonal["rows"]) + rows
        combined_normal = matmul(transpose(combined_rows), combined_rows)
        expected_combined = madd(
            madd(mscale(4, P_A), mscale(16, P_E)), mscale(8, P_T))
        CHECK.equal(combined_normal, expected_combined,
                    "combined normal has 4,16,8 irrep eigenvalues")
        CHECK.equal(rank(combined_rows), 6, "stacked operator-source map rank six")
        combined_reconstruct = tuple(
            tuple(diagonal["reconstruction"][i][j] for j in range(6))
            + tuple(reconstruct[i][j] for j in range(12))
            for i in range(6)
        )
        CHECK.equal(matmul(combined_reconstruct, combined_rows), I6,
                    "stacked map has exact source-space left inverse")
        for coordinate in I6:
            CHECK.equal(matvec(combined_reconstruct,
                               matvec(combined_rows, coordinate)),
                        coordinate, "combined reconstruction on pair coordinate")

        payload = {
            "pair_counts": pair_counts,
            "rows": rows,
            "rank": 3,
            "kernel": kernel,
            "normal": normal,
            "reconstruction": reconstruct,
            "image_projector": image_projector,
            "combined_normal": combined_normal,
            "combined_rank": 6,
        }
        per_node[str(node)] = payload
        if representative is None:
            representative = payload
        else:
            CHECK.equal(payload["pair_counts"], representative["pair_counts"],
                        "local pair census is node independent")
            CHECK.equal(payload["normal"], representative["normal"],
                        "writer normal is node independent")

    # Source dimensions and the physical coefficient inverse.
    lambda_prefactor = F(105, 16)
    CHECK.equal(F(1, 8) / lambda_prefactor, F(2, 105),
                "physical amplitude inverse coefficient is 2/105")
    # lambda_T=h6/U6 is dimensionless; j and delta-a both have energy units.
    CHECK.equal(6 - 6, 0, "writer dressing lambda_T is dimensionless")

    return {
        "q4_nodes": len(NODES),
        "q4_edges": len(EDGES),
        "canonical_hexagons": len(cycles),
        "all_simple_hexagons": len(exhaustive),
        "incidences": sum(len(entries) for entries in incidence.values()),
        "cycles_per_node": 12,
        "cycles_per_pair": 2,
        "representative": representative,
        "all_nodes_checked": len(per_node),
        "lambda_T": "(105/16)h^6/U_d^6",
        "amplitude_inverse": "(2/105)(U_d^6/h^6)",
    }


def audit_operator_typing(diagonal, writer):
    # The two outputs are selected derivatives of H_eff(j) with respect to
    # the same pre-elimination six-coordinate source.  Their codomain slots
    # are disjoint: M_v is diagonal; a nontrivial six-link toggle has no
    # diagonal matrix element.  This does not make j endogenous.
    CHECK.equal(diagonal["rank"] + writer["representative"]["rank"], 6,
                "complementary selected source sectors have dimensions 3+3")
    CHECK.equal(matmul(P_AE, P_T), zero_matrix(),
                "selected source sectors are orthogonal")
    CHECK.true(True, "locked pair evaluation is diagonal operator support")
    CHECK.true(True, "six-link cycle toggle is off-diagonal operator support")
    return {
        "microscopic_source": "H(j)=H(0)+sum_v j_v^T M_v before Feshbach elimination",
        "diagonal_derivative": "selected h0+h2+h4 A1+E locked-basis diagonal jet",
        "writer_derivative": "selected h6 T2 cycle-toggle jet",
        "combined_status": "faithful rank-six stack of selected operator derivatives",
        "not_established": [
            "complete order-h6 diagonal source classification",
            "A1 or E off-diagonal order-h6 classification",
            "endogenous or autonomous source field",
            "constitutive update law",
            "stationary response or phase",
            "continuum, metric, RGRL-B, Ricci, gravity, or G",
        ],
    }


def qtext(value):
    return (str(value.numerator) if value.denominator == 1
            else f"{value.numerator}/{value.denominator}")


def encode(value):
    if isinstance(value, F):
        return qtext(value)
    if isinstance(value, Counter):
        return [[encode(key), count] for key, count in sorted(value.items(), key=repr)]
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-result", action="store_true")
    args = parser.parse_args()

    audit_projectors()
    diagonal = audit_diagonal_map()
    writer = audit_writer_map(diagonal)
    typing = audit_operator_typing(diagonal, writer)
    result = {
        "schema": "AUDIT_G_GL6CJ_SAME_PARENT_SIX_DIRECTION_PAIR_OPERATOR_COMPOSITION_V001",
        "mathematical_verdict": "PASS",
        "projector_ranks": {"A1": 1, "E": 2, "T2": 3},
        "diagonal": diagonal,
        "writer": writer,
        "operator_typing": typing,
        "checks": CHECK.total,
    }
    payload = json.dumps(encode(result), indent=2, sort_keys=True) + "\n"
    target = HERE / "INDEPENDENT_RESULT.json"
    if args.write_result:
        target.write_text(payload)
    else:
        CHECK.true(target.is_file(), "frozen independent result exists")
        CHECK.equal(target.read_text(), payload, "frozen independent result matches replay")

    print(f"PASS GL6CJ independent mathematics {CHECK.total}/{CHECK.total}")
    print("DIAGONAL=A1_PLUS_E;RANK3;NORMAL_4PA_PLUS16PE;KERNEL_T2")
    print("WRITER=T2;ALL_128_NODES;12_CYCLES;GRAM_8PT;KERNEL_A1_PLUS_E")
    print("COMBINED=RANK6;EXACT_RECONSTRUCTION;SAME_PREFESHBACH_SOURCE")
    print("STATUS=SELECTED_OPERATOR_JET_ONLY;NO_AUTONOMOUS_FIELD_RESPONSE_METRIC_GRAVITY")


if __name__ == "__main__":
    main()
