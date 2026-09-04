#!/usr/bin/env python3
"""Exact GL6CJ same-parent six-direction pair-operator composition.

The replay constructs the locked diagonal pair read and the GL6CH order-h6
off-diagonal tensor writer from one pre-Feshbach six-pair source.  It proves
their complementary ranks and right inverses at every node of the declared
Q4 incidence quotient.  All arithmetic is exact and only the standard
library is used.
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
A = tuple(map(F, (1, 1, 1, 1, 1, 1)))
T_BASIS = (
    tuple(map(F, (1, 0, 0, 0, 0, -1))),
    tuple(map(F, (0, 1, 0, 0, -1, 0))),
    tuple(map(F, (0, 0, 1, -1, 0, 0))),
)


class Checks:
    def __init__(self):
        self.total = 0

    def equal(self, got, want, label):
        self.total += 1
        if got != want:
            raise AssertionError(f"{label}: got {got!r}, want {want!r}")

    def true(self, condition, label):
        self.total += 1
        if not condition:
            raise AssertionError(label)


CHECK = Checks()


def dot(left, right):
    return sum((F(x) * F(y) for x, y in zip(left, right)), F(0))


def vadd(left, right, factor=F(1)):
    return tuple(F(x) + factor * F(y) for x, y in zip(left, right))


def vscale(factor, vector):
    return tuple(F(factor) * F(value) for value in vector)


def zero_matrix(rows, columns):
    return tuple(tuple(F(0) for _ in range(columns)) for _ in range(rows))


def identity(size):
    return tuple(tuple(F(i == j) for j in range(size)) for i in range(size))


def transpose(matrix):
    return tuple(tuple(matrix[i][j] for i in range(len(matrix)))
                 for j in range(len(matrix[0])))


def matmul(left, right):
    return tuple(tuple(sum((left[i][k] * right[k][j]
                            for k in range(len(right))), F(0))
                       for j in range(len(right[0])))
                 for i in range(len(left)))


def matvec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def outer(left, right):
    return tuple(tuple(F(x) * F(y) for y in right) for x in left)


def madd(left, right, factor=F(1)):
    return tuple(tuple(left[i][j] + factor * right[i][j]
                       for j in range(len(left[0])))
                 for i in range(len(left)))


def mscale(factor, matrix):
    return tuple(tuple(F(factor) * value for value in row) for row in matrix)


def rank(matrix):
    work = [list(map(F, row)) for row in matrix]
    if not work:
        return 0
    rows, columns = len(work), len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][column]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            value = work[row][column]
            work[row] = [work[row][j] - value * work[pivot_row][j]
                         for j in range(columns)]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


P_A = mscale(F(1, 6), outer(A, A))
P_T = zero_matrix(6, 6)
for direction in T_BASIS:
    P_T = madd(P_T, mscale(F(1, 2), outer(direction, direction)))
P_E = madd(madd(identity(6), P_A, F(-1)), P_T, F(-1))
P_AE = madd(P_A, P_E)

CHECK.equal(matmul(P_A, P_A), P_A, "P_A is a projector")
CHECK.equal(matmul(P_E, P_E), P_E, "P_E is a projector")
CHECK.equal(matmul(P_T, P_T), P_T, "P_T is a projector")
CHECK.equal(matmul(P_A, P_E), zero_matrix(6, 6), "A and E are orthogonal")
CHECK.equal(matmul(P_A, P_T), zero_matrix(6, 6), "A and T are orthogonal")
CHECK.equal(matmul(P_E, P_T), zero_matrix(6, 6), "E and T are orthogonal")
CHECK.equal(rank(P_A), 1, "A rank one")
CHECK.equal(rank(P_E), 2, "E rank two")
CHECK.equal(rank(P_T), 3, "T rank three")
CHECK.equal(madd(madd(P_A, P_E), P_T), identity(6),
            "A+E+T resolves six-pair identity")


def pair_memory(bits):
    z = tuple(F(1 - 2 * bit) for bit in bits)
    return tuple(z[a] * z[b] for a, b in PAIR_ORDER)


def diagonal_read_census():
    locked_words = tuple(bits for bits in product((0, 1), repeat=4)
                         if sum(bits) == 2)
    rows = tuple(pair_memory(bits) for bits in locked_words)
    histogram = Counter(rows)
    CHECK.equal(len(locked_words), 6, "six degree-two locked local words")
    CHECK.equal(len(histogram), 3, "three complement-identified pair words")
    CHECK.equal(sorted(histogram.values()), [2, 2, 2],
                "each pair word has two complement-related locked words")
    CHECK.true(all(matvec(P_T, row) == (F(0),) * 6 for row in rows),
               "every locked pair word lies in A+E")
    CHECK.equal(rank(rows), 3, "locked diagonal read has rank three")

    normal = zero_matrix(6, 6)
    for row in rows:
        normal = madd(normal, outer(row, row))
    expected_normal = madd(mscale(F(4), P_A), mscale(F(16), P_E))
    CHECK.equal(normal, expected_normal,
                "sum over six locked words is 4 P_A + 16 P_E")

    # R_D maps compatible six locked-word scores back to the A+E part of j.
    read_inverse_columns = []
    for row in rows:
        column = vadd(vscale(F(1, 4), matvec(P_A, row)),
                      vscale(F(1, 16), matvec(P_E, row)))
        read_inverse_columns.append(column)
    read_inverse = transpose(tuple(read_inverse_columns))
    read_matrix = rows
    CHECK.equal(matmul(read_inverse, read_matrix), P_AE,
                "diagonal-read right inverse reconstructs P_A+P_E")
    target_projector = matmul(read_matrix, read_inverse)
    CHECK.equal(rank(target_projector), 3,
                "compatible diagonal-read target projector has rank three")
    CHECK.equal(matmul(target_projector, target_projector), target_projector,
                "diagonal-read target map is a projector")

    for coordinate in identity(6):
        scores = matvec(read_matrix, coordinate)
        recovered = matvec(read_inverse, scores)
        CHECK.equal(recovered, matvec(P_AE, coordinate),
                    "explicit diagonal reconstruction of a pair coordinate")

    return {
        "locked_words": locked_words,
        "read_rows": rows,
        "unique_pair_words": histogram,
        "rank": 3,
        "normal": normal,
        "right_inverse": read_inverse,
        "target_projector": target_projector,
    }


def complement_pair(pair):
    return tuple(sorted(set(range(4)) - set(pair)))


def theta(pair):
    pair = tuple(sorted(pair))
    answer = [F(0)] * 6
    answer[PAIR_INDEX[pair]] = F(1)
    answer[PAIR_INDEX[complement_pair(pair)]] = F(-1)
    return tuple(answer)


# -------------------------------------------------------------- Q4 incidence
PERIOD = 4
STEPS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
CELLS = tuple(product(range(PERIOD), repeat=3))
EDGES = tuple((cell, port) for cell in CELLS for port in range(4))
NODES = tuple((kind, cell) for kind in ("P", "C") for cell in CELLS)


def qadd(left, right):
    return tuple((left[index] + right[index]) % PERIOD for index in range(3))


def qsub(left, right):
    return tuple((left[index] - right[index]) % PERIOD for index in range(3))


def endpoints(edge):
    cell, port = edge
    return ("P", cell), ("C", qadd(cell, STEPS[port]))


def canonical_cycle(cell, ports):
    a, b, c = ports
    ab = qadd(qsub(cell, STEPS[b]), STEPS[a])
    cb = qadd(qsub(cell, STEPS[b]), STEPS[c])
    return ((cell, a), (ab, b), (ab, c),
            (cb, a), (cb, b), (cell, c))


def writer_census(read):
    cycles = []
    seen = set()
    by_node = defaultdict(list)
    edge_counts = Counter()
    for cell in CELLS:
        for ports in combinations(range(4), 3):
            cycle = canonical_cycle(cell, ports)
            mask = frozenset(cycle)
            CHECK.true(mask not in seen, "Q4 canonical cycle is unique")
            seen.add(mask)
            CHECK.equal(len(mask), 6, "Q4 canonical cycle has six links")
            local = defaultdict(list)
            for edge in cycle:
                edge_counts[edge] += 1
                for node in endpoints(edge):
                    local[node].append(edge[1])
            CHECK.equal(len(local), 6, "Q4 canonical cycle has six nodes")
            CHECK.true(all(len(ports_at_node) == 2
                           for ports_at_node in local.values()),
                       "two cycle ports meet at every cycle node")
            for node, ports_at_node in local.items():
                pair = tuple(sorted(ports_at_node))
                by_node[node].append((len(cycles), pair, theta(pair)))
            cycles.append((cycle, ports))

    CHECK.equal(len(cycles), 256, "Q4 has 256 elementary hexagons")
    CHECK.true(all(edge_counts[edge] == 6 for edge in EDGES),
               "each Q4 link lies in six hexagons")
    CHECK.equal(sum(len(rows) for rows in by_node.values()), 256 * 6,
                "Q4 node-cycle incidence count is 1536")
    CHECK.equal(set(by_node), set(NODES), "every Q4 node is reached")

    reference_writer = None
    local_payload = {}
    for node in NODES:
        entries = by_node[node]
        CHECK.equal(len(entries), 12, "exactly twelve hexagons meet each node")
        pair_counts = Counter(pair for _index, pair, _vector in entries)
        CHECK.equal(pair_counts, Counter({pair: 2 for pair in PAIR_ORDER}),
                    "each local port pair occurs on exactly two hexagons")
        rows = tuple(vector for _index, _pair, vector in entries)
        CHECK.true(all(matvec(P_T, row) == row for row in rows),
                   "every local cycle tensor lies in T")
        CHECK.true(all(dot(memory, vector) == 0
                       for memory in read["read_rows"] for vector in rows),
                   "locked diagonal words are orthogonal to writer tensors")
        normal = zero_matrix(6, 6)
        for row in rows:
            normal = madd(normal, outer(row, row))
        CHECK.equal(normal, mscale(F(8), P_T),
                    "local cycle Gram is 8 P_T")
        CHECK.equal(rank(rows), 3, "local tensor writer has rank three")

        writer_inverse = mscale(F(1, 8), transpose(rows))
        CHECK.equal(matmul(writer_inverse, rows), P_T,
                    "writer right inverse reconstructs P_T")
        target_projector = matmul(rows, writer_inverse)
        CHECK.equal(rank(target_projector), 3,
                    "writer target projector has rank three")
        CHECK.equal(matmul(target_projector, target_projector), target_projector,
                    "writer target map is a projector")
        for coordinate in identity(6):
            coefficients = matvec(rows, coordinate)
            recovered = matvec(writer_inverse, coefficients)
            CHECK.equal(recovered, matvec(P_T, coordinate),
                        "explicit writer reconstruction of a pair coordinate")

        combined = tuple(read["read_rows"]) + rows
        CHECK.equal(rank(combined), 6,
                    "diagonal plus off-diagonal source map has rank six")
        combined_normal = madd(read["normal"], normal)
        expected_combined = madd(
            madd(mscale(F(4), P_A), mscale(F(16), P_E)),
            mscale(F(8), P_T))
        CHECK.equal(combined_normal, expected_combined,
                    "combined normal is 4 P_A + 16 P_E + 8 P_T")
        for coordinate in identity(6):
            read_scores = matvec(read["read_rows"], coordinate)
            writer_scores = matvec(rows, coordinate)
            reconstructed = vadd(
                matvec(read["right_inverse"], read_scores),
                matvec(writer_inverse, writer_scores))
            CHECK.equal(reconstructed, coordinate,
                        "combined right inverse reconstructs full pair source")

        signature = (pair_counts, normal, writer_inverse, target_projector)
        if reference_writer is None:
            reference_writer = signature
            local_payload = {
                "representative_node": node,
                "rows": rows,
                "pair_counts": pair_counts,
                "normal": normal,
                "right_inverse": writer_inverse,
                "target_projector": target_projector,
                "combined_normal": combined_normal,
            }
        else:
            # Row order depends on cycle enumeration, so compare all
            # translation-invariant algebraic data except the literal rows.
            CHECK.equal(signature[0], reference_writer[0],
                        "node pair census is uniform")
            CHECK.equal(signature[1], reference_writer[1],
                        "node writer Gram is uniform")

    return {
        "q4_cells": len(CELLS),
        "q4_nodes": len(NODES),
        "q4_links": len(EDGES),
        "q4_hexagons": len(cycles),
        "total_node_cycle_incidences": sum(len(rows) for rows in by_node.values()),
        "cycles_per_node": 12,
        "cycles_per_local_pair": 2,
        "rank": 3,
        **local_payload,
    }


def qtext(value):
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


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
    parser.add_argument("--write-ledger", action="store_true")
    args = parser.parse_args()

    read = diagonal_read_census()
    writer = writer_census(read)
    ledger = {
        "lane": "GL6CJ",
        "single_source": "H(j)=H(0)+sum_v j_v.M_v before Feshbach elimination",
        "pair_order": PAIR_ORDER,
        "projectors": {"P_A": P_A, "P_E": P_E, "P_T": P_T},
        "diagonal_locked_read": read,
        "h6_tensor_writer": writer,
        "combined": {
            "rank": 6,
            "normal": madd(madd(mscale(F(4), P_A), mscale(F(16), P_E)),
                            mscale(F(8), P_T)),
            "reconstruction": "j=R_D d + (1/8) sum_{c incident v} w_c Theta_{v,c}",
            "dressed_writer_scale": "lambda_T=(105/16)h^6/U_d^6",
            "physical_amplitude_inverse": "j_T=(2/105)(U_d^6/h^6) sum_c delta_a_c Theta_{v,c}",
            "operator_support": "locked read diagonal; ring writer off-diagonal",
        },
        "scope": "finite-range perturbative operator/source-access theorem; no autonomous field, response, metric, Ricci, gravity, or G",
    }
    payload = json.dumps(encode(ledger), indent=2, sort_keys=True) + "\n"
    target = HERE / "EXACT_LEDGER.json"
    if args.write_ledger:
        target.write_text(payload)
    else:
        CHECK.true(target.is_file(), "frozen exact ledger exists")
        CHECK.equal(target.read_text(), payload, "frozen exact ledger matches replay")

    print(f"PASS__GL6CJ_SAME_PARENT_PAIR_COMPOSITION__{CHECK.total}/{CHECK.total}")
    print("DIAGONAL_LOCKED_READ=RANK3_A1_PLUS_E;NORMAL=4PA_PLUS16PE")
    print("H6_TENSOR_WRITER=RANK3_T2;LOCAL_GRAM=8PT;HEXAGONS_PER_NODE=12")
    print("COMBINED_SINGLE_SOURCE=RANK6;EXACT_COMPLEMENTARY_RIGHT_INVERSES")
    print("TYPED_SPLIT=OPERATOR_JET_CLOSED;AUTONOMOUS_CONSTITUTIVE_RESPONSE_OPEN")
    print("NO_STATIONARY_RESPONSE_METRIC_RGRLB_RICCI_GRAVITY_G")


if __name__ == "__main__":
    main()
