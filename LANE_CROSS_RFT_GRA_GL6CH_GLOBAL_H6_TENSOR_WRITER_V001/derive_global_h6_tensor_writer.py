#!/usr/bin/env python3
"""Exact global order-h6 tensor-writer derivation for GL6CH.

This replay is deliberately independent of the GL6AO/BX executables.  It
enumerates all 720 orders of an alternating hexagon, differentiates the
canonical Hermitian endpoint-symmetrized amplitude, projects the result into
the six-pair T2 sector, and reconstructs all 256 hexagons of the declared Q4
incidence quotient.  It also proves locally that the global h2 and h4 first
vertices have no T2 component, so lower-order folds cannot manufacture the
off-diagonal h6 result.

Only the Python standard library is used; all arithmetic is exact.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
import argparse
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAIR_ORDER = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIR_ORDER)}
A1 = (F(1),) * 6
EA = tuple(map(F, (1, 1, -2, -2, 1, 1)))
EB = tuple(map(F, (1, -1, 0, 0, -1, 1)))
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
            raise AssertionError(f"{label}: got {got!r}, want {want!r}")


CHECK = Checks()


def dot(left, right):
    return sum((F(x) * F(y) for x, y in zip(left, right)), F(0))


def vadd(left, right, factor=F(1)):
    return tuple(F(x) + factor * F(y) for x, y in zip(left, right))


def vscale(factor, vector):
    return tuple(F(factor) * F(value) for value in vector)


def outer(left, right):
    return tuple(tuple(F(x) * F(y) for y in right) for x in left)


def madd(left, right, factor=F(1)):
    return tuple(tuple(left[i][j] + factor * right[i][j]
                       for j in range(len(left[i])))
                 for i in range(len(left)))


def zero_matrix(rows, columns):
    return tuple(tuple(F(0) for _ in range(columns)) for _ in range(rows))


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


def project_t(vector):
    answer = (F(0),) * 6
    for direction in T_BASIS:
        CHECK.equal(dot(direction, direction), F(2),
                    "T basis vector has norm squared two")
        answer = vadd(answer, direction, dot(direction, vector) / 2)
    return answer


def complement_pair(pair):
    return tuple(sorted(set(range(4)) - set(pair)))


def theta(pair):
    pair = tuple(sorted(pair))
    complement = complement_pair(pair)
    answer = [F(0)] * 6
    answer[PAIR_INDEX[pair]] = F(1)
    answer[PAIR_INDEX[complement]] = F(-1)
    return tuple(answer)


def pair_products(z):
    return tuple(F(z[a]) * F(z[b]) for a, b in PAIR_ORDER)


def cycle_energy(selected, signs):
    """Defect energy of a partial alternating six-cycle toggle."""
    charges = [0] * 6
    for edge in selected:
        charges[edge] += signs[edge]
        charges[(edge + 1) % 6] += signs[edge]
    return sum(value * value for value in charges)


def local_words(vertex, signs, selected, pair, pair_orientation,
                complement_orientation):
    """Return initial, intermediate, final pair words at one cycle vertex."""
    previous_edge = (vertex - 1) % 6
    next_edge = vertex
    pair_ports = tuple(pair if pair_orientation == 0 else reversed(pair))
    complement = complement_pair(pair)
    complement_ports = tuple(complement if complement_orientation == 0
                             else reversed(complement))
    z0 = [None] * 4
    z0[pair_ports[0]] = signs[previous_edge]
    z0[pair_ports[1]] = signs[next_edge]
    z0[complement_ports[0]] = 1
    z0[complement_ports[1]] = -1
    CHECK.true(z0.count(1) == 2 and z0.count(-1) == 2,
               "local initial word is degree two")
    z = list(z0)
    if previous_edge in selected:
        z[pair_ports[0]] *= -1
    if next_edge in selected:
        z[pair_ports[1]] *= -1
    zf = list(z0)
    zf[pair_ports[0]] *= -1
    zf[pair_ports[1]] *= -1
    return pair_products(z0), pair_products(z), pair_products(zf), pair_ports


def direct_hexagon_census():
    """All direct Q-only paths and their first microscopic source jet."""
    energy_histograms = []
    base_amplitudes = []
    canonical_gradients = Counter()
    tempting_differences = Counter()
    context_count = 0

    for phase in (1, -1):
        signs = tuple(phase * (1 if edge % 2 == 0 else -1)
                      for edge in range(6))
        amplitude = F(0)
        profile_histogram = Counter()
        prefix_weights = Counter()
        for order in permutations(range(6)):
            selected = set()
            energies = []
            weight = F(1)
            for edge in order[:-1]:
                selected.add(edge)
                energy = cycle_energy(selected, signs)
                CHECK.true(energy > 0,
                           "proper nonempty cycle prefix remains outside lock")
                energies.append(energy)
                weight /= energy
            profile_histogram[tuple(energies)] += 1
            amplitude -= weight
            selected = set()
            for step, edge in enumerate(order[:-1]):
                selected.add(edge)
                prefix_weights[frozenset(selected)] += weight / energies[step]
        CHECK.equal(sum(profile_histogram.values()), 720,
                    "all 720 hexagon orders retained")
        CHECK.equal(amplitude, F(-63, 8),
                    "alternating-hexagon direct amplitude")
        energy_histograms.append(profile_histogram)
        base_amplitudes.append(amplitude)

        for vertex in range(6):
            for pair in PAIR_ORDER:
                direction = theta(pair)
                CHECK.equal(dot(direction, direction), F(2),
                            "Theta has norm squared two")
                CHECK.equal(project_t(direction), direction,
                            "Theta is pure T2")
                CHECK.equal(dot(A1, direction), F(0),
                            "Theta is trace free")
                CHECK.equal(dot(EA, direction), F(0),
                            "Theta is orthogonal to Ea")
                CHECK.equal(dot(EB, direction), F(0),
                            "Theta is orthogonal to Eb")
                for pair_orientation in (0, 1):
                    for complement_orientation in (0, 1):
                        context_count += 1
                        gradient = [F(0)] * 6
                        first_delta = None
                        second_delta = None
                        for selected, weighted_denominator in prefix_weights.items():
                            before, middle, final, pair_ports = local_words(
                                vertex, signs, selected, pair,
                                pair_orientation, complement_orientation)
                            midpoint = tuple((before[index] + final[index]) / 2
                                             for index in range(6))
                            score = tuple(middle[index] - midpoint[index]
                                          for index in range(6))
                            # Pure-T endpoints vanish, so the canonical
                            # midpoint score equals the ordinary defect score.
                            theta_score = dot(direction, score)
                            defective = ((vertex - 1) % 6 in selected) ^ (
                                vertex in selected)
                            CHECK.equal(theta_score,
                                        F(2) if defective else F(0),
                                        "Theta score is two exactly at a defect")
                            for component in range(6):
                                gradient[component] += (
                                    weighted_denominator * score[component])

                        # Single-edge changes define a tempting but incomplete
                        # full-vector shortcut.  Canonical endpoint
                        # symmetrization agrees with it only in A1 and T2.
                        empty = frozenset()
                        _, initial, _, pair_ports = local_words(
                            vertex, signs, empty, pair,
                            pair_orientation, complement_orientation)
                        _, after_first, _, _ = local_words(
                            vertex, signs, {((vertex - 1) % 6)}, pair,
                            pair_orientation, complement_orientation)
                        _, after_second, _, _ = local_words(
                            vertex, signs, {vertex}, pair,
                            pair_orientation, complement_orientation)
                        first_delta = vadd(after_first, initial, F(-1))
                        second_delta = vadd(after_second, initial, F(-1))
                        tempting = vscale(F(105, 32),
                                          vadd(first_delta, second_delta))
                        gradient = tuple(gradient)
                        expected_full = tuple(
                            F(105, 8) if index == PAIR_INDEX[tuple(sorted(pair))]
                            else F(0) for index in range(6))
                        CHECK.equal(gradient, expected_full,
                                    "canonical full gradient is (105/8)e_ab")
                        CHECK.equal(dot(A1, gradient), F(105, 8),
                                    "universal A1 derivative per vertex")
                        CHECK.equal(project_t(gradient),
                                    vscale(F(105, 16), direction),
                                    "canonical T2 gradient is (105/16)Theta")
                        CHECK.equal(dot(direction, gradient), F(105, 8),
                                    "j=Theta derivative is 105/8")
                        CHECK.equal(dot(A1, tempting), dot(A1, gradient),
                                    "tempting shortcut shares A1 projection")
                        CHECK.equal(project_t(tempting), project_t(gradient),
                                    "tempting shortcut shares T2 projection")
                        difference = vadd(tempting, gradient, F(-1))
                        CHECK.equal(project_t(difference), (F(0),) * 6,
                                    "shortcut error has no T2 component")
                        CHECK.equal(dot(A1, difference), F(0),
                                    "shortcut error has no A1 component")
                        canonical_gradients[gradient] += 1
                        tempting_differences[difference] += 1

    CHECK.equal(energy_histograms[0], energy_histograms[1],
                "energy profiles are independent of alternating orientation")
    CHECK.equal(base_amplitudes, [F(-63, 8), F(-63, 8)],
                "both directions have the same source-free amplitude")
    CHECK.equal(context_count, 2 * 6 * 6 * 2 * 2,
                "all local port/orientation contexts exhausted")
    CHECK.true(any(any(value for value in row)
                   for row in tempting_differences),
               "full-vector shortcut genuinely differs in E2")
    return {
        "amplitude": F(-63, 8),
        "theta_derivative": F(105, 8),
        "tensor_vector_coefficient": F(105, 16),
        "canonical_pair_gradient_coefficient": F(105, 8),
        "context_count": context_count,
        "energy_profile_histogram": energy_histograms[0],
        "canonical_gradient_histogram": canonical_gradients,
        "shortcut_error_histogram": tempting_differences,
    }


# ---------------------------------------------------------------- lower orders
def central_delta(bits, flips):
    after = list(bits)
    for port in flips:
        after[port] ^= 1
    return vadd(pair_products(tuple(1 - 2 * bit for bit in after)),
                pair_products(tuple(1 - 2 * bit for bit in bits)), F(-1))


def central_v4(bits, neighbor_externals):
    """Complete owner-once h4 first-source coefficient at one locked node."""
    singles = [central_delta(bits, (port,)) for port in range(4)]
    answer = [F(-3, 16) * sum(row[pair] for row in singles)
              for pair in range(6)]
    for left, right in PAIR_ORDER:
        sign_left = 1 if bits[left] == 0 else -1
        sign_right = 1 if bits[right] == 0 else -1
        energy = F((sign_left + sign_right) ** 2 + 2)
        pair_score = central_delta(bits, (left, right))
        for component in range(6):
            answer[component] += (
                (F(1, 2 * energy) - F(3, 16))
                * (singles[left][component] + singles[right][component])
                + pair_score[component] / energy ** 2)
    for port, external_bits in enumerate(neighbor_externals):
        sign_shared = 1 if bits[port] == 0 else -1
        for external_bit in external_bits:
            sign_external = 1 if external_bit == 0 else -1
            energy = F((sign_shared + sign_external) ** 2 + 2)
            coefficient = (F(1, 2 * energy) - F(3, 16)
                           + F(1, energy ** 2))
            for component in range(6):
                answer[component] += coefficient * singles[port][component]
    return tuple(answer)


def lower_order_t_census():
    locked_words = tuple(bits for bits in product((0, 1), repeat=4)
                         if sum(bits) == 2)
    neighborhood_count = 0
    v4_histogram = Counter()
    for bits in locked_words:
        memory = pair_products(tuple(1 - 2 * bit for bit in bits))
        CHECK.equal(project_t(memory), (F(0),) * 6,
                    "every locked pair word has zero T2 projection")
        singles = [central_delta(bits, (port,)) for port in range(4)]
        v2 = vscale(F(1, 4), tuple(sum(row[index] for row in singles)
                                   for index in range(6)))
        CHECK.equal(v2, vscale(F(-1), memory),
                    "owner-once h2 first vertex is -M")
        CHECK.equal(project_t(v2), (F(0),) * 6,
                    "owner-once h2 first vertex has no T2")
        neighbor_words = {
            shared: tuple(row for row in product((0, 1), repeat=3)
                          if sum(row) == 2 - shared)
            for shared in (0, 1)
        }
        choices = [neighbor_words[bits[port]] for port in range(4)]
        for external in product(*choices):
            neighborhood_count += 1
            v4 = central_v4(bits, external)
            expected = tuple(F(-4, 9) - F(37, 12) * value
                             for value in memory)
            CHECK.equal(v4, expected,
                        "owner-once h4 local identity")
            CHECK.equal(project_t(v4), (F(0),) * 6,
                        "owner-once h4 first vertex has no T2")
            v4_histogram[v4] += 1
    CHECK.equal(neighborhood_count, 486,
                "all simple locked radius-one neighborhoods exhausted")
    CHECK.equal(sorted(v4_histogram.values()), [162, 162, 162],
                "h4 identity has three equal locked-pair classes")
    return {
        "locked_local_words": len(locked_words),
        "radius_one_neighborhoods": neighborhood_count,
        "h2_identity": "V2_v=-M_v",
        "h4_identity": "V4_v=-(4/9)1_6-(37/12)M_v",
        "t2_projection": "zero at h0, h2, and h4",
        "h4_histogram": v4_histogram,
    }


# --------------------------------------------------------------- Q4 geometry
PERIOD = 4
STEPS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
CELLS = tuple(product(range(PERIOD), repeat=3))
EDGES = tuple((cell, port) for cell in CELLS for port in range(4))


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


def orientation_and_global_census():
    # A four-cycle in a simple bipartite incidence graph would be a pair of
    # parent cells sharing two distinct child neighbors.  The declared Q4
    # quotient has none, so no distinct locked configuration is reachable in
    # four flips.  This check is deliberately global: the diagonal local h4
    # stencil alone would not exclude an alternating four-cycle in a generic
    # simple 4-regular bipartite graph.
    child_neighbors = {
        cell: frozenset(qadd(cell, step) for step in STEPS)
        for cell in CELLS
    }
    four_cycle_count = sum(
        len(set(child_neighbors[left]) & set(child_neighbors[right]))
        * (len(set(child_neighbors[left]) & set(child_neighbors[right])) - 1) // 2
        for left, right in combinations(CELLS, 2)
    )
    CHECK.equal(four_cycle_count, 0,
                "declared Q4 incidence has no four-cycle")

    cycles = []
    masks = set()
    edge_counts = Counter()
    orientation_counts = Counter()
    u_by_missing = {}
    for cell in CELLS:
        for ports in combinations(range(4), 3):
            cycle = canonical_cycle(cell, ports)
            mask = frozenset(cycle)
            CHECK.true(mask not in masks, "canonical Q4 hexagon is unique")
            masks.add(mask)
            CHECK.equal(len(mask), 6, "canonical Q4 cycle has six links")
            nodes = {}
            for edge in cycle:
                edge_counts[edge] += 1
                for node in endpoints(edge):
                    nodes.setdefault(node, []).append(edge[1])
            CHECK.equal(len(nodes), 6, "canonical Q4 cycle has six vertices")
            CHECK.true(all(len(local_ports) == 2
                           for local_ports in nodes.values()),
                       "cycle has two ports at every vertex")
            vertex_thetas = []
            for local_ports in nodes.values():
                pair = tuple(sorted(local_ports))
                vector = theta(pair)
                CHECK.equal(project_t(vector), vector,
                            "every geometric cycle Theta is pure T2")
                vertex_thetas.append(vector)
            u = tuple(sum(vector[index] for vector in vertex_thetas)
                      for index in range(6))
            missing = next(iter(set(range(4)) - set(ports)))
            expected = vscale(F(2), tuple(sum(theta(pair)[index]
                                                   for pair in combinations(ports, 2))
                                           for index in range(6)))
            CHECK.equal(u, expected,
                        "six vertex Thetas equal twice the orientation triangle")
            if missing in u_by_missing:
                CHECK.equal(u, u_by_missing[missing],
                            "orientation vector is translation independent")
            else:
                u_by_missing[missing] = u
            orientation_counts[missing] += 1
            cycles.append((cycle, ports, missing))

    CHECK.equal(len(cycles), 256, "Q4 has 256 canonical hexagons")
    CHECK.true(all(edge_counts[edge] == 6 for edge in EDGES),
               "each Q4 link lies on six canonical hexagons")
    CHECK.equal(orientation_counts, Counter({0: 64, 1: 64, 2: 64, 3: 64}),
                "four hexagon orientations occur equally")

    expected_coordinates = {
        3: (F(2), F(2), F(-2)),
        2: (F(2), F(-2), F(2)),
        1: (F(-2), F(2), F(2)),
        0: (F(-2), F(-2), F(-2)),
    }
    coordinates = {}
    for missing, vector in u_by_missing.items():
        coordinates[missing] = tuple(dot(vector, basis) / 2
                                     for basis in T_BASIS)
        CHECK.equal(coordinates[missing], expected_coordinates[missing],
                    "orientation vector has tetrahedral T coordinates")
        CHECK.equal(dot(vector, vector), F(24),
                    "orientation vector norm squared is 24")
    CHECK.equal(tuple(sum(u_by_missing[d][index] for d in range(4))
                      for index in range(6)), (F(0),) * 6,
                "four orientation vectors sum to zero")
    for left, right in combinations(range(4), 2):
        CHECK.equal(dot(u_by_missing[left], u_by_missing[right]), F(-8),
                    "distinct orientation-vector inner product is -8")
    CHECK.equal(rank([coordinates[d] for d in range(4)]), 3,
                "orientation vectors span all three T2 directions")

    projector = zero_matrix(6, 6)
    for basis in T_BASIS:
        projector = madd(projector, outer(basis, basis), F(1, 2))
    gram = zero_matrix(6, 6)
    for vector in u_by_missing.values():
        gram = madd(gram, outer(vector, vector))
    CHECK.equal(gram,
                tuple(tuple(F(32) * projector[i][j] for j in range(6))
                      for i in range(6)),
                "orientation sum u_d u_d^T equals 32 P_T")

    # Toggle typing: every mask is real, involutive, and its two alternating
    # directions have the same exact coefficient from the direct census.
    cycle_mask = frozenset(cycles[0][0])
    witness = frozenset(cycles[0][0][::2])
    toggled = witness.symmetric_difference(cycle_mask)
    CHECK.equal(toggled, frozenset(cycles[0][0][1::2]),
                "cycle toggle exchanges the two alternating words")
    CHECK.equal(toggled.symmetric_difference(cycle_mask), witness,
                "formal cycle toggle is an involution")

    return {
        "q4_cells": len(CELLS),
        "q4_links": len(EDGES),
        "q4_hexagons": len(cycles),
        "q4_four_cycles": four_cycle_count,
        "hexagons_per_link": 6,
        "orientation_counts": orientation_counts,
        "u_by_missing": u_by_missing,
        "u_coordinates": coordinates,
        "orientation_rank": 3,
        "orientation_gram": gram,
        "t_projector": projector,
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

    direct = direct_hexagon_census()
    lower = lower_order_t_census()
    geometry = orientation_and_global_census()
    ledger = {
        "lane": "GL6CH",
        "claim": "global order-h6 off-diagonal T2 candidate-field writer",
        "direct_hexagon": direct,
        "lower_order_exclusion": lower,
        "global_geometry": geometry,
        "operator": {
            "source_free": "-(63/8)(h^6/U_d^5) sum_c T_c",
            "first_T_source": "+(105/16)(h^6/U_d^6) sum_c T_c sum_{v in c} j_v.Theta_{v,c}",
            "scope": "candidate-field-dependent future writer; not phase, Ricci, gravity, or G",
        },
    }
    payload = json.dumps(encode(ledger), indent=2, sort_keys=True) + "\n"
    target = HERE / "EXACT_LEDGER.json"
    if args.write_ledger:
        target.write_text(payload)
    else:
        CHECK.true(target.is_file(), "frozen exact ledger exists")
        CHECK.equal(target.read_text(), payload,
                    "frozen exact ledger matches replay")

    print(f"PASS__GL6CH_GLOBAL_H6_TENSOR_WRITER__{CHECK.total}/{CHECK.total}")
    print("HEXAGON_DIRECT=-63/8;CANONICAL_PAIR_GRADIENT=105/8_E_AB")
    print("T2_WRITER=+105/16_THETA;THETA_NORM2=2;BOTH_DIRECTIONS=YES")
    print("LOWER_T_FIRST_VERTICES=H0_H2_H4_ZERO;H6_FIRST_SURVIVOR")
    print("ORIENTATION_TETRAHEDRON=SUM0_NORM24_INNER-8_GRAM32PT_RANK3")
    print("CLAIM=CANDIDATE_FIELD_DEPENDENT_FUTURE_WRITER;NO_PHASE_RICCI_GRAVITY_G")


if __name__ == "__main__":
    main()
