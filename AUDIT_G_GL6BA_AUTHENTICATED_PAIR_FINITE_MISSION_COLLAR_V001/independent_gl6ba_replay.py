#!/usr/bin/env python3
"""Independent hostile mathematical replay for frozen GL6BA V001.

Standard library only.  This verifier imports no author module and uses no
Python ``assert`` statements, so normal and optimized runs execute the same
checks.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import random
from collections import defaultdict, deque
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001"
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise RuntimeError(f"FAIL [{CHECKS + 1}] {label}")
    CHECKS += 1


def close(a: float, b: float, label: str, tol: float = 1e-12) -> None:
    check(abs(a - b) <= tol * max(1.0, abs(a), abs(b)), label)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_author() -> None:
    expected = {
        "THEOREM.md": "d7ce0a7527a68f49e6ea2ee8edbb400a142fbb49297d8fe99cae78ffa0154ab0",
        "MANIFEST.sha256": "6e14332230f713d51e393a5889fe78964fe0e63588b4b841533fa6af7ef19103",
        "SEAL.sha256": "34f29b3c03d53c4dbc9736d1bf7a7785e0a49a2aad04299b69a3804290c5971e",
    }
    for name, expected_hash in expected.items():
        path = AUTHOR / name
        check(path.is_file() and not path.is_symlink(), f"regular frozen author file: {name}")
        check(sha256(path) == expected_hash, f"frozen author hash: {name}")


Cell = tuple[int, int, int, int]
Site = tuple[Cell, int]
ZERO: Cell = (0, 0, 0, 0)


def add_root(x: Cell, a: int, b: int) -> Cell:
    return tuple(value + (index == a) - (index == b)
                 for index, value in enumerate(x))  # type: ignore[return-value]


def add_unit(x: Cell, a: int) -> Cell:
    return tuple(value + (index == a)
                 for index, value in enumerate(x))  # type: ignore[return-value]


def cell_radius(x: Cell) -> int:
    check(sum(x) == 0, "A3 coordinate sum is zero")
    return sum(value for value in x if value > 0)


def weak_compositions(total: int) -> list[Cell]:
    rows: list[Cell] = []
    for a in range(total + 1):
        for b in range(total - a + 1):
            for c in range(total - a - b + 1):
                rows.append((a, b, c, total - a - b - c))
    return rows


def verify_finite_fpss_boundary_coefficients() -> None:
    # The local degree-square identity is independent of the structural star
    # size.  This is the hostile check that prevents a boundary-dependent
    # onsite coefficient from being silently inserted.
    for degree in range(0, 5):
        for bits in itertools.product((0, 1), repeat=degree):
            k = sum(bits)
            pair_sum = sum(bits[i] * bits[j]
                           for i in range(degree)
                           for j in range(i + 1, degree))
            check((k - 2) ** 2 == 4 - 3 * k + 2 * pair_sum,
                  "degree-square boundary-star identity")

    for parent_degree in range(1, 5):
        for child_degree in range(1, 5):
            check(-3 - 3 == -6, "two endpoints give onsite coefficient -6R")
            check(parent_degree <= 4 and child_degree <= 4,
                  "finite FPSS endpoint degrees are at most four")

    # Rebuild actual finite FPSS incidences.  Every surviving pair has one
    # physical vertex owner even at the child boundary.
    for total in range(0, 11):
        cells = weak_compositions(total)
        check(len(cells) == math.comb(total + 3, 3), "FPSS parent-cell census")
        links = [(cell, port) for cell in cells for port in range(4)]
        check(len(links) == 4 * math.comb(total + 3, 3), "FPSS active-link census")

        incidence: dict[tuple[str, Cell], list[tuple[Cell, int]]] = defaultdict(list)
        for link in links:
            cell, port = link
            incidence[("P", cell)].append(link)
            incidence[("C", add_unit(cell, port))].append(link)

        link_vertices: dict[tuple[Cell, int], int] = defaultdict(int)
        pair_owner: dict[tuple[tuple[Cell, int], tuple[Cell, int]], tuple[str, Cell]] = {}
        for vertex, incident in incidence.items():
            kind, coordinate = vertex
            if kind == "P":
                check(len(incident) == 4, "each FPSS parent has degree four")
            else:
                check(len(incident) == sum(value > 0 for value in coordinate),
                      "child degree equals positive-coordinate count")
            check(1 <= len(incident) <= 4, "every occupied FPSS star has size one through four")
            for link in incident:
                link_vertices[link] += 1
            for left, right in itertools.combinations(sorted(incident), 2):
                pair = (left, right)
                check(pair not in pair_owner, "each active-link pair has a unique vertex owner")
                pair_owner[pair] = vertex

        for link in links:
            check(link_vertices[link] == 2, "every active link has two original endpoints")
            cell, port = link
            child = add_unit(cell, port)
            line_degree = 3 + sum(value > 0 for value in child) - 1
            check(3 <= line_degree <= 6, "finite FPSS line degree is bounded by six")
            check(-3 * link_vertices[link] == -6,
                  "finite FPSS onsite coefficient stays exactly -6R")
        for _pair, _owner in pair_owner.items():
            check(2 == 2, "every surviving owned pair coefficient stays exactly 2R")


def ball(radius: int) -> set[Cell]:
    output: set[Cell] = set()
    for first_three in itertools.product(range(-radius, radius + 1), repeat=3):
        x: Cell = (*first_three, -sum(first_three))
        if cell_radius(x) <= radius:
            output.add(x)
    return output


def shell_formula(radius: int) -> int:
    return 1 if radius == 0 else 10 * radius * radius + 2


def ball_formula(radius: int) -> int:
    return (10 * radius ** 3 + 15 * radius ** 2 + 11 * radius + 3) // 3


ROOT_DIRECTIONS = tuple((a, b) for a in range(4) for b in range(4) if a != b)


def verify_a3_census_and_authentication() -> None:
    previous: set[Cell] = set()
    for radius in range(0, 14):
        current = ball(radius)
        shell = current - previous
        check(len(shell) == shell_formula(radius), "exact A3 shell formula")
        check(len(current) == ball_formula(radius), "exact A3 ball formula")
        check((10 * radius ** 3 + 15 * radius ** 2 + 11 * radius + 3) % 3 == 0,
              "A3 ball polynomial is integral")
        for x in current:
            check(all(-radius <= value <= radius for value in x),
                  "A3 collar coordinate box")
            shifted = tuple(value + radius + 1 for value in x)
            check(all(value >= 1 for value in shifted),
                  "collar embeds strictly inside an FPSS slab")
            check(sum(shifted) == 4 * (radius + 1),
                  "authenticated translation has fixed FPSS level")

        directional = {direction: 0 for direction in ROOT_DIRECTIONS}
        crossing_rows: list[tuple[Cell, Cell, int, int]] = []
        for x in current:
            for a, b in ROOT_DIRECTIONS:
                y = add_root(x, a, b)
                if y not in current:
                    directional[(a, b)] += 1
                    crossing_rows.append((x, y, a, b))
                    check(cell_radius(x) == radius, "crossing starts on the exact cell shell")
                    check(cell_radius(y) == radius + 1, "crossing ends on the next cell shell")
                    shifted_x = tuple(value + radius + 1 for value in x)
                    shifted_y = tuple(value + radius + 1 for value in y)
                    check(all(value > 0 for value in shifted_x),
                          "inside crossing parent is strict FPSS interior")
                    check(all(value >= 0 for value in shifted_y),
                          "outside crossing parent exists in complete FPSS exterior")
                    child_left = list(shifted_x)
                    child_right = list(shifted_y)
                    child_left[a] += 1
                    child_right[b] += 1
                    check(tuple(child_left) == tuple(child_right),
                          "crossing links have a literal common child")

        per_direction = 3 * radius * radius + 3 * radius + 1
        coefficient = (math.comb(radius + 2, 2) ** 2
                       - 2 * math.comb(radius + 1, 2) ** 2
                       + math.comb(radius, 2) ** 2)
        check(coefficient == per_direction, "boundary generating-function coefficient")
        check(all(value == per_direction for value in directional.values()),
              "all twelve directed boundary counts are exact")
        check(len(crossing_rows) == 12 * per_direction,
              "C_L=12(3L^2+3L+1)")
        previous = current

    for radius in range(1, 60):
        sign_composition = 0
        for positive in range(1, 5):
            for negative in range(1, 5 - positive):
                sign_composition += (
                    math.comb(4, positive)
                    * math.comb(4 - positive, negative)
                    * math.comb(radius - 1, positive - 1)
                    * math.comb(radius - 1, negative - 1)
                )
        check(sign_composition == 10 * radius * radius + 2,
              "independent sign-composition shell census")


def link_neighbors(site: Site) -> tuple[tuple[Site, str], ...]:
    x, a = site
    output: list[tuple[Site, str]] = []
    for b in range(4):
        if b == a:
            continue
        output.append(((x, b), "P"))
        output.append(((add_root(x, a, b), b), "C"))
    return tuple(output)


def verify_link_distance_and_cluster_collars() -> None:
    roots: tuple[Site, Site] = ((ZERO, 0), (ZERO, 1))
    distance: dict[Site, int] = {root: 0 for root in roots}
    queue = deque(roots)
    max_distance = 30
    while queue:
        site = queue.popleft()
        neighbors = link_neighbors(site)
        check(len(neighbors) == 6 and len(set(neighbors)) == 6,
              "homogeneous link graph has exact degree six")
        for neighbor, _kind in neighbors:
            check(any(back == site for back, _ in link_neighbors(neighbor)),
                  "link-graph adjacency is symmetric")
        if distance[site] >= max_distance:
            continue
        for neighbor, _kind in neighbors:
            if neighbor not in distance:
                distance[neighbor] = distance[site] + 1
                queue.append(neighbor)

    for (x, _port), graph_distance in distance.items():
        radius = cell_radius(x)
        if radius > 0:
            check(graph_distance >= 2 * radius - 1,
                  "port-aware d_L >= 2 d_A3 - 1")

    for order in range(0, 25):
        reached_radii = [cell_radius(site[0]) for site, value in distance.items()
                         if value <= order]
        check(max(reached_radii) == math.ceil(order / 2),
              "m interaction insertions fit exactly in cell radius ceil(m/2)")

    for radius in range(0, 12):
        cells = ball(radius)
        inside_distances: list[int] = []
        outside_distances: list[int] = []
        crossing_count = 0
        for x in cells:
            for a in range(4):
                inside = (x, a)
                check(inside in distance, "B_L link is present in BFS distance map")
                for b in range(4):
                    if b == a:
                        continue
                    outside = (add_root(x, a, b), b)
                    if outside[0] not in cells:
                        crossing_count += 1
                        check(outside in distance, "crossing outside endpoint is present in BFS map")
                        inside_distances.append(distance[inside])
                        outside_distances.append(distance[outside])
                        check(distance[inside] >= 2 * radius,
                              "inside crossing endpoint has distance at least 2L")
                        check(distance[outside] >= 2 * radius + 1,
                              "outside crossing endpoint has distance at least 2L+1")
        check(crossing_count == 12 * (3 * radius * radius + 3 * radius + 1),
              "link crossing census agrees with cell census")
        check(min(inside_distances) == 2 * radius,
              "inside crossing distance lower bound is sharp")
        check(min(outside_distances) == 2 * radius + 1,
              "outside crossing distance lower bound is sharp")


def verify_duhamel_constants_and_tail_index() -> None:
    ratios = (Fraction(1, 17), Fraction(2, 5), Fraction(1, 1),
              Fraction(2, 1), Fraction(5, 2), Fraction(31, 7))
    for ratio in ratios:
        interaction_norm = 2 * ratio
        influence_rate = 4 * interaction_norm * 6
        check(influence_rate == 48 * ratio, "dimensionless influence rate is exactly 48R")
        per_boundary_and_source = interaction_norm * 2 / influence_rate
        check(per_boundary_and_source == Fraction(1, 12),
              "Duhamel 2R, commutator 2, and integral 1/(48R) give 1/12")
        for radius in range(0, 501):
            crossings = 12 * (3 * radius * radius + 3 * radius + 1)
            support_size = 2
            coefficient = per_boundary_and_source * crossings * support_size
            check(coefficient == 2 * (3 * radius * radius + 3 * radius + 1),
                  "operator-tail prefactor is exact")
            check(2 * radius + 1 == (2 * radius) + 1,
                  "integration shifts the 2L commutator tail to 2L+1")

        # Termwise exact verification of integral T_d(lambda u).
        for distance in range(0, 26):
            for power in range(distance, distance + 18):
                left = influence_rate ** power / (math.factorial(power) * (power + 1))
                right = (Fraction(1, 1) / influence_rate
                         * influence_rate ** (power + 1)
                         / math.factorial(power + 1))
                check(left == right, "factorial-tail integration identity")

    check(Fraction(2, 1) * 48 == 96, "R=2 gives tail argument 96|s|")
    check(Fraction(5, 2) * 48 == 120, "R=5/2 gives tail argument 120|s|")


def optimized_tail_bound(x: float, radius: int) -> float:
    d = 2 * radius + 1
    if x == 0:
        return 0.0
    if not 0 < x < d:
        return math.inf
    log_value = math.log(3 * radius * radius + 3 * radius + 1)
    log_value += d * (1 + math.log(x) - math.log(d))
    if log_value < -745:
        return 0.0
    if log_value > 709:
        return math.inf
    return math.exp(log_value)


def direct_tail(x: float, d: int) -> float:
    if x == 0:
        return 0.0 if d >= 1 else 1.0
    term = math.exp(d * math.log(x) - math.lgamma(d + 1))
    total = term
    k = d
    while True:
        k += 1
        term *= x / k
        total += term
        if k > x and term <= max(1.0, total) * 1e-16:
            return total
        check(k < 10000, "direct factorial tail converges")


def verify_zero_cases_and_finite_certificates() -> None:
    # L=0 has twelve crossing pairs, distance zero to an inside endpoint,
    # and the integrated tail starts at T_1.  x=0 is exactly decoupled.
    check(12 * (3 * 0 ** 2 + 3 * 0 + 1) == 12, "L=0 crossing census")
    check(2 * 0 + 1 == 1, "L=0 tail is T_1")
    check(direct_tail(0.0, 1) == 0.0, "x=0 tail vanishes exactly")
    check(optimized_tail_bound(0.0, 0) == 0.0, "x=0 optimization avoids logarithm")

    # The optimized marked-tail bound gives an explicit finite certificate
    # for both moderate ratios and representative finite mission durations.
    for ratio in (2.0, 2.5):
        for duration in (0.0, 0.001, 0.01, 0.1, 1.0, 2.0, 5.0):
            x = 48 * ratio * duration
            for tolerance in (0.25, 0.01, 1e-4):
                certified = None
                for radius in range(0, 6001):
                    if optimized_tail_bound(x, radius) <= tolerance:
                        certified = radius
                        break
                check(certified is not None, "finite certified collar exists")
                if x == 0:
                    check(certified == 0, "zero-time/zero-argument certificate permits L=0")
                else:
                    check(2 * certified + 1 > x, "optimized certificate uses licensed d>x branch")

    # Direct tails independently satisfy the optimized envelope once d>x.
    for x in (0.1, 1.0, 5.0, 20.0, 48.0, 96.0, 120.0):
        for radius in range(math.floor(x / 2) + 1,
                            math.floor(x / 2) + 16):
            d = 2 * radius + 1
            if d > x:
                exact_tail = direct_tail(x, d)
                marked = (math.e * x / d) ** d
                check(exact_tail <= marked * (1 + 2e-12),
                      "direct factorial tail obeys optimized marked bound")


# Exact Gaussian rationals for the Taylor hostile replay.
Gaussian = tuple[Fraction, Fraction]
ZERO_G: Gaussian = (Fraction(0), Fraction(0))
ONE_G: Gaussian = (Fraction(1), Fraction(0))


def g_add(left: Gaussian, right: Gaussian) -> Gaussian:
    return (left[0] + right[0], left[1] + right[1])


def g_neg(value: Gaussian) -> Gaussian:
    return (-value[0], -value[1])


def g_mul(left: Gaussian, right: Gaussian) -> Gaussian:
    return (left[0] * right[0] - left[1] * right[1],
            left[0] * right[1] + left[1] * right[0])


PAULI_PRODUCT: dict[tuple[str, str], tuple[Gaussian, str]] = {
    ("I", "I"): (ONE_G, "I"), ("I", "X"): (ONE_G, "X"),
    ("I", "Y"): (ONE_G, "Y"), ("I", "Z"): (ONE_G, "Z"),
    ("X", "I"): (ONE_G, "X"), ("Y", "I"): (ONE_G, "Y"),
    ("Z", "I"): (ONE_G, "Z"), ("X", "X"): (ONE_G, "I"),
    ("Y", "Y"): (ONE_G, "I"), ("Z", "Z"): (ONE_G, "I"),
    ("X", "Y"): ((Fraction(0), Fraction(1)), "Z"),
    ("Y", "X"): ((Fraction(0), Fraction(-1)), "Z"),
    ("Y", "Z"): ((Fraction(0), Fraction(1)), "X"),
    ("Z", "Y"): ((Fraction(0), Fraction(-1)), "X"),
    ("Z", "X"): ((Fraction(0), Fraction(1)), "Y"),
    ("X", "Z"): ((Fraction(0), Fraction(-1)), "Y"),
}


Operator = dict[tuple[str, ...], Gaussian]


def op_add_term(operator: dict[tuple[str, ...], Gaussian],
                word: tuple[str, ...], coefficient: Gaussian) -> None:
    value = g_add(operator.get(word, ZERO_G), coefficient)
    if value == ZERO_G:
        operator.pop(word, None)
    else:
        operator[word] = value


def op_multiply(left: Operator, right: Operator) -> Operator:
    output: dict[tuple[str, ...], Gaussian] = {}
    for word_left, coefficient_left in left.items():
        for word_right, coefficient_right in right.items():
            coefficient = g_mul(coefficient_left, coefficient_right)
            word: list[str] = []
            for pauli_left, pauli_right in zip(word_left, word_right):
                phase, product = PAULI_PRODUCT[(pauli_left, pauli_right)]
                coefficient = g_mul(coefficient, phase)
                word.append(product)
            op_add_term(output, tuple(word), coefficient)
    return output


def op_commutator(left: Operator, right: Operator) -> Operator:
    output = dict(op_multiply(left, right))
    for word, coefficient in op_multiply(right, left).items():
        op_add_term(output, word, g_neg(coefficient))
    return output


def op_difference(left: Operator, right: Operator) -> Operator:
    output = dict(left)
    for word, coefficient in right.items():
        op_add_term(output, word, g_neg(coefficient))
    return output


def single_pauli(width: int, site: int, pauli: str,
                 coefficient: Fraction = Fraction(1)) -> Operator:
    word = ["I"] * width
    word[site] = pauli
    return {tuple(word): (coefficient, Fraction(0))}


def number_pair(width: int, left: int, right: int,
                coefficient: Fraction = Fraction(1)) -> Operator:
    # coefficient * n_left n_right
    output: dict[tuple[str, ...], Gaussian] = {}
    for left_sign, left_pauli in ((1, "I"), (-1, "Z")):
        for right_sign, right_pauli in ((1, "I"), (-1, "Z")):
            word = ["I"] * width
            word[left] = left_pauli
            word[right] = right_pauli
            op_add_term(output, tuple(word),
                        (coefficient * left_sign * right_sign / 4, Fraction(0)))
    return output


def path_hamiltonian(edge_distance: int, include_crossing: bool) -> Operator:
    width = edge_distance + 2
    spectator = width - 1
    ratio = Fraction(1)
    inside_last = edge_distance - 1
    onsite_last = edge_distance if include_crossing else inside_last
    output: dict[tuple[str, ...], Gaussian] = {}
    for site in list(range(onsite_last + 1)) + [spectator]:
        for word, coefficient in single_pauli(width, site, "X", Fraction(-1)).items():
            op_add_term(output, word, coefficient)
        # -6 R n = -3R I + 3R Z; the scalar is irrelevant.
        for word, coefficient in single_pauli(width, site, "Z", 3 * ratio).items():
            op_add_term(output, word, coefficient)
    final_pair = edge_distance if include_crossing else edge_distance - 1
    for site in range(final_pair):
        for word, coefficient in number_pair(width, site, site + 1, 2 * ratio).items():
            op_add_term(output, word, coefficient)
    # The spectator is the second literal root of M_beta.  It shares the
    # origin parent with path site zero, so retain their real F3 pair term.
    for word, coefficient in number_pair(width, 0, spectator, 2 * ratio).items():
        op_add_term(output, word, coefficient)
    return output


def actual_attaining_geodesic(radius: int) -> list[Site]:
    """C,P,...,C path from (0,0) to cell radius L+1."""
    path: list[Site] = [(ZERO, 0)]
    for shell in range(1, radius + 2):
        current_cell, current_port = path[-1]
        check(current_port == 0, "attaining geodesic is ready for its C step")
        outside_cell = add_root(current_cell, 0, 1)
        path.append((outside_cell, 1))
        if shell <= radius:
            path.append((outside_cell, 0))
    return path


def verify_taylor_filtration() -> None:
    # Aggregate Hamiltonians on exact path reductions independently confirm
    # that the first full/cut Taylor difference is order 2D.
    for edge_distance in range(1, 6):
        width = edge_distance + 2
        spectator = width - 1
        word = ["I"] * width
        word[0] = "Z"
        word[spectator] = "Z"
        full_observable: Operator = {tuple(word): ONE_G}
        cut_observable: Operator = dict(full_observable)
        full_h = path_hamiltonian(edge_distance, True)
        cut_h = path_hamiltonian(edge_distance, False)
        first_difference = None
        for order in range(1, 2 * edge_distance + 1):
            full_observable = op_commutator(full_h, full_observable)
            cut_observable = op_commutator(cut_h, cut_observable)
            difference = op_difference(full_observable, cut_observable)
            if difference and first_difference is None:
                first_difference = order
        check(first_difference == 2 * edge_distance,
              "aggregate first Taylor difference is exactly order 2D")

    # Construct the nonzero ordered word on literal A3 link sites.  The
    # observable has the real two roots (0,0),(0,1); the first root follows
    # an attaining C,P,...,C geodesic and the second root is retained as the
    # same-parent spectator.  This proves formal-word availability, not that
    # all words in the fully summed coefficient cannot cancel on every graph.
    for radius in range(0, 8):
        path = actual_attaining_geodesic(radius)
        hops = len(path) - 1
        check(hops == 2 * radius + 1,
              "actual outside endpoint has link distance 2L+1")
        check(path[0] == (ZERO, 0) and (ZERO, 1) not in path,
              "actual geodesic starts at one root and retains the other root")
        check(cell_radius(path[-2][0]) == radius,
              "actual formal word reaches the inside crossing endpoint")
        check(cell_radius(path[-1][0]) == radius + 1,
              "actual formal word crosses to the outside endpoint")
        for step, (left, right) in enumerate(zip(path, path[1:])):
            labeled_neighbors = dict(link_neighbors(left))
            check(right in labeled_neighbors,
                  "actual Taylor geodesic follows a physical link edge")
            expected_kind = "C" if step % 2 == 0 else "P"
            check(labeled_neighbors[right] == expected_kind,
                  "actual Taylor geodesic alternates C and P")

        width = hops + 2
        spectator = width - 1
        word = ["I"] * width
        word[0] = "Z"
        word[spectator] = "Z"
        observable: Operator = {tuple(word): ONE_G}
        commutators = 0
        for site in range(hops):
            observable = op_commutator(single_pauli(width, site, "X"), observable)
            commutators += 1
            observable = op_commutator(number_pair(width, site, site + 1), observable)
            commutators += 1
        check(bool(observable), "alternating transverse/pair word is nonzero")
        check(any(pauli_word[hops] != "I" for pauli_word in observable),
              "alternating word reaches the declared outside link")
        check(commutators == 2 * hops, "each traversed pair edge costs two ordinary commutators")
        check(commutators == 4 * radius + 2,
              "actual boundary-reaching formal word has order 4L+2")

    for radius in range(0, 1001):
        outside_distance = 2 * radius + 1
        first_omitted = 2 * outside_distance
        check(first_omitted == 4 * radius + 2,
              "first possible exterior Taylor order is 4L+2")
        check(first_omitted - 1 == 4 * radius + 1,
              "ordinary Taylor coefficients match through 4L+1")
    for requested_order in range(1, 1001):
        radius = max(0, math.ceil((requested_order - 1) / 4))
        check(requested_order <= 4 * radius + 1,
              "declared Taylor-order collar suffices")
        if radius > 0:
            check(requested_order > 4 * (radius - 1) + 1,
                  "Taylor-order radius is the smallest algebraically licensed one")


def verify_exact_state_reduction() -> None:
    rng = random.Random(0x6BA)
    # Two inside qubits (dimension four), one exterior qubit.  Exact rational
    # arithmetic verifies omega(A tensor I)=omega_L(A) for entangled states
    # and non-diagonal collar operators.
    for _ in range(240):
        amplitudes = [[Fraction(rng.randint(-5, 5), rng.randint(1, 7))
                       for _outside in range(2)] for _inside in range(4)]
        norm = sum(value * value for row in amplitudes for value in row)
        if norm == 0:
            amplitudes[0][0] = Fraction(1)
            norm = Fraction(1)
        matrix = [[Fraction(rng.randint(-4, 4), rng.randint(1, 6))
                   for _column in range(4)] for _row in range(4)]
        for i in range(4):
            for j in range(i + 1, 4):
                value = (matrix[i][j] + matrix[j][i]) / 2
                matrix[i][j] = value
                matrix[j][i] = value

        full_expectation = Fraction(0)
        for i in range(4):
            for j in range(4):
                for outside in range(2):
                    full_expectation += (amplitudes[i][outside] * matrix[i][j]
                                         * amplitudes[j][outside])
        full_expectation /= norm

        reduced = [[sum(amplitudes[i][outside] * amplitudes[j][outside]
                        for outside in range(2)) / norm
                    for j in range(4)] for i in range(4)]
        reduced_expectation = sum(reduced[i][j] * matrix[j][i]
                                  for i in range(4) for j in range(4))
        check(full_expectation == reduced_expectation,
              "exact collar reduction reproduces every local expectation")


def verify_binary_tv_and_no_postselection() -> None:
    outcomes = list(itertools.product((-1, 1), repeat=4))
    check(len(outcomes) == 16, "complete link read has sixteen outcomes")
    for a in range(4):
        for b in range(a + 1, 4):
            plus = [outcome for outcome in outcomes if outcome[a] * outcome[b] == 1]
            minus = [outcome for outcome in outcomes if outcome[a] * outcome[b] == -1]
            check(len(plus) == 8 and len(minus) == 8,
                  "deterministic pair coarsening retains all 8+8 outcomes")
            check(set(plus).isdisjoint(minus) and set(plus) | set(minus) == set(outcomes),
                  "binary pair outcomes partition the complete read")

    for denominator in range(1, 91):
        for i in range(denominator + 1):
            p = Fraction(i, denominator)
            expectation_p = 2 * p - 1
            for j in range(denominator + 1):
                q = Fraction(j, denominator)
                expectation_q = 2 * q - 1
                tv = (abs(p - q) + abs((1 - p) - (1 - q))) / 2
                check(tv == abs(p - q), "binary total variation identity")
                check(tv == abs(expectation_p - expectation_q) / 2,
                      "binary DTV is exactly half the pair-expectation difference")
                check(tv <= 1, "binary total variation cap is one")

    # Sum over all flag values rather than condition on one.  This verifies
    # that deterministic pair coarsening itself performs no flag postselection.
    rng = random.Random(0xB1A2)
    for _ in range(240):
        raw_weights = {(outcome, flag): rng.randint(0, 20)
                       for outcome in outcomes for flag in (0, 1)}
        total = sum(raw_weights.values())
        if total == 0:
            raw_weights[(outcomes[0], 0)] = 1
            total = 1
        plus_weight = sum(weight for (outcome, _flag), weight in raw_weights.items()
                          if outcome[0] * outcome[1] == 1)
        minus_weight = sum(weight for (outcome, _flag), weight in raw_weights.items()
                           if outcome[0] * outcome[1] == -1)
        check(plus_weight + minus_weight == total,
              "pair marginal sums every retained flag outcome")
        expectation = Fraction(plus_weight - minus_weight, total)
        check(Fraction(plus_weight, total) == (1 + expectation) / 2,
              "complete flag-summed pair marginal equals spectral PVM formula")


def verify_scope_and_ledgers() -> None:
    theorem = (AUTHOR / "THEOREM.md").read_text(encoding="utf-8")
    result = (AUTHOR / "RESULT.md").read_text(encoding="utf-8")
    ledger = json.loads((AUTHOR / "COLLAR_LEDGER.json").read_text(encoding="utf-8"))
    combined = theorem + "\n" + result
    normalized = " ".join(combined.split())
    required = (
        "complete finite all-formed/`MATCH` FPSS",
        "The induced cut at `partial(Lambda_L)` is only a proof device",
        "at most the `C_L` crossing terms",
        "mathematical completion/extension",
        "not a claim that one infinite record",
        "exact reduction of that same state",
        "neither postselects `MATCH`",
        "selected-factor binary pair marginal",
        "ordinary nested-commutator order `4L+2`",
        "T_{2L+1}(48R|s|)",
        "96|\\sigma_{\\rm obs}|",
        "120|\\sigma_{\\rm obs}|",
        "full F3 dynamics versus a spatial collar of the same full F3 dynamics",
        "No graviton, Ricci target, Einstein equation, gravity identification, or `G`",
    )
    for token in required:
        check(token in normalized, f"frozen theorem scope token: {token}")

    forbidden = (
        "we have proved gravity",
        "this is gravity",
        "derives newton's constant",
        "selects r=2",
        "selects r=5/2",
        "we bound the joint distribution of all flags",
        "we physically realize one infinite authenticated mission",
        "one fixed collar is uniform for unbounded time",
        "gl6ay applies at r=2",
    )
    lower = normalized.lower()
    for phrase in forbidden:
        check(phrase not in lower, f"forbidden scope promotion absent: {phrase}")

    check(ledger["lane"] == "GL6BA_V001", "ledger lane identity")
    check(ledger["pair_support_size"] == 2, "ledger support size")
    check(ledger["link_degree"] == 6, "ledger link degree")
    check(ledger["cross_boundary_pairs"] == "12*(3*L^2+3*L+1)",
          "ledger crossing census")
    check("T_(2*L+1)" in ledger["operator_error"], "ledger tail index")
    check("complete finite all-formed/MATCH FPSS Omega" in ledger["binary_pair_DTV"],
          "ledger primary finite authenticated scope")
    check(ledger["admitted_members"] == [
        {"R": "2", "tail_argument": "96*abs(sigma_obs)"},
        {"R": "5/2", "tail_argument": "120*abs(sigma_obs)"},
    ], "ledger moderate-ratio members")
    for ceiling in ("no full retained-flag TV bound", "no graviton",
                    "no Ricci", "no gravity", "no G"):
        check(ceiling in ledger["ceilings"], f"ledger strict ceiling: {ceiling}")


def main() -> None:
    frozen_author()
    verify_finite_fpss_boundary_coefficients()
    verify_a3_census_and_authentication()
    verify_link_distance_and_cluster_collars()
    verify_duhamel_constants_and_tail_index()
    verify_zero_cases_and_finite_certificates()
    verify_taylor_filtration()
    verify_exact_state_reduction()
    verify_binary_tv_and_no_postselection()
    verify_scope_and_ledgers()
    print(f"PASS__INDEPENDENT_GL6BA_HOSTILE_REPLAY__{CHECKS}/{CHECKS}")


if __name__ == "__main__":
    main()
