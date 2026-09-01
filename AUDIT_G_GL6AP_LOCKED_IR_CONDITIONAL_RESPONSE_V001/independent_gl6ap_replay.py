#!/usr/bin/env python3
"""Independent exact replay for the GL6AP hostile audit.

This file imports neither the GL6AP nor GL6AN author verifier.
"""

from __future__ import annotations

import cmath
import hashlib
import itertools
import math
from collections import Counter, deque
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


# 1. Independent class-table reconstruction of the S4 statements.
class_sizes = (1, 6, 3, 8, 6)  # 1, (12), (12)(34), (123), (1234)
characters = {
    "A1": (1, 1, 1, 1, 1),
    "A2": (1, -1, 1, 1, -1),
    "E": (2, 0, 2, -1, 0),
    "T1": (3, -1, -1, 0, 1),
    "T2": (3, 1, -1, 0, -1),
}
square_class = (0, 0, 0, 3, 2)
cube_class = (0, 1, 2, 0, 4)


def inner(left, right):
    return Fraction(sum(class_sizes[i] * left[i] * right[i]
                        for i in range(5)), 24)


def sym2(character):
    return tuple((character[i] ** 2 + character[square_class[i]]) // 2
                 for i in range(5))


def sym3(character):
    return tuple((character[i] ** 3
                  + 3 * character[i] * character[square_class[i]]
                  + 2 * character[cube_class[i]]) // 6
                 for i in range(5))


for name, character in characters.items():
    check(inner(character, character) == 1, f"irreducible {name}")
for left, right in itertools.combinations(characters, 2):
    check(inner(characters[left], characters[right]) == 0,
          f"orthogonal irreps {left}/{right}")

chi_e = characters["E"]
chi_t = characters["T2"]
chi_a1 = characters["A1"]
chi_a2 = characters["A2"]
end_e = tuple(value * value for value in chi_e)
sym2e = sym2(chi_e)
sym2t = sym2(chi_t)
sym3e = sym3(chi_e)
sym3t = sym3(chi_t)
check(sym2e == tuple(chi_a1[i] + chi_e[i] for i in range(5)),
      "Sym2(E)=A1+E")
check(sym2t == tuple(chi_a1[i] + chi_e[i] + chi_t[i] for i in range(5)),
      "Sym2(T2)=A1+E+T2")
check(end_e == tuple(chi_a1[i] + chi_a2[i] + chi_e[i] for i in range(5)),
      "End(E)=A1+A2+E")
check(inner(sym2e, sym2t) == 2, "exactly two quadratic contractions")
check(inner(chi_t, end_e) == 0, "no matrix covariant linear in T2")
check(inner(sym3t, chi_a2) == 0, "no reciprocal cubic A2 covariant")
check(inner(sym3e, chi_a1) == 1, "E cubic invariant allowed")
check(inner(end_e, chi_a1) == 1, "E scalar mass allowed")


# 2. Locked pair plane and trivial-character link kernel.
pairs = tuple(itertools.combinations(range(4), 2))
incidence = [[Fraction(a in pair) for pair in pairs] for a in range(4)]


def rank(matrix):
    work = [list(row) for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next((r for r in range(row, len(work))
                      if work[r][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        scale = work[row][column]
        work[row] = [value / scale for value in work[row]]
        for r in range(len(work)):
            if r != row and work[r][column]:
                factor = work[r][column]
                work[r] = [work[r][c] - factor * work[row][c]
                           for c in range(len(work[0]))]
        row += 1
    return row


check(rank(incidence) == 4, "pair incidence rank four")
pair_basis = (
    (1, 0, -1, -1, 0, 1),
    (0, 1, -1, -1, 1, 0),
)
for vector in pair_basis:
    check(all(sum(incidence[a][j] * vector[j] for j in range(6)) == 0
              for a in range(4)), "pair E basis in kernel")
check(rank([list(vector) for vector in pair_basis]) == 2,
      "pair E basis independent")

locked_spins = tuple(z for z in itertools.product((-1, 1), repeat=4)
                     if sum(z) == 0)
opposite = (((0, 1), (2, 3)), ((0, 2), (1, 3)),
            ((0, 3), (1, 2)))
triples = set()
for z in locked_spins:
    values = {(a, b): z[a] * z[b] for a, b in pairs}
    triple = tuple(sum(values[tuple(sorted(pair))] for pair in partition)
                   for partition in opposite)
    check(sum(triple) == -2, "locked opposite-pair affine plane")
    triples.add(triple)
check(len(locked_spins) == 6 and len(triples) == 3,
      "six assignments collapse to three pair-E vertices")


def matmul(left, right):
    return [[sum(left[i][k] * right[k][j] for k in range(len(right)))
             for j in range(len(right[0]))] for i in range(len(left))]


def dagger(matrix):
    return [[matrix[i][j].conjugate() for i in range(len(matrix))]
            for j in range(len(matrix[0]))]


def max_error(left, right):
    return max(abs(left[i][j] - right[i][j])
               for i in range(len(left)) for j in range(len(left[0])))


phase_sets = (
    tuple(cmath.exp(1j * angle) for angle in (0.0, 0.31, 1.17, -0.74)),
    (1 + 0j, -1 + 0j, 1j, -1j),
    tuple(cmath.exp(1j * angle) for angle in (0.2, -0.4, 0.9, -1.1)),
)
for phases in phase_sets:
    s = sum(phases)
    determinant = 16 - abs(s) ** 2
    check(determinant > 1e-10, "generic symbol rank two")
    B = [[1 + 0j] * 4, list(phases)]
    inverse = [[4 / determinant, -s.conjugate() / determinant],
               [-s / determinant, 4 / determinant]]
    correction = matmul(matmul(dagger(B), inverse), B)
    projector = [[complex(i == j) - correction[i][j] for j in range(4)]
                 for i in range(4)]
    check(max_error(matmul(B, projector), [[0j] * 4, [0j] * 4]) < 1e-11,
          "projector is annihilated by B")
    check(max_error(matmul(projector, projector), projector) < 1e-11,
          "projector idempotent")
    check(max_error(projector, dagger(projector)) < 1e-11,
          "projector Hermitian")
    check(abs(sum(projector[i][i] for i in range(4)) - 2) < 1e-11,
          "projector trace two")

trivial_projector = [[complex(i == j) - 0.25 for j in range(4)]
                     for i in range(4)]
check(abs(sum(trivial_projector[i][i] for i in range(4)) - 3) < 1e-12,
      "trivial link kernel rank three")
check(max_error(matmul([[1 + 0j] * 4], trivial_projector), [[0j] * 4]) < 1e-12,
      "trivial link kernel equals centered T2")

theta_directions = ((1, -1, 2, -2), (3, -2, -4, 3), (1, 1, -1, -1))
for theta in theta_directions:
    check(sum(theta) == 0, "centered character direction")
    norm2 = sum(value * value for value in theta)
    for scale in (1e-3, 3e-4):
        s = sum(cmath.exp(1j * scale * value) for value in theta)
        ratio = (4 - abs(s)) / (scale * scale)
        check(abs(ratio - norm2 / 2) < 2e-5,
              "Gram eigenvalue quadratic coefficient")
        singular_ratio = math.sqrt(4 - abs(s)) / scale
        check(abs(singular_ratio - math.sqrt(norm2 / 2)) < 2e-5,
              "singular value linear coefficient")


# 3. Rebuild Q4 and check every S4 affine automorphism independently.
L = 4
cells = tuple(itertools.product(range(L), repeat=3))
directions = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
permutations = tuple(itertools.permutations(range(4)))


def add_mod(left, right):
    return tuple((left[i] + right[i]) % L for i in range(3))


def child(cell, port):
    return add_mod(cell, directions[port])


edges = tuple((cell, port) for cell in cells for port in range(4))
endpoints = {edge: (("P", edge[0]), ("C", child(*edge))) for edge in edges}
incident = {}
for edge, ends in endpoints.items():
    for vertex in ends:
        incident.setdefault(vertex, []).append(edge)
check(len(cells) == 64 and len(edges) == 256 and len(incident) == 128,
      "Q4 census")
check(all(len(row) == 4 for row in incident.values()), "Q4 degree four")
check(all(len(set(ends)) == 2 for ends in endpoints.values()), "Q4 no loops")


def determinant3(matrix):
    return (matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
            - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
            + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]))


def port_matrix(permutation):
    columns = [tuple(directions[permutation[j]][i] - directions[permutation[3]][i]
                     for i in range(3)) for j in range(3)]
    return [[columns[j][i] for j in range(3)] for i in range(3)]


def apply_matrix(matrix, vector):
    return tuple(sum(matrix[i][j] * vector[j] for j in range(3)) % L
                 for i in range(3))


determinants = Counter()
for permutation in permutations:
    matrix = port_matrix(permutation)
    determinant = determinant3(matrix)
    determinants[determinant] += 1
    check(abs(determinant) == 1, "unimodular port map")
    shift = directions[permutation[3]]
    parent_images = {apply_matrix(matrix, cell) for cell in cells}
    child_images = {add_mod(apply_matrix(matrix, cell), shift) for cell in cells}
    check(len(parent_images) == 64 and len(child_images) == 64,
          "Q4 vertex bijection")
    edge_images = set()
    for edge in edges:
        cell, port = edge
        image_edge = (apply_matrix(matrix, cell), permutation[port])
        edge_images.add(image_edge)
        old_child = child(cell, port)
        mapped_child = add_mod(apply_matrix(matrix, old_child), shift)
        check(mapped_child == child(*image_edge), "incidence covariance")
    check(len(edge_images) == 256, "Q4 edge bijection")
check(determinants == Counter({1: 12, -1: 12}), "all 24 determinant signs")


# 4. Independent Edmonds--Karp completion of the locked hexagon collar.
cycle = (
    ((0, 0, 0), 0),
    ((1, 3, 0), 1),
    ((1, 3, 0), 2),
    ((0, 3, 1), 0),
    ((0, 3, 1), 1),
    ((0, 0, 0), 2),
)
cycle_vertices = []
for edge in cycle:
    for vertex in endpoints[edge]:
        if vertex not in cycle_vertices:
            cycle_vertices.append(vertex)
check(len(cycle_vertices) == 6, "target is a six-cycle")
fixed = {edge: int(index % 2 == 0) for index, edge in enumerate(cycle)}
for vertex in cycle_vertices:
    external = sorted(edge for edge in incident[vertex] if edge not in cycle)
    check(len(external) == 2, "two external collar links")
    fixed[external[0]] = 1
    fixed[external[1]] = 0

source, sink = ("SOURCE",), ("SINK",)
capacity = {}


def add_arc(left, right, value):
    capacity.setdefault(left, {})[right] = value
    capacity.setdefault(right, {}).setdefault(left, 0)


parent_fixed = Counter()
child_fixed = Counter()
for edge, value in fixed.items():
    if value:
        parent_fixed[edge[0]] += 1
        child_fixed[child(*edge)] += 1
for cell in cells:
    check(parent_fixed[cell] <= 2 and child_fixed[cell] <= 2,
          "fixed collar locally feasible")
    add_arc(source, ("P", cell), 2 - parent_fixed[cell])
    add_arc(("C", cell), sink, 2 - child_fixed[cell])
for edge in edges:
    if edge not in fixed:
        add_arc(("P", edge[0]), ("C", child(*edge)), 1)

residual = {node: dict(row) for node, row in capacity.items()}
flow = 0
while True:
    predecessor = {source: None}
    queue = deque((source,))
    while queue and sink not in predecessor:
        node = queue.popleft()
        for target, value in residual[node].items():
            if value > 0 and target not in predecessor:
                predecessor[target] = node
                queue.append(target)
    if sink not in predecessor:
        break
    amount = 10**9
    node = sink
    while predecessor[node] is not None:
        amount = min(amount, residual[predecessor[node]][node])
        node = predecessor[node]
    node = sink
    while predecessor[node] is not None:
        prior = predecessor[node]
        residual[prior][node] -= amount
        residual[node][prior] += amount
        node = prior
    flow += amount

required = sum(2 - parent_fixed[cell] for cell in cells)
check(flow == required, "collar extends to global degree-two configuration")
occupied = {edge for edge, value in fixed.items() if value}
for edge in edges:
    if edge not in fixed:
        parent = ("P", edge[0])
        child_node = ("C", child(*edge))
        if residual[child_node].get(parent, 0) == 1:
            occupied.add(edge)
check(len(occupied) == 128, "Q4 has 128 occupied links")
for vertex in incident:
    check(sum(edge in occupied for edge in incident[vertex]) == 2,
          "global degree-two lock")
check(tuple(int(edge in occupied) for edge in cycle) == (1, 0, 1, 0, 1, 0),
      "target cycle alternating")


def local_pair_type(vertex, occupation):
    occupied_ports = frozenset(edge[1] for edge in incident[vertex]
                               if edge in occupation)
    for index, partition in enumerate(opposite):
        if occupied_ports in tuple(map(frozenset, partition)):
            return index
    raise AssertionError("locked vertex lacks pair type")


toggled = occupied.symmetric_difference(cycle)
before = Counter(local_pair_type(vertex, occupied) for vertex in cycle_vertices)
after = Counter(local_pair_type(vertex, toggled) for vertex in cycle_vertices)
before_tuple = tuple(before[index] for index in range(3))
after_tuple = tuple(after[index] for index in range(3))
delta = tuple(after_tuple[index] - before_tuple[index] for index in range(3))
check(before_tuple == (1, 2, 3), "initial E-type counts")
check(after_tuple == (3, 2, 1), "final E-type counts")
check(delta == (2, 0, -2), "nonzero centered E count change")


def subset_energy(indices):
    charge = Counter()
    for index in indices:
        change = -1 if cycle[index] in occupied else 1
        for vertex in endpoints[cycle[index]]:
            charge[vertex] += change
    return sum(value * value for value in charge.values())


coefficient = Fraction(0)
profiles = Counter()
for ordering in itertools.permutations(range(6)):
    contribution = Fraction(1)
    energies = []
    for length in range(1, 6):
        energy = subset_energy(ordering[:length])
        check(energy > 0, "proper cycle prefix exits lock")
        energies.append(energy)
        contribution *= Fraction(-1, energy)
    profiles[tuple(energies)] += 1
    coefficient += contribution
check(coefficient == Fraction(-63, 8), "native hexagon coefficient -63/8")
check(sum(profiles.values()) == math.factorial(6), "all 720 orderings classified")

partitions = tuple(frozenset(map(frozenset, partition)) for partition in opposite)
delta_orbit = set()
for permutation in permutations:
    moved = []
    for partition in partitions:
        image = frozenset(frozenset(permutation[index] for index in pair)
                          for pair in partition)
        moved.append(delta[partitions.index(image)])
    delta_orbit.add(tuple(moved))
check(len(delta_orbit) == 6, "full six-vector E orbit")
check(all(sum(vector) == 0 for vector in delta_orbit), "orbit centered")
check(rank([list(vector) for vector in delta_orbit]) == 2,
      "loop-count orbit spans E")


# 5. Positive-frequency quotient and residue/threshold distinctions.
for support in itertools.combinations(range(1, 7), 3):
    for weights in itertools.product(range(1, 4), repeat=3):
        total = sum(weights)
        moment = sum(Fraction(weight * energy) for weight, energy
                     in zip(weights, support))
        quotient = moment / total
        check(Fraction(min(support)) <= quotient,
              "positive-frequency single-mode upper bound")

positive_gap = Fraction(2)
elastic_weight = Fraction(100)
inelastic_weight = Fraction(1)
full_quotient = positive_gap * inelastic_weight / (elastic_weight + inelastic_weight)
positive_quotient = positive_gap * inelastic_weight / inelastic_weight
check(full_quotient < positive_gap, "elastic atom spoils positive-gap bound")
check(positive_quotient == positive_gap, "positive-frequency quotient restores bound")
pole_sequence = [(Fraction(1, n), Fraction(1, n)) for n in range(1, 20)]
check(all(residue > 0 and frequency > 0 for frequency, residue in pole_sequence),
      "positive atom at every sequence point")
check(pole_sequence[-1][0] < pole_sequence[0][0], "pole energy can approach zero")
check(pole_sequence[-1][1] < pole_sequence[0][1], "residue can vanish independently")


# 6. Fail-closed direct custody and textual ceiling checks.
dependency_rows = [line for line in
                   (ROOT / "LANE_CROSS_RFT_GRA_GL6AP_LOCKED_IR_CONDITIONAL_RESPONSE_V001"
                    / "DEPENDENCIES.sha256").read_text().splitlines()
                   if line.strip()]
check(len(dependency_rows) == 7, "exact GL6AP dependency count")
for row in dependency_rows:
    expected, relative = row.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), "dependency exists")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          "dependency hash exact")
    check("GL6AN" in relative, "dependency confined to GL6AN custody")
    check(all(marker not in relative for marker in ("GL6AL", "GL6AO", "GL6AQ")),
          "mutable and parallel lanes excluded")

theorem = (ROOT / "LANE_CROSS_RFT_GRA_GL6AP_LOCKED_IR_CONDITIONAL_RESPONSE_V001"
           / "THEOREM.md").read_text()
for required_text in (
    "scalar mass `r_E` is symmetry allowed",
    "has a discrete character set and by itself has no literal",
    "a pole sequence can mathematically have residues tending to",
    "physical embedding `X:A3 -> R^3`",
    "no all-orders effective",
    "Nothing in GL6AP assumes or derives a gauge phase",
):
    check(required_text in theorem, "required conditional ceiling")
for forbidden_text in (
    "the E-sector pole is proved",
    "character is physical momentum",
    "symmetry forbids the mass",
    "derives gravity",
):
    check(forbidden_text.lower() not in theorem.lower(), "forbidden promotion absent")

print(f"PASS__GL6AP_INDEPENDENT_REPLAY__{checks}/{checks}")
print("REPRESENTATIONS: KER_B1=T2_DIM3; LOCKED_PAIR=E_DIM2; HOM_S4_ZERO")
print("Q4: ALL_24_AFFINE_AUTOMORPHISMS; DEGREE2_WITNESS; E_ORBIT_SPANS")
print("LOOP: COUNTS_(1,2,3)_TO_(3,2,1); COEFFICIENT_MINUS63_OVER8")
print("RESPONSE: MASS_ALLOWED; TWO_QUADRATIC_COVARIANTS; RECIPROCAL_O_THETA4")
print("SPECTRAL: POSITIVE_FREQUENCY_BOUND; THRESHOLD_ATOM_RESIDUE_DISTINCT")
print("CEILING: NO_POLE_MOMENTUM_CONE_GAUGE_GRAVITY_OR_G_PROMOTION")
