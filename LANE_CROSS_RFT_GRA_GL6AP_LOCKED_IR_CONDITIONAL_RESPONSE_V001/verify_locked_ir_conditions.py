#!/usr/bin/env python3
"""Exact checks for the GL6AP locked-sector infrared conditional theorem."""

from __future__ import annotations

import cmath
import hashlib
import itertools
import math
import re
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


def inner_product(left, right):
    return Fraction(sum(a * b for a, b in zip(left, right)), len(left))


# 1. Exact S4 character algebra.
permutations = tuple(itertools.permutations(range(4)))
partitions = (
    frozenset((frozenset((0, 1)), frozenset((2, 3)))),
    frozenset((frozenset((0, 2)), frozenset((1, 3)))),
    frozenset((frozenset((0, 3)), frozenset((1, 2)))),
)


def power(permutation, exponent):
    result = tuple(range(4))
    for _ in range(exponent):
        result = tuple(permutation[result[index]] for index in range(4))
    return result


def fixed_ports(permutation):
    return sum(permutation[index] == index for index in range(4))


def apply_partition(permutation, partition):
    return frozenset(frozenset(permutation[index] for index in pair)
                     for pair in partition)


def fixed_partitions(permutation):
    return sum(apply_partition(permutation, partition) == partition
               for partition in partitions)


trivial = tuple(1 for _ in permutations)
chi_t = tuple(fixed_ports(p) - 1 for p in permutations)
chi_e = tuple(fixed_partitions(p) - 1 for p in permutations)


def permutation_sign(permutation):
    inversions = sum(permutation[i] > permutation[j]
                     for i in range(4) for j in range(i + 1, 4))
    return -1 if inversions % 2 else 1


chi_a2 = tuple(permutation_sign(p) for p in permutations)
check(inner_product(chi_t, chi_t) == 1, "T2 irreducible character")
check(inner_product(chi_e, chi_e) == 1, "E irreducible character")
check(inner_product(chi_t, chi_e) == 0, "Hom_S4(T2,E)=0")
check(inner_product(trivial, chi_t) == 0, "T2 excludes A1")
check(inner_product(trivial, chi_e) == 0, "E excludes A1")


def symmetric_square(character):
    out = []
    for index, p in enumerate(permutations):
        p2 = power(p, 2)
        j = permutations.index(p2)
        out.append((character[index] ** 2 + character[j]) // 2)
    return tuple(out)


def symmetric_cube(character):
    out = []
    for index, p in enumerate(permutations):
        p2 = permutations.index(power(p, 2))
        p3 = permutations.index(power(p, 3))
        value = (character[index] ** 3
                 + 3 * character[index] * character[p2]
                 + 2 * character[p3]) // 6
        out.append(value)
    return tuple(out)


sym2e = symmetric_square(chi_e)
sym2t = symmetric_square(chi_t)
sym3e = symmetric_cube(chi_e)
sym3t = symmetric_cube(chi_t)
end_e = tuple(value * value for value in chi_e)
check(inner_product(sym2e, trivial) == 1, "Sym2(E) has one A1")
check(inner_product(sym2e, chi_e) == 1, "Sym2(E) has one E")
check(inner_product(sym2e, chi_t) == 0, "Sym2(E) has no T2")
check(inner_product(sym2t, trivial) == 1, "Sym2(T2) has one A1")
check(inner_product(sym2t, chi_e) == 1, "Sym2(T2) has one E")
check(inner_product(sym2t, chi_t) == 1, "Sym2(T2) has one T2")
check(inner_product(sym2e, sym2t) == 2, "two quadratic spatial invariants")
check(inner_product(sym3e, trivial) == 1, "one E cubic invariant")
check(inner_product(end_e, trivial) == 1, "End(E) has one A1")
check(inner_product(end_e, chi_a2) == 1, "End(E) has one A2")
check(inner_product(end_e, chi_e) == 1, "End(E) has one E")
check(inner_product(end_e, chi_t) == 0, "End(E) has no T2 linear invariant")
check(inner_product(sym3t, chi_a2) == 0,
      "reciprocal cubic antisymmetric invariant absent")


# 2. Incidence symbol ranks, projector support, and scaling type.
def matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def dagger(a):
    return [list(row) for row in zip(*[[value.conjugate() for value in row]
                                      for row in a])]


def identity(n):
    return [[complex(i == j) for j in range(n)] for i in range(n)]


def max_error(a, b):
    return max(abs(a[i][j] - b[i][j])
               for i in range(len(a)) for j in range(len(a[0])))


phase_samples = (
    (1 + 0j, 1 + 0j, -1 + 0j, -1 + 0j),
    (1 + 0j, 1j, -1 + 0j, -1j),
    (1 + 0j, 1 + 0j, 1 + 0j, -1 + 0j),
)
for phases in phase_samples:
    s = sum(phases)
    B = [[1 + 0j] * 4, list(phases)]
    gram_rows = [[4 + 0j, s.conjugate()], [s, 4 + 0j]]
    determinant = 16 - abs(s) ** 2
    inverse = [[4 / determinant, -s.conjugate() / determinant],
               [-s / determinant, 4 / determinant]]
    projector = identity(4)
    correction = matmul(matmul(dagger(B), inverse), B)
    projector = [[projector[i][j] - correction[i][j] for j in range(4)]
                 for i in range(4)]
    check(max_error(matmul(B, projector), [[0j] * 4, [0j] * 4]) < 1e-12,
          "symbol projector lies in kernel")
    check(max_error(matmul(projector, projector), projector) < 1e-12,
          "symbol projector idempotent")
    check(abs(sum(projector[i][i] for i in range(4)) - 2) < 1e-12,
          "generic kernel rank two")
    check(abs((4 - abs(s)) + (4 + abs(s)) - 8) < 1e-12,
          "Gram eigenvalue trace")
    check(abs((4 - abs(s)) * (4 + abs(s)) - determinant) < 1e-12,
          "Gram eigenvalue determinant")

trivial_projector = [[complex(i == j) - Fraction(1, 4) for j in range(4)]
                     for i in range(4)]
check(abs(sum(trivial_projector[i][i] for i in range(4)) - 3) < 1e-12,
      "trivial kernel rank three")
ones = [[1 + 0j] * 4]
check(max_error(matmul(ones, trivial_projector), [[0j] * 4]) < 1e-12,
      "trivial kernel removes A1")

centered = (
    (1, -1, 0, 0),
    (1, 2, -1, -2),
    (Fraction(3, 2), Fraction(-1, 2), -2, 1),
)
for theta in centered:
    check(sum(theta) == 0, "centered character tangent")
    norm2 = sum(value * value for value in theta)
    pair2 = sum((theta[a] - theta[b]) ** 2
                for a, b in itertools.combinations(range(4), 2))
    check(pair2 == 4 * norm2, "centered pair-square identity")
    t = 1e-5
    s = sum(cmath.exp(1j * t * float(value)) for value in theta)
    lam = 4 - abs(s)
    check(abs(lam / t**2 - float(norm2) / 2) < 2e-5,
          "Gram eigenvalue quadratic")
    check(abs(math.sqrt(lam) / t - math.sqrt(float(norm2) / 2)) < 2e-5,
          "singular value linear")


# 3. Locked pair E plane and opposite-pair harmonics.
pairs = tuple(itertools.combinations(range(4), 2))
locked = tuple(z for z in itertools.product((-1, 1), repeat=4) if sum(z) == 0)


def pair_type_from_z(z):
    occupied = frozenset(index for index, value in enumerate(z) if value == -1)
    for index, partition in enumerate(partitions):
        if occupied in partition:
            return index
    raise AssertionError("locked occupation lacks opposite-pair type")


type_vectors = set()
type_count = Counter()
for z in locked:
    vector = tuple(z[a] * z[b] for a, b in pairs)
    p = (vector[0] + vector[5], vector[1] + vector[4],
         vector[2] + vector[3])
    check(sum(p) == -2, "locked opposite-pair affine plane")
    type_vectors.add(p)
    type_count[pair_type_from_z(z)] += 1
check(len(locked) == 6, "six locked spin assignments")
check(len(type_vectors) == 3, "three particle-hole-even pair types")
check(tuple(type_count[index] for index in range(3)) == (2, 2, 2),
      "two complements per E type")

# The quadratic character harmonics transform through the same three
# opposite-pair permutation representation.
sample_theta = (2, -1, 3, -4)


def r_triple(theta):
    return (theta[0] * theta[1] + theta[2] * theta[3],
            theta[0] * theta[2] + theta[1] * theta[3],
            theta[0] * theta[3] + theta[1] * theta[2])


for permutation in permutations:
    moved_theta = tuple(sample_theta[permutation[index]] for index in range(4))
    moved_r = r_triple(moved_theta)
    original_r = r_triple(sample_theta)
    expected = []
    for partition in partitions:
        inverse_image = frozenset(
            frozenset(permutation[index] for index in pair)
            for pair in partition
        )
        expected.append(original_r[partitions.index(inverse_image)])
    check(tuple(expected) == moved_r, "quadratic harmonic transforms as E+A1")


# 4. Q4 locked background and native-loop nonconservation of uniform E.
L = 4
cells = tuple(itertools.product(range(L), repeat=3))
links = tuple((x, port) for x in cells for port in range(4))


def child(x, port):
    y = list(x)
    if port < 3:
        y[port] = (y[port] + 1) % L
    return tuple(y)


directions = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))


def port_linear_map(permutation, x):
    base = directions[permutation[3]]
    columns = tuple(
        tuple((directions[permutation[j]][i] - base[i]) % L
              for i in range(3))
        for j in range(3)
    )
    return tuple(sum(columns[j][i] * x[j] for j in range(3)) % L
                 for i in range(3))


# Exact S4 covariance of the declared period-four quotient.  Parents map by
# A_sigma and children by A_sigma followed by d_{sigma(4)}.
for permutation in permutations:
    parent_images = {port_linear_map(permutation, x) for x in cells}
    check(len(parent_images) == len(cells), "Q4 port map is invertible")
    child_shift = directions[permutation[3]]
    covariance = True
    for x in cells:
        moved_parent = port_linear_map(permutation, x)
        for port in range(4):
            old_child = child(x, port)
            moved_child = tuple(
                (port_linear_map(permutation, old_child)[i] + child_shift[i]) % L
                for i in range(3)
            )
            if moved_child != child(moved_parent, permutation[port]):
                covariance = False
    check(covariance, "Q4 incidence is S4 covariant")


endpoints = {edge: (("P", edge[0]), ("C", child(*edge))) for edge in links}
incident = {}
for edge, ends in endpoints.items():
    for vertex in ends:
        incident.setdefault(vertex, []).append(edge)

cycle = (
    ((0, 0, 0), 0),
    ((1, 3, 0), 1),
    ((1, 3, 0), 2),
    ((0, 3, 1), 0),
    ((0, 3, 1), 1),
    ((0, 0, 0), 2),
)
fixed = {edge: int(index % 2 == 0) for index, edge in enumerate(cycle)}
cycle_vertices = []
for edge in cycle:
    for vertex in endpoints[edge]:
        if vertex not in cycle_vertices:
            cycle_vertices.append(vertex)
check(len(cycle_vertices) == 6, "six target hexagon vertices")

# Fix one external occupied and one external empty link at every cycle vertex.
# This chooses a complete local assignment; the max flow proves it extends to
# a global degree-two Q4 background.
for vertex in cycle_vertices:
    external = sorted(edge for edge in incident[vertex] if edge not in cycle)
    check(len(external) == 2, "two external links at hexagon vertex")
    fixed[external[0]] = 1
    fixed[external[1]] = 0

source, sink = ("S",), ("T",)
capacity = {}


def add_edge(u, v, cap):
    capacity.setdefault(u, {})[v] = cap
    capacity.setdefault(v, {}).setdefault(u, 0)


p_fixed, c_fixed = Counter(), Counter()
for edge, value in fixed.items():
    if value:
        p_fixed[edge[0]] += 1
        c_fixed[child(*edge)] += 1
for x in cells:
    check(p_fixed[x] <= 2 and c_fixed[x] <= 2, "fixed assignment locally feasible")
    add_edge(source, ("P", x), 2 - p_fixed[x])
    add_edge(("C", x), sink, 2 - c_fixed[x])
for edge in links:
    if edge not in fixed:
        add_edge(("P", edge[0]), ("C", child(*edge)), 1)

# Dinic flow, distinct from GL6AN's author background construction.
residual = {u: dict(row) for u, row in capacity.items()}


def levels():
    level = {source: 0}
    queue = deque((source,))
    while queue:
        u = queue.popleft()
        for v, cap in residual[u].items():
            if cap and v not in level:
                level[v] = level[u] + 1
                queue.append(v)
    return level


def send(u, amount, level, cursor):
    if u == sink:
        return amount
    keys = tuple(residual[u])
    while cursor[u] < len(keys):
        v = keys[cursor[u]]
        if residual[u][v] and level.get(v) == level[u] + 1:
            pushed = send(v, min(amount, residual[u][v]), level, cursor)
            if pushed:
                residual[u][v] -= pushed
                residual[v][u] += pushed
                return pushed
        cursor[u] += 1
    return 0


flow = 0
while True:
    level = levels()
    if sink not in level:
        break
    cursor = {u: 0 for u in residual}
    while True:
        pushed = send(source, 10**9, level, cursor)
        if not pushed:
            break
        flow += pushed
required = sum(2 - p_fixed[x] for x in cells)
check(flow == required, "fixed hexagon assignment extends to locked Q4")

occupied = {edge for edge, value in fixed.items() if value}
for edge in links:
    if edge in fixed:
        continue
    pnode = ("P", edge[0])
    cnode = ("C", child(*edge))
    if residual[cnode].get(pnode, 0) == 1:
        occupied.add(edge)
check(len(occupied) == 128, "two occupied links per Q4 parent")
for vertex, edges in incident.items():
    check(sum(edge in occupied for edge in edges) == 2, "global degree-two lock")
check(tuple(int(edge in occupied) for edge in cycle) == (1, 0, 1, 0, 1, 0),
      "target hexagon alternating")


def vertex_type(vertex, occupation):
    occupied_ports = frozenset(edge[1] for edge in incident[vertex]
                               if edge in occupation)
    for index, partition in enumerate(partitions):
        if occupied_ports in partition:
            return index
    raise AssertionError("locked vertex has no E type")


toggled = occupied.symmetric_difference(cycle)
before = Counter(vertex_type(vertex, occupied) for vertex in cycle_vertices)
after = Counter(vertex_type(vertex, toggled) for vertex in cycle_vertices)
before_tuple = tuple(before[index] for index in range(3))
after_tuple = tuple(after[index] for index in range(3))
check(before_tuple == (1, 2, 3), "initial local E-type counts")
check(after_tuple == (3, 2, 1), "toggled local E-type counts")
check(tuple(after_tuple[i] - before_tuple[i] for i in range(3)) == (2, 0, -2),
      "native hexagon changes uniform E count")

# The conserved-covector kernel is S4 invariant.  The orbit of this nonzero
# centered E vector spans the full centered plane, so irreducibility excludes
# every nonzero conserved uniform E combination.
delta = tuple(after_tuple[i] - before_tuple[i] for i in range(3))
delta_orbit = set()
for permutation in permutations:
    moved = []
    for partition in partitions:
        image = apply_partition(permutation, partition)
        moved.append(delta[partitions.index(image)])
    delta_orbit.add(tuple(moved))
check(all(sum(vector) == 0 for vector in delta_orbit),
      "loop-count orbit remains in E plane")
check(any(delta[0] * vector[1] != delta[1] * vector[0]
          for vector in delta_orbit),
      "loop-count orbit spans both E directions")

# Independently replay the nonzero hexagon coefficient needed in the
# commutator witness.
initial_cycle = tuple(int(edge in occupied) for edge in cycle)


def subset_energy(subset):
    charge = Counter()
    for index in subset:
        delta = 1 if initial_cycle[index] == 0 else -1
        for vertex in endpoints[cycle[index]]:
            charge[vertex] += delta
    return sum(value * value for value in charge.values())


all_indices = frozenset(range(6))
amplitude = {frozenset(): Fraction(1)}
for size in range(1, 6):
    for values in itertools.combinations(range(6), size):
        subset = frozenset(values)
        energy = subset_energy(subset)
        check(energy > 0, "proper hexagon subset leaves lock")
        amplitude[subset] = -sum(amplitude[subset - {edge}] for edge in subset) / energy
coefficient = sum(amplitude[all_indices - {edge}] for edge in all_indices)
check(coefficient == Fraction(-63, 8), "nonzero native hexagon coefficient")


# 5. Custody and strict claim ceilings.
dependency_rows = [line for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
                   if line.strip()]
check(len(dependency_rows) == 7, "exact dependency count")
for line in dependency_rows:
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), "dependency exists: " + relative)
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          "dependency hash: " + relative)

raw_theorem = (HERE / "THEOREM.md").read_text()
theorem = " ".join(raw_theorem.replace("**", "").split())
check(re.search(r"(?<!\\)\b(?:quad|qquad)\b", raw_theorem) is None,
      "no malformed TeX spacing commands")
for token in (
    "Hom}_{S_4}(T_2,E)=0",
    "(N_1,N_2,N_3)_i=(1,2,3)",
    "(N_1,N_2,N_3)_f=(3,2,1)",
    "scalar mass `r_E` is symmetry allowed",
    "Sym}^3(E)\\supset A_1",
    "no nonzero uniform linear combination of pair-`E` is conserved",
    "a pole sequence can mathematically have residues tending to zero",
    "The restriction to `(0,infinity)` removes any elastic zero-frequency atom",
    "checks all 24 graph automorphisms",
    "has a discrete character set and by itself has no literal",
    "does not insert a global locked projector",
    "no all-orders effective Hamiltonian is being asserted here",
    "A vanishing threshold alone can be a continuum and is not a pole",
    "physical embedding `X:A3 -> R^3`",
    "Nothing in GL6AP assumes or derives a gauge phase",
):
    check(token in theorem, "required theorem scope: " + token)
for forbidden in (
    "the E-sector pole is proved",
    "character is physical momentum",
    "quadratically soft singular value",
    "six-record operation",
    "[H_{\\rm eff},\\mathbf N_E]\\ne0",
):
    check(forbidden.lower() not in theorem.lower(),
          "forbidden promotion absent: " + forbidden)

print(f"GL6AP exact verification: PASS ({checks}/{checks})")
print("REPRESENTATIONS: LINK_TRIVIAL=T2_DIM3; PAIR_LOCKED=E_DIM2; HOM_S4=0")
print("LOOP: NATIVE_HEXAGON_CHANGES_UNIFORM_E_COUNTS_BY_(2,0,-2); NO_E_WARD")
print("QUADRATIC: MASS_ALLOWED; TWO_SPATIAL_INVARIANTS; CUBIC_E_INVARIANT_ALLOWED")
print("SPECTRAL: GAP_GAPLESS_POLE_REQUIRE_SELECTED_STATE_AND_COMPLETE_HEFF")
print("CEILING: NO_GAUGE_PHASE_PHOTON_GRAVITON_MOMENTUM_CONE_GRAVITY_G")
