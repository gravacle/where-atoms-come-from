#!/usr/bin/env python3
"""Independent exact checks for GL6AS native hexagon collective response."""

from __future__ import annotations

import hashlib
import itertools
import math
import re
from collections import Counter
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


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def matmul(left, right):
    return [[sum(left[i][k] * right[k][j] for k in range(len(right)))
             for j in range(len(right[0]))] for i in range(len(left))]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def rank_exact(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [work[row][j] - scale * work[pivot_row][j]
                         for j in range(columns)]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def nullspace(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    rows = len(work)
    columns = len(work[0])
    pivot_columns = []
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row:
                continue
            scale = work[row][column]
            if scale:
                work[row] = [work[row][j] - scale * work[pivot_row][j]
                             for j in range(columns)]
        pivot_columns.append(column)
        pivot_row += 1
    free_columns = [column for column in range(columns)
                    if column not in pivot_columns]
    basis = []
    for free in free_columns:
        vector = [Fraction(0) for _ in range(columns)]
        vector[free] = 1
        for row, pivot in enumerate(pivot_columns):
            vector[pivot] = -work[row][free]
        basis.append(tuple(vector))
    return tuple(basis)


# 1. S4 representation products and selection rules.
permutations = tuple(itertools.permutations(range(4)))
partitions = (
    frozenset((frozenset((0, 1)), frozenset((2, 3)))),
    frozenset((frozenset((0, 2)), frozenset((1, 3)))),
    frozenset((frozenset((0, 3)), frozenset((1, 2)))),
)


def compose_power(permutation, exponent):
    result = tuple(range(4))
    for _ in range(exponent):
        result = tuple(permutation[result[index]] for index in range(4))
    return result


def fixed_partitions(permutation):
    count = 0
    for partition in partitions:
        moved = frozenset(
            frozenset(permutation[index] for index in pair)
            for pair in partition
        )
        count += moved == partition
    return count


def inner(left, right):
    return Fraction(sum(a * b for a, b in zip(left, right)), 24)


trivial = tuple(1 for _ in permutations)
chi_ports = tuple(sum(permutation[i] == i for i in range(4))
                  for permutation in permutations)
chi_t2 = tuple(value - 1 for value in chi_ports)
chi_e = tuple(fixed_partitions(permutation) - 1
              for permutation in permutations)


def symmetric_square(character):
    values = []
    for index, permutation in enumerate(permutations):
        squared = permutations.index(compose_power(permutation, 2))
        values.append((character[index] ** 2 + character[squared]) // 2)
    return tuple(values)


sym2_t2 = symmetric_square(chi_t2)
check(inner(chi_t2, chi_t2) == 1, "T2 irreducible")
check(inner(chi_e, chi_e) == 1, "E irreducible")
check(inner(chi_t2, chi_e) == 0, "Hom(T2,E)=0")
check(inner(chi_ports, chi_e) == 0, "four-source representation has no E")
check(inner(sym2_t2, trivial) == 1, "Sym2(T2) has A1")
check(inner(sym2_t2, chi_e) == 1, "Sym2(T2) has E")
check(inner(sym2_t2, chi_t2) == 1, "Sym2(T2) has T2")
check(inner(tuple(value * value for value in chi_t2), chi_e) == 1,
      "theta times T2 source has one E")


# 2. Exact Q4 hexagon census and port-charge conservation.
L = 4
cells = tuple(itertools.product(range(L), repeat=3))
links = tuple((x, port) for x in cells for port in range(4))


def child(x, port):
    value = list(x)
    if port < 3:
        value[port] = (value[port] + 1) % L
    return tuple(value)


endpoints = {
    edge: ((0, edge[0]), (1, child(*edge)))
    for edge in links
}
adjacency = {}
for edge, (left, right) in endpoints.items():
    adjacency.setdefault(left, []).append((right, edge))
    adjacency.setdefault(right, []).append((left, edge))
for vertex in adjacency:
    adjacency[vertex].sort()

cycles = {}


def walk_cycle(start, current, visited, path_edges):
    if len(path_edges) == 6:
        return
    for neighbor, edge in adjacency[current]:
        if neighbor == start:
            if len(path_edges) == 5:
                ordered = tuple(path_edges + (edge,))
                cycles.setdefault(frozenset(ordered), ordered)
            continue
        if neighbor in visited or len(path_edges) == 5:
            continue
        walk_cycle(start, neighbor, visited | {neighbor}, path_edges + (edge,))


for start in sorted(adjacency):
    walk_cycle(start, start, {start}, tuple())

check(len(cycles) == 256, "Q4 has 256 undirected hexagons")
link_cycle_count = Counter()
missing_count = Counter()
for key, ordered in cycles.items():
    check(len(key) == 6, "simple six-edge cycle")
    ports = [edge[1] for edge in ordered]
    counts = Counter(ports)
    check(sorted(counts.values()) == [2, 2, 2],
          "hexagon uses three port labels twice")
    used = set(counts)
    missing = next(port for port in range(4) if port not in used)
    missing_count[missing] += 1
    for port in used:
        parities = {index % 2 for index, value in enumerate(ports)
                    if value == port}
        check(parities == {0, 1}, "each port occurs once on each alternating half")
    for edge in key:
        link_cycle_count[edge] += 1
check(tuple(missing_count[port] for port in range(4)) == (64, 64, 64, 64),
      "four hexagon orientations have equal count")
check(all(link_cycle_count[edge] == 6 for edge in links),
      "every link lies on six hexagons")


# 3. Exact cycle symbol, kernel, rank, and leading norm.
triples = tuple(itertools.combinations(range(4), 3))


def cycle_symbol(z):
    columns = []
    for a, b, c in triples:
        column = [z[0] * 0 for _ in range(4)]
        column[a] = z[b] - z[c]
        column[b] = z[c] - z[a]
        column[c] = z[a] - z[b]
        columns.append(tuple(column))
    return [list(row) for row in zip(*columns)]


for z in (
    tuple(map(Fraction, (1, 2, 3, 4))),
    tuple(map(Fraction, (2, -1, 5, 3))),
    tuple(map(Fraction, (1, 1, -1, -1))),
    tuple(map(Fraction, (1, 1, 1, -1))),
):
    C = cycle_symbol(z)
    check(all(sum(C[row][column] for row in range(4)) == 0
              for column in range(4)), "cycle symbol has zero A1 boundary")
    check(all(sum(z[row] * C[row][column] for row in range(4)) == 0
              for column in range(4)), "cycle symbol satisfies child constraint")
    check(rank_exact(C) == 2, "nontrivial cycle symbol rank two")
    B = [[Fraction(1) for _ in range(4)], list(z)]
    check(all(value == 0 for row in matmul(B, C) for value in row),
          "B times C vanishes exactly")


def leading_cycle(theta):
    return cycle_symbol(tuple(map(Fraction, theta)))


centered_theta = (
    (1, -1, 0, 0),
    (1, 2, -1, -2),
    (3, -1, -4, 2),
    (5, -2, 1, -4),
)
for theta_values in centered_theta:
    theta = tuple(map(Fraction, theta_values))
    check(sum(theta) == 0, "centered theta")
    C1 = leading_cycle(theta)
    lhs = matmul(C1, transpose(C1))
    norm2 = dot(theta, theta)
    projector_t = [
        [Fraction(int(i == j)) - Fraction(1, 4) for j in range(4)]
        for i in range(4)
    ]
    rhs = [
        [4 * (norm2 * projector_t[i][j] - theta[i] * theta[j])
         for j in range(4)]
        for i in range(4)
    ]
    check(lhs == rhs, "C1 C1T exact transverse identity")
    check(rank_exact(C1) == 2, "leading cycle symbol rank two")
    check(all(sum(row) == 0 for row in lhs), "leading symbol removes A1")
    check(all(dot(row, theta) == 0 for row in lhs),
          "leading symbol removes longitudinal theta")


# 4. Strict local pair E composite and exact overlap.
pairs = tuple(itertools.combinations(range(4), 2))
R = [[Fraction(int(port in pair)) for pair in pairs] for port in range(4)]
e_basis = nullspace(R)
check(rank_exact(R) == 4, "pair incidence rank four")
check(len(e_basis) == 2, "pair E kernel dimension two")

locked_z = tuple(z for z in itertools.product((-1, 1), repeat=4) if sum(z) == 0)
pair_vectors = []
for z in locked_z:
    M = tuple(Fraction(z[a] * z[b]) for a, b in pairs)
    pair_vectors.append(M)
    check(tuple(dot(row, M) for row in R) == (-1, -1, -1, -1),
          "locked affine pair identity")
    e = tuple(Fraction(-value, 2) for value in z)
    check(all(M[index] == 4 * e[a] * e[b]
              for index, (a, b) in enumerate(pairs)),
          "authenticated pair read is exact quadratic composite")
    for partition in partitions:
        first, second = tuple(partition)
        first_index = pairs.index(tuple(sorted(first)))
        second_index = pairs.index(tuple(sorted(second)))
        check(M[first_index] == M[second_index],
              "strict local pair T2 component vanishes")

for alpha, beta in itertools.product(range(-3, 4), repeat=2):
    if alpha == beta == 0:
        continue
    c = tuple(alpha * e_basis[0][index] + beta * e_basis[1][index]
              for index in range(6))
    values = tuple(dot(c, M) for M in pair_vectors)
    check(sum(values) == 0, "locked E read has zero uniform mean")
    mean_square = sum(value * value for value in values) / 6
    check(mean_square == Fraction(8, 3) * dot(c, c),
          "locked E read has exact 8/3 norm")

# Every one-link exchange between degree-two configurations gives the AQ
# local E displacement and isolated-doublet overlap.
for initial in itertools.combinations(range(4), 2):
    initial_set = set(initial)
    for removed in initial:
        for inserted in set(range(4)) - initial_set:
            final_set = (initial_set - {removed}) | {inserted}
            zi = tuple(-1 if port in initial_set else 1 for port in range(4))
            zf = tuple(-1 if port in final_set else 1 for port in range(4))
            Mi = tuple(Fraction(zi[a] * zi[b]) for a, b in pairs)
            Mf = tuple(Fraction(zf[a] * zf[b]) for a, b in pairs)
            delta = tuple(Mf[index] - Mi[index] for index in range(6))
            check(all(dot(row, delta) == 0 for row in R),
                  "loop-vertex displacement lies in E")
            check(dot(delta, delta) == 16, "loop-vertex E displacement norm")
            check(dot(delta, Mf) - dot(delta, Mi) == 16,
                  "matched E read distinguishes endpoints by 16")
            check(abs(Fraction(dot(delta, Mi) - dot(delta, Mf), 2)) == 8,
                  "isolated doublet E transition amplitude eight")


# 5. Retained-lineage maps and complement parity.
retained_map = [
    [Fraction(0 if row == column else 2) for column in range(4)]
    for row in range(4)
]
ones = (Fraction(1),) * 4
check(tuple(dot(row, ones) for row in retained_map) == (6, 6, 6, 6),
      "retained map A1 eigenvalue six")
for vector in (
    (1, -1, 0, 0),
    (1, 0, -1, 0),
    (1, 0, 0, -1),
):
    moved = tuple(dot(row, vector) for row in retained_map)
    check(moved == tuple(-2 * value for value in vector),
          "retained map T2 eigenvalue minus two")

for support in itertools.product((0, 1), repeat=6):
    product = math.prod(support)
    check(product == int(all(support)), "six-support product gate")

for z, M in zip(locked_z, pair_vectors):
    complemented = tuple(-value for value in z)
    Mc = tuple(Fraction(complemented[a] * complemented[b]) for a, b in pairs)
    check(Mc == M, "pair E read is complement even")
    e = tuple(Fraction(-value, 2) for value in z)
    ec = tuple(Fraction(-value, 2) for value in complemented)
    check(ec == tuple(-value for value in e), "link density is complement odd")


# 6. Exact two-state oscillator-strength and spectral-channel factors.
T = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
H = [[-value for value in row] for row in T]  # J=1
D = [[Fraction(0), Fraction(0)], [Fraction(0), Fraction(1)]]


def subtract(left, right):
    return [[left[i][j] - right[i][j] for j in range(len(left[0]))]
            for i in range(len(left))]


def commutator(left, right):
    return subtract(matmul(left, right), matmul(right, left))


double = commutator(D, commutator(H, D))
plus_expectation = sum(double[i][j] for i in range(2) for j in range(2)) / 2
f_value = plus_expectation / 2
check(f_value == Fraction(1, 2),
      "double-commutator factor J times t times delta squared over two")

O = [[Fraction(-8), Fraction(0)], [Fraction(0), Fraction(8)]]
minus = (Fraction(1), Fraction(-1))
plus = (Fraction(1), Fraction(1))
transition = sum(plus[i] * O[i][j] * minus[j]
                 for i in range(2) for j in range(2)) / 2
check(abs(transition) == 8, "isolated doublet pair-read overlap eight")
loop_transition = sum(plus[i] * T[i][j] * minus[j]
                      for i in range(2) for j in range(2)) / 2
check(loop_transition == 0, "own loop-amplitude source is doublet diagonal")


# 7. Fail-closed ancestry and theorem ceilings.
dependency_rows = [
    line for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
    if line.strip()
]
check(len(dependency_rows) == 17, "exact dependency count")
dependency_paths = set()
for line in dependency_rows:
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), "dependency exists: " + relative)
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          "dependency hash: " + relative)
    dependency_paths.add(relative)
check(len(dependency_paths) == 17, "dependency paths unique")
check(sum("GL6AO" in path for path in dependency_paths) == 6,
      "six AO custody objects")
check(sum("GL6AP" in path for path in dependency_paths) == 6,
      "six AP custody objects")
check(sum("GL6AQ" in path for path in dependency_paths) == 5,
      "five AQ custody objects")
check(not any(path.endswith(
    "LANE_CROSS_RFT_GRA_GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION_V001/SEAL.sha256"
) for path in dependency_paths), "no invented AQ author seal")

raw_theorem = (HERE / "THEOREM.md").read_text()
theorem = " ".join(raw_theorem.replace("**", "").split())
check(re.search(r"(?<!\\)\b(?:quad|qquad)\b", raw_theorem) is None,
      "no malformed TeX spacing commands")
check("\\hexagon" not in raw_theorem, "no undefined hexagon TeX command")
for token in (
    "[H_{\\rm hex},N_a]=0",
    "\\operatorname{im}C(\\chi)=\\ker B(\\chi)",
    "C_1(\\theta)C_1(\\theta)^T",
    "only charge/continuity-supported soft candidate",
    "conditional diagnostic, not a dispersion theorem",
    "\\omega_1^2=\\omega_2^2",
    "cannot be set equal to the bare",
    "exact quadratic composite",
    "\\langle {\\rm one}\\text{-}T_2|O(c)|0\\rangle=0",
    "F_c^{\\lambda\\mu}(q,k-q)",
    "not calibrated physical momenta",
    "support thresholds only if the form factor",
    "K_{E\\leftarrow\\mathrm{loop}}(\\omega,1)=0",
    "\\sum_c\\tau_c=-H_{\\rm hex}/J",
    "\\operatorname{Sym}^2_0(T_2)=E\\oplus T_2",
    "algebraic tensor-like composite channel",
    "leading Fock overlap is therefore a two-mode channel",
    "normalized product trace remains an exact universal counterexample",
    "Nothing here assumes a conventional gauge phase",
):
    check(token in theorem, "required theorem scope: " + token)

aggregate = " ".join(
    (HERE / name).read_text().lower()
    for name in ("THEOREM.md", "README.md", "RESULT.md", "SELF_AUDIT.md")
)
for forbidden in (
    "the t2 mode is proved gapless",
    "the e pole is proved",
    "is a photon",
    "is a graviton",
    "is gravity",
    "derives a physical cone",
    "derives newton's constant",
    "fixed q4 has an infrared limit",
    "g=j",
):
    check(forbidden not in aggregate, "forbidden promotion absent: " + forbidden)

print(f"PASS__GL6AS_NATIVE_HEXAGON_COLLECTIVE__{checks}/{checks}")
print("COLLECTIVE: CENTERED_PORT_CHARGE=T2; GENERIC_LOCKED_BUNDLE=TRANSVERSE_RANK2")
print("CYCLE_SYMBOL: IMAGE_EQUALS_KERNEL; C1C1T=4(THETA2_PT-THETA_THETAT)")
print("SOFTNESS: OSCILLATOR_STRENGTH_O(THETA2); STRUCTURE_FACTOR_AND_STATE_REQUIRED")
print("PAIR_READ: EXACT_E_COMPOSITE; PH_EVEN; GAUSSIAN_TWO_T2_CHANNEL")
print("RETAINED: A1_PLUS_T2; ZERO_CHARACTER_E_CROSS=0; SIX_SUPPORT_PRODUCT")
print("TENSOR: SYM2_0_T2=E_PLUS_T2_ALGEBRA_ONLY; NO_GRAVITY_PROMOTION")
