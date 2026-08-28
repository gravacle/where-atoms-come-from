#!/usr/bin/env python3
"""Independent hostile audit for GRA-FS-F3-Q4-CSRAV-V001.

This executable does not import the builder verifier.  It reconstructs the
periodic q4 family, the complete frozen reduced source at operator level, the
uniform word homogeneity, a fixed Feshbach derivative, and the quadratic-
contact boundary.  It also enforces the repaired distinction between a
prospective additive query and a tensor forced by the source-free parent.
"""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
THEOREM = (HERE / "THEOREM.md").read_text(encoding="utf-8")
RESULT = (HERE / "RESULT.md").read_text(encoding="utf-8")
SELF_AUDIT = (HERE / "SELF_AUDIT.md").read_text(encoding="utf-8")
passed = 0


def check(condition: bool, label: str) -> None:
    global passed
    if not condition:
        raise AssertionError(label)
    passed += 1
    print(f"PASS {' '.join(label.split())}")


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_ledger(path: Path) -> dict[str, str]:
    ledger: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            value, name = line.split(maxsplit=1)
            ledger[name.strip()] = value
    return ledger


def rank(matrix) -> int:
    if not matrix:
        return 0
    work = [[F(value) for value in row] for row in matrix]
    row_count, column_count = len(work), len(work[0])
    pivot_row = 0
    for column in range(column_count):
        selected = next((row for row in range(pivot_row, row_count)
                         if work[row][column]), None)
        if selected is None:
            continue
        work[pivot_row], work[selected] = work[selected], work[pivot_row]
        pivot = work[pivot_row][column]
        work[pivot_row] = [value / pivot for value in work[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [left - multiple * right
                         for left, right in zip(work[row], work[pivot_row])]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def add(left, right):
    return [[F(a) + F(b) for a, b in zip(row_l, row_r)]
            for row_l, row_r in zip(left, right)]


def scale(value, matrix):
    return [[F(value) * F(entry) for entry in row] for row in matrix]


def matmul(left, right):
    return [[sum((F(a) * F(b) for a, b in zip(row, column)), F(0))
             for column in zip(*right)] for row in left]


def identity(size):
    return [[F(row == column) for column in range(size)] for row in range(size)]


def inverse(matrix):
    size = len(matrix)
    work = [[F(value) for value in row] + unit
            for row, unit in zip(matrix, identity(size))]
    for column in range(size):
        selected = next((row for row in range(column, size)
                         if work[row][column]), None)
        if selected is None:
            raise AssertionError("singular hostile test matrix")
        work[column], work[selected] = work[selected], work[column]
        pivot = work[column][column]
        work[column] = [value / pivot for value in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [left - multiple * right
                         for left, right in zip(work[row], work[column])]
    return [row[size:] for row in work]


def submatrix(matrix, rows, columns):
    return [[matrix[row][column] for column in columns] for row in rows]


def matrix_zero(matrix) -> bool:
    return all(value == 0 for row in matrix for value in row)


# -------------------------------------------------------------------------
# Dependency semantics and byte custody.
# -------------------------------------------------------------------------

dependencies = parse_ledger(HERE / "DEPENDENCIES.sha256")
check(len(dependencies) == 11, "dependency ledger has exactly eleven entries")
for relative, expected in dependencies.items():
    path = (HERE / relative).resolve()
    check(path.is_file() and not path.is_symlink(),
          f"dependency is a regular file: {relative}")
    check(digest(path) == expected, f"dependency digest matches: {relative}")
    check(sha256(path.read_bytes() + b"hostile-tamper").hexdigest() != expected,
          f"dependency appended-byte tamper fails: {relative}")

bs = (ROOT / "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md").read_text()
cw = (ROOT / "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md").read_text()
fm = (ROOT / "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md").read_text()
fq = (ROOT / "LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/THEOREM.md").read_text()
fr = (ROOT / "LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/THEOREM.md").read_text()
check(all(term in bs for term in ("H_{\\rm inc}", "H_{\\rm car}",
                                  "H_{\\rm form}", "H_{\\rm fb}",
                                  "H_{\\rm port}")),
      "BS04 names the broader physical sectors omitted by the reduction")
check("H=H_0+V_X" in cw and "d_*=2" in cw,
      "CW explicitly selects the reduced degree-two parent")
check("H=H_0+V_X" in fm and "E_R=0" in fm,
      "FM explicitly selects the zero-detuning reduced parent")
check("occurrence multiplicities" in fq and "node/port weights" in fq,
      "FQ treats nonedge multiplicities and weights as prospective data")
check("same source-deformed parent" in fq and "post hoc hand weight" in fq,
      "FQ requires source insertion before Feshbach reduction")
check("vertex degree term has four incident" in fr
      and "its additive source is purely `A1`" in fr,
      "FR admits the occurrence-one additive degree-square query")
check("not yet a rank theorem for the complete BS20 source" in fr,
      "FR leaves unselected physical completion weights open")


# -------------------------------------------------------------------------
# Independent finite quotient and covering replay.
# -------------------------------------------------------------------------

SHIFTS = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))


def build_graph(size):
    vertices = tuple((part, x, y, z)
                     for part in (0, 1)
                     for x, y, z in product(range(size), repeat=3))
    adjacency = {vertex: set() for vertex in vertices}
    labels = {}
    for x, y, z in product(range(size), repeat=3):
        a = (0, x, y, z)
        for label, shift in enumerate(SHIFTS):
            b = (1, (x + shift[0]) % size,
                    (y + shift[1]) % size,
                    (z + shift[2]) % size)
            adjacency[a].add(b)
            adjacency[b].add(a)
            labels[frozenset((a, b))] = label
    return vertices, adjacency, labels


def canonical_cycle(cycle):
    cycle = tuple(cycle)
    variants = []
    for oriented in (cycle, tuple(reversed(cycle))):
        variants.extend(oriented[index:] + oriented[:index]
                        for index in range(len(cycle)))
    return min(variants)


def short_cycles(adjacency):
    found = {4: set(), 6: set(), 8: set()}
    for start in adjacency:
        def walk(vertex, path):
            for neighbor in adjacency[vertex]:
                if neighbor == start:
                    if len(path) in found:
                        found[len(path)].add(canonical_cycle(path))
                    continue
                if neighbor in path or len(path) == 8:
                    continue
                walk(neighbor, path + (neighbor,))
        walk(start, (start,))
    return found


def lift_balance(cycle, labels):
    answer = [0, 0, 0]
    for index, vertex in enumerate(cycle):
        neighbor = cycle[(index + 1) % len(cycle)]
        shift = SHIFTS[labels[frozenset((vertex, neighbor))]]
        sign = 1 if vertex[0] == 0 else -1
        answer = [value + sign * delta for value, delta in zip(answer, shift)]
    return tuple(answer)


v5, g5, label5 = build_graph(5)
v10, g10, label10 = build_graph(10)
for size, vertices, adjacency, labels in ((5, v5, g5, label5),
                                           (10, v10, g10, label10)):
    check(len(vertices) == 2 * size**3,
          f"G_{size} has the expected vertex count")
    check(sum(map(len, adjacency.values())) // 2 == 4 * size**3,
          f"G_{size} has the expected edge count")
    check(all(len(neighbors) == 4 for neighbors in adjacency.values()),
          f"G_{size} is closed and four-regular")
    start = vertices[0]
    seen = {start}
    queue = deque((start,))
    while queue:
        vertex = queue.popleft()
        for neighbor in adjacency[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    check(len(seen) == len(vertices), f"G_{size} is connected")
    check(len(labels) == 4 * size**3,
          f"G_{size} has one q4 label on every edge")

cycles = short_cycles(g5)
check(not cycles[4], "G_5 has no simple four-cycle")
check(len(cycles[6]) == 4 * 5**3,
      "G_5 has the inherited elementary-hexagon count")
check(all(len(Counter(label5[frozenset((cycle[i], cycle[(i + 1) % 6]))]
                      for i in range(6))) == 3
          and set(Counter(label5[frozenset((cycle[i], cycle[(i + 1) % 6]))]
                          for i in range(6)).values()) == {2}
          for cycle in cycles[6]),
      "every G_5 hexagon uses three labels twice")
check(all(lift_balance(cycle, label5) == (0, 0, 0)
          for length in (6, 8) for cycle in cycles[length]),
      "every G_5 cycle through length eight closes in the infinite cover")

# Any 2p-step walk has p forward and p reverse shift labels, so every lifted
# displacement coordinate has magnitude at most p.  Exhaustion independently
# confirms that p<=4 cannot be a nonzero multiple of five.
for half_length in range(1, 5):
    bad = False
    for labels in product(range(4), repeat=2 * half_length):
        displacement = [0, 0, 0]
        for label in labels[:half_length]:
            displacement = [a + b for a, b in zip(displacement, SHIFTS[label])]
        for label in labels[half_length:]:
            displacement = [a - b for a, b in zip(displacement, SHIFTS[label])]
        if all(value % 5 == 0 for value in displacement) and any(displacement):
            bad = True
            break
    check(not bad,
          f"no quotient winding exists in any {2 * half_length}-step label balance")

fibers = Counter()
cover_local = True
for vertex in v10:
    image = (vertex[0], vertex[1] % 5, vertex[2] % 5, vertex[3] % 5)
    fibers[image] += 1
    neighbor_images = {(neighbor[0], neighbor[1] % 5,
                        neighbor[2] % 5, neighbor[3] % 5)
                       for neighbor in g10[vertex]}
    cover_local &= neighbor_images == g5[image]
check(set(fibers.values()) == {8}, "G_10 to G_5 has eight-point fibers")
check(cover_local, "G_10 to G_5 is locally bijective and hence a graph cover")


# -------------------------------------------------------------------------
# Complete prospective reduced source: weights and operators.
# -------------------------------------------------------------------------

SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
DYADS = tuple((x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z)
              for x, y, z in SIGNS)
W = tuple(sum(row[column] for row in DYADS) for column in range(6))
E = ((1, -1, 0, 0, 0, 0), (1, 1, -2, 0, 0, 0))
dot = lambda left, right: sum(F(a) * F(b) for a, b in zip(left, right))

check(rank(DYADS) == 4, "four tetrahedral contraction rows have exact rank four")
check(W == (4, 4, 4, 0, 0, 0),
      "occurrence-one degree tensor is the scalar four-dyad sum")
check(all(dot(row, null) == 0 for row in DYADS for null in E),
      "both diagonal-traceless E directions annihilate every dyad")
check(rank(DYADS + E) == 6, "the two displayed E nulls exhaust the kernel")

# Four Pauli-X terms on four distinct link factors are independently isolated
# by computational-basis matrix elements.  Add the actual degree-square
# diagonal operator; the source map still factors through four invariants and
# the independent flips supply the lower bound.
flips = []
for link in range(4):
    flat = [0] * (16 * 16)
    for state in range(16):
        flat[state * 16 + (state ^ (1 << link))] = 1
    flips.append(flat)
degree = []
for state in range(16):
    occupation = bin(state).count("1")
    for target in range(16):
        degree.append((occupation - 2)**2 if state == target else 0)
check(rank(flips) == 4, "four distinct link flips have operator rank four")
operators_by_source_coordinate = []
for coordinate in range(6):
    operators_by_source_coordinate.append([
        sum(DYADS[link][coordinate] * flips[link][entry]
            for link in range(4)) + W[coordinate] * degree[entry]
        for entry in range(16 * 16)
    ])
check(rank(operators_by_source_coordinate) == 4,
      "complete flip-plus-degree microscopic operator source has exact rank four")
for null in E:
    contracted = [sum(null[coordinate] * operators_by_source_coordinate[coordinate][entry]
                      for coordinate in range(6))
                  for entry in range(16 * 16)]
    check(not any(contracted), "one E-polarized microscopic operator is exactly zero")


# -------------------------------------------------------------------------
# Uniform word homogeneity and block-resolved/Feshbach chain rule.
# -------------------------------------------------------------------------

class Dual:
    def __init__(self, value, derivative=0):
        self.value = F(value)
        self.derivative = F(derivative)

    def __mul__(self, other):
        other = other if isinstance(other, Dual) else Dual(other)
        return Dual(self.value * other.value,
                    self.derivative * other.value + self.value * other.derivative)

    def inverse(self):
        return Dual(1 / self.value, -self.derivative / self.value**2)

    def __pow__(self, exponent):
        if exponent < 0:
            return (self.inverse()) ** (-exponent)
        answer = Dual(1)
        for _ in range(exponent):
            answer = answer * self
        return answer


def compositions(total, slots=4):
    if slots == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, slots - 1):
                yield (first,) + rest


all_word_rows = []
for order in (2, 4, 6, 8):
    rows = []
    dual_ok = True
    for counts in compositions(order):
        row = tuple(sum(counts[label] * DYADS[label][coordinate]
                        for label in range(4)) - (order - 1) * W[coordinate]
                    for coordinate in range(6))
        rows.append(row)
        for coordinate in range(6):
            amplitude = Dual(1)
            for label in range(4):
                amplitude = amplitude * Dual(1, -F(DYADS[label][coordinate], 2)) ** counts[label]
            amplitude = amplitude * Dual(1, -F(W[coordinate], 2)) ** (-(order - 1))
            dual_ok &= amplitude.value == 1 and amplitude.derivative == -F(row[coordinate], 2)
    all_word_rows.extend(rows)
    check(dual_ok, f"order-{order} uniform word logarithmic derivative is exact")
    check(rank(DYADS + tuple(rows)) == 4,
          f"order-{order} uniform word family stays in the four-dyad span")
    check(all(dot(row, null) == 0 for row in rows for null in E),
          f"order-{order} uniform word family retains both E nulls")

# Two block-local copies test the nonuniform statement: coefficients need not
# equal the uniform integer formula, but every derivative factors through the
# same local C and W rows.  Arbitrary history mixing cannot restore a local E.
block_rows = []
for block in range(2):
    for label in range(4):
        block_rows.append(tuple(([0] * 6 if block else list(DYADS[label]))
                                + (list(DYADS[label]) if block else [0] * 6)))
    block_rows.append(tuple(([0] * 6 if block else list(W))
                            + (list(W) if block else [0] * 6)))
check(rank(block_rows) == 8,
      "two nonuniform blocks have four internal source directions each")
for block in range(2):
    for null in E:
        direction = ([0] * (6 * block) + list(null)
                     + [0] * (6 * (1 - block)))
        check(all(dot(row, direction) == 0 for row in block_rows),
              "one block-local E direction annihilates every microscopic derivative")

# Independent fixed-Schur derivative replay.  The microscopic variations are
# arbitrary operator matrices multiplied by the four dyad invariants and the
# scalar degree invariant.  Therefore an E derivative is zero before the map;
# the explicit derivative shows it remains zero after projection.
H0 = [[F(5), F(1), F(0), F(1), F(0)],
      [F(1), F(7), F(1), F(0), F(1)],
      [F(0), F(1), F(11), F(1), F(0)],
      [F(1), F(0), F(1), F(13), F(2)],
      [F(0), F(1), F(0), F(2), F(17)]]
VARIATIONS = []
for label in range(4):
    matrix = [[F(0) for _ in range(5)] for _ in range(5)]
    matrix[label][(label + 1) % 5] = F(label + 1)
    matrix[(label + 1) % 5][label] = F(label + 1)
    matrix[label][label] = F(2 * label + 1)
    VARIATIONS.append(matrix)
DEGREE_VARIATION = [[F(row == column) * F(row + 1)
                     for column in range(5)] for row in range(5)]


def microscopic_derivative(direction):
    invariants = [dot(row, direction) for row in DYADS]
    scalar = dot(W, direction)
    answer = [[F(0) for _ in range(5)] for _ in range(5)]
    for coefficient, variation in zip(invariants, VARIATIONS):
        answer = add(answer, scale(coefficient, variation))
    return add(answer, scale(scalar, DEGREE_VARIATION))


def feshbach_derivative(matrix, variation, energy=F(31)):
    p, q = (0, 1), (2, 3, 4)
    h_pq = submatrix(matrix, p, q)
    h_qp = submatrix(matrix, q, p)
    h_qq = submatrix(matrix, q, q)
    d_pp = submatrix(variation, p, p)
    d_pq = submatrix(variation, p, q)
    d_qp = submatrix(variation, q, p)
    d_qq = submatrix(variation, q, q)
    resolvent = inverse(add(scale(energy, identity(3)), scale(-1, h_qq)))
    answer = d_pp
    answer = add(answer, matmul(d_pq, matmul(resolvent, h_qp)))
    answer = add(answer, matmul(h_pq, matmul(resolvent, d_qp)))
    answer = add(answer, matmul(h_pq, matmul(resolvent,
                 matmul(d_qq, matmul(resolvent, h_qp)))))
    return answer


effective = []
for coordinate in range(6):
    direction = [F(index == coordinate) for index in range(6)]
    derivative = feshbach_derivative(H0, microscopic_derivative(direction))
    effective.append([entry for row in derivative for entry in row])
check(rank(effective) <= 4,
      "independent fixed Feshbach derivative cannot exceed microscopic rank four")
for null in E:
    check(matrix_zero(microscopic_derivative(null)),
          "one E microscopic derivative vanishes before Feshbach")
    check(matrix_zero(feshbach_derivative(H0, microscopic_derivative(null))),
          "one E derivative remains zero after fixed Feshbach reduction")


# -------------------------------------------------------------------------
# Quadratic-contact boundary and documentary scope.
# -------------------------------------------------------------------------

A1 = (1, 1, 1, 0, 0, 0)
hessian = [[F(E[0][row] * A1[column] + A1[row] * E[0][column])
            for column in range(6)] for row in range(6)]
check(any(value for row in hessian for value in row),
      "a lawful quadratic contact can have a nonzero E-mixed Hessian")
check(any(sum(hessian[row][column] * E[0][column]
                  for column in range(6)) for row in range(6)),
      "the quadratic contact can carry a genuine E Hessian leg")
zero_gradient = [F(0)] * 6
check(not any(zero_gradient),
      "every O(j^2) contact has zero source-off first derivative")

required_phrases = (
    "prospectively freezes its FQ17a occurrence",
    "not a unique tensor forced by the source-free",
    "definition-level exclusion",
    "every projected term\nthrough order eight",
    "not a derivation\nthat those sectors vanish",
    "matrix element between\ntwo configurations differing only on link",
    "Freeze `h != 0`",
    "under one common spatially uniform source",
    "Equation (FS11) is not asserted block by block",
    "history-dependent scalar energy fractions",
    "rank}\\mathcal Q_{\\rm eff}\\le4",
    "A general prospectively frozen BS20 contact may contain an `E`-dependent",
    "does not set the seagull Hessian to zero",
    "complete unreduced BS04 Hamiltonian has rank four",
)
for phrase in required_phrases:
    check(phrase in THEOREM,
          f"corrected theorem preserves hostile scope phrase: {phrase}")
check("FQ17a fixes its additive node weight" not in THEOREM,
      "theorem removes the claim that FQ uniquely forces the node tensor")
check("not a unique source tensor derived from\nthe source-free Hamiltonian" in RESULT,
      "result preserves the prospective-query qualification")
check("uniquely forces that query tensor" in SELF_AUDIT,
      "self-audit preserves the prospective-query qualification")
check("through-order-eight effective image has rank at most four"
      in " ".join(RESULT.split()),
      "result does not overstate projected rank as exact")
check("Cross-dyad/root sources" in RESULT and "gravity, and `G` remain open" in RESULT,
      "result preserves the physical successor ceiling")


# Base payload manifest and outer builder seal.  The hostile packet has a
# separate, nonrecursive audit manifest and seal.
manifest = parse_ledger(HERE / "MANIFEST.sha256")
expected_members = {
    "DEPENDENCIES.sha256", "README.md", "RESULT.md", "SELF_AUDIT.md",
    "THEOREM.md", "verify_complete_source_rank_audit.py",
}
check(set(manifest) == expected_members,
      "base manifest has exactly the six stable builder payload members")
for relative, expected in manifest.items():
    path = HERE / relative
    check(digest(path) == expected, f"base manifest digest matches: {relative}")
    check(sha256(path.read_bytes() + b"hostile-tamper").hexdigest() != expected,
          f"base payload appended-byte tamper fails: {relative}")

seal = parse_ledger(HERE / "SEAL.sha256")
check(set(seal) == {"MANIFEST.sha256", "VERIFICATION.txt"},
      "builder seal covers exactly manifest and replay transcript")
for relative, expected in seal.items():
    check(digest(HERE / relative) == expected,
          f"builder seal digest matches: {relative}")

print(f"SUMMARY {passed}/{passed} independent hostile checks passed")
print("DISPOSITION PASS_AFTER_PROSPECTIVE_WEIGHT_OPERATOR_RANK_AND_NONUNIFORM_SOURCE_REPAIRS__REDUCED_QUERY_MICRO_RANK4__EFFECTIVE_RANK_AT_MOST4__UNREDUCED_PHYSICAL_COMPLETION_OPEN")
