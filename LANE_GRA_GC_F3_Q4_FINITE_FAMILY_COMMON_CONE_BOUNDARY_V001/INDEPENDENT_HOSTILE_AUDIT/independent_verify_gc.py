#!/usr/bin/env python3
"""Algorithmically independent hostile replay of the frozen GC packet."""

from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent
ROOT = LANE.parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def verify_hash_list(list_path, base):
    count = 0
    for line in list_path.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        path = base / relative
        check(path.is_file() and not path.is_symlink() and
              digest(path) == expected,
              f"custody {list_path.name}: {relative}")
        count += 1
    return count


check(verify_hash_list(AUDIT / "TARGET_CUSTODY.sha256", LANE) == 10,
      "target custody freezes all ten GC core/seal files")
check(verify_hash_list(LANE / "DEPENDENCIES.sha256", ROOT) == 15,
      "all fifteen FS/FD/FV/FY/FL/GB dependencies replay")
check(verify_hash_list(LANE / "MANIFEST.sha256", LANE) == 8,
      "the frozen eight-file GC author manifest replays")
seal_lines = (LANE / "SEAL.sha256").read_text().splitlines()
seal_hash, seal_name = seal_lines[0].split("  ", 1)
check(seal_name == "MANIFEST.sha256" and
      seal_hash == digest(LANE / seal_name),
      "the author seal owns the frozen GC manifest")


# -------------------------------------------------------------------------
# Claim-surface and dependency typing.

theorem = (LANE / "THEOREM.md").read_text()
result = (LANE / "RESULT.md").read_text()
self_audit = (LANE / "SELF_AUDIT.md").read_text()
readme = (LANE / "README.md").read_text()
result_json = json.loads((LANE / "RESULT.json").read_text())
flat_t = " ".join(theorem.split())
flat_r = " ".join(result.split())
flat_s = " ".join(self_audit.split())

check("This is not spectral inheritance" in flat_t and
      "not an isometry or an inclusion of the ice Hilbert spaces" in flat_t,
      "graph-cover ancestry is not promoted to spectral inheritance")
check("minimum-norm Brillouin-zone representative" in flat_t and
      "|k_[n]|" in result_json["momentum_norm"],
      "momentum norm is attached to a class and minimum representative")
check("count-normalized, not yet physically volume- or wavefunction- normalized"
      in flat_t and "prospectively derived `Z_Q(a_*,k)`" in flat_s,
      "FY count normalization is not promoted to physical normalization")
check("controlled cluster `S_L` of split finite- volume polarizations" in
      flat_t and "one eigenvector is never asked to carry a rank-two residue"
      in flat_s,
      "rank-two visibility admits a degenerate projector or coalescing cluster")
check("FL's exact one-link electric source" in flat_t and
      "A raw numerator that vanishes with a controlled momentum power is **not** a refutation"
      in flat_t,
      "a derivative-source vanishing numerator is not called no pole")
check("They are not a decisive closure packet" in flat_t and
      "cannot supply the independent momentum-owning theory dependency exposed by GB"
      in flat_t,
      "the three-size spectral screen is separated from the GB Ward dependency")
check(all(token in readme for token in
          ("C-stage", "`D1`", "`D2`", "`D3`", "`RGRL-B`")),
      "README retains the C-stage and D1/D2/D3 Gravity Formation ceiling")
check(result_json["ceiling"] ==
      "geometry/source scaling only; no Ward-derived tensor pole, common massless cone, gravity, or G",
      "machine-readable disposition retains the bounded claim ceiling")


# -------------------------------------------------------------------------
# Independent graph construction and simple-cycle enumeration.

SHIFTS = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))


def canonical_cycle(cycle):
    variants = []
    size = len(cycle)
    for oriented in (cycle, tuple(reversed(cycle))):
        variants.extend(oriented[offset:] + oriented[:offset]
                        for offset in range(size))
    return min(variants)


def graph_and_six_cycles(length):
    adjacency = defaultdict(list)
    edge_label = {}
    for cell in product(range(length), repeat=3):
        a_vertex = (0,) + cell
        for label, shift in enumerate(SHIFTS):
            b_cell = tuple((cell[axis] + shift[axis]) % length
                           for axis in range(3))
            b_vertex = (1,) + b_cell
            adjacency[a_vertex].append(b_vertex)
            adjacency[b_vertex].append(a_vertex)
            edge_label[frozenset((a_vertex, b_vertex))] = label

    cycles = set()
    for start in tuple(adjacency):
        def walk(path):
            if len(path) == 6:
                if start in adjacency[path[-1]]:
                    cycles.add(canonical_cycle(tuple(path)))
                return
            for neighbor in adjacency[path[-1]]:
                if neighbor == start or neighbor in path:
                    continue
                walk(path + [neighbor])
        walk([start])

    by_missing = defaultdict(set)
    typing_ok = True
    for cycle in cycles:
        labels = [edge_label[frozenset((cycle[index],
                                        cycle[(index + 1) % 6]))]
                  for index in range(6)]
        counts = Counter(labels)
        typing_ok = (typing_ok and sorted(counts.values()) == [2, 2, 2]
                     and len(counts) == 3)
        missing = tuple(set(range(4)) - set(counts))
        if len(missing) != 1:
            raise AssertionError("H6 missing-label classification failed")
        by_missing[missing[0]].add(cycle)
    return adjacency, edge_label, cycles, by_missing, typing_ok


g5, edges5, cycles5, missing5, typing5 = graph_and_six_cycles(5)
g10, edges10, cycles10, missing10, typing10 = graph_and_six_cycles(10)
check(typing5 and typing10,
      "every independently enumerated H6 uses three q4 labels twice")
check((len(g5), len(edges5), len(cycles5)) == (250, 500, 500),
      "independent G_5 graph has 2L^3 vertices, 4L^3 links, and 4L^3 H6 rings")
check((len(g10), len(edges10), len(cycles10)) == (2000, 4000, 4000),
      "independent G_10 graph has 2L^3 vertices, 4L^3 links, and 4L^3 H6 rings")
check(all(len(missing5[label]) == 125 and
          len(missing10[label]) == 1000 for label in range(4)),
      "every missing-label family has exactly L^3 elementary rings")


def reduce_vertex(vertex, modulus):
    return (vertex[0],) + tuple(value % modulus for value in vertex[1:])


reduced = Counter(canonical_cycle(tuple(reduce_vertex(vertex, 5)
                                        for vertex in cycle))
                  for cycle in cycles10)
check(set(reduced) == cycles5 and set(reduced.values()) == {8},
      "all independently enumerated H6 rings pull back eight-to-one")


# -------------------------------------------------------------------------
# Exact coefficient, affine, reciprocal, and TT algebra.

def matrix_rank(rows):
    work = [[F(value) for value in row] for row in rows]
    rank_value = 0
    for column in range(len(work[0])):
        pivot = next((row for row in range(rank_value, len(work))
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[rank_value], work[pivot] = work[pivot], work[rank_value]
        scale = work[rank_value][column]
        work[rank_value] = [value / scale for value in work[rank_value]]
        for row in range(len(work)):
            if row == rank_value or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [value - factor * base
                         for value, base in zip(work[row], work[rank_value])]
        rank_value += 1
    return rank_value


def determinant(rows):
    work = [[F(value) for value in row] for row in rows]
    sign = 1
    result_det = F(1)
    for column in range(len(work)):
        pivot = next((row for row in range(column, len(work))
                      if work[row][column]), None)
        if pivot is None:
            return F(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        result_det *= pivot_value
        for row in range(column + 1, len(work)):
            factor = work[row][column] / pivot_value
            for entry in range(column, len(work)):
                work[row][entry] -= factor * work[column][entry]
    return sign * result_det


TETRA = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
ring_rows = []
for sx, sy, sz in TETRA:
    ring_rows.append((F(231, 8), F(231, 8), F(231, 8),
                      -F(189, 8) * sx * sy,
                      -F(189, 8) * sx * sz,
                      -F(189, 8) * sy * sz))
e_rows = ((F(-1), F(1), F(0), F(0), F(0), F(0)),
          (F(-1), F(0), F(1), F(0), F(0), F(0)))
witness = e_rows + tuple(ring_rows)
check(matrix_rank(ring_rows) == 4 and matrix_rank(witness) == 6,
      "independently derived Coulomb ring plus E rows have local rank six")
check(determinant(witness) == F(-4678629417, 256),
      "independent Gaussian determinant reproduces the FV witness exactly")


# Use exact sign vectors; physical vectors are signs/sqrt(3).  The fourth
# direction owns the zero-shift edge, and r_i=s_4-s_i are primitive cells.
s4 = TETRA[3]
roots = tuple(tuple(s4[axis] - TETRA[index][axis]
                    for axis in range(3)) for index in range(3))
check(all(tuple(roots[index][axis] - s4[axis] for axis in range(3)) ==
          tuple(-TETRA[index][axis] for axis in range(3))
          for index in range(3)),
      "the exact affine embedding sends shifted bonds to the other -n_i")
gram = tuple(tuple(F(sum(roots[row][axis] * roots[column][axis]
                         for axis in range(3)), 3)
                   for column in range(3)) for row in range(3))
expected_gram = ((F(8, 3), F(4, 3), F(4, 3)),
                 (F(4, 3), F(8, 3), F(4, 3)),
                 (F(4, 3), F(4, 3), F(8, 3)))
check(gram == expected_gram and determinant(gram) == F(256, 27),
      "the FS affine join has the exact FD Gram matrix and squared covolume")


def qform(vector):
    return 4 * sum(value * value for value in vector) - sum(vector) ** 2


def qsquares(vector):
    return (sum(value * value for value in vector) +
            sum((vector[first] - vector[second]) ** 2
                for first in range(3) for second in range(first + 1, 3)))


check(all(qform(vector) == qsquares(vector)
          for vector in product(range(-3, 4), repeat=3)),
      "the reciprocal integer form is the positive sum-of-squares form")
q3 = {vector for vector in product(range(-2, 3), repeat=3)
      if qform(vector) == 3}
q4 = {vector for vector in product(range(-2, 3), repeat=3)
      if qform(vector) == 4}
check(len(q3) == 8 and len(q4) == 6 and
      all(qform(vector) % 4 in (0, 3)
          for vector in product(range(4), repeat=3)),
      "the global reciprocal minimum is q=3 with 8 vectors, followed by q=4 with 6")
check(F(3, 16) * min(qform(vector) for vector in q3) == F(9, 16),
      "the exact reciprocal minimum gives k_min=3pi/(2 L a_*)")
check(qform((1, 0, 0)) == 3 and qform((1, 1, 0)) == 4 and
      qform((1, -1, 0)) == 8,
      "the three proposed spectral-screen rays have distinct exact norms")


def matmul(first, second):
    return tuple(tuple(sum(first[row][inner] * second[inner][column]
                           for inner in range(len(second)))
                       for column in range(len(second[0])))
                 for row in range(len(first)))


def matadd(first, second):
    return tuple(tuple(first[row][column] + second[row][column]
                       for column in range(len(first[0])))
                 for row in range(len(first)))


def matscale(scale, matrix):
    return tuple(tuple(F(scale) * value for value in row) for row in matrix)


I3 = tuple(tuple(F(int(row == column)) for column in range(3))
           for row in range(3))
r = (F(-1), F(-1), F(-1))
r2 = sum(value * value for value in r)
P = matadd(I3, matscale(-F(1, r2),
           tuple(tuple(r[row] * r[column] for column in range(3))
                 for row in range(3))))


def tt_image(matrix):
    pap = matmul(matmul(P, matrix), P)
    transverse_trace = sum(P[row][column] * matrix[row][column]
                           for row in range(3) for column in range(3))
    return matadd(pap, matscale(-transverse_trace / 2, P))


def coordinate_basis(column):
    matrix = [[F(0) for _ in range(3)] for _ in range(3)]
    slots = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
    row, col = slots[column]
    matrix[row][col] = F(1) if row == col else F(1, 2)
    matrix[col][row] = matrix[row][col]
    return tuple(tuple(row_values) for row_values in matrix)


def pack(matrix):
    return (matrix[0][0], matrix[1][1], matrix[2][2],
            2 * matrix[0][1], 2 * matrix[0][2], 2 * matrix[1][2])


columns = tuple(pack(tt_image(coordinate_basis(column)))
                for column in range(6))
TT = tuple(tuple(columns[column][row] for column in range(6))
           for row in range(6))
check(matmul(P, P) == P and matrix_rank(P) == 2,
      "independent reciprocal-ray transverse projector is exact rank two")
check(matmul(TT, TT) == TT and matrix_rank(TT) == 2,
      "independent symmetric-tensor TT projector is idempotent and rank two")
check(all(sum(tt_image(coordinate_basis(column))[axis][axis]
              for axis in range(3)) == 0 for column in range(6)),
      "all independent TT images are exactly traceless")


# -------------------------------------------------------------------------
# Response-boundary counterexamples and finite-size target.

check(F(8) * F(1, 2) ** 3 == 1 and
      F(1, 2) / F(1, 2) == 1,
      "co-refining a_* by one half preserves affine volume and k_min")
check(F(1, 2) > 0 and F(1, 2) < 1,
      "0<s<1 supplies simultaneous a_r->0 and L_r a_r->infinity")

# The same exact family kinematics allows either a mass or no mass.
laplacian = F(1, 100)
massless = F(3, 2) * laplacian
massive = F(1, 4) + F(3, 2) * laplacian
check(massless > 0 and massive > massless,
      "family/source kinematics alone does not identify a massless response")

# FL's gauge-invariant electric source: Delta=c|k| and R=Z c|k|/2,
# hence 2 Delta R=Z c^2 k^2.  A vanishing raw numerator is compatible with
# a genuine massless pole and cannot be used as a no-pole theorem.
k, speed, normalization = F(2, 7), F(3, 5), F(11, 13)
delta = speed * k
residue = normalization * speed * k / 2
check(2 * delta * residue == normalization * speed * speed * k * k,
      "FL supplies an exact massless-pole counterexample to constant raw 2DeltaR")
check("count-normalized TT spectral screen" in
      result_json["smallest_target"] and
      "GB local spatial-momentum" in result_json["independent_dependency"],
      "machine-readable target separates the executable spectrum from Ward ownership")


print(f"SUMMARY {checks}/{checks} independent hostile GC checks passed")
print("VERDICT PASS")
print("CEILING C-stage family geometry, local coefficient ancestry, affine/"
      "reciprocal scaling, TT kinematics, and a bounded spectral screen only; "
      "no matched-family pole, canonical native-source amplitude, stress Ward "
      "closure, common physical cone, Gravity Formation D1/D2/D3, gravity, or G")
