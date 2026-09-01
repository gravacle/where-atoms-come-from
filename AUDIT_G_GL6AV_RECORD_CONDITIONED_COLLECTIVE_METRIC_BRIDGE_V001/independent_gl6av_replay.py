#!/usr/bin/env python3
"""Independent exact replay for GL6AV.

This file deliberately imports no GL6AV or upstream verifier.
"""

from fractions import Fraction as F
from itertools import permutations, product


checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def dot(u, v):
    return sum((x * y for x, y in zip(u, v)), F(0))


def rank(rows):
    matrix = [[F(x) for x in row] for row in rows]
    if not matrix:
        return 0
    row = 0
    for col in range(len(matrix[0])):
        pivot = next((i for i in range(row, len(matrix)) if matrix[i][col]), None)
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        scale = matrix[row][col]
        matrix[row] = [x / scale for x in matrix[row]]
        for i in range(len(matrix)):
            if i != row and matrix[i][col]:
                scale = matrix[i][col]
                matrix[i] = [x - scale * y
                             for x, y in zip(matrix[i], matrix[row])]
        row += 1
    return row


def determinant(matrix):
    n = len(matrix)
    total = F(0)
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = F(-1 if inversions % 2 else 1)
        for i in range(n):
            term *= matrix[i][p[i]]
        total += term
    return total


def matmul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return [list(row) for row in zip(*a)]


def cycle_type(p):
    seen = set()
    cycles = []
    for i in range(len(p)):
        if i in seen:
            continue
        j = i
        length = 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = p[j]
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))


# Binary six-support law and its unique product extension.
for bits in product((0, 1), repeat=6):
    coefficient = 1
    for bit in bits:
        coefficient *= bit
    check(coefficient == int(all(bits)), f"binary six-factor AND {bits}")
for omitted in range(6):
    bits = [1] * 6
    bits[omitted] = 0
    coefficient = 1
    for bit in bits:
        coefficient *= bit
    check(coefficient == 0, f"one omitted support kills loop {omitted}")

# The product is affine in each coordinate, hence multilinear.
sample = [F(2, 3), F(-3, 5), F(7, 11), F(5, 4), F(-2, 9), F(13, 8)]
for coordinate in range(6):
    def value(x):
        entries = sample[:]
        entries[coordinate] = x
        answer = F(1)
        for entry in entries:
            answer *= entry
        return answer
    x0, x1, lam = F(-2), F(3), F(4, 7)
    check(value((1 - lam) * x0 + lam * x1)
          == (1 - lam) * value(x0) + lam * value(x1),
          f"multilinearity coordinate {coordinate}")

# Homogeneous clock law and sample spectral pushforward.
for q in (F(-3), F(-1), F(0), F(1, 2), F(1), F(3, 2), F(2)):
    coefficient = F(1)
    for _ in range(6):
        coefficient *= q
    check(coefficient == q**6, f"homogeneous q sixth power {q}")

energies = [F(-3), F(-1), F(0), F(2), F(5)]
q = F(3, 2)
scaled = [q**6 * energy for energy in energies]
base_gaps = [energy - energies[0] for energy in energies]
scaled_gaps = [energy - scaled[0] for energy in scaled]
check(scaled_gaps == [q**6 * gap for gap in base_gaps], "finite spectral scaling")
window = (F(-2), F(3))
base_window = (window[0] / q**6, window[1] / q**6)
for energy, energy_q in zip(energies, scaled):
    check((window[0] <= energy_q <= window[1])
          == (base_window[0] <= energy <= base_window[1]),
          f"spectral projector rescaling {energy}")
check(q**-6 * q**6 == 1, "Fourier Jacobian q^-6")

# Piecewise-constant prescribed q(t): time parameters add and depend on q^6.
segments = [(F(2), F(1, 3)), (F(3, 2), F(2, 5)), (F(1, 4), F(7, 3))]
pieces = [duration * amplitude**6 for duration, amplitude in segments]
sigma = sum(pieces, F(0))
check(sigma == sum(duration * amplitude**6 for duration, amplitude in segments),
      "prescribed time integral")
check(sigma == pieces[0] + (pieces[1] + pieces[2]), "time composition")

# Positive binary log chart has only the all-one authenticated point.
positive_binary = [bits for bits in product((0, 1), repeat=4) if all(x > 0 for x in bits)]
check(positive_binary == [(1, 1, 1, 1)], "positive binary chart is trivial")

# A=2(11^T-I): exact spectrum, determinant, and two-sided inverse.
A = [[F(0) if i == j else F(2) for j in range(4)] for i in range(4)]
one = [F(1)] * 4
check([dot(row, one) for row in A] == [F(6)] * 4, "A1 eigenvalue 6")
centered = ([F(1), F(-1), F(0), F(0)],
            [F(0), F(1), F(-1), F(0)],
            [F(0), F(0), F(1), F(-1)])
for i, vector in enumerate(centered):
    check([dot(row, vector) for row in A] == [F(-2) * x for x in vector],
          f"T2 eigenvalue -2 basis {i}")
check(rank(A) == 4, "orientation log map rank four")
check(determinant(A) == F(-48), "orientation log map determinant -48")
B = [[F(1, 6) - (F(1, 2) if i == j else F(0))
      for j in range(4)] for i in range(4)]
identity4 = [[F(int(i == j)) for j in range(4)] for i in range(4)]
check(matmul(A, B) == identity4, "orientation inverse right")
check(matmul(B, A) == identity4, "orientation inverse left")
rho = [F(2, 5), F(-1, 3), F(7, 11), F(5, 13)]
j = [dot(row, rho) for row in A]
recovered = [sum(j, F(0)) / 6 - entry / 2 for entry in j]
check(recovered == rho, "AV11 inverse factors")

# Tetrahedral Gram, resolution, evaluation rank, kernel, and AV16 inverse.
t = [
    [F(1, 2), F(1, 2), F(1, 2)],
    [F(1, 2), F(-1, 2), F(-1, 2)],
    [F(-1, 2), F(1, 2), F(-1, 2)],
    [F(-1, 2), F(-1, 2), F(1, 2)],
]
for a in range(4):
    for b in range(4):
        expected = F(3, 4) if a == b else F(-1, 4)
        check(dot(t[a], t[b]) == expected, f"tetrahedral Gram {a},{b}")
check([sum((t[a][i] for a in range(4)), F(0)) for i in range(3)] == [F(0)] * 3,
      "tetrahedral vectors centered")
resolution = [[sum((t[a][i] * t[a][j] for a in range(4)), F(0))
               for j in range(3)] for i in range(3)]
identity3 = [[F(int(i == j)) for j in range(3)] for i in range(3)]
check(resolution == identity3, "tetrahedral resolution identity")

# Coordinate order: xx, yy, zz, xy, xz, yz.
evaluation = []
for x, y, z in t:
    evaluation.append([x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z])
check(rank(evaluation) == 4, "tetrahedral evaluation rank four")
e_basis = ([F(1), F(-1), F(0), F(0), F(0), F(0)],
           [F(1), F(0), F(-1), F(0), F(0), F(0)])
for i, vector in enumerate(e_basis):
    check([dot(row, vector) for row in evaluation] == [F(0)] * 4,
          f"traceless diagonal kernel {i}")
check(6 - rank(evaluation) == len(e_basis), "kernel dimension exactly two")

orientation_data = [F(2, 7), F(-3, 5), F(11, 13), F(17, 19)]
trace = sum(orientation_data, F(0))
sxy = (orientation_data[0] - orientation_data[1]
       - orientation_data[2] + orientation_data[3]) / 2
sxz = (orientation_data[0] - orientation_data[1]
       + orientation_data[2] - orientation_data[3]) / 2
syz = (orientation_data[0] + orientation_data[1]
       - orientation_data[2] - orientation_data[3]) / 2
representative = [trace / 3, trace / 3, trace / 3, sxy, sxz, syz]
check([dot(row, representative) for row in evaluation] == orientation_data,
      "AV16 exact right inverse")
for basis_index in range(4):
    data = [F(int(i == basis_index)) for i in range(4)]
    trace = sum(data, F(0))
    coeffs = [
        trace / 3, trace / 3, trace / 3,
        (data[0] - data[1] - data[2] + data[3]) / 2,
        (data[0] - data[1] + data[2] - data[3]) / 2,
        (data[0] + data[1] - data[2] - data[3]) / 2,
    ]
    check([dot(row, coeffs) for row in evaluation] == data,
          f"AV16 basis right inverse {basis_index}")

# Equivariance under all port permutations.
sym_basis = [
    [[F(1), F(0), F(0)], [F(0), F(0), F(0)], [F(0), F(0), F(0)]],
    [[F(0), F(0), F(0)], [F(0), F(1), F(0)], [F(0), F(0), F(0)]],
    [[F(0), F(0), F(0)], [F(0), F(0), F(0)], [F(0), F(0), F(1)]],
    [[F(0), F(1), F(0)], [F(1), F(0), F(0)], [F(0), F(0), F(0)]],
    [[F(0), F(0), F(1)], [F(0), F(0), F(0)], [F(1), F(0), F(0)]],
    [[F(0), F(0), F(0)], [F(0), F(0), F(1)], [F(0), F(1), F(0)]],
]


def evaluate_tensor(S):
    return [dot(t[a], [dot(S[i], t[a]) for i in range(3)]) for a in range(4)]


for p in permutations(range(4)):
    O = [[sum((t[p[a]][i] * t[a][j] for a in range(4)), F(0))
          for j in range(3)] for i in range(3)]
    check(matmul(transpose(O), O) == identity3, f"port action orthogonal {p}")
    equivariant = True
    for S in sym_basis:
        transformed = matmul(matmul(O, S), transpose(O))
        old = evaluate_tensor(S)
        new = evaluate_tensor(transformed)
        equivariant &= all(new[p[a]] == old[a] for a in range(4))
    check(equivariant, f"evaluation equivariant {p}")

# S4 character check: R4=A1+T2 and Sym2(T2)=A1+E+T2.
class_sizes = {(1, 1, 1, 1): 1, (2, 1, 1): 6, (2, 2): 3,
               (3, 1): 8, (4,): 6}
expected_t2 = {(1, 1, 1, 1): 3, (2, 1, 1): 1, (2, 2): -1,
               (3, 1): 0, (4,): -1}
expected_e = {(1, 1, 1, 1): 2, (2, 1, 1): 0, (2, 2): 2,
              (3, 1): -1, (4,): 0}
observed_sizes = {kind: 0 for kind in class_sizes}
observed_t2 = {kind: set() for kind in class_sizes}
observed_e = {kind: set() for kind in class_sizes}
for p in permutations(range(4)):
    kind = cycle_type(p)
    observed_sizes[kind] += 1
    fixed = sum(p[i] == i for i in range(4))
    chi_t2 = fixed - 1
    p2 = tuple(p[p[i]] for i in range(4))
    chi_t2_squared_element = sum(p2[i] == i for i in range(4)) - 1
    chi_sym2 = (chi_t2**2 + chi_t2_squared_element) // 2
    observed_t2[kind].add(chi_t2)
    observed_e[kind].add(chi_sym2 - 1 - chi_t2)
for kind in class_sizes:
    check(observed_sizes[kind] == class_sizes[kind], f"S4 class size {kind}")
    check(observed_t2[kind] == {expected_t2[kind]}, f"T2 character {kind}")
    check(observed_e[kind] == {expected_e[kind]}, f"E character {kind}")
inner = sum(class_sizes[kind] * expected_t2[kind] * expected_e[kind]
            for kind in class_sizes)
check(inner == 0, "E and T2 character orthogonality")

# AV19: a speed q^6 v has squared temporal coefficient q^12 v^2.
v1 = F(7, 5)
for q in (F(1, 3), F(1), F(5, 4), F(2)):
    vq = q**6 * v1
    check(vq**2 == q**12 * v1**2, f"metric temporal q^12 factor {q}")

# Two different mode speeds produce two distinct cones, not one cone.
k = F(1)
speed_one = F(1)
speed_two = F(2)
null_frequencies = {speed_one * k, -speed_one * k,
                    speed_two * k, -speed_two * k}
check(len(null_frequencies) == 4, "two-speed characteristic has four roots")
check(speed_one != speed_two, "multiple-cone counterexample is nondegenerate")

print(f"PASS__GL6AV_INDEPENDENT_REPLAY__{checks}/{checks}")
