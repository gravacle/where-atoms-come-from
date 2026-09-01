#!/usr/bin/env python3
"""Independent exact replay of the frozen GL6U degree-interaction theorem.

This script uses only the Python standard library.  It does not import or
execute the GL6U author verifier.  The small-time series is rebuilt on the
complete 16-state active-star Hilbert space with symbolic polynomial
coefficients in (h, Delta, U_d, d_star).
"""

from __future__ import annotations

from fractions import Fraction as F
from math import comb, factorial


checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


# ---------------------------------------------------------------------------
# A tiny exact polynomial ring Q[h, Delta, U, p].
# ---------------------------------------------------------------------------

NVAR = 4
ZERO_MON = (0,) * NVAR


class Poly:
    def __init__(self, terms=None):
        raw = {} if terms is None else dict(terms)
        self.terms = {m: F(c) for m, c in raw.items() if c}

    @staticmethod
    def const(value):
        value = F(value)
        return Poly({ZERO_MON: value}) if value else Poly()

    @staticmethod
    def var(index):
        mon = [0] * NVAR
        mon[index] = 1
        return Poly({tuple(mon): F(1)})

    def __add__(self, other):
        other = as_poly(other)
        out = dict(self.terms)
        for mon, coefficient in other.terms.items():
            out[mon] = out.get(mon, F(0)) + coefficient
            if not out[mon]:
                del out[mon]
        return Poly(out)

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-as_poly(other))

    def __rsub__(self, other):
        return as_poly(other) - self

    def __mul__(self, other):
        other = as_poly(other)
        out = {}
        for left_mon, left_coefficient in self.terms.items():
            for right_mon, right_coefficient in other.terms.items():
                mon = tuple(a + b for a, b in zip(left_mon, right_mon))
                out[mon] = out.get(mon, F(0)) + left_coefficient * right_coefficient
        return Poly(out)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        check(exponent >= 0, "nonnegative polynomial exponent")
        out = Poly.const(1)
        base = self
        n = exponent
        while n:
            if n & 1:
                out = out * base
            base = base * base
            n //= 2
        return out

    def __truediv__(self, denominator):
        denominator = F(denominator)
        return Poly({m: c / denominator for m, c in self.terms.items()})

    def __eq__(self, other):
        return self.terms == as_poly(other).terms

    def evaluate(self, values):
        total = F(0)
        for mon, coefficient in self.terms.items():
            term = coefficient
            for value, exponent in zip(values, mon):
                term *= F(value) ** exponent
            total += term
        return total

    def is_zero(self):
        return not self.terms


def as_poly(value):
    return value if isinstance(value, Poly) else Poly.const(value)


h, Delta, U, p = (Poly.var(index) for index in range(NVAR))


def gap(r):
    """Exact active-star diagonal gap relative to the all-blank scalar."""
    return r * Delta + U * r * (r + 1 - 4 * p)


delta1 = gap(1)
delta2 = gap(2)
check(2 * delta1 - delta2 == -2 * U, "exact gap cancellation")
check(2 * delta2 - delta1 == 3 * Delta + U * (10 - 12 * p),
      "exact E-sector gap combination")


# ---------------------------------------------------------------------------
# Exact FPSS N=0 census and degree ownership.
# ---------------------------------------------------------------------------

N = 0
a0 = comb(N + 3, 3)
b0 = comb(N + 4, 3)
M = b0
guards = b0 - a0
raw_links = M * M
active_links = 4 * a0
blank_nonedges = raw_links - active_links
check((a0, b0, M, guards) == (1, 4, 4, 3), "N=0 site census")
check((raw_links, active_links, blank_nonedges) == (16, 4, 12),
      "N=0 link census")

for state in range(16):
    bits = [(state >> site) & 1 for site in range(4)]
    degree = sum(bits)
    # One active parent, four active children, and three guard parents.
    full_penalty = (as_poly(degree) - p) ** 2
    full_penalty += sum((as_poly(bit) - p) ** 2 for bit in bits)
    full_penalty += 3 * p**2
    blank_penalty = 8 * p**2
    check(full_penalty - blank_penalty == degree * (degree + 1 - 4 * p),
          f"degree ownership state {state}")


# ---------------------------------------------------------------------------
# Exact Dicke combinatorics, without diagonalizing or using the author rows.
# ---------------------------------------------------------------------------

for r in range(4):
    left_count = comb(4, r)
    right_count = comb(4, r + 1)
    total_flip_edges = left_count * (4 - r)
    check(total_flip_edges == right_count * (r + 1),
          f"Dicke flip-edge double count r={r}")
    # Square of <r+1|sum X|r>.
    squared_hop = F(total_flip_edges**2, left_count * right_count)
    check(squared_hop == (r + 1) * (4 - r),
          f"Dicke hopping square r={r}")

    # One-site X matrix element: count configurations of the other 3 sites.
    one_site_numerator = comb(3, r)
    squared_x = F(one_site_numerator**2, left_count * right_count)
    check(squared_x == F((r + 1) * (4 - r), 16),
          f"Dicke X-row square r={r}")

    # One-site X times Z on two other sites: exact signed transition count.
    signed = 0
    for other_state in range(8):
        other_bits = [(other_state >> q) & 1 for q in range(3)]
        if sum(other_bits) == r:
            signed += (-1) ** (other_bits[0] + other_bits[1])
    expected_sign = (1, -1, -1, 1)[r]
    expected_square = (F(1, 4), F(1, 24), F(1, 24), F(1, 4))[r]
    check((signed > 0) - (signed < 0) == expected_sign,
          f"Dicke XZZ sign r={r}")
    check(F(signed**2, left_count * right_count) == expected_square,
          f"Dicke XZZ square r={r}")


# ---------------------------------------------------------------------------
# Complete 16-state symbolic time series through order five.
# ---------------------------------------------------------------------------

DIM = 16


def apply_hamiltonian(vector):
    out = [Poly() for _ in range(DIM)]
    for state in range(DIM):
        r = bin(state).count("1")
        out[state] = out[state] + gap(r) * vector[state]
        for site in range(4):
            out[state] = out[state] - h * vector[state ^ (1 << site)]
    return out


v0 = [Poly() for _ in range(DIM)]
v0[0] = Poly.const(1)
powers = [v0]
for _ in range(5):
    powers.append(apply_hamiltonian(powers[-1]))


def gzero():
    return (Poly(), Poly())


def gadd(left, right):
    return (left[0] + right[0], left[1] + right[1])


def gscale(value, scalar):
    return (value[0] * scalar, value[1] * scalar)


def gconj(value):
    return (value[0], -value[1])


def gmul(left, right):
    return (left[0] * right[0] - left[1] * right[1],
            left[0] * right[1] + left[1] * right[0])


minus_i = ((F(1), F(0)), (F(0), F(-1)),
           (F(-1), F(0)), (F(0), F(1)))
coefficients = []
for order, vector in enumerate(powers):
    real_factor, imag_factor = minus_i[order % 4]
    coefficients.append([
        (entry * real_factor / factorial(order),
         entry * imag_factor / factorial(order))
        for entry in vector
    ])


def expectation_series(operator_action, max_order=5):
    result = []
    for total_order in range(max_order + 1):
        value = gzero()
        for left_order in range(total_order + 1):
            right_order = total_order - left_order
            for source in range(DIM):
                for target, weight in operator_action(source):
                    value = gadd(
                        value,
                        gmul(gconj(coefficients[left_order][target]),
                             gscale(coefficients[right_order][source], weight)),
                    )
        check(value[1].is_zero(), f"real expectation order {total_order}")
        result.append(value[0])
    return result


def action_x0(state):
    return [(state ^ 1, 1)]


def action_z0(state):
    return [(state, 1 if (state & 1) == 0 else -1)]


def action_x0z1z2(state):
    sign = -1 if (((state >> 1) & 1) + ((state >> 2) & 1)) % 2 else 1
    return [(state ^ 1, sign)]


x_series = expectation_series(action_x0)
z_series = expectation_series(action_z0)
y_series = expectation_series(action_x0z1z2)

for order in (1, 3, 5):
    check(x_series[order].is_zero(), f"x odd order {order}")
    check(z_series[order].is_zero(), f"z odd order {order}")
    check(y_series[order].is_zero(), f"y odd order {order}")

check(x_series[0].is_zero(), "x blank value")
check(z_series[0] == 1, "z blank value")
check(y_series[0].is_zero(), "y blank value")
check(x_series[2] == h * delta1, "x order-two coefficient")
check(z_series[2] == -2 * h**2, "z order-two coefficient")


def convolve(left, right, max_order=5):
    return [sum((left[q] * right[n - q] for q in range(n + 1)), Poly())
            for n in range(max_order + 1)]


x_minus_y = [x_series[n] - y_series[n] for n in range(6)]
z_squared = convolve(z_series, z_series)
xz_squared = convolve(x_series, z_squared)
factor_defect = [y_series[n] - xz_squared[n] for n in range(6)]

for order in range(4):
    check(x_minus_y[order].is_zero(), f"x-y lower order {order}")
    check(factor_defect[order].is_zero(), f"factor defect lower order {order}")
check(x_minus_y[4] == F(4, 3) * h**3 * (2 * delta2 - delta1),
      "exact E-sector order-four coefficient")
check(factor_defect[4] == -F(16, 3) * h**3 * U,
      "exact interaction factorization defect")
check((-4 * h * factor_defect[4]) == F(64, 3) * h**4 * U,
      "exact shared-response correction")
check(factor_defect[5].is_zero(), "factor defect order-five zero")

# Sector leading coefficients reconstructed from D=-8 h x I-4 h y A_L.
a1_order2 = -8 * h * (x_series[2] + 2 * y_series[2])
e2_order4 = -8 * h * x_minus_y[4]
t2_order2 = -8 * h * x_series[2]
check(a1_order2 == -24 * h**2 * delta1, "A1 leading coefficient")
check(e2_order4 == -F(32, 3) * h**4 * (2 * delta2 - delta1),
      "E2 leading coefficient")
check(t2_order2 == -8 * h**2 * delta1, "T2 leading coefficient")

witness = (F(1), F(13), F(1), F(2))
check(delta1.evaluate(witness) == 7, "witness delta1")
check(delta2.evaluate(witness) == 16, "witness delta2")
check(factor_defect[4].evaluate(witness) == -F(16, 3),
      "witness factor defect")
check(a1_order2.evaluate(witness) == -168, "witness A1")
check(t2_order2.evaluate(witness) == -56, "witness T2")
check(e2_order4.evaluate(witness) == -F(800, 3), "witness E2")


# ---------------------------------------------------------------------------
# Independent exact Pauli double-commutator and L(K4) sector replay.
# ---------------------------------------------------------------------------


def zero_matrix():
    return [[0 for _ in range(DIM)] for _ in range(DIM)]


def matrix_add(left, right, scale_right=1):
    return [[left[i][j] + scale_right * right[i][j]
             for j in range(DIM)] for i in range(DIM)]


def matrix_mul(left, right):
    out = zero_matrix()
    for i in range(DIM):
        for k in range(DIM):
            if left[i][k]:
                for j in range(DIM):
                    if right[k][j]:
                        out[i][j] += left[i][k] * right[k][j]
    return out


def commutator(left, right):
    return matrix_add(matrix_mul(left, right), matrix_mul(right, left), -1)


def x_matrix(site):
    out = zero_matrix()
    for state in range(DIM):
        out[state ^ (1 << site)][state] = 1
    return out


def z_matrix(site):
    out = zero_matrix()
    for state in range(DIM):
        out[state][state] = 1 if ((state >> site) & 1) == 0 else -1
    return out


def scaled(matrix, scalar):
    return [[scalar * value for value in row] for row in matrix]


X = [x_matrix(site) for site in range(4)]
Z = [z_matrix(site) for site in range(4)]
H_transverse = zero_matrix()
for matrix in X:
    H_transverse = matrix_add(H_transverse, matrix, -1)

edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
pair_matrices = [matrix_mul(Z[a], Z[b]) for a, b in edges]

for row, edge_b in enumerate(edges):
    for column, edge_a in enumerate(edges):
        actual = commutator(commutator(H_transverse, pair_matrices[row]),
                            pair_matrices[column])
        common = set(edge_b) & set(edge_a)
        if edge_b == edge_a:
            expected = scaled(matrix_add(X[edge_b[0]], X[edge_b[1]]), -4)
        elif len(common) == 1:
            shared = next(iter(common))
            other_b = next(q for q in edge_b if q != shared)
            other_a = next(q for q in edge_a if q != shared)
            expected = scaled(matrix_mul(matrix_mul(X[shared], Z[other_b]),
                                         Z[other_a]), -4)
        else:
            expected = zero_matrix()
        check(actual == expected, f"double commutator {edge_b},{edge_a}")

# The diagonal BREAK Hamiltonian commutes with every pair query.
H_break = zero_matrix()
for state in range(DIM):
    H_break[state][state] = bin(state).count("1")  # any diagonal degree law suffices
for pair, matrix in zip(edges, pair_matrices):
    check(commutator(H_break, matrix) == zero_matrix(), f"BREAK pair {pair}")

line = [[int(i != j and len(set(edge_i) & set(edge_j)) == 1)
         for j, edge_j in enumerate(edges)] for i, edge_i in enumerate(edges)]


def matvec(matrix, vector):
    return [sum(entry * value for entry, value in zip(row, vector))
            for row in matrix]


sector_vectors = {
    4: [[1, 1, 1, 1, 1, 1]],
    -2: [[1, -1, 0, 0, -1, 1], [1, 0, -1, -1, 0, 1]],
    0: [[1, 0, 0, 0, 0, -1], [0, 1, 0, 0, -1, 0],
        [0, 0, 1, -1, 0, 0]],
}
all_vectors = []
for eigenvalue, vectors in sector_vectors.items():
    for vector in vectors:
        check(matvec(line, vector) == [eigenvalue * q for q in vector],
              f"line-graph eigenvector {eigenvalue},{vector}")
        all_vectors.append(vector)

# Exact rational Gaussian elimination confirms the six vectors are a basis.
matrix = [[F(all_vectors[column][row]) for column in range(6)] for row in range(6)]
rank = 0
for column in range(6):
    pivot = next((row for row in range(rank, 6) if matrix[row][column]), None)
    if pivot is None:
        continue
    matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
    pivot_value = matrix[rank][column]
    matrix[rank] = [entry / pivot_value for entry in matrix[rank]]
    for row in range(6):
        if row != rank and matrix[row][column]:
            factor = matrix[row][column]
            matrix[row] = [a - factor * b for a, b in zip(matrix[row], matrix[rank])]
    rank += 1
check(rank == 6, "complete A1+E2+T2 sector basis")

print(f"N0_COUNTS parents={a0}+{guards} children={b0} active={active_links} nonedges={blank_nonedges}")
print("SYMBOLIC_DEFECT -(16/3)*h^3*U_d*s^4")
print("WITNESS_SECTORS A1=-168*s^2 E2=-(800/3)*s^4 T2=-56*s^2")
print(f"PASS__INDEPENDENT_GL6U_REPLAY__{checks}/{checks}")
