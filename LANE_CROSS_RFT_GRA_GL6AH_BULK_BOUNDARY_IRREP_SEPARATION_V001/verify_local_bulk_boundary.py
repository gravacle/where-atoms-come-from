#!/usr/bin/env python3
"""Exact local replays for the mutable GL6AH bulk/boundary separation.

The Hamiltonians are the GL6AB/GL6AG integer witness at delta=0.  Dense
integer vectors suffice for the two-cell (2^8) and three-cell (2^12)
blocks.  A separate sparse polynomial replay retains h, U_d, and a
dimensionless tag g on the direct connector through order six.
"""

from fractions import Fraction
from math import factorial


PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
E_ROWS = ((1, 1), (-1, 0), (0, -1), (0, -1), (-1, 0), (1, 1))
Q12 = 63371264
Q16 = 123422773248
checks = 0


def require(condition, message):
    global checks
    if not condition:
        raise AssertionError(message)
    checks += 1


def choose(n, k):
    out = 1
    for j in range(1, k + 1):
        out = out * (n - k + j) // j
    return out


def e_contract(pair_vector):
    return tuple(sum(E_ROWS[row][column] * pair_vector[row]
                     for row in range(6))
                 for column in range(2))


def vector_sub(left, right):
    return tuple(a - b for a, b in zip(left, right))


def vector_add(*vectors):
    return tuple(sum(entries) for entries in zip(*vectors))


class DenseReplay:
    """Exact delta=0 witness on complete four-link cells."""

    def __init__(self, cells, source_mask, bridges, max_order):
        self.cells = cells
        self.links = 4 * cells
        self.dimension = 1 << self.links
        self.formed = [False] * self.links
        for port in range(4):
            self.formed[port] = bool((source_mask >> port) & 1)
        for link in range(4, self.links):
            self.formed[link] = True

        interactions = []
        for cell in range(cells):
            for left, right in PAIRS:
                interactions.append((4 * cell + left, 4 * cell + right))
        interactions.extend(bridges)
        self.diagonal = [
            sum(2 * ((word >> left) & 1) * ((word >> right) & 1)
                for left, right in interactions)
            for word in range(self.dimension)
        ]

        blank = [0] * self.dimension
        blank[0] = 1
        self.powers = [blank]
        for _ in range(max_order):
            self.powers.append(self.apply_h(self.powers[-1]))

    def apply_h(self, incoming):
        outgoing = [0] * self.dimension
        for word, amplitude in enumerate(incoming):
            if amplitude == 0:
                continue
            outgoing[word] += self.diagonal[word] * amplitude
            for link, formed in enumerate(self.formed):
                if formed:
                    outgoing[word ^ (1 << link)] -= amplitude
        return outgoing

    def pair_raw(self, cell, order):
        answer = [0] * 6
        for split in range(order + 1):
            matrix = [0] * 6
            bra = self.powers[order - split]
            ket = self.powers[split]
            for word in range(self.dimension):
                if bra[word] == 0 or ket[word] == 0:
                    continue
                product = bra[word] * ket[word]
                z = [1 - 2 * ((word >> (4 * cell + port)) & 1)
                     for port in range(4)]
                for row, (left, right) in enumerate(PAIRS):
                    matrix[row] += product * z[left] * z[right]
            weight = (-1 if split & 1 else 1) * choose(order, split)
            for row in range(6):
                answer[row] += weight * matrix[row]
        return tuple(answer)

    def e_raw(self, cell, order):
        return e_contract(self.pair_raw(cell, order))


# Sparse polynomials use exponent triples (power(h), power(U_d), power(g)).
def poly_add_term(poly, exponent, coefficient):
    if coefficient == 0:
        return
    poly[exponent] = poly.get(exponent, 0) + coefficient
    if poly[exponent] == 0:
        del poly[exponent]


def poly_add_scaled(target, source, coefficient):
    for exponent, value in source.items():
        poly_add_term(target, exponent, coefficient * value)


def poly_times_monomial(poly, exponent, coefficient):
    out = {}
    for old, value in poly.items():
        new = tuple(old[j] + exponent[j] for j in range(3))
        poly_add_term(out, new, coefficient * value)
    return out


def poly_times_poly(left, right):
    out = {}
    for a, x in left.items():
        for b, y in right.items():
            exponent = tuple(a[j] + b[j] for j in range(3))
            poly_add_term(out, exponent, x * y)
    return out


def polynomial_powers(source_formed):
    links = 8
    dimension = 1 << links
    formed = [False] * links
    formed[0] = source_formed
    for link in range(4, 8):
        formed[link] = True
    within = [(4 * cell + left, 4 * cell + right)
              for cell in range(2) for left, right in PAIRS]

    powers = [[{} for _ in range(dimension)] for _ in range(7)]
    powers[0][0] = {(0, 0, 0): 1}
    for order in range(6):
        incoming = powers[order]
        outgoing = powers[order + 1]
        for word, poly in enumerate(incoming):
            if not poly:
                continue
            within_count = sum(((word >> left) & 1) *
                               ((word >> right) & 1)
                               for left, right in within)
            if within_count:
                poly_add_scaled(
                    outgoing[word],
                    poly_times_monomial(poly, (0, 1, 0),
                                        2 * within_count),
                    1,
                )
            direct_occupied = ((word >> 0) & 1) * ((word >> 4) & 1)
            if direct_occupied:
                poly_add_scaled(
                    outgoing[word],
                    poly_times_monomial(poly, (0, 1, 1), 2),
                    1,
                )
            for link, is_formed in enumerate(formed):
                if is_formed:
                    poly_add_scaled(
                        outgoing[word ^ (1 << link)],
                        poly_times_monomial(poly, (1, 0, 0), -1),
                        1,
                    )
    return powers


def polynomial_pair_raw(powers, pair, order):
    answer = {}
    left, right = pair
    for split in range(order + 1):
        matrix = {}
        for word in range(256):
            bra = powers[order - split][word]
            ket = powers[split][word]
            if not bra or not ket:
                continue
            sign = ((1 - 2 * ((word >> (4 + left)) & 1)) *
                    (1 - 2 * ((word >> (4 + right)) & 1)))
            poly_add_scaled(matrix, poly_times_poly(bra, ket), sign)
        poly_add_scaled(answer, matrix,
                        (-1 if split & 1 else 1) * choose(order, split))
    return answer


def polynomial_sub(left, right):
    out = dict(left)
    poly_add_scaled(out, right, -1)
    return out


# Representation checks: W = A1 + T2 + E in the inherited pair order.
incidence = [[int(port in pair) for port in range(4)] for pair in PAIRS]
for column in range(2):
    for port in range(4):
        require(sum(E_ROWS[row][column] * incidence[row][port]
                    for row in range(6)) == 0,
                "E is the pair-incidence null")

one = (1, 1, 1, 1, 1, 1)
for b in range(4):
    u_in = tuple(row[b] for row in incidence)
    s_b = tuple(2 * value - 1 for value in u_in)
    require(all(Fraction(u_in[j]) ==
                Fraction(one[j] + s_b[j], 2) for j in range(6)),
            "u_in=(one+s_b)/2")
    require(sum(s_b) == 0, "T2 star is orthogonal to A1")
    require(sum(value * value for value in s_b) == 6,
            "T2 star normalization")
    require(e_contract(u_in) == (0, 0), "one-port star has exact E null")
    u_out = tuple(1 - value for value in u_in)
    require(e_contract(u_out) == (0, 0),
            "one-port complement has exact E null")

for b in range(4):
    fixed_rows = [row for row, pair in enumerate(PAIRS) if b in pair]
    other_rows = [row for row, pair in enumerate(PAIRS) if b not in pair]
    require(len(fixed_rows) == 3 and len(other_rows) == 3,
            "H_b has the two expected pair orbits")

# Symbolic two-cell q6: delta=0, direct connector tagged by g.
poly_sham = polynomial_powers(False)
poly_source = polynomial_powers(True)
expected_incident = {(4, 2, 2): -128}
for row, pair in enumerate(PAIRS):
    matched = polynomial_sub(polynomial_pair_raw(poly_source, pair, 6),
                             polynomial_pair_raw(poly_sham, pair, 6))
    require(matched == (expected_incident if 0 in pair else {}),
            "symbolic q6 equals -128*h^4*U_d^2*g^2*u_0^in")
for order in range(6):
    for pair in PAIRS:
        matched = polynomial_sub(
            polynomial_pair_raw(poly_source, pair, order),
            polynomial_pair_raw(poly_sham, pair, order),
        )
        require(matched == {}, "symbolic direct coefficient below q6")

# Dense two-cell replay through q16, including the exact all-orders symmetry
# pattern in this finite jet and the bridge-off factorization diagnostic.
direct_bridge = ((0, 4),)
two_sham = DenseReplay(2, 0, direct_bridge, 16)
two_source = DenseReplay(2, 1, direct_bridge, 16)
for order in range(6):
    require(vector_sub(two_source.pair_raw(1, order),
                       two_sham.pair_raw(1, order)) == (0,) * 6,
            "two-cell matched coefficient below q6")
direct_q6 = vector_sub(two_source.pair_raw(1, 6),
                       two_sham.pair_raw(1, 6))
require(direct_q6 == (-128, -128, -128, 0, 0, 0),
        "two-cell full six-pair q6")
require(direct_q6 == vector_add(tuple(-64 for _ in range(6)),
                                tuple(-64 * x for x in
                                      (1, 1, 1, -1, -1, -1))),
        "equal A1 and T2 coordinate coefficients")
for order in range(17):
    direct_difference = vector_sub(two_source.pair_raw(1, order),
                                   two_sham.pair_raw(1, order))
    require(e_contract(direct_difference) == (0, 0),
            "direct one-port E null through replay horizon")

off_sham = DenseReplay(2, 0, (), 16)
off_source = DenseReplay(2, 1, (), 16)
for order in range(17):
    require(vector_sub(off_source.pair_raw(1, order),
                       off_sham.pair_raw(1, order)) == (0,) * 6,
            "two-cell bridge-off factorization jet")

require(Fraction((-1) * (-128), factorial(6)) == Fraction(8, 45),
        "i^6 sign and q6 factorial")
require(Fraction(Q12, factorial(12)) == Fraction(5626, 42525),
        "q12 factorial")
require(Fraction(Q16, factorial(16)) ==
        Fraction(1116019, 189189000), "q16 factorial")

# Three-cell chain: source port a=1 -> receiver port b=0, followed by
# receiver port d=2 -> helper port 1.
direct = (1, 4)
helper = (4 + 2, 8 + 1)
chain_bridges = (direct, helper)
chain_sham = DenseReplay(3, 0, chain_bridges, 16)
chain_source = DenseReplay(3, 1 << 1, chain_bridges, 16)
for order in range(12):
    require(vector_sub(chain_source.e_raw(1, order),
                       chain_sham.e_raw(1, order)) == (0, 0),
            "three-cell matched E coefficient below q12")
chain_q12 = vector_sub(chain_source.e_raw(1, 12),
                       chain_sham.e_raw(1, 12))
require(chain_q12 == (-Q12, 0), "chain q12 is +Q12*w_02")

# Two cells (or a disconnected helper) remain exactly E-null at q12.
direct_only_sham = DenseReplay(3, 0, (direct,), 12)
direct_only_source = DenseReplay(3, 1 << 1, (direct,), 12)
helper_only_sham = DenseReplay(3, 0, (helper,), 12)
helper_only_source = DenseReplay(3, 1 << 1, (helper,), 12)
none_sham = DenseReplay(3, 0, (), 12)
none_source = DenseReplay(3, 1 << 1, (), 12)
q12_support = vector_add(
    chain_q12,
    tuple(-x for x in vector_sub(direct_only_source.e_raw(1, 12),
                                 direct_only_sham.e_raw(1, 12))),
    tuple(-x for x in vector_sub(helper_only_source.e_raw(1, 12),
                                 helper_only_sham.e_raw(1, 12))),
    vector_sub(none_source.e_raw(1, 12), none_sham.e_raw(1, 12)),
)
require(q12_support == (-Q12, 0),
        "two-connector Mobius support owns chain q12")

# A two-connector path which omits the direct source--receiver bridge must
# traverse two distinct ports of its middle cell.  That requires a further
# within-cell diagonal pair.  The exact bounded lemma is zero in all six
# receiver coordinates below q14 and opens at q14 only in the receiver's
# one-port star, hence remains E-null.  This excludes the distance-two wedge
# from the complete q12 E ledger without extrapolating beyond this motif.
alternate_sham = DenseReplay(3, 0, chain_bridges, 14)
alternate_source = DenseReplay(3, 1 << 1, chain_bridges, 14)
for order in range(14):
    alternate = vector_sub(alternate_source.pair_raw(2, order),
                           alternate_sham.pair_raw(2, order))
    require(alternate == (0,) * 6,
            "alternate source-helper-receiver path below q14")
alternate_q14 = vector_sub(alternate_source.pair_raw(2, 14),
                           alternate_sham.pair_raw(2, 14))
require(alternate_q14 == (-14721024, 0, 0, -14721024, -14721024, 0),
        "alternate distance-two path q14 one-port star")
require(e_contract(alternate_q14) == (0, 0),
        "alternate distance-two q14 remains E-null")

# Source pair {x=0,a=1}; the pair includes the direct source port a.
chain_x = DenseReplay(3, 1 << 0, chain_bridges, 16)
chain_ax = DenseReplay(3, (1 << 0) | (1 << 1), chain_bridges, 16)
for order in range(16):
    mobius = vector_add(
        chain_ax.e_raw(1, order),
        tuple(-x for x in chain_source.e_raw(1, order)),
        tuple(-x for x in chain_x.e_raw(1, order)),
        chain_sham.e_raw(1, order),
    )
    require(mobius == (0, 0), "source-pair remote Mobius below q16")
chain_q16 = vector_add(
    chain_ax.e_raw(1, 16),
    tuple(-x for x in chain_source.e_raw(1, 16)),
    tuple(-x for x in chain_x.e_raw(1, 16)),
    chain_sham.e_raw(1, 16),
)
require(chain_q16 == (Q16, 0), "chain q16 is -Q16*w_02")

# Endpoint incidence identity and the N=1 boundary specialization.
for b in range(4):
    incident_rows = [PAIRS.index(tuple(sorted((b, d))))
                     for d in range(4) if d != b]
    require(tuple(sum(E_ROWS[row][coordinate] for row in incident_rows)
                  for coordinate in range(2)) == (0, 0),
            "sum_{d!=b} w_bd=0")

# Receiver n=e_1 with direct port b=0: r=(1,0,1,1), and d=b is excluded.
r = (1, 0, 1, 1)
eta = [0, 0]
for d in range(1, 4):
    row = PAIRS.index((0, d))
    for coordinate in range(2):
        eta[coordinate] += r[d] * E_ROWS[row][coordinate]
require(tuple(eta) == tuple(-x for x in E_ROWS[0]),
        "N=1 endpoint eta is -w_01")
require(vector_add(tuple(Q12 * x for x in eta)) == (-Q12, -Q12),
        "N=1 q12 helper sum")
require(tuple(-Q16 * x for x in eta) == (Q16, Q16),
        "N=1 q16 helper sum")
for b in range(4):
    homogeneous_eta = [0, 0]
    for d in range(4):
        if d == b:
            continue
        row = PAIRS.index(tuple(sorted((b, d))))
        for coordinate in range(2):
            homogeneous_eta[coordinate] += 3 * E_ROWS[row][coordinate]
    require(tuple(homogeneous_eta) == (0, 0),
            "homogeneous r_nd=3 endpoint cancellation")

print(f"PASS__GL6AH_LOCAL_BULK_BOUNDARY__{checks}/{checks}")
print("DIRECT_Q6=-128_H4_UD2_U_B_IN;EXPECTATION=+8/45")
print("TAG_POLYNOMIAL=-128_G2_AT_H_EQ_UD_EQ_1_DELTA_EQ_0")
print("DIRECT_IRREPS=EQUAL_A1_T2;E_EXACT_NULL")
print("CHAIN_Q12=+63371264_W_BD;CHAIN_Q16=-123422773248_W_BD")
print("ALTERNATE_NO_DIRECT=ZERO_BELOW_Q14;Q14=-14721024_U_1_IN;E_NULL")
print("ENDPOINT=ETA_SUM_R_ND_W_BD;HOMOGENEOUS_ETA_ZERO")
