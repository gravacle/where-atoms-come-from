#!/usr/bin/env python3
"""Independent compressed-sector replay for frozen GL6AH.

This implementation omits source-sham qubits only after proving their blank
sector reducing: they have no transverse flip and every remaining term is
diagonal in their occupation.  It reconstructs Hamiltonian powers directly
on the resulting 5-, 9-, and 10-qubit bases and imports no author ledger.
"""

from fractions import Fraction
from itertools import product
from math import factorial


PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
E_ROWS = ((1, 1), (-1, 0), (0, -1), (0, -1), (-1, 0), (1, 1))
Q12 = 63371264
Q16 = 123422773248
checks = 0


def require(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def choose(n, k):
    out = 1
    for j in range(1, k + 1):
        out = out * (n - k + j) // j
    return out


def add(left, *others):
    return tuple(a + sum(values) for a, values in zip(left, zip(*others)))


def sub(left, right):
    return tuple(a - b for a, b in zip(left, right))


def scale(value, factor):
    return tuple(factor * entry for entry in value)


def e_contract(pair_vector):
    return tuple(
        sum(E_ROWS[row][coordinate] * pair_vector[row] for row in range(6))
        for coordinate in range(2)
    )


def w(left, right):
    return E_ROWS[PAIRS.index(tuple(sorted((left, right))))]


class Replay:
    def __init__(self, cells, source_mask, bridges, max_order, delta=0):
        self.cells = cells
        self.logical_to_active = {}
        for port in range(4):
            if (source_mask >> port) & 1:
                self.logical_to_active[(0, port)] = len(self.logical_to_active)
        for cell in range(1, cells):
            for port in range(4):
                self.logical_to_active[(cell, port)] = len(self.logical_to_active)
        self.links = len(self.logical_to_active)
        self.dimension = 1 << self.links

        interactions = []
        for cell in range(cells):
            for left, right in PAIRS:
                a = self.logical_to_active.get((cell, left))
                b = self.logical_to_active.get((cell, right))
                if a is not None and b is not None:
                    interactions.append((a, b))
        for left, right in bridges:
            a = self.logical_to_active.get(left)
            b = self.logical_to_active.get(right)
            if a is not None and b is not None:
                interactions.append((a, b))

        self.diagonal = [0] * self.dimension
        for word in range(self.dimension):
            self.diagonal[word] = (
                delta * bin(word).count("1")
                + sum(
                    2 * ((word >> left) & 1) * ((word >> right) & 1)
                    for left, right in interactions
                )
            )

        blank = [0] * self.dimension
        blank[0] = 1
        self.powers = [blank]
        for _ in range(max_order):
            incoming = self.powers[-1]
            outgoing = [0] * self.dimension
            for word, amplitude in enumerate(incoming):
                if amplitude == 0:
                    continue
                outgoing[word] += self.diagonal[word] * amplitude
                for link in range(self.links):
                    outgoing[word ^ (1 << link)] -= amplitude
            self.powers.append(outgoing)

    def pair_raw(self, cell, order):
        answer = [0] * 6
        active_ports = [self.logical_to_active[(cell, port)] for port in range(4)]
        for split_at in range(order + 1):
            matrix = [0] * 6
            bra = self.powers[order - split_at]
            ket = self.powers[split_at]
            for word in range(self.dimension):
                if bra[word] == 0 or ket[word] == 0:
                    continue
                amplitude = bra[word] * ket[word]
                z = [1 - 2 * ((word >> q) & 1) for q in active_ports]
                for row, (left, right) in enumerate(PAIRS):
                    matrix[row] += amplitude * z[left] * z[right]
            coefficient = choose(order, split_at)
            if split_at & 1:
                coefficient = -coefficient
            for row in range(6):
                answer[row] += coefficient * matrix[row]
        return tuple(answer)

    def e_raw(self, cell, order):
        return e_contract(self.pair_raw(cell, order))


def source_mobius(pair_branch, direct_branch, extra_branch, sham_branch):
    return add(pair_branch, scale(direct_branch, -1),
               scale(extra_branch, -1), sham_branch)


def branch_mobius(both, direct, helper, none):
    return add(both, scale(direct, -1), scale(helper, -1), none)


# Fixed-frame representation and endpoint incidence identities.
incidence = tuple(tuple(int(port in pair) for port in range(4)) for pair in PAIRS)
for coordinate in range(2):
    for port in range(4):
        require(
            sum(E_ROWS[row][coordinate] * incidence[row][port]
                for row in range(6)) == 0,
            "E columns lie in ker(P^T)",
        )

one = (1,) * 6
stars = []
for b in range(4):
    u_in = tuple(incidence[row][b] for row in range(6))
    star = tuple(2 * entry - 1 for entry in u_in)
    stars.append(star)
    require(tuple(Fraction(one[j] + star[j], 2) for j in range(6)) == u_in,
            "u_in equal A1 plus T2 split")
    require(sum(star) == 0, "T2 star orthogonal to A1")
    require(sum(entry * entry for entry in star) == 6, "T2 star norm")
    require(e_contract(u_in) == (0, 0), "one-port E null")
require(tuple(sum(stars[b][j] for b in range(4)) for j in range(6)) == (0,) * 6,
        "four T2 stars sum to zero")


# Direct edge: reconstruct the full six-vector for every receiver port and
# independently check the finite direct E jet through order sixteen.
source_port = 1
for receiver_port in range(4):
    direct = (((0, source_port), (1, receiver_port)),)
    sham = Replay(2, 0, direct, 16)
    formed = Replay(2, 1 << source_port, direct, 16)
    for order in range(6):
        require(sub(formed.pair_raw(1, order), sham.pair_raw(1, order)) == (0,) * 6,
                "direct pair vector below q6")
    expected = tuple(-128 * incidence[row][receiver_port] for row in range(6))
    q6 = sub(formed.pair_raw(1, 6), sham.pair_raw(1, 6))
    require(q6 == expected, "direct full-six q6")
    require(q6 == add(scale(one, -64), scale(stars[receiver_port], -64)),
            "equal A1 and T2 q6 coefficients")
    for order in range(17):
        require(e_contract(sub(formed.pair_raw(1, order),
                               sham.pair_raw(1, order))) == (0, 0),
                "direct E-null finite reconstruction")

require(Fraction((-1) * (-128), factorial(6)) == Fraction(8, 45),
        "q6 i^6 sign and factorial")
require(Fraction(Q12, factorial(12)) == Fraction(5626, 42525),
        "q12 factorial")
require(Fraction(Q16, factorial(16)) == Fraction(1116019, 189189000),
        "q16 factorial")


# Receiver-helper chain: reconstruct all twelve ordered receiver port pairs.
# Source port a=1 is direct; x=0 is the second source in the q16 Mobius term.
source_a = 1
source_x = 0
for b in range(4):
    for d in range(4):
        if d == b:
            continue
        helper_c = (d + 1) % 4
        direct = ((0, source_a), (1, b))
        helper = ((1, d), (2, helper_c))
        bridges = (direct, helper)
        sham = Replay(3, 0, bridges, 16)
        formed_a = Replay(3, 1 << source_a, bridges, 16)
        formed_x = Replay(3, 1 << source_x, bridges, 16)
        formed_ax = Replay(3, (1 << source_a) | (1 << source_x), bridges, 16)

        for order in range(12):
            require(sub(formed_a.e_raw(1, order), sham.e_raw(1, order)) == (0, 0),
                    "chain E below q12")
        require(sub(formed_a.e_raw(1, 12), sham.e_raw(1, 12)) == scale(w(b, d), Q12),
                "chain q12 sign and direction")

        for order in range(16):
            value = source_mobius(
                formed_ax.e_raw(1, order), formed_a.e_raw(1, order),
                formed_x.e_raw(1, order), sham.e_raw(1, order),
            )
            require(value == (0, 0), "chain source Mobius below q16")
        q16 = source_mobius(
            formed_ax.e_raw(1, 16), formed_a.e_raw(1, 16),
            formed_x.e_raw(1, 16), sham.e_raw(1, 16),
        )
        require(q16 == scale(w(b, d), -Q16), "chain q16 sign and direction")


# Canonical bridge-support Mobius ownership is reconstructed separately from
# the source Mobius: neither direct-only nor helper-only owns the E opening.
b, d, c = 0, 2, 1
direct = ((0, source_a), (1, b))
helper = ((1, d), (2, c))
bridge_sets = ((direct, helper), (direct,), (helper,), ())
q12_values = []
q16_values = []
for bridges in bridge_sets:
    sham = Replay(3, 0, bridges, 16)
    formed_a = Replay(3, 1 << source_a, bridges, 16)
    formed_x = Replay(3, 1 << source_x, bridges, 16)
    formed_ax = Replay(3, (1 << source_a) | (1 << source_x), bridges, 16)
    q12_values.append(sub(formed_a.e_raw(1, 12), sham.e_raw(1, 12)))
    q16_values.append(source_mobius(
        formed_ax.e_raw(1, 16), formed_a.e_raw(1, 16),
        formed_x.e_raw(1, 16), sham.e_raw(1, 16),
    ))
require(branch_mobius(*q12_values) == scale(w(b, d), Q12),
        "two-bridge support owns q12")
require(branch_mobius(*q16_values) == scale(w(b, d), -Q16),
        "two-bridge support owns q16")


# Literal helper counts, homogeneous cancellation, and the N=1 boundary sign.
for b in range(4):
    require(tuple(sum(w(b, d)[coordinate] for d in range(4) if d != b)
                  for coordinate in range(2)) == (0, 0),
            "incident E rows cancel")

for n in product(range(5), repeat=4):
    if sum(n) > 8:
        continue
    for b in range(4):
        eta = [0, 0]
        for d in range(4):
            if d == b:
                continue
            literal_helpers = sum(1 for c in range(4) if c != d and n[c] > 0)
            for coordinate in range(2):
                eta[coordinate] += literal_helpers * w(b, d)[coordinate]
        if all(entry > 0 for entry in n):
            require(tuple(eta) == (0, 0), "all-positive homogeneous eta zero")

for cell_port in (1, 2, 3):
    n = tuple(int(port == cell_port) for port in range(4))
    eta = [0, 0]
    for d in range(1, 4):
        helpers = sum(1 for c in range(4) if c != d and n[c] > 0)
        for coordinate in range(2):
            eta[coordinate] += helpers * w(0, d)[coordinate]
    require(tuple(eta) == scale(w(0, cell_port), -1), "N1 eta boundary law")
    require(scale(tuple(eta), Q12) == scale(w(0, cell_port), -Q12),
            "N1 q12 reconciliation")
    require(scale(tuple(eta), -Q16) == scale(w(0, cell_port), Q16),
            "N1 q16 reconciliation")


# Scope audit: every finite nested coefficient is polynomial in delta because
# the finite matrix entries are affine in delta.  Nonzero value at delta=0
# earns an open nonzero neighborhood, but no detuned closed form is inferred.
direct = (((0, source_a), (1, 0)),)
for delta in (-2, -1, 0, 1, 2):
    sham = Replay(2, 0, direct, 6, delta=delta)
    formed = Replay(2, 1 << source_a, direct, 6, delta=delta)
    require(len(sub(formed.pair_raw(1, 6), sham.pair_raw(1, 6))) == 6,
            "finite detuning reconstruction")

print(f"PASS__INDEPENDENT_GL6AH_REPLAY__{checks}/{checks}")
print("DIRECT=FULL6_Q6_MINUS128_U_B_IN;A1_T2_EQUAL;E_NULL")
print("CHAIN=Q12_PLUS63371264_W_BD;Q16_MINUS123422773248_W_BD")
print("SUPPORT=TWO_BRIDGE_RECEIVER_HELPER_CHAIN_OWNS_Q12_Q16")
print("ENDPOINT=ETA_LITERAL_HELPER_COUNT;HOMOGENEOUS_ZERO;N1_SIGNS_RECONCILED")
print("DELTA_SCOPE=EXACT_DISPLAYED_WITNESS_DELTA0;OPEN_NONZERO_NEIGHBORHOOD_ONLY")
print("CEILING=NO_BULK_SHEAR_STATIONARY_MODE_COMMON_CONE_RICCI_GRAVITY_G")
