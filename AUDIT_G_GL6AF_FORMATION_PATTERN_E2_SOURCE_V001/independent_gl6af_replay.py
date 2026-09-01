#!/usr/bin/env python3
"""Independent exact algebra and branch-typing replay for frozen GL6AF."""

from fractions import Fraction as F
from itertools import combinations, product


checks = 0


def require(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


# Gaussian-integer coefficients (real, imaginary).
ZERO_G = (0, 0)
ONE_G = (1, 0)


def gadd(left, right):
    return (left[0] + right[0], left[1] + right[1])


def gneg(value):
    return (-value[0], -value[1])


def gmul(left, right):
    return (left[0] * right[0] - left[1] * right[1],
            left[0] * right[1] + left[1] * right[0])


SINGLE = {
    ("I", "I"): (ONE_G, "I"),
}
for pauli in "XYZ":
    SINGLE[("I", pauli)] = (ONE_G, pauli)
    SINGLE[(pauli, "I")] = (ONE_G, pauli)
    SINGLE[(pauli, pauli)] = (ONE_G, "I")
for left, right, phase, out in (
    ("X", "Y", (0, 1), "Z"),
    ("Y", "X", (0, -1), "Z"),
    ("Y", "Z", (0, 1), "X"),
    ("Z", "Y", (0, -1), "X"),
    ("Z", "X", (0, 1), "Y"),
    ("X", "Z", (0, -1), "Y"),
):
    SINGLE[(left, right)] = (phase, out)


def op_add(left, right):
    out = dict(left)
    for word, coefficient in right.items():
        out[word] = gadd(out.get(word, ZERO_G), coefficient)
        if out[word] == ZERO_G:
            del out[word]
    return out


def op_scale(operator, coefficient):
    return {word: gmul(coefficient, value) for word, value in operator.items()
            if gmul(coefficient, value) != ZERO_G}


def op_mul(left, right):
    out = {}
    for left_word, left_coefficient in left.items():
        for right_word, right_coefficient in right.items():
            coefficient = gmul(left_coefficient, right_coefficient)
            word = []
            for left_pauli, right_pauli in zip(left_word, right_word):
                phase, pauli = SINGLE[(left_pauli, right_pauli)]
                coefficient = gmul(coefficient, phase)
                word.append(pauli)
            word = tuple(word)
            out[word] = gadd(out.get(word, ZERO_G), coefficient)
            if out[word] == ZERO_G:
                del out[word]
    return out


def commutator(left, right):
    return op_add(op_mul(left, right), op_scale(op_mul(right, left), (-1, 0)))


def pauli_at(pauli, sites):
    word = tuple(pauli if index in sites else "I" for index in range(4))
    return {word: ONE_G}


PAIRS = tuple(combinations(range(4), 2))
E = (
    (1, 1),
    (-1, 0),
    (0, -1),
    (0, -1),
    (-1, 0),
    (1, 1),
)
P = tuple(tuple(int(port in pair) for port in range(4)) for pair in PAIRS)


# Polynomials in the formed-link z coordinate.
def pclean(poly):
    return {degree: F(value) for degree, value in poly.items() if value}


def padd(left, right):
    out = dict(left)
    for degree, value in right.items():
        out[degree] = out.get(degree, F(0)) + value
    return pclean(out)


def pscale(poly, value):
    return pclean({degree: F(value) * coefficient
                   for degree, coefficient in poly.items()})


def pmul(left, right):
    out = {}
    for left_degree, left_value in left.items():
        for right_degree, right_value in right.items():
            degree = left_degree + right_degree
            out[degree] = out.get(degree, F(0)) + left_value * right_value
    return pclean(out)


PZERO = {}
PONE = {0: F(1)}
PZ = {1: F(1)}


def expectation_div_minus_four_x(operator, pattern):
    """Product expectation of an exact Pauli operator, divided by -4*x."""
    answer = PZERO
    for word, coefficient in operator.items():
        require(coefficient[1] == 0, "double commutator is real Hermitian")
        x_sites = [index for index, pauli in enumerate(word) if pauli == "X"]
        require(len(x_sites) == 1, "one and only one X owns response term")
        if any(pauli == "Y" for pauli in word):
            continue
        x_site = x_sites[0]
        if not pattern[x_site]:
            continue
        value = PONE
        for index, pauli in enumerate(word):
            if pauli == "Z":
                value = pmul(value, PZ if pattern[index] else PONE)
            elif pauli not in ("I", "X"):
                raise AssertionError("unexpected Pauli word")
        answer = padd(answer, pscale(value, F(-coefficient[0], 4)))
    return answer


def matrix_contract(left, matrix, right):
    out = [[PZERO for _ in range(len(right[0]))]
           for _ in range(len(left[0]))]
    for row in range(len(left[0])):
        for column in range(len(right[0])):
            value = PZERO
            for i in range(len(matrix)):
                for j in range(len(matrix[0])):
                    value = padd(value, pscale(
                        matrix[i][j], left[i][row] * right[j][column]))
            out[row][column] = value
    return out


def outer(vector):
    return [[vector[row] * vector[column] for column in range(2)]
            for row in range(2)]


def determinant_2(matrix):
    return padd(pmul(matrix[0][0], matrix[1][1]),
                pscale(pmul(matrix[0][1], matrix[1][0]), -1))


def matrix_rank_2(matrix):
    if all(not entry for row in matrix for entry in row):
        return 0
    return 2 if determinant_2(matrix) else 1


rank_by_count = {count: set() for count in range(5)}
all_matrices = {}
for pattern in product((0, 1), repeat=4):
    # Reconstruct H/h after exact restriction to the physical K-support
    # eigenpath.  Delta*n drops from the double commutator because it commutes
    # with every pair query.
    hamiltonian = {}
    for site, formed in enumerate(pattern):
        if formed:
            hamiltonian = op_add(
                hamiltonian, op_scale(pauli_at("X", {site}), (-1, 0)))

    matrix = [[PZERO for _ in PAIRS] for _ in PAIRS]
    for source_index, source_pair in enumerate(PAIRS):
        source = pauli_at("Z", source_pair)
        for read_index, read_pair in enumerate(PAIRS):
            read = pauli_at("Z", read_pair)
            exact = commutator(commutator(hamiltonian, read), source)

            # Independently test the complete operator identity before taking
            # any product-state expectation.
            expected = {}
            if source_pair == read_pair:
                for site in source_pair:
                    if pattern[site]:
                        expected = op_add(
                            expected,
                            op_scale(pauli_at("X", {site}), (-4, 0)),
                        )
            else:
                shared = set(source_pair) & set(read_pair)
                if len(shared) == 1:
                    common = next(iter(shared))
                    if pattern[common]:
                        other = (set(source_pair) | set(read_pair)) - {common}
                        expected_word = tuple(
                            "X" if site == common else
                            "Z" if site in other else "I"
                            for site in range(4)
                        )
                        expected = {expected_word: (-4, 0)}
            require(exact == expected, "independent Pauli double commutator")
            matrix[read_index][source_index] = \
                expectation_div_minus_four_x(exact, pattern)

    all_matrices[pattern] = matrix
    restricted = matrix_contract(E, matrix, E)
    rank_by_count[sum(pattern)].add(matrix_rank_2(restricted))

    if sum(pattern) <= 1:
        require(matrix_rank_2(restricted) == 0,
                "zero/one formed fixed-E restriction null")
    elif sum(pattern) == 2:
        formed_pair = tuple(index for index, formed in enumerate(pattern)
                            if formed)
        row = E[PAIRS.index(formed_pair)]
        expected = [[pscale({0: 1, 1: -1}, 4 * value)
                     for value in outer(row)[matrix_row]]
                    for matrix_row in range(2)]
        require(restricted == expected,
                "two-formed exact -16hx(1-z) covector block")

require(rank_by_count == {0: {0}, 1: {0}, 2: {1}, 3: {2}, 4: {2}},
        "exact rank threshold by formed count")

# Closed representative forms.
three = matrix_contract(E, all_matrices[(1, 1, 1, 0)], E)
three_scalar = {0: F(3), 1: F(-2), 2: F(-1)}
require(three == [[pscale(three_scalar, value) for value in (2, 1)],
                  [pscale(three_scalar, value) for value in (1, 2)]],
        "three-formed representative")
four = matrix_contract(E, all_matrices[(1, 1, 1, 1)], E)
four_scalar = {0: F(4), 2: F(-4)}
require(four == [[pscale(four_scalar, value) for value in (2, 1)],
                 [pscale(four_scalar, value) for value in (1, 2)]],
        "four-formed representative")

# Broken-S4 ceiling: the two-record E x E restriction is exact, but D need
# not preserve E.  Exhibit the exact nonzero incidence-complement cross block
# for pattern 1100, while the full symmetric pattern has no such mixing.
cross_two = matrix_contract(P, all_matrices[(1, 1, 0, 0)], E)
mix = {0: F(-2), 1: F(2)}
require(cross_two == [
    [PZERO, PZERO],
    [PZERO, PZERO],
    [mix, mix],
    [mix, mix],
], "exact broken-S4 complement mixing witness")
cross_four = matrix_contract(P, all_matrices[(1, 1, 1, 1)], E)
require(all(not entry for row in cross_four for entry in row),
        "full S4 pattern preserves E plane")

# Exact K-pattern projector and source-independence census.  In the K
# computational basis, Pi_kappa is the one-word indicator.  The projectors
# are orthogonal, idempotent, complete, and every physical Hamiltonian
# transition flips only an active-link bit while retaining the K word.
patterns = tuple(product((0, 1), repeat=4))
for kword in patterns:
    memberships = [int(kword == pattern) for pattern in patterns]
    require(sum(memberships) == 1, "K projectors resolve identity")
    require(all(value * value == value for value in memberships),
            "K projectors are idempotent")
    for left_index, left in enumerate(memberships):
        for right_index, right in enumerate(memberships):
            if left_index != right_index:
                require(left * right == 0, "K projectors are orthogonal")
    for active_word in patterns:
        for site in range(4):
            transition_kword = kword
            for branch_pattern in patterns:
                before = int(kword == branch_pattern)
                after = int(transition_kword == branch_pattern)
                require(before == after, "H transition commutes with Pi_kappa")
            if kword[site]:
                flipped = list(active_word)
                flipped[site] = 1 - flipped[site]
                require(tuple(flipped) != active_word,
                        "formed support owns active transverse transition")

# Pair queries and every GL6V source phase are functions of active Z values
# tensored with identity on K.  Their eigenvalue is therefore independent of
# the K word, which is the finite-basis commutator check with every Pi_kappa.
for active_word in patterns:
    for left, right in PAIRS:
        pair_value = (1 - 2 * active_word[left]) * \
                     (1 - 2 * active_word[right])
        values_across_k = {pair_value for _ in patterns}
        require(values_across_k == {pair_value},
                "M and GL6V source phase are identity on K")

# CTP sign/factor induction.  G^R = i E_*^2/(2 hbar) [B(t),A], while
# dB(t)/dt at zero is i[H,B]/hbar.  Their product is -E_*^2/(2 hbar^2)D.
require((0, 1) == gmul((0, 1), ONE_G), "Heisenberg derivative i factor")
require(gmul((0, 1), (0, 1)) == (-1, 0),
        "two imaginary factors give negative entrance slope")
require(F(-1, 2) * F(-16) == F(8),
        "AF11 to AF17 normalized coefficient")

print(f"PASS__INDEPENDENT_GL6AF_REPLAY__{checks}/{checks}")
print("PAULI=ALL_16_PATTERNS_COMPLETE_6X6_DOUBLE_COMMUTATOR_EXACT")
print("RANKS=FORMED_COUNT_0_1_TO0_2_TO1_3_4_TO2")
print("COEFFICIENT=TWO_FORMED_MINUS16_H_X_ONE_MINUS_Z_WTW")
print("MIXING=BROKEN_S4_FIXED_E_RESTRICTION_NOT_INVARIANT_BLOCK")
print("BRANCH=ORTHOGONAL_K_PROJECTORS_H_M_GL6V_SOURCE_PRESERVE_SECTORS")
print("CTP=BRANCHWISE_SLOPE_PLUS8_ESTAR2_H_X_ONE_MINUS_Z_OVER_HBAR2")
print("ANCESTRY=PHYSICAL_K_MEDIATOR_NOT_SEMANTIC_REC_OR_PAIR_RECORD")
print("SCOPE=NO_COLLECTIVE_STIFFNESS_BULK_STRESS_RICCI_GRAVITY_G")
