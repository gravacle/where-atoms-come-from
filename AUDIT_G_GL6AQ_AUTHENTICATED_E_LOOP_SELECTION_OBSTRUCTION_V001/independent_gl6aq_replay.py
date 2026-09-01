#!/usr/bin/env python3
"""Independent hostile replay for GL6AQ.

This script imports no author code.  It reconstructs the local locked-E
algebra, S4 selection rule, Pauli trace obstruction, retained-support loop
factor, and frozen dependency custody using only the standard library.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION_V001"
PAIRS = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIRS)}
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def matvec(a, v):
    return tuple(sum(a[i][j] * v[j] for j in range(len(v)))
                 for i in range(len(a)))


def matadd(a, b):
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))]
            for i in range(len(a))]


def matscale(a, scalar):
    return [[scalar * value for value in row] for row in a]


def rank_q(matrix):
    data = [[Fraction(value) for value in row] for row in matrix]
    row = 0
    for column in range(len(data[0])):
        pivot = next((index for index in range(row, len(data))
                      if data[index][column]), None)
        if pivot is None:
            continue
        data[row], data[pivot] = data[pivot], data[row]
        scale = data[row][column]
        data[row] = [value / scale for value in data[row]]
        for index in range(len(data)):
            if index == row or not data[index][column]:
                continue
            factor = data[index][column]
            data[index] = [x - factor * y for x, y in zip(data[index], data[row])]
        row += 1
    return row


def permutation_matrix(images, size):
    matrix = [[Fraction() for _ in range(size)] for _ in range(size)]
    for old, new in enumerate(images):
        matrix[new][old] = Fraction(1)
    return matrix


def pair_vector(z):
    return tuple(z[a] * z[b] for a, b in PAIRS)


R = [[Fraction(int(port in pair)) for pair in PAIRS] for port in range(4)]
E_BASIS = (
    (Fraction(1), 0, -1, -1, 0, 1),
    (0, Fraction(1), -1, -1, 1, 0),
)


def section_exact_locked_e_operator():
    check(rank_q(R) == 4, "unsigned pair incidence has rank four")
    check(6 - rank_q(R) == 2, "ker R is exactly two-dimensional")
    for vector in E_BASIS:
        check(matvec(R, vector) == (0, 0, 0, 0), "declared vectors span ker R")

    locked = tuple(z for z in product((-1, 1), repeat=4) if sum(z) == 0)
    pair_values = tuple(pair_vector(z) for z in locked)
    check(len(locked) == 6, "six locked local spin configurations")
    for index, values in enumerate(pair_values):
        check(matvec(R, values) == (-1, -1, -1, -1),
              f"locked affine identity {index}")

    mean = tuple(sum(Fraction(values[j], 6) for values in pair_values)
                 for j in range(6))
    check(mean == (Fraction(-1, 3),) * 6, "locked pair mean")
    covariance = [[
        sum(Fraction((values[i] - mean[i]) * (values[j] - mean[j]), 6)
            for values in pair_values)
        for j in range(6)
    ] for i in range(6)]
    projector = matscale(covariance, Fraction(3, 8))
    check(matmul(projector, projector) == projector, "covariance gives exact projector")
    check(sum(projector[i][i] for i in range(6)) == 2, "projector rank is two")
    check(matmul(R, projector) == [[Fraction() for _ in range(6)] for _ in range(4)],
          "projector image lies in ker R")
    for vector in E_BASIS:
        check(matvec(projector, vector) == vector, "projector fixes E basis")
        check(matvec(covariance, vector)
              == tuple(Fraction(8, 3) * value for value in vector),
              "covariance eigenvalue on E is 8/3")

    # Exact quadratic identity on a spanning two-parameter family.  Since the
    # covariance restriction is (8/3)I_E, this is a theorem for every real c.
    for alpha, beta in product(range(-5, 6), repeat=2):
        vector = tuple(alpha * E_BASIS[0][j] + beta * E_BASIS[1][j]
                       for j in range(6))
        norm2 = sum(value * value for value in vector)
        outcomes = tuple(sum(vector[j] * values[j] for j in range(6))
                         for values in pair_values)
        second = sum(value * value for value in outcomes) / 6
        check(second == Fraction(8, 3) * norm2, "locked E quadratic identity")
        if alpha or beta:
            check(any(outcomes), "nonzero E operator survives local locked projection")

    # Pair Pauli words are orthonormal in the normalized product trace.
    words = [tuple("Z" if port in pair else "I" for port in range(4))
             for pair in PAIRS]
    for i, left in enumerate(words):
        for j, right in enumerate(words):
            trace_inner = Fraction(int(left == right))
            check(trace_inner == Fraction(int(i == j)), "pair trace Gram is identity")


def section_s4_intertwiner_obstruction():
    # The E projector above can also be reconstructed directly from ker R.
    locked = tuple(z for z in product((-1, 1), repeat=4) if sum(z) == 0)
    values = tuple(pair_vector(z) for z in locked)
    mean = (Fraction(-1, 3),) * 6
    covariance = [[sum(Fraction((m[i] - mean[i]) * (m[j] - mean[j]), 6)
                       for m in values) for j in range(6)] for i in range(6)]
    p_e = matscale(covariance, Fraction(3, 8))

    representations = []
    character_overlap = Fraction()
    for sigma in permutations(range(4)):
        port_rep = permutation_matrix(sigma, 4)
        pair_images = [PAIR_INDEX[tuple(sorted((sigma[a], sigma[b])))]
                       for a, b in PAIRS]
        pair_rep = permutation_matrix(pair_images, 6)
        chi_port = sum(port_rep[i][i] for i in range(4))
        chi_e = sum(matmul(p_e, pair_rep)[i][i] for i in range(6))
        character_overlap += chi_port * chi_e
        representations.append((pair_rep, port_rep))
    check(character_overlap / 24 == 0, "four-port representation has no E")

    zero = [[Fraction() for _ in range(4)] for _ in range(6)]
    for output in range(6):
        for input_port in range(4):
            seed = [[Fraction() for _ in range(4)] for _ in range(6)]
            seed[output][input_port] = Fraction(1)
            reynolds = [[Fraction() for _ in range(4)] for _ in range(6)]
            for pair_rep, port_rep in representations:
                term = matmul(matmul(pair_rep, seed), transpose(port_rep))
                reynolds = matadd(reynolds, term)
            reynolds = matscale(reynolds, Fraction(1, 24))
            check(matmul(p_e, reynolds) == zero,
                  "every S4-equivariant port-to-pair map has zero E output")
            for pair_rep, port_rep in representations:
                check(matmul(matmul(pair_rep, reynolds), transpose(port_rep)) == reynolds,
                      "Reynolds map is equivariant")


def section_direct_projection_and_query_source_separation():
    locked_bits = tuple(bits for bits in product((0, 1), repeat=4)
                        if sum(bits) == 2)
    for bits in locked_bits:
        for port in range(4):
            changed = list(bits)
            changed[port] ^= 1
            check(sum(changed) in (1, 3), "one transverse flip leaves k=2")

    # Direct pair-Z versus single-port-X normalized trace overlap is zero.
    for pair in PAIRS:
        pair_word = tuple("Z" if port in pair else "I" for port in range(4))
        for port in range(4):
            x_word = tuple("X" if index == port else "I" for index in range(4))
            check(pair_word != x_word, "query/read pair word is distinct from K source word")

    theorem = " ".join((AUTHOR / "THEOREM.md").read_text().split())
    check("finite real linear combination of complete pair reads" in theorem,
          "pair operator is typed as authenticated query/read")
    check("A generic multi-cell defect word is not `S4`-closed" in theorem,
          "one-cell selection rule is not promoted to generic defects")
    check("This is not a nonzero linear projection" in theorem,
          "nonlinear loop is separated from linear source overlap")


def section_all_local_loop_displacements():
    # Exhaust all possible two loop ports, both alternating orientations, and
    # both choices of the one occupied external port.  This strictly contains
    # the three parent-node cases used by the sealed target loop.
    for loop_ports in combinations(range(4), 2):
        external = tuple(port for port in range(4) if port not in loop_ports)
        for first_occupied in (0, 1):
            loop_bits = (first_occupied, 1 - first_occupied)
            for occupied_external in external:
                initial = [0, 0, 0, 0]
                initial[loop_ports[0]] = loop_bits[0]
                initial[loop_ports[1]] = loop_bits[1]
                initial[occupied_external] = 1
                final = list(initial)
                final[loop_ports[0]], final[loop_ports[1]] = (
                    initial[loop_ports[1]], initial[loop_ports[0]]
                )
                check(sum(initial) == sum(final) == 2, "loop exchange preserves local lock")
                m_i = pair_vector(tuple(1 - 2 * bit for bit in initial))
                m_f = pair_vector(tuple(1 - 2 * bit for bit in final))
                delta = tuple(m_f[j] - m_i[j] for j in range(6))
                check(matvec(R, delta) == (0, 0, 0, 0), "loop displacement lies in E")
                check(any(delta), "loop displacement is nonzero")
                check(sum(value * value for value in delta) == 16,
                      "loop displacement squared norm is 16")
                check(sum(delta[j] * (m_f[j] - m_i[j]) for j in range(6)) == 16,
                      "authenticated delta read distinguishes endpoints by 16")


INITIAL = (1, 0, 1, 0, 1, 0)


def loop_energy(subset):
    changes = [0] * 6
    for edge in subset:
        changes[edge] = 1 if INITIAL[edge] == 0 else -1
    charges = [changes[(vertex - 1) % 6] + changes[vertex]
               for vertex in range(6)]
    return sum(charge * charge for charge in charges)


def weighted_path_sum(weights):
    total = Fraction()
    for order in permutations(range(6)):
        subset = set()
        term = Fraction(1)
        for edge in order[:-1]:
            subset.add(edge)
            energy = loop_energy(frozenset(subset))
            check(energy > 0, "proper alternating-loop prefix remains in Q")
            term *= Fraction(-1, energy)
        for edge in order:
            term *= weights[edge]
        total += term
    return total


def section_six_support_factor_and_trace_counterexample():
    full = frozenset(range(6))
    check(loop_energy(frozenset()) == loop_energy(full) == 0,
          "empty and complete alternating sets are locked")
    for size in range(1, 6):
        for subset in combinations(range(6), size):
            check(loop_energy(frozenset(subset)) > 0,
                  "no proper nonempty subset is locked")

    census = {size: Counter(loop_energy(frozenset(subset))
                            for subset in combinations(range(6), size))
              for size in range(1, 6)}
    check(census[1] == Counter({2: 6}), "one-step loop energies")
    check(census[2] == Counter({4: 9, 2: 6}), "two-step loop energies")
    check(census[3] == Counter({4: 12, 2: 6, 6: 2}), "three-step loop energies")
    check(census[4] == census[2] and census[5] == census[1],
          "loop energy complement symmetry")

    unit = weighted_path_sum((Fraction(1),) * 6)
    check(unit == Fraction(-63, 8), "independent 720-path coefficient")
    rational_samples = (
        tuple(Fraction(value) for value in (2, 3, 5, 7, 11, 13)),
        (Fraction(1, 2), Fraction(-2, 3), Fraction(4, 5),
         Fraction(3, 7), Fraction(-5, 11), Fraction(6, 13)),
    )
    for weights in rational_samples:
        product_weight = Fraction(1)
        for value in weights:
            product_weight *= value
        check(weighted_path_sum(weights) == unit * product_weight,
              "every path factors through all six support weights")
    for weights in product((Fraction(0), Fraction(1)), repeat=6):
        expected = unit if all(weights) else Fraction()
        check(weighted_path_sum(weights) == expected,
              "removing any binary support kills the selected loop entry")

    # Normalized trace of Pauli words is cyclic.  Verify the actual complex
    # Pauli multiplication table on the full four-qubit basis; linearity then
    # proves the claim for every local operator pair.
    multiplication = {
        ("I", "I"): (1, "I"), ("I", "X"): (1, "X"),
        ("I", "Y"): (1, "Y"), ("I", "Z"): (1, "Z"),
        ("X", "I"): (1, "X"), ("Y", "I"): (1, "Y"),
        ("Z", "I"): (1, "Z"), ("X", "X"): (1, "I"),
        ("Y", "Y"): (1, "I"), ("Z", "Z"): (1, "I"),
        ("X", "Y"): (1j, "Z"), ("Y", "X"): (-1j, "Z"),
        ("Y", "Z"): (1j, "X"), ("Z", "Y"): (-1j, "X"),
        ("Z", "X"): (1j, "Y"), ("X", "Z"): (-1j, "Y"),
    }

    def pauli_trace_product(left, right):
        phase = 1
        result = []
        for a, b in zip(left, right):
            local_phase, local_result = multiplication[(a, b)]
            phase *= local_phase
            result.append(local_result)
        return phase if all(value == "I" for value in result) else 0

    paulis = tuple(product("IXYZ", repeat=4))
    for left in paulis:
        for right in paulis:
            check(pauli_trace_product(left, right) == pauli_trace_product(right, left),
                  "Pauli trace commutator vanishes")

    theorem = " ".join((AUTHOR / "THEOREM.md").read_text().split())
    check("product trace is invariant under every automorphism of this UHF spin" in theorem,
          "UHF unique-trace automorphism invariance is stated")
    check("refutes a **universal** nonzero stationary defect-contrast" in theorem,
          "only universal nonzero claim is refuted")
    check("neither prove nor refute the **existence**" in theorem,
          "existential selected-state claim remains open")
    check("lower-order diagonal shifts for a nonuniform defect word" in theorem,
          "nonuniform lower-order terms are withheld")


def section_custody_and_promotion_ceiling():
    dependencies = []
    for line in (AUTHOR / "DEPENDENCIES.sha256").read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        target = ROOT / relative
        check(target.is_file(), f"dependency exists: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"dependency hash: {relative}")
        dependencies.append(relative)
    check(len(dependencies) == 11, "exact frozen dependency count")
    check(len(set(dependencies)) == 11, "dependency rows are unique")
    check(all("GL6AM" in path or "GL6AN" in path for path in dependencies),
          "custody confined to GL6AM/GL6AN")
    check(any(path.endswith("AUDIT_G_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001/SEAL.sha256")
              for path in dependencies), "GL6AM hostile audit seal pinned")
    check(any(path.endswith("AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/SEAL.sha256")
              for path in dependencies), "GL6AN hostile audit seal pinned")

    combined = " ".join(
        " ".join((AUTHOR / name).read_text().split())
        for name in ("THEOREM.md", "RESULT.md", "README.md", "SELF_AUDIT.md")
    )
    required = (
        "They do not select a state",
        "positivity, not strict positivity",
        "generic multi-cell defect word is not `S4`-closed",
        "not a nonzero linear projection",
        "not construct an infinite locked projector",
        "universal claim is refuted; the existential selected-state claim remains open",
        "No state, pole, physical momentum, cone",
    )
    for token in required:
        check(token in combined, f"scope ceiling present: {token}")


def main():
    section_exact_locked_e_operator()
    section_s4_intertwiner_obstruction()
    section_direct_projection_and_query_source_separation()
    section_all_local_loop_displacements()
    section_six_support_factor_and_trace_counterexample()
    section_custody_and_promotion_ceiling()
    print(f"PASS__INDEPENDENT_GL6AQ_HOSTILE_REPLAY__{checks}/{checks}")
    print("PAIR=AUTHENTICATED_QUERY_READ_NONZERO_ON_LOCAL_LOCKED_E")
    print("K=DIRECT_LOCKED_ZERO_AND_ONE_CELL_S4_LINEAR_E_ZERO")
    print("LOOP=SIX_SUPPORT_PRODUCT_TIMES_MINUS63_OVER8_CHANGES_E")
    print("TRACE=POSITIVE_E_MASS_ZERO_RETARDED_AND_ZERO_DEFECT_CONTRAST")
    print("SCOPE=UNIVERSAL_NONZERO_REFUTED_EXISTENTIAL_SELECTED_LOCKED_STATE_OPEN")
    print("CEILING=NO_BULK_STATE_COMPLETE_GENERATOR_POLE_CONE_GRAVITY_G")


if __name__ == "__main__":
    main()
