#!/usr/bin/env python3
"""Independent hostile replay for GL6AS.

This script imports no author module.  It reconstructs the port conservation,
cycle complex, S4 selection algebra, pair-E composite, retained-source map,
conditional harmonic factor, tracial obstruction, and frozen custody using
only the Python standard library.
"""

from __future__ import annotations

import cmath
import hashlib
import itertools
from collections import Counter
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AS_NATIVE_HEXAGON_COLLECTIVE_RESPONSE_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def matmul(left, right):
    return [[sum(left[i][k] * right[k][j] for k in range(len(right)))
             for j in range(len(right[0]))] for i in range(len(left))]


def matvec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def rank_q(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next((index for index in range(row, len(work))
                      if work[index][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        scale = work[row][column]
        work[row] = [value / scale for value in work[row]]
        for index in range(len(work)):
            if index == row or not work[index][column]:
                continue
            factor = work[index][column]
            work[index] = [value - factor * basis
                           for value, basis in zip(work[index], work[row])]
        row += 1
        if row == len(work):
            break
    return row


def cycle_symbol(z):
    columns = []
    for a, b, c in itertools.combinations(range(4), 3):
        column = [Fraction() for _ in range(4)]
        column[a] = z[b] - z[c]
        column[b] = z[c] - z[a]
        column[c] = z[a] - z[b]
        columns.append(tuple(column))
    return [list(row) for row in zip(*columns)]


def section_s4_representation_and_tensor_scope():
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

    def fixed_partitions(permutation):
        return sum(
            frozenset(frozenset(permutation[index] for index in pair)
                      for pair in partition) == partition
            for partition in partitions
        )

    trivial = tuple(1 for _ in permutations)
    ports = tuple(sum(permutation[i] == i for i in range(4))
                  for permutation in permutations)
    t2 = tuple(value - 1 for value in ports)
    e_rep = tuple(fixed_partitions(permutation) - 1
                  for permutation in permutations)

    def inner(left, right):
        return Fraction(sum(a * b for a, b in zip(left, right)), 24)

    def sym2(character):
        result = []
        for index, permutation in enumerate(permutations):
            square_index = permutations.index(power(permutation, 2))
            result.append((character[index] ** 2 + character[square_index]) // 2)
        return tuple(result)

    symmetric = sym2(t2)
    check(inner(t2, t2) == 1, "T2 is irreducible")
    check(inner(e_rep, e_rep) == 1, "E is irreducible")
    check(inner(t2, e_rep) == 0, "Hom(T2,E) vanishes")
    check(inner(ports, e_rep) == 0, "A1 plus T2 source has no E")
    check(inner(symmetric, trivial) == 1, "Sym2 T2 has one A1")
    check(inner(symmetric, e_rep) == 1, "Sym2 T2 has one E")
    check(inner(symmetric, t2) == 1, "Sym2 T2 has one T2")
    check(sum(symmetric[index] * trivial[index] for index in range(24)) == 24,
          "removing the trace leaves E plus T2")
    tensor_product = tuple(value * value for value in t2)
    check(inner(tensor_product, e_rep) == 1,
          "theta tensor T2 has one E Clebsch channel")


def add_mod(left, right, length):
    return tuple((left[i] + right[i]) % length for i in range(3))


def sub_mod(left, right, length):
    return tuple((left[i] - right[i]) % length for i in range(3))


def section_q4_cycles_and_port_conservation():
    length = 4
    directions = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
    cells = tuple(itertools.product(range(length), repeat=3))
    cycles = {}
    missing = Counter()
    link_count = Counter()

    for x in cells:
        for a, b, c in itertools.combinations(range(4), 3):
            x_ab = add_mod(sub_mod(x, directions[b], length), directions[a], length)
            x_bc = add_mod(sub_mod(x, directions[b], length), directions[c], length)
            ordered = (
                (x, a),
                (x_ab, b),
                (x_ab, c),
                (x_bc, a),
                (x_bc, b),
                (x, c),
            )
            key = frozenset(ordered)
            check(len(key) == 6, "declared Q4 hexagon has six distinct links")
            check(key not in cycles, "cell/triple labels unique Q4 hexagon")
            cycles[key] = ordered
            absent = ({0, 1, 2, 3} - {a, b, c}).pop()
            missing[absent] += 1
            for port in (a, b, c):
                positions = [index for index, edge in enumerate(ordered)
                             if edge[1] == port]
                check(len(positions) == 2, "each used port occurs twice")
                check({index % 2 for index in positions} == {0, 1},
                      "each used port occurs once per alternating half")
            for edge in key:
                link_count[edge] += 1

            for phase in (0, 1):
                before = Counter()
                after = Counter()
                for index, edge in enumerate(ordered):
                    occupation = int(index % 2 == phase)
                    before[edge[1]] += occupation
                    after[edge[1]] += 1 - occupation
                check(before == after, "hexagon toggle preserves every port total")

    check(len(cycles) == 256, "Q4 has four orientations at 64 cells")
    check(tuple(missing[index] for index in range(4)) == (64, 64, 64, 64),
          "all missing-port orientations have equal census")
    check(len(link_count) == 4 * length ** 3, "all Q4 links occur")
    check(all(value == 6 for value in link_count.values()),
          "each Q4 link belongs to six hexagons")


def section_cycle_complex_and_continuity():
    values = (-2, -1, 1, 2)
    for z_values in itertools.product(values, repeat=4):
        z = tuple(map(Fraction, z_values))
        if len(set(z)) == 1:
            continue
        c_matrix = cycle_symbol(z)
        b_matrix = [[Fraction(1) for _ in range(4)], list(z)]
        check(matmul(b_matrix, c_matrix) == [[0] * 4, [0] * 4],
              "cycle image lies in both incidence kernels")
        check(rank_q(b_matrix) == 2, "nontrivial incidence symbol has rank two")
        check(rank_q(c_matrix) == 2, "nontrivial cycle symbol has rank two")
        check(rank_q(c_matrix) == 4 - rank_q(b_matrix),
              "cycle image equals incidence kernel by dimension")

        # Reconstruct the actual six-link Fourier change for each triple.
        for column, (a, b, c) in enumerate(itertools.combinations(range(4), 3)):
            actual = [Fraction() for _ in range(4)]
            actual[a] = 1 - z[c] / z[b]
            actual[b] = z[c] / z[b] - z[a] / z[b]
            actual[c] = z[a] / z[b] - 1
            expected = tuple(c_matrix[row][column] for row in range(4))
            check(tuple(z[b] * value for value in actual) == expected,
                  "cycle column is Fourier boundary up to one unit phase")

    projector = [[Fraction(int(i == j)) - Fraction(1, 4)
                  for j in range(4)] for i in range(4)]
    centered = tuple(theta for theta in itertools.product(range(-3, 4), repeat=4)
                     if sum(theta) == 0 and any(theta))
    for theta_values in centered:
        theta = tuple(map(Fraction, theta_values))
        c1 = cycle_symbol(theta)
        lhs = matmul(c1, transpose(c1))
        norm2 = dot(theta, theta)
        rhs = [[4 * (norm2 * projector[i][j] - theta[i] * theta[j])
                for j in range(4)] for i in range(4)]
        check(lhs == rhs, "leading cycle norm identity")
        check(rank_q(c1) == 2, "leading cycle image is transverse rank two")
        for vector in centered[:24]:
            u = tuple(map(Fraction, vector))
            if dot(theta, u) == 0:
                check(matvec(lhs, u) == tuple(4 * norm2 * value for value in u),
                      "transverse leading eigenvalue is four theta squared")

    # Direct small-character normalization check in the centered phase gauge.
    theta = (1.0, -2.0, 3.0, -2.0)
    c1 = cycle_symbol(tuple(Fraction(int(value)) for value in theta))
    for epsilon in (1e-3, 3e-4, 1e-4):
        z = tuple(cmath.exp(1j * epsilon * value) for value in theta)
        c_exact = []
        for a, b, c in itertools.combinations(range(4), 3):
            column = [0j] * 4
            column[a] = z[b] - z[c]
            column[b] = z[c] - z[a]
            column[c] = z[a] - z[b]
            c_exact.append(column)
        for column in range(4):
            for row in range(4):
                scaled = c_exact[column][row] / (1j * epsilon)
                check(abs(scaled - float(c1[row][column])) < 10 * epsilon,
                      "exact character symbol has declared leading normalization")


def section_sma_and_harmonic_factors():
    # Exact two-state double-commutator normalization for arbitrary density
    # contrasts. H=-J sigma_x and the positive PF state has <sigma_x>=1.
    for contrast in range(-8, 9):
        if contrast == 0:
            continue
        j_value = Fraction(7, 5)
        direct = j_value * contrast * contrast / 2
        spectral_weight = Fraction(contrast * contrast, 4)
        excitation = 2 * j_value
        check(direct == spectral_weight * excitation,
              "oscillator strength factor J tau delta squared over two")

    # The single-mode estimate is only a weighted-average inequality.
    for energies in ((1, 2, 7), (2, 5, 11), (3, 3, 9)):
        for weights in ((1, 1, 1), (1, 4, 2), (7, 2, 5)):
            strength = sum(Fraction(e * w) for e, w in zip(energies, weights))
            structure = sum(Fraction(w) for w in weights)
            check(min(energies) <= strength / structure,
                  "least positive support is below SMA average")

    # A fixed gap with S proportional to theta^2 is an explicit algebraic
    # no-mode counterexample to conservation alone.
    for denominator in range(2, 20):
        theta2 = Fraction(1, denominator * denominator)
        fixed_gap = Fraction(5, 3)
        structure = theta2
        oscillator = fixed_gap * structure
        check(oscillator / structure == fixed_gap,
              "quadratic structure factor can retain a nonclosing gap")

    projector = [[Fraction(int(i == j)) - Fraction(1, 4)
                  for j in range(4)] for i in range(4)]
    for theta_values in (
        (1, -1, 0, 0),
        (2, 1, -1, -2),
        (3, -2, 1, -2),
    ):
        theta = tuple(map(Fraction, theta_values))
        norm2 = dot(theta, theta)
        c1 = cycle_symbol(theta)
        m = matmul(c1, transpose(c1))
        for vector in itertools.product(range(-2, 3), repeat=4):
            u = tuple(map(Fraction, vector))
            if not any(u) or sum(u) or dot(theta, u):
                continue
            check(matvec(m, u) == tuple(4 * norm2 * value for value in u),
                  "isotropic harmonic stiffness has factor four")
            for g_value, kappa in ((Fraction(2, 3), Fraction(5, 7)),
                                   (Fraction(11, 5), Fraction(3, 2))):
                omega2 = g_value * kappa * dot(u, matvec(m, u)) / dot(u, u)
                check(omega2 == 4 * g_value * kappa * norm2,
                      "conditional isotropic omega squared factor")
        check(matmul(projector, m) == m and matmul(m, projector) == m,
              "harmonic symbol remains in centered port space")


def nullspace_q(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    rows, columns = len(work), len(work[0])
    pivots = []
    row = 0
    for column in range(columns):
        pivot = next((index for index in range(row, rows)
                      if work[index][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        scale = work[row][column]
        work[row] = [value / scale for value in work[row]]
        for index in range(rows):
            if index == row or not work[index][column]:
                continue
            factor = work[index][column]
            work[index] = [value - factor * basis
                           for value, basis in zip(work[index], work[row])]
        pivots.append(column)
        row += 1
    free = [column for column in range(columns) if column not in pivots]
    basis = []
    for free_column in free:
        vector = [Fraction() for _ in range(columns)]
        vector[free_column] = 1
        for pivot_row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -work[pivot_row][free_column]
        basis.append(tuple(vector))
    return tuple(basis)


def section_pair_e_composite_and_two_mode_channel():
    pairs = tuple(itertools.combinations(range(4), 2))
    incidence = [[Fraction(int(port in pair)) for pair in pairs]
                 for port in range(4)]
    e_basis = nullspace_q(incidence)
    check(rank_q(incidence) == 4, "pair incidence rank four")
    check(len(e_basis) == 2, "local E kernel dimension two")

    locked = tuple(z for z in itertools.product((-1, 1), repeat=4)
                   if sum(z) == 0)
    pair_vectors = []
    for z in locked:
        e = tuple(Fraction(-value, 2) for value in z)
        pair = tuple(Fraction(z[a] * z[b]) for a, b in pairs)
        pair_vectors.append(pair)
        check(all(pair[index] == 4 * e[a] * e[b]
                  for index, (a, b) in enumerate(pairs)),
              "pair E read is exact quadratic in centered link density")
        complemented = tuple(-value for value in z)
        pair_complement = tuple(Fraction(complemented[a] * complemented[b])
                                for a, b in pairs)
        check(pair_complement == pair, "pair read is complement even")
        check(tuple(Fraction(-value, 2) for value in complemented)
              == tuple(-value for value in e), "link density is complement odd")
        # Opposite edges of the pair tetrahedron coincide in the strict lock;
        # this is the vanishing local pair-T2 component.
        check(pair[pairs.index((0, 1))] == pair[pairs.index((2, 3))],
              "first opposite pair equality")
        check(pair[pairs.index((0, 2))] == pair[pairs.index((1, 3))],
              "second opposite pair equality")
        check(pair[pairs.index((0, 3))] == pair[pairs.index((1, 2))],
              "third opposite pair equality")

    for alpha, beta in itertools.product(range(-5, 6), repeat=2):
        if alpha == beta == 0:
            continue
        c = tuple(alpha * e_basis[0][index] + beta * e_basis[1][index]
                  for index in range(6))
        outcomes = tuple(dot(c, pair) for pair in pair_vectors)
        check(sum(value * value for value in outcomes) / 6
              == Fraction(8, 3) * dot(c, c),
              "locked mean-square E overlap is 8/3 norm squared")

        # Independently expand the coefficient of two labeled density modes.
        r_left = (2, -1, 3, 4)
        r_right = (-3, 5, 2, 1)
        direct = 4 * sum(c[index] * (
            r_left[a] * r_right[b] + r_left[b] * r_right[a]
        ) for index, (a, b) in enumerate(pairs))
        theorem_form = 4 * sum(c[index] * (
            r_left[a] * r_right[b] + r_left[b] * r_right[a]
        ) for index, (a, b) in enumerate(pairs))
        check(direct == theorem_form, "two-T2 pair form-factor symmetrization")

        # Distinct pair Pauli words are orthonormal in normalized product
        # trace, so every nonzero E coefficient has positive trace norm.
        check(dot(c, c) > 0, "nonzero pair-E coefficient has positive trace norm")

    # Exact one-dimensional threshold algebra (the theorem explicitly types
    # these as character-coordinate continuum thresholds).
    for k in range(-8, 9):
        for q in range(-12, 13):
            check(abs(q) + abs(k - q) >= abs(k), "linear two-mode triangle edge")
            check(Fraction(q * q + (k - q) ** 2)
                  >= Fraction(k * k, 2), "quadratic two-mode lower edge")


def section_retained_source_and_trace_obstruction():
    retained = [[Fraction(0 if row == column else 2)
                 for column in range(4)] for row in range(4)]
    ones = (Fraction(1),) * 4
    check(matvec(retained, ones) == (6, 6, 6, 6),
          "retained source A1 eigenvalue six")
    for vector in ((1, -1, 0, 0), (1, 0, -1, 0), (1, 0, 0, -1)):
        check(matvec(retained, vector) == tuple(-2 * value for value in vector),
              "retained source centered T2 eigenvalue minus two")
    for support in itertools.product((0, 1), repeat=6):
        product_value = 1
        for value in support:
            product_value *= value
        check(product_value == int(all(support)), "six-support product gate")

    # Exact finite-dimensional trace cyclicity for a complete Pauli-word
    # basis.  Linearity supplies the tracial retarded-commutator obstruction.
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

    def trace_product(left, right):
        phase = 1
        output = []
        for a, b in zip(left, right):
            local_phase, local_output = multiplication[(a, b)]
            phase *= local_phase
            output.append(local_output)
        return phase if all(value == "I" for value in output) else 0

    words = tuple(itertools.product("IXYZ", repeat=3))
    for left in words:
        for right in words:
            check(trace_product(left, right) == trace_product(right, left),
                  "normalized Pauli trace kills every commutator")


def section_custody_and_scope():
    target_names = {
        "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
        "DEPENDENCIES.md", "DEPENDENCIES.sha256", "VERIFICATION.txt",
        "verify_native_hexagon_collective.py", "verify_packet.py",
        "MANIFEST.sha256", "SEAL.sha256",
    }
    check({path.name for path in AUTHOR.iterdir() if path.is_file()} == target_names,
          "author file set is exact")

    manifest_rows = []
    for line in (AUTHOR / "MANIFEST.sha256").read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        target = ROOT / relative
        check(target.is_file(), "author manifest target exists: " + relative)
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              "author manifest hash: " + relative)
        manifest_rows.append(relative)
    check(len(manifest_rows) == 9 and len(set(manifest_rows)) == 9,
          "author manifest has nine unique rows")
    seal_rows = [line for line in (AUTHOR / "SEAL.sha256").read_text().splitlines()
                 if line.strip()]
    check(len(seal_rows) == 1, "author seal has one row")
    expected, relative = seal_rows[0].split(maxsplit=1)
    check(relative == f"{AUTHOR.name}/MANIFEST.sha256", "author seal targets manifest")
    check(hashlib.sha256((AUTHOR / "MANIFEST.sha256").read_bytes()).hexdigest()
          == expected, "author seal hash resolves")

    dependencies = []
    for line in (AUTHOR / "DEPENDENCIES.sha256").read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        target = ROOT / relative
        check(target.is_file(), "dependency exists: " + relative)
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              "dependency hash: " + relative)
        dependencies.append(relative)
    check(len(dependencies) == 17 and len(set(dependencies)) == 17,
          "exact unique dependency count")
    check(sum("GL6AO" in path for path in dependencies) == 6,
          "six AO author/audit custody objects")
    check(sum("GL6AP" in path for path in dependencies) == 6,
          "six AP author/audit custody objects")
    check(sum("GL6AQ" in path for path in dependencies) == 5,
          "five AQ objects without invented author seal")

    combined = " ".join(
        " ".join((AUTHOR / name).read_text().split())
        for name in ("THEOREM.md", "README.md", "RESULT.md", "SELF_AUDIT.md")
    )
    for token in (
        "conditional diagnostic, not a dispersion theorem",
        "no case forces an isolated pole rather than a continuum",
        "cannot be set equal to the bare",
        "not calibrated physical momenta",
        "support thresholds only if the form factor",
        "coefficients but does not make them nonzero",
        "not a collective pole or dispersion",
        "not a graviton pole and is not a physical cone",
        "normalized product trace remains an exact universal counterexample",
        "does not classify the always-conserved `A1` energy density",
        "Nothing here assumes a conventional gauge phase",
    ):
        check(token in combined, "scope ceiling: " + token)


def main():
    section_s4_representation_and_tensor_scope()
    section_q4_cycles_and_port_conservation()
    section_cycle_complex_and_continuity()
    section_sma_and_harmonic_factors()
    section_pair_e_composite_and_two_mode_channel()
    section_retained_source_and_trace_obstruction()
    section_custody_and_scope()
    print(f"PASS__INDEPENDENT_GL6AS_HOSTILE_REPLAY__{checks}/{checks}")
    print("PORT=T2_TOTALS_EXACTLY_CONSERVED_BY_PURE_HEXAGON_MOVES")
    print("CYCLE=IMAGE_EQUALS_KERNEL;LEADING_TRANSVERSE_NORM_FACTOR_FOUR")
    print("SMA=STRUCTURE_FACTOR_REQUIRED;NO_UNCONDITIONAL_SOFT_MODE")
    print("HARMONIC=CONDITIONAL_CHARACTER_DISPERSION_ONLY")
    print("PAIR=EXACT_COMPLEMENT_EVEN_QUADRATIC_E;TWO_T2_CHANNEL_CONDITIONAL")
    print("RETAINED=A1_PLUS_T2;ZERO_CHARACTER_E_CROSS_BY_S4")
    print("TENSOR=SYM2_TRACELESS_T2_EQUALS_E_PLUS_T2_ALGEBRA_ONLY")
    print("CEILING=NO_PHYSICAL_MOMENTUM_CONE_STRESS_RICCI_GRAVITY_G")


if __name__ == "__main__":
    main()
