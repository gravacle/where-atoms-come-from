#!/usr/bin/env python3
"""Independent hostile replay for frozen GL6BC V001.

Standard library only.  This program imports no author module and uses no
Python ``assert`` statements, so normal and optimized runs execute the same
checks.
"""

from __future__ import annotations

import cmath
import hashlib
import itertools
import math
import random
import re
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6BC_FINITE_COLLAR_RECIPROCAL_LINEAGE_OBSTRUCTION_V001"
CHECKS = 0
TOL = 3.0e-11


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise RuntimeError(f"FAIL [{CHECKS + 1}] {label}")
    CHECKS += 1


def close(left: complex | float, right: complex | float, label: str,
          tolerance: float = TOL) -> None:
    scale = max(1.0, abs(left), abs(right))
    check(abs(left - right) <= tolerance * scale, label)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


TARGET_HASHES = {
    "DEPENDENCIES.md": "688198d9b7f151ba9f4aafe24b6b3590869156ce45ed1e12be9e47332275a08a",
    "DEPENDENCIES.sha256": "45dd8b87edbebe8f3be9930377722001b2afbdd0c7c7e16a48602eb6d63f1ea5",
    "MANIFEST.sha256": "614570080fd9ce3ebac3edfef4c74ede52766c96948659eba643e2b3b286ab5a",
    "README.md": "edd5c3b2b963c1ecc723bbbbcdb9ad7e0c76803b65c263cdc0104f82c3284819",
    "RECIPROCITY_LEDGER.json": "1bd68d5c739bbf56ce24350860f5fccddf993c91896fb7c3293999afdca6f053",
    "RESULT.md": "d9f57f46007057d8d6ea99e9738bb3b995d84bf607d8817b5a8873a8ab92b05e",
    "SEAL.sha256": "c6ccb51c0796ef889b471fde9c691407cdda16d23239deb917cc4d6b0405cd03",
    "SELF_AUDIT.md": "7282e51a4555627022efdc4cd7b705aaa72c0a3e6c3c94f3dac83edc9c968e88",
    "THEOREM.md": "6355e6e1dda470e363122f5e3342c01346dc81e8992d05b1436f30b899041ea6",
    "VERIFICATION.txt": "f4f4ff69f2902c9dbe5e906a580acc4ebabee8c58f48a34037b33dd534fb21d5",
    "verify_finite_collar_reciprocal_lineage_obstruction.py": "fa619989713935f339a443d4ed21a983db2793028f4bfb6092bb62d41a56afa7",
    "verify_packet.py": "1ca37c32825463389a0b87a1d818f7cf4cd58073d46b8f946c1237ad7f7a44ce",
}


DEPENDENCY_HASHES = {
    "GRAVITY_RECORD_FIRST_WORKING_THEORY_CLOSURE_V001.md": "cf9229586268f054b473b1641085ebafc3bca01fa0691a191cbda923ae1fa7f2",
    "GRAVITY_RGRL_LOCAL_RESPONSE_G_IDENTIFIABILITY_CEILING_V002.md": "140eb379f22c369b0442b5585a4828122ac4f5858b54ae3f2bf6391866ac84a3",
    "GRAVITY_RGRL_ONSHELL_OFFSHELL_CLARIFICATION_ADOPTION_V001.md": "4959f99898b216edc7da3e212ce2e26422287899fcf8f3b41cd34ef5d8bb3ff8",
    "AUDIT_G_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/AUDIT.md": "03bda2dba369211542dfef1af065490e21033483443cd6a67a21b06bf74e0bc9",
    "AUDIT_G_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001/AUDIT.md": "5734ad57122c64e3174aa7706b0e7aa86102b3a18a3b868aca20af0997ab462a",
    "AUDIT_G_GL6T_F3_LINEAGE_GATED_Q4_E2_RESPONSE_V001/AUDIT.md": "cefa812da0cc0b75fd689097326f74f30d679571ff5001d0fa00d09632f944b1",
    "AUDIT_G_GL6V_Q4_PAIR_PHASE_SOURCE_READ_DILATION_V001/AUDIT.md": "1c217195d61cfb503f49fb4950d7cb4e8da9873e4214d0e64841d4cf3d284bbd",
    "LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/THEOREM.md": "d7ce0a7527a68f49e6ea2ee8edbb400a142fbb49297d8fe99cae78ffa0154ab0",
    "LANE_CROSS_RFT_GRA_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001/THEOREM.md": "a51e802f6ba148e5f9848e95f41a80073795b24b7eaf87e36c0766b0856aa494",
    "LANE_CROSS_RFT_GRA_GL6T_F3_LINEAGE_GATED_Q4_E2_RESPONSE_V001/THEOREM.md": "c610f0dc4092e8ca738c1b1a32a5e034f1ccbaafb032d3c1b5afdf3efaca539b",
    "LANE_CROSS_RFT_GRA_GL6V_Q4_PAIR_PHASE_SOURCE_READ_DILATION_V001/THEOREM.md": "3f5b52aa066d4d6f56f75a06a1f6623d49d988c531d0b9ef82e590dd92aec51d",
    "LANE_GRA_DQ_F3_RECIPROCAL_BACKPROPAGATION_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md": "7b2758e4254d95501bc3658ff8ab44c6f18b35135ab5c251ee4975d8b2f88920",
    "LANE_GRA_DQ_F3_RECIPROCAL_BACKPROPAGATION_SCREEN_V001/THEOREM.md": "7fcf779ae58e7f80d03a1f664a4371380134ae870e98b6a38407b9adf21ecc13",
}


def parse_hashes(path: Path) -> dict[str, str]:
    rows: dict[str, str] = {}
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        pieces = raw.split("  ", 1)
        check(len(pieces) == 2, f"two-column hash row {path.name}:{line_number}")
        digest, relative = pieces
        check(bool(re.fullmatch(r"[0-9a-f]{64}", digest)),
              f"valid SHA-256 {path.name}:{line_number}")
        check(relative not in rows, f"unique hash target {path.name}:{line_number}")
        rows[relative] = digest
    return rows


def verify_frozen_custody() -> None:
    check(AUTHOR.is_dir() and not AUTHOR.is_symlink(), "frozen author directory exists")
    actual_names = {path.name for path in AUTHOR.iterdir() if path.is_file()}
    check(actual_names == set(TARGET_HASHES), "exact frozen author file inventory")
    for name, expected in sorted(TARGET_HASHES.items()):
        path = AUTHOR / name
        check(path.is_file() and not path.is_symlink(), f"regular author target {name}")
        check(sha256(path) == expected, f"frozen author digest {name}")

    manifest = parse_hashes(AUTHOR / "MANIFEST.sha256")
    expected_manifest = set(TARGET_HASHES) - {"MANIFEST.sha256", "SEAL.sha256"}
    check(set(manifest) == expected_manifest, "author manifest has exact pre-seal inventory")
    for relative, expected in sorted(manifest.items()):
        check(sha256(AUTHOR / relative) == expected,
              f"author manifest replay {relative}")
    seal = parse_hashes(AUTHOR / "SEAL.sha256")
    check(seal == {"MANIFEST.sha256": TARGET_HASHES["MANIFEST.sha256"]},
          "author seal pins the exact manifest")

    ledger = parse_hashes(AUTHOR / "DEPENDENCIES.sha256")
    check(ledger == DEPENDENCY_HASHES, "dependency ledger equals independently frozen map")
    for relative, expected in sorted(DEPENDENCY_HASHES.items()):
        path = ROOT / relative
        check(path.is_file() and not path.is_symlink(), f"regular dependency {relative}")
        check(sha256(path) == expected, f"dependency digest {relative}")


Matrix = list[list[complex]]


def zeros(rows: int, columns: int) -> Matrix:
    return [[0j for _ in range(columns)] for _ in range(rows)]


def eye(size: int) -> Matrix:
    output = zeros(size, size)
    for index in range(size):
        output[index][index] = 1
    return output


def add(*matrices: Matrix) -> Matrix:
    output = zeros(len(matrices[0]), len(matrices[0][0]))
    for matrix in matrices:
        for row in range(len(matrix)):
            for column in range(len(matrix[0])):
                output[row][column] += matrix[row][column]
    return output


def scale(coefficient: complex, matrix: Matrix) -> Matrix:
    return [[coefficient * value for value in row] for row in matrix]


def mul(left: Matrix, right: Matrix) -> Matrix:
    output = zeros(len(left), len(right[0]))
    for row in range(len(left)):
        for middle in range(len(right)):
            value = left[row][middle]
            if value == 0:
                continue
            for column in range(len(right[0])):
                output[row][column] += value * right[middle][column]
    return output


def dagger(matrix: Matrix) -> Matrix:
    return [[matrix[column][row].conjugate()
             for column in range(len(matrix))]
            for row in range(len(matrix[0]))]


def kron(left: Matrix, right: Matrix) -> Matrix:
    output = zeros(len(left) * len(right), len(left[0]) * len(right[0]))
    for i, row in enumerate(left):
        for j, value in enumerate(row):
            if value == 0:
                continue
            for k, right_row in enumerate(right):
                for ell, right_value in enumerate(right_row):
                    output[i * len(right) + k][j * len(right[0]) + ell] = value * right_value
    return output


def tensor(*matrices: Matrix) -> Matrix:
    output = [[1 + 0j]]
    for matrix in matrices:
        output = kron(output, matrix)
    return output


def comm(left: Matrix, right: Matrix) -> Matrix:
    return add(mul(left, right), scale(-1, mul(right, left)))


def max_abs(matrix: Matrix) -> float:
    return max(abs(value) for row in matrix for value in row)


def mat_vec(matrix: Matrix, vector: list[complex]) -> list[complex]:
    return [sum(value * vector[column] for column, value in enumerate(row))
            for row in matrix]


def inner(left: list[complex], right: list[complex]) -> complex:
    return sum(a.conjugate() * b for a, b in zip(left, right))


def expectation(vector: list[complex], matrix: Matrix) -> complex:
    return inner(vector, mat_vec(matrix, vector))


def single(operator: Matrix, slot: int, factors: int) -> Matrix:
    identity = [[1, 0], [0, 1]]
    return tensor(*(operator if index == slot else identity for index in range(factors)))


def permutation_unitary(factors: int, mapping) -> Matrix:
    dimension = 1 << factors
    output = zeros(dimension, dimension)
    for source in range(dimension):
        bits = [(source >> (factors - 1 - slot)) & 1 for slot in range(factors)]
        target_bits = mapping(bits[:])
        target = 0
        for bit in target_bits:
            target = (target << 1) | bit
        output[target][source] = 1
    return output


def cnot(factors: int, control: int, target: int) -> Matrix:
    def mapping(bits: list[int]) -> list[int]:
        if bits[control] == 1:
            bits[target] ^= 1
        return bits
    return permutation_unitary(factors, mapping)


def swap(factors: int, first: int, second: int) -> Matrix:
    def mapping(bits: list[int]) -> list[int]:
        bits[first], bits[second] = bits[second], bits[first]
        return bits
    return permutation_unitary(factors, mapping)


def normalized_random_state(dimension: int, generator: random.Random) -> list[complex]:
    vector = [complex(generator.uniform(-1, 1), generator.uniform(-1, 1))
              for _ in range(dimension)]
    norm = math.sqrt(inner(vector, vector).real)
    return [value / norm for value in vector]


def verify_support_word_algebra_and_schedule() -> None:
    identity = [[1, 0], [0, 1]]
    x = [[0, 1], [1, 0]]
    y = [[0, -1j], [1j, 0]]
    z = [[1, 0], [0, -1]]
    p = [[0, 0], [0, 1]]
    q = [[1, 0], [0, 0]]
    factors = 5  # K1,K2,link1,link2,source-query
    dimension = 1 << factors
    identity_full = eye(dimension)

    p1 = single(p, 0, factors)
    p2 = single(p, 1, factors)
    q1 = single(q, 0, factors)
    q2 = single(q, 1, factors)
    xk1 = single(x, 0, factors)
    x1 = single(x, 2, factors)
    x2 = single(x, 3, factors)
    n1 = single(p, 2, factors)
    n2 = single(p, 3, factors)
    z1 = single(z, 2, factors)
    z2 = single(z, 3, factors)
    z_query = single(z, 4, factors)
    pair = mul(z1, z2)

    projectors: dict[tuple[int, int], Matrix] = {}
    for beta1, beta2 in itertools.product((0, 1), repeat=2):
        left = p1 if beta1 else q1
        right = p2 if beta2 else q2
        projector = mul(left, right)
        projectors[(beta1, beta2)] = projector
        close(max_abs(add(mul(projector, projector), scale(-1, projector))), 0,
              f"Pi_{beta1}{beta2} is idempotent")
        close(max_abs(add(projector, scale(-1, dagger(projector)))), 0,
              f"Pi_{beta1}{beta2} is Hermitian")
    close(max_abs(add(*projectors.values(), scale(-1, identity_full))), 0,
          "route-word projectors resolve the identity")
    for left_word, right_word in itertools.combinations(projectors, 2):
        close(max_abs(mul(projectors[left_word], projectors[right_word])), 0,
              "distinct route-word projectors are orthogonal")

    ratio = Fraction(5, 2)
    h_nd = scale(Fraction(7, 19), z_query)
    h_full = add(
        scale(-1, mul(p1, x1)),
        scale(-1, mul(p2, x2)),
        scale(-6 * ratio, add(n1, n2)),
        scale(2 * ratio, mul(n1, n2)),
        h_nd,
    )
    for word, projector in projectors.items():
        close(max_abs(comm(h_full, projector)), 0,
              f"route-lifted Hamiltonian conserves Pi_{word}")
    close(max_abs(comm(mul(p1, x1), p1)), 0,
          "existing gated transverse term does not write K1")
    close(max_abs(comm(mul(p2, x2), p2)), 0,
          "existing gated transverse term does not write K2")
    close(max_abs(comm(h_nd, projectors[(1, 0)])), 0,
          "disjoint spectator Hamiltonian preserves route words")

    # GL6V's literal pair-copy/phase/uncompute schedule, including its source
    # ancilla, acts as identity on both K factors.
    copy = mul(cnot(factors, 3, 4), cnot(factors, 2, 4))
    close(max_abs(add(mul(copy, copy), scale(-1, identity_full))), 0,
          "pair copy is involutive")
    angle = 0.713
    phase = add(scale(math.cos(angle / 2), identity_full),
                scale(1j * math.sin(angle / 2), z_query))
    source = mul(copy, mul(phase, copy))
    for word, projector in projectors.items():
        close(max_abs(comm(copy, projector)), 0,
              f"pair copy preserves Pi_{word}")
        close(max_abs(comm(source, projector)), 0,
              f"phase-source sandwich preserves Pi_{word}")
    close(max_abs(comm(source, pair)), 0,
          "phase-source sandwich preserves scored pair")

    read_projectors: list[Matrix] = []
    for bit1, bit2 in itertools.product((0, 1), repeat=2):
        read_projector = mul(n1 if bit1 else add(identity_full, scale(-1, n1)),
                             n2 if bit2 else add(identity_full, scale(-1, n2)))
        read_projectors.append(read_projector)
        for word, projector in projectors.items():
            close(max_abs(comm(read_projector, projector)), 0,
                  f"complete link-read outcome preserves Pi_{word}")
    close(max_abs(add(*read_projectors, scale(-1, identity_full))), 0,
          "complete link read has no postselection gap")

    # A nontrivial piecewise propagator containing controlled transverse
    # rotations, diagonal interaction/spectator evolution, and a GL6V pulse.
    def controlled_rotation(projector: Matrix, link_x: Matrix, rotation: float) -> Matrix:
        correction = add(
            scale(math.cos(rotation) - 1, identity_full),
            scale(-1j * math.sin(rotation), link_x),
        )
        return add(identity_full, mul(projector, correction))

    control1 = controlled_rotation(p1, x1, 0.319)
    control2 = controlled_rotation(p2, x2, -0.227)
    pair_phase = add(identity_full,
                     scale(cmath.exp(-0.173j) - 1, mul(n1, n2)))
    spectator_phase = add(scale(math.cos(0.109), identity_full),
                          scale(-1j * math.sin(0.109), z_query))
    propagator = mul(source, mul(spectator_phase,
                     mul(pair_phase, mul(control2, control1))))
    close(max_abs(add(mul(dagger(propagator), propagator), scale(-1, identity_full))), 0,
          "piecewise current-parent propagator is unitary")
    for word, projector in projectors.items():
        close(max_abs(comm(propagator, projector)), 0,
              f"complete current schedule preserves Pi_{word}")

    # Exhaust a basis of the allowed algebra: I/Z on retained K and arbitrary
    # Pauli operators on link/link/source-read factors.
    paulis = (identity, x, y, z)
    allowed_basis: list[Matrix] = []
    for k_left, k_right in itertools.product((identity, z), repeat=2):
        for tail in itertools.product(paulis, repeat=3):
            observable = tensor(k_left, k_right, *tail)
            allowed_basis.append(observable)
            for word, projector in projectors.items():
                close(max_abs(comm(observable, projector)), 0,
                      f"allowed source/read basis commutes with Pi_{word}")
    check(len(allowed_basis) == 256, "allowed algebra basis census")

    # Both retarded orientations remain zero after nontrivial evolution.
    for index in range(0, len(allowed_basis), 13):
        observable = allowed_basis[index]
        evolved = mul(dagger(propagator), mul(observable, propagator))
        for word, projector in projectors.items():
            close(max_abs(comm(projector, observable)), 0,
                  "Pi(s)-then-B retarded commutator is zero")
            close(max_abs(comm(evolved, projector)), 0,
                  "B(s)-then-Pi retarded commutator is zero")

    generator = random.Random(60901)
    for sample in range(64):
        state = normalized_random_state(dimension, generator)
        evolved = mat_vec(propagator, state)
        for word, projector in projectors.items():
            before = expectation(state, projector)
            after = expectation(evolved, projector)
            close(after, before, f"entangled route probability conserved sample {sample} word {word}")

    # In all-BREAK, every link term and every allowed source/read operation
    # commutes with M.  This is state-independent, not a blank-state accident.
    break_projector = projectors[(0, 0)]
    h_break = mul(break_projector, mul(h_full, break_projector))
    close(max_abs(comm(h_break, pair)), 0, "all-BREAK full Hamiltonian leaves M static")
    close(max_abs(comm(copy, pair)), 0, "all-BREAK copy leaves M static")
    close(max_abs(comm(phase, pair)), 0, "all-BREAK phase leaves M static")
    close(max_abs(comm(source, pair)), 0, "all-BREAK uncompute leaves M static")
    for sample in range(32):
        state = normalized_random_state(dimension, generator)
        projected = mat_vec(break_projector, state)
        norm = math.sqrt(max(0.0, inner(projected, projected).real))
        check(norm > 1.0e-8, "random state has nonzero all-BREAK component")
        state_break = [value / norm for value in projected]
        evolved_break = mat_vec(propagator, state_break)
        close(expectation(evolved_break, pair), expectation(state_break, pair),
              "all-BREAK pair statistic is unchanged for arbitrary state")

    # Writer/route operations are the precise negative controls: they do not
    # commute with K projectors and are absent from the declared free/source/read schedule.
    check(max_abs(comm(xk1, p1)) > 0.5, "K-changing writer would break conservation")
    route_swap = swap(factors, 0, 2)
    check(max_abs(comm(route_swap, p1)) > 0.5,
          "K/link route swap would break conservation if left on")

    # Complete-read conservation does not forbid Bayesian updating after
    # conditioning on a correlated read.  This negative control guards the
    # theorem's unconditioned dynamical meaning.
    state = [0j] * dimension
    def basis_index(bits: tuple[int, ...]) -> int:
        output = 0
        for bit in bits:
            output = (output << 1) | bit
        return output
    state[basis_index((0, 0, 0, 0, 0))] = 1 / math.sqrt(2)  # Pi_00, M=+1
    state[basis_index((1, 1, 0, 1, 0))] = 1 / math.sqrt(2)  # Pi_11, M=-1
    pair_plus = scale(Fraction(1, 2), add(identity_full, pair))
    pair_minus = scale(Fraction(1, 2), add(identity_full, scale(-1, pair)))
    prior = expectation(state, break_projector).real
    plus_probability = expectation(state, pair_plus).real
    joint = expectation(state, mul(break_projector, pair_plus)).real
    close(prior, 0.5, "mixed-route prior")
    close(plus_probability, 0.5, "correlated read probability")
    close(joint / plus_probability, 1.0,
          "postselection can update a route posterior without back-reaction")
    nonselective = (expectation(state, mul(pair_plus, mul(break_projector, pair_plus))).real
                    + expectation(state, mul(pair_minus, mul(break_projector, pair_minus))).real)
    close(nonselective, prior, "complete nonselective read preserves route marginal")


def poisson_tail(distance: int, argument: float) -> float:
    check(distance >= 1 and argument >= 0, "Poisson-tail domain")
    if argument == 0:
        return 0.0
    logarithm = distance * math.log(argument) - math.lgamma(distance + 1)
    term = 0.0 if logarithm < -745 else math.exp(logarithm)
    total = term
    index = distance
    while term > 0:
        index += 1
        term *= argument / index
        total += term
        if term <= 1.0e-16 * max(1.0, total):
            break
        check(index < 100000, "Poisson-tail recurrence terminates")
    return total


def collar_delta(radius: int, ratio: float, time: float) -> float:
    if ratio == 0 or time == 0:
        return 0.0
    coefficient = 3 * radius * radius + 3 * radius + 1
    return min(1.0, coefficient * poisson_tail(2 * radius + 1,
                                                48 * ratio * abs(time)))


def verify_one_tail_collar_transfer() -> None:
    check(48 * Fraction(2, 1) == 96, "R=2 substitution is exact")
    check(48 * Fraction(5, 2) == 120, "R=5/2 substitution is exact")
    for radius in range(0, 80):
        for ratio in (0.0, 2.0, 2.5, 17.0):
            close(collar_delta(radius, ratio, 0.0), 0,
                  "zero-time collar error is exact")
            close(collar_delta(radius, 0.0, 4.1), 0,
                  "R=0 collar error is exact")

    # Factorial decay defeats the quadratic shell for each fixed finite input.
    for ratio, time in ((2.0, 0.01), (2.5, 0.01),
                        (2.0, 0.1), (2.5, 0.1), (7.0, 0.05)):
        values = [collar_delta(radius, ratio, time) for radius in range(1, 181)]
        check(values[-1] < 1.0e-30,
              f"collar tail converges for finite R={ratio},s={time}")
        check(all(values[index + 1] <= values[index] + 1.0e-15
                  for index in range(100, len(values) - 1)),
              f"collar tail is eventually decreasing R={ratio},s={time}")

    # Exact exhaustive binary reverse-triangle replay.  The shared BREAK
    # probability is the algebraic content of p^(Omega,0)=p^(L,0).
    grid = [Fraction(index, 24) for index in range(25)]
    for keep_full in grid:
        for keep_collar in grid:
            error = abs(keep_full - keep_collar)
            for break_shared in grid:
                distance_full = abs(keep_full - break_shared)
                distance_collar = abs(keep_collar - break_shared)
                check(abs(distance_full - distance_collar) <= error,
                      "one-tail reverse triangle is exact")
                check(max(Fraction(0), distance_collar - error) <= distance_full,
                      "one-tail lower interval is exact")
                check(distance_full <= min(Fraction(1), distance_collar + error),
                      "one-tail upper interval is exact")
                if distance_collar > error:
                    check(distance_full > 0,
                          "strict collar margin certifies exterior positivity")

    # Hostile negative control: if BREAK also has a tail, paying only one tail
    # is false, while the generic two-tail bound is sharp.
    keep_collar = Fraction(1)
    break_collar = Fraction(0)
    keep_full = Fraction(4, 5)
    break_full = Fraction(1, 5)
    delta = Fraction(1, 5)
    distance_collar = abs(keep_collar - break_collar)
    distance_full = abs(keep_full - break_full)
    check(distance_full < distance_collar - delta,
          "one-tail bound fails when BREAK is not exact")
    check(distance_full == distance_collar - 2 * delta,
          "generic two-tail lower bound is attained")


def verify_forward_positive_control() -> None:
    identity = [[1, 0], [0, 1]]
    x = [[0, 1], [1, 0]]
    z = [[1, 0], [0, -1]]
    pair = tensor(z, z)
    blank = [1 + 0j, 0j, 0j, 0j]
    for step in range(65):
        time = math.pi * step / 128
        rotation = add(scale(math.cos(time), identity),
                       scale(1j * math.sin(time), x))
        keep = mat_vec(tensor(rotation, rotation), blank)
        mean_keep = expectation(keep, pair).real
        mean_break = expectation(blank, pair).real
        expected_mean = math.cos(2 * time) ** 2
        distance = abs(mean_keep - mean_break) / 2
        close(mean_keep, expected_mean, "positive-control KEEP mean")
        close(mean_break, 1, "positive-control BREAK mean")
        close(distance, 0.5 * math.sin(2 * time) ** 2,
              "positive-control binary total variation")
    close(0.5 * math.sin(2 * (math.pi / 4)) ** 2, 0.5,
          "positive control reaches one half")


def verify_future_writer_boundary() -> None:
    identity = [[1, 0], [0, 1]]
    x = [[0, 1], [1, 0]]
    p = [[0, 0], [0, 1]]
    plus_y = [1 / math.sqrt(2), 1j / math.sqrt(2)]
    check(max_abs(comm(x, p)) > 0.5,
          "a K-changing Hamiltonian writer evades the current commutant")

    def writer_probability(parameter: float) -> float:
        writer = add(scale(math.cos(parameter / 2), identity),
                     scale(-1j * math.sin(parameter / 2), x))
        return expectation(mat_vec(writer, plus_y), p).real

    step = 1.0e-6
    derivative = (writer_probability(step) - writer_probability(-step)) / (2 * step)
    check(abs(derivative) > 0.49,
          "an explicit source-dependent K writer can give nonzero first response")
    identity_derivative = (expectation(plus_y, p).real
                           - expectation(plus_y, p).real) / (2 * step)
    close(identity_derivative, 0,
          "identity-spectator fresh carrier has zero response")


def verify_scope_text() -> None:
    theorem = (AUTHOR / "THEOREM.md").read_text(encoding="utf-8")
    result = (AUTHOR / "RESULT.md").read_text(encoding="utf-8")
    required = [
        "finite route word",
        "explicit premise",
        "formation, writer, route, and terminal-query pulses are off",
        "p^{\\Omega,0}=p^{L,0}",
        "If both comparator arms have nontrivial transverse route words",
        "including states entangled across support and link",
        "entire `K` density matrix is fixed",
        "quasi-local limit for which the inherited",
        "future writer/formation channel",
        "undefined",
        "merely appended as an identity",
        "not a no-go",
        "does not prove autonomous whole-exterior formation or routing",
    ]
    for token in required:
        check(token in theorem, f"frozen theorem scope token: {token}")
    for token in ("one-way physical-support-gated response",
                  "does not yet supply reciprocal",
                  "not a metric", "gravity, or `G` theorem"):
        check(token in result, f"frozen result ceiling token: {token}")


def main() -> None:
    verify_frozen_custody()
    verify_support_word_algebra_and_schedule()
    verify_one_tail_collar_transfer()
    verify_forward_positive_control()
    verify_future_writer_boundary()
    verify_scope_text()
    print(f"PASS__INDEPENDENT_GL6BC_HOSTILE_REPLAY__{CHECKS}/{CHECKS}")
    print("COLLAR=ALL_BREAK_EXACT;ONE_KEEP_TAIL;STRICT_MARGIN_CERTIFICATE")
    print("SUPPORT=FINITE_WORD_PROJECTORS_CONSERVED_BY_FULL_DECLARED_SCHEDULE")
    print("KERNEL=BOTH_RETARDED_ORIENTATIONS_ZERO_FOR_COMPLETE_COMMUTANT_ALGEBRA")
    print("SCOPE=UNCONDITIONED_COMPLETE_READ;FINITE_CYLINDERS_UNDER_EXHAUSTION")
    print("MISSING=ONE_SOURCE_DEPENDENT_FRESH_WRITER_INSTRUMENT_AT_OPERATIONAL_TYPE_LEVEL")


if __name__ == "__main__":
    main()
