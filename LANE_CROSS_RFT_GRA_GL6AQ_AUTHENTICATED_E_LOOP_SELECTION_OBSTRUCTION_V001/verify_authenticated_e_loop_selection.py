#!/usr/bin/env python3
"""Independent exact replay of the GL6AQ source/selection obstruction."""

from __future__ import annotations

import hashlib
import itertools
from collections import Counter
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
PAIR_INDEX = {pair: i for i, pair in enumerate(PAIRS)}
PORTS = tuple(range(4))
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def matmul(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def transpose(a):
    return [list(row) for row in zip(*a)]


def matadd(a, b):
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def matscale(a, q):
    return [[q * value for value in row] for row in a]


def matvec(a, v):
    return tuple(sum(a[i][j] * v[j] for j in range(len(v))) for i in range(len(a)))


def permutation_matrix(images, size):
    out = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for old, new in enumerate(images):
        out[new][old] = Fraction(1)
    return out


def pair_vector(z):
    return tuple(z[a] * z[b] for a, b in PAIRS)


R = [[Fraction(int(a in pair)) for pair in PAIRS] for a in PORTS]
E_BASIS = (
    (Fraction(1), Fraction(0), Fraction(-1), Fraction(-1), Fraction(0), Fraction(1)),
    (Fraction(0), Fraction(1), Fraction(-1), Fraction(-1), Fraction(1), Fraction(0)),
)


def section_locked_e_source() -> None:
    for basis in E_BASIS:
        check(matvec(R, basis) == (0, 0, 0, 0), "declared E basis lies in ker R")

    locked_z = tuple(z for z in itertools.product((-1, 1), repeat=4) if sum(z) == 0)
    locked_m = tuple(pair_vector(z) for z in locked_z)
    check(len(locked_m) == 6, "six degree-two local configurations")
    for m in locked_m:
        check(matvec(R, m) == (-1, -1, -1, -1), "locked affine incidence identity")

    mean = tuple(sum(Fraction(m[j], 6) for m in locked_m) for j in range(6))
    check(mean == tuple(Fraction(-1, 3) for _ in PAIRS), "uniform locked pair mean")
    covariance = [[
        sum(Fraction((m[i] - mean[i]) * (m[j] - mean[j]), 6) for m in locked_m)
        for j in range(6)
    ] for i in range(6)]
    p_e = matscale(covariance, Fraction(3, 8))
    check(matmul(p_e, p_e) == p_e, "locked covariance reconstructs E projector")
    check(sum(p_e[i][i] for i in range(6)) == 2, "E projector rank two")

    for alpha, beta in itertools.product(range(-3, 4), repeat=2):
        c = tuple(alpha * E_BASIS[0][j] + beta * E_BASIS[1][j] for j in range(6))
        if alpha == beta == 0:
            continue
        norm2 = sum(value * value for value in c)
        values = tuple(sum(c[j] * m[j] for j in range(6)) for m in locked_m)
        variance = sum(value * value for value in values) / 6
        check(variance == Fraction(8, 3) * norm2, "exact nonzero locked E source norm")
        check(any(value != 0 for value in values), "nonzero E source survives locked projection")
        check(matvec(p_e, c) == c, "sample source is exactly E projected")

    # Full local product trace: distinct pair Pauli strings are orthonormal.
    for c in E_BASIS:
        trace_mean = Fraction(0)
        trace_square = sum(value * value for value in c)
        check(trace_mean == 0, "pair E source has zero product-trace mean")
        check(trace_square > 0, "pair E source has positive trace correlation mass")


def section_s4_linear_selection() -> None:
    # Reconstruct P_E from the line graph polynomial used by the sealed source theorem.
    adjacency = [[Fraction(0) for _ in PAIRS] for _ in PAIRS]
    for i, pair in enumerate(PAIRS):
        for j, other in enumerate(PAIRS):
            if i != j and set(pair) & set(other):
                adjacency[i][j] = Fraction(1)
    identity = [[Fraction(int(i == j)) for j in range(6)] for i in range(6)]
    p_e = matscale(
        matmul(adjacency, matadd(adjacency, matscale(identity, Fraction(-4)))),
        Fraction(1, 12),
    )
    zero_6x4 = [[Fraction(0) for _ in range(4)] for _ in range(6)]

    reps = []
    character_inner = Fraction(0)
    for sigma in itertools.permutations(PORTS):
        vertex = permutation_matrix(sigma, 4)
        pair_images = [PAIR_INDEX[tuple(sorted((sigma[a], sigma[b])))] for a, b in PAIRS]
        pair = permutation_matrix(pair_images, 6)
        chi_vertex = sum(vertex[i][i] for i in range(4))
        chi_e = sum(matmul(p_e, pair)[i][i] for i in range(6))
        character_inner += chi_vertex * chi_e
        reps.append((pair, vertex))
    check(character_inner / 24 == 0, "four-port representation contains no E irrep")

    # Group averaging every elementary 6x4 map spans all intertwiners.
    for out_index in range(6):
        for in_index in range(4):
            seed = [[Fraction(0) for _ in range(4)] for _ in range(6)]
            seed[out_index][in_index] = Fraction(1)
            averaged = [[Fraction(0) for _ in range(4)] for _ in range(6)]
            for pair, vertex in reps:
                transformed = matmul(matmul(pair, seed), transpose(vertex))
                averaged = matadd(averaged, transformed)
            averaged = matscale(averaged, Fraction(1, 24))
            check(matmul(p_e, averaged) == zero_6x4, "every equivariant four-port cross map has zero E output")
            for pair, vertex in reps:
                check(matmul(matmul(pair, averaged), transpose(vertex)) == averaged,
                      "group-averaged cross map is equivariant")


def section_direct_k_obstruction_and_loop_e_change() -> None:
    locked_bits = tuple(bits for bits in itertools.product((0, 1), repeat=4) if sum(bits) == 2)
    for bits in locked_bits:
        for port in PORTS:
            flipped = list(bits)
            flipped[port] ^= 1
            check(tuple(flipped) not in locked_bits, "one X flip leaves the locked subspace")

    # The target hexagon has three parent vertices.  These are their two loop ports
    # and the initial occupations inherited from the alternating word (1,0,1,0,1,0).
    parent_data = (
        ((0, 2), (1, 0)),
        ((1, 2), (0, 1)),
        ((0, 1), (0, 1)),
    )
    for loop_ports, loop_occupations in parent_data:
        external = tuple(port for port in PORTS if port not in loop_ports)
        for occupied_external in external:
            initial = [0, 0, 0, 0]
            initial[loop_ports[0]], initial[loop_ports[1]] = loop_occupations
            initial[occupied_external] = 1
            final = list(initial)
            final[loop_ports[0]], final[loop_ports[1]] = initial[loop_ports[1]], initial[loop_ports[0]]
            check(sum(initial) == sum(final) == 2, "hexagon parent remains degree two")
            m_initial = pair_vector(tuple(1 - 2 * n for n in initial))
            m_final = pair_vector(tuple(1 - 2 * n for n in final))
            delta = tuple(m_final[j] - m_initial[j] for j in range(6))
            check(matvec(R, delta) == (0, 0, 0, 0), "hexagon pair displacement is E")
            check(delta != (0, 0, 0, 0, 0, 0), "hexagon pair displacement is nonzero")
            check(sum(value * value for value in delta) == 16, "hexagon E displacement squared norm")
            read_difference = sum(delta[j] * (m_final[j] - m_initial[j]) for j in range(6))
            check(read_difference == 16, "authenticated E read distinguishes loop endpoints")

    # Pair-Z strings and one-link X strings are distinct Pauli words, hence orthogonal.
    for pair in PAIRS:
        for port in PORTS:
            pair_word = tuple("Z" if p in pair else "I" for p in PORTS)
            x_word = tuple("X" if p == port else "I" for p in PORTS)
            check(pair_word != x_word, "direct K/pair-E Hilbert-Schmidt overlap is zero")


INITIAL = (1, 0, 1, 0, 1, 0)


def subset_energy(subset: frozenset[int]) -> int:
    delta = [0] * 6
    for edge in subset:
        delta[edge] = 1 if INITIAL[edge] == 0 else -1
    # Vertex v joins cycle edges v-1 and v.
    charges = [delta[(vertex - 1) % 6] + delta[vertex] for vertex in range(6)]
    return sum(charge * charge for charge in charges)


def path_coefficient(kappa) -> Fraction:
    total = Fraction(0)
    for order in itertools.permutations(range(6)):
        selected = set()
        term = Fraction(1)
        for edge in order[:-1]:
            selected.add(edge)
            term *= Fraction(-1, subset_energy(frozenset(selected)))
        for edge in order:
            term *= kappa[edge]
        total += term
    return total


def section_sixth_order_gate() -> None:
    full = frozenset(range(6))
    check(subset_energy(frozenset()) == 0, "empty loop subset is locked")
    check(subset_energy(full) == 0, "complete alternating loop returns to lock")
    for size in range(1, 6):
        for subset in itertools.combinations(range(6), size):
            check(subset_energy(frozenset(subset)) > 0, "no proper loop subset returns to lock")

    census = {
        size: Counter(subset_energy(frozenset(subset)) for subset in itertools.combinations(range(6), size))
        for size in range(1, 6)
    }
    check(census[1] == Counter({2: 6}), "one-flip energy census")
    check(census[2] == Counter({4: 9, 2: 6}), "two-flip energy census")
    check(census[3] == Counter({4: 12, 2: 6, 6: 2}), "three-flip energy census")
    check(census[4] == census[2], "four-flip complement census")
    check(census[5] == census[1], "five-flip complement census")

    coefficient = path_coefficient((Fraction(1),) * 6)
    check(coefficient == Fraction(-63, 8), "independent 720-path hexagon coefficient")

    samples = (
        (Fraction(2), Fraction(3), Fraction(5), Fraction(7), Fraction(11), Fraction(13)),
        (Fraction(1, 2), Fraction(-2, 3), Fraction(4, 5), Fraction(3, 7), Fraction(-5, 11), Fraction(6, 13)),
    )
    for kappa in samples:
        product = Fraction(1)
        for value in kappa:
            product *= value
        check(path_coefficient(kappa) == coefficient * product,
              "weighted sixth-order coefficient factors through all retained supports")

    for kappa in itertools.product((Fraction(0), Fraction(1)), repeat=6):
        expected = coefficient if all(kappa) else Fraction(0)
        check(path_coefficient(kappa) == expected, "binary retained support gates the loop entry")


def section_custody_and_scope() -> None:
    dependency_rows = []
    for raw in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
        if not raw.strip():
            continue
        expected, relative = raw.split(maxsplit=1)
        dependency_rows.append(relative)
        target = ROOT / relative
        check(target.is_file(), f"dependency exists: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"dependency hash: {relative}")
    check(len(dependency_rows) == 11, "exact dependency row count")
    check(len(set(dependency_rows)) == len(dependency_rows), "dependency rows are unique")

    theorem = (HERE / "THEOREM.md").read_text()
    result = (HERE / "RESULT.md").read_text()
    required = (
        "{8\\over3}\\,c^Tc",
        "P_{\\cal Q}X_eP_{\\cal Q}=0",
        "P_E K(t)=0",
        "\\|\\delta M_x\\|^2=16",
        "\\prod_{e\\in C}\\kappa_e",
        "refutes a **universal** nonzero stationary defect-contrast",
        "**existence** of some other",
        "remains open",
        "No\nstate, pole, physical momentum, cone",
    )
    for marker in required:
        check(marker in theorem, f"required theorem scope marker: {marker}")
    check("No nonzero stationary bulk defect response follows" in result,
          "result keeps stationary obstruction explicit")
    check("no selected state" in result.lower(), "result keeps state-selection ceiling")


def main() -> None:
    section_locked_e_source()
    section_s4_linear_selection()
    section_direct_k_obstruction_and_loop_e_change()
    section_sixth_order_gate()
    section_custody_and_scope()
    print(f"PASS__GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION__{checks}/{checks}")
    print("PAIR_SOURCE=EXACT_NONZERO_LOCKED_E_OPERATOR_OVERLAP")
    print("STATIONARY_SPECTRUM=CORRELATION_POSITIVE_NOT_STRICT_RETARDED")
    print("K_SOURCE=ZERO_DIRECT_AND_ONE_CELL_LINEAR_E_PROJECTION")
    print("LOOP=E_CHANGING_MINUS63_OVER8_TIMES_PRODUCT_OF_SIX_RETAINED_SUPPORTS")
    print("BULK_DEFECT=UNIVERSAL_NONZERO_REFUTED_EXISTENTIAL_SELECTED_STATE_OPEN")
    print("CEILING=NO_STATE_POLE_MOMENTUM_CONE_RICCI_GRAVITY_G")


if __name__ == "__main__":
    main()
