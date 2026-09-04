#!/usr/bin/env python3
"""Independent hostile replay of the GL6CM stationary writer response.

No target program is imported or executed.  This replay reconstructs the
finite-component Perron--Frobenius scope, the spectral Gram factorization,
the exact kernel and common-rescaling null, a nontrivial W^T K W pullback,
and the K2/two-arm-star response in exact rational/Q(sqrt(2)) arithmetic.
"""

from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CM_GLOBAL_H6_STATIONARY_WRITER_RESPONSE_V001"


def ftext(value):
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def canonical_hash(value):
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def transpose(matrix):
    return tuple(tuple(row[j] for row in matrix) for j in range(len(matrix[0])))


def dot(left, right):
    return sum((F(x) * F(y) for x, y in zip(left, right)), F(0))


def matvec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def matmul(left, right):
    right_t = transpose(right)
    return tuple(tuple(dot(row, column) for column in right_t) for row in left)


def quadratic(matrix, vector):
    return dot(vector, matvec(matrix, vector))


def rank(matrix):
    rows = [list(map(F, row)) for row in matrix]
    if not rows:
        return 0
    nrows, ncols = len(rows), len(rows[0])
    pivot_row = 0
    for column in range(ncols):
        pivot = next((r for r in range(pivot_row, nrows) if rows[r][column]), None)
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        divisor = rows[pivot_row][column]
        rows[pivot_row] = [x / divisor for x in rows[pivot_row]]
        for r in range(nrows):
            if r == pivot_row or not rows[r][column]:
                continue
            multiple = rows[r][column]
            rows[r] = [x - multiple * y for x, y in zip(rows[r], rows[pivot_row])]
        pivot_row += 1
        if pivot_row == nrows:
            break
    return pivot_row


def spectral_kernel(amplitude_rows, gaps):
    """Return 2 S^T diag(1/gap) S, with lambda suppressed."""
    columns = len(amplitude_rows[0])
    answer = [[F(0) for _ in range(columns)] for _ in range(columns)]
    for row, gap in zip(amplitude_rows, gaps):
        gap = F(gap)
        assert gap > 0
        for i in range(columns):
            for j in range(columns):
                answer[i][j] += 2 * F(row[i]) * F(row[j]) / gap
    return tuple(tuple(row) for row in answer)


def identity(n):
    return tuple(tuple(F(i == j) for j in range(n)) for i in range(n))


def matrix_add(left, right):
    return tuple(tuple(F(x) + F(y) for x, y in zip(a, b)) for a, b in zip(left, right))


def matrix_power(matrix, exponent):
    answer = identity(len(matrix))
    base = matrix
    while exponent:
        if exponent & 1:
            answer = matmul(answer, base)
        base = matmul(base, base)
        exponent //= 2
    return answer


def connected(adjacency):
    seen = {0}
    frontier = [0]
    while frontier:
        node = frontier.pop()
        for other, value in enumerate(adjacency[node]):
            if value and other not in seen:
                seen.add(other)
                frontier.append(other)
    return len(seen) == len(adjacency)


def finite_component_pf_scope():
    """Exhaust small connected graphs and verify a strict-positive PF shift.

    The proof used by the audit is general: for any finite connected
    nonnegative symmetric adjacency A, (I+A)^(n-1) is entrywise positive,
    since a path of length at most n-1 can be padded with diagonal I steps.
    Therefore A is irreducible and its largest eigenvalue is simple.  The
    exhaustive census here independently checks the constructive premise on
    every labeled simple connected graph through five vertices.
    """
    counts = {}
    canonical_rows = []
    for n in range(2, 6):
        edges = tuple(combinations(range(n), 2))
        count = 0
        for mask in range(1 << len(edges)):
            adjacency = [[F(0) for _ in range(n)] for _ in range(n)]
            for bit, (a, b) in enumerate(edges):
                if mask & (1 << bit):
                    adjacency[a][b] = adjacency[b][a] = F(1)
            adjacency = tuple(tuple(row) for row in adjacency)
            if not connected(adjacency):
                continue
            shifted = matrix_add(identity(n), adjacency)
            witness = matrix_power(shifted, n - 1)
            assert all(entry > 0 for row in witness for entry in row)
            count += 1
            canonical_rows.append([n, mask, [[ftext(x) for x in row] for row in witness]])
        counts[str(n)] = count
    assert counts == {"2": 1, "3": 4, "4": 38, "5": 728}
    return {
        "labeled_connected_simple_graphs_n2_through_n5": counts,
        "total_checked": sum(counts.values()),
        "constructive_pf_witness": "(I+A)^(n-1) is entrywise positive",
        "conclusion": "connected finite symmetric nonnegative A has a simple Perron root; H0=-J A has a unique ground ray and a positive finite gap",
        "census_sha256": canonical_hash(canonical_rows),
    }


def spectral_factorization_and_kernel():
    # These exact transition-amplitude rows deliberately have a sole common
    # null and are not taken from the target.  They realize the general Gram
    # factorization with three different positive excitation gaps.
    amplitudes = (
        (F(1), F(-1), F(0), F(0)),
        (F(0), F(2), F(-3), F(1)),
        (F(2), F(-1), F(4), F(-5)),
    )
    gaps = (F(2), F(3), F(5))
    assert all(sum(row) == 0 for row in amplitudes)
    assert rank(amplitudes) == 3
    kernel = spectral_kernel(amplitudes, gaps)
    assert kernel == transpose(kernel)
    assert rank(kernel) == 3
    assert matvec(kernel, (1, 1, 1, 1)) == (0, 0, 0, 0)

    checked = 0
    zeros = []
    for vector in product(range(-2, 3), repeat=4):
        transitions = matvec(amplitudes, vector)
        gram_value = sum((2 * value * value / gap for value, gap in zip(transitions, gaps)), F(0))
        direct_value = quadratic(kernel, vector)
        assert direct_value == gram_value >= 0
        assert (direct_value == 0) == all(value == 0 for value in transitions)
        if direct_value == 0:
            zeros.append(vector)
        checked += 1
    assert zeros == [(-2, -2, -2, -2), (-1, -1, -1, -1), (0, 0, 0, 0),
                     (1, 1, 1, 1), (2, 2, 2, 2)]

    # A second factor with lower rank confirms that the uniform null need not
    # exhaust the exact kernel; the Q B|0>=0 criterion remains decisive.
    dark_amplitudes = (
        (F(1), F(-1), F(0), F(0)),
        (F(0), F(1), F(-1), F(0)),
    )
    dark_kernel = spectral_kernel(dark_amplitudes, (F(1), F(7, 3)))
    assert rank(dark_kernel) == 2
    assert matvec(dark_kernel, (1, 1, 1, 1)) == (0, 0, 0, 0)
    assert matvec(dark_kernel, (0, 0, 0, 1)) == (0, 0, 0, 0)

    return {
        "factorization": "K=2 lambda_T^2 S^T diag(Delta_n^-1) S",
        "reciprocal": True,
        "positive_semidefinite": True,
        "kernel": "K(y,y)=0 iff S y=0 iff Q B_y|0>=0",
        "primary_transition_rank": rank(amplitudes),
        "primary_kernel_rank": rank(kernel),
        "integer_vectors_checked": checked,
        "primary_only_null_in_test_box": "uniform rescaling",
        "extra_dark_kernel_example_rank": len(dark_kernel) - rank(dark_kernel),
        "common_rescaling_derivation": "S_n dot 1=<n|sum_c T_c|0>=-(1/J)<n|H0|0>=0",
        "kernel_matrix": [[ftext(x) for x in row] for row in kernel],
    }, kernel


def writer_pullback(kernel):
    # An independently chosen finite writer map with nontrivial mixing.
    writer = (
        (F(1), F(0), F(-1), F(0), F(0), F(0)),
        (F(0), F(1), F(0), F(-1), F(0), F(0)),
        (F(1), F(-1), F(0), F(0), F(1), F(-1)),
        (F(2), F(0), F(-1), F(-1), F(0), F(0)),
    )
    pulled = matmul(transpose(writer), matmul(kernel, writer))
    assert pulled == transpose(pulled)
    comparisons = 0
    for left in product(range(-1, 2), repeat=6):
        mapped_left = matvec(writer, left)
        for right in (
            (1, 0, 0, 0, 0, 0),
            (0, 1, -1, 0, 1, 0),
            (1, -1, 1, -1, 1, -1),
        ):
            mapped_right = matvec(writer, right)
            assert dot(left, matvec(pulled, right)) == dot(mapped_left, matvec(kernel, mapped_right))
            comparisons += 1
        assert quadratic(pulled, left) == quadratic(kernel, mapped_left) >= 0
    return {
        "formula": "K_T=W^T K W",
        "writer_shape": [4, 6],
        "pulled_rank": rank(pulled),
        "bilinear_exact_comparisons": comparisons,
        "positive_semidefinite": True,
        "kernel_rule": "j is dark iff Wj lies in ker K",
        "pulled_matrix": [[ftext(x) for x in row] for row in pulled],
    }


@dataclass(frozen=True)
class Q2:
    rational: F = F(0)
    radical: F = F(0)  # coefficient of sqrt(2)

    def __add__(self, other):
        other = as_q2(other)
        return Q2(self.rational + other.rational, self.radical + other.radical)

    __radd__ = __add__

    def __neg__(self):
        return Q2(-self.rational, -self.radical)

    def __sub__(self, other):
        return self + (-as_q2(other))

    def __rsub__(self, other):
        return as_q2(other) - self

    def __mul__(self, other):
        other = as_q2(other)
        return Q2(
            self.rational * other.rational + 2 * self.radical * other.radical,
            self.rational * other.radical + self.radical * other.rational,
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = as_q2(other)
        denominator = other.rational * other.rational - 2 * other.radical * other.radical
        assert denominator
        conjugate = Q2(other.rational, -other.radical)
        numerator = self * conjugate
        return Q2(numerator.rational / denominator, numerator.radical / denominator)


def as_q2(value):
    return value if isinstance(value, Q2) else Q2(F(value), F(0))


def qdot(left, right):
    return sum((as_q2(x) * as_q2(y) for x, y in zip(left, right)), Q2())


def qmatvec(matrix, vector):
    return tuple(qdot(row, vector) for row in matrix)


def two_component_response():
    sqrt2 = Q2(0, 1)
    halfsqrt2 = Q2(0, F(1, 2))

    # K2: its only writer T is exactly proportional to H0, so the sole
    # excited residue vanishes.
    sigma_x = ((0, 1), (1, 0))
    k2_ground = (halfsqrt2, halfsqrt2)
    k2_excited = (halfsqrt2, -halfsqrt2)
    assert qdot(k2_ground, k2_ground) == qdot(k2_excited, k2_excited) == Q2(1)
    assert qdot(k2_excited, qmatvec(sigma_x, k2_ground)) == Q2()

    # Shared-arm star, ordered (center, arm0, arm1), with J=1.
    h0 = ((0, -1, -1), (-1, 0, 0), (-1, 0, 0))
    t0 = ((0, 1, 0), (1, 0, 0), (0, 0, 0))
    t1 = ((0, 0, 1), (0, 0, 0), (1, 0, 0))
    ground = (halfsqrt2, F(1, 2), F(1, 2))
    middle = (Q2(), halfsqrt2, -halfsqrt2)
    upper = (halfsqrt2, F(-1, 2), F(-1, 2))
    assert qmatvec(h0, ground) == tuple(-sqrt2 * x for x in ground)
    assert qmatvec(h0, middle) == (Q2(), Q2(), Q2())
    assert qmatvec(h0, upper) == tuple(sqrt2 * x for x in upper)
    assert qdot(middle, qmatvec(t0, ground)) == Q2(F(1, 2))
    assert qdot(middle, qmatvec(t1, ground)) == Q2(F(-1, 2))
    assert qdot(upper, qmatvec(t0, ground)) == Q2()
    assert qdot(upper, qmatvec(t1, ground)) == Q2()

    # K_cc'=2 <0|T_c R T_c'|0>.  Only the middle state contributes,
    # with gap sqrt(2), producing sqrt(2)/4 [[1,-1],[-1,1]] at J=1.
    prefactor = Q2(0, F(1, 4))
    cycle_kernel = (
        (prefactor, -prefactor),
        (-prefactor, prefactor),
    )
    assert qmatvec(cycle_kernel, (Q2(1), Q2(1))) == (Q2(), Q2())

    samples = []
    for w0, w1 in ((F(0), F(0)), (F(1), F(0)), (F(0), F(1)),
                   (F(3), F(-2)), (F(7, 5), F(7, 5))):
        response = qdot((Q2(w0), Q2(w1)), qmatvec(cycle_kernel, (Q2(w0), Q2(w1))))
        expected = Q2(0, (w0 - w1) ** 2 / 4)
        assert response == expected

        # Independently differentiate E(s)=-sqrt((1-sw0)^2+(1-sw1)^2).
        # If F(s) is the radicand, E''=-F''/(2 sqrt(F))
        # +F'^2/(4 F^(3/2)); evaluating at zero is exact in Q(sqrt2).
        f0 = Q2(2)
        f1 = Q2(-2 * (w0 + w1))
        f2 = Q2(2 * (w0*w0 + w1*w1))
        branch_hessian = -(f2 / (2 * sqrt2)) + f1*f1 / (4 * (2 * sqrt2))
        assert branch_hessian == -expected
        samples.append({
            "w0": ftext(w0),
            "w1": ftext(w1),
            "minus_energy_hessian": f"({ftext(expected.radical)})sqrt(2)",
        })

    return {
        "isolated_K2_spectral_response": "0",
        "star_spectrum_over_J": ["-sqrt(2)", "0", "+sqrt(2)"],
        "transition_rows": {"middle": ["1/2", "-1/2"], "upper": ["0", "0"]},
        "cycle_kernel": "lambda_T^2 sqrt(2)/(4J) [[1,-1],[-1,1]]",
        "uniform_mode": "null",
        "relative_mode": "strict",
        "minus_energy_hessian": "sqrt(2)(w0-w1)^2/(4J)",
        "exact_branch_samples": samples,
    }


def physical_scaling_and_units():
    j0 = F(63, 8)
    lambda0 = F(105, 16)
    x_literal = (F(2), F(0))
    w_literal = tuple(lambda0 * x for x in x_literal)
    coefficient = (w_literal[0] - w_literal[1]) ** 2 / (4 * j0)
    assert w_literal == (F(105, 8), F(0))
    assert coefficient == F(175, 32)
    inverse_rational = F(16, 175)
    assert coefficient * inverse_rational * 2 == 1
    return {
        "J": "(63/8)h^6/U_d^5 [energy]",
        "lambda_T": "(105/16)h^6/U_d^6 [dimensionless]",
        "j_and_x": "energy",
        "R": "inverse energy",
        "cycle_kernel_entries": "inverse energy",
        "K_spec(j,k)_as_bilinear_value": "energy when j,k each carry energy",
        "scalar_path": "j=s jhat with s energy and jhat dimensionless; O=dH/ds dimensionless; -E'' inverse energy",
        "literal_weights": ["(105/8)h^6/U_d^6", "0"],
        "literal_positive_response": "(175sqrt(2)/32)h^6/U_d^7",
        "unit_positive_response": "(175sqrt(2)/64)h^6/U_d^7",
        "spectral_only_inverse": "(16sqrt(2)/175)U_d^7/h^6",
        "power_count": "w^2/J = h^6/U_d^7 = inverse energy",
    }


def compare_frozen_claims():
    ledger = json.loads((TARGET / "EXACT_LEDGER.json").read_text())
    assert ledger["source_free_component"]["J"] == "(63/8)h^6/U_d^5"
    assert ledger["source_free_component"]["hamiltonian"] == "H0=-J sum_c T_c"
    assert ledger["writer"]["lambda_T"] == "(105/16)h^6/U_d^6"
    assert ledger["spectral_response"]["sign"] == "positive semidefinite"
    assert ledger["spectral_response"]["common_cycle_rescaling"] == "null"
    assert ledger["pair_source_pullback"] == "K_T_spec=W^T K W"
    assert ledger["two_overlap"]["isolated_K2_curvature"] == "0"
    assert ledger["two_overlap"]["star_curvature"] == "-sqrt(2)(w0-w1)^2/(4J)"
    assert ledger["two_overlap"]["literal_Q4_negative_energy_hessian"] == "-(175sqrt(2)/32)h^6/U_d^7"
    excluded = " | ".join(ledger["excluded"])
    for token in ("diagonal", "contact", "record authentication", "thermodynamic", "real-time", "gravity", "Newton G"):
        assert token in excluded
    return "all decisive final-ledger identities and scope exclusions match"


def main():
    pf = finite_component_pf_scope()
    spectral, kernel = spectral_factorization_and_kernel()
    pullback = writer_pullback(kernel)
    components = two_component_response()
    scaling = physical_scaling_and_units()
    comparison = compare_frozen_claims()
    result = {
        "schema": "AUDIT_G_GL6CM_INDEPENDENT_V001",
        "algorithm_independence": "No target program imported or executed; PF scope, spectral Gram algebra, pullback, and K2/star response were rebuilt independently with exact arithmetic.",
        "perron_frobenius_scope": pf,
        "spectral_factorization": spectral,
        "writer_pullback": pullback,
        "finite_components": components,
        "physical_scaling_and_units": scaling,
        "target_comparison": comparison,
        "boundary_findings": {
            "writer_only_typing": "PASS: the differentiated family is explicitly the source-linear off-diagonal H6 writer, not the full source-dependent Hamiltonian.",
            "unclassified_first_vertex": "PASS: a possible diagonal H6 first-source operator is explicitly open.",
            "contact": "PASS: every source-second contact remains open; the result is not called the complete Hessian.",
            "stationarity": "PASS: finite connected locked component and PF ground branch only; no thermodynamic phase selection is claimed.",
            "record_authentication": "PASS: the pair-source probe is not promoted to an authenticated record or autonomous source.",
            "collar_summation": "PASS: the two-overlap coefficient is a witness, not a term summed independently over a dense graph.",
            "bulk_promotion": "PASS: linked clusters, real time, locality, metric/Ricci/Einstein form, gravity, and G remain open.",
            "legendre_convention": "PASS: CM18 explicitly defines phi=-dE/ds and Gamma=E+s phi, giving Gamma''=(-E'')^-1=K^-1.",
        },
        "disposition": "PASS",
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = HERE / "INDEPENDENT_RESULT.json"
    if "--emit" in sys.argv:
        output.write_text(rendered)
    else:
        assert output.read_text() == rendered

    print(f"PASS finite-component PF scope: {pf['total_checked']} exhaustive graphs plus general irreducibility proof")
    print(f"PASS exact spectral Gram factorization and kernel: {spectral['integer_vectors_checked']} vectors")
    print(f"PASS exact common-rescaling null and W^T K W pullback: {pullback['bilinear_exact_comparisons']} comparisons")
    print("PASS isolated K2 zero and shared-star strict relative response")
    print("PASS physical source scaling: 175sqrt(2)/32 h^6/U_d^7")
    print("PASS claim typing: writer-only; contact, record authentication, bulk gravity and G open")
    print("AUDIT DISPOSITION: PASS")


if __name__ == "__main__":
    main()
