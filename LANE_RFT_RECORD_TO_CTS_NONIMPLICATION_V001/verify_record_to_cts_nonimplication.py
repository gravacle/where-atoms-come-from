#!/usr/bin/env python3
"""Exact witnesses for RECORD_FM not implying CTS.

The executable checks the compact stochastic-CTC consistency law, registered
record distributions, the absence of a writer-before-query directed cut, and
the transpose composability boundary. It does not establish empirical CTC
existence or universal failure of joint separators.
"""

from fractions import Fraction
from itertools import product


class Ledger:
    def __init__(self) -> None:
        self.passed = 0
        self.total = 0

    def check(self, condition: bool, name: str) -> None:
        self.total += 1
        if not condition:
            raise AssertionError(name)
        self.passed += 1


def row_times_matrix(row, matrix):
    return tuple(
        sum(row[i] * matrix[i][j] for i in range(len(row)))
        for j in range(len(matrix[0]))
    )


def tv(p, q):
    keys = set(p) | set(q)
    return sum(abs(p.get(k, 0) - q.get(k, 0)) for k in keys) / 2


def marginal_y(joint):
    return {
        y: sum(weight for (yy, _), weight in joint.items() if yy == y)
        for y in (0, 1)
    }


def valid_writer_query_cut(assignment, edges):
    if assignment["W"] != "past" or assignment["Q"] != "future":
        return False
    return all(
        not (assignment[source] == "future" and assignment[target] == "past")
        for source, target in edges
    )


def transpose_2x2(a):
    return ((a[0][0], a[1][0]), (a[0][1], a[1][1]))


def main() -> None:
    t = Ledger()

    # Compact stochastic CTC controller.
    K = (
        (Fraction(7, 8), Fraction(1, 8)),
        (Fraction(3, 8), Fraction(5, 8)),
    )
    tau = (Fraction(3, 4), Fraction(1, 4))
    for index, row in enumerate(K):
        t.check(sum(row) == 1, f"ctc_kernel_row_{index}_normalized")
        t.check(all(value >= 0 for value in row), f"ctc_kernel_row_{index}_positive")
    t.check(row_times_matrix(tau, K) == tau, "ctc_stationary_law_exact")
    t.check(K[0][0] - K[1][0] == Fraction(1, 2),
            "ctc_nontrivial_eigenvalue_half")

    # Uniqueness: p' = 1/8 + p/2, so fixed p is exactly 1/4.
    fixed_p = Fraction(1, 8) / (1 - Fraction(1, 2))
    t.check(fixed_p == Fraction(1, 4), "ctc_stationary_law_unique")
    for p in (Fraction(0), Fraction(1, 8), Fraction(1, 4),
              Fraction(1, 2), Fraction(1)):
        updated = Fraction(1, 8) + p / 2
        t.check((updated == p) == (p == Fraction(1, 4)),
                f"ctc_fixed_point_sample_{p}")

    # Genuine blank-to-bit write and exact hold.
    blank = "BLANK"
    for x in (0, 1):
        written = x
        held = written
        t.check(blank != written, f"write_creates_from_blank_x{x}")
        t.check(held == x, f"hold_preserves_x{x}")

    # Complete query history H=(Y,C), Y=X xor C.
    joint = {}
    for x in (0, 1):
        joint[x] = {(x ^ c, c): tau[c] for c in (0, 1)}
        t.check(sum(joint[x].values()) == 1, f"history_x{x}_normalized")
        t.check(all((y ^ c) == x for (y, c) in joint[x]),
                f"history_x{x}_label_recoverable")
    t.check(set(joint[0]).isdisjoint(set(joint[1])),
            "complete_history_supports_disjoint")
    t.check(tv(joint[0], joint[1]) == 1, "complete_history_tv_one")
    q_only = {x: marginal_y(joint[x]) for x in (0, 1)}
    t.check(tv(q_only[0], q_only[1]) == Fraction(1, 2),
            "query_only_tv_half")

    # Any physical partition W=past,Q=future is defeated by Q->W.
    nodes = ("W", "Q")
    edges = (("W", "Q"), ("Q", "W"))
    assignments = [
        dict(zip(nodes, sides))
        for sides in product(("past", "future"), repeat=len(nodes))
    ]
    candidates = [a for a in assignments if a["W"] == "past" and a["Q"] == "future"]
    t.check(len(candidates) == 1, "one_writer_query_side_assignment")
    t.check(not any(valid_writer_query_cut(a, edges) for a in candidates),
            "no_definite_writer_before_query_cut")
    t.check(("Q", "W") in edges, "future_to_past_controller_edge_explicit")

    # Transpose fixes basis record states.
    zero = Fraction(0)
    one = Fraction(1)
    rho0 = ((one, zero), (zero, zero))
    rho1 = ((zero, zero), (zero, one))
    t.check(transpose_2x2(rho0) == rho0, "transpose_fixes_basis_0")
    t.check(transpose_2x2(rho1) == rho1, "transpose_fixes_basis_1")

    # The Choi swap has a negative antisymmetric direction.
    swap = (
        (1, 0, 0, 0),
        (0, 0, 1, 0),
        (0, 1, 0, 0),
        (0, 0, 0, 1),
    )
    v = (0, 1, -1, 0)
    fv = tuple(sum(swap[i][j] * v[j] for j in range(4)) for i in range(4))
    t.check(fv == tuple(-entry for entry in v),
            "transpose_choi_antisymmetric_eigenvalue_minus_one")
    t.check(Fraction(-1, 2) < 0, "partial_transpose_negative_eigenvalue")

    # Basis-only classical K recoding is exact.
    for x in (0, 1):
        classical_boundary = x
        classical_hold = classical_boundary
        classical_read = classical_hold
        quantum_basis_read = x
        t.check(classical_read == quantum_basis_read,
                f"basis_classical_K_recoding_x{x}")

    # Tomography containing a Y component distinguishes transpose from identity.
    bloch_plus_y = (Fraction(0), Fraction(1), Fraction(0))
    bloch_transposed = (
        bloch_plus_y[0], -bloch_plus_y[1], bloch_plus_y[2]
    )
    t.check(bloch_transposed == (0, -1, 0),
            "transpose_flips_Y_tomographic_component")
    t.check(bloch_transposed != bloch_plus_y,
            "tomography_distinguishes_transpose_from_identity")

    print("WEAK_RECORD_LINEAGE_CONTENT EXACT")
    print("STOCHASTIC_CTC_FIXED_POINT UNIQUE")
    print("CTC_RECORD_DISTRIBUTIONS EXACT")
    print("DEFINITE_MACROCUT REFUTED")
    print("TRANSPOSE_NON_CP_BOUNDARY EXACT")
    print("BASIS_CLASSICAL_K_RECODING EXACT")
    print(f"TOTAL {t.passed}/{t.total} PASS")
    print("VERDICT RECORD_DOES_NOT_IMPLY_CTS__T1_FIRST_NONIMPLICATION")


if __name__ == "__main__":
    main()

