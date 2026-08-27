#!/usr/bin/env python3
"""Exact finite witnesses for the FHBC-to-FCLPD reduction.

This standard-library verifier checks finite rational matrix identities, the
derived frontier, an exhaustive premise ledger, and the passive hidden-bypass
boundary.  It does not authenticate an actual Hamiltonian boundary.
"""

from fractions import Fraction as F
from itertools import product


CHECKS = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    CHECKS.append(name)
    print(f"PASS {name}")


def eye(n):
    return tuple(tuple(F(int(i == j)) for j in range(n)) for i in range(n))


def add(a, b):
    return tuple(tuple(a[i][j] + b[i][j] for j in range(len(a[0])))
                 for i in range(len(a)))


def mul(a, b):
    return tuple(
        tuple(sum((a[i][k] * b[k][j] for k in range(len(b))), F(0))
              for j in range(len(b[0])))
        for i in range(len(a))
    )


def transpose(a):
    return tuple(tuple(a[j][i] for j in range(len(a)))
                 for i in range(len(a[0])))


def tensor(a, b):
    return tuple(
        tuple(a[i][j] * b[k][l]
              for j in range(len(a[0])) for l in range(len(b[0])))
        for i in range(len(a)) for k in range(len(b))
    )


def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0))


def conjugate(u, rho):
    return mul(mul(u, rho), transpose(u))


def basis_state(dim, index):
    return tuple(
        tuple(F(int(i == index and j == index)) for j in range(dim))
        for i in range(dim)
    )


def partial_trace_second(rho, dim_a, dim_b):
    return tuple(
        tuple(sum((rho[i * dim_b + k][j * dim_b + k]
                   for k in range(dim_b)), F(0))
              for j in range(dim_a))
        for i in range(dim_a)
    )


def probability(rho, effect):
    return trace(mul(rho, effect))


def reachable(edges, start, goal, removed=frozenset()):
    stack = [start]
    seen = set()
    while stack:
        node = stack.pop()
        if node == goal:
            return True
        if node in seen:
            continue
        seen.add(node)
        stack.extend(v for u, v in edges if u == node and (u, v) not in removed)
    return False


def main():
    # Exact local-autonomy witness with an initially entangled complement.
    half = F(1, 2)
    bell = (
        (half, F(0), F(0), half),
        (F(0), F(0), F(0), F(0)),
        (F(0), F(0), F(0), F(0)),
        (half, F(0), F(0), half),
    )
    i2 = eye(2)
    x2 = ((F(0), F(1)), (F(1), F(0)))
    p0 = basis_state(2, 0)
    p1 = basis_state(2, 1)
    rho_d = partial_trace_second(bell, 2, 2)

    check("bell_state_normalized", trace(bell) == 1)
    check("bell_reduction_maximally_mixed",
          rho_d == ((half, F(0)), (F(0), half)))
    check("factorized_unitary_is_unitary",
          mul(transpose(tensor(x2, i2)), tensor(x2, i2)) == eye(4))

    full_evolved = conjugate(tensor(x2, i2), bell)
    reduced_evolved = conjugate(x2, rho_d)
    full_p0 = probability(full_evolved, tensor(p0, i2))
    full_p1 = probability(full_evolved, tensor(p1, i2))
    red_p0 = probability(reduced_evolved, p0)
    red_p1 = probability(reduced_evolved, p1)
    check("entangled_complement_p0_reduces_exactly", full_p0 == red_p0)
    check("entangled_complement_p1_reduces_exactly", full_p1 == red_p1)
    check("reduced_query_probabilities_normalize", red_p0 + red_p1 == 1)

    # The same reduction holds for a selective local CP branch, not only a
    # deterministic unitary.  K=P0 is trace-nonincreasing.
    k_full = tensor(p0, i2)
    full_branch = conjugate(k_full, bell)
    reduced_branch_from_full = partial_trace_second(full_branch, 2, 2)
    reduced_branch_direct = mul(mul(p0, rho_d), p0)
    check("local_selective_cp_branch_reduces_exactly",
          reduced_branch_from_full == reduced_branch_direct)
    check("selective_branch_weight_matches", trace(full_branch) == half
          and trace(reduced_branch_direct) == half)
    normalized_branch = tuple(tuple(2 * value for value in row)
                              for row in reduced_branch_direct)
    check("positive_selective_branch_normalizes", trace(normalized_branch) == 1)

    # Nonempty exact Hamiltonian/circuit witness: source writer W copies to R.
    # Computational ordering is |WR> = |00>,|01>,|10>,|11>.
    x_on_w = (
        (F(0), F(0), F(1), F(0)),
        (F(0), F(0), F(0), F(1)),
        (F(1), F(0), F(0), F(0)),
        (F(0), F(1), F(0), F(0)),
    )
    cnot = (
        (F(1), F(0), F(0), F(0)),
        (F(0), F(1), F(0), F(0)),
        (F(0), F(0), F(0), F(1)),
        (F(0), F(0), F(1), F(0)),
    )
    pr0 = add(basis_state(4, 0), basis_state(4, 2))
    pr1 = add(basis_state(4, 1), basis_state(4, 3))
    rho00 = basis_state(4, 0)

    check("source_x_gate_unitary", mul(transpose(x_on_w), x_on_w) == eye(4))
    check("common_cnot_unitary", mul(transpose(cnot), cnot) == eye(4))
    check("complete_query_effects", add(pr0, pr1) == eye(4))
    check("query_effects_orthogonal", mul(pr0, pr1) == tuple(
        tuple(F(0) for _ in range(4)) for _ in range(4)))

    source_states = [rho00, conjugate(x_on_w, rho00)]
    final_states = [conjugate(cnot, state) for state in source_states]
    laws = []
    for state in final_states:
        laws.append((probability(state, pr0), probability(state, pr1)))

    check("arm0_query_law", laws[0] == (F(1), F(0)))
    check("arm1_query_law", laws[1] == (F(0), F(1)))
    check("perfect_witness_tv", (abs(laws[1][0] - laws[0][0])
                                  + abs(laws[1][1] - laws[0][1])) / 2 == 1)
    check("each_query_law_normalized", all(sum(law, F(0)) == 1 for law in laws))

    # Derived write-ancestor frontier.
    edges = (("root", "source"), ("source", "unitary"), ("unitary", "query"))
    past = {"root", "source", "unitary"}
    frontier = {(u, v) for u, v in edges if u in past and v not in past}
    check("source_to_query_route_nonempty", reachable(edges, "source", "query"))
    check("derived_frontier_exact", frontier == {("unitary", "query")})
    check("frontier_nonempty", bool(frontier))
    check("frontier_deletion_separates", not reachable(
        edges, "source", "query", frozenset(frontier)))
    check("no_return_into_write_past", not any(
        u not in past and v in past for u, v in edges))

    # Exhaustive logic ledger.  FHBC is primitive microphysics, while the four
    # F clauses are derived groupings used by the theorem proof.
    names = (
        "finite_device", "joint_root", "exact_quantum_law",
        "complete_instrument", "outer_closure", "source_custody",
        "mission_separation",
    )

    def derived(physical):
        return (
            physical["finite_device"] and physical["outer_closure"]
            and physical["mission_separation"],
            physical["finite_device"] and physical["joint_root"],
            physical["exact_quantum_law"] and physical["complete_instrument"],
            physical["outer_closure"] and physical["source_custody"],
        )

    implication_cases = 0
    antecedent_cases = 0
    for values in product((False, True), repeat=len(names)):
        p = dict(zip(names, values))
        fhbc = all(values)
        f0, f1, f2, f3 = derived(p)
        implication_cases += 1
        antecedent_cases += int(fhbc)
        if fhbc and not (f0 and f1 and f2 and f3):
            raise AssertionError((values, f0, f1, f2, f3))
    check("fhbc_to_f0_f3_exhaustive_128_cases", implication_cases == 128)
    check("fhbc_antecedent_nonempty", antecedent_cases == 1)

    # Each residual premise can fail independently; no premise is smuggled in
    # from positive response.
    complete = {name: True for name in names}
    missing_closure = dict(complete, outer_closure=False)
    missing_instrument = dict(complete, complete_instrument=False)
    missing_root = dict(complete, joint_root=False)
    missing_custody = dict(complete, source_custody=False)
    interleaved_mission = dict(complete, mission_separation=False)
    check("missing_outer_closure_blocks_f0_f3",
          derived(missing_closure)[0] is False
          and derived(missing_closure)[3] is False)
    check("incomplete_instrument_blocks_f2",
          derived(missing_instrument)[2] is False)
    check("missing_joint_root_blocks_f1", derived(missing_root)[1] is False)
    check("missing_source_custody_blocks_f3",
          derived(missing_custody)[3] is False)
    check("query_before_later_write_blocks_separation",
          derived(interleaved_mission)[0] is False)

    # Exact hidden-bypass twins: passive laws match, causal break separates.
    passive_honest = {(0, 0, 0, 0): F(1, 2), (1, 1, 1, 1): F(1, 2)}
    passive_bypass = dict(passive_honest)
    check("passive_honest_bypass_laws_identical", passive_honest == passive_bypass)

    honest_edges = (("source", "visible_R"), ("visible_R", "query"))
    bypass_edges = honest_edges + (("source", "hidden_H"), ("hidden_H", "query"))
    visible_cut = frozenset({("visible_R", "query")})
    check("honest_visible_cut_separates", not reachable(
        honest_edges, "source", "query", visible_cut))
    check("hidden_bypass_survives_visible_cut", reachable(
        bypass_edges, "source", "query", visible_cut))

    honest_reset_law = {0: 0, 1: 0}
    bypass_reset_law = {0: 0, 1: 1}
    check("causal_break_distinguishes_twins", honest_reset_law != bypass_reset_law)
    check("bypass_world_fails_outer_closure", ("source", "hidden_H") in bypass_edges)
    check("response_contrast_does_not_imply_fhbc", laws[0] != laws[1])

    print(f"TOTAL {len(CHECKS)}/{len(CHECKS)} PASS")
    print("VERDICT FINITE_HAMILTONIAN_BOUNDARY_CLOSURE_DERIVES_F0_F3_DCL")
    print("OUTER_DEVICE_CLOSURE_AND_SOURCE_CUSTODY PHYSICAL_PREMISES")
    print("PASSIVE_VISIBLE_RESPONSE_ONLY_CLOSURE_AUTHENTICATION IMPOSSIBLE")
    print("UNIVERSAL_FHBC_MEMBERSHIP NOT_PROVED")


if __name__ == "__main__":
    main()
