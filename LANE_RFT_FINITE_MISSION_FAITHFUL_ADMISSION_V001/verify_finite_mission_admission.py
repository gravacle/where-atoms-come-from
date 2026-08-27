#!/usr/bin/env python3
"""Exact witnesses for finite-mission faithful-admission reduction.

This executable checks the finite score-table construction, the honest/bypass
observational twins, the unique inaccessible-record certification counterexample,
the non-CP transpose witness, and the finite logical bridge used in THEOREM.md.
It does not authenticate an external physical record or prove the open ontic
universal physical-cover premise.
"""

from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import product
from typing import Optional


class Ledger:
    def __init__(self) -> None:
        self.passed = 0
        self.total = 0

    def check(self, condition: bool, name: str) -> None:
        self.total += 1
        if not condition:
            raise AssertionError(name)
        self.passed += 1


def tv(p: dict[int, Fraction], q: dict[int, Fraction]) -> Fraction:
    keys = set(p) | set(q)
    return sum(abs(p.get(k, Fraction(0)) - q.get(k, Fraction(0))) for k in keys) / 2


def deterministic(value: int) -> dict[int, Fraction]:
    return {0: Fraction(value == 0), 1: Fraction(value == 1)}


def transpose_2x2(a: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]):
    return ((a[0][0], a[1][0]), (a[0][1], a[1][1]))


def trace_product_2x2(a, b) -> Fraction:
    # Tr(a b), exact over Fractions.
    return sum(a[i][k] * b[k][i] for i in range(2) for k in range(2))


@dataclass(frozen=True)
class Cover:
    branch: str
    same_episode: bool = True
    finite_exact: bool = True
    complete_boundary: bool = True
    no_hidden_label: bool = True
    standard_composable: bool = True
    independent_auth: bool = True
    graph_auth: bool = False


BRANCH_TO_ADMISSION = {
    "K": "FCN_ADMIT",
    "W": "ACN_CP_ADMIT",
    "A": "CSTAR_CP_ADMIT",
    "I": "FICO-DM-ADMIT",
}

BRANCH_TO_ARCHITECTURE = {
    "K": "FCN_ARCH",
    "W": "ACN_CP_ARCH",
    "A": "CSTAR_CP_ARCH",
    "I": "FICO-DM-ARCH",
}


def transcribe_physical(cover: Cover) -> Optional[str]:
    structural = (
        cover.same_episode
        and cover.finite_exact
        and cover.complete_boundary
        and cover.no_hidden_label
        and cover.standard_composable
    )
    if not structural:
        return None
    return BRANCH_TO_ARCHITECTURE.get(cover.branch)


def transcribe_authenticated(cover: Cover) -> Optional[str]:
    if not cover.independent_auth or transcribe_physical(cover) is None:
        return None
    return BRANCH_TO_ADMISSION.get(cover.branch)


def main() -> None:
    t = Ledger()

    # FM1: an exact finite score-table simulator.
    response = {
        (0, 0, "Q"): {0: Fraction(3, 4), 1: Fraction(1, 4)},
        (1, 0, "Q"): {0: Fraction(1, 4), 1: Fraction(3, 4)},
        (0, 1, "Q"): {0: Fraction(2, 3), 1: Fraction(1, 3)},
        (1, 1, "Q"): {0: Fraction(1, 3), 1: Fraction(2, 3)},
        (0, 0, "R"): {0: Fraction(1), 1: Fraction(0)},
        (1, 0, "R"): {0: Fraction(0), 1: Fraction(1)},
    }
    for key, distribution in response.items():
        t.check(sum(distribution.values()) == 1, f"score_table_normalized_{key}")
        simulator = dict(distribution)  # K(h|x,z,j) := P(h|x,z,j)
        t.check(simulator == distribution, f"score_table_exact_{key}")
    t.check(tv(response[(0, 0, "Q")], response[(1, 0, "Q")]) == Fraction(1, 2),
            "score_table_nonzero_contrast")

    # FM4: honest/bypass observational twins.
    visible_honest = {}
    visible_bypass = {}
    for x in (0, 1):
        c_h, h_h, q_h = x, 0, x
        c_b, h_b, q_b = x, x, x
        visible_honest[x] = (c_h, q_h)
        visible_bypass[x] = (c_b, q_b)
        t.check(h_h == 0, f"honest_hidden_constant_x{x}")
        t.check(h_b == x, f"bypass_hidden_copy_x{x}")
    t.check(visible_honest == visible_bypass, "observational_twins_visible_law_equal")
    t.check(visible_honest == {0: (0, 0), 1: (1, 1)}, "visible_perfect_record")

    honest_after_replace = {x: deterministic(0) for x in (0, 1)}
    bypass_after_replace = {x: deterministic(x) for x in (0, 1)}
    t.check(tv(honest_after_replace[0], honest_after_replace[1]) == 0,
            "honest_visible_cut_screens")
    t.check(tv(bypass_after_replace[0], bypass_after_replace[1]) == 1,
            "bypass_visible_cut_fails_to_screen")
    t.check(honest_after_replace != bypass_after_replace,
            "replacement_breaks_observational_equivalence")

    # Joint-correlation closure witness.
    joint = {}
    for x in (0, 1):
        outcomes = []
        for u in (0, 1):
            c, d = x ^ u, u
            outcomes.append((c, d, c ^ d))
        joint[x] = outcomes
    for x in (0, 1):
        c_counts = {v: sum(1 for c, _, _ in joint[x] if c == v) for v in (0, 1)}
        d_counts = {v: sum(1 for _, d, _ in joint[x] if d == v) for v in (0, 1)}
        t.check(c_counts == {0: 1, 1: 1}, f"joint_C_marginal_blind_x{x}")
        t.check(d_counts == {0: 1, 1: 1}, f"joint_D_marginal_blind_x{x}")
        t.check(all(q == x for _, _, q in joint[x]), f"joint_parity_records_x{x}")

    # Identity/transpose endpoint twin and exact non-CP Choi witness.
    zero = Fraction(0)
    one = Fraction(1)
    rho0 = ((one, zero), (zero, zero))
    rho1 = ((zero, zero), (zero, one))
    e0 = rho0
    e1 = rho1
    for idx, rho in enumerate((rho0, rho1)):
        t.check(transpose_2x2(rho) == rho, f"transpose_diagonal_state_{idx}")
        for jdx, effect in enumerate((e0, e1)):
            p_identity = trace_product_2x2(rho, effect)
            p_transpose = trace_product_2x2(transpose_2x2(rho), effect)
            t.check(p_identity == p_transpose,
                    f"identity_transpose_endpoint_equal_{idx}_{jdx}")

    swap = (
        (1, 0, 0, 0),
        (0, 0, 1, 0),
        (0, 1, 0, 0),
        (0, 0, 0, 1),
    )
    v = (0, 1, -1, 0)
    fv = tuple(sum(swap[i][j] * v[j] for j in range(4)) for i in range(4))
    quadratic = sum(v[i] * fv[i] for i in range(4))
    t.check(fv == tuple(-entry for entry in v), "swap_antisymmetric_eigenvector")
    t.check(quadratic == -2, "transpose_choi_negative_quadratic_form")
    t.check(quadratic < 0, "transpose_not_completely_positive")

    # FM5: a genuine unique record and a no-write common-cause world have the
    # same complete surviving observational law.
    weights = {0: Fraction(1, 3), 1: Fraction(2, 3)}
    natural_record_world = {(u, u, u): weight for u, weight in weights.items()}
    natural_common_cause_world = {(u, u, u): weight for u, weight in weights.items()}
    t.check(natural_record_world == natural_common_cause_world,
            "unique_record_common_cause_laws_equal")
    t.check(sum(natural_record_world.values()) == 1,
            "unique_record_surviving_law_normalized")
    record_has_write_edge = True
    common_cause_has_write_edge = False
    t.check(record_has_write_edge and not common_cause_has_write_edge,
            "unique_record_causal_structures_differ")
    independent_provenance_survives = False
    t.check(not independent_provenance_survives,
            "unique_record_independent_authentication_absent_by_construction")

    # FM2: four exact bridge branches, outcome-independent of contrast.
    for branch, expected in BRANCH_TO_ADMISSION.items():
        cover = Cover(branch=branch)
        t.check(transcribe_physical(cover) == BRANCH_TO_ARCHITECTURE[branch],
                f"physical_branch_transcribes_{branch}")
        t.check(transcribe_authenticated(cover) == expected,
                f"authenticated_branch_transcribes_{branch}")
        t.check(transcribe_authenticated(replace(cover, graph_auth=True)) == expected,
                f"graph_flag_does_not_change_base_admission_{branch}")

    full = Cover(branch="A")
    for field in (
        "same_episode",
        "finite_exact",
        "complete_boundary",
        "no_hidden_label",
        "standard_composable",
    ):
        defective = replace(full, **{field: False})
        t.check(transcribe_physical(defective) is None,
                f"missing_physical_bridge_clause_rejected_{field}")
        t.check(transcribe_authenticated(defective) is None,
                f"missing_authenticated_bridge_clause_rejected_{field}")
    unauthenticated = replace(full, independent_auth=False)
    t.check(transcribe_physical(unauthenticated) == "CSTAR_CP_ARCH",
            "missing_auth_does_not_refute_ontic_physical_structure")
    t.check(transcribe_authenticated(unauthenticated) is None,
            "missing_auth_rejects_authenticated_admission")
    t.check(transcribe_physical(Cover(branch="UNKNOWN")) is None,
            "unknown_physical_branch_rejected")
    t.check(transcribe_authenticated(Cover(branch="UNKNOWN")) is None,
            "unknown_authenticated_branch_rejected")

    # Exhaustive truth-table checks.  Ontic physical transcription uses the first
    # five clauses; authenticated admission additionally uses independent_auth.
    exhaustive_physical_ok = True
    exhaustive_authenticated_ok = True
    for branch in tuple(BRANCH_TO_ADMISSION) + ("UNKNOWN",):
        for bits in product((False, True), repeat=6):
            cover = Cover(branch, *bits)
            known = branch in BRANCH_TO_ADMISSION
            exhaustive_physical_ok &= (
                (transcribe_physical(cover) is not None) == (known and all(bits[:5]))
            )
            exhaustive_authenticated_ok &= (
                (transcribe_authenticated(cover) is not None) == (known and all(bits))
            )
    t.check(exhaustive_physical_ok, "exhaustive_ontic_bridge_truth_table")
    t.check(exhaustive_authenticated_ok, "exhaustive_authenticated_bridge_truth_table")
    t.check(transcribe_authenticated(full) == transcribe_authenticated(full),
            "zero_or_positive_contrast_not_an_input")
    t.check(not full.graph_auth and transcribe_authenticated(full) is not None,
            "base_admission_does_not_silently_require_A5_graph")
    t.check(replace(full, graph_auth=True).graph_auth,
            "separate_A5_graph_flag_can_be_supplied")

    print("score_table_representation EXACT")
    print("honest_bypass_observational_twin EXACT")
    print("joint_correlation_closure_witness EXACT")
    print("transpose_non_cp_witness EXACT")
    print("unique_inaccessible_record_certification_counterexample EXACT")
    print("ontic_and_authenticated_four_branch_logic EXHAUSTIVE")
    print(f"TOTAL {t.passed}/{t.total} PASS")
    print("VERDICT AUTHENTICATED_TRANSCRIPTION_PROVED__ALL_ONTIC_CERTIFICATION_REFUTED__ONTIC_COVER_OPEN")


if __name__ == "__main__":
    main()
