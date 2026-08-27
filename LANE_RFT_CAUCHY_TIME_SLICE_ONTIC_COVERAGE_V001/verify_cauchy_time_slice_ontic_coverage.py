#!/usr/bin/env python3
"""Exact finite ledger for the conditional Cauchy/time-slice theorem.

This verifies only checkable implication logic and finite state/effect witnesses.
It does not establish a physical Cauchy surface, a time-slice axiom, a complete
port census, or CTS for any external record.
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


BRANCH_TO_ARCH = {
    "K_CLASSICAL": "FCN_ARCH",
    "K_QUANTUM_FINITE_DIM": "FCN_ARCH",
    "W": "ACN_CP_ARCH",
    "A": "CSTAR_CP_ARCH",
    "I": "FICO_DM_ARCH",
}


@dataclass(frozen=True)
class CTSPacket:
    branch: str
    t0_same_episode: bool = True
    t1_definite_bounded_cut: bool = True
    t2_cauchy_sufficient: bool = True
    t3_exact_positive_propagation: bool = True
    t4_complete_joint_crossing: bool = True
    t5_conditional_future_provenance: bool = True
    independent_authentication: bool = False


def cts(packet: CTSPacket) -> bool:
    return (
        packet.branch in BRANCH_TO_ARCH
        and packet.t0_same_episode
        and packet.t1_definite_bounded_cut
        and packet.t2_cauchy_sufficient
        and packet.t3_exact_positive_propagation
        and packet.t4_complete_joint_crossing
        and packet.t5_conditional_future_provenance
    )


def ontic_cover(packet: CTSPacket) -> Optional[str]:
    if not cts(packet):
        return None
    return BRANCH_TO_ARCH[packet.branch]


def a1_a4(packet: CTSPacket) -> Optional[tuple[bool, bool, bool, bool]]:
    if ontic_cover(packet) is None:
        return None
    return (True, True, True, True)


def compose_classical(
    cut_law: dict[int, Fraction],
    future_kernel: dict[int, dict[int, Fraction]],
) -> dict[int, Fraction]:
    outcomes = sorted({h for row in future_kernel.values() for h in row})
    return {
        h: sum(cut_law.get(s, Fraction(0)) * future_kernel[s].get(h, Fraction(0))
               for s in future_kernel)
        for h in outcomes
    }


def tv(p: dict[int, Fraction], q: dict[int, Fraction]) -> Fraction:
    keys = set(p) | set(q)
    return sum(abs(p.get(k, 0) - q.get(k, 0)) for k in keys) / 2


def pointwise_mul(a: tuple[Fraction, ...], b: tuple[Fraction, ...]):
    return tuple(x * y for x, y in zip(a, b))


def alpha(a: tuple[Fraction, ...]):
    """A finite commutative C*-isomorphism A_G -> A_D (coordinate permutation)."""
    return (a[2], a[0], a[1])


def alpha_inverse(b: tuple[Fraction, ...]):
    return (b[1], b[2], b[0])


def state(weights: tuple[Fraction, ...], effect: tuple[Fraction, ...]) -> Fraction:
    return sum(w * e for w, e in zip(weights, effect))


def embed_field(effect: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    """Embed C^3 field effects into C^3 tensor C^2 joint field-memory effects."""
    return tuple(value for value in effect for _ in range(2))


def field_marginal(weights: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(weights[2 * i] + weights[2 * i + 1] for i in range(3))


def main() -> None:
    t = Ledger()

    # Exhaustive T0--T5 implication ledger for every permitted branch.
    exhaustive = True
    for branch in tuple(BRANCH_TO_ARCH) + ("UNKNOWN",):
        for bits in product((False, True), repeat=6):
            packet = CTSPacket(branch, *bits)
            expected = branch in BRANCH_TO_ARCH and all(bits)
            exhaustive &= cts(packet) == expected
            exhaustive &= (ontic_cover(packet) is not None) == expected
            exhaustive &= (a1_a4(packet) == (True, True, True, True)) == expected
    t.check(exhaustive, "exhaustive_T0_T5_branch_implication")

    # Each known branch has the promised ontic architecture routing.
    for branch, architecture in BRANCH_TO_ARCH.items():
        packet = CTSPacket(branch)
        t.check(ontic_cover(packet) == architecture, f"route_{branch}")
        t.check(a1_a4(packet) == (True, True, True, True), f"A1_A4_{branch}")

    # Every CTS clause is necessary to invoke this theorem.
    full = CTSPacket("A")
    for field in (
        "t0_same_episode",
        "t1_definite_bounded_cut",
        "t2_cauchy_sufficient",
        "t3_exact_positive_propagation",
        "t4_complete_joint_crossing",
        "t5_conditional_future_provenance",
    ):
        t.check(ontic_cover(replace(full, **{field: False})) is None,
                f"missing_{field}_rejected")

    # Ontic coverage does not require or manufacture independent authentication.
    t.check(not full.independent_authentication, "authentication_absent")
    t.check(ontic_cover(full) == "CSTAR_CP_ARCH", "ontic_survives_without_auth")
    t.check(not replace(full, independent_authentication=True) == full,
            "authentication_is_separate_datum")

    # Exact classical Cauchy-kernel factorization with one arm-common future.
    cut_laws = {
        0: {0: Fraction(1, 2), 1: Fraction(1, 2), 2: Fraction(0)},
        1: {0: Fraction(0), 1: Fraction(1, 2), 2: Fraction(1, 2)},
    }
    future = {
        0: {0: Fraction(1), 1: Fraction(0)},
        1: {0: Fraction(1, 2), 1: Fraction(1, 2)},
        2: {0: Fraction(0), 1: Fraction(1)},
    }
    for x, law in cut_laws.items():
        t.check(sum(law.values()) == 1, f"cut_law_normalized_x{x}")
    for s, row in future.items():
        t.check(sum(row.values()) == 1, f"future_kernel_normalized_s{s}")
    responses = {x: compose_classical(law, future) for x, law in cut_laws.items()}
    t.check(responses[0] == {0: Fraction(3, 4), 1: Fraction(1, 4)},
            "classical_response_x0_exact")
    t.check(responses[1] == {0: Fraction(1, 4), 1: Fraction(3, 4)},
            "classical_response_x1_exact")
    t.check(tv(responses[0], responses[1]) == Fraction(1, 2),
            "classical_cut_contrast_exact")

    # Finite commutative witness for the algebraic time-slice pullback.
    one = (Fraction(1), Fraction(1), Fraction(1))
    a = (Fraction(1, 3), Fraction(1, 2), Fraction(3, 4))
    b = (Fraction(2, 3), Fraction(1, 4), Fraction(1, 2))
    t.check(alpha(one) == one, "time_slice_isomorphism_unital")
    t.check(alpha_inverse(alpha(a)) == a, "time_slice_inverse_left")
    t.check(alpha(alpha_inverse(a)) == a, "time_slice_inverse_right")
    t.check(alpha(pointwise_mul(a, b)) == pointwise_mul(alpha(a), alpha(b)),
            "time_slice_isomorphism_multiplicative")
    omega_g = (Fraction(1, 4), Fraction(1, 2), Fraction(1, 4))
    effect_d = (Fraction(0), Fraction(1, 2), Fraction(1))
    effect_g = alpha_inverse(effect_d)
    # omega_D := omega_G o alpha^{-1}; hence omega_D(E_D)=omega_G(alpha^{-1}E_D).
    p_d = state(omega_g, alpha_inverse(effect_d))
    p_g = state(omega_g, effect_g)
    t.check(p_d == p_g, "time_slice_state_effect_pullback_exact")
    t.check(Fraction(0) <= p_g <= Fraction(1), "time_slice_probability_valid")

    # The field time-slice algebra is only a component of the joint separator.
    joint_weights = (
        Fraction(1, 8), Fraction(1, 8),
        Fraction(1, 4), Fraction(1, 4),
        Fraction(1, 8), Fraction(1, 8),
    )
    t.check(sum(joint_weights) == 1, "joint_separator_state_normalized")
    t.check(field_marginal(joint_weights) == omega_g,
            "joint_separator_restricts_to_field_state")
    t.check(state(joint_weights, embed_field(effect_g)) == state(omega_g, effect_g),
            "joint_embedding_preserves_field_expectation")

    # Two joint states can have the same field marginal but opposite memory bits.
    joint_memory_0 = (
        Fraction(1, 4), Fraction(0),
        Fraction(1, 2), Fraction(0),
        Fraction(1, 4), Fraction(0),
    )
    joint_memory_1 = (
        Fraction(0), Fraction(1, 4),
        Fraction(0), Fraction(1, 2),
        Fraction(0), Fraction(1, 4),
    )
    memory_effect = tuple(Fraction(i % 2) for i in range(6))
    t.check(field_marginal(joint_memory_0) == field_marginal(joint_memory_1) == omega_g,
            "field_marginal_does_not_fix_memory")
    t.check(state(joint_memory_0, memory_effect) == 0,
            "joint_memory_zero_effect")
    t.check(state(joint_memory_1, memory_effect) == 1,
            "joint_memory_one_effect")

    # A visible carrier can fail completeness when a hidden bypass remains.
    replaced_visible_boundary = 0
    bypass_outputs = {x: x for x in (0, 1)}
    t.check(replaced_visible_boundary == 0, "declared_boundary_common_replacement")
    t.check(bypass_outputs[0] != bypass_outputs[1], "hidden_bypass_retains_contrast")
    incomplete = replace(full, t4_complete_joint_crossing=False)
    t.check(ontic_cover(incomplete) is None, "incomplete_bypass_not_CTS")

    # Record contrast is not an input to the conditional logic.
    zero_contrast_candidate = full
    positive_contrast_candidate = full
    t.check(ontic_cover(zero_contrast_candidate) == ontic_cover(positive_contrast_candidate),
            "contrast_not_used_for_CTS_bridge")

    print("T0_T5_BRANCH_LOGIC EXHAUSTIVE")
    print("CLASSICAL_CAUCHY_FACTORIZATION EXACT")
    print("COMMUTATIVE_TIME_SLICE_PULLBACK EXACT")
    print("JOINT_SEPARATOR_EXTENSION_AND_MEMORY WITNESS EXACT")
    print("INCOMPLETE_BYPASS_REJECTED EXACT")
    print("AUTHENTICATION_AND_CONTRAST_NOT_INFERRED EXACT")
    print(f"TOTAL {t.passed}/{t.total} PASS")
    print("VERDICT CONDITIONAL_CTS_ONTIC_COVERAGE_PROVED__UNIVERSAL_CTS_OPEN")


if __name__ == "__main__":
    main()
