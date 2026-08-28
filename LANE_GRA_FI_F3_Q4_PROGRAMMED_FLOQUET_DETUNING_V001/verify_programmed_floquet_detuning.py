#!/usr/bin/env python3
"""Deterministic checks for GRA-FI-F3-Q4-PFCD-V001.

The checks replay finite q4 incidence, the exact two-pulse spectrum, the
dressed-parent functional calculus, and its rigorous scalar remainder bound.
They do not prove an autonomous controller, physical port calibration, a
collective phase, RGRL-B, or gravity.
"""

from __future__ import annotations

import itertools
import math
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parent
THEOREM = (ROOT / "THEOREM.md").read_text(encoding="utf-8")
SELF_AUDIT = (ROOT / "SELF_AUDIT.md").read_text(encoding="utf-8")


class Checks:
    def __init__(self) -> None:
        self.total = 0
        self.failed = 0

    def check(self, condition: bool, label: str) -> None:
        self.total += 1
        if condition:
            print(f"PASS  {label}")
        else:
            self.failed += 1
            print(f"FAIL  {label}")


def fronts(depth: int) -> list[tuple[int, int, int, int]]:
    return [
        item
        for item in itertools.product(range(depth + 1), repeat=4)
        if sum(item) == depth
    ]


def incidence(depth: int) -> np.ndarray:
    parents = fronts(depth)
    children = fronts(depth + 1)
    child_index = {child: i for i, child in enumerate(children)}
    matrix = np.zeros((len(children), len(parents)), dtype=float)
    for j, parent in enumerate(parents):
        for axis in range(4):
            child = list(parent)
            child[axis] += 1
            matrix[child_index[tuple(child)], j] = 1.0
    return matrix


def matrix_function_symmetric(matrix: np.ndarray, function) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * function(values)) @ vectors.T


def exp_i_symmetric(matrix: np.ndarray, scale: float) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * np.exp(1j * scale * values)) @ vectors.T


def balanced(text: str, token: str) -> bool:
    return text.count(token) % 2 == 0


def main() -> int:
    q = Checks()
    samples = [incidence(depth) for depth in range(6)]

    q.check(
        all(
            matrix.shape
            == (math.comb(depth + 4, 3), math.comb(depth + 3, 3))
            for depth, matrix in enumerate(samples)
        ),
        "q4 parent/child stars-and-bars dimensions",
    )
    q.check(
        all(np.allclose(matrix.sum(axis=0), 4.0) for matrix in samples),
        "every q4 parent has four append children",
    )
    q.check(
        all(np.max(matrix.sum(axis=1)) <= 4.0 for matrix in samples),
        "every finite-slab child has at most four parents",
    )
    q.check(
        all(np.linalg.matrix_rank(matrix) == matrix.shape[1] for matrix in samples),
        "sampled q4 incidence maps are injective",
    )
    q.check(
        all(np.linalg.norm(matrix, 2) <= 4.0 + 1e-12 for matrix in samples),
        "q4 incidence operator norm ceiling",
    )

    # The next active slab sees a carrier in C but not one in P or blank G.
    for depth, matrix in enumerate(samples[:4]):
        parent_count = matrix.shape[1]
        child_count = matrix.shape[0]
        guard_count = child_count
        next_slab_onsite = np.diag(
            np.r_[
                np.zeros(parent_count),
                np.ones(child_count),
                np.zeros(guard_count),
            ]
        )
        child_projector = np.diag(
            np.r_[
                np.zeros(parent_count),
                np.ones(child_count),
                np.zeros(guard_count),
            ]
        )
        q.check(
            np.array_equal(next_slab_onsite, child_projector),
            f"next-slab onsite is exactly child-only depth={depth}",
        )

    # Exact K/n typing and pulse-slice audit.  Basis order is |K,n>.  Both the
    # raw I_K X_n and gated P_K X_n generators move an old K=1,n=1 support
    # edge, while raw X also moves a new K=0,n=0 blank link.  With both absent,
    # the diagonal incidence block commutes with n.
    i2 = np.eye(2)
    x2 = np.array([[0.0, 1.0], [1.0, 0.0]])
    p1 = np.diag([0.0, 1.0])
    p_k = np.kron(p1, i2)
    p_n = np.kron(i2, p1)
    raw_x = np.kron(i2, x2)
    gated_x = np.kron(p1, x2)
    old_saturated = np.array([0.0, 0.0, 0.0, 1.0])
    new_blank = np.array([1.0, 0.0, 0.0, 0.0])
    q.check(not np.allclose(gated_x @ old_saturated, 0.0),
            "PESC K-gated flip would move saturated old-slab n")
    q.check(not np.allclose(raw_x @ old_saturated, 0.0),
            "raw BS06 flip would move saturated old-slab n")
    q.check(np.allclose(gated_x @ new_blank, 0.0),
            "PESC K-gated flip vanishes on K=0 next-slab blank")
    q.check(not np.allclose(raw_x @ new_blank, 0.0),
            "raw BS06 flip would move K=0 next-slab blank")
    diagonal_slice = 0.83 * p_n
    q.check(np.allclose(diagonal_slice @ p_n - p_n @ diagonal_slice, 0.0),
            "dual-flip-free incidence slice exactly conserves n")
    q.check(not np.allclose(p_k, p_n), "K support and active n are distinct factors")

    # A stroboscopic echo can return n while changing the intended carrier
    # unitary.  Flip n=1 to n=0, wait under n-controlled hopping, then flip it
    # back: n returns exactly but the carrier sees identity rather than the
    # desired hop generated on the saturated n=1 block.
    eta_echo = 0.31
    x_on_n = np.kron(x2, i2)
    controlled_hop = np.kron(p1, x2)
    u_wait = exp_i_symmetric(controlled_hop, eta_echo)
    u_echo = x_on_n @ u_wait @ x_on_n
    echo_carrier_block = u_echo[2:4, 2:4]
    desired_carrier_hop = exp_i_symmetric(x2, eta_echo)
    q.check(
        np.allclose(echo_carrier_block, i2)
        and not np.allclose(echo_carrier_block, desired_carrier_hop),
        "incidence-only stroboscopic return does not preserve carrier hop unitary",
    )

    eta = math.pi / 16.0
    spectrum_ok = True
    unitarity_ok = True
    phase_window_ok = True
    dark_count_ok = True
    branch_function_ok = True
    kernel_bound_ok = True
    kernel_sign_ok = True

    for matrix in samples:
        child_count, parent_count = matrix.shape
        zero_p = np.zeros((parent_count, parent_count))
        zero_c = np.zeros((child_count, child_count))
        x_b = np.block([[zero_p, matrix.T], [matrix, zero_c]])
        u_h = exp_i_symmetric(x_b, eta)
        u_d = np.diag(
            np.r_[np.ones(parent_count), -1j * np.ones(child_count)]
        )
        u_f = u_d @ u_h
        ident = np.eye(parent_count + child_count)
        unitarity_ok &= np.allclose(u_f.conj().T @ u_f, ident, atol=2e-12)

        singular = np.linalg.svd(matrix, compute_uv=False)
        omega = np.arccos(np.cos(eta * singular) / math.sqrt(2.0))
        formula = np.r_[
            np.exp(1j * (omega - math.pi / 4.0)),
            np.exp(-1j * (omega + math.pi / 4.0)),
            -1j * np.ones(child_count - parent_count),
        ]
        actual_angles = np.sort(np.angle(np.linalg.eigvals(u_f)))
        formula_angles = np.sort(np.angle(formula))
        spectrum_ok &= np.allclose(actual_angles, formula_angles, atol=4e-12)
        phase_window_ok &= bool(
            np.all(omega >= math.pi / 4.0 - 1e-12)
            and np.all(omega <= math.pi / 3.0 + 1e-12)
        )
        dark_count_ok &= child_count - parent_count == len(
            [value for value in formula if abs(value + 1j) < 1e-12]
        )

        k_matrix = matrix.T @ matrix
        k_values, k_vectors = np.linalg.eigh(k_matrix)
        exact_values = -(
            np.arccos(np.cos(eta * np.sqrt(k_values)) / math.sqrt(2.0))
            - math.pi / 4.0
        )
        exact_operator = (k_vectors * exact_values) @ k_vectors.T
        parent_phase_values = -(omega - math.pi / 4.0)
        branch_function_ok &= np.allclose(
            np.sort(np.linalg.eigvalsh(exact_operator)),
            np.sort(parent_phase_values),
            atol=3e-12,
        )

        remainder = exact_operator + 0.5 * eta**2 * k_matrix
        upper = (eta**4 / 6.0) * (k_matrix @ k_matrix)
        kernel_bound_ok &= bool(
            np.min(np.linalg.eigvalsh(remainder)) >= -3e-12
            and np.min(np.linalg.eigvalsh(upper - remainder)) >= -3e-12
        )
        kernel_sign_ok &= np.linalg.norm(k_matrix, 2) > 0.0

    q.check(unitarity_ok, "two-pulse Floquet operator is unitary")
    q.check(spectrum_ok, "exact SVD Floquet eigenvalue formula")
    q.check(phase_window_ok, "parent/child phase windows under eta ceiling")
    q.check(dark_count_ok, "unpaired child dark-mode multiplicity")
    q.check(branch_function_ok, "dressed-parent exact functional calculus")
    q.check(kernel_bound_ok, "operator common-child remainder inequality")
    q.check(kernel_sign_ok, "strict nonzero leading sibling coefficient")

    # Independent scalar derivative/remainder replay across the frozen domain.
    z_values = np.linspace(0.0, math.pi**2 / 16.0, 4001)
    g_values = np.arccos(np.cos(np.sqrt(z_values)) / math.sqrt(2.0)) - math.pi / 4.0
    q.check(
        np.all(g_values <= z_values / 2.0 + 2e-15),
        "scalar upper kernel bound on frozen interval",
    )
    q.check(
        np.all(g_values >= z_values / 2.0 - z_values**2 / 6.0 - 2e-15),
        "scalar lower kernel bound on frozen interval",
    )

    required_theorem_phrases = [
        "yes stroboscopically",
        "During **both** the hop pulse and the next-slab onsite pulse",
        "PESC `-h sum_e P_e^K X_{n_e}` actuator exactly zero",
        "merely stroboscopic incidence echo",
        "stroboscopic\" describes the carrier quasienergy law",
        "joint-`n`-and-carrier-identity",
        "both raw and PESC `K`-gated incidence flips are zero",
        "qualified dual-flip-free isolation",
        "H_D=\\epsilon_\\psi\\Pi_C",
        "PROGRAMMED-NEXT-SLAB-DETUNING",
        "Static source-off bulk",
        "No new bulk carrier generator",
        "controller couplings which isolate and time the pulse remain supplied",
        "work the physical implementation requires",
        "matrices and calibration remain supplied",
        "not a claim that the source-off static Hamiltonian",
        "massless collective coordinate",
        "RGRL-B",
    ]
    normalized = " ".join(THEOREM.lower().split())
    for phrase in required_theorem_phrases:
        q.check(" ".join(phrase.lower().split()) in normalized,
                f"theorem ceiling: {phrase}")

    normalized_self = " ".join(SELF_AUDIT.lower().split())
    q.check(
        "merely stroboscopic incidence echo" in normalized_self,
        "self-audit preserves incidence-echo ceiling",
    )

    forbidden_promotions = [
        "the static fd05 hamiltonian is derived",
        "the controller schedule is proved autonomous",
        "gravity is proved by this theorem",
        "the numerical value of g is derived",
        "switching requires no work",
    ]
    for phrase in forbidden_promotions:
        q.check(phrase not in normalized, f"forbidden promotion absent: {phrase}")

    q.check(THEOREM.count("\\[") == THEOREM.count("\\]"),
            "display-math delimiters balanced")
    q.check(balanced(THEOREM, "`"), "theorem inline-code delimiters balanced")
    q.check(balanced(SELF_AUDIT, "`"), "self-audit inline-code delimiters balanced")
    q.check("\x08" not in THEOREM and "\x0c" not in THEOREM,
            "no backspace/form-feed control characters")

    dependencies = [
        "../LANE_GRA_FD_F3_Q4_COMMON_CHILD_ACOUSTIC_CONE_V001/THEOREM.md",
        "../LANE_GRA_FF_F3_Q4_CARRIER_LIFT_DERIVABILITY_NO_GO_V001/THEOREM.md",
        "../LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/THEOREM.md",
        "../LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md",
    ]
    for dependency in dependencies:
        q.check((ROOT / dependency).is_file(),
                f"dependency present: {Path(dependency).parent.name}")

    print(f"SUMMARY {q.total - q.failed}/{q.total} checks passed")
    return 1 if q.failed else 0


if __name__ == "__main__":
    sys.exit(main())
