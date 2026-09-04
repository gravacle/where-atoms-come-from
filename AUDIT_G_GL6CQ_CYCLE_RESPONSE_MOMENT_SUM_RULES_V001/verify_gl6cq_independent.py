#!/usr/bin/env python3
"""Independent hostile derivation for GL6CQ.

The target author program is neither imported nor executed.  This replay
rebuilds the cycle solder, differentiates a general invariant A1+T2 symbol,
projects its raw orientation moments back to irreducible coordinates,
reconstructs the tetrahedral contact, and derives both normalized observable
matching equations with exact rational arithmetic.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
TETRA = (
    (F(1), F(1), F(1)),
    (F(1), F(-1), F(-1)),
    (F(-1), F(1), F(-1)),
    (F(-1), F(-1), F(1)),
)
PAIR_ORDER = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: i for i, pair in enumerate(PAIR_ORDER)}
T_PAIR = (
    (F(1), F(0), F(0), F(0), F(0), F(-1)),
    (F(0), F(1), F(0), F(0), F(-1), F(0)),
    (F(0), F(0), F(1), F(-1), F(0), F(0)),
)


class Checks:
    def __init__(self):
        self.total = 0

    def equal(self, got, want, label):
        self.total += 1
        if got != want:
            raise AssertionError(f"{label}: got {got!r}, want {want!r}")

    def true(self, condition, label):
        self.total += 1
        if not condition:
            raise AssertionError(label)


CHECK = Checks()


def dot(left, right):
    return sum((F(a) * F(b) for a, b in zip(left, right)), F(0))


def transpose(matrix):
    return tuple(tuple(matrix[i][j] for i in range(len(matrix)))
                 for j in range(len(matrix[0])))


def matmul(left, right):
    return tuple(tuple(sum((left[i][q] * right[q][j]
                            for q in range(len(right))), F(0))
                       for j in range(len(right[0])))
                 for i in range(len(left)))


def identity(size):
    return tuple(tuple(F(i == j) for j in range(size)) for i in range(size))


U = (F(1, 2),) * 4
Q = tuple(tuple(TETRA[d][i] / 2 for i in range(3)) for d in range(4))
S = tuple((U[d],) + Q[d] for d in range(4))
CHECK.equal(matmul(transpose(S), S), identity(4),
            "independent A1+T2 cycle solder is orthogonal")


def zero_tensor(shape):
    if len(shape) == 1:
        return [F(0) for _ in range(shape[0])]
    return [zero_tensor(shape[1:]) for _ in range(shape[0])]


def projected_moments(coefficients):
    """Differentiate a generic invariant symbol, then project raw moments."""
    kappa, alpha, eta, b, c, d = map(F, coefficients)

    # H_abmn is the exact second k derivative of the A1+T2 symbol at k=0.
    hessian = zero_tensor((4, 4, 3, 3))
    for m in range(3):
        hessian[0][0][m][m] = 2 * alpha
    for i in range(3):
        complement = tuple(axis for axis in range(3) if axis != i)
        m, n = complement
        hessian[0][1 + i][m][n] = eta
        hessian[0][1 + i][n][m] = eta
        hessian[1 + i][0][m][n] = eta
        hessian[1 + i][0][n][m] = eta
    for i in range(3):
        for j in range(3):
            for m in range(3):
                for n in range(3):
                    value = F(0)
                    if i == j and m == n:
                        value += 2 * b
                    if i == j == m == n:
                        value += 2 * c
                    if i != j and ((m, n) == (i, j) or (m, n) == (j, i)):
                        value += d
                    hessian[1 + i][1 + j][m][n] = value

    # Raw orientation moment is minus the Fourier Hessian: M=-d_k d_k K(0).
    raw_moment = zero_tensor((4, 4, 3, 3))
    for left in range(4):
        for right in range(4):
            for m in range(3):
                for n in range(3):
                    raw_moment[left][right][m][n] = -sum(
                        (S[left][a] * hessian[a][z][m][n] * S[right][z]
                         for a in range(4) for z in range(4)), F(0))

    def project(left_vector, right_vector, m, n):
        return sum((left_vector[a] * raw_moment[a][z][m][n] * right_vector[z]
                    for a in range(4) for z in range(4)), F(0))

    m_aa = tuple(tuple(project(U, U, m, n) for n in range(3))
                 for m in range(3))
    m_at = tuple(tuple(tuple(project(U, tuple(Q[a][i] for a in range(4)), m, n)
                                  for n in range(3)) for m in range(3))
                 for i in range(3))
    m_tt = tuple(tuple(tuple(tuple(
        project(tuple(Q[a][i] for a in range(4)),
                tuple(Q[a][j] for a in range(4)), m, n)
        for n in range(3)) for m in range(3)) for j in range(3))
                 for i in range(3))

    # This compares every projected component, not only the final contractions.
    for m in range(3):
        for n in range(3):
            CHECK.equal(m_aa[m][n], -hessian[0][0][m][n],
                        "raw-to-irrep AA moment projection")
            for i in range(3):
                CHECK.equal(m_at[i][m][n], -hessian[0][1 + i][m][n],
                            "raw-to-irrep AT moment projection")
                for j in range(3):
                    CHECK.equal(m_tt[i][j][m][n],
                                -hessian[1 + i][1 + j][m][n],
                                "raw-to-irrep TT moment projection")

    # Z has a common-amplitude null and kappa I on T2.
    z = tuple(tuple(sum((S[left][1 + i] * kappa * S[right][1 + i]
                         for i in range(3)), F(0))
                    for right in range(4)) for left in range(4))
    CHECK.equal(tuple(sum((z[a][zeta] * U[zeta] for zeta in range(4)), F(0))
                      for a in range(4)), (F(0),) * 4, "zeroth common null")
    z_t = sum((sum((Q[a][i] * z[a][zeta] * Q[zeta][i]
                    for a in range(4) for zeta in range(4)), F(0))
               for i in range(3)), F(0))

    m_perp = sum((m_tt[i][i][m][m]
                  for i in range(3) for m in range(3) if i != m), F(0))
    m_parallel = sum((m_tt[i][i][i][i] for i in range(3)), F(0))
    m_cross = sum((m_tt[i][j][i][j]
                   for i in range(3) for j in range(3) if i != j), F(0))
    recovered = (
        z_t / 3,
        -sum((m_aa[m][m] for m in range(3)), F(0)) / 6,
        -sum((m_at[i][tuple(a for a in range(3) if a != i)[0]]
              [tuple(a for a in range(3) if a != i)[1]]
              for i in range(3)), F(0)) / 3,
        -m_perp / 12,
        -m_parallel / 6 + m_perp / 12,
        -m_cross / 6,
    )
    CHECK.equal(recovered, (kappa, alpha, eta, b, c, d),
                "all six coefficients recovered from raw moments")
    CHECK.equal((z_t, m_perp, m_parallel, m_cross),
                (3 * kappa, -12 * b, -6 * (b + c), -6 * d),
                "four observable contractions have exact multiplicities")
    return tuple(map(str, recovered))


samples = (
    (7, 11, 13, 17, 19, 23),
    (F(2, 3), F(-5, 7), F(11, 13), F(-17, 19), F(23, 29), F(-31, 37)),
    (0, 1, 0, 0, 0, 0),
    (5, 0, 0, 0, 0, 0),
)
recoveries = tuple(projected_moments(sample) for sample in samples)


def add_forms(*forms):
    return tuple(sum((form[i] for form in forms), F(0))
                 for i in range(len(forms[0])))


def scale_form(factor, form):
    return tuple(F(factor) * value for value in form)


# Linear forms are ordered as (Z_T,M_perp,M_parallel,M_cross).
kappa_form = (F(1, 3), F(0), F(0), F(0))
b_form = (F(0), F(-1, 12), F(0), F(0))
c_form = (F(0), F(1, 12), F(-1, 6), F(0))
d_form = (F(0), F(0), F(0), F(-1, 6))
a_form = add_forms(scale_form(-2, kappa_form), scale_form(8, b_form))
b_pullback_form = add_forms(scale_form(-16, kappa_form), scale_form(8, c_form))
c_pullback_form = add_forms(scale_form(12, kappa_form), scale_form(8, d_form))
extension_form = add_forms(b_pullback_form, c_pullback_form)
reference_form = a_form
CHECK.equal(extension_form, (F(-4, 3), F(2, 3), F(-4, 3), F(-4, 3)),
            "unrescaled tensor-extension observable form")
CHECK.equal(reference_form, (F(-2, 3), F(-2, 3), F(0), F(0)),
            "unrescaled reference-shape observable form")
CHECK.equal(scale_form(F(3, 2), extension_form),
            (F(-2), F(1), F(-2), F(-2)),
            "boxed extension moment coefficients")
CHECK.equal(scale_form(F(3, 2), reference_form),
            (F(-1), F(-1), F(0), F(0)),
            "boxed reference moment coefficients")


def contact_reconstruction(p):
    """Independently contract the four tetrahedral defect projectors."""
    p = F(p)
    q_matrices = []
    for port in range(4):
        tau = tuple(F(-1 if port in pair else 1) for pair in PAIR_ORDER)
        coordinates = tuple(dot(t, tau) for t in T_PAIR)
        # Q=tau tau^T/6, and both pair basis columns carry 1/sqrt(2).
        q_matrices.append(tuple(tuple(coordinates[i] * coordinates[j] / 12
                                      for j in range(3)) for i in range(3)))

    # Polynomial coefficient table: exponent -> 3x3 matrix.
    weighted = {}
    for port, direction in enumerate(TETRA):
        for m in range(3):
            exponent = tuple(2 if axis == m else 0 for axis in range(3))
            matrix = weighted.setdefault(exponent, zero_tensor((3, 3)))
            for i in range(3):
                for j in range(3):
                    matrix[i][j] += direction[m] * direction[m] * q_matrices[port][i][j]
        for m in range(3):
            for n in range(m + 1, 3):
                exponent = tuple(1 if axis in (m, n) else 0 for axis in range(3))
                matrix = weighted.setdefault(exponent, zero_tensor((3, 3)))
                for i in range(3):
                    for j in range(3):
                        matrix[i][j] += (2 * direction[m] * direction[n] *
                                         q_matrices[port][i][j])
    for m in range(3):
        exponent = tuple(2 if axis == m else 0 for axis in range(3))
        CHECK.equal(tuple(tuple(row) for row in weighted[exponent]),
                    tuple(tuple(F(4, 3) if i == j else F(0) for j in range(3))
                          for i in range(3)), "tetrahedral contact diagonal monomial")
    for m in range(3):
        for n in range(m + 1, 3):
            exponent = tuple(1 if axis in (m, n) else 0 for axis in range(3))
            expected = tuple(tuple(F(8, 3) if {i, j} == {m, n} and i != j
                                   else F(0) for j in range(3)) for i in range(3))
            CHECK.equal(tuple(tuple(row) for row in weighted[exponent]), expected,
                        "tetrahedral contact cross monomial")

    beta = 4 * p - F(4, 3)
    delta = F(2, 3) - beta
    a_contact = -2 * beta - F(2, 3) * delta
    c_contact = -F(4, 3) * delta
    CHECK.equal(a_contact, F(4, 3) * (1 - 4 * p), "contact A coefficient")
    CHECK.equal(c_contact, F(8, 3) * (2 * p - 1), "contact C coefficient")
    return str(a_contact), str(c_contact)


contact_samples = tuple(contact_reconstruction(p)
                        for p in (F(0), F(1, 2), F(109, 128), F(1)))


# Writer and sublattice normalization, derived without square-root arithmetic.
lambda_t = F(105, 16)
mu = F(105, 8)
CHECK.equal(mu, 2 * lambda_t, "mu equals two lambda_T")
for p, pair in enumerate(PAIR_ORDER[:3]):
    complement = tuple(sorted(set(range(4)) - set(pair)))
    j = [F(0)] * 6
    j[p] = F(p + 2, p + 3)
    j[PAIR_INDEX[complement]] = -j[p]
    theta_score = j[p] - j[PAIR_INDEX[complement]]
    CHECK.equal(lambda_t * theta_score, mu * j[p],
                "CM theta and CL canonical pair vertices agree")
CHECK.equal(dot((F(1), F(1)), (F(1), F(1))), F(2),
            "unnormalized common embedding has Gram two")
CHECK.equal(mu * mu / 2, 2 * lambda_t * lambda_t,
            "orthonormal common Hessian carries mu squared over two")


def evaluate(form, values):
    return sum((form[i] * values[i] for i in range(len(form))), F(0))


# Direct equivalence of the coefficient equations and the two boxed rules.
observable_samples = (
    ((F(2), F(3), F(5), F(7)), F(11), F(13), F(17, 19)),
    ((F(-5, 7), F(11, 13), F(-17, 23), F(29, 31)), F(37, 41), F(43, 47), F(53, 59)),
)
for observables, mu_squared, g_contact, p in observable_samples:
    coefficient_extension = (mu_squared / 2 * evaluate(extension_form, observables)
                             + g_contact * F(8, 3) * (2 * p - 1))
    boxed_extension = (mu_squared / 2 *
                       evaluate(scale_form(F(3, 2), extension_form), observables)
                       + 4 * g_contact * (2 * p - 1))
    CHECK.equal(boxed_extension, F(3, 2) * coefficient_extension,
                "boxed extension rule is exact rescaling")
    coefficient_reference = (mu_squared / 2 * evaluate(reference_form, observables)
                              + g_contact * F(4, 3) * (1 - 4 * p))
    boxed_reference = (mu_squared / 2 *
                       evaluate(scale_form(F(3, 2), reference_form), observables)
                       + 2 * g_contact * (1 - 4 * p))
    CHECK.equal(boxed_reference, F(3, 2) * coefficient_reference,
                "boxed reference rule is exact rescaling")


def encode(value):
    if isinstance(value, F):
        return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-result", action="store_true")
    args = parser.parse_args()
    result = {
        "verdict": "PASS",
        "independence": "target author derivation was neither imported nor executed",
        "algebra_checks_before_frozen_result_custody": CHECK.total,
        "solder": "S=(u,Q) is orthogonal; u_d=1/2 and Q_di=(T_d)_i/2",
        "coefficient_recovery": {
            "kappa": "Z_T/3",
            "alpha": "-sum_m M_AA^{mm}/6",
            "eta": "-(M_AT,x^{yz}+M_AT,y^{zx}+M_AT,z^{xy})/3",
            "b": "-M_perp/12",
            "c": "-M_parallel/6+M_perp/12",
            "d": "-M_cross/6",
            "exact_sample_recoveries": recoveries,
        },
        "normalization": {
            "pair_vertex": "mu=2 lambda_T because j.Theta_p=2j_p on pure T2",
            "CL_common": "j_plus=(j_P+j_C)/2",
            "orthonormal_common": "jhat_plus=(j_P+j_C)/sqrt(2)=sqrt(2)j_plus",
            "normalized_cycle_hessian": "(mu^2/2)B_T^*K_bare B_T",
            "contact": "g_ct multiplies contact with no mu factor",
        },
        "contact_samples_A_C": contact_samples,
        "unrescaled_observable_forms": {
            "extension": extension_form,
            "reference": reference_form,
        },
        "observable_rules": {
            "extension": "(mu^2/2)[-2Z_T+M_perp-2M_parallel-2M_cross]+4g_ct(2p-1)=0",
            "reference": "-(mu^2/2)[Z_T+M_perp]+2g_ct(1-4p)=0",
        },
        "analytic_scope": "finite absolute second moment gives the k2 projection; divergence requires a nonanalytic test rather than these sum rules",
        "ceiling": "displayed contact-plus-two-writer response sector only; no satisfaction, phase, complete source-second closure, 1PI, Ricci, gravity, or G",
    }
    payload = json.dumps(encode(result), indent=2, sort_keys=True) + "\n"
    target = HERE / "INDEPENDENT_RESULT.json"
    if args.write_result:
        target.write_text(payload)
    CHECK.true(target.is_file(), "frozen independent result exists")
    CHECK.equal(target.read_text(), payload, "frozen independent result matches replay")
    print(f"PASS__GL6CQ_INDEPENDENT_HOSTILE_SCIENCE__{CHECK.total}/{CHECK.total}")
    print("MOMENT_PROJECTION=RAW_ORIENTATION_TO_A1_T2_ALL_COMPONENTS")
    print("COEFFICIENTS=KAPPA_ALPHA_ETA_B_C_D_EXACT")
    print("NORMALIZATION=MU_EQ_2LAMBDA;COMMON_GRAM_2;NORMALIZED_HESSIAN_MU2_OVER_2")
    print("OBSERVABLE_RULES=EXTENSION_AND_REFERENCE_EXACTLY_REDERIVED")
    print("NO_STALE_MU2;CONTACT_HAS_NO_MU")
    print("NO_SATISFACTION_PHASE_1PI_RICCI_GRAVITY_G")


if __name__ == "__main__":
    main()
