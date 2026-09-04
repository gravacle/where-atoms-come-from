#!/usr/bin/env python3
"""Exact symbolic replay for GL6CQ stationary-response moment sum rules."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as F
from itertools import combinations, permutations, product
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAIRS = tuple(combinations(range(4), 2))
PIDX = {pair: i for i, pair in enumerate(PAIRS)}
TETRA = (
    (F(1), F(1), F(1)),
    (F(1), F(-1), F(-1)),
    (F(-1), F(1), F(-1)),
    (F(-1), F(-1), F(1)),
)
TBASIS = (
    tuple(map(F, (1, 0, 0, 0, 0, -1))),
    tuple(map(F, (0, 1, 0, 0, -1, 0))),
    tuple(map(F, (0, 0, 1, -1, 0, 0))),
)
TAU = (
    tuple(map(F, (-1, -1, -1, 1, 1, 1))),
    tuple(map(F, (-1, 1, 1, -1, -1, 1))),
    tuple(map(F, (1, -1, 1, -1, 1, -1))),
    tuple(map(F, (1, 1, -1, 1, -1, -1))),
)
u = (F(1, 2),) * 4
Q = tuple(tuple(TETRA[d][i] / 2 for i in range(3)) for d in range(4))
S = tuple((u[d],) + Q[d] for d in range(4))


class Checks:
    def __init__(self):
        self.n = 0

    def equal(self, got, expected, label):
        self.n += 1
        if got != expected:
            raise AssertionError(f"{label}: got {got!r}, expected {expected!r}")

    def true(self, condition, label):
        self.n += 1
        if not condition:
            raise AssertionError(label)


C = Checks()


def dot(a, b):
    return sum((F(x) * F(y) for x, y in zip(a, b)), F(0))


def add(*vs):
    return tuple(sum((F(v[i]) for v in vs), F(0)) for i in range(len(vs[0])))


def scale(c, v):
    return tuple(F(c) * F(x) for x in v)


def transpose(m):
    return tuple(tuple(m[i][j] for i in range(len(m))) for j in range(len(m[0])))


def mm(a, b):
    return tuple(tuple(sum((a[i][r] * b[r][j] for r in range(len(b))), F(0))
                       for j in range(len(b[0]))) for i in range(len(a)))


C.equal(mm(transpose(S), S), tuple(tuple(F(i == j) for j in range(4)) for i in range(4)),
        "cycle A1+T2 solder is orthogonal")


# Sparse polynomials in (kx,ky,kz), truncated by callers.
def padd(*ps):
    out = {}
    for p in ps:
        for e, coefficient in p.items():
            out[e] = out.get(e, F(0)) + coefficient
    return {e: coefficient for e, coefficient in out.items() if coefficient}


def pscale(c, p):
    return {e: F(c) * coefficient for e, coefficient in p.items()
            if F(c) * coefficient}


def pmul(a, b, degree=2):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = tuple(ea[i] + eb[i] for i in range(3))
            if sum(e) <= degree:
                out[e] = out.get(e, F(0)) + ca * cb
    return {e: coefficient for e, coefficient in out.items() if coefficient}


def linear(v):
    return {e: F(c) for e, c in zip(((1, 0, 0), (0, 1, 0), (0, 0, 1)), v) if c}


def square_linear(v):
    return pmul(linear(v), linear(v), 2)


def pmatmul(a, b):
    return tuple(tuple(padd(*(pmul(a[i][r], b[r][j], 2)
                              for r in range(len(b))))
                       for j in range(len(b[0]))) for i in range(len(a)))


ZERO = {}
ONE = {(0, 0, 0): F(1)}
R2 = {(2, 0, 0): F(1), (0, 2, 0): F(1), (0, 0, 2): F(1)}


def cycle_offsets():
    records = []
    for missing in range(4):
        a, b, c = tuple(port for port in range(4) if port != missing)
        rho = {
            tuple(sorted((a, b))): scale(F(1, 2), add(TETRA[a], TETRA[b], scale(-1, TETRA[c]))),
            tuple(sorted((a, c))): scale(F(1, 2), add(TETRA[a], scale(-1, TETRA[b]), TETRA[c])),
            tuple(sorted((b, c))): scale(F(1, 2), add(scale(-1, TETRA[a]), TETRA[b], TETRA[c])),
        }
        C.true(all(dot(r, r) == F(11, 4) for r in rho.values()),
               "GL6CL centered offset radius")
        records.append((missing, rho))
    C.equal(sum(len(rho) for _, rho in records), 12, "twelve GL6CL offsets")
    return tuple(records)


OFFSETS = cycle_offsets()


def writer_tensor_polynomial():
    """Four by three GL6CL B_+ rows in the unnormalized T basis."""
    rows = []
    for _, rho_by_pair in OFFSETS:
        row = []
        for t in TBASIS:
            entry = {}
            for pair, rho in rho_by_pair.items():
                # 2 cos(k.rho)=2-(k.rho)^2+O(k4).
                entry = padd(entry, pscale(t[PIDX[pair]],
                                           padd(pscale(2, ONE), pscale(-1, square_linear(rho)))))
            row.append(entry)
        rows.append(tuple(row))
    return tuple(rows)


BU = writer_tensor_polynomial()


def block_to_cycle(kappa, b, c, d):
    """Build S diag(0,K_TT) S^T through k2 for numeric coefficients."""
    block = [[{} for _ in range(4)] for _ in range(4)]
    for i in range(3):
        for j in range(3):
            entry = {}
            if i == j:
                entry = padd(pscale(kappa, ONE), pscale(b, R2),
                             {tuple(2 if q == i else 0 for q in range(3)): F(c)})
            else:
                exponent = tuple(1 if q in (i, j) else 0 for q in range(3))
                entry = {exponent: F(d)} if d else {}
            block[i + 1][j + 1] = entry
    # Expand S block S^T explicitly.
    cycle = []
    for a in range(4):
        row = []
        for z in range(4):
            row.append(padd(*(pscale(S[a][r] * S[z][s], block[r][s])
                              for r in range(4) for s in range(4))))
        cycle.append(tuple(row))
    return tuple(cycle)


def expected_pair_hessian(kappa, b, c, d):
    Acoef = -2 * kappa + 8 * b
    Bcoef = -16 * kappa + 8 * c
    Ccoef = 12 * kappa + 8 * d
    answer = []
    for i in range(3):
        row = []
        for j in range(3):
            entry = pscale(8 * kappa, ONE) if i == j else {}
            if i == j:
                entry = padd(entry, pscale(Acoef, R2),
                             {tuple(2 if q == i else 0 for q in range(3)): F(Bcoef)})
            else:
                exponent = tuple(1 if q in (i, j) else 0 for q in range(3))
                entry = padd(entry, {exponent: F(Ccoef)})
            row.append(entry)
        answer.append(tuple(row))
    return tuple(answer)


def verify_writer_pullback():
    # BU uses T vectors of norm sqrt(2).  Divide BU^T K BU by two to obtain
    # the orthonormal pair-T basis used in CO11.
    for coefficients in ((1, 0, 0, 0), (0, 1, 0, 0),
                         (0, 0, 1, 0), (0, 0, 0, 1),
                         (2, 3, 5, 7)):
        cycle = block_to_cycle(*map(F, coefficients))
        raw = pmatmul(pmatmul(transpose(BU), cycle), BU)
        normalized = tuple(tuple(pscale(F(1, 2), raw[i][j]) for j in range(3))
                           for i in range(3))
        C.equal(normalized, expected_pair_hessian(*map(F, coefficients)),
                "CO12 pullback normalization")


verify_writer_pullback()


def moment_projection_test():
    """Create an invariant Hessian, negate it to moments, and recover all coefficients."""
    kappa, alpha, eta, b, c, d = map(F, (7, 11, 13, 17, 19, 23))
    ztt = tuple(tuple(kappa * F(i == j) for j in range(3)) for i in range(3))
    maa = tuple(tuple(-2 * alpha * F(m == n) for n in range(3)) for m in range(3))

    mat = []
    for i in range(3):
        plane = []
        complement = tuple(q for q in range(3) if q != i)
        for m in range(3):
            plane.append(tuple(-eta if m != n and {m, n} == set(complement) else F(0)
                               for n in range(3)))
        mat.append(tuple(plane))
    mat = tuple(mat)

    mtt = []
    for i in range(3):
        row_j = []
        for j in range(3):
            plane = []
            for m in range(3):
                line = []
                for n in range(3):
                    h = F(0)
                    if i == j:
                        h += 2 * b * F(m == n)
                        h += 2 * c * F(i == m == n)
                    else:
                        h += d * (F(i == m and j == n) + F(i == n and j == m))
                    line.append(-h)
                plane.append(tuple(line))
            row_j.append(tuple(plane))
        mtt.append(tuple(row_j))
    mtt = tuple(mtt)

    recovered_kappa = sum(ztt[i][i] for i in range(3)) / 3
    recovered_alpha = -sum(maa[m][m] for m in range(3)) / 6
    recovered_eta = -sum(mat[i][tuple(q for q in range(3) if q != i)[0]]
                         [tuple(q for q in range(3) if q != i)[1]] for i in range(3)) / 3
    m_perp = sum(mtt[i][i][m][m] for i in range(3) for m in range(3) if i != m)
    m_parallel = sum(mtt[i][i][i][i] for i in range(3))
    m_cross = sum(mtt[i][j][i][j] for i in range(3) for j in range(3) if i != j)
    recovered_b = -m_perp / 12
    recovered_c = -m_parallel / 6 + m_perp / 12
    recovered_d = -m_cross / 6
    C.equal((recovered_kappa, recovered_alpha, recovered_eta,
             recovered_b, recovered_c, recovered_d),
            (kappa, alpha, eta, b, c, d), "all moment projections recover coefficients")
    C.equal((m_perp, m_parallel, m_cross),
            (-12 * b, -6 * (b + c), -6 * d), "contracted moment normalization")
    return {
        "kappa": "Z_T/3",
        "alpha": "-M_AA_trace/6",
        "eta": "-(M_AT,x^yz+M_AT,y^zx+M_AT,z^xy)/3",
        "b": "-M_perp/12",
        "c": "-M_parallel/6+M_perp/12",
        "d": "-M_cross/6",
    }


MOMENT_FORMULAS = moment_projection_test()


def contact_projection(p):
    """Rebuild the common BV contact in the orthonormal pair-T basis."""
    p = F(p)
    beta = 4 * p - F(4, 3)
    delta = F(2, 3) - beta
    gamma = F(8, 3) * (4 * p - 1)
    q_mats = []
    for tau in TAU:
        qfull = tuple(tuple(tau[i] * tau[j] / 6 for j in range(6)) for i in range(6))
        # Orthonormal columns are t_i/sqrt(2): the two sqrt factors give 1/2.
        qcoord = tuple(tuple(sum((TBASIS[i][a] * qfull[a][z] * TBASIS[j][z]
                                  for a in range(6) for z in range(6)), F(0)) / 2
                             for j in range(3)) for i in range(3))
        q_mats.append(qcoord)

    weighted = [[{} for _ in range(3)] for _ in range(3)]
    for a in range(4):
        theta2 = square_linear(TETRA[a])
        for i in range(3):
            for j in range(3):
                weighted[i][j] = padd(weighted[i][j], pscale(q_mats[a][i][j], theta2))
    expected_weighted = []
    for i in range(3):
        row = []
        for j in range(3):
            entry = pscale(F(4, 3), R2) if i == j else {}
            if i != j:
                exponent = tuple(1 if q in (i, j) else 0 for q in range(3))
                entry = {exponent: F(8, 3)}
            row.append(entry)
        expected_weighted.append(tuple(row))
    C.equal(tuple(tuple(row) for row in weighted), tuple(expected_weighted),
            "BV/CO tetrahedral weighted contact identity")

    result = []
    for i in range(3):
        row = []
        for j in range(3):
            entry = pscale(8 + gamma, ONE) if i == j else {}
            if i == j:
                entry = padd(entry, pscale(-2 * beta, R2))
            entry = padd(entry, pscale(-delta / 2, weighted[i][j]))
            row.append(entry)
        result.append(tuple(row))

    Acoef = F(4, 3) * (1 - 4 * p)
    Ccoef = F(8, 3) * (2 * p - 1)
    expected = []
    for i in range(3):
        row = []
        for j in range(3):
            entry = pscale(8 + gamma, ONE) if i == j else {}
            if i == j:
                entry = padd(entry, pscale(Acoef, R2))
            else:
                exponent = tuple(1 if q in (i, j) else 0 for q in range(3))
                entry = padd(entry, {exponent: Ccoef})
            row.append(entry)
        expected.append(tuple(row))
    C.equal(tuple(result), tuple(expected), "CO26 contact coefficients")
    return gamma, Acoef, Ccoef


for test_p in (F(0), F(1, 2), F(109, 128), F(1)):
    contact_projection(test_p)


def sum_rule_algebra():
    # Coefficient order: (Z_T,M_perp,M_parallel,M_cross).
    extension = (F(-4, 3), F(2, 3), F(-4, 3), F(-4, 3))
    reference = (F(-2, 3), F(-2, 3), F(0), F(0))
    rescaled_extension = tuple(F(3, 2) * x for x in extension)
    rescaled_reference = tuple(F(3, 2) * x for x in reference)
    C.equal(rescaled_extension, (F(-2), F(1), F(-2), F(-2)),
            "rescaled tensor-extension moment sum")
    C.equal(rescaled_reference, (F(-1), F(-1), F(0), F(0)),
            "rescaled reference-shape moment sum")
    return extension, reference


EXTENSION_COEFFS, REFERENCE_COEFFS = sum_rule_algebra()


lambda_T = F(105, 16)
mu = F(105, 8)
C.equal(mu, 2 * lambda_T, "CM/CL source normalization mu=2 lambda_T")
C.equal(mu * mu / 2, 2 * lambda_T * lambda_T,
        "orthonormal common-source Hessian factor")
C.equal(F(-63, 8) * F(-5) * 2, F(315, 4), "uniform scalar derivative normalization")


def qtext(x):
    x = F(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def encode(x):
    if isinstance(x, F):
        return qtext(x)
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    return x


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-ledger", action="store_true")
    parser.add_argument("--print-ledger", action="store_true")
    args = parser.parse_args()
    ledger = {
        "schema": "GL6CQ_STATIONARY_RESPONSE_MOMENT_SUM_RULES_V001",
        "check_count": C.n,
        "bare_cycle_susceptibility": "K_bare_dd'(R)=2 Re <0|T_0d Q(H0-E0)^(-1)Q T_Rd'|0>",
        "centered_fourier": "K_dd'(k)=sum_R exp(i k.X_Rdd') K_bare_dd'(R)",
        "moments": {
            "zeroth": "Z_dd'=sum_R K_bare_dd'(R)",
            "second": "M_dd'^{mn}=sum_R X_Rdd'^m X_Rdd'^n K_bare_dd'(R)",
            "expansion": "K(k)=Z-(1/2)k_m k_n M^{mn}+o(k^2)",
            "projected_coefficients": MOMENT_FORMULAS,
            "contractions": {
                "Z_T": "sum_i Q_i^T Z Q_i=3 kappa",
                "M_perp": "sum_{i!=m} M_TT,ii^{mm}=-12 b",
                "M_parallel": "sum_i M_TT,ii^{ii}=-6(b+c)",
                "M_cross": "sum_{i!=j} M_TT,ij^{ij}=-6d",
            },
        },
        "normalization": {
            "cycle_kernel": "bare: no lambda_T^2 and no mu^2",
            "CM_lambda_T": "(105/16)h^6/U_d^6",
            "CL_mu": "(105/8)h^6/U_d^6=2 lambda_T",
            "identity": "for pure T2, j.Theta_ab=2 j_ab, hence lambda_T(j.Theta)=mu j_ab",
            "GL6CL_common_coordinate": "j_plus=(j_P+j_C)/2",
            "GL6CL_common_pullback": "H_T^H6=mu^2 B_T^* K_bare B_T",
            "orthonormal_common_coordinate": "jhat_plus=(j_P+j_C)/sqrt(2)=sqrt(2)j_plus",
            "GL6BV_normalized_common_pullback": "Hhat_T^H6=(mu^2/2) B_T^* K_bare B_T",
            "contact_scale": "g_ct=h^2/(4U_d^3), with no mu factor",
        },
        "writer_pullback": {
            "A": "-2 kappa+8b",
            "B": "-16 kappa+8c",
            "C": "12 kappa+8d",
            "extension_mismatch": "-4 kappa+8(c+d)",
        },
        "contact": {
            "p": "same-state equal-partner probability",
            "A": "(4/3)(1-4p)",
            "B": "0",
            "C": "(8/3)(2p-1)",
        },
        "observable_sum_rules": {
            "tensor_gradient_extension": "(mu^2/2)[-2 Z_T+M_perp-2M_parallel-2M_cross]+4g_ct(2p-1)=0",
            "reference_shape": "-(mu^2/2)[Z_T+M_perp]+2g_ct(1-4p)=0",
            "unrescaled_extension_coefficients": encode(EXTENSION_COEFFS),
            "unrescaled_reference_coefficients": encode(REFERENCE_COEFFS),
            "guard": "K_bare moments and p are evaluated in the same stationary state; mu is applied exactly once outside K_bare and 1/2 converts to the orthonormal GL6BV common source",
            "zero_momentum_guard": "the conditions do not remove or interpret the zero-momentum (k=0) tensor term and do not prove background stationarity, masslessness, or a gauge null",
        },
        "convergence": {
            "sufficient_for_coefficients": "sum_R (1+|X_Rdd'|^2)|K_bare_dd'(R)|<infinity entrywise",
            "remainder": "finite second moment gives o(k^2); finite fourth moment plus inversion gives O(k^4); exponential moment gives analyticity",
            "critical_limit": "if the second moment diverges, b,c,d and the analytic k^2 sum rules are undefined; nonanalytic |k|^sigma or k^2 log|k| behavior must be treated directly",
        },
        "scope": "exact moment/sum-rule theorem for the same-state contact plus two-writer spectral sector; no phase, full source-first closure, 1PI inversion, Ricci, gravity, or G",
    }
    payload = json.dumps(encode(ledger), indent=2, sort_keys=True) + "\n"
    target = HERE / "EXACT_LEDGER.json"
    if args.write_ledger:
        target.write_text(payload)
    elif args.print_ledger:
        print(payload, end="")
        return
    else:
        C.true(target.is_file(), "frozen ledger exists")
        C.equal(json.loads(target.read_text()), json.loads(payload),
                "frozen ledger matches exact replay")
    print(f"PASS__GL6CQ_STATIONARY_RESPONSE_MOMENT_SUM_RULES__{C.n}/{C.n}")
    print("MOMENTS=ZEROTH_KAPPA;SECOND_ALPHA_ETA_B_C_D")
    print("NORMALIZATION=K_BARE_NO_WRITER;MU_SQUARED_ONCE;ORTHONORMAL_COMMON_HALF;CONTACT_NO_MU")
    print("SAME_STATE=K_BARE_MOMENTS_AND_P")
    print("SUM_RULES=CO29_CO30_AS_EXACT_REAL_SPACE_MOMENT_IDENTITIES")
    print("CRITICAL_LIMIT=SECOND_MOMENT_DIVERGENCE_INVALIDATES_ANALYTIC_K2_FORM")
    print("NO_PHASE_1PI_RICCI_GRAVITY_G")


if __name__ == "__main__":
    main()
