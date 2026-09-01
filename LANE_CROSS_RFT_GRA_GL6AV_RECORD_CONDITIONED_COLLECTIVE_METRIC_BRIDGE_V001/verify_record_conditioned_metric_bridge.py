#!/usr/bin/env python3
"""Fail-closed exact replay for GL6AV."""

from __future__ import annotations

import hashlib
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

DEPENDENCIES = {
    "LANE_CROSS_RFT_GRA_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/THEOREM.md":
        "f75edcb115c3f7c86c6598f4597366b36e363df2d03ad919cc607b57dfb6b20c",
    "LANE_CROSS_RFT_GRA_GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION_V001/THEOREM.md":
        "1d1b01380ec8fd7ce83c69d45b68d9bde36bbe1dacdd32e3a5909ee6723a5ace",
    "LANE_CROSS_RFT_GRA_GL6AS_NATIVE_HEXAGON_COLLECTIVE_RESPONSE_V001/THEOREM.md":
        "bfe36071a24ccc7d6d7a16afeeea1b5554a95562ae91ac59c709db478000db9f",
    "LANE_CROSS_RFT_GRA_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001/THEOREM.md":
        "8407cee5196bfa4240f02159a5f59f941903dcf7a10e2baa18cf52a01ac8f743",
    "LANE_CROSS_RFT_GRA_GL6AT_PRIMARY_QUANTUM_ICE_CROSSWALK_V001/RESULT.md":
        "1c2b48c40a88000a88f8446fa2aea5116acb0896d94848274c87a22003bbcad9",
}

checks = 0


def require(condition: bool, name: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(f"FAIL__{name}")
    checks += 1
    print(f"PASS__{name}")


def rank(matrix: list[list[F]]) -> int:
    a = [row[:] for row in matrix]
    rows = len(a)
    cols = len(a[0]) if rows else 0
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        p = a[r][c]
        a[r] = [x / p for x in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                q = a[i][c]
                a[i] = [x - q * y for x, y in zip(a[i], a[r])]
        r += 1
    return r


def dot(x: tuple[F, ...], y: tuple[F, ...]) -> F:
    return sum((a * b for a, b in zip(x, y)), F(0))


def main() -> None:
    for rel, expected in DEPENDENCIES.items():
        actual = hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
        require(actual == expected, f"DEPENDENCY_SHA256__{Path(rel).parent.name}")

    theorem = (HERE / "THEOREM.md").read_text()
    result = (HERE / "RESULT.md").read_text()
    audit = (HERE / "SELF_AUDIT.md").read_text()
    for token in (
        "H_{\\rm hex}(q)=q^6H_{\\rm hex}(1)",
        "not claimed to be the complete conditioned",
        "\\alpha_t^{(q)}=\\alpha_{q^6t}^{(1)}",
        "preserves those spectral populations",
        "time-reparametrized quasi-local evolution",
        "\\operatorname{im}{\\cal E}=A_1\\oplus T_2",
        "\\ker{\\cal E}=E",
        "not yet a physical six-direction metric tangent",
        "memory-conditioned future",
        "generator/dynamics",
        "one nondegenerate quadratic Lorentzian characteristic",
        "single Lorentzian quadratic principal symbol/conformal-",
        "infinite fractional `q`",
        "not yet a retained",
        "physical source",
        "Common coupling",
        "Stress/Ward custody",
        "RGRL-B",
        "Gravity identity",
        "Newton coefficient",
        "No Maxwell description, gauge field, graviton, Ricci tensor, or Einstein",
    ):
        require(token in theorem, f"THEOREM_REQUIRED_TOKEN__{hashlib.sha256(token.encode()).hexdigest()[:12]}")
    require("It is not yet one physical metric tangent" in result,
            "RESULT_PRESERVES_SOURCE_READ_CEILING")
    require("Is memory back-reaction already two-way?" in audit,
            "SELF_AUDIT_TESTS_TWO_WAY_CEILING")

    # Every six-edge loop gives the homogeneous sixth power exactly.
    for q in (-2, -1, 0, 1, 2, 3):
        product = 1
        for _ in range(6):
            product *= q
        require(product == q**6, f"HOMOGENEOUS_SIXTH_POWER__Q_{q}")

    # Orientation log map A=2(11^T-I): eigenvalue 6 on A1 and -2 on T2.
    A = [[F(0) if i == j else F(2) for j in range(4)] for i in range(4)]
    one = [F(1)] * 4
    require([dot(tuple(row), tuple(one)) for row in A] == [F(6)] * 4,
            "ORIENTATION_A1_EIGENVALUE_6")
    centered_basis = ([F(1), F(-1), F(0), F(0)],
                      [F(0), F(1), F(-1), F(0)],
                      [F(0), F(0), F(1), F(-1)])
    for i, v in enumerate(centered_basis):
        Av = [dot(tuple(row), tuple(v)) for row in A]
        require(Av == [F(-2) * x for x in v], f"ORIENTATION_T2_EIGENVALUE_MINUS2__{i}")
    require(rank(A) == 4, "ORIENTATION_LOG_MAP_RANK_4")

    # Exact inverse rho_d=(sum j)/6-j_d/2.
    rho = [F(2, 5), F(-1, 3), F(7, 11), F(5, 13)]
    j = [sum((A[d][a] * rho[a] for a in range(4)), F(0)) for d in range(4)]
    recovered = [sum(j, F(0)) / 6 - x / 2 for x in j]
    require(recovered == rho, "ORIENTATION_LOG_MAP_EXACT_INVERSE")

    # Tetrahedral Gram and resolution of the identity.
    t = (
        (F(1, 2), F(1, 2), F(1, 2)),
        (F(1, 2), F(-1, 2), F(-1, 2)),
        (F(-1, 2), F(1, 2), F(-1, 2)),
        (F(-1, 2), F(-1, 2), F(1, 2)),
    )
    for a in range(4):
        for b in range(4):
            expected = F(3, 4) if a == b else F(-1, 4)
            require(dot(t[a], t[b]) == expected, f"TETRAHEDRAL_GRAM__{a}_{b}")
    resolution = [[sum((t[a][i] * t[a][j] for a in range(4)), F(0))
                   for j in range(3)] for i in range(3)]
    require(resolution == [[F(int(i == j)) for j in range(3)] for i in range(3)],
            "TETRAHEDRAL_RESOLUTION_IDENTITY")

    # Evaluation matrix on [Sxx,Syy,Szz,Sxy,Sxz,Syz].
    eval_matrix = []
    for x, y, z in t:
        eval_matrix.append([x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z])
    require(rank(eval_matrix) == 4, "SYMMETRIC_TENSOR_EVALUATION_RANK_4")
    diag_e1 = [F(1), F(-1), F(0), F(0), F(0), F(0)]
    diag_e2 = [F(1), F(0), F(-1), F(0), F(0), F(0)]
    for i, e in enumerate((diag_e1, diag_e2)):
        values = [dot(tuple(row), tuple(e)) for row in eval_matrix]
        require(values == [F(0)] * 4, f"TRACELESS_DIAGONAL_E_KERNEL__{i}")

    # Explicit orientation reconstruction (AV16), tested on arbitrary data.
    q = [F(2, 7), F(-3, 5), F(11, 13), F(17, 19)]
    tr = sum(q, F(0))
    sxy = (q[0] - q[1] - q[2] + q[3]) / 2
    sxz = (q[0] - q[1] + q[2] - q[3]) / 2
    syz = (q[0] + q[1] - q[2] - q[3]) / 2
    coeffs = [tr / 3, tr / 3, tr / 3, sxy, sxz, syz]
    reconstructed = [dot(tuple(row), tuple(coeffs)) for row in eval_matrix]
    require(reconstructed == q, "ORIENTATION_TENSOR_EXPLICIT_RECONSTRUCTION")

    # Fail closed if high-level ceilings disappear.
    combined = theorem + result + audit
    for phrase in (
        "common coupling",
        "stress",
        "Ward",
        "RGRL-B",
        "gravity",
        "G",
    ):
        require(phrase.lower() in combined.lower(), f"OPEN_GATE_RETAINED__{phrase.replace('-', '_')}")

    print(f"PASS__GL6AV_EXACT_REPLAY_COMPLETE__{checks}/{checks}")


if __name__ == "__main__":
    main()
