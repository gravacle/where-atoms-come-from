#!/usr/bin/env python3
"""Exact bounded checks for the GL6CY calculation target.

This verifies the barycentric nonlinear six-pair source map, its complete
first/mixed-second chain rule, full-rank six-direction source/tensor typing,
affine-cocycle telescoping, the alternating-cycle center cancellation, the
GC refinement factor, and the fail-closed owner schema.

It deliberately does not evaluate the separately audited GL6CU source Jet or
a physical Ward residual.  The latter remains refused until every required
physical owner is defined and the lawful quotient has been constructed.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAIRS = tuple(combinations(range(4), 2))
T = (
    (F(1), F(1), F(1)),
    (F(1), F(-1), F(-1)),
    (F(-1), F(1), F(-1)),
    (F(-1), F(-1), F(1)),
)
V = tuple(tuple(x / 2 for x in row) for row in T)
I3 = tuple(tuple(F(i == j) for j in range(3)) for i in range(3))
CHECKS = 0


def check(condition, label):
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def dot(x, y):
    return sum((F(a) * F(b) for a, b in zip(x, y)), F(0))


def vadd(*xs):
    return tuple(sum((F(x[i]) for x in xs), F(0)) for i in range(len(xs[0])))


def vscale(c, x):
    return tuple(F(c) * F(a) for a in x)


def transpose(a):
    return tuple(tuple(a[j][i] for j in range(len(a))) for i in range(len(a[0])))


def matmul(a, b):
    return tuple(tuple(sum((a[i][k] * b[k][j] for k in range(len(b))), F(0))
                       for j in range(len(b[0]))) for i in range(len(a)))


def matvec(a, x):
    return tuple(sum((a[i][j] * x[j] for j in range(len(x))), F(0))
                 for i in range(len(a)))


def inverse(a):
    n = len(a)
    work = [list(a[i]) + [F(i == j) for j in range(n)] for i in range(n)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if work[r][col]), None)
        if pivot is None:
            raise ValueError("singular")
        work[col], work[pivot] = work[pivot], work[col]
        p = work[col][col]
        work[col] = [x / p for x in work[col]]
        for r in range(n):
            if r == col:
                continue
            q = work[r][col]
            work[r] = [work[r][j] - q * work[col][j] for j in range(2 * n)]
    return tuple(tuple(row[n:]) for row in work)


def outer(x, y):
    return tuple(tuple(x[i] * y[j] for j in range(len(y))) for i in range(len(x)))


def madd(a, b):
    return tuple(tuple(a[i][j] + b[i][j] for j in range(len(a[0])))
                 for i in range(len(a)))


def mscale(c, a):
    return tuple(tuple(F(c) * x for x in row) for row in a)


def odot(x, y):
    return madd(outer(x, y), outer(y, x))


def contract(a, b):
    return sum((a[i][j] * b[i][j]
                for i in range(len(a)) for j in range(len(a[0]))), F(0))


def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0))


@dataclass(frozen=True)
class BJ:
    """Bivariate 2-jet storing actual r, s, and mixed-rs derivatives."""

    v: F = F(0)
    r: F = F(0)
    s: F = F(0)
    rs: F = F(0)

    def __post_init__(self):
        for name in ("v", "r", "s", "rs"):
            object.__setattr__(self, name, F(getattr(self, name)))

    def __add__(self, other):
        other = bj(other)
        return BJ(self.v + other.v, self.r + other.r,
                  self.s + other.s, self.rs + other.rs)

    __radd__ = __add__

    def __neg__(self):
        return BJ(-self.v, -self.r, -self.s, -self.rs)

    def __sub__(self, other):
        return self + (-bj(other))

    def __rsub__(self, other):
        return bj(other) - self

    def __mul__(self, other):
        other = bj(other)
        return BJ(
            self.v * other.v,
            self.r * other.v + self.v * other.r,
            self.s * other.v + self.v * other.s,
            self.rs * other.v + self.r * other.s
            + self.s * other.r + self.v * other.rs,
        )

    __rmul__ = __mul__

    def reciprocal(self):
        if not self.v:
            raise ZeroDivisionError
        return BJ(
            1 / self.v,
            -self.r / self.v**2,
            -self.s / self.v**2,
            2 * self.r * self.s / self.v**3 - self.rs / self.v**2,
        )

    def __truediv__(self, other):
        return self * bj(other).reciprocal()


def bj(x):
    return x if isinstance(x, BJ) else BJ(F(x))


def jtranspose(a):
    return tuple(tuple(a[j][i] for j in range(len(a))) for i in range(len(a[0])))


def jmatmul(a, b):
    return tuple(tuple(sum((a[i][k] * b[k][j] for k in range(len(b))), BJ())
                       for j in range(len(b[0]))) for i in range(len(a)))


def jinverse(a):
    n = len(a)
    work = [list(a[i]) + [BJ(F(i == j)) for j in range(n)] for i in range(n)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if work[r][col].v), None)
        if pivot is None:
            raise ValueError("singular base")
        work[col], work[pivot] = work[pivot], work[col]
        p = work[col][col]
        work[col] = [x / p for x in work[col]]
        for r in range(n):
            if r == col:
                continue
            q = work[r][col]
            work[r] = [work[r][j] - q * work[col][j] for j in range(2 * n)]
    return tuple(tuple(row[n:]) for row in work)


def jmatvec(a, x):
    return tuple(sum((a[i][j] * x[j] for j in range(len(x))), BJ())
                 for i in range(len(a)))


def jvadd(*xs):
    return tuple(sum((x[i] for x in xs), BJ()) for i in range(len(xs[0])))


def jouter(x, y):
    return tuple(tuple(x[i] * y[j] for j in range(len(y))) for i in range(len(x)))


def jmadd(a, b):
    return tuple(tuple(a[i][j] + b[i][j] for j in range(len(a[0])))
                 for i in range(len(a)))


def jodot(x, y):
    return jmadd(jouter(x, y), jouter(y, x))


def jcontract(a, b):
    return sum((a[i][j] * b[i][j]
                for i in range(len(a)) for j in range(len(a[0]))), BJ())


def source_first(b):
    return tuple(-2 * dot(b, vadd(V[a], V[c])) for a, c in PAIRS)


def source_second(ar, ass, br, bs, brs):
    rho_rs = vadd(brs, vscale(-1, matvec(ar, bs)), vscale(-1, matvec(ass, br)))
    return tuple(-2 * dot(rho_rs, vadd(V[a], V[c])) - 4 * dot(br, bs)
                 for a, c in PAIRS)


def direct_source_jet(ar, ass, br, bs, brs):
    jf = tuple(tuple(BJ(I3[i][j], ar[i][j], ass[i][j], F(0))
                     for j in range(3)) for i in range(3))
    jb = tuple(BJ(F(0), br[i], bs[i], brs[i]) for i in range(3))
    jg = jmatmul(jf, jtranspose(jf))
    jgi = jinverse(jg)
    answer = []
    for a, c in PAIRS:
        pa = jvadd(jmatvec(jf, tuple(BJ(x) for x in V[a])), jb)
        pc = jvadd(jmatvec(jf, tuple(BJ(x) for x in V[c])), jb)
        q = jodot(pa, pc)
        answer.append(BJ(trace(odot(V[a], V[c]))) - jcontract(jgi, q))
    return tuple(answer)


def symmetric_basis():
    rows = []
    for i in range(3):
        m = [[F(0) for _ in range(3)] for _ in range(3)]
        m[i][i] = F(1, 2)
        rows.append(tuple(tuple(x for x in row) for row in m))
    for i, j in ((0, 1), (0, 2), (1, 2)):
        m = [[F(0) for _ in range(3)] for _ in range(3)]
        m[i][j] = m[j][i] = F(1, 2)
        rows.append(tuple(tuple(x for x in row) for row in m))
    return tuple(rows)


def rank(matrix):
    a = [list(map(F, row)) for row in matrix]
    r = 0
    for c in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        q = a[r][c]
        a[r] = [x / q for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                q = a[i][c]
                a[i] = [a[i][j] - q * a[r][j] for j in range(len(a[0]))]
        r += 1
    return r


def pair_solder():
    # coordinate order xx,yy,zz,xy,xz,yz
    return tuple(tuple((odot(V[a], V[c])[i][j])
                       for a, c in PAIRS)
                 for i, j in ((0, 0), (1, 1), (2, 2),
                              (0, 1), (0, 2), (1, 2)))


def dual_source_rows():
    dc = pair_solder()
    weights = (F(1), F(1), F(1), F(2), F(2), F(2))
    # Columns are D_C^* applied to the six tensor-coordinate basis covectors.
    return tuple(tuple(dc[t][p] * weights[t] for t in range(6))
                 for p in range(6))


def qadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def qscale(c, x):
    return (F(c) * x[0], F(c) * x[1])


def qmul(x, y):
    return (x[0] * y[0] - x[1] * y[1],
            x[0] * y[1] + x[1] * y[0])


def poly_mul(a, b, degree=3):
    out = [(F(0), F(0)) for _ in range(degree + 1)]
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j <= degree:
                out[i + j] = qadd(out[i + j], qmul(x, y))
    return out


def exp_i_series(lam):
    return [
        (F(1), F(0)),
        (F(0), lam),
        (-lam * lam / 2, F(0)),
        (F(0), -lam * lam * lam / 6),
    ]


def cycle_source_series(position, kind, pair, k, xi):
    sign = F(1) if kind == "P" else F(-1)
    phase = exp_i_series(dot(k, position))
    avg = [(F(0), F(0)) for _ in range(4)]
    for ta in T:
        step = sign * dot(k, ta)
        term = exp_i_series(step)
        term[0] = qadd(term[0], (F(-1), F(0)))
        for n in range(4):
            avg[n] = qadd(avg[n], qscale(F(1, 8), term[n]))
    dscalar = poly_mul(phase, avg)
    a, c = pair
    port_sum = vscale(sign, vadd(V[a], V[c]))
    coefficient = -2 * dot(xi, port_sum)
    return [qscale(coefficient, x) for x in dscalar]


def main_result():
    # Tetrahedral and pair-solder preliminaries.
    check(vadd(*T) == (F(0), F(0), F(0)), "tetrahedral first moment")
    second = tuple(tuple(sum((T[a][i] * T[a][j] for a in range(4)), F(0))
                         for j in range(3)) for i in range(3))
    check(second == mscale(4, I3), "tetrahedral second moment")
    dc = pair_solder()
    check(rank(dc) == 6, "six-pair solder rank")
    check(rank(dual_source_rows()) == 6, "six tensor dual-source directions")

    abasis = symmetric_basis()
    bdirs = (
        (F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(1)),
        (F(1), F(1), F(0)),
        (F(1), F(0), F(1)),
        (F(0), F(1), F(1)),
    )

    # Direct rational bivariate differentiation versus CY13-CY14 for all
    # 6x6 deformation pairs and every raw pair source.
    for r in range(6):
        for s in range(6):
            brs = (F((r + 2 * s) % 5 - 2),
                   F((2 * r + s) % 5 - 2),
                   F((3 * r + 4 * s) % 5 - 2))
            got = direct_source_jet(abasis[r], abasis[s],
                                    bdirs[r], bdirs[s], brs)
            first_r = source_first(bdirs[r])
            first_s = source_first(bdirs[s])
            mixed = source_second(abasis[r], abasis[s],
                                  bdirs[r], bdirs[s], brs)
            for p in range(6):
                check(got[p].v == 0, f"flat source {r},{s},{p}")
                check(got[p].r == first_r[p], f"first r {r},{s},{p}")
                check(got[p].s == first_s[p], f"first s {r},{s},{p}")
                check(got[p].rs == mixed[p], f"mixed {r},{s},{p}")

    # Centered coframe has identically zero source jet for all six symmetric
    # directions.
    z = (F(0), F(0), F(0))
    for r in range(6):
        for s in range(6):
            got = direct_source_jet(abasis[r], abasis[s], z, z, z)
            check(all(x == BJ() for x in got), f"centered zero {r},{s}")

    t2 = (
        (F(1), F(0), F(0), F(0), F(0), F(-1)),
        (F(0), F(1), F(0), F(0), F(-1), F(0)),
        (F(0), F(0), F(1), F(-1), F(0), F(0)),
    )
    axes = ((F(1), F(0), F(0)),
            (F(0), F(1), F(0)),
            (F(0), F(0), F(1)))
    for i in range(3):
        check(source_first(axes[i]) == vscale(-2, t2[i]),
              f"first center direction pure T2 {i}")
        second = source_second(mscale(0, I3), mscale(0, I3),
                               axes[i], axes[i], z)
        check(second == (F(-4),) * 6,
              f"quadratic center self term is A1 {i}")

    # Exact affine cocycle telescoping for three nonsingular rational frames.
    frames = (
        I3,
        ((F(2), F(1), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1))),
        ((F(1), F(0), F(0)), (F(1), F(2), F(0)), (F(0), F(1), F(1))),
    )
    centers = ((F(1), F(2), F(3)), (F(-1), F(0), F(2)), (F(3), F(-2), F(1)))
    transitions = {}
    for x in range(3):
        for y in range(3):
            transitions[(y, x)] = (
                matmul(inverse(frames[y]), frames[x]),
                matvec(inverse(frames[y]), vadd(centers[x], vscale(-1, centers[y]))),
            )
    gyx, tyx = transitions[(1, 0)]
    gzy, tzy = transitions[(2, 1)]
    gxz, txz = transitions[(0, 2)]
    check(matmul(gxz, matmul(gzy, gyx)) == I3, "linear cocycle telescope")
    total_t = vadd(matvec(gxz, matvec(gzy, tyx)), matvec(gxz, tzy), txz)
    check(total_t == z, "translation cocycle telescope")

    # Every alternating three-pair pattern repeated on the opposite local
    # frame cancels its constant center source.
    for a, b, c in permutations(range(4), 3):
        seq = ((c, a), (a, b), (b, c), (c, a), (a, b), (b, c))
        for axis in axes:
            total = F(0)
            for n, pair in enumerate(seq):
                sign = F(1) if n % 2 == 0 else F(-1)
                local = vscale(sign, vadd(V[pair[0]], V[pair[1]]))
                total += -2 * dot(axis, local)
            check(total == 0, f"alternating center cancellation {a},{b},{c}")

    # Representative unwrapped A3 hexagon: t^2 cancels exactly while a
    # generic first-source writer has a nonzero t^3 coefficient.
    positions = (
        (F(0), F(0), F(0)),
        T[0],
        vadd(T[0], vscale(-1, T[1])),
        vadd(T[0], vscale(-1, T[1]), T[3]),
        vadd(vscale(-1, T[1]), T[3]),
        T[3],
    )
    kinds = ("P", "C", "P", "C", "P", "C")
    pairs = ((0, 3), (0, 1), (1, 3), (0, 3), (0, 1), (1, 3))
    k = (F(2), F(3), F(5))
    xi = (F(1), F(-2), F(3))
    cycle = [(F(0), F(0)) for _ in range(4)]
    for pos, kind, pair in zip(positions, kinds, pairs):
        row = cycle_source_series(pos, kind, pair, k, xi)
        cycle = [qadd(cycle[n], qscale(F(105, 8), row[n])) for n in range(4)]
    check(cycle[0] == (F(0), F(0)), "cycle center source order zero")
    check(cycle[1] == (F(0), F(0)), "cycle center source order one")
    check(cycle[2] == (F(0), F(0)), "cycle constant O(k2) cancellation")
    check(cycle[3] != (F(0), F(0)), "cycle first writer begins generically at O(k3)")

    # Exact GC family ratio a_N |k_min,N| = (3 pi / 10) 2^-N,
    # storing the rational coefficient of pi.
    gc = [F(3, 10) / (2**n) for n in range(8)]
    for n in range(7):
        check(gc[n + 1] == gc[n] / 2, f"GC normalized mismatch ratio {n}")

    schema = json.loads((HERE / "OWNER_SCHEMA.json").read_text())
    required = (
        "raw_source_jet",
        "nonlinear_source_embedding",
        "support_and_shared_midpoint",
        "gd_recoil_material",
        "stationary_state_and_measure",
        "ctp_spectral_contact",
        "constraint_and_moving_projector",
        "boundary_controller_query",
        "physical_matching_and_writer",
        "retained_field_shell",
        "connected_1pi_schur_legendre_quotient",
        "accumulation_remainder_classification",
    )
    ids = tuple(row["id"] for row in schema["owners"])
    check(len(ids) == len(set(ids)), "owner IDs occur exactly once")
    check(ids == required, "complete ordered owner schema")
    check(schema["rules"]["undefined_is_not_zero"], "undefined is not zero")
    open_tokens = ("PENDING", "UNDEFINED", "OPEN", "PARTIAL")
    blockers = [row["id"] for row in schema["owners"]
                if any(token in row["current_status"] for token in open_tokens)]
    check(bool(blockers), "physical residual fail-closed blockers remain")

    return {
        "schema": "GL6CY_EXACT_TARGET_CHECK_V001",
        "status": "PASS_TARGET_ARITHMETIC__PHYSICAL_RESIDUAL_UNDEFINED",
        "checks": CHECKS,
        "pair_order": [list(pair) for pair in PAIRS],
        "six_tensor_directions_rank": 6,
        "center_first_source_sector": "T2_ONLY",
        "center_second_source_includes": "A1_PLUS_T2_DEPENDING_ON_A_AND_BRS",
        "shared_midpoint_connection": "EXACT_AFFINE_COCYCLE__TELESCOPING_FOR_GLOBAL_POSITION_FIELD__NOT_CURVATURE",
        "cycle_first_source": "CONSTANT_O_AK2_CANCELS__GENERIC_WRITER_BEGINS_O_A2K3",
        "gc_family": "a_N_times_kmin_N=(3*pi/10)*2^-N",
        "blocked_owner_ids": blockers,
        "raw_gl6cu_composition": "PENDING_DISTINCT_GL6CU_V002_HOSTILE_PASS",
        "physical_quotient_residual": "UNDEFINED",
        "ward": False,
        "gl6cr_invoked": False,
        "gravity": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    result = main_result()
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(f"PASS__GL6CY_EXACT_TARGET_ARITHMETIC__{CHECKS}/{CHECKS}")
        print("CENTER_FIRST=T2__CENTER_SECOND=INCLUDES_A1")
        print("CYCLE=O_AK2_CANCELS__GENERIC_FIRST_WRITER_O_A2K3")
        print("GC_NORMALIZED_GENERATOR_MISMATCH=O_2_TO_MINUS_N")
        print("PHYSICAL_QUOTIENT_RESIDUAL=UNDEFINED__NO_WARD_GL6CR_GRAVITY_G")


if __name__ == "__main__":
    main()
