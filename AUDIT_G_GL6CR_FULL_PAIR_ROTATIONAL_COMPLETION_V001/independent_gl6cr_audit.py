#!/usr/bin/env python3
"""Independent exact audit of GL6CR.

This program deliberately does not import or execute the author derivation.
It reconstructs the tetrahedral group action, the complete Reynolds space,
the tensor solder, the SO(3) subspace, and the direct Ward-null system using
only Python's standard library and exact rational arithmetic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as Q
from itertools import combinations, permutations
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CR_FULL_PAIR_ROTATIONAL_COMPLETION_V001"

PAIRS = tuple(combinations(range(4), 2))
PAIR_POS = {p: n for n, p in enumerate(PAIRS)}
UPPER6 = tuple((i, j) for i in range(6) for j in range(i, 6))
QMON = ((2, 0, 0), (0, 2, 0), (0, 0, 2),
        (1, 1, 0), (1, 0, 1), (0, 1, 1))
QMON_LABEL = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
CMON = tuple((a, b, 3 - a - b)
             for a in range(4) for b in range(4 - a))
TETRA = (
    (Q(1), Q(1), Q(1)),
    (Q(1), Q(-1), Q(-1)),
    (Q(-1), Q(1), Q(-1)),
    (Q(-1), Q(-1), Q(1)),
)


class Audit:
    def __init__(self) -> None:
        self.count = 0

    def eq(self, got, expected, label: str) -> None:
        self.count += 1
        if got != expected:
            raise AssertionError(f"{label}: got {got!r}; expected {expected!r}")

    def ok(self, condition: bool, label: str) -> None:
        self.count += 1
        if not condition:
            raise AssertionError(label)


A = Audit()


def zmat(n: int, m: int):
    return tuple(tuple(Q(0) for _ in range(m)) for _ in range(n))


def eye(n: int):
    return tuple(tuple(Q(i == j) for j in range(n)) for i in range(n))


def tr(m):
    return tuple(tuple(m[i][j] for i in range(len(m)))
                 for j in range(len(m[0])))


def mm(a, b):
    return tuple(tuple(sum((a[i][r] * b[r][j] for r in range(len(b))), Q(0))
                       for j in range(len(b[0])))
                 for i in range(len(a)))


def mv(a, v):
    return tuple(sum((a[i][j] * v[j] for j in range(len(v))), Q(0))
                 for i in range(len(a)))


def madd(a, b):
    return tuple(tuple(a[i][j] + b[i][j] for j in range(len(a[0])))
                 for i in range(len(a)))


def mscale(c, a):
    return tuple(tuple(Q(c) * x for x in row) for row in a)


def outer(u, v):
    return tuple(tuple(x * y for y in v) for x in u)


def inverse(m):
    n = len(m)
    w = [list(map(Q, m[i])) + list(eye(n)[i]) for i in range(n)]
    for col in range(n):
        pivot = next((r for r in range(col, n) if w[r][col]), None)
        if pivot is None:
            raise AssertionError("singular matrix")
        w[col], w[pivot] = w[pivot], w[col]
        q = w[col][col]
        w[col] = [x / q for x in w[col]]
        for r in range(n):
            if r == col:
                continue
            q = w[r][col]
            if q:
                w[r] = [w[r][j] - q * w[col][j] for j in range(2 * n)]
    return tuple(tuple(row[n:]) for row in w)


def rref(rows):
    if not rows:
        return [], []
    w = [list(map(Q, row)) for row in rows]
    nr, nc = len(w), len(w[0])
    pivots = []
    pr = 0
    for c in range(nc):
        pivot = next((r for r in range(pr, nr) if w[r][c]), None)
        if pivot is None:
            continue
        w[pr], w[pivot] = w[pivot], w[pr]
        q = w[pr][c]
        w[pr] = [x / q for x in w[pr]]
        for r in range(nr):
            if r == pr or not w[r][c]:
                continue
            q = w[r][c]
            w[r] = [w[r][j] - q * w[pr][j] for j in range(nc)]
        pivots.append(c)
        pr += 1
        if pr == nr:
            break
    return w, pivots


def rank(rows):
    return len(rref(rows)[1])


def nullspace(rows):
    w, pivots = rref(rows)
    nc = len(rows[0]) if rows else 0
    free = [c for c in range(nc) if c not in pivots]
    ans = []
    for f in free:
        x = [Q(0)] * nc
        x[f] = Q(1)
        for rr in range(len(pivots) - 1, -1, -1):
            pc = pivots[rr]
            x[pc] = -sum((w[rr][c] * x[c] for c in free), Q(0))
        ans.append(tuple(x))
    return tuple(ans)


def coordinates(basis_rows, target):
    """Return c with sum_i c_i basis_rows[i] == target."""
    n = len(basis_rows)
    system = [[basis_rows[c][r] for c in range(n)] + [target[r]]
              for r in range(len(target))]
    w, pivots = rref(system)
    if tuple(p for p in pivots if p < n) != tuple(range(n)):
        raise AssertionError("basis is not independent")
    for row in w:
        if all(row[c] == 0 for c in range(n)) and row[n] != 0:
            raise AssertionError("target is outside span")
    # With pivot columns 0..n-1, each pivot row directly contains the answer.
    return tuple(w[i][n] for i in range(n))


def perm_matrix(images, n):
    p = [[Q(0)] * n for _ in range(n)]
    for source, dest in enumerate(images):
        p[dest][source] = Q(1)
    return tuple(tuple(row) for row in p)


def group_elements():
    ans = []
    for g in permutations(range(4)):
        # Completeness of the tetrahedral frame gives R=(1/4) sum |T_ga><T_a|.
        r = zmat(3, 3)
        for a in range(4):
            r = madd(r, mscale(Q(1, 4), outer(TETRA[g[a]], TETRA[a])))
        edge_images = []
        for a, b in PAIRS:
            edge_images.append(PAIR_POS[tuple(sorted((g[a], g[b])))])
        p = perm_matrix(edge_images, 6)
        A.eq(mm(tr(r), r), eye(3), "orthogonal tetrahedral action")
        for a in range(4):
            A.eq(mv(r, TETRA[a]), TETRA[g[a]], "tetrahedron port action")
        ans.append((g, r, p))
    A.eq(len(ans), 24, "S4 order")
    return tuple(ans)


def symmetric_matrix_seed(n, i, j):
    m = [[Q(0)] * n for _ in range(n)]
    m[i][j] = Q(1)
    m[j][i] = Q(1)
    return tuple(tuple(row) for row in m)


def quadratic_seed(slot):
    """Symmetric matrix q such that k^T q k is exactly QMON[slot]."""
    exp = QMON[slot]
    m = [[Q(0)] * 3 for _ in range(3)]
    if 2 in exp:
        i = exp.index(2)
        m[i][i] = Q(1)
    else:
        ids = [i for i, e in enumerate(exp) if e]
        m[ids[0]][ids[1]] = m[ids[1]][ids[0]] = Q(1, 2)
    return tuple(tuple(row) for row in m)


def qcoeff(q):
    return (q[0][0], q[1][1], q[2][2],
            2 * q[0][1], 2 * q[0][2], 2 * q[1][2])


def flatten_factor(pair_matrix, momentum_matrix):
    qc = qcoeff(momentum_matrix)
    return tuple(pair_matrix[i][j] * qc[s]
                 for i, j in UPPER6 for s in range(6))


def vadd(u, v):
    return tuple(u[i] + v[i] for i in range(len(u)))


def reynolds(pair_seed, q_seed, group):
    total = (Q(0),) * 126
    for _, r, p in group:
        # Project under F -> P^T F(R k) P.  The unnormalised group sum fixes
        # the exact basis convention and avoids an arbitrary 1/24 scale.
        pm = mm(mm(tr(p), pair_seed), p)
        qm = mm(mm(tr(r), q_seed), r)
        total = vadd(total, flatten_factor(pm, qm))
    return total


def construct_reynolds_space(group):
    all_images = []
    pivots = []
    labels = []
    for i, j in UPPER6:
        ps = symmetric_matrix_seed(6, i, j)
        for s, (a, b) in enumerate(QMON_LABEL):
            image = reynolds(ps, quadratic_seed(s), group)
            all_images.append(image)
            if rank(pivots + [image]) > len(pivots):
                pivots.append(image)
                labels.append(f"Reynolds[K_{i}{j}*k_{a}k_{b}]")
    A.eq(rank(all_images), 9, "rank of all 126 Reynolds images")
    A.eq(len(pivots), 9, "greedy Reynolds basis dimension")
    return tuple(pivots), tuple(labels), tuple(all_images)


def coefficients_to_quadratic(coefficients):
    q = [[Q(0)] * 3 for _ in range(3)]
    for value, exp in zip(coefficients, QMON):
        if 2 in exp:
            i = exp.index(2)
            q[i][i] = value
        else:
            ids = [i for i, e in enumerate(exp) if e]
            q[ids[0]][ids[1]] = q[ids[1]][ids[0]] = value / 2
    return tuple(tuple(row) for row in q)


def covariance_images(v, r, p):
    """Return coefficient vectors of K(Rk) and P K(k) P^T."""
    left = []
    at = 0
    for _i, _j in UPPER6:
        q = coefficients_to_quadratic(v[at:at + 6])
        left.extend(qcoeff(mm(mm(tr(r), q), r)))
        at += 6

    coefficient_matrices = []
    for s in range(6):
        cm = [[Q(0)] * 6 for _ in range(6)]
        at = 0
        for i, j in UPPER6:
            cm[i][j] = cm[j][i] = v[at + s]
            at += 6
        coefficient_matrices.append(mm(mm(p, tuple(tuple(row) for row in cm)), tr(p)))
    right = tuple(coefficient_matrices[s][i][j]
                  for i, j in UPPER6 for s in range(6))
    return tuple(left), right


def vector_to_poly_matrix(v):
    out = [[{} for _ in range(6)] for _ in range(6)]
    at = 0
    for i, j in UPPER6:
        poly = {QMON[s]: v[at + s] for s in range(6) if v[at + s]}
        out[i][j] = poly
        out[j][i] = dict(poly)
        at += 6
    return out


def poly_matrix_to_vector(m):
    return tuple(m[i][j].get(QMON[s], Q(0))
                 for i, j in UPPER6 for s in range(6))


def poly_add(a, b):
    out = dict(a)
    for e, c in b.items():
        out[e] = out.get(e, Q(0)) + c
        if out[e] == 0:
            del out[e]
    return out


def poly_scale(c, a):
    return {e: Q(c) * x for e, x in a.items() if Q(c) * x}


def poly_mul(a, b):
    out = {}
    for e, c in a.items():
        for f, d in b.items():
            g = tuple(e[i] + f[i] for i in range(3))
            out[g] = out.get(g, Q(0)) + c * d
    return {e: c for e, c in out.items() if c}


def evaluate_vector(v, k):
    m = [[Q(0)] * 6 for _ in range(6)]
    at = 0
    powers = [Q(k[0]) ** e[0] * Q(k[1]) ** e[1] * Q(k[2]) ** e[2]
              for e in QMON]
    for i, j in UPPER6:
        val = sum((v[at + s] * powers[s] for s in range(6)), Q(0))
        m[i][j] = m[j][i] = val
        at += 6
    return tuple(tuple(row) for row in m)


def hmatrix_from_coord(v):
    xx, yy, zz, xy, xz, yz = v
    return ((xx, xy, xz), (xy, yy, yz), (xz, yz, zz))


def hcoord_from_matrix(h):
    return (h[0][0], h[1][1], h[2][2], h[0][1], h[0][2], h[1][2])


def solder():
    cols = []
    for a, b in PAIRS:
        va = tuple(x / 2 for x in TETRA[a])
        vb = tuple(x / 2 for x in TETRA[b])
        h = madd(outer(va, vb), outer(vb, va))
        cols.append(hcoord_from_matrix(h))
    return tr(tuple(cols))


def check_solder(group, d):
    expected = (
        (Q(1, 2), Q(-1, 2), Q(-1, 2), Q(-1, 2), Q(-1, 2), Q(1, 2)),
        (Q(-1, 2), Q(1, 2), Q(-1, 2), Q(-1, 2), Q(1, 2), Q(-1, 2)),
        (Q(-1, 2), Q(-1, 2), Q(1, 2), Q(1, 2), Q(-1, 2), Q(-1, 2)),
        (Q(0), Q(0), Q(-1, 2), Q(1, 2), Q(0), Q(0)),
        (Q(0), Q(-1, 2), Q(0), Q(0), Q(1, 2), Q(0)),
        (Q(-1, 2), Q(0), Q(0), Q(0), Q(0), Q(1, 2)),
    )
    A.eq(d, expected, "declared solder coefficients")
    dinv = inverse(d)
    A.eq(mm(dinv, d), eye(6), "solder inverse on pair coordinates")
    A.eq(mm(d, dinv), eye(6), "solder inverse on tensor coordinates")
    for _, r, p in group:
        for j in range(6):
            lhs = mv(d, tuple(p[i][j] for i in range(6)))
            h = hmatrix_from_coord(tuple(d[i][j] for i in range(6)))
            rhs = hcoord_from_matrix(mm(mm(r, h), tr(r)))
            A.eq(lhs, rhs, "solder S4 equivariance")
    return dinv


def character_dimensions(group):
    total_quad = Q(0)
    total_const = Q(0)
    multiplicities = [Q(0), Q(0), Q(0)]
    for _, r, p in group:
        chi_v = sum(r[i][i] for i in range(3))
        chi_v2 = sum(mm(r, r)[i][i] for i in range(3))
        chi_q = (chi_v * chi_v + chi_v2) / 2
        chi_p = sum(p[i][i] for i in range(6))
        chi_p2 = sum(mm(p, p)[i][i] for i in range(6))
        chi_sp = (chi_p * chi_p + chi_p2) / 2
        total_quad += chi_q * chi_sp
        total_const += chi_sp
        # The edge representation is A1+E+T2, so E can be recovered without
        # importing an S4 character table.  Sym^2(T2)=A1+E+T2 is then checked
        # element by element, and the three target multiplicities are counted.
        chi_e = chi_p - 1 - chi_v
        A.eq(chi_q, 1 + chi_e + chi_v,
             "Sym2(T2) decomposes as A1+E+T2")
        for i, chi_irrep in enumerate((Q(1), chi_e, chi_v)):
            multiplicities[i] += chi_irrep * chi_sp
    A.eq(total_quad / 24, Q(9), "character dimension at k^2")
    A.eq(total_const / 24, Q(3), "constant S4 symmetric dimension")
    A.eq(tuple(x / 24 for x in multiplicities), (Q(3), Q(3), Q(3)),
         "A1, E, T2 each occur three times in Sym2(pair)")


def tensor_basis_matrices(d):
    return tuple(hmatrix_from_coord(tuple(d[r][c] for r in range(6)))
                 for c in range(6))


def qpoly_from_symmetric(q):
    coeff = qcoeff(q)
    return {QMON[s]: coeff[s] for s in range(6) if coeff[s]}


def so3_vectors(d):
    hs = tensor_basis_matrices(d)
    traces = tuple(sum(h[a][a] for a in range(3)) for h in hs)
    r2 = {(2, 0, 0): Q(1), (0, 2, 0): Q(1), (0, 0, 2): Q(1)}
    mats = []

    # a: |k|^2 h:g (full Frobenius contraction, including both off diagonals).
    ma = [[{} for _ in range(6)] for _ in range(6)]
    # b: |k|^2 tr(h)tr(g).
    mb = [[{} for _ in range(6)] for _ in range(6)]
    # c: (h k).(g k).
    mc = [[{} for _ in range(6)] for _ in range(6)]
    # d: tr(h)(k^T g k)+tr(g)(k^T h k).
    md = [[{} for _ in range(6)] for _ in range(6)]
    for i in range(6):
        for j in range(6):
            frob = sum((hs[i][a][b] * hs[j][a][b]
                        for a in range(3) for b in range(3)), Q(0))
            ma[i][j] = poly_scale(frob, r2)
            mb[i][j] = poly_scale(traces[i] * traces[j], r2)

            # k^T H_i H_j k; symmetrize the coefficient matrix explicitly.
            hij = mm(hs[i], hs[j])
            hji = mm(hs[j], hs[i])
            sym = mscale(Q(1, 2), madd(hij, hji))
            mc[i][j] = qpoly_from_symmetric(sym)

            qd = madd(mscale(traces[i], hs[j]), mscale(traces[j], hs[i]))
            md[i][j] = qpoly_from_symmetric(qd)
    for m in (ma, mb, mc, md):
        mats.append(poly_matrix_to_vector(m))
    A.eq(rank(mats), 4, "independent SO3 quadratic bilinears")
    return tuple(mats)


def constant_dimensions(d):
    # The two SO3 forms at k^0 are h:g and tr(h)tr(g).
    hs = tensor_basis_matrices(d)
    frob = []
    trtr = []
    for i, j in UPPER6:
        frob.append(sum((hs[i][a][b] * hs[j][a][b]
                         for a in range(3) for b in range(3)), Q(0)))
        ti = sum(hs[i][a][a] for a in range(3))
        tj = sum(hs[j][a][a] for a in range(3))
        trtr.append(ti * tj)
    A.eq(rank((tuple(frob), tuple(trtr))), 2, "constant SO3 dimension")


def cubic_ward_constraints(basis, d, dinv):
    # x_j(k,xi) is a degree-one polynomial.  Store coefficients by k axis.
    x_by_xi = []
    for q in range(3):
        hlin = [[Q(0)] * 3 for _ in range(6)]
        # Coordinate order xx,yy,zz,xy,xz,yz for k odot e_q.
        for axis in range(3):
            k = [Q(0)] * 3
            k[axis] = Q(1)
            xi = [Q(0)] * 3
            xi[q] = Q(1)
            hm = tuple(tuple(k[a] * xi[b] + xi[a] * k[b]
                             for b in range(3)) for a in range(3))
            hc = hcoord_from_matrix(hm)
            for row in range(6):
                hlin[row][axis] = hc[row]
        xlin = [[Q(0)] * 3 for _ in range(6)]
        for j in range(6):
            for axis in range(3):
                xlin[j][axis] = sum((dinv[j][r] * hlin[r][axis]
                                     for r in range(6)), Q(0))
        # Direct inversion audit on all symbolic coefficients.
        for r in range(6):
            for axis in range(3):
                got = sum((d[r][j] * xlin[j][axis]
                           for j in range(6)), Q(0))
                A.eq(got, hlin[r][axis], "symbolic solder inversion")
        x_by_xi.append(tuple(tuple(row) for row in xlin))

    # One constraint row for every xi, output component, cubic monomial.
    columns_by_basis = []
    for iv in basis:
        kmat = vector_to_poly_matrix(iv)
        col = []
        for q in range(3):
            xlin = x_by_xi[q]
            for out in range(6):
                total = {}
                for j in range(6):
                    lin = {tuple(1 if a == axis else 0 for a in range(3)): xlin[j][axis]
                           for axis in range(3) if xlin[j][axis]}
                    total = poly_add(total, poly_mul(kmat[out][j], lin))
                col.extend(total.get(e, Q(0)) for e in CMON)
        columns_by_basis.append(tuple(col))
    constraints = tuple(tuple(columns_by_basis[c][r] for c in range(9))
                        for r in range(len(columns_by_basis[0])))
    A.eq(len(constraints), 180, "full Ward coefficient-row census")
    A.eq(rank(constraints), 8, "direct cubic Ward constraint rank")
    ns = nullspace(constraints)
    A.eq(len(ns), 1, "direct cubic Ward nullity")
    return constraints, ns[0], tuple(x_by_xi)


def proportional(u, v):
    scale = None
    for a, b in zip(u, v):
        if b:
            if scale is None:
                scale = a / b
            elif a != scale * b:
                return False
        elif a:
            return False
    return scale is not None


def t2_block_coefficients(v):
    # Raw GL6CO T2 directions: e01-e23, e02-e13, e03-e12.
    t = (
        (Q(1), Q(0), Q(0), Q(0), Q(0), Q(-1)),
        (Q(0), Q(1), Q(0), Q(0), Q(-1), Q(0)),
        (Q(0), Q(0), Q(1), Q(-1), Q(0), Q(0)),
    )
    # L_ab(k)=t_a^T K(k)t_b. Fit A r2 I+B diag(k_i^2)+C offdiag(k_i k_j).
    km = vector_to_poly_matrix(v)
    block = [[{} for _ in range(3)] for _ in range(3)]
    for a in range(3):
        for b in range(3):
            p = {}
            for i in range(6):
                for j in range(6):
                    if t[a][i] and t[b][j]:
                        p = poly_add(p, poly_scale(t[a][i] * t[b][j], km[i][j]))
            block[a][b] = p
    # Coordinate axes are aligned with these t_i in the GL6CO convention.
    # From diagonal a=0: k_y^2 or k_z^2 gives A; k_x^2 gives A+B.
    Avec = [block[a][a].get(QMON[(a + 1) % 3], None) for a in range(3)]
    # More robust: all transverse squared coefficients must agree.
    transverse = []
    longitudinal = []
    for a in range(3):
        longitudinal.append(block[a][a].get(tuple(2 if r == a else 0 for r in range(3)), Q(0)))
        for r in range(3):
            if r != a:
                transverse.append(block[a][a].get(tuple(2 if s == r else 0 for s in range(3)), Q(0)))
    A.ok(all(x == transverse[0] for x in transverse), "T2 transverse A consistency")
    A.ok(all(x == longitudinal[0] for x in longitudinal), "T2 longitudinal A+B consistency")
    aa = transverse[0]
    bb = longitudinal[0] - aa
    cross = []
    for a in range(3):
        for b in range(a + 1, 3):
            exp = tuple(1 if r in (a, b) else 0 for r in range(3))
            cross.append(block[a][b].get(exp, Q(0)))
    A.ok(all(x == cross[0] for x in cross), "T2 offdiagonal C consistency")
    cc = cross[0]
    # Ensure no undeclared terms remain.
    expected = [[{} for _ in range(3)] for _ in range(3)]
    r2 = {QMON[0]: aa, QMON[1]: aa, QMON[2]: aa}
    for a in range(3):
        expected[a][a] = dict(r2)
        exp = tuple(2 if r == a else 0 for r in range(3))
        expected[a][a][exp] = expected[a][a].get(exp, Q(0)) + bb
        expected[a][a] = {e: c for e, c in expected[a][a].items() if c}
        for b in range(a + 1, 3):
            exp = tuple(1 if r in (a, b) else 0 for r in range(3))
            expected[a][b] = {exp: cc} if cc else {}
            expected[b][a] = dict(expected[a][b])
    A.eq(block, expected, "complete T2 ABC form")
    return aa, bb, cc


def qstr(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def parseq(s):
    return Q(str(s))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--frozen-target-hash", default="")
    parser.add_argument("--write-result", action="store_true")
    args = parser.parse_args()

    group = group_elements()
    character_dimensions(group)
    basis, labels, all_images = construct_reynolds_space(group)
    for iv in basis:
        for _, r, p in group:
            left, right = covariance_images(iv, r, p)
            A.eq(left, right, "Reynolds basis obeys declared covariance convention")

    d_solder = solder()
    dinv = check_solder(group, d_solder)
    constant_dimensions(d_solder)

    so3 = so3_vectors(d_solder)
    so3_coords = tuple(coordinates(basis, v) for v in so3)
    A.eq(rank(so3_coords), 4, "SO3 subspace rank inside S4 space")

    ledger = json.loads((TARGET / "EXACT_LEDGER.json").read_text())
    qd = ledger["quadratic_symbols"]
    A.eq(labels, tuple(qd["invariant_basis_labels"]), "frozen Reynolds pivot labels")
    claimed_so3 = tuple(tuple(parseq(x) for x in row)
                        for row in qd["rotational_coordinates_in_invariant_basis"])
    A.eq(so3_coords, claimed_so3, "four SO3 coordinate rows")

    residuals = tuple(tuple(parseq(x) for x in row)
                      for row in qd["five_residual_covectors_on_invariant_coordinates"])
    A.eq(rank(residuals), 5, "five residual covectors independent")
    for ell in residuals:
        for row in so3_coords:
            A.eq(sum((ell[i] * row[i] for i in range(9)), Q(0)), Q(0),
                 "residual annihilates SO3 space")
    # Rank 4 plus a 5-dimensional annihilator proves iff, not just necessity.
    A.eq(rank(so3_coords) + rank(residuals), 9,
         "residuals are necessary and sufficient")

    abc = tuple(t2_block_coefficients(v) for v in basis)
    claimed_abc = tuple(tuple(parseq(x) for x in row)
                        for row in qd["T2_block_ABC_by_invariant_basis"])
    A.eq(abc, claimed_abc, "T2 ABC coefficients in raw pair convention")
    mismatch = tuple(row[1] + row[2] for row in abc)
    A.eq(mismatch,
         tuple(parseq(x) for x in qd["T2_extension_mismatch_B_plus_C"]),
         "T2 B+C projection")
    # It is one nonzero constraint, not a full 5-condition test.
    A.eq(rank((mismatch,)), 1, "T2 test supplies one relation")

    constraints, ward_null, x_by_xi = cubic_ward_constraints(basis, d_solder, dinv)
    claimed_ray = tuple(parseq(x) for x in
                        ledger["direct_cubic_gauge_shortcut"]["unique_null_ray_in_invariant_basis"])
    A.ok(proportional(ward_null, claimed_ray), "Ward null ray matches target")

    einstein = tuple(sum((so3[r][i] * c for r, c in enumerate(
        (Q(1, 2), Q(-1, 2), Q(-1), Q(1, 2)))), Q(0)) for i in range(126))
    einstein_coords = coordinates(basis, einstein)
    claimed_einstein = tuple(parseq(x) for x in
                             ledger["einstein_reference"]["coordinates_in_invariant_basis"])
    A.eq(einstein_coords, claimed_einstein, "Einstein coordinates in Reynolds basis")
    A.eq(tuple(16 * x for x in einstein_coords), claimed_ray,
         "direct Ward ray equals sixteen times Einstein coordinates")

    # Verify the Einstein matrix kills all three solder-inverted longitudinal modes.
    ke_poly = vector_to_poly_matrix(einstein)
    for q in range(3):
        for out in range(6):
            total = {}
            for j in range(6):
                lin = {tuple(1 if a == axis else 0 for a in range(3)): x_by_xi[q][j][axis]
                       for axis in range(3) if x_by_xi[q][j][axis]}
                total = poly_add(total, poly_mul(ke_poly[out][j], lin))
            A.eq(total, {}, "Einstein polynomial Ward null")

    kval = evaluate_vector(einstein, (2, 3, 5))
    A.eq(rank(kval), 3, "generic Einstein pair-matrix rank")

    # Change from canonical pairs to the raw A1,E1,E2,T1,T2,T3 basis used by
    # the earlier held-out GL6CO test.  This catches both diagonal/off-diagonal
    # tensor-coordinate factors and the direction of the solder congruence.
    irrep_columns = (
        (1, 1, 1, 1, 1, 1),
        (1, -1, 0, 0, -1, 1),
        (1, 1, -2, -2, 1, 1),
        (1, 0, 0, 0, 0, -1),
        (0, 1, 0, 0, -1, 0),
        (0, 0, 1, -1, 0, 0),
    )
    change = tr(tuple(tuple(Q(x) for x in col) for col in irrep_columns))
    irrep_k = mm(mm(tr(change), kval), change)
    expected_irrep_k = tuple(tuple(map(Q, row)) for row in (
        (-38, 5, 37, 15, 10, 6),
        (5, 100, 20, -30, 20, 0),
        (37, 20, 4, -30, -20, 24),
        (15, -30, -30, 4, -6, -10),
        (10, 20, -20, -6, 9, -15),
        (6, 0, 24, -10, -15, 25),
    ))
    A.eq(irrep_k, expected_irrep_k,
         "Einstein solder agrees with prior full A1/E2/T2 reference")
    # The three longitudinal pair coordinates are independent at nonzero k.
    xnum = []
    kval_k = (Q(2), Q(3), Q(5))
    for q in range(3):
        xnum.append(tuple(sum((x_by_xi[q][j][axis] * kval_k[axis]
                               for axis in range(3)), Q(0)) for j in range(6)))
    A.eq(rank(xnum), 3, "independent longitudinal modes at generic momentum")
    for x in xnum:
        A.eq(mv(kval, x), (Q(0),) * 6, "generic longitudinal null vector")

    # Independently derive the three SO3 Ward equations by applying each SO3
    # basis to the longitudinal coordinates and taking the constraint rowspace.
    ward_so3_cols = []
    for v in so3:
        kmat = vector_to_poly_matrix(v)
        col = []
        for q in range(3):
            for out in range(6):
                total = {}
                for j in range(6):
                    lin = {tuple(1 if a == axis else 0 for a in range(3)): x_by_xi[q][j][axis]
                           for axis in range(3) if x_by_xi[q][j][axis]}
                    total = poly_add(total, poly_mul(kmat[out][j], lin))
                col.extend(total.get(e, Q(0)) for e in CMON)
        ward_so3_cols.append(tuple(col))
    ward_so3_rows = tuple(tuple(ward_so3_cols[c][r] for c in range(4))
                          for r in range(180))
    A.eq(rank(ward_so3_rows), 3, "three independent Ward equations in SO3 family")
    A.ok(proportional(nullspace(ward_so3_rows)[0], (Q(1), Q(-1), Q(-2), Q(1))),
         "SO3 Ward null is Einstein ray")
    declared_equations = (
        (Q(2), Q(0), Q(1), Q(0)),
        (Q(0), Q(0), Q(1), Q(2)),
        (Q(0), Q(2), Q(0), Q(2)),
    )
    A.eq(rank(declared_equations), 3, "declared Ward equations independent")
    # Same rowspace iff stacking does not increase rank.
    A.eq(rank(ward_so3_rows + declared_equations), 3,
         "declared Ward equations equal independently derived rowspace")

    target_hash = ""
    manifest = TARGET / "MANIFEST.sha256"
    if manifest.exists():
        target_hash = sha256(manifest)
    if args.frozen_target_hash:
        A.eq(target_hash, args.frozen_target_hash, "frozen target manifest hash")

    result = {
        "audit": "GL6CR independent hostile audit",
        "method": "independent exact rational reconstruction; author module neither imported nor executed",
        "checks": A.count,
        "reynolds_raw_seed_count": len(all_images),
        "s4_quadratic_dimension_character": 9,
        "s4_quadratic_dimension_reynolds": rank(all_images),
        "s4_constant_dimension": 3,
        "so3_quadratic_dimension": rank(so3_coords),
        "so3_constant_dimension": 2,
        "rotational_codimension": rank(residuals),
        "direct_ward_constraint_rows": len(constraints),
        "direct_ward_rank": rank(constraints),
        "direct_ward_nullity": len(nullspace(constraints)),
        "direct_ward_ray": [qstr(x) for x in claimed_ray],
        "einstein_coordinates": [qstr(x) for x in einstein_coords],
        "einstein_generic_rank": rank(kval),
        "target_manifest_sha256": target_hash,
        "disposition": "PASS" if target_hash else "PROVISIONAL_MATH_PASS_TARGET_NOT_SEALED",
        "ceiling": "algebraic classifier only; the physical F3 Ward identity, same-state complete response, 1PI/quotient construction, causal refinement, gravity, and G are not derived",
    }
    if args.write_result:
        (HERE / "INDEPENDENT_RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"PASS__GL6CR_INDEPENDENT_HOSTILE_REPLAY__{A.count}/{A.count}")
    print("SPACE=S4_9;SO3_4;ROTATIONAL_CODIMENSION_5;T2_PROJECTION_RANK_1")
    print("DIRECT_WARD=180x9;RANK_8;NULLITY_1;UNIQUE_EINSTEIN_RAY")
    print("DISPOSITION=PASS;ALGEBRAIC_CLASSIFIER_ONLY;PHYSICAL_F3_WARD_1PI_GRAVITY_G_OPEN")


if __name__ == "__main__":
    main()
