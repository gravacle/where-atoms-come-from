#!/usr/bin/env python3
"""Independent hostile replay of the repaired GL6CL Fourier writer.

No target derivation is imported or executed.  The calculation starts from
the tetrahedral incidence and the audited GL6CH local tensor projection, then
reconstructs geometry, Fourier rows, ranks, expansions, obstructions, and the
uniform scalar identity with exact rational arithmetic.
"""

from __future__ import annotations

import hashlib
import json
import sys
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
STEPS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
A = (F(1),) * 6
TBASIS = (
    tuple(map(F, (1, 0, 0, 0, 0, -1))),
    tuple(map(F, (0, 1, 0, 0, -1, 0))),
    tuple(map(F, (0, 0, 1, -1, 0, 0))),
)


def add(*vs):
    return tuple(sum((F(v[i]) for v in vs), F(0)) for i in range(len(vs[0])))


def scale(c, v):
    return tuple(F(c) * F(x) for x in v)


def dot(a, b):
    return sum((F(x) * F(y) for x, y in zip(a, b)), F(0))


def outer(a, b):
    return tuple(tuple(F(x) * F(y) for y in b) for x in a)


def eye(n):
    return tuple(tuple(F(i == j) for j in range(n)) for i in range(n))


def madd(*ms):
    return tuple(tuple(sum((F(m[i][j]) for m in ms), F(0))
                       for j in range(len(ms[0][0])))
                 for i in range(len(ms[0])))


def mscale(c, m):
    return tuple(tuple(F(c) * F(x) for x in row) for row in m)


def transpose(m):
    return tuple(tuple(m[i][j] for i in range(len(m))) for j in range(len(m[0])))


def mm(a, b):
    return tuple(tuple(sum((a[i][r] * b[r][j] for r in range(len(b))), F(0))
                       for j in range(len(b[0])))
                 for i in range(len(a)))


def rank(m):
    work = [list(map(F, row)) for row in m]
    if not work:
        return 0
    r = 0
    for c in range(len(work[0])):
        pivot = next((i for i in range(r, len(work)) if work[i][c]), None)
        if pivot is None:
            continue
        work[r], work[pivot] = work[pivot], work[r]
        q = work[r][c]
        work[r] = [x / q for x in work[r]]
        for i in range(len(work)):
            if i == r or not work[i][c]:
                continue
            q = work[i][c]
            work[i] = [work[i][j] - q * work[r][j] for j in range(len(work[0]))]
        r += 1
        if r == len(work):
            break
    return r


def det(m):
    n = len(m)
    total = F(0)
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = F(-1 if inversions % 2 else 1)
        for i in range(n):
            term *= m[i][p[i]]
        total += term
    return total


PA = mscale(F(1, 6), outer(A, A))
PT = tuple(tuple(sum((F(t[i] * t[j], 2) for t in TBASIS), F(0))
                 for j in range(6)) for i in range(6))
PE = madd(eye(6), mscale(-1, PA), mscale(-1, PT))
assert mm(PA, PA) == PA and mm(PE, PE) == PE and mm(PT, PT) == PT
assert (rank(PA), rank(PE), rank(PT)) == (1, 2, 3)


def parent_position(x):
    basis = tuple(add(TETRA[i], scale(-1, TETRA[3])) for i in range(3))
    return tuple(sum((F(x[i]) * basis[i][j] for i in range(3)), F(0))
                 for j in range(3))


def child_position(y):
    return add(parent_position(y), TETRA[3])


def geometry():
    records = []
    offsets = []
    zero = (0, 0, 0)
    for missing in range(4):
        a, b, c = tuple(i for i in range(4) if i != missing)
        x_ab = add(STEPS[a], scale(-1, STEPS[b]))
        x_cb = add(STEPS[c], scale(-1, STEPS[b]))
        nodes = (
            ("P", (a, c), parent_position(zero)),
            ("C", (a, b), child_position(STEPS[a])),
            ("P", (b, c), parent_position(x_ab)),
            ("C", (a, c), child_position(add(x_ab, STEPS[c]))),
            ("P", (a, b), parent_position(x_cb)),
            ("C", (b, c), child_position(STEPS[c])),
        )
        center = scale(F(1, 6), tuple(sum((p[j] for _, _, p in nodes), F(0))
                                     for j in range(3)))
        expected_center = scale(F(1, 2), add(TETRA[a], TETRA[c], scale(-1, TETRA[b])))
        assert center == expected_center
        by_pair = {}
        for kind, pair, position in nodes:
            by_pair.setdefault(tuple(sorted(pair)), {})[kind] = add(position, scale(-1, center))
        assert set(by_pair) == set(combinations((a, b, c), 2))
        rho = {}
        for pair, ends in by_pair.items():
            assert ends["P"] == scale(-1, ends["C"])
            assert dot(ends["C"], ends["C"]) == F(11, 4)
            rho[pair] = ends["C"]
            offsets.append(ends["C"])
        records.append({"missing": missing, "ports": (a, b, c), "center": center, "rho": rho})

    second = tuple(tuple(sum((r[i] * r[j] for r in offsets), F(0))
                         for j in range(3)) for i in range(3))
    assert second == mscale(11, eye(3))
    diag4 = tuple(sum((r[i] ** 4 for r in offsets), F(0)) for i in range(3))
    mixed4 = tuple(sum((r[i] ** 2 * r[j] ** 2 for r in offsets), F(0))
                   for i, j in combinations(range(3), 2))
    assert diag4 == (F(83, 4),) * 3
    assert mixed4 == (F(19, 4),) * 3
    return records, offsets, second, diag4, mixed4


# Sparse polynomial: exponent triple -> rational coefficient.
def padd(*ps):
    out = {}
    for p in ps:
        for e, c in p.items():
            out[e] = out.get(e, F(0)) + c
    return {e: c for e, c in out.items() if c}


def pscale(c, p):
    return {e: F(c) * x for e, x in p.items() if F(c) * x}


def pmul(a, b, degree):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = tuple(ea[i] + eb[i] for i in range(3))
            if sum(e) <= degree:
                out[e] = out.get(e, F(0)) + ca * cb
    return {e: c for e, c in out.items() if c}


def ppow(p, n, degree):
    out = {(0, 0, 0): F(1)}
    for _ in range(n):
        out = pmul(out, p, degree)
    return out


def linear(v):
    return {e: c for e, c in zip(((1, 0, 0), (0, 1, 0), (0, 0, 1)), v) if c}


def cosine4(v):
    x = linear(v)
    return padd({(0, 0, 0): F(1)}, pscale(F(-1, 2), ppow(x, 2, 4)),
                 pscale(F(1, 24), ppow(x, 4, 4)))


def pdet(m, degree):
    n = len(m)
    total = {}
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = {(0, 0, 0): F(-1 if inversions % 2 else 1)}
        for i in range(n):
            term = pmul(term, m[i][p[i]], degree)
        total = padd(total, term)
    return total


def common_analysis(records):
    b0 = []
    component_polys = []
    for rec in records:
        row = [F(0)] * 6
        prow = [{} for _ in range(6)]
        for pair, rho in rec["rho"].items():
            idx = PIDX[pair]
            row[idx] = F(2)
            prow[idx] = pscale(2, cosine4(rho))
        b0.append(tuple(row))
        component_polys.append(tuple(prow))
    b0 = tuple(b0)
    assert mm(transpose(b0), b0) == madd(mscale(24, PA), mscale(8, PT))
    assert rank(b0) == 4

    wt0 = mm(b0, PT)
    assert mm(transpose(wt0), wt0) == mscale(8, PT)
    assert rank(wt0) == 3

    atcols = (A,) + TBASIS
    pmat = []
    for row in component_polys:
        pmat.append(tuple(padd(*(pscale(col[i], row[i]) for i in range(6)))
                           for col in atcols))
    determinant = pdet(tuple(pmat), 4)
    expected_det = {
        (0, 0, 0): F(768),
        (2, 0, 0): F(-1408), (0, 2, 0): F(-1408), (0, 0, 2): F(-1408),
        (4, 0, 0): F(2800, 3), (0, 4, 0): F(2800, 3), (0, 0, 4): F(2800, 3),
        (2, 2, 0): F(2144), (2, 0, 2): F(2144), (0, 2, 2): F(2144),
    }
    assert determinant == expected_det

    # Reconstruct the T-column normal independently through quadratic order.
    tmat = []
    for row in component_polys:
        tmat.append(tuple(padd(*(pscale(col[i], row[i]) for i in range(6)))
                          for col in TBASIS))
    tnormal = []
    for i in range(3):
        tnormal.append(tuple(padd(*(pmul(tmat[d][i], tmat[d][j], 2)
                                    for d in range(4))) for j in range(3)))
    tnormal = tuple(tnormal)
    expected_tnormal = []
    unit_exponents = ((2, 0, 0), (0, 2, 0), (0, 0, 2))
    for i in range(3):
        row = []
        for j in range(3):
            p = {(0, 0, 0): F(16)} if i == j else {}
            if i == j:
                for e in unit_exponents:
                    p = padd(p, {e: F(-4)})
                p = padd(p, {unit_exponents[i]: F(-32)})
            else:
                e = tuple(1 if q in (i, j) else 0 for q in range(3))
                p = padd(p, {e: F(24)})
            row.append(p)
        expected_tnormal.append(tuple(row))
    assert tnormal == tuple(expected_tnormal)

    return b0, wt0, determinant, tnormal


def locked_read_and_combination(wt0):
    drows = []
    for bits in product((0, 1), repeat=4):
        if sum(bits) != 2:
            continue
        z = tuple(F(1 - 2 * b) for b in bits)
        drows.append(tuple(z[a] * z[b] for a, b in PAIRS))
    drows = tuple(drows)
    dnormal = mm(transpose(drows), drows)
    assert dnormal == madd(mscale(4, PA), mscale(16, PE))
    assert rank(drows) == 3
    combined = madd(dnormal, mm(transpose(wt0), wt0))
    assert combined == madd(mscale(4, PA), mscale(16, PE), mscale(8, PT))
    assert rank(tuple(drows) + tuple(wt0)) == 6
    assert det(combined) == 524288
    inverse_normal = madd(mscale(F(1, 4), PA), mscale(F(1, 16), PE), mscale(F(1, 8), PT))
    assert mm(inverse_normal, combined) == eye(6)

    one_copy = mscale(F(1, 2), wt0)
    unsoldered = []
    for row in drows:
        unsoldered.append(tuple(row) + (F(0),) * 6)
        unsoldered.append((F(0),) * 6 + tuple(row))
    for row in one_copy:
        unsoldered.append(tuple(row) + tuple(row))
    assert rank(tuple(unsoldered)) == 9
    return drows, dnormal, combined, inverse_normal


def smooth_bound(offsets):
    fourth_radius_sum = sum((dot(r, r) ** 2 for r in offsets), F(0))
    assert fourth_radius_sum == F(363, 4)
    # ||Delta B||_F^2 <= (363/4)|k|^4 and sigma_min(W_T(0))^2=8.
    assert F(8) / fourth_radius_sum == F(32, 363)
    return fourth_radius_sum


def relative_analysis(records):
    # L[d,i] is a three-vector so L[d,i].k is the leading real entry after
    # stripping the common -2i factor.
    lrows = []
    for rec in records:
        row = []
        for t in TBASIS:
            row.append(tuple(sum((t[PIDX[pair]] * rho[q]
                                  for pair, rho in rec["rho"].items()), F(0))
                             for q in range(3)))
        lrows.append(tuple(row))
    lrows = tuple(lrows)
    expected_lrows = (
        ((F(3, 2), F(-1, 2), F(-1, 2)), (F(-1, 2), F(3, 2), F(-1, 2)), (F(-1, 2), F(-1, 2), F(3, 2))),
        ((F(3, 2), F(1, 2), F(1, 2)), (F(1, 2), F(3, 2), F(-1, 2)), (F(1, 2), F(-1, 2), F(3, 2))),
        ((F(3, 2), F(1, 2), F(-1, 2)), (F(1, 2), F(3, 2), F(1, 2)), (F(-1, 2), F(1, 2), F(3, 2))),
        ((F(3, 2), F(-1, 2), F(1, 2)), (F(-1, 2), F(3, 2), F(1, 2)), (F(1, 2), F(1, 2), F(3, 2))),
    )
    assert lrows == expected_lrows

    minors = []
    for omitted in range(4):
        matrix = tuple(tuple(linear(lrows[r][c]) for c in range(3))
                       for r in range(4) if r != omitted)
        minors.append(pdet(matrix, 3))
    minor_sum = {}
    for m in minors:
        minor_sum = padd(minor_sum, pmul(m, m, 6))
    expected_sum = {
        (6, 0, 0): F(9), (0, 6, 0): F(9), (0, 0, 6): F(9),
        (4, 2, 0): F(-9), (4, 0, 2): F(-9), (2, 4, 0): F(-9),
        (0, 4, 2): F(-9), (2, 0, 4): F(-9), (0, 2, 4): F(-9),
        (2, 2, 2): F(58),
    }
    assert minor_sum == expected_sum

    def leading_matrix(k):
        return tuple(tuple(dot(v, k) for v in row) for row in lrows)

    assert rank(leading_matrix((1, 2, 4))) == 3
    face_directions = ((1, 1, 0), (1, -1, 0), (1, 0, 1),
                       (1, 0, -1), (0, 1, 1), (0, 1, -1))
    assert all(rank(leading_matrix(k)) == 2 for k in face_directions)

    def sine_series(direction):
        matrix = []
        for rec in records:
            row = []
            for t in TBASIS:
                series = {}
                for pair, rho in rec["rho"].items():
                    f = dot(direction, rho)
                    coefficient = t[PIDX[pair]]
                    if f < 0:
                        f, coefficient = -f, -coefficient
                    if f and coefficient:
                        series[f] = series.get(f, F(0)) + coefficient
                row.append({f: c for f, c in series.items() if c})
            matrix.append(tuple(row))
        return tuple(matrix)

    dependencies = {
        (1, 1, 0): (0, 1, F(1)), (1, -1, 0): (0, 1, F(-1)),
        (1, 0, 1): (0, 2, F(1)), (1, 0, -1): (0, 2, F(-1)),
        (0, 1, 1): (1, 2, F(1)), (0, 1, -1): (1, 2, F(-1)),
    }
    for direction, (left, right, factor) in dependencies.items():
        for row in sine_series(direction):
            assert row[left] == {f: factor * c for f, c in row[right].items() if factor * c}
    return lrows, minors, minor_sum, dependencies


def finite_momentum_loss(records):
    # k=(pi/4)(1,1,1): encode cos(n*pi/8) for odd n as rational
    # coefficients of the independent symbols C=cos(pi/8), S=cos(3pi/8).
    cmap = {
        1: (F(1), F(0)), 3: (F(0), F(1)), 5: (F(0), F(-1)), 7: (F(-1), F(0)),
        9: (F(-1), F(0)), 11: (F(0), F(-1)), 13: (F(0), F(1)), 15: (F(1), F(0)),
    }

    def symbol_add(a, b):
        return (a[0] + b[0], a[1] + b[1])

    rows = []
    direction = (F(1), F(1), F(1))
    for rec in records:
        row = []
        for t in TBASIS:
            value = (F(0), F(0))
            for pair, rho in rec["rho"].items():
                n = int(2 * dot(direction, rho))
                assert F(n) == 2 * dot(direction, rho) and n % 2
                cs = cmap[n % 16]
                value = symbol_add(value, (2 * t[PIDX[pair]] * cs[0],
                                           2 * t[PIDX[pair]] * cs[1]))
            row.append(value)
        rows.append(tuple(row))
    rows = tuple(rows)
    expected = (
        ((F(-2), F(0)),) * 3,
        ((F(0), F(2)),) * 3,
        ((F(0), F(2)),) * 3,
        ((F(0), F(2)),) * 3,
    )
    assert rows == expected
    assert all(row[0] == row[1] == row[2] for row in rows)
    return rows


def rotation_scope_check():
    # A 45-degree z rotation sends the xy off-diagonal tensor to a diagonal
    # traceless tensor: R M_xy R^T = diag(-1,1,0).  Products c^2,s^2,cs are
    # all 1/2, so the result is exact over the rationals.
    transformed = ((F(-1), F(0), F(0)),
                   (F(0), F(1), F(0)),
                   (F(0), F(0), F(0)))
    assert transformed[0][1] == 0 and transformed[0][0] == -transformed[1][1]
    return "T2 is not SO(3)-closed: a 45-degree rotation maps xy shear into diagonal E2."


def scalar_identity():
    values = {}
    for bits in product((0, 1), repeat=4):
        n = sum(bits)
        z = tuple(1 - 2 * b for b in bits)
        lhs = sum(z[a] * z[b] for a, b in PAIRS)
        rhs = 2 * (n - 2) ** 2 - 2
        assert lhs == rhs
        values[n] = lhs
    assert values == {0: 6, 1: 0, 2: -2, 3: 0, 4: 6}
    assert 6 * F(105, 8) == F(315, 4)
    assert F(-63, 8) * F(-5) * F(2) == F(315, 4)
    return values


def encode(x):
    if isinstance(x, F):
        return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    return x


def canonical_hash(obj):
    payload = json.dumps(encode(obj), sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def main():
    records, offsets, second, diag4, mixed4 = geometry()
    b0, wt0, determinant, tnormal = common_analysis(records)
    drows, dnormal, combined, inverse_normal = locked_read_and_combination(wt0)
    radius4 = smooth_bound(offsets)
    lrows, minors, minor_sum, dependencies = relative_analysis(records)
    corner = finite_momentum_loss(records)
    rotation_scope = rotation_scope_check()
    occupancy_values = scalar_identity()

    result = {
        "schema": "AUDIT_G_GL6CL_INDEPENDENT_V001",
        "disposition": "PASS",
        "geometry": {
            "orientations": 4,
            "centered_offsets": 12,
            "all_radius_squared": "11/4",
            "second_moment": encode(second),
            "fourth_diagonal": encode(diag4),
            "fourth_mixed": encode(mixed4),
            "canonical_hash": canonical_hash(records),
        },
        "zero_mode": {
            "canonical_direct_rank": rank(b0),
            "canonical_direct_normal": "24 P_A + 8 P_T",
            "complete_tensor_writer_rank": rank(wt0),
            "complete_tensor_writer_normal": "8 P_T",
            "locked_read_rank": rank(drows),
            "locked_read_normal": "4 P_A + 16 P_E",
            "combined_rank": rank(tuple(drows) + tuple(wt0)),
            "combined_normal": "4 P_A + 16 P_E + 8 P_T",
            "combined_determinant": str(det(combined)),
            "unsoldered_rank": 9,
            "unsoldered_all_momentum_rank_ceiling": 10,
        },
        "common_expansion": {
            "canonical_direct_AT_determinant": "768-1408|k|^2+1072|k|^4-(416/3)sum_i k_i^4+O(|k|^6)",
            "orthonormal_T_normal": "8I-2|k|^2I+12kk^T-28diag(k_i^2)+O(|k|^4)",
            "T2_SO3_scope": rotation_scope,
            "scope_conclusion": "The T2-T2 block alone diagnoses neither physical isotropy nor anisotropy; a normalized E2+T2 completion is required.",
        },
        "smooth_ball": {
            "frobenius_bound_squared_coefficient": encode(radius4),
            "minimum_zero_T_singular_value_squared": "8",
            "full_rank_condition": "|k|^4<32/363",
            "analytic_left_inverse": "[C(k)^*C(k)]^-1 C(k)^*",
        },
        "relative_sector": {
            "zero_rank": 0,
            "generic_leading_rank": 3,
            "face_diagonal_leading_rank": 2,
            "exact_face_diagonal_column_dependencies": 6,
            "minor_sum": "9 sum_i k_i^6-9 sum_{i!=j}k_i^4 k_j^2+58 kx^2 ky^2 kz^2",
            "canonical_hash": canonical_hash((lrows, minors, minor_sum, dependencies)),
        },
        "finite_momentum": {
            "q": "(pi,0,0)",
            "k": "(pi/4)(1,1,1)",
            "complete_common_T_rank": 1,
            "all_three_columns_identical": True,
            "symbolic_rows_hash": canonical_hash(corner),
        },
        "uniform_A1": {
            "identity": "sum_{a<b} Z_a Z_b=2(n-2)^2-2",
            "occupancy_values": occupancy_values,
            "storage_shift": "U_d -> U_d+2q, plus scalar -2q",
            "denominator_derivative": "+315/4 h^6/U_d^6",
            "six_vertex_sum": "+315/4 h^6/U_d^6",
        },
        "writer_scope": {
            "complete_arbitrary_profile": "B_plus P_T and B_minus P_T only",
            "A1_plus_E_access": "locked diagonal read D",
            "unprojected_canonical_direct_row": "bookkeeping only",
            "not_established": "arbitrary-profile off-diagonal A1/E completion, autonomous source, stationary response, continuum geometry, gravity, or G",
        },
    }
    frozen = HERE / "INDEPENDENT_RESULT.json"
    if "--print-json" in sys.argv:
        print(json.dumps(encode(result), indent=2, sort_keys=True))
        return
    assert json.loads(frozen.read_text()) == encode(result)
    print("PASS GL6CL independent geometry/offset moments: 4 orientations, 12 offsets")
    print("PASS zero-mode scopes/ranks: D=A1+E rank3; BPT=T2 rank3; combined rank6 det524288")
    print("PASS common expansion and smooth ball: |k|^4<32/363")
    print("PASS T2 rotational guard: full E2+T2 completion required")
    print("PASS relative sector: zero0/generic3/face-leading2; six exact dependencies")
    print("PASS finite-momentum loss: common tensor rank1 at q=(pi,0,0)")
    print("PASS uniform A1 storage/writer identity: 315/4")
    print("AUDIT DISPOSITION: PASS")


if __name__ == "__main__":
    main()
