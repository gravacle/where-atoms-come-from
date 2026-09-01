#!/usr/bin/env python3
"""Exact/combinatorial replay for the mutable GL6AK packet."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent
CHECKS = 0


def check(condition: bool, message: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(message)
    CHECKS += 1


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


E = tuple(tuple(int(i == a) for i in range(4)) for a in range(4))
LABELS = range(4)
PAIRS = tuple(itertools.combinations(LABELS, 2))


def dist(x, y=(0, 0, 0, 0)):
    delta = sub(x, y)
    check(sum(delta) == 0, "A3 distance received a non-A3 displacement")
    return sum(abs(v) for v in delta) // 2


def internal_edge(x, a, b):
    return frozenset(((x, a), (x, b)))


def shared_edge(x, a, b):
    return frozenset(((x, a), (add(x, sub(E[a], E[b])), b)))


def incident_edges(x, a):
    result = set()
    for b in LABELS:
        if b == a:
            continue
        result.add(internal_edge(x, a, b))
        result.add(frozenset(((x, a), (add(x, sub(E[a], E[b])), b))))
    return result


def permute_x(x, sigma):
    out = [0] * 4
    for old, new in enumerate(sigma):
        out[new] = x[old]
    return tuple(out)


def permute_site(site, sigma):
    x, a = site
    return permute_x(x, sigma), sigma[a]


def matmul(a, b):
    n = len(a)
    return [[sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def matsub(a, b):
    return [[a[i][j] - b[i][j] for j in range(len(a))] for i in range(len(a))]


def eye(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def scale(a, c):
    return [[c * v for v in row] for row in a]


def tail(d, x):
    # Stable positive evaluation of T_d(x)=sum_{k>=d}x^k/k!.
    if x == 0:
        return 1.0 if d == 0 else 0.0
    log_term = d * math.log(x) - math.lgamma(d + 1)
    term = 0.0 if log_term < -745 else math.exp(log_term)
    if term == 0.0:
        return 0.0
    total = term
    k = d
    while term > max(1e-18, abs(total) * 1e-16):
        k += 1
        term *= x / k
        total += term
        if k > 10000:
            raise RuntimeError("tail did not converge")
    return total


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


# 1. Infinite graph census and unique oriented edge representatives.
origins = [
    (0, 0, 0, 0),
    (1, -1, 0, 0),
    (2, -1, -1, 0),
    (-1, -1, 1, 1),
]
for x in origins:
    check(sum(x) == 0, "sample origin must lie in A3")
    oriented = [internal_edge(x, a, b) for a, b in PAIRS]
    oriented += [shared_edge(x, a, b) for a, b in PAIRS]
    check(len(oriented) == 12, "each cell must anchor twelve pair terms")
    check(len(set(oriented)) == 12, "oriented pair representatives must be unique")
    containing = set()
    for a in LABELS:
        containing.update(incident_edges(x, a))
    check(len(containing) == 18, "a cell must touch six internal and twelve shared terms")
    for a in LABELS:
        inc = incident_edges(x, a)
        check(len(inc) == 6, "every infinite site must have degree six")
        same = sum(all(site[0] == x for site in edge) for edge in inc)
        cross = len(inc) - same
        check((same, cross) == (3, 3), "degree must split as three plus three")


# 2. Exact translation and S4 covariance of both edge families.
translations = [(0, 0, 0, 0), (1, -1, 0, 0), (2, -1, 0, -1)]
permutations = list(itertools.permutations(LABELS))
for x in origins:
    for a, b in PAIRS:
        for z in translations:
            lhs = frozenset((add(site[0], z), site[1]) for site in internal_edge(x, a, b))
            rhs = internal_edge(add(x, z), a, b)
            check(lhs == rhs, "internal edge translation covariance failed")
            lhs = frozenset((add(site[0], z), site[1]) for site in shared_edge(x, a, b))
            rhs = shared_edge(add(x, z), a, b)
            check(lhs == rhs, "shared edge translation covariance failed")
        for sigma in permutations:
            got_i = frozenset(permute_site(site, sigma) for site in internal_edge(x, a, b))
            want_i = internal_edge(permute_x(x, sigma), sigma[a], sigma[b])
            check(got_i == want_i, "internal edge S4 covariance failed")
            got_s = frozenset(permute_site(site, sigma) for site in shared_edge(x, a, b))
            xa = permute_x(x, sigma)
            aa, bb = sigma[a], sigma[b]
            want_s = frozenset(((xa, aa), (add(xa, sub(E[aa], E[bb])), bb)))
            check(got_s == want_s, "shared edge S4 covariance failed")


# 3. Finite-patch strict-interior embeddings.
patches = [
    {(0, 0, 0, 0)},
    {(0, 0, 0, 0), (1, -1, 0, 0), (0, 1, -1, 0)},
    {(2, -1, -1, 0), (-1, 2, 0, -1), (0, -1, 2, -1)},
]
for patch in patches:
    mins = [min(x[i] for x in patch) for i in LABELS]
    m = tuple(max(1, 1 - mins[i]) for i in LABELS)
    n = sum(m)
    for x in patch:
        image = add(m, x)
        check(all(v >= 1 for v in image), "embedded parent must be strict interior")
        check(sum(image) == n, "embedded parent must lie in one S_N layer")
        for a in LABELS:
            child = add(image, E[a])
            check(all(v >= 1 for v in child), "embedded child must remain positive")
            check(sum(child) == n + 1, "embedded child must lie in S_(N+1)")
        for a, b in PAIRS:
            y = add(x, sub(E[a], E[b]))
            left_child = add(add(m, x), E[a])
            right_child = add(add(m, y), E[b])
            check(left_child == right_child, "shared-child embedding failed")


# 4. Cell balls and interaction-shell census bounds.
for r in range(0, 7):
    points = []
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            for c in range(-r, r + 1):
                x = (a, b, c, -a - b - c)
                if dist(x) <= r:
                    points.append(x)
    check(len(points) <= (2 * r + 1) ** 3, "A3 ball cubic bound failed")
    check(len(points) == len(set(points)), "A3 ball enumeration duplicated a point")


# 5. Exact constant reduction and factorial boundary-tail behavior.
for j in (Fraction(1, 3), Fraction(2), Fraction(17, 5)):
    hbar = Fraction(7, 4)
    lam = 4 * j * 6 / hbar
    check(lam == 24 * j / hbar, "lambda=24J/hbar reduction failed")
    check((72 * j / hbar) / lam == 3, "boundary prefactor must reduce to three")

for x in (0.1, 0.75, 2.5):
    for d in range(0, 10):
        # Fundamental theorem identity checked by independent quadrature.
        steps = 20000
        t = 1.0
        acc = 0.0
        for k in range(steps + 1):
            u = t * k / steps
            weight = 0.5 if k in (0, steps) else 1.0
            acc += weight * tail(d, x * u)
        integral = acc * t / steps
        target = tail(d + 1, x) / x
        check(abs(integral - target) < 2e-8, "integrated factorial-tail identity failed")

def boundary_series(R, r_x, x):
    total = 0.0
    # Terms after this cutoff are far below double precision for tested x.
    for r in range(R, R + 220):
        total += (2 * r + 1) ** 3 * tail(r - r_x + 1, x)
    return total


for x in (0.2, 1.0, 4.0):
    vals = [boundary_series(r, 2, x) for r in (5, 10, 20, 40, 60)]
    check(all(vals[i + 1] < vals[i] for i in range(len(vals) - 1)), "tail must decrease")
    check(vals[-1] < 1e-20, "boundary series must vanish at large radius")


# 6. Concrete Folner sequence diagnostic under fixed Z3 shifts.
def folner_box(r):
    return {
        (a, b, c, -a - b - c)
        for a in range(-r, r + 1)
        for b in range(-r, r + 1)
        for c in range(-r, r + 1)
    }


for shift in ((1, -1, 0, 0), (2, -1, -1, 0)):
    ratios = []
    for r in (3, 7, 15):
        box = folner_box(r)
        moved = {add(x, shift) for x in box}
        ratios.append(len(box.symmetric_difference(moved)) / len(box))
    check(ratios[2] < ratios[1] < ratios[0], "Folner symmetric-difference ratio must decrease")


# 7. Exact multiplicity-free A1+E+T2 projectors in six-pair space.
adj = [[Fraction(0) for _ in PAIRS] for _ in PAIRS]
for i, p in enumerate(PAIRS):
    for j, q in enumerate(PAIRS):
        if i != j and set(p).intersection(q):
            adj[i][j] = Fraction(1)
i6 = eye(6)
a_minus_4i = matsub(adj, scale(i6, Fraction(4)))
a_plus_2i = [[adj[i][j] + 2 * i6[i][j] for j in range(6)] for i in range(6)]
p_a1 = [[Fraction(1, 6) for _ in range(6)] for _ in range(6)]
p_e = scale(matmul(adj, a_minus_4i), Fraction(1, 12))
p_t2 = scale(matmul(a_minus_4i, a_plus_2i), Fraction(-1, 8))
zero = [[Fraction(0) for _ in range(6)] for _ in range(6)]
for p, rank in ((p_a1, 1), (p_e, 2), (p_t2, 3)):
    check(matmul(p, p) == p, "irrep projector must be idempotent")
    check(sum(p[i][i] for i in range(6)) == rank, "irrep projector rank failed")
check(matmul(p_a1, p_e) == zero, "A1 and E projectors must be orthogonal")
check(matmul(p_a1, p_t2) == zero, "A1 and T2 projectors must be orthogonal")
check(matmul(p_e, p_t2) == zero, "E and T2 projectors must be orthogonal")
check(
    [[p_a1[i][j] + p_e[i][j] + p_t2[i][j] for j in range(6)] for i in range(6)] == i6,
    "A1+E+T2 projectors must resolve identity",
)


# 8. Frozen dependency custody.
deps = {
    "LANE_CROSS_RFT_GRA_GL6Y_FPSS_SHARED_CHILD_CONNECTOR_V001/THEOREM.md":
        "ee98b5eb31832a9e3c9bf7874c9e9d5530f3b6ef6c6e9023c7dbcb89b1153019",
    "LANE_CROSS_RFT_GRA_GL6Y_FPSS_SHARED_CHILD_CONNECTOR_V001/MANIFEST.sha256":
        "c1fc74d0c8b974d76420676f60c8584b3abcb4cfff5906ce210f2aa1e95b9d0e",
    "AUDIT_G_GL6Y_FPSS_SHARED_CHILD_CONNECTOR_V001/AUDIT.md":
        "b33b749b775b7ab5a3100f927845c2fef6b10c1718a9d65ffdfe542fb00d8983",
    "LANE_CROSS_RFT_GRA_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001/THEOREM.md":
        "faea49e3dcd5f2b4d5b3bab9026432d741192339dd789a939ef2318236848c0e",
    "LANE_CROSS_RFT_GRA_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001/MANIFEST.sha256":
        "9bebed96f2b864738639571f870ae7a34dc440dea5ab4418bc4e7cc9c2eb2a63",
    "AUDIT_G_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001/AUDIT.md":
        "cc34f4f8ea3824f6788209db4e3b9a1e03f034d38048b1098b221f5d741b0e0a",
    "LANE_CROSS_RFT_GRA_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001/THEOREM.md":
        "a51e802f6ba148e5f9848e95f41a80073795b24b7eaf87e36c0766b0856aa494",
    "LANE_CROSS_RFT_GRA_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001/MANIFEST.sha256":
        "fc50cad54dca00aab1c30d7c12ef07147df1242f94483f63955185695073f706",
    "AUDIT_G_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001/AUDIT.md":
        "5734ad57122c64e3174aa7706b0e7aa86102b3a18a3b868aca20af0997ab462a",
}
for rel, expected in deps.items():
    path = ROOT / rel
    check(path.is_file(), f"missing dependency {rel}")
    check(sha(path) == expected, f"dependency hash drift: {rel}")


# 9. Ledger and theorem ceiling checks.
ledger = json.loads((HERE / "BULK_DYNAMICS_LEDGER.json").read_text())
check(ledger["edges"]["site_degree"] == 6, "ledger site degree drift")
check("48 |Ud|/hbar" in ledger["interaction"]["lambda_F3"], "ledger lambda drift")
theorem = (HERE / "THEOREM.md").read_text()
for token in (
    "AUTHOR_FROZEN_POSTFREEZE_AUDIT_REQUIRED",
    "Not claimed",
    "physical momentum",
    "No artificial wraparound",
    "gapless, pole-dominated",
    "No finite-volume state limit is assumed",
    "no new coupling or free weight",
):
    check(token.lower() in theorem.lower(), f"missing theorem ceiling token: {token}")


print(f"PASS {CHECKS}/{CHECKS}")
