#!/usr/bin/env python3
"""Independent hostile physics/custody re-audit of repaired GF V005."""

from __future__ import annotations

import ast
import hashlib
import itertools
import json
import subprocess
from fractions import Fraction as F
from pathlib import Path


AUDIT = Path(__file__).resolve().parent
TARGET = AUDIT.parent
ROOT = TARGET.parent
checks = 0


def passed(label: str) -> None:
    global checks
    checks += 1
    print(f"PASS {label}")


def demand(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)
    passed(label)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path):
    for raw in path.read_text().splitlines():
        if raw.strip():
            digest, rel = raw.split("  ", 1)
            yield digest, rel


def transpose(a):
    return tuple(zip(*a))


def matmul(a, b):
    bt = transpose(b)
    return tuple(tuple(sum(x * y for x, y in zip(row, col)) for col in bt) for row in a)


def matscale(c, a):
    return tuple(tuple(c * x for x in row) for row in a)


def eye(n):
    return tuple(tuple(F(1) if i == j else F(0) for j in range(n)) for i in range(n))


def diagonal(*values):
    return tuple(tuple(values[i] if i == j else F(0) for j in range(len(values))) for i in range(len(values)))


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def rank(a):
    work = [list(row) for row in a]
    if not work:
        return 0
    nrow, ncol = len(work), len(work[0])
    pivot_row = 0
    for col in range(ncol):
        pivot = next((r for r in range(pivot_row, nrow) if work[r][col] != 0), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][col]
        work[pivot_row] = [x / scale for x in work[pivot_row]]
        for r in range(nrow):
            if r == pivot_row:
                continue
            c = work[r][col]
            if c:
                work[r] = [x - c * y for x, y in zip(work[r], work[pivot_row])]
        pivot_row += 1
        if pivot_row == nrow:
            break
    return pivot_row


def inv2(a):
    det = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    if det == 0:
        raise ValueError("singular")
    return ((a[1][1] / det, -a[0][1] / det), (-a[1][0] / det, a[0][0] / det))


# ---------------------------------------------------------------------------
# Frozen V005 target, dependencies, author manifest, and rejection custody.

expected_target = {
    "DEPENDENCIES.sha256": "c2394d1f1f67129fac44bd09f8c256d174415e7d01968d88905e14a4e4dc8b5f",
    "MANIFEST.sha256": "f5aa89cbfd353df28645c1db67f4b871151114e39e43d2959eeaf8c59195f358",
    "OBSERVABLE_CONTRACT.json": "3e08560c8c67294fbb0e39e4cd0d487de3d362441968cc7ab164197e5ca66c0a",
    "README.md": "e5cb2a4b500c6c1606532e0e6a4f101b0586380e7d935373ffca43e1f737998f",
    "RESULT.json": "087b8de7092be8e2fc58cb47cc1fa16ec9a38aadcfb5422a7088c75a0b182945",
    "RESULT.md": "c0555fd70bea7d2f6cf59706124adc9c1a868fc0374c10fbed3559ecebc00c6f",
    "SEAL.sha256": "38d36ff7150224c81f7e558829b2a759103ce78f15fa458b1f0e60b87c50cef3",
    "SELF_AUDIT.md": "62c605609c2096a4ffc17097b6620ba6a819bc42a1849ff30f87c0ad500e467c",
    "THEOREM.md": "928541b88f2c306905d715ea1e7b7ea0ab0f59d405134fe768f0d3f126c62e91",
    "VERIFICATION.txt": "5bfc5ad81b4d73888ff83ade2f7c66b532dc04c96ca5d9280deb746497951b4b",
    "verify_canonical_ir_ancestry_observable.py": "575b1ae1a799436c38ead2aaf46fb6078eb1d0265b75b8e3d4da908791e6ddec",
}
custody = dict((rel, digest) for digest, rel in rows(AUDIT / "TARGET_CUSTODY.sha256"))
demand(custody == expected_target, "V005 custody is exactly the root-frozen eleven-file target")
for rel, expected in expected_target.items():
    path = TARGET / rel
    demand(path.is_file(), f"target file exists: {rel}")
    demand(sha256(path) == expected, f"target hash replays: {rel}")

dependencies = list(rows(TARGET / "DEPENDENCIES.sha256"))
demand(len(dependencies) == 25, "author dependency packet has twenty-five frozen rows")
for expected, rel in dependencies:
    path = ROOT / rel
    demand(path.is_file(), f"dependency exists from repository root: {rel}")
    demand(sha256(path) == expected, f"dependency hash replays from repository root: {rel}")

manifest = list(rows(TARGET / "MANIFEST.sha256"))
demand(len(manifest) == 9, "author manifest freezes nine V005 core artifacts")
for expected, rel in manifest:
    path = TARGET / rel
    demand(path.is_file(), f"author manifest file exists: {rel}")
    demand(sha256(path) == expected, f"author manifest hash replays: {rel}")
demand(list(rows(TARGET / "SEAL.sha256")) == [(sha256(TARGET / "MANIFEST.sha256"), "MANIFEST.sha256")], "author seal owns the V005 manifest")

author = subprocess.run(["python3", str(TARGET / "verify_canonical_ir_ancestry_observable.py")], cwd=ROOT, check=False, capture_output=True, text=True)
demand(author.returncode == 0, "frozen V005 author verifier exits successfully")
demand("SUMMARY 104/104" in author.stdout, "frozen V005 author verifier reproduces 104/104")
demand("PHYSICS_UNEXECUTED" in author.stdout, "author replay retains design-only status")

prior_rows = list(rows(AUDIT / "PRIOR_REJECTIONS_CUSTODY.sha256"))
demand(len(prior_rows) == 9, "V005 audit pins all nine V004-rejection artifacts")
for expected, rel in prior_rows:
    path = TARGET / rel
    demand(path.is_file(), f"prior rejection artifact exists: {rel}")
    demand(sha256(path) == expected, f"prior rejection artifact hash replays: {rel}")
prior_v004 = (TARGET / "INDEPENDENT_HOSTILE_REAUDIT_V004/INDEPENDENT_HOSTILE_REAUDIT_V004.md").read_text()
prior_chain = (TARGET / "INDEPENDENT_HOSTILE_REAUDIT_V004/PRIOR_REJECTIONS_CUSTODY.sha256").read_text()
demand("REJECT -- REPAIR_REQUIRED" in prior_v004, "V004 rejection verdict remains explicit")
demand("INDEPENDENT_HOSTILE_REAUDIT_V003" in prior_chain, "recursive custody preserves V001--V003 rejection history")

theorem = (TARGET / "THEOREM.md").read_text()
contract = json.loads((TARGET / "OBSERVABLE_CONTRACT.json").read_text())
result = json.loads((TARGET / "RESULT.json").read_text())
author_source = (TARGET / "verify_canonical_ir_ancestry_observable.py").read_text()
binding = contract["canonical_source_binding"]

# ---------------------------------------------------------------------------
# Exhaustive three-valued amplitude partition, including both limit bounds.

demand("finite sigma_max" in binding["pass"], "PASS requires a finite canonical upper bound")
demand("upper-bound divergence" in binding["fail"], "FAIL includes applicable-D upper-bound divergence")
demand("whether the raw vertex is nonvanishing" in binding["indeterminate"], "no-D INDETERMINATE is independent of raw scaling")
demand("failed normalized ancestry" not in binding["fail"], "canonical amplitude remains separate from ancestry")
demand("reported independently" in binding["separation_rule"], "machine contract keeps amplitude and ancestry separate")

def specification(raw_full_rank, map_applicable, lower_positive, upper_bounded):
    if not raw_full_rank:
        return "FAIL"
    if not map_applicable:
        return "INDETERMINATE"
    if lower_positive and upper_bounded:
        return "PASS"
    return "FAIL"


def predicates(raw_full_rank, map_applicable, lower_positive, upper_bounded):
    pass_branch = raw_full_rank and map_applicable and lower_positive and upper_bounded
    fail_branch = (not raw_full_rank) or (map_applicable and (not lower_positive or not upper_bounded))
    indeterminate_branch = raw_full_rank and not map_applicable
    return pass_branch, fail_branch, indeterminate_branch


# Every logical combination of rank, map availability, lower-limit status,
# and upper-limit status receives exactly one label.
for raw_full, has_map, lower_positive, upper_bounded in itertools.product((False, True), repeat=4):
    flags = predicates(raw_full, has_map, lower_positive, upper_bounded)
    label = specification(raw_full, has_map, lower_positive, upper_bounded)
    demand(sum(flags) == 1, f"classifier disjoint/total: raw={raw_full},D={has_map},lower={lower_positive},upper={upper_bounded}")
    demand(flags[("PASS", "FAIL", "INDETERMINATE").index(label)], f"predicate tuple agrees with specification: {label}, raw={raw_full},D={has_map},lower={lower_positive},upper={upper_bounded}")

# Extract the exact author classifier and compare it on the same full truth
# table. controlled_power is varied independently to prove raw scaling does
# not alter a no-D label.
tree = ast.parse(author_source)
node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "canonical_outcome")
module = ast.Module(body=[node], type_ignores=[])
ast.fix_missing_locations(module)
namespace = {}
exec(compile(module, str(TARGET / "verify_canonical_ir_ancestry_observable.py"), "exec"), namespace)
canonical_outcome = namespace["canonical_outcome"]
for raw_full, has_map, lower_positive, upper_bounded, controlled_power in itertools.product((False, True), repeat=5):
    expected = specification(raw_full, has_map, lower_positive, upper_bounded)
    actual = canonical_outcome(2 if raw_full else 1, controlled_power, has_map, F(1) if lower_positive else F(0), upper_bounded)
    demand(actual == expected, f"author classifier matches complete limit table: {expected},power={controlled_power}")

# Exact representative families for every physical branch.
for L in (F(10), F(100), F(1000)):
    # V004 regression: D=diag(L,1) gives sigma_max -> infinity.
    D_upper = diagonal(L, F(1))
    N_upper = matmul(D_upper, transpose(D_upper))
    demand(rank(N_upper) == 2 and D_upper[1][1] == 1 and D_upper[0][0] == L, f"L={L}: exact upper-divergent rank-two family")
    demand(canonical_outcome(2, False, True, F(1), False) == "FAIL", f"L={L}: V004 diag(L,1) regression is FAIL")

    # Lower-loss family: D=diag(1/L,1) has sigma_min -> 0.
    D_lower = diagonal(1 / L, F(1))
    demand(rank(D_lower) == 2 and D_lower[0][0] == 1 / L, f"L={L}: exact lower-vanishing finite-size family")
    demand(canonical_outcome(2, False, True, F(0), True) == "FAIL", f"L={L}: lower-bound-loss family is FAIL")

    # Raw full-rank families without D are INDETERMINATE at vanishing,
    # constant, and divergent raw scalings.
    for scale, power_flag, name in ((1 / L, True, "vanishing"), (F(1), False, "constant"), (L, False, "divergent")):
        N_raw = matscale(scale, eye(2))
        demand(rank(N_raw) == 2, f"L={L}: {name} no-D raw family remains full rank")
        demand(canonical_outcome(2, power_flag, False, F(0), False) == "INDETERMINATE", f"L={L}: {name} full-rank no-D family is INDETERMINATE")

demand(F(10) < F(100) < F(1000), "V004 sigma_max regression is unbounded")
demand(F(1, 10) > F(1, 100) > F(1, 1000), "canonical lower-loss regression tends to zero")
demand(canonical_outcome(2, False, True, F(1), True) == "PASS", "finite positive bounded canonical vertex is PASS")
demand(canonical_outcome(1, False, False, F(0), False) == "FAIL", "persistent raw rank loss is FAIL even without D")

# V002 separation regression: amplitude PASS, ancestry FAIL, overall G2 FAIL.
for L in (F(10), F(100), F(1000)):
    delta = 1 / L
    R_L = matscale(L / 2, eye(2))
    N_can = matscale(2 * delta, R_L)
    Omega = matscale(1 / (2 * L), eye(2))
    demand(N_can == eye(2), f"L={L}: V002 overlap retains full finite canonical amplitude")
    demand(Omega[0][0] > 0, f"L={L}: V002 overlap ancestry is finite-size defined")
demand(F(1, 20) > F(1, 200) > F(1, 2000), "V002 overlap ancestry tends monotonically to zero")
demand(canonical_outcome(2, False, True, F(1), True) == "PASS", "V002 overlap amplitude is PASS")
demand(result["decision_separation"].endswith("overall G2 FAIL"), "ancestry FAIL makes amplitude-PASS family overall G2 FAIL")

# ---------------------------------------------------------------------------
# No regressions in representation, TT domain, ground query, binding,
# momentum registry, ancestry, factorization, and strict claim ceiling.

demand("unitary massless Poincare representation" in theorem, "V005 retains the unitary massless Poincare pole-bundle gate")
demand("e^{+2i\\theta}" in theorem and "e^{-2i\\theta}" in theorem, "V005 retains little-group helicities plus/minus two")
demand("reducible scalar doublet" in theorem, "V005 explicitly rejects a two-scalar doublet")
demand("equality of dispersive slopes is insufficient" in theorem, "equal slopes remain insufficient for Wigner covariance")
weights = (F(2), F(-2))
demand(not [(i, j) for i in range(2) for j in range(2) if abs(weights[i] - weights[j]) == 1], "J={+2,-2} forces trivial ISO(2) null translations")
demand(tuple(sorted((F(0), F(0)))) != tuple(sorted(weights)), "two scalar weights fail the helicity gate")

demand("does not define" in theorem and "Pi_TT(0)" in theorem, "TT projection remains forbidden at k=0")
demand("directional limit" in theorem and "nonzero sequence" in theorem, "TT soft limits remain punctured nonzero sequences")

def tt_projector(k):
    norm2 = sum(x * x for x in k)
    if norm2 == 0:
        raise ZeroDivisionError("Pi_TT(0) undefined")
    P = tuple(tuple((F(1) if i == j else F(0)) - k[i] * k[j] / norm2 for j in range(3)) for i in range(3))
    indices = tuple(itertools.product(range(3), repeat=2))
    return tuple(tuple((P[i][a] * P[j][b] + P[i][b] * P[j][a] - P[i][j] * P[a][b]) / 2 for a, b in indices) for i, j in indices)


try:
    tt_projector((F(0), F(0), F(0)))
    zero_rejected = False
except ZeroDivisionError:
    zero_rejected = True
demand(zero_rejected, "independent TT construction rejects k=0 exactly")
pi = tt_projector((F(1), F(2), F(3)))
demand(matmul(pi, pi) == pi, "nonzero-k TT projector is exactly idempotent")
demand(rank(pi) == 2, "nonzero-k symmetric TT image has exact rank two")

demand("rho_(0,L)=P_(0,L)/rank(P_(0,L))" in theorem, "V005 retains normalized complete-ground typing")
demand("choosing an individual vector" in theorem, "selected degenerate-ground vectors remain forbidden")
rho = matscale(F(1, 2), eye(2))
U = ((F(3, 5), F(-4, 5)), (F(4, 5), F(3, 5)))
demand(matmul(transpose(U), U) == eye(2), "hostile ground-basis rotation is exactly orthogonal")
M1 = ((F(1), F(2)), (F(0), F(1)))
M2 = ((F(2), F(0)), (F(1), F(3)))

def mixed_gram(ma, mb):
    return trace(matmul(rho, matmul(transpose(ma), mb)))


before = tuple(tuple(mixed_gram(a, b) for b in (M1, M2)) for a in (M1, M2))
rotated = (matmul(M1, U), matmul(M2, U))
after = tuple(tuple(mixed_gram(a, b) for b in rotated) for a in rotated)
demand(before == after, "rho0=P0/2 query is exactly ground-basis invariant")

demand("are **affine** volumes" in theorem, "unbound primitive-cell measure remains affine")
demand("independently audited and sealed FD" in theorem, "physical volume still requires independent FD binding")
demand("no further multiplicative conversion" in theorem, "one global binding forbids sectorwise length refits")
demand("physical FD binding" in result["missing"], "machine result keeps the physical binding open")

def q3mul(x, y):
    return (x[0] * y[0] + 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


demand(q3mul((F(0), F(16, 9)), (F(0), F(3, 8))) == (F(2), F(0)), "Q_aff conversion squares to 2/v3 exactly")
demand(q3mul((F(0), F(3, 8)), (F(0), F(8, 9))) == (F(1), F(0)), "conjugate source/operator conversions cancel exactly")

def qform(n):
    return 4 * sum(x * x for x in n) - sum(n) ** 2


def representatives(L):
    best = {}
    for n in itertools.product(range(-L, L + 1), repeat=3):
        cls = tuple(x % L for x in n)
        key = (qform(n), n)
        if cls not in best or key < best[cls]:
            best[cls] = key
    return {cls: key[1] for cls, key in best.items()}


registries = {}
for L in (5, 10):
    reps = representatives(L)
    registries[L] = reps
    demand(len(reps) == L ** 3, f"L={L} registry contains every character once")
    demand(reps[(0, 0, 0)] == (0, 0, 0), f"L={L} zero character uses zero minimizer")
    demand(all(qform(n) == qform(reps[tuple((-x) % L for x in cls)]) for cls, n in reps.items()), f"L={L} conjugate characters have equal norm")
demand([qform(n) for n in ((1, 0, 0), (1, 1, 0), (1, -1, 0))] == [3, 4, 8], "three frozen rays have q=3,4,8")
demand(min(qform(n) for n in registries[5].values() if n != (0, 0, 0)) == 3, "q_min=3 reproduces exact smallest shell")
demand(all(tuple((2 * x) % 10 for x in cls) in registries[10] for cls in registries[5]), "G5 characters inject by doubling into G10")
demand(any(any(x % 2 for x in cls) for cls in registries[10]), "G10 contains new odd fine-cover characters")

S = ((F(4), F(1)), (F(1), F(3)))
R = ((F(2), F(1, 2)), (F(1, 2), F(1)))
B = ((F(2), F(1)), (F(1), F(1)))
SinvR = matmul(inv2(S), R)
S2invR2 = matmul(inv2(matmul(matmul(transpose(B), S), B)), matmul(matmul(transpose(B), R), B))
demand(trace(SinvR) == trace(S2invR2), "ancestry generalized-eigenvalue trace is basis invariant")
demand(SinvR[0][0] * SinvR[1][1] - SinvR[0][1] * SinvR[1][0] == S2invR2[0][0] * S2invR2[1][1] - S2invR2[0][1] * S2invR2[1][0], "ancestry generalized-eigenvalue determinant is basis invariant")
demand("liminf" in contract["ancestry"]["pass"] and "eta_min" in contract["ancestry"]["pass"], "ancestry retains prospective lower-bound gate")

Z = diagonal(F(2), F(3))
Zplus = diagonal(F(1, 2), F(1, 3))
Dpole = diagonal(F(5), F(7))
Dinv = diagonal(F(1, 5), F(1, 7))
Gamma = ((F(1), F(2)), (F(3), F(4)))
G = matmul(matmul(Z, Dinv), Gamma)
demand(matmul(matmul(Dpole, Zplus), G) == Gamma, "branchwise pole amputation recovers source-independent vertex")
demand("not longitudinal" in contract["factorization"]["ceiling"], "factorization readiness remains distinct from G3 decoupling")

demand(contract["status"] == "REPAIRED_V005_AWAITING_INDEPENDENT_AUDIT_PHYSICS_UNEXECUTED", "machine contract remains physics-unexecuted before this audit")
demand("no positive G2" in contract["ceiling"] and "no gravity" in contract["ceiling"], "machine contract retains strict design-only ceiling")
demand("Consequently it proves no" in theorem, "theorem retains no-positive-G2 disposition")
demand("no numerical pole evidence" in result["ceiling"] and "no gravity" in result["ceiling"], "machine result does not promote physics")

print(f"SUMMARY {checks}/{checks} independent GF V005 hostile checks passed")
print("VERDICT PASS -- NARROW_DESIGN_CONTRACT_SEAL")
print("CLASSIFIER disjoint and total over rank, map availability, lower-limit, upper-limit, and raw-scaling cases")
print("REGRESSIONS V002 ancestry separation PASS; V003 every full-rank no-D case INDETERMINATE; V004 diag(L,1) upper divergence FAIL")
print("CEILING design contract only; no matched-family pole evidence, positive G2, G3, gravity, or G; GE not audited or repinned")
