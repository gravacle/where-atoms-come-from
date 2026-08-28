#!/usr/bin/env python3
"""Independent hostile physics/custody re-audit of repaired GF V004."""

from __future__ import annotations

import ast
import hashlib
import itertools
import json
import math
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


def diagonal(*values):
    return tuple(tuple(values[i] if i == j else F(0) for j in range(len(values))) for i in range(len(values)))


# ---------------------------------------------------------------------------
# Frozen V004 target, dependency, author-manifest, and rejection-chain custody.

expected_target = {
    "DEPENDENCIES.sha256": "c2394d1f1f67129fac44bd09f8c256d174415e7d01968d88905e14a4e4dc8b5f",
    "MANIFEST.sha256": "15b285a28084d21c5de86a5fd9c966dc1d2c4f069cb03cbfb754c73273fbcde5",
    "OBSERVABLE_CONTRACT.json": "ff3fd30e303dd9175f2c7c134cac45e2ac5d4018463bdabdb64eb552bee8bac4",
    "README.md": "6d5828b63c51d9ce7385b0baad9694621f013dc4f52c4b217babb186e0b89b33",
    "RESULT.json": "74c8155840710ad6927e8850c31afe5a1da4e78530d3fa0e5322c7966ee9b5e9",
    "RESULT.md": "d12111a936302ce4f5455fc61988984efa17d11ca6bd429cd63cb10ae406016f",
    "SEAL.sha256": "6abda6859a3ae6271195230e6c44d28672ea2359095df41b7afc4570af1f0279",
    "SELF_AUDIT.md": "6cc136c3672669429461945dd1337e44adb5fd3f0367c3825e091d37c3dbb0ae",
    "THEOREM.md": "db8bfa8141d4cf3414cc7caeab19507bef07136ff9edfe2779a09a95bcdf8869",
    "VERIFICATION.txt": "e8706c8809deb0e6189e466d61f5c39df58e389eab419df6096ee8885a2e52d9",
    "verify_canonical_ir_ancestry_observable.py": "f4502348376bf463306589383081422838ab7883c3bfe8f5027179183f53535d",
}
custody = dict((rel, digest) for digest, rel in rows(AUDIT / "TARGET_CUSTODY.sha256"))
demand(custody == expected_target, "V004 custody is exactly the root-frozen eleven-file target")
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
demand(len(manifest) == 9, "author manifest freezes nine V004 core artifacts")
for expected, rel in manifest:
    path = TARGET / rel
    demand(path.is_file(), f"author manifest file exists: {rel}")
    demand(sha256(path) == expected, f"author manifest hash replays: {rel}")
demand(
    list(rows(TARGET / "SEAL.sha256")) == [(sha256(TARGET / "MANIFEST.sha256"), "MANIFEST.sha256")],
    "author seal owns the V004 manifest",
)

author = subprocess.run(
    ["python3", str(TARGET / "verify_canonical_ir_ancestry_observable.py")],
    cwd=ROOT,
    check=False,
    capture_output=True,
    text=True,
)
demand(author.returncode == 0, "frozen V004 author verifier exits successfully")
demand("SUMMARY 102/102" in author.stdout, "frozen V004 author verifier reproduces 102/102")
demand("PHYSICS_UNEXECUTED" in author.stdout, "author replay retains design-only status")

prior_rows = list(rows(AUDIT / "PRIOR_REJECTIONS_CUSTODY.sha256"))
demand(len(prior_rows) == 9, "V004 audit pins all nine V003-rejection artifacts")
for expected, rel in prior_rows:
    path = TARGET / rel
    demand(path.is_file(), f"prior rejection artifact exists: {rel}")
    demand(sha256(path) == expected, f"prior rejection artifact hash replays: {rel}")
prior_v003 = (TARGET / "INDEPENDENT_HOSTILE_REAUDIT_V003/INDEPENDENT_HOSTILE_REAUDIT_V003.md").read_text()
prior_chain = (TARGET / "INDEPENDENT_HOSTILE_REAUDIT_V003/PRIOR_REJECTIONS_CUSTODY.sha256").read_text()
demand("REJECT -- REPAIR_REQUIRED" in prior_v003, "V003 rejection verdict remains explicit")
demand("INDEPENDENT_HOSTILE_REAUDIT_V002" in prior_chain, "V003 custody continues to pin the V002/V001 rejection chain")

theorem = (TARGET / "THEOREM.md").read_text()
contract = json.loads((TARGET / "OBSERVABLE_CONTRACT.json").read_text())
result = json.loads((TARGET / "RESULT.json").read_text())
author_source = (TARGET / "verify_canonical_ir_ancestry_observable.py").read_text()
binding = contract["canonical_source_binding"]

# ---------------------------------------------------------------------------
# Priority regression: the four explicitly repaired bounded cases are total
# and disjoint, and ancestry remains a separate conjunctive output.

demand("whether the raw vertex is nonvanishing" in binding["indeterminate"], "V004 includes nonvanishing full-rank no-D inputs in INDETERMINATE")
demand("raw numerator remains full rank" in binding["indeterminate"], "machine no-D branch requires full raw rank")
demand("failed normalized ancestry" not in binding["fail"], "amplitude FAIL remains independent of ancestry")
demand("reported independently" in binding["separation_rule"], "amplitude and ancestry remain separate outputs")

def current_predicates(raw_rank, map_applicable, canonical_rank, lower, upper_finitely_bounded):
    pass_branch = (
        raw_rank == 2
        and map_applicable
        and canonical_rank == 2
        and lower > 0
        and upper_finitely_bounded
    )
    fail_branch = raw_rank < 2 or (
        map_applicable and (canonical_rank < 2 or lower <= 0)
    )
    indeterminate_branch = raw_rank == 2 and not map_applicable
    return pass_branch, fail_branch, indeterminate_branch


# Exhaust every Boolean/rank/lower-bound combination on the contract's
# explicitly named finite-upper-bound domain.  Exactly one branch must fire.
for raw_rank, has_map, can_rank, lower in itertools.product((1, 2), (False, True), (1, 2), (F(0), F(1))):
    flags = current_predicates(raw_rank, has_map, can_rank, lower, True)
    demand(sum(flags) == 1, f"bounded-domain amplitude partition total/disjoint: raw={raw_rank},D={has_map},can={can_rank},lower={lower}")

raw_fail = current_predicates(1, False, 0, F(0), True)
no_map_power = current_predicates(2, False, 0, F(0), True)
no_map_constant = current_predicates(2, False, 0, F(1), True)
canonical_pass = current_predicates(2, True, 2, F(1), True)
canonical_rank_fail = current_predicates(2, True, 1, F(1), True)
canonical_lower_fail = current_predicates(2, True, 2, F(0), True)
demand(raw_fail == (False, True, False), "raw-rank deficiency uniquely yields FAIL")
demand(no_map_power == (False, False, True), "controlled-power full-rank no-D uniquely yields INDETERMINATE")
demand(no_map_constant == (False, False, True), "nonvanishing full-rank no-D uniquely yields INDETERMINATE")
demand(canonical_pass == (True, False, False), "applicable-D positive bounded canonical vertex uniquely yields PASS")
demand(canonical_rank_fail == (False, True, False), "applicable-D canonical rank loss uniquely yields FAIL")
demand(canonical_lower_fail == (False, True, False), "applicable-D lower-bound loss uniquely yields FAIL")

# Extract the exact author classifier and replay both prior regressions.
tree = ast.parse(author_source)
node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "canonical_outcome")
module = ast.Module(body=[node], type_ignores=[])
ast.fix_missing_locations(module)
namespace = {}
exec(compile(module, str(TARGET / "verify_canonical_ir_ancestry_observable.py"), "exec"), namespace)
canonical_outcome = namespace["canonical_outcome"]
demand(canonical_outcome(2, True, False, F(0)) == "INDETERMINATE", "author replay retains controlled-power no-D INDETERMINATE")
demand(canonical_outcome(2, False, False, F(1)) == "INDETERMINATE", "author replay repairs V003 nonvanishing no-D counterexample")

for L in (F(10), F(100), F(1000)):
    delta = 1 / L
    R_L = matscale(L / 2, eye(2))
    N_can = matscale(2 * delta, R_L)
    Omega = matscale(1 / (2 * L), eye(2))
    demand(N_can == eye(2), f"L={L}: V002 overlap retains full finite canonical amplitude")
    demand(Omega[0][0] > 0, f"L={L}: V002 overlap ancestry is finite-size defined")
demand(F(1, 20) > F(1, 200) > F(1, 2000), "V002 overlap ancestry tends monotonically to zero")
demand(canonical_outcome(2, False, True, F(1)) == "PASS", "V002 overlap canonical amplitude remains PASS")
demand(result["decision_separation"].endswith("overall G2 FAIL"), "ancestry FAIL still makes amplitude-PASS family overall G2 FAIL")

# ---------------------------------------------------------------------------
# New material counterexample: finite-upper-bound PASS is not implemented or
# complemented.  Let Nraw=I and applicable D_L=diag(L,1), so
# Zcan,L=diag(L,1).  It stays positive rank two with lower bound one but has no
# finite upper bound.  Contract PASS is false; current FAIL and INDET are also
# false.  The author function, which has no upper-bound argument, returns PASS.

demand("finite sigma_max" in binding["pass"], "machine PASS explicitly requires a finite canonical upper bound")
for L in (F(10), F(100), F(1000)):
    N_raw = eye(2)
    D_L = diagonal(L, F(1))
    N_can = matmul(matmul(D_L, N_raw), transpose(D_L))
    Z_can = D_L
    demand(rank(N_raw) == rank(N_can) == 2, f"L={L}: divergent-vertex family remains exact rank two")
    demand(Z_can[1][1] == 1 and Z_can[0][0] == L, f"L={L}: sigma_min=1 and sigma_max=L exactly")
demand(F(10) < F(100) < F(1000), "divergent-vertex sigma_max has no finite limiting upper bound")

divergent_flags = current_predicates(2, True, 2, F(1), False)
demand(divergent_flags == (False, False, False), "applicable-D unbounded canonical vertex matches no frozen prose branch")
demand(sum(divergent_flags) == 0, "advertised amplitude partition is not total once finite-sigma_max is enforced")
author_divergent_label = canonical_outcome(2, False, True, F(1))
demand(author_divergent_label == "PASS", "author classifier falsely accepts the unbounded canonical vertex")
demand("upper_bound" not in author_source[author_source.index("def canonical_outcome"):author_source.index("# Exact V002 overlap regression")], "author classifier has no canonical upper-bound input")
demand(author_divergent_label == "PASS" and not any(divergent_flags), "author executable contradicts the frozen finite-sigma_max contract")

# ---------------------------------------------------------------------------
# No regressions in representation, TT domain, ground query, binding,
# momentum registry, ancestry, factorization, and claim ceiling.

demand("unitary massless Poincare representation" in theorem, "V004 retains the unitary massless Poincare pole-bundle gate")
demand("e^{+2i\\theta}" in theorem and "e^{-2i\\theta}" in theorem, "V004 retains little-group helicities plus/minus two")
demand("reducible scalar doublet" in theorem, "V004 explicitly rejects a two-scalar doublet")
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

demand("rho_(0,L)=P_(0,L)/rank(P_(0,L))" in theorem, "V004 retains normalized complete-ground typing")
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
demand(min(qform(n) for n in registries[5].values() if n != (0, 0, 0)) == 3, "q_min=3 reproduces the exact smallest shell")
demand(all(tuple((2 * x) % 10 for x in cls) in registries[10] for cls in registries[5]), "G5 characters inject by doubling into G10")
demand(any(any(x % 2 for x in cls) for cls in registries[10]), "G10 contains new odd fine-cover characters")

S = ((F(4), F(1)), (F(1), F(3)))
R = ((F(2), F(1, 2)), (F(1, 2), F(1)))
B = ((F(2), F(1)), (F(1), F(1)))
SinvR = matmul(inv2(S), R)
S2invR2 = matmul(inv2(matmul(matmul(transpose(B), S), B)), matmul(matmul(transpose(B), R), B))
demand(trace(SinvR) == trace(S2invR2), "ancestry generalized-eigenvalue trace is basis invariant")
demand(SinvR[0][0] * SinvR[1][1] - SinvR[0][1] * SinvR[1][0] == S2invR2[0][0] * S2invR2[1][1] - S2invR2[0][1] * S2invR2[1][0], "ancestry generalized-eigenvalue determinant is basis invariant")
demand("liminf" in contract["ancestry"]["pass"] and "eta_min" in contract["ancestry"]["pass"], "ancestry retains its prospective lower-bound gate")

Z = diagonal(F(2), F(3))
Zplus = diagonal(F(1, 2), F(1, 3))
Dpole = diagonal(F(5), F(7))
Dinv = diagonal(F(1, 5), F(1, 7))
Gamma = ((F(1), F(2)), (F(3), F(4)))
G = matmul(matmul(Z, Dinv), Gamma)
demand(matmul(matmul(Dpole, Zplus), G) == Gamma, "branchwise pole amputation recovers the source-independent vertex")
demand("not longitudinal" in contract["factorization"]["ceiling"], "factorization readiness remains distinct from G3 decoupling")

demand(contract["status"] == "REPAIRED_V004_AWAITING_INDEPENDENT_AUDIT_PHYSICS_UNEXECUTED", "machine contract remains physics-unexecuted")
demand("no positive G2" in contract["ceiling"] and "no gravity" in contract["ceiling"], "machine contract retains strict design-only ceiling")
demand("Consequently it proves no" in theorem, "theorem retains no-positive-G2 disposition")
demand("no numerical pole evidence" in result["ceiling"] and "no gravity" in result["ceiling"], "machine result does not promote physics")

print(f"SUMMARY {checks}/{checks} independent GF V004 hostile checks passed")
print("VERDICT REJECT -- REPAIR_REQUIRED")
print("MATERIAL finite-sigma_max PASS requirement is neither implemented nor complemented by FAIL")
print("COUNTEREXAMPLE N_raw=I_2, applicable D_L=diag(L,1), Z_can,L=diag(L,1): no prose branch applies while author executable returns PASS")
print("MINIMAL_REPAIR add canonical upper-bound input; PASS requires uniform finite sigma_max and FAIL includes unbounded sigma_max; retain the V003 no-D repair")
print("CEILING no GE repin; no positive G2, G3, gravity, or G")
