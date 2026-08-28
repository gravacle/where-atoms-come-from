#!/usr/bin/env python3
"""Independent hostile physics/custody re-audit of repaired GF V003."""

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


def diag(*values):
    return tuple(tuple(values[i] if i == j else F(0) for j in range(len(values))) for i in range(len(values)))


# ---------------------------------------------------------------------------
# Frozen V003 target, dependency, author-manifest, and prior-rejection custody.

expected_target = {
    "DEPENDENCIES.sha256": "c2394d1f1f67129fac44bd09f8c256d174415e7d01968d88905e14a4e4dc8b5f",
    "MANIFEST.sha256": "2a2fa56eb2c4d20cb1f31f3f80af7342485ca48ecc62bb310f7c9bf053edd91b",
    "OBSERVABLE_CONTRACT.json": "5d09c0b68ec2089b77085b9edf70b39cf4fe8ab285c8063c800cef1c60d5ee39",
    "README.md": "20210d378ba03253386f7e6ba698453640648d2e63c71d68eae7a1d758b18905",
    "RESULT.json": "ce2245a321ec4bcff0d46ee680899e874b75668d9160b7f6bd8431c83f08bc21",
    "RESULT.md": "cc5e07b9ee29a03e17c36d87e536c044da9e25e0555d85cdc6a9595b007e9c79",
    "SEAL.sha256": "64deb94c9d2539bdc600d66127c2e096c94c965cc035ea3c176cf921f4e7767b",
    "SELF_AUDIT.md": "a16190136451db25f63aea7a145212c55230401095a0d4e5b981d046964068b0",
    "THEOREM.md": "8bd228871e7c64eadfd009e0abb9f7dbfd95ac7e4fdd67d343f69a4ea18002e1",
    "VERIFICATION.txt": "c61739a6c89070ff67e414a5aa4f3a68b542abd32f611a051220c60dccd3bdd6",
    "verify_canonical_ir_ancestry_observable.py": "7c4b5877e81c1dbda9612ef954b43211901711d1e6e57a5586d17a6a11664b69",
}
custody = dict((rel, digest) for digest, rel in rows(AUDIT / "TARGET_CUSTODY.sha256"))
demand(custody == expected_target, "V003 custody is exactly the root-frozen eleven-file target")
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
demand(len(manifest) == 9, "author manifest freezes nine V003 core artifacts")
for expected, rel in manifest:
    path = TARGET / rel
    demand(path.is_file(), f"author manifest file exists: {rel}")
    demand(sha256(path) == expected, f"author manifest hash replays: {rel}")

seal = list(rows(TARGET / "SEAL.sha256"))
demand(seal == [(sha256(TARGET / "MANIFEST.sha256"), "MANIFEST.sha256")], "author seal owns the V003 manifest")

author = subprocess.run(
    ["python3", str(TARGET / "verify_canonical_ir_ancestry_observable.py")],
    cwd=ROOT,
    check=False,
    capture_output=True,
    text=True,
)
demand(author.returncode == 0, "frozen V003 author verifier exits successfully")
demand("SUMMARY 101/101" in author.stdout, "frozen V003 author verifier reproduces 101/101")
demand("PHYSICS_UNEXECUTED" in author.stdout, "author replay retains the design-only status")

prior_rows = list(rows(AUDIT / "PRIOR_REJECTIONS_CUSTODY.sha256"))
demand(len(prior_rows) == 9, "V003 audit pins all nine prior V002-rejection artifacts")
for expected, rel in prior_rows:
    path = TARGET / rel
    demand(path.is_file(), f"prior rejection artifact exists: {rel}")
    demand(sha256(path) == expected, f"prior rejection artifact hash replays: {rel}")
prior_v002 = (TARGET / "INDEPENDENT_HOSTILE_REAUDIT_V002/INDEPENDENT_HOSTILE_REAUDIT_V002.md").read_text()
prior_v001 = (TARGET / "INDEPENDENT_HOSTILE_REAUDIT_V002/PRIOR_V001_REJECTION.md").read_text()
demand("REJECT -- REPAIR_REQUIRED" in prior_v002, "V002 rejection verdict remains explicit")
demand("V001 hostile audit rejected" in prior_v001, "V001 rejection remains preserved inside V002 custody")

theorem = (TARGET / "THEOREM.md").read_text()
contract = json.loads((TARGET / "OBSERVABLE_CONTRACT.json").read_text())
result = json.loads((TARGET / "RESULT.json").read_text())
author_source = (TARGET / "verify_canonical_ir_ancestry_observable.py").read_text()

# ---------------------------------------------------------------------------
# Required regression: V002 amplitude/ancestry overlap is genuinely repaired.

binding = contract["canonical_source_binding"]
demand("failed normalized ancestry" not in binding["fail"], "V003 removes ancestry failure from canonical-amplitude FAIL")
demand("reported independently" in binding["separation_rule"], "contract keeps amplitude and ancestry as separate outputs")
demand("canonical amplitude only" in theorem, "theorem explicitly types the three-valued label as amplitude only")
demand("amplitude PASS and ancestry FAIL" in theorem, "theorem states the separated V002 regression outcome")

# Exact V002 regression: Delta=1/L, R=(L/2)I gives Ncan=I, while
# S=L^2 I gives Omega=(1/(2L))I -> 0.  The amplitude label is PASS; ancestry
# and hence the outer positive-G2 conjunction fail.
for L in (F(10), F(100), F(1000)):
    delta = 1 / L
    R_L = matscale(L / 2, eye(2))
    N_can = matscale(2 * delta, R_L)
    Omega = matscale(1 / (2 * L), eye(2))
    demand(N_can == eye(2), f"L={L}: V002 regression retains full finite canonical amplitude")
    demand(Omega[0][0] > 0, f"L={L}: V002 regression ancestry is finite-size well-defined")
demand(F(1, 20) > F(1, 200) > F(1, 2000), "V002 regression ancestry tends monotonically to zero")
amplitude_pass = True
amplitude_fail = False
ancestry_pass = False
overall_g2 = amplitude_pass and ancestry_pass
demand(amplitude_pass and not amplitude_fail, "V002 regression has one canonical-amplitude label: PASS")
demand(not ancestry_pass and not overall_g2, "separate ancestry gate makes the V002 regression overall G2 FAIL")
demand(result["decision_separation"].endswith("overall G2 FAIL"), "machine result preserves the corrected outer conjunction")

# ---------------------------------------------------------------------------
# New hostile regression: the advertised amplitude partition is not total.

def frozen_predicates(raw_rank, raw_vanishes_with_controlled_power, map_applicable, canonical_rank, lower, upper):
    pass_branch = (
        map_applicable
        and canonical_rank == 2
        and lower > 0
        and math.isfinite(float(lower))
        and math.isfinite(float(upper))
    )
    fail_branch = raw_rank < 2 or (
        map_applicable and (canonical_rank < 2 or lower <= 0)
    )
    indeterminate_branch = (
        raw_rank == 2
        and raw_vanishes_with_controlled_power
        and not map_applicable
    )
    return pass_branch, fail_branch, indeterminate_branch


ordinary_cases = (
    (1, False, False, 0, F(0), F(0), "raw-rank FAIL"),
    (2, True, False, 0, F(0), F(0), "controlled-power INDETERMINATE"),
    (2, False, True, 2, F(1, 3), F(2), "canonical PASS"),
    (2, False, True, 1, F(0), F(2), "canonical-rank FAIL"),
    (2, False, True, 2, F(0), F(2), "canonical-lower-bound FAIL"),
)
for raw_rank, vanishes, has_map, can_rank, lower, upper, label in ordinary_cases:
    flags = frozen_predicates(raw_rank, vanishes, has_map, can_rank, lower, upper)
    demand(sum(flags) == 1, f"frozen amplitude branches are disjoint on {label}")

# Exact uncovered family: N_raw,L = I_2 at every nonzero IR momentum.  It is
# full rank and does not vanish.  No same-parent D has been derived.  Because
# the contract itself says D=I is raw reporting only, the canonical amplitude
# is unknown.  None of the three frozen prose predicates applies.
N_raw = eye(2)
demand(rank(N_raw) == 2, "uncovered raw numerator family has exact full rank two")
demand(N_raw == eye(2), "uncovered raw numerator is finite and nonvanishing on every size")
demand(binding["identity_rule"].startswith("D=I is the raw reporting convention"), "D=I cannot silently convert raw reporting into a canonical claim")
hole_flags = frozen_predicates(2, False, False, 0, F(0), F(0))
demand(hole_flags == (False, False, False), "full-rank nonvanishing raw source without D matches no frozen amplitude branch")
demand(sum(hole_flags) == 0, "frozen PASS/FAIL/INDETERMINATE amplitude partition is not total")

# Extract and execute the author's exact classifier without importing/running
# the rest of the author module.  Its fallback assigns the uncovered input to
# FAIL, contradicting the theorem/contract rather than closing the contract.
tree = ast.parse(author_source)
node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "canonical_outcome")
module = ast.Module(body=[node], type_ignores=[])
ast.fix_missing_locations(module)
namespace = {}
exec(compile(module, str(TARGET / "verify_canonical_ir_ancestry_observable.py"), "exec"), namespace)
author_label = namespace["canonical_outcome"](2, False, False, F(1))
demand(author_label == "FAIL", "author executable silently labels the uncovered no-D constant family FAIL")
demand("vanishes with a controlled" in binding["indeterminate"], "contract INDETERMINATE prose excludes the nonvanishing no-D family")
demand("vanishes with a controlled" in theorem, "theorem INDETERMINATE prose excludes the nonvanishing no-D family")
demand(author_label == "FAIL" and sum(hole_flags) == 0, "author executable and frozen decision prose disagree on the same exact input")

# The outer G2 conjunction cannot repair a missing/contradictory component
# label: it can prevent positive promotion, but a prospective observable
# contract still must distinguish physical FAIL from unavailable canonical
# normalization.  The minimal total extension is every full-rank/no-D input
# -> INDETERMINATE, with raw scaling reported unchanged.
proposed_label = "INDETERMINATE" if rank(N_raw) == 2 and not False else "FAIL"
demand(proposed_label == "INDETERMINATE", "minimal total repair assigns every full-rank no-D family INDETERMINATE")
demand(not overall_g2, "outer G2 conjunction remains conservative while the classifier is repaired")

# ---------------------------------------------------------------------------
# No regressions in representation, TT domain, ground query, binding,
# momentum registry, ancestry, factorization, and claim ceiling.

demand("unitary massless Poincare representation" in theorem, "V003 retains the unitary massless Poincare pole-bundle gate")
demand("e^{+2i\\theta}" in theorem and "e^{-2i\\theta}" in theorem, "V003 retains little-group helicities plus/minus two")
demand("reducible scalar doublet" in theorem, "V003 explicitly rejects a two-scalar doublet")
demand("equality of dispersive slopes is insufficient" in theorem, "equal slopes remain insufficient for Wigner covariance")

weights = (F(2), F(-2))
allowed_plus = [(i, j) for i in range(2) for j in range(2) if weights[i] - weights[j] == 1]
allowed_minus = [(i, j) for i in range(2) for j in range(2) if weights[i] - weights[j] == -1]
demand(allowed_plus == [] and allowed_minus == [], "J={+2,-2} forces trivial ISO(2) null translations on the two-state fiber")
demand(tuple(sorted(weights)) == (F(-2), F(2)), "pole fiber has exactly the required Wigner weights")
demand(tuple(sorted((F(0), F(0)))) != (F(-2), F(2)), "two scalar weights fail the helicity gate")

demand("does not define" in theorem and "Pi_TT(0)" in theorem, "TT projection remains forbidden at k=0")
demand("directional limit" in theorem and "nonzero sequence" in theorem, "TT soft limits remain punctured nonzero sequences")

def tt_projector(k):
    norm2 = sum(x * x for x in k)
    if norm2 == 0:
        raise ZeroDivisionError("Pi_TT(0) undefined")
    P = tuple(tuple((F(1) if i == j else F(0)) - k[i] * k[j] / norm2 for j in range(3)) for i in range(3))
    indices = tuple(itertools.product(range(3), repeat=2))
    return tuple(
        tuple((P[i][a] * P[j][b] + P[i][b] * P[j][a] - P[i][j] * P[a][b]) / 2 for a, b in indices)
        for i, j in indices
    )


try:
    tt_projector((F(0), F(0), F(0)))
    zero_rejected = False
except ZeroDivisionError:
    zero_rejected = True
demand(zero_rejected, "independent TT construction rejects k=0 exactly")
pi = tt_projector((F(1), F(2), F(3)))
demand(matmul(pi, pi) == pi, "nonzero-k TT projector is exactly idempotent")
demand(rank(pi) == 2, "nonzero-k symmetric TT image has exact rank two")

demand("rho_(0,L)=P_(0,L)/rank(P_(0,L))" in theorem, "V003 retains normalized complete-ground typing")
demand("choosing an individual vector" in theorem, "selected vectors in a degenerate ground space remain forbidden")
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
e1 = ((F(1),), (F(0),))
Ue1 = matmul(U, e1)
selected_before = sum(x[0] * x[0] for x in matmul(M1, e1))
selected_after = sum(x[0] * x[0] for x in matmul(M1, Ue1))
demand(selected_before != selected_after, "selected degenerate-ground vector is basis dependent")

demand("are **affine** volumes" in theorem, "unbound primitive-cell measure remains affine")
demand("independently audited and sealed FD" in theorem, "physical volume still requires independent FD binding")
demand("no further multiplicative conversion" in theorem, "one global binding forbids sectorwise length refits")
demand("physical FD binding" in result["missing"], "machine result keeps the physical binding open")

def q3mul(x, y):
    return (x[0] * y[0] + 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


v3 = (F(0), F(16, 9))
two_over_v3 = (F(0), F(3, 8))
v3_over_two = (F(0), F(8, 9))
demand(q3mul(v3, two_over_v3) == (F(2), F(0)), "Q_aff conversion squares to 2/v3 exactly")
demand(q3mul(two_over_v3, v3_over_two) == (F(1), F(0)), "conjugate source/operator conversions cancel exactly")

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
    demand(reps[(0, 0, 0)] == (0, 0, 0), f"L={L} zero character uses the zero minimizer")
    demand(all(qform(n) == qform(reps[tuple((-x) % L for x in cls)]) for cls, n in reps.items()), f"L={L} conjugate characters have equal norm")
demand([qform(n) for n in ((1, 0, 0), (1, 1, 0), (1, -1, 0))] == [3, 4, 8], "three frozen rays have q=3,4,8")
demand(min(qform(n) for n in registries[5].values() if n != (0, 0, 0)) == 3, "q_min=3 reproduces the exact smallest shell")
demand(all(tuple((2 * x) % 10 for x in cls) in registries[10] for cls in registries[5]), "G5 characters inject by doubling into G10")
demand(any(any(x % 2 for x in cls) for cls in registries[10]), "G10 contains new odd fine-cover characters")

S = ((F(4), F(1)), (F(1), F(3)))
R = ((F(2), F(1, 2)), (F(1, 2), F(1)))
B = ((F(2), F(1)), (F(1), F(1)))
SinvR = matmul(inv2(S), R)
S2 = matmul(matmul(transpose(B), S), B)
R2 = matmul(matmul(transpose(B), R), B)
S2invR2 = matmul(inv2(S2), R2)
demand(trace(SinvR) == trace(S2invR2), "ancestry generalized-eigenvalue trace is basis invariant")
demand(
    SinvR[0][0] * SinvR[1][1] - SinvR[0][1] * SinvR[1][0]
    == S2invR2[0][0] * S2invR2[1][1] - S2invR2[0][1] * S2invR2[1][0],
    "ancestry generalized-eigenvalue determinant is basis invariant",
)
demand("liminf" in contract["ancestry"]["pass"] and "eta_min" in contract["ancestry"]["pass"], "ancestry retains its prospective positive lower-bound gate")

Z = diag(F(2), F(3))
Zplus = diag(F(1, 2), F(1, 3))
Dpole = diag(F(5), F(7))
Dinv = diag(F(1, 5), F(1, 7))
Gamma = ((F(1), F(2)), (F(3), F(4)))
G = matmul(matmul(Z, Dinv), Gamma)
demand(matmul(matmul(Dpole, Zplus), G) == Gamma, "branchwise pole amputation recovers the source-independent vertex")
demand("not longitudinal" in contract["factorization"]["ceiling"], "factorization readiness remains distinct from G3 decoupling")

demand(contract["status"] == "REPAIRED_V003_AWAITING_INDEPENDENT_AUDIT_PHYSICS_UNEXECUTED", "machine contract remains physics-unexecuted")
demand("no positive G2" in contract["ceiling"] and "no gravity" in contract["ceiling"], "machine contract retains the strict design-only ceiling")
demand("Consequently it proves no" in theorem, "theorem retains the no-positive-G2 disposition")
demand("no numerical pole evidence" in result["ceiling"] and "no gravity" in result["ceiling"], "machine result does not promote physics")

print(f"SUMMARY {checks}/{checks} independent GF V003 hostile checks passed")
print("VERDICT REJECT -- REPAIR_REQUIRED")
print("MATERIAL canonical-amplitude PASS/FAIL/INDETERMINATE is disjoint but not total")
print("COUNTEREXAMPLE N_raw,L=I_2, full rank and nonvanishing, with no applicable frozen D: no prose branch applies while the author executable returns FAIL")
print("MINIMAL_REPAIR classify every full-rank raw source with no applicable same-parent D as INDETERMINATE and report its raw scaling; keep ancestry separate")
print("CEILING no GE repin; no positive G2, G3, gravity, or G")
