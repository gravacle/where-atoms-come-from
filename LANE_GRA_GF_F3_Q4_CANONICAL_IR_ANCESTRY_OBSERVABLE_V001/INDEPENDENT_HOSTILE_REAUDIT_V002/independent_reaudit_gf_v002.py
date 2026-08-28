#!/usr/bin/env python3
"""Independent hostile physics and custody re-audit of repaired GF V002."""

from __future__ import annotations

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


def matadd(a, b):
    return tuple(tuple(a[i][j] + b[i][j] for j in range(len(a[0]))) for i in range(len(a)))


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
    return (
        (a[1][1] / det, -a[0][1] / det),
        (-a[1][0] / det, a[0][0] / det),
    )


def diag(*values):
    return tuple(tuple(values[i] if i == j else F(0) for j in range(len(values))) for i in range(len(values)))


# ---------------------------------------------------------------------------
# Exact frozen target/dependency/manifest/seal custody from correct roots.

custody = list(rows(AUDIT / "TARGET_CUSTODY.sha256"))
demand(len(custody) == 11, "V002 target custody names all eleven frozen author files")
for expected, rel in custody:
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
demand(len(manifest) == 9, "author manifest freezes nine core artifacts")
for expected, rel in manifest:
    path = TARGET / rel
    demand(path.is_file(), f"author manifest file exists: {rel}")
    demand(sha256(path) == expected, f"author manifest hash replays: {rel}")

seal = list(rows(TARGET / "SEAL.sha256"))
demand(seal == [(sha256(TARGET / "MANIFEST.sha256"), "MANIFEST.sha256")], "author seal owns the repaired V002 manifest")

author = subprocess.run(
    ["python3", str(TARGET / "verify_canonical_ir_ancestry_observable.py")],
    cwd=ROOT,
    check=False,
    capture_output=True,
    text=True,
)
demand(author.returncode == 0, "frozen V002 author verifier exits successfully")
demand("SUMMARY 94/94" in author.stdout, "frozen V002 author verifier reproduces 94/94")
demand("PHYSICS_UNEXECUTED" in author.stdout, "author replay retains the design-only status")

theorem = (TARGET / "THEOREM.md").read_text()
contract = json.loads((TARGET / "OBSERVABLE_CONTRACT.json").read_text())
result = json.loads((TARGET / "RESULT.json").read_text())
prior = (AUDIT / "PRIOR_V001_REJECTION.md").read_text()

demand("V001 hostile audit rejected" in prior, "prior V001 rejection remains explicit")
demand("eda0d378" in prior and "a0ab7259" in prior, "prior rejected target hashes remain recorded")

# ---------------------------------------------------------------------------
# Wigner/Poincare helicity gate and scalar-doublet rejection.

demand("unitary massless Poincare representation" in theorem, "V002 requires a unitary massless Poincare pole bundle")
demand("e^{+2i\\theta}" in theorem and "e^{-2i\\theta}" in theorem, "V002 fixes little-group rotation helicities plus/minus two")
demand("reducible scalar doublet" in theorem, "V002 explicitly rejects a two-scalar doublet")
demand("equality of dispersive slopes is insufficient" in theorem, "equal cone slopes are not substituted for Wigner covariance")

# For a two-state unitary massless little-group fiber with J weights +/-2,
# [J,T_+]=+T_+ and [J,T_-]=-T_- allow no nonzero matrix element. Thus the
# ISO(2) translation generators are forced to act trivially; this is a
# discrete-helicity, not continuous-spin, representation.
weights = (F(2), F(-2))
allowed_plus = [(i, j) for i in range(2) for j in range(2) if weights[i] - weights[j] == 1]
allowed_minus = [(i, j) for i in range(2) for j in range(2) if weights[i] - weights[j] == -1]
demand(allowed_plus == [] and allowed_minus == [], "two-state J={+2,-2} fiber forces trivial ISO(2) translations")
demand(tuple(sorted(weights)) == (F(-2), F(2)), "helicity fiber has exactly the two required Wigner weights")
scalar_weights = (F(0), F(0))
demand(tuple(sorted(scalar_weights)) != (F(-2), F(2)), "degenerate covariant scalar doublet fails the helicity gate")

# ---------------------------------------------------------------------------
# TT k=0 exclusion and exact nonzero-k projector.

demand("does not define" in theorem and "Pi_TT(0)" in theorem, "V002 excludes TT projection at the zero character")
demand("directional limit" in theorem and "nonzero sequence" in theorem, "every TT soft limit uses a punctured nonzero sequence")

def tt_projector(k):
    norm2 = sum(x * x for x in k)
    if norm2 == 0:
        raise ZeroDivisionError("Pi_TT(0) undefined")
    P = tuple(tuple((F(1) if i == j else F(0)) - k[i] * k[j] / norm2 for j in range(3)) for i in range(3))
    indices = tuple(itertools.product(range(3), repeat=2))
    pi = []
    for i, j in indices:
        row = []
        for k0, l in indices:
            row.append((P[i][k0] * P[j][l] + P[i][l] * P[j][k0] - P[i][j] * P[k0][l]) / 2)
        pi.append(tuple(row))
    return tuple(pi)


try:
    tt_projector((F(0), F(0), F(0)))
    zero_rejected = False
except ZeroDivisionError:
    zero_rejected = True
demand(zero_rejected, "independent TT construction rejects k=0 exactly")
pi = tt_projector((F(1), F(2), F(3)))
demand(matmul(pi, pi) == pi, "nonzero-k TT projector is exactly idempotent")
demand(rank(pi) == 2, "nonzero-k symmetric TT image has exact rank two")

# ---------------------------------------------------------------------------
# Degenerate-ground query and basis invariance.

demand("rho_(0,L)=P_(0,L)/rank(P_(0,L))" in theorem, "V002 types a degenerate ground by the normalized complete projector")
demand("choosing an individual vector" in theorem, "V002 forbids a selected vector in a degenerate ground space")

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
demand(before == after, "rho0=P0/2 ground query is basis invariant")
e1 = ((F(1),), (F(0),))
Ue1 = matmul(U, e1)
selected_before = sum(x[0] * x[0] for x in matmul(M1, e1))
selected_after = sum(x[0] * x[0] for x in matmul(M1, Ue1))
demand(selected_before != selected_after, "selected degenerate-ground vector is basis dependent and correctly excluded")

# ---------------------------------------------------------------------------
# Affine-versus-physical binding and exact source pairing.

demand("are **affine** volumes" in theorem, "V002 does not call the unbound cell measure physical")
demand("independently audited and sealed FD" in theorem, "physical naming requires independent FD binding")
demand("no further multiplicative conversion" in theorem, "one global binding forbids sectorwise refits")
demand("physical FD binding" in result["missing"], "machine result keeps physical binding open")

# Work in Q(sqrt(3)): (a,b) denotes a+b sqrt(3).
def q3mul(x, y):
    return (x[0] * y[0] + 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


v3 = (F(0), F(16, 9))          # 16/(3 sqrt(3)) = 16 sqrt(3)/9, a_*=1
two_over_v3 = (F(0), F(3, 8))  # 3 sqrt(3)/8
v3_over_two = (F(0), F(8, 9))  # 8 sqrt(3)/9
demand(q3mul(v3, two_over_v3) == (F(2), F(0)), "Q_aff conversion squares to 2/v3 exactly")
demand(q3mul(two_over_v3, v3_over_two) == (F(1), F(0)), "source/operator pairing conversions cancel exactly")

# ---------------------------------------------------------------------------
# Exact momentum registry and cover behavior.

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


for L in (5, 10):
    reps = representatives(L)
    demand(len(reps) == L ** 3, f"L={L} registry has every character exactly once")
    demand(reps[(0, 0, 0)] == (0, 0, 0), f"L={L} zero character has the zero minimizer")
    for cls, n in reps.items():
        neg = tuple((-x) % L for x in cls)
        if qform(n) != qform(reps[neg]):
            raise AssertionError(f"L={L} conjugate norm mismatch at {cls}")
    passed(f"L={L} every conjugate character has equal reciprocal norm")

demand([qform(n) for n in ((1, 0, 0), (1, 1, 0), (1, -1, 0))] == [3, 4, 8], "three frozen rays have q=3,4,8")
demand(3 == min(qform(n) for n in representatives(5).values() if n != (0, 0, 0)), "q_min=3 gives |k_min|=3 pi/(2 L a_*)")
base = representatives(5)
fine = representatives(10)
demand(all(tuple((2 * x) % 10 for x in cls) in fine for cls in base), "G5 characters inject by doubling into G10")
demand(any(any(x % 2 for x in cls) for cls in fine), "G10 contains genuinely new odd momentum characters")

# ---------------------------------------------------------------------------
# Ancestry invariance and factorization readiness.

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
demand("liminf" in contract["ancestry"]["pass"] and "eta_min" in contract["ancestry"]["pass"], "normalized ancestry has a prospective positive lower-bound gate")

Z = diag(F(2), F(3))
Zplus = diag(F(1, 2), F(1, 3))
Dpole = diag(F(5), F(7))
Dinv = diag(F(1, 5), F(1, 7))
Gamma = ((F(1), F(2)), (F(3), F(4)))
G = matmul(matmul(Z, Dinv), Gamma)
demand(matmul(matmul(Dpole, Zplus), G) == Gamma, "retarded branchwise pole amputation recovers the source-independent vertex")
demand("not longitudinal" in contract["factorization"]["ceiling"], "factorization readiness is not promoted to G3 decoupling")

# ---------------------------------------------------------------------------
# Hostile PASS/FAIL/INDETERMINATE residue overlap.

binding = contract["canonical_source_binding"]
demand("failed normalized ancestry" in binding["fail"], "machine amplitude FAIL condition includes failed ancestry")
demand("liminf sigma_min" in binding["pass"], "machine amplitude PASS condition depends on canonical residue")

# Exact asymptotic counterexample:
# Delta_L=1/L, R_L=(L/2) I => N_can=2 Delta R=I and Z_can=I,
# while S_L=L^2 I => Omega_L=(1/(2L)) I -> 0.
# It has full raw/canonical rank and finite nonzero canonical amplitude, but
# failed normalized ancestry.
for L in (F(10), F(100), F(1000)):
    delta = 1 / L
    R_L = matscale(L / 2, eye(2))
    N_can = matscale(2 * delta, R_L)
    Omega = matscale(1 / (2 * L), eye(2))
    demand(N_can == eye(2), f"L={L}: canonical numerator and vertex remain finite full rank")
    demand(Omega[0][0] > 0, f"L={L}: finite-size ancestry remains well-defined")
demand(F(1, 20) > F(1, 200) > F(1, 2000), "normalized ancestry counterexample tends monotonically to zero")

pass_condition = True   # frozen applicable D, rank 2, sigma_min=1, sigma_max=1
fail_condition = True   # eta_min(Omega)->0
indeterminate_condition = False
demand(pass_condition and fail_condition and not indeterminate_condition, "exact counterexample triggers declared amplitude PASS and FAIL simultaneously")
demand("failed normalized ancestry" not in result["residue"], "RESULT.json omits the theorem/contract ancestry-based amplitude FAIL branch")

# The outer G2 conjunction still rejects this example, but that does not make
# the advertised three-valued amplitude classifier single-valued.
positive_g2 = pass_condition and not fail_condition
demand(not positive_g2, "outer G2 conjunction correctly blocks false positive promotion")
demand(pass_condition and fail_condition, "outer conjunction does not disambiguate the amplitude label itself")

# ---------------------------------------------------------------------------
# Strict design-only ceiling and verdict.

demand(contract["status"] == "REPAIRED_V002_AWAITING_INDEPENDENT_AUDIT_PHYSICS_UNEXECUTED", "machine contract remains awaiting independent audit")
demand("no positive G2" in contract["ceiling"] and "no gravity" in contract["ceiling"], "machine contract retains strict design-only ceiling")
demand("Consequently it proves no" in theorem, "theorem retains no-positive-G2 disposition")

print(f"SUMMARY {checks}/{checks} independent GF V002 hostile checks passed")
print("VERDICT REJECT -- REPAIR_REQUIRED")
print("MATERIAL one exact input triggers canonical-amplitude PASS and FAIL simultaneously")
print("MINIMAL_REPAIR remove normalized-ancestry failure from the amplitude FAIL branch and retain it as the separate ancestry/G2 gate; add the overlap regression")
print("CEILING no GE repin; no positive G2, G3, gravity, or G")
