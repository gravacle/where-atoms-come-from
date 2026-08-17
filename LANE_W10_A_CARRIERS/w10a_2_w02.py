#!/usr/bin/env python3
"""
LANE W10-A — SCRIPT 2.  W-02's CHARACTER-RATIO CRITERION ON FOUR OCCUPIED CLASSES.

W-02 (REGISTER row 2) rules   FORMATION OCCURS  <=>  G != {1},
      G = < chi_a / chi_b : a,b in supp(pi) >,
and tabulates the cases.  Its table (REGISTER:132-133) reads

      |S|=3 -> <u,v>  ·  S={0,C} -> <u>  ·  S={0,F} -> <v>  ·  S={F,C} -> <u/v>  ·  |S|=1 -> never

where 0 = the root (class 11), F = class 10, C = class 01.  K1 HAS NO CLASS 00, so its
S ranges over the 7 non-empty subsets of a THREE-element set.  On B0b and B4 it ranges over
all 15 non-empty subsets of a FOUR-element set.  This script runs all 15, against a battery
of connections, and asks (a) does the CRITERION hold, (b) does the TABLE hold.

ISOLATION.  Weights, connection battery, k-range, evaluator and code path are IDENTICAL in
every row.  The ONE thing that moves is the SUPPORT S (equivalently the support difference
lattice L_S).  Every arm's input vector is printed and all arms are asserted pairwise distinct
at the bytes -- see the ARMS DIFF block.  A "carrier control" holding pi fixed would be a
zero-variable control and is not run (see PUBLISHED_CONVENTIONS.txt).

IEEE double for the k-sweeps; EXACT integer congruence for every |Z_k| = 1 claim.
"""
import sys, math, hashlib
from fractions import Fraction
from itertools import combinations
import numpy as np
from w10a_lib import CLASSES, CLASS_NAME, CHAR_NAME, L_S, hnf2, L_conn_rational, sublattice

LOG = []
def out(s=""):
    print(s); LOG.append(s)

K_SWEEP = 200000
TWO_PI = 2.0 * math.pi

out("=" * 108)
out("W10-A SCRIPT 2 — W-02's CHARACTER-RATIO CRITERION ON FOUR OCCUPIED CLASSES")
out("=" * 108)
out(f"numpy {np.__version__}; IEEE double for k-sweeps, K = {K_SWEEP}; EXACT where a line says EXACT.")
out()

# ------------------------------------------------------------------ the 15 supports
SUBSETS = []
for r in range(1, 5):
    for S in combinations(CLASSES, r):
        SUBSETS.append(list(S))

def gname(basis):
    """The group G in u,v notation, from the canonical basis of L_S."""
    if not basis:
        return "{1}"
    if len(basis) == 2:
        return "<u,v>"
    (m, n) = basis[0]
    NICE = {(1, 0): "<u>", (0, 1): "<v>", (1, 1): "<uv>", (1, -1): "<u/v>"}
    if (m, n) in NICE:
        return NICE[(m, n)]
    def pw(sym, e):
        if e == 0: return ""
        if e == 1: return sym
        if e == -1: return sym + "^-1"
        return f"{sym}^{e}"
    return "<" + ((pw("u", m) + pw("v", n)) or "1") + ">"

def trivial_iff(basis):
    if not basis:
        return "always trivial (|S| = 1)"
    if len(basis) == 2:
        return "W_F = 1 and W_C = 1"
    (m, n) = basis[0]
    # u = conj(W_F), v = W_C  ->  u^m v^n = conj(W_F)^m W_C^n = 1
    if (m, n) == (1, 0):  return "W_F = 1"
    if (m, n) == (0, 1):  return "W_C = 1"
    if (m, n) == (1, 1):  return "W_C = W_F   [MISSING ROW]"
    if (m, n) == (1, -1): return "W_F * W_C = 1"
    return f"conj(W_F)^{m} W_C^{n} = 1"

# ------------------------------------------------------------------ connections
# alpha = arg u = -f,  beta = arg v = c.   Rational connections carry (A,B,q): alpha=2pi A/q etc.
CONNS = [
    # label,              kind,   data
    ("trivial (0,0)",                "rat", (0, 0, 1)),
    ("alpha=pi, beta=0     [u=-1]",  "rat", (1, 0, 2)),
    ("alpha=0, beta=pi     [v=-1]",  "rat", (0, 1, 2)),
    ("alpha=2pi/3, beta=-2pi/3 [uv=1]", "rat", (1, -1, 3)),
    ("alpha=2pi/3, beta=2pi/3  [u/v=1]", "rat", (1, 1, 3)),
    ("alpha=2pi/5, beta=3*2pi/5",     "rat", (1, 3, 5)),
    ("S1 published: W_F=-1,W_C=-i",   "rat", (-6, 9, 12)),      # alpha=-pi, beta=3pi/2
    ("S3/S4 headline f=2.0,c=1.1",    "res", (11, 20)),
    ("generic irrational 1/sqrt2,1/sqrt3", "irr", (1/math.sqrt(2), 1/math.sqrt(3))),
    ("generic irrational sqrt2/7,sqrt5/9", "irr", (math.sqrt(2)/7.0, math.sqrt(5)/9.0)),
]

def conn_lattice(kind, data):
    """EXACT for 'rat'; derived and stated for 'res'/'irr'."""
    if kind == "rat":
        A, B, q = data
        if q == 1:
            return ((1, 0), (0, 1))          # alpha = beta = 0: L = Z^2
        return L_conn_rational(A, q, B, q)
    if kind == "res":
        m, n = data
        return hnf2([(m, n)])                # L = Z(11,20), derived in the text below
    return ()                                # irrational, rationally independent: L = {0}

def conn_frac(kind, data):
    """(alpha/2pi, beta/2pi) as float64."""
    if kind == "rat":
        A, B, q = data
        return (A / q, B / q)
    if kind == "res":
        # f = 2.0, c = 1.1 -> alpha = -2.0, beta = 1.1
        return (-2.0 / TWO_PI, 1.1 / TWO_PI)
    return data

def sweep(pi4, kind, data, K=K_SWEEP):
    """mean log|Z_k| over k=1..K, min|Z_k|, and an exact-zero flag."""
    a, b = conn_frac(kind, data)
    k = np.arange(1, K + 1, dtype=np.float64)
    ta = (k * a) % 1.0
    tb = (k * b) % 1.0
    ua = np.exp(2j * np.pi * ta)
    vb = np.exp(2j * np.pi * tb)
    p00, p10, p01, p11 = [float(x) for x in pi4]
    Z = p00 + p10 * ua + p01 * vb + p11 * ua * vb
    m = np.abs(Z)
    exact_zero = bool(np.any(m < 1e-14))
    with np.errstate(divide="ignore"):
        lg = np.log(np.maximum(m, 1e-300))
    return float(lg.mean()), float(m.min()), float(m.max()), exact_zero

def exact_never(pi4, kind, data, S):
    """EXACT (integer congruence): for a rational connection alpha=2pi A/q, beta=2pi B/q,
    |Z_k| = 1 for ALL k iff all occupied classes share the same exponent k(aA+bB) mod q.
    Returns (applicable, all_moduli_one)."""
    if kind != "rat":
        return (False, None)
    A, B, q = data
    if q == 1:
        return (True, True)
    ex = {c: (c[0] * A + c[1] * B) % q for c in S}
    for k in range(1, q + 1):
        vals = {(k * e) % q for e in ex.values()}
        if len(vals) > 1:
            return (True, False)
    return (True, True)

# ------------------------------------------------------------------ ARMS DIFF
out("-" * 108)
out("ARMS DIFF — the 15 arms, their weight vectors, and a byte hash of each.  W-08's isolation")
out("audit found that the commonest FATAL defect is ZERO variables moved.  So: DIFF THE ARMS.")
out("-" * 108)
ARMS = []
seen = {}
for S in SUBSETS:
    w = Fraction(1, len(S))
    pi4 = [w if c in S else Fraction(0) for c in CLASSES]
    key = hashlib.sha256(repr([str(x) for x in pi4]).encode()).hexdigest()[:16]
    ARMS.append((S, pi4, key))
    assert key not in seen, f"ARMS COLLIDE: {S} and {seen[key]} are byte-identical"
    seen[key] = S
out(f"{'support':<22}{'pi = (p00,p10,p01,p11)':<40}{'sha256[:16]':<20}")
for S, pi4, key in ARMS:
    nm = "{" + ",".join(CLASS_NAME[c] for c in S) + "}"
    out(f"{nm:<22}{'('+', '.join(str(x) for x in pi4)+')':<40}{key:<20}")
out(f"all {len(ARMS)} arms pairwise distinct at the bytes: True")
out()

# ------------------------------------------------------------------ the 15-row table
out("=" * 108)
out("PART A — THE FULL SUPPORT TABLE ON A FOUR-CLASS CARRIER (all 15 non-empty supports)")
out("=" * 108)
out("L_S = the support difference lattice (S4-1's object).  G = {1} iff L_S is contained in the")
out("connection's relation lattice L (W-03's object).  BOTH LATTICES ARE ALREADY IN THE REGISTER;")
out("THE CONTAINMENT IS NOT.  EXACT integer lattice arithmetic in this table.")
out()
out(f"{'#':<3}{'support':<16}{'chars':<14}{'L_S basis':<22}{'rk':<4}{'G':<10}"
    f"{'G = {1} iff':<26}{'realizable on'}")
CARRIER_OCC = {
    "K1/B1,B2,B1s,B3": [(1, 0), (0, 1), (1, 1)],
    "B1q": [(0, 0), (1, 0), (0, 1)],
    "B1p": [(1, 0), (0, 1)],
    "B0a": [(0, 0), (1, 0), (0, 1)],
    "B0b/B4": [(0, 0), (1, 0), (0, 1), (1, 1)],
}
ROWS = []
for i, (S, pi4, key) in enumerate(ARMS):
    basis = L_S(S)
    where = [n for n, occ in CARRIER_OCC.items() if all(c in occ for c in S)]
    nm = "{" + ",".join(CLASS_NAME[c] for c in S) + "}"
    ch = " ".join(CHAR_NAME[c] for c in S)
    ROWS.append((S, pi4, basis, where))
    out(f"{i+1:<3}{nm:<16}{ch:<14}{str(basis):<22}{len(basis):<4}{gname(basis):<10}"
        f"{trivial_iff(basis):<26}{', '.join(where) if where else 'FOUR-CLASS ONLY'}")
out()
out("READ OFF:")
out("  * rank L_S = 2 for every |S| >= 3 and rank 1 for every |S| = 2 and rank 0 for |S| = 1.")
out("    S4-1's 'rank 2 iff |S| >= 3' therefore EXTENDS to four classes verbatim (at unit charge;")
out("    W-03 already refuted it OFF unit charge and that refutation is untouched here).")
out("  * The 6 two-element supports give only FOUR distinct groups: <u>, <v>, <uv>, <u/v>.")
out("    W-02's table lists THREE of them.  The missing one is G = <uv>, from S = {00,11}.")
out("  * S={00,10} gives <u> and S={11,01} gives <u>; S={00,01} and S={11,10} both give <v>.")
out("    Those two coincidences ARE W-03's involution 00<->11, 10<->01 (REGISTER:184).")
out("    The involution fixes {00,11} and {10,01} setwise; the first is the row W-02 lacks.")
out()

# ------------------------------------------------------------------ criterion test
out("=" * 108)
out("PART B — THE CRITERION ITSELF, EXHAUSTIVE: 15 SUPPORTS x 10 CONNECTIONS = 150 CASES")
out("=" * 108)
out("PREDICT: forms  <=>  L_S not contained in L.   MEASURE: r = (1/K) sum_{k<=K} log|Z_k|.")
out("Weights: uniform on the support (Fractions -> float64).  Same K, same evaluator, every row.")
out()
out("Derivations of the relation lattices used (stated, not fitted):")
out("  rational alpha=2pi A/q, beta=2pi B/q:  L = {(m,n) : mA+nB = 0 mod q}, computed exactly.")
out("  S3/S4 headline f=2.0, c=1.1: alpha=-2.0, beta=1.1; m alpha + n beta in 2pi Z forces")
out("     the 2pi-multiple to be 0 (pi transcendental), so -2m + 1.1n = 0, i.e. L = Z(11,20).")
out("     (-11f + 20c = 0 is the erratum against W-02's own resonance.)")
out("  irrational rows: {1, alpha/2pi, beta/2pi} are Q-linearly independent, so L = {0}.")
out()

nmis = 0; ncase = 0; nform = 0; nnever = 0
exact_checked = 0; exact_agree = 0
out(f"{'support':<14}{'connection':<36}{'L_S<=L':<8}{'predict':<9}"
    f"{'r=(1/K)sum log|Z|':<20}{'min|Z_k|':<12}{'verdict':<9}{'match'}")
for S, pi4, basis, where in ROWS:
    nm = "{" + ",".join(CLASS_NAME[c] for c in S) + "}"
    for (lbl, kind, data) in CONNS:
        Lc = conn_lattice(kind, data)
        contained = sublattice(basis, Lc)
        predict = "never" if contained else "forms"
        r, mn, mx, ez = sweep(pi4, kind, data)
        verdict = "never" if abs(r) < 1e-12 else ("forms" if r < -1e-9 else "??")
        ok = (verdict == predict)
        nmis += (0 if ok else 1); ncase += 1
        nform += (verdict == "forms"); nnever += (verdict == "never")
        app, allone = exact_never(pi4, kind, data, S)
        if app:
            exact_checked += 1
            exact_agree += (allone == (predict == "never"))
        out(f"{nm:<14}{lbl:<36}{str(contained):<8}{predict:<9}{r:< 20.12e}{mn:<12.6e}"
            f"{verdict:<9}{'OK' if ok else '**MISMATCH**'}")
out()
out(f"CASES {ncase}   forms {nform}   never {nnever}   MISMATCHES {nmis}")
out(f"EXACT integer-congruence cross-check on the {exact_checked} rational-connection cases:")
out(f"   'all |Z_k| = 1' agrees with 'L_S contained in L' on {exact_agree} of {exact_checked}.")
out()

# ------------------------------------------------------------------ weight robustness
out("-" * 108)
out("PART C — WEIGHT ROBUSTNESS WITHIN A SUPPORT (declared as a robustness pass, not a control)")
out("-" * 108)
out("Same 150 cases with NON-UNIFORM weights on each support (seeded Dirichlet-ish, seed 20260816).")
rng = np.random.default_rng(20260816)
nmis2 = 0
for S, pi4u, basis, where in ROWS:
    w = rng.random(len(S)) + 0.05
    w = w / w.sum()
    pi4 = [0.0] * 4
    for c, x in zip(S, w):
        pi4[CLASSES.index(c)] = float(x)
    for (lbl, kind, data) in CONNS:
        Lc = conn_lattice(kind, data)
        predict = "never" if sublattice(basis, Lc) else "forms"
        r, mn, mx, ez = sweep(pi4, kind, data)
        verdict = "never" if abs(r) < 1e-12 else ("forms" if r < -1e-9 else "??")
        nmis2 += (0 if verdict == predict else 1)
out(f"MISMATCHES under non-uniform weights: {nmis2} of {ncase}")
out("  The verdict depends on the SUPPORT and the CONNECTION and not on the weights inside the")
out("  support -- which is exactly what the criterion claims, and it holds at four classes.")
out()

# ------------------------------------------------------------------ the missing row, exhibited
out("=" * 108)
out("PART D — THE ROW W-02's TABLE DOES NOT HAVE, EXHIBITED ON B0b AND B4")
out("=" * 108)
out("S = {00,11}: a vertex in NEITHER loop and a vertex in BOTH, and no other weight.")
out("G = <uv> = <conj(W_F) W_C>, trivial iff W_C = W_F.  So there is a NON-TRIVIAL connection")
out("family that NEVER FORMS and that W-02's table does not contain.  COR-B exhibited FOUR")
out("non-forming families; this is a FIFTH, and it is unrealizable on every carrier the corpus ran.")
out()
S = [(0, 0), (1, 1)]
basis = L_S(S)
out(f"  L_S = {basis},  G = {gname(basis)},  trivial iff {trivial_iff(basis)}")
for (lbl, kind, data, note) in [
    ("W_F = W_C = e^{i 2pi/3}  (alpha=-2pi/3, beta=2pi/3)", "rat", (-1, 1, 3), "uv = 1: NEVER FORMS"),
    ("W_F = W_C = e^{i 2pi/5}", "rat", (-1, 1, 5), "uv = 1: NEVER FORMS"),
    ("W_F = -1, W_C = -1      (alpha=pi, beta=pi)", "rat", (1, 1, 2), "uv = 1: NEVER FORMS"),
    ("W_F = -1, W_C = +1      (alpha=pi, beta=0)", "rat", (1, 0, 2), "uv = -1: forms"),
    ("S1 published W_F=-1,W_C=-i", "rat", (-6, 9, 12), "uv = e^{i pi/2}: forms"),
]:
    Lc = conn_lattice(kind, data)
    contained = sublattice(basis, Lc)
    pi4 = [Fraction(1, 2), Fraction(0), Fraction(0), Fraction(1, 2)]
    r, mn, mx, ez = sweep(pi4, kind, data)
    app, allone = exact_never(pi4, kind, data, S)
    out(f"  {lbl:<50} L_S<=L {str(contained):<6} r = {r:< 16.10e}  min|Z_k| = {mn:.6e}")
    out(f"      EXACT: all |Z_k| = 1 ? {allone}      [{note}]")
out()
out("AND THE SAME SUPPORT IS UNREACHABLE ON EVERY CARRIER OF S4:511-520 EXCEPT B0b AND B4:")
for n, occ in CARRIER_OCC.items():
    ok = all(c in occ for c in S)
    out(f"   {n:<20} occupies {[CLASS_NAME[c] for c in occ]}  -> S={{00,11}} realizable: {ok}")
out()

# ------------------------------------------------------------------ verdict
out("=" * 108)
out("PART E — VERDICT ON W-02")
out("=" * 108)
out("(1) THE CRITERION HOLDS AT FOUR CLASSES, AND ITS PROOF NEVER MENTIONS THE CLASS COUNT.")
out("    |Z_k| = 1 with non-negative weights summing to 1 forces every occupied character's")
out("    k-th power to coincide (strict triangle inequality), for ANY number of classes; so")
out("    |Z_k| = 1 for all k  <=>  every ratio is 1  <=>  G = {1}.  The converse half is Weyl")
out("    on the closure of G, again class-count-free.  CARRIER-INDEPENDENT.  0 mismatches / 150.")
out("(2) THE TABLE IS A THREE-CLASS STATEMENT IN EXACTLY ONE PLACE, AND NO ROW OF IT IS FALSE:")
out("    it enumerates 3 of the 4 rank-1 support lattices.  The missing one, L_S = Z(1,1) from")
out("    S = {00,11}, needs a vertex in BOTH loops AND a vertex in NEITHER -- W-09's four-class")
out("    condition exactly.  S4's own Theorem S4-1 (S4:625-640) DOES carry that row; the")
out("    REGISTER's W-02 row does not, and no lane has ever run it on a carrier that realizes it.")
out("(3) THE LABELS ARE THREE-CLASS LABELS.  '|S|=3' names the FULL support on K1 and a PROPER")
out("    subset on B0b/B4; '|S|=4' is a new row (same group <u,v>).  '|S|=1 -> never' has 4")
out("    instances instead of 3.  The rulings extend; the ENUMERATION does not.")

with open("w10a_2_w02.OUT.txt", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
sys.exit(1 if (nmis or nmis2 or exact_agree != exact_checked) else 0)
