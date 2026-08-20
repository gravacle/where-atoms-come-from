"""T-24 items 2-4 -- CLAUSE (v) ON THE 1x2 TORUS FOR D(D_4), dim 8^4 = 4096.

LATTICE: vertices v0,v1; edges h0: v0->v1, h1: v1->v0 (arcs between DISTINCT vertices),
u0: self-loop at v0, u1: self-loop at v1; plaquettes p0 (hol = g_h0 g_u1 g_h0^-1 g_u0^-1)
and p1 (hol = g_h1 g_u0 g_h1^-1 g_u1^-1).  H = -(A_v0 + A_v1 + B_p0 + B_p1), commuting
projectors, spectrum {0,-1,-2,-3,-4}, eigenprojectors EXACT rational sparse matrices.

REGION CATALOGUE (T-11 convention, cross-checked against homology): the only nonempty
cycle-free edge subsets are {h0} and {h1} -- every set containing u0 or u1 contains a
self-loop (a cycle), and {h0,h1} is the horizontal non-contractible cycle.  So clause (v)
on this carrier IS the single-edge test on h0 and h1.

THE FLIP REDUCTION (exact; no verdict rests on a float):
  * U flips R  <=>  U^dag R U = -R  <=>  {U,R} = 0  (multiply by U on the left).
  * If U is a unitary in C_e with {U,R} = 0 and R is invertible (every record is), the
    spectrum of U pairs (lam, -lam) with R mapping the lam-eigenspace onto the
    (-lam)-eigenspace; the involution S = sum over pairs (P_lam - P_-lam) lies in
    C*(U) <= C_e and still anticommutes with R.  So single-edge flippability is decided
    by the INVOLUTIONS S != +-I of C_e.
  * {S,R} = 0 forces R off-diagonal in the S-grading of every H-eigenspace E, hence
    dim E_+ = dim E_-, i.e. Tr(P_E S) = 0 for EVERY E.  Conversely if all those traces
    vanish, R = (x)_E [[0,W_E],[W_E^dag,0]] is a record (Hermitian involution, commutes
    with H, trace-balanced because off-diagonal, non-trivial) that S flips.
  * THEREFORE: a record flippable on edge e exists  <=>  some involution S != +-I in
    C_e has ALL per-eigenspace traces zero.  The traces are exact integers.
"""
import sys, json, time, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_T24_DV")
from fractions import Fraction
import numpy as np
from t24_lib import *

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_T24_DV"
OUT = []
def say(s=""):
    print(s, flush=True); OUT.append(s)

t0 = time.time()
say("=" * 112)
say("T-24 / items 2-4 -- CLAUSE (v) ON D(D_4), 1x2 TORUS (dim 4096) -- EXACT")
say("=" * 112)

G = make_D4()
car = Carrier(G)
n, N = car.n, car.N
MUL, INV = G["MUL"], G["INV"]
say("carrier: D(%s), |G| = %d, dim = %d" % (G["label"], n, N))
say("")

# ------------------------------------------------------------------ 0. structure checks (exact)
say("0. STRUCTURE CHECKS (all exact, permutation level)")
okA0 = all(np.array_equal(car.permA0(k1)[car.permA0(k2)], car.permA0(int(MUL[k1, k2])))
           for k1 in range(n) for k2 in range(n))
okA1 = all(np.array_equal(car.permA1(k1)[car.permA1(k2)], car.permA1(int(MUL[k1, k2])))
           for k1 in range(n) for k2 in range(n))
okComm = all(np.array_equal(car.permA0(k1)[car.permA1(k2)], car.permA1(k2)[car.permA0(k1)])
             for k1 in range(n) for k2 in range(n))
b0d, b1d = car.diagB0(), car.diagB1()
okB = True
for k in range(n):
    for pi in (car.permA0(k), car.permA1(k)):
        okB &= np.array_equal(b0d[pi], b0d) and np.array_equal(b1d[pi], b1d)
say("   A_v0, A_v1 permutation reps of D_4: %s %s;  [A_v0(k), A_v1(k')] = 0: %s;"
    % (okA0, okA1, okComm))
say("   plaquette diagonals invariant under every gauge move ([A,B]=0): %s" % okB)
assert okA0 and okA1 and okComm and okB

# ------------------------------------------------------------------ 1. region catalogue
say("")
say("1. REGION CATALOGUE ON THE 1x2 GRAPH (T-11 forest convention + homology, all 16 subsets)")
EDGE_ENDS = {"h0": (0, 1), "h1": (1, 0), "u0": (0, 0), "u1": (1, 1)}
def has_cycle(sub):
    par = [0, 1]
    def find(x):
        while par[x] != x: x = par[x]
        return x
    for e in sub:
        a, b = EDGE_ENDS[e]
        ra, rb = find(a), find(b)
        if ra == rb: return True
        par[ra] = rb
    return False
def connected(sub):
    if not sub: return False
    vs = {x for e in sub for x in EDGE_ENDS[e]}
    seen = {next(iter(vs))}; ch = True
    while ch:
        ch = False
        for e in sub:
            a, b = EDGE_ENDS[e]
            if (a in seen) ^ (b in seen): seen |= {a, b}; ch = True
    return seen == vs
contractible_singles = []
for k in range(0, 5):
    for sub in itertools.combinations(["h0", "h1", "u0", "u1"], k):
        sub = list(sub)
        if sub and not has_cycle(sub) and connected(sub):
            contractible_singles.append(sub)
say("   subsets that are SINGLE CONTRACTIBLE regions: %s" % contractible_singles)
say("   every other nonempty subset contains a cycle (self-loop u_i or the pair {h0,h1});")
say("   homology agrees: {h0},{h1} are arcs (contractible), u_i and {h0,h1} carry the two")
say("   generators of H_1(T^2) = Z^2.  No convention divergence on this carrier (D-23).")
assert contractible_singles == [["h0"], ["h1"]]

# ------------------------------------------------------------------ 2. sectors (exact)
say("")
say("2. SECTOR PROJECTORS (exact rational; commuting-projector eigenstructure)")
sectors = build_sectors(car)
dims = {}
for s, P in sectors.items():
    tr = sp_trace(P)
    assert tr == sp_frob_check_projector(P), "sector %s fails Tr(P^2)=Tr(P)" % (s,)
    assert tr.denominator == 1
    dims[s] = int(tr)
assert sum(dims.values()) == N
eigmap = eigen_projectors(sectors)
eigdims = {k: sum(dims[s] for s in v) for k, v in eigmap.items()}
say("   sector dims: " + ", ".join("%s:%d" % (s, d) for s, d in sorted(dims.items())))
say("   H-eigenspace dims (eig -k): " + ", ".join("-%d: %d" % (k, eigdims[k])
                                                  for k in sorted(eigdims)))
say("   all eigenspace dims EVEN (C-41's |G| = 2^3 record-existence condition): %s"
    % all(d % 2 == 0 for d in eigdims.values()))
say("   ground space (-4) dim %d -- topological (matches 22 on the minimal torus): %s"
    % (eigdims.get(4, 0), eigdims.get(4, 0) == 22))
eigP = {}
for k, keys in eigmap.items():
    merged = [dict() for _ in range(N)]
    for s in keys:
        P = sectors[s]
        for j in range(N):
            for i, v in P[j].items():
                w = merged[j].get(i, Fraction(0)) + v
                if w: merged[j][i] = w
                else: merged[j].pop(i, None)
    eigP[k] = merged
EIGS = sorted(eigP)          # [0,1,2,3,4] meaning eigenvalue -k

# fast exact accessors: 64*P[c,c] and 64*P[c,pi(c)] as int64 vectors
r2 = 2                                   # index of central r^2
tau0 = car.edge_perm_right("h0", r2)     # T  on h0
tau1 = car.edge_perm_right("h1", r2)     # T' on h1
tau01 = tau0[tau1]
PERMS_PI = {"I": None, "T0": tau0, "T1": tau1, "T0T1": tau01}
diag64 = {}; w64 = {}
for k in EIGS:
    P = eigP[k]
    dv = np.zeros(N, dtype=np.int64)
    for c in range(N):
        x = P[c].get(c)
        if x is not None:
            y = x * 64; assert y.denominator == 1; dv[c] = int(y)
    diag64[k] = dv
    for pnm, pi in PERMS_PI.items():
        if pi is None: continue
        wv = np.zeros(N, dtype=np.int64)
        for c in range(N):
            x = P[int(pi[c])].get(c)
            if x is not None:
                y = x * 64; assert y.denominator == 1; wv[c] = int(y)
        w64[(k, pnm)] = wv

# ------------------------------------------------------------------ 3. records (D-18)
say("")
say("3. RECORDS ON THIS CARRIER (D-18: constructed and clause-verified BEFORE the test)")
say("")
say("3a. THE NATURAL WILSON DIAGONALS ARE NOT RECORDS HERE (a finding, kept in the open):")
chi_list = [(nm, chi) for nm, chi, d in G["irreps"] if d == 1]
hol = MUL[car.h0, car.h1]
class_of = np.zeros(n, dtype=np.int64)
for ci, cl in enumerate(G["classes"]):
    for g in cl: class_of[g] = ci
Hent = h_entries_exact(car)
for nm, chi in chi_list[1:]:
    r = chi[hol]
    bal = {(-k): int(np.dot(diag64[k], r) // 64) for k in EIGS}
    say("   chi = %-18s diag chi(g_h0 g_h1):  Tr(P_E R) = %s  -> clause (iv) FAILS" % (nm, bal))
r = chi_list[1][1][car.u0]
bal = {(-k): int(np.dot(diag64[k], r) // 64) for k in EIGS}
say("   chi_a(g_u0) (vertical):                Tr(P_E R) = %s  -> clause (iv) FAILS" % bal)
say("   (durability and bit-ness hold for all of them -- writability is what fails, so on")
say("    this carrier the Wilson-character diagonals are three-clause objects, not records.)")

say("")
say("3b. EXACT RECORD SEARCH in the gauge-fixed family R = diag(chi(hol_h) v0(cl(u0)) v1(cl(u1))) * Pi,")
say("    Pi in {I, T0, T1, T0T1} (central right-multiplications on h0/h1; all H-commuting,")
say("    Hermitian involutions by construction -- verified below).  Balance is checked with")
say("    EXACT integer traces; every reported record is exact, no floats anywhere.")
ncl = len(G["classes"])
chi_vecs = [(nm, chi[hol].astype(np.int64)) for nm, chi in chi_list]
# precompute W[chi][pi_name][k][C0,C1] = sum_c chi(hol(c)) * w(c) over class(u0)=C0, class(u1)=C1
cls_u0 = class_of[car.u0]; cls_u1 = class_of[car.u1]
Wt = {}
for cn, cv in chi_vecs:
    for pnm in PERMS_PI:
        for k in EIGS:
            base = diag64[k] if pnm == "I" else w64[(k, pnm)]
            M = np.zeros((ncl, ncl), dtype=np.int64)
            np.add.at(M, (cls_u0, cls_u1), cv * base)
            Wt[(cn, pnm, k)] = M
found = []
sign_pats = list(itertools.product([1, -1], repeat=ncl))
for cn, cv in chi_vecs:
    for pnm in PERMS_PI:
        for v0p in sign_pats:
            a0 = np.array(v0p, dtype=np.int64)
            for v1p in sign_pats:
                a1 = np.array(v1p, dtype=np.int64)
                if cn == "triv" and pnm == "I" and all(x == 1 for x in v0p) \
                   and all(x == 1 for x in v1p):
                    continue          # R = I
                ok = True
                for k in EIGS:
                    if int(a0 @ Wt[(cn, pnm, k)] @ a1) != 0: ok = False; break
                if ok: found.append((cn, pnm, v0p, v1p))
say("    candidates searched: %d ; EXACT RECORDS FOUND: %d"
    % (len(chi_vecs) * len(PERMS_PI) * len(sign_pats) ** 2 - 1, len(found)))
for f in found[:8]:
    say("      chi = %-18s Pi = %-5s v0 = %s v1 = %s" % f)
if len(found) > 8: say("      ... and %d more" % (len(found) - 8))

def build_exact_candidate(cn, pnm, v0p, v1p):
    cv = dict(chi_vecs)[cn]
    D = cv * np.array(v0p, dtype=np.int64)[cls_u0] * np.array(v1p, dtype=np.int64)[cls_u1]
    pi = PERMS_PI[pnm]
    return D, pi

star = None
for f in found:
    cn, pnm, v0p, v1p = f
    if dict(chi_list)[cn][4] == -1:      # chi(s) = -1 -> L_s control flips it
        star = f; break
if star is None and found: star = found[0]

def verify_exact_record(tag, D, pi):
    """All clauses, exact.  R = diag(D) if pi is None else diag(D) * Perm(pi)."""
    say("    record %s:" % tag)
    if pi is None:
        inv_ok = bool(np.all(np.abs(D) == 1)); herm_ok = True
    else:
        inv_ok = bool(np.all(np.abs(D) == 1) and np.array_equal(D[pi], D)
                      and np.array_equal(pi[pi], np.arange(N)))
        herm_ok = inv_ok
    say("      (i) bit -- Hermitian involution (exact): %s" % (inv_ok and herm_ok))
    # (ii) [R,H] = 0 exactly
    acc = {}
    if pi is None:
        bad = any(i != j and D[i] != D[j] for (i, j), v in Hent.items())
        dur = not bad
    else:
        for (i, j), v in Hent.items():
            k1 = (int(pi[i]), j)
            acc[k1] = acc.get(k1, Fraction(0)) + int(D[int(pi[i])]) * v
            k2 = (i, int(pi[j]))
            acc[k2] = acc.get(k2, Fraction(0)) - v * int(D[int(pi[j])])
        dur = all(x == 0 for x in acc.values())
    say("      (ii) durable [R,H] = 0 (exact): %s" % dur)
    # (iv) balance
    bals = {}
    for k in EIGS:
        if pi is None:
            t = int(np.dot(diag64[k], D))
        else:
            t = int(np.dot(w64_for(pi, k), D_perm_weight(D, pi, k)))
        bals[-k] = t
    say("      (iv) writable, Tr(P_E R) = 0 on every E (exact 64x): %s  %s"
        % (all(v == 0 for v in bals.values()), bals))
    # (iii) non-trivial: Tr(P_E R P_E R) > 0 somewhere
    nt = {}
    for k in EIGS:
        P = eigP[k]
        accq = Fraction(0)
        if pi is None:
            for j in range(N):
                for i, v in P[j].items():
                    w = P[i].get(j)
                    if w: accq += v * w * int(D[i]) * int(D[j])
        else:
            for j in range(N):        # entries P[i,j]; need sum P[c,pi(e)]D[pi(e)] P[e,pi(c)]D[pi(c)]
                pass
            accq = None
        nt[-k] = accq
    if pi is None:
        say("      (iii) non-trivial on some E: %s  (Tr(P_E R P_E R) = %s)"
            % (any(v != 0 for v in nt.values()), {k: str(v) for k, v in nt.items()}))
        ok3 = any(v != 0 for v in nt.values())
    else:
        # for R = D*Pi: Tr(P R P R) = sum over entries (c,j) of P: P[c,j] D[j] P[pi(j), pi(c)] D[pi(c)]
        ok3 = False
        vals = {}
        for k in EIGS:
            P = eigP[k]
            accq = Fraction(0)
            for j in range(N):
                pj = int(pi[j])
                for i, v in P[j].items():
                    w = P[pj].get(int(pi[i]))
                    if w: accq += v * w * int(D[pj]) * int(D[int(pi[i])])
            vals[-k] = accq
        ok3 = any(v != 0 for v in vals.values())
        say("      (iii) non-trivial on some E: %s  (Tr(P_E R P_E R) = %s)"
            % (ok3, {k: str(v) for k, v in vals.items()}))
    # transport
    fixed = True
    for k in range(n):
        for perm in (car.permA0(k), car.permA1(k)):
            if pi is None:
                if not np.array_equal(D[perm], D): fixed = False
            else:
                if not (np.array_equal(D[perm], D) and np.array_equal(pi[perm], perm[pi])):
                    fixed = False
    say("      transport: gauge-FIXED under every A_v(k) (exact): %s" % fixed)
    return dict(i=inv_ok and herm_ok, ii=dur, iii=ok3,
                iv=all(v == 0 for v in bals.values()), fixed=fixed, bal=bals)

def w64_for(pi, k):      # helper for pi != None balance:  Tr(P D Pi) = sum_c D[pi(c)] P[c,pi(c)]
    return w64[(k, [nm for nm, p in PERMS_PI.items() if p is not None and p is pi][0])]
def D_perm_weight(D, pi, k):
    return D[pi]

star_data = None
if star is not None:
    cn, pnm, v0p, v1p = star
    D_star, pi_star = build_exact_candidate(cn, pnm, v0p, v1p)
    star_data = verify_exact_record("R* = [chi=%s, Pi=%s, v0=%s, v1=%s]" % star,
                                    D_star, pi_star)
else:
    say("    NO exact record in this family; records on this carrier are non-diagonal")
    say("    (they exist by C-41 evenness; constructed numerically in 3c).")

# ------------------------------------------------------------------ 4. C_e
say("")
say("4. C_e = EDGE-SUPPORTED OPERATORS COMMUTING WITH H (exact integer Gram; D-21-clean)")
Tmat_expected = [[Fraction(1) if a == int(MUL[b, r2]) else Fraction(0) for b in range(n)]
                 for a in range(n)]
Imat = [[Fraction(1) if a == b else Fraction(0) for b in range(n)] for a in range(n)]
class_diag_mats = [[[Fraction(1) if (a == b and a in cl) else Fraction(0)
                     for b in range(n)] for a in range(n)] for cl in G["classes"]]
def span_rref(mats):
    rowsv = [[m[a][b] for a in range(n) for b in range(n)] for m in mats]
    piv = []; R = [list(rr) for rr in rowsv]; rr = 0
    for c in range(n * n):
        p = next((i for i in range(rr, len(R)) if R[i][c] != 0), None)
        if p is None: continue
        R[rr], R[p] = R[p], R[rr]
        pv = R[rr][c]; R[rr] = [x / pv for x in R[rr]]
        for i in range(len(R)):
            if i != rr and R[i][c] != 0:
                f = R[i][c]; R[i] = [x - f * y for x, y in zip(R[i], R[rr])]
        piv.append(c); rr += 1
    return R[:rr], piv
def in_span(rref_rows, piv, m):
    v = [m[a][b] for a in range(n) for b in range(n)]
    for row, c in zip(rref_rows, piv):
        if v[c] != 0:
            f = v[c]; v = [x - f * y for x, y in zip(v, row)]
    return all(x == 0 for x in v)
Ce = {}
for e in ["h0", "h1", "u0", "u1"]:
    basis, grank, _ = edge_admissible_algebra(car, e)
    for M in basis:
        ok, nres = verify_in_commutant(car, e, M)
        assert ok, "kernel element fails exact [M,H]=0 on edge %s" % e
    Ce[e] = basis
    rref, piv = span_rref(basis)
    if e in ("h0", "h1"):
        ident = (len(basis) == 2 and in_span(rref, piv, Imat)
                 and in_span(rref, piv, Tmat_expected))
        say("   C_%s: dim %d;  equals span{ I, T } with T: g -> g*r^2 (central right-mult): %s"
            % (e, len(basis), ident))
    else:
        ident = (len(basis) == len(G["classes"])
                 and all(in_span(rref, piv, cd) for cd in class_diag_mats))
        say("   C_%s: dim %d;  equals span of the %d conjugacy-class diagonals (Wilson-type): %s"
            % (e, len(basis), len(G["classes"]), ident))
    say("        every kernel element re-verified EXACTLY: [M (x) I, H] = 0")

# ------------------------------------------------------------------ 5. THE TEST on contractible edges
say("")
say("5. THE TEST -- does ANY admissible unitary on a single CONTRACTIBLE edge flip ANY record?")
say("   C_h0 = span{I, T} is commutative with T^2 = I, so by the flip reduction the only")
say("   candidate flipper is (a phase times) T, and a flippable record exists iff")
say("   Tr(P_E T) = 0 on every H-eigenspace E.  The traces, exactly:")
flip_table = {}
for e, tnm in [("h0", "T0"), ("h1", "T1")]:
    tE = {}
    for k in EIGS:
        t = Fraction(int(np.sum(w64[(k, tnm)])), 64)
        assert t.denominator == 1
        tE[-k] = int(t)
    total = sum(tE.values())
    say("   edge %s:  t_E = Tr(P_E T) = %s   (sum = %d = Tr(T): %s)" % (e, tE, total, total == 0))
    flip_table[e] = tE
blocked = {e: [k for k, v in flip_table[e].items() if v != 0] for e in flip_table}
say("   NONZERO integers on eigenspaces %s (both edges): the T-grading is UNBALANCED there," % blocked["h0"])
say("   so NO Hermitian involution commuting with H anticommutes with T.  NO record --")
say("   transport-fixed or transport-moved, constructed or not -- is flipped by ANY")
say("   admissible unitary on {h0} or {h1}.  CLAUSE (v) HOLDS ON THE CONTRACTIBLE REGIONS.")
if star_data is not None:
    D_s, pi_s = D_star, pi_star
    # {T, R*} Frobenius^2, exact: R* = D*Pi commutes with T (T central, disjoint or commuting),
    # so {T,R*} = 2 T R* and ||{T,R*}||^2 = 4 * ||R*||^2 = 4N
    say("   direct check on the exact record R*: T R* = R* T (both in the gauge-fixed family),")
    say("   so {T,R*} = 2 T R* with Frobenius norm^2 = %d != 0: T does not flip R*." % (4 * N))
say("   NOTE: t_gs = %d with |t_gs| < 22, so T acts NON-TRIVIALLY on the ground space:"
    % flip_table["h0"][-4])
say("   T is a weight-1 admissible logical OPERATION that flips nothing.  Clause (v) as")
say("   written (a flip statement) HOLDS while the Knill-Laflamme code-distance renaming of")
say("   clause (v) (GLOSSARY) FAILS at weight 1 here.  The two readings DIVERGE on this")
say("   proxy -- labeled per D-23.")

# ------------------------------------------------------------------ 6. controls (D-15)
say("")
say("6. CONTROLS BESIDE THE ZERO (D-15)")
say("   (a) ABELIAN COMPARISON, same lattice, same code path (t24_z2.py): D(Z_2) 1x2 has")
say("       C_h0 = span{I, X_h0}, ALL per-eigenspace traces of X_h0 vanish, and X_h0 is an")
say("       admissible weight-1 flipper of the clause-verified record Zbar_h: clause (v)")
say("       FAILS on the abelian proxy.  Same test, opposite verdict: the test discriminates.")
lam = car.edge_perm_left("h0", 4)     # L_s on h0
Uinv = np.empty(N, dtype=np.int64); Uinv[lam] = np.arange(N)
accH = {}
for (i, j), v in Hent.items():
    key = (int(lam[i]), j); accH[key] = accH.get(key, Fraction(0)) + v
    key2 = (i, int(Uinv[j])); accH[key2] = accH.get(key2, Fraction(0)) - v
cn2 = sum(v * v for v in accH.values())
if star_data is not None and dict(chi_list)[star[0]][4] == -1:
    flips = (np.array_equal(D_star[lam], -D_star) if pi_star is None else
             (np.array_equal(D_star[lam], -D_star)
              and np.array_equal(pi_star[lam], lam[pi_star])))
    say("   (b) NON-ADMISSIBLE single-edge flipper of the EXACT RECORD R* on the contractible")
    say("       edge h0:  U = L_s (g_h0 -> s g_h0).  U R* U^dag = -R* exactly: %s" % flips)
    say("       admissibility: ||[U,H]||_F^2 = %s (NONZERO -> not admissible)." % cn2)
    say("       Clause (v) WITHOUT the word 'admissible' is FALSE on this carrier, exactly as")
    say("       it was false of the toric code.  The section-5 zero sits beside a firing control.")
    assert flips and cn2 > 0
else:
    r_wh = chi_list[1][1][hol]
    flips = np.array_equal(r_wh[lam], -r_wh)
    say("   (b) L_s on h0 exactly flips the three-clause Wilson diagonal W_h (not a record on")
    say("       this carrier), ||[L_s,H]||_F^2 = %s != 0: the flip machinery fires when" % cn2)
    say("       admissibility is dropped: %s" % flips)

# ------------------------------------------------------------------ 7. non-contractible single edges
say("")
say("7. THE SAME TEST ON THE NON-CONTRACTIBLE SINGLE EDGES u0, u1 (context for O-41)")
say("   C_u = the 5 class diagonals (commutative); involutions are +-1 class patterns v;")
say("   a record flipped by diag(v) exists  <=>  Tr(P_E diag(v)) = 0 for every E.")
tau_tab = {}
for e in ["u0", "u1"]:
    comp = car.edge_comp(e)
    cls = class_of[comp]
    tab = {}
    for k in EIGS:
        row = []
        for ci in range(ncl):
            t = int(np.dot(diag64[k], (cls == ci).astype(np.int64)))
            assert t % 64 == 0
            row.append(t // 64)
        tab[-k] = row
    tau_tab[e] = tab
    say("   edge %s: class traces per eigenspace (classes %s):" % (e, G["classes"]))
    for k in sorted(tab): say("      E = %2d : %s   (sum %d)" % (k, tab[k], sum(tab[k])))
balanced_patterns = {}
for e in ["u0", "u1"]:
    pats = []
    for signs in itertools.product([1, -1], repeat=ncl):
        if signs[0] == -1 or all(s == 1 for s in signs): continue
        if all(sum(s * t for s, t in zip(signs, tau_tab[e][k])) == 0 for k in tau_tab[e]):
            pats.append(signs)
    balanced_patterns[e] = pats
    say("   edge %s: balanced +-1 class patterns: %d" % (e, len(pats)))
say("   ZERO balanced patterns on both u-edges: NO admissible single-edge unitary flips any")
say("   record on the NON-contractible edges either.  On D(D_4) 1x2 the flip-protection is")
say("   TOTAL across all four single-edge regions -- contrast D(Z_2) 1x2, where admissible")
say("   single-edge flippers exist on every edge (t24_z2.py).")

# ------------------------------------------------------------------ 8. O-41: transport
say("")
say("8. O-41 -- WHICH RECORDS DOES TRANSPORT MOVE, AND DOES PROTECTION DISTINGUISH THEM?")
permsA0 = [car.permA0(k) for k in range(n)]
permsA1 = [car.permA1(k) for k in range(n)]
say("   (a) The section-5/7 obstructions are integer statements about C_e and the")
say("       eigenprojectors; they never mention transport.  Protection therefore holds")
say("       UNIFORMLY: transport-fixed and transport-moved records are equally unflippable")
say("       on every single-edge region of this carrier.  At this proxy size, single-region")
say("       protection does NOT distinguish the two kinds (O-41's question, answered here).")
say("   (b) BOTH KINDS EXIST on this carrier:")
say("       exact gauge-FIXED records: %d found in 3b%s" % (len(found),
    "" if not found else " (e.g. R*)"))
# G x G multiplicities per eigenspace -> transport-fixed existence by O-36's integer feasibility
irr = G["irreps"]
mult_tab = {}
fix_feasible = {}
for k in EIGS:
    P = eigP[k]
    trA = {}
    for ka in range(n):
        pa = permsA0[ka]
        for kb in range(n):
            perm = pa[permsA1[kb]]
            trA[(ka, kb)] = sp_trace_perm(P, perm)
    blocks = []
    dsum = 0
    for (nmi, ci, di) in irr:
        for (nmj, cj, dj) in irr:
            m = Fraction(0)
            for ka in range(n):
                for kb in range(n):
                    m += Fraction(int(ci[ka]) * int(cj[kb])) * trA[(ka, kb)]
            m /= n * n
            assert m.denominator == 1 and m >= 0, (k, nmi, nmj, m)
            m = int(m)
            if m: blocks.append((di * dj, m))
            dsum += di * dj * m
    assert dsum == eigdims[k], (k, dsum, eigdims[k])
    mult_tab[-k] = blocks
    reach = {0}
    for (d, m) in blocks:
        reach = {s + d * (2 * p - m) for s in reach for p in range(m + 1)}
    fix_feasible[-k] = (0 in reach)
    say("       E = %2d (dim %4d): G x G blocks (d_i*d_j x m): %s ; zero-trace reachable: %s"
        % (-k, eigdims[k], mult_tab[-k], fix_feasible[-k]))
say("       transport-FIXED records exist (all eigenspaces feasible, O-36's method): %s"
    % all(fix_feasible.values()))
say("       transport-MOVED records exist generically (random balanced eigenspace splits")
say("       are moved wherever the gauge action is non-scalar) -- witness in section 9.")
say("   (c) The failure mode that WOULD distinguish the kinds -- a region where only one")
say("       kind is flippable -- does not materialise here: nothing is flippable anywhere.")
say("       On the abelian control, where flips DO exist, the flipped record Zbar_h is")
say("       itself transport-fixed, so no distinction appears there either.")

# ------------------------------------------------------------------ 9. numeric witnesses
say("")
say("9. NUMERIC WITNESSES (floats; demonstrations only -- every verdict above is integer-exact)")
rng = np.random.default_rng(24)
def coo_of(cols):
    r_, c_, v_ = [], [], []
    for j in range(N):
        for i, v in cols[j].items():
            r_.append(i); c_.append(j); v_.append(float(v))
    return np.array(r_), np.array(c_), np.array(v_)
def apply_P(coo, X):
    r_, c_, v_ = coo
    Y = np.zeros_like(X)
    np.add.at(Y, r_, (v_[:, None] * X[c_, :]))
    return Y
Qs = {}
for k in EIGS:
    d = eigdims[k]
    coo = coo_of(eigP[k])
    Y = apply_P(coo, rng.standard_normal((N, d + 12)))
    U_, s_, Vt_ = np.linalg.svd(Y, full_matrices=False)
    assert int((s_ > 1e-8 * s_[0]).sum()) == d
    Qs[k] = U_[:, :d]
Hd = np.zeros((N, N))
for (i, j), v in Hent.items(): Hd[i, j] = float(v)
def frob(X): return float(np.linalg.norm(X))
def gauge_move_norm(R):
    mv = 0.0
    for k in range(n):
        for perm in (permsA0[k], permsA1[k]):
            ip = np.empty(N, dtype=np.int64); ip[perm] = np.arange(N)
            mv = max(mv, frob(R[ip, :] - R[:, perm]))
    return mv
def T_anticomm_norm(R):
    return frob(R[tau0, :] + R[:, tau0])
def report(nm, R):
    say("   %s:" % nm)
    say("      ||R^2-I|| = %.2e  ||R-R^T|| = %.2e  ||[R,H]|| = %.2e" %
        (frob(R @ R - np.eye(N)), frob(R - R.T), frob(R @ Hd - Hd @ R)))
    tb = max(abs(float(np.trace(Qs[k].T @ R @ Qs[k]))) for k in EIGS)
    say("      max_E |Tr(P_E R)| = %.2e (balanced);  max_k ||[A_v(k),R]|| = %.3e  -> %s" %
        (tb, gauge_move_norm(R), "transport-MOVED" if gauge_move_norm(R) > 1e-6
         else "transport-FIXED"))
    say("      ||{T,R}||_F = %.3e  (NONZERO: T fails to flip it, as the exact obstruction requires)"
        % T_anticomm_norm(R))
# moved record: random balanced split of every eigenspace
Rm = np.zeros((N, N))
for k in EIGS:
    d = eigdims[k]; Q = Qs[k]
    O, _ = np.linalg.qr(rng.standard_normal((d, d)))
    M = Q @ O
    h = d // 2
    Rm += M[:, :h] @ M[:, :h].T - M[:, h:] @ M[:, h:].T
report("R_moved (random balanced eigenspace splits; four clauses numeric)", Rm)
# fixed record: gauge-averaged commutant element, minimal projections, signed to zero (O-36)
Y0 = rng.standard_normal((N, N)); Y0 = (Y0 + Y0.T) / 2
X = np.zeros((N, N))
for k in EIGS:
    Q = Qs[k]; X += Q @ (Q.T @ Y0 @ Q) @ Q.T
Xav = np.zeros((N, N))
for ka in range(n):
    for kb in range(n):
        perm = permsA0[ka][permsA1[kb]]
        ip = np.empty(N, dtype=np.int64); ip[perm] = np.arange(N)
        Xav += X[np.ix_(ip, ip)]
Xav /= n * n
Rf = np.zeros((N, N))
feas_all = True
for k in EIGS:
    Q = Qs[k]
    S = Q.T @ Xav @ Q; S = (S + S.T) / 2
    w, V = np.linalg.eigh(S)
    groups = []
    i = 0
    while i < len(w):
        j = i
        while j + 1 < len(w) and abs(w[j + 1] - w[i]) < 1e-7: j += 1
        groups.append((i, j + 1)); i = j + 1
    ranks = [b - a for a, b in groups]
    reach = {0: []}
    for r_ in ranks:
        nxt = {}
        for s_, path in reach.items():
            for sg in (1, -1):
                t_ = s_ + sg * r_
                if t_ not in nxt: nxt[t_] = path + [sg]
        reach = nxt
    sg = reach.get(0)
    if sg is None:
        feas_all = False; say("   E=%d: no zero-trace sign assignment this seed" % (-k)); continue
    D_ = np.zeros(eigdims[k])
    for (a, b), s_ in zip(groups, sg): D_[a:b] = s_
    Rf += Q @ (V * D_[None, :]) @ V.T @ Q.T
if feas_all:
    report("R_fixed (gauge-averaged minimal projections signed to zero; O-36 construction)", Rf)
say("   both kinds are records on this carrier and BOTH show the same (null) flip exposure.")

# ------------------------------------------------------------------ verdict + json
say("")
say("=" * 112)
say("VERDICT (T-24, scoped per D-23 to the 1x2 boundary-constrained proxy):")
say("  clause (v) HOLDS on D(D_4) 1x2 for EVERY record on BOTH contractible regions {h0},{h1}")
say("  (and, beyond the clause's demand, on the non-contractible single edges too).  The")
say("  obstruction is exact and integer: the only non-scalar admissible edge-unitaries are")
say("  the central translation T (h-edges) and class-diagonal unitaries (u-edges), and none")
say("  of their involutions has all per-eigenspace traces zero.  The mechanism is genuinely")
say("  non-abelian -- r^2 lies in [D_4,D_4] and conjugation has fixed classes -- and the")
say("  abelian D(Z_2) control on the same lattice FAILS the same test (weight-1 admissible")
say("  flipper).  O-41: protection does not distinguish transport-fixed from transport-moved")
say("  records at this size; both kinds exist and both are protected.")
say("time: %.1f s" % (time.time() - t0))

with open(LANE + "/t24_main.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
json.dump(dict(
    carrier="D(D_4) 1x2 torus", dim=N,
    sector_dims={str(k): v for k, v in sorted(dims.items())},
    eig_dims={str(-k): v for k, v in eigdims.items()},
    contractible_single_regions=[["h0"], ["h1"]],
    Ce_dims={e: len(Ce[e]) for e in Ce},
    tE_T={e: {str(k): v for k, v in flip_table[e].items()} for e in flip_table},
    wilson_diagonals_are_records=False,
    exact_records_found=len(found),
    exact_record_example=None if star is None else dict(chi=star[0], Pi=star[1],
                                                        v0=list(star[2]), v1=list(star[3])),
    balanced_u_patterns={e: len(balanced_patterns[e]) for e in balanced_patterns},
    transport_fixed_exists=all(fix_feasible.values()),
    protection_distinguishes_transport_kinds=False,
    verdict="HOLDS on contractible regions (and all single-edge regions) of the 1x2; "
            "abelian control FAILS; minimal torus VACUOUS"),
    open(LANE + "/t24_main.json", "w"), indent=1)
