"""STEP 7 -- THE CROWDING RELATION MATRIX, and whether IT has geometric structure.

DEFINITION.  Two records at a time.  Each couples to the bath with the SAME total strength and
with a SITE PROFILE read off its physical support (physical qubit q -> bath site q mod nqh).
    CROWDED : both records couple into bath region A (sites 0..nqh-1)
    SPREAD  : record i into region A, record j into region B (sites nqh..2nqh-1) --
              identical profile shapes, identical strengths, disjoint sites
    C_ij = 1 - <chi_i>_crowded / <chi_i>_spread          (25 times in [1,13])
D-16 IS OBSERVED: the denominator is the SPREAD control, never an "alone" value.

AN EXACT SECTOR DECOMPOSITION MAKES THIS CHEAP.  From a maximally mixed code state and with
only two record-couplings, the system content is the algebra the two records generate:
4-dimensional if they commute (four bath sectors), 2-dimensional if they anticommute.
SC-12 validates it against the full 2^k code space, which step 4 tied to the full 2^n space.

CONTROLS IN THE SAME TABLE
  FREE   k unentangled qubits, H = 0, records = the 2k single-qubit Paulis (support {q}).
  RAND   random symmetric matrix, same size, same value range.
D-17    lam in {0.4,0.8,1.2} and region size nqh in {2,3} are both swept.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC")
from lib_scaleb import *
from lib_operational import *
import numpy as np

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

TIMES = np.linspace(1.0, 13.0, 25)
rng = np.random.default_rng(5150)
Hd = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)

def vn(r):
    e = np.linalg.eigvalsh((r + r.conj().T) / 2); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

def profile(supp, nqh):
    w = np.zeros(nqh)
    for q in supp: w[q % nqh] += 1
    return w / w.sum()

def bath_ops(env, wi, wj, nqh, crowded):
    Bi = sum(wi[s] * env.site[s] for s in range(nqh))
    Bj = sum(wj[s] * env.site[(s if crowded else nqh + s)] for s in range(nqh))
    return Bi, Bj

# ---------------------------------------------------------------- exact sector chi
def chi_pair_fast(anti, wi, wj, crowded, lam, env, nqh, times=TIMES):
    Bi, Bj = bath_ops(env, wi, wj, nqh, crowded)
    rth = env.thermal(); nB = env.dim
    if not anti:
        sect = {}
        for r in (1, -1):
            for s in (1, -1):
                w, U = np.linalg.eigh(env.HB + lam * (r * Bi + s * Bj))
                sect[(r, s)] = (w, U, U.conj().T @ rth @ U)
        ai = aj = 0.0
        for t in times:
            rb = {}
            for key, (w, U, C) in sect.items():
                ph = np.exp(-1j * w * t)
                rb[key] = U @ (ph[:, None] * C * ph.conj()[None, :]) @ U.conj().T
            avg = sum(rb.values()) / 4.0
            Sav = vn(avg)
            ai += Sav - 0.5 * (vn((rb[(1, 1)] + rb[(1, -1)]) / 2) + vn((rb[(-1, 1)] + rb[(-1, -1)]) / 2))
            aj += Sav - 0.5 * (vn((rb[(1, 1)] + rb[(-1, 1)]) / 2) + vn((rb[(1, -1)] + rb[(-1, -1)]) / 2))
        return max(ai / len(times), 0.0), max(aj / len(times), 0.0)
    Ht = np.kron(I2, env.HB) + lam * (np.kron(Zm, Bi) + np.kron(Xm, Bj))
    w, U = np.linalg.eigh(Ht)
    r0 = np.kron(I2 / 2, rth)
    C = U.conj().T @ r0 @ U
    ai = aj = 0.0
    for t in times:
        ph = np.exp(-1j * w * t)
        rho = U @ (ph[:, None] * C * ph.conj()[None, :]) @ U.conj().T
        for which in (0, 1):
            rr = rho
            if which == 1:
                T = rr.reshape(2, nB, 2, nB)
                T = np.einsum('ab,bicj,cd->aidj', Hd, T, Hd.conj().T)
                rr = T.reshape(2 * nB, 2 * nB)
            T = rr.reshape(2, nB, 2, nB)
            outs = []
            for s in (0, 1):
                rB = T[s, :, s, :]
                p = float(np.real(np.trace(rB)))
                if p > 1e-12: outs.append((p, rB / p))
            if len(outs) < 2: continue
            av = sum(p * r for p, r in outs)
            v = max(vn(av) - sum(p * vn(r) for p, r in outs), 0.0)
            if which == 0: ai += v
            else: aj += v
    return ai / len(times), aj / len(times)

def chi_pair_full(ri, rj, wi, wj, crowded, lam, env, nqh, Hs, st0):
    """the same quantity computed in the FULL code space, for validation only"""
    nS = ri.shape[0]; nB = env.dim
    Bi, Bj = bath_ops(env, wi, wj, nqh, crowded)
    HINT = np.kron(ri, Bi) + np.kron(rj, Bj)
    ai = aj = 0.0
    for rho in evolve_cached(Hs, env, HINT, lam, TIMES, st0):
        for R, which in ((ri, 0), (rj, 1)):
            outs = []
            for s in (+1, -1):
                Pr = np.kron((np.eye(nS) + s * R) / 2, np.eye(nB))
                blk = Pr @ rho @ Pr
                p = float(np.real(np.trace(blk)))
                if p < 1e-12: continue
                outs.append((p, (blk / p).reshape(nS, nB, nS, nB).trace(axis1=0, axis2=2)))
            if len(outs) < 2: continue
            av = sum(p * r for p, r in outs)
            v = max(vn(av) - sum(p * vn(r) for p, r in outs), 0.0)
            if which == 0: ai += v
            else: aj += v
    return ai / len(TIMES), aj / len(TIMES)

NQH = 3
env6 = Environment(nq=6, energies=(1.0, 1.4, 0.7, 1.1, 0.9, 1.3), beta=2.0)
env4 = Environment(nq=4, energies=(1.0, 1.4, 0.7, 1.1), beta=2.0)

P("=" * 118)
P("LANE_SCALE_B_METRIC  STEP 7 -- CROWDING RELATION MATRIX")
P("=" * 118)
P("")
P("SC-12  exact sector decomposition vs the FULL 2^k code space (n = 6, lam = 0.8, nqh = 3)")
P("  %-10s %-7s %-8s %-15s %-15s %-11s %-8s" %
  ("pair", "anti", "mode", "chi_i FULL", "chi_i SECTOR", "|diff|", "verdict"))
P("  " + "-" * 84)
n = 6
stab, pairs = carrier(n); k = len(pairs); vs, lab = record_vectors(pairs, n)
Heff = -2 * np.eye(2 ** k, dtype=complex); st0k = np.eye(2 ** k, dtype=complex) / 2 ** k
Rred = [std_pauli(k, i % k, 'X' if i < k else 'Z') for i in range(2 * k)]
worst = 0.0
for (i, j) in [(0, 1), (0, k), (1, k + 2), (2, 3)]:
    anti = bool(sp_form(vs[i], vs[j], n))
    wi = profile(sorted(support(vs[i], n)), NQH); wj = profile(sorted(support(vs[j], n)), NQH)
    for crowded in (True, False):
        fi, fj = chi_pair_full(Rred[i], Rred[j], wi, wj, crowded, 0.8, env6, NQH, Heff, st0k)
        si, sj = chi_pair_fast(anti, wi, wj, crowded, 0.8, env6, NQH)
        dd = max(abs(fi - si), abs(fj - sj)); worst = max(worst, dd)
        P("  %-10s %-7s %-8s %-15.9f %-15.9f %-11.2e %-8s" %
          (lab[i] + "," + lab[j], anti, "crowd" if crowded else "spread", fi, si, dd,
           "PASS" if dd < 1e-9 else "FAIL"))
P("")
P("  SC-12 worst |diff| = %.3e  %s" % (worst, "PASS" if worst < 1e-9 else "FAIL"))
if worst >= 1e-9:
    P("  SETUP BROKEN -- no crowding result is reported.")
    open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/07_crowding_relation.txt", "w").write("\n".join(OUT) + "\n")
    sys.exit(0)

# ---------------------------------------------------------------- the matrix
CACHE = {}
def cpair(anti, wi, wj, crowded, lam, env, nqh):
    key = (anti, tuple(np.round(wi, 9)), tuple(np.round(wj, 9)), crowded, lam, id(env), nqh)
    if key not in CACHE:
        CACHE[key] = chi_pair_fast(anti, wi, wj, crowded, lam, env, nqh)
    return CACHE[key]

def crowding_matrix(supps, antif, m, nqh, lam, env):
    C = np.full((m, m), np.nan)
    profs = [profile(s, nqh) for s in supps]
    for i in range(m):
        for j in range(i + 1, m):
            a = antif(i, j)
            ci_c, cj_c = cpair(a, profs[i], profs[j], True, lam, env, nqh)
            ci_s, cj_s = cpair(a, profs[i], profs[j], False, lam, env, nqh)
            C[i, j] = 1.0 - ci_c / ci_s if ci_s > 1e-12 else np.nan
            C[j, i] = 1.0 - cj_c / cj_s if cj_s > 1e-12 else np.nan
    np.fill_diagonal(C, 0.0)
    return C

def fast_dc(D):
    D2 = D ** 2
    return -0.5 * (D2 - D2.mean(1, keepdims=True) - D2.mean(0, keepdims=True) + D2.mean())

def dimsD(D):
    B = fast_dc(D); w = np.linalg.eigvalsh((B + B.T) / 2)[::-1]
    pos = w[w > 1e-10]
    if not len(pos): return 0, 0.0, 0.0, w
    c = np.cumsum(pos) / pos.sum()
    neg = float(-w[w < -1e-10].sum())
    return int(np.searchsorted(c, 0.90) + 1), float(pos.sum() ** 2 / (pos ** 2).sum()), \
           neg / (neg + pos.sum()), w

def geom(C):
    m = C.shape[0]
    Ms = (C + C.T) / 2
    off = ~np.eye(m, dtype=bool)
    D = np.sqrt(np.maximum(Ms.max() - Ms, 0.0)); np.fill_diagonal(D, 0.0)
    d90, dpr, negf, w = dimsD(D)
    nv, wv = triangle_violations(D) if m <= 40 else (-1, -1.0)
    return dict(cmin=C[off].min(), cmax=C[off].max(), asym=float(np.abs(C - C.T).max()),
                d90=d90, dpr=dpr, negf=negf, tri=nv, cdim=corr_dim(D), evals=w)

P("")
P("=" * 118)
P("7B.  THE CROWDING MATRIX AND ITS GEOMETRY.  lam = 0.8, region size nqh = 3.")
P("     tri = ordered triples violating the triangle inequality after d = sqrt(Cmax - C_sym).")
P("")
P("  %-4s %-4s | %-40s | %-40s | %-13s" %
  ("n", "2k", "CODE   Cmin    Cmax    asym   d90 dPR  tri cdim",
   "FREE   Cmin    Cmax    asym   d90 dPR  tri cdim", "RAND d90 dPR"))
P("  " + "-" * 112)
NS7 = [4, 6, 8, 10, 12, 14, 16, 18, 20]
store = {}
for n in NS7:
    stab, pairs = carrier(n); k = len(pairs); vs, lab = record_vectors(pairs, n); m = 2 * k
    supps = [sorted(support(v, n)) for v in vs]
    C = crowding_matrix(supps, lambda i, j: bool(sp_form(vs[i], vs[j], n)), m, NQH, 0.8, env6)
    F = crowding_matrix([[i % k] for i in range(m)], lambda i, j: abs(i - j) == k, m, NQH, 0.8, env6)
    off = ~np.eye(m, dtype=bool)
    lo, hi = C[off].min(), C[off].max()
    Mr = random_control(m, 1.0, rng, vals=None) * (hi - lo) + lo
    Mr = (Mr + Mr.T) / 2; np.fill_diagonal(Mr, 0.0)
    Dr = np.sqrt(np.maximum(Mr.max() - Mr, 0.0)); np.fill_diagonal(Dr, 0.0)
    d90r, dprr, _, _ = dimsD(Dr)
    gc, gf = geom(C), geom(F)
    store[n] = (C, F, gc, gf)
    P("  %-4d %-4d | %7.4f %7.4f %7.1e %3d %5.2f %4d %5.2f | %7.4f %7.4f %7.1e %3d %5.2f %4d %5.2f | %3d %5.2f" %
      (n, m, gc["cmin"], gc["cmax"], gc["asym"], gc["d90"], gc["dpr"], gc["tri"], gc["cdim"],
       gf["cmin"], gf["cmax"], gf["asym"], gf["d90"], gf["dpr"], gf["tri"], gf["cdim"],
       d90r, dprr))

P("")
P("  Eigenvalue spectrum of the CODE crowding matrix (top 8) at each n:")
for n in NS7:
    C = store[n][0]
    w = np.linalg.eigvalsh((C + C.T) / 2)[::-1]
    P("   n=%-3d " % n + "[" + " ".join("%+.4f" % x for x in w[:8]) + (" ...]" if len(w) > 8 else "]"))

P("")
P("  The n = 6 CODE crowding matrix in full:")
n = 6
stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n)
C = store[6][0]
P("        " + " ".join("%8s" % l for l in lab))
for i, l in enumerate(lab):
    P("  %-6s" % l + " ".join("%8.4f" % C[i, j] for j in range(len(lab))))

# ---------------------------------------------------------------- D-17
P("")
P("=" * 118)
P("7C.  D-17 -- VARY THE VENUE'S OWN SCALE before calling anything new.")
P("  %-4s %-6s %-6s | %-9s %-9s %-6s %-7s %-7s" %
  ("n", "lam", "nqh", "Cmin", "Cmax", "d90", "dPR", "cdim"))
P("  " + "-" * 66)
for n in [8, 12, 16]:
    stab, pairs = carrier(n); k = len(pairs); vs, lab = record_vectors(pairs, n); m = 2 * k
    supps = [sorted(support(v, n)) for v in vs]
    af = lambda i, j: bool(sp_form(vs[i], vs[j], n))
    for lam in [0.4, 0.8, 1.2]:
        for nqh, ev in [(2, env4), (3, env6)]:
            C = crowding_matrix(supps, af, m, nqh, lam, ev)
            g = geom(C)
            P("  %-4d %-6.1f %-6d | %-9.5f %-9.5f %-6d %-7.2f %-7.2f" %
              (n, lam, nqh, g["cmin"], g["cmax"], g["d90"], g["dpr"], g["cdim"]))

# ---------------------------------------------------------------- what carries the relation
P("")
P("=" * 118)
P("7D.  WHAT DOES C_ij ACTUALLY DEPEND ON?  Least squares of C_ij on (1) the anticommutation")
P("     bit alone, (2) the site-profile inner product w_i . w_j alone, (3) both.")
P("     R^2 = 1 on 'both' means the operational relation carries NO information beyond the")
P("     symplectic bit and the support profile -- i.e. nothing that (a) and (d) do not already have.")
P("")
P("  %-4s %-4s %-16s %-16s %-16s" % ("n", "2k", "R2 anti only", "R2 profile only", "R2 both"))
P("  " + "-" * 62)
for n in NS7:
    stab, pairs = carrier(n); k = len(pairs); vs, lab = record_vectors(pairs, n); m = 2 * k
    supps = [sorted(support(v, n)) for v in vs]
    profs = [profile(s, NQH) for s in supps]
    C = store[n][0]
    y, xa, xp = [], [], []
    for i in range(m):
        for j in range(m):
            if i == j: continue
            y.append(C[i, j]); xa.append(float(sp_form(vs[i], vs[j], n)))
            xp.append(float(np.dot(profs[i], profs[j])))
    y = np.array(y)
    def r2(cols):
        A = np.column_stack([np.ones(len(y))] + cols)
        b, *_ = np.linalg.lstsq(A, y, rcond=None)
        return 1.0 - (y - A @ b).var() / y.var() if y.var() > 0 else 1.0
    P("  %-4d %-4d %-16.4f %-16.4f %-16.4f" %
      (n, m, r2([np.array(xa)]), r2([np.array(xp)]), r2([np.array(xa), np.array(xp)])))

P("")
P("7E.  DRIFT.  d90 of the crowding geometry against 2k.")
P("  %-6s %-6s %-8s %-8s %-10s" % ("n", "2k", "CODE d90", "FREE d90", "CODE/(2k-1)"))
P("  " + "-" * 44)
xs, ys = [], []
for n in NS7:
    m = store[n][0].shape[0]
    P("  %-6d %-6d %-8d %-8d %-10.3f" % (n, m, store[n][2]["d90"], store[n][3]["d90"],
                                         store[n][2]["d90"] / (m - 1)))
    xs.append(m); ys.append(store[n][2]["d90"])
sl = np.polyfit(np.array(xs, float), np.array(ys, float), 1)[0]
P("")
P("  d90 vs 2k slope = %+.4f  -> %s" %
  (sl, "SATURATES" if abs(sl) < 0.05 else "GROWS WITH n -- no fixed intrinsic dimension"))

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/07_crowding_relation.txt", "w").write("\n".join(OUT) + "\n")
