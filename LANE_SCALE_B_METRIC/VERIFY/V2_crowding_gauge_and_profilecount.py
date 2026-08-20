"""ADVERSARIAL CHECK V2 -- THE CONTROL THE LANE RAN ON (d) AND (e) BUT NOT ON (b2).

The lane's step 2D gauge test covers relations (a),(c),(d),(e) ONLY.  The crowding relation
(b2) -- which the lane calls "the one operational relation with content" and whose RANGE it
reports as the study's headline SATURATION -- was never gauge-tested and never basis-tested.

But the lane's own step 8 establishes that C_ij is an exact function of
(anticommutation bit, support profile_i, support profile_j), and the support profile is read
off the PHYSICAL SUPPORT, which the lane's own step 2D shows changes by the MAXIMUM POSSIBLE
AMOUNT (1.0000) under stabiliser multiplication.  So (b2) should be maximally gauge-dependent.

TWO CHECKS, both run:
  V2a  recompute the crowding matrix under 10 random stabiliser (gauge) draws per n.
       Report max|C_gauge - C_base|, Cmin, Cmax, asym, d90.  Controls in the same table.
  V2b  is the crowding d90 just the COUNT OF DISTINCT SUPPORT PROFILES?  If rank and d90
       track the profile count, the "growth" is the Gram-Schmidt weight staircase again.
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
Hd = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
NQH = 3
env6 = Environment(nq=6, energies=(1.0, 1.4, 0.7, 1.1, 0.9, 1.3), beta=2.0)

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

def chi_pair_fast(anti, wi, wj, crowded, lam, env, nqh, times=TIMES):
    """VERBATIM the lane's own sector routine (07_crowding_relation.py)."""
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
    C = U.conj().T @ np.kron(I2 / 2, rth) @ U
    ai = aj = 0.0
    nvalid = 0
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
            nvalid += 1
            av = sum(p * r for p, r in outs)
            v = max(vn(av) - sum(p * vn(r) for p, r in outs), 0.0)
            if which == 0: ai += v
            else: aj += v
    return ai / len(times), aj / len(times)

CACHE = {}
def C_entry(anti, wi, wj, lam):
    key = (anti, tuple(np.round(wi, 9)), tuple(np.round(wj, 9)), lam)
    if key in CACHE: return CACHE[key]
    ci_c, _ = chi_pair_fast(anti, wi, wj, True, lam, env6, NQH)
    ci_s, _ = chi_pair_fast(anti, wi, wj, False, lam, env6, NQH)
    v = 1.0 - ci_c / ci_s if abs(ci_s) > 1e-12 else float('nan')
    CACHE[key] = v
    return v

def crowd_matrix(vecs, n, lam=0.8):
    m = len(vecs)
    profs = [profile(sorted(support(v, n)), NQH) for v in vecs]
    M = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            if i == j: continue
            M[i, j] = C_entry(bool(sp_form(vecs[i], vecs[j], n)), profs[i], profs[j], lam)
    return M, profs

def geom(M):
    Ms = (M + M.T) / 2
    D = np.sqrt(np.maximum(M.max() - Ms, 0.0)); np.fill_diagonal(D, 0.0)
    st = dim_stats(double_centre(D))
    return st["d_frac"], st["d_pr"]

P("=" * 116)
P("V2a  GAUGE TEST ON THE CROWDING RELATION (b2) -- the control the lane ran on (d) and (e)")
P("     but NOT on (b2).  10 random stabiliser draws per n, lam = 0.8, nqh = 3, 25 times in [1,13].")
P("     CONTROL IN THE SAME TABLE: relation (a), which the lane already showed is gauge-INVARIANT")
P("     (max-change 0.0000) -- if (b2) also returns 0.0000 the test is dead and proves nothing.")
P("=" * 116)
P("")
P("  %-4s %-4s | %-9s %-9s %-9s %-9s %-9s | %-9s %-9s | %-9s" %
  ("n", "2k", "b2 maxChg", "Cmin base", "Cmin rng", "Cmax rng", "b2 d90 rng", "(a) maxChg", "(a) d90 rng", "#profs"))
P("  " + "-" * 112)
rng = np.random.default_rng(31337)
for n in [4, 6, 8, 10, 12, 14, 16]:
    stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n); m = len(vs)
    Mb, profs = crowd_matrix(vs, n)
    off = ~np.eye(m, dtype=bool)
    d0, _ = geom(Mb)
    Ab = M_symplectic(vs, n)
    chg = 0.0; achg = 0.0
    cmins = [Mb[off].min()]; cmaxs = [Mb[off].max()]; d90s = [d0]; ad90s = []
    Da = 1.0 - Ab; np.fill_diagonal(Da, 0.0)
    ad90s.append(dim_stats(double_centre(Da))["d_frac"])
    S = stab_group(n)
    for _ in range(10):
        gvs = [pauli_mul(v, S[rng.integers(0, 4)], n) for v in vs]
        Mg, _ = crowd_matrix(gvs, n)
        chg = max(chg, float(np.max(np.abs(Mg - Mb))))
        cmins.append(Mg[off].min()); cmaxs.append(Mg[off].max())
        d90s.append(geom(Mg)[0])
        Ag = M_symplectic(gvs, n)
        achg = max(achg, float(np.max(np.abs(Ag - Ab))))
        Dg = 1.0 - Ag; np.fill_diagonal(Dg, 0.0)
        ad90s.append(dim_stats(double_centre(Dg))["d_frac"])
    nprof = len(set(tuple(np.round(p, 9)) for p in profs))
    P("  %-4d %-4d | %-9.4f %-9.4f %-9s %-9s %-9s | %-9.4f %-9s | %-9d" %
      (n, m, chg, Mb[off].min(), "%.3f-%.3f" % (min(cmins), max(cmins)),
       "%.3f-%.3f" % (min(cmaxs), max(cmaxs)), "%d-%d" % (min(d90s), max(d90s)),
       achg, "%d-%d" % (min(ad90s), max(ad90s)), nprof))

P("")
P("=" * 116)
P("V2b  IS THE CROWDING d90 JUST THE COUNT OF DISTINCT SUPPORT PROFILES?")
P("     If numerical rank and d90 track #profiles, the (b2) 'growth' is the Gram-Schmidt")
P("     weight staircase -- the SAME artifact the lane already diagnosed for (d) and (e).")
P("=" * 116)
P("")
P("  %-4s %-4s | %-8s %-10s %-10s %-8s" % ("n", "2k", "#profs", "rank(C)", "d90(b2)", "2k"))
P("  " + "-" * 60)
xs, ys, zs = [], [], []
for n in [4, 6, 8, 10, 12, 14, 16, 18, 20]:
    stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n); m = len(vs)
    Mb, profs = crowd_matrix(vs, n)
    nprof = len(set(tuple(np.round(p, 9)) for p in profs))
    r = int(np.linalg.matrix_rank(Mb, tol=1e-8))
    d, _ = geom(Mb)
    P("  %-4d %-4d | %-8d %-10d %-10d %-8d" % (n, m, nprof, r, d, m))
    xs.append(nprof); ys.append(d); zs.append(m)
xs = np.array(xs, float); ys = np.array(ys, float); zs = np.array(zs, float)
P("")
P("  corr(d90, #distinct profiles) = %.4f     slope d90 vs #profs = %+.4f" %
  (np.corrcoef(xs, ys)[0, 1], np.polyfit(xs, ys, 1)[0]))
P("  corr(d90, 2k)                 = %.4f     slope d90 vs 2k     = %+.4f" %
  (np.corrcoef(zs, ys)[0, 1], np.polyfit(zs, ys, 1)[0]))
P("  #distinct profiles vs 2k slope = %+.4f" % np.polyfit(zs, xs, 1)[0])

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/VERIFY/V2_crowding_gauge_and_profilecount.txt",
     "w").write("\n".join(OUT) + "\n")
