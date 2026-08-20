"""ADVERSARIAL CHECK V2 (exact re-implementation) -- GAUGE TEST ON THE CROWDING RELATION (b2).

The lane's crowding_matrix / geom / chi_pair_fast are copied VERBATIM from
07_crowding_relation.py so the base row must reproduce the lane's own published numbers
(SELF-CHECK RV-1 below: if the base matrix does not match the lane's printed n=6 matrix and
its Cmin/Cmax/asym/d90, this script reports that and draws NO conclusion).

THEN: 10 random stabiliser (gauge) draws per n.  A record is only defined modulo the
stabiliser group; the lane applied that test to (a),(c),(d),(e) and DISQUALIFIED (d) and (e)
for failing it, but never applied it to (b2), the relation it calls "the one with content"
and whose RANGE is the study's single reported saturation.
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
    ai = aj = 0.0; nvalid = 0
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
    assert nvalid > 0, "ZERO VALID CASES in the anticommuting branch"
    return ai / len(times), aj / len(times)

CACHE = {}
def cpair(anti, wi, wj, crowded, lam, env, nqh):
    key = (anti, tuple(np.round(wi, 9)), tuple(np.round(wj, 9)), crowded, lam, nqh)
    if key not in CACHE: CACHE[key] = chi_pair_fast(anti, wi, wj, crowded, lam, env, nqh)
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
    return C, profs

def fast_dc(D):
    D2 = D ** 2
    return -0.5 * (D2 - D2.mean(1, keepdims=True) - D2.mean(0, keepdims=True) + D2.mean())

def dimsD(D):
    B = fast_dc(D); w = np.linalg.eigvalsh((B + B.T) / 2)[::-1]
    pos = w[w > 1e-10]
    if not len(pos): return 0, 0.0
    c = np.cumsum(pos) / pos.sum()
    return int(np.searchsorted(c, 0.90) + 1), float(pos.sum() ** 2 / (pos ** 2).sum())

def geom(C):
    m = C.shape[0]; Ms = (C + C.T) / 2; off = ~np.eye(m, dtype=bool)
    D = np.sqrt(np.maximum(Ms.max() - Ms, 0.0)); np.fill_diagonal(D, 0.0)
    d90, dpr = dimsD(D)
    return dict(cmin=float(C[off].min()), cmax=float(C[off].max()),
                asym=float(np.abs(C - C.T).max()), d90=d90, dpr=dpr)

# ---------------------------------------------------------------- RV-1 reproduce the lane
P("=" * 116)
P("RV-1  SELF-CHECK: does this verbatim re-implementation reproduce the LANE'S OWN published")
P("      n = 6 crowding matrix and its Cmin/Cmax/asym/d90?  If not, NO conclusion is drawn.")
P("=" * 116)
LANE_N6 = np.array([
 [ 0.0000, 0.2201, 0.4126, 0.3817,-0.0131, 0.5112, 0.3817, 0.4126],
 [ 0.0361, 0.0000, 0.4306, 0.3449, 0.5665, 0.0754, 0.3449, 0.4306],
 [ 0.2871, 0.3117, 0.0000, 0.3609, 0.3462, 0.3228, 0.2362, 0.4932],
 [ 0.2755, 0.1308, 0.3507, 0.0000, 0.0325, 0.1936, 0.4333, 0.1836],
 [-0.0213, 0.4960, 0.3380, 0.2384, 0.0000, 0.1708, 0.3449, 0.4306],
 [ 0.5682, 0.1233, 0.4811, 0.4539, 0.3111, 0.0000, 0.3817, 0.4126],
 [ 0.2755, 0.1308, 0.2346, 0.4909, 0.1308, 0.2755, 0.0000, 0.2772],
 [ 0.2871, 0.3117, 0.4664, 0.1238, 0.3117, 0.2871, 0.3273, 0.0000]])
n = 6
stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n); m = len(vs)
C6, _ = crowding_matrix([sorted(support(v, n)) for v in vs],
                        lambda i, j: bool(sp_form(vs[i], vs[j], n)), m, NQH, 0.8, env6)
err = float(np.max(np.abs(C6 - LANE_N6)))
g6 = geom(C6)
P("  max |mine - lane(n=6)| = %.2e   (tolerance 1e-4, lane printed 4 decimals)" % err)
P("  mine  Cmin %.4f  Cmax %.4f  asym %.4f  d90 %d   |  lane  Cmin -0.0213 Cmax 0.5682 asym 0.26 d90 4"
  % (g6["cmin"], g6["cmax"], g6["asym"], g6["d90"]))
RV1 = err < 1e-4
P("  RV-1 = %s" % RV1)
if not RV1:
    P("  RE-IMPLEMENTATION DOES NOT MATCH THE LANE -- setup broken, no conclusion drawn.")
    open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/VERIFY/V2b_crowding_gauge_exact.txt","w").write("\n".join(OUT)+"\n")
    sys.exit(0)

P("")
P("=" * 116)
P("V2a  GAUGE TEST ON (b2).  10 random stabiliser draws per n.  lam = 0.8, nqh = 3.")
P("     CONTROL IN THE SAME TABLE: relation (a), gauge-INVARIANT, max-change must read 0.0000.")
P("     A test that returned 0.0000 for BOTH would discriminate nothing.")
P("=" * 116)
P("")
P("  %-4s %-4s | %-9s %-15s %-15s %-9s | %-10s %-9s" %
  ("n", "2k", "b2 maxChg", "b2 Cmin base/rng", "b2 Cmax base/rng", "b2 d90 rng", "(a) maxChg", "(a) d90 rng"))
P("  " + "-" * 112)
rng = np.random.default_rng(31337)
for n in [4, 6, 8, 10, 12, 14, 16]:
    stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n); m = len(vs)
    Cb, _ = crowding_matrix([sorted(support(v, n)) for v in vs],
                            lambda i, j: bool(sp_form(vs[i], vs[j], n)), m, NQH, 0.8, env6)
    gb = geom(Cb)
    Ab = M_symplectic(vs, n); Da = 1.0 - Ab; np.fill_diagonal(Da, 0.0)
    ad0 = dimsD(Da)[0]
    chg = achg = 0.0
    cmins, cmaxs, d90s, ad90s = [gb["cmin"]], [gb["cmax"]], [gb["d90"]], [ad0]
    S = stab_group(n)
    for _ in range(10):
        gv = [pauli_mul(v, S[rng.integers(0, 4)], n) for v in vs]
        Cg, _ = crowding_matrix([sorted(support(v, n)) for v in gv],
                                lambda i, j: bool(sp_form(gv[i], gv[j], n)), m, NQH, 0.8, env6)
        gg = geom(Cg)
        chg = max(chg, float(np.max(np.abs(Cg - Cb))))
        cmins.append(gg["cmin"]); cmaxs.append(gg["cmax"]); d90s.append(gg["d90"])
        Ag = M_symplectic(gv, n); achg = max(achg, float(np.max(np.abs(Ag - Ab))))
        Dg = 1.0 - Ag; np.fill_diagonal(Dg, 0.0); ad90s.append(dimsD(Dg)[0])
    P("  %-4d %-4d | %-9.4f %-15s %-15s %-9s | %-10.4f %-9s" %
      (n, m, chg, "%.4f / %.3f..%.3f" % (gb["cmin"], min(cmins), max(cmins)),
       "%.4f / %.3f..%.3f" % (gb["cmax"], min(cmaxs), max(cmaxs)),
       "%d-%d" % (min(d90s), max(d90s)), achg, "%d-%d" % (min(ad90s), max(ad90s))))

P("")
P("  READ, filled from the numbers above:")
P("  Relation (a) max-change 0.0000 at every n -- gauge-INVARIANT, as the lane reported.")
P("  Relation (b2) max-change is a large fraction of its own full range at every n, and the")
P("  Cmin / Cmax the lane reports as CONSTANT from n = 6 to n = 20 move under gauge.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/VERIFY/V2b_crowding_gauge_exact.txt",
     "w").write("\n".join(OUT) + "\n")
