"""S4 -- DOES THE ENVIRONMENT HOLD RELATIONS BETWEEN RECORDS, AND ONLY AT SCALE?

THE OBJECT.  The k commuting records R_1..R_k generate an abelian group of 2^k elements.  EVERY
non-identity element g_S = prod_{i in S} R_i is itself a record on the same carrier (clauses
(i)-(iv) checked below, matrix route at n <= 8 and F_2 route at every n).  So the question
'does the environment hold relations between records, or only the records themselves' is not a
question about derived statistics: it is the question of WHICH ELEMENT OF THE RECORD GROUP the
environment knows best.  Call |S| the DEPTH of g_S.  Depth is well posed here because the
COUPLING singles out R_1..R_k: record i couples to site j with weight W[i,j].

THE MEASUREMENTS, all time-averaged over 25 times in [1,13]:
  chi(g_S : F)            for every non-identity S and every fragment F
  chi_depth1(F)           = max over |S| = 1
  chi_star(F)             = max over all S
  EXCESS(F)               = chi_star - chi_depth1   >= 0
  depth*(F)               = |argmax S|
  SYNERGY(i,j;F)          = chi(R_i,R_j : F) - chi(R_i : F) - chi(R_j : F)
  C3(i,j,l;F)             = chi_ijl - sum chi_pairs + sum chi_singles     (connected 3-body)
  MARGIN(i,j;F)           = chi(R_iR_j : F) - max(chi(R_i:F), chi(R_j:F))

CONTROL IN THE SAME TABLE (D-15).  'separate' -- site j couples to record (j mod k) alone at the
same total coupling per site.  That state is exactly k independent carriers on k disjoint baths
(confirmed through the FULL model in s2), so no fragment can hold a relation.  Every relational
column must read 0 there while the ordinary chi columns stay large.
"""
import sys, io, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY")
from common import *

OUT = io.StringIO()
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.write(s + "\n")

LAM = 0.8
EPS = 0.01          # bits: 'holds essentially nothing'


# ------------------------------------------------------------------ parity matrix / masks
def parity_matrix(k):
    """PV[p, m] = product of r_i over i in the subset m, for the sign pattern p.  Walsh matrix."""
    p = np.arange(2 ** k)[:, None]
    m = np.arange(2 ** k)[None, :]
    pc = np.zeros((2 ** k, 2 ** k), dtype=np.int64)
    x = p & m
    for b in range(k):
        pc += (x >> b) & 1
    return (1 - 2 * (pc % 2)).astype(float)


def mask_rows(PV, specs):
    """specs: list of list of (subset_mask, sign).  mask = prod_x (1 + s_x PV[:,x]) / nP."""
    nP = PV.shape[0]
    rows = np.empty((len(specs), nP))
    for q, sp in enumerate(specs):
        v = np.ones(nP)
        for m, s in sp: v *= (1.0 + s * PV[:, m])
        rows[q] = v / nP
    return rows


def entropy_stack(B, sites, ti, rows):
    rhoF = kron_sites(B.rho, sites, ti)
    return entropies(averages(rhoF, rows))


def bits(m):
    return [b for b in range(20) if (m >> b) & 1]


# ------------------------------------------------------------------ are all group elements records?
P("=" * 132)
P("S4  RELATIONS BETWEEN RECORDS")
P("=" * 132)
P("")
P("PRE-CHECK (D-18): is every non-identity element of the record group itself a record?")
P(f"{'n':>3} {'k':>3} | {'group elts':>10} {'F2: not in <S>':>15} | "
  f"{'M (i)':>6} {'M (ii)':>7} {'M (iii)':>8} {'M (iv)':>7}")
P("-" * 72)
for n in (4, 6, 8):
    car = carrier(n); k = car['k']
    stabspan = set()
    for a in (0, 1):
        for b in (0, 1):
            stabspan.add(tuple((a * car['stab'][0][t] + b * car['stab'][1][t]) % 2 for t in range(2 * n)))
    nf2 = 0
    for m in range(1, 2 ** k):
        v = [0] * (2 * n)
        for i in bits(m):
            v = [(v[t] + car['recs_xz'][i][t]) % 2 for t in range(2 * n)]
        nf2 += (tuple(v) not in stabspan)
    H = code_hamiltonian(n); Rs = record_matrices(car); es = eigenspaces(H)
    ci = cii = ciii = civ = 0
    for m in range(1, 2 ** k):
        G = np.eye(2 ** n, dtype=complex)
        for i in bits(m): G = G @ Rs[i]
        ci += (np.linalg.norm(G - G.conj().T) < 1e-9 and np.linalg.norm(G @ G - np.eye(2 ** n)) < 1e-9)
        cii += (np.linalg.norm(H @ G - G @ H) < 1e-9)
        ciii += clause_iii(G, es); civ += clause_iv(G, es)
    P(f"{n:>3} {k:>3} | {2**k-1:>10} {nf2:>15} | {ci:>6} {cii:>7} {ciii:>8} {civ:>7}")
    del H, Rs
P("READ: all 2^k - 1 non-identity group elements satisfy (i)-(iv).  The 'relations' measured")
P("      below are records, not derived statistics.")

# ------------------------------------------------------------------ TABLE C: pairwise structure
P("")
P("=" * 132)
P("TABLE C  PAIRWISE STRUCTURE, nq = 6, lam = 0.8.  'crd' = crowded, 'sep' = SPREAD CONTROL")
P("          (k independent carriers on disjoint baths -- no relation can exist by construction).")
P("=" * 132)
P(f"{'k':>3} {'frag':>6} | {'chi_i crd':>9} {'chi_i sep':>9} | {'SYN mean crd':>12} {'SYN mean sep':>12} "
  f"{'SYN max crd':>11} {'SYN max sep':>11} | {'chiPAR crd':>10} {'chiPAR sep':>10} | "
  f"{'MARGIN crd':>10} {'MARGIN sep':>10} | {'#pure crd':>9} {'#pure sep':>9}")
P("-" * 168)
NQ = 6
KS = [2, 4, 6, 8, 10]
tabC = {}
for k in KS:
    PV = parity_matrix(k)
    single = [1 << i for i in range(k)]
    prs = list(itertools.combinations(range(k), 2))
    specs = [[]]                                          # uniform
    idx_single = {}
    for i in range(k):
        idx_single[i] = (len(specs), len(specs) + 1)
        specs += [[(single[i], +1)], [(single[i], -1)]]
    idx_joint, idx_par = {}, {}
    for (i, j) in prs:
        idx_joint[(i, j)] = tuple(range(len(specs), len(specs) + 4))
        specs += [[(single[i], a), (single[j], b)] for a in (+1, -1) for b in (+1, -1)]
        idx_par[(i, j)] = (len(specs), len(specs) + 1)
        specs += [[(single[i] | single[j], +1)], [(single[i] | single[j], -1)]]
    rows = mask_rows(PV, specs)
    for kind in ('crowded', 'separate'):
        W = weights(kind, k, NQ)
        B = Broadcast(k, NQ, W, LAM)
        for fname, frags in (('site', [[j] for j in range(NQ)]), ('whole', [list(range(NQ))])):
            acc = None
            for F in frags:
                a = np.zeros(len(specs))
                for ti in range(len(TIMES)):
                    a += entropy_stack(B, F, ti, rows)
                a /= len(TIMES)
                acc = a if acc is None else np.vstack([acc, a]) if acc.ndim > 1 else np.vstack([acc, a])
            acc = np.atleast_2d(acc)
            # per-fragment chi, then average over fragments
            chi_s = np.zeros((acc.shape[0], k))
            syn = np.zeros((acc.shape[0], len(prs)))
            par = np.zeros((acc.shape[0], len(prs)))
            mar = np.zeros((acc.shape[0], len(prs)))
            pure = np.zeros(acc.shape[0])
            for fi in range(acc.shape[0]):
                S = acc[fi]; S0 = S[0]
                for i in range(k):
                    p, q = idx_single[i]; chi_s[fi, i] = max(S0 - 0.5 * (S[p] + S[q]), 0.0)
                for pi_, (i, j) in enumerate(prs):
                    q4 = idx_joint[(i, j)]
                    cij = max(S0 - 0.25 * sum(S[q] for q in q4), 0.0)
                    a2, b2 = idx_par[(i, j)]
                    cp = max(S0 - 0.5 * (S[a2] + S[b2]), 0.0)
                    syn[fi, pi_] = cij - chi_s[fi, i] - chi_s[fi, j]
                    par[fi, pi_] = cp
                    mar[fi, pi_] = cp - max(chi_s[fi, i], chi_s[fi, j])
                    if chi_s[fi, i] < EPS and chi_s[fi, j] < EPS and cp > EPS: pure[fi] += 1
            tabC[(k, kind, fname)] = dict(chi=chi_s, syn=syn, par=par, mar=mar, pure=pure)
    for fname in ('site', 'whole'):
        c = tabC[(k, 'crowded', fname)]; s = tabC[(k, 'separate', fname)]
        P(f"{k:>3} {fname:>6} | {c['chi'].mean():>9.5f} {s['chi'].mean():>9.5f} | "
          f"{c['syn'].mean():>12.6f} {s['syn'].mean():>12.2e} "
          f"{np.abs(c['syn']).max():>11.6f} {np.abs(s['syn']).max():>11.2e} | "
          f"{c['par'].mean():>10.5f} {s['par'].mean():>10.5f} | "
          f"{c['mar'].max():>10.6f} {s['mar'].max():>10.6f} | "
          f"{c['pure'].sum():>9.0f} {s['pure'].sum():>9.0f}")
P("-" * 168)
P("SYN = chi(R_i,R_j:F) - chi(R_i:F) - chi(R_j:F).  chiPAR = chi of the parity record R_iR_j.")
P("MARGIN = chiPAR - max(chi_i, chi_j): positive means the fragment knows the RELATION better")
P("than either record it relates.  #pure = pairs where chi_i,chi_j < 0.01 bit but chiPAR > 0.01,")
P("summed over fragments -- the fragment holds the relation and neither relatum.")

# ------------------------------------------------------------------ TABLE D: the group scan
P("")
P("=" * 132)
P("TABLE D  THE GROUP SCAN.  chi of EVERY one of the 2^k - 1 records in the group, by DEPTH.")
P("          nq = 6, lam = 0.8.  Control column 'sep' in the same table.")
P("=" * 132)


def group_scan(k, nq, kind, lam, frags, seed=7, times=TIMES):
    PV = parity_matrix(k)
    specs = [[]]
    for m in range(1, 2 ** k):
        specs += [[(m, +1)], [(m, -1)]]
    rows = mask_rows(PV, specs)
    W = weights(kind, k, nq, seed=seed)
    B = Broadcast(k, nq, W, lam, times=times)
    depth = np.array([bin(m).count('1') for m in range(2 ** k)])
    out = []
    for F in frags:
        a = np.zeros(len(specs))
        for ti in range(len(times)):
            a += entropy_stack(B, F, ti, rows)
        a /= len(times)
        chi = np.zeros(2 ** k)
        for m in range(1, 2 ** k):
            chi[m] = max(a[0] - 0.5 * (a[2 * m - 1] + a[2 * m]), 0.0)
        out.append(chi)
    return np.array(out), depth


P(f"{'k':>3} {'geom':>9} {'frag':>6} | {'chi_depth1':>10} {'chi_star':>9} {'EXCESS':>9} "
  f"{'depth*':>7} {'frac depth*>1':>13} | " + " ".join(f"{'d'+str(d):>8}" for d in range(1, 7)))
P("-" * 150)
GD = {}
for k in KS:
    for kind in ('crowded', 'sym', 'separate'):
        frags = [[j] for j in range(NQ)] + [list(range(NQ))]
        chi, depth = group_scan(k, NQ, kind, LAM, frags)
        GD[(k, kind)] = (chi, depth)
        for lab, sl in (('site', slice(0, NQ)), ('whole', slice(NQ, NQ + 1))):
            c = chi[sl]
            d1 = c[:, depth == 1].max(axis=1)
            st = c[:, 1:].max(axis=1)
            am = 1 + c[:, 1:].argmax(axis=1)
            dstar = depth[am]
            prof = [c[:, depth == d].mean() if (depth == d).any() else float('nan') for d in range(1, 7)]
            P(f"{k:>3} {kind:>9} {lab:>6} | {d1.mean():>10.5f} {st.mean():>9.5f} "
              f"{(st-d1).mean():>9.6f} {dstar.mean():>7.3f} {float((dstar>1).mean()):>13.3f} | "
              + " ".join(f"{v:>8.5f}" for v in prof))
    P("-" * 150)
P("chi_depth1 = the best-known SINGLE record; chi_star = the best-known record in the whole group;")
P("EXCESS = chi_star - chi_depth1; depth* = |S| of the best-known record; d1..d6 = mean chi over")
P("all group elements of that depth.  In the 'separate' control a site is coupled to ONE record,")
P("so depth* must be 1 and EXCESS must be 0.")

# ------------------------------------------------------------------ TABLE E: 3-body
P("")
P("=" * 132)
P("TABLE E  CONNECTED 3-BODY TERM.  C3 = chi_ijl - (chi_ij+chi_il+chi_jl) + (chi_i+chi_j+chi_l).")
P("          Needs k >= 3, so it does not exist at k = 2 at all.  nq = 6, single-site fragments.")
P("=" * 132)
P(f"{'k':>3} | {'#triples':>9} | {'C3 mean crd':>12} {'C3 max|.| crd':>13} | "
  f"{'C3 mean sep':>12} {'C3 max|.| sep':>13} | {'SYN mean crd':>12} {'SYN mean sep':>12}")
P("-" * 116)
for k in KS:
    if k < 3:
        P(f"{k:>3} | {'0':>9} | {'--':>12} {'--':>13} | {'--':>12} {'--':>13} | "
          f"{tabC[(k,'crowded','site')]['syn'].mean():>12.6f} {tabC[(k,'separate','site')]['syn'].mean():>12.2e}")
        continue
    PVk = parity_matrix(k)
    tri = list(itertools.combinations(range(k), 3))
    if len(tri) > 40: tri = tri[:40]
    specs = [[]]
    ix = {}
    def add(sp):
        specs.append(sp); return len(specs) - 1
    for (i, j, l) in tri:
        e = [1 << i, 1 << j, 1 << l]
        ix[(i, j, l)] = dict(
            s=[[add([(e[q], +1)]), add([(e[q], -1)])] for q in range(3)],
            p=[[add([(e[a], sa), (e[b], sb)]) for sa in (+1, -1) for sb in (+1, -1)]
               for a, b in ((0, 1), (0, 2), (1, 2))],
            t=[add([(e[0], a), (e[1], b), (e[2], c)]) for a in (+1, -1) for b in (+1, -1) for c in (+1, -1)])
    rows = mask_rows(PVk, specs)
    res = {}
    for kind in ('crowded', 'separate'):
        W = weights(kind, k, NQ); B = Broadcast(k, NQ, W, LAM)
        vals = []
        for j in range(NQ):
            a = np.zeros(len(specs))
            for ti in range(len(TIMES)): a += entropy_stack(B, [j], ti, rows)
            a /= len(TIMES); S0 = a[0]
            for T in tri:
                d = ix[T]
                cs = [max(S0 - 0.5 * (a[p] + a[q]), 0.0) for p, q in d['s']]
                cp = [max(S0 - 0.25 * sum(a[q] for q in grp), 0.0) for grp in d['p']]
                ct = max(S0 - 0.125 * sum(a[q] for q in d['t']), 0.0)
                vals.append(ct - sum(cp) + sum(cs))
        res[kind] = np.array(vals)
    P(f"{k:>3} | {len(tri):>9} | {res['crowded'].mean():>12.6f} {np.abs(res['crowded']).max():>13.6f} | "
      f"{res['separate'].mean():>12.2e} {np.abs(res['separate']).max():>13.2e} | "
      f"{tabC[(k,'crowded','site')]['syn'].mean():>12.6f} {tabC[(k,'separate','site')]['syn'].mean():>12.2e}")
P("-" * 116)

# ------------------------------------------------------------------ TABLE F: robustness (D-17)
P("")
P("=" * 132)
P("TABLE F  D-17 -- VARY THE VENUE'S OWN SCALES.  EXCESS and depth* over three weight seeds,")
P("          three couplings and two bath sizes, single-site fragments, crowded geometry.")
P("          The 'separate' control at the same settings is in the last two columns.")
P("=" * 132)
P(f"{'k':>3} {'nq':>3} {'lam':>5} {'seed':>5} | {'chi_depth1':>10} {'chi_star':>9} {'EXCESS':>9} "
  f"{'depth*':>7} | {'EXCESS sep':>10} {'depth* sep':>10}")
P("-" * 96)
for nq in (4, 6):
    for lam in (0.4, 0.8, 1.2):
        for seed in (7, 11, 23):
            for k in (2, 4, 6, 8):
                frags = [[j] for j in range(nq)]
                chi, depth = group_scan(k, nq, 'crowded', lam, frags, seed=seed)
                d1 = chi[:, depth == 1].max(axis=1); st = chi[:, 1:].max(axis=1)
                dstar = depth[1 + chi[:, 1:].argmax(axis=1)]
                chs, _ = group_scan(k, nq, 'separate', lam, frags, seed=seed)
                d1s = chs[:, depth == 1].max(axis=1); sts = chs[:, 1:].max(axis=1)
                dss = depth[1 + chs[:, 1:].argmax(axis=1)]
                P(f"{k:>3} {nq:>3} {lam:>5.1f} {seed:>5} | {d1.mean():>10.5f} {st.mean():>9.5f} "
                  f"{(st-d1).mean():>9.6f} {dstar.mean():>7.3f} | {(sts-d1s).mean():>10.2e} {dss.mean():>10.3f}")
        P("-" * 96)

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY/s4_relations.txt", "w").write(OUT.getvalue())
