"""S6 -- THE MASTER TABLE.  One k-series, every quantity, control column in the same table.

FINE k SERIES.  [[n,n-2,2]] gives k = n-2, always EVEN.  Odd k is reached by coupling only k of
the n-2 records on the carrier [[k+2 or k+3, ., 2]] -- the uncoupled records contribute nothing to
the bath Hamiltonian, so the broadcast state is EXACTLY the k-coupled-record state and every chi
about an uncoupled record, or about any group element involving one, is exactly zero.  That is
itself a reading: an uncoupled record is invisible however many coupled ones surround it.

THRESHOLD-FREE REDUNDANCY.  At bath sizes 4-6 the whole-bath chi is nowhere near the classical
plateau, so any (1-delta) threshold measure sits pinned near R = 1 and cannot resolve k.  The
threshold-free companion carried here is

    E(k) = nq * mean_j chi(R_i : site j) / chi(R_i : whole bath)

-- the number of single-site copies the whole-bath information is worth.  E = nq means every site
holds the whole thing (perfect objectivity); E -> 0 means the record is only legible to the bath
as a whole.  E is reported beside R_delta, never instead of it.
"""
import sys, io, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY")
from common import *
from s4_helpers import group_scan_chi, group_scan_chi_fast

OUT = io.StringIO()
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.write(s + "\n")

LAM = 0.8
KS = list(range(1, 11))


def carrier_for(k):
    n = k + 2 if (k + 2) % 2 == 0 else k + 3
    return n


def full_profile(k, nq, kind, lam, seed=7, times=TIMES, do_group=True):
    W = weights(kind, k, nq, seed=seed)
    B = Broadcast(k, nq, W, lam, times=times)
    allsub = [tuple(c) for f in range(1, nq + 1) for c in itertools.combinations(range(nq), f)]
    specs = [(f"R{i}", [(0.5, {i: +1}), (0.5, {i: -1})]) for i in range(k)]
    Stau = thermal_entropies(nq)
    chi = {s: np.zeros(k + 1) for s in allsub}
    for ti in range(len(times)):
        for s in allsub:
            r = chi_batch(B, list(s), ti, specs)
            chi[s] += np.array([r[f"R{i}"] for i in range(k)] +
                               [max(r['_Savg'] - float(Stau[list(s)].sum()), 0.0)])
    for s in allsub: chi[s] /= len(times)
    whole = tuple(range(nq))
    cw = chi[whole][:k]; csite = np.array([chi[(j,)][:k] for j in range(nq)])
    out = dict(chi_whole=cw, chi_site=csite, chi_all=chi[whole][k],
               chi_all_site=np.array([chi[(j,)][k] for j in range(nq)]))
    out['E'] = float(np.mean([nq * csite[:, i].mean() / cw[i] for i in range(k) if cw[i] > 1e-12]))
    for d in (0.10, 0.25):
        thr = (1 - d) * cw
        out[f'N1_{d}'] = float((csite >= thr[None, :]).sum(axis=0).mean())
        fd = []
        for i in range(k):
            if cw[i] < 1e-12: continue
            got = nq
            for f in range(1, nq + 1):
                if float(np.mean([chi[s][i] for s in allsub if len(s) == f])) >= thr[i]:
                    got = f; break
            fd.append(got)
        out[f'R_{d}'] = float(np.mean([nq / q for q in fd])) if fd else float('nan')
    if do_group:
        frags = [[j] for j in range(nq)] + [list(range(nq))]
        g = group_scan_chi_fast(k, nq, W, lam, frags, times=times)
        depth = np.array([bin(m).count('1') for m in range(2 ** k)])
        d1 = g[:, depth == 1].max(axis=1); st = g[:, 1:].max(axis=1)
        ds = depth[1 + g[:, 1:].argmax(axis=1)]
        nb = np.array([int((g[q, 1:] > d1[q] + 1e-12).sum()) for q in range(g.shape[0])])
        out.update(chi_d1_site=float(d1[:nq].mean()), chi_star_site=float(st[:nq].mean()),
                   exc_site=float((st - d1)[:nq].mean()), dstar_site=float(ds[:nq].mean()),
                   nbetter_site=float(nb[:nq].mean()),
                   chi_d1_whole=float(d1[nq]), chi_star_whole=float(st[nq]),
                   exc_whole=float(st[nq] - d1[nq]), dstar_whole=float(ds[nq]),
                   nbetter_whole=int(nb[nq]))
    return out


# SELF-CHECK: the Walsh route must reproduce the mask route exactly, or every group-scan
# number below is worthless.
_d = 0.0
for _k, _nq in ((2, 4), (3, 5), (4, 4)):
    _W = weights('crowded', _k, _nq); _f = [[0], [1], list(range(_nq))]
    _d = max(_d, float(np.abs(group_scan_chi(_k, _nq, _W, LAM, _f, times=TIMES[:6])
                              - group_scan_chi_fast(_k, _nq, _W, LAM, _f, times=TIMES[:6])).max()))
P(f"SELF-CHECK  Walsh group scan vs explicit-mask group scan: max |diff| = {_d:.3e} -> "
  + ("PASSED" if _d < 1e-10 else "FAILED, DRAW NO CONCLUSION"))
P("")
P("=" * 156)
P("S6  MASTER TABLE -- one k-series from k = 1 to k = 10, nq = 6, lam = 0.8")
P("=" * 156)
P("")
P("carrier used per k:  " + ", ".join(f"k={k}->[[{carrier_for(k)},{carrier_for(k)-2},2]]" for k in KS))
P("")
P(f"{'k':>3} {'geom':>9} | {'chi_wh':>8} {'chi_site':>8} {'E=nq*rat':>9} {'N1(.1)':>7} {'R(.1)':>6} "
  f"{'R(.25)':>7} | {'chi_ALL':>8} {'sum_i':>8} | {'exc site':>9} {'d* site':>8} {'Nbet site':>9} | "
  f"{'chi_d1 wh':>9} {'chi_* wh':>9} {'exc wh':>8} {'d* wh':>6} {'ratio wh':>9}")
P("-" * 156)
MASTER = {}
for kind in ('crowded', 'sym', 'separate'):
    for k in KS:
        r = full_profile(k, 6, kind, LAM)
        MASTER[(k, kind)] = r
        rat = r['chi_star_whole'] / r['chi_d1_whole'] if r['chi_d1_whole'] > 1e-12 else float('nan')
        P(f"{k:>3} {kind:>9} | {r['chi_whole'].mean():>8.5f} {r['chi_site'].mean():>8.5f} "
          f"{r['E']:>9.4f} {r['N1_0.1']:>7.3f} {r['R_0.1']:>6.3f} {r['R_0.25']:>7.3f} | "
          f"{r['chi_all']:>8.5f} {float(r['chi_whole'].sum()):>8.5f} | "
          f"{r['exc_site']:>9.6f} {r['dstar_site']:>8.3f} {r['nbetter_site']:>9.2f} | "
          f"{r['chi_d1_whole']:>9.5f} {r['chi_star_whole']:>9.5f} {r['exc_whole']:>8.5f} "
          f"{r['dstar_whole']:>6.1f} {rat:>9.4f}")
    P("-" * 156)

# ------------------------------------------------------------------ constant sites per record
P("")
P("=" * 156)
P("TABLE J  BATH GROWN IN PROPORTION TO k AT CONSTANT SITES-PER-RECORD c = nq/k.")
P("          This is the D-16 separation: c fixed removes the crowding of the bath, leaving only")
P("          the record count.  Control column 'separate' in the same table.")
P("=" * 156)
P(f"{'c':>4} {'k':>3} {'nq':>3} | {'chi_wh crd':>10} {'chi_st crd':>10} {'E crd':>8} {'R(.1) crd':>9} "
  f"{'R(.25) crd':>10} | {'chi_wh sep':>10} {'chi_st sep':>10} {'E sep':>8} {'R(.1) sep':>9} {'R(.25) sep':>10}")
P("-" * 132)
for c, pairs in ((1, [(2, 2), (4, 4), (6, 6), (8, 8)]),
                 (2, [(2, 4), (3, 6), (4, 8)]),
                 (3, [(1, 3), (2, 6)])):
    for k, nq in pairs:
        a = full_profile(k, nq, 'crowded', LAM, do_group=False)
        b = full_profile(k, nq, 'separate', LAM, do_group=False)
        P(f"{c:>4} {k:>3} {nq:>3} | {a['chi_whole'].mean():>10.5f} {a['chi_site'].mean():>10.5f} "
          f"{a['E']:>8.4f} {a['R_0.1']:>9.3f} {a['R_0.25']:>10.3f} | "
          f"{b['chi_whole'].mean():>10.5f} {b['chi_site'].mean():>10.5f} {b['E']:>8.4f} "
          f"{b['R_0.1']:>9.3f} {b['R_0.25']:>10.3f}")
    P("-" * 132)

# ------------------------------------------------------------------ the pointed verdict
P("")
P("=" * 156)
P("TABLE K  THE POINTED QUESTION -- is any quantity ZERO or TRIVIAL at small k and NON-ZERO at")
P("          larger k?  Each candidate is printed across the whole k-series with its control.")
P("=" * 156)
CANDS = [
    ("EXCESS, single site, crowded", 'exc_site', 'crowded'),
    ("EXCESS, single site, sym", 'exc_site', 'sym'),
    ("EXCESS, single site, SEPARATE (control)", 'exc_site', 'separate'),
    ("EXCESS, whole bath, crowded", 'exc_whole', 'crowded'),
    ("EXCESS, whole bath, sym", 'exc_whole', 'sym'),
    ("EXCESS, whole bath, SEPARATE (control)", 'exc_whole', 'separate'),
    ("depth* of best-known record, whole, sym", 'dstar_whole', 'sym'),
    ("depth* of best-known record, whole, crowded", 'dstar_whole', 'crowded'),
    ("depth* of best-known record, whole, SEPARATE", 'dstar_whole', 'separate'),
    ("N_better, single site, sym", 'nbetter_site', 'sym'),
    ("N_better, single site, crowded", 'nbetter_site', 'crowded'),
    ("N_better, single site, SEPARATE (control)", 'nbetter_site', 'separate'),
    ("chi(best single record), whole, sym", 'chi_d1_whole', 'sym'),
    ("chi(best record in group), whole, sym", 'chi_star_whole', 'sym'),
    ("E = nq*chi_site/chi_whole, crowded", 'E', 'crowded'),
    ("E = nq*chi_site/chi_whole, SEPARATE (control)", 'E', 'separate'),
]
P(f"{'quantity':>46} | " + " ".join(f"{'k='+str(k):>9}" for k in KS))
P("-" * 156)
for lab, key, kind in CANDS:
    P(f"{lab:>46} | " + " ".join(f"{MASTER[(k,kind)][key]:>9.5f}" for k in KS))
P("-" * 156)
P("")
P("VERDICT LINES -- filled from the numbers above, not in advance:")
for lab, key, kind in CANDS:
    v = np.array([MASTER[(k, kind)][key] for k in KS])
    small = float(np.max(np.abs(v[:2])))         # k = 1, 2
    large = float(np.max(np.abs(v[6:])))         # k = 7..10
    if small < 1e-8 and large > 1e-6:
        verd = "ZERO at small k, NON-ZERO at large k"
    elif small < 1e-8 and large < 1e-8:
        verd = "zero everywhere"
    elif large > 2 * small + 1e-9:
        verd = "non-zero at small k, GROWS with k"
    elif large < 0.5 * small:
        verd = "non-zero at small k, SHRINKS with k"
    else:
        verd = "non-zero at small k, roughly flat"
    P(f"   {lab:>46} : {verd}")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY/s6_master.txt", "w").write(OUT.getvalue())
