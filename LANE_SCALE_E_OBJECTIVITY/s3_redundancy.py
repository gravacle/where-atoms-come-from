"""S3 -- REDUNDANCY OF EACH RECORD, AND HOW IT SCALES WITH k.

chi held by the WHOLE bath and by EVERY fragment, time-averaged over 25 times in [1,13].
R_delta is reported two ways, both from the time-averaged chi:

   N1(delta)  = how many DISJOINT single-qubit fragments each hold at least (1-delta) of the
                whole-bath chi.  This is the literal reading of the brief.
   f(delta)   = the smallest fragment SIZE whose AVERAGE chi reaches (1-delta) of the whole-bath
                chi;  R_delta = nq / f(delta).  This is the standard quantum-Darwinism measure.

CONTROL IN THE SAME TABLE (D-15, D-16): 'separate' -- site j couples to record (j mod k) ALONE,
at the SAME total coupling per site.  That state is exactly k independent carriers on k disjoint
baths (verified through the full model in s2), so it is the SPREAD control against which crowding
is measured -- never against an 'alone' value.

D-17: the venue's own scales are varied -- lam in {0.4, 0.8, 1.2}, bath nq in {4, 6}, and the
degenerate 'sym' geometry is carried alongside the generic one.
"""
import sys, io, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY")
from common import *

OUT = io.StringIO()
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.write(s + "\n")

DELTAS = (0.10, 0.25)


def redundancy_profile(k, nq, kind, lam, times=TIMES, seed=7):
    """time-averaged chi_i(F) for EVERY fragment F, plus derived redundancy measures."""
    W = weights(kind, k, nq, seed=seed)
    B = Broadcast(k, nq, W, lam, times=times)
    nT = len(times)
    specs = [(f"R{i}", [(0.5, {i: +1}), (0.5, {i: -1})]) for i in range(k)]
    Stau = thermal_entropies(nq)
    allsub = [tuple(c) for f in range(1, nq + 1) for c in itertools.combinations(range(nq), f)]
    chi = {s: np.zeros(k + 1) for s in allsub}
    for ti in range(nT):
        for s in allsub:
            r = chi_batch(B, list(s), ti, specs)
            call = max(r['_Savg'] - float(Stau[list(s)].sum()), 0.0)
            chi[s] += np.array([r[f"R{i}"] for i in range(k)] + [call])
    for s in allsub: chi[s] /= nT
    whole = tuple(range(nq))
    cw = chi[whole][:k]                                    # per-record whole-bath chi
    csite = np.array([chi[(j,)][:k] for j in range(nq)])    # (nq, k)
    out = dict(k=k, nq=nq, kind=kind, lam=lam, W=W,
               chi_whole=cw, chi_site=csite, chi_all_whole=chi[whole][k],
               chi_all_site=np.array([chi[(j,)][k] for j in range(nq)]))
    for d in DELTAS:
        thr = (1 - d) * cw
        out[f"N1_{d}"] = (csite >= thr[None, :]).sum(axis=0)
        fdel = np.zeros(k)
        for i in range(k):
            if cw[i] < 1e-9: fdel[i] = np.nan; continue
            got = nq
            for f in range(1, nq + 1):
                subs = [s for s in allsub if len(s) == f]
                av = float(np.mean([chi[s][i] for s in subs]))
                if av >= thr[i]: got = f; break
            fdel[i] = got
        out[f"f_{d}"] = fdel
        out[f"R_{d}"] = nq / fdel
    return out


def mean_nan(a):
    a = np.asarray(a, dtype=float)
    return float(np.nanmean(a)) if a.size else float('nan')


P("=" * 130)
P("S3  REDUNDANCY vs NUMBER OF RECORDS -- crowded and its SPREAD control in the same table")
P("=" * 130)
P("")
P(f"times: 25 points in [1,13];  beta = {BETA};  bath energies {ENERGIES[:6]}")
P("'crowded'  = every record couples to every bath site, generic weights, sum_i W[i,j]^2 = 1")
P("'separate' = site j couples to record (j mod k) alone, same total coupling per site")
P("             (= k independent carriers on disjoint baths; for k > nq the surplus records are")
P("             unhosted and hold chi = 0, which is itself the finding for that geometry)")
P("")

KS = [2, 4, 6, 8, 10]
LAM = 0.8

# ---------------- SELF-CHECK: the analytic chi_ALL identity ----------------
# chi_ALL(F) = S(rho-bar_F) - sum_{j in F} S(tau_j) is used everywhere below because the 2^k-fold
# brute force is unaffordable at k = 10.  It rests on every conditional bath state being a
# unitary image of the thermal state.  If this check fails, NOTHING below is interpretable.
_bad = 0.0
for _k, _nq, _kind, _lam in ((2, 4, 'crowded', 0.8), (4, 4, 'crowded', 0.8),
                             (3, 5, 'crowded', 1.2), (4, 4, 'separate', 0.8)):
    _W = weights(_kind, _k, _nq); _B = Broadcast(_k, _nq, _W, _lam, times=(2.0, 7.0))
    _St = thermal_entropies(_nq)
    for _ti in (0, 1):
        for _F in ([0], [1, 2], list(range(_nq))):
            _brute = _B.chi_all(_F, _ti)
            _r = chi_batch(_B, _F, _ti, [("R0", [(0.5, {0: +1}), (0.5, {0: -1})])])
            _ana = max(_r['_Savg'] - float(_St[_F].sum()), 0.0)
            _bad = max(_bad, abs(_brute - _ana))
P(f"SELF-CHECK  analytic chi_ALL vs brute-force 2^k averaging: max |diff| = {_bad:.3e}  -> "
  + ("PASSED" if _bad < 1e-10 else "FAILED, DRAW NO CONCLUSION"))
P("")

# ------------------------------------------------------------------ TABLE A: FIXED bath
for nq in (6, 4):
    P("=" * 130)
    P(f"TABLE A(nq={nq})  FIXED BATH SIZE nq = {nq}, lam = {LAM}.  Columns for 'crowded' and the "
      f"'separate' SPREAD CONTROL side by side.")
    P("=" * 130)
    P(f"{'k':>3} {'nq':>3} | {'chi_wh crd':>10} {'chi_wh sep':>10} {'crd/sep':>8} | "
      f"{'chi_st crd':>10} {'chi_st sep':>10} | {'N1(.1) crd':>10} {'N1(.1) sep':>10} | "
      f"{'f(.1) crd':>9} {'f(.1) sep':>9} | {'R(.1) crd':>9} {'R(.1) sep':>9} | "
      f"{'chiALL crd':>10} {'sum_i crd':>9} {'chiALL sep':>10} {'sum_i sep':>9}")
    P("-" * 172)
    store = {}
    for k in KS:
        c = redundancy_profile(k, nq, 'crowded', LAM)
        s = redundancy_profile(k, nq, 'separate', LAM)
        store[k] = (c, s)
        ratio = mean_nan(c['chi_whole']) / mean_nan(s['chi_whole']) if mean_nan(s['chi_whole']) > 1e-12 else float('nan')
        P(f"{k:>3} {nq:>3} | {mean_nan(c['chi_whole']):>10.5f} {mean_nan(s['chi_whole']):>10.5f} "
          f"{ratio:>8.3f} | {mean_nan(c['chi_site']):>10.5f} {mean_nan(s['chi_site']):>10.5f} | "
          f"{mean_nan(c['N1_0.1']):>10.3f} {mean_nan(s['N1_0.1']):>10.3f} | "
          f"{mean_nan(c['f_0.1']):>9.3f} {mean_nan(s['f_0.1']):>9.3f} | "
          f"{mean_nan(c['R_0.1']):>9.3f} {mean_nan(s['R_0.1']):>9.3f} | "
          f"{c['chi_all_whole']:>10.5f} {float(np.sum(c['chi_whole'])):>9.5f} "
          f"{s['chi_all_whole']:>10.5f} {float(np.sum(s['chi_whole'])):>9.5f}")
    P("-" * 172)
    P("chi_wh = whole-bath chi per record (mean over the k records);  chi_st = single-site chi")
    P("(mean over records and sites);  N1 = # single sites holding >= 0.9 of the whole-bath chi;")
    P("f = smallest fragment size whose AVERAGE chi reaches 0.9 of whole-bath; R = nq/f.")
    P("chiALL = Holevo of the whole k-record register; sum_i = sum of the individual chi.")
    P("")
    if nq == 6:
        P("PER-RECORD DETAIL at nq=6, crowded (shows the spread across records, not just the mean):")
        for k in KS:
            c = store[k][0]
            P(f"   k={k:>2}  chi_whole per record: " +
              " ".join(f"{v:6.4f}" for v in c['chi_whole']))
            P(f"        N1(0.10) per record: " + " ".join(f"{int(v):6d}" for v in c['N1_0.1']) +
              f"   N1(0.25): " + " ".join(f"{int(v):6d}" for v in c['N1_0.25']))
        P("")
        P("PER-RECORD DETAIL at nq=6, separate (the SPREAD control):")
        for k in KS:
            s = store[k][1]
            P(f"   k={k:>2}  chi_whole per record: " +
              " ".join(f"{v:6.4f}" for v in s['chi_whole']) +
              f"   hosted records: {int((s['chi_whole']>1e-9).sum())}/{k}")
        P("")

# ------------------------------------------------------------------ TABLE B: PROPORTIONAL bath
P("=" * 130)
P("TABLE B  BATH GROWN IN PROPORTION TO k:  nq = k  (one bath site per record), lam = 0.8.")
P("          The 'separate' control here gives every record exactly one private site.")
P("=" * 130)
P(f"{'k':>3} {'nq':>3} | {'chi_wh crd':>10} {'chi_wh sep':>10} {'crd/sep':>8} | "
  f"{'chi_st crd':>10} {'chi_st sep':>10} | {'N1(.1) crd':>10} {'N1(.1) sep':>10} | "
  f"{'f(.1) crd':>9} {'f(.1) sep':>9} | {'R(.1) crd':>9} {'R(.1) sep':>9} | "
  f"{'chiALL crd':>10} {'sum_i crd':>9}")
P("-" * 152)
for k in [2, 4, 6, 8]:
    nq = k
    c = redundancy_profile(k, nq, 'crowded', LAM)
    s = redundancy_profile(k, nq, 'separate', LAM)
    ratio = mean_nan(c['chi_whole']) / mean_nan(s['chi_whole'])
    P(f"{k:>3} {nq:>3} | {mean_nan(c['chi_whole']):>10.5f} {mean_nan(s['chi_whole']):>10.5f} "
      f"{ratio:>8.3f} | {mean_nan(c['chi_site']):>10.5f} {mean_nan(s['chi_site']):>10.5f} | "
      f"{mean_nan(c['N1_0.1']):>10.3f} {mean_nan(s['N1_0.1']):>10.3f} | "
      f"{mean_nan(c['f_0.1']):>9.3f} {mean_nan(s['f_0.1']):>9.3f} | "
      f"{mean_nan(c['R_0.1']):>9.3f} {mean_nan(s['R_0.1']):>9.3f} | "
      f"{c['chi_all_whole']:>10.5f} {float(np.sum(c['chi_whole'])):>9.5f}")
P("-" * 152)

# ------------------------------------------------------------------ D-17: vary lam, and the sym venue
P("")
P("=" * 130)
P("D-17  VARY THE VENUE'S OWN SCALES.  Same quantities at three couplings and in the DEGENERATE")
P("      'sym' geometry (all weights equal), nq = 6.")
P("=" * 130)
P(f"{'geom':>9} {'lam':>5} {'k':>3} | {'chi_wh':>8} {'chi_st':>8} {'N1(.1)':>7} {'f(.1)':>6} "
  f"{'R(.1)':>6} {'chiALL':>8} {'sum_i':>8} {'chiALL/sum':>10}")
P("-" * 96)
for kind in ('crowded', 'sym', 'separate'):
    for lam in (0.4, 0.8, 1.2):
        for k in (2, 4, 6, 8, 10):
            r = redundancy_profile(k, 6, kind, lam)
            si = float(np.sum(r['chi_whole']))
            P(f"{kind:>9} {lam:>5.1f} {k:>3} | {mean_nan(r['chi_whole']):>8.5f} "
              f"{mean_nan(r['chi_site']):>8.5f} {mean_nan(r['N1_0.1']):>7.3f} "
              f"{mean_nan(r['f_0.1']):>6.3f} {mean_nan(r['R_0.1']):>6.3f} "
              f"{r['chi_all_whole']:>8.5f} {si:>8.5f} "
              f"{(r['chi_all_whole']/si if si>1e-12 else float('nan')):>10.4f}")
    P("-" * 96)

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY/s3_redundancy.txt", "w").write(OUT.getvalue())
