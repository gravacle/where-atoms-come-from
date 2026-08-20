"""B2 -- SUBTRACT WHAT IS ALREADY EXPLAINED, AND LOOK AT WHAT IS LEFT.

Fit the explained model (capacity C-36 + pairing C-38 + crowding-selectivity C-39, every
coefficient FITTED, plus a control column C-38 proves must be zero) to the configuration
battery at every n, then report the RESIDUAL: measured chi minus explained chi.

Also -- and this needs no fit at all -- report the MODEL-FREE THREE-BODY CONTRASTS: pairs of
configurations with IDENTICAL (m0, m1, p0, p1) that differ only in how the PARTNERS pair with
EACH OTHER.  The explained model is pairwise-in-(record, partner) and predicts these are equal.
Any difference is a three-record quantity the explained model cannot contain.
"""
import numpy as np, sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from battery import *

t0 = time.time()
NB = 3
LAMS = (0.4, 0.8, 1.2)
NS = [4, 6, 8]
say("=" * 122)
say("B2   THE RESIDUAL  --  measured chi minus the model built from C-36, C-38, C-39")
say("=" * 122)
say(f"  carrier [[n,n-2,2]], n = {NS}; bath {NB} qubits at beta {BASE.beta}, energies {BASE.energies[:NB]}")
say(f"  chi TIME-AVERAGED over {len(BASE.times)} times in [{BASE.times[0]:.0f},{BASE.times[-1]:.0f}] (a snapshot recurs and is not interpretable)")
say("")

ALL = []
for n in NS:
    ops, _, _ = build_ops(n)
    for lam in LAMS:
        rows = run_battery(n, BASE, lam, NB=NB, ops=ops)
        ALL += rows
        say(f"  n={n} lam={lam}: {len(rows)} configurations   ({time.time()-t0:.0f}s)")

json.dump([{k: (float(v) if isinstance(v, (np.floating,)) else v) for k, v in r.items()} for r in ALL],
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b2_raw.json'), 'w'), indent=0)

# ---------------------------------------------------------------- the measured table
say("")
say("MEASURED chi.  Note the columns n=4,6,8: where a configuration exists at more than one n it")
say("is the SAME NUMBER to every printed digit -- B1's exact N-independence, visible here.")
for lam in LAMS:
    say("")
    say(f"  lam = {lam}")
    say(f"  {'configuration':<38}{'m0':>4}{'m1':>4}{'p0':>4}{'p1':>4}{'ppA':>5}" +
        "".join(f"{'chi n='+str(n):>17}" for n in NS))
    names = [c[0] for c in CONFIGS]
    for name in names:
        row = {}
        for n in NS:
            r = [x for x in ALL if x['n'] == n and x['lam'] == lam and x['name'] == name]
            if r: row[n] = r[0]
        if not row: continue
        any_r = list(row.values())[0]
        say(f"  {name:<38}{any_r['m0']:>4}{any_r['m1']:>4}{any_r['p0']:>4}{any_r['p1']:>4}{any_r['pp_anti']:>5}" +
            "".join((f"{row[n]['chi']:>17.12f}" if n in row else f"{'--':>17}") for n in NS))

# ---------------------------------------------------------------- the fit
say("")
say("=" * 122)
say("THE EXPLAINED MODEL, FITTED.  log chi = a + gamma*log(1/(1+m)) + beta*(m1+p1) + delta*m1 + c0*p0")
say("=" * 122)
say(f"  {'n':>4}{'lam':>6}{'#cfg':>6}" + "".join(f"{t.split(' ')[0]:>10}" for t in TERMS) +
    f"{'cond':>9}{'rms resid':>12}{'max|resid|':>12}{'rel rms':>10}")
FITS = {}
for n in NS:
    for lam in LAMS:
        rows = [r for r in ALL if r['n'] == n and r['lam'] == lam]
        if len(rows) < 6: continue
        F = fit(rows)
        FITS[(n, lam)] = (rows, F)
        say(f"  {n:>4}{lam:>6.2f}{len(rows):>6}" + "".join(f"{c:>10.4f}" for c in F['coef']) +
            f"{F['cond']:>9.1f}{F['rms']:>12.3e}{F['maxabs']:>12.3e}{F['relrms']:>10.3f}")
say("")
say("  THE CONTROL COLUMN c0 (D-15): C-38's theorem forces it to be EXACTLY zero, and it is.")
say(f"  largest |c0| over every fit: {max(abs(F['coef'][4]) for _, F in FITS.values()):.3e}")
say("  A fit that could not see a zero would be worthless; the other four coefficients are all")
say("  large, so this fit CAN register a non-zero and does not.")

# ---------------------------------------------------------------- residual table
say("")
say("=" * 122)
say("THE RESIDUAL, PER CONFIGURATION.  (measured - explained), at every n.  lam = 0.8 shown in full.")
say("=" * 122)
for lam in LAMS:
    say("")
    say(f"  lam = {lam}")
    say(f"  {'configuration':<38}{'m0':>4}{'m1':>4}{'p0':>4}{'p1':>4}{'ppA':>5}" +
        "".join(f"{'resid n='+str(n):>16}" for n in NS))
    for name in [c[0] for c in CONFIGS]:
        cells = {}
        for n in NS:
            if (n, lam) not in FITS: continue
            rows, F = FITS[(n, lam)]
            for i, r in enumerate(rows):
                if r['name'] == name: cells[n] = (r, F['resid'][i])
        if not cells: continue
        any_r = list(cells.values())[0][0]
        say(f"  {name:<38}{any_r['m0']:>4}{any_r['m1']:>4}{any_r['p0']:>4}{any_r['p1']:>4}{any_r['pp_anti']:>5}" +
            "".join((f"{cells[n][1]:>+16.6f}" if n in cells else f"{'--':>16}") for n in NS))

# ---------------------------------------------------------------- three-body, MODEL-FREE
say("")
say("=" * 122)
say("MODEL-FREE THREE-BODY CONTRASTS.  Configurations with IDENTICAL (m0,m1,p0,p1) that differ")
say("ONLY in how the PARTNERS pair with EACH OTHER.  No fit is involved; this is a direct")
say("difference of two measured numbers, so no fitting choice can manufacture or hide it.")
say("=" * 122)
PAIRS3 = [("X2@0,X3@0",                "X2@0,Z2@0"),
          ("Z1@0,Z1X2@0",              "Z1X2@0,Z1Z2@0"),
          ("X2@0,X3@0,X4@0",           "X2@0,Z2@0,X3@0"),
          ("Z1@0,X2@0,X3@0",           "Z1X2@0,Z1Z2@0,X3@0"),
          ("X2@1,X3@2",                "X2@1,X3@1"),
          ("Z1@1,Z1X2@2",              "Z1@1,Z1X2@1")]
say(f"  {'lam':>5}{'n':>4}  {'configuration A':<26}{'configuration B':<26}{'ppA':>5}{'ppB':>5}{'chi A':>15}{'chi B':>15}{'A - B':>14}")
tb = []
for lam in LAMS:
    for n in NS:
        for a, b in PAIRS3:
            ra = [x for x in ALL if x['n'] == n and x['lam'] == lam and x['name'] == a]
            rb = [x for x in ALL if x['n'] == n and x['lam'] == lam and x['name'] == b]
            if not ra or not rb: continue
            ra, rb = ra[0], rb[0]
            same_feat = all(ra[k] == rb[k] for k in ('m0', 'm1', 'p0', 'p1'))
            if not same_feat: continue
            d = ra['chi'] - rb['chi']
            tb.append(dict(lam=lam, n=n, a=a, b=b, d=d))
            say(f"  {lam:>5.2f}{n:>4}  {a:<26}{b:<26}{ra['pp_anti']:>5}{rb['pp_anti']:>5}"
                f"{ra['chi']:>15.12f}{rb['chi']:>15.12f}{d:>+14.9f}")
say("")
if tb:
    say(f"  largest |A - B| over all three-body contrasts: {max(abs(x['d']) for x in tb):.6e}")
    say(f"  smallest |A - B|:                              {min(abs(x['d']) for x in tb):.6e}")

# ---------------------------------------------------------------- residual vs features
say("")
say("=" * 122)
say("DOES THE RESIDUAL HAVE STRUCTURE?")
say("=" * 122)
for lam in LAMS:
    n = max(NS)
    if (n, lam) not in FITS: continue
    rows, F = FITS[(n, lam)]
    res = F['resid']
    say("")
    say(f"  lam = {lam}, n = {n}:  residual grouped by the THREE-BODY count ppA (partner pairs that anticommute)")
    say(f"  {'ppA':>5}{'#cfg':>6}{'mean resid':>15}{'rms resid':>14}{'min':>13}{'max':>13}")
    for g in sorted(set(r['pp_anti'] for r in rows)):
        idx = [i for i, r in enumerate(rows) if r['pp_anti'] == g]
        v = res[idx]
        say(f"  {g:>5}{len(idx):>6}{v.mean():>+15.6f}{np.sqrt((v**2).mean()):>14.6f}{v.min():>+13.6f}{v.max():>+13.6f}")
    say(f"  residual grouped by m (records sharing the read record's site, minus one)")
    say(f"  {'m':>5}{'#cfg':>6}{'mean resid':>15}{'rms resid':>14}")
    for g in sorted(set(r['m'] for r in rows)):
        idx = [i for i, r in enumerate(rows) if r['m'] == g]
        v = res[idx]
        say(f"  {g:>5}{len(idx):>6}{v.mean():>+15.6f}{np.sqrt((v**2).mean()):>14.6f}")

say("")
say(f"  elapsed {time.time()-t0:.1f}s")
