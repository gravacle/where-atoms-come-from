"""B6 -- THE MODEL HIERARCHY AND THE UPPER BOUND ON WHAT COULD STILL BE HIDING.

Nested fits on the pooled battery (CONFIGS + EXTRA_CONFIGS) at n = 8:

  M0  capacity only            log chi = a + gamma*log(1/(1+m))
  M1  THE EXPLAINED MODEL      M0 + beta*(m1+p1) + delta*m1 + c0*p0        [C-36 + C-38 + C-39]
  M2  + THREE-RECORD TERMS     M1 + counts over PAIRS OF PARTNERS
  M3  + one four-record term   M2 + (pairs on the read site)^2

Reported against FLOOR-M (float64) and FLOOR-V (the venue floor of B3).  Then: for a candidate
extra term, the smallest coefficient the data could still tell from zero.

THE CONTROL COLUMN c0 IS THE CALIBRATION OF THE FITTING MACHINERY.  C-38's theorem makes it
EXACTLY zero.  Whatever the fit returns for it is the size of a coefficient this design can
invent out of misspecification alone, and no fitted coefficient smaller than that means anything.
"""
import numpy as np, sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from battery import *

t0 = time.time()
NB, N = 3, 8
LAMS = (0.4, 0.8, 1.2)
ALLCFG = CONFIGS + EXTRA_CONFIGS
say("=" * 126)
say("B6   THE MODEL HIERARCHY, AND HOW LARGE A HIDDEN TERM COULD BE")
say("=" * 126)
ops, _, _ = build_ops(N)

DATA = {}
for lam in LAMS:
    rows = run_battery(N, BASE, lam, NB=NB, ops=ops, configs=ALLCFG)
    for r in rows:
        cfg = dict(ALLCFG)[r['name']]
        r.update(three_body_counts(cfg, ops, N))
    DATA[lam] = rows
    say(f"  lam={lam}: {len(rows)} configurations   ({time.time()-t0:.0f}s)")

FEAT = {
 'M0': [('const', lambda r: 1.0), ('gamma cap', lambda r: np.log(1.0 / (1.0 + r['m'])))],
}
FEAT['M1'] = FEAT['M0'] + [('beta pair', lambda r: float(r['m1'] + r['p1'])),
                           ('delta crowd', lambda r: float(r['m1'])),
                           ('c0 CONTROL', lambda r: float(r['p0']))]
FEAT['M2'] = FEAT['M1'] + [('A_read', lambda r: float(r['A_read'])),
                           ('A_colo', lambda r: float(r['A_colo'])),
                           ('A_split', lambda r: float(r['A_split'])),
                           ('C_colo_11', lambda r: float(r['C_colo_11'])),
                           ('C_colo_10', lambda r: float(r['C_colo_10'])),
                           ('C_colo_00', lambda r: float(r['C_colo_00']))]
FEAT['M3'] = FEAT['M2'] + [('A_read^2', lambda r: float(r['A_read'] ** 2)),
                           ('m*A_read', lambda r: float(r['m'] * r['A_read']))]

def dofit(rows, feats):
    A = np.array([[f(r) for _, f in feats] for r in rows])
    y = np.log(np.array([r['chi'] for r in rows]))
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    lres = y - A @ coef
    meas = np.array([r['chi'] for r in rows]); pred = np.exp(A @ coef)
    dof = max(len(rows) - np.linalg.matrix_rank(A), 1)
    s2 = float(lres @ lres) / dof
    return dict(A=A, coef=coef, s=float(np.sqrt(s2)), resid=meas - pred, pred=pred, meas=meas,
                rms=float(np.sqrt(np.mean((meas - pred) ** 2))),
                maxabs=float(np.max(np.abs(meas - pred))),
                relrms=float(np.sqrt(np.mean(((meas - pred) / meas) ** 2))))

try:
    VJ = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b3_venues.json')))
    allv = np.array([VJ[k]['resid'] for k in VJ])
    FLOOR_V = float(np.median(allv.max(axis=0) - allv.min(axis=0)))
except Exception as e:
    FLOOR_V = float('nan')
FLOOR_M = 1.0e-15

say("")
say("THE HIERARCHY.  rms and max are in chi; 'rel rms' is the residual as a fraction of chi.")
say(f"  FLOOR-M (float64, B1/B3) = {FLOOR_M:.1e} in chi.   FLOOR-V (venue floor, B3) = {FLOOR_V:.5f} in chi.")
say("")
say(f"  {'lam':>5}{'model':>6}{'#par':>6}{'rms resid':>12}{'max|resid|':>12}{'rel rms':>10}"
    f"{'rms / FLOOR-V':>15}{'rms / FLOOR-M':>15}")
FITS = {}
for lam in LAMS:
    rows = DATA[lam]
    for m in ('M0', 'M1', 'M2', 'M3'):
        F = dofit(rows, FEAT[m]); FITS[(lam, m)] = F
        say(f"  {lam:>5.2f}{m:>6}{len(FEAT[m]):>6}{F['rms']:>12.4e}{F['maxabs']:>12.4e}{F['relrms']:>10.3f}"
            f"{F['rms']/FLOOR_V:>15.1f}{F['rms']/FLOOR_M:>15.2e}")
    say("")

say("COEFFICIENTS.")
for m in ('M1', 'M2'):
    say("")
    say(f"  {m}")
    say(f"  {'lam':>5}" + "".join(f"{nm:>13}" for nm, _ in FEAT[m]))
    for lam in LAMS:
        say(f"  {lam:>5.2f}" + "".join(f"{c:>13.4f}" for c in FITS[(lam, m)]['coef']))
say("")
say("  THE CONTROL c0 IS EXACTLY ZERO BY THEOREM.  In M1 the fit returns it as large as")
say(f"  {max(abs(FITS[(l,'M1')]['coef'][4]) for l in LAMS):.4f}; in M2 as large as "
    f"{max(abs(FITS[(l,'M2')]['coef'][4]) for l in LAMS):.4f}.")
say("  THAT IS THE CALIBRATION: a coefficient of this size can be produced by misspecification")
say("  alone.  Nothing smaller than it, in any fitted column, means anything.")

say("")
say("HOW MUCH OF THE RESIDUAL DOES EACH THREE-RECORD CHANNEL CARRY?  Drop one column at a time")
say("from M2 and see how far the residual moves back toward M1.")
say(f"  {'lam':>5}  {'column dropped':<14}{'rms resid':>12}{'increase over M2':>19}{'as % of the M1->M2 gain':>26}")
for lam in LAMS:
    F1, F2 = FITS[(lam, 'M1')], FITS[(lam, 'M2')]
    gain = F1['rms'] - F2['rms']
    for i in range(5, len(FEAT['M2'])):
        sub = [f for j, f in enumerate(FEAT['M2']) if j != i]
        Fd = dofit(DATA[lam], sub)
        say(f"  {lam:>5.2f}  {FEAT['M2'][i][0]:<14}{Fd['rms']:>12.4e}{Fd['rms']-F2['rms']:>19.4e}"
            f"{100*(Fd['rms']-F2['rms'])/gain:>26.1f}")
    say("")

say("THE UPPER BOUND ON A FURTHER TERM.  For a candidate feature f, project out everything M2")
say("already fits; the remaining norm |f_perp| sets how small a coefficient the data could still")
say("resolve.  Error scale: the M2 residual, and separately the venue floor of B3.")
say("")
CANDS = {
 "extensive: total records written":     lambda r: float(r['npart'] + 1),
 "records on OTHER sites":               lambda r: float(r['p0'] + r['p1']),
 "site occupancy squared":               lambda r: float(r['m'] ** 2),
 "pairing partners squared":             lambda r: float((r['m1'] + r['p1']) ** 2),
 "four-record: (A_read choose 2)":       lambda r: float(r['A_read'] * (r['A_read'] - 1) / 2),
 "m * (pairing partners)":               lambda r: float(r['m'] * (r['m1'] + r['p1'])),
}
say(f"  {'candidate extra term f':<38}{'lam':>5}{'|f_perp|':>11}{'bound |coef|, fit scale':>26}{'bound |coef|, venue floor':>27}")
for name, fn in CANDS.items():
    for lam in LAMS:
        rows = DATA[lam]; F = FITS[(lam, 'M2')]
        f = np.array([fn(r) for r in rows])
        c, *_ = np.linalg.lstsq(F['A'], f, rcond=None)
        fp = f - F['A'] @ c
        nrm = float(np.linalg.norm(fp))
        chibar = float(np.mean([r['chi'] for r in rows]))
        if nrm < 1e-9:
            say(f"  {name:<38}{lam:>5.2f}{nrm:>11.2e}{'COLLINEAR with M2 -- no bound':>53}")
            continue
        b_fit = 2.0 * F['s'] / nrm
        b_ven = 2.0 * (FLOOR_V / chibar) * np.sqrt(len(rows)) / nrm
        say(f"  {name:<38}{lam:>5.2f}{nrm:>11.4f}{b_fit:>26.5f}{b_ven:>27.5f}")
say("")
say("  For comparison the FITTED leading coefficients are of order: gamma ~ 1.2-1.6,")
say("  delta ~ 0.3-0.5, beta ~ 0.1-0.2, and the three-record coefficients are listed above.")

json.dump({str(l): dict(names=[r['name'] for r in DATA[l]],
                        chi=[float(r['chi']) for r in DATA[l]],
                        M1=list(map(float, FITS[(l,'M1')]['resid'])),
                        M2=list(map(float, FITS[(l,'M2')]['resid'])))
           for l in LAMS},
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b6_hierarchy.json'), 'w'), indent=0)
say("")
say(f"  elapsed {time.time()-t0:.1f}s")
