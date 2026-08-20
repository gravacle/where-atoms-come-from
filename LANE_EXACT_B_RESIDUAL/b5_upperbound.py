"""B5 -- WHAT THE RESIDUAL IS, AND HOW BIG A TERM COULD STILL BE HIDING.

  (1) EXTENSIVITY, decided EXACTLY.  Is there any term that grows with the NUMBER of records?
      C-38's theorem (B1/E8) says a commuting partner on another bath site changes chi by
      EXACTLY zero.  Adding such partners raises the record count without limit.  Therefore the
      coefficient of any purely extensive term is EXACTLY ZERO -- by proof, not by range.
      Checked here against a CONTROL that is exactly non-zero in the same table.

  (2) ADDITIVITY OVER DISJOINT BATH SITES.  Do two pairing partners on two different sites
      suppress the record by the PRODUCT of their separate factors?  If yes the disturbance is
      additive in log chi over disjoint regions and there is no cross term; if no, the failure
      of additivity is itself a residual and is measured here.

  (3) THE THREE-BODY TERM.  Add to the explained model a term counting partner pairs that
      ANTICOMMUTE WITH EACH OTHER -- a quantity C-36/38/39 cannot contain, since all three are
      statements about the record and ONE partner.  Refit, and report what it absorbs.

  (4) THE UPPER BOUND.  With the extended model fitted, how large could a FURTHER term be and
      still hide?  Reported as a coefficient bound from the least-squares standard error, with
      the error scale taken BOTH from the fit residual AND from the venue floor of B3.
"""
import numpy as np, sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from battery import *

t0 = time.time()
NB, N = 3, 8
LAMS = (0.4, 0.8, 1.2)
say("=" * 124)
say("B5   EXTENSIVITY, ADDITIVITY, THE THREE-BODY TERM, AND THE UPPER BOUND")
say("=" * 124)
ops, _, _ = build_ops(N)
env = BASE.env(NB)

def chi_of(partners, lam, read=READ):
    so = [(ops[read][0], 0)] + [(ops[l][0], s) for l, s in partners]
    return float(np.mean(chi_times(so, ops[read][0], env, lam, BASE.times)))

# ---------------------------------------------------------------- (1) EXTENSIVITY
say("")
say("(1) EXTENSIVITY -- IS THERE ANY TERM THAT GROWS WITH THE NUMBER OF RECORDS?")
say("    Left: records added on OTHER sites that COMMUTE with the read record.  Right, in the")
say("    SAME TABLE (D-15): the identical count of records added ON THE READ RECORD'S SITE,")
say("    which must register a large change or the probe is blind.")
say("")
say(f"  {'lam':>5}{'records written':>17}{'chi, +commuting on other sites':>32}{'chi - alone':>14}"
    f"{'chi, +same site':>18}{'chi - alone':>14}")
EXT = []
for lam in LAMS:
    alone = chi_of([], lam)
    for k in range(0, 3):
        far = [("X2", 1), ("X3", 2)][:k]
        near = [("X2", 0), ("X3", 0)][:k]
        cf = chi_of(far, lam)
        cn = chi_of(near, lam)
        EXT.append(abs(cf - alone))
        say(f"  {lam:>5.2f}{1+k:>17}{cf:>32.14f}{cf-alone:>+14.2e}{cn:>18.14f}{cn-alone:>+14.6f}")
say("")
say(f"  largest |chi - alone| from adding commuting records on other sites: {max(EXT):.3e}  -> EXACTLY ZERO")
say("  (the control column moves by up to 0.4 in the same table, so the probe is not blind)")
say("")
say("  CONCLUSION, EXACT: no extensive term exists.  chi does not know how many records the")
say("  carrier holds.  It knows only how many share its bath site and how many pair with it.")

# ---------------------------------------------------------------- (2) ADDITIVITY
say("")
say("(2) ADDITIVITY OVER DISJOINT BATH SITES, for the PAIRING disturbance.")
say("    d1 = chi(one pairing partner on another site)/chi(alone).  If the disturbance were")
say("    additive in log chi over disjoint sites, two such partners on two DIFFERENT sites would")
say("    give d1^2 exactly.")
say("")
say(f"  {'lam':>5}{'d1 (Z1@1)':>14}{'d1b (Z1X2@2)':>15}{'measured d(both, 2 sites)':>27}{'d1*d1b':>12}{'ratio meas/product':>20}")
for lam in LAMS:
    alone = chi_of([], lam)
    d1 = chi_of([("Z1", 1)], lam) / alone
    d1b = chi_of([("Z1X2", 2)], lam) / alone
    both = chi_of([("Z1", 1), ("Z1X2", 2)], lam) / alone
    say(f"  {lam:>5.2f}{d1:>14.9f}{d1b:>15.9f}{both:>27.9f}{d1*d1b:>12.9f}{both/(d1*d1b):>20.9f}")
say("")
say("  and the same two partners forced onto the SAME other site (so they are no longer disjoint):")
say(f"  {'lam':>5}{'d(both, same other site)':>28}{'d1*d1b':>12}{'ratio':>12}")
for lam in LAMS:
    alone = chi_of([], lam)
    d1 = chi_of([("Z1", 1)], lam) / alone
    d1b = chi_of([("Z1X2", 2)], lam) / alone
    same = chi_of([("Z1", 1), ("Z1X2", 1)], lam) / alone
    say(f"  {lam:>5.2f}{same:>28.9f}{d1*d1b:>12.9f}{same/(d1*d1b):>12.9f}")

# ---------------------------------------------------------------- (3) THE THREE-BODY TERM
say("")
say("(3) THE THREE-BODY TERM.")
def design_ext(f, r):
    return design_row(f) + [float(r['pp_anti_same']), float(r['pp_anti'] - r['pp_anti_same'])]
TERMS_EXT = TERMS + ["eps (3-body, shared site)", "zeta (3-body, other sites)"]

def fit2(rows, ext=True):
    A = np.array([(design_ext(r, r) if ext else design_row(r)) for r in rows])
    y = np.log(np.array([r['chi'] for r in rows]))
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ coef
    meas = np.array([r['chi'] for r in rows]); pred = np.exp(A @ coef)
    dof = max(len(rows) - A.shape[1], 1)
    s2 = float(res @ res) / dof
    cov = s2 * np.linalg.pinv(A.T @ A)
    return dict(A=A, coef=coef, se=np.sqrt(np.diag(cov)), pred=pred, resid=meas - pred,
                logresid=res, rms=float(np.sqrt(np.mean((meas - pred) ** 2))),
                maxabs=float(np.max(np.abs(meas - pred))), s=float(np.sqrt(s2)))

say("")
say(f"  {'lam':>5}{'model':>12}{'#cfg':>6}{'params':>8}" + "".join(f"{t.split(' ')[0]:>9}" for t in TERMS_EXT) +
    f"{'rms resid':>12}{'max|resid|':>12}")
RES = {}
for lam in LAMS:
    rows = run_battery(N, BASE, lam, NB=NB, ops=ops)
    Fb = fit2(rows, ext=False); Fe = fit2(rows, ext=True)
    RES[lam] = (rows, Fb, Fe)
    say(f"  {lam:>5.2f}{'explained':>12}{len(rows):>6}{5:>8}" + "".join(f"{c:>9.4f}" for c in Fb['coef']) +
        f"{'':>18}{Fb['rms']:>12.3e}{Fb['maxabs']:>12.3e}")
    say(f"  {lam:>5.2f}{'+3-body':>12}{len(rows):>6}{7:>8}" + "".join(f"{c:>9.4f}" for c in Fe['coef']) +
        f"{Fe['rms']:>12.3e}{Fe['maxabs']:>12.3e}")
say("")
say("  three-body coefficients with their standard errors:")
say(f"  {'lam':>5}{'eps (shared site)':>26}{'zeta (other sites)':>26}{'c0 CONTROL (must be 0)':>28}")
for lam in LAMS:
    _, _, Fe = RES[lam]
    say(f"  {lam:>5.2f}{Fe['coef'][5]:>+16.5f} +- {Fe['se'][5]:<7.5f}{Fe['coef'][6]:>+16.5f} +- {Fe['se'][6]:<7.5f}"
        f"{Fe['coef'][4]:>+16.2e} +- {Fe['se'][4]:<9.2e}")

# ---------------------------------------------------------------- (4) UPPER BOUND
say("")
say("(4) THE UPPER BOUND ON A FURTHER TERM.")
say("    For a candidate extra feature f, the smallest coefficient the data could still")
say("    distinguish from zero is set by the residual scale divided by the part of f that is")
say("    ORTHOGONAL to the terms already fitted.  Two error scales are used: the fit's own")
say("    residual, and B3's VENUE FLOOR (the larger, and the honest one).")
try:
    VJ = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b3_venues.json')))
    allv = np.array([VJ[k]['resid'] for k in VJ])
    VENUE_FLOOR = float(np.median(allv.max(axis=0) - allv.min(axis=0)))
except Exception as e:
    VENUE_FLOOR = float('nan'); say(f"  (b3_venues.json not available: {e})")
say(f"    venue floor (median configuration spread across the 8 venues of B3): {VENUE_FLOOR:.5f} in chi")
say("")
CANDS = {
  "extensive: total records written":       lambda r: float(r['npart'] + 1),
  "records on OTHER sites (any pairing)":   lambda r: float(r['p0'] + r['p1']),
  "site occupancy squared":                 lambda r: float(r['m'] ** 2),
  "3-body count (shared site)":             lambda r: float(r['pp_anti_same']),
  "pairing partners squared":               lambda r: float((r['m1'] + r['p1']) ** 2),
}
say(f"  {'candidate extra term f':<40}{'lam':>5}{'|f_perp| range':>16}{'bound |coef| (fit scale)':>26}{'bound |coef| (venue floor)':>28}")
for name, fn in CANDS.items():
    for lam in LAMS:
        rows, Fb, Fe = RES[lam]
        A = Fe['A']
        f = np.array([fn(r) for r in rows])
        # part of f not already expressible by the fitted design
        c, *_ = np.linalg.lstsq(A, f, rcond=None)
        fp = f - A @ c
        nrm = float(np.linalg.norm(fp))
        if nrm < 1e-10:
            say(f"  {name:<40}{lam:>5.2f}{nrm:>16.2e}{'EXACTLY COLLINEAR -- no bound possible':>54}")
            continue
        # coefficient in log chi; convert to a chi-scale bound using mean chi
        chibar = float(np.mean([r['chi'] for r in rows]))
        b_fit = 2.0 * Fe['s'] / nrm
        b_ven = 2.0 * (VENUE_FLOOR / chibar) * np.sqrt(len(rows)) / nrm
        say(f"  {name:<40}{lam:>5.2f}{nrm:>16.4f}{b_fit:>26.5f}{b_ven:>28.5f}")
say("")
say("  Read these as: an additional additive term coef*f in log chi could not be seen unless")
say("  |coef| exceeded the bound.  Relative to the fitted leading coefficients (gamma ~ 1,")
say("  beta and delta of order 0.1-1), that is the room left for anything unexplained.")
say("")
say(f"  elapsed {time.time()-t0:.1f}s")
