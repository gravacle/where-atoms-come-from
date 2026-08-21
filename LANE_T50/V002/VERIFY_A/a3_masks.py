#!/usr/bin/env python3
"""REFUTER A2 -- ATTACK 2: DENSITY-SHAPE CONFOUND IN THE UNCHECKED WINDOW.

The C1 condition pins the MEDIAN per-block f-hat ladder -- but only on N >= 64
(DENS_MED_I0 = 2), while the exponent fit starts at the guard point, usually
N = 16. The two low grid points are constrained by NOTHING except the void
guard at the single first-admissible point. New declared literal shapes:

  PERIODIC ISLAND MASKS. Period P divides 64 -> every window of N >= P carries
  EXACTLY the same programmed count (deterministic), so the checked ladder
  (N >= 64) is EXACTLY flat: slope 0.0000, infinite margin. For P = 128 the
  N = 64 median is forced to f by the pairing symmetry c(o) + c(o+64) = c_128.
  The medians at N = 16 (and 32 for P = 128) are suppressed by the island
  geometry -> the first one/two FIT points sit low -> fitted beta_WU rises.

  The guard bounds the suppression at the FIRST admissible point only
  (guard_start checks a single i0); the forced-median symmetry at P/2 protects
  midpoints; so the reachable excursion is capped -- MEASURE THE CAP, and
  whether a certified (state OK) read lands ABOVE the band's upper edge:
  point_certificate counts only surrogates fitting BELOW BAND_ACC[0] -- the
  upper edge is uncertified.

All masks are declared literals (deterministic 0/1), one-carrier, written
identically to every W sector -- exactly the legal-pattern class of A-D1.
"""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import v2_pipeline as V

MASTER = np.random.SeedSequence(777003)
SECT = V.SECT

def island_mask(period, L, bg=()):  # deterministic literal
    m = np.zeros(SECT, bool)
    for start in range(0, SECT, period):
        m[start:start + L] = True
        for b in bg:
            m[start + b] = True
    return m

def hole_mask(period, L):           # complement: dense sea with holes
    return ~island_mask(period, L)

MASKS = [
    ("island P=64 L=18",  island_mask(64, 18)),
    ("island P=64 L=20",  island_mask(64, 20)),
    ("island P=128 L=50", island_mask(128, 50)),
    ("island P=128 L=52", island_mask(128, 52)),
    ("island P=128 L=54", island_mask(128, 54)),
    ("island P=128 L=60", island_mask(128, 60)),
    ("hole   P=64 L=17",  hole_mask(64, 17)),
]

print("REFUTER A2 -- ATTACK 2: periodic island literals vs the C1 window")
print("deterministic median ladder (V.dens_med_ladder), then 30 seeds through")
print("the sealed measure_occ; DENS_TOL_MED=%.2f on N>=%d" %
      (V.DENS_TOL_MED, V.GRID[V.DENS_MED_I0]))
print("")
roleW0 = V.occ_roles()
for name, mask in MASKS:
    progm = np.zeros((V.NSECT, SECT)); progm[roleW0] = mask.astype(float)
    fmed = V.dens_med_ladder(progm, roleW0)
    ds, _ = V.fit_lin(V.GRID[V.DENS_MED_I0:], fmed[V.DENS_MED_I0:])
    print("  %-18s f=%.4f  fmed ladder %s" %
          (name, mask.mean(), " ".join("%.4f" % x for x in fmed)))
    print("  %-18s checked-window slope %+.5f (%s)" %
          ("", ds, "PASSES" if abs(ds) <= V.DENS_TOL_MED else "VOID"))
    sts = {}; betas = []; okb = []; fires = 0; b1 = 0; cbad = []
    for ss in MASTER.spawn(30):
        rng = V.rng_from(ss)
        vals, prog, roleW, _, _ = V.build_occ(rng, mask=mask)
        m = V.measure_occ(vals, prog, roleW, rng)
        sts[m["state"]] = sts.get(m["state"], 0) + 1
        if m["state"] in ("OK", "SEAM"):
            betas.append(m["bWU"]); fires += m["fire_a"]; b1 += m["B1"]
            if m["state"] == "OK":
                okb.append(m["bWU"]); cbad.append(m["cert_bad"])
    print("  %-18s states %s | fire_a %d/30 | B1 %d | beta med %s | "
          "certified range [%s, %s]" %
          ("", sts, fires, b1,
           ("%+.4f" % float(np.median(betas))) if betas else "--",
           ("%+.4f" % min(okb)) if okb else "--",
           ("%+.4f" % max(okb)) if okb else "--"))
    print("")

# ---- clause (b): the strongest mask against uniform f=0.5 on the SAME part ----
print("CLAUSE (b): island literal (W1) vs uniform f=0.5 (W2), same part, same")
print("read, sealed machinery (measure_occ with wsel/shared as measure_occ2):")
best = island_mask(128, 52)
def build_two_class(rng, mask, f2=0.5):
    roleW = V.occ_roles(V.NSECT)
    vals = np.empty((V.NSECT, SECT)); prog = np.zeros((V.NSECT, SECT))
    wclass = np.empty(V.NSECT, dtype=object); wclass[:] = None
    widx = 0
    for s in range(V.NSECT):
        r = V.residual(rng, "uniform", SECT)
        if roleW[s]:
            wclass[s] = "W1" if widx % 2 == 0 else "W2"
            p = mask.copy() if wclass[s] == "W1" else (rng.random(SECT) < f2)
            vals[s] = np.where(p, -float(V.N_E), r); prog[s] = p.astype(float)
            widx += 1
        else:
            vals[s] = r
    return vals, prog, roleW, wclass

fb = 0; dstats = []; st12 = {}
for ss in MASTER.spawn(30):
    rng = V.rng_from(ss)
    vals, prog, roleW, wclass = build_two_class(rng, best)
    v = V.detrend_addr(vals, roleW)
    cs = V.seccs(v)
    uidx = list(np.where(~roleW)[0])
    shared = V.uu_pool(cs, V.same_role_pairs(uidx), rng)
    res = {}
    for wcl in ("W1", "W2"):
        res[wcl] = V.measure_occ(vals, prog, roleW, rng, wsel=wcl,
                                 wclass=wclass, shared=shared)
    key = (res["W1"]["state"], res["W2"]["state"])
    st12[key] = st12.get(key, 0) + 1
    if all(res[w]["state"] == "OK" for w in res):
        d = abs(res["W1"]["bWU"] - res["W2"]["bWU"])
        se = math.hypot(res["W1"]["seWU"], res["W2"]["seWU"])
        dstats.append((d, d - 2 * se))
        fb += (d - 2 * se) > 0.2
print("  states (W1,W2): %s" % st12)
if dstats:
    ds_ = np.array(dstats)
    print("  |delta beta| med %.4f  fire margin med %+.4f  clause (b) FIRES %d/%d" %
          (float(np.median(ds_[:, 0])), float(np.median(ds_[:, 1])), fb, len(dstats)))
print("")
print("END ATTACK 2")
