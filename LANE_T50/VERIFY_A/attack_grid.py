#!/usr/bin/env python3
"""REFUTER A -- final escalation of ATTACK A.

The registered text leaves the GRID choice to the reader ("geometric N-grid >= 8 points
spanning >= 1.5 decades") and the SECTOR LENGTH to the reader's preparation.  The
median->mean pinch at N = sector length is what saved clause (b) at the full grid: the
top grid point is pinned to the sector MEAN density, flattening the fitted slope.  A
reader whose grid stops at N = 1024 (16..1024: 13 sqrt(2)-spaced points, 1.8 decades,
all >= 8 points and >= 1.5 decades) never sees the pinch.  Same part, same read, same
registered rules.  Patterns: cascade literals at hi = 1.6 and 1.7.
"""
import math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import sys; sys.path.insert(0, HERE)
from verify_a import (SECT, DENS_TOL, BAND_ACC, BAND_CTL,
                      fit_loglog, fit_lin_vs_logN, gen_uniform, gen_cascade_fixed,
                      build_occ, adj_pairs, same_pairs, csum, guard_i0, med_iqr)
from attack_precision import block_tables_med

GRID13 = [16, 23, 32, 45, 64, 91, 128, 181, 256, 362, 512, 724, 1024]
MASTER = np.random.SeedSequence(31337)

def measure_pair(vals, prog, classes, rng, K, grid):
    cs = csum(vals); pcs = csum(prog)
    uuA, _, _ = block_tables_med(cs, same_pairs(classes, "U"), rng, K, grid)
    res = {}
    for w in ("W1", "W2"):
        wuA, fhm, fhd = block_tables_med(cs, adj_pairs(classes, w, "U"), rng, K, grid,
                                         prog_cs=pcs)
        i0 = guard_i0(wuA, uuA, grid, min_pts=6)
        if i0 is None:
            return dict(state="READS_UNWRITTEN")
        g = slice(i0, None); Ns = grid[i0:]
        b, se = fit_loglog(Ns, wuA[g]); bU, seU = fit_loglog(Ns, uuA[g])
        if b is None or bU is None:
            return dict(state="INCONCLUSIVE_ZERO_A")
        ds_mean, _ = fit_lin_vs_logN(Ns, fhm[g])
        if abs(ds_mean) > DENS_TOL:
            return dict(state="VOID_DENSITY")
        if not (BAND_CTL[0] <= bU <= BAND_CTL[1]):
            return dict(state="INCONCLUSIVE_CONTROL")
        res[w] = (b, se)
    b1, s1 = res["W1"]; b2, s2 = res["W2"]
    d = abs(b1 - b2); sed = math.hypot(s1, s2)
    return dict(state="OK", b1=b1, b2=b2, delta=d, sed=sed,
                fire_b=(d - 2 * sed) > 0.2, above1=b1 > BAND_ACC[1])

def main():
    out = []
    def P(s=""):
        out.append(s); print(s, flush=True)
    P("REFUTER A -- GRID-CHOICE ESCALATION.  Registered-compliant reader: 1024-sector")
    P("part, K=256 pairs, 13-point sqrt(2) grid N=16..1024 (1.8 decades -- the registered")
    P("text requires only >=8 points and >=1.5 decades and leaves the grid to the reader).")
    P("Cascade literal vs uniform f=0.5, SAME part SAME read; all registered non-fire")
    P("rules applied.  master seed SeedSequence(31337)")
    P("")
    prng = np.random.default_rng(np.random.SeedSequence(777))
    cas16, f16 = gen_cascade_fixed(prng, levels=8, hi=1.6, f0=0.25)
    cas17, f17 = gen_cascade_fixed(prng, levels=8, hi=1.7, f0=0.25)
    P("declared literal patterns: cascade(1.6/0.4) f=%.4f | cascade(1.7/0.3) f=%.4f" %
      (f16, f17))
    P("")
    for tag, gpat, nseeds in (("cascade 1.6/0.4 vs uniform", cas16, 30),
                              ("cascade 1.7/0.3 vs uniform", cas17, 30)):
        fires = okc = above = 0; deltas, margins, b1s = [], [], []
        states = {}
        for ss in MASTER.spawn(nseeds):
            rng = np.random.default_rng(ss)
            v, p, c = build_occ(rng, {"W1": gpat, "W2": gen_uniform(0.5)},
                                ["W1", "U", "W2", "U"], nsect=1024)
            m = measure_pair(v, p, c, rng, K=256, grid=GRID13)
            states[m["state"]] = states.get(m["state"], 0) + 1
            if m["state"] != "OK": continue
            okc += 1; deltas.append(m["delta"]); margins.append(m["delta"] - 2 * m["sed"] - 0.2)
            b1s.append(m["b1"])
            if m["fire_b"]: fires += 1
            if m["above1"]: above += 1
        P("%s:" % tag)
        P("   REGISTERED CLAUSE (b) FIRES %d/%d OK reads (of %d seeds); states %s" %
          (fires, okc, nseeds, states))
        if deltas:
            dm, dl, dh = med_iqr(deltas); mm, ml, mh = med_iqr(margins)
            bm, bl, bh = med_iqr(b1s)
            P("   |delta beta| med %+0.4f IQR[%+0.4f,%+0.4f]" % (dm, dl, dh))
            P("   fire margin  med %+0.4f IQR[%+0.4f,%+0.4f]" % (mm, ml, mh))
            P("   skew beta_WU med %+0.4f IQR[%+0.4f,%+0.4f]; band-above %d/%d" %
              (bm, bl, bh, above, okc))
        P("")
    P("END OF GRID ESCALATION")
    with open(os.path.join(HERE, "attack_grid_run.txt"), "w") as fh:
        fh.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
