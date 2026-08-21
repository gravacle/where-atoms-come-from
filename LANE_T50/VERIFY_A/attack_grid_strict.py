#!/usr/bin/env python3
"""REFUTER A -- ATTACK A, strict-compliance re-run.

Self-audit: on a sqrt(2)-spaced grid the registered guard's '>= 1.5 decades' clause is
the binding admissibility rule (the '>= 6 points' parenthetical presumes the 2x grid).
Here the guard admits N_min ONLY if the surviving span log10(N_max/N_min) >= 1.5, and
the read is declared unwritten otherwise.  The actual N_min used is reported per read.
Grids tried: 16..1024 (1.81 dec) and 16..1448 (1.96 dec), both sqrt(2)-spaced, both
registered-compliant.  Cascade literals hi = 1.65, 1.7 vs uniform f=0.5, same part,
same read; every registered non-fire rule applied.
"""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from verify_a import (SECT, DENS_TOL, BAND_ACC, BAND_CTL,
                      fit_loglog, fit_lin_vs_logN, gen_uniform, gen_cascade_fixed,
                      build_occ, adj_pairs, same_pairs, csum, med_iqr)
from attack_precision import block_tables_med

MASTER = np.random.SeedSequence(60609)

def guard_i0_strict(A_wu, A_uu, grid, kappa=8.0, decades=1.5):
    for i0 in range(len(grid)):
        if math.log10(grid[-1] / grid[i0]) < decades:
            return None
        if A_wu[i0] > kappa * A_uu[i0]:
            return i0
    return None

def measure_pair_strict(vals, prog, classes, rng, K, grid):
    cs = csum(vals); pcs = csum(prog)
    uuA, _, _ = block_tables_med(cs, same_pairs(classes, "U"), rng, K, grid)
    res = {}
    for w in ("W1", "W2"):
        wuA, fhm, fhd = block_tables_med(cs, adj_pairs(classes, w, "U"), rng, K, grid,
                                         prog_cs=pcs)
        i0 = guard_i0_strict(wuA, uuA, grid)
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
        res[w] = (b, se, grid[i0])
    b1, s1, n1 = res["W1"]; b2, s2, n2 = res["W2"]
    d = abs(b1 - b2); sed = math.hypot(s1, s2)
    return dict(state="OK", b1=b1, b2=b2, delta=d, sed=sed, Nmin1=n1, Nmin2=n2,
                fire_b=(d - 2 * sed) > 0.2, above1=b1 > BAND_ACC[1])

def sqrt2_grid(lo, hi):
    g = []; x = lo
    while x <= hi + 0.5:
        g.append(int(round(x))); x *= math.sqrt(2)
    return g

def main():
    out = []
    def P(s=""):
        out.append(s); print(s, flush=True)
    P("REFUTER A -- STRICT-COMPLIANCE RE-RUN OF THE GRID ESCALATION")
    P("guard admissibility: N_min only where log10(N_max/N_min) >= 1.5 (the registered")
    P("decades rule, binding on sqrt(2) grids); reads failing it are UNWRITTEN.")
    P("1024-sector part, K=256 of 512 disjoint pairs.  master seed SeedSequence(60609)")
    P("")
    prng = np.random.default_rng(np.random.SeedSequence(777))
    _cas16, _ = gen_cascade_fixed(prng, levels=8, hi=1.6, f0=0.25)   # keep draw order
    _cas17, _ = gen_cascade_fixed(prng, levels=8, hi=1.7, f0=0.25)
    prng2 = np.random.default_rng(np.random.SeedSequence(778))
    cas165, f165 = gen_cascade_fixed(prng2, levels=8, hi=1.65, f0=0.25)
    cas17b, f17b = gen_cascade_fixed(prng2, levels=8, hi=1.7, f0=0.25)
    P("declared literals: cascade(1.65/0.35) f=%.4f | cascade(1.7/0.3) f=%.4f" %
      (f165, f17b))
    P("")
    for glo, ghi in ((16, 1024), (16, 1448)):
        grid = sqrt2_grid(glo, ghi)
        P("GRID %s (%d points, %.2f decades):" %
          (grid, len(grid), math.log10(grid[-1] / grid[0])))
        for tag, gpat in (("cascade 1.65/0.35 vs uniform", cas165),
                          ("cascade 1.7/0.3  vs uniform", cas17b)):
            fires = okc = above = 0; deltas, margins, b1s, nmins = [], [], [], []
            states = {}
            for ss in MASTER.spawn(30):
                rng = np.random.default_rng(ss)
                v, p, c = build_occ(rng, {"W1": gpat, "W2": gen_uniform(0.5)},
                                    ["W1", "U", "W2", "U"], nsect=1024)
                m = measure_pair_strict(v, p, c, rng, K=256, grid=grid)
                states[m["state"]] = states.get(m["state"], 0) + 1
                if m["state"] != "OK": continue
                okc += 1; deltas.append(m["delta"])
                margins.append(m["delta"] - 2 * m["sed"] - 0.2)
                b1s.append(m["b1"]); nmins.append(m["Nmin1"])
                if m["fire_b"]: fires += 1
                if m["above1"]: above += 1
            P("  %-28s CLAUSE (b) FIRES %d/%d OK (of 30); states %s" %
              (tag, fires, okc, states))
            if deltas:
                dm, dl, dh = med_iqr(deltas); mm, ml, mh = med_iqr(margins)
                bm, bl, bh = med_iqr(b1s)
                P("     |delta| med %+0.4f IQR[%+0.4f,%+0.4f]; margin med %+0.4f "
                  "IQR[%+0.4f,%+0.4f]" % (dm, dl, dh, mm, ml, mh))
                P("     skew beta_WU med %+0.4f; band-above %d/%d; N_min used: %s" %
                  (bm, above, okc, sorted(set(nmins))))
        P("")
    P("END OF STRICT RE-RUN")
    with open(os.path.join(HERE, "attack_grid_strict_run.txt"), "w") as fh:
        fh.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
