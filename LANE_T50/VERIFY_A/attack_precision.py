#!/usr/bin/env python3
"""REFUTER A -- escalation of ATTACK A: a registered-COMPLIANT precision reader.

The registered text declares only MINIMUMS (K >= 8 pairs, >= 8 points, >= 1.5 decades)
and leaves part size, K, and grid density to the reader.  A reader with a 1024-sector
part, K = 256 disjoint adjacent pairs, and a 17-point sqrt(2)-spaced grid violates
nothing.  Against that reader we run the declared cascade(1.6/0.4) literal pattern and
uniform f=0.5 on the SAME part in the SAME read, and count registered clause (b) fires.
Also: prediction band-above rate, and the repair diagnostic (MEDIAN per-block f-hat
trend) on both skew and honest patterns.
"""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from verify_a import (SECT, NE, KAPPA, DENS_TOL, BAND_ACC, BAND_CTL, XI_MIN,
                      fit_loglog, fit_lin_vs_logN, gen_uniform, gen_cascade_fixed,
                      gen_twolevel_fixed, build_occ, adj_pairs, same_pairs, csum,
                      guard_i0, med_iqr)

GRID17 = [16, 23, 32, 45, 64, 91, 128, 181, 256, 362, 512, 724, 1024, 1448, 2048, 2896, 4096]
MASTER = np.random.SeedSequence(9191)

def block_tables_med(cs, pairs, rng, K, grid, prog_cs=None):
    A, FHmean, FHmed = [], [], []
    for N in grid:
        sel = rng.choice(len(pairs), size=min(K, len(pairs)), replace=False)
        ds, fs = [], []
        for i in sel:
            sA, sB = pairs[i]
            oA = int(rng.integers(0, SECT - N + 1)); oB = int(rng.integers(0, SECT - N + 1))
            ds.append((cs[sA, oA + N] - cs[sA, oA]) - (cs[sB, oB + N] - cs[sB, oB]))
            if prog_cs is not None:
                fs.append((prog_cs[sA, oA + N] - prog_cs[sA, oA]) / N)
        ds = np.asarray(ds, float)
        A.append(float(np.median(np.abs(ds))))
        FHmean.append(float(np.mean(fs)) if fs else float("nan"))
        FHmed.append(float(np.median(fs)) if fs else float("nan"))
    return np.array(A), np.array(FHmean), np.array(FHmed)

def measure_pair_precision(vals, prog, classes, rng, K, grid):
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
        ds_med, _ = fit_lin_vs_logN(Ns, fhd[g])
        if abs(ds_mean) > DENS_TOL:
            return dict(state="VOID_DENSITY")
        if not (BAND_CTL[0] <= bU <= BAND_CTL[1]):
            return dict(state="INCONCLUSIVE_CONTROL")
        res[w] = (b, se, ds_mean, ds_med)
    b1, s1, dm1, dd1 = res["W1"]; b2, s2, dm2, dd2 = res["W2"]
    d = abs(b1 - b2); sed = math.hypot(s1, s2)
    return dict(state="OK", b1=b1, se1=s1, b2=b2, se2=s2, delta=d, sed=sed,
                fire_b=(d - 2 * sed) > 0.2, dmean1=dm1, dmed1=dd1, dmean2=dm2, dmed2=dd2,
                above1=b1 > BAND_ACC[1])

def main():
    out = []
    def P(s=""):
        out.append(s); print(s, flush=True)
    P("REFUTER A -- PRECISION-READER ESCALATION (registered-compliant: K=256 of the")
    P("512 disjoint W1-U pairs on a 1024-sector part, 17-point sqrt(2) grid 16..4096,")
    P("2.4 decades; every declared analysis rule of the registered text applied,")
    P("including the control-band and density non-fire rules).")
    P("master seed SeedSequence(9191)")
    P("")
    prng = np.random.default_rng(np.random.SeedSequence(777))
    cas15, _ = gen_cascade_fixed(prng, levels=8, hi=1.5, f0=0.25)
    cas16, _ = gen_cascade_fixed(prng, levels=8, hi=1.6, f0=0.25)
    two1, _ = gen_twolevel_fixed(prng, 0.03, 0.95, 0.25)

    for tag, gpat, nseeds in (("cascade 1.6/0.4 vs uniform f=0.5", cas16, 30),
                              ("cascade 1.5/0.5 vs uniform f=0.5", cas15, 30),
                              ("two-level 0.03/0.95 vs uniform", two1, 30)):
        fires = okc = above = 0; deltas, margins, b1s = [], [], []
        states = {}
        dmeds1 = []
        for ss in MASTER.spawn(nseeds):
            rng = np.random.default_rng(ss)
            v, p, c = build_occ(rng, {"W1": gpat, "W2": gen_uniform(0.5)},
                                ["W1", "U", "W2", "U"], nsect=1024)
            m = measure_pair_precision(v, p, c, rng, K=256, grid=GRID17)
            states[m["state"]] = states.get(m["state"], 0) + 1
            if m["state"] != "OK": continue
            okc += 1; deltas.append(m["delta"]); margins.append(m["delta"] - 2 * m["sed"] - 0.2)
            b1s.append(m["b1"]); dmeds1.append(m["dmed1"])
            if m["fire_b"]: fires += 1
            if m["above1"]: above += 1
        P("%s:" % tag)
        P("   registered clause (b) FIRES %d/%d OK reads (of %d seeds); states %s" %
          (fires, okc, nseeds, states))
        if deltas:
            dm, dl, dh = med_iqr(deltas)
            mm, ml, mh = med_iqr(margins)
            bm, bl, bh = med_iqr(b1s)
            P("   |delta beta| med %+0.4f IQR[%+0.4f,%+0.4f]; fire margin med %+0.4f "
              "IQR[%+0.4f,%+0.4f]" % (dm, dl, dh, mm, ml, mh))
            P("   skew-pattern beta_WU med %+0.4f IQR[%+0.4f,%+0.4f]; prediction band-above"
              " %d/%d" % (bm, bl, bh, above, okc))
            P("   [repair diagnostic] MEDIAN per-block f-hat slope on log10 N, skew class:"
              " med %+0.4f" % float(np.median(dmeds1)))
        P("")

    P("repair-diagnostic control -- honest patterns, median f-hat slope (should pass a")
    P("0.05 tolerance):")
    for tag, g in (("uniform f=0.5", gen_uniform(0.5)), ("uniform f=0.1", gen_uniform(0.1))):
        sl = []
        for ss in MASTER.spawn(20):
            rng = np.random.default_rng(ss)
            v, p, c = build_occ(rng, {"W1": g, "W2": gen_uniform(0.5)},
                                ["W1", "U", "W2", "U"], nsect=1024)
            cs = csum(v); pcs = csum(p)
            wuA, fhm, fhd = block_tables_med(cs, adj_pairs(c, "W1", "U"), rng, 256, GRID17,
                                             prog_cs=pcs)
            sl.append(fit_lin_vs_logN(GRID17, fhd)[0])
        P("   %-16s median f-hat slope med %+0.4f  max|.| %0.4f" %
          (tag, float(np.median(sl)), float(np.max(np.abs(sl)))))
    P("")
    P("END OF PRECISION ESCALATION")
    with open(os.path.join(HERE, "attack_precision_run.txt"), "w") as fh:
        fh.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
