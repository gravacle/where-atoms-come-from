#!/usr/bin/env python3
"""REFUTER A2 -- ATTACK 1: CLAUSE (c) UNDER ORDINARY BIASED DATA.

The orientation prediction (v2_design.md section 3) predicts screening for
'random AND DC-free-coded data'. Clause (c) as registered (section 4) and as
implemented (measure_ori) carries NO condition on the declared DATA pattern's
balance: the guard-scope table lists only the DC positive control, the U-U
control band, and the DCF voider. The C4 text names 'the declared write
pattern itself, verified from the reader's own write record' as the
orientation analogue of the density condition -- but registers no condition,
no tolerance, and no code path reads the declared pattern.

A physicist's own data is ordinarily biased (real payloads carry DC content;
the program's own C-72 cell says orientation encoding 'accumulate[s] only
under all-one-way writing' -- i.e. partial one-way bias accumulates in
proportion). Under the model's own physics a biased DATA class accumulates:
that is the ENCODING WORKING AS CLAIMED. If clause (c) fires there, the
falsifier fires on correct physics.

Sweep: data_p1 in {0.50, 0.52, 0.55, 0.60, 0.65, 0.75}, 40 seeds each,
through the SEALED pipeline's own measure_ori (one implementation), all
guards as registered. Also the same sweep through my own independent
orientation rebuild (below) to show it is physics, not a pipeline artifact.
"""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # LANE_T50/V002
import v2_pipeline as V

MASTER = np.random.SeedSequence(777001)

print("REFUTER A2 -- ATTACK 1: clause (c) vs ordinary biased data")
print("sealed pipeline measure_ori, all registered guards; 40 seeds per bias")
print("bias b = 2*p1 - 1 (fraction of one-way excess in the reader's own data)")
print("")
for p1 in (0.50, 0.52, 0.55, 0.60, 0.65, 0.75):
    sts = {}; fires = 0; xis = []; bdus = []
    for ss in MASTER.spawn(40):
        rng = V.rng_from(ss)
        vals, classes, _ = V.build_ori(rng, data_p1=p1)
        m = V.measure_ori(vals, classes, rng)
        sts[m["state"]] = sts.get(m["state"], 0) + 1
        fires += m.get("fire_c", False)
        if m["state"] == "OK":
            xis.append(m["xi"]); bdus.append(m["bDU"])
    print("  p1=%.2f (b=%.2f): states %s | fire_c %d/40 | xi med %s | bDU med %s" %
          (p1, 2 * p1 - 1, sts, fires,
           ("%+.4f" % float(np.median(xis))) if xis else "--",
           ("%+.4f" % float(np.median(bdus))) if bdus else "--"))
print("")
print("Independent rebuild (my own Kerr model and fits, no lane code):")

SECT = 4096; GRID = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
GR = 64; PSFS = 0.7; RN = 0.05; MARGIN = 4
CYC = ["DATA", "U", "DCF", "U", "DC", "U"]; NS = 252

def psf():
    xs = np.arange(-8, 9, dtype=float); k = np.exp(-0.5 * (xs / PSFS) ** 2)
    return k / k.sum()

def my_build(rng, p1):
    cls = [CYC[s % 6] for s in range(NS)]
    M = np.empty((NS, SECT))
    for s, c in enumerate(cls):
        if c == "DATA":
            M[s] = ((rng.random(SECT) < p1) * 2 - 1).astype(float)
        elif c == "DCF":
            M[s] = np.tile([1.0, -1.0], SECT // 2)
        elif c == "DC":
            M[s] = 1.0
        else:
            M[s] = (2.0 * rng.binomial(GR, 0.5, SECT) - GR) / GR
    flat = np.convolve(M.reshape(-1), psf(), mode="same")
    v = flat.reshape(NS, SECT) + rng.normal(0, RN, (NS, SECT))
    return v, cls

def fitll(Ns, A):
    x = np.log10(np.array(Ns, float)); y = np.log10(np.array(A, float))
    xm, ym = x.mean(), y.mean(); Sxx = ((x - xm) ** 2).sum()
    b = ((x - xm) * (y - ym)).sum() / Sxx
    r = y - (ym + b * (x - xm))
    return b, math.sqrt(max((r ** 2).sum(), 0) / max(len(x) - 2, 1) / Sxx)

def pairs_of(cls, kind):
    used = set(); out = []
    for i in range(len(cls) - 1):
        if i in used or i + 1 in used:
            continue
        a, b = cls[i], cls[i + 1]
        if a == kind and b == "U":
            out.append((i, i + 1)); used |= {i, i + 1}
        elif b == kind and a == "U":
            out.append((i + 1, i)); used |= {i, i + 1}
    return out

def upairs(cls):
    ui = [i for i, c in enumerate(cls) if c == "U"]
    return [(ui[j], ui[j + 1]) for j in range(0, len(ui) - 1, 2)]

def ladder(cs, prs, rng):
    A = []
    for N in GRID:
        lo, hi = MARGIN, SECT - N - MARGIN
        if hi < lo:
            lo, hi = 0, 0
        ds = []
        for (sa, sb) in prs:
            for _ in range(4):
                oa = int(rng.integers(lo, hi + 1)); ob = int(rng.integers(lo, hi + 1))
                ds.append((cs[sa, oa + N] - cs[sa, oa]) - (cs[sb, ob + N] - cs[sb, ob]))
        A.append(float(np.median(np.abs(ds))))
    return np.array(A)

M2 = np.random.SeedSequence(777002)
for p1 in (0.55, 0.60, 0.75):
    fires = 0; xis = []
    for ss in M2.spawn(20):
        rng = np.random.default_rng(ss)
        v, cls = my_build(rng, p1)
        cs = np.concatenate([np.zeros((NS, 1)), np.cumsum(v, axis=1)], axis=1)
        A_uu = ladder(cs, upairs(cls), rng)
        A_du = ladder(cs, pairs_of(cls, "DATA"), rng)
        A_fu = ladder(cs, pairs_of(cls, "DCF"), rng)
        A_cu = ladder(cs, pairs_of(cls, "DC"), rng)
        bDU, seDU = fitll(GRID, A_du); bUU, seUU = fitll(GRID, A_uu)
        bFU, _ = fitll(GRID, A_fu)
        # guard on DC control as the pipeline does
        i0 = None
        for j in range(0, len(GRID) - 6 + 1):
            if A_cu[j] > 8.0 * A_uu[j]:
                i0 = j; break
        if i0 is None:
            continue
        bCU, seCU = fitll(GRID[i0:], A_cu[i0:])
        if bCU < 0.9 or not (0.35 <= bUU <= 0.65) or not (0.35 <= bFU <= 0.65):
            continue
        xi = bDU - bUU; sexi = math.hypot(seDU, seUU)
        xis.append(xi)
        fires += (xi - 2 * sexi >= 0.25) and (bCU - 2 * seCU >= 0.9)
    print("  p1=%.2f: my fire_c %d/20 | xi med %s" %
          (p1, fires, ("%+.4f" % float(np.median(xis))) if xis else "--"))
print("")
print("END ATTACK 1")
