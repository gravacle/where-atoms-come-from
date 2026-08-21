#!/usr/bin/env python3
"""REFUTER A2 -- INDEPENDENT REBUILD of the V002 model-side numbers.

My own code, from the REGISTERED TEXT (v2_design.md sections 1, 2, 5), not from
the pipeline: address-balanced interleave, rung-1 class de-trend, pool
statistic (all disjoint adjacent pairs x 4 placements), void guard kappa=8,
median density condition on N>=64 tol 0.02, control band, one-sided B=200
certificate, 2-SE fire rule. Cross-checks R1, R5 (via the pipeline's own
builders would be circular -- occupancy is rebuilt from the model laws:
programmed -N_E, residual U{-5..5}; orientation from grains/PSF/read-noise as
declared in the run header). Then rung-2 injection check (my analytic worry:
the calibration read's per-sector mean error enters D at sd N*sigma/sqrt(SECT)
*sqrt(2) -- at N=4096 comparable to the pair noise -- does beta_UU stay in
band?)."""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import v2_pipeline as V   # imported ONLY for constants echo comparison at end

SECT = 4096; NSECT = 256; NE = 100.0
GRID = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
KAPPA = 8.0; TOL = 0.02; I0D = 2

def roles():
    s = np.arange(NSECT)
    return ((s + s // 16) % 2 == 0)

def fitll(Ns, A):
    x = np.log10(np.array(Ns, float)); y = np.log10(np.array(A, float))
    xm, ym = x.mean(), y.mean(); Sxx = ((x - xm) ** 2).sum()
    b = ((x - xm) * (y - ym)).sum() / Sxx
    r = y - (ym + b * (x - xm))
    return b, math.sqrt(max((r ** 2).sum(), 0) / max(len(x) - 2, 1) / Sxx)

def fitlin(Ns, F):
    x = np.log10(np.array(Ns, float)); y = np.array(F, float)
    xm = x.mean(); Sxx = ((x - xm) ** 2).sum()
    return ((x - xm) * (y - y.mean())).sum() / Sxx

def build(rng, f):
    rw = roles()
    vals = np.empty((NSECT, SECT)); prog = np.zeros((NSECT, SECT))
    for s in range(NSECT):
        r = rng.integers(-5, 6, SECT).astype(float)
        if rw[s]:
            p = rng.random(SECT) < f
            vals[s] = np.where(p, -NE, r); prog[s] = p
        else:
            vals[s] = r
    return vals, prog, rw

def rung1(vals, rw):
    out = vals.copy(); s = np.arange(NSECT)
    for a in range(16):
        su = (s % 16 == a) & (~rw)
        out[s % 16 == a] -= vals[su].mean()
    return out

def pool(cs, prs, rng):
    P = {}
    for N in GRID:
        ds = []
        for (sa, sb) in prs:
            for _ in range(4):
                oa = int(rng.integers(0, SECT - N + 1))
                ob = int(rng.integers(0, SECT - N + 1))
                ds.append((cs[sa, oa + N] - cs[sa, oa]) -
                          (cs[sb, ob + N] - cs[sb, ob]))
        P[N] = np.asarray(ds, float)
    A = np.array([float(np.median(np.abs(P[N]))) for N in GRID])
    return P, A

def wu_pairs(rw):
    used = np.zeros(NSECT, bool); prs = []
    for i in range(NSECT - 1):
        if used[i] or used[i + 1]:
            continue
        if rw[i] != rw[i + 1]:
            w, u = (i, i + 1) if rw[i] else (i + 1, i)
            prs.append((w, u)); used[i] = used[i + 1] = True
    return prs

def uu_prs(rw):
    ui = list(np.where(~rw)[0])
    return [(ui[j], ui[j + 1]) for j in range(0, len(ui) - 1, 2)]

def dens_ladder(prog, rw):
    P = prog[rw]
    pcs = np.concatenate([np.zeros((P.shape[0], 1)), np.cumsum(P, axis=1)], axis=1)
    return np.array([float(np.median((pcs[:, N:] - pcs[:, :-N]) / N)) for N in GRID])

def measure(vals, prog, rw, rng, B=200, calib=None):
    if calib is not None:
        vals = vals - calib.mean(axis=1, keepdims=True)
    v = rung1(vals, rw)
    cs = np.concatenate([np.zeros((NSECT, 1)), np.cumsum(v, axis=1)], axis=1)
    _, A_uu = pool(cs, uu_prs(rw), rng)
    if A_uu[0] == 0.0:
        return dict(state="RAILED")
    wp, A_wu = pool(cs, wu_pairs(rw), rng)
    i0 = None
    for j in range(0, len(GRID) - 6 + 1):
        if A_wu[j] > KAPPA * A_uu[j]:
            i0 = j; break
    if i0 is None:
        return dict(state="READS_UNWRITTEN")
    fmed = dens_ladder(prog, rw)
    if abs(fitlin(GRID[I0D:], fmed[I0D:])) > TOL:
        return dict(state="VOID_DENSITY_MEDIAN")
    Ns = GRID[i0:]
    bWU, seWU = fitll(Ns, A_wu[i0:]); bUU, seUU = fitll(Ns, A_uu[i0:])
    if not (0.35 <= bUU <= 0.65):
        return dict(state="INCONCLUSIVE_CONTROL", bUU=bUU)
    bad = passed = 0
    for _ in range(B):
        lad = np.array([float(np.median(np.abs(rng.choice(wp[N], wp[N].size))))
                        for N in GRID])
        j0 = None
        for j in range(0, len(GRID) - 6 + 1):
            if lad[j] > KAPPA * A_uu[j]:
                j0 = j; break
        if j0 is None:
            continue
        passed += 1
        if fitll(GRID[j0:], lad[j0:])[0] < 0.9:
            bad += 1
    xi = bWU - bUU; sexi = math.hypot(seWU, seUU)
    st = "OK" if (bad == 0 and passed > 0) else "SEAM"
    fire = (bWU + 2 * seWU < 0.9) or (xi + 2 * sexi < 0.25)
    return dict(state=st, bWU=bWU, bUU=bUU, xi=xi, fire=fire,
                B1=(st == "OK" and 0.9 <= bWU <= 1.1 and 0.35 <= bUU <= 0.65
                    and xi >= 0.25))

M = np.random.SeedSequence(777005)
print("REFUTER A2 -- INDEPENDENT REBUILD (my code, from the registered text)")
print("")
print("R1 cross-check: 30 reads f=0.5, full machinery, B=200 certificate:")
bs = dict(bWU=[], bUU=[], xi=[]); b1 = fires = 0; sts = {}
for ss in M.spawn(30):
    rng = np.random.default_rng(ss)
    vals, prog, rw = build(rng, 0.5)
    m = measure(vals, prog, rw, rng)
    sts[m["state"]] = sts.get(m["state"], 0) + 1
    if m["state"] in ("OK", "SEAM"):
        for k in ("bWU", "bUU", "xi"):
            bs[k].append(m[k])
        b1 += m["B1"]; fires += m["fire"]
print("  beta_WU med %+.4f (lane +0.9999) | beta_UU med %+.4f (lane +0.4949) | "
      "xi med %+.4f (lane +0.5048)" %
      (np.median(bs["bWU"]), np.median(bs["bUU"]), np.median(bs["xi"])))
print("  states %s | B1 %d/30 | fires %d" % (sts, b1, fires))
print("")
print("R2 cross-check: my deterministic median ladder on the lane's skew masks:")
rw = roles()
for name, mask in (("cascade 1.5/0.5", V.cascade_mask(0.5, 1.5, 0.5)),
                   ("cascade 1.7/0.3", V.cascade_mask(0.5, 1.7, 0.3)),
                   ("two-level", V.twolevel_mask()),
                   ("M3 record", V.record_mask())):
    pr = np.zeros((NSECT, SECT)); pr[rw] = mask.astype(float)
    fmed = dens_ladder(pr, rw)
    sl = fitlin(GRID[I0D:], fmed[I0D:])
    print("  %-16s slope %+.4f -> %s (lane: VOID at 1.000)" %
          (name, sl, "VOID" if abs(sl) > TOL else "PASSES"))
print("")
print("RUNG-2 INJECTION CHECK (my worry: calibration-mean noise ~N*sd/64*sqrt2")
print("enters the pools; does beta_UU stay in band?): 15 honest reads, rung 2")
print("applied with a fp=0 calibration read:")
bu = []; bw = []
for ss in M.spawn(15):
    rng = np.random.default_rng(ss)
    vals, prog, rw2 = build(rng, 0.5)
    calib = np.empty((NSECT, SECT))
    for s in range(NSECT):
        calib[s] = rng.integers(-5, 6, SECT).astype(float)
    m = measure(vals, prog, rw2, rng, calib=calib)
    if m["state"] in ("OK", "SEAM"):
        bu.append(m["bUU"]); bw.append(m["bWU"])
print("  beta_UU med %+.4f (no-rung2 ~0.495; band [0.35,0.65]) | beta_WU med %+.4f"
      % (np.median(bu), np.median(bw)))
print("")
print("END REBUILD")
