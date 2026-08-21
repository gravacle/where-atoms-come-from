#!/usr/bin/env python3
"""REFUTER A -- COMPUTATION -- T-50 DESIGN ONE (THE SAME-READ CONTRAST EXPONENT).

Independent rebuild of the design's model-side numbers with refuter-owned code, then
attacks: baseline offsets (common-mode AND sector-differential), density confounds,
nominal-vs-actual errors, ensemble-vs-single-read statistics, and falsifier clauses
fired against correct physics.  Model laws taken from the sealed formation layer
(model/geometry.py): programmed NAND cell = -100 e, over-erase residual = uniform
integers {-5..+5} e per cell, orientation grain = +-1 (DC-saturated all +1).

This code implements the REGISTERED text, including the guarded non-fire rules the
lane's own fire_a omits (U-U control outside [0.35,0.65] => INCONCLUSIVE, never a fire),
so every fire reported here is unimpeachable under the design's own registered rules.

Writes nothing outside LANE_T50/VERIFY_A. No git, no reproduce.sh, no r3.sh.
"""
import math, os, sys, time
import numpy as np

HERE  = os.path.dirname(os.path.abspath(__file__))
REPO  = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, "model"))
import geometry as GE   # read-only, anchor only

SECT      = 4096
GRID9     = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
GRID8MIN  = [16, 32, 64, 128, 256, 512, 1024, 2048]      # registered minimum: 8 pts, >=1.5 dec
NE        = 100.0
KAPPA     = 8.0
DENS_TOL  = 0.05
MIN_PTS   = 6
BAND_ACC  = (0.9, 1.1)
BAND_CTL  = (0.35, 0.65)
XI_MIN    = 0.25

MASTER = np.random.SeedSequence(4242)   # refuter's own seed, distinct from the lane's

# ---------------------------------------------------------------- fitting
def fit_loglog(Ns, A):
    A = np.asarray(A, float)
    if np.any(A <= 0):
        return None, None            # undefined -> caller marks inconclusive
    x = np.log10(np.asarray(Ns, float)); y = np.log10(A)
    n = len(x); xm, ym = x.mean(), y.mean()
    Sxx = ((x - xm) ** 2).sum()
    b = ((x - xm) * (y - ym)).sum() / Sxx
    r = y - (ym + b * (x - xm))
    se = math.sqrt(max((r ** 2).sum(), 0.0) / max(n - 2, 1) / Sxx)
    return b, se

def fit_lin_vs_logN(Ns, F):
    x = np.log10(np.asarray(Ns, float)); y = np.asarray(F, float)
    n = len(x); xm, ym = x.mean(), y.mean()
    Sxx = ((x - xm) ** 2).sum()
    b = ((x - xm) * (y - ym)).sum() / Sxx
    r = y - (ym + b * (x - xm))
    se = math.sqrt(max((r ** 2).sum(), 0.0) / max(n - 2, 1) / Sxx)
    return b, se

def med_iqr(v):
    v = np.asarray(v, float)
    return float(np.median(v)), float(np.percentile(v, 25)), float(np.percentile(v, 75))

# ---------------------------------------------------------------- pattern generators
def gen_uniform(f):
    def g(rng): return rng.random(SECT) < f
    return g

def gen_alternating():
    def g(rng): return np.arange(SECT) % 2 == 0
    return g

def gen_fixed_record(f, cells=512):
    def g(rng):
        p = np.zeros(SECT, bool); p[:cells] = rng.random(cells) < f
        return p
    return g

def gen_twolevel(f_lo, f_hi, dense_frac):
    nd = int(SECT * dense_frac)
    def g(rng):
        p = np.full(SECT, f_lo); p[SECT - nd:] = f_hi
        return rng.random(SECT) < p
    return g

def cascade_probs(rng, levels, hi, f0):
    """Dyadic multiplicative cascade of programming probabilities over one sector.
    At each level every parent segment splits its density: one half x hi, other x (2-hi),
    which half chosen at random.  Mean density preserved (before capping at 1)."""
    p = np.full(SECT, f0)
    seg = SECT
    for _ in range(levels):
        seg //= 2
        npar = SECT // (seg * 2)
        pick = rng.integers(0, 2, npar)
        halves = np.where(pick[:, None] == 1, [hi, 2.0 - hi], [2.0 - hi, hi])
        p *= np.repeat(halves.reshape(-1), seg)
    return np.minimum(p, 1.0)

def gen_cascade_fixed(rng_pat, levels=8, hi=1.5, f0=0.25):
    """ONE declared literal pattern (drawn once), written identically to every W sector."""
    probs = cascade_probs(rng_pat, levels, hi, f0)
    mask = rng_pat.random(SECT) < probs
    def g(rng): return mask.copy()
    return g, float(mask.mean())

def gen_twolevel_fixed(rng_pat, f_lo, f_hi, dense_frac):
    nd = int(SECT * dense_frac)
    p = np.full(SECT, f_lo); p[SECT - nd:] = f_hi
    mask = rng_pat.random(SECT) < p
    def g(rng): return mask.copy()
    return g, float(mask.mean())

# ---------------------------------------------------------------- part builders
def build_occ(rng, wgens, cycle, nsect=256, mu=0.0, two_signed=False,
              offset=0.0, drift=0.0, sector_off_sd=0.0):
    """wgens: dict classname -> mask generator.  'U' handled internally.
    Returns vals (nsect,SECT), prog (nsect,SECT), classes list."""
    vals = np.empty((nsect, SECT)); prog = np.zeros((nsect, SECT)); classes = []
    for s in range(nsect):
        c = cycle[s % len(cycle)]; classes.append(c)
        r = rng.integers(-5, 6, SECT).astype(float) + mu
        if c == "U":
            vals[s] = r
        else:
            m = wgens[c](rng)
            if two_signed:
                sign = (rng.integers(0, 2, SECT) * 2 - 1).astype(float)
            else:
                sign = -np.ones(SECT)
            vals[s] = np.where(m, sign * NE, r); prog[s] = m
    if offset: vals += offset
    if drift:
        vals += np.linspace(-0.5, 0.5, nsect * SECT).reshape(nsect, SECT) * drift
    if sector_off_sd:
        vals += rng.normal(0.0, sector_off_sd, nsect)[:, None]
    return vals, prog, classes

def build_ori(rng, data="random", nsect=256, erase_bias=0.0, offset=0.0,
              tilt_deg=0.0, sector_off_sd=0.0, class_offset=0.0):
    """cycle [DATA, U, DC, U].  class_offset: additive per-cell read offset applied to
    the WRITTEN classes only (DATA and DC) -- the sector-differential/class-correlated
    baseline attack."""
    vals = np.empty((nsect, SECT)); classes = []
    cyc = ["DATA", "U", "DC", "U"]
    for s in range(nsect):
        c = cyc[s % 4]; classes.append(c)
        if c == "DATA":
            if data == "random":
                v = (rng.integers(0, 2, SECT) * 2 - 1).astype(float)
            else:                                    # dcfree
                v = np.tile([1.0, -1.0], SECT // 2)
        elif c == "DC":
            v = np.ones(SECT)
        else:
            v = (rng.integers(0, 2, SECT) * 2 - 1).astype(float) + erase_bias
        if tilt_deg:
            v = v * np.cos(rng.normal(0.0, math.radians(tilt_deg), SECT))
        if class_offset and c in ("DATA", "DC"):
            v = v + class_offset
        vals[s] = v
    if offset: vals += offset
    if sector_off_sd:
        vals += rng.normal(0.0, sector_off_sd, nsect)[:, None]
    return vals, classes

# ---------------------------------------------------------------- pairing + tables
def adj_pairs(classes, a, b):
    return [(i, i + 1) for i in range(len(classes) - 1)
            if classes[i] == a and classes[i + 1] == b]

def same_pairs(classes, a):
    idx = [i for i, c in enumerate(classes) if c == a]
    return [(idx[j], idx[j + 1]) for j in range(0, len(idx) - 1, 2)]

def block_tables(cs, pairs, rng, K, grid, prog_cs=None, aligned=False):
    A, M2, VAR, FH, RAW = [], [], [], [], []
    for N in grid:
        sel = rng.choice(len(pairs), size=min(K, len(pairs)), replace=False)
        ds, fs, raws = [], [], []
        for i in sel:
            sA, sB = pairs[i]
            if aligned:
                oA = oB = 0
            else:
                oA = int(rng.integers(0, SECT - N + 1)); oB = int(rng.integers(0, SECT - N + 1))
            a = cs[sA, oA + N] - cs[sA, oA]; b = cs[sB, oB + N] - cs[sB, oB]
            ds.append(a - b); raws.append(a)
            if prog_cs is not None:
                fs.append((prog_cs[sA, oA + N] - prog_cs[sA, oA]) / N)
        ds = np.asarray(ds, float)
        A.append(float(np.median(np.abs(ds)))); M2.append(float(np.median(ds ** 2)))
        VAR.append(float(np.var(ds, ddof=1))); RAW.append(float(np.median(np.abs(raws))))
        FH.append(float(np.mean(fs)) if fs else float("nan"))
    return dict(A=np.array(A), M2=np.array(M2), VAR=np.array(VAR),
                FH=np.array(FH), RAW=np.array(RAW))

def csum(vals):
    return np.concatenate([np.zeros((vals.shape[0], 1)), np.cumsum(vals, axis=1)], axis=1)

def guard_i0(A_wu, A_uu, grid, kappa=KAPPA, min_pts=MIN_PTS):
    for i0 in range(0, len(grid) - min_pts + 1):
        if A_wu[i0] > kappa * A_uu[i0]:
            return i0
    return None

# ---------------------------------------------------------------- occupancy measurement
def measure_occ(vals, prog, classes, rng, wname="W", K=16, grid=GRID9,
                kappa=KAPPA, aligned=False):
    cs = csum(vals); pcs = csum(prog)
    wu = block_tables(cs, adj_pairs(classes, wname, "U"), rng, K, grid, prog_cs=pcs, aligned=aligned)
    uu = block_tables(cs, same_pairs(classes, "U"), rng, K, grid, aligned=aligned)
    ww = block_tables(cs, same_pairs(classes, wname), rng, K, grid, aligned=aligned)
    i0 = guard_i0(wu["A"], uu["A"], grid, kappa)
    out = dict(wu=wu, uu=uu, ww=ww, i0=i0)
    if i0 is None:
        out.update(state="READS_UNWRITTEN"); return out
    Ns = grid[i0:]; g = slice(i0, None)
    bWU, seWU = fit_loglog(Ns, wu["A"][g]); bUU, seUU = fit_loglog(Ns, uu["A"][g])
    if bWU is None or bUU is None:
        out.update(state="INCONCLUSIVE_ZERO_A"); return out
    bWW = fit_loglog(Ns, ww["A"][g])[0]
    ds, dse = fit_lin_vs_logN(Ns, wu["FH"][g])
    xi = bWU - bUU; sexi = math.hypot(seWU, seUU)
    density_ok = abs(ds) <= DENS_TOL
    ctl_ok = BAND_CTL[0] <= bUU <= BAND_CTL[1]
    state = ("OK" if density_ok and ctl_ok else
             "VOID_DENSITY" if not density_ok else "INCONCLUSIVE_CONTROL")
    out.update(state=state, Nmin=grid[i0], bWU=bWU, seWU=seWU, bUU=bUU, seUU=seUU,
               bWW=bWW, xi=xi, sexi=sexi, dslope=ds,
               density_ok=density_ok, ctl_ok=ctl_ok,
               B1=density_ok and ctl_ok and BAND_ACC[0] <= bWU <= BAND_ACC[1] and xi >= XI_MIN,
               band_above=density_ok and ctl_ok and bWU > BAND_ACC[1],
               band_below=density_ok and ctl_ok and bWU < BAND_ACC[0],
               # registered clause (a): guard-passing, density-trend-free, control in band
               fire_a=density_ok and ctl_ok and
                      ((bWU + 2 * seWU < BAND_ACC[0]) or (xi + 2 * sexi < XI_MIN)),
               # lane's own (unregistered-guarded) version, for comparison:
               fire_a_lane=density_ok and
                      ((bWU + 2 * seWU < BAND_ACC[0]) or (xi + 2 * sexi < XI_MIN)))
    return out

def measure_occ_pair(vals, prog, classes, rng, K=16, grid=GRID9):
    """Two written classes W1,W2 against shared U -- registered clause (b)."""
    cs = csum(vals); pcs = csum(prog)
    uu = block_tables(cs, same_pairs(classes, "U"), rng, K, grid)
    res = {}
    for w in ("W1", "W2"):
        wu = block_tables(cs, adj_pairs(classes, w, "U"), rng, K, grid, prog_cs=pcs)
        i0 = guard_i0(wu["A"], uu["A"], grid)
        if i0 is None:
            return dict(state="READS_UNWRITTEN", which=w)
        Ns = grid[i0:]; g = slice(i0, None)
        b, se = fit_loglog(Ns, wu["A"][g])
        bU, seU = fit_loglog(Ns, uu["A"][g])
        if b is None or bU is None:
            return dict(state="INCONCLUSIVE_ZERO_A", which=w)
        ds, _ = fit_lin_vs_logN(Ns, wu["FH"][g])
        if abs(ds) > DENS_TOL:
            return dict(state="VOID_DENSITY", which=w)
        if not (BAND_CTL[0] <= bU <= BAND_CTL[1]):
            return dict(state="INCONCLUSIVE_CONTROL", which=w)
        res[w] = (b, se)
    (b1, s1), (b2, s2) = res["W1"], res["W2"]
    d = abs(b1 - b2); sed = math.hypot(s1, s2)
    return dict(state="OK", b1=b1, se1=s1, b2=b2, se2=s2, delta=d, sed=sed,
                fire_b=(d - 2 * sed) > 0.2)

# ---------------------------------------------------------------- orientation measurement
def measure_ori(vals, classes, rng, K=16, grid=GRID9):
    cs = csum(vals)
    dU = block_tables(cs, adj_pairs(classes, "DATA", "U"), rng, K, grid)
    cU = block_tables(cs, adj_pairs(classes, "DC", "U"), rng, K, grid)
    uu = block_tables(cs, same_pairs(classes, "U"), rng, K, grid)
    bDU, seDU = fit_loglog(grid, dU["A"])
    bUU, seUU = fit_loglog(grid, uu["A"])
    if bDU is None or bUU is None:
        return dict(state="INCONCLUSIVE_ZERO_A")
    i0 = guard_i0(cU["A"], uu["A"], grid)
    if i0 is None:
        bCU, seCU = None, None
    else:
        bCU, seCU = fit_loglog(grid[i0:], cU["A"][i0:])
    xi = bDU - bUU; sexi = math.hypot(seDU, seUU)
    ctl_ok = BAND_CTL[0] <= bUU <= BAND_CTL[1]
    dc_ok = (bCU is not None) and (bCU >= BAND_ACC[0])
    return dict(state="OK" if ctl_ok else "INCONCLUSIVE_CONTROL",
                bDU=bDU, seDU=seDU, bUU=bUU, seUU=seUU, bCU=bCU, seCU=seCU,
                xi=xi, sexi=sexi, ctl_ok=ctl_ok, dc_ok=dc_ok,
                B3=ctl_ok and BAND_CTL[0] <= bDU <= BAND_CTL[1] and xi < XI_MIN and dc_ok,
                fire_c=ctl_ok and (xi - 2 * sexi >= XI_MIN) and
                       (bCU is not None) and (seCU is not None) and
                       (bCU - 2 * seCU >= BAND_ACC[0]))

# ================================================================ runs
OUT = []
def P(s=""):
    OUT.append(s); print(s, flush=True)

def spawn(n):
    return MASTER.spawn(n)

def main():
    t0 = time.time()
    P("REFUTER A (COMPUTATION) -- T-50 DESIGN ONE -- INDEPENDENT VERIFICATION RUN")
    P("date 2026-08-21 | refuter master seed SeedSequence(4242) | numpy %s" % np.__version__)
    P("all pipeline code refuter-owned; model laws: programmed=-100 e, residual U{-5..5} e,")
    P("grains +-1, DC=+1 (model/geometry.py, read-only anchor).")
    P("")

    # ---------- S0 anchor ----------
    pat = GE.occupancy_patterns()
    err = pat["unwritten_e"].astype(float)
    rs = [abs(np.where(w == 1, -NE, err).sum()) / np.abs(np.where(w == 1, -NE, err)).sum()
          for w in pat["written"]]
    P("S0 ANCHOR: sealed-page rho min %.15f max %.15f (brief: 0.968134 / 0.979799)"
      % (min(rs), max(rs)))
    P("")

    # ---------- S1 rebuild baseline ----------
    P("S1 REBUILD -- 50 independent occupancy reads (f=0.5) and 50 orientation reads,")
    P("   refuter's own pipeline, K=16, 9-point grid (the lane's settings).")
    occ = dict(bWU=[], bUU=[], bWW=[], xi=[], ds=[], B1=[], states=[])
    for ss in spawn(50):
        rng = np.random.default_rng(ss)
        v, p, c = build_occ(rng, {"W": gen_uniform(0.5)}, ["W", "U"])
        m = measure_occ(v, p, c, rng)
        occ["states"].append(m["state"])
        if "bWU" in m:
            for k in ("bWU", "bUU", "bWW", "xi"):
                occ[k].append(m[k])
            occ["ds"].append(m["dslope"]); occ["B1"].append(m["B1"])
    for k, lab in (("bWU", "beta_WU"), ("bUU", "beta_UU"), ("bWW", "beta_WW"),
                   ("xi", "xi"), ("ds", "density slope")):
        m_, lo, hi = med_iqr(occ[k])
        P("   occ %-14s median %+8.4f  IQR [%+8.4f, %+8.4f]   (lane: %s)" %
          (lab, m_, lo, hi, dict(bWU="+0.9998", bUU="+0.5003", bWW="+0.5152",
                                 xi="+0.4993", ds="+0.0001")[k]))
    P("   occ B1 true %d/50 (lane 49/50); read states: %s" %
      (sum(occ["B1"]), {s: occ["states"].count(s) for s in set(occ["states"])}))
    ori = dict(bDU=[], bUU=[], bCU=[], xi=[], B3=[], states=[])
    for ss in spawn(50):
        rng = np.random.default_rng(ss)
        v, c = build_ori(rng)
        m = measure_ori(v, c, rng)
        ori["states"].append(m["state"])
        if "bDU" in m:
            ori["bDU"].append(m["bDU"]); ori["bUU"].append(m["bUU"])
            ori["xi"].append(m["xi"]); ori["B3"].append(m["B3"])
            if m["bCU"] is not None: ori["bCU"].append(m["bCU"])
    for k, lab in (("bDU", "beta_DATA-U"), ("bUU", "beta_UU"), ("bCU", "beta_DC-U"),
                   ("xi", "xi")):
        m_, lo, hi = med_iqr(ori[k])
        P("   ori %-14s median %+8.4f  IQR [%+8.4f, %+8.4f]   (lane: %s)" %
          (lab, m_, lo, hi, dict(bDU="+0.4968", bUU="+0.5008", bCU="+0.9993",
                                 xi="+0.0075")[k]))
    P("   ori B3 true %d/50 (lane 49/50)" % sum(ori["B3"]))
    P("")

    # ---------- S2 constraint-2 table ----------
    P("S2 CONSTRAINT-2 -- moment choices, one read each (lane: 0.9934/1.9867/1.0412 and")
    P("   0.4696/0.9314/0.9306):")
    rng = np.random.default_rng(spawn(1)[0])
    v, p, c = build_occ(rng, {"W": gen_uniform(0.5)}, ["W", "U"])
    tw = block_tables(csum(v), adj_pairs(c, "W", "U"), rng, 16, GRID9)
    bA = fit_loglog(GRID9, tw["A"])[0]; bM = fit_loglog(GRID9, tw["M2"])[0]
    bV = fit_loglog(GRID9, tw["VAR"])[0]
    P("   occ  W-U   |D| %+7.4f | uncentred D^2 %+7.4f | centred var %+7.4f" % (bA, bM, bV))
    P("   identity check: 2 x |D| slope = %+7.4f vs D^2 slope %+7.4f (exact: median of a" %
      (2 * bA, bM))
    P("   monotone transform)")
    rng = np.random.default_rng(spawn(1)[0])
    v, c = build_ori(rng)
    td = block_tables(csum(v), adj_pairs(c, "DATA", "U"), rng, 16, GRID9)
    bA2 = fit_loglog(GRID9, td["A"])[0]; bM2 = fit_loglog(GRID9, td["M2"])[0]
    bV2 = fit_loglog(GRID9, td["VAR"])[0]
    P("   ori DATA-U |D| %+7.4f | uncentred D^2 %+7.4f | centred var %+7.4f" % (bA2, bM2, bV2))
    P("   centred variance ~1 for BOTH encodings: CONFIRMED non-discriminating")
    P("")

    # ---------- S3 offsets and drift ----------
    P("S3 COMMON-MODE OFFSET / DRIFT SWEEP (10 reps each, medians; constraint 1):")
    for name, kw in (("offset 0", {}), ("offset +0.5 e", dict(offset=0.5)),
                     ("offset -0.5 e", dict(offset=-0.5)), ("offset +5 e", dict(offset=5.0)),
                     ("offset +50 e", dict(offset=50.0)), ("drift 50 e", dict(drift=50.0))):
        bc, br = [], []
        for ss in spawn(10):
            rng = np.random.default_rng(ss)
            v, p, c = build_occ(rng, {"W": gen_uniform(0.5)}, ["W", "U"], **kw)
            cs = csum(v)
            wu = block_tables(cs, adj_pairs(c, "W", "U"), rng, 16, GRID9)
            uu = block_tables(cs, same_pairs(c, "U"), rng, 16, GRID9)
            i0 = guard_i0(wu["A"], uu["A"], GRID9) or 0
            bc.append(fit_loglog(GRID9[i0:], wu["A"][i0:])[0])
            br.append(fit_loglog(GRID9, wu["RAW"])[0])
        P("   occ %-14s contrast %+8.4f  raw %+8.4f" %
          (name, float(np.median(bc)), float(np.median(br))))
    for name, off in (("offset 0", 0.0), ("offset +0.05", 0.05), ("offset +0.5", 0.5)):
        bc, br = [], []
        for ss in spawn(10):
            rng = np.random.default_rng(ss)
            v, c = build_ori(rng, offset=off)
            cs = csum(v)
            du = block_tables(cs, adj_pairs(c, "DATA", "U"), rng, 16, GRID9)
            bc.append(fit_loglog(GRID9, du["A"])[0])
            br.append(fit_loglog(GRID9, du["RAW"])[0])
        P("   ori %-14s contrast %+8.4f  raw %+8.4f" %
          (name, float(np.median(bc)), float(np.median(br))))
    P("")

    # ---------- S4 mutations ----------
    P("S4 MUTATIONS (50 seeds each):")
    flips = fires = 0; bws = []
    for ss in spawn(50):
        rng = np.random.default_rng(ss)
        v, p, c = build_occ(rng, {"W": gen_uniform(0.5)}, ["W", "U"], two_signed=True)
        m = measure_occ(v, p, c, rng)
        if m["state"] == "READS_UNWRITTEN" or not m.get("B1", False): flips += 1
        if m.get("fire_a", False): fires += 1
        if "bWU" in m: bws.append(m["bWU"])
    P("   M1 two-signed write: B1 flips %d/50, registered clause (a) fires %d/50," % (flips, fires))
    P("      guard-passing beta_WU median %+0.4f (lane: flips 50/50, fires 50/50, %+0.4f)" %
      (float(np.median(bws)) if bws else float("nan"), 0.4859))
    for mu in (2.0, 20.0):
        ok = 0; bwl = []
        for ss in spawn(50):
            rng = np.random.default_rng(ss)
            v, p, c = build_occ(rng, {"W": gen_uniform(0.5)}, ["W", "U"], mu=mu)
            m = measure_occ(v, p, c, rng)
            if m.get("B1", False): ok += 1
            if "bWU" in m: bwl.append(m["bWU"])
        P("   M2 residual mu=%+5.1f e: B1 %d/50, beta_WU median %+0.4f" %
          (mu, ok, float(np.median(bwl))))
    for mode, aligned in (("aligned", True), ("random-offset", False)):
        paths = {"READS_UNWRITTEN": 0, "VOID_DENSITY": 0, "INCONCLUSIVE_CONTROL": 0,
                 "INCONCLUSIVE_ZERO_A": 0, "OK": 0}
        fired = 0
        for ss in spawn(50):
            rng = np.random.default_rng(ss)
            v, p, c = build_occ(rng, {"W": gen_fixed_record(0.5)}, ["W", "U"])
            m = measure_occ(v, p, c, rng, aligned=aligned)
            paths[m["state"]] += 1
            if m.get("fire_a", False): fired += 1
        P("   M3 fixed-record %-14s state counts %s; clause (a) fires %d/50" %
          (mode, {k: v_ for k, v_ in paths.items() if v_}, fired))
    P("")

    # ---------- S5 honest false-fire at lane settings ----------
    P("S5 HONEST FALSE-FIRE, lane settings (K=16, 9 pts) -- registered clause rules:")
    fa = 0; n_ok = 0
    for ss in spawn(200):
        rng = np.random.default_rng(ss)
        v, p, c = build_occ(rng, {"W": gen_uniform(0.5)}, ["W", "U"])
        m = measure_occ(v, p, c, rng)
        if m["state"] == "OK": n_ok += 1
        if m.get("fire_a", False): fa += 1
    P("   clause (a): %d/200 fires (%d OK reads)" % (fa, n_ok))
    fb = 0; n_ok = 0
    for ss in spawn(200):
        rng = np.random.default_rng(ss)
        v, p, c = build_occ(rng, {"W1": gen_uniform(0.35), "W2": gen_uniform(0.65)},
                            ["W1", "U", "W2", "U"])
        m = measure_occ_pair(v, p, c, rng)
        if m["state"] == "OK": n_ok += 1
        if m.get("fire_b", False): fb += 1
    P("   clause (b) f=0.35 vs 0.65: %d/200 fires (%d OK reads)" % (fb, n_ok))
    fc = 0
    for ss in spawn(200):
        rng = np.random.default_rng(ss)
        v, c = build_ori(rng)
        if measure_ori(v, c, rng).get("fire_c", False): fc += 1
    P("   clause (c): %d/200 fires" % fc)
    P("")

    # ---------- S6 seam hammer ----------
    P("S6 SEAM HAMMER -- low-density occupancy, 150 seeds per f, kappa=8 (the design's")
    P("   guard); counting clause-(a) fires and guard-passing reads below the 0.9 band:")
    for f in (0.012, 0.015, 0.02, 0.03):
        unw = fires = below = okc = 0; worst = float("inf")
        for ss in spawn(150):
            rng = np.random.default_rng(ss)
            v, p, c = build_occ(rng, {"W": gen_uniform(f)}, ["W", "U"])
            m = measure_occ(v, p, c, rng)
            if m["state"] == "READS_UNWRITTEN": unw += 1; continue
            if m["state"] != "OK": continue
            okc += 1
            worst = min(worst, m["bWU"])
            if m["bWU"] < 0.9: below += 1
            if m["fire_a"]: fires += 1
        P("   f=%.3f  unwritten %3d/150  OK %3d  fires %d  below-band %d  worst beta %+0.4f" %
          (f, unw, okc, fires, below, worst if okc else float("nan")))
    P("")

    # ---------- S7 registered minimum ----------
    P("S7 REGISTERED-MINIMUM READER (K=8, 8-point grid 16..2048) -- the floor the text")
    P("   permits; honest surfaces, 200 seeds each:")
    fa = inc = 0; states = {}
    for ss in spawn(200):
        rng = np.random.default_rng(ss)
        v, p, c = build_occ(rng, {"W": gen_uniform(0.5)}, ["W", "U"])
        m = measure_occ(v, p, c, rng, K=8, grid=GRID8MIN)
        states[m["state"]] = states.get(m["state"], 0) + 1
        if m.get("fire_a", False): fa += 1
    P("   occ clause (a): %d/200 fires; states %s" % (fa, states))
    fc = 0; states = {}; zero_a = 0
    for ss in spawn(200):
        rng = np.random.default_rng(ss)
        v, c = build_ori(rng)
        m = measure_ori(v, c, rng, K=8, grid=GRID8MIN)
        states[m["state"]] = states.get(m["state"], 0) + 1
        if m["state"] == "INCONCLUSIVE_ZERO_A": zero_a += 1
        if m.get("fire_c", False): fc += 1
    P("   ori clause (c): %d/200 fires; states %s; zero-median-A incidents %d" %
      (fc, states, zero_a))
    P("")

    # ---------- S8 ATTACK A: density-shape confound ----------
    P("=" * 78)
    P("S8 ATTACK A -- WITHIN-SECTOR DENSITY SKEW (constraint-3 in a new costume).")
    P("   Legal one-carrier patterns whose MEAN density is flat in N (the registered")
    P("   density check passes) but whose MEDIAN block density rises with N, because")
    P("   the registered statistic median|D_k| reads the TYPICAL block, not the mean.")
    P("   Each pattern is ONE declared literal mask, written identically to every W")
    P("   sector; guard, density check, control band all applied as registered.")
    P("")

    def attack_pattern(tag, genfn, K, nsect=256, nseeds=40):
        stats = dict(bWU=[], se=[], states={}, above=0, ok=0, dsl=[])
        for ss in spawn(nseeds):
            rng = np.random.default_rng(ss)
            v, p, c = build_occ(rng, {"W": genfn}, ["W", "U"], nsect=nsect)
            m = measure_occ(v, p, c, rng, K=K)
            stats["states"][m["state"]] = stats["states"].get(m["state"], 0) + 1
            if m["state"] != "OK": continue
            stats["ok"] += 1
            stats["bWU"].append(m["bWU"]); stats["se"].append(m["seWU"])
            stats["dsl"].append(m["dslope"])
            if m["band_above"]: stats["above"] += 1
        if stats["bWU"]:
            bm, lo, hi = med_iqr(stats["bWU"])
            P("   %-34s K=%-3d  OK %2d/%d  beta_WU med %+7.4f IQR[%+7.4f,%+7.4f]" %
              (tag, K, stats["ok"], nseeds, bm, lo, hi))
            P("   %-34s        SE med %.4f  dens-slope med %+0.4f  band-above %d/%d" %
              ("", float(np.median(stats["se"])), float(np.median(stats["dsl"])),
               stats["above"], stats["ok"]))
        else:
            P("   %-34s K=%-3d  no OK reads; states %s" % (tag, K, stats["states"]))
        if stats["states"].get("OK", 0) < nseeds:
            P("   %-34s        non-OK states: %s" %
              ("", {k: v_ for k, v_ in stats["states"].items() if k != "OK"}))
        return stats

    # fixed literal patterns, one per family (drawn from a dedicated pattern seed)
    prng = np.random.default_rng(np.random.SeedSequence(777))
    cas15, f15 = gen_cascade_fixed(prng, levels=8, hi=1.5, f0=0.25)
    cas16, f16 = gen_cascade_fixed(prng, levels=8, hi=1.6, f0=0.25)
    two1, ft1 = gen_twolevel_fixed(prng, 0.03, 0.95, 0.25)
    P("   declared literal patterns: cascade(1.5/0.5,8lvl) f=%.4f | cascade(1.6/0.4,8lvl)"
      % f15)
    P("   f=%.4f | two-level(0.03/0.95, dense 1/4) f=%.4f" % (f16, ft1))
    P("")
    st_c15_16 = attack_pattern("cascade 1.5/0.5", cas15, K=16)
    st_c15_64 = attack_pattern("cascade 1.5/0.5", cas15, K=64)
    st_c16_16 = attack_pattern("cascade 1.6/0.4", cas16, K=16)
    st_c16_64 = attack_pattern("cascade 1.6/0.4", cas16, K=64)
    st_tl_16  = attack_pattern("two-level 0.03/0.95 1/4", two1, K=16)
    st_tl_64  = attack_pattern("two-level 0.03/0.95 1/4", two1, K=64)
    P("")
    P("   one cascade(1.6/0.4) read, A_WU(N) table (the shape of the confound):")
    rng = np.random.default_rng(spawn(1)[0])
    v, p, c = build_occ(rng, {"W": cas16}, ["W", "U"])
    cs = csum(v)
    wu = block_tables(cs, adj_pairs(c, "W", "U"), rng, 64, GRID9, prog_cs=csum(p))
    uu = block_tables(cs, same_pairs(c, "U"), rng, 64, GRID9)
    for i, N in enumerate(GRID9):
        P("      N=%5d  A_WU %12.1f e   A_UU %9.1f e   mean f-hat %.4f" %
          (N, wu["A"][i], uu["A"][i], wu["FH"][i]))
    P("")
    P("   REGISTERED CLAUSE (b): the skew pattern paired with uniform f=0.5 on the SAME")
    P("   part, same read (cycle [W1,U,W2,U]); both classes must be guard-passing,")
    P("   density-trend-free, control in band; fire iff |delta beta| - 2 SE_delta > 0.2:")
    for tag, genfn, K, nsect, nseeds in (
            ("cascade 1.5/0.5 vs uniform", cas15, 16, 256, 40),
            ("cascade 1.5/0.5 vs uniform", cas15, 64, 512, 40),
            ("cascade 1.6/0.4 vs uniform", cas16, 16, 256, 40),
            ("cascade 1.6/0.4 vs uniform", cas16, 64, 512, 40),
            ("cascade 1.6/0.4 vs uniform", cas16, 128, 512, 25),
            ("two-level vs uniform", two1, 64, 512, 40)):
        fb = okc = 0; deltas = []; margins = []
        states = {}
        for ss in spawn(nseeds):
            rng = np.random.default_rng(ss)
            v, p, c = build_occ(rng, {"W1": genfn, "W2": gen_uniform(0.5)},
                                ["W1", "U", "W2", "U"], nsect=nsect)
            m = measure_occ_pair(v, p, c, rng, K=K)
            states[m["state"]] = states.get(m["state"], 0) + 1
            if m["state"] != "OK": continue
            okc += 1; deltas.append(m["delta"]); margins.append(m["delta"] - 2 * m["sed"])
            if m["fire_b"]: fb += 1
        P("   %-30s K=%-3d nsect=%d: FIRES %d/%d OK reads (of %d seeds)" %
          (tag, K, nsect, fb, okc, nseeds))
        if deltas:
            P("      |delta beta| med %+0.4f  fire margin (|d|-2SE-0.2) med %+0.4f  states %s" %
              (float(np.median(deltas)), float(np.median(margins)) - 0.2,
               {k: v_ for k, v_ in states.items() if k != "OK"} or "{all OK}"))
    P("")

    # ---------- S9 ATTACK B: sector-differential offsets ----------
    P("=" * 78)
    P("S9 ATTACK B -- SECTOR-DIFFERENTIAL BASELINE OFFSETS (the half of constraint 1 the")
    P("   design does not cancel: D cancels only offsets COMMON to the two sectors of a")
    P("   pair; per-sector and per-class response offsets survive the subtraction).")
    P("")
    P("   B-1 iid per-sector offsets (instrument response varying sector to sector):")
    for surf, sds in (("occ", (0.5, 2.0, 5.0)), ("ori", (0.02, 0.05, 0.2))):
        for sd in sds:
            fires = 0; states = {}
            for ss in spawn(40):
                rng = np.random.default_rng(ss)
                if surf == "occ":
                    v, p, c = build_occ(rng, {"W": gen_uniform(0.5)}, ["W", "U"],
                                        sector_off_sd=sd)
                    m = measure_occ(v, p, c, rng)
                    fires += bool(m.get("fire_a", False))
                else:
                    v, c = build_ori(rng, sector_off_sd=sd)
                    m = measure_ori(v, c, rng)
                    fires += bool(m.get("fire_c", False))
                states[m["state"]] = states.get(m["state"], 0) + 1
            P("   %s sector-offset sd=%-5g: fires %d/40, states %s" % (surf, sd, fires, states))
    P("")
    P("   B-2 CLASS-CORRELATED read offset on the magnetic surface: the mapper reads")
    P("   WRITTEN sectors (DATA and DC alike) with a small per-cell additive offset that")
    P("   is not magnetization (topography/charging crosstalk).  Medium physics correct")
    P("   (data screens); registered clause (c) fire rate:")
    for co in (0.02, 0.05, 0.1, 0.2):
        fires = 0; xis = []; states = {}
        for ss in spawn(40):
            rng = np.random.default_rng(ss)
            v, c = build_ori(rng, class_offset=co)
            m = measure_ori(v, c, rng)
            states[m["state"]] = states.get(m["state"], 0) + 1
            if "xi" in m: xis.append(m["xi"])
            if m.get("fire_c", False): fires += 1
        P("   class offset %5.2f /grain: clause (c) FIRES %d/40, xi med %+0.4f, states %s" %
          (co, fires, float(np.median(xis)) if xis else float("nan"), states))
    P("")

    # ---------- S10 supplementary rebuilds ----------
    P("S10 SUPPLEMENTARY REBUILDS:")
    bds = []
    for ss in spawn(20):
        rng = np.random.default_rng(ss)
        v, c = build_ori(rng, data="dcfree")
        m = measure_ori(v, c, rng)
        if "bDU" in m: bds.append(m["bDU"])
    P("   R6 dc-free coded: beta_DATA-U median %+0.4f (lane +0.5266)" % float(np.median(bds)))
    for tag, gen in (("f=0.10", gen_uniform(0.10)), ("f=0.90", gen_uniform(0.90)),
                     ("alternating", gen_alternating())):
        bl = []
        for ss in spawn(20):
            rng = np.random.default_rng(ss)
            v, p, c = build_occ(rng, {"W": gen}, ["W", "U"])
            m = measure_occ(v, p, c, rng)
            if "bWU" in m: bl.append(m["bWU"])
        P("   R7 %-12s beta_WU median %+0.4f" % (tag, float(np.median(bl))))
    bd, bc2 = [], []
    for ss in spawn(20):
        rng = np.random.default_rng(ss)
        v, c = build_ori(rng, tilt_deg=30.0)
        m = measure_ori(v, c, rng)
        if "bDU" in m:
            bd.append(m["bDU"])
            if m["bCU"] is not None: bc2.append(m["bCU"])
    P("   R8 tilt 30deg: beta_DATA-U %+0.4f, beta_DC-U %+0.4f (lane +0.4999/+1.0015)" %
      (float(np.median(bd)), float(np.median(bc2))))
    bias_fire = 0
    for ss in spawn(50):
        rng = np.random.default_rng(ss)
        v, c = build_ori(rng, erase_bias=0.25)
        if measure_ori(v, c, rng).get("fire_c", False): bias_fire += 1
    P("   M2-ori biased erase +0.25: clause (c) fires %d/50 if scope ignored (lane 34/50)"
      % bias_fire)
    P("")
    P("elapsed %.1f s" % (time.time() - t0))
    P("END OF REFUTER RUN")

if __name__ == "__main__":
    main()
    with open(os.path.join(HERE, "verify_a_run.txt"), "w") as fh:
        fh.write("\n".join(OUT) + "\n")
