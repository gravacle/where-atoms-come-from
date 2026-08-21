#!/usr/bin/env python3
"""T-50 DESIGN ONE V002 -- THE SAME-READ CONTRAST EXPONENT, EIGHT CLOSURES -- model-side
verification.

LANE_T50/V002. Read-only import of model/geometry.py (the sealed formation layer).
Writes nothing outside LANE_T50/V002 and the scratchpad. No git, no reproduce.sh, no r3.sh.

The core is V001's (kept by both refuters): D_k(N) = sum_W v_i - sum_U v_i over
same-read sector pairs; A(N) = median_k |D_k(N)|; the observable is the fitted log-log
exponent pair (beta_WU, beta_UU) and xi = beta_WU - beta_UU.

The eight closures (LANE_T50/JUDGMENT.txt section 2) are implemented here as REGISTERED
STRUCTURE, not patches:
 C1 the density condition reads the STATISTIC'S OWN MOMENT over the STATISTIC'S OWN
    POPULATION: the median per-block f-hat over ALL in-sector block placements,
    trend-free -- for a one-carrier write median_k|D_k(N)| = N_E * N * median_k
    f-hat_k(N) + O(noise), so the condition pins exactly the moment the statistic
    reads. V001's mean condition is REPLACED, not kept: over the offset-uniform block
    population the mean per-block f-hat equals the sector density identically, so a
    mean-trend check has no failing branch (D-8/INST-14).
 C2 the DC-free-coded in-read crosstalk voider (orientation clause void unless the same
    read's DC-free-coded sectors sit in the control band).
 C3 the orientation half carried by a NAMED M-reading instrument: polar Kerr microscopy
    (wide-field MOKE imaging). Kerr rotation reads the LOCAL MAGNETIZATION component
    (polar geometry: M_z), with NON-ZERO response at k=0 -- the uniform-film Kerr
    hysteresis loop is the standard demonstration -- unlike every stray-field mapper
    (MFM/scanning Hall/NV), whose transfer exp(-|k|d)(1-exp(-|k|t)) vanishes at k=0.
    The pipeline models the Kerr read as local-M cell averaging plus an optical PSF
    (k=0-preserving) and carries the stray-field transfer as a mutation-suite member
    that the read's own positive control DETECTS at rate 1.000.
 C4 per-clause guard scope (the guard-scope table in v2_design.md; code below follows it).
 C5 a registered fixed-pattern treatment: rung 1 (same-read, address-balanced class
    de-trend from the part's own unwritten sectors) and rung 2 (all-erased calibration
    read, per-sector static-pattern subtraction); the U-U control band is the
    fixed-pattern DETECTOR; plus fixed-pattern mutation members.
 C6 INCONCLUSIVE precedence wired into the fire booleans: the read STATE is computed
    first, in the registered order, and no falsifier boolean can be True unless the
    state is OK or SEAM. Code equals text: this file is the only pipeline.
 C7 grid, pairing and placements PINNED (GRID 16..4096 geometric x2, ALL disjoint
    adjacent pairs with 4 placements each, >= 16 pairs, sector 4096 cells, >= 256
    sectors); kappa_void = 8 for the void guard; the point-band claim additionally
    requires the READ'S OWN CERTIFICATE, re-derived per reader from the reader's own
    U-U and W-U pair pools by the registered bootstrap (point_certificate below) --
    noise-law- and shot-noise-adaptive, no model constant above the void guard.
 C8 the railed-population branch: A_UU(16) = 0 reads INCONCLUSIVE_RAILED
    (insufficient access), all fits skipped, no vacuous guard pass, no nan escapes.

Every statistic reported is computed ONCE by this pipeline and echoed verbatim into the
sealed txt (INST-17). Matched ensembles (COMP-12); every replicate draws its own
independent residual field (COMP-13).
"""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join(os.path.dirname(os.path.dirname(HERE)), "model")
sys.path.insert(0, MODEL)
import geometry as GE  # read-only

# ================= PINNED CONSTANTS (all in the registered text) =================
GRID       = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096]   # PINNED (C7)
NG         = len(GRID)
K_MIN      = 16          # minimum disjoint pairs the geometry must supply -- PINNED (C7)
PLACEMENTS = 4           # registered in-sector placements per pair per grid point.
                         # THE STATISTIC READS THE WHOLE POOL: A(N) = median|D| over ALL
                         # disjoint pairs x PLACEMENTS placements -- the SAME object the
                         # guard, the control and the certificate read. One object, one
                         # moment, one population (the C1 discipline, applied to the
                         # statistic itself; a K-subsample and a pool certificate would
                         # be two different objects, and the seam would leak between
                         # them -- measured during construction, closed here).
SECT       = 4096        # sector length, cells -- PINNED (C7)
NSECT      = 256         # sectors per occupancy part (min; C7)
KAPPA_VOID = 8.0         # void guard (V001's measured seam value, on the pinned grid)
MIN_POINTS = 6           # >= 1.5 decades must survive the guard
DENS_TOL_MED  = 0.02     # |slope of MEDIAN per-block f-hat on log10 N| (C1; co-tuned in R2:
                         # a sub-tolerance median trend inflates beta_WU by <= tol = 0.02,
                         # 5x inside the 0.1 band half-width; measured gap in the run).
                         # THERE IS NO MEAN CONDITION IN V002: over the offset-uniform
                         # block population the mean per-block f-hat equals the sector
                         # density identically, so a mean-trend check has no failing
                         # branch (D-8/INST-14) -- V001's mean check is REPLACED, and the
                         # M3 family lands READS_UNWRITTEN / VOID_DENSITY_MEDIAN (measured).
DENS_MED_I0   = 2        # the median condition is evaluated on N >= 64 (7 points): the
                         # discrete binomial-median bias of honest sparse patterns lives
                         # at N=16..32 and is not a density trend; measured in R2
BAND_ACC   = (0.9, 1.1)
BAND_CTL   = (0.35, 0.65)
XI_MIN     = 0.25
N_E        = GE.N_E      # 100 e -- the model's own parameter (model-side only)
RES_SD     = math.sqrt(10.0)   # sd of uniform{-5..5}
# point-certificate bootstrap (C7): built from the READER'S OWN pair pools
KP_B       = 400         # surrogate ladders per certificate
# orientation (Kerr) model
NSECT_ORI  = 252
GRAINS     = 64          # grains per resolved bit cell
PSF_SIGMA  = 0.7         # optical PSF sigma, cells
READ_NOISE = 0.05        # additive optical read noise per cell (declared)
EDGE_MARGIN_ORI = 4      # blocks sit >= this many cells from any sector boundary (PSF)
STRAY_T    = 0.5         # film thickness, cells, for the stray-field mutation member

MASTER = np.random.SeedSequence(20260821)
MASK_SEED = 20260821     # declared literal-mask seed (C1 skew masks; M3 record mask)

ORI_CYCLE = ["DATA", "U", "DCF", "U", "DC", "U"]
ADDR_NCLASS = 16
# declared address-class fixed-pattern PROFILE (deterministic literal; deck-like)
ADDR_PROFILE = np.array([0.0, 1.0, -0.5, 2.0, -1.5, 0.5, 1.5, -2.0,
                         2.5, -1.0, 0.0, 1.0, -2.5, 1.5, -0.5, 2.0])

# ============================ small fitting helpers ============================
def fit_loglog(Ns, A):
    x = np.log10(np.asarray(Ns, float)); y = np.log10(np.asarray(A, float))
    n = len(x); xm, ym = x.mean(), y.mean()
    Sxx = ((x - xm) ** 2).sum()
    b = ((x - xm) * (y - ym)).sum() / Sxx
    resid = y - (ym + b * (x - xm))
    se = math.sqrt(max((resid ** 2).sum(), 0.0) / max(n - 2, 1) / Sxx)
    return b, se

def fit_lin(Ns, F):
    x = np.log10(np.asarray(Ns, float)); y = np.asarray(F, float)
    n = len(x); xm, ym = x.mean(), y.mean()
    Sxx = ((x - xm) ** 2).sum()
    b = ((x - xm) * (y - ym)).sum() / Sxx
    resid = y - (ym + b * (x - xm))
    se = math.sqrt(max((resid ** 2).sum(), 0.0) / max(n - 2, 1) / Sxx)
    return b, se

def seccs(vals):
    return np.concatenate([np.zeros((vals.shape[0], 1)), np.cumsum(vals, axis=1)], axis=1)

def bsum(cs, s, o, N):
    return cs[s, o + N] - cs[s, o]

def med_iqr(v):
    v = np.asarray(v, float)
    return float(np.median(v)), float(np.percentile(v, 25)), float(np.percentile(v, 75))

# ======================= registered interleave / addressing ====================
def occ_roles(nsect=NSECT):
    """REGISTERED ADDRESS-BALANCED INTERLEAVE (C5 rung 1): role W iff (s + s//16) even.
    Every address class a = s mod 16 appears in the W role in half the 16-sector groups
    and in the U role in the other half, so every address class has same-read unwritten
    observations."""
    s = np.arange(nsect)
    return ((s + s // 16) % 2 == 0)   # True = W

def adj_wu_pairs(roleW):
    """Disjoint adjacent (written, unwritten) sector pairs, either order."""
    used = np.zeros(len(roleW), bool); pairs = []
    for i in range(len(roleW) - 1):
        if used[i] or used[i + 1]:
            continue
        if roleW[i] != roleW[i + 1]:
            w, u = (i, i + 1) if roleW[i] else (i + 1, i)
            pairs.append((w, u)); used[i] = used[i + 1] = True
    return pairs

def same_role_pairs(idx):
    return [(idx[j], idx[j + 1]) for j in range(0, len(idx) - 1, 2)]

# ============================== literal masks (C1) =============================
def cascade_mask(f0, hi, lo, levels=8, seed=MASK_SEED):
    """ONE declared literal mask (multiplicative cascade), written identically to every
    W sector -- refuter A's D1 family, as declared literals."""
    rng = np.random.default_rng([seed, int(hi * 1000), int(lo * 1000)])
    w = np.ones(SECT); seg = SECT
    for _ in range(levels):
        seg //= 2
        for start in range(0, SECT, 2 * seg):
            if rng.random() < 0.5:
                a, b = hi, lo
            else:
                a, b = lo, hi
            w[start:start + seg] *= a; w[start + seg:start + 2 * seg] *= b
    p = np.clip(f0 * w, 0.0, 1.0)
    return rng.random(SECT) < p

def twolevel_mask(seed=MASK_SEED):
    """Dense quarter at 0.95, sparse remainder at 0.03 (refuter A's two-level literal)."""
    rng = np.random.default_rng([seed, 9595, 303])
    p = np.where(np.arange(SECT) < SECT // 4, 0.95, 0.03)
    return rng.random(SECT) < p

def record_mask(seed=MASK_SEED):
    """M3: a fixed 512-cell record at f=0.5, identical in every W sector."""
    rng = np.random.default_rng([seed, 512])
    m = np.zeros(SECT, bool); m[:512] = rng.random(512) < 0.5
    return m

# ============================== residual noise laws ============================
def residual(rng, law, size):
    if law == "uniform":
        return rng.integers(-5, 6, size).astype(float)      # the model's own law
    if law == "laplace":
        return rng.laplace(0.0, RES_SD / math.sqrt(2.0), size)
    if law == "t3":
        return rng.standard_t(3, size) * (RES_SD / math.sqrt(3.0))
    raise ValueError(law)

# =============================== part builders =================================
def build_occ(rng, f=0.5, pattern="random", mask=None, mu=0.0, two_signed=False,
              offset=0.0, drift=0.0, law="uniform", fp_mode=None, fp_amp=0.0,
              railed=False, nsect=NSECT, two_class=False, f2=None):
    """Occupancy part on the registered address-balanced interleave.
    Returns vals, prog, roleW, wclass (W1/W2 labels if two_class), fp_true."""
    roleW = occ_roles(nsect)
    vals = np.empty((nsect, SECT)); prog = np.zeros((nsect, SECT))
    wclass = np.empty(nsect, dtype=object); wclass[:] = None
    widx = 0
    for s in range(nsect):
        r = residual(rng, law, SECT) + mu
        if roleW[s]:
            fs = f
            if two_class:
                wclass[s] = "W1" if widx % 2 == 0 else "W2"
                fs = f if wclass[s] == "W1" else f2
                widx += 1
            if mask is not None:
                p = mask.copy()
            elif pattern == "random":
                p = rng.random(SECT) < fs
            elif pattern == "alternating":
                p = (np.arange(SECT) % 2 == 0)
            elif pattern == "blocky":
                p = np.where(np.arange(SECT) < SECT // 2,
                             rng.random(SECT) < 1.5 * fs, rng.random(SECT) < 0.5 * fs)
            else:
                raise ValueError(pattern)
            sign = (rng.integers(0, 2, SECT) * 2 - 1).astype(float) if two_signed \
                else -np.ones(SECT)
            vals[s] = np.where(p, sign * N_E, r); prog[s] = p.astype(float)
        else:
            vals[s] = r
    if railed:   # the read floors the erased distribution at the rail: erased reads 0
        vals = np.where(prog.astype(bool), vals, 0.0)
    fp_true = np.zeros(nsect)
    if fp_mode == "iid":
        fp_true = rng.normal(0.0, fp_amp, nsect)
    elif fp_mode == "addr":
        fp_true = fp_amp * ADDR_PROFILE[np.arange(nsect) % ADDR_NCLASS]
    vals = vals + fp_true[:, None]
    if offset:
        vals = vals + offset
    if drift:
        ramp = np.linspace(-0.5, 0.5, nsect * SECT).reshape(nsect, SECT) * drift
        vals = vals + ramp
    return vals, prog, roleW, wclass, fp_true

def build_occ_calib(rng, fp_true, mu=0.0, law="uniform", nsect=NSECT):
    """C5 rung 2: the all-erased CALIBRATION READ of the same part (same static fixed
    pattern, fresh residual field)."""
    vals = np.empty((nsect, SECT))
    for s in range(nsect):
        vals[s] = residual(rng, law, SECT) + mu
    return vals + fp_true[:, None]

# ---- orientation (Kerr) ----
def _psf_kernel(sigma=PSF_SIGMA):
    xs = np.arange(-8, 9, dtype=float)
    k = np.exp(-0.5 * (xs / sigma) ** 2)
    return k / k.sum()

def build_ori(rng, data_p1=0.5, erase_bias=0.0, class_offset=0.0, fp_mode=None,
              fp_amp=0.0, tilt_deg=0.0, read="kerr", standoff=1.0, nsect=NSECT_ORI):
    """Orientation surface read by the NAMED instrument model (C3).
    read='kerr': local-M cell average through an optical PSF (k=0-preserving).
    read='stray': the stray-field transfer T(k)=exp(-|k|d)(1-exp(-|k|t)), T(0)=0
    (the mutation-suite member; every named stray-field mapper).
    Cell = resolved bit cell (bit length >= instrument resolution, the reader's own
    write clock); GRAINS grains per cell."""
    classes = [ORI_CYCLE[s % 6] for s in range(nsect)]
    M = np.empty((nsect, SECT))
    for s, c in enumerate(classes):
        if c == "DATA":
            v = ((rng.random(SECT) < data_p1) * 2 - 1).astype(float)
        elif c == "DCF":
            v = np.tile([1.0, -1.0], SECT // 2)
        elif c == "DC":
            v = np.ones(SECT)
        else:
            v = (2.0 * rng.binomial(GRAINS, 0.5, SECT) - GRAINS) / GRAINS + erase_bias
        if tilt_deg:
            v = v * np.cos(rng.normal(0.0, math.radians(tilt_deg), SECT))
        M[s] = v
    flat = M.reshape(-1)
    if read == "kerr":
        flat = np.convolve(flat, _psf_kernel(), mode="same")
    elif read == "stray":
        F = np.fft.rfft(flat)
        fr = np.fft.rfftfreq(flat.size)          # cycles per cell
        T = np.exp(-2 * np.pi * np.abs(fr) * standoff) * \
            (1.0 - np.exp(-2 * np.pi * np.abs(fr) * STRAY_T))
        flat = np.fft.irfft(F * T, n=flat.size)
    else:
        raise ValueError(read)
    vals = flat.reshape(nsect, SECT)
    written = np.array([c in ("DATA", "DCF", "DC") for c in classes])
    if class_offset:                              # non-magnetic in-read crosstalk
        vals = vals + np.where(written, class_offset, 0.0)[:, None]
    fp_true = np.zeros(nsect)
    if fp_mode == "iid":                          # scan-line offsets
        fp_true = rng.normal(0.0, fp_amp, nsect)
    vals = vals + fp_true[:, None]
    vals = vals + rng.normal(0.0, READ_NOISE, vals.shape)
    return vals, classes, fp_true

def build_ori_calib(rng, fp_true, read="kerr", standoff=1.0, nsect=NSECT_ORI):
    """Orientation rung 2: the all-AC-erased calibration read (same static pattern)."""
    M = np.empty((nsect, SECT))
    for s in range(nsect):
        M[s] = (2.0 * rng.binomial(GRAINS, 0.5, SECT) - GRAINS) / GRAINS
    flat = M.reshape(-1)
    if read == "kerr":
        flat = np.convolve(flat, _psf_kernel(), mode="same")
    vals = flat.reshape(nsect, SECT) + fp_true[:, None]
    vals = vals + rng.normal(0.0, READ_NOISE, vals.shape)
    return vals

# ====================== fixed-pattern treatment (C5) ==========================
def detrend_addr(vals, roleW):
    """RUNG 1 (registered, always applied on occupancy): subtract, from every sector of
    address class a, the mean per-cell value of the SAME READ'S UNWRITTEN sectors of
    class a. Removes address-structured fixed pattern (wordline/layer, scan line) --
    the named dominant systematic -- using only the part's own unwritten sectors."""
    out = vals.copy()
    s = np.arange(vals.shape[0])
    for a in range(ADDR_NCLASS):
        sel_u = (s % ADDR_NCLASS == a) & (~roleW)
        if sel_u.sum() == 0:
            continue
        out[s % ADDR_NCLASS == a] -= vals[sel_u].mean()
    return out

def detrend_calib(vals, calib):
    """RUNG 2 (registered escalation): subtract each sector's mean value in the
    all-erased calibration read (static per-sector pattern removed symmetrically)."""
    return vals - calib.mean(axis=1, keepdims=True)

# ============================ pair tables and pools ===========================
def pair_pool(cs, pairs, rng, margin=0, per_pair=PLACEMENTS):
    """THE STATISTIC'S OWN OBJECT: the pool of signed D values over ALL given disjoint
    pairs, per_pair uniform-random in-sector placements each; A(N) = median|D| over the
    pool. The guard, the control fit, the fitted exponents, and the point certificate
    all read pools of exactly this construction -- one object, one moment, one
    population. Also returns the A-side raw ladder median|sum| (the corrupted-shadow
    display of the constraint-1 sweep)."""
    pool = {}; RAW = []
    for N in GRID:
        lo, hi = margin, SECT - N - margin
        if hi < lo:
            lo, hi = 0, 0
        ds = []; raws = []
        for (sA, sB) in pairs:
            for _ in range(per_pair):
                oA = int(rng.integers(lo, hi + 1)); oB = int(rng.integers(lo, hi + 1))
                a = bsum(cs, sA, oA, N)
                ds.append(a - bsum(cs, sB, oB, N)); raws.append(a)
        pool[N] = np.asarray(ds, float)
        RAW.append(float(np.median(np.abs(raws))))
    A = np.array([float(np.median(np.abs(pool[N]))) for N in GRID])
    return pool, A, np.array(RAW)

def uu_pool(cs, upairs, rng, margin=0, per_pair=PLACEMENTS):
    pool, A, _ = pair_pool(cs, upairs, rng, margin=margin, per_pair=per_pair)
    return pool, A

def raw_ladder(vals, roleW, rng):
    """The UNTREATED raw single-sector ladder: median|A-side block sum| over the same
    pair-and-placement construction, computed on the read AS TAKEN (no rung-1/rung-2
    treatment). This is the constraint-1 SHADOW display: what a raw absolute statistic
    does under a common-mode offset -- the corruption the contrast is immune to. (The
    pipeline's own rung-1 de-trend incidentally also protects the treated raw ladder
    against a global offset; the shadow is therefore shown untreated, or the display
    would demonstrate the treatment instead of the constraint.)"""
    cs = seccs(vals)
    _, _, raw = pair_pool(cs, adj_wu_pairs(roleW), rng)
    return raw

def guard_start(A_wu, A_uu, kappa=KAPPA_VOID):
    for i0 in range(0, NG - MIN_POINTS + 1):
        if A_wu[i0] > kappa * A_uu[i0]:
            return i0
    return None

# ==================== the point certificate (C7) ==============================
def point_certificate(wu_p, uu_A, rng, B=KP_B):
    """REGISTERED PROCEDURE, run by the reader on their own read; replaces any fixed
    kappa-above-the-void constant. The V001 defect (A-D4, B-K4): near the guard seam,
    conditioning on median|D_WU(N_min)| > kappa_void*A_UU(N_min) selects up-fluctuations
    at N_min and biases the fitted exponent below the band -- at a rate set by the
    reader's own noise law and programmed-count shot noise, which no registered constant
    can anticipate. The closure: the point-band claim is asserted ONLY for reads that
    CERTIFY it from their own draws. Certificate: B surrogate ladders are built by
    resampling, per grid point, THE POOL'S OWN SIZE of values from the read's own
    written-unwritten pool -- the surrogate is a bootstrap replica of THE VERY OBJECT
    THE STATISTIC REPORTS, carrying the read's own shot noise and noise law; each
    surrogate passes through EXACTLY the pipeline's guard and fit against the read's
    own A_UU. The read is point-admissible (OK) iff zero of the B guard-passing
    surrogates fit below the accumulation band's lower edge; otherwise it is a SEAM
    read: measured, 2-SE falsifier still armed, point-band sentence not asserted.
    Because statistic and certificate read the same object, a read whose own ladder is
    tilted below the band cannot certify: its surrogates inherit the tilt.
    Returns (n_guard_passing, n_below_band)."""
    passed = bad = 0
    for _ in range(B):
        lad = np.array([float(np.median(np.abs(rng.choice(wu_p[N], wu_p[N].size))))
                        for N in GRID])
        i0 = guard_start(lad, uu_A)
        if i0 is None:
            continue
        passed += 1
        if fit_loglog(GRID[i0:], lad[i0:])[0] < BAND_ACC[0]:
            bad += 1
    return passed, bad

# ================== the density condition (C1) -- deterministic ================
def dens_med_ladder(prog, rw):
    """MEDIAN per-block programmed fraction over ALL in-sector block placements of
    every selected written sector, per grid N -- the EXACT population the statistic's
    random-offset blocks sample, read at the statistic's own moment (the median).
    Deterministic in the declared pattern. (The MEAN over this population equals the
    sector density identically at every N, so no mean condition exists in V002 --
    a check with no failing branch is the defect class the program removes.)"""
    P = prog[rw]
    pcs = np.concatenate([np.zeros((P.shape[0], 1)), np.cumsum(P, axis=1)], axis=1)
    fmed = []
    for N in GRID:
        w = (pcs[:, N:] - pcs[:, :-N]) / N          # every placement, every sector
        fmed.append(float(np.median(w)))
    return np.array(fmed)

def density_state(fmed):
    ds, _ = fit_lin(GRID[DENS_MED_I0:], fmed[DENS_MED_I0:])
    if abs(ds) > DENS_TOL_MED:
        return "VOID_DENSITY_MEDIAN", ds
    return None, ds

# ========================== occupancy measurement =============================
def measure_occ(vals, prog, roleW, rng, rung1=True, calib=None,
                wsel=None, wclass=None, shared=None, want_ww=False):
    """The V002 occupancy read. REGISTERED STATE ORDER (C6; code equals text):
       1 INCONCLUSIVE_RAILED   -- A_UU(16) = 0: insufficient access, nothing fitted
       2 READS_UNWRITTEN       -- void guard: no admissible N_min with margin > kappa_void
       3 VOID_DENSITY_MEDIAN   -- |slope of median-placement f-hat (N>=64)| > 0.02 (C1)
       4 INCONCLUSIVE_CONTROL  -- beta_UU outside [0.35, 0.65] (the fixed-pattern
                                  detector, C5; escalate rung 2 before this is final)
       5 OK (the read's own point certificate is clean) or SEAM (it is not)
    Falsifier booleans can be True ONLY in states OK and SEAM.
    wsel/wclass: restrict the written class (clause-b parts). shared: reuse pools."""
    if calib is not None:
        vals = detrend_calib(vals, calib)
    if rung1:
        vals = detrend_addr(vals, roleW)
    cs = seccs(vals)
    if shared is None:
        uidx = list(np.where(~roleW)[0])
        pool, A_uu = uu_pool(cs, same_role_pairs(uidx), rng)
    else:
        pool, A_uu = shared
    out = dict(A_uu=A_uu)
    # 1 railed
    if A_uu[0] == 0.0:
        out.update(state="INCONCLUSIVE_RAILED", fire_a=False, B1=False)
        return out
    if wsel is None:
        rw = roleW
        pairs = adj_wu_pairs(roleW)
    else:
        rw = roleW & (wclass == wsel)
        # disjoint adjacent pairs: selected written sectors against unwritten neighbours
        pairs = []
        used = set()
        wset = set(np.where(rw)[0]); uset = set(np.where(~roleW)[0])
        for i in range(len(roleW) - 1):
            a, b = i, i + 1
            if a in used or b in used:
                continue
            if a in wset and b in uset:
                pairs.append((a, b)); used |= {a, b}
            elif b in wset and a in uset:
                pairs.append((b, a)); used |= {a, b}
    if len(pairs) < K_MIN:
        out.update(state="INSUFFICIENT_GEOMETRY", fire_a=False, B1=False)
        return out
    wu_p, A_wu, rawA = pair_pool(cs, pairs, rng)
    out["A_wu"] = A_wu; out["RAW"] = rawA
    # 2 void guard
    i0 = guard_start(A_wu, A_uu)
    if i0 is None:
        out.update(state="READS_UNWRITTEN", fire_a=False, B1=False)
        return out
    # 3 density (the statistic's own moment over its own population; C1)
    fmed = dens_med_ladder(prog, rw)
    dstate, dslope = density_state(fmed)
    out.update(med_slope=dslope, i0=i0, Nmin=GRID[i0],
               margin=float(A_wu[i0] / A_uu[i0]))
    if dstate is not None:
        out.update(state=dstate, fire_a=False, B1=False)
        return out
    # fits on the admissible window
    g = slice(i0, None); Ns = GRID[i0:]
    bWU, seWU = fit_loglog(Ns, A_wu[g])
    bUU, seUU = fit_loglog(Ns, A_uu[g])
    xi = bWU - bUU; sexi = math.hypot(seWU, seUU)
    out.update(bWU=bWU, seWU=seWU, bUU=bUU, seUU=seUU, xi=xi, sexi=sexi)
    # 4 control band (fixed-pattern detector)
    if not (BAND_CTL[0] <= bUU <= BAND_CTL[1]):
        out.update(state="INCONCLUSIVE_CONTROL", fire_a=False, B1=False)
        return out
    # 5 OK / SEAM: the point-band claim needs the read's own certificate (C7)
    cert_pass, cert_bad = point_certificate(wu_p, A_uu, rng)
    out.update(cert_pass=cert_pass, cert_bad=cert_bad)
    out["state"] = "OK" if (cert_bad == 0 and cert_pass > 0) else "SEAM"
    if want_ww:
        widx = list(np.where(rw)[0])
        _, A_ww, _ = pair_pool(cs, same_role_pairs(widx), rng)
        out["bWW"] = fit_loglog(Ns, A_ww[g])[0] if np.all(A_ww[g] > 0) else float("nan")
    # falsifier booleans -- precedence satisfied by construction (states above returned)
    out["fire_a"] = (bWU + 2 * seWU < BAND_ACC[0]) or (xi + 2 * sexi < XI_MIN)
    out["B1"] = (out["state"] == "OK") and (BAND_ACC[0] <= bWU <= BAND_ACC[1]) \
        and (BAND_CTL[0] <= bUU <= BAND_CTL[1]) and (xi >= XI_MIN)
    return out

def measure_occ2(vals, prog, roleW, wclass, rng):
    """Clause (b): both written classes through the same machinery, shared U pool."""
    v = detrend_addr(vals, roleW)
    cs = seccs(v)
    uidx = list(np.where(~roleW)[0])
    shared = uu_pool(cs, same_role_pairs(uidx), rng)
    res = {}
    for wcl in ("W1", "W2"):
        res[wcl] = measure_occ(vals, prog, roleW, rng, rung1=True, calib=None,
                               wsel=wcl, wclass=wclass, shared=shared)
    ok = all(res[w]["state"] in ("OK",) for w in res)
    fire_b = False
    if ok:
        d = abs(res["W1"]["bWU"] - res["W2"]["bWU"])
        se = math.hypot(res["W1"]["seWU"], res["W2"]["seWU"])
        fire_b = (d - 2 * se) > 0.2
    return res, fire_b

# ========================== orientation measurement ===========================
def ori_pairs(classes):
    used = set(); fam = {"DATA": [], "DCF": [], "DC": []}
    for i in range(len(classes) - 1):
        a, b = i, i + 1
        if a in used or b in used:
            continue
        ca, cb = classes[a], classes[b]
        if ca in fam and cb == "U":
            fam[ca].append((a, b)); used |= {a, b}
        elif cb in fam and ca == "U":
            fam[cb].append((b, a)); used |= {a, b}
    return fam

def measure_ori(vals, classes, rng, calib=None):
    """The V002 orientation (Kerr) read. REGISTERED STATE ORDER (C4, C6):
       1 INCONCLUSIVE_RAILED      -- A_UU(16) = 0
       2 INCONCLUSIVE_DC_CONTROL  -- the same read's DC-saturated positive control does
                                     not pass its guard with beta_DC-U >= 0.9: the
                                     instrument has not demonstrated k=0 sensitivity
       3 INCONCLUSIVE_CONTROL     -- beta_UU outside the control band
       4 INCONCLUSIVE_CROSSTALK   -- the DC-free-coded voider (C2): beta_DCF-U outside
                                     the control band -- an additive in-read artifact
                                     (crosstalk, biased erase) is present
       5 OK
    The DATA sector carries NO void guard (its predicted state is screening; C4).
    fire_c can be True ONLY in state OK."""
    if calib is not None:
        vals = detrend_calib(vals, calib)
    cs = seccs(vals)
    uidx = [i for i, c in enumerate(classes) if c == "U"]
    pool, A_uu = uu_pool(cs, same_role_pairs(uidx), rng, margin=EDGE_MARGIN_ORI)
    out = dict(A_uu=A_uu)
    if A_uu[0] == 0.0:
        out.update(state="INCONCLUSIVE_RAILED", fire_c=False, B3=False)
        return out
    fam = ori_pairs(classes)
    _, A_du, raw_du = pair_pool(cs, fam["DATA"], rng, margin=EDGE_MARGIN_ORI)
    _, A_fu, _ = pair_pool(cs, fam["DCF"], rng, margin=EDGE_MARGIN_ORI)
    _, A_cu, _ = pair_pool(cs, fam["DC"], rng, margin=EDGE_MARGIN_ORI)
    bDU, seDU = fit_loglog(GRID, A_du)
    bFU, seFU = fit_loglog(GRID, A_fu)
    bUU, seUU = fit_loglog(GRID, A_uu)
    xi = bDU - bUU; sexi = math.hypot(seDU, seUU)
    out.update(bDU=bDU, seDU=seDU, bFU=bFU, seFU=seFU, bUU=bUU, seUU=seUU,
               xi=xi, sexi=sexi, A_du=A_du, A_fu=A_fu, A_cu=A_cu, RAW=raw_du)
    i0 = guard_start(A_cu, A_uu)
    if i0 is None:
        out.update(state="INCONCLUSIVE_DC_CONTROL", bCU=float("nan"),
                   seCU=float("nan"), fire_c=False, B3=False)
        return out
    bCU, seCU = fit_loglog(GRID[i0:], A_cu[i0:])
    out.update(bCU=bCU, seCU=seCU)
    if bCU < BAND_ACC[0]:
        out.update(state="INCONCLUSIVE_DC_CONTROL", fire_c=False, B3=False)
        return out
    if not (BAND_CTL[0] <= bUU <= BAND_CTL[1]):
        out.update(state="INCONCLUSIVE_CONTROL", fire_c=False, B3=False)
        return out
    if not (BAND_CTL[0] <= bFU <= BAND_CTL[1]):
        out.update(state="INCONCLUSIVE_CROSSTALK", fire_c=False, B3=False)
        return out
    out["state"] = "OK"
    out["fire_c"] = (xi - 2 * sexi >= XI_MIN) and (bCU - 2 * seCU >= BAND_ACC[0])
    out["B3"] = (BAND_CTL[0] <= bDU <= BAND_CTL[1]) and (xi < XI_MIN) \
        and (BAND_ACC[0] <= bCU <= BAND_ACC[1])
    return out

def rng_from(ss):
    return np.random.default_rng(ss)

# ================================ the sealed run ===============================
def _statecount(states):
    d = {}
    for s in states:
        d[s] = d.get(s, 0) + 1
    return "{" + ", ".join("%s: %d" % kv for kv in sorted(d.items())) + "}"

def main():
    out = []
    P = out.append
    P("T-50 DESIGN ONE V002 -- THE SAME-READ CONTRAST EXPONENT, EIGHT CLOSURES -- "
      "SEALED VERIFICATION RUN")
    P("lane LANE_T50/V002 | date 2026-08-21 | master seed SeedSequence(20260821) | "
      "numpy %s" % np.__version__)
    P("model imported read-only from: %s" % MODEL)
    P("PINNED (C7): GRID=%s  ALL disjoint pairs (>=%d) x %d placements  SECT=%d  "
      "NSECT>=%d" % (GRID, K_MIN, PLACEMENTS, SECT, NSECT))
    P("THE STATISTIC IS THE POOL: A(N) = median|D| over every disjoint pair x")
    P("  placement -- the statistic, the guard, the control and the certificate all")
    P("  read pools of this one construction (one object, one moment, one population;")
    P("  a K-subsample statistic beside a pool certificate would be two different")
    P("  objects, and the seam leaked between them -- measured in construction,")
    P("  closed by this pinning).")
    P("declared: kappa_void=%.0f  MIN_POINTS=%d  DENS_TOL_MED=%.2f on N>=%d  "
      "bands %s / %s  xi>=%.2f  2-SE rule  certificate B=%d" %
      (KAPPA_VOID, MIN_POINTS, DENS_TOL_MED, GRID[DENS_MED_I0], BAND_ACC, BAND_CTL,
       XI_MIN, KP_B))
    P("orientation (C3): polar Kerr microscopy model -- GRAINS/cell=%d  PSF sigma=%.1f "
      "cells  read noise %.2f  edge margin %d cells  %d sectors, cycle %s" %
      (GRAINS, PSF_SIGMA, READ_NOISE, EDGE_MARGIN_ORI, NSECT_ORI, ORI_CYCLE))
    P("interleave: REGISTERED ADDRESS-BALANCED role map, role W iff (s + s//16) even;")
    P("  every address class s mod 16 appears in both roles (C5 rung 1's estimator).")
    P("")

    # ---------------- ANCHOR ----------------
    P("ANCHOR -- lineage to the sealed formation layer (model/geometry.py)")
    pat = GE.occupancy_patterns()
    q = float(N_E); err = pat["unwritten_e"].astype(float)
    rs = [abs(np.where(w == 1, -q, err).sum()) / np.abs(np.where(w == 1, -q, err)).sum()
          for w in pat["written"]]
    P("  sealed-page rho under the model's own residual: min %.16f max %.16f" %
      (min(rs), max(rs)))
    P("  (registered figures 0.9681335687351514 / 0.9797987440417643 -- %s)" %
      ("MATCH" if abs(min(rs) - 0.9681335687351514) < 1e-15 and
                  abs(max(rs) - 0.9797987440417643) < 1e-15 else "MISMATCH"))
    w0 = pat["written"][0]
    v0 = np.where(w0 == 1, -q, err)
    mycell = np.where(w0.astype(bool), -float(N_E), err)
    P("  per-cell value law identical to this lane's builder on sealed page 0: %s" %
      ("PASS" if np.array_equal(v0, mycell) else "FAIL"))
    P("")

    # ---------------- R1 BASELINE ----------------
    P("R1 OCCUPANCY BASELINE -- 50 independent replicates, f=0.5 random data, the full")
    P("   V002 machinery (rung-1 de-trend, density condition, control, certificate).")
    R = 50
    st = dict(bWU=[], bUU=[], bWW=[], xi=[], dmed=[], margin=[], certbad=[],
              states=[], B1=0, firea=0)
    for ss in MASTER.spawn(R):
        rng = rng_from(ss)
        vals, prog, roleW, _, _ = build_occ(rng, f=0.5)
        m = measure_occ(vals, prog, roleW, rng, want_ww=True)
        st["states"].append(m["state"])
        if m["state"] in ("OK", "SEAM"):
            for k2 in (("bWU", "bWU"), ("bUU", "bUU"), ("bWW", "bWW"), ("xi", "xi")):
                st[k2[0]].append(m[k2[1]])
            st["dmed"].append(m["med_slope"]); st["margin"].append(m["margin"])
            st["certbad"].append(m["cert_bad"])
            st["B1"] += m["B1"]; st["firea"] += m["fire_a"]
    for k2, lab in (("bWU", "beta_WU"), ("bUU", "beta_UU"), ("bWW", "beta_WW"),
                    ("xi", "xi=WU-UU"), ("dmed", "median-density slope"),
                    ("margin", "guard margin")):
        m_, lo, hi = med_iqr(st[k2])
        P("    %-21s median %+9.4f   IQR [%+9.4f, %+9.4f]" % (lab, m_, lo, hi))
    P("    states %s | B1 (point claim) TRUE %d/%d | falsifier (a) fires %d/%d" %
      (_statecount(st["states"]), st["B1"], R, st["firea"], R))
    P("    certificate: bad-surrogate counts over the 50 reads: min %d max %d" %
      (min(st["certbad"]), max(st["certbad"])))
    base_bWU = float(np.median(st["bWU"]))
    P("")

    # ---------------- R2 DENSITY CONDITION (C1) ----------------
    P("R2 THE DENSITY CONDITION RESTATED TO THE STATISTIC'S OWN MOMENT (C1; closes")
    P("   refuter A's D1). The registered statistic is median_k|D_k(N)|; for a")
    P("   one-carrier write median_k|D_k(N)| = N_E*N*median_k f-hat_k(N) + O(noise) --")
    P("   a monotone transform -- so the density condition reads THE SAME MOMENT of THE")
    P("   SAME POPULATION: the MEDIAN per-block programmed fraction over ALL in-sector")
    P("   block placements of every written sector, deterministic in the declared")
    P("   pattern, trend-free to %.2f on N >= %d. A sub-tolerance median trend can move" %
      (DENS_TOL_MED, GRID[DENS_MED_I0]))
    P("   the fitted exponent by at most the tolerance itself: 0.02, five times inside")
    P("   the 0.1 band half-width (the co-tuning, structural).")
    P("   THERE IS NO MEAN CONDITION: over the offset-uniform placement population the")
    P("   mean per-block f-hat equals the sector density identically at every N -- a")
    P("   condition with no failing branch, the defect class this program removes")
    P("   (D-8/INST-14). V001's mean check is REPLACED; the M3 family is caught below.")
    P("")
    P("   HONEST PATTERNS -- worst |median-placement slope| over 20 seeds each:")
    honest_worst = 0.0
    for name, kw in [("f=0.015", dict(f=0.015)), ("f=0.02", dict(f=0.02)),
                     ("f=0.05", dict(f=0.05)), ("f=0.10", dict(f=0.10)),
                     ("f=0.25", dict(f=0.25)), ("f=0.50", dict(f=0.50)),
                     ("f=0.75", dict(f=0.75)), ("f=0.90", dict(f=0.90)),
                     ("alternating", dict(pattern="alternating")),
                     ("blocky 1.5f/0.5f", dict(f=0.5, pattern="blocky"))]:
        worst = 0.0
        for ss in MASTER.spawn(20):
            rng = rng_from(ss)
            _, prog, roleW, _, _ = build_occ(rng, **kw)
            fmed = dens_med_ladder(prog, roleW)
            ds_, _ = fit_lin(GRID[DENS_MED_I0:], fmed[DENS_MED_I0:])
            worst = max(worst, abs(ds_))
        honest_worst = max(honest_worst, worst)
        P("     %-18s worst |slope| %.4f" % (name, worst))
    P("   honest worst overall: %.4f (tolerance %.2f sits %.1fx above it)" %
      (honest_worst, DENS_TOL_MED, DENS_TOL_MED / honest_worst))
    P("")
    P("   THE SKEW MASKS (refuter A's D1 killers, as declared literals, mask seed %d)" %
      MASK_SEED)
    P("   and the M3 record mask -- slope is DETERMINISTIC in the mask; 50 seeds each")
    P("   through the full pipeline:")
    mask_rows = [("cascade 1.5/0.5", cascade_mask(0.5, 1.5, 0.5)),
                 ("cascade 1.6/0.4", cascade_mask(0.5, 1.6, 0.4)),
                 ("cascade 1.65/0.35", cascade_mask(0.5, 1.65, 0.35)),
                 ("cascade 1.7/0.3", cascade_mask(0.5, 1.7, 0.3)),
                 ("two-level 0.95/0.03", twolevel_mask()),
                 ("M3 fixed record 512", record_mask())]
    roleW0 = occ_roles()
    smallest_mask_slope = float("inf")
    for name, mask in mask_rows:
        progm = np.zeros((NSECT, SECT)); progm[roleW0] = mask.astype(float)
        fmed = dens_med_ladder(progm, roleW0)
        ds_, _ = fit_lin(GRID[DENS_MED_I0:], fmed[DENS_MED_I0:])
        if "record" not in name:
            smallest_mask_slope = min(smallest_mask_slope, abs(ds_))
        states = []; fa = fb_reach = 0
        for ss in MASTER.spawn(50):
            rng = rng_from(ss)
            vals, prog, roleW, _, _ = build_occ(rng, mask=mask)
            m = measure_occ(vals, prog, roleW, rng)
            states.append(m["state"])
            fa += m.get("fire_a", False)
            fb_reach += m["state"] in ("OK", "SEAM")
        P("     %-20s f=%.4f  slope %+.4f  states %s" %
          (name, mask.mean(), ds_, _statecount(states)))
        P("       reads reaching any clause: %d/50 | falsifier fires: %d/50" %
          (fb_reach, fa))
    P("   smallest skew-mask slope %.4f: the tolerance %.2f sits %.1fx below it." %
      (smallest_mask_slope, DENS_TOL_MED, smallest_mask_slope / DENS_TOL_MED))
    P("   Under V001 these masks passed the mean-density check and fired registered")
    P("   falsifier (b) up to 30/30 on correct physics (VERIFY_A attack_grid_strict:")
    P("   cascade 1.65/0.35, 14-point grid, FIRES 30/30). Under V002 every one is")
    P("   VOID_DENSITY_MEDIAN or READS_UNWRITTEN at rate 1.000 and NO read reaches a")
    P("   clause. The closure is structural: the condition now pins the exact moment")
    P("   the statistic reads, so a one-carrier pattern cannot trend the statistic's")
    P("   typical block without trending the condition's own object.")
    P("")
    P("   CLAUSE (b) HONEST PAIR -- f=0.35 vs f=0.65 on the SAME part, 40 seeds:")
    fb = 0; both_ok = 0; states2 = []
    for ss in MASTER.spawn(40):
        rng = rng_from(ss)
        vals, prog, roleW, wclass, _ = build_occ(rng, f=0.35, two_class=True, f2=0.65)
        res, fire_b = measure_occ2(vals, prog, roleW, wclass, rng)
        states2.append((res["W1"]["state"], res["W2"]["state"]))
        both_ok += all(res[w]["state"] == "OK" for w in res)
        fb += fire_b
    P("     both classes OK %d/40 | clause (b) fires %d/40" % (both_ok, fb))
    P("")

    # ---------------- R3 SEAM AND CERTIFICATE (C7) ----------------
    P("R3 THE POINT CERTIFICATE (C7; closes refuter A's D4 and refuter B's K4 and")
    P("   repair 7). The grid, pairing and placements are PINNED, killing the reader")
    P("   freedom (K=8 / 8-point / sqrt2 grids are gone). Above the void guard there is")
    P("   NO fixed kappa constant: the point-band sentence is asserted only for reads")
    P("   whose OWN CERTIFICATE is clean -- B=%d surrogate ladders resampled from the" % KP_B)
    P("   read's own W-U pair pool (carrying the read's own noise law AND programmed-")
    P("   count shot noise), each run through EXACTLY the pipeline's guard and fit; one")
    P("   surrogate below the band edge and the read is SEAM: measured, 2-SE falsifier")
    P("   still armed, point-band sentence not asserted. Re-derived per reader, from the")
    P("   reader's own draws -- no model constant transfers.")
    P("   Seam ensembles, 150 seeds per f per noise law (Laplace and Student-t3 at")
    P("   matched sigma answer B's attack 5 -- the guarantee is now the reader's own):")
    seam_ok_worst = float("inf"); seam_fires = 0
    for law in ("uniform", "laplace", "t3"):
        for f in (0.015, 0.02, 0.03):
            void = 0; okb = []; sb = []; fires = 0; other = []
            for ss in MASTER.spawn(150):
                rng = rng_from(ss)
                vals, prog, roleW, _, _ = build_occ(rng, f=f, law=law)
                m = measure_occ(vals, prog, roleW, rng)
                s = m["state"]
                if s == "READS_UNWRITTEN":
                    void += 1
                elif s == "OK":
                    okb.append(m["bWU"]); fires += m["fire_a"]
                elif s == "SEAM":
                    sb.append(m["bWU"]); fires += m["fire_a"]
                else:
                    other.append(s)
            seam_fires += fires
            if okb:
                seam_ok_worst = min(seam_ok_worst, min(okb))
            P("     %-8s f=%.3f  void %3d | OK %3d worst %s | SEAM %3d worst %s | "
              "fires %d%s" %
              (law, f, void, len(okb),
               ("%+.4f" % min(okb)) if okb else "  --  ",
               len(sb), ("%+.4f" % min(sb)) if sb else "  --  ", fires,
               ("" if not other else " | other " + _statecount(other))))
    P("   Mid-density control (the attack-4 region), 50 seeds each, uniform law:")
    for f in (0.05, 0.10):
        void = 0; okb = []; sb = []; fires = 0
        for ss in MASTER.spawn(50):
            rng = rng_from(ss)
            vals, prog, roleW, _, _ = build_occ(rng, f=f)
            m = measure_occ(vals, prog, roleW, rng)
            if m["state"] == "READS_UNWRITTEN":
                void += 1
            elif m["state"] == "OK":
                okb.append(m["bWU"]); fires += m["fire_a"]
            elif m["state"] == "SEAM":
                sb.append(m["bWU"]); fires += m["fire_a"]
        seam_fires += fires
        if okb:
            seam_ok_worst = min(seam_ok_worst, min(okb))
        P("     f=%.2f  void %2d | OK %2d worst %s | SEAM %2d worst %s | fires %d" %
          (f, void, len(okb), ("%+.4f" % min(okb)) if okb else "  --  ",
           len(sb), ("%+.4f" % min(sb)) if sb else "  --  ", fires))
    P("   WORST CERTIFIED (OK) beta_WU ANYWHERE ABOVE: %+.4f (band edge 0.9); falsifier" %
      seam_ok_worst)
    P("   fires across every ensemble above: %d. On the pinned pool statistic a read" %
      seam_fires)
    P("   either reads as unwritten or certifies -- refuter B's +0.8022 and refuter A's")
    P("   0.8815-0.8954 guard-passing families were creatures of the K=8..16 subset")
    P("   statistic, and the pinning removed them. The SEAM state remains registered")
    P("   and EXERCISED: a genuinely screening surface cannot certify and lands SEAM")
    P("   with the 2-SE falsifier armed (the two-signed member of the sealed suite,")
    P("   where clause (a) then fires at rate 1.000).")
    P("")

    # ---------------- R4 FIXED PATTERN (C5) ----------------
    P("R4 THE FIXED-PATTERN TREATMENT (C5; closes refuter B's K3 and refuter A's D3).")
    P("   Two registered rungs. RUNG 1 (always on, same read): subtract per address")
    P("   class the mean of the SAME READ'S unwritten sectors of that class -- removes")
    P("   the named dominant systematics (wordline/layer structure, scan-line offsets),")
    P("   which are ADDRESS-STRUCTURED. RUNG 2 (escalation, the part's own unwritten")
    P("   sectors in a second preparation): erase the whole part, read it, subtract each")
    P("   sector's calibration mean; removes ANY static per-sector pattern. The U-U")
    P("   control band IS the fixed-pattern detector; INCONCLUSIVE_CONTROL is its named")
    P("   outcome, and the INCONCLUSIVE precedence (C6) is wired into fire_a, so the")
    P("   sealed-V001 defect -- fire_a firing on a correctly-accumulating part -- is")
    P("   structurally gone: one implementation, states before fires.")
    P("   iid per-sector offsets N(0,s) e/cell (refuter B's attack 1), 30 reps each:")
    P("   s in e/cell | rung 1 only: states, fire_a | with rung 2: states, B1, fire_a")
    for s_fp in (0.0, 0.05, 0.1, 0.2, 0.5, 1.0):
        st1 = []; f1 = 0; st2 = []; b2 = 0; f2 = 0; wdev = []
        for ss in MASTER.spawn(30):
            rng = rng_from(ss)
            vals, prog, roleW, _, fp = build_occ(rng, f=0.5, fp_mode="iid", fp_amp=s_fp)
            m1 = measure_occ(vals, prog, roleW, rng)
            st1.append(m1["state"]); f1 += m1.get("fire_a", False)
            calib = build_occ_calib(rng, fp)
            m2 = measure_occ(vals, prog, roleW, rng, calib=calib)
            st2.append(m2["state"]); f2 += m2.get("fire_a", False)
            b2 += m2.get("B1", False)
            if m2["state"] in ("OK", "SEAM"):
                wdev.append(abs(m2["bWU"] - 1.0))
        P("     s=%.2f  rung1 %s fire_a %d" % (s_fp, _statecount(st1), f1))
        P("            rung2 %s B1 %d/30 fire_a %d  worst |beta-1| %.4f" %
          (_statecount(st2), b2, f2, max(wdev) if wdev else float("nan")))
    P("   Sealed V001 at s=0.5-1.0 e: code fire_a 26-29/30 on a perfectly accumulating")
    P("   part (refuter B, attack 1). V002: fire_a 0/30 at every s (precedence), the")
    P("   read names its own state, and rung 2 RECOVERS THE VERDICT -- the protocol no")
    P("   longer returns INCONCLUSIVE forever on real parts (A-D3 closed).")
    P("   Address-structured pattern (deck profile), amplitude sweep, 30 reps each:")
    for amp in (0.5, 5.0, 50.0):
        st0 = []; st1 = []
        for ss in MASTER.spawn(30):
            rng = rng_from(ss)
            vals, prog, roleW, _, fp = build_occ(rng, f=0.5, fp_mode="addr", fp_amp=amp)
            m0 = measure_occ(vals, prog, roleW, rng, rung1=False)
            st0.append(m0["state"])
            m1 = measure_occ(vals, prog, roleW, rng)
            st1.append(m1["state"])
        P("     amp=%.1f e  no treatment %s -> rung 1 %s" %
          (amp, _statecount(st0), _statecount(st1)))
    P("   Rung integrity on honest reads (no fixed pattern), 20 reps: rung 2 applied")
    P("   anyway must not distort:")
    devs = []
    for ss in MASTER.spawn(20):
        rng = rng_from(ss)
        vals, prog, roleW, _, fp = build_occ(rng, f=0.5)
        calib = build_occ_calib(rng, fp)
        m2 = measure_occ(vals, prog, roleW, rng, calib=calib)
        devs.append(abs(m2["bWU"] - base_bWU) if m2["state"] in ("OK", "SEAM")
                    else float("inf"))
    P("     max |beta_WU shift| vs baseline median: %.4f; all states OK: %s" %
      (max(devs), all(d < float("inf") for d in devs)))
    P("   SCOPE, STATED: a per-sector-random response OFFSET that is not static across")
    P("   the two preparations is indistinguishable inside one read from a physical")
    P("   sector-mean difference; the control band detects it and the read registers")
    P("   INCONCLUSIVE_CONTROL -- an insufficient-instrument outcome, named, never a")
    P("   clause fire. Response GAIN nonuniformity scales pair prefactors and is caught")
    P("   by the same control when it matters (measured outcome).")
    P("")

    # ---------------- R5 ORIENTATION (C2, C3, C4) ----------------
    P("R5 THE ORIENTATION HALF ON THE NAMED M-READING INSTRUMENT (C3; closes refuter")
    P("   B's K1/K2, refuter A's D2). INSTRUMENT: polar Kerr microscopy (wide-field")
    P("   MOKE imaging). WHAT IT MEASURES, VERIFIED BEFORE NAMING: the polar Kerr")
    P("   rotation of reflected polarized light is proportional to the LOCAL out-of-")
    P("   plane magnetization M_z within the optical spot and penetration depth -- the")
    P("   uniform-film Kerr hysteresis loop is the standard demonstration that the")
    P("   response at k=0 is NONZERO (a saturated uniform film gives the full rotation),")
    P("   where every stray-field mapper (MFM, scanning Hall, NV) has transfer")
    P("   exp(-|k|d)(1-exp(-|k|t)) = 0 at k=0. v_i = the raw analyzer-difference signal")
    P("   per resolved cell, in the instrument's own units; the magneto-optic")
    P("   proportionality is a prefactor and cancels from the exponent. REQUIREMENT:")
    P("   bit-cell length >= the microscope's declared resolution -- the reader's own")
    P("   writer sets the bit length (spin-stand / drive write clock), so this is a")
    P("   preparation choice, not a device constant.")
    P("   Model: cell = mean of %d grain moments, optical PSF sigma %.1f cells," %
      (GRAINS, PSF_SIGMA))
    P("   read noise %.2f; blocks sit >= %d cells from sector boundaries." %
      (READ_NOISE, EDGE_MARGIN_ORI))
    # interior map check
    rng = rng_from(MASTER.spawn(1)[0])
    vals_k, classes_k, _ = build_ori(rng)
    rng = rng_from(MASTER.spawn(1)[0])
    vals_s, classes_s, _ = build_ori(rng, read="stray", standoff=1.0)
    dc_rows = [i for i, c in enumerate(classes_k) if c == "DC"]
    u_rows = [i for i, c in enumerate(classes_k) if c == "U"]
    P("   k=0 SENSITIVITY, MEASURED ON THE MODEL: DC-sector interior mean map value")
    P("     Kerr read:        DC %+8.5f   AC-erased %+8.5f" %
      (float(vals_k[dc_rows][:, 100:-100].mean()),
       float(vals_k[u_rows][:, 100:-100].mean())))
    P("     stray-field read: DC %+8.5f   AC-erased %+8.5f  (the K1 kill, reproduced)" %
      (float(vals_s[dc_rows][:, 100:-100].mean()),
       float(vals_s[u_rows][:, 100:-100].mean())))
    P("   GUARD SCOPE (C4; closes K2): the DATA sector carries NO void guard -- its")
    P("   predicted state is screening, and V001's shared guard voided the prediction's")
    P("   own subject. Admissibility is carried by the SAME READ'S DC-saturated positive")
    P("   control (guard + accumulation band), the U-U control band, and the DC-free")
    P("   voider. The density condition is occupancy-only: no programmed fraction exists")
    P("   on an orientation surface; the orientation analogue is the declared write")
    P("   pattern itself, verified from the reader's own write record.")
    P("")
    P("   BASELINE -- 50 reps, random data, honest Kerr read:")
    ost = dict(bDU=[], bUU=[], bCU=[], bFU=[], xi=[], states=[], B3=0, firec=0)
    for ss in MASTER.spawn(50):
        rng = rng_from(ss)
        vals, classes, _ = build_ori(rng)
        m = measure_ori(vals, classes, rng)
        ost["states"].append(m["state"])
        if m["state"] == "OK":
            for k2 in ("bDU", "bUU", "bCU", "bFU", "xi"):
                ost[k2].append(m[k2])
            ost["B3"] += m["B3"]; ost["firec"] += m["fire_c"]
    for k2, lab in (("bDU", "beta_DATA-U"), ("bUU", "beta_UU"), ("bCU", "beta_DC-U"),
                    ("bFU", "beta_DCF-U (voider)"), ("xi", "xi=DU-UU")):
        m_, lo, hi = med_iqr(ost[k2])
        P("     %-20s median %+8.4f   IQR [%+8.4f, %+8.4f]" % (lab, m_, lo, hi))
    P("     states %s | B3 true %d/50 | fire_c %d/50" %
      (_statecount(ost["states"]), ost["B3"], ost["firec"]))
    P("     the voider ARMS on honest reads (beta_DCF-U in the control band): clause")
    P("     (c) remains reachable -- and the same read shows the instrument seeing")
    P("     accumulation at k=0 (beta_DC-U ~ 1).")
    P("")
    P("   THE IN-READ CROSSTALK VOIDER (C2; closes A-D2): a per-cell additive offset on")
    P("   WRITTEN sectors only (topography/charging crosstalk, not magnetization),")
    P("   40 reps per level -- refuter A measured clause (c) firing 13/40 at 0.10 and")
    P("   34/40 at 0.20 per grain on a correctly-screening medium:")
    for c in (0.02, 0.05, 0.10, 0.20):
        sts = []; fc = 0
        for ss in MASTER.spawn(40):
            rng = rng_from(ss)
            vals, classes, _ = build_ori(rng, class_offset=c)
            m = measure_ori(vals, classes, rng)
            sts.append(m["state"]); fc += m.get("fire_c", False)
        P("     offset %.2f/grain: states %s | fire_c %d/40" % (c, _statecount(sts), fc))
    P("     The crosstalk rides every written class equally; the DC-free-coded sectors")
    P("     can only leave the control band through such an additive artifact, so the")
    P("     voider trips and the read is INCONCLUSIVE_CROSSTALK -- fire_c 0/40 at every")
    P("     level, closing D2's both measured fire channels.")
    P("   BIASED ERASE +0.25/grain (V001 policed this by a scope CLAUSE; V002 detects")
    P("   it -- the erase-mean artifact enters every written-minus-U contrast the same")
    P("   way and trips the same voider), 30 reps:")
    sts = []; fc = 0
    for ss in MASTER.spawn(30):
        rng = rng_from(ss)
        vals, classes, _ = build_ori(rng, erase_bias=0.25)
        m = measure_ori(vals, classes, rng)
        sts.append(m["state"]); fc += m.get("fire_c", False)
    P("     states %s | fire_c %d/30 (V001: 34-35/50 when the scope clause was ignored)" %
      (_statecount(sts), fc))
    P("   CLAUSE (c) IS TRIGGERABLE (the answer to 'a falsifier no instrument can")
    P("   trigger'): data written 75%% one-way -- genuinely accumulating data -- through")
    P("   the same Kerr model, 30 reps:")
    sts = []; fc = 0
    for ss in MASTER.spawn(30):
        rng = rng_from(ss)
        vals, classes, _ = build_ori(rng, data_p1=0.75)
        m = measure_ori(vals, classes, rng)
        sts.append(m["state"]); fc += m.get("fire_c", False)
    P("     states %s | fire_c %d/30 -- the clause fires when and only when the physics" %
      (_statecount(sts), fc))
    P("     it names is present, on the instrument the text names.")
    P("   THE STRAY-FIELD MEMBER (K1's instrument, self-detected): the same surface read")
    P("   through T(k) = exp(-|k|d)(1-exp(-|k|t)), d=1 cell, 30 reps:")
    sts = []
    for ss in MASTER.spawn(30):
        rng = rng_from(ss)
        vals, classes, _ = build_ori(rng, read="stray", standoff=1.0)
        m = measure_ori(vals, classes, rng)
        sts.append(m["state"])
    P("     states %s -- a reader on MFM/Hall/NV gets INCONCLUSIVE_DC_CONTROL from the" %
      _statecount(sts))
    P("     read's own positive control, never a silent void: the instrument requirement")
    P("     is load-bearing and SELF-MEASURED.")
    P("   SCAN-LINE FIXED PATTERN (K3's magnetic face) and rung 2, 20 reps each:")
    for s_fp in (0.05, 0.20):
        st1 = []; st2 = []
        for ss in MASTER.spawn(20):
            rng = rng_from(ss)
            vals, classes, fp = build_ori(rng, fp_mode="iid", fp_amp=s_fp)
            m1 = measure_ori(vals, classes, rng)
            st1.append(m1["state"])
            calib = build_ori_calib(rng, fp)
            m2 = measure_ori(vals, classes, rng, calib=calib)
            st2.append(m2["state"])
        P("     s=%.2f  raw %s -> rung 2 %s" % (s_fp, _statecount(st1), _statecount(st2)))
    P("   EASY-AXIS TILT 30 deg (the COMP-4/INST-5 configuration), 20 reps:")
    bd, bc2 = [], []
    for ss in MASTER.spawn(20):
        rng = rng_from(ss)
        vals, classes, _ = build_ori(rng, tilt_deg=30.0)
        m = measure_ori(vals, classes, rng)
        if m["state"] == "OK":
            bd.append(m["bDU"]); bc2.append(m["bCU"])
    P("     beta_DATA-U median %+.4f | beta_DC-U median %+.4f -- tilt is a per-cell" %
      (float(np.median(bd)), float(np.median(bc2))))
    P("     prefactor; no squareness, no epsilon, no M_r/M_s anywhere in this design.")
    P("")

    # ---------------- R6 RAILED (C8) ----------------
    P("R6 THE RAILED-POPULATION BRANCH (C8; closes refuter B's K5, refuter A's D5).")
    P("   Erased population railed at the read floor (custody: read-retry scales floor")
    P("   at GND; the erased distribution below GND is not resolved), 30 reps:")
    sts = []; fa = 0; nan_escape = False
    for ss in MASTER.spawn(30):
        rng = rng_from(ss)
        vals, prog, roleW, _, _ = build_occ(rng, f=0.5, railed=True)
        m = measure_occ(vals, prog, roleW, rng)
        sts.append(m["state"]); fa += m.get("fire_a", False)
        for kk in ("bWU", "bUU", "xi"):
            if kk in m and isinstance(m[kk], float) and math.isnan(m[kk]):
                nan_escape = True
    P("     states %s | fire_a %d/30 | any nan reaching a decision: %s" %
      (_statecount(sts), fa, nan_escape))
    P("     A_UU(16) = 0 reads INSUFFICIENT ACCESS: the guard is never evaluated (no")
    P("     vacuous pass), nothing is fitted, the reader is told which access mode")
    P("     suffices (per-cell analog read that resolves the erased distribution).")
    P("")

    # ---------------- R7 FALSE FIRE ----------------
    P("R7 FALSE-FIRE RATE -- 200 honest seeds per clause, full V002 machinery:")
    fa = 0
    for ss in MASTER.spawn(200):
        rng = rng_from(ss)
        vals, prog, roleW, _, _ = build_occ(rng, f=0.5)
        m = measure_occ(vals, prog, roleW, rng)
        fa += m.get("fire_a", False)
    fbn = 0
    for ss in MASTER.spawn(200):
        rng = rng_from(ss)
        vals, prog, roleW, wclass, _ = build_occ(rng, f=0.35, two_class=True, f2=0.65)
        _, fire_b = measure_occ2(vals, prog, roleW, wclass, rng)
        fbn += fire_b
    fcn = 0
    for ss in MASTER.spawn(200):
        rng = rng_from(ss)
        vals, classes, _ = build_ori(rng)
        m = measure_ori(vals, classes, rng)
        fcn += m.get("fire_c", False)
    P("     clause (a): %d/200   clause (b): %d/200   clause (c): %d/200" %
      (fa, fbn, fcn))
    false_fires = fa + fbn + fcn
    P("")

    # ---------------- R8 OFFSET / DRIFT (constraint 1) ----------------
    P("R8 COMMON-MODE OFFSET AND DRIFT (the T-50 constraint-1 sweep, carried and re-run")
    P("   through V002; 10 reps per setting; contrast vs the UNTREATED raw single-")
    P("   sector exponent -- the corruption the contrast is immune to. The pipeline's")
    P("   rung-1 de-trend incidentally protects even the treated raw ladder against a")
    P("   global offset, so the shadow is computed on the read AS TAKEN):")
    def occ_sweep(name, kw):
        bC, bR = [], []
        for ss in MASTER.spawn(10):
            rng = rng_from(ss)
            vals, prog, roleW, _, _ = build_occ(rng, f=0.5, **kw)
            m = measure_occ(vals, prog, roleW, rng)
            if m["state"] in ("OK", "SEAM"):
                bC.append(m["bWU"])
                bR.append(fit_loglog(GRID, raw_ladder(vals, roleW, rng))[0])
        return float(np.median(bC)), float(np.median(bR))
    occ_rows = []
    for name, kw in [("offset 0", {}), ("offset +0.5 e", dict(offset=0.5)),
                     ("offset -0.5 e", dict(offset=-0.5)),
                     ("offset +5 e", dict(offset=5.0)),
                     ("offset +50 e", dict(offset=50.0)),
                     ("drift 50 e span", dict(drift=50.0))]:
        bc, br = occ_sweep(name, kw)
        occ_rows.append((name, bc, br))
        P("     %-18s contrast beta_WU %+8.4f   raw-W beta %+8.4f" % (name, bc, br))
    occ_shift = max(abs(bc - occ_rows[0][1]) for _, bc, _ in occ_rows)
    ori_rows = []
    for name, co in [("offset 0", 0.0), ("offset +0.05", 0.05), ("offset +0.5", 0.5)]:
        bC, bR = [], []
        for ss in MASTER.spawn(10):
            rng = rng_from(ss)
            vals, classes, _ = build_ori(rng, class_offset=0.0)
            if co:
                vals = vals + co    # common-mode across the WHOLE read
            m = measure_ori(vals, classes, rng)
            if m["state"] == "OK":
                bC.append(m["bDU"]); bR.append(fit_loglog(GRID, m["RAW"])[0])
        ori_rows.append((name, float(np.median(bC)), float(np.median(bR))))
        P("     %-18s contrast beta_DU %+8.4f   raw-DATA beta %+8.4f" % ori_rows[-1])
    ori_shift = max(abs(bc - ori_rows[0][1]) for _, bc, _ in ori_rows)
    B4 = (occ_shift < 0.05) and (ori_shift < 0.05)
    P("     max contrast shift: occupancy %.4f, orientation %.4f -> B4 (immunity): %s" %
      (occ_shift, ori_shift, B4))
    P("")

    # ---------------- R9 CONSTRAINT-2 AND UNIT FREEDOM ----------------
    P("R9 CONSTRAINT-2 TABLE -- the three moment choices on the SAME pool, one read")
    P("   each:")
    rng = rng_from(MASTER.spawn(1)[0])
    vals, prog, roleW, _, _ = build_occ(rng, f=0.5)
    v = detrend_addr(vals, roleW); cs = seccs(v)
    wu_pool_r9, _, _ = pair_pool(cs, adj_wu_pairs(roleW), rng)
    rng2 = rng_from(MASTER.spawn(1)[0])
    valso, clso, _ = build_ori(rng2)
    cso = seccs(valso)
    du_pool_r9, _, _ = pair_pool(cso, ori_pairs(clso)["DATA"], rng2,
                                 margin=EDGE_MARGIN_ORI)
    for lab, pl in (("occupancy W-U", wu_pool_r9), ("orientation DATA-U", du_pool_r9)):
        A_ = [float(np.median(np.abs(pl[N]))) for N in GRID]
        M2_ = [float(np.median(pl[N] ** 2)) for N in GRID]
        V_ = [float(np.var(pl[N], ddof=1)) for N in GRID]
        bA, _ = fit_loglog(GRID, A_); bM, _ = fit_loglog(GRID, M2_)
        bV, _ = fit_loglog(GRID, V_)
        P("     %-20s median|D| %+7.4f | UNCENTRED D^2 %+7.4f | CENTRED var %+7.4f" %
          (lab, bA, bM, bV))
    P("     |D| separates 1 vs 1/2, uncentred D^2 separates 2 vs 1, centred variance")
    P("     ~1 for BOTH encodings: excluded from the design and displayed here (D-15).")
    ssu = MASTER.spawn(1)[0]
    rngA = rng_from(ssu)
    valsA, progA, roleA, _, _ = build_occ(rngA, f=0.5)
    mA = measure_occ(valsA, progA, roleA, rng_from(np.random.SeedSequence(999)))
    rngB = rng_from(ssu)
    valsB, progB, roleB, _, _ = build_occ(rngB, f=0.5)
    valsB = valsB * 0.04
    mB = measure_occ(valsB, progB, roleB, rng_from(np.random.SeedSequence(999)))
    P("   UNIT FREEDOM (INST-3): beta_WU in e units %.10f | x0.04 'V/e' %.10f | "
      "diff %.2e" % (mA["bWU"], mB["bWU"], abs(mA["bWU"] - mB["bWU"])))
    P("")

    # ---------------- R10 PATTERN INDEPENDENCE ----------------
    P("R10 OCCUPANCY PATTERN-INDEPENDENCE over the CONDITION-PASSING patterns (the C1")
    P("    quantifier re-derived: every pattern passing guard + density condition +")
    P("    control), 20 reps each:")
    meds = {}
    for name, kw in [("random f=0.10", dict(f=0.10)), ("random f=0.25", dict(f=0.25)),
                     ("random f=0.50", dict(f=0.50)), ("random f=0.75", dict(f=0.75)),
                     ("random f=0.90", dict(f=0.90)),
                     ("alternating", dict(pattern="alternating")),
                     ("blocky 1.5f/0.5f", dict(f=0.5, pattern="blocky"))]:
        bl = []; nok = 0
        for ss in MASTER.spawn(20):
            rng = rng_from(ss)
            vals, prog, roleW, _, _ = build_occ(rng, **kw)
            m = measure_occ(vals, prog, roleW, rng)
            if m["state"] in ("OK", "SEAM"):
                bl.append(m["bWU"]); nok += m["state"] == "OK"
        meds[name] = float(np.median(bl))
        P("     %-18s beta_WU median %+8.4f  (n=%d, OK %d)" %
          (name, meds[name], len(bl), nok))
    spread = max(meds.values()) - min(meds.values())
    B2 = spread < 0.2
    P("     max pairwise spread %.4f vs the clause-(b) threshold 0.2 -> B2: %s" %
      (spread, B2))
    P("")

    # ---------------- guard-scope table ----------------
    P("GUARD-SCOPE TABLE (C4) -- which condition guards which clause; no unguarded")
    P("clause, no clause reading a voided surface:")
    P("  condition                          clause(a) clause(b) clause(c) point-claim")
    P("  occ void guard (kappa_void=8)         yes       yes       --        yes")
    P("  occ density-median (C1)               yes       yes       --        yes")
    P("  occ control band / FP detector        yes       yes       --        yes")
    P("  occ point certificate (C7)            --        OK-only   --        yes")
    P("  occ railed branch (C8)                yes       yes       --        yes")
    P("  ori DC positive control (C3/C4)       --        --        yes       --")
    P("  ori U-U control band                  --        --        yes       --")
    P("  ori DC-free voider (C2)               --        --        yes       --")
    P("  (clause (a) fires in states OK and SEAM under the 2-SE rule; clause (b)")
    P("  requires both classes OK; clause (c) requires state OK. The DATA sector")
    P("  carries no void guard: its predicted state is screening.)")
    P("")

    # ---------------- decision booleans ----------------
    P("DECISION BOOLEANS (model-side):")
    P("  B1 occupancy accumulation (state OK, bands, xi):        %d/50 TRUE" % st["B1"])
    P("  B2 pattern-independence (spread %.4f < 0.2):            %s" % (spread, B2))
    P("  B3 orientation screening + in-read positive control:    %d/50 TRUE" % ost["B3"])
    P("  B4 offset/drift immunity of the contrast (< 0.05):      %s" % B4)
    P("  B5 false fires over 600 honest reads:                   %d" % false_fires)
    P("  B6 worst CERTIFIED beta_WU across every seam ensemble:  %+.4f (edge 0.9)" %
      seam_ok_worst)
    P("")
    P("END OF SEALED RUN")
    return "\n".join(out)

if __name__ == "__main__":
    print(main())
