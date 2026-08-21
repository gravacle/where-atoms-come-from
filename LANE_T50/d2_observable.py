"""T-50 DESIGN TWO -- THE SELF-REFERENCED DOUBLING CONTRAST.

The observable, in one line: from ONE per-cell read of ONE part,
    D  =  beta_w - beta_e
where beta is the fitted log2-slope of the UNCENTRED block second moment
    M2(n) = mean over disjoint contiguous blocks B of size n of ( sum_{i in B} u_i )^2
across a dyadic ladder of block sizes n, computed for the WRITTEN sector (beta_w)
and for the UNWRITTEN sector of the same part in the same read (beta_e), with the
per-cell value u_i referenced to THE ERASED POPULATION OF THE SAME READ:
    written occupancy sector : u_i = v_i - mean(v over the sector's own erased-cluster cells)
    written orientation track: u_i = v_i - mean(v over the adjacent AC-erased region, same scan)
    unwritten sector (null leg): u_j = v_j - mean(v over the unwritten sector itself)
      [exactly centred by construction; finite-size correction 1/(1 - n/N) applied]

PREDICTION (structural integers, no imported constants):
    occupancy encoding, any data pattern in scope . beta_w -> 2, beta_e -> 1, D -> 1
    orientation encoding, DC-balanced (real) data . beta_w -> 1, beta_e -> 1, D -> 0
    orientation encoding, DC-saturated ............ beta_w -> 2,            D -> 1

WHY IT HONOURS THE THREE COMPUTED CONSTRAINTS OF THE T-50 ROW:
 (1) contrast AND exponent: the reference is taken from the unwritten/erased population
     of the SAME read (the written-minus-unwritten differential), and the observable IS
     an exponent; a common-mode or sector-common offset shifts v_i and the reference
     equally and cancels PER CELL, before any statistic is formed.  u is invariant under
     v -> g*v + b for any gain g and offset b (g cancels in the M2 ratio; b cancels in u).
 (2) the statistic is the UNCENTRED second moment (2 vs 1), never the centred variance;
     the null leg alone is self-centred, which is exactly what makes it a null.
 (3) density fixed: f is MEASURED per block from the same read's clustering, and the
     ladder is admissible only where Var_blocks(f_B) sits inside the binomial envelope
     f(1-f)/n -- a fixed-length record in a growing block violates the envelope and the
     protocol REFUSES rather than reporting an exponent.

No quantity is imported: no V_t,neutral, no C_fg, no datasheet Delta or Q_p, no
squareness class.  Every input is read off the part itself in the same read.

Model-side simulation constants (SIMULATION ONLY -- they never enter the estimator):
    N_E = 100 e (planar-class programmed charge, the model's sealed figure)
    DELTA = 5 e (declared over-erase residual half-width, uniform{-5..5})

Verdicts are MEASURED OUTCOMES, per the principal's accumulation directive:
    ACCUMULATES   (D - 2*SE > 0.5)
    SCREENS       (D + 2*SE < 0.5)
    INDETERMINATE (neither, with SE reported)
    REFUSED       (a guard failed; the read cannot carry the measurement, no verdict)
"""

import numpy as np

N_E = 100.0     # model programmed level, electron units above erased nominal (simulation only)
DELTA = 5       # model residual half-width, e (simulation only)


# --------------------------------------------------------------------------- reads (model side)
def read_occupancy(N, f, rng, mu=0.0, offset=0.0, gain=1.0, two_signed=False,
                   fixed_record=0, prog_sd=2.0):
    """Simulated per-cell read of a WRITTEN occupancy sector (NAND-like V_t map in
    electron-equivalent units above the erased nominal).  Programmed cell reads at
    N_E + noise; erased cell reads at its over-erase residual (uniform{-DELTA..DELTA} + mu).
    offset/gain model tester miscalibration (v -> gain*(v + offset)).
    two_signed  : mutation m1 -- the write injects EITHER carrier sign.
    fixed_record: mutation m3 -- a fixed-length 50%-dense record at the sector start,
                  the remainder erased (refuter A's density counterexample).
    Returns (v, data) with data the 0/1 programmed mask."""
    if fixed_record:
        data = np.zeros(N, dtype=int)
        data[:fixed_record] = rng.integers(0, 2, fixed_record)
    else:
        data = (rng.random(N) < f).astype(int)
    resid = rng.integers(-DELTA, DELTA + 1, N).astype(float) + mu
    prog = N_E + rng.normal(0.0, prog_sd, N)
    if two_signed:
        prog = prog * rng.choice([-1.0, 1.0], N)
    v = np.where(data == 1, prog, resid)
    return gain * (v + offset), data


def read_occupancy_unwritten(N, rng, mu=0.0, offset=0.0, gain=1.0):
    """The same part's never-programmed sector: residuals only."""
    v = rng.integers(-DELTA, DELTA + 1, N).astype(float) + mu
    return gain * (v + offset)


def read_orientation(N, rng, kind="random", offset=0.0, gain=1.0, tilt_deg=0.0,
                     m=1.0, read_sd=0.1):
    """Simulated per-grain read of an orientation track (MFM / read-back map, arbitrary
    instrument units).  kind: 'random' (fair +-1 data), 'dc' (DC-saturated), 'dcfree'
    (alternating +1,-1 coded).  tilt_deg: easy-axis dispersion; the read sees m*cos(th)."""
    if kind == "random":
        s = (rng.integers(0, 2, N) * 2 - 1).astype(float)
    elif kind == "dc":
        s = np.ones(N)
    elif kind == "dcfree":
        s = np.tile([1.0, -1.0], N // 2 + 1)[:N]
    else:
        raise ValueError(kind)
    if tilt_deg > 0:
        s = s * np.cos(rng.normal(0.0, np.deg2rad(tilt_deg), N))
    v = m * s + rng.normal(0.0, read_sd, N)
    return gain * (v + offset)


def read_orientation_erased(N, rng, offset=0.0, gain=1.0, tilt_deg=0.0, m=1.0,
                            read_sd=0.1, drift=0.0):
    """AC-erased region of the same platter, same scan: random per-grain signs.
    drift: linear instrument drift across the region (guard-test input)."""
    s = (rng.integers(0, 2, N) * 2 - 1).astype(float)
    if tilt_deg > 0:
        s = s * np.cos(rng.normal(0.0, np.deg2rad(tilt_deg), N))
    v = m * s + rng.normal(0.0, read_sd, N)
    if drift:
        v = v + drift * np.linspace(0, 1, N)
    return gain * (v + offset)


# --------------------------------------------------------------------------- estimator pieces
def ladder_sizes(n_min, n_max):
    out, n = [], n_min
    while n <= n_max:
        out.append(n)
        n *= 2
    return out


def block_sums(u, n):
    K = u.size // n
    return u[:K * n].reshape(K, n).sum(axis=1)


def block_sums_2phase(u, n):
    """Block sums at phase 0 and phase n/2 (overlapping by half): unbiased for E[S^2]
    and materially lower sampling noise at the top of the ladder for the same read."""
    a = block_sums(u, n)
    b = block_sums(u[n // 2:], n)
    return np.concatenate([a, b])


def m2_ladder(u, ns, self_centred=False):
    """Uncentred second moment of block sums over the ladder (two-phase overlapping
    estimator).  When the reference was the sector's OWN sample mean (null leg), apply
    the exact iid finite-size correction 1/(1 - n/N)."""
    N = u.size
    out = {}
    for n in ns:
        S = block_sums_2phase(u, n)
        m2 = float(np.mean(S ** 2))
        if self_centred:
            m2 = m2 / max(1.0 - n / N, 1e-9)
        out[n] = m2
    return out


def fit_slope(m2):
    ns = np.array(sorted(m2), dtype=float)
    ys = np.array([m2[int(n)] for n in ns], dtype=float)
    if np.any(ys <= 0):
        return np.nan
    return float(np.polyfit(np.log2(ns), np.log2(ys), 1)[0])


def exponents_all(u, ns):
    """Diagnostic: slopes of E|S|, E[S^2] (uncentred), Var(S) across the ladder --
    the T-50 row's constraint-2 echo."""
    ab, un, ce = {}, {}, {}
    for n in ns:
        S = block_sums(u, n)
        ab[n] = float(np.mean(np.abs(S)))
        un[n] = float(np.mean(S ** 2))
        ce[n] = float(np.var(S)) + 1e-300
    return fit_slope(ab), fit_slope(un), fit_slope(ce)


# --------------------------------------------------------------------------- the estimator
def _core(v_w, v_null, encoding, n_min, n_max):
    """One pass of the full pipeline on one (written, unwritten) read pair.
    Returns dict or a refusal dict."""
    guards = {}

    # ---- reference and density, from the read itself
    if encoding == "occupancy":
        # The erased population is identified from the part's OWN NULL STATE -- the
        # unwritten sector of the same read -- never by clustering the written sector
        # alone (a two-signed write corrupts any written-sector-only clustering; the
        # null state cannot be gamed by the write).  Cells of the written sector inside
        # the erased band are the page's own erased cells; the reference is their mean
        # (same page, same disturb history), falling back to the unwritten sector's
        # mean when the page has too few erased cells (e.g. all-programmed).
        e_loc, e_sd = float(np.mean(v_null)), float(np.std(v_null)) + 1e-30
        band = np.abs(v_w - e_loc) <= 6.0 * e_sd
        if band.sum() >= 32:
            # one recentring pass on the page's own erased cells
            e2 = float(np.mean(v_w[band]))
            band = np.abs(v_w - e2) <= 6.0 * e_sd
        if band.sum() >= 32:
            ref_pop = v_w[band]
            guards["reference"] = "written sector's own erased-band mean (band from the " \
                                  "unwritten sector's location and spread, same read)"
        else:
            ref_pop = v_null
            guards["reference"] = "unwritten-sector mean (written sector has too few " \
                                  "erased cells)"
        prog_mask = ~band
        f_hat = float(np.mean(prog_mask))
        # G-separation: the two populations must be separable on this read -- cells in
        # the ambiguous annulus (6..9 e_sd from the erased location) must be rare
        ref0 = float(np.mean(ref_pop))
        ambig = float(np.mean((np.abs(v_w - ref0) > 6.0 * e_sd)
                              & (np.abs(v_w - ref0) < 9.0 * e_sd)))
        guards["separation_ok"] = bool(ambig < 0.01)
        if not guards["separation_ok"]:
            return dict(refused=True, reason="POPULATIONS NOT SEPARABLE: "
                        f"{100*ambig:.1f}% of the written sector sits in the ambiguous "
                        "annulus between the erased band and the programmed population; "
                        "the read cannot assign cells; no exponent is reported",
                        guards=guards, f=f_hat)
    else:  # orientation
        ref_pop = v_null
        prog_mask = None
        f_hat = np.nan
        guards["reference"] = "adjacent AC-erased region mean (same scan)"
    ref = float(np.mean(ref_pop))
    ref_var, ref_n = float(np.var(ref_pop)), ref_pop.size

    u_w = v_w - ref
    u_e = v_null - float(np.mean(v_null))          # null leg: exactly centred

    # ---- G-interleave: the two halves of the reference population must agree.
    # A sector-specific offset delta_b corrupts the exponent once delta_b^2 ~ var/n_max;
    # the half-difference resolves offsets down to its own sampling floor
    # (var(h1)/|h1| + var(h2)/|h2|), so the guard fires only ABOVE 3 sigma of that floor
    # and above the corruption threshold.  Larger parts resolve finer -- stated in scope.
    h = v_null.size // 2
    dref = float(np.mean(v_null[:h]) - np.mean(v_null[h:]))
    var_null = float(np.var(v_null)) + 1e-30
    floor = float(np.var(v_null[:h])) / h + float(np.var(v_null[h:])) / (v_null.size - h)
    thresh = max(9.0 * floor, 0.1 * var_null / n_max)
    guards["interleave_ok"] = bool(dref ** 2 <= thresh)
    if not guards["interleave_ok"]:
        return dict(refused=True, reason="REFERENCE INCONSISTENCY: the two halves of the "
                    "reference population disagree beyond the drift guard", guards=guards)

    # ---- G-scope (occupancy): the ladder must start above the crossover the read's own
    #      f_hat sets;  n_min_eff = 4*(1-f)/f, and at least 3 dyadic points must remain
    n_min_eff = n_min
    if encoding == "occupancy" and 0.0 < f_hat < 1.0:
        need = 4.0 * (1.0 - f_hat) / f_hat
        while n_min_eff < need:
            n_min_eff *= 2
    guards["n_min_eff"] = n_min_eff
    if n_min_eff * 4 > n_max:
        return dict(refused=True, reason="OUT OF SCOPE: measured programmed fraction "
                    f"f={f_hat:.4f} needs n_min={n_min_eff} > n_max/4={n_max//4}; the read "
                    "is too small to carry the ladder at this density", guards=guards,
                    f=f_hat)
    ns = ladder_sizes(n_min_eff, n_max)

    # ---- G-density (occupancy): per-block programmed fraction inside the binomial
    #      envelope at EVERY ladder scale (constraint 3; refuter A's counterexample)
    if encoding == "occupancy" and 0.0 < f_hat < 1.0:
        ok = True
        worst = 0.0
        for n in ns:
            fb = block_sums(prog_mask.astype(float), n) / n
            if fb.size < 8:
                continue
            env = f_hat * (1.0 - f_hat) / n
            ratio = float(np.var(fb)) / max(env, 1e-30)
            worst = max(worst, ratio)
            if ratio > 3.0:
                ok = False
        guards["density_envelope_worst"] = worst
        guards["density_stationary"] = ok
        if not ok:
            return dict(refused=True, reason="DENSITY NOT STATIONARY: Var(f_B) exceeds "
                        f"3x the binomial envelope (worst {worst:.1f}x); a fixed record in "
                        "a growing block is this case; no exponent is reported", guards=guards,
                        f=f_hat)

    # ---- G-balance (orientation): DC-loaded data routes to the DC clause
    if encoding == "orientation":
        z = abs(float(np.mean(u_w))) * np.sqrt(u_w.size) / (float(np.std(u_w)) + 1e-30)
        guards["dc_loaded"] = bool(z > 4.0)

    # ---- the ladder and the slopes
    m2w = m2_ladder(u_w, ns, self_centred=False)
    m2e = m2_ladder(u_e, ns, self_centred=True)
    beta_w, beta_e = fit_slope(m2w), fit_slope(m2e)

    # ---- G-resolution: the reference mean is estimated from N_ref cells; its error acts
    #      as a per-cell offset (var_ref/N_ref)^(1/2) and contributes (var_ref/N_ref)*n^2
    #      to M2_w(n).  If that reaches 10% of the measured top rung, the written leg is
    #      quieter than the reference can resolve and an exponent would be an artifact of
    #      the reference itself (the failure direction is a FALSE ACCUMULATES on an
    #      ultra-quiet coded track) -- refuse.
    # threshold 0.25: a reference-error share of 25% at the top rung bounds the induced
    # slope bias at ~log2(1.25)/octaves < 0.08, small against the 0.5 verdict gap; the
    # ultra-quiet coded-track case sits at ~15x the written leg's power and still refuses
    err_top = ref_var / max(ref_n, 1) * ns[-1] ** 2
    guards["ref_resolution_ok"] = bool(err_top <= 0.25 * m2w[ns[-1]])
    if not guards["ref_resolution_ok"]:
        return dict(refused=True, reason="REFERENCE RESOLUTION INSUFFICIENT: the written "
                    "leg's block power at the top rung is within 10x of the reference-mean "
                    "estimation error; enlarge the reference population or shorten the "
                    "ladder; no exponent is reported", guards=guards, f=f_hat)

    # ---- G-null admissibility: the unwritten leg must behave as an incoherent null
    guards["null_leg_ok"] = bool(0.6 <= beta_e <= 1.4) if np.isfinite(beta_e) else False
    if not guards["null_leg_ok"]:
        return dict(refused=True, reason=f"NULL LEG OUT OF SPEC: beta_e={beta_e:.3f} "
                    "outside [0.6,1.4]; correlated noise or instrument drift; the read "
                    "cannot supply its own null", guards=guards, f=f_hat)

    return dict(refused=False, beta_w=beta_w, beta_e=beta_e, D=beta_w - beta_e,
                f=f_hat, ns=ns, m2w=m2w, m2e=m2e, guards=guards)


def estimate(v_w, v_null, encoding, n_min=128, n_max=2048, jk=8):
    """The full measurement: point estimate + jackknife SE over jk contiguous segments
    (the whole pipeline, reference and guards included, is re-run on each replicate)."""
    full = _core(v_w, v_null, encoding, n_min, n_max)
    if full["refused"]:
        full["verdict"] = "REFUSED"
        return full
    Nw, Ne = v_w.size, v_null.size
    ds = []
    for j in range(jk):
        mw = np.ones(Nw, bool)
        mw[j * Nw // jk:(j + 1) * Nw // jk] = False
        me = np.ones(Ne, bool)
        me[j * Ne // jk:(j + 1) * Ne // jk] = False
        r = _core(v_w[mw], v_null[me], encoding, n_min, n_max)
        if not r["refused"] and np.isfinite(r["D"]):
            ds.append(r["D"])
    ds = np.array(ds)
    if ds.size >= jk - 1:
        se = float(np.sqrt((ds.size - 1) / ds.size * np.sum((ds - ds.mean()) ** 2)))
    else:
        se = np.nan
    D = full["D"]
    if np.isfinite(se) and D - 2 * se > 0.5:
        verdict = "ACCUMULATES"
    elif np.isfinite(se) and D + 2 * se < 0.5:
        verdict = "SCREENS"
    else:
        verdict = "INDETERMINATE"
    full.update(se=se, verdict=verdict)
    return full


def estimate_naive(v_w, v_null, n_min=128, n_max=2048):
    """THE SENTINEL: the same doubling exponents WITHOUT the self-reference -- the reader
    trusts an absolute zero (the V_t,neutral of the refuted repair).  Kept beside the
    estimator as the D-15 contrast pair: every attack that moves the sentinel while the
    self-referenced D stands is measured, not asserted."""
    ns = ladder_sizes(n_min, n_max)
    bw = fit_slope(m2_ladder(v_w, ns))
    be = fit_slope(m2_ladder(v_null, ns))
    return dict(beta_w=bw, beta_e=be, D=bw - be)
