"""REFUTER A (computation) -- INDEPENDENT REIMPLEMENTATION of the T-50 DESIGN TWO
observable, written from the design's PREDICTION/FALSIFIER text before consulting
d2_observable.py's internals for anything but the declared constants (ladder 128..2048,
jackknife 8, guard thresholds as stated in the sealed text).  Used to rebuild the
design's model-side numbers with code the designer never saw.

The estimator, per the sealed text:
  D = beta_w - beta_e over a dyadic ladder of block sizes; M2(n) = mean over two-phase
  contiguous blocks of (sum u_i)^2; occupancy reference = mean of the written sector's
  own erased-band cells (band located from the unwritten sector's mean/sd, 6 sd, one
  recentring pass, fallback = unwritten mean); orientation reference = mean of the
  adjacent AC-erased region; null leg self-centred with 1/(1-n/N).
Guards, per the sealed text: separation annulus (6..9 sd) < 1%; interleave (halves of
the reference population); scope n_min >= 4(1-f)/f with >= 3 dyadic points; density
envelope Var(f_B) <= 3 f(1-f)/n at every scale; DC balance z > 4 routes; resolution
(ref_var/N_ref) n_max^2 <= 0.25 M2_w(n_max); null slope in [0.6, 1.4].
Verdicts: ACCUMULATES (D-2SE>0.5), SCREENS (D+2SE<0.5), else INDETERMINATE; REFUSED.
"""
import numpy as np


def _ladder(n_min, n_max):
    ns = []
    n = n_min
    while n <= n_max:
        ns.append(n)
        n *= 2
    return ns


def _bsums(u, n):
    k = u.size // n
    a = u[:k * n].reshape(k, n).sum(1)
    v = u[n // 2:]
    k2 = v.size // n
    b = v[:k2 * n].reshape(k2, n).sum(1)
    return np.r_[a, b]


def _m2(u, ns, self_centred=False):
    N = u.size
    out = {}
    for n in ns:
        s = _bsums(u, n)
        m = float(np.mean(s * s))
        if self_centred:
            m /= max(1.0 - n / N, 1e-9)
        out[n] = m
    return out


def _slope(m2):
    x = np.log2(np.array(sorted(m2), float))
    y = np.log2(np.array([m2[int(2 ** xi + 0.5)] for xi in x], float))
    x = x - x.mean()
    return float(np.dot(x, y - y.mean()) / np.dot(x, x))


def core(v_w, v_null, encoding, n_min=128, n_max=2048):
    g = {}
    if encoding == "occupancy":
        e_loc, e_sd = float(np.mean(v_null)), float(np.std(v_null)) + 1e-30
        band = np.abs(v_w - e_loc) <= 6 * e_sd
        if band.sum() >= 32:
            band = np.abs(v_w - float(np.mean(v_w[band]))) <= 6 * e_sd
        if band.sum() >= 32:
            ref_pop = v_w[band]
        else:
            ref_pop = v_null
            band = np.zeros(v_w.size, bool)
        f = float(np.mean(~band))
        r0 = float(np.mean(ref_pop))
        amb = float(np.mean((np.abs(v_w - r0) > 6 * e_sd) & (np.abs(v_w - r0) < 9 * e_sd)))
        g["ambig"] = amb
        if amb >= 0.01:
            return dict(refused=True, why="separation", g=g, f=f)
        prog = ~band
    else:
        ref_pop, prog, f = v_null, None, np.nan
    ref = float(np.mean(ref_pop))
    u_w = v_w - ref
    u_e = v_null - float(np.mean(v_null))

    h = v_null.size // 2
    dref = float(np.mean(v_null[:h]) - np.mean(v_null[h:]))
    floor = float(np.var(v_null[:h])) / h + float(np.var(v_null[h:])) / (v_null.size - h)
    thr = max(9 * floor, 0.1 * (float(np.var(v_null)) + 1e-30) / n_max)
    if dref * dref > thr:
        return dict(refused=True, why="interleave", g=g, f=f)

    nme = n_min
    if encoding == "occupancy" and 0 < f < 1:
        need = 4 * (1 - f) / f
        while nme < need:
            nme *= 2
    if nme * 4 > n_max:
        return dict(refused=True, why="scope", g=g, f=f)
    ns = _ladder(nme, n_max)

    if encoding == "occupancy" and 0 < f < 1:
        for n in ns:
            k = prog.size // n
            fb = prog[:k * n].reshape(k, n).mean(1)
            v2 = prog[n // 2:]
            k2 = v2.size // n
            fb = np.r_[fb, v2[:k2 * n].reshape(k2, n).mean(1)] if k2 else fb
            if fb.size < 8:
                continue
            if float(np.var(fb)) > 3 * f * (1 - f) / n:
                return dict(refused=True, why="density", g=g, f=f)

    if encoding == "orientation":
        z = abs(float(np.mean(u_w))) * np.sqrt(u_w.size) / (float(np.std(u_w)) + 1e-30)
        g["dc_loaded"] = bool(z > 4.0)
        g["z"] = z

    m2w = _m2(u_w, ns)
    m2e = _m2(u_e, ns, self_centred=True)
    bw, be = _slope(m2w), _slope(m2e)

    err_top = (float(np.var(ref_pop)) / max(ref_pop.size, 1)) * ns[-1] ** 2
    if err_top > 0.25 * m2w[ns[-1]]:
        return dict(refused=True, why="resolution", g=g, f=f)
    if not (0.6 <= be <= 1.4):
        return dict(refused=True, why="null_slope", g=g, f=f, beta_e=be)
    return dict(refused=False, beta_w=bw, beta_e=be, D=bw - be, f=f, g=g)


def run(v_w, v_null, encoding, n_min=128, n_max=2048, jk=8):
    full = core(v_w, v_null, encoding, n_min, n_max)
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
        r = core(v_w[mw], v_null[me], encoding, n_min, n_max)
        if not r["refused"] and np.isfinite(r["D"]):
            ds.append(r["D"])
    ds = np.array(ds)
    se = float(np.sqrt((ds.size - 1) / ds.size * np.sum((ds - ds.mean()) ** 2))) \
        if ds.size >= jk - 1 else np.nan
    D = full["D"]
    if np.isfinite(se) and D - 2 * se > 0.5:
        v = "ACCUMULATES"
    elif np.isfinite(se) and D + 2 * se < 0.5:
        v = "SCREENS"
    else:
        v = "INDETERMINATE"
    full.update(se=se, verdict=v)
    return full


# ---- my own read models (independent of d2_observable's), same declared physics
NE, DELTA = 100.0, 5


def occ(N, f, rng, mu=0.0, off_w=0.0, off_e=0.0, gain_w=1.0, gain_e=1.0,
        two_signed=False, prog_sd=2.0):
    data = (rng.random(N) < f).astype(int)
    resid = rng.integers(-DELTA, DELTA + 1, N).astype(float) + mu
    prog = NE + rng.normal(0, prog_sd, N)
    if two_signed:
        prog *= rng.choice([-1.0, 1.0], N)
    v_w = gain_w * (np.where(data == 1, prog, resid) + off_w)
    v_e = gain_e * (rng.integers(-DELTA, DELTA + 1, N).astype(float) + mu + off_e)
    return v_w, v_e


def ori(N, rng, kind="random", off_w=0.0, off_e=0.0, m=1.0, sd=0.1,
        drift_w=0.0, drift_e=0.0, cluster=1, tilt=0.0):
    if kind == "random":
        if cluster > 1:
            nb = N // cluster + 1
            s = np.repeat(rng.integers(0, 2, nb) * 2 - 1, cluster)[:N].astype(float)
        else:
            s = (rng.integers(0, 2, N) * 2 - 1).astype(float)
    elif kind == "dc":
        s = np.ones(N)
    elif kind == "dcfree":
        s = np.tile([1.0, -1.0], N // 2 + 1)[:N]
    if tilt > 0:
        s = s * np.cos(rng.normal(0, np.deg2rad(tilt), N))
    v_w = m * s + rng.normal(0, sd, N) + off_w
    if drift_w:
        v_w = v_w + drift_w * (np.linspace(0, 1, N) - 0.5)   # mean-centred ramp
    se_ = (rng.integers(0, 2, N) * 2 - 1).astype(float)
    if tilt > 0:
        se_ = se_ * np.cos(rng.normal(0, np.deg2rad(tilt), N))
    v_e = m * se_ + rng.normal(0, sd, N) + off_e
    if drift_e:
        v_e = v_e + drift_e * (np.linspace(0, 1, N) - 0.5)
    return v_w, v_e
