"""V1 -- ADVERSARIAL RE-FIT OF THE LANE'S chi SCALING.

The lane reports  chi_total ~ N^a  with a = -0.4638 +- 0.0212 (distributed, nq=3) and
a = -0.2879 +- 0.0090 (shared, nq=3), and calls the CONTROL a = +1.0010 +- 0.0012 "LINEAR".

Three attacks, all run rather than argued:
  A1  COMPETING ASYMPTOTIC FORMS on the lane's own numbers: power law vs exponential vs
      inverse-log vs saturate-to-a-positive-constant vs stretched exponential.  Same data,
      same y-space, rms + AIC printed side by side.
  A2  IS THE EXPONENT A LAW OR A LOCAL SLOPE?  The data are deterministic, so the OLS
      sigma(a) is not an uncertainty -- it is the misfit of a wrong model divided by sqrt(n).
      The honest uncertainty is the DRIFT of the local slope across the range.  Printed.
  A3  IS THE EXPONENT A PROPERTY OF THE SYSTEM OR OF THE 25-POINT TIME GRID?  D-17 applied to
      the lane's own venue scale: re-run the identical engine with other time grids, other
      lam, other beta, and re-fit.  If the exponent moves far outside +-0.02, the reported
      exponent is a property of the grid, not of the record count.

The engine is COPIED VERBATIM from the lane's s3_chi_scaling.py (functions chi_distributed /
chi_shared) so that this is a re-fit of THEIR numbers, not of a different computation.  A
reproduction check against their s3_chi_scaling.json is printed first; if it fails, nothing
is concluded.
"""
import sys, json, math
import numpy as np

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_A_ZERO_OR_GROWING"
OUT = []
def p(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

I2 = np.eye(2); Xm = np.array([[0, 1], [1, 0]], complex); Zm = np.array([[1, 0], [0, -1]], complex)
ENERGY_POOL = (1.0, 1.4, 0.7, 1.2, 0.9, 1.6, 0.8, 1.1)
def energies(nq): return tuple(ENERGY_POOL[j % len(ENERGY_POOL)] for j in range(nq))

def vN(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

def uprop(H, t):
    w, V = np.linalg.eigh(H)
    return (V * np.exp(-1j * w * t)) @ V.conj().T

def binom_pmf(m):
    return {2 * a - m: math.comb(m, a) / 2 ** m for a in range(m + 1)}

def chi_distributed(k, nq, t, lam=0.8, beta=2.0):
    E = energies(nq); per = [0.0] * k
    for j in range(nq):
        idx = [i for i in range(k) if i % nq == j]
        m = len(idx)
        if m == 0: continue
        Hb = E[j] * Zm
        w = np.exp(-beta * np.array([E[j], -E[j]]))
        rth = np.diag(w / w.sum()).astype(complex)
        def rho(c):
            U = uprop(Hb + lam * c * Xm, t)
            return U @ rth @ U.conj().T
        pm = binom_pmf(m)
        rbar = sum(pr * rho(c) for c, pr in pm.items())
        pm1 = binom_pmf(m - 1)
        cond = {sgn: sum(pr * rho(c + sgn) for c, pr in pm1.items()) for sgn in (+1, -1)}
        chi = vN(rbar) - 0.5 * (vN(cond[+1]) + vN(cond[-1]))
        for i in idx: per[i] = max(chi, 0.0)
    return per, float(sum(per))

def chi_shared(k, nq, t, lam=0.8, beta=2.0):
    E = energies(nq)
    def bop(jj, P):
        M = np.array([[1]], complex)
        for q in range(nq): M = np.kron(M, P if q == jj else I2)
        return M
    HB = sum(E[j] * bop(j, Zm) for j in range(nq))
    probe = sum(bop(j, Xm) for j in range(nq))
    ww, VV = np.linalg.eigh(HB)
    pth = np.exp(-beta * ww); pth /= pth.sum()
    rth = (VV * pth) @ VV.conj().T
    cache = {}
    def rho(C):
        if C not in cache:
            U = uprop(HB + lam * C * probe, t)
            cache[C] = U @ rth @ U.conj().T
        return cache[C]
    pm = binom_pmf(k)
    rbar = sum(pr * rho(C) for C, pr in pm.items())
    pm1 = binom_pmf(k - 1)
    cond = {sgn: sum(pr * rho(C + sgn) for C, pr in pm1.items()) for sgn in (+1, -1)}
    chi = max(vN(rbar) - 0.5 * (vN(cond[+1]) + vN(cond[-1])), 0.0)
    return [chi] * k, float(k * chi)

def series(engine, nq, KS, times, lam=0.8, beta=2.0):
    out = []
    for k in KS:
        tot = 0.0
        for t in times:
            tot += engine(k, nq, float(t), lam, beta)[1]
        out.append(tot / len(times))
    return out

KS = [2, 4, 6, 8, 10, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]
TIMES = np.linspace(1.0, 13.0, 25)

p("=" * 118)
p("V1  ADVERSARIAL RE-FIT OF THE chi SCALING CLAIMS  (engine copied verbatim from the lane's s3)")
p("=" * 118)
p("")

# ---------- reproduction check against the lane's own json
lane = json.load(open(LANE + "/s3_chi_scaling.json"))
mine = {}
for mode, eng in (("distributed", chi_distributed), ("shared", chi_shared)):
    mine[mode] = series(eng, 3, KS, TIMES)
dev = 0.0
for mode in ("distributed", "shared"):
    for i, k in enumerate(KS):
        dev = max(dev, abs(mine[mode][i] - lane["total"]["%s|3|%d" % (mode, k)]))
ctrl = [np.mean([chi_distributed(k, k, float(t))[1] for t in TIMES]) for k in KS]
for i, k in enumerate(KS):
    dev = max(dev, abs(ctrl[i] - lane["control_bath_scaled"][str(k)]))
p("REPRODUCTION CHECK vs the lane's s3_chi_scaling.json, max abs deviation: %.3e  (must be < 1e-9)" % dev)
if dev >= 1e-9:
    p("REPRODUCTION FAILED -- CONCLUDING NOTHING."); sys.exit(1)
p("reproduction ok: what follows re-fits THEIR numbers.")
p("")

# ---------- A1 competing forms
def fit_forms(ks, ys):
    ks = np.array(ks, float); ys = np.array(ys, float)
    lk = np.log(ks); ly = np.log(ys)
    out = {}
    def rec(name, pred, npar):
        r = ys - pred
        rms = float(np.sqrt(np.mean(r ** 2)))
        n = len(ys)
        aic = n * math.log(max(rms ** 2, 1e-300)) + 2 * npar
        out[name] = (rms, aic, float(np.max(np.abs(r))))
    # power law  y = A k^a   (fit in log space, the lane's own fit)
    A = np.vstack([lk, np.ones_like(lk)]).T
    c = np.linalg.lstsq(A, ly, rcond=None)[0]
    rec("power law  A*k^a", np.exp(A @ c), 2); out["power_a"] = float(c[0])
    # exponential  y = A exp(-b k)
    A2 = np.vstack([ks, np.ones_like(ks)]).T
    c2 = np.linalg.lstsq(A2, ly, rcond=None)[0]
    rec("exponential A*exp(-b k)", np.exp(A2 @ c2), 2)
    # inverse-log  y = A/(b + ln k)  -> 1/y = (b + ln k)/A linear in ln k
    A3 = np.vstack([lk, np.ones_like(lk)]).T
    c3 = np.linalg.lstsq(A3, 1.0 / ys, rcond=None)[0]
    rec("inverse-log A/(b+ln k)", 1.0 / (A3 @ c3), 2)
    # saturating to a POSITIVE constant  y = c0 + A k^a  (grid over c0 and a)
    best = None
    for c0 in np.linspace(0.0, float(ys.min()) * 0.999, 60):
        for a in np.linspace(-2.5, 0.0, 120):
            b = ks ** a
            A4 = np.vstack([b]).T
            amp = np.linalg.lstsq(A4, ys - c0, rcond=None)[0]
            pred = c0 + A4 @ amp
            r = float(np.sqrt(np.mean((ys - pred) ** 2)))
            if best is None or r < best[0]: best = (r, c0, a, pred)
    rec("saturate c0 + A*k^a", best[3], 3); out["sat_c0"] = best[1]; out["sat_a"] = best[2]
    # stretched exponential y = A exp(-b k^q)
    best2 = None
    for q in np.linspace(0.05, 1.0, 40):
        A5 = np.vstack([ks ** q, np.ones_like(ks)]).T
        c5 = np.linalg.lstsq(A5, ly, rcond=None)[0]
        pred = np.exp(A5 @ c5)
        r = float(np.sqrt(np.mean((ys - pred) ** 2)))
        if best2 is None or r < best2[0]: best2 = (r, q, pred)
    rec("stretched A*exp(-b k^q)", best2[2], 3); out["stretch_q"] = best2[1]
    return out

p("A1  COMPETING ASYMPTOTIC FORMS, fitted to the SAME numbers, residuals in y-space (chi bits).")
p("    A model is only 'supported by the data' if it beats the alternatives.  Lower rms / AIC is better.")
p("")
for label, ys in (("distributed fixed bath nq=3", mine["distributed"]),
                  ("shared      fixed bath nq=3", mine["shared"]),
                  ("CONTROL bath scaled nq=k   ", ctrl)):
    f = fit_forms(KS, ys)
    p("  %s   (lane calls this: %s)" % (label, "decaying power law" if "CONTROL" not in label else "LINEAR"))
    p("     %-26s %12s %12s %12s" % ("form", "rms(y)", "AIC", "max|resid|"))
    for nm in ("power law  A*k^a", "exponential A*exp(-b k)", "inverse-log A/(b+ln k)",
               "saturate c0 + A*k^a", "stretched A*exp(-b k^q)"):
        rms, aic, mx = f[nm]
        p("     %-26s %12.5f %12.2f %12.5f" % (nm, rms, aic, mx))
    p("     power-law exponent a = %+.4f ;  best saturating fit has c0 = %.4f (a positive floor), a = %+.3f ;"
      % (f["power_a"], f["sat_c0"], f["sat_a"]))
    p("     best stretched exponent q = %.3f" % f["stretch_q"])
    p("")

# ---------- A2 local slopes
p("A2  IS THE EXPONENT A LAW, OR A LOCAL SLOPE THAT DRIFTS?")
p("    local slope a_i = ln(y_{i+1}/y_i) / ln(k_{i+1}/k_i).  A true power law has a CONSTANT column.")
p("")
p("    k range        distributed nq=3     shared nq=3      CONTROL nq=k")
for i in range(len(KS) - 1):
    row = []
    for ys in (mine["distributed"], mine["shared"], ctrl):
        row.append(math.log(ys[i + 1] / ys[i]) / math.log(KS[i + 1] / KS[i]))
    p("    %4d -> %4d   %16.4f %16.4f %16.4f" % (KS[i], KS[i + 1], row[0], row[1], row[2]))
for label, ys in (("distributed nq=3", mine["distributed"]), ("shared nq=3", mine["shared"]),
                  ("CONTROL nq=k", ctrl)):
    sl = [math.log(ys[i + 1] / ys[i]) / math.log(KS[i + 1] / KS[i]) for i in range(len(KS) - 1)]
    tail = sl[-5:]
    p("    %-18s local slopes span [%+.4f, %+.4f]  (spread %.4f);  last five mean %+.4f, spread %.4f"
      % (label, min(sl), max(sl), max(sl) - min(sl), float(np.mean(tail)), max(tail) - min(tail)))
p("")

# ---------- A2b window-dependence of the fitted exponent
def loglog_a(ks, ys):
    x = np.log(np.array(ks, float)); y = np.log(np.array(ys, float))
    A = np.vstack([x, np.ones_like(x)]).T
    c, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ c; resid = y - pred
    dof = max(len(x) - 2, 1); s2 = float((resid ** 2).sum() / dof)
    cov = s2 * np.linalg.inv(A.T @ A)
    return float(c[0]), float(np.sqrt(cov[0, 0])), float(np.sqrt(s2))

p("A2b FITTING WINDOW DEPENDENCE (the lane fitted the whole range, k = 2..256).")
p("    window            distributed a +- s     shared a +- s        CONTROL a +- s")
for lo in (2, 8, 16, 32, 64):
    sel = [i for i, k in enumerate(KS) if k >= lo]
    if len(sel) < 3: continue
    ks = [KS[i] for i in sel]
    r = []
    for ys in (mine["distributed"], mine["shared"], ctrl):
        r.append(loglog_a(ks, [ys[i] for i in sel]))
    p("    k >= %-4d (%2d pts)  %+.4f +- %.4f    %+.4f +- %.4f    %+.4f +- %.4f"
      % (lo, len(ks), r[0][0], r[0][1], r[1][0], r[1][1], r[2][0], r[2][1]))
p("")

# ---------- A3 venue-scale sweep on the fit
p("A3  D-17 ON THE FIT ITSELF: does the exponent belong to the system or to the 25-point time grid?")
p("    Same engine, same k list; only the time grid / lam / beta change.  Fitted over the full range.")
p("")
p("    setting                                  distributed a +- s   shared a +- s   CONTROL a +- s")
settings = [
    ("lane's own: 25 t in [1,13], lam .8 b 2", np.linspace(1, 13, 25), 0.8, 2.0),
    ("101 t in [1,13]                       ", np.linspace(1, 13, 101), 0.8, 2.0),
    ("25 t in [1.5,13.5]                     ", np.linspace(1.5, 13.5, 25), 0.8, 2.0),
    ("25 t in [0.5,30]                       ", np.linspace(0.5, 30, 25), 0.8, 2.0),
    ("25 t in [20,40]                        ", np.linspace(20, 40, 25), 0.8, 2.0),
    ("25 t in [1,13], lam = 0.4              ", np.linspace(1, 13, 25), 0.4, 2.0),
    ("25 t in [1,13], lam = 1.6              ", np.linspace(1, 13, 25), 1.6, 2.0),
    ("25 t in [1,13], beta = 0.5             ", np.linspace(1, 13, 25), 0.8, 0.5),
]
sweep = {}
for name, tt, lam, beta in settings:
    d = series(chi_distributed, 3, KS, tt, lam, beta)
    s = series(chi_shared, 3, KS, tt, lam, beta)
    c = [np.mean([chi_distributed(k, k, float(t), lam, beta)[1] for t in tt]) for k in KS]
    rd, rs, rc = loglog_a(KS, d), loglog_a(KS, s), loglog_a(KS, c)
    sweep[name.strip()] = dict(dist=rd, shared=rs, ctrl=rc,
                               dist_first=d[0], dist_last=d[-1], ctrl_last=c[-1])
    p("    %-40s %+.4f +- %.4f  %+.4f +- %.4f  %+.4f +- %.4f"
      % (name, rd[0], rd[1], rs[0], rs[1], rc[0], rc[1]))
p("")
ds = [sweep[n.strip()]["dist"][0] for n, _, _, _ in settings]
ss = [sweep[n.strip()]["shared"][0] for n, _, _, _ in settings]
cs = [sweep[n.strip()]["ctrl"][0] for n, _, _, _ in settings]
p("    distributed exponent across settings: min %+.4f  max %+.4f  SPREAD %.4f   (lane's quoted sigma 0.0212)"
  % (min(ds), max(ds), max(ds) - min(ds)))
p("    shared      exponent across settings: min %+.4f  max %+.4f  SPREAD %.4f   (lane's quoted sigma 0.0090)"
  % (min(ss), max(ss), max(ss) - min(ss)))
p("    CONTROL     exponent across settings: min %+.4f  max %+.4f  SPREAD %.4f   (lane's quoted sigma 0.0012)"
  % (min(cs), max(cs), max(cs) - min(cs)))
p("")

# ---------- the exact bound, re-checked
p("EXACT-BOUND RE-CHECK (this is the part of the lane's (S) verdict that does NOT rest on a fit):")
worst = 0.0; worstlab = ""
for mode, eng in (("distributed", chi_distributed), ("shared", chi_shared)):
    for nq in (1, 2, 3, 4):
        for k in KS:
            v = np.mean([eng(k, nq, float(t))[1] for t in TIMES])
            if v / nq > worst: worst, worstlab = v / nq, "%s nq=%d k=%d chi=%.4f" % (mode, nq, k, v)
p("    max over all (mode, nq, k) of chi_total / nq = %.4f   (bound is 1.0);  attained at %s" % (worst, worstlab))
p("    -> total chi never exceeds log2 dim(bath) = nq.  Standard Holevo bound; the lane's chain-rule")
p("       argument is correct and is an ORDINARY result, not a new one.")

with open(LANE + "/VERIFY/v1_refit_chi.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
json.dump(dict(sweep={a: dict(dist=b["dist"], shared=b["shared"], ctrl=b["ctrl"]) for a, b in sweep.items()}),
          open(LANE + "/VERIFY/v1_refit_chi.json", "w"), indent=1)
