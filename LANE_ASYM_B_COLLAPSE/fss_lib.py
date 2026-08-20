"""fss_lib -- the finite-size-scaling engine of S4, factored out so S6 can import it
without re-running S4.  Same bodies; the two can never drift apart."""

import sys, numpy as np
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE"

OUT = []
def say(s=""):
    print(s); OUT.append(s)

# ------------------------------------------------------------------ fitting primitives
def _wls(X, y, sig):
    W = 1.0/np.asarray(sig)
    A = X*W[:, None]; b = y*W
    beta, *_ = np.linalg.lstsq(A, b, rcond=None)
    r = y - X@beta
    chi2 = float(np.sum((r/sig)**2))
    cov = np.linalg.pinv(A.T@A)
    return beta, chi2, r, cov

def fit_form(name, N, Q, sig):
    N = np.asarray(N, float); Q = np.asarray(Q, float); sig = np.asarray(sig, float)
    if name == "LIN":  X = np.c_[N, np.ones_like(N)]
    elif name == "LOG": X = np.c_[np.log(N), np.ones_like(N)]
    elif name == "SAT1": X = np.c_[np.ones_like(N), -1.0/N]
    elif name == "SAT2":
        best = None
        for xi in np.exp(np.linspace(np.log(0.2), np.log(200.0), 400)):
            Xx = np.c_[np.ones_like(N), -np.exp(-N/xi)]
            beta, chi2, r, cov = _wls(Xx, Q, sig)
            if best is None or chi2 < best[1]: best = (beta, chi2, r, cov, xi)
        beta, chi2, r, cov, xi = best
        return dict(name=name, p=3, beta=list(beta)+[xi], chi2=chi2, resid=r, cov=cov)
    else: raise ValueError(name)
    beta, chi2, r, cov = _wls(X, Q, sig)
    return dict(name=name, p=2, beta=list(beta), chi2=chi2, resid=r, cov=cov)

def aicc(chi2, m, p):
    if m - p - 1 <= 0: return np.inf
    return chi2 + 2*p + 2*p*(p+1)/(m-p-1)

def classify(N, Q, sig, label):
    """Returns a dict with everything the master table needs."""
    N = np.asarray(N, float); Q = np.asarray(Q, float); sig = np.asarray(sig, float)
    m = len(N); scale = max(abs(Q).max(), 1e-300)
    res = dict(label=label, m=m, Nmin=N.min(), Nmax=N.max(), Qmin=Q.min(), Qmax=Q.max())
    # --- exact categories first
    if abs(Q).max() < 1e-12:
        res.update(category="IDENTICALLY ZERO", best="-", dAICc=np.inf, alpha=np.nan,
                   cv=np.nan, Q0=0.0, Q0e=0.0, dbl=np.nan, rms1N=0.0, extensive=False,
                   expo=np.nan, expo_e=np.nan, note="exact zero at every N tested")
        return res
    if (Q.max()-Q.min())/scale < 1e-12:
        res.update(category="CONSTANT", best="-", dAICc=np.inf, alpha=0.0, cv=0.0,
                   Q0=float(Q.mean()), Q0e=0.0, dbl=1.0, rms1N=0.0, extensive=False,
                   expo=0.0, expo_e=0.0, note="exact constant at every N tested")
        return res
    # --- model selection
    fits = {nm: fit_form(nm, N, Q, sig) for nm in ("LIN", "LOG", "SAT1", "SAT2")}
    sc = sorted(((aicc(f["chi2"], m, f["p"]), nm) for nm, f in fits.items()))
    best, second = sc[0], sc[1]
    d = second[0]-best[0]
    # --- collapse
    als = np.linspace(-1.5, 2.5, 4001)
    cvs = [np.std(Q/N**a, ddof=0)/abs(np.mean(Q/N**a)) for a in als]
    ia = int(np.argmin(cvs)); alpha, cv = float(als[ia]), float(cvs[ia])
    # --- power-law exponent with uncertainty (only where Q > 0 everywhere)
    if Q.min() > 0:
        X = np.c_[np.log(N), np.ones_like(N)]
        rel = sig/Q
        beta, chi2, r, cov = _wls(X, np.log(Q), np.maximum(rel, 1e-12))
        dof = max(m-2, 1)
        se = float(np.sqrt(cov[0, 0])*max(1.0, np.sqrt(chi2/dof)))
        expo, expo_e = float(beta[0]), se
    else:
        expo, expo_e = np.nan, np.nan
    # --- 1/N extrapolation on the largest-N half
    h = max(4, m//2); Nh, Qh, sh = N[-h:], Q[-h:], sig[-h:]
    X = np.c_[np.ones_like(Nh), 1.0/Nh, 1.0/Nh**2]
    beta, chi2, r, cov = _wls(X, Qh, sh)
    dof = max(len(Nh)-3, 1)
    Q0 = float(beta[0]); Q0e = float(np.sqrt(cov[0, 0])*max(1.0, np.sqrt(chi2/dof)))
    # --- doubling
    dbls = []
    lut = {int(round(x)): q for x, q in zip(N, Q)}
    for x in N:
        a, b = int(round(x)), int(round(2*x))
        if a in lut and b in lut and lut[a] != 0: dbls.append(lut[b]/lut[a])
    dbl = float(dbls[-1]) if dbls else np.nan
    # --- 1/N fit quality: for a DIVERGENT quantity the intercept is meaningless and the
    #     residual of the 1/N fit is what says so.
    rms1N = float(np.sqrt(np.mean(r**2)))
    # --- gravity's requirement (a), tested directly and not through a fit
    ext = (d >= 4.0 and best[1] == "LIN" and dbl == dbl and 1.6 <= dbl <= 2.6
           and expo == expo and abs(expo-1.0) < 0.25)
    # --- verdict
    decaying = Q[-1] < Q[0]
    if decaying and d >= 4.0:
        cat = "DECAYING"
    elif decaying:
        cat = "DECAYING (form undetermined)"
    elif d < 4.0:
        cat = "CANNOT DISTINGUISH"
    elif best[1] == "LIN":
        cat = "GROWING (linear)"
    elif best[1] == "LOG":
        cat = "GROWING (logarithmic, SUB-extensive)"
    else:
        cat = "SATURATING"
    if best[1] == "LIN" and expo == expo and expo > 1.5 and d >= 4:
        cat = "GROWING (super-linear)"
    res.update(category=cat, best=best[1], dAICc=float(d), alpha=alpha, cv=cv,
               Q0=Q0, Q0e=Q0e, dbl=dbl, expo=expo, expo_e=expo_e, rms1N=rms1N, extensive=bool(ext),
               fits=fits, rank=[(nm, float(a)) for a, nm in sc],
               rss=float(np.sum(fits[best[1]]["resid"]**2)),
               note="")
    return res

def ascii_plot(N, Q, title, width=64, height=12):
    N = np.asarray(N, float); Q = np.asarray(Q, float)
    lines = ["    " + title]
    lo, hi = Q.min(), Q.max()
    if hi - lo < 1e-15: hi = lo + 1.0
    grid = [[" "]*width for _ in range(height)]
    for x, y in zip(N, Q):
        c = int(round((x-N.min())/(N.max()-N.min()+1e-30)*(width-1)))
        r = height-1-int(round((y-lo)/(hi-lo)*(height-1)))
        grid[r][c] = "o"
    for i, row in enumerate(grid):
        v = hi - (hi-lo)*i/(height-1)
        lines.append("   %10.4f |%s" % (v, "".join(row)))
    lines.append("              +" + "-"*width)
    lines.append("              %-*s%s" % (width-8, "N=%g" % N.min(), "N=%g" % N.max()))
    return lines

