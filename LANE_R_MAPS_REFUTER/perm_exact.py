"""High-precision lambda: split at crossovers and at zeros; mpmath tanh-sinh."""
import itertools, numpy as np, mpmath as mp
mp.mp.dps = 30

def lam_exact(pi):
    a0,a1,b0,b1 = [mp.mpf(float(x)) for x in pi]
    def A2(t): return a0**2 + a1**2 + 2*a0*a1*mp.cos(t)
    def B2(t): return b0**2 + b1**2 + 2*b0*b1*mp.cos(t)
    c0 = a0**2+a1**2-b0**2-b1**2
    c1 = 2*(a0*a1 - b0*b1)
    pts = [mp.mpf(0), mp.pi]
    if c1 != 0:
        x = -c0/c1
        if -1 < x < 1: pts.append(mp.acos(x))
    # zeros of either branch occur at t=pi when a0==a1 or b0==b1 (log singularity)
    pts = sorted(set(pts))
    def f(t):
        m = mp.sqrt(max(A2(t), B2(t)))
        return mp.log(m) if m > 0 else mp.mpf('-inf')
    tot = mp.mpf(0)
    for lo, hi in zip(pts[:-1], pts[1:]):
        tot += mp.quad(f, [lo, hi])
    return tot/mp.pi

def spread(pi):
    vals = [lam_exact(np.array(pi)[list(p)]) for p in itertools.permutations(range(4))]
    return max(vals)-min(vals), sorted({mp.nstr(v,14) for v in vals})

tests = [(0.1,0.2,0.3,0.4), (0.1,0.4,0.4,0.1), (0.4,0.05,0.3,0.25),
         (0.5,0.02,0.28,0.20), (0.6,0.01,0.2,0.19), (0.34,0.33,0.32,0.01),
         (0.55,0.30,0.10,0.05), (0.7,0.1,0.1,0.1), (0.25,0.25,0.25,0.25),
         (0.46,0.04,0.26,0.24), (0.9,0.04,0.03,0.03)]
for pi in tests:
    s, vals = spread(pi)
    print(f"pi={pi}  spread={mp.nstr(s,6):>14s}  distinct={len(vals)}  {vals if len(vals)<5 else vals[:4]}")
