import numpy as np
# is the 3.8e-11 in Block 9 a real deviation or floating-point u**k for large k?
rng = np.random.default_rng(2718281828)
for kmax in (10, 100, 10**3, 10**5):
    w1 = w2 = 0.0
    r = np.random.default_rng(2718281828)
    for _ in range(4000):
        pi = r.dirichlet(np.ones(4)); f, c = r.uniform(0, 2*np.pi, 2)
        u, v = np.exp(-1j*f), np.exp(1j*c); ks = r.integers(1, kmax+1, 30)
        Z  = pi[0]+pi[1]*u**ks+pi[2]*v**ks+pi[3]*(u*v)**ks
        pa = pi[[3,2,1,0]]
        Za = pa[0]+pa[1]*u**ks+pa[2]*v**ks+pa[3]*(u*v)**ks
        w1 = max(w1, float(np.abs(np.abs(Z)-np.abs(Za)).max()))
        pb = pi[[1,0,3,2]]; ub = 1/u
        Zb = pb[0]+pb[1]*ub**ks+pb[2]*v**ks+pb[3]*(ub*v)**ks
        w2 = max(w2, float(np.abs(np.abs(Z)-np.abs(Zb)).max()))
    print(f"k <= {kmax:>6d}:  (00 11)(10 01) dev = {w1:.2e}   (00 10)(01 11) dev = {w2:.2e}")
