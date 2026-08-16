import itertools, numpy as np
M = 400001
t = 2*np.pi*(np.arange(M)+0.5)/M
X = np.exp(1j*t)
def lam(pi):
    A = np.abs(pi[0]+pi[1]*X); B = np.abs(pi[2]+pi[3]*X)
    return float(np.mean(np.log(np.maximum(A,B))))

rng = np.random.default_rng(97531)
worst = 0.0; argw = None
for _ in range(300):
    pi = rng.dirichlet(np.ones(4)*0.4)
    vals = [lam(pi[list(p)]) for p in itertools.permutations(range(4))]
    sp = max(vals)-min(vals)
    if sp > worst: worst, argw = sp, pi.copy()
print("seed 97531, 300 random pi (Dirichlet a=0.4), all 24 perms each")
print("worst spread =", worst, "at pi =", argw)
for pi in (np.array([0.30,0.30,0.39,0.01]), np.array([0.25,0.25,0.49,0.01]),
           np.array([0.20,0.20,0.59,0.01]), np.array([0.1,0.4,0.4,0.1]),
           np.array([0.05,0.45,0.45,0.05]), np.array([0.02,0.48,0.30,0.20])):
    vals = sorted({round(lam(pi[list(p)]),10) for p in itertools.permutations(range(4))})
    print(f"pi={pi} -> {len(vals)} distinct: {vals}")
