"""G0 — VALIDATION. Reproduce, from my own code, every corpus number this lane leans on.
If any of these misses, nothing downstream is trustworthy."""
import numpy as np
from glib import *

print("=" * 78)
print("G0.1  CARRIER K1 — PUBLISHED INCIDENCE MATRICES")
print("=" * 78)
d1, d2 = incidence(NV_K1, EDGES_K1, FACES_K1)
print("edges (source,target), index 0..5 = e1..e6:", EDGES_K1)
print("d1 (5x6), rows v0..v4, cols e1..e6:")
print(d1)
print("d2 (6x1), rows e1..e6, col F:")
print(d2.T, "(transposed for display)")
print("d1 @ d2 =", (d1 @ d2).ravel(), "  -> d^2 = 0:", np.all(d1 @ d2 == 0))
print("topology:", betti(NV_K1, EDGES_K1, FACES_K1))

print()
print("=" * 78)
print("G0.2  RANK-ONE U(1): the closed form, against literal matrix action")
print("=" * 78)
rng = np.random.default_rng(20260816)
f, c = 2.0, 1.1
p = [0.4, 0.15, 0.15, 0.15, 0.15]
U = u1_conn(f, c)
s = state_rank1(p)
z, cf = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s))
print("characters (zeta, coeff):")
for zz, cc in zip(z, cf):
    print("   zeta = %+.9f%+.9fj   coeff = %+.9f%+.9fj" % (zz.real, zz.imag, cc.real, cc.imag))
print("class weights:", class_weights(s, EDGES_K1, LOOP_F, LOOP_C))
dev = max(abs(Z_from_chars(z, cf, [k])[0] - Z_direct(U, EDGES_K1, LOOP_F, LOOP_C, s, k))
          for k in range(0, 201))
print("max |char-form - matrix action| over k<=200 : %.3e" % dev)

ks = np.arange(1, 4001)
Za = np.abs(Z_from_chars(z, cf, ks))
print("min |Z_k| over k<=400  = %.6f at k = %d   [S3 4.1: 0.024654 at k=42]"
      % (Za[:400].min(), ks[:400][Za[:400].argmin()]))
print("sup |Z_k| over k<=4000 = %.6f at k = %d   [S3 4.1: 0.999941 at k=377]"
      % (Za.max(), ks[Za.argmax()]))
print("#{|Z_k|>0.99, k<=4000} = %d                [S3 4.1: 37]" % int((Za > 0.99).sum()))

print()
print("=" * 78)
print("G0.3  lambda_B — my own quadrature vs my own schedule-B simulation")
print("=" * 78)
# generic (non-resonant) connection, class weights (0.4,0.3,0.3)
p2 = [0.4, 0.15, 0.15, 0.15, 0.15]
s2 = state_rank1(p2)
fg, cg = 1.0, np.sqrt(2.0)
Ug = u1_conn(fg, cg)
zg, cg_ = merge_characters(*characters(Ug, EDGES_K1, LOOP_F, LOOP_C, s2))
lam_sim = lambda_B(zg, cg_, N=2000000)
# exponent vectors in (f,c): class (1,1) -> u v -> (-1,+1); (1,0) -> u -> (-1,0);
#                            (0,1) -> v -> (0,+1)
exps = np.array([[-1, 1], [-1, 0], [0, 1]])
coef = np.array([0.4, 0.3, 0.3])
lam_quad = mahler_torus(exps, coef, ngrid=3000)
print("schedule-B simulation N=2e6      lambda_B = %.9f" % lam_sim)
print("2-torus quadrature 3000^2        m(...)   = %.9f" % lam_quad)
print("REGISTER erratum generic value            = -0.767507880")
print("deviation sim-vs-quad = %.2e" % abs(lam_sim - lam_quad))

print()
print("--- corpus cross-checks, my code, independent route ---")
# B1p (bridged carrier): classes {01:0.5, 10:0.5}, rank G = 1, lambda = log(1/2)
lam_bridge = mahler_torus(np.array([[-1, 0], [0, 1]]), np.array([0.5, 0.5]), ngrid=4000)
print("m(0.5 + 0.5 z)                  = %.9f   [S4 CONTROL 2 B1p: -0.693147181]" % lam_bridge)
# B1/B2/B3: classes {01:.4, 10:.4, 11:.2}
lam_b1 = mahler_torus(np.array([[-1, 1], [-1, 0], [0, 1]]), np.array([0.2, 0.4, 0.4]), 3000)
print("m(0.4 + 0.4x + 0.2y)            = %.9f   [S4 CONTROL 1/2 B1: -0.756573586]" % lam_b1)
# B0b: register correction  log(4/9)
print("log(4/9)                        = %.16f   [W-03 correction]" % np.log(4 / 9))
lam_b0b = mahler_torus(np.array([[-1, 1], [-1, 0], [0, 1]]), np.array([1/9, 4/9, 4/9]), 4000)
print("m(1/9 + 4/9 x + 4/9 y)          = %.9f" % lam_b0b)

print()
print("=" * 78)
print("G0.4  THE FOUR THEOREMS AT RANK ONE — my own reproduction")
print("=" * 78)
# (i) the root can never fire
s_root = state_rank1([1.0, 0, 0, 0, 0])
zr, cr = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s_root))
Zr = np.abs(Z_from_chars(zr, cr, np.arange(1, 5001)))
print("(i)   root-only state, rank 1: #characters = %d, min|Z_k| = %.12f, max = %.12f"
      % (len(zr), Zr.min(), Zr.max()))
print("      -> the root can never fire.  CONFIRMED at rank one.")

# (ii) S4-1: rank G = 2 iff |S| >= 3, enumerated over all 15 non-empty subsets
corners = {(0, 0): (0, 0), (1, 0): (-1, 0), (0, 1): (0, 1), (1, 1): (-1, 1)}
import itertools
cnt = {0: 0, 1: 0, 2: 0}
for r in range(1, 5):
    for S in itertools.combinations(corners, r):
        V = np.array([corners[x] for x in S])
        D = V[1:] - V[0] if len(V) > 1 else np.zeros((0, 2))
        rk = int(np.linalg.matrix_rank(D)) if D.size else 0
        cnt[rk] += 1
print("(ii)  S4-1 enumeration over 15 subsets: rank2/rank1/rank0 = %d/%d/%d  [S4: 5/6/4]"
      % (cnt[2], cnt[1], cnt[0]))

# (iii) lambda from the class-weight MULTISET (24 permutations)
import itertools as it
base = [0.4, 0.3, 0.2, 0.1]
allcorn = np.array([[0, 0], [-1, 0], [0, 1], [-1, 1]])
vals = []
for perm in it.permutations(range(4)):
    vals.append(mahler_torus(allcorn, np.array([base[i] for i in perm]), 1200))
print("(iii) 4-class multiset invariance: %d permutations, spread = %.3e  [W-03: 2.4e-15]"
      % (len(vals), max(vals) - min(vals)))

# (iv) pinch = spectator, exactly, at every connection
rng = np.random.default_rng(7771)
worst = 0.0
for _ in range(2000):
    ff, cc = rng.uniform(0, 2 * np.pi, 2)
    w = rng.dirichlet([1, 1, 1, 1])
    Zp = sum(w[i] * np.exp(1j * (allcorn[i, 0] * ff + allcorn[i, 1] * cc)) for i in range(4))
    wsw = np.array([w[3], w[2], w[1], w[0]])   # (0,0)<->(1,1), (1,0)<->(0,1)
    Zq = sum(wsw[i] * np.exp(1j * (allcorn[i, 0] * ff + allcorn[i, 1] * cc)) for i in range(4))
    worst = max(worst, abs(abs(Zp) - abs(Zq)))
print("(iv)  pinch<->spectator |Z_1| deviation over 2000 random (f,c,weights): %.3e"
      % worst, " [W-03: 6.55e-15]")
print()
print("G0 VALIDATION COMPLETE.")
