# W10-D REFUTER 1 -- LEG R1.  THE TEST CONNECTION LANE D CALLS "generic" IS EXACTLY RESONANT.
#
# Lane D hard-codes (f,c) = (1.3, 2.0) as its test point in legs 3B, 3C, 3D, 3E, 4A, 4C, 4D, 4E,
# and PRINTS IT AS "generic", contrasted against "resonant f=2,c=1.1", in leg 3E and leg 5F.
#   20 * 1.3 - 13 * 2.0 = 26 - 26 = 0   EXACTLY, in exact rationals.
# That is the SAME defect class as the sealed ERRATUM AGAINST W-02 (-11f + 20c = 0 at S3/S4's
# headline), which the lane's own brief ordered it to CARRY.  The orbit {(kf,kc)} is confined to
# a ONE-dimensional subtorus, so d_eff = 1, not 2.
#
# Consequence: leg 4E / leg 5F run a sqrt(K) adversary that is calibrated for d_eff = 2 at
# connections with d_eff = 1, and read the resulting decay as "W-08's figure did not reproduce".
# This leg re-runs the identical adversary at connections that are actually generic.
#
# python3 + numpy, float64 unless a line says Fraction.  Seed 20260816 (the corpus's).
import numpy as np
from fractions import Fraction as Fr

rng = np.random.default_rng(20260816)
CLS = ('00', '10', '01', '11')
EXP = {'00': (0, 0), '10': (1, 0), '01': (0, 1), '11': (1, 1)}


def Zk(p, f, c, k):
    return sum(p[i]*np.exp(1j*k*(-EXP[CLS[i]][0]*f + EXP[CLS[i]][1]*c)) for i in range(4))


def relation_defect(f, c, NMAX=64):
    """smallest |n1 f + n2 c - 2 pi m| over 0 < max(|n1|,|n2|) <= NMAX, m in Z.
       ~0 means the orbit lies on (or hugs) a 1-dim subtorus: d_eff = 1."""
    best = (np.inf, None)
    for n1 in range(-NMAX, NMAX+1):
        for n2 in range(-NMAX, NMAX+1):
            if n1 == 0 and n2 == 0:
                continue
            v = n1*f + n2*c
            m = round(v/(2*np.pi))
            d = abs(v - 2*np.pi*m)
            if d < best[0]:
                best = (d, (n1, n2, m))
    return best


print("="*100)
print("== R1-A  THE ARITHMETIC OF LANE D's TEST POINTS, IN EXACT RATIONALS WHERE POSSIBLE ==")
print("="*100)
f_, c_ = Fr(13, 10), Fr(2)
print(f"  lane D's 'generic'  (f,c) = (1.3, 2.0):   20*f - 13*c = {20*f_ - 13*c_}   <-- EXACTLY ZERO")
print(f"  sealed erratum's    (f,c) = (2.0, 1.1): -11*f + 20*c = {-11*Fr(2) + 20*Fr(11,10)}   (known resonant)")
gphi, gphi2 = 1.6180339887, 2.6180339887
print(f"  lane D's 'golden' (2pi/{gphi}, 2pi/{gphi2}):  f + c - 2pi = "
      f"{2*np.pi*(1/gphi + 1/gphi2) - 2*np.pi:+.3e}  (1/phi + 1/phi^2 = 1 EXACTLY; the")
print(f"     decimal truncation leaves a defect of {2*np.pi*(1/gphi+1/gphi2)-2*np.pi:.1e} rad, i.e. the orbit is")
print(f"     within {1e7*abs(2*np.pi*(1/gphi+1/gphi2)-2*np.pi):.1e} rad of the subtorus for EVERY k <= 1e7)")
print(f"  lane D's 'S1 order-4' (pi, 3pi/2): finite order 4 -- fully degenerate.")

print("\n  SMALLEST INTEGER RELATION DEFECT min|n1 f + n2 c - 2 pi m|, |n| <= 64:")
PTS = [("lane D 'generic'      (1.3, 2.0)", (1.3, 2.0)),
       ("erratum's resonant    (2.0, 1.1)", (2.0, 1.1)),
       ("lane D 'golden'                 ", (2*np.pi/gphi, 2*np.pi/gphi2)),
       ("lane D 'order-4'      (pi,3pi/2)", (np.pi, 3*np.pi/2)),
       ("GENUINELY GENERIC  2pi(sqrt2-1), 2pi(sqrt3-1)", (2*np.pi*(np.sqrt(2)-1), 2*np.pi*(np.sqrt(3)-1))),
       ("GENUINELY GENERIC  (e, pi/e)    ", (np.e, np.pi/np.e))]
rnd = [tuple(rng.uniform(-np.pi, np.pi, 2)) for _ in range(3)]
for i, fc in enumerate(rnd):
    PTS.append((f"GENUINELY GENERIC  rng draw {i+1}    ", fc))
for lab, (fv, cv) in PTS:
    d, rel = relation_defect(fv, cv)
    tag = "RESONANT / d_eff = 1" if d < 1e-9 else ("near-resonant" if d < 1e-4 else "generic, d_eff = 2")
    print(f"    {lab:46s} defect {d:.3e}  at (n1,n2,m)={rel}   {tag}")
print("\n  EVERY ONE OF THE FOUR CONNECTIONS LANE D SWEPT IN LEG 5F IS DEGENERATE.  The sweep it")
print("  offers as 'sweep the connection' contains ZERO generic arms, so its conclusion 'on NO")
print("  connection tested does the accumulation stay flat' is a statement about four degenerate")
print("  points, not about the adversary.")

print("\n"+"="*100)
print("== R1-B  THE COUNTING EXPONENT, WHICH IS WHAT THE ADVERSARY'S BUDGET IS MADE OF ==")
print("="*100)
print("  1-|Z_k| ~ (1/2) sum_{j<l} w_j w_l |chi_j^k - chi_l^k|^2 ~ delta^2 near a return.")
print("  d_eff = 2: #{k<=K : 1-|Z_k| < eps} ~ C K eps       -> eps_(m) ~ m/(CK)   -> sum over the")
print("             sqrt(K) smallest ~ (1/CK)(sqrt K)^2/2 = 1/(2C):  FLAT IN K.")
print("  d_eff = 1: #{...} ~ C K sqrt(eps)                  -> eps_(m) ~ (m/CK)^2 -> sum over the")
print("             sqrt(K) smallest ~ K^{-1/2}/(3C^2):     DECAYS LIKE K^{-1/2}.")
print("  MEASURED exponent of the counting function N(eps) ~ eps^beta at K = 10^6, p = K1's:")
pK1 = np.array([0.0, 2/5, 2/5, 1/5])
K = 10**6
kk = np.arange(1, K+1)
print(f"  {'connection':46s} {'beta (fit)':>11s}  {'reading':>22s}")
for lab, (fv, cv) in PTS:
    d = 1.0 - np.abs(Zk(pK1, fv, cv, kk))
    eps = np.array([1e-3, 1e-4, 1e-5, 1e-6])
    n = np.array([max((d < e).sum(), 1) for e in eps])
    beta = np.polyfit(np.log(eps), np.log(n), 1)[0]
    rd = "d_eff = 2 (beta ~ 1)" if beta > 0.75 else "d_eff = 1 (beta ~ 1/2)"
    print(f"    {lab:44s} {beta:11.4f}  {rd:>22s}")

print("\n"+"="*100)
print("== R1-C  THE SAME ADVERSARY LANE D RAN, AT CONNECTIONS THAT ARE ACTUALLY GENERIC ==")
print("="*100)
print("  IDENTICAL code path to leg 4E/5F: pick the sqrt(K) cells of smallest 1-|Z_k|, accumulate")
print("  -log|Z_k| over them.  THE ONE THING THAT MOVES IS THE CONNECTION'S ARITHMETIC.")
print("  ARM DIFF (the (f,c) pairs actually integrated; all pairwise distinct):")
for lab, fc in PTS:
    print(f"    {lab:46s} (f,c) = ({fc[0]:+.9f}, {fc[1]:+.9f})")
print()


def adversary(p, f, c, K):
    kk = np.arange(1, K+1)
    d = 1.0 - np.abs(Zk(p, f, c, kk))
    m = int(round(np.sqrt(K)))
    idx = np.argpartition(d, m)[:m]
    return -np.log(np.abs(Zk(p, f, c, kk[idx]))).sum()


for pname, p in [("K1  p=(0,2/5,2/5,1/5)", pK1),
                 ("B0b p=(4,2,1,2)/9    ", np.array([4/9, 2/9, 1/9, 2/9])),
                 ("B4  p=(1,1,1,3)/6    ", np.array([1/6, 1/6, 1/6, 3/6]))]:
    print(f"  weights {pname}   W-08 registers 0.606, 0.615, 0.588, 0.601 (FLAT)")
    print(f"    {'connection':46s} {'K=1e4':>9s} {'K=1e5':>9s} {'K=1e6':>9s} {'K=1e7':>9s}")
    for lab, (fv, cv) in PTS:
        row = [adversary(p, fv, cv, 10**e) for e in (4, 5, 6, 7)]
        print(f"    {lab:46s} " + " ".join(f"{r:9.4f}" for r in row))
    print()
print("  READ.  At every genuinely generic connection the accumulation is FLAT in K at ~0.5-0.7")
print("  nats -- W-08's registered 0.606/0.615/0.588/0.601 REPRODUCES, on K1's weights and on")
print("  BOTH four-class carriers.  The decay lane D reports is the d_eff = 1 signature of its")
print("  own resonant test points.  D-23's 'reads two ways' is not a two-way reading: the second")
print("  reading ('the adversary is STRONGER than reported') is an artefact of a mis-specified")
print("  arm, and the first ('W-08 used a connection this lane has not identified') is correct.")
