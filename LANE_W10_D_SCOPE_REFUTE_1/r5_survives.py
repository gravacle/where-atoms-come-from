# W10-D REFUTER 1 -- LEG R5.  WHAT I TRIED THAT COULD HAVE BROKEN THE LANE AND DID NOT.
#
# R5-A  D-01's two firing regions, on a DETERMINISTIC GRID instead of lane D's random draws
#       (different sampling scheme, same closed forms), plus the exact area argument.
# R5-B  D-03/D-04's torus-zero claim in EXACT RATIONALS, and min|Z_1| in CLOSED FORM
#       (lane D reported 1.11e-01 and 2.72e-01 as grid minima; they are 1/9 and sqrt(6)/9).
# R5-C  D-11/D-12's rates by direct simulation at connections that are NOT resonant.
# R5-D  D-19/D-20's dressed restoration on **B4**, a carrier lane D never built.  Lane D
#       exhibited it on B0b only; B4 has a different V, E, F, chi, b1, b2, tree and class
#       multiset, so this is a genuinely second arm, not a reprint.
# R5-E  D-22's monotonicity at generic connections, and the character identity in EXACT
#       Gaussian rationals (lane D labelled its float64 check "exact rationals"; it is not).
import numpy as np
from fractions import Fraction as Fr
from itertools import combinations
import collections

rng = np.random.default_rng(20260816)
CLS = ('00', '10', '01', '11')
EXP = {'00': (0, 0), '10': (1, 0), '01': (0, 1), '11': (1, 1)}
GEN = [("2pi(sqrt2-1), 2pi(sqrt3-1)", (2*np.pi*(np.sqrt(2)-1), 2*np.pi*(np.sqrt(3)-1))),
       ("(e, pi/e)                 ", (np.e, np.pi/np.e))]
for i in range(2):
    GEN.append((f"rng connection {i+1}          ", tuple(rng.uniform(-np.pi, np.pi, 2))))


def Zk(p, f, c, k):
    return sum(p[i]*np.exp(1j*k*(-EXP[CLS[i]][0]*f + EXP[CLS[i]][1]*c)) for i in range(4))


print("="*100)
print("== R5-A  D-01's FIRING REGIONS ON A DETERMINISTIC GRID (lane D used random draws) ==")
print("="*100)
n = 1201
g = (np.arange(n)+0.5)*2*np.pi/n - np.pi
F, C = np.meshgrid(g, g, indexing='ij')
F, C = F.ravel(), C.ravel()


def hull_maxgap(pts):
    A = np.sort(np.angle(np.stack(pts, 0)), axis=0)
    gap = np.diff(np.concatenate([A, A[:1]+2*np.pi], 0), axis=0)
    return gap.max(axis=0) <= np.pi + 1e-12


CH = {'00': np.ones_like(F, dtype=complex), '10': np.exp(-1j*F),
      '01': np.exp(1j*C), '11': np.exp(1j*(C-F))}
print(f"  grid {n} x {n} = {n*n} points, midpoint rule on (-pi,pi]^2 -- no seed, no draws")
for lab, occ, exact in [("three classes {10,01,11}  (B1,B2,B1s,B3)", ('10', '01', '11'), 0.25),
                        ("three classes {00,10,01}  (B1q)         ", ('00', '10', '01'), 0.25),
                        ("two   classes {10,01}     (B1p)         ", ('10', '01'), 0.0),
                        ("FOUR  classes             (B0b, B4)     ", ('00', '10', '01', '11'), 0.5)]:
    h = hull_maxgap([CH[o] for o in occ])
    print(f"  {lab}  fires {h.mean():.6f}   exact {exact}   |dev| {abs(h.mean()-exact):.2e}")
cf = (np.cos(F) + np.cos(C) <= 0)
h4 = hull_maxgap([CH[o] for o in ('00', '10', '01', '11')])
print(f"  four-class closed form cos f + cos c <= 0 agrees on {int((h4 == cf).sum())} of {n*n}")
print("  D-01 REPRODUCED under a different sampling scheme.  Both values exact.  SURVIVES.")

print("\n"+"="*100)
print("== R5-B  D-03/D-04 IN EXACT RATIONALS, AND THE TWO MINIMA IN CLOSED FORM ==")
print("="*100)
print("  P has a torus zero iff  alpha + beta cos t = 0  has a root in [-1,1], with")
print("  alpha = p00^2+p01^2-p10^2-p11^2 and beta = 2(p00 p01 - p10 p11).   ALL IN Fraction:")
for nm, pf in [("B0b (4,2,1,2)/9", [Fr(4, 9), Fr(2, 9), Fr(1, 9), Fr(2, 9)]),
               ("B4  (1,1,1,3)/6", [Fr(1, 6), Fr(1, 6), Fr(1, 6), Fr(3, 6)]),
               ("B1  K1 (0,2,2,1)/5", [Fr(0), Fr(2, 5), Fr(2, 5), Fr(1, 5)]),
               ("SENSE-C (1,1,1,1)/4", [Fr(1, 4)]*4)]:
    al = pf[0]**2 + pf[2]**2 - pf[1]**2 - pf[3]**2
    be = 2*(pf[0]*pf[2] - pf[1]*pf[3])
    if be == 0:
        verdict = "beta = 0 and alpha = 0 -> ZERO CURVE" if al == 0 else f"beta = 0, alpha = {al} != 0 -> NO ZERO"
    else:
        r = -al/be
        verdict = f"cos t = {r} -> {'ROOT' if abs(r) <= 1 else 'NO ROOT'}"
    print(f"  {nm:20s} alpha = {str(al):>8s}  beta = {str(be):>8s}   {verdict}")
print("  EXACT.  Both four-class carriers are torus-zero-FREE at their own SENSE-U weights and")
print("  SENSE-C has a zero CURVE.  D-03 and D-04 SURVIVE, in exact arithmetic.")
print("\n  AND THE MINIMA, WHICH LANE D REPORTED AS GRID NUMBERS -- they have closed forms:")
print("    |Z_1| >= | A(t) - B(t) | with A^2 - B^2 = alpha + beta cos t, so min |Z_1| =")
print("    min_t |alpha + beta cos t| / (A + B).")
print(f"    B0b: alpha = 1/9, beta = 0, max(A+B) = 5/9 + 4/9 = 1 at cos t = 1")
print(f"         -> min |Z_1| = 1/9 = {1/9:.12f}     lane D's grid value 1.111146e-01")
cstar = Fr(-2, 3)
val = (np.sqrt(10 + 6*float(cstar)) - np.sqrt(2 + 2*float(cstar)))/6
print(f"    B4 : stationary at cos t = {cstar}, giving min |Z_1| = sqrt(6)/9 = {np.sqrt(6)/9:.12f}")
print(f"         (check from the branch difference: {val:.12f})   lane D's grid value 2.721663e-01")
kk = np.arange(1, 200001)
for nm, p in [("B0b", np.array([4/9, 2/9, 1/9, 2/9])), ("B4 ", np.array([1/6, 1/6, 1/6, 3/6]))]:
    m = min(np.abs(Zk(p, fv, cv, kk)).min() for _, (fv, cv) in GEN)
    print(f"    {nm} min_k |Z_k| over k <= 2e5 at 4 generic connections: {m:.12f}")

print("\n"+"="*100)
print("== R5-C  D-11 / D-12's RATES, DIRECT, AT CONNECTIONS THAT ARE NOT RESONANT ==")
print("="*100)
print(f"  {'carrier':10s} {'connection':28s} {'(1/N)log|Omega_N|, N=2e6':>26s} {'closed form':>16s} {'dev':>10s}")
kk = np.arange(1, 2000001)
for nm, p, exact in [("B0b", np.array([4/9, 2/9, 1/9, 2/9]), np.log(4/9)),
                     ("B4 ", np.array([1/6, 1/6, 1/6, 3/6]), np.log(0.5)),
                     ("B1 ", np.array([0, 2/5, 2/5, 1/5]), -0.756573585640)]:
    for lab, (fv, cv) in GEN:
        r = np.log(np.abs(Zk(p, fv, cv, kk))).mean()
        print(f"  {nm:10s} {lab:28s} {r:26.9f} {exact:16.9f} {abs(r-exact):10.2e}")
print("  D-11 and D-12's ARITHMETIC SURVIVES at non-resonant connections.  log(4/9) and log(1/2)")
print("  are the rates; the lane's identification of N1 off K1 holds.")

print("\n"+"="*100)
print("== R5-D  D-19 / D-20 ON **B4** -- A CARRIER LANE D NEVER BUILT ==")
print("="*100)
Vn = 6
edges = [(0, 2), (2, 1), (1, 3), (3, 0), (0, 4), (4, 1), (1, 5), (5, 0)]
E = len(edges)
eidx = {e: k for k, e in enumerate(edges)}


def chain(seq):
    ch = np.zeros(E)
    for (s, t) in seq:
        if (s, t) in eidx:
            ch[eidx[(s, t)]] += 1
        else:
            ch[eidx[(t, s)]] -= 1
    return ch


gF = chain([(0, 2), (2, 1), (1, 3), (3, 0)])
gC = chain([(0, 2), (2, 1), (1, 4), (4, 0)])
FV = {v for v in range(Vn) if any(gF[k] != 0 and v in edges[k] for k in range(E))}
CV = {v for v in range(Vn) if any(gC[k] != 0 and v in edges[k] for k in range(E))}
print(f"  ARM DIFF vs lane D's leg 5 (which ran B0b):")
print(f"    lane D  B0b : V=9 E=18 F=9 chi=0 b1=2 b2=1  classes {{00:4,01:1,10:2,11:2}}")
print(f"    here    B4  : V={Vn} E={E} F=4 chi=2 b1=1 b2=2  classes "
      f"{{00:1,01:1,10:1,11:3}}  gamma_F {sorted(FV)}  gamma_C {sorted(CV)}")
adj = collections.defaultdict(list)
for k, (s, t) in enumerate(edges):
    adj[s].append((t, k, +1)); adj[t].append((s, k, -1))
tree = {0: []}
dq = collections.deque([0])
while dq:
    x = dq.popleft()
    for (y, k, sg) in adj[x]:
        if y not in tree:
            tree[y] = tree[x] + [(k, sg)]
            dq.append(y)
print(f"    spanning tree from root 0: " + ", ".join(f"v{v}:{[k for k,_ in tree[v]]}" for v in range(Vn)))


def holo(a, ch):
    return np.exp(1j*sum(sg*a[k] for k, sg in ch))


def loop_holo(a, g):
    return np.exp(1j*float(np.dot(g, a)))


def dress(a, s):
    return np.array([np.conjugate(holo(a, tree[v]))*s[v] for v in range(Vn)])


def A(a, s):
    t = dress(a, s)
    return np.outer(t.conjugate(), t)


a0 = rng.uniform(0, 2*np.pi, E)
s0 = rng.normal(size=Vn) + 1j*rng.normal(size=Vn)
s0 /= np.linalg.norm(s0)
worst_g = 0.0
for _ in range(200):
    th = rng.uniform(0, 2*np.pi, Vn)
    a1 = np.array([a0[k] + th[t] - th[s] for k, (s, t) in enumerate(edges)])
    s1 = np.exp(1j*th)*s0
    worst_g = max(worst_g, np.abs(np.abs(A(a1, s1)) - np.abs(A(a0, s0))).max())
print(f"  FULL gauge action (connection AND section), 200 random gauges: max |A_uv| dev = {worst_g:.3e}")
WF, WC = loop_holo(a0, gF), loop_holo(a0, gC)
sF = np.array([s0[v]*(WF if v in FV else 1) for v in range(Vn)])
sC = np.array([s0[v]*(WC if v in CV else 1) for v in range(Vn)])
diag = np.abs(np.abs(sF)**2 - np.abs(sC)**2).max()
AF, AC = A(a0, sF), A(a0, sC)
dressed = np.abs(AF-AC)
A0 = A(a0, s0)
pred = np.zeros((Vn, Vn))
for u in range(Vn):
    for v in range(Vn):
        da = (1 if v in FV else 0) - (1 if u in FV else 0)
        db = (1 if v in CV else 0) - (1 if u in CV else 0)
        pred[u, v] = abs(A0[u, v])*abs(WF**da - WC**db)
print(f"  S3's diagonal observables separate the branches by            {diag:.3e}")
print(f"  the DRESSED observable separates them by                      {dressed.max():.12f}")
print(f"  closed form |A_uv[s]| . |W_F^(da) - W_C^(db)| max dev over {Vn*Vn} pairs: "
      f"{np.abs(pred-dressed).max():.3e}")
print("  D-19 HOLDS ON A SECOND FOUR-CLASS CARRIER WITH DIFFERENT TOPOLOGY.  SURVIVES, and is")
print("  now exhibited on TWO arms rather than one.")
print("\n  AND THE RECURRENCE, D_k = max_(u,v) |A_uv[M_dF^k s] - A_uv[M_c^k s]|, k = 1..4000:")
K = 4000
kk4 = np.arange(1, K+1)
for lab, (wF, wC) in [("generic (random a_e)          ", (WF, WC)),
                      ("order 4: W_F=-1, W_C=-i       ", (-1+0j, -1j)),
                      ("order 2: W_F=-1, W_C=-1       ", (-1+0j, -1+0j)),
                      ("irrational e^{i sqrt2}, sqrt3 ", (np.exp(1j*np.sqrt(2)), np.exp(1j*np.sqrt(3))))]:
    D = np.zeros(K)
    for u in range(Vn):
        for v in range(Vn):
            da = (1 if v in FV else 0) - (1 if u in FV else 0)
            db = (1 if v in CV else 0) - (1 if u in CV else 0)
            D = np.maximum(D, abs(A0[u, v])*np.abs(wF**(da*kk4) - wC**(db*kk4)))
    print(f"    {lab} cells below 1e-9: {int((D < 1e-9).sum()):5d} of {K}   min {D.min():.3e}")
print("  1000 of 4000 at order 4 and 2000 of 4000 at order 2, ON B4 AS ON B0b AND AS ON K1.")
print("  D-20 SURVIVES on a second arm: it is the connection's arithmetic, not the carrier.")

print("\n"+"="*100)
print("== R5-E  D-22, AND THE ONE PRECISION LABEL THAT IS WRONG ==")
print("="*100)
K = 10**6
kk = np.arange(1, K+1)
for nm, p in [("B0b", np.array([4/9, 2/9, 1/9, 2/9])), ("B4 ", np.array([1/6, 1/6, 1/6, 3/6])),
              ("B1 ", np.array([0, 2/5, 2/5, 1/5]))]:
    w = max((np.abs(Zk(p, fv, cv, kk))-1).max() for _, (fv, cv) in GEN)
    print(f"  {nm}  max_k (|Z_k| - 1) over 4 GENERIC connections, k <= 1e6: {w:+.3e}")
print("  D-22's monotonicity SURVIVES (it is a triangle inequality; it could not have failed,")
print("  which voids the control and leaves the theorem standing -- the brief's own rule).")
print("\n  THE CHARACTER IDENTITY, NOW ACTUALLY IN EXACT ARITHMETIC.  Lane D's leg 3E prints")
print("  'Checked at FOUR classes in EXACT RATIONALS' and then evaluates np.exp(2j*pi*...) in")
print("  float64, reporting residual 3.331e-16 -- a float residual, not an exact check.  Redone")
print("  in Z[i] with q = 4 (characters in {1,i,-1,-i}), all (a,b), k = 1..8, weights in Fraction:")
pf = [Fr(4, 9), Fr(2, 9), Fr(1, 9), Fr(2, 9)]
POW = [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(-1), Fr(0)), (Fr(0), Fr(-1))]   # i^0..i^3 as (re,im)
worst = None
for a in range(4):
    for b in range(4):
        for k in range(1, 9):
            re = im = Fr(0)
            chs = []
            for i, cl in enumerate(CLS):
                ea, eb = EXP[cl]
                e = (k*(-ea*a + eb*b)) % 4
                cr, ci = POW[e]
                chs.append((cr, ci))
                re += pf[i]*cr; im += pf[i]*ci
            lhs = re*re + im*im
            rhs = Fr(1)
            for i in range(4):
                for j in range(i+1, 4):
                    dr = chs[i][0]-chs[j][0]; di = chs[i][1]-chs[j][1]
                    rhs -= pf[i]*pf[j]*(dr*dr + di*di)
            d = abs(lhs-rhs)
            worst = d if worst is None or d > worst else worst
print(f"    max |LHS - RHS| over q=4, all (a,b), k=1..8, in Fraction arithmetic: {worst}"
      f"   ({'EXACTLY ZERO' if worst == 0 else 'NONZERO'})")
print("  D-22's identity SURVIVES and is now exact.  The lane's arithmetic is right; its")
print("  precision LABEL was wrong, which under its own conventions page is a defect.")
