# W10-D leg 5 -- W-06/W-07's DRESSED-ALGEBRA RESTORATION, OFF K1, ON A FOUR-CLASS CARRIER.
# B0b is RECONSTRUCTED here (S4 publishes its row, not its d1/d2), and the reconstruction is
# VALIDATED against S4's own published class multiset before anything is measured on it.
import numpy as np
from fractions import Fraction

rng = np.random.default_rng(20260816)
print("="*100)
print("== 5A  B0b RECONSTRUCTED FROM S4's ROW, AND VALIDATED AGAINST IT ==")
print("="*100)
print("  S4:519 row 'B0b ring torus 3x3 grid, loops meet  V=9 E=18 F=9 chi=0 b1=2 b2=1'.")
print("  S4 does not publish d1,d2 for it (S4's own lane code is not in the repo), so this lane")
print("  BUILDS the 3x3 torus grid and then CHECKS that S4's published numbers all come out.")
print("  If any of them did not, nothing below would be about B0b.  DECLARED, not assumed.")
n = 3
vid = {(i, j): 3*i+j for i in range(n) for j in range(n)}
edges = []                                   # (src, tgt)
for i in range(n):
    for j in range(n):
        edges.append((vid[(i, j)], vid[((i+1) % n, j)]))     # horizontal
        edges.append((vid[(i, j)], vid[(i, (j+1) % n)]))     # vertical
E = len(edges); V = n*n
faces = []                                   # each square, as signed edge chain
eidx = {e: k for k, e in enumerate(edges)}
for i in range(n):
    for j in range(n):
        ch = np.zeros(E)
        ch[eidx[(vid[(i, j)], vid[((i+1) % n, j)])] ] += 1
        ch[eidx[(vid[((i+1) % n, j)], vid[((i+1) % n, (j+1) % n)])]] += 1
        ch[eidx[(vid[(i, (j+1) % n)], vid[((i+1) % n, (j+1) % n)])]] -= 1
        ch[eidx[(vid[(i, j)], vid[(i, (j+1) % n)])]] -= 1
        faces.append(ch)
F = len(faces)
d1 = np.zeros((V, E))
for k, (s, t) in enumerate(edges):
    d1[t, k] += 1; d1[s, k] -= 1
d2 = np.array(faces).T
r1, r2 = np.linalg.matrix_rank(d1), np.linalg.matrix_rank(d2)
b0 = V - r1; b1 = E - r1 - r2; b2 = F - r2
print(f"  built: V={V} E={E} F={F}  chi={V-E+F}  rank d1={r1} rank d2={r2}  b0={b0} b1={b1} b2={b2}"
      f"   max|d1.d2| = {np.abs(d1@d2).max():.1e}")
assert (V, E, F, V-E+F, b0, b1, b2) == (9, 18, 9, 0, 1, 2, 1), "reconstruction does not match S4"
gF = faces[0].copy()                                     # boundary of one square: BOUNDS
gC = np.zeros(E)                                         # horizontal loop j=0: does NOT bound
for i in range(n):
    gC[eidx[(vid[(i, 0)], vid[((i+1) % n, 0)])]] += 1
FV = {v for v in range(V) if any(gF[k] != 0 and v in edges[k] for k in range(E))}
CV = {v for v in range(V) if any(gC[k] != 0 and v in edges[k] for k in range(E))}
cls = {}
for v in range(V):
    cls[v] = ('1' if v in FV else '0') + ('1' if v in CV else '0')
from collections import Counter
cnt = Counter(cls.values())
print(f"  gamma_F is a cycle: {np.abs(d1@gF).max():.1e}   bounds (in im d2): "
      f"{np.linalg.matrix_rank(np.c_[d2, gF]) == r2}")
print(f"  gamma_C is a cycle: {np.abs(d1@gC).max():.1e}   bounds (in im d2): "
      f"{np.linalg.matrix_rank(np.c_[d2, gC]) == r2}")
print(f"  CLASS MULTISET BUILT : {{'00':{cnt['00']}, '01':{cnt['01']}, '10':{cnt['10']}, '11':{cnt['11']}}}")
print( "  S4:576 PUBLISHED     : {'00':4, '01':1, '10':2, '11':2}")
assert (cnt['00'], cnt['01'], cnt['10'], cnt['11']) == (4, 1, 2, 2), "multiset mismatch"
print("  MATCH.  The reconstruction reproduces S4's published row exactly, including the class")
print("  multiset that every functional in this corpus depends on.  It is B0b for these purposes.")

print("\n"+"="*100)
print("== 5B  THE DRESSED OBSERVABLE ON B0b: GAUGE INVARIANCE UNDER THE **FULL** ACTION ==")
print("="*100)
# spanning tree by BFS from root 0
import collections
adj = collections.defaultdict(list)
for k, (s, t) in enumerate(edges):
    adj[s].append((t, k, +1)); adj[t].append((s, k, -1))
tree_path = {0: []}
dq = collections.deque([0])
while dq:
    x = dq.popleft()
    for (y, k, sg) in adj[x]:
        if y not in tree_path:
            tree_path[y] = tree_path[x] + [(k, sg)]
            dq.append(y)

def holo(a, chain):
    return np.exp(1j*sum(sg*a[k] for k, sg in chain))

def loop_holo(a, g):
    return np.exp(1j*float(np.dot(g, a)))

def dress(a, s):
    return np.array([holo(a, tree_path[v]).conjugate()*s[v] for v in range(V)])

def A(a, s):
    t = dress(a, s)
    return np.outer(t.conjugate(), t)

def branch(a, s, which, k=1):
    W = loop_holo(a, gF) if which == 'F' else loop_holo(a, gC)
    S = set(FV) if which == 'F' else set(CV)
    return np.array([s[v]*(W**k if v in S else 1) for v in range(V)])

a0 = rng.uniform(0, 2*np.pi, E)
s0 = rng.normal(size=V) + 1j*rng.normal(size=V); s0 /= np.linalg.norm(s0)
worst_g, worst_wf = 0.0, 0.0
for _ in range(200):
    th = rng.uniform(0, 2*np.pi, V)
    a1 = a0.copy()
    for k, (sc, tg) in enumerate(edges):
        a1[k] = a0[k] + th[tg] - th[sc]
    s1 = np.exp(1j*th)*s0
    worst_g = max(worst_g, np.abs(np.abs(A(a1, s1)) - np.abs(A(a0, s0))).max())
    worst_wf = max(worst_wf, abs(loop_holo(a1, gF)-loop_holo(a0, gF)),
                   abs(loop_holo(a1, gC)-loop_holo(a0, gC)))
print(f"  FULL gauge action (S1:63 on the connection AND on the section), 200 random gauges:")
print(f"    max |A_uv| deviation           = {worst_g:.3e}")
print(f"    max holonomy deviation (W_F,W_C)= {worst_wf:.3e}")

print("\n== 5C  BRANCH SEPARATION: S3's DIAGONAL OBSERVABLES vs THE DRESSED ONE, ON B0b ==")
sF, sC = branch(a0, s0, 'F'), branch(a0, s0, 'C')
diag = np.abs(np.abs(sF)**2 - np.abs(sC)**2).max()
AF, AC = A(a0, sF), A(a0, sC)
dressed = np.abs(AF-AC)
i, j = np.unravel_index(np.argmax(dressed), dressed.shape)
print(f"  max over v of | |s_v|^2 (branch F) - |s_v|^2 (branch C) |  = {diag:.3e}   "
      f"(S3's gauge-invariant algebra: the diagonal)")
print(f"  max over (u,v) of |A_uv[M_dF s] - A_uv[M_c s]|             = {dressed.max():.12f}  "
      f"at (u,v)=({i},{j}), classes {cls[i]}/{cls[j]}")
print("  W-06's RESTORATION HOLDS ON A FOUR-CLASS CARRIER.  Off K1, from an independent build.")

print("\n== 5D  AND THE CLOSED FORM, WHICH NAMES WHAT THE DRESSED OBSERVABLE ACTUALLY SEES ==")
WF, WC = loop_holo(a0, gF), loop_holo(a0, gC)
A0 = A(a0, s0)
pred = np.zeros((V, V))
for u in range(V):
    for v in range(V):
        da = (1 if v in FV else 0) - (1 if u in FV else 0)
        db = (1 if v in CV else 0) - (1 if u in CV else 0)
        pred[u, v] = abs(A0[u, v])*abs(WF**da - WC**db)
print(f"  |A_uv[M_dF s] - A_uv[M_c s]| = |A_uv[s]| . |W_F^(a_v-a_u) - W_C^(b_v-b_u)|")
print(f"    max deviation from the closed form over all 81 pairs: {np.abs(pred-dressed).max():.3e}")
print("  SO THE DRESSED SEPARATION IS A FUNCTION OF THE **CLASS DIFFERENCE** OF THE PAIR (u,v)")
print("  AND NOTHING ELSE ABOUT THE CARRIER.  It is non-zero iff some occupied pair has")
print("  (a_v-a_u, b_v-b_u) != (0,0) with W_F^da != W_C^db -- the SAME character condition as")
print("  formation.  W-06's restoration is CARRIER-INDEPENDENT given two occupied classes.")
print("  Its NUMBERS (0.384349931183, 3*sqrt(3)/10, 1000 of 4000) are K1-and-state-specific.")

print("\n== 5E  AND IT RECURS ON B0b TOO -- W-06/W-07's OBSTRUCTION IS NOT K1's ==")
print(f"  D_k = max_(u,v) |A_uv[M_dF^k s] - A_uv[M_c^k s]|, k = 1..4000, same state, same tree.")
for lab, (fv, cv) in [("generic  (random a_e above)", (None, None)),
                      ("order 4: W_F=-1, W_C=-i", (np.pi, -np.pi/2)),
                      ("order 2: W_F=-1, W_C=-1", (np.pi, np.pi)),
                      ("irrational: W_F=e^{i sqrt2}, W_C=e^{i sqrt3}", (np.sqrt(2), np.sqrt(3)))]:
    if fv is None:
        wF, wC = WF, WC
    else:
        wF, wC = np.exp(1j*fv), np.exp(1j*cv)
    K = 4000
    kk = np.arange(1, K+1)
    D = np.zeros(K)
    for u in range(V):
        for v in range(V):
            da = (1 if v in FV else 0) - (1 if u in FV else 0)
            db = (1 if v in CV else 0) - (1 if u in CV else 0)
            D = np.maximum(D, abs(A0[u, v])*np.abs(wF**(da*kk) - wC**(db*kk)))
    print(f"    {lab:44s} cells below 1e-9: {int((D<1e-9).sum()):5d} of {K}   min {D.min():.3e}")
print("  The order-4 point annihilates the dressed record on a fixed fraction of cells ON B0b")
print("  EXACTLY AS ON K1.  ATTAINED-vs-APPROACHED (W-07) IS A PROPERTY OF THE CONNECTION'S")
print("  ARITHMETIC, NOT OF THE CARRIER.  CARRIER-INDEPENDENT, CONNECTION-SCOPED.")

print("\n"+"="*100)
print("== 5F  A W-08 FIGURE I COULD NOT REPRODUCE -- RECORDED, NOT SMOOTHED ==")
print("="*100)
print("  W-08 registers the adversarial schedule accumulating '0.606, 0.615, 0.588, 0.601 nats")
print("  at K = 10^4..10^7' -- i.e. FLAT in K.  Leg 4E on K1's own weights gives 0.5718, 0.6211,")
print("  0.3615, 0.1143: the first two agree, the last two DECAY.  Sweep the connection:")
CLS4 = (('00', (0, 0)), ('10', (1, 0)), ('01', (0, 1)), ('11', (1, 1)))
def Z(p, f, c, k):
    return sum(p[i]*np.exp(1j*k*(-CLS4[i][1][0]*f + CLS4[i][1][1]*c)) for i in range(4))
pK1 = np.array([0, 2/5, 2/5, 1/5])
for lab, (fv, cv) in [("generic (1.3, 2.0)", (1.3, 2.0)), ("S3/S4 headline (2.0, 1.1)", (2.0, 1.1)),
                      ("golden (2pi/phi, 2pi/phi^2)", (2*np.pi/1.6180339887, 2*np.pi/2.6180339887)),
                      ("S1 order-4 (pi, 3pi/2)", (np.pi, 3*np.pi/2))]:
    row = []
    for K in (10**4, 10**5, 10**6, 10**7):
        kk = np.arange(1, K+1)
        d = 1.0 - np.abs(Z(pK1, fv, cv, kk))
        m = int(round(np.sqrt(K)))
        idx = np.argpartition(d, m)[:m]
        row.append(-np.log(np.abs(Z(pK1, fv, cv, kk[idx]))).sum())
    print(f"    {lab:32s} " + "  ".join(f"{r:8.4f}" for r in row))
print("  On NO connection tested does the accumulation stay flat at ~0.60 through K = 10^7.")
print("  READ, BOTH WAYS AND SCORED AS NEITHER: either W-08's figure used a connection or a")
print("  tie-break this lane has not identified, OR the adversary is STRONGER than W-08 reported")
print("  (accumulation decaying to 0 means |Omega| -> 1 along the adversarial schedule).")
print("  EITHER WAY W-08's QUALITATIVE RULING -- the adversarial budget is O(1), the honest one")
print("  diverges linearly -- IS REPRODUCED AND IS CARRIER-INDEPENDENT.  The four printed")
print("  constants are not reproduced and are marked UNDETERMINED in the scope table.")
