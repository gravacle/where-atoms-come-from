"""SECOND-LINE ADVERSARIAL PROBE for T-44-C (external verifier, independent lineage).

The lane's own VERIFY reimplemented the transfer-matrix series.  This probe attacks with
DIFFERENT MACHINERY entirely:
  P1  Dense float64 LINEAR SOLVE  G = (I - mu.Lam.A)^{-1} e_S  (numpy) -- no series,
      no transfer matrix -- checked against the sealed OUT's printed decimals.
      (Valid because the sealed tail bounds put |series - limit| far below 1e-12 at
      the printed values' scale; a disagreement at ~1e-9 relative refutes the numbers.)
  P2  Exact integer per-order counts to w=30 by an independently written iteration:
      R1 additivity identity, mirror count equality, w_min=d, parity, shadow onset.
  P3  TRUNCATION HONESTY, exact: one run to W=820, partial(500) vs partial(820);
      remainder must be > 0 (bounds not display zeros) and <= the sealed tail formula
      K*(mu*B)^(W+1)/(1-mu*B) with B recomputed here (different iteration depth).
  P4  FRESH-MU ATTACK with rows NEITHER the lane NOR its verify used:
      3-D mu = 1/9 (between the lane's rows) and mu = 9/50 = 0.18 -- ABOVE the
      infinite-lattice landmark 1/6 yet still certified convergent on the finite
      venue (mu*B < 1): if the composition classes were secretly tied to the
      borrowed 'critical' landmark, a super-landmark row should break them.
      2-D mu = 13/50 = 0.26, ABOVE the 2-D landmark 1/4, certified convergent
      (a first try at 27/100 FAILED the certificate, mu*B > 1, and the resolvent
      went unphysical -- kept as a note: the certificate does real work).  Classes checked: R2 defect < 0, partial
      dipole screening, occupancy positivity+additivity, mirror zero, shadow R<1
      with axial minimum.
  P5  Cross-check B against the exact grid spectral radius 2*sum_i cos(pi/(n_i+1))
      (separable grid Laplacian spectrum -- standard): B must be >= rho and close.
Every check prints a computed boolean; exit 0 iff all OK.
"""
import sys
import math
import itertools
from fractions import Fraction
from collections import deque
import numpy as np

CH = []
def chk(tag, ok, msg):
    CH.append((tag, bool(ok)))
    print("  [%s] %s -- %s" % ("OK " if ok else "BAD", tag, msg))

def venue(dims):
    coords = list(itertools.product(*[range(d) for d in dims]))
    idx = {c: i for i, c in enumerate(coords)}
    A = np.zeros((len(coords), len(coords)))
    nbrs = [[] for _ in coords]
    for c in coords:
        for ax in range(len(dims)):
            q = list(c); q[ax] += 1; q = tuple(q)
            if q in idx:
                A[idx[c], idx[q]] = 1; A[idx[q], idx[c]] = 1
                nbrs[idx[c]].append(idx[q]); nbrs[idx[q]].append(idx[c])
    return coords, idx, nbrs, A

def bfsd(nbrs, srcs, blocked=frozenset()):
    dist = {s: 0 for s in srcs}
    dq = deque(srcs)
    while dq:
        u = dq.popleft()
        for v in nbrs[u]:
            if v not in dist and v not in blocked:
                dist[v] = dist[u] + 1; dq.append(v)
    return dist

def gsolve(A, mu, src, lam0=(), lamhalf=()):
    """G(src -> .) = sum_{w>=0} mu^w y_w with y_w = (Lam A)^T-applied... implemented as
       the linear solve (I - mu * Lam @ A) x = e_src summed over src; Lam applies on
       ENTERING a site (start unweighted): y_w[v] = lam_v * sum_u A[u,v] y_{w-1}[u].
       So x = e_S + mu*Lam A^T x ... with symmetric A: x = (I - mu*diag(lam)@A)^{-1} e_S."""
    n = A.shape[0]
    lam = np.ones(n)
    for s in lamhalf: lam[s] = 0.5
    for s in lam0: lam[s] = 0.0
    e = np.zeros(n)
    for s in src: e[s] = 1.0
    M = np.eye(n) - mu * (lam[:, None] * A)   # row v gets lam_v * sum_u A[v,u] y[u]
    return np.linalg.solve(M, e)

def close(a, b, rel=1e-9):
    return abs(a - b) <= rel * max(abs(a), abs(b), 1e-300)

print("=" * 78)
print("T-44-C SECOND-LINE ADVERSARIAL PROBE -- linear-solve + exact counts, fresh mu")
print("=" * 78)

co3, ix3, nb3, A3 = venue((11, 5, 5))
co2, ix2, nb2, A2 = venue((13, 7))

# ---- P5 spectral: exact grid rho vs sealed B
rho3 = 2 * (math.cos(math.pi / 12) + math.cos(math.pi / 6) + math.cos(math.pi / 6))
rho2 = 2 * (math.cos(math.pi / 14) + math.cos(math.pi / 8))
ev3 = max(abs(np.linalg.eigvalsh(A3)))
ev2 = max(abs(np.linalg.eigvalsh(A2)))
B3_sealed, B2_sealed = 5.396145, 3.797765  # printed (truncated) in sealed OUT

# independent Collatz-Wielandt at a DIFFERENT depth (120 iters, exact rationals)
def cw(nbrs, iters=120):
    n = len(nbrs)
    phi = [1] * n
    for _ in range(iters):
        phi = [phi[v] + sum(phi[u] for u in nbrs[v]) for v in range(n)]
    B = max(Fraction(sum(phi[u] for u in nbrs[v]), phi[v]) for v in range(n))
    return B, phi

B3, phi3 = cw(nb3)
B2, phi2 = cw(nb2)
chk("P5.a", abs(ev3 - rho3) < 1e-9 and abs(ev2 - rho2) < 1e-9,
    "numeric spectral radius matches the separable-grid formula 2*sum cos(pi/(n+1)) "
    "(rho3=%.6f rho2=%.6f)" % (rho3, rho2))
chk("P5.b", float(B3) >= ev3 and float(B3) - ev3 < 1e-4 and
            float(B2) >= ev2 and float(B2) - ev2 < 1e-4,
    "independent CW bound (120 iters) is a genuine upper bound, tight to <1e-4 "
    "(B3=%.6f B2=%.6f)" % (float(B3), float(B2)))
chk("P5.c", B3_sealed >= ev3 - 1e-9 and B3_sealed - ev3 < 1e-3 and
            B2_sealed >= ev2 - 1e-9 and B2_sealed - ev2 < 1e-3,
    "sealed printed B (80 iters: 5.396145 / 3.797765) is a GENUINE upper bound on "
    "the true rho (5.395953 / 3.797615), looser than my 120-iter bound as CW "
    "iteration predicts -- conservative, i.e. sealed tail bounds only over-cover")

# ---- geometry (transcribed from sealed OUT text, not from the model source)
A1 = [ix3[(1, 2, 2)]]; C3s = [ix3[(10, 2, 2)]]
Bsep = {s: [ix3[(1 + s, 2, 2)]] for s in (2, 3, 4, 6)}
Am = [ix3[(1, 2, 2)]]; Bm = [ix3[(9, 2, 2)]]
plane5 = [ix3[(5, y, z)] for y in range(5) for z in range(5)]
Dp = [ix3[(2, 2, 2)]]; Dm = [ix3[(3, 2, 2)]]
occ = [ix3[(2 + i, 2, 2)] for i in range(4)]
slab = [ix3[(5, y, z)] for y in range(5) for z in range(5)]
disk = [ix3[(5, y, z)] for y in (1, 2, 3) for z in (1, 2, 3)]
center = [ix3[(5, 2, 2)]]
slab1 = [s for s in slab if s not in center]
ring = [s for s in slab if s not in disk]

# ---- P2 exact per-order counts (independent iteration, weights as Fractions)
def orders(nbrs, src, Wmax, lam0=frozenset(), lamhalf=frozenset()):
    n = len(nbrs)
    lam = {}
    for s in lamhalf: lam[s] = Fraction(1, 2)
    for s in lam0: lam[s] = Fraction(0)
    y = [Fraction(0)] * n
    for s in src: y[s] += 1
    out = [y[:]]
    for w in range(Wmax):
        yn = [Fraction(0)] * n
        for v in range(n):
            lv = lam.get(v, 1)
            if lv:
                acc = sum(y[u] for u in nbrs[v])
                if acc: yn[v] = lv * acc
        y = yn
        out.append(y[:])
    return out

W = 30
oA = orders(nb3, A1, W)
d_AC = min(bfsd(nb3, A1)[c] for c in C3s)
cnt = lambda snaps, w, T: sum(snaps[w][t] for t in T)
chk("P2.a", d_AC == 9 and all(cnt(oA, w, C3s) == 0 for w in range(9)) and
    cnt(oA, 9, C3s) == 1 and all(cnt(oA, w, C3s) == 0 for w in range(9, W + 1) if (w - 9) % 2),
    "w_min = d = 9 with N_9 = 1 and bipartite parity (independent count)")

ok_r1 = True
for s in (2, 3, 4, 6):
    oB = orders(nb3, Bsep[s], W)
    oU = orders(nb3, A1 + Bsep[s], W)
    ok_r1 &= all(cnt(oU, w, C3s) == cnt(oA, w, C3s) + cnt(oB, w, C3s)
                 for w in range(W + 1))
chk("P2.b", ok_r1, "R1 per-order additivity identity, seps 2/3/4/6, independent counts")

oMp = orders(nb3, Am, W); oMm = orders(nb3, Bm, W)
chk("P2.c", all(cnt(oMp, w, [c]) == cnt(oMm, w, [c]) for w in range(W + 1)
                for c in plane5),
    "mirror-pair count equality at every mid-plane site, every order (exact zero Phi)")

oD = orders(nb3, A1, W, lam0=frozenset(disk))
d_punc = min(bfsd(nb3, A1, blocked=frozenset(disk))[b] for b in Bm)
fo = next((w for w in range(W + 1) if cnt(oD, w, Bm) != 0), None)
d_free = min(bfsd(nb3, A1)[b] for b in Bm)
chk("P2.d", fo == 12 and d_punc == 12 and d_free == 8,
    "shadow onset %s = punctured BFS %d > free %d (independent)" % (fo, d_punc, d_free))

oS = orders(nb3, A1, W, lam0=frozenset(slab))
cut = all(b not in bfsd(nb3, A1, blocked=frozenset(slab)) for b in Bm)
chk("P2.e", cut and all(cnt(oS, w, Bm) == 0 for w in range(W + 1)),
    "full slab: vertex cut and zero count at every order (independent)")

# ---- P1 linear-solve reproduction of printed decimals (mu rows of the sealed OUT)
mu12, mu6 = 1.0 / 12.0, 1.0 / 6.0
g6 = {}
g6['free'] = gsolve(A3, mu6, A1)
g6['cen'] = gsolve(A3, mu6, A1, lam0=center)
g6['disk'] = gsolve(A3, mu6, A1, lam0=disk)
g6['slab1'] = gsolve(A3, mu6, A1, lam0=slab1)
g6['ring'] = gsolve(A3, mu6, A1, lam0=ring)
b = Bm[0]
seal = {'free': 0.0012367408720414, 'cen': 0.0007977198503738,
        'disk': 0.0001313921462753, 'slab1': 0.0000554135496265,
        'ring': 0.0006292395119519}
okP1 = all(close(g6[k][b], v, 1e-8) for k, v in seal.items())
chk("P1.a", okP1, "shadow-chain values reproduce by LINEAR SOLVE (no series): " +
    ", ".join("%s=%.10e" % (k, g6[k][b]) for k in seal))

r_ax = g6['disk'][ix3[(9, 2, 2)]] / g6['free'][ix3[(9, 2, 2)]]
r_co = g6['disk'][ix3[(9, 0, 0)]] / g6['free'][ix3[(9, 0, 0)]]
chk("P1.b", close(r_ax, 0.10624064, 1e-6) and close(r_co, 0.15402869, 1e-6),
    "far-plane ratios axial %.8f / corner %.8f reproduce (sealed 0.10624064/0.15402869)"
    % (r_ax, r_co))

# R2 defect table by linear solve
okd = True
vals = {}
seal_d = {(2, mu12): -0.008174820104, (3, mu12): -0.000784482164,
          (4, mu12): -0.000077104643, (6, mu12): -0.000000808188,
          (2, mu6): -0.088532997166, (3, mu6): -0.038627719073,
          (4, mu6): -0.017807581541, (6, mu6): -0.004065283309}
for s in (2, 3, 4, 6):
    for mu in (mu12, mu6):
        U = A1 + Bsep[s]
        gU = gsolve(A3, mu, U, lam0=U)
        gA = gsolve(A3, mu, A1, lam0=A1)
        gB = gsolve(A3, mu, Bsep[s], lam0=Bsep[s])
        c = C3s[0]
        rel = (gU[c] - gA[c] - gB[c]) / (gA[c] + gB[c])
        vals[(s, mu)] = rel
        okd &= close(rel, seal_d[(s, mu)], 5e-7)
chk("P1.c", okd, "all 8 R2 relative defects reproduce by linear solve "
    "(e.g. sep3 mu=1/6: %.12f vs sealed -0.038627719073)" % vals[(3, mu6)])

# dipole/monopole and occupancy by linear solve
gp12 = gsolve(A3, mu12, Dp); gm12 = gsolve(A3, mu12, Dm)
gp6 = gsolve(A3, mu6, Dp);  gm6 = gsolve(A3, mu6, Dm)
c5 = ix3[(5, 2, 2)]; c10 = ix3[(10, 2, 2)]
rr = [abs(gp12[c5] - gm12[c5]) / gm12[c5], abs(gp12[c10] - gm12[c10]) / gm12[c10],
      abs(gp6[c5] - gm6[c5]) / gm6[c5],   abs(gp6[c10] - gm6[c10]) / gm6[c10]]
chk("P1.d", close(rr[0], 0.90474560, 1e-6) and close(rr[1], 0.89183281, 1e-6) and
            close(rr[2], 0.58758338, 1e-6) and close(rr[3], 0.51901413, 1e-6),
    "dipole/monopole ratios reproduce: %.8f %.8f %.8f %.8f" % tuple(rr))

seal_k = [0.001044979827, 0.003217558908, 0.007721839900, 0.017143768301]
okk = True
for k in range(1, 5):
    gk = gsolve(A3, mu6, occ[:k])
    okk &= close(gk[C3s[0]], seal_k[k - 1], 1e-8)
    okk &= all(gk[i] > 0 for i in range(len(co3)))
chk("P1.e", okk, "occupancy k=1..4 values reproduce by linear solve; Phi > 0 at every site")

# 2-D discriminator
mu4 = 0.25
A12 = [ix2[(1, 3)]]; B12 = [ix2[(11, 3)]]
disk2 = [ix2[(6, y)] for y in (2, 3, 4)]
gf2 = gsolve(A2, mu4, A12); gd2 = gsolve(A2, mu4, A12, lam0=disk2)
r2ax = gd2[ix2[(10, 3)]] / gf2[ix2[(10, 3)]]
o2 = orders(nb2, A12, 30, lam0=frozenset(disk2))
fo2 = next((w for w in range(31) if cnt(o2, w, B12) != 0), None)
d2p = min(bfsd(nb2, A12, blocked=frozenset(disk2))[b] for b in B12)
d2f = min(bfsd(nb2, A12)[b] for b in B12)
chk("P1.f", close(r2ax, 0.099030, 2e-5) and fo2 == 14 and d2p == 14 and d2f == 10,
    "2-D: axial ratio %.6f (sealed 0.099030), onset %s = punctured BFS %d > free %d"
    % (r2ax, fo2, d2p, d2f))

# ---- P3 truncation honesty, exact big-int (independent code)
def partials(nbrs, src, p, q, Ws):
    """exact rational partial sums of sum mu^w N_w(src->target) at the depths in Ws;
       returns dict W -> Fraction, target = Bm site; integer Horner accumulation."""
    n = len(nbrs)
    y = [0] * n
    for s in src: y[s] += 1
    accs = [0] * n
    for v in range(n): accs[v] = y[v]
    P = 1
    out = {}
    Wmax = max(Ws)
    for w in range(1, Wmax + 1):
        yn = [0] * n
        for v in range(n):
            sacc = 0
            for u in nbrs[v]: sacc += y[u]
            yn[v] = sacc
        y = yn
        P *= p
        for v in range(n): accs[v] = accs[v] * q + y[v] * P
        if w in Ws:
            out[w] = Fraction(accs[Bm[0]], q ** w)
    return out

for (p, q) in ((1, 12), (1, 6)):
    ps = partials(nb3, A1, p, q, (500, 820))
    rem = ps[820] - ps[500]
    K = Fraction(sum(phi3[s] for s in A1), min(phi3))
    x = Fraction(p, q) * B3
    tb = K * x ** 501 / (1 - x)
    ok = rem > 0 and rem <= tb
    def l10(fr):
        return math.log10(fr.numerator) - math.log10(fr.denominator)
    chk("P3.%s" % q, ok,
        "mu=%d/%d: remainder(500->820) = 10^%.1f > 0 and <= tail bound 10^%.1f "
        "(slack ~10^%.1f) -- bounds are real, not display zeros"
        % (p, q, l10(rem) if rem > 0 else float('nan'),
           l10(tb), l10(tb / rem) if rem > 0 else float('nan')))

# ---- P4 fresh-mu attack (rows never used by lane OR its verify)
print("\n  P4 FRESH-MU ATTACK: 3-D mu = 1/9 and 9/50 (9/50 = 0.18 > landmark 1/6);")
print("     2-D mu = 13/50 (> landmark 1/4).  Certificates first.")
for (mu_f, B, tag) in ((Fraction(1, 9), B3, "3D-1/9"), (Fraction(9, 50), B3, "3D-9/50"),
                       (Fraction(13, 50), B2, "2D-13/50")):
    chk("P4.cert-%s" % tag, mu_f * B < 1,
        "mu*B = %.6f < 1 certified on the finite venue" % float(mu_f * B))

for mu in (1.0 / 9.0, 0.18):
    tag = "1/9" if mu < 0.15 else "9/50"
    # R2 defect
    s = 3
    U = A1 + Bsep[s]
    gU = gsolve(A3, mu, U, lam0=U); gA = gsolve(A3, mu, A1, lam0=A1)
    gB = gsolve(A3, mu, Bsep[s], lam0=Bsep[s])
    c = C3s[0]
    dneg = (gU[c] - gA[c] - gB[c]) < 0
    # dipole screening
    gp = gsolve(A3, mu, Dp); gm = gsolve(A3, mu, Dm)
    scr = all(abs(gp[ix3[(x, 2, 2)]] - gm[ix3[(x, 2, 2)]]) < gm[ix3[(x, 2, 2)]]
              for x in range(5, 11))
    # occupancy positivity + additivity (float, 1e-10)
    g4 = gsolve(A3, mu, occ)
    gsum = sum(gsolve(A3, mu, [s0]) for s0 in occ)
    addf = np.max(np.abs(g4 - gsum)) < 1e-10 * np.max(np.abs(g4))
    posf = np.all(g4 > 0)
    # shadow
    gf = gsolve(A3, mu, A1); gd = gsolve(A3, mu, A1, lam0=disk)
    plane = [(y, z) for y in range(5) for z in range(5)]
    Rv = {yz: gd[ix3[(9,) + yz]] / gf[ix3[(9,) + yz]] for yz in plane}
    shad = all(v < 1 for v in Rv.values()) and all(Rv[(2, 2)] <= v for v in Rv.values())
    chk("P4.%s" % tag, dneg and scr and addf and posf and shad,
        "mu=%s: R2 defect<0, dipole screened, occupancy additive+positive, shadow "
        "R<1 axial-min (rel defect %.6f, x=5 ratio %.6f, axial R %.6f)"
        % (tag, (gU[c] - gA[c] - gB[c]) / (gA[c] + gB[c]),
           abs(gp[c5] - gm[c5]) / gm[c5], Rv[(2, 2)]))

# 2-D fresh super-landmark row
mu = 0.26
gf2 = gsolve(A2, mu, A12); gd2 = gsolve(A2, mu, A12, lam0=disk2)
Rv2 = {y: gd2[ix2[(10, y)]] / gf2[ix2[(10, y)]] for y in range(7)}
U2 = A12 + [ix2[(4, 3)]]
gU2 = gsolve(A2, mu, U2, lam0=U2); gA2 = gsolve(A2, mu, A12, lam0=A12)
gB2 = gsolve(A2, mu, [ix2[(4, 3)]], lam0=[ix2[(4, 3)]])
c2 = ix2[(11, 3)]
chk("P4.2D", (gU2[c2] - gA2[c2] - gB2[c2]) < 0 and all(0 < v < 1 for v in Rv2.values())
    and all(Rv2[3] <= v for v in Rv2.values()),
    "2-D mu=13/50 (ABOVE the 1/4 landmark, certified): R2 defect<0, shadow 0<R<1 "
    "axial-min (axial R %.6f) -- classes hold above the borrowed landmark" % Rv2[3])

# mirror zero is mu-independent (P2.c per-order equality) -- restate as covered.
nbad = sum(1 for (_, ok) in CH if not ok)
print("=" * 78)
print("SECOND-LINE PROBE: %d/%d OK%s" % (len(CH) - nbad, len(CH),
      "" if nbad == 0 else " -- FAILURES: " + ", ".join(t for t, ok in CH if not ok)))
print("=" * 78)
sys.exit(1 if nbad else 0)
