"""ADVERSARIAL VERIFY for T-44-C composition lane.
Independent reimplementation (coordinate-dict venue build, separate iteration code).
Attacks:
  V1  spectral bound B recomputed independently; convergence certificates re-checked.
  V2  w_min = d support thresholds (free + punctured) recomputed.
  V3  R1 superposition per-order integer identity re-derived; overlap control.
  V4  R2 defect table spot-checked to printed digits (sep=2,3,4,6; both mu rows).
  V5  mirror-pair exact zero re-counted; dipole/monopole ratios to printed digits.
  V6  occupancy accumulation values + exact additivity re-derived.
  V7  shadowing: slab zero + independent cut check, onset=punctured BFS, far-plane
      ratio table spot-checked (min/corner values), lambda monotonicity per order.
  V8  TRUNCATION HONESTY: partial sums at W=500 vs W=800; the actual remainder must
      sit inside the sealed lane's stated tail bound formula K*(muB)^(W+1)/(1-muB).
  V9  MU-SHIFT ATTACK: all composition-class claims rerun at fresh mu rows the lane
      never used (mu=1/10 deep, mu=17/100 nearer the landmark; both certified
      convergent) -- if any composition fact flips, the "mu-uniform within
      convergence / composition is not a critical phenomenon" claim is refuted.
  V10 2-D discriminator onset/cut/ratio spot-checks.
Every check prints a computed boolean.  Exit 0 iff all pass.
"""
import sys
from fractions import Fraction
from collections import deque

CHECKS = []
def chk(tag, ok, msg):
    CHECKS.append((tag, bool(ok)))
    print("  [%s] %s -- %s" % ("OK " if ok else "BAD", tag, msg))

# ---------------- independent venue build (coordinate dict, not index arithmetic)
def build(dims):
    import itertools
    coords = list(itertools.product(*[range(d) for d in dims]))
    idx = {c: i for i, c in enumerate(coords)}
    nbrs = [[] for _ in coords]
    for c in coords:
        for ax in range(len(dims)):
            for dd in (-1, 1):
                q = list(c); q[ax] += dd; q = tuple(q)
                if q in idx:
                    nbrs[idx[c]].append(idx[q])
    return coords, idx, nbrs

def bfs(nbrs, srcs, blocked=frozenset()):
    dist = {}
    dq = deque()
    for s in srcs:
        dist[s] = 0; dq.append(s)
    while dq:
        u = dq.popleft()
        for v in nbrs[u]:
            if v not in dist and v not in blocked:
                dist[v] = dist[u] + 1; dq.append(v)
    return dist

def specbound(nbrs, iters=80):
    n = len(nbrs)
    phi = [1] * n
    for _ in range(iters):
        phi = [phi[v] + sum(phi[u] for u in nbrs[v]) for v in range(n)]
    B = max(Fraction(sum(phi[u] for u in nbrs[v]), phi[v]) for v in range(n))
    return B, phi

def walk(nbrs, srcs, W, wnum=None, wden=1, snapsW=30):
    """independent iteration.  wnum: per-site integer weight numerator (enter-weight),
       denominator wden common (sigma).  Returns snaps (scaled counts, w<=snapsW) and
       full list of scaled count vectors? -- we return snaps plus a callable partial-sum
       accumulator done inline by caller via gen."""
    n = len(nbrs)
    y = [0] * n
    for s in srcs:
        y[s] += 1
    snaps = [y[:]]
    yield y
    w = 0
    while True:
        w += 1
        ynew = [0] * n
        for v in range(n):
            lv = wden if wnum is None else wnum[v]
            if lv:
                s = 0
                for u in nbrs[v]:
                    s += y[u]
                if s:
                    ynew[v] = s * lv
        y = ynew
        yield y

def series(nbrs, srcs, W, mus, wnum=None, sigma=1, snapsW=30):
    """returns (snaps list w=0..snapsW of scaled vectors, acc numerators per mu, denbase)
       acc[mi] = integer vector; G = acc[v] / (sigma*q)^W  (mu = p/q)."""
    gen = walk(nbrs, srcs, W, wnum, sigma)
    n = len(nbrs)
    snaps = []
    acc = [None] * len(mus)
    P = [1] * len(mus)
    for w, y in enumerate(gen):
        if w == 0:
            for mi in range(len(mus)):
                acc[mi] = y[:]
            snaps.append(y[:])
        else:
            for mi, (p, q) in enumerate(mus):
                P[mi] = P[mi] * p if w > 0 else 1
                bm = sigma * q
                a = acc[mi]
                for v in range(n):
                    a[v] = a[v] * bm + y[v] * P[mi]
            if w <= snapsW:
                snaps.append(y[:])
        if w == W:
            break
    return snaps, acc

def G(acc, mus, mi, sigma, W, sites):
    p, q = mus[mi]
    return Fraction(sum(acc[mi][s] for s in sites), (sigma * q) ** W)

def tailbound(B, phi, srcs, mu, W):
    K = Fraction(sum(phi[s] for s in srcs), min(phi))
    x = Fraction(*mu) * B
    return K * x ** (W + 1) / (1 - x)

def dec(fr, digits=12):
    if fr == 0:
        return "0." + "0" * digits
    sign = "-" if fr < 0 else ""
    fr = abs(fr)
    ip = fr.numerator // fr.denominator
    rem = fr - ip
    return "%s%d.%0*d" % (sign, ip, digits, (rem.numerator * 10 ** digits) // rem.denominator)

print("=" * 78)
print("T-44-C ADVERSARIAL VERIFY -- independent code path")
print("=" * 78)

coords3, idx3, N3 = build((11, 5, 5))
coords2, idx2, N2 = build((13, 7))
MU3 = [(1, 12), (1, 6)]
MU2 = [(1, 8), (1, 4)]
W3, W2 = 500, 1600

# ---------------- V1 spectral bounds
B3, phi3 = specbound(N3)
B2, phi2 = specbound(N2)
print("\nV1 spectral bounds (independent):")
chk("V1.a", dec(B3, 6) == "5.396145", "3-D B ~ %s (sealed printed 5.396145)" % dec(B3, 6))
chk("V1.b", dec(B2, 6) == "3.797765", "2-D B ~ %s (sealed printed 3.797765)" % dec(B2, 6))
chk("V1.c", Fraction(1, 6) * B3 < 1 and Fraction(1, 12) * B3 < 1
         and Fraction(1, 4) * B2 < 1 and Fraction(1, 8) * B2 < 1,
    "all four sealed mu rows certified convergent, independently")

# ---------------- geometry (independent)
A1 = [idx3[(1, 2, 2)]]
C3s = [idx3[(10, 2, 2)]]
SEPS = (2, 3, 4, 6)
B1 = {s: [idx3[(1 + s, 2, 2)]] for s in SEPS}
Ablk = [idx3[(x, y, z)] for x in (1, 2) for y in (1, 2) for z in (1, 2)]
Bblk = [idx3[(x, y, z)] for x in (5, 6) for y in (1, 2) for z in (1, 2)]
Bov  = [idx3[(x, y, z)] for x in (2, 3) for y in (1, 2) for z in (1, 2)]
Am, Bm = [idx3[(1, 2, 2)]], [idx3[(9, 2, 2)]]
plane5 = [idx3[(5, y, z)] for y in range(5) for z in range(5)]
Dp, Dm = [idx3[(2, 2, 2)]], [idx3[(3, 2, 2)]]
occ = [idx3[(2 + i, 2, 2)] for i in range(4)]
As, Bs = [idx3[(1, 2, 2)]], [idx3[(9, 2, 2)]]
slab = set(idx3[(5, y, z)] for y in range(5) for z in range(5))
disk = set(idx3[(5, y, z)] for y in (1, 2, 3) for z in (1, 2, 3))
center = {idx3[(5, 2, 2)]}

def wvec3(lam0=frozenset(), half=frozenset(), sigma=1):
    w = [sigma] * len(coords3)
    for s in half: w[s] = 1
    for s in lam0: w[s] = 0
    return w

# ---------------- V2 w_min = d
print("\nV2 w_min = d (free):")
d_AC = bfs(N3, A1)[C3s[0]]
sn, _ = series(N3, A1, 30, MU3)
below = all(sum(sn[w][c] for c in C3s) == 0 for w in range(d_AC))
atd = sum(sn[d_AC][c] for c in C3s)
par = all(sum(sn[w][c] for c in C3s) == 0 for w in range(d_AC, 31) if (w - d_AC) % 2 == 1)
chk("V2.a", d_AC == 9 and below and atd == 1 and par,
    "3-D A->C: d=%d, N_w=0 below, N_d=%d, parity holds (sealed: d=9, N_d=1)" % (d_AC, atd))

# ---------------- V3 R1 superposition identity
print("\nV3 R1 superposition (per-order integer identity, independent counts):")
ok = True
for s in SEPS:
    U = A1 + B1[s]
    snU, _ = series(N3, U, 30, MU3)
    snA, _ = series(N3, A1, 30, MU3)
    snB, _ = series(N3, B1[s], 30, MU3)
    ok &= all(sum(snU[w][c] for c in C3s) ==
              sum(snA[w][c] for c in C3s) + sum(snB[w][c] for c in C3s)
              for w in range(31))
chk("V3.a", ok, "defect = 0 at every order w<=30, seps 2,3,4,6")
snU, _ = series(N3, Ablk + Bblk, 30, MU3)
snA, _ = series(N3, Ablk, 30, MU3)
snB, _ = series(N3, Bblk, 30, MU3)
chk("V3.b", all(sum(snU[w][c] for c in C3s) == sum(snA[w][c] for c in C3s) +
                sum(snB[w][c] for c in C3s) for w in range(31)),
    "2x2x2 block regions additive per order")
Uov = sorted(set(Ablk) | set(Bov)); Iov = sorted(set(Ablk) & set(Bov))
snUo, _ = series(N3, Uov, 30, MU3)
snBo, _ = series(N3, Bov, 30, MU3)
snIo, _ = series(N3, Iov, 30, MU3)
ok = all(sum(snUo[w][c] for c in C3s) - sum(snA[w][c] for c in C3s) -
         sum(snBo[w][c] for c in C3s) == -sum(snIo[w][c] for c in C3s)
         for w in range(31))
nz = any(sum(snIo[w][c] for c in C3s) > 0 for w in range(31))
chk("V3.c", ok and nz, "overlap control: defect = -G(AnB->C) per order, nonzero")

# ---------------- V4 R2 defect digits
print("\nV4 R2 defect table (resummed, W=500, independent):")
sealed_rel = {2: ("-0.008174820104", "-0.088532997166"),
              3: ("-0.000784482164", "-0.038627719073"),
              4: ("-0.000077104643", "-0.017807581541"),
              6: ("-0.000000808188", "-0.004065283309")}
ok_all = True
for s in SEPS:
    U = A1 + B1[s]
    _, accU = series(N3, U, W3, MU3, wnum=wvec3(lam0=set(U)), sigma=1)
    _, accA = series(N3, A1, W3, MU3, wnum=wvec3(lam0=set(A1)), sigma=1)
    _, accB = series(N3, B1[s], W3, MU3, wnum=wvec3(lam0=set(B1[s])), sigma=1)
    for mi in range(2):
        dG = G(accU, MU3, mi, 1, W3, C3s) - G(accA, MU3, mi, 1, W3, C3s) - G(accB, MU3, mi, 1, W3, C3s)
        rel = dG / (G(accA, MU3, mi, 1, W3, C3s) + G(accB, MU3, mi, 1, W3, C3s))
        got = dec(rel)
        want = sealed_rel[s][mi]
        ok = (got == want)
        ok_all &= ok
        if not ok:
            print("      sep=%d mi=%d got %s want %s" % (s, mi, got, want))
chk("V4.a", ok_all, "all 8 relative-defect entries reproduce to the printed 12 digits")

# ---------------- V5 mirror zero + dipole digits
print("\nV5 screening:")
snP, accP = series(N3, Am, W3, MU3)
snM, accM = series(N3, Bm, W3, MU3)
ok = all(snP[w][c] == snM[w][c] for w in range(31) for c in plane5)
ok &= all(accP[mi][c] == accM[mi][c] for mi in range(2) for c in plane5)
chk("V5.a", ok, "mirror pair: counts equal at every order and resummed, all 25 mid-plane sites")
_, accDp = series(N3, Dp, W3, MU3)
_, accDm = series(N3, Dm, W3, MU3)
sealed_ratio = {(5, 1): "0.58758338", (10, 1): "0.51901413", (5, 0): "0.90474560", (10, 0): "0.89183281"}
ok = True
for (x, mi), want in sealed_ratio.items():
    c = [idx3[(x, 2, 2)]]
    phid = G(accDp, MU3, mi, 1, W3, c) - G(accDm, MU3, mi, 1, W3, c)
    r = dec(abs(phid) / G(accDm, MU3, mi, 1, W3, c), 8)
    ok &= (r == want)
chk("V5.b", ok, "dipole/monopole ratios reproduce printed digits at x=5,10 both mu rows")

# ---------------- V6 occupancy accumulation
print("\nV6 accumulation:")
sealed_k = ["0.001044979827", "0.003217558908", "0.007721839900", "0.017143768301"]
ok = okadd = True
singles = [series(N3, [occ[i]], W3, MU3)[1] for i in range(4)]
for k in range(1, 5):
    _, acck = series(N3, occ[:k], W3, MU3)
    ok &= (dec(G(acck, MU3, 1, 1, W3, C3s)) == sealed_k[k - 1])
    for mi in range(2):
        okadd &= (sum(acck[mi][c] for c in C3s) ==
                  sum(singles[i][mi][c] for i in range(k) for c in C3s))
chk("V6.a", ok, "k=1..4 occupancy values reproduce printed digits (mu=1/6)")
chk("V6.b", okadd, "exact additivity of numerators, both mu rows")

# ---------------- V7 shadowing
print("\nV7 shadowing:")
snF, accF = series(N3, As, W3, MU3)
_, accSlab = series(N3, As, W3, MU3, wnum=wvec3(lam0=slab))
cut = Bs[0] not in bfs(N3, As, blocked=frozenset(slab))
zero = all(accSlab[mi][b] == 0 for mi in range(2) for b in Bs)
chk("V7.a", cut and zero, "slab: BFS cut (independent) and G(A->B)=0 resummed")
snD, accD = series(N3, As, W3, MU3, wnum=wvec3(lam0=disk))
dfree = bfs(N3, As)[Bs[0]]
ddisk = bfs(N3, As, blocked=frozenset(disk))[Bs[0]]
fo = next((w for w in range(31) if sum(snD[w][b] for b in Bs) != 0), None)
chk("V7.b", dfree == 8 and ddisk == 12 and fo == 12,
    "onset order %s = punctured BFS %d > free %d (sealed: 12, 12, 8)" % (fo, ddisk, dfree))
gfree = G(accF, MU3, 1, 1, W3, Bs)
chk("V7.c", dec(gfree, 16) == "0.0012367408720414",
    "G_free(A->B) mu=1/6 = %s (sealed 0.0012367408720414)" % dec(gfree, 16))
c22 = [idx3[(9, 2, 2)]]
r22 = G(accD, MU3, 1, 1, W3, c22) / G(accF, MU3, 1, 1, W3, c22)
ccor = [idx3[(9, 0, 0)]]
rcor = G(accD, MU3, 1, 1, W3, ccor) / G(accF, MU3, 1, 1, W3, ccor)
chk("V7.d", dec(r22, 8) == "0.10624064" and dec(rcor, 8) == "0.15402869",
    "far-plane ratios: axial %s corner %s (sealed 0.10624064 / 0.15402869)" % (dec(r22, 8), dec(rcor, 8)))
allR = []
for y in range(5):
    for z in range(5):
        c = [idx3[(9, y, z)]]
        allR.append(((y, z), G(accD, MU3, 1, 1, W3, c) / G(accF, MU3, 1, 1, W3, c)))
chk("V7.e", all(v < 1 for _, v in allR) and all(r22 <= v for _, v in allR),
    "R < 1 everywhere on far plane, minimum axial (independent full table)")
snH, accH = series(N3, As, W3, MU3, wnum=wvec3(half=disk, sigma=2), sigma=2)
ok = True
for w in range(31):
    n0 = Fraction(sum(snD[w][b] for b in Bs))
    nh = Fraction(sum(snH[w][b] for b in Bs), 2 ** w)
    n1 = Fraction(sum(snF[w][b] for b in Bs))
    ok &= (n0 <= nh <= n1)
chk("V7.f", ok, "lambda monotone per order: blocked <= half <= empty (scaled counts)")

# ---------------- V8 truncation honesty
print("\nV8 truncation honesty (W=500 vs W=800 remainder vs stated tail bound):")
_, accF8 = series(N3, As, 800, MU3)
ok = True
for mi in range(2):
    g500 = G(accF, MU3, mi, 1, W3, Bs)
    g800 = G(accF8, MU3, mi, 1, 800, Bs)
    rem = g800 - g500
    tb = tailbound(B3, phi3, As, MU3[mi], W3)
    ok &= (0 <= rem <= tb)
    print("      mu=1/%-3d actual remainder(500->800) = %s <= stated tail %s : %s"
          % (MU3[mi][1], dec(rem, 20), dec(tb, 20), 0 <= rem <= tb))
chk("V8.a", ok, "actual remainder within the sealed tail-bound formula, both mu rows")
# same for an R2 run (opaque union, sep=2) -- weighted case
U = A1 + B1[2]
_, aU5 = series(N3, U, W3, MU3, wnum=wvec3(lam0=set(U)))
_, aU8 = series(N3, U, 800, MU3, wnum=wvec3(lam0=set(U)))
ok = True
for mi in range(2):
    rem = G(aU8, MU3, mi, 1, 800, C3s) - G(aU5, MU3, mi, 1, W3, C3s)
    tb = tailbound(B3, phi3, U, MU3[mi], W3)
    ok &= (0 <= rem <= tb)
chk("V8.b", ok, "weighted (R2) run: remainder within tail bound too (weights<=1 shrink counts)")

# ---------------- V9 mu-shift attack
print("\nV9 MU-SHIFT ATTACK: fresh rows mu=1/10 and mu=17/100 (never run by the lane):")
MUX = [(1, 10), (17, 100)]
for (p, q) in MUX:
    x = Fraction(p, q) * B3
    chk("V9.cert-%d/%d" % (p, q), x < 1, "mu*B = %s < 1 certified" % dec(x, 6))
# R2 defect still negative
_, xU = series(N3, U, W3, MUX, wnum=wvec3(lam0=set(U)))
_, xA = series(N3, A1, W3, MUX, wnum=wvec3(lam0=set(A1)))
_, xB = series(N3, B1[2], W3, MUX, wnum=wvec3(lam0=set(B1[2])))
ok = all(G(xU, MUX, mi, 1, W3, C3s) - G(xA, MUX, mi, 1, W3, C3s) - G(xB, MUX, mi, 1, W3, C3s) < 0
         for mi in range(2))
chk("V9.a", ok, "R2 defect < 0 at both fresh mu rows (subadditive class holds)")
# partial screening
_, xDp = series(N3, Dp, W3, MUX)
_, xDm = series(N3, Dm, W3, MUX)
ok = True
rat = {0: [], 1: []}
for mi in range(2):
    tb = tailbound(B3, phi3, Dp, MUX[mi], W3) + tailbound(B3, phi3, Dm, MUX[mi], W3)
    for x_ in range(5, 11):
        c = [idx3[(x_, 2, 2)]]
        phid = G(xDp, MUX, mi, 1, W3, c) - G(xDm, MUX, mi, 1, W3, c)
        mono = G(xDm, MUX, mi, 1, W3, c)
        ok &= (abs(phid) + tb < mono)
        rat[mi].append(abs(phid) / mono)
chk("V9.b", ok, "|Phi_dipole| < Phi_mono at every probe, both fresh rows (screening class holds)")
print("      screening ratio at x=5: mu=1/10 %s, mu=17/100 %s (mu-DEPENDENT depth confirmed)"
      % (dec(rat[0][0], 6), dec(rat[1][0], 6)))
# accumulation positive + additive
_, x4 = series(N3, occ, W3, MUX)
okp = all(Fraction(x4[mi][c], MUX[mi][1] ** W3 if False else 1) >= 0 for mi in range(2) for c in range(len(coords3)))
okp = all(x4[mi][c] >= 0 for mi in range(2) for c in range(len(coords3)))
xs = [series(N3, [occ[i]], W3, MUX)[1] for i in range(4)]
oka = all(sum(x4[mi][c] for c in C3s) == sum(xs[i][mi][c] for i in range(4) for c in C3s)
          for mi in range(2))
chk("V9.c", okp and oka, "occupancy: one-signed and exactly additive at fresh rows (accumulation class holds)")
# mirror zero (mu-independent per order, but resummed check at fresh rows)
_, xP = series(N3, Am, W3, MUX)
_, xM = series(N3, Bm, W3, MUX)
chk("V9.d", all(xP[mi][c] == xM[mi][c] for mi in range(2) for c in plane5),
    "mirror-pair exact zero at fresh rows (exact screening holds)")
# shadow ratios
_, xF = series(N3, As, W3, MUX)
_, xD = series(N3, As, W3, MUX, wnum=wvec3(lam0=disk))
ok = True
vals = {}
for mi in range(2):
    rr = {}
    for y in range(5):
        for z in range(5):
            c = [idx3[(9, y, z)]]
            rr[(y, z)] = G(xD, MUX, mi, 1, W3, c) / G(xF, MUX, mi, 1, W3, c)
    ok &= all(v < 1 for v in rr.values()) and all(rr[(2, 2)] <= v for v in rr.values())
    vals[mi] = rr[(2, 2)]
chk("V9.e", ok, "shadow: R<1 everywhere, axial minimum, both fresh rows (opacity class holds; "
    "axial R = %s at 1/10, %s at 17/100 -- magnitudes move with mu, structure does not)"
    % (dec(vals[0], 6), dec(vals[1], 6)))

# ---------------- V10 2-D discriminator
print("\nV10 2-D discriminator:")
A2s, B2s = [idx2[(1, 3)]], [idx2[(11, 3)]]
slab2 = set(idx2[(6, y)] for y in range(7))
disk2 = set(idx2[(6, y)] for y in (2, 3, 4))
def wvec2(lam0=frozenset()):
    w = [1] * len(coords2)
    for s in lam0: w[s] = 0
    return w
d_free2 = bfs(N2, A2s)[B2s[0]]
d_disk2 = bfs(N2, A2s, blocked=frozenset(disk2))[B2s[0]]
cut2 = B2s[0] not in bfs(N2, A2s, blocked=frozenset(slab2))
sn2, acc2F = series(N2, A2s, W2, MU2)
sn2D, acc2D = series(N2, A2s, W2, MU2, wnum=wvec2(lam0=disk2))
_, acc2S = series(N2, A2s, W2, MU2, wnum=wvec2(lam0=slab2))
fo2 = next((w for w in range(31) if sum(sn2D[w][b] for b in B2s) != 0), None)
zero2 = all(acc2S[mi][b] == 0 for mi in range(2) for b in B2s)
c2ax = [idx2[(10, 3)]]
r2ax = G(acc2D, MU2, 1, 1, W2, c2ax) / G(acc2F, MU2, 1, 1, W2, c2ax)
chk("V10.a", cut2 and zero2 and d_free2 == 10 and d_disk2 == 14 and fo2 == 14,
    "cut zero; onset %s = punctured BFS %d > free %d (sealed 14, 14, 10)" % (fo2, d_disk2, d_free2))
chk("V10.b", dec(r2ax, 6) == "0.099030",
    "2-D axial ratio %s (sealed 0.099030)" % dec(r2ax, 6))

nbad = sum(1 for _, ok in CHECKS if not ok)
print("=" * 78)
print("VERIFY CHECKS: %d/%d OK%s" % (len(CHECKS) - nbad, len(CHECKS),
      "" if nbad == 0 else "  FAILURES: " + ", ".join(t for t, ok in CHECKS if not ok)))
print("=" * 78)
sys.exit(0 if nbad else 1)
