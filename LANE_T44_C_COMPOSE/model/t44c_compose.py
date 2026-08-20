"""T-44-C  THREE-REGION COMPOSITION OF THE MEDIATED COUPLING -- exact, the venue's own lattice.

THE QUESTION (assignment): with the T-44 weighted instrument -- the mediated coupling
G(S->T) = sum over admissible strings from S to T, per-link amplitude mu (the coupling
tier's ONE declared parameter), a string of weight w contributing mu^w -- does the coupling
COMPOSE like a source should?
  (1) SUPERPOSITION   G(A u B -> C)  vs  G(A->C) + G(B->C), disjoint A,B: exact defect.
  (2) SCREENING vs ACCUMULATION with record VALUES on A and B under the C-72 encoding split
      (occupancy vs orientation): does the encoding-level sign structure survive mediation?
  (3) SHADOWING: an intervening region D between A and B -- blocked or reweighted paths --
      the surface's own opacity, computed not asserted.  Controls: empty D; fully-blocking D.
This lane feeds the judge COMPOSITION FACTS; it does not locate criticality (that is the
exponent lane's assignment) and earns no exponent.

THE INSTRUMENT (declared once, used for all three probes):
  Venue lattices (the venue's own, not imported):
    WORLD tier (earned dimension 3, C-78): the GR1/T43-B grain block -- grains at integer
      points, adjacency = shares a face; here an 11 x 5 x 5 open block.
    CORNER-CLASS discriminator (earned dimension 2): the 2-D lattice the corner tier's
      coupling-writer strings live on; here a 13 x 7 open block.  (D-15: a venue where the
      walk class is KNOWN to differ, D=2 vs D=3.)
  ADMISSIBLE STRING from S to T: a walk on the adjacency graph, start in S, end in T,
    weight w = number of links.  C-80's one genuine separation law, w_min = d exactly,
    appears in the instrument as: the count N_w(S->T) is 0 for all w < d(S,T) and positive
    at w = d(S,T) -- GATED below, not assumed.
  Per-order counts N_w(S->T): exact integers (exact scaled integers when a per-visit weight
    lambda is declared), by sparse integer transfer-matrix iteration.
  Resummed G_mu(S->T) = sum_w mu^w N_w(S->T): exact rational partial sum to W_res plus an
    exact rational TAIL BOUND from a Collatz-Wielandt certificate  A phi <= B phi  (phi > 0
    an integer power-iteration vector), giving N_w <= K_S B^w and
    tail(W) <= K_S (mu B)^{W+1} / (1 - mu B), all exact rationals.  Every resummed
    comparison below is gated on INTERVALS [partial, partial + tail]; every equality is a
    per-order integer identity, exact with no truncation.
  ADMISSIBILITY RULES for source regions (both declared; the split is a finding, not a bug):
    R1 ENDPOINT RULE (transparent sources): any walk with start in S, end in T.
    R2 OPAQUE-SOURCE RULE: additionally, the walk's interior may not visit any PRESENT
      source site (records are matter to the mediator).  Implemented exactly: per-visit
      weight 0 on all present source sites (start unweighted; entering weighted).

mu ROWS (declared): world 3-D  mu in {1/12, 1/6};  2-D  mu in {1/8, 1/4}.  The values 1/6
  and 1/4 are the infinite-lattice walk-critical amplitudes 1/(2D) -- BORROWED landmarks
  (Polya 1921; lattice Green's functions: Watson 1939, Montroll-Weiss 1965, Spitzer 1964)
  used ONLY as declared benchmark rows; on these finite venues both rows are certified
  convergent by the computed spectral bound (G0.3).  Locating criticality is the exponent
  lane's job, not exercised here.

OWNERS (D-24, everything standard named): walk generating functions / transfer matrices
  (Feller; Stanley, Enumerative Combinatorics); lattice Green's functions and recurrence
  (Polya; Watson; Montroll-Weiss; Spitzer, potential theory); Perron-Frobenius /
  Collatz-Wielandt bound (Collatz 1942, Wielandt 1950; Varga, Matrix Iterative Analysis);
  kernel positivity as sign-preservation is the positivity half of potential theory's
  maximum principle territory (Doob); Menger/BFS separation for the blocking control
  (Menger 1927).  OURS: the composition facts of THIS instrument on THIS venue -- the
  R1-exact / R2-defect superposition split, the C-72 sign-survival result, the computed
  opacity taxonomy, and the earned-separation-lengthening reading of shadow onset.

D-1: no gravitational form in any construction step; Newton is named ONLY in the final
  comparison section.  D-15: a control beside every zero and every flat line.  D-24 audit
  in D24_AUDIT.txt.  No literal verdicts: every PASS/FAIL prints from a computed boolean.
"""
import sys
from fractions import Fraction
from collections import deque

# ---------------------------------------------------------------- venue: open grid block
class Venue:
    def __init__(self, dims, name):
        self.dims = dims
        self.name = name
        self.D = len(dims)
        n = 1
        for d in dims:
            n *= d
        self.n = n
        self.nbrs = [[] for _ in range(n)]
        for site in self._sites():
            i = self.idx(site)
            for ax in range(self.D):
                for dd in (-1, 1):
                    q = list(site)
                    q[ax] += dd
                    if 0 <= q[ax] < dims[ax]:
                        self.nbrs[i].append(self.idx(tuple(q)))

    def _sites(self):
        def rec(prefix, ax):
            if ax == self.D:
                yield tuple(prefix)
                return
            for v in range(self.dims[ax]):
                yield from rec(prefix + [v], ax + 1)
        yield from rec([], 0)

    def idx(self, site):
        i = 0
        for ax in range(self.D):
            i = i * self.dims[ax] + site[ax]
        return i

    def bfs_dist(self, srcs, blocked=frozenset()):
        """exact hop distance from the SET srcs to every site; blocked sites impassable
           (a blocked site gets no distance unless it is a source, which we disallow)."""
        INF = -1
        dist = [INF] * self.n
        dq = deque()
        for s in srcs:
            dist[s] = 0
            dq.append(s)
        while dq:
            u = dq.popleft()
            for v in self.nbrs[u]:
                if dist[v] == INF and v not in blocked:
                    dist[v] = dist[u] + 1
                    dq.append(v)
        return dist

    def spectral_bound(self, iters=80):
        """Collatz-Wielandt upper bound on the adjacency spectral radius:
           for any phi > 0,  lambda_max <= max_v (A phi)_v / phi_v.  phi = (A+I)^k 1,
           exact integers; the bound is an exact rational.  BORROWED (Collatz, Wielandt)."""
        phi = [1] * self.n
        for _ in range(iters):
            phi = [phi[v] + sum(phi[u] for u in self.nbrs[v]) for v in range(self.n)]
        B = Fraction(0)
        for v in range(self.n):
            r = Fraction(sum(phi[u] for u in self.nbrs[v]), phi[v])
            if r > B:
                B = r
        self.phi = phi
        self.phimin = min(phi)
        self.B = B
        return B

    def K_of(self, srcs):
        """N_w(S->T) <= K_S B^w for every target set T (proof in header: A phi <= B phi,
           1_T <= phi/phi_min)."""
        return Fraction(sum(self.phi[s] for s in srcs), self.phimin)

# ---------------------------------------------------------------- the weighted string sum
W_ORD = 30  # per-order (exact-integer) table depth

def run(V, src_sites, W_res, mus, lam0=frozenset(), lamhalf=frozenset(),
        force_sigma=None):
    """Transfer-matrix iteration of the weighted string sum from the region src_sites.
       Per-visit weight: 0 on lam0 sites, 1/2 on lamhalf sites, 1 elsewhere; weights apply
       on ENTERING a site (the start endpoint is unweighted).  sigma = 2 iff lamhalf used
       (force_sigma exercises the scaled code path with all-unit weights: the empty-D
       control's independent route); scaled counts ctilde_w = sigma^w * N_w are integers.
       Returns dict with:
         snaps  : list of scaled count vectors, w = 0..W_ORD (exact)
         acc    : per mu row, per site, integer numerator of the partial sum over the
                  common denominator (sigma*q)^W_res  (mu = p/q)
         sigma, W_res, mus, K (tail prefactor for this source set)."""
    n = V.n
    sigma = force_sigma if force_sigma else (2 if lamhalf else 1)
    wvec = [sigma] * n
    for s in lamhalf:
        wvec[s] = 1
    for s in lam0:
        wvec[s] = 0
    y = [0] * n
    for s in src_sites:
        y[s] += 1
    snaps = [y[:]]
    nm = len(mus)
    acc = [y[:] for _ in range(nm)]
    P = [1] * nm
    base = [sigma * q for (p, q) in mus]
    for w in range(1, W_res + 1):
        ynew = [0] * n
        for v in range(n):
            lv = wvec[v]
            if lv:
                s = 0
                for u in V.nbrs[v]:
                    s += y[u]
                if s:
                    ynew[v] = s * lv
        y = ynew
        if w <= W_ORD:
            snaps.append(y[:])
        for mi in range(nm):
            p = mus[mi][0]
            P[mi] *= p
            bm = base[mi]
            a = acc[mi]
            for v in range(n):
                a[v] = a[v] * bm + y[v] * P[mi]
    return dict(snaps=snaps, acc=acc, sigma=sigma, W=W_res, mus=mus,
                K=V.K_of(src_sites), V=V)

def gval(R, mi, sites):
    """exact rational partial sum of G_mu(src -> sites)."""
    p, q = R['mus'][mi]
    den = (R['sigma'] * q) ** R['W']
    num = sum(R['acc'][mi][s] for s in sites)
    return Fraction(num, den)

def gtail(R, mi):
    """exact rational tail bound for ONE target site; for a target set multiply by |T|."""
    p, q = R['mus'][mi]
    V = R['V']
    x = Fraction(p, q) * V.B  # per-step bound (weighted matrix <= sigma*A entrywise,
                              # true weight <= mu*B per step regardless of sigma)
    assert x < 1
    return R['K'] * x ** (R['W'] + 1) / (1 - x)

def ordcount(R, w, sites):
    """exact per-order scaled count ctilde_w summed over target sites (integer)."""
    return sum(R['snaps'][w][s] for s in sites)

def dec(fr, digits=12):
    """exact-rational -> fixed decimal string (truncated, exact arithmetic)."""
    if fr == 0:
        return "0." + "0" * digits
    sign = "-" if fr < 0 else ""
    fr = abs(fr)
    ip = fr.numerator // fr.denominator
    rem = fr - ip
    frac_digits = (rem.numerator * 10 ** digits) // rem.denominator
    return "%s%d.%0*d" % (sign, ip, digits, frac_digits)

GATES = []
def gate(tag, ok, msg):
    GATES.append((tag, bool(ok)))
    print("  [%s] %s -- %s" % ("PASS" if ok else "FAIL", tag, msg))

# ================================================================ build venues
print("=" * 78)
print("T-44-C  THREE-REGION COMPOSITION -- exact weighted string sums, both tiers")
print("=" * 78)

V3 = Venue((11, 5, 5), "WORLD 3-D 11x5x5 grain block (face adjacency)")
V2 = Venue((13, 7), "2-D 13x7 block (corner-class string lattice; D-15 discriminator)")
MU3 = [(1, 12), (1, 6)]   # declared; 1/6 = 1/(2*3) infinite-lattice landmark (BORROWED)
MU2 = [(1, 8), (1, 4)]    # declared; 1/4 = 1/(2*2) infinite-lattice landmark (BORROWED)
WRES3, WRES2 = 500, 1600

print("\nSECTION 0 -- venue and instrument self-checks")
print("-" * 78)
for V in (V3, V2):
    print("venue: %s   sites=%d  max degree=%d" %
          (V.name, V.n, max(len(nb) for nb in V.nbrs)))
    sym = all(v in V.nbrs[u] for u in range(V.n) for v in V.nbrs[u])
    gate("G0.1-%dD" % V.D, sym, "adjacency symmetric (undirected venue)")
    B = V.spectral_bound()
    print("    Collatz-Wielandt bound on spectral radius: B = %s ~ %s" %
          (str(B)[:40], dec(B, 6)))
for V, MUS, W in ((V3, MU3, WRES3), (V2, MU2, WRES2)):
    for (p, q) in MUS:
        x = Fraction(p, q) * V.B
        gate("G0.3-%dD-mu%d/%d" % (V.D, p, q), x < 1,
             "mu*B = %s < 1: resummation certified convergent on this venue (W=%d)"
             % (dec(x, 6), W))

# geometry -- 3-D
def S3(*sites):
    return frozenset(V3.idx(s) for s in sites)
A1 = S3((1, 2, 2))
C3 = S3((10, 2, 2))
SEPS = (2, 3, 4, 6)
B1 = {s: S3((1 + s, 2, 2)) for s in SEPS}
Ablk = frozenset(V3.idx((x, y, z)) for x in (1, 2) for y in (1, 2) for z in (1, 2))
Bblk = frozenset(V3.idx((x, y, z)) for x in (5, 6) for y in (1, 2) for z in (1, 2))
Bov  = frozenset(V3.idx((x, y, z)) for x in (2, 3) for y in (1, 2) for z in (1, 2))
Am, Bm = S3((1, 2, 2)), S3((9, 2, 2))
plane5 = frozenset(V3.idx((5, y, z)) for y in range(5) for z in range(5))
Dp, Dm = S3((2, 2, 2)), S3((3, 2, 2))
occ_sites = [V3.idx((2 + i, 2, 2)) for i in range(4)]
As, Bs = S3((1, 2, 2)), S3((9, 2, 2))
slab  = frozenset(V3.idx((5, y, z)) for y in range(5) for z in range(5))
disk  = frozenset(V3.idx((5, y, z)) for y in (1, 2, 3) for z in (1, 2, 3))
center = S3((5, 2, 2))
ring  = slab - disk
slab1 = slab - center
plane9 = [(y, z) for y in range(5) for z in range(5)]

# geometry -- 2-D
def S2(*sites):
    return frozenset(V2.idx(s) for s in sites)
A1_2 = S2((1, 3))
C2 = S2((11, 3))
B1_2 = {s: S2((1 + s, 3)) for s in SEPS}
Am2, Bm2 = S2((1, 3)), S2((11, 3))
plane6_2 = frozenset(V2.idx((6, y)) for y in range(7))
slab2 = frozenset(V2.idx((6, y)) for y in range(7))
disk2 = frozenset(V2.idx((6, y)) for y in (2, 3, 4))
center2 = S2((6, 3))
ring2 = slab2 - disk2
plane10_2 = [y for y in range(7)]

# w_min = d contact with C-80 (per pair used; parity too -- bipartite venue)
print("\n  w_min = d check (the C-80 growing law as the mediator's support threshold):")
runs_cache = {}
def RUN(V, src, W, mus, lam0=frozenset(), lamhalf=frozenset(), force_sigma=None):
    key = (id(V), tuple(sorted(src)), W, tuple(mus),
           tuple(sorted(lam0)), tuple(sorted(lamhalf)), force_sigma)
    if key not in runs_cache:
        runs_cache[key] = run(V, sorted(src), W, mus, lam0, lamhalf, force_sigma)
    return runs_cache[key]

r_A1 = RUN(V3, A1, WRES3, MU3)
r_A12 = RUN(V2, A1_2, WRES2, MU2)
for (V, r, src, tgt, lbl) in ((V3, r_A1, A1, C3, "3D A->C"),
                              (V2, r_A12, A1_2, C2, "2D A->C")):
    d = min(V.bfs_dist(src)[t] for t in tgt)
    below = all(ordcount(r, w, tgt) == 0 for w in range(d))
    at = ordcount(r, d, tgt) > 0
    par = all(ordcount(r, w, tgt) == 0
              for w in range(d, W_ORD + 1) if (w - d) % 2 == 1)
    gate("G0.2-%s" % lbl.replace(" ", ""), below and at and par,
         "%s: N_w=0 for w<d=%d, N_d=%d>0, odd-parity orders vanish (bipartite)"
         % (lbl, d, ordcount(r, d, tgt)))

# kernel positivity: G_mu(a->c) > 0 for EVERY site c (both mu rows, both venues)
for (V, r, lbl, W) in ((V3, r_A1, "3D", WRES3), (V2, r_A12, "2D", WRES2)):
    okall = True
    for mi in range(2):
        t = gtail(r, mi)
        p, q = r['mus'][mi]
        den = (r['sigma'] * q) ** r['W']
        for c in range(V.n):
            if Fraction(r['acc'][mi][c], den) <= t:
                okall = False
    gate("G0.4-%s" % lbl, okall,
         "kernel positivity: G(a->c) > tail bound > 0 at every site, both mu rows "
         "(mediation carries no sign of its own)")

# ================================================================ SECTION 1
print("\nSECTION 1 -- SUPERPOSITION: G(AuB->C) vs G(A->C)+G(B->C), exact defect")
print("-" * 78)
print("  R1 ENDPOINT RULE (transparent sources):")
ok_cert = ok_ord = ok_res = True
for (V, W, mus, Aset, Bmap, Cset, lbl) in (
        (V3, WRES3, MU3, A1, B1, C3, "3D"),
        (V2, WRES2, MU2, A1_2, B1_2, C2, "2D")):
    for s in SEPS:
        Bset = Bmap[s]
        U = Aset | Bset
        # indicator-additivity certificate: 1_AuB - 1_A - 1_B = 0 (disjoint)
        cert = (len(Aset & Bset) == 0) and (len(U) == len(Aset) + len(Bset))
        rA, rB, rU = (RUN(V, Aset, W, mus), RUN(V, Bset, W, mus), RUN(V, U, W, mus))
        po = all(ordcount(rU, w, Cset) ==
                 ordcount(rA, w, Cset) + ordcount(rB, w, Cset)
                 for w in range(W_ORD + 1))
        rs = all(sum(rU['acc'][mi][c] for c in Cset) ==
                 sum(rA['acc'][mi][c] for c in Cset) +
                 sum(rB['acc'][mi][c] for c in Cset) for mi in range(2))
        ok_cert &= cert; ok_ord &= po; ok_res &= rs
gate("G1.1a", ok_cert, "disjointness certificate holds at every separation, both venues")
gate("G1.1b", ok_ord, "per-order defect = 0 at EVERY order w<=%d, every separation, "
     "both venues (integer identity)" % W_ORD)
gate("G1.1c", ok_res, "resummed partial sums exactly additive at all mu rows "
     "(equal numerators over the common denominator)")

# block-region row (3-D)
rAb, rBb, rUb = (RUN(V3, Ablk, WRES3, MU3), RUN(V3, Bblk, WRES3, MU3),
                 RUN(V3, Ablk | Bblk, WRES3, MU3))
pob = all(ordcount(rUb, w, C3) == ordcount(rAb, w, C3) + ordcount(rBb, w, C3)
          for w in range(W_ORD + 1))
gate("G1.1d", pob, "2x2x2 BLOCK regions: same exact additivity, every order")

# D-15 control beside the zero defect: overlapping regions register the inclusion-
# exclusion defect exactly (the zero above measures disjointness, not blindness)
rBo = RUN(V3, Bov, WRES3, MU3)
rUo = RUN(V3, Ablk | Bov, WRES3, MU3)
rIo = RUN(V3, Ablk & Bov, WRES3, MU3)
poo = all(ordcount(rUo, w, C3) - ordcount(rAb, w, C3) - ordcount(rBo, w, C3)
          == -ordcount(rIo, w, C3) for w in range(W_ORD + 1))
some_neg = any(ordcount(rIo, w, C3) > 0 for w in range(W_ORD + 1))
gate("G1.2", poo and some_neg,
     "OVERLAP control: defect = -G(AnB->C) exactly at every order, and it is nonzero "
     "(|AnB|=%d sites) -- the R1 zero is a measurement" % len(Ablk & Bov))

print("\n  R2 OPAQUE-SOURCE RULE (records are matter to the mediator):")
print("  defect(sep) = G_R2(AuB->C) - G_R2(A->C alone) - G_R2(B->C alone), exact")
ok_sign = True
rows = []
for s in SEPS:
    Bset = B1[s]
    U = A1 | Bset
    rUa = RUN(V3, U, WRES3, MU3, lam0=U)
    rAa = RUN(V3, A1, WRES3, MU3, lam0=A1)
    rBa = RUN(V3, Bset, WRES3, MU3, lam0=Bset)
    ords = [ordcount(rUa, w, C3) - ordcount(rAa, w, C3) - ordcount(rBa, w, C3)
            for w in range(W_ORD + 1)]
    ok_sign &= all(d <= 0 for d in ords) and any(d < 0 for d in ords)
    dvals = []
    for mi in range(2):
        dG = gval(rUa, mi, C3) - gval(rAa, mi, C3) - gval(rBa, mi, C3)
        Gm = gval(rAa, mi, C3) + gval(rBa, mi, C3)
        tb = gtail(rUa, mi) + gtail(rAa, mi) + gtail(rBa, mi)
        dvals.append((dG, Gm, tb))
    rows.append((s, ords, dvals))
gate("G1.3a", ok_sign, "R2 per-order defect <= 0 at every order and < 0 somewhere, "
     "every separation (self-shadowing: walks through the other source are lost)")
print("    sep | first nonzero defect order | defect/G_sum at mu=1/12 | at mu=1/6 "
      "(exact rationals; |tail| bounds printed)")
mags = []
for (s, ords, dvals) in rows:
    fo = next((w for w, d in enumerate(ords) if d != 0), None)
    rels = []
    for (dG, Gm, tb) in dvals:
        rels.append(dG / Gm)
    mags.append((s, rels))
    print("    %3d | %3s | %s | %s   (tails <= %s, %s)"
          % (s, str(fo), dec(rels[0]), dec(rels[1]),
             dec(dvals[0][2], 20), dec(dvals[1][2], 20)))
mono = all(abs(mags[i][1][mi]) > abs(mags[i + 1][1][mi])
           for i in range(len(mags) - 1) for mi in range(2))
gate("G1.3b", mono, "|relative defect| strictly decreasing in the A-B separation at "
     "both mu rows (magnitudes tabulated; NO exponent fitted or claimed)")
# 2-D check row
s = 3
U2 = A1_2 | B1_2[s]
rUa2 = RUN(V2, U2, WRES2, MU2, lam0=U2)
rAa2 = RUN(V2, A1_2, WRES2, MU2, lam0=A1_2)
rBa2 = RUN(V2, B1_2[s], WRES2, MU2, lam0=B1_2[s])
ords2 = [ordcount(rUa2, w, C2) - ordcount(rAa2, w, C2) - ordcount(rBa2, w, C2)
         for w in range(W_ORD + 1)]
gate("G1.3c", all(d <= 0 for d in ords2) and any(d < 0 for d in ords2),
     "2-D discriminator: same R2 defect sign structure (composition split is not a "
     "dimension artifact)")

# ================================================================ SECTION 2
print("\nSECTION 2 -- SCREENING vs ACCUMULATION under the C-72 encoding split")
print("-" * 78)
print("  composed coupling at probe c:  Phi(c) = sum_i q_i * G(x_i -> c);")
print("  ORIENTATION encoding: q_i in {+1,-1}.  OCCUPANCY encoding: q_i = 1 if present.")
print("  G0.4 already gated: every kernel G(x_i->c) > 0 -- the mediator is sign-blind;")
print("  whatever sign structure the encoding has is EXACTLY what survives mediation.")

# exact mirror cancellation (orientation, +1/-1, symmetric pair)
rm_p = RUN(V3, Am, WRES3, MU3)
rm_m = RUN(V3, Bm, WRES3, MU3)
po = all(ordcount(rm_p, w, [c]) == ordcount(rm_m, w, [c])
         for w in range(W_ORD + 1) for c in plane5)
rs = all(rm_p['acc'][mi][c] == rm_m['acc'][mi][c] for mi in range(2) for c in plane5)
gate("G2.2", po and rs, "EXACT SCREENING: mirror pair (+1 at x=1, -1 at x=9), all 25 "
     "probes on the mid-plane x=5: Phi = 0 at every order and in both resummed rows "
     "(count equality computed, not assumed from symmetry)")

# dipole vs monopole (orientation, adjacent +1/-1)
rd_p = RUN(V3, Dp, WRES3, MU3)
rd_m = RUN(V3, Dm, WRES3, MU3)
print("\n    dipole (+ at x=2, - at x=3) vs monopole (+ at x=3 alone), probes (x,2,2):")
print("    probe x | mu    | Phi_dipole (exact) | Phi_mono | |Phi_dip|/Phi_mono")
ok_screen = True
ratseq = {0: [], 1: []}
for mi in range(2):
    t3 = gtail(rd_p, mi) + gtail(rd_m, mi)
    for x in range(5, 11):
        c = [V3.idx((x, 2, 2))]
        phid = gval(rd_p, mi, c) - gval(rd_m, mi, c)
        mono = gval(rd_m, mi, c)
        ok_screen &= (abs(phid) + t3 < mono)
        ratseq[mi].append(abs(phid) / mono)
        print("    %7d | 1/%-3d | %s | %s | %s"
              % (x, MU3[mi][1], dec(phid), dec(mono), dec(abs(phid) / mono, 8)))
gate("G2.3a", ok_screen, "PARTIAL SCREENING: |Phi_dipole| < Phi_monopole at every probe, "
     "both mu rows (interval-gated with tail bounds)")
for mi in range(2):
    dec_flags = [ratseq[mi][i] > ratseq[mi][i + 1] for i in range(len(ratseq[mi]) - 1)]
    print("    mu=1/%d: screening ratio strictly decreasing along probes: %s "
          "(REPORTED AS FOUND -- depth of screening is a mu-dependent fact, no law claimed)"
          % (MU3[mi][1], all(dec_flags)))

# occupancy accumulation
print("\n    OCCUPANCY encoding, k one-signed sources on the axis, probe (10,2,2):")
rs_occ = [RUN(V3, frozenset(occ_sites[:k]), WRES3, MU3) for k in range(1, 5)]
ok_inc = ok_add = ok_pos = True
for mi in range(2):
    prev = None
    for k in range(1, 5):
        Rk = rs_occ[k - 1]
        Gk = gval(Rk, mi, C3)
        add = sum(RUN(V3, frozenset([occ_sites[i]]), WRES3, MU3)['acc'][mi][c]
                  for i in range(k) for c in C3) == sum(Rk['acc'][mi][c] for c in C3)
        ok_add &= add
        if prev is not None:
            ok_inc &= (Gk > prev + gtail(Rk, mi))
        prev = Gk
        den = (Rk['sigma'] * MU3[mi][1]) ** Rk['W']
        tb = gtail(Rk, mi)
        ok_pos &= all(Fraction(Rk['acc'][mi][c], den) > tb for c in range(V3.n))
        if mi == 1:
            print("    k=%d  Phi(C) at mu=1/6 = %s" % (k, dec(Gk)))
gate("G2.4a", ok_add, "R1 accumulation is EXACTLY additive: Phi_k = sum of the k "
     "single-source couplings (integer identity on numerators)")
gate("G2.4b", ok_inc, "Phi strictly increasing in k, both mu rows")
gate("G2.4c", ok_pos, "one-signed everywhere: Phi_k(c) > 0 at EVERY site c "
     "(occupancy cannot screen -- there is nothing to cancel with)")

# R2 cross: opacity affects magnitude, never sign; mirror still cancels
Uocc = frozenset(occ_sites)
rUocc = RUN(V3, Uocc, WRES3, MU3, lam0=Uocc)
sub = True
for mi in range(2):
    lhs = gval(rUocc, mi, C3)
    rhs = sum(gval(RUN(V3, frozenset([s]), WRES3, MU3, lam0=frozenset([s])), mi, C3)
              for s in occ_sites)
    tb = gtail(rUocc, mi) + 4 * gtail(rs_occ[0], mi)
    sub &= (lhs + tb < rhs)
den0 = MU3[0][1] ** rUocc['W']
posR2 = all(ordcount(rUocc, w, [c]) >= 0 for w in range(W_ORD + 1) for c in range(V3.n))
rmm_p = RUN(V3, Am, WRES3, MU3, lam0=Am | Bm)
rmm_m = RUN(V3, Bm, WRES3, MU3, lam0=Am | Bm)
mirR2 = all(ordcount(rmm_p, w, [c]) == ordcount(rmm_m, w, [c])
            for w in range(W_ORD + 1) for c in plane5)
gate("G2.5a", sub, "R2 occupancy: subadditive (opacity defect < 0) yet still one-signed")
gate("G2.5b", posR2, "R2 occupancy: every per-order count >= 0 -- opacity reduces, "
     "NEVER flips a sign (screening-by-sign and shadowing-by-opacity are different "
     "mechanisms, computed apart)")
gate("G2.5c", mirR2, "R2 orientation mirror pair: still cancels exactly on the mid-plane")
print("\n  C-72 CONNECTION (computed, both venues' facts above): the encoding-level sign")
print("  structure SURVIVES mediation unchanged -- orientation encodings can screen")
print("  (G2.2 exact zero, G2.3 partial), occupancy encodings can only accumulate")
print("  (G2.4, G2.5); the mediator contributes magnitude only (G0.4 positivity).")

# ================================================================ SECTION 3
print("\nSECTION 3 -- SHADOWING: intervening region D between A and B")
print("-" * 78)
r_free = RUN(V3, As, WRES3, MU3)
# empty-D control -- INDEPENDENT code path: the scaled-weight machinery (sigma=2) with
# every weight = 1; must reproduce the free run exactly after unscaling.
r_empty = RUN(V3, As, WRES3, MU3, force_sigma=2)
po = all(Fraction(ordcount(r_empty, w, [c]), 2 ** w) == ordcount(r_free, w, [c])
         for w in range(W_ORD + 1) for c in Bs | C3)
rs_eq = all(gval(r_empty, mi, Bs) == gval(r_free, mi, Bs) for mi in range(2))
gate("G3.1", po and rs_eq,
     "EMPTY-D control: the weighted-obstacle machinery at lambda=1 everywhere "
     "(independent scaled code path, sigma=2) reproduces the free coupling EXACTLY")

# full slab, lambda=0: total opacity iff D is a vertex cut
r_slab0 = RUN(V3, As, WRES3, MU3, lam0=slab)
dcut = V3.bfs_dist(As, blocked=slab)
nopath = all(dcut[b] == -1 for b in Bs)
zero = all(ordcount(r_slab0, w, Bs) == 0 for w in range(W_ORD + 1)) and \
       all(sum(r_slab0['acc'][mi][b] for b in Bs) == 0 for mi in range(2))
gate("G3.2a", nopath and zero, "FULLY-BLOCKING control: 25-site slab at x=5, lambda=0: "
     "G(A->B) = 0 at every order and resummed; BFS confirms the slab is a vertex cut "
     "(total opacity is the geometry's, verified two ways)")
r_slabh = RUN(V3, As, WRES3, MU3, lamhalf=slab)
trans = all(gval(r_slabh, mi, Bs) > gtail(r_slabh, mi) for mi in range(2))
gate("G3.2b", trans, "beside the zero: the SAME slab at lambda=1/2 transmits (G > 0) -- "
     "the zero is the blocking, not the instrument")

# nested obstacle monotonicity (coverage sweep)
r_cen0 = RUN(V3, As, WRES3, MU3, lam0=frozenset(center))
r_disk0 = RUN(V3, As, WRES3, MU3, lam0=disk)
r_slab1 = RUN(V3, As, WRES3, MU3, lam0=slab1)
r_ring0 = RUN(V3, As, WRES3, MU3, lam0=ring)
chain = [("empty", r_free), ("center 1", r_cen0), ("3x3 disk", r_disk0),
         ("slab-1", r_slab1), ("full slab", r_slab0)]
po_mono = True
for i in range(len(chain) - 1):
    for w in range(W_ORD + 1):
        if ordcount(chain[i + 1][1], w, Bs) > ordcount(chain[i][1], w, Bs):
            po_mono = False
gate("G3.3a", po_mono, "NESTED COVERAGE empty c {center} c disk c slab-1 c slab: "
     "per-order counts monotone non-increasing at every order")
print("    G(A->B) at mu=1/6 by obstacle (exact partial sums; tails bounded):")
for (lbl, r) in chain + [("ring (slab-disk)", r_ring0)]:
    print("      %-16s %s" % (lbl, dec(gval(r, 1, Bs), 16)))
strict = all(gval(chain[i][1], mi, Bs) >
             gval(chain[i + 1][1], mi, Bs) + gtail(chain[i][1], mi) +
             gtail(chain[i + 1][1], mi)
             for mi in range(2) for i in range(len(chain) - 2))
gate("G3.3b", strict, "resummed G strictly decreasing along the nested chain "
     "(down to slab-1; the last step lands on exact zero)")

# first-order onset = earned separation in the punctured venue
d_free = min(V3.bfs_dist(As)[b] for b in Bs)
d_disk = min(V3.bfs_dist(As, blocked=disk)[b] for b in Bs)
fo = next((w for w in range(W_ORD + 1) if ordcount(r_disk0, w, Bs) != 0), None)
gate("G3.4", fo == d_disk and d_disk > d_free,
     "SHADOW ONSET = earned separation in the punctured venue: first nonzero order "
     "%s = BFS distance around the disk %d (> free distance %d) -- opacity LENGTHENS "
     "the earned separation; w_min = d survives obstacles" % (fo, d_disk, d_free))

# lambda monotonicity on the disk
r_diskh = RUN(V3, As, WRES3, MU3, lamhalf=disk)
lam_mono_ord = True
for w in range(W_ORD + 1):
    n0 = Fraction(ordcount(r_disk0, w, Bs))                      # sigma=1
    nh = Fraction(ordcount(r_diskh, w, Bs), 2 ** w)              # sigma=2 scaled
    n1 = Fraction(ordcount(r_free, w, Bs))
    if not (n0 <= nh <= n1):
        lam_mono_ord = False
lam_mono_res = all(gval(r_disk0, mi, Bs) + gtail(r_disk0, mi) + gtail(r_diskh, mi)
                   < gval(r_diskh, mi, Bs) and
                   gval(r_diskh, mi, Bs) + gtail(r_diskh, mi) + gtail(r_free, mi)
                   < gval(r_free, mi, Bs) for mi in range(2))
gate("G3.5", lam_mono_ord and lam_mono_res,
     "REWEIGHTED D: per-order and resummed coupling monotone in lambda "
     "(0 < 1/2 < 1) -- opacity is graded, the venue's own dial")

# shadow pattern on the far plane
print("\n    SHADOW PATTERN: source (1,2,2), 3x3 disk at x=5, probes = plane x=9;")
print("    R(y,z) = G_disk / G_free at mu=1/6 (exact ratio, 8 digits):")
Rmap = {}
for (y, z) in plane9:
    c = [V3.idx((9, y, z))]
    Rmap[(y, z)] = gval(r_disk0, 1, c) / gval(r_free, 1, c)
for y in range(5):
    print("      " + "  ".join(dec(Rmap[(y, z)], 8) for z in range(5)))
allless = all(v < 1 for v in Rmap.values())
minc = all(Rmap[(2, 2)] <= v for v in Rmap.values())
corners = [Rmap[(y, z)] for y in (0, 4) for z in (0, 4)]
maxcorn = all(any(abs(cv - max(Rmap.values())) == 0 for cv in corners)
              for _ in [0]) and max(Rmap.values()) in corners
gate("G3.6a", allless, "R < 1 at every far-plane probe (the obstacle dims everything)")
gate("G3.6b", minc, "R minimal at the axial probe (2,2) -- the geometric shadow, computed")
gate("G3.6c", maxcorn and len(set(corners)) == 1,
     "R maximal and equal at the four far-plane corners (venue symmetry on the numbers)")

# 2-D discriminator venue: same shadow facts
r_free2 = RUN(V2, Am2, WRES2, MU2)
r_slab02 = RUN(V2, Am2, WRES2, MU2, lam0=slab2)
r_disk02 = RUN(V2, Am2, WRES2, MU2, lam0=disk2)
nopath2 = all(V2.bfs_dist(Am2, blocked=slab2)[b] == -1 for b in Bm2)
zero2 = all(ordcount(r_slab02, w, Bm2) == 0 for w in range(W_ORD + 1))
d_free2 = min(V2.bfs_dist(Am2)[b] for b in Bm2)
d_disk2 = min(V2.bfs_dist(Am2, blocked=disk2)[b] for b in Bm2)
fo2 = next((w for w in range(W_ORD + 1) if ordcount(r_disk02, w, Bm2) != 0), None)
Rmap2 = {}
for y in plane10_2:
    c = [V2.idx((10, y))]
    Rmap2[y] = gval(r_disk02, 1, c) / gval(r_free2, 1, c)
gate("G3.7", nopath2 and zero2 and fo2 == d_disk2 and d_disk2 > d_free2 and
     all(v < 1 for v in Rmap2.values()) and
     all(Rmap2[3] <= v for v in Rmap2.values()),
     "2-D DISCRIMINATOR: full-column cut gives exact zero; disk onset %s = punctured "
     "BFS %d > free %d; R < 1 with the minimum on-axis -- every composition/shadow "
     "fact reproduces in D=2 (facts are venue-structural, magnitudes differ)"
     % (fo2, d_disk2, d_free2))
print("    2-D far-column ratios R(y), mu=1/4: " +
      "  ".join(dec(Rmap2[y], 6) for y in range(7)))

# ================================================================ SECTION 4
print("\nSECTION 4 -- COMPOSITION FACTS FOR THE JUDGE, and the one named comparison (D-1)")
print("-" * 78)
print("""  What the mediated coupling DOES, computed exactly on the venue's own lattices:
  F1  SUPERPOSITION is a property of the ADMISSIBILITY RULE, not of the walk sum:
      exact (defect identically zero, per order, an integer identity) under the
      transparent endpoint rule R1; strictly subadditive under the opaque-source rule
      R2, with a negative defect that shrinks as the sources separate (tabulated).
  F2  The mediator kernel is POSITIVE: sign lives in the source values only.  The C-72
      encoding split survives mediation intact -- ORIENTATION encodings screen (exact
      mirror zero; dipole partially screened, mu-dependent depth), OCCUPANCY encodings
      accumulate (exact additivity under R1, one-signed at every probe, no possible
      cancellation).  Opacity (R2) changes magnitudes, never signs.
  F3  SHADOWING exists under R2/lambda and is the venue's own opacity: exact zero
      through a separating slab, graded in lambda, monotone in coverage, directional
      (axial minimum on the far plane), and equal to a LENGTHENING of the earned
      separation (shadow onset order = BFS distance in the punctured venue -- C-80's
      w_min = d, surviving obstacles).  Under R1 with transparent sources there is NO
      shadowing at all.
  THE COMPARISON (named here only): a Newtonian source term requires exact
  superposition, one-signed accumulation, no screening, and no shadowing (owners:
  Newton; the no-screening/no-shadowing character of gravity vs the two-signed,
  screenable electromagnetic coupling -- standard classical field taxonomy, Maxwell/
  Heaviside territory).  The instrument REACHES that composition profile in exactly one
  place: the ENDPOINT RULE R1 with OCCUPANCY-encoded sources -- where all four facts
  hold exactly (F1-R1, F2-occupancy).  Under the opaque rule, or with orientation
  encoding, the composition profile is electromagnetic-like (screening) or absorptive
  (shadowing) instead.  WHICH admissibility rule and encoding the record surface itself
  enforces is NOT decided by this lane; that is the judge's question, and these are the
  facts it needs.  No falloff exponent is used, fitted, or earned anywhere above.""")

nfail = sum(1 for (_, ok) in GATES if not ok)
print("=" * 78)
print("GATES: %d/%d PASS%s" % (len(GATES) - nfail, len(GATES),
      "" if nfail == 0 else "  -- FAILURES: " +
      ", ".join(t for (t, ok) in GATES if not ok)))
print("=" * 78)
sys.exit(0 if nfail == 0 else 1)
