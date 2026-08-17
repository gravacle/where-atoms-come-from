# rule_main2.py -- W-19 RULING, part 2.  The three tests that decide the ruling.
import numpy as np, time
from rule_verify import *

t0 = time.time()
P("=" * 118)
P("W-19 RULING -- PART 2.  Fragment-rule dependence, algebra dependence, and falsifiability.")
P("=" * 118)

def three(vec, L, l, F):
    vX = hadamard_all(vec, L); pX = np.abs(vX) ** 2
    return (channel_EXT(vec, L, l, F), channel_CHI(vX, L, l, F), channel_CL(pX, L, l, F))

def HS_HX(vec, L, l):
    HS = S_ax(vec, L, AX(L, [l]))
    pX = np.abs(hadamard_all(vec, L)) ** 2
    pm = pX.reshape([2] * L).sum(axis=tuple(a for a in range(L) if a != L - 1 - l))
    HX = float(-(pm[pm > 1e-15] * np.log2(pm[pm > 1e-15])).sum())
    return HS, HX

def curve(name, car, vec, l, frags):
    L = car.L
    HS, HX = HS_HX(vec, L, l)
    vX = hadamard_all(vec, L); pX = np.abs(vX) ** 2
    E = [channel_EXT(vec, L, l, F) for F in frags]
    Cc = [channel_CHI(vX, L, l, F) for F in frags]
    K = [channel_CL(pX, L, l, F) for F in frags]
    nz = lambda x, h: (x / h if h > 1e-12 else 0.0)
    P("  %-34s H_EXT(S)=%.9f  H_elec(S)=%.9f" % (name, HS, HX))
    P("      |F| :" + "".join("%10d" % len(F) for F in frags))
    P("      EXT :" + "".join("%10.6f" % nz(x, HS) for x in E))
    P("      CHI :" + "".join("%10.6f" % nz(x, HX) for x in Cc))
    P("      CL  :" + "".join("%10.6f" % nz(x, HX) for x in K))
    pts = lambda vs: sum(1 for v in vs if abs(v - 1.0) <= 0.10)
    P("      points: EXT=%d CHI=%d CL=%d" % (pts([nz(x, HS) for x in E]),
                                             pts([nz(x, HX) for x in Cc]),
                                             pts([nz(x, HX) for x in K])))
    return HS, HX, E, Cc, K

# ================================================================== 6  lane B READING 1 geometry
P("")
P("[6] LANE B READING 1, RE-RUN THROUGH ALL THREE ALGEBRA CHANNELS.")
P("    theta_L, S = link 0, NESTED-BY-INDEX fragments {1},{1,2},...  -- lane B's own fragment rule,")
P("    not lane A's BFS rule (on theta, d = 1, so lane A's rule has only one fragment).")
for Lk in (6, 8):
    car = Carrier("theta_%d" % Lk, *theta(Lk))
    frags = [list(range(1, k + 1)) for k in range(1, Lk)]
    wp, psi = car.ground(0.50); gs = car.lift(psi)
    curve("theta_%d ground g2=0.50" % Lk, car, gs, 0, frags)
    ghz = np.zeros(1 << Lk); ghz[0] = ghz[(1 << Lk) - 1] = 1 / np.sqrt(2)
    P("      [magnetic GHZ gauss residual %.1e]" % gauss_residual(car, ghz))
    curve("theta_%d MAGNETIC GHZ (g2->0)" % Lk, car, ghz, 0, frags)
    # electric GHZ = (|+..+> + |-..->)/sqrt2, written in the Z basis
    plus = np.ones(1 << Lk) / 2 ** (Lk / 2.0)
    sgn = 1.0 - 2.0 * (np.bitwise_count(np.arange(1 << Lk)) & 1)
    minus = sgn / 2 ** (Lk / 2.0)
    eghz = (plus + minus) / np.linalg.norm(plus + minus)
    P("      [electric GHZ gauss residual %.1e]" % gauss_residual(car, eghz))
    curve("theta_%d ELECTRIC GHZ" % Lk, car, eghz, 0, frags)
    for s in (11, 22):
        curve("theta_%d HAAR physical seed=%d" % (Lk, s), car, car.lift(haar_physical(car, s)), 0, frags)

# ================================================================== 7  partition dependence
P("")
P("[7] IS THE PLATEAU A PROPERTY OF THE STATE OR OF THE CHOSEN PARTITION?")
P("    Lane A reads R_delta off the d BFS-LEVEL CUTS.  Here the SAME state, SAME carrier, SAME")
P("    algebra and SAME number of disjoint fragments, with the PARTITION RANDOMISED.")
P("    Channel = CL (gauge-invariant on both sides).  delta = 0.10.")
rng = np.random.default_rng(20260817)
for nm, (V, E) in [("dbl_chain9", mg_chain(5)), ("tri_chain12", tri_chain12()), ("heawood", heawood())]:
    car = Carrier(nm, V, E); L = car.L
    wp, psi = car.ground(0.50); vec = car.lift(psi)
    pX = np.abs(hadamard_all(vec, L)) ** 2
    HS, HX = HS_HX(vec, L, 0)
    # lane A's rule-C cuts: BFS level cuts
    u, v = car.edges[0]
    dist = bfs_dist(V, E, u, 0); d = dist[v]
    cuts = []
    for k in range(d):
        cuts.append([i for i, (a, b) in enumerate(E) if i != 0 and min(dist[a], dist[b]) == k and max(dist[a], dist[b]) == k + 1])
    cuts = [c for c in cuts if c]
    okc = sum(1 for c in cuts if channel_CL(pX, L, 0, c) >= 0.9 * HX)
    P("  %-12s L=%2d  H_elec(S)=%.9f   BFS-CUT partition sizes %s -> R_delta = %d of %d"
      % (nm, L, HX, [len(c) for c in cuts], okc, len(cuts)))
    env = [i for i in range(L) if i != 0]
    k = len(cuts)
    tally = []
    for _ in range(200):
        pm = list(env); rng.shuffle(pm)
        parts = [sorted(pm[j::k]) for j in range(k)]
        tally.append(sum(1 for pt in parts if channel_CL(pX, L, 0, pt) >= 0.9 * HX))
    tally = np.array(tally)
    P("      200 RANDOM equal-size disjoint partitions into %d parts: R_delta mean %.3f  max %d  "
      "distribution %s" % (k, tally.mean(), tally.max(),
                           {int(x): int((tally == x).sum()) for x in np.unique(tally)}))

# ================================================================== 8  falsifiability of the gauge-invariant channel
P("")
P("[8] CAN THE GAUGE-INVARIANT PLATEAU FAIL?  Two regimes, same criterion, same channel (CL).")
P("    (a) d >= 4 carriers, BFS-cut fragments: ground state vs HAAR PHYSICAL states.")
for nm, (V, E) in [("tri_chain12", tri_chain12()), ("heawood", heawood())]:
    car = Carrier(nm, V, E); L = car.L
    frags, d = rule_A_fragments(V, E, 0)
    for tag, vec in [("ground g2=0.50", car.lift(car.ground(0.50)[1])),
                     ("ground g2=3.00", car.lift(car.ground(3.00)[1])),
                     ("HAAR seed 7", car.lift(haar_physical(car, 7))),
                     ("HAAR seed 8", car.lift(haar_physical(car, 8)))]:
        pX = np.abs(hadamard_all(vec, L)) ** 2
        HS, HX = HS_HX(vec, L, 0)
        vals = [channel_CL(pX, L, 0, F) / HX if HX > 1e-12 else 0.0 for F in frags]
        P("      %-12s %-16s H_elec(S)=%.9f  CL/H = %s  points=%d"
          % (nm, tag, HX, " ".join("%.6f" % x for x in vals),
             sum(1 for x in vals if abs(x - 1) <= 0.10)))
P("    (b) d = 1 carrier (theta_6), 5 DISJOINT single-link fragments: state-carried, not forced.")
car = Carrier("theta_6", *theta(6)); L = 6
frags = [[j] for j in range(1, 6)]
plus = np.ones(1 << L) / 2 ** (L / 2.0)
sgn = 1.0 - 2.0 * (np.bitwise_count(np.arange(1 << L)) & 1)
eghz = (plus + sgn / 2 ** (L / 2.0)); eghz /= np.linalg.norm(eghz)
for tag, vec in [("ELECTRIC GHZ", eghz),
                 ("ground g2=0.50", car.lift(car.ground(0.50)[1])),
                 ("ground g2=0.10", car.lift(car.ground(0.10)[1])),
                 ("ground g2=3.00", car.lift(car.ground(3.00)[1])),
                 ("HAAR seed 7", car.lift(haar_physical(car, 7))),
                 ("HAAR seed 8", car.lift(haar_physical(car, 8)))]:
    pX = np.abs(hadamard_all(vec, L)) ** 2
    HS, HX = HS_HX(vec, L, 0)
    vals = [channel_CL(pX, L, 0, F) / HX if HX > 1e-12 else 0.0 for F in frags]
    P("      theta_6 %-16s H_elec(S)=%.9f  CL/H on 5 disjoint fragments = %s  R_delta=%d"
      % (tag, HX, " ".join("%.6f" % x for x in vals), sum(1 for x in vals if x >= 0.9)))
P("    (c) the same on theta_8 with 7 disjoint single-link fragments.")
car = Carrier("theta_8", *theta(8)); L = 8
frags = [[j] for j in range(1, 8)]
plus = np.ones(1 << L) / 2 ** (L / 2.0)
sgn = 1.0 - 2.0 * (np.bitwise_count(np.arange(1 << L)) & 1)
eghz = (plus + sgn / 2 ** (L / 2.0)); eghz /= np.linalg.norm(eghz)
ghz = np.zeros(1 << L); ghz[0] = ghz[(1 << L) - 1] = 1 / np.sqrt(2)
for tag, vec in [("ELECTRIC GHZ", eghz), ("MAGNETIC GHZ", ghz),
                 ("ground g2=0.50", car.lift(car.ground(0.50)[1])),
                 ("HAAR seed 7", car.lift(haar_physical(car, 7)))]:
    pX = np.abs(hadamard_all(vec, L)) ** 2
    HS, HX = HS_HX(vec, L, 0)
    vals = [channel_CL(pX, L, 0, F) / HX if HX > 1e-12 else 0.0 for F in frags]
    P("      theta_8 %-16s H_elec(S)=%.9f  CL/H = %s  R_delta=%d"
      % (tag, HX, " ".join("%.4f" % x for x in vals), sum(1 for x in vals if x >= 0.9)))

# ================================================================== 9  is the electric GHZ gauge content zero?
P("")
P("[9] ZERO-VARIABLE CHECK ON THE (b)/(c) EXHIBIT.  Does the Gauss law enter the theta_6 electric-GHZ")
P("    plateau at all?  Compare with the SAME state on L BARE QUBITS with NO Gauss law and no")
P("    plaquettes: identical vector, identical estimator.  If the curves are byte-equal, the")
P("    exhibit has ZERO VARIABLES MOVED with respect to gauge invariance.")
for L in (6, 8):
    plus = np.ones(1 << L) / 2 ** (L / 2.0)
    sgn = 1.0 - 2.0 * (np.bitwise_count(np.arange(1 << L)) & 1)
    e = (plus + sgn / 2 ** (L / 2.0)); e /= np.linalg.norm(e)
    pX = np.abs(hadamard_all(e, L)) ** 2
    vals = [channel_CL(pX, L, 0, [j]) for j in range(1, L)]
    P("      L=%d  electric GHZ, CL(X_0:X_j) = %s   (the state is a classical repetition code in the"
      " X basis; the Gauss law only fixes the PARITY, which the fragments never use)"
      % (L, " ".join("%.6f" % x for x in vals)))

P("")
P("elapsed %.1f s" % (time.time() - t0))
open("OUT_rule_main2.txt", "w").write("\n".join(LOG) + "\n")
