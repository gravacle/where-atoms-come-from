#!/usr/bin/env python3
# LANE W08 / M3 — script 4
# (3) THE THREE-CLASS vs FOUR-CLASS QUESTION, and the SCOPE of W-03's multiset theorem.
# (4) m(P) at K1's published ready state, against the registered generic and resonant values.
# Seed 20260816.  Double precision unless a line says EXACT.
import numpy as np
from itertools import permutations

rng = np.random.default_rng(20260816)
L = []
def out(s=""):
    print(s); L.append(s)

out("=" * 96)
out("M3-4  THREE CLASSES vs FOUR, THE SCOPE OF THE MULTISET THEOREM, AND THE MAHLER VALUES")
out("=" * 96)
out("numpy %s ; IEEE double unless a line says EXACT." % np.__version__)
out()

# ---------------------------------------------------------------- p00 = 0 from incidence
EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
FACE_V, CYC_V = {0, 1, 2}, {0, 3, 4}
classes = {v: (int(v in FACE_V), int(v in CYC_V)) for v in range(5)}
counts = {}
for v, cl in classes.items():
    counts[cl] = counts.get(cl, 0) + 1
out("(a) p00 = 0 ON K1 IS FORCED BY INCIDENCE, NOT CHOSEN.")
out("    FACE_V u CYC_V = %s  =  all five vertices." % sorted(FACE_V | CYC_V))
out("    class counts on K1: %s" % {("%d%d" % c): n for c, n in sorted(counts.items())})
out("    class 00 is EMPTY, so p00 = 0 for EVERY ready state on K1.  Any statement about a")
out("    weight vector with p00 > 0 is a statement about some OTHER carrier.")
out("    S4's own carrier list contains three such carriers (S4_THE_MEASUREMENT_V001.md:582,")
out("    :519, :546): B0b {00:4,01:1,10:2,11:2}, B4 spindle {00:1,01:1,10:1,11:3},")
out("    B1q K1-bridged+spectator {00:1,01:3,10:3}.  So the p00 > 0 regime is NOT vacuous")
out("    for the corpus -- only for K1.")
out()

# ---------------------------------------------------------------- permutation scope
out("(b) HOW MUCH OF W-03's 24-of-24 PERMUTATION INVARIANCE IS ABOUT K1 AT ALL?")
out("    A permutation sigma of the four classes maps a K1 weight vector (0,a,b,c) to a vector")
out("    realisable on K1 iff it leaves 0 in slot 00, i.e. iff sigma^{-1}(00) carries weight 0.")


def in_family(p, sigma):
    q = tuple(p[i] for i in sigma)
    return q[0] == 0.0


for lab, p in [("generic K1 state (0, .3, .3, .4)", (0.0, .3, .3, .4)),
               ("K1 PUBLISHED state (0, 0, .5, .5)", (0.0, 0.0, .5, .5)),
               ("K1 uniform state  (0, .4, .4, .2)", (0.0, .4, .4, .2))]:
    n = sum(in_family(p, s) for s in permutations(range(4)))
    out("      %-36s : %2d of 24 permutations stay inside K1's family" % (lab, n))
out("    So at a generic K1 ready state EIGHTEEN of the twenty-four permutations move weight")
out("    INTO the empty class 00 and leave K1 entirely.  The surviving six are exactly S_3 on")
out("    {10,01,11} -- and those six are the ones that matter for K1, because they are the ones")
out("    that say the loop-incidence LABELS (curvature loop vs flat loop vs root) are invisible.")
out("    VERDICT: three quarters of the multiset theorem's permutation content is a statement")
out("    about carriers K1 does not contain.  It is a true and proved statement (the proof is")
out("    in the corpus, LANE_S5_CHARGE_CODE/s5_B_sweep.OUT.txt B5), it is non-vacuous for the")
out("    corpus's OTHER carriers, and on K1 itself it says exactly six things, not twenty-four.")
out()

# ---------------------------------------------------------------- which results break at p00>0
out("(c) DOES ANY REGISTERED RESULT SILENTLY USE THE FOUR-CLASS FORM WHERE p00 = 0 IS DOING")
out("    UNDECLARED WORK?  Four candidates, each tested.")
out()


def m_quad(p, n=400000):
    a, b, c, d = p
    t = 2 * np.pi * (np.arange(n) + 0.5) / n
    x = np.exp(1j * t)
    return float(np.mean(np.log(np.maximum(np.abs(a + b * x), np.abs(c + d * x)))))


def orbit_avg(p, f, c, N):
    a, b, cc, d = p
    k = np.arange(1, N + 1, dtype=np.float64)
    z = a + b * np.exp(-1j * f * k) + cc * np.exp(1j * c * k) + d * np.exp(1j * (c - f) * k)
    return float(np.mean(np.log(np.abs(z))))


out("    (c1) N1, lambda = m(p00+p10 x+p01 y+p11 xy).  SAFE.  Its proof (Jensen in y, then")
out("         Weyl on the rank-0 orbit) never uses p00 = 0.  Checked on a p00 > 0 carrier state:")
f1, c1 = 2 * np.pi * (np.sqrt(5) - 1) / 2, 2 * np.pi * np.sqrt(2)
for lab, p in [("B4 spindle uniform (1/6,1/6,1/6,1/2)", (1 / 6, 1 / 6, 1 / 6, 0.5)),
               ("B0b uniform        (4/9,1/9,2/9,2/9)", (4 / 9, 1 / 9, 2 / 9, 2 / 9)),
               ("balanced           (1/4,1/4,1/4,1/4)", (.25, .25, .25, .25))]:
    mq, oa = m_quad(p), orbit_avg(p, f1, c1, 2 * 10 ** 6)
    out("         %-38s m(P) = %+.9f   orbit avg N=2e6 = %+.9f   dev %.2e"
        % (lab, mq, oa, abs(mq - oa)))
out("         (B0b uniform reproduces the register's log(4/9) = -0.810930216216 exactly:")
out("          m = %.12f, log(4/9) = %.12f)" % (m_quad((4/9, 1/9, 2/9, 2/9)), np.log(4/9)))
out("         NOTE ON THE THIRD ROW, stated rather than buried: the balanced state has")
out("         P = (1+x)(1+y)/4, whose zero set is TWO CIRCLES, so both the quadrature and the")
out("         orbit average converge slowly around an integrable log singularity.  The exact")
out("         value is log(1/4) = %.12f; quadrature n=4e5 gives %.9f and the N=2e6 orbit gives"
    % (np.log(0.25), m_quad((.25, .25, .25, .25))))
out("         %.9f.  Agreement to 2.4e-06 is the CONVERGENCE RATE, not the accuracy of N1."
    % orbit_avg((.25, .25, .25, .25), f1, c1, 2 * 10 ** 6))
out()

out("    (c2) W-03's MULTISET THEOREM.  SAFE AS A RANK-0 THEOREM, and it does NOT need p00 = 0.")
out("         24 permutations of a p00 > 0 state, quadrature n = 4e5:")
p = (0.4, 0.3, 0.2, 0.1)
vals = [m_quad(tuple(p[i] for i in s)) for s in permutations(range(4))]
out("         p = (0.4,0.3,0.2,0.1): spread over 24 permutations = %.3e ; value = %.12f"
    % (max(vals) - min(vals), vals[0]))
out("         BUT IT IS A RANK-0 STATEMENT ONLY.  At a TORSION connection (rank L = 2) the")
out("         average is over a finite orbit and permutation invariance dies.  S1's own")
out("         published connection, order 4, is such a point.  EXACT ARITHMETIC (the orbit is")
out("         u = -1, v = -i, uv = i, period 4), so no float can hide an exact zero:")
out("           Z_1 = (p00-p10) + i(p11-p01)      Z_3 = conj(Z_1)")
out("           Z_2 = p00+p10-p01-p11             Z_4 = 1")
out("           lambda = (1/4)[ 2 log|Z_1| + log|Z_2| ]")
from fractions import Fraction as Fr
pf = (Fr(2, 5), Fr(3, 10), Fr(1, 5), Fr(1, 10))


def lam_order4_exact(q):
    z1sq = (q[0] - q[1]) ** 2 + (q[3] - q[2]) ** 2      # EXACT rational
    z2 = abs(q[0] + q[1] - q[2] - q[3])                 # EXACT rational
    if z1sq == 0 or z2 == 0:
        return float("-inf")
    # lambda = (1/4)[log|Z_1| + log|Z_2| + log|Z_3| + log|Z_4|], |Z_3|=|Z_1|, |Z_4|=1
    return (np.log(float(z1sq)) + np.log(float(z2))) / 4.0


vals2 = sorted(set(round(lam_order4_exact(tuple(pf[i] for i in s)), 9)
                   for s in permutations(range(4))))
nz = sum(1 for s in permutations(range(4)) if lam_order4_exact(tuple(pf[i] for i in s)) == float("-inf"))
out("         p = (2/5,3/10,1/5,1/10) at (f,c) = (pi, 3pi/2), all 24 permutations:")
out("           distinct lambda values = %s"
    % ["-inf" if not np.isfinite(x) else "%.9f" % x for x in vals2])
out("           %d of 24 permutations give lambda = -infinity EXACTLY (Z_2 = 0, since" % nz)
out("           2/5 + 1/10 = 1/2 = 3/10 + 1/5 -- a pairing coincidence in this multiset).")
out("         SPREAD = INFINITE -> multiset invariance FAILS at rank two, and fails hardest.")
out("         *** CONFOUND RECORDED: my first pass computed this in double over 4000 cells and")
out("         *** got a finite spread of 7.102087, because an exactly-zero |Z_k| came back as")
out("         *** ~1e-17.  Same float trap as m3_3(e).  The exact answer is -infinity.")
out("         This is the SAME defect LANE_G_GROUP_REFUTER/g4b_multiset.OUT.txt found for the")
out("         refined weights (spread 0.2792).  It is here reproduced for the FOUR CLASS")
out("         weights on K1's OWN published connection, which is where it bites this corpus.")
out("         DISCOUNT: this is a re-finding, not a new finding.  Credit is theirs.")
out()

out("    (c3) W-03's PINCH = SPECTATOR INVOLUTION (00<->11, 10<->01).  SAFE, and stronger than")
out("         the multiset theorem: it is exact POINTWISE in k and at every connection, not")
out("         only after integration.  Verified, and contrasted with a non-involution swap:")
worst_inv = worst_other = 0.0
for _ in range(2000):
    q = rng.dirichlet([1, 1, 1, 1])
    f, c = rng.uniform(0, 2 * np.pi, 2)
    k = rng.integers(1, 50)
    u, v = np.exp(-1j * f * k), np.exp(1j * c * k)
    base = abs(q[0] + q[1] * u + q[2] * v + q[3] * u * v)
    inv = (q[3], q[2], q[1], q[0])            # 00<->11, 10<->01
    oth = (q[1], q[0], q[2], q[3])            # 00<->10 alone: NOT a pointwise symmetry
    zi = abs(inv[0] + inv[1] * u + inv[2] * v + inv[3] * u * v)
    zo = abs(oth[0] + oth[1] * u + oth[2] * v + oth[3] * u * v)
    worst_inv = max(worst_inv, abs(base - zi))
    worst_other = max(worst_other, abs(base - zo))
out("         2000 random (p,f,c,k): max | |Z_k| - |Z_k| under 00<->11,10<->01 | = %.3e" % worst_inv)
out("                                max | |Z_k| - |Z_k| under 00<->10 only     | = %.3e" % worst_other)
out("         The first is an identity; the second is not.  The multiset theorem's other 22")
out("         permutations hold ONLY after averaging over the full torus.")
out()

out("    (c4) W-01's CONVEX-HULL CRITERION, READ ON THE STATE SIDE.  ***THIS IS THE ONE THAT")
out("         BREAKS.***  'Fires iff the weights obey the triangle/polygon inequality' is TRUE")
out("         at three classes (M3-1) and FALSE at four (M3-2).  On K1 nobody can notice,")
out("         because p00 = 0 makes the two coincide.  Quantified in m3_2: the polygon reading")
out("         calls exactly 1/2 of the four-class simplex a firer; the truth is exactly 1/4.")
out("         The corpus's own B4 spindle at its own uniform state (1/6,1/6,1/6,1/2) satisfies")
out("         the polygon inequality WITH EQUALITY (1/2 = 1/6+1/6+1/6) and has NO torus zero:")
out("           min_{T^2} |P| = %.6f   (not 0)" % (
    min(abs(abs(1/6 + 1/6 * np.exp(1j * t)) - abs(1/6 + 0.5 * np.exp(1j * t)))
        for t in 2 * np.pi * np.arange(4000) / 4000)))
out("         and S4's own closed form for it, lambda = log(1/2), is exactly the 'no zero'")
out("         value m(P) = log p_max.  The counterexample was already sitting in the corpus,")
out("         computed and published, with nobody having asked the question it answers.")
out("         NOTE THE READING IS TWO-WAY: this does NOT show any registered SENTENCE is false.")
out("         The register states the criterion with THREE coefficients (REGISTER_V001.md:43)")
out("         and S4 states it for (p0,q,r) (S4:229-231).  Both are correct as written and")
out("         both are K1-only.  What is established is that the obvious four-class")
out("         generalisation -- the one N1's own polynomial invites -- is FALSE, and that no")
out("         page in the corpus says so.")
out()

# ---------------------------------------------------------------- (4) the Mahler values
out("=" * 96)
out("(4)  m(P) AT K1's PUBLISHED READY STATE, AGAINST THE REGISTERED VALUES")
out("=" * 96)
P_PUB = (0.0, 0.0, 0.5, 0.5)
P_S3 = (0.0, 0.3, 0.3, 0.4)
out("    K1 published p_v=(1/2,0,0,1/4,1/4) -> (p00,p10,p01,p11) = (0, 0, 1/2, 1/2)")
out("      P = (y/2)(1+x),  m(P) = log(1/2) + m(y) + m(1+x) = -log 2   EXACTLY")
out("      -log 2                       = %.12f" % (-np.log(2)))
out("      quadrature n = 4e5           = %.12f" % m_quad(P_PUB))
out("      quadrature n = 4e6           = %.12f" % m_quad(P_PUB, 4 * 10 ** 6))
out()
out("    S3/S4 sense-C state (0, .3, .3, .4)  [ = S3's p=(0.4,.15,.15,.15,.15) ]")
out("      m(P) quadrature n = 4e6      = %.12f" % m_quad(P_S3, 4 * 10 ** 6))
out("      registered generic value     = -0.767507880   (Cassaigne-Maillot, S4:912)")
out("      m(0.4+0.3x+0.3y) same state, other Jensen grouping = %.12f"
    % m_quad((0.4, 0.3, 0.3, 0.0), 4 * 10 ** 6))
out()
out("    THE RESONANT POINT OF THE ERRATUM AGAINST W-02: f = 2.0, c = 1.1, -11f + 20c = 0.")
out("      check: -11*2.0 + 20*1.1 = %.1e" % (-11 * 2.0 + 20 * 1.1))
for lab, p, reg in [("S3 state (0,.3,.3,.4)", P_S3, -0.767014993),
                    ("K1 PUBLISHED (0,0,.5,.5)", P_PUB, None)]:
    o = orbit_avg(p, 2.0, 1.1, 10 ** 7)
    s = "      %-26s orbit average N=1e7 on the resonant connection = %+.9f" % (lab, o)
    if reg is not None:
        s += "   (registered corrected value %+.9f, dev %.1e)" % (reg, abs(o - reg))
    out(s)
out("      AT THE PUBLISHED STATE THE RESONANCE IS INVISIBLE: c does not enter |Z_k| at all,")
out("      so the subtorus average equals the full-torus average equals -log 2 = %.9f."
    % (-np.log(2)))
out("      The erratum against W-02 corrects a number that, at K1's own published ready state,")
out("      does not exist as a separate number.")
out()
out("    THE FOUR VALUES SIDE BY SIDE")
out("      -0.767507880   registered GENERIC (rank L = 0), S3/S4 sense-C state")
out("      -0.767014993   registered RESONANT subtorus value, same state, f=2.0 c=1.1")
out("      -0.693147181   m(P) at K1's OWN PUBLISHED ready state          [= -log 2, EXACT]")
out("      -infinity      the ACTUAL per-cell rate at K1's own published state AND connection")
out("    differences:  -log2 - (-0.767507880) = %+0.9f" % (-np.log(2) + 0.767507880))
out("                  -log2 - (-0.767014993) = %+0.9f" % (-np.log(2) + 0.767014993))
out("    The published state is SLOWER than generic by 0.0744 in the Mahler reading and")
out("    INFINITELY FASTER than generic in the actual-orbit reading.  Both are correct; they")
out("    are answers to different questions, and only the second is a rate of anything.")
out()
out("    CROSS-CHECK AGAINST THE CORPUS'S OWN TABLE (S4:785, :582).  -log 2 = log(1/2) is")
out("    already in S4 twice, as B1p's and B4's uniform-state rate.  K1's published ready")
out("    state is the same degenerate two-class point reached by a different route:")
out("      S4 B1p K1-bridged, uniform, classes {01:3,10:3}: log max(1/2,1/2) = %.12f"
    % np.log(0.5))
out("      M3 K1, PUBLISHED state, classes {11:1/2, 01:1/2}: -log 2      = %.12f" % (-np.log(2)))
out()
out("DONE.")

open("m3_4_threeclass_vs_fourclass.OUT.txt", "w").write("\n".join(L) + "\n")
