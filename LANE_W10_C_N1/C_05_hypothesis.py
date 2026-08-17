#!/usr/bin/env python3
"""
LANE W-10 / C — STEP 5.  THE HYPOTHESIS N1 MUST STATE.  THE LANE'S LOAD-BEARING ITEM.

The brief: "log|P| is not Riemann-integrable when P has zeros on T^2, and a four-class P HAS
zeros on half the connection space (W-09).  So the Birkhoff convergence needs more than Weyl.
Name the theorem.  Lawton 1983 is the candidate and W-03 recorded it MISSING from the IMPORT
AUDIT.  Say what N1 must state as a hypothesis to be publishable."

FIVE LEGS.
  1  THE PREMISE IS CORRECTED FIRST.  Neither four-class carrier the corpus owns has a zero on
     T^2 under SENSE U (C_03 sec 4).  The set of four-class ready states whose P has a torus
     zero has measure EXACTLY 1/4 — the same as for three classes.  W-09's 1/2 is a different
     set (connection space, weights free).  So four-class occupancy is NOT the variable that
     turns the hypothesis on.  THE VARIABLE IS THE READY STATE, through one inequality.
  2  B0b* — THE SAME COMPLEX, ONE DESIGNATED LOOP MOVED.  A four-class carrier the corpus's own
     3x3 torus supplies, whose P DOES have torus zeros.  Built from incidence.
  3  THE EXACT FAILURE ALREADY IN THE PUBLISHED TABLE.  S4's four-class SENSE C row
     (lambda = -1.386294361120, both B0b and B4) at S1's own published connection: Z_k = 0
     EXACTLY unless k = 0 mod 4 (3 of every 4 cells), so lambda = -infinity for every
     N >= 1.  Shown in exact rational arithmetic, and
     shown to be INVISIBLE to float64, which reports a finite -28.09.
  4  THE FINITE-N FAILURE, SWEPT.  A four-class, NON-FACTORING P with an EXACT Mahler measure
     log(2/5) and a torus zero at (-1,-1), on a connection whose orbit is dense in T^2
     (L = {0}, proved).  The average sits arbitrarily far from m(P) for arbitrarily long.
  5  THE LIMIT FAILURE, RIGOROUSLY.  A Liouville alpha with L = {0} for which
     liminf (1/N) sum log|Z_k| = -infinity on a FOUR-CLASS ready state, by a bound that needs
     no summation: every term is <= 0, so the average at N = q is at most log|Z_q| / q.

Precision: exact Fractions / mpmath where marked EXACT; float64 elsewhere.
Seed: master 20260816, offset +5.
"""
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
import sys
from C_01_carriers import Carrier, build_B0b, rank_Q, in_span

mp.mp.dps = 60

# ============================================================ LEG 2 : B0b*, from incidence
def build_B0b_star():
    """The SAME 3x3 torus complex as B0b.  ONE THING MOVES: gamma_C is the (1,1)-diagonal
    cycle instead of the row-0 horizontal cycle.  Same gamma_F, same d1, same d2."""
    b = build_B0b()
    H, W = b.H, b.W
    gC = {H[(0, 0)]: 1, W[(0, 1)]: 1, H[(1, 1)]: 1, W[(1, 2)]: 1, H[(2, 2)]: 1, W[(2, 0)]: 1}
    return Carrier("B0b* ring torus 3x3, gamma_C = the (1,1) DIAGONAL", 9,
                   b.edges, b.faces, dict(b.gF), gC)

def push_U(car):
    VF = car.loop_vertices(car.gF); VC = car.loop_vertices(car.gC)
    cnt = {(0,0):0,(1,0):0,(0,1):0,(1,1):0}
    for v in range(car.V):
        cnt[(1 if v in VF else 0, 1 if v in VC else 0)] += 1
    V = car.V
    return (Fr(cnt[(0,0)],V), Fr(cnt[(1,0)],V), Fr(cnt[(0,1)],V), Fr(cnt[(1,1)],V)), cnt

def has_zero(p):
    hi = (p[0] + p[1]) - (p[2] + p[3]); lo = abs(p[0] - p[1]) - abs(p[2] - p[3])
    return hi * lo <= 0, hi, lo

def factors(p):
    return p[0] * p[3] == p[1] * p[2]

def m_jensen_mp(p, n=1 << 16, splits=None):
    """Jensen reduction in mpmath.  If `splits` (branch-crossing angles) are given, integrate
    each smooth piece separately with Gauss-Legendre — spectral accuracy despite the kink."""
    p00, p10, p01, p11 = [mp.mpf(q.numerator) / q.denominator if isinstance(q, Fr) else mp.mpf(q)
                          for q in p]
    def f(t):
        c = mp.cos(t)
        A2 = p00 ** 2 + p10 ** 2 + 2 * p00 * p10 * c
        B2 = p01 ** 2 + p11 ** 2 + 2 * p01 * p11 * c
        return mp.log(mp.sqrt(A2 if A2 > B2 else B2))
    if splits:
        pts = [mp.mpf(0)] + [mp.mpf(s) for s in splits] + [mp.pi]
        tot = mp.mpf(0)
        for a, b in zip(pts[:-1], pts[1:]):
            tot += mp.quad(f, [a, b])
        return tot / mp.pi
    return mp.quad(f, [0, mp.pi]) / mp.pi

def crossing_angles(p):
    """cos t where |A|^2 = |B|^2 ; affine in cos t so at most one root."""
    p00, p10, p01, p11 = [float(q) for q in p]
    A0 = p00 ** 2 + p10 ** 2 - p01 ** 2 - p11 ** 2
    B0 = 2 * (p00 * p10 - p01 * p11)
    if B0 == 0:
        return []
    c = -A0 / B0
    return [float(np.arccos(c))] if -1 <= c <= 1 else []

# ============================================================ orbit machinery
def running_average_half_plus(p, da, db, Ns):
    """(1/N) sum_{k<=N} log|Z_k| for alpha = 1/2 + da, beta = 1/2 + db with da,db TINY.

    PRECISION DEFECT FOUND AND FIXED IN THIS LANE, RECORDED RATHER THAN PATCHED SILENTLY:
    computing np.mod(k*(0.5+da), 1.0) in float64 destroys the answer for da < 1e-10, because
    k*(0.5+da) has magnitude ~5e6 at k = 1e7, whose float64 ulp is ~1e-9 — larger than the
    quantity k*da being measured.  The first run of this table was float-noise for the last two
    rows.  The fixed form splits the phase EXACTLY:  k*alpha mod 1 = (k mod 2)/2 + k*da,
    where k is an exact integer and k*da < 1e-7 never wraps."""
    p00, p10, p01, p11 = [float(q) for q in p]
    Nmax = max(Ns); out = {}
    tot = 0.0; done = 0; CH = 10 ** 6
    while done < Nmax:
        n = min(CH, Nmax - done)
        ki = np.arange(done + 1, done + n + 1, dtype=np.int64)
        k = ki.astype(np.float64)
        fa = np.mod((ki % 2) * 0.5 + k * da, 1.0)
        fb = np.mod((ki % 2) * 0.5 + k * db, 1.0)
        x = np.exp(2j * np.pi * fa); y = np.exp(2j * np.pi * fb)
        a = np.abs(p00 + p10 * x + p01 * y + p11 * x * y)
        cs = np.cumsum(np.log(np.maximum(a, 1e-323)))
        for N in Ns:
            if done < N <= done + n:
                out[N] = (tot + cs[N - done - 1]) / N
        tot += float(cs[-1]); done += n
    return out

# ============================================================ main
if __name__ == "__main__":
    print("=" * 104)
    print("LEG 1 — THE BRIEF'S PREMISE, CORRECTED BEFORE IT IS USED.")
    print("=" * 104)
    print("""  'A four-class P HAS zeros on T^2 (W-09: the firing region is exactly 1/2)' conflates:
     (i)  {(f,c) : SOME non-negative p annihilates Z_1}   — W-09's set, measure 1/2 for four
          classes and 1/4 for three;  and
     (ii) {(x,y) in T^2 : P_p(x,y) = 0} for a FIXED p     — the set the integrability question
          is about, which is measure ZERO for every p (2 points, or one circle).
  (i) is the UNION over p of (ii).  A union of null sets is not null.  The question 'does the
  hypothesis bite?' is answered by (ii), i.e. by ONE INEQUALITY ON THE READY STATE:
        P has a zero on T^2  <=>  ((p00+p10)-(p01+p11)) * (|p00-p10|-|p01-p11|) <= 0.
  Measure of that set in the ready-state simplex: EXACTLY 1/4 for four classes and EXACTLY 1/4
  for three (C_03 sec 4).  FOUR-CLASS OCCUPANCY IS NOT THE OPERATIVE VARIABLE HERE.""")

    print("\n" + "=" * 104)
    print("LEG 2 — B0b*, THE SAME COMPLEX WITH ONE DESIGNATED LOOP MOVED.")
    print("=" * 104)
    b, bs = build_B0b(), build_B0b_star()
    rb, rs = b.report(), bs.report()
    print(f"  ARMS.  gamma_F identical in both: {sorted(b.gF.items()) == sorted(bs.gF.items())}")
    print(f"         gamma_C differs:  B0b {sorted(b.gF.keys())}->{sorted(b.gC.keys())}  "
          f"B0b* {sorted(bs.gC.keys())}   (byte-identical arms? "
          f"{sorted(b.gC.items()) == sorted(bs.gC.items())})")
    assert sorted(b.gC.items()) != sorted(bs.gC.items()), "ARMS BYTE-IDENTICAL — CONTROL VOID"
    for car, r in ((b, rb), (bs, rs)):
        p, cnt = push_U(car)
        z, hi, lo = has_zero(p)
        print(f"\n  {car.name}")
        print(f"     V={r['V']} E={r['E']} F={r['F']} chi={r['chi']} b1={r['b1']} b2={r['b2']}"
              f"   gF bounds {r['gF_bounds']}  gC bounds {r['gC_bounds']}  indep {r['independent']}")
        print(f"     class counts {{00:{cnt[(0,0)]}, 10:{cnt[(1,0)]}, 01:{cnt[(0,1)]}, 11:{cnt[(1,1)]}}}"
              f"   FOUR CLASSES OCCUPIED: {all(cnt[k] > 0 for k in cnt)}")
        print(f"     pi = {tuple(str(q) for q in p)}   factors: {factors(p)}"
              f"   ZEROS ON T^2: {z}   (hi={hi}, lo={lo})")
        sp = crossing_angles(p)
        mP = m_jensen_mp(p, splits=sp)
        print(f"     m(P) [mpmath, kink-split at t = {['%.6f'%s for s in sp]}] = {mp.nstr(mP, 20)}")
    print("""
  ONE VARIABLE MOVED: which non-bounding cycle is designated gamma_C on the SAME complex, with
  the SAME gamma_F.  Topology, incidence, faces, gauge count and four-class occupancy are all
  unchanged.  The zero set of P appears.  THE HYPOTHESIS IS TURNED ON BY THE DESIGNATED LOOP,
  NOT BY THE CARRIER AND NOT BY CLASS OCCUPANCY.""")

    print("\n" + "=" * 104)
    print("LEG 3 — A FAILURE ALREADY SITTING IN S4's PUBLISHED TABLE.  EXACT ARITHMETIC.")
    print("=" * 104)
    pC = (Fr(1,4), Fr(1,4), Fr(1,4), Fr(1,4))
    print("  S4:582 publishes lambda(SENSE C) = -1.386294361120 for BOTH four-class carriers,")
    print("  with the identity '(1+x+y+xy)/4 = (1+x)(1+y)/4, so lambda = log(1/4) EXACTLY'.")
    print("  At S1 sec6's OWN published connection W_F = -1, W_C = -i  (u = -1, v = -i):")
    from fractions import Fraction as F2
    for k in range(1, 9):
        # exact Gaussian-integer arithmetic: u^k in {1,-1}, v^k in {1,-i,-1,i}
        ue = (-1) ** k
        ver, vei = [(1, 0), (0, -1), (-1, 0), (0, 1)][k % 4]
        # Z_k = (1+u^k)(1+v^k)/4
        ar, ai = 1 + ue, 0
        br, bi = 1 + ver, vei
        zr = F2(ar * br - ai * bi, 4); zi = F2(ar * bi + ai * br, 4)
        print(f"     k={k}:  u^k={ue:+d}  v^k={ver:+d}{vei:+d}i   "
              f"Z_k = {zr} {'+' if zi>=0 else '-'} {abs(zi)}i   "
              f"|Z_k|^2 = {zr*zr+zi*zi}  {'ZERO' if zr==0 and zi==0 else ''}")
    print("""     Z_k = 0 EXACTLY unless k = 0 mod 4 — odd k kill it through (1+u^k) = 0, and
     k = 2 mod 4 through (1+v^k) = 0.  Three of every four cells annihilate.
     So Omega_N = 0 for every N >= 1 and (1/N) sum log|Z_k| = -infinity for every N >= 1.
     THE PUBLISHED FOUR-CLASS SENSE C RATE IS NOT THE RATE AT THE CORPUS'S OWN CONNECTION.
     float64 does not see it: C_04 reports a finite -28.088521944 because |Z_2| evaluates to
     4.13e-17 instead of 0.  A lane trusting float64 would publish a finite wrong number.""")

    print("\n" + "=" * 104)
    print("LEG 4 — THE FINITE-N FAILURE, ON A NON-FACTORING FOUR-CLASS P WITH AN EXACT m(P).")
    print("=" * 104)
    pX = (Fr(4,10), Fr(2,10), Fr(3,10), Fr(1,10))
    z, hi, lo = has_zero(pX)
    print(f"  pi = (4,2,3,1)/10   four classes occupied, class-level ready state (the sense S4")
    print(f"  uses for its SENSE C column).   factors: {factors(pX)}   zeros on T^2: {z}"
          f"  (hi={hi}, lo={lo})")
    print("""  ZERO LOCATED EXACTLY: at x = -1, A = p00-p10 = 1/5 and B = p01-p11 = 1/5, so
  y = -A/B = -1.  P(-1,-1) = 4/10 - 2/10 - 3/10 + 1/10 = 0.   ISOLATED torus zero at (-1,-1).
  AND m(P) IS EXACT ANYWAY:  |0.4+0.2x|^2 - |0.3+0.1x|^2 = 0.1(1+cos t) >= 0, so the FIRST
  Jensen branch dominates everywhere and m(P) = m(2/5 + (1/5)x) = log(2/5) EXACTLY.
  So 'has a torus zero' and 'has a closed form' are INDEPENDENT — a second, independent blow
  to S4:599's inference that a non-factoring four-term P forces quadrature.""")
    mX = mp.log(mp.mpf(2) / 5)
    print(f"  m(P) = log(2/5) = {mp.nstr(mX, 20)}    Jensen check (mpmath) = "
          f"{mp.nstr(m_jensen_mp(pX), 20)}")
    print("""
  THE CONNECTION.  alpha = 1/2 + delta*sqrt2, beta = 1/2 + delta*sqrt3.
  L = {(m,n) : m alpha + n beta in Z} = {0}:  m alpha + n beta = (m+n)/2 + delta(m sqrt2 +
  n sqrt3), and 1, sqrt2, sqrt3 are linearly independent over Q, so the only solution is
  m = n = 0.  THE ORBIT IS DENSE AND EQUIDISTRIBUTED IN T^2 — every hypothesis of the corpus's
  'orbit dense in T^2 therefore the average converges to m(P)' clause is satisfied.
  For every ODD k, (u^k,v^k) = (-1,-1) * (tiny rotation), i.e. the orbit sits at distance
  ~2 pi k delta from the zero.  ONE VARIABLE MOVES DOWN THE TABLE: delta.""")
    Ns = [10 ** 2, 4 * 10 ** 3, 10 ** 4, 10 ** 6, 10 ** 7]
    print(f"\n  {'delta':>10s} | " + " | ".join(f"N={N:>9d}" for N in Ns) + " |  m(P)")
    print("  " + "-" * 100)
    for dexp in (3, 6, 9, 12, 15):
        delta = 10.0 ** (-dexp)
        out = running_average_half_plus(pX, delta * np.sqrt(2.0), delta * np.sqrt(3.0), Ns)
        print(f"  1e-{dexp:<8d} | " + " | ".join(f"{out[N]:11.6f}" for N in Ns)
              + f" | {float(mX):.6f}")
    print("  " + "-" * 100)
    print("  CROSSOVER CHECK: the N at which the average first comes within 0.01 of m(P)")
    for dexp in (3, 4, 5, 6, 7):
        delta = 10.0 ** (-dexp)
        grid = [int(10 ** (e / 4)) for e in range(8, 33)]
        out = running_average_half_plus(pX, delta * np.sqrt(2.0), delta * np.sqrt(3.0), grid)
        first = next((N for N in grid if abs(out[N] - float(mX)) < 0.01), None)
        print(f"     delta = 1e-{dexp}   first N within 0.01 of m(P): {first}"
              f"   (1/delta = {int(1/delta)})   ratio N/(1/delta) = "
              f"{first*delta if first else float('nan'):.3f}")
    print("""
  READ THE 'N = 4000' COLUMN.  4000 is S3's own measurement window (the register's erratum
  against W-02 is about exactly that window).  At delta = 1e-9 the honest answer is log(2/5) =
  -0.916291 and the measurement returns -6.154476, a factor of 6.72.  The deficit grows without
  bound as delta -> 0 AT EVERY FIXED N.  There is therefore NO N at which lambda can be read
  off without a hypothesis on the connection — and no amount of extra sampling fixes it,
  because the failure is not statistical.""")

    print("\n" + "=" * 104)
    print("LEG 5 — THE LIMIT FAILURE, RIGOROUSLY, ON A FOUR-CLASS READY STATE.")
    print("=" * 104)
    print("""  CONSTRUCTION.  gamma = 2^-3 + 2^-b2 + 2^-b3 + ... with b2 = 256 and b_{j+1} = 2^(2^b_j);
  alpha = gamma/2; beta = 4^(1/3) mod 1.
     * alpha is LIOUVILLE, beta is algebraic of degree 3.  L = {0}: m alpha + n beta in Z with
       n != 0 would make beta Liouville, contradicting Roth's theorem; with n = 0 it would make
       alpha rational.  SO THE ORBIT IS DENSE AND EQUIDISTRIBUTED IN T^2.
     * q_j = 2^(b_j) satisfies || q_j alpha - 1/2 || ~ 2^(-b_{j+1}).
  READY STATE: SENSE C four-class (1,1,1,1)/4, so log|Z_k| = log|1+u^k| + log|1+v^k| - log 4.
  EVERY TERM IS <= 0 (|Z_k| <= 1, W-08), so for every j
        (1/q_j) sum_{k<=q_j} log|Z_k|   <=   log|Z_{q_j}| / q_j ,
  and no summation is needed to bound the average from above.""")
    # || q_j alpha - 1/2 || = eps_j ~ 2^{b_j - b_{j+1} - 1};  |1 + u^{q_j}| = 2 |sin(pi eps_j)|
    # computed as 2 sin(pi eps) directly — forming exp(2 pi i (1/2 + eps)) and adding 1 would
    # cancel to nothing at any working precision.
    for j, (bj, bnext) in enumerate([(3, 256), (256, None)]):
        q = 2 ** bj
        if bnext is None:
            print(f"     j=2: q_2 = 2^256 ~ 1.16e77,  b_3 = 2^(2^256),  eps_2 ~ 2^-(2^256).")
            print(f"          log|Z_{{q_2}}| <= log(2 pi eps_2) + log(1/2) ~ -(2^256) log 2, so")
            print(f"          (1/q_2) sum log|Z_k|  <=  -(2^256) log 2 / 2^256 ... and with")
            print(f"          b_3 = 2^(2^256) the bound is -2^(2^256) * 0.693 / 1.16e77:")
            print(f"          ASTRONOMICALLY NEGATIVE, and it is a BOUND, not a measurement.")
            continue
        eps = mp.mpf(2) ** (bj - bnext - 1)
        z1 = 2 * mp.sin(mp.pi * eps)                          # = |1 + u^q|
        bound = mp.log(z1 * 2 / 4) / q                        # log|Z_q| <= log(|1+u^q| * 2/4)
        print(f"     j={j+1}: q_1 = 2^{bj} = {q}   ||q alpha - 1/2|| = eps_1 ~ 2^-{bnext+1-bj}"
              f" = {mp.nstr(eps, 6)}")
        print(f"          |1 + u^q| = 2 sin(pi eps_1) = {mp.nstr(z1, 6)}"
              f"   log|Z_q| <= {mp.nstr(mp.log(z1/2), 10)}")
        print(f"          => (1/q) sum_{{k<=q}} log|Z_k|  <=  {mp.nstr(bound, 10)}"
              f"      while m(P) = {mp.nstr(-mp.log(4), 10)}"
              f"   (factor {mp.nstr(bound/(-mp.log(4)), 4)})")
    print("""
  THE FIRST BOUND ALREADY EXCEEDS m(P) BY A FACTOR OF ~16, AT N = 8, ON A DENSE ORBIT.
  The second is beyond any floating-point representation and needs none: it is a bound, not a
  measurement.  liminf_N (1/N) sum log|Z_k| = -infinity, while m(P) = -log 4 is finite.
  SO 'THE ORBIT IS DENSE IN T^2, THEREFORE THE AVERAGE CONVERGES TO m(P)' IS A NON-SEQUITUR ON
  A FOUR-CLASS READY STATE, exactly as M1 showed it is on K1's three-class one.""")
    sys.exit(0)
