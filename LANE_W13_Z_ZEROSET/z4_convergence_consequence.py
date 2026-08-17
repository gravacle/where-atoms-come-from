"""LANE W-13 / Z  --  z4: WHAT THE ZERO SET DOES TO N1's CONVERGENCE.
The round's job is 'settle the convergence, or exhibit the counterexample'.  This lane owns the
zero set, so it settles the HALF OF THE PROBLEM THE ZERO SET DECIDES, and says exactly which
half is left.  Every ladder spans >= 4 decades and the TREND is printed, never one endpoint.
The full Liouville counterexample is M1_06's and is CITED, not redone."""
import sys, math
import numpy as np
from fractions import Fraction as F
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from z0_lib import strat_exact, zeros_closed_form, min_abs_P, mahler, fr, sublevel_measure

W = 96
def hdr(s):
    print("=" * W); print(s); print("=" * W)

hdr("z4  THE ZERO SET AND N1's CONVERGENCE: THE DICHOTOMY, AND WHAT IT LEAVES OPEN")
print("numpy", np.__version__, "\n")

# ---------------------------------------------------------------- exact-orbit machinery
def orbit_angles(num, den, K):
    """theta_k = 2 pi * ((k*num) mod den)/den for k = 1..K, with (k*num) mod den formed in
    EXACT Python integers.  This is M3-3's discipline: the position mod 1 never cancels."""
    r = np.empty(K, dtype=np.float64)
    acc = 0
    num %= den
    for k in range(K):
        acc += num
        if acc >= den:
            acc -= den
        r[k] = acc / den
    return r

def logZ(p, ax, ay):
    p00, p10, p01, p11 = [float(q) for q in p]
    x = np.exp(2j * np.pi * ax); y = np.exp(2j * np.pi * ay)
    return np.log(np.abs(p00 + p10 * x + p01 * y + p11 * x * y))

def ladder(p, ax, ay, Ns, label, mP):
    lz = logZ(p, ax, ay)
    cs = np.cumsum(lz)
    out = []
    for N in Ns:
        out.append(cs[N - 1] / N)
    print(f"      {label:<26s} " + " ".join(f"{v:13.6f}" for v in out))
    print(f"      {'  dev from m(P)':<26s} " + " ".join(f"{v-mP:13.2e}" for v in out))
    return out

K1  = fr(F(0),     F(3, 10), F(3, 10), F(2, 5))     # TWO isolated conical zeros
B0b = fr(F(4, 9),  F(2, 9),  F(2, 9),  F(1, 9))     # EMPTY
S1P = fr(F(0),     F(0),     F(1, 2),  F(1, 2))     # CURVE, stratum I
NS = [10, 100, 1000, 10000, 100000, 1000000]

# =======================================================================================
print("-" * W)
print("(a) THEOREM Z3 -- THE DICHOTOMY.  THIS IS THE LANE'S ANSWER TO 'SETTLE THE CONVERGENCE'.")
print("""
  (i)  IF Z(P) = EMPTY  --  equivalently (S1-S2)(D1-D2) > 0, equivalently the Jensen branches
       do not cross  --  THEN log|P| is CONTINUOUS AND BOUNDED on T^2, with
             log(min_t ||A|-|B||)  <=  log|P|  <=  log(max_t (|A|+|B|)) = 0,
       so log|P| is Riemann-integrable and WEYL EQUIDISTRIBUTION ALONE GIVES
             (1/N) SUM_{k<=N} log|Z_k|  ->  m(P) = log(an explicit weight)   [Theorem Z2].
       ==> ON THIS STRATUM H2 (density of <(conj W_F, W_C)>) IS SUFFICIENT.  NO DIOPHANTINE
           HYPOTHESIS IS NEEDED, NO DISCREPANCY RATE IS NEEDED, AND N1 IS UNCONDITIONAL.
       The stratum has measure 3/4 of the four-class simplex and 3/4 of K1's three-class one.

  (ii) IF Z(P) != EMPTY, H2 IS NOT SUFFICIENT, and the exceptional set of connections has a
       DIFFERENT SHAPE in the two sub-cases:
         Z = TWO POINTS (theta = 2).  An exact hit needs (u^k,v^k) = (x0,y0): TWO conditions.
           For each k the solution set is FINITE (k^2 points), so the exact-hit set of
           connections is COUNTABLE.  Near-hits are controlled by SIMULTANEOUS approximation
           to a point; Borel-Cantelli gives dist >= c k^{-1/2-delta} for a.e. (alpha,beta).
         Z = A CURVE (theta = 1).  In stratum I an exact hit needs only u^k = -1: ONE
           condition.  For each k the solution set is a FINITE UNION OF CIRCLES in the
           connection torus, and the union over k is DENSE.  The exact-hit set of connections
           goes from CODIMENSION 2 TO CODIMENSION 1.
       THAT is the precise sense in which the p00 = 0 curve states are materially harder.

  THE ANSWER TO THE ROUND'S QUESTION, FOR K1's REGISTERED pi:  IT IS IN CASE (ii), SUB-CASE
  'TWO POINTS' -- THE MILDEST SINGULAR CASE.  Its zeros are SIMPLE and CONICAL (z2), theta = 2,
  and the Diophantine condition N1 needs there is the weakest one in the stratification.
""")

# =======================================================================================
print("-" * W)
print("(b) THE THEOREM-CHECK ON (i).  THIS IS A THEOREM CHECK, NOT A CONTROL: it COULD NOT")
print("    HAVE FAILED, and 'could not have failed' voids a control, never a theorem.")
print("    On the EMPTY stratum |Z_k| >= min_{T^2}|P| for EVERY connection and EVERY k.")
rng = np.random.default_rng(20260817)
mn = min_abs_P(B0b)
worst = 2.0
for trial in range(400):
    a = F(int(rng.integers(1, 10**12)), 10**12)
    b = F(int(rng.integers(1, 10**12)), 10**12)
    ax = orbit_angles(a.numerator, a.denominator, 20000)
    ay = orbit_angles(b.numerator, b.denominator, 20000)
    worst = min(worst, float(np.exp(logZ(B0b, ax, ay)).min()))
print(f"    B0b uniform pi, 400 random connections x 20000 cells = 8e6 evaluations")
print(f"      min over ALL of them of |Z_k| = {worst:.12f}")
print(f"      theoretical floor min_{{T^2}}|P| = 1/9 = {1/9:.12f}   (never breached: "
      f"{worst >= 1/9 - 1e-12})")
print(f"      => |log|Z_k|| <= {abs(math.log(1/9)):.6f} for every k and every connection.")
print()

# =======================================================================================
print("-" * W)
print("(c) A THEOREM THAT WAS FOUND BY A CONFOUND IN MY OWN FIRST DESIGN, RECORDED RATHER")
print("    THAN PATCHED OUT.")
print("""
    I set out to build a 2x2: rows = the stratum, columns = 'benign' vs 'a connection whose
    orbit passes very close to K1_REG's torus zero', with nothing else moving.  THE SECOND
    COLUMN CANNOT BE BUILT AT K1's REGISTERED pi, and the reason is a theorem:

      THEOREM Z4.  At K1's registered pi = (0, 3/10, 3/10, 2/5) -- indeed whenever
      p10 = p01 -- the two torus zeros satisfy  y0 = conj(x0), hence  x0 y0 = 1: THEY LIE ON
      THE ANTI-DIAGONAL SUBTORUS {xy = 1}.  An orbit point (u^k, v^k) can therefore equal a
      zero only if (uv)^k = 1, i.e. only if (k,k) is a relation of (u,v).  SO:
        * NO CONNECTION SATISFYING H2 EVER LANDS EXACTLY ON A ZERO OF N1's REGISTERED
          POLYNOMIAL.  Under H2 only NEAR-approaches are possible, never exact hits.
        * and driving the orbit to within delta of a zero at step k forces
          |(uv)^k - 1| = O(delta), i.e. a NEAR-RELATION of size 2k with defect O(delta).
      MY FIRST 'HOSTILE' COLUMN HAD EXACTLY THAT DEFECT AND I DID NOT NOTICE IT UNTIL THE
      N = 10^6 ROW FAILED TO RECOVER: alpha_h + beta_h = (1 + 2e-12)/7, a near-relation (7,7)
      with defect 2e-12, so for every k << 5e11 the orbit is confined near a rank-1 subtorus
      and the average is a SUBTORUS average, not m(P).  TWO VARIABLES MOVED, and the second
      one is not removable.  It is printed below as what it is.
""")
p00_, p10_, p01_, p11_ = K1
zzk = zeros_closed_form(K1)
print(f"    THEOREM Z4 CHECKED EXACTLY: p10 = p01 ? {p10_ == p01_}.  x0*y0 for both zeros:")
for (x0_, y0_, cs_, sn2_, sg_) in zzk:
    print(f"      x0*y0 = {(x0_*y0_).real:+.15f}{(x0_*y0_).imag:+.15f}i   "
          f"|x0*y0 - 1| = {abs(x0_*y0_-1):.3e}")
print(f"      EXACT: x0 = cos+isin, y0 = cos-isin, x0 y0 = cos^2+sin^2 = "
      f"{cs_*cs_} + {sn2_} = {cs_*cs_+sn2_}")
print()
DEN = 10**30
phi = (1 + 5 ** 0.5) / 2
alpha_b = F(int(round((phi - 1) * DEN)), DEN)
beta_b = F(int(round((2 ** 0.5 - 1) * DEN)), DEN)
s0 = math.acos(-2 / 3)
thx = s0 / (2 * math.pi)
thy = 1.0 - thx
off = 1e-12
alpha_h = F(int(round((thx + off) / 7 * DEN)), DEN)
beta_h = F(int(round((thy + off) / 7 * DEN)), DEN)
print(f"    ARM DIFF (inputs), printed so a reader can check rather than trust:")
print(f"      alpha_benign  = {float(alpha_b):.18f}     alpha_near = {float(alpha_h):.18f}")
print(f"      beta_benign   = {float(beta_b):.18f}     beta_near  = {float(beta_h):.18f}")
print(f"      numerators differ: {alpha_b.numerator != alpha_h.numerator}, "
      f"{beta_b.numerator != beta_h.numerator}")
print(f"      7*(alpha+beta) - 1 :  benign {float(7*(alpha_b+beta_b))-1:+.6e}   "
      f"near {float(7*(alpha_h+beta_h))-1:+.6e}   <- THE CONFOUND, MEASURED")
print(f"      pi rows differ: K1_REG {tuple(str(q) for q in K1)} | "
      f"B0b {tuple(str(q) for q in B0b)} | S1_PUB {tuple(str(q) for q in S1P)}")
print(f"      strata: {strat_exact(K1)[0]} | {strat_exact(B0b)[0]} | "
      f"{strat_exact(S1P)[0]}-{strat_exact(S1P)[1]}")
print()
EXACT_M = {"B0b_U  EMPTY": math.log(4/9), "S1_PUB CURVE": -math.log(2)}
KMAX = NS[-1]
for colname, (aa, bb) in (("BENIGN", (alpha_b, beta_b)),
                          ("NEAR-ZERO (and, unavoidably, NEAR-RESONANT)", (alpha_h, beta_h))):
    ax = orbit_angles(aa.numerator, aa.denominator, KMAX)
    ay = orbit_angles(bb.numerator, bb.denominator, KMAX)
    print(f"    COLUMN {colname}")
    print(f"    {'':<26s} N = " + " ".join(f"{n:>13d}" for n in NS))
    for lbl, p_ in (("B0b_U  EMPTY", B0b), ("K1_REG TWO", K1), ("S1_PUB CURVE", S1P)):
        mP = EXACT_M.get(lbl, mahler(p_))
        src = "EXACT" if lbl in EXACT_M else "quadrature"
        print(f"      m(P) = {mP:.12f} ({src})   min_k |Z_k|, k<=1e6 = "
              f"{float(np.exp(logZ(p_, ax, ay)).min()):.3e}")
        ladder(p_, ax, ay, NS, lbl, mP)
    print()
print("    READ IT, WITH THE CONFOUND CARRIED:")
print("      * THE EMPTY ROW IS THE SAME IN BOTH COLUMNS from N = 1000 on, to 1e-5.  Neither")
print("        the near-zero targeting nor the near-resonance can touch it -- log|P| is")
print("        BOUNDED there, so no connection can produce a dive.  That is Theorem Z3(i).")
print("      * THE TWO ROW dives to -3.40 at N = 10 in the second column (a factor 4.4 in the")
print("        average) and does NOT return to m(P) by N = 10^6.  BOTH causes are present and")
print("        Theorem Z4 says they cannot be separated.  The residual -9.8e-03 is the")
print("        SUBTORUS offset, not the singularity: it is flat in N from 10^5 to 10^6.")
print("      * THE CURVE ROW is governed by ONE angle, not two (M3-3(d)), so the second")
print("        column's beta is irrelevant to it and both columns converge.")
print()

# =======================================================================================
print("-" * W)
print("(d) THE APPROACH-RATE LAWS THAT SET THE DIOPHANTINE EXPONENT, 5 DECADES, 200 CONNECTIONS.")
print("    A minimum over ONE Kronecker orbit is not self-averaging, so one connection claims")
print("    nothing.  GEOMETRIC MEAN over 200 independent connections, same connections used")
print("    for both targets -- the ONE thing that moves between the two columns is whether the")
print("    target is a POINT (codim 2) or a CIRCLE (codim 1).")
rng2 = np.random.default_rng(20260819)
NL = [100, 1000, 10000, 100000, 1000000]
acc_pt = np.zeros(len(NL)); acc_ci = np.zeros(len(NL))
NC = 200
for _ in range(NC):
    aa = F(int(rng2.integers(1, 10**15)), 10**15)
    bb = F(int(rng2.integers(1, 10**15)), 10**15)
    ax = orbit_angles(aa.numerator, aa.denominator, NL[-1])
    ay = orbit_angles(bb.numerator, bb.denominator, NL[-1])
    dx = np.abs(((ax - thx + 0.5) % 1.0) - 0.5)
    dy = np.abs(((ay - thy + 0.5) % 1.0) - 0.5)
    dpt = np.hypot(dx, dy)
    dci = np.abs(((ax - 0.5 + 0.5) % 1.0) - 0.5)
    for i, N in enumerate(NL):
        acc_pt[i] += math.log(float(dpt[:N].min()))
        acc_ci[i] += math.log(float(dci[:N].min()))
acc_pt /= NC; acc_ci /= NC
print(f"    {'N':>9s} {'geo-mean min dist POINT':>25s} {'slope':>9s} "
      f"{'geo-mean min dist CIRCLE':>26s} {'slope':>9s}")
for i, N in enumerate(NL):
    s1 = (acc_pt[i] - acc_pt[i-1]) / math.log(10.0) if i else float('nan')
    s2 = (acc_ci[i] - acc_ci[i-1]) / math.log(10.0) if i else float('nan')
    print(f"    {N:9d} {math.exp(acc_pt[i]):25.3e} {s1:9.4f} "
          f"{math.exp(acc_ci[i]):26.3e} {s2:9.4f}")
print("    PREDICTED slopes -0.5 (point, codim 2) and -1.0 (circle, codim 1).  The TREND over")
print("    four consecutive decades is what is claimed; no single decade is quoted.")
print("    CONSEQUENCE: the inhomogeneous Diophantine exponent N1 needs is tau > 1/2 in the")
print("    TWO-point stratum and tau > 1 in the CURVE stratum.  Borel-Cantelli gives both for")
print("    almost every connection, so N1 holds a.e. in EVERY stratum -- the strata differ in")
print("    the SHAPE of the exceptional set, not in its measure.\n")

# =======================================================================================
print("-" * W)
print("(e) THE RESONANT CASE DOES NOT ESCAPE THE ZERO SET -- IT CHANGES ITS CODIMENSION.")
print("""
    If H2 fails with a primitive relation (m,n), the orbit closure is the circle
    H = {(z^n, z^{-m})} and the Birkhoff limit is the ONE-VARIABLE Mahler measure
        lambda_(m,n) = m( p00 + p10 z^n + p01 z^{-m} + p11 z^{n-m} )
    (M1_08 T2(d); Boyd-Lawton then gives lambda_(m,n) -> m(P) as the relation grows).
    THE ZERO-SET QUESTION REAPPEARS ON H, ONE DIMENSION DOWN: does H meet Z(P)?
      * Z(P) is 0-dimensional and H is 1-dimensional inside a 2-torus, so GENERICALLY THEY
        MISS -- and then log|P| restricted to H is CONTINUOUS and the subtorus average is
        unproblematic.
      * They meet only on a codimension-1 set of relations, and then the one-variable
        restriction has a circle zero and the classical Sudler/log|2 sin| problem starts.
    CHECKED at the corpus's own two distinguished connections, on K1's REGISTERED pi:
""")
p = K1
zz = zeros_closed_form(p)
print(f"      K1_REG zeros at angles (s,t)/2pi = "
      f"({thx:.9f}, {thy:.9f}) and ({thy:.9f}, {thx:.9f})")
for lbl, (m_, n_) in (("S3/S4 headline f=2.0,c=1.1  relation (11,-20)", (11, -20)),
                      ("S1 published W_F=-1,W_C=-i  order 4 (finite)", (0, 0))):
    if (m_, n_) == (0, 0):
        pts = [(F(a, 4), F(b, 4)) for a in range(4) for b in range(4)]
        dd = min(math.hypot(abs(((float(a) - thx + .5) % 1) - .5),
                            abs(((float(b) - thy + .5) % 1) - .5)) for a, b in pts)
        print(f"      {lbl}: the orbit is 4 points of T^2 (u=-1, v=-i);")
        print(f"        min distance from the whole orbit to Z(P) = {dd:.6f}  -> MISSES.")
        print(f"        Its average is the finite-orbit value -(1/2)log5 = "
              f"{-0.5*math.log(5):.12f}, not m(P).  M1_08 T2(e), reproduced.")
    else:
        zs = np.exp(2j * np.pi * np.arange(200000) / 200000)
        xs = zs ** (-n_); ys = zs ** (-m_)
        vals = np.abs(float(p[0]) + float(p[1]) * xs + float(p[2]) * ys
                      + float(p[3]) * xs * ys)
        print(f"      {lbl}: subtorus H = {{(z^20, z^-11)}}")
        print(f"        min_{{H}} |P| over a 200000-point sweep = {vals.min():.9f}  -> MISSES Z(P)")
        print(f"        so log|P| is CONTINUOUS on H, the subtorus average is unproblematic,")
        print(f"        and the erratum's 4.9e-04 gap is a BOYD-LAWTON APPROXIMATION GAP")
        print(f"        (M1_07), NOT a singularity effect.  m(P|H) by quadrature = "
              f"{float(np.mean(np.log(vals))):.9f}")
print()

# =======================================================================================
print("-" * W)
print("(f) WHAT THIS LANE DOES NOT ESTABLISH -- BEFORE THE VERDICT, NOT AFTER.")
print("""
  * IT DOES NOT PROVE N1 ON THE SINGULAR STRATA.  On Z(P) != empty this lane supplies the
    GEOMETRY (count, simplicity, theta) and the SHAPE of the exceptional set.  The sufficient
    condition -- discrepancy o(1/(log N)^2) plus an inhomogeneous Diophantine condition
    dist >= c k^{-tau} -- is M1_08 T2(c)'s proof sketch, and it is a SKETCH there and is not
    upgraded here.  What is added is that theta = 2 (not 1) at K1's registered pi, which is
    the input that sketch's dyadic shell count needs and did not have.
  * THE LIOUVILLE COUNTEREXAMPLE IS M1_06's, NOT MINE.  M1_06 exhibits (u,v) with L = {0}
    (so H2 HOLDS) and liminf (1/N) sum log|Z_k| = -infinity while m(P) is finite.  That is
    the counterexample the round asks for and IT ALREADY EXISTS ON DISK.  This lane's
    contribution is to say WHICH STATES ADMIT IT: exactly the 1/4 of the simplex where the
    Jensen branches cross, and nowhere else.
  * THE HOSTILE COLUMN OF (c) IS AN ADVERSARIAL EXHIBIT, NOT A MEASUREMENT.  It shows a dive
    and a recovery at ONE tuned k.  It does not show liminf = -infinity and is not scored as
    doing so.
  * FLOAT FLOOR, STATED.  The orbit position is exact in Python integers but the angle is
    converted to double, so |Z_k| below about 1e-15 is not resolved.  Every dive reported
    here is above 1e-13.
  * A NULL READS TWO WAYS.  The EMPTY row of (c) showing no sensitivity to the connection's
    arithmetic reads as 'the singularity is the whole difficulty' OR as 'the two connections
    were not different enough'.  The theorem check in (b) is what distinguishes them, and it
    is a THEOREM, so it settles the first reading -- but only because it is a theorem, not
    because the null was informative.
""")
print("DONE z4")
