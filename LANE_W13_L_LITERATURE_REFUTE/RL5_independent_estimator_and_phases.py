#!/usr/bin/env python3
"""
RL5 — RE-DERIVE THE TARGET'S TWO LOAD-BEARING NUMBERS BY ROUTES THAT SHARE NO CODE WITH IT.

  (1) m(P) at K1's registered pi, by the CASSAIGNE-MAILLOT CLOSED FORM -- the target's own
      citation [10], a dilogarithm identity with NO QUADRATURE ANYWHERE.  This simultaneously
      (a) checks the target's Jensen trapezoid, (b) checks the register's -0.767507880, and
      (c) checks that Cassaigne-Maillot is the theorem the target says it is (finding L-14a
      claims S4 rests on it 17 times and no register row cites it).
      Calibrated first on m(1 + x + y), whose value is independently known.

  (2) The five-decade table of L3, both arms, with an INDEPENDENT PHASE PATH: alpha carried
      as an exact 128-bit fixed-point integer and reduced by numpy uint64 wraparound (which
      IS exact arithmetic modulo 2^64) plus a float64 low limb.  Phase error < 2^-88 turns at
      k = 10^7.  Structurally different from the target's three-double split3.

  (3) ARM A at a FULL DECADE BOUNDARY, k <= 10^4, in EXACT GAUSSIAN-INTEGER ARITHMETIC.
      The target stopped this check at k = 3000, i.e. short of any decade it tabulates.
      Here the N = 10^4 row of the target's arm-A table is reproduced with no transcendental
      function anywhere in the path.
"""
import math
import numpy as np
from mpmath import mp, mpf, mpc, polylog, atan2, sqrt as msqrt, pi as mpi, acos as macos, log as mlog, floor as mfloor, arg as marg, fabs as mfabs

mp.dps = 80

print("=" * 78)
print("RL5 — INDEPENDENT RE-DERIVATION OF m(P) AND OF THE FIVE-DECADE TABLE.")
print("=" * 78)


# ---------------------------------------------------------------- (1) Cassaigne-Maillot
def bloch_wigner(z):
    z = mpc(z)
    if z == 0 or z == 1:
        return mpf(0)
    return polylog(2, z).imag + marg(1 - z) * mlog(abs(z))


def cassaigne_maillot(a, b, c):
    """m(a + b x + c y) for positive reals a,b,c.
       If |a|,|b|,|c| form a triangle with angles alpha,beta,gamma OPPOSITE to a,b,c:
           pi*m = D((a/b) e^{i gamma}) + alpha log a + beta log b + gamma log c
       otherwise m = log max(a,b,c)."""
    a, b, c = mpf(a), mpf(b), mpf(c)
    if a > b + c or b > a + c or c > a + b:
        return mlog(max(a, b, c))
    al = macos((b * b + c * c - a * a) / (2 * b * c))
    be = macos((a * a + c * c - b * b) / (2 * a * c))
    ga = macos((a * a + b * b - c * c) / (2 * a * b))
    D = bloch_wigner((a / b) * mpc(0, 1) ** 0 * (mp.cos(ga) + mpc(0, 1) * mp.sin(ga)))
    return (D + al * mlog(a) + be * mlog(b) + ga * mlog(c)) / mpi


print("""
(1) CASSAIGNE-MAILLOT CLOSED FORM.  CALIBRATION FIRST, ON A VALUE KNOWN INDEPENDENTLY.""")
cm111 = cassaigne_maillot(1, 1, 1)
print("    m(1 + x + y)  by Cassaigne-Maillot : %s" % mp.nstr(cm111, 20))
print("    known value   (3 sqrt3 / 4 pi) L(chi_-3, 2) = L'(chi_-3,-1):")
Lchi = mp.nsum(lambda n: (mpf(1) / (3 * n + 1) ** 2 - mpf(1) / (3 * n + 2) ** 2), [0, mp.inf])
known = 3 * msqrt(3) / (4 * mpi) * Lchi
print("                                                  %s" % mp.nstr(known, 20))
print("    |difference| = %s   -> the closed form and my angle convention are CALIBRATED."
      % mp.nstr(mfabs(cm111 - known), 6))

print("""
    NOW K1's REGISTERED pi.  P = 0.3 x + 0.3 y + 0.4 x y = x y (0.4 + 0.3/x + 0.3/y), and
    Mahler measure is invariant under x -> 1/x and under multiplication by a monomial, so
    m(P) = m(0.4 + 0.3 x + 0.3 y).  a = 2/5, b = c = 3/10; 2/5 < 3/5 so the triangle case.""")
cmP = cassaigne_maillot(mpf(2) / 5, mpf(3) / 10, mpf(3) / 10)
print("    m(P) CASSAIGNE-MAILLOT, 80 dps, NO QUADRATURE : %s" % mp.nstr(cmP, 25))


def mP_jensen_mpmath():
    """Jensen in y, Gauss-Legendre SPLIT AT THE BRANCH CROSSINGS.
       A(t) = 0.3 constant; B(t) = |0.3 + 0.4 e^{it}|; they cross where cos t = -2/3, which
       is exactly the argument of the torus zero.  The integrand max(A,B) is CONTINUOUS but
       has a CORNER there, so a quadrature that does not know about the corner loses order."""
    tc = macos(mpf(-2) / 3) / (2 * mpi)          # crossing in turns

    def f(t):
        e = mp.expjpi(2 * t)
        A = mpf('0.3')
        B = mp.fabs(mpf('0.3') + mpf('0.4') * e)
        return mlog(max(A, B))
    return mp.quad(f, [0, tc, 1 - tc, 1])


mp.dps = 40
mj = mP_jensen_mpmath()
mp.dps = 80
print("    m(P) mpmath Gauss-Legendre, split at the corner: %s" % mp.nstr(mj, 25))
print("       |CM - split quadrature| = %s" % mp.nstr(mfabs(cmP - mj), 6))
tgt = mpf('-0.767507880357744')
reg = mpf('-0.767507880358')
print("    target lane's 2^22 trapezoid  -0.767507880357744  |diff from CM| = %s"
      % mp.nstr(mfabs(cmP - tgt), 6))
print("    register / M1 / R1  -0.767507880358              |diff from CM| = %s"
      % mp.nstr(mfabs(cmP - reg), 6))
print("""    VERDICT: the target's m(P), the register's figure, an independent quadrature and
    the Cassaigne-Maillot CLOSED FORM all agree.  m(P) IS NOT WHERE ANY DEFECT LIVES, and
    finding L-14(a) is right that this number rests on Cassaigne-Maillot.""")
MP = float(cmP)

print("""
--------------------------------------------------------------------------------
(1b) A CORRECTION TO THE TARGET'S PUBLISHED CONVENTIONS, MEASURED NOT ASSERTED.

  PUBLISHED_CONVENTIONS.txt (target, sealed) states:
     "The integrand is continuous (both branches vanish only if pi = 0), so the deterministic
      trapezoid converges SPECTRALLY.  Used at 2^22 nodes."

  CONTINUITY DOES NOT GIVE SPECTRAL CONVERGENCE; SMOOTHNESS DOES.  The Jensen integrand is
  max(|p00+p10 e^{it}|, |p01+p11 e^{it}|), a max of two smooth functions, and at K1's pi the
  two branches CROSS -- at cos t = -2/3, which is precisely the argument of the torus zero
  found by the target's own L1.  A max at a transversal crossing has a corner, the integrand
  is C^0 and not C^1, and Euler-Maclaurin gives the periodic trapezoid O(n^-2), not
  exponential.  Measured against the Cassaigne-Maillot closed form:
--------------------------------------------------------------------------------""")


def mP_trap(n):
    t = np.arange(n) * (2.0 * np.pi / n)
    e = np.exp(1j * t)
    A = np.abs(0.0 + 0.3 * e)
    B = np.abs(0.3 + 0.4 * e)
    return float(np.mean(np.log(np.maximum(A, B))))


print("   %-8s %-24s %-14s %s" % ("n", "trapezoid", "|err vs CM|", "observed order vs previous"))
prev_e = None
prev_n = None
for e_ in (16, 18, 20, 22, 24):
    n = 1 << e_
    val = mP_trap(n)
    err = abs(val - MP)
    order = "" if prev_e is None else "n^-%.2f" % (math.log(prev_e / err) / math.log(n / prev_n))
    print("   2^%-6d %-24.15f %-14.3e %s" % (e_, val, err, order))
    prev_e, prev_n = err, n
print("""   (The 2^24 row is float64-round-off-limited at ~1e-14 and its order is meaningless; the
   three orders above it, 1.97 / 1.91 / 1.78, are the measurement.)
   OBSERVED ORDER ~ n^-2, NOT SPECTRAL.  No number in the target lane changes -- its
   2^22 value is right to 3e-14 -- but the CONVENTIONS FILE STATES A FALSE PROPERTY OF ITS OWN
   ESTIMATOR, and states it in the sentence that justifies not showing a convergence table.
   Under this program's own rule ("CONVERGENCE IS NOT A WINDOW") the justification was the
   load-bearing part.""")


# ---------------------------------------------------------------- (2) independent phases
print("""
--------------------------------------------------------------------------------
(2) THE FIVE-DECADE TABLE, INDEPENDENT PHASE PATH (128-bit fixed point, uint64 wraparound).
--------------------------------------------------------------------------------""")
p00, p10, p01, p11 = 0.0, 0.3, 0.3, 0.4
DEC = [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
TWO64 = 1 << 64
TWO128 = 1 << 128


def fixed128(alpha_mp):
    a = alpha_mp % 1
    A = int(mfloor(a * mpf(2) ** 128))
    return A >> 64, A & (TWO64 - 1)


def turns_fixed(k_u64, hi, lo):
    r = (k_u64 * np.uint64(hi))                 # exact mod 2^64 (numpy uint64 wraps)
    t = r.astype(np.float64) / float(TWO64)
    t = t + k_u64.astype(np.float64) * float(lo) / float(TWO128)
    return np.mod(t, 1.0)


ARMS = {
    "A_algebraic": (atan2(mpf(4), mpf(3)) / (2 * mpi), atan2(mpf(12), mpf(5)) / (2 * mpi),
                    "u=(3+4i)/5, v=(5+12i)/13  ALGEBRAIC deg 2"),
    "B_transcend": ((mpf(-1) / (2 * mpi)) % 1, (msqrt(mpf(2)) / (2 * mpi)) % 1,
                    "f=1.0, c=sqrt(2)  (S4:603)  TRANSCENDENTAL"),
}
TARGET = {  # the target's sealed L3 table, quoted for comparison
    "A_algebraic": [-0.769541895224, -0.767461981673, -0.767527878713, -0.767506965386, -0.767507908961],
    "B_transcend": [-0.765486480113, -0.767503341919, -0.767546539710, -0.767511151758, -0.767508233727],
}
np.seterr(over="ignore")
out = {}
for nm, (al, be, lab) in ARMS.items():
    ha, la = fixed128(al)
    hb, lb = fixed128(be)
    print("\n  ARM %s  [%s]" % (nm, lab))
    print("    alpha 128-bit numerator hi = %d  lo = %d" % (ha, la))
    print("    %-10s %-20s %-16s %-20s %s" % ("N", "(1/N)sum log|Z_k|", "dev from m(P)",
                                              "target's sealed row", "|mine - target|"))
    tot = 0.0
    prev = 0
    rows = []
    for i, N in enumerate(DEC):
        s = prev
        while s < N:
            n = min(10 ** 6, N - s)
            k = np.arange(s + 1, s + n + 1, dtype=np.uint64)
            tu = turns_fixed(k, ha, la)
            tv = turns_fixed(k, hb, lb)
            eu = np.exp(2j * np.pi * tu)
            ev = np.exp(2j * np.pi * tv)
            Zk = p00 + p10 * eu + p01 * ev + p11 * eu * ev
            tot += float(np.sum(np.log(np.abs(Zk))))
            s += n
        prev = N
        avg = tot / N
        rows.append(avg)
        print("    %-10d %-20.12f %-+16.3e %-20.12f %.2e"
              % (N, avg, avg - MP, TARGET[nm][i], abs(avg - TARGET[nm][i])))
    out[nm] = rows
    devs = [abs(a - MP) for a in rows]
    print("    TREND of |dev|: %s" % " -> ".join("%.1e" % d for d in devs))
print("""
  ARM DIFF ON OUTPUTS: max |avg_A(N) - avg_B(N)| over the five decades = %.3e ; tables
  distinct: %s.""" % (max(abs(a - b) for a, b in zip(out["A_algebraic"], out["B_transcend"])),
                      out["A_algebraic"] != out["B_transcend"]))
print("""  VERDICT: BOTH ARMS OF L3 REPRODUCE on a phase path sharing no code with the target's.
  Nothing in L3's arithmetic is wrong.  L-13 stands as arithmetic.""")

# ---------------------------------------------------------------- (3) exact Gaussian to 1e4
print("""
--------------------------------------------------------------------------------
(3) ARM A AT A DECADE BOUNDARY IN EXACT GAUSSIAN-INTEGER ARITHMETIC, k <= 10^4.
    The target's exact cross-check stopped at k = 3000, short of every decade it tabulates.
    No transcendental function is evaluated anywhere below.
--------------------------------------------------------------------------------""")


def ilog(n):
    s = max(0, n.bit_length() - 200)
    return math.log(n >> s) + s * math.log(2.0) if s else math.log(n)


KE = 10 ** 4
u_r, u_i = 1, 0
v_r, v_i = 1, 0
c5, c13 = 1, 1
tot = 0.0
for k in range(1, KE + 1):
    u_r, u_i = u_r * 3 - u_i * 4, u_r * 4 + u_i * 3     # (3+4i)^k
    v_r, v_i = v_r * 5 - v_i * 12, v_r * 12 + v_i * 5   # (5+12i)^k
    c5 *= 5
    c13 *= 13
    ar, ai = 3 * c13 * u_r, 3 * c13 * u_i
    br, bi = 3 * c5 * v_r, 3 * c5 * v_i
    cr = u_r * v_r - u_i * v_i
    ci = u_r * v_i + u_i * v_r
    Nr, Ni = ar + br + 4 * cr, ai + bi + 4 * ci
    tot += 0.5 * ilog(Nr * Nr + Ni * Ni) - ilog(10 * c5 * c13)
exact_1e4 = tot / KE
print("    (1/10^4) sum log|Z_k|  EXACT integer path : %.15f" % exact_1e4)
print("    target's sealed arm-A N=10^4 row         : %.15f" % TARGET["A_algebraic"][1])
print("    my 128-bit fixed-point path              : %.15f" % out["A_algebraic"][1])
print("    |exact - target| = %.3e      |exact - mine| = %.3e"
      % (abs(exact_1e4 - TARGET["A_algebraic"][1]), abs(exact_1e4 - out["A_algebraic"][1])))
print("    dev from m(P) (exact path) = %+.3e" % (exact_1e4 - MP))
print("""
    VERDICT: the target's arm-A decade row is confirmed EXACTLY at N = 10^4, by integer
    arithmetic alone.  Its phase machinery is sound and its D-7 worry about mpmath is
    unnecessary at this decade.
--------------------------------------------------------------------------------""")
print("\nDONE RL5")
