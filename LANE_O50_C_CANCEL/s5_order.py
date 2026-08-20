"""
O-50-C  STEP 5.   PART 2.  WHAT WOULD MAKE IT ACCUMULATE?

A sum of m terms with independent signs grows as sqrt(m).  A sum with correlated signs can grow
as m.  Gravity's source accumulates because mass has ONE SIGN.  This step measures the coherence
of a record-configuration functional as a function of the configuration's ORDER, by THREE
INDEPENDENT INSTRUMENTS, all exact:

  I-1  DEFECT SWEEP -- deterministic.  Start from the fully ordered configuration and flip q
       records.  Order parameter: the defect fraction q/m.  No measure, no ensemble: this is the
       value AT a configuration, which is the object Part 1 showed matters.
  I-2  BIASED MEASURE -- a product measure with P(record = +1) = p.  Order parameter:
       lambda = 2p-1 = E[s_i], a SYMMETRY-BREAKING bias.  Exact rational binomial sums.
  I-3  CORRELATION WITHOUT A FIELD -- a nearest-neighbour Ising measure on the records at coupling
       u = e^(2 beta), kept RATIONAL so the whole computation stays exact; zero field.  Order
       parameter: the nearest-neighbour correlation t = (u-1)/(u+1).  E[s_i] = 0 at every u.

I-3 IS THE INSTRUMENT THAT DECIDES THE QUESTION, and it is not a null control -- it is a live
instrument that produces arbitrarily strong correlations and still fails to accumulate.  If
"correlated signs" were sufficient, I-3 would show accumulation.  It does not.

INSERTED vs INDUCED.  All three order parameters are INSERTED BY HAND -- that is the point of the
step: it asks what WOULD make the quantity accumulate, so the ordering must be put in to see what
it buys.  Step 6 asks whether anything INDUCES it.  Nothing here is a property of the carrier.
"""
import sys, os, math
from fractions import Fraction
from decimal import Decimal, getcontext
from math import comb
getcontext().prec = 60

OUT = []
def P(s=""):
    OUT.append(str(s)); print(s)
BAR = "=" * 104
bar = "-" * 104
PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494")
SQ2PI = (Decimal(2) / PI).sqrt()

def dec(fr):
    return Decimal(fr.numerator) / Decimal(fr.denominator)

def coherence_typical(m):
    """EXACT coherence of the uniform (maximally disordered) measure, as a rational."""
    return Fraction(comb(m - 1, (m - 1) // 2), 1 << (m - 1))

def coherence_biased_exact(m, a, b):
    """EXACT E|sum s_i| / m under the product measure with P(+1) = a/(a+b), a,b positive integers.
       Pure INTEGER arithmetic with RUNNING powers and a running binomial, so the loop is O(m)
       big-integer multiplications rather than O(m) exponentiations."""
    if b == 0: return Fraction(1)
    num = 0
    ak = 1                      # a^k
    bk = b ** m                 # b^(m-k)
    c = 1                       # C(m,k)
    for k in range(m + 1):
        d = 2 * k - m
        if d: num += c * ak * bk * (d if d > 0 else -d)
        if k < m:
            ak *= a
            bk //= b
            c = c * (m - k) // (k + 1)
    return Fraction(num, (a + b) ** m * m)

def coherence_biased_decimal(m, lam):
    """the SAME quantity, evaluated at 60 digits, for m beyond the exact route's reach.
       Cross-checked against coherence_biased_exact at the largest common m."""
    p = (Decimal(1) + lam) / 2
    q = 1 - p
    mu = Decimal(m) * (2 * p - 1)
    sd = (Decimal(m) * 4 * p * q).sqrt() if p * q > 0 else Decimal(0)
    lo = max(0, int(float(m) * float(p) - 60 * float(sd) / 2 - 10))
    hi = min(m, int(float(m) * float(p) + 60 * float(sd) / 2 + 10))
    tot = Decimal(0)
    lp, lq = p.ln(), (q.ln() if q > 0 else Decimal(0))
    lgm = Decimal(math.lgamma(m + 1))
    for k in range(lo, hi + 1):
        d = abs(2 * k - m)
        if d == 0: continue
        lw = lgm - Decimal(math.lgamma(k + 1)) - Decimal(math.lgamma(m - k + 1)) \
             + Decimal(k) * lp + (Decimal(m - k) * lq if q > 0 else Decimal(0))
        tot += lw.exp() * Decimal(d)
    return tot / Decimal(m)

P(BAR)
P("O-50-C  STEP 5.   PART 2: WHAT WOULD MAKE IT ACCUMULATE?   THREE EXACT INSTRUMENTS.")
P(BAR)
P()
P("D-23 SCOPE: TORUS.  Carrier, records and writer group as established in steps 1-2.")
P()

# ================================================================ I-1
P(bar)
P("  I-1.  THE DEFECT SWEEP.  DETERMINISTIC, EXACT, NO ENSEMBLE.")
P(bar)
P("  F(s) = s_1 + ... + s_m evaluated AT a configuration with q records flipped away from ordered.")
P("  coherence(config) = |F| / m = |1 - 2q/m|, exactly.  The order parameter is the defect fraction.")
P()
fracs = (0.0, 0.05, 0.10, 0.25, 0.40, 0.45, 0.49, 0.50)
P(f"  {'m':>6} " + " ".join(f"{('q/m=' + f'{x:.2f}'):>12}" for x in fracs))
for m in (16, 64, 256, 1024, 4096, 65536):
    row = []
    for x in fracs:
        q = round(x * m)
        row.append(f"{float(abs(Fraction(m - 2 * q, m))):>12.6f}")
    P(f"  {m:>6} " + " ".join(row))
P()
P("  and the TYPICAL configuration for comparison, exactly:")
P(f"  {'m':>7} {'typical defect fraction':>26} {'EXACT typical coherence':>26} {'sqrt(2/(pi m))':>17}")
for m in (16, 64, 256, 1024, 4096, 65536):
    c = coherence_typical(m)
    P(f"  {m:>7} {('1/2 +- ' + f'{0.5 / (m ** 0.5):.6f}'):>26} {dec(c):>26.16f} "
      f"{(SQ2PI / Decimal(m).sqrt()):>17.12f}")
P()
P("  READ: coherence IS the order parameter -- at any configuration it equals |1 - 2q/m| exactly,")
P("  with NO m-dependence and no asymptotics involved.  A configuration held at a FIXED defect")
P("  fraction below 1/2 has coherence bounded away from zero at every m, so F is EXTENSIVE there:")
P("  at q/m = 0.45 the coherence is 0.125000 at m = 16 and 0.100006 at m = 65536 -- the difference is")
P("  entirely the rounding of q = 0.45m to an integer (7 of 16), not an m-dependence of the law; the")
P("  q/m = 0.25 column is 0.500000 at every m because 0.25m is an integer there.  The typical config")
P("  sits at defect fraction 1/2 to within m^(-1/2), which is exactly why its coherence is m^(-1/2).")
P()
P("  ==> THE m^(-1/2) LAW IS A STATEMENT ABOUT WHERE THE TYPICAL CONFIGURATION SITS, NOT ABOUT WHAT")
P("      THE FUNCTIONAL CAN DO.  The same functional on the same carrier is fully extensive at an")
P("      ordered configuration.  Nothing about the carrier changed between those two columns.")

# ================================================================ I-2
P()
P(bar)
P("  I-2.  THE BIASED MEASURE.  ORDER PARAMETER lambda = 2p-1 = E[s_i].  EXACT RATIONALS.")
P(bar)
P("  A symmetry-breaking bias is INSERTED: each record independently reads +1 with probability p.")
P("  lambda = 0 is the unique writer-invariant measure of step 3; lambda = 1 is fully ordered.")
P()
lams = [(0, 1), (1, 64), (1, 32), (1, 16), (1, 8), (1, 4), (1, 2), (3, 4), (1, 1)]
P(f"  {'m':>6} " + " ".join(f"{('lam=' + f'{ln}/{ld}'):>11}" for ln, ld in lams))
for m in (16, 64, 256, 1024, 4096):
    row = []
    for ln, ld in lams:
        a, b = ld + ln, ld - ln            # p = a/(a+b) = (1+lam)/2
        c = coherence_biased_exact(m, a, b) if b else Fraction(1)
        row.append(f"{float(c):>11.7f}")
    P(f"  {m:>6} " + " ".join(row))
P()
P("  THE CROSSOVER, EXACTLY.  The bias wins over the fluctuation when lambda exceeds the disordered")
P("  coherence sqrt(2/(pi m)).")
P()
P(f"  {'m':>7} {'floor=sqrt(2/(pi m))':>22} {'coh at lam=floor/4':>20} {'coh at lam=floor':>18}"
  f" {'coh at lam=4*floor':>20} {'(4*floor)':>12}")
for m in (16, 64, 256, 1024, 4096):
    floor = SQ2PI / Decimal(m).sqrt()
    def rat(x, D=10 ** 6):
        n = int(x * D)
        return (D + n, D - n)
    r4, r1, r0 = rat(floor / 4), rat(floor), rat(floor * 4)
    P(f"  {m:>7} {floor:>22.12f} {float(coherence_biased_exact(m, *r4)):>20.8f}"
      f" {float(coherence_biased_exact(m, *r1)):>18.8f}"
      f" {float(coherence_biased_exact(m, *r0)):>20.8f} {float(floor * 4):>12.8f}")
P()
P("  READ: at lambda a quarter of the floor the coherence is still essentially the floor; at four")
P("  times the floor the coherence has become lambda itself.  THE CROSSOVER IS AT lambda ~ m^(-1/2).")
P("  An arbitrarily SMALL but m-INDEPENDENT bias converts the sqrt(m) residual into an extensive")
P("  one.  Accumulation does not require strong ordering -- it requires ANY ordering that does not")
P("  shrink with m.")
P()
P("  OUT-OF-SAMPLE (D-20).  Law read off m <= 256 only: coherence(m,lambda) = sqrt(lambda^2 +")
P("  2/(pi m)).  PREDICTED, then compared to the exact value at m far outside that range.")
P()
P(f"  {'m':>8} {'lambda':>10} {'lam/floor':>10} {'PREDICTED':>24} {'MEASURED':>24} {'rel. error':>12}"
  f" {'route':>18}")
for m in (4096, 16384, 262144):
    for ln, ld in [(0, 1), (1, 128), (1, 16), (1, 2)]:
        lam = Decimal(ln) / Decimal(ld)
        pred = (lam * lam + Decimal(2) / (PI * Decimal(m))).sqrt()
        if m <= 4096:
            a, b = ld + ln, ld - ln
            meas = dec(coherence_biased_exact(m, a, b)); route = "exact rational"
        else:
            meas = coherence_biased_decimal(m, lam); route = "60-digit decimal"
        floor = SQ2PI / Decimal(m).sqrt()
        P(f"  {m:>8} {f'{ln}/{ld}':>10} {float(lam / floor):>10.3f} {pred:>24.14f} {meas:>24.14f}"
          f" {float(abs(pred - meas) / meas):>12.2e} {route:>18}")
P()
# cross-check the two routes at a common m
for ln, ld in [(0, 1), (1, 16)]:
    a, b = ld + ln, ld - ln
    ex = dec(coherence_biased_exact(4096, a, b))
    dc = coherence_biased_decimal(4096, Decimal(ln) / Decimal(ld))
    P(f"  ROUTE CROSS-CHECK at m = 4096, lambda = {ln}/{ld}:  exact rational {ex:.16f}   "
      f"60-digit decimal {dc:.16f}   agree to {abs(ex - dc):.2e}")
P()
P("  READ, FROM THE TABLE AND NOT AROUND IT.  The relative error tracks lambda/floor and nothing")
P("  else, falling monotonically as lambda/floor moves away from 1 in either direction -- read off")
P("  the table: ratio 1.25 -> 9.7e-2 ;  5.0 -> 2.0e-2 ;  10.0 -> 5.0e-3 ;  40.1 -> 3.1e-4 ;")
P("  80.2 -> 7.8e-5 ;  320.8 -> 4.9e-6 ;  and ratio 0 (lambda = 0) -> 1e-6 to 6e-5.")
P("  So sqrt(lambda^2 + 2/(pi m)) is EXACT in the two limits and no better than a ~10% interpolation")
P("  AT the crossing itself.  D-20: reported, not hidden.  Nothing in the Part-2 conclusion depends")
P("  on the interpolation -- only on the two limits, which are exact.")
P("  The two independent evaluation routes agree to better than 1e-13 where they overlap.")

# ================================================================ I-3
P()
P(bar)
P("  I-3.  CORRELATION WITHOUT SYMMETRY BREAKING.  THE INSTRUMENT THAT DECIDES THE QUESTION.")
P(bar)
P("  Nearest-neighbour Ising measure on the records: weight proportional to u^(aligned bonds),")
P("  u = e^(2 beta) RATIONAL, zero field.  The measure is invariant under the GLOBAL flip, so")
P("  E[s_i] = 0 at every u: the signs are strongly CORRELATED and the symmetry is NOT broken.")
P()

def ising_hist(m):
    """histogram over configurations of (aligned bonds, total spin).  Built once per m."""
    h = {}
    for x in range(1 << m):
        s = [1 - 2 * ((x >> i) & 1) for i in range(m)]
        al = sum(1 for i in range(m - 1) if s[i] == s[i + 1])
        S = sum(s)
        h[(al, S)] = h.get((al, S), 0) + 1
    return h

def ising_exact(h, u):
    Z = 0; e1 = 0; e2 = 0; ea = 0
    for (al, S), c in h.items():
        w = c * u ** al
        Z += w; e1 += w * S; e2 += w * S * S; ea += w * abs(S)
    return Fraction(e1, Z), Fraction(e2, Z), Fraction(ea, Z)

def ising_var_sum(m, t):
    """EXACT Var(sum s) for the OPEN zero-field chain, as the literal double sum:
       <s_i s_j> = t^|i-j|, so Var = m + 2 * sum_{d=1}^{m-1} (m-d) t^d.  Exact rational.
       Only usable at small m: t^d has denominator den(t)^d."""
    v = Fraction(m); td = Fraction(1)
    for d in range(1, m):
        td *= t
        v += 2 * (m - d) * td
    return v

def ising_var_closed(m, t):
    """the SAME sum in closed form, in 60-digit Decimal, so large m is reachable:
       Var = m(1+t)/(1-t) - 2t(1 - t^m)/(1-t)^2 .   VERIFIED against ising_var_sum below."""
    td = Decimal(t.numerator) / Decimal(t.denominator)
    if td == 0: return Decimal(m)
    return Decimal(m) * (1 + td) / (1 - td) - 2 * td * (1 - td ** m) / (1 - td) ** 2

P(f"  {'m':>4} {'u=e^2b':>8} {'t=corr':>10} {'xi=-1/ln t':>11} {'E[sum s] EXACT':>15}"
  f" {'Var EXACT':>13} {'Var/m':>9} {'E|sum s| EXACT':>15} {'coherence':>11} {'coh*sqrt(m)':>12}"
  f" {'closed==brute':>14}")
for m in (16, 18):
    h = ising_hist(m)
    for u in (Fraction(1), Fraction(2), Fraction(4), Fraction(10), Fraction(100), Fraction(1000)):
        t = Fraction(u - 1, u + 1)
        e1, e2, ea = ising_exact(h, u)
        var = e2 - e1 * e1
        vt = ising_var_sum(m, t)
        xi = "inf" if t == 1 else ("0" if t == 0 else f"{-1 / math.log(float(t)):.3f}")
        P(f"  {m:>4} {str(u):>8} {float(t):>10.6f} {xi:>11} {str(e1):>15}"
          f" {float(var):>13.5f} {float(var) / m:>9.4f} {float(ea):>15.6f}"
          f" {float(ea) / m:>11.6f} {float(ea) / m * (m ** 0.5):>12.6f}"
          f" {str(var == vt):>14}")
P()
P("  The transfer-matrix variance and the brute enumeration of all 2^m configurations agree EXACTLY")
P("  as rationals in every row, so the large-m table below rests on a verified formula.")
P()
P("  THE CLOSED FORM FOR Var, VERIFIED AGAINST THE EXACT RATIONAL DOUBLE SUM BEFORE IT IS USED:")
P(f"  {'m':>5} {'t':>10} {'Var exact rational (double sum)':>34} {'Var closed form (60 digits)':>30} {'agree?':>8}")
for m in (16, 32, 64):
    for t in (Fraction(0), Fraction(1, 2), Fraction(9, 11), Fraction(99, 101)):
        vs = ising_var_sum(m, t); vc = ising_var_closed(m, t)
        ok = abs(dec(vs) - vc) < Decimal("1e-30") * (abs(vc) + 1)
        P(f"  {m:>5} {float(t):>10.6f} {float(vs):>34.14f} {vc:>30.14f} {str(ok):>8}")
P()
P("  LARGE m, from the VERIFIED closed form.  Entry is coherence * sqrt(m); a FLAT row means the")
P("  m^(-1/2) law still holds and only the constant changed.")
P()
P(f"  {'t = corr':>10} {'xi':>10} " + " ".join(f"{('m=' + str(m)):>13}" for m in
                                              (64, 256, 1024, 4096, 16384, 262144))
  + f" {'sqrt(2/pi)*sqrt((1+t)/(1-t))':>30}")
for t in (Fraction(0), Fraction(1, 2), Fraction(9, 11), Fraction(99, 101), Fraction(999, 1001)):
    row = []
    for m in (64, 256, 1024, 4096, 16384, 262144):
        var = ising_var_closed(m, t)
        ea = SQ2PI * var.sqrt()
        row.append(f"{float(ea / Decimal(m) * Decimal(m).sqrt()):>13.6f}")
    lim = SQ2PI * ((Decimal(1) + dec(t)) / (Decimal(1) - dec(t))).sqrt()
    xi = "0" if t == 0 else f"{-1 / math.log(float(t)):.2f}"
    P(f"  {float(t):>10.6f} {xi:>10} " + " ".join(row) + f" {float(lim):>30.6f}")
P()
P("  THE GAUSSIAN STEP USED IN THAT TABLE (E|X| = sqrt(2/pi)*sqrt(Var) for mean-zero X), CHECKED")
P("  AGAINST THE EXACT ENUMERATION -- INCLUDING WHERE IT FAILS:")
P(f"  {'m':>4} {'t':>10} {'E|sum s| EXACT':>16} {'sqrt(2/pi)*sqrt(Var)':>22} {'ratio':>10}")
for m in (16, 18):
    h = ising_hist(m)
    for u in (Fraction(1), Fraction(4), Fraction(100)):
        t = Fraction(u - 1, u + 1)
        e1, e2, ea = ising_exact(h, u)
        var = e2 - e1 * e1
        g = SQ2PI * (Decimal(var.numerator) / Decimal(var.denominator)).sqrt()
        P(f"  {m:>4} {float(t):>10.6f} {float(ea):>16.6f} {float(g):>22.6f} {float(dec(ea) / g):>10.5f}")
P()
P("  READ THE GAUSSIAN CHECK HONESTLY: the step is good to 1.5-3% when m is well above the")
P("  correlation length xi, and is off by 22% at t = 0.980 where xi = 50 EXCEEDS m = 16.  That is")
P("  the strongly-correlated small-m regime, and it is exactly where the large-m table's entries")
P("  are still rising rather than flat.  THE PART-2 CONCLUSION DOES NOT REST ON THIS STEP: it rests")
P("  on Var, which is exact at every m and every t (closed form verified against the double sum")
P("  above), and Var/m converges to the CONSTANT (1+t)/(1-t) -- printed in the Var/m column, which")
P("  is exact.  The Gaussian step only converts a variance into a mean absolute value for display.")
P()
P("  READ, AND IT IS THE FINDING OF PART 2:")
P("  E[sum s] is the EXACT RATIONAL 0 in every row, at every correlation strength.  The variance")
P("  grows to (1+t)/(1-t) times m, so coherence stays proportional to m^(-1/2) and only the")
P("  CONSTANT in front changes -- each row of the large-m table is FLAT in m once m exceeds the")
P("  correlation length xi, and converges to the predicted sqrt(2/pi)*sqrt((1+t)/(1-t)) in the last")
P("  column.  The rows that are NOT yet flat are exactly the rows where xi is still comparable to m.")
P()
P("  ==> CORRELATED SIGNS ARE NOT ENOUGH.  A measure can correlate the records arbitrarily strongly")
P("      and the signed sum still fails to accumulate, because correlation renormalises the")
P("      COEFFICIENT of sqrt(m), not the POWER.  What converts sqrt(m) into m is SYMMETRY BREAKING")
P("      -- a non-zero E[s_i] -- which is what I-2 inserts and I-3 refuses to.  'Mass has one sign'")
P("      is a statement about a BROKEN SYMMETRY, not about a correlation.")

# ================================================================ the two numbers asked for
P()
P(bar)
P("  THE TWO NUMBERS THE PROBE ASKS FOR, EXACTLY, AT SEVERAL m, WITH AN OUT-OF-SAMPLE CHECK")
P(bar)
P(f"  {'m':>8} {'ORDERED coherence':>20} {'TYPICAL coherence (exact rational)':>36}"
  f" {'ratio ordered/typical':>22} {'sqrt(pi m/2) predicted':>24}")
for m in (16, 64, 256, 1024, 4096, 16384, 65536):
    ct = coherence_typical(m)
    P(f"  {m:>8} {'1 exactly':>20} {dec(ct):>36.20f}"
      f" {float(1 / dec(ct)):>22.4f} {float((Decimal(m) * PI / 2).sqrt()):>24.4f}")
P()
P("  ORDERED coherence is EXACTLY 1 at every m -- not asymptotically, exactly, because |sum s_i| = m")
P("  at an ordered configuration.  TYPICAL coherence is the exact rational C(m,m/2)/2^m.  Their")
P("  RATIO is the exact factor by which ordering amplifies the functional, and it GROWS as")
P("  sqrt(pi m/2): the larger the carrier, the more ordering is worth.")
P()
P("  OUT-OF-SAMPLE: the ratio law sqrt(pi m/2) was read off m <= 256 and reproduces the exact ratio")
P("  at m = 65536 (last two columns), which was not used to fix it.")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "s5_order.txt"), "w") as f:
    f.write("\n".join(OUT) + "\n")
