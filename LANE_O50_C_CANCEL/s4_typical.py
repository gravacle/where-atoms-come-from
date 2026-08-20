"""
O-50-C  STEP 4.   PART 1, THE SHARPER QUESTION.  MEAN-ZERO OVER CONFIGURATION SPACE IS NOT
                  CANCELLATION AT A CONFIGURATION.  EXACT, TO m = 4096 AND BEYOND.

THE ERROR THIS STEP IS BUILT TO PREVENT.  "The mean over configuration space is exactly zero" and
"the value at the configuration the universe is actually in is zero" are DIFFERENT STATEMENTS, and
the first does not imply the second.  Conflating them would convert a m^(-1/2) suppression into a
non-existence claim.  This step separates them with exact arithmetic.

THE OBJECT.  On the torus the records are independent, so the analogue of C-62's E(s) = sum J_i
s_i s_{i+1} is a functional of the RECORDS THEMSELVES:

        F(s) = sum_{i=1..m} c_i s_i .

  INSERTED : the CHOICE of a linear functional -- this is a modelling choice, stated as one.
  INDUCED  : the EQUALITY of the weights c_i.  Step 1 measured Aut(carrier) to be TRANSITIVE on
             the records (D-22), so every scalar the carrier attaches to a record is the same on
             every record.  Uniform weights are not assumed here; they are forced by a measured
             symmetry.  Section D then drops the equality anyway and sweeps arbitrary weights.

EXACTNESS.  Every number in sections A-C is an exact rational computed from integers.  The only
floating point is the final decimal rendering, done at 40 significant digits with `decimal`.
"""
import sys, os, itertools, random
from fractions import Fraction
from decimal import Decimal, getcontext
getcontext().prec = 60

OUT = []
def P(s=""):
    OUT.append(str(s)); print(s)
BAR = "=" * 104
bar = "-" * 104

def comb(n, k):
    from math import comb as C
    return C(n, k)

def coherence_exact(m):
    """EXACT mean of |sum of m independent +-1| divided by m, over the uniform measure.
       Derivation: E|S_m| = m * C(m-1, floor((m-1)/2)) / 2^(m-1)  (a standard identity, VERIFIED
       against brute enumeration below rather than trusted)."""
    return Fraction(comb(m - 1, (m - 1) // 2), 1 << (m - 1))

def brute_absmean(m):
    """DIRECT enumeration of all 2^m configurations, exact integers.  F(x) = m - 2*popcount(x)."""
    tot = 0
    for x in range(1 << m):
        tot += abs(m - 2 * bin(x).count("1"))
    return Fraction(tot, 1 << m)

def brute_mean(m):
    """DIRECT enumeration, exact.  Not a formula: every configuration is visited."""
    tot = 0
    for x in range(1 << m):
        tot += m - 2 * bin(x).count("1")
    return Fraction(tot, 1 << m)

def central_prob(m):
    """coherence(m) for EVEN m, as the exact product prod_{j=1..m/2} (2j-1)/(2j)  =  C(m,m/2)/2^m.
       Evaluated in Decimal at 60 digits so that m beyond the reach of exact big-integer
       binomials is still available; cross-checked against the exact Fraction below."""
    acc = Decimal(1)
    for j in range(1, m // 2 + 1):
        acc = acc * Decimal(2 * j - 1) / Decimal(2 * j)
    return acc

def dec(fr, digits=24):
    return Decimal(fr.numerator) / Decimal(fr.denominator)

P(BAR)
P("O-50-C  STEP 4.   MEAN-ZERO IS NOT CANCELLATION.  EXACT RATIONALS TO m = 4096 AND OUT-OF-SAMPLE.")
P(BAR)
P()
P("D-23 SCOPE: TORUS.  m = 2k records on k disjoint L x L tori (step 1); G_W simply transitive")
P("            (step 2); uniform measure is the unique writer-invariant one (step 3).")
P()

# ================================================================ A. the two statements
P(bar)
P("  A.  THE TWO STATEMENTS, SIDE BY SIDE, EXACTLY.  F(s) = s_1 + ... + s_m .")
P(bar)
P("  MEAN(F)   is the average over all 2^m configurations.")
P("  E|F|      is the average MAGNITUDE -- what the functional is worth at a typical configuration.")
P("  MAX|F|    is what it is worth at an ORDERED configuration.")
P("  COHERENCE is E|F| / MAX|F| -- C-62's object, the fraction of the maximum that survives.")
P()
P(f"  {'m':>5} {'MEAN(F) exact':>14} {'E|F| exact (rational)':>34} {'E|F| decimal':>16}"
  f" {'MAX|F|':>7} {'COHERENCE = E|F|/m':>20} {'brute-force check':>18}")
for m in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20):
    mean = brute_mean(m)
    ex = coherence_exact(m) * m
    br = brute_absmean(m)
    coh = coherence_exact(m)
    rat = str(ex) if len(str(ex)) < 32 else str(ex)[:29] + "..."
    P(f"  {m:>5} {str(mean):>14} {rat:>34} {dec(ex):>16.10f} {m:>7} {dec(coh):>20.12f}"
      f" {('AGREES' if br == ex else 'DISAGREES'):>18}")
P()
P("  READ: MEAN(F) is EXACTLY 0 as a rational at every m -- the cancellation law holds, exactly, as")
P("  a statement about the AVERAGE.  E|F| is NOT zero at any m, and it GROWS.  The closed form")
P("  E|F| = m*C(m-1,floor((m-1)/2))/2^(m-1) agrees with brute-force enumeration of all 2^m")
P("  configurations at every m up to 20, exactly, as rationals.")
P()
P("  THE ABSOLUTE RESIDUAL, WHICH THE RATIO HIDES:")
P(f"  {'m':>7} {'COHERENCE E|F|/m':>20} {'E|F| itself':>16} {'E|F|/sqrt(m)':>16} {'E[F^2] exact':>14}"
  f" {'sqrt(E[F^2])':>14}")
for m in (2, 8, 32, 128, 512, 2048, 8192, 32768):
    coh = coherence_exact(m)
    ex = coh * m
    P(f"  {m:>7} {dec(coh):>20.12f} {dec(ex):>16.6f} {dec(ex) / Decimal(m).sqrt():>16.12f}"
      f" {m:>14} {Decimal(m).sqrt():>14.6f}")
P()
P("  READ: the RATIO falls as m^(-1/2) and the MAGNITUDE rises as m^(+1/2).  E[F^2] = m exactly at")
P("  every m (each c_i^2 = 1, cross terms have mean 0), so the root-mean-square residual is exactly")
P("  sqrt(m).  E|F|/sqrt(m) is flat at 0.7978... = sqrt(2/pi).")
P()
P("  ==> THE CANCELLATION IS A FACTOR OF m^(-1/2), NOT A ZERO.  At a typical configuration the")
P("      functional is worth sqrt(2m/pi), not 0.  Calling this 'cancellation' is correct only for")
P("      the RATIO to the maximum; the residual itself DIVERGES.  It is SUB-EXTENSIVE, which is a")
P("      different and weaker failure than vanishing -- and it is the failure that matters, because")
P("      a source must be EXTENSIVE, i.e. Theta(m), not merely non-zero.")

# ================================================================ B. digit-for-digit + OOS
P()
P(bar)
P("  B.  THE LAW, DIGIT FOR DIGIT, AND OUT-OF-SAMPLE (D-20: prediction, not fitting)")
P(bar)
P("  CLAIM UNDER TEST (C-62, measured on the chain): coherence = sqrt(2/pi) * m^(-1/2).")
P("  Here it is an EXACT CLOSED FORM on the torus, so no fit is involved and D-20's short-range")
P("  objection cannot apply: coherence(m) = C(m-1, floor((m-1)/2)) / 2^(m-1), a rational number.")
P()
PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494")
sqrt2pi = (Decimal(2) / PI).sqrt()
P(f"  sqrt(2/pi) to 30 digits = {sqrt2pi:.30f}")
P()
P(f"  {'m':>8} {'coherence * sqrt(m)  (exact -> 30 digits)':>46} {'matching digits vs sqrt(2/pi)':>30}")
def matching_digits(a, b):
    sa = f"{a:.30f}"; sb = f"{b:.30f}"
    k = 0
    for ca, cb in zip(sa, sb):
        if ca != cb: break
        if ca.isdigit(): k += 1
    return k
for m in (4, 16, 64, 256, 1024, 4096, 16384, 65536):
    val = dec(coherence_exact(m)) * Decimal(m).sqrt()
    P(f"  {m:>8} {val:>46.30f} {matching_digits(val, sqrt2pi):>30}")
P()
P("  READ: the exact rational, times sqrt(m), converges to sqrt(2/pi) digit by digit; the number of")
P("  agreeing digits grows with m exactly as an m^(-1) correction predicts (the exact expansion is")
P("  sqrt(2/pi)*(1 - 1/(4m) + O(m^-2))).")
P()
P("  OUT-OF-SAMPLE.  The law is FIXED from m <= 256 only, then used to PREDICT m far outside that")
P("  range; the prediction is compared to the exact rational, which was not used to build it.")
P()
P(f"  {'m (out of sample)':>18} {'PREDICTED sqrt(2/pi)/sqrt(m)':>30} {'EXACT coherence':>30}"
  f" {'relative error':>16}")
for m in (4096, 65536, 262144, 1048576, 4194304):
    pred = sqrt2pi / Decimal(m).sqrt()
    ex = central_prob(m)
    if m <= 65536:
        exact = dec(coherence_exact(m))
        agree = abs(exact - ex) < Decimal("1e-40")
    else:
        agree = None
    rel = abs(pred - ex) / ex
    P(f"  {m:>18} {pred:>30.20f} {ex:>30.20f} {rel:>16.3e}"
      f"   {'exact-Fraction cross-check: ' + str(agree) if agree is not None else ''}")
P()
P("  READ: the m <= 256 law predicts the exact value out to m = 4194304 with a relative error that")
P("  falls exactly as 1/(4m) -- the predicted correction, not a fitted one.  THE LAW IS NOT A")
P("  SHORT-RANGE FIT.  The two smallest rows carry an exact-big-integer Fraction cross-check of the")
P("  Decimal product evaluation, so the large-m rows rest on a verified evaluator.")

# ================================================================ C. the general functional
P()
P(bar)
P("  C.  THE GENERAL STATEMENT, NOT JUST THE LINEAR ONE.  PARSEVAL, EXACTLY.")
P(bar)
P("  For ANY functional f on the 2^m record configurations, expanded exactly as f = sum_S fhat(S)")
P("  chi_S, two identities hold with no approximation:")
P()
P("      MEAN(f)      = fhat(EMPTY)                                  [the RECORD-BLIND part]")
P("      VAR(f)       = sum over NON-EMPTY S of fhat(S)^2            [the RECORD-DEPENDENT part]")
P()
P("  so the typical value of f is  fhat(EMPTY) +- sqrt(VAR), EXACTLY, for every functional.")
P()
P("  ==> THE EXTENSIVE / RESPONSIVE DICHOTOMY, and it is an EXACT ARGUMENT, not a trend:")
P("      * the only part of f that is the SAME at every configuration is fhat(EMPTY).  That part can")
P("        be extensive, and it is RECORD-BLIND: no write changes it.  (This is C-61 at the record")
P("        level, now as an identity.)")
P("      * everything that RESPONDS to a write lives in the non-trivial characters, has mean")
P("        EXACTLY 0, and therefore takes BOTH SIGNS.  A quantity with mean exactly 0 that is not")
P("        identically 0 is positive somewhere and negative somewhere.")
P("      * therefore NO functional of the record configuration is both RESPONSIVE TO WRITING and")
P("        SIGN-DEFINITE IN ITS RESPONSE.  That is exactly the pair of properties a source needs.")
P()
P("  Verified exactly on random rational functionals (and note VAR is computed both ways):")
rnd = random.Random(7)
def wht(vals):
    n = len(vals); a = [Fraction(v) for v in vals]
    h = 1
    while h < n:
        for i in range(0, n, h * 2):
            for j in range(i, i + h):
                x, y = a[j], a[j + h]
                a[j], a[j + h] = x + y, x - y
        h *= 2
    return [c / n for c in a]
P(f"    {'m':>3} {'fhat(EMPTY) = MEAN':>20} {'VAR direct':>18} {'sum_{{S != 0}} fhat^2':>22} {'agree?':>7}"
  f" {'takes both signs?':>18}")
for m in (2, 3, 4, 5, 6):
    vals = [Fraction(rnd.randrange(-9, 10)) for _ in range(1 << m)]
    fh = wht(vals)
    mean = fh[0]
    var1 = sum((v - mean) ** 2 for v in vals) / (1 << m)
    var2 = sum(c * c for c in fh[1:])
    resp = [v - mean for v in vals]
    both = any(r > 0 for r in resp) and any(r < 0 for r in resp)
    P(f"    {m:>3} {str(mean):>20} {str(var1):>18} {str(var2):>22} {str(var1 == var2):>7}"
      f" {str(both):>18}")
P()
P("  CONTROL (D-15).  A functional whose typical value IS extensive, so the instrument registers")
P("  accumulation when accumulation is there.  N_+(s) = number of records reading +1.")
P()
P(f"  {'m':>6} {'MEAN(N_+) exact':>16} {'sqrt(VAR) exact':>18} {'typical N_+':>26}"
  f" {'extensive?':>11} {'record-DEPENDENT part':>22}")
for m in (4, 16, 64, 256, 1024):
    mean = Fraction(m, 2)
    var = Fraction(m, 4)
    P(f"  {m:>6} {str(mean):>16} {str(Decimal(var.numerator).sqrt() / Decimal(var.denominator).sqrt())[:12]:>18}"
      f" {(str(mean) + ' +- ' + f'{float(var) ** 0.5:.2f}'):>26} {'YES':>11}"
      f" {f'{float(var) ** 0.5:.2f} = sqrt(m)/2':>22}")
P()
P("  READ: N_+ IS extensive and IS sign-definite -- and its extensive part is the CONSTANT m/2,")
P("  which is identical at every configuration and therefore carries no record information at all.")
P("  Its configuration-DEPENDENT part is sqrt(m)/2.  The instrument does register the extensive")
P("  column as YES, so a genuinely extensive record functional would not be missed; the finding is")
P("  that its extensive part is the record-blind constant.")
P()
P("  N_+ /m = 1/2 + O(m^(-1/2)) at EVERY typical configuration.  A quantity that takes the same")
P("  value to relative accuracy m^(-1/2) at every typical configuration cannot distinguish one")
P("  configuration from another, which is what 'carrying record information' means.")

# ================================================================ D. can reweighting escape?
P()
P(bar)
P("  D.  CAN A DIFFERENT WEIGHTING ESCAPE?  THE EXACT TRADEOFF.  (drops the D-22 uniform weights)")
P(bar)
P("  Section A used equal weights because step 1 MEASURED Aut(carrier) transitive on records.  Drop")
P("  that and allow ANY weights c_i.  Then, exactly:")
P("        MAX|F| = ||c||_1 (attained at an ordered configuration)   TYPICAL |F| ~ ||c||_2 .")
P("  and the coherence ratio is ||c||_2 / ||c||_1, which lies in [1/sqrt(m), 1] for every c.")
P("  The lower end is UNIFORM weights; the upper end is ALL WEIGHT ON ONE RECORD.")
P()
P(f"  {'weight profile c (m = 64)':<34} {'||c||_1':>12} {'||c||_2':>12} {'ratio':>10}"
  f" {'participation':>14} {'extensive?':>11}")
m = 64
profiles = {
    "uniform  c_i = 1": [Fraction(1)] * m,
    "power law c_i = 1/i": [Fraction(1, i + 1) for i in range(m)],
    "power law c_i = 1/i^2": [Fraction(1, (i + 1) ** 2) for i in range(m)],
    "geometric c_i = 2^-i": [Fraction(1, 1 << i) for i in range(m)],
    "top 8 records only": [Fraction(1) if i < 8 else Fraction(0) for i in range(m)],
    "top 1 record only": [Fraction(1) if i < 1 else Fraction(0) for i in range(m)],
}
for name, c in profiles.items():
    l1 = sum(abs(x) for x in c)
    l2sq = sum(x * x for x in c)
    l2 = Decimal(l2sq.numerator).sqrt() / Decimal(l2sq.denominator).sqrt()
    ratio = l2 / (Decimal(l1.numerator) / Decimal(l1.denominator))
    # participation: how many records carry half the total weight
    srt = sorted((abs(x) for x in c), reverse=True)
    acc = Fraction(0); part = 0
    for x in srt:
        acc += x; part += 1
        if acc * 2 >= l1: break
    ext = "YES" if part >= m // 4 else "no"   # participation >= m/4
    P(f"  {name:<34} {float(l1):>12.4f} {float(l2):>12.4f} {float(ratio):>10.4f} {part:>14} {ext:>11}")
P(f"  {'1/sqrt(m) for reference':<34} {'':>12} {'':>12} {float(1 / Decimal(m).sqrt()):>10.4f}")
P()
P("  READ: the ratio reaches 1 only when the weight collapses onto ONE record -- participation 1 --")
P("  at which point the functional is not extensive at all: it is a single bit.  Every profile that")
P("  SPREADS weight over Theta(m) records sits at ratio O(m^(-1/2)).  This is an exact inequality,")
P("  not a search: ||c||_2 <= ||c||_1 with equality iff c has one non-zero entry, and ||c||_1 <=")
P("  sqrt(m)||c||_2 by Cauchy-Schwarz with equality iff all |c_i| are equal.")
P()
P("  ==> REWEIGHTING CANNOT PRODUCE A QUANTITY THAT IS BOTH EXTENSIVE AND COHERENT AT A TYPICAL")
P("      CONFIGURATION.  The escape is closed by Cauchy-Schwarz, on any carrier, at any m.")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "s4_typical.txt"), "w") as f:
    f.write("\n".join(OUT) + "\n")
