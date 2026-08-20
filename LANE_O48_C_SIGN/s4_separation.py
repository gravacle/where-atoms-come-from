"""
O-48-C  STEP 4.   CLAUSE (e): SEPARATION DEPENDENCE, AND WHO PUT IT THERE.

Clause (e) of the standard wants a SEPARATION-DEPENDENT quantity with a POWER-LAW falloff.
Three questions, each answered exactly:

  Q1  Does the correlation between two records fall off with their separation at all?
  Q2  Does the ENERGY COST of changing that correlation depend on the separation?
  Q3  If a power law does appear anywhere, was it INSERTED or INDUCED?

The arithmetic that makes Q2 exact:  Z_i Z_j = product of t_k over the window k in [i, j-1].
Changing the value of Z_i Z_j therefore means flipping an ODD NUMBER of bonds in that window
and nothing outside it.  Flipping bond k costs dE_k = -2 J_k t_k.  The minimum cost is the
minimum of sum_{k in S} dE_k over odd-sized subsets S of the window, which is computed here
in closed form in exact integers (sort the costs; take all the negative ones; if that count is
even, repair it with the single cheapest available swap).  Nothing is sampled.
"""
import itertools, math, random
from fractions import Fraction

OUT = []
def P(s=""):
    OUT.append(s); print(s)

def couplings(name, m, seed=0):
    rnd = random.Random(seed)
    if name == "uniform":   return [1] * m
    if name == "randpos":   return [rnd.randrange(1, 61) for _ in range(m)]
    if name == "randsign":  return [rnd.choice((1, -1)) * rnd.randrange(1, 61) for _ in range(m)]
    raise ValueError(name)

def min_odd_subset_sum(costs):
    """EXACT minimum of sum over odd-sized subsets of the integer list `costs`."""
    negs = [c for c in costs if c < 0]
    poss = [c for c in costs if c >= 0]
    base = sum(negs)
    if len(negs) % 2 == 1:
        best = base
    else:
        cand = []
        if poss: cand.append(base + min(poss))          # add the cheapest non-negative
        if negs: cand.append(base - max(negs))          # drop the least useful negative
        best = min(cand) if cand else None
    # brute-force cross-check on short lists -- the closed form must not be trusted alone
    if len(costs) <= 14:
        bf = None
        for r in range(1, len(costs) + 1, 2):
            for S in itertools.combinations(costs, r):
                v = sum(S)
                if bf is None or v < bf: bf = v
        assert bf == best, f"closed form {best} != brute force {bf} on {costs}"
    return best

P("=" * 104)
P("O-48-C  STEP 4.   SEPARATION DEPENDENCE, AND WHETHER IT WAS INSERTED")
P("=" * 104)
P()

# ------------------------------------------------------------------------ Q1: the correlation
P("-" * 104)
P("  Q1.  DOES THE RECORD-PAIR CORRELATION FALL OFF WITH SEPARATION?")
P("  <Z_i Z_j> is evaluated exactly on three ensembles: the GROUND eigenspace, a single generic")
P("  eigenspace, and the uniform ensemble over all 2^n configurations.")
P("-" * 104)
P(f"  {'family':<10} {'n':>4} {'ensemble':<22} " + " ".join(f"{'d='+str(d):>9}" for d in (1, 2, 4, 8, 16)))
for fam in ("uniform", "randpos", "randsign"):
    n = 33
    J = couplings(fam, n - 1)
    # ground: t_k = -sign(J_k)
    tg = [-1 if j > 0 else 1 for j in J]
    rnd = random.Random(11)
    tr = [rnd.choice((1, -1)) for _ in J]
    for label, t in (("GROUND eigenspace", tg), ("a generic eigenspace", tr)):
        cells = []
        for d in (1, 2, 4, 8, 16):
            v = 1
            for k in range(0, d): v *= t[k]
            cells.append(f"{v:>9}")
        P(f"  {fam:<10} {n:>4} {label:<22} " + " ".join(cells))
    # uniform ensemble over all configurations: exact
    cells = []
    for d in (1, 2, 4, 8, 16):
        # mean over all 2^(n-1) bond patterns of prod_{k<d} t_k is exactly 0 for d >= 1
        cells.append(f"{0:>9}")
    P(f"  {fam:<10} {n:>4} {'uniform over all 2^n':<22} " + " ".join(cells))
    P()
P("  READ: on ANY single eigenspace |<Z_iZ_j>| = 1 at EVERY separation -- perfectly FLAT, no decay")
P("  of any kind, neither power law nor exponential.  Averaged over the whole configuration space")
P("  it is exactly 0 at every separation -- also flat.  Neither is a falloff.  A 1D chain at fixed")
P("  energy is long-range ordered, so there is no length scale in the correlation to measure.")
P("  Clause (e) asks for separation DEPENDENCE; this quantity has none at all.")
P()

# ------------------------------------------------------- Q2: cost of changing the correlation
P("-" * 104)
P("  Q2.  DOES THE ENERGY COST OF CHANGING THE CORRELATION DEPEND ON THE SEPARATION?")
P("  cost(d) = min over odd subsets of the d-bond window of the exact integer energy change.")
P("  Evaluated from the GROUND configuration (every bond satisfied, so every single flip costs).")
P("-" * 104)
P(f"  {'family':<10} {'n':>4} " + " ".join(f"{'cost d='+str(d):>11}" for d in (1, 2, 4, 8, 16, 32))
  + f"   {'d-dependent?':>13}")
for fam in ("uniform", "randpos", "randsign"):
    for n in (33, 65):
        J = couplings(fam, n - 1)
        tg = [-1 if j > 0 else 1 for j in J]
        cells, vals = [], []
        for d in (1, 2, 4, 8, 16, 32):
            if d > len(J): cells.append(f"{'-':>11}"); continue
            costs = [-2 * J[k] * tg[k] for k in range(0, d)]
            c = min_odd_subset_sum(costs)
            cells.append(f"{c:>11}"); vals.append(c)
        P(f"  {fam:<10} {n:>4} " + " ".join(cells) +
          f"   {('YES' if len(set(vals)) > 1 else 'NO -- FLAT'):>13}")
    P()
P("  READ THE 'uniform' ROWS FIRST.  With every J_k equal the cost is the SAME NUMBER at every")
P("  separation -- exactly flat, no d-dependence whatsoever.  That is the control that decides")
P("  the meaning of the disordered rows: any d-dependence seen there comes from the SPREAD OF THE")
P("  COUPLINGS, which was INSERTED into the venue, and not from separation as such.")
P()

# ------------------------------------------------- what shape does the disordered case take?
P("-" * 104)
P("  WHAT SHAPE DOES THE DISORDERED CASE TAKE, AND IS IT A LAW?")
P("  From the ground configuration cost(d) = 2 * min_{k<d} |J_k|, a MINIMUM OVER d SAMPLES of the")
P("  inserted coupling distribution.  For J drawn uniformly on {1..K} the exact expectation of")
P("  that minimum is a closed form; it is computed below and compared to the measured mean.")
P("-" * 104)
K = 60
P(f"  {'d':>5} {'measured mean cost (300 draws)':>31} {'exact E[2*min of d draws]':>27} "
  f"{'ratio':>8} {'cost*d':>10}   {'CONTROL uniform J':>18}")
for d in (1, 2, 4, 8, 16, 32, 64, 128):
    rnd = random.Random(2024)
    tot = 0
    for _ in range(300):
        tot += 2 * min(rnd.randrange(1, K + 1) for _ in range(d))
    meas = Fraction(tot, 300)
    # exact: E[min] = sum_{v=1}^{K} P(min >= v) = sum_v ((K-v+1)/K)^d
    exact = 2 * sum(Fraction((K - v + 1) ** d, K ** d) for v in range(1, K + 1))
    P(f"  {d:>5} {float(meas):>31.4f} {float(exact):>27.6f} "
      f"{float(meas / exact) if exact else 0:>8.4f} {float(exact) * d:>10.4f}   {2:>18}")
P()
P("  READ: E[cost(d)] -> 2 as d grows (the floor of the inserted distribution), and cost*d is NOT")
P("  constant, so this is not a 1/d law either.  It is the extreme-value statistic of a coupling")
P("  distribution that was PUT IN BY HAND, and it saturates at the smallest coupling present.")
P("  It carries no information about the records: min_{k<d} |J_k| is a function of H alone --")
P("  the same record-blindness found in step 3, now in the separation channel.")
P()

# ------------------------------------------------------------ Q3: insert a power law, get one
P("-" * 104)
P("  Q3.  THE HONESTY CONTROL.  Insert a long-range coupling J_ij = round(A / |i-j|^p) and measure")
P("  the two-body energy at separation d.  If the measured exponent comes back equal to p, the")
P("  measurement was of the INPUT.  This is run so that no later reading of clause (e) on this")
P("  construction can mistake an inserted exponent for an induced one.")
P("-" * 104)
P(f"  {'inserted p':>11} {'measured exponent (log-log, exact ints)':>41} {'|recovered - inserted|':>24}")
A = 10 ** 9
devs = []
for p in (1.0, 1.5, 2.0, 3.0):
    ds = [1, 2, 4, 8, 16, 32]
    Js = [max(1, int(A / d ** p)) for d in ds]
    # slope of log J vs log d, from the two extreme points -- no fitting freedom
    slope = (math.log(Js[-1]) - math.log(Js[0])) / (math.log(ds[-1]) - math.log(ds[0]))
    devs.append(abs(-slope - p))
    P(f"  {p:>11.1f} {-slope:>41.6f} {abs(-slope - p):>24.2e}")
P()
P(f"  READ: the recovered exponent reproduces the inserted exponent to within {max(devs):.1e} in every")
P("  row -- the residual is integer rounding of A/d^p, nothing else.  This control is DELIBERATELY")
P("  trivial: that is the point.  A power law in this construction is exactly as strong as the")
P("  power law written into H, and no stronger.")
P("  NOTHING in the nearest-neighbour chain INDUCES a power law: its energy is contact-only, its")
P("  correlations are flat, and the only separation dependence anywhere in the lane came from the")
P("  extreme-value statistics of couplings that were inserted by hand and saturates.")
P()

# ------------------------------------------------------------------------- contact-or-nothing
P("-" * 104)
P("  THE CONTACT STATEMENT, EXACTLY.   Energy attributable to the pair (i, j) at separation d,")
P("  defined as the coefficient of Z_i Z_j in H.  No fitting, no ensemble, just read off H.")
P("-" * 104)
P(f"  {'family':<10} {'n':>4} " + " ".join(f"{'d='+str(d):>10}" for d in (1, 2, 3, 4, 8, 16))
  + f"   {'CONTROL: sum over all d':>24}")
for fam in ("uniform", "randpos"):
    n = 33
    J = couplings(fam, n - 1)
    cells = []
    for d in (1, 2, 3, 4, 8, 16):
        coeff = J[0] if d == 1 else 0
        cells.append(f"{coeff:>10}")
    P(f"  {fam:<10} {n:>4} " + " ".join(cells) + f"   {sum(abs(j) for j in J):>24}")
P()
P("  READ: the two-body energy is exactly zero at every separation beyond 1, while the same table's")
P("  control column shows the total is large -- the instrument is reading a real H.  The energy of")
P("  this construction is CONTACT-OR-NOTHING, which is exactly the verdict C-47 reached for")
P("  spread(Phi).  Clause (e) is not merely failed by an exponential here; there is no falloff to")
P("  measure, because there is nothing at range at all.")
P()

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_C_SIGN/s4_separation.txt", "w").write("\n".join(OUT) + "\n")
