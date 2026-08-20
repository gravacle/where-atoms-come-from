"""
VERIFY 1.  Independent recomputation of every combinatorial number the lane reports in
Parts 1 and 2, WITHOUT importing any lane code and WITHOUT any carrier at all.

THE POINT OF DOING IT CARRIER-FREE: if the identical numbers fall out of the bare hypercube
{+-1}^m with (Z_2)^m acting by translation, then NOTHING in Parts 1 and 2 is a torus result.
"""
from fractions import Fraction
from itertools import combinations
import math

def line(s=""): print(s, flush=True)

line("=" * 100)
line("V1.  PARTS 1 AND 2 RECOMPUTED FROM THE BARE HYPERCUBE.  NO CARRIER, NO TORIC CODE, NO LANE CODE.")
line("=" * 100)

# ---------------------------------------------------------------- A. complete functional enumeration
line()
line("A.  COMPLETE ENUMERATION OF ALL +-1-VALUED FUNCTIONALS  (lane step 3A)")
line(f"  {'m':>3} {'#f':>8} {'invariant':>10} {'non-inv':>8} {'writer-ODD':>11} {'mean=0':>8} {'non-inv&mean!=0':>16} {'lane says':>28}")
lane = {2: (16, 2, 14, 6, 6, 8), 3: (256, 2, 254, 70, 70, 184), 4: (65536, 2, 65534, 2790, 12870, 52664)}
for m in (2, 3, 4):
    N = 1 << m
    tot = 1 << N
    inv = nodd = mz = 0
    ninv_mnz = 0
    gens = [1 << j for j in range(m)]          # translations by e_j
    allg = list(range(1, N))                    # every non-identity group element
    for code in range(tot):
        vals = [1 if (code >> x) & 1 else -1 for x in range(N)]
        is_inv = all(vals[x ^ g] == vals[x] for g in gens for x in range(N))
        is_odd = any(all(vals[x ^ g] == -vals[x] for x in range(N)) for g in allg)
        s = sum(vals)
        if is_inv: inv += 1
        if is_odd: nodd += 1
        if s == 0: mz += 1
        if (not is_inv) and s != 0: ninv_mnz += 1
    got = (tot, inv, tot - inv, nodd, mz, ninv_mnz)
    line(f"  {m:>3} {tot:>8} {inv:>10} {tot-inv:>8} {nodd:>11} {mz:>8} {ninv_mnz:>16} "
         f"{str(lane[m]):>28}  MATCH={got == lane[m]}")
line("  READ: reproduced exactly, with no carrier in the calculation.  Also note the identity")
line("  #(non-invariant AND mean != 0) = 2^(2^m) - 2 - C(2^m, 2^(m-1)) -- pure counting:")
for m in (2, 3, 4):
    N = 1 << m
    line(f"    m={m}:  2^{N} - 2 - C({N},{N//2}) = {(1<<N) - 2 - math.comb(N, N//2)}")

# ---------------------------------------------------------------- B. is writer-odd carrier-dependent?
line()
line("B.  IS 'writer-odd' A TORUS FACT?  Count it for the CHAIN-SHAPED group (one global flip) too.")
line(f"  {'m':>3} {'group':<28} {'writer-ODD':>11} {'mean=0':>8} {'invariant':>10}")
for m in (2, 3, 4):
    N = 1 << m
    for gname, G in (("(Z_2)^m  (torus)", list(range(1, N))),
                     ("Z_2 global flip (chain)", [N - 1])):
        gens = [1 << j for j in range(m)] if "torus" in gname else [N - 1]
        inv = nodd = mz = 0
        for code in range(1 << N):
            vals = [1 if (code >> x) & 1 else -1 for x in range(N)]
            if all(vals[x ^ g] == vals[x] for g in gens for x in range(N)): inv += 1
            if any(all(vals[x ^ g] == -vals[x] for x in range(N)) for g in G): nodd += 1
            if sum(vals) == 0: mz += 1
        line(f"  {m:>3} {gname:<28} {nodd:>11} {mz:>8} {inv:>10}")
line("  READ: 'mean=0' count is C(2^m,2^(m-1)) regardless of the group -- it is not a group fact at")
line("  all.  The lane's decisive 52664 is 2^16 - 2 - C(16,8) and would be identical on the chain.")

# ---------------------------------------------------------------- C. the coherence law
line()
line("C.  THE COHERENCE LAW, EXACT RATIONALS  (lane steps 4A/4B, 5)")
def EabsF_over_m(m):
    """E|s_1+...+s_m| / m, exact."""
    tot = 0
    for q in range(m + 1):
        tot += math.comb(m, q) * abs(m - 2 * q)
    return Fraction(tot, m * (1 << m))
line(f"  {'m':>6} {'E|F|/m exact':>28} {'C(m,m/2)/2^m':>28} {'equal?':>7} {'lane value':>22}")
lane_coh = {16: "0.19638061523437500000", 256: "0.04981910993614015124",
            4096: "0.01246618536376025958", 65536: "0.00311672467625241587"}
for m in (2, 4, 16, 256, 4096, 65536):
    a = EabsF_over_m(m)
    b = Fraction(math.comb(m, m // 2), 1 << m)
    la = lane_coh.get(m, "")
    line(f"  {m:>6} {float(a):>28.20f} {float(b):>28.20f} {str(a==b):>7} {la:>22}")
line("  READ: E|F|/m is EXACTLY the central binomial C(m,m/2)/2^m at every m, whose Stirling")
line("  asymptotic is sqrt(2/(pi m)).  The 'law' sqrt(2/pi)*m^(-1/2) is the textbook asymptotic of")
line("  the central binomial coefficient; the 1/(4m) correction the lane calls 'predicted' is the")
line("  standard next Stirling term.  Reproduced here without any carrier.")
line(f"  E|F| at m=65536 = {float(EabsF_over_m(65536)*65536):.4f}   (lane: 204.2577)   "
     f"sqrt(2*65536/pi) = {math.sqrt(2*65536/math.pi):.4f}")

# ---------------------------------------------------------------- D. Ising variance
line()
line("D.  I-3 ISING: Var of sum s on a 1D RING with nn correlation t, exact, and the limit.")
def var_ring_exact(m, t):
    """<s_i s_j> = (t^d + t^(m-d))/(1+t^m) on a ring of m sites, d = |i-j| cyclic."""
    tot = Fraction(0)
    tm = t ** m
    for d in range(m):
        tot += m * (t ** d + t ** (m - d)) / (1 + tm)
    return tot
line(f"  {'m':>5} {'t':>10} {'Var exact (ring)':>24} {'lane Var':>16} {'Var/m':>10} {'(1+t)/(1-t)':>12}")
lane_var = {(16, Fraction(1,3)): "30.50000", (16, Fraction(3,5)): "56.50212",
            (16, Fraction(9,11)): "112.49619", (16, Fraction(99,101)): "230.84338"}
for t in (Fraction(1,3), Fraction(3,5), Fraction(9,11), Fraction(99,101)):
    v = var_ring_exact(16, t)
    line(f"  {16:>5} {float(t):>10.6f} {float(v):>24.5f} {lane_var[(16,t)]:>16} "
         f"{float(v)/16:>10.4f} {float((1+t)/(1-t)):>12.4f}")
line("  (large-m rows use the closed form below)")
line(f"  {'m':>8} {'t':>10} {'coh*sqrt(m) = sqrt(2/pi)*sqrt(Var/m)':>38} {'lane':>12}")
lane_i3 = {(64, Fraction(1,2)): 1.367505, (256, Fraction(1,2)): 1.378373,
           (4096, Fraction(1,2)): 1.381752, (262144, Fraction(1,2)): 1.381973,
           (262144, Fraction(999,1001)): 25.207251, (4096, Fraction(999,1001)): 23.641679}
for (m, t), lv in sorted(lane_i3.items()):
    # closed form for the ring: Var/m = sum_d (t^d + t^(m-d))/(1+t^m)
    tf = float(t)
    tm = tf ** m if m * math.log(abs(tf) if tf else 1e-300) > -700 else 0.0
    s = 0.0
    for d in range(m):
        s += (tf ** d + (tf ** (m - d) if (m - d) * (-math.log(tf)) < 700 else 0.0))
    s /= (1 + tm)
    val = math.sqrt(2 / math.pi) * math.sqrt(s)
    line(f"  {m:>8} {float(t):>10.6f} {val:>38.6f} {lv:>12.6f}")
line("  READ: Var/m -> (1+t)/(1-t) is the textbook 1D-Ising susceptibility sum "
     "sum_d t^|d| = (1+t)/(1-t).")

# ---------------------------------------------------------------- E. biased measure crossover
line()
line("E.  I-2 BIASED MEASURE.  coherence = E|sum s|/m with iid P(+1)=(1+lam)/2, exact.")
def coh_biased(m, lam):
    p = (1 + lam) / 2
    tot = Fraction(0)
    num = Fraction(0)
    for q in range(m + 1):
        w = Fraction(math.comb(m, q)) * (Fraction(p) ** (m - q)) * ((1 - Fraction(p)) ** q)
        num += w * abs(m - 2 * q)
    return num / m
line(f"  {'m':>5} {'lam':>8} {'coherence exact':>20} {'lane':>12} {'sqrt(lam^2+2/(pi m))':>22}")
lane_i2 = {(16, Fraction(0)): 0.1963806, (16, Fraction(1,64)): 0.1967641, (16, Fraction(1,16)): 0.2024897,
           (16, Fraction(1,4)): 0.2879194, (16, Fraction(1,2)): 0.5023604,
           (256, Fraction(0)): 0.0498191, (256, Fraction(1,16)): 0.0728554, (256, Fraction(1,4)): 0.2500006}
for (m, lam), lv in sorted(lane_i2.items()):
    c = coh_biased(m, lam)
    line(f"  {m:>5} {float(lam):>8.5f} {float(c):>20.7f} {lv:>12.7f} "
         f"{math.sqrt(float(lam)**2 + 2/(math.pi*m)):>22.7f}")
line("  READ: reproduced.  The crossover at lam ~ m^(-1/2) is the ordinary law-of-large-numbers")
line("  crossover between a bias lam*m and a fluctuation sqrt(m).")

# ---------------------------------------------------------------- F. the 'both signs' claim
line()
line("F.  THE [C5](b) 'RESPONSE TAKES BOTH SIGNS' CLAIM: how much of it is about records?")
line("  Claim: for any bijection g of a finite set and any f, sum_s [f(g.s)-f(s)] = 0.")
import random
random.seed(7)
ok = True
for trial in range(200):
    N = random.choice([5, 6, 7, 11])            # NOT a power of 2, NOT a group, NOT a carrier
    perm = list(range(N)); random.shuffle(perm)
    f = [random.randint(-50, 50) for _ in range(N)]
    if sum(f[perm[s]] - f[s] for s in range(N)) != 0: ok = False
line(f"  200 random permutations of arbitrary finite sets (sizes 5,6,7,11): sum of deltas always 0? {ok}")
line("  READ: the lane's own caveat 6 says this; it is true of ANY permutation of ANY finite set and")
line("  has no content about records, the five clauses, or the toric code.")
