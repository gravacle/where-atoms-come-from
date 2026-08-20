"""
O-50-C  STEP 3.   PART 1.  IS MEAN-ZERO FORCED, AND BY WHAT?

The theorem candidate as written says:

    IF records are independently writable, G_W acts simply transitively on record configurations;
    THEN every G_W-invariant functional is CONSTANT,
    AND  every non-constant functional is G_W-ODD, hence has mean exactly zero, hence CANCELS.

Step 2 proved the first line and the second line on the torus, exactly.  This step tests the
THIRD line by COMPLETE ENUMERATION -- every functional in an exactly enumerable class, not a
sample -- and finds the third line FALSE.  It then states, and verifies, the exact condition that
replaces it, and separates the role of the GROUP from the role of the MEASURE.

EVERYTHING HERE IS EXACT.  Values are integers or fractions.Fraction; no floating point is used
in any classification or count (D-19).

CARRIER.  Step 2's torus carrier: k disjoint L x L tori, m = 2k records, G_W = (Z_2)^m acting
simply transitively by translation.  Because the action is FIXED by step 2, this step needs only
the action, and the counts below are therefore counts on the torus.
"""
import sys, os, itertools, random
from fractions import Fraction

OUT = []
def P(s=""):
    OUT.append(str(s)); print(s)
BAR = "=" * 104
bar = "-" * 104

def configs(m):
    """all 2^m record configurations as +-1 tuples; index x has s_i = +1 iff bit i of x is 0"""
    return [tuple(1 - 2 * ((x >> i) & 1) for i in range(m)) for x in range(1 << m)]

def wht(vals):
    """EXACT Walsh-Hadamard transform.  Returns fhat[S] = 2^-m sum_s f(s) chi_S(s), as Fractions.
       fhat[0] is the MEAN of f over the uniform measure on configurations."""
    n = len(vals); a = [Fraction(v) for v in vals]
    h = 1
    while h < n:
        for i in range(0, n, h * 2):
            for j in range(i, i + h):
                x, y = a[j], a[j + h]
                a[j], a[j + h] = x + y, x - y
        h *= 2
    return [c / n for c in a]

P(BAR)
P("O-50-C  STEP 3.   PART 1: IS MEAN-ZERO FORCED, AND BY WHAT?   COMPLETE ENUMERATION, EXACT.")
P(BAR)
P()
P("D-23 SCOPE: TORUS (the action is step 2's, established on k disjoint L x L tori).")
P()

# ================================================================ A. complete enumeration
P(bar)
P("  A.  EVERY +-1-VALUED FUNCTIONAL, COUNTED.  NOT A SAMPLE -- THE WHOLE CLASS.")
P(bar)
P("  A functional f : {+-1}^m -> {+1,-1} is a vector of 2^m signs, so there are 2^(2^m) of them:")
P("  16 at m=2, 256 at m=3, 65536 at m=4.  Every one is enumerated and classified on three axes:")
P("     INVARIANT      f(g.s) = f(s) for all g in G_W")
P("     MEAN ZERO      sum over all configurations of f(s) is exactly 0")
P("     WRITER-ODD     f(g.s) = -f(s) for SOME g in G_W  (the theorem candidate's hypothesis)")
P()
P(f"  {'m':>3} {'#functionals':>13} {'invariant':>10} {'non-invariant':>14} {'writer-ODD':>11}"
  f" {'mean = 0':>9} {'NON-INVARIANT AND MEAN != 0':>28} {'fraction':>9}")
tableA = {}
for m in (2, 3, 4):
    N = 1 << m
    gens = list(range(m))
    allg = list(range(N))                       # every group element as an xor-mask
    inv = odd = mz = bad = oddbad = 0
    total = 1 << N
    for code in range(total):
        f = [1 if (code >> x) & 1 else -1 for x in range(N)]
        is_inv = all(f[x ^ g] == f[x] for g in allg for x in range(N))
        is_mz = (sum(f) == 0)
        is_odd = any(all(f[x ^ g] == -f[x] for x in range(N)) for g in allg if g)
        inv += is_inv; mz += is_mz; odd += is_odd
        if (not is_inv) and (not is_mz): bad += 1
        if is_odd and not is_mz: oddbad += 1
    assert oddbad == 0, "a writer-odd functional with non-zero mean would refute the one true implication"
    tableA[m] = (total, inv, total - inv, odd, mz, bad)
    P(f"  {m:>3} {total:>13} {inv:>10} {total - inv:>14} {odd:>11} {mz:>9} {bad:>28}"
      f" {Fraction(bad, total)!s:>9}")
P()
tot4, inv4, non4, odd4, mz4, bad4 = tableA[4]
P("  READ, and it is the decisive number of Part 1:")
P(f"  at m = 4 there are {non4} NON-INVARIANT functionals; only {mz4} of them have mean zero.")
P(f"  {bad4} of them -- {float(Fraction(bad4, non4)) * 100:.1f}% -- are NON-INVARIANT AND HAVE"
  f" NON-ZERO MEAN.")
P()
P("  ==> THE THEOREM CANDIDATE'S THIRD LINE IS FALSE.  NON-INVARIANCE DOES NOT IMPLY WRITER-ODD,")
P("      AND IT DOES NOT IMPLY MEAN ZERO.  This is a REFUTATION BY COMPLETE ENUMERATION, not by a")
P("      failed search: the class was exhausted.")
P()
P(f"  And WRITER-ODD is a strictly stronger condition than mean-zero: at m = 4, {odd4} functionals")
P(f"  are writer-odd, {mz4} have mean zero.  Every writer-odd one has mean zero (the implication")
P(f"  that IS true), but the converse fails on {mz4 - odd4} functionals.")
P()
P("  THE SMALLEST EXPLICIT COUNTEREXAMPLES, exhibited rather than counted:")

def show(m, name, fn):
    cs = configs(m)
    vals = [fn(s) for s in cs]
    N = 1 << m
    is_inv = all(vals[x ^ g] == vals[x] for g in range(N) for x in range(N))
    is_odd = any(all(vals[x ^ g] == -vals[x] for x in range(N)) for g in range(1, N))
    fh = wht(vals)
    P(f"    m={m}  {name:<26} values {str(vals):<44} invariant {str(is_inv):<6} writer-odd "
      f"{str(is_odd):<6} EXACT mean {str(fh[0]):<8}")

show(2, "f = (1 + s1)/2", lambda s: Fraction(1 + s[0], 2))
show(2, "f = s1", lambda s: s[0])
show(2, "f = 1", lambda s: 1)
show(2, "f = s1 + s2 + s1s2", lambda s: s[0] + s[1] + s[0] * s[1])
show(3, "f = s1 + s2 + s1s2", lambda s: s[0] + s[1] + s[0] * s[1])
P()
P("  f = (1+s1)/2 -- 'record 1 reads +1' -- is NON-INVARIANT (the writer of record 1 changes it),")
P("  is NOT writer-odd, and has EXACT mean 1/2.  It does NOT cancel.  One line kills the candidate.")
P("  f = s1 + s2 + s1s2 has EXACT mean 0 and is odd under NO writer: mean-zero does not require")
P("  writer-oddness either.  The two conditions are genuinely distinct, in both directions.")

# ================================================================ B. the exact condition
P()
P(bar)
P("  B.  THE EXACT CONDITION THAT REPLACES IT, AND ITS PROOF")
P(bar)
P("  Because G_W = (Z_2)^m acts simply transitively, the space of functionals is the REGULAR")
P("  representation.  It decomposes into the 2^m one-dimensional characters chi_S, and every")
P("  functional has a unique exact expansion  f = sum_S fhat(S) chi_S.  Then, identically,")
P()
P("        mean of f over the uniform measure  =  fhat(EMPTY SET) .")
P()
P("  EXACT CONDITION.   mean(f) = 0  <=>  fhat(EMPTY) = 0  <=>  f is ORTHOGONAL TO THE CONSTANTS")
P("                                    <=>  f lies in the sum of the NON-TRIVIAL isotypic components.")
P()
P("  Proof: chi_EMPTY = 1, and the characters are orthogonal, so fhat(EMPTY) = <f,1> = mean(f).")
P("  The group action enters ONLY through the fact that chi_S for S non-empty has mean 0, which")
P("  holds because such a chi_S is negated by any writer that flips an odd number of the records")
P("  in S -- and simple transitivity guarantees such a writer EXISTS.")
P()
P("  RELATIONS AMONG THE THREE CONDITIONS, verified by the complete enumeration above:")
P("        writer-odd  ==>  mean zero  ==>  non-invariant (for non-constant f)")
P(f"        mean zero   =/=>  writer-odd            ({mz4 - odd4} counterexamples at m=4)")
P(f"        non-invariant =/=> mean zero            ({bad4} counterexamples at m=4)")
P()
P("  VERIFICATION that fhat(EMPTY) is the mean, on exact random rational functionals:")
rnd = random.Random(11)
P(f"    {'m':>3} {'trial':>6} {'sum_s f(s) / 2^m (exact)':>26} {'fhat(EMPTY) (exact)':>22} {'equal?':>7}")
for m in (2, 3, 4, 5):
    for t in range(2):
        vals = [Fraction(rnd.randrange(-9, 10), rnd.randrange(1, 5)) for _ in range(1 << m)]
        direct = sum(vals) / (1 << m)
        fh = wht(vals)
        P(f"    {m:>3} {t:>6} {str(direct):>26} {str(fh[0]):>22} {str(direct == fh[0]):>7}")
P()
P("  Exact rational agreement in every case.  The transform is the mean; no floating point involved.")

# ============================================== C. group or measure?
P()
P(bar)
P("  C.  IS MEAN-ZERO A CONSEQUENCE OF THE GROUP, OR AN ARTIFACT OF THE UNIFORM MEASURE?")
P(bar)
P("  BOTH ARE NEEDED, and the exact statement separating them is this:")
P()
P("    (1) On a G-TORSOR -- a set on which G acts simply transitively -- the UNIFORM measure is the")
P("        UNIQUE G-invariant probability measure.  (If mu(g.s) = mu(s) for all g and the action is")
P("        transitive, mu is constant; normalisation fixes the constant.)")
P("    (2) Under ANY G-invariant measure, the average of f equals the average of its INVARIANT part.")
P("        On a torsor the invariant part is the constant fhat(EMPTY), so the mean is fhat(EMPTY).")
P()
P("  So mean-zero for the non-trivial isotypic components is an EXACT CONSEQUENCE OF THE GROUP")
P("  ACTION, PROVIDED THE MEASURE IS REQUIRED TO BE WRITER-INVARIANT -- and on the torus that")
P("  requirement pins the measure down to exactly one, the uniform one.  Drop writer-invariance of")
P("  the measure and NOTHING survives: the mean of chi_S becomes whatever the measure says.")
P()
P("  Verified exactly.  Uniqueness of the invariant measure, by exact linear algebra on the actual")
P("  action (D-15: the same solver on the two non-transitive controls, which have MANY invariant")
P("  measures):")
P()

def invariant_measure_dim(m, masks):
    """dimension of the affine space of G-invariant probability measures = (#orbits - 1)."""
    N = 1 << m
    parent = list(range(N))
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for x in range(N):
        for g in masks:
            ra, rb = find(x), find(x ^ g)
            if ra != rb: parent[ra] = rb
    return len({find(x) for x in range(N)})

P(f"  {'carrier / writer group':<44} {'m':>3} {'#orbits':>8} {'#invariant probability measures':>32}"
  f" {'unique?':>8}")
for m in (2, 4, 6):
    full = [1 << j for j in range(m)]
    o = invariant_measure_dim(m, full)
    P(f"  {'TORUS  G_W = (Z_2)^m, simply transitive':<44} {m:>3} {o:>8} "
      f"{('a ' + str(o - 1) + '-parameter family' if o > 1 else 'EXACTLY ONE (uniform)'):>32}"
      f" {str(o == 1):>8}")
    o1 = invariant_measure_dim(m, [(1 << m) - 1])
    P(f"  {'CTRL-1 chain shape: one global flip':<44} {m:>3} {o1:>8} "
      f"{('a ' + str(o1 - 1) + '-parameter family'):>32} {str(o1 == 1):>8}")
    o2 = invariant_measure_dim(m, [1 << j for j in range(m - 1)])
    P(f"  {'CTRL-2 writers reach only m-1 records':<44} {m:>3} {o2:>8} "
      f"{('a ' + str(o2 - 1) + '-parameter family'):>32} {str(o2 == 1):>8}")
P()
P("  READ: on the torus there is EXACTLY ONE writer-invariant measure at every m, so 'uniform' is")
P("  not a modelling choice -- it is the unique measure compatible with the writer symmetry.  On")
P("  CTRL-1 there are 2^{m-1}-1 free parameters and on CTRL-2 there is 1: the solver returns a")
P("  family when a family is there, so the uniqueness on the torus is a measured result.")
P()
P("  AND THE OTHER SIDE, WHICH IS THE ONE THAT MATTERS FOR PART 2.  Exact means of the SAME")
P("  functional under measures that are NOT writer-invariant:")
P()
m = 4
cs = configs(m)
f_sum = lambda s: sum(s)
P(f"  {'measure on record configurations':<46} {'writer-invariant?':>18} {'EXACT mean of s1':>18}"
  f" {'EXACT mean of (s1+..+s4)':>25}")
def rep(name, wts, invflag):
    tot = sum(wts)
    m1 = sum(w * s[0] for w, s in zip(wts, cs)) / Fraction(tot)
    ms = sum(w * sum(s) for w, s in zip(wts, cs)) / Fraction(tot)
    P(f"  {name:<46} {str(invflag):>18} {str(m1):>18} {str(ms):>25}")
rep("uniform (the unique invariant one)", [Fraction(1)] * len(cs), True)
for p in (Fraction(1, 2), Fraction(3, 5), Fraction(3, 4), Fraction(1)):
    wts = []
    for s in cs:
        w = Fraction(1)
        for si in s: w *= (p if si == 1 else 1 - p)
        wts.append(w)
    if p == 1:
        wts = [Fraction(1) if all(si == 1 for si in s) else Fraction(0) for s in cs]
    rep(f"product, P(record = +1) = {p}", wts, p == Fraction(1, 2))
rep("delta at the single configuration (+,+,+,+)",
    [Fraction(1) if all(si == 1 for si in s) else Fraction(0) for s in cs], False)
rep("delta at the single configuration (+,-,+,-)",
    [Fraction(1) if s == (1, -1, 1, -1) else Fraction(0) for s in cs], False)
P()
P("  READ: the moment the measure stops being writer-invariant, the mean of s1 and of the sum are")
P("  NON-ZERO and grow to their maximum.  The cancellation law is a statement ABOUT THE MEASURE as")
P("  much as about the group.  Nothing in the five clauses names a measure.")
P()
P("  THE LAST LINE OF THIS TABLE IS THE PIVOT OF THE WHOLE PROBE.  The universe is in ONE")
P("  configuration -- a delta measure -- which is never writer-invariant unless m = 0.  The next")
P("  step asks what the functional is WORTH at one configuration, which is a different question")
P("  from what its mean is.")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "s3_means.txt"), "w") as f:
    f.write("\n".join(OUT) + "\n")
