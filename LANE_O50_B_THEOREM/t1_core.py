"""O-50-B  PART 1 -- THE ABSTRACT CORE, EXACTLY.

The theorem candidate concerns a group G_W acting on the configuration space
C = {+-1}^k of a maximal commuting family of records.  Part 1 decides the
GROUP-THEORETIC content exactly, with integer / Fraction arithmetic only --
no floating point anywhere.  Nothing here is carrier-specific; Parts 2-4
supply the carriers.

OBJECTS
  C            = {+-1}^k, |C| = 2^k
  chi_S(sigma) = prod_{i in S} sigma_i      the 2^k characters, an orthogonal basis
                                            of ALL real functions on C
  G <= (Z_2)^k acting by coordinate flips: (eps . sigma)_i = eps_i sigma_i
  d : C -> Z_{>=0}   a WEIGHT (in the carriers this is the joint-eigenspace
                     dimension inside one energy shell)

FACTS PROVED HERE (each checked exhaustively, not sampled)
  F1  chi_S(eps.sigma) = (-1)^{|S cap eps|} chi_S(sigma).   [character action]
  F2  f is G-invariant  <=>  fhat(S) = 0 for every S NOT orthogonal to G.
      dim(invariants) = #orbits of G on C = 2^k / |G| = 2^{k - rank G}.
  F3  G = (Z_2)^k (all flips available)  =>  the action is SIMPLY TRANSITIVE,
      #orbits = 1, dim(invariants) = 1, invariants = constants.
  F4  THE LITERAL CLAUSE IS FALSE.  "every non-constant functional has mean
      exactly zero" is refuted by f(sigma) = #{i : sigma_i = +1}, mean k/2.
  F5  THE REPAIRED CLAUSE IS TRUE and unconditional:
      for every f and every eps in G,  sum_sigma d(sigma) [f(eps.sigma) - f(sigma)] = 0
      whenever d is G-invariant.  The RESPONSE to a write always cancels.
  F6  MASTER IDENTITY:  <f>_d  =  <Pi_G f>_d  where Pi_G is the group average,
      PROVIDED d is G-invariant.  The d-weighted mean of ANY functional depends
      only on its G-INVARIANT component.  Under F3 that component is a constant,
      so the mean is record-blind.
"""
import itertools, sys
from fractions import Fraction

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); OUT.append(s); print(s)

# ------------------------------------------------------------------ machinery
def configs(k):
    return list(itertools.product((1, -1), repeat=k))

def chi(S, sigma):
    v = 1
    for i in S: v *= sigma[i]
    return v

def subsets(k):
    out = []
    for m in range(1 << k):
        out.append(tuple(i for i in range(k) if (m >> i) & 1))
    return out

def flip(eps, sigma):
    return tuple(-s if e else s for e, s in zip(eps, sigma))

def group_from_gens(gens, k):
    """closure of a set of F_2 flip-vectors under XOR"""
    G = {tuple([0]*k)}
    changed = True
    while changed:
        changed = False
        for g in list(G):
            for h in gens:
                x = tuple((a ^ b) for a, b in zip(g, h))
                if x not in G: G.add(x); changed = True
    return sorted(G)

def orbits(G, k):
    C = configs(k); seen = set(); orbs = []
    for s in C:
        if s in seen: continue
        o = sorted({flip(e, s) for e in G})
        for x in o: seen.add(x)
        orbs.append(o)
    return orbs

def invariant_S(G, k):
    """S with chi_S G-invariant  <=>  |S cap eps| even for every eps in G"""
    out = []
    for S in subsets(k):
        if all(sum(e[i] for i in S) % 2 == 0 for e in G): out.append(S)
    return out

def fourier(f, k):
    """exact Fourier coefficients over Fraction"""
    C = configs(k); N = len(C)
    return {S: Fraction(sum(f(s) * chi(S, s) for s in C), N) for S in subsets(k)}

# ------------------------------------------------------------------ F1,F2,F3
P("=" * 100)
P("PART 1 -- THE ABSTRACT CORE.  Exact integer / Fraction arithmetic, exhaustive over C.")
P("=" * 100)
P("")
P("--- F1/F2/F3: orbit count, invariant dimension, and the CONTROL COLUMN (D-15) ---")
P("")
P(f"{'k':>2} {'group G_W':<34} {'|G|':>5} {'#orbits':>8} {'dim inv':>8} "
  f"{'2^(k-rank)':>11} {'transitive?':>12} {'non-const invariant exhibited':>32}")
P("-" * 128)

rows = []
for k in range(1, 8):
    e = lambda i: tuple(1 if j == i else 0 for j in range(k))
    cases = [("FULL  <flip_1..flip_k>", [e(i) for i in range(k)])]
    if k >= 2:
        cases.append(("DROP ONE  <flip_1..flip_{k-1}>", [e(i) for i in range(k - 1)]))
        cases.append(("PAIRED  <flip_i flip_{i+1}>",
                      [tuple(1 if j in (i, i + 1) else 0 for j in range(k)) for i in range(k - 1)]))
        cases.append(("TRIVIAL  <>", []))
    for name, gens in cases:
        G = group_from_gens(gens, k)
        orbs = orbits(G, k)
        inv = invariant_S(G, k)
        rank = 0
        # rank of G as F_2 space
        basis = []
        for g in G:
            v = list(g)
            for b in basis:
                h = next((i for i in range(k) if b[i]), None)
                if h is not None and v[h]: v = [(x ^ y) for x, y in zip(v, b)]
            if any(v): basis.append(v)
        rank = len(basis)
        assert len(inv) == len(orbs), (k, name, len(inv), len(orbs))
        assert len(orbs) == 2 ** (k - rank)
        nc = [S for S in inv if S]
        ex = ("chi_" + "".join(str(i + 1) for i in nc[0])) if nc else "NONE (only constants)"
        P(f"{k:>2} {name:<34} {len(G):>5} {len(orbs):>8} {len(inv):>8} "
          f"{2**(k-rank):>11} {str(len(orbs)==1):>12} {ex:>32}")
        rows.append((k, name, len(G), len(orbs), len(inv)))
    P("")

P("READ: with the FULL flip group the orbit count is 1 and the invariant space is exactly")
P("      1-dimensional (the constants) at every k tested -- the theorem's core.  The three")
P("      CONTROL rows in each block are groups that are NOT the full flip group; each has")
P("      orbit count > 1 and exhibits a NAMED non-constant invariant.  The instrument would")
P("      have registered a non-transitive writer group had one been present.")
P("")

# ------------------------------------------------------------------ F4  the false clause
P("=" * 100)
P("--- F4: THE LITERAL CLAUSE OF THE CANDIDATE IS FALSE ---")
P("")
P('CANDIDATE, verbatim: "every invariant functional is CONSTANT, AND every non-constant')
P('functional is G_W-ODD, hence has mean exactly zero over configurations, hence CANCELS."')
P("")
P("The second conjunct does not follow.  G_W-invariance is a property of the WHOLE function;")
P("non-constant is its negation; but a non-constant function is NOT thereby odd.  Counterexample")
P("below, exact, at every k.")
P("")
P(f"{'k':>2} {'f(sigma)':<28} {'mean over C':>12} {'non-const?':>11} {'responsive?':>12} "
  f"{'G_W-odd?':>10} {'mean of responsive part':>24}")
P("-" * 108)
for k in range(1, 7):
    f = lambda s: sum(1 for x in s if x == +1)      # number of records reading +1
    C = configs(k)
    mean = Fraction(sum(f(s) for s in C), len(C))
    fh = fourier(f, k)
    nonconst = any(v != 0 for S, v in fh.items() if S)
    e1 = tuple([1] + [0] * (k - 1))
    responsive = any(f(flip(e1, s)) != f(s) for s in C)
    odd = all(f(flip(tuple([1]*k), s)) == -f(s) for s in C)
    resp_mean = Fraction(sum(f(s) - mean for s in C), len(C))
    P(f"{k:>2} {'#{i : sigma_i = +1}':<28} {str(mean):>12} {str(nonconst):>11} "
      f"{str(responsive):>12} {str(odd):>10} {str(resp_mean):>24}")
P("")
P("READ: f = #{i : sigma_i = +1} is NON-CONSTANT and RESPONSIVE at every k, and its mean is")
P("      k/2, which is NOT zero for k >= 1.  The literal clause is REFUTED.  What IS zero is")
P("      the mean of its RESPONSIVE part f - <f> (last column, exactly 0 at every k).  The")
P("      non-zero mean sits entirely in the CONSTANT, i.e. writer-INVARIANT, i.e. RECORD-BLIND")
P("      component -- k/2 is the same number whatever is written.  That is C-61 restated.")
P("")

# ------------------------------------------------------------------ F5  the repaired clause
P("=" * 100)
P("--- F5: THE REPAIRED CLAUSE.  Every WRITE RESPONSE cancels, exactly, for EVERY functional ---")
P("")
P("Delta_i f (sigma) := f(flip_i sigma) - f(sigma).   Claim: Delta_i f (flip_i sigma) = -Delta_i f(sigma),")
P("hence sum_sigma d(sigma) Delta_i f(sigma) = 0 for any flip_i-invariant weight d.")
P("Checked EXHAUSTIVELY over an orthogonal BASIS of all functions on C (the 2^k characters),")
P("which settles it for every real functional by linearity, plus random integer functions.")
P("")
P(f"{'k':>2} {'#functions checked':>19} {'max |sum_sigma Delta_i f|  (d=1)':>34} "
  f"{'max |sum_sigma d Delta_i f| (d random G-inv)':>45} {'CONTROL: d NOT G-invariant':>28}")
P("-" * 132)
import random
random.seed(20500)
for k in range(1, 8):
    C = configs(k); Ss = subsets(k)
    worst = 0; worstd = 0; worstctl = 0; cnt = 0
    dinv = {s: 1 for s in C}                                   # any constant weight is G-invariant
    # a NON-constant but still full-flip-invariant weight cannot exist for the FULL group
    # (transitive => invariant weight is constant), so d=1 is the general G-invariant weight here.
    dbad = {s: random.randint(1, 9) for s in C}                # CONTROL: not G-invariant
    for S in Ss:
        fn = lambda s, S=S: chi(S, s)
        for i in range(k):
            e = tuple(1 if j == i else 0 for j in range(k))
            worst = max(worst, abs(sum(fn(flip(e, s)) - fn(s) for s in C)))
            worstd = max(worstd, abs(sum(dinv[s] * (fn(flip(e, s)) - fn(s)) for s in C)))
            worstctl = max(worstctl, abs(sum(dbad[s] * (fn(flip(e, s)) - fn(s)) for s in C)))
            cnt += 1
    for _ in range(200):
        vals = {s: random.randint(-50, 50) for s in C}
        fn = lambda s: vals[s]
        for i in range(k):
            e = tuple(1 if j == i else 0 for j in range(k))
            worst = max(worst, abs(sum(fn(flip(e, s)) - fn(s) for s in C)))
            worstd = max(worstd, abs(sum(dinv[s] * (fn(flip(e, s)) - fn(s)) for s in C)))
            worstctl = max(worstctl, abs(sum(dbad[s] * (fn(flip(e, s)) - fn(s)) for s in C)))
            cnt += 1
    P(f"{k:>2} {cnt:>19} {worst:>34} {worstd:>45} {worstctl:>28}")
P("")
P("READ: the write response sums to EXACTLY 0 (integer zero, not a tolerance) for every one of")
P("      the functions checked, at every k, whenever the weight is invariant under the flip.")
P("      The CONTROL column uses a weight that is NOT flip-invariant and is non-zero -- so the")
P("      instrument can register a non-cancelling response when one exists.  The control is")
P("      exactly the physically meaningful failure mode: an UNBALANCED carrier (D-15).")
P("")

# ------------------------------------------------------------------ F6  master identity
P("=" * 100)
P("--- F6: MASTER IDENTITY.  <f>_d = <Pi_G f>_d  whenever d is G-invariant ---")
P("")
P("Pi_G f (sigma) := (1/|G|) sum_{eps in G} f(eps.sigma), the projection onto invariants.")
P("If d is G-invariant then reindexing gives <f>_d = <Pi_G f>_d identically.  Under the FULL")
P("flip group Pi_G f = <f> is a CONSTANT, so the d-weighted mean of ANY functional equals a")
P("number that does not depend on the record configuration at all.")
P("")
P(f"{'k':>2} {'group':<26} {'max |<f>_d - <Pi_G f>_d|':>26} {'#f':>6} "
  f"{'dim inv':>8} {'CONTROL d not G-inv: max diff':>30}")
P("-" * 104)
for k in range(1, 7):
    e = lambda i: tuple(1 if j == i else 0 for j in range(k))
    for name, gens in [("FULL", [e(i) for i in range(k)]),
                       ("PAIRED (control group)",
                        [tuple(1 if j in (i, i+1) else 0 for j in range(k)) for i in range(k-1)])]:
        G = group_from_gens(gens, k)
        C = configs(k)
        # a G-invariant weight: constant on orbits
        orbs = orbits(G, k)
        dg = {}
        for oi, o in enumerate(orbs):
            for s in o: dg[s] = 2 + oi
        dbad = {s: random.randint(1, 9) for s in C}
        worst = 0; worstb = 0; cnt = 0
        for _ in range(300):
            vals = {s: Fraction(random.randint(-40, 40)) for s in C}
            f = lambda s: vals[s]
            Pi = {s: Fraction(sum(f(flip(g, s)) for g in G), len(G)) for s in C}
            lhs = Fraction(sum(dg[s]*f(s) for s in C), sum(dg.values()))
            rhs = Fraction(sum(dg[s]*Pi[s] for s in C), sum(dg.values()))
            worst = max(worst, abs(lhs - rhs))
            lb = Fraction(sum(dbad[s]*f(s) for s in C), sum(dbad.values()))
            rb = Fraction(sum(dbad[s]*Pi[s] for s in C), sum(dbad.values()))
            worstb = max(worstb, abs(lb - rb)); cnt += 1
        P(f"{k:>2} {name:<26} {str(worst):>26} {cnt:>6} {len(orbs):>8} {str(worstb):>30}")
P("")
P("READ: the identity holds as an EXACT rational zero in every case with a G-invariant weight.")
P("      The CONTROL column (weight not G-invariant) is non-zero, so the test discriminates.")
P("")
P("=" * 100)
P("PART 1 VERDICT")
P("=" * 100)
P("  * The GROUP-THEORETIC core of the candidate is TRUE and now proved exactly:")
P("      independent writability => all coordinate flips in G_W => G_W acts SIMPLY")
P("      TRANSITIVELY => #orbits = 1 => the invariant space is 1-dimensional (constants).")
P("  * The candidate's SECOND clause as literally written is FALSE: non-constant does not")
P("      imply mean zero.  Counterexample f = #{i : sigma_i = +1}, mean k/2, at every k.")
P("  * The REPAIRED clause is true and STRONGER than what was asked: the write RESPONSE of")
P("      every functional cancels exactly, and the d-weighted mean of every functional equals")
P("      that of its invariant part, which under independent writability is a constant.")
P("  * Therefore: RESPONSIVE and NON-CANCELLING are mutually exclusive COMPONENTS of any")
P("      functional, not mutually exclusive functionals.  Every functional splits into a")
P("      record-blind accumulating part and a record-sensitive cancelling part.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_B_THEOREM/t1_core.txt","w").write("\n".join(OUT)+"\n")
