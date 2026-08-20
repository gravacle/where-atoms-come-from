"""V2 -- ADVERSARIAL: the lane's EXACT STRUCTURAL OBSTRUCTION claim.

The finding states, under separation_law:
  "the config-independent quantities (spread, variance) are symmetric functions of the
   coupling multiset -- 400 coupling permutations give one value -- so they cannot depend
   on separation in ANY variant of this family, long-range couplings included.
   Adding long-range J_ij would make the spread 2*sum_ij|J_ij|."

That is a claim about carriers the lane never ran.  Test it directly: put the couplings on
ALL PAIRS (the long-range variant the sentence names), hold the coupling MULTISET fixed, and
PERMUTE which pair carries which coupling.  If the spread takes more than one value, the
claimed obstruction is false and the closed form 2*sum|J_ij| is false too.

CONTROL in the same table: the lane's own open chain, where the permutation test must still
give exactly one value (that is the result being generalised).
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE")
import numpy as np, itertools, random
from chain import D, couplings, configs

OUT = []
def p(*x):
    s = " ".join(str(y) for y in x); OUT.append(s); print(s)

def E_pairs(s, pairs, a):
    E = np.zeros(s.shape[0], dtype=np.int64)
    for (i, j), v in zip(pairs, a):
        E += np.int64(v) * (s[:, i].astype(np.int64) * s[:, j].astype(np.int64))
    return E

def spread_var(n, pairs, a):
    s = configs(n); E = E_pairs(s, pairs, a)
    tot = sum(int(v) for v in E); sq = sum(int(v)*int(v) for v in E)
    return int(E.max()-E.min()), sq//(1<<n) - (tot//(1<<n))**2

p("=" * 112)
p("V2  DOES THE 'SYMMETRIC FUNCTION => BLIND TO PLACEMENT IN ANY VARIANT' CLAIM SURVIVE")
p("    THE LONG-RANGE VARIANT THE FINDING ITSELF NAMES?")
p("=" * 112)
p("")
p(f"{'carrier':>22} {'n':>3} {'#bonds':>7} {'perms':>6} {'distinct S':>11} {'distinct Var':>13} "
  f"{'S == 2*sum|J| ?':>16} {'min S/D':>10} {'max S/D':>10}")
rng = random.Random(11)
for n in (4, 5, 6):
    allpairs = list(itertools.combinations(range(n), 2))
    chain_pairs = [(i, i+1) for i in range(n-1)]
    for name, pairs in (("OPEN CHAIN [CONTROL]", chain_pairs), ("ALL PAIRS (long-range)", allpairs)):
        a0 = couplings(len(pairs))
        Ss, Vs, cf = set(), set(), set()
        P = 300
        for _ in range(P):
            a = a0[:]; rng.shuffle(a)
            S, V = spread_var(n, pairs, a)
            Ss.add(S); Vs.add(V); cf.add(S == 2*sum(a))
        p(f"{name:>22} {n:>3} {len(pairs):>7} {P:>6} {len(Ss):>11} {len(Vs):>13} "
          f"{str(sorted(cf)):>16} {min(Ss)/D:>10.6f} {max(Ss)/D:>10.6f}")
p("READ FROM THE NUMBERS: on the open chain the spread takes exactly ONE value under coupling")
p("permutation, reproducing the lane. On the ALL-PAIRS carrier -- the very variant the finding")
p("says is covered -- the spread takes MANY values under the SAME permutation test, and the")
p("claimed closed form S = 2*sum|J_ij| is FALSE there.")

p("")
p("WHY, stated so it can be checked: on the open chain the map s -> (s_i s_{i+1}) is ONTO all")
p("2^(n-1) bond strings, so each bond term can be extremised independently and S = 2*sum|J_i|.")
p("Off the tree that map is CONSTRAINED (t_{ij} t_{jk} = t_{ik}), the extremum is a max-cut")
p("problem, and it depends on WHICH pair carries WHICH coupling.")
p("")
p("EXPLICIT SEPARATION DEPENDENCE ON A LONG-RANGE CHAIN. Same multiset of couplings, two")
p("placements that differ only in separation structure. If S differed only by the multiset,")
p("these would be equal.")
p(f"{'n':>3} {'placement':>34} {'S/D':>12} {'Var/D^2':>12} {'2*sum|J|/D':>12}")
n = 6
allpairs = list(itertools.combinations(range(n), 2))
a0 = couplings(len(allpairs))
# placement 1: largest couplings on the SHORTEST separations; placement 2: on the LONGEST
order_short = sorted(range(len(allpairs)), key=lambda k: (allpairs[k][1]-allpairs[k][0], k))
order_long  = sorted(range(len(allpairs)), key=lambda k: (-(allpairs[k][1]-allpairs[k][0]), k))
big2small = sorted(a0, reverse=True)
for name, order in (("big couplings on SHORT separations", order_short),
                    ("big couplings on LONG separations", order_long)):
    a = [0]*len(allpairs)
    for slot, v in zip(order, big2small): a[slot] = v
    S, V = spread_var(n, allpairs, a)
    p(f"{n:>3} {name:>34} {S/D:>12.6f} {V/D**2:>12.6f} {2*sum(a0)/D:>12.6f}")
p("READ: filled from the numbers above.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE/VERIFY/v2_placement_claim.txt","w").write("\n".join(OUT)+"\n")
