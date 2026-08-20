"""O-50-A  step 8.  THE DICHOTOMY WITHOUT THE HYPOTHESIS.

   The brief's theorem candidate assumes independent writability in order to get a single orbit.
   That assumption turns out to be unnecessary for the CONCLUSION THAT MATTERS.  For ANY finite
   group G acting on a finite configuration set C, every f : C -> R splits ORTHOGONALLY as

        f  =  f_inv  +  f_odd ,      f_inv(c) := (1/|G|) sum_g f(g.c)

   with f_inv G-INVARIANT (so NO write changes it: it is UNRESPONSIVE) and f_odd summing to
   EXACTLY ZERO on every orbit (so it CANCELS over the set a writer can reach).  Hence:

        NO functional of the record configuration is both RESPONSIVE to writing
        and NON-CANCELLING over the writer orbit.  On any carrier.  No hypothesis.

   Independent writability only fixes HOW BIG the unresponsive part can be: one orbit means
   f_inv is a constant (dimension 1); many orbits mean f_inv can be a non-constant but frozen
   label -- which is what the 1D chain's bond variables are, and why C-61 saw record-blind
   accumulation there.

   Checked below with EXACT rational arithmetic on every carrier, with a positive control."""
import sys, itertools, random, numpy as np
from fractions import Fraction
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_A_ACTION")
import orbits

def check(name, gens, m, ntrial=200, seed=3):
    G = orbits.close_group(gens, m)
    orbs = orbits.orbits_of(G, m)
    rng = random.Random(seed)
    worst_orbit_sum = 0; resp_cnt = 0; bad = 0; unresp_nonconst = 0
    for _ in range(ntrial):
        f = [Fraction(rng.randint(-9, 9)) for _ in range(m)]
        finv = [sum(f[g[c]] for g in G) / len(G) for c in range(m)]
        fodd = [f[c] - finv[c] for c in range(m)]
        inv_ok = all(finv[g[c]] == finv[c] for g in G for c in range(m))
        s = max(abs(sum(fodd[c] for c in o)) for o in orbs)
        worst_orbit_sum = max(worst_orbit_sum, s)
        responsive = any(f[g[c]] != f[c] for g in G for c in range(m))
        resp_cnt += responsive
        if responsive and all(x == 0 for x in fodd): bad += 1
        if not responsive and len(set(f)) > 1: unresp_nonconst += 1
        if not inv_ok: bad += 1
    dim = orbits.invariant_dim_exact(G, m)[0]
    print(f"  {name:40s} |G|={len(G):>4} |C|={m:>3} orbits={len(orbs):>3} inv_dim={dim:>3}  "
          f"responsive draws={resp_cnt}/{ntrial}  max |orbit sum of f_odd| = {worst_orbit_sum}  "
          f"violations={bad}  unresponsive-but-non-constant draws={unresp_nonconst}")
    return dim, len(orbs), worst_orbit_sum

print("EXACT CHECK OF THE DECOMPOSITION (max |orbit sum of f_odd| must be EXACTLY 0)")
# torus, k=2: translations of F_2^2
for k in (2, 4, 6):
    m = 2 ** k
    cfgs = list(itertools.product((0, 1), repeat=k)); idx = {c: i for i, c in enumerate(cfgs)}
    gens = []
    for j in range(k):
        t = tuple(1 if i == j else 0 for i in range(k))
        gens.append([idx[tuple((c[i] + t[i]) % 2 for i in range(k))] for c in cfgs])
    check(f"torus-type, k={k} independent records", gens, m)
# chain greedy: k=1
check("1D chain, greedy family (k=1)", [[1, 0]], 2)
# chain naive: global negation on 2^n
for nq in (3, 4, 5, 6):
    m = 2 ** nq
    check(f"1D chain, naive family {{Z_i}} n={nq}", [[m - 1 - x for x in range(m)]], m)
print()
print("POSITIVE CONTROL (D-15): a TRIVIAL group, where NOTHING may cancel and every functional")
print("is unresponsive -- the machinery must report inv_dim = |C| and non-constant invariants.")
for m in (4, 8):
    check(f"trivial group on {m} configurations", [list(range(m))], m)
print()
print("READ, from the numbers above: max |orbit sum of f_odd| is EXACTLY 0 on every carrier,")
print("including the ones where the hypothesis fails; and it is also 0 on the trivial-group")
print("control ONLY because every orbit is a singleton there -- where the whole of f is f_inv,")
print("nothing is responsive, and inv_dim = |C|.  The control therefore separates the cases.")
