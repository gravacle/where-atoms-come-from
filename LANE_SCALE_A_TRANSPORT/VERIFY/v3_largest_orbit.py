"""V3 -- THE 'LARGEST ORBIT' COLUMN: ASSERTED OR MEASURED?

s3_orbits.py computes the generic stabiliser ONLY for |G| <= 16 ("if G.n <= 16: ... else:
gen_triv = None") and prints '?' in the largest-orbit column for every larger carrier.  Its
own READ block then prints s['Gb'] = |G/Z| UNCONDITIONALLY under the heading 'largest orbit
(records)', and the finding's MASTER TABLE and ARM B carry that number for |G| = 32 and 64
with the gloss "(generic stabiliser verified trivial)".  ARM B -- the only arm that varies
transport at fixed record count -- is D_8xZ_2 vs ES_2^(1+4), BOTH order 32, i.e. both in the
band the script left unverified.

I run the check the lane skipped.  Generic stabiliser in Gbar = G/Z is trivial iff for EVERY
MINIMAL subgroup Z < L <= G (image of order 2 in Gbar) the fixed set is a proper subvariety:
    dim_L = sum_E sum_{sigma in Irr(L)} m_sigma(E)^2   <   dim_all = sum_E dim(E)^2.
Fixed sets are nested, so a non-trivial stabiliser contains a minimal one: minimal L is
necessary AND sufficient.  L is abelian of order 2|Z|, so this is cheap at |G| = 64 too.

POSITIVE CONTROL, SAME TABLE (D-15): abelian carriers must return NO minimal L and largest
orbit 1; and the code path can return 'NON-trivial' -- a carrier with dim_L == dim_all would
print it and refute the column.  Both outcomes are reachable.
"""
import sys, numpy as np
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT"
sys.path.insert(0, LANE)
import glib
from carriers import census, isotypic, phi
def say(*a): print(*a); sys.stdout.flush()

def subgroup_grp(G, K): return glib.Grp("K", sorted(K), lambda a, b: int(G.mt[a, b]))

def dimL(G, ce, K):
    Kg = subgroup_grp(G, K)
    cl, chi, d, cls_of = Kg.chars()
    tot = 0
    for v in (-2, -1, 0):
        if ce['dims'][v] == 0: continue
        cf = np.array([ce['chis'][v][Kg.el[c[0]]] for c in cl])
        m, dd = Kg.decompose(cf)
        assert np.max(np.abs(m - np.round(m))) < 1e-6, "non-integer multiplicity"
        m = np.round(m).astype(int)
        assert int(np.sum(dd * m)) == ce['dims'][v], "SELF-CHECK sum d*m != dim E"
        tot += int(np.sum(m ** 2))
    return tot

say("="*136)
say("V3   GENERIC STABILISER AND LARGEST ORBIT, COMPUTED ON EVERY CARRIER INCLUDING |G| = 32 AND 64")
say("="*136)
say("")
say(f"  {'carrier':<13}{'|G|':>5}{'abel':>6}{'|Z|':>5}{'|G/Z|':>7}{'#minimal L':>12}{'dim_all':>10}"
    f"{'max dim_L':>11}{'max dim_L/dim_all':>19}{'generic stab':>14}{'largest orbit':>15}{'lane printed':>14}")
lane_says = {}   # what the finding's master table asserts
for G in glib.ladder(64):
    ce = census(G)
    Z = set(G.centre); n = G.n
    dall = sum(ce['dims'][v] ** 2 for v in (-2, -1, 0))
    mins = set()
    for g in range(n):
        if g in Z: continue
        if int(G.mt[g, g]) not in Z: continue          # image in G/Z must have order 2
        L = frozenset(Z | {int(G.mt[z, g]) for z in Z})
        assert len(L) == 2 * len(Z), f"{G.name}: <Z,g> not of order 2|Z|"
        mins.add(L)
    best = 0
    for L in mins:
        best = max(best, dimL(G, ce, L))
    triv = all(dimL(G, ce, L) < dall for L in mins) if mins else True
    GZ = n // len(Z)
    largest = GZ if triv else "NON-TRIVIAL->smaller"
    say(f"  {G.name:<13}{n:>5}{str(G.abelian):>6}{len(Z):>5}{GZ:>7}{len(mins):>12}{dall:>10}"
        f"{best:>11}{(best/dall if dall else 0):>19.6f}"
        f"{('trivial' if triv else 'NON-trivial'):>14}{str(largest):>15}{GZ:>14}")
say("")
say("  'lane printed' is the column the finding reports as 'largest orbit'.  It equals |G/Z| by")
say("  construction in s3_orbits.py's READ block, whether or not the stabiliser was checked.")
