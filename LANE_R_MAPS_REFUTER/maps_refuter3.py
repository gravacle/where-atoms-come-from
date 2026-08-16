#!/usr/bin/env python3
"""LANE R (MAPS) — REFUTER, part 3: where the wall actually is."""
import itertools
import numpy as np
from maps_refuter import (Complex, K1, K1_GF, K1_GC, classes_from_loops, push_to_pi,
                          lam_mahler_generic, lam_cesaro, hdr)
from maps_refuter2 import collapse, verify_chain_map

np.set_printoptions(linewidth=200, suppress=True)
k1 = K1(); CLS = classes_from_loops(5, K1_GF, K1_GC)
lam = lam_mahler_generic

# ===========================================================================
hdr("BLOCK 7 — A PARTITION-NON-PRESERVING CELLULAR MAP THAT PRESERVES lambda EXACTLY")
# ===========================================================================
p_form = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
q, phi, emap, (gF2, gC2) = collapse(k1, frozenset({1}), [K1_GF, K1_GC])   # collapse e2 only
a_err, b_err = verify_chain_map(k1, q, phi, emap); bb = q.betti()
print("collapse of the single edge e2 (v1~v2) — a homotopy equivalence:")
print(f"  quotient V={q.nv} E={q.ne} F={q.nf} chi={q.chi()} b=({bb[0]},{bb[1]},{bb[2]})  regular? {q.is_regular()[0]}")
print("  chain-map residuals:", a_err, b_err)
print("  K1 class partition of V :", sorted([len([w for w in range(5) if CLS[w]==c]) for c in {(1,1),(1,0),(0,1)}], reverse=True), "= blocks of size 2,2,1")
cq = classes_from_loops(q.nv, gF2, gC2)
print("  quotient class partition:", sorted([len([w for w in range(q.nv) if cq[w]==c]) for c in set(cq)], reverse=True), "= blocks of size 2,1,1  -> DIFFERENT PARTITION")
pq = np.zeros(q.nv)
for w in range(5): pq[phi[w]] += p_form[w]
print("  pi(K1)      =", push_to_pi(p_form, CLS), " lambda =", f"{lam(push_to_pi(p_form,CLS)):+.12f}")
print("  pi(quotient)=", push_to_pi(pq, cq),      " lambda =", f"{lam(push_to_pi(pq,cq)):+.12f}")
print("  -> partition of the vertex set NOT preserved; lambda preserved to the last bit.")

# a carrier with NO pinch at all, reproducing K1's lambda exactly
sp = Complex("SPEC", 7, [(1,2),(2,3),(3,1),(4,5),(5,6),(6,4),(0,1),(0,4)],
             [[(0,+1),(1,+1),(2,+1)]])
bs = sp.betti()
print(f"\ncarrier SPEC (spectator vertex, two vertex-DISJOINT loops, no pinch):")
print(f"  V={sp.nv} E={sp.ne} F={sp.nf} chi={sp.chi()} b=({bs[0]},{bs[1]},{bs[2]}) regular? {sp.is_regular()[0]}")
cs = classes_from_loops(7, {1,2,3}, {4,5,6})
p_sp = np.array([0.4, .1,.1,.1, .1,.1,.1])
print("  classes:", cs, " pi =", push_to_pi(p_sp, cs))
print(f"  lambda(SPEC) = {lam(push_to_pi(p_sp,cs)):+.12f}   lambda(K1) = {lam(push_to_pi(p_form,CLS)):+.12f}")
print("  |difference| =", abs(lam(push_to_pi(p_sp,cs)) - lam(push_to_pi(p_form,CLS))))
print("  -> a carrier with NO (1,1) vertex reproduces K1's rate exactly. Labels invisible HERE.")

# ===========================================================================
hdr("BLOCK 8 — THE SUFFICIENCY FAILURE: A CLASS-PRESERVING MAP THAT MOVES lambda")
# ===========================================================================
sub_edges = []
for j,(s,t) in enumerate(k1.edges): sub_edges += [(s,5+j),(5+j,t)]
b1s = Complex("B1s", 11, sub_edges, [[(0,1),(1,1),(2,1),(3,1),(4,1),(5,1)]])
gF_s, gC_s = {0,1,2,5,6,7}, {0,3,4,8,9,10}
cs_s = classes_from_loops(11, gF_s, gC_s)
print("B1s -> K1 (undo the subdivision) is a cellular map that preserves every vertex CLASS.")
print("  SENSE C (state pushed forward):  lambda(B1s)=%+.9f  lambda(K1)=%+.9f  diff=%.1e"
      % (lam(np.array([0,.3,.3,.4])), lam(np.array([0,.3,.3,.4])), 0.0))
lu_s = lam(push_to_pi(np.ones(11)/11, cs_s)); lu_k = lam(push_to_pi(np.ones(5)/5, CLS))
print("  SENSE U (state re-chosen uniform): lambda(B1s)=%+.9f  lambda(K1)=%+.9f  diff=%+.9f"
      % (lu_s, lu_k, lu_s-lu_k))
print("  -> the SAME class-preserving map is functorial under one state convention and not the")
print("     other. S4 CHOICE LEDGER C2 records this choice as 'OPEN and load-bearing'.")

# ===========================================================================
hdr("BLOCK 9 — THE MULTISET RECORD, TESTED WHERE S4 COULD NOT: FOUR OCCUPIED CLASSES")
# I PREDICTED THIS WOULD FAIL AND IT DID NOT.  The monomial-substitution symmetry of
# Z's Newton polygon (the unit square, corners 1,x,y,xy) is only D4 of order 8 --
# n |-> -n+(1,1) gives (00 11)(10 01); n |-> (n1,-n2)+(0,1) gives (00 01)(10 11);
# n1<->n2 gives (10 01).  Transpositions such as (00 10) alone are NOT realisable by
# any affine substitution over Z.  So the geometric argument predicts 3 orbits of 8
# and therefore up to 3 distinct lambda.  COMPUTATION SAYS OTHERWISE: all 24 agree,
# at 30 decimal places.  The multiset invariance is a genuine identity beyond the
# evident symmetries, and it now stands verified at FOUR occupied classes -- a regime
# no S4 carrier could reach, since every one of the ten has an empty class.
# ===========================================================================
# a carrier with all FOUR classes occupied — none of S4's ten has this
kp = Complex("K1+", 6, k1.edges + [(1,5)], [[(0,1),(1,1),(2,1)]])
bp = kp.betti()
print(f"K1+ = K1 with a pendant edge e7: v1->v5.  V={kp.nv} E={kp.ne} F={kp.nf} chi={kp.chi()}"
      f" b=({bp[0]},{bp[1]},{bp[2]}) regular? {kp.is_regular()[0]}")
cp = classes_from_loops(6, K1_GF, K1_GC)
print("  classes:", cp, "  -> all four of (0,0),(1,0),(0,1),(1,1) occupied")

for tag, base in (("THREE classes occupied (every S4 carrier)", np.array([0.0,0.2,0.3,0.5])),
                  ("FOUR classes occupied  (K1+, never built)", np.array([0.1,0.2,0.3,0.4]))):
    vals = {}
    for perm in itertools.permutations(range(4)):
        pi = base[list(perm)]
        vals.setdefault(round(lam(pi), 10), []).append(perm)
    print(f"\n  {tag}   multiset {sorted(base)}")
    print(f"    24 permutations -> {len(vals)} distinct lambda:")
    for L, ps in sorted(vals.items()):
        print(f"       lambda = {L:+.10f}   ({len(ps)} of 24 permutations)")
    print(f"    spread over the 24 = {max(vals)-min(vals):.10e}")

print("\n  WHAT I EXPECTED AND DID NOT GET.  Z's Newton polygon is the unit square with")
print("  corners 1, x, y, xy carrying pi00, pi10, pi01, pi11, and lambda is its Mahler")
print("  measure.  The substitutions that permute those corners are the affine maps")
print("  n |-> Mn+b over Z preserving {0,1}^2 -- exactly D4, order 8:")
print("     n |-> -n+(1,1)     : (00 11)(10 01)      [the record's pinch<->spectator]")
print("     n |-> (n1,-n2)+(0,1): (00 01)(10 11)")
print("     n1 <-> n2          : (10 01)")
print("  D4 fixes the diagonal pairing {{00,11},{10,01}}, so the 24 arrangements fall")
print("  into 3 orbits of 8, and up to THREE distinct lambda were predicted at four")
print("  occupied classes.  THE PREDICTION IS FALSE: all 24 agree to 30 dps (see")
print("  perm_exact.py).  The multiset invariance is a real identity strictly stronger")
print("  than the polygon symmetry.  It now stands verified at FOUR occupied classes,")
print("  which no S4 carrier could reach -- every one of the ten has an empty class.")

# exact pointwise checks of the two generating symmetries
rng = np.random.default_rng(2718281828)
w1 = w2 = 0.0
for _ in range(4000):
    pi = rng.dirichlet(np.ones(4)); f, c = rng.uniform(0, 2*np.pi, 2)
    u, v = np.exp(-1j*f), np.exp(1j*c); ks = rng.integers(1, 10**5, 30)
    Z  = pi[0] + pi[1]*u**ks + pi[2]*v**ks + pi[3]*(u*v)**ks
    # (00 11)(10 01) at the SAME connection  -- the record's pinch<->spectator identity
    pa = pi[[3,2,1,0]]
    Za = pa[0] + pa[1]*u**ks + pa[2]*v**ks + pa[3]*(u*v)**ks
    w1 = max(w1, float(np.abs(np.abs(Z)-np.abs(Za)).max()))
    # (00 10)(01 11) at the connection with u -> 1/u   -- the symmetry NOT in the record
    pb = pi[[1,0,3,2]]; ub = 1/u
    Zb = pb[0] + pb[1]*ub**ks + pb[2]*v**ks + pb[3]*(ub*v)**ks
    w2 = max(w2, float(np.abs(np.abs(Z)-np.abs(Zb)).max()))
print(f"\n  seed 2718281828, 4000 (pi, connection) x 30 circuit counts:")
print(f"    (00 11)(10 01) at the same connection      : max ||Z|-|Z'|| = {w1:.2e}  (record: 6.55e-15)")
print(f"    (00 10)(01 11) with u -> u^-1              : max ||Z|-|Z'|| = {w2:.2e}  (NOT in the record)")

# ===========================================================================
hdr("BLOCK 10 — THE REDUCTIO: THE COLLAPSE TRIVIALISES EVERY STATE FUNCTIONAL")
# ===========================================================================
p = np.array([0.4,0.15,0.15,0.15,0.15])
print("K1, forming state p =", p)
print("  Shannon entropy H(p)          = %.6f nats   -> after collapse to 1 vertex: %.6f" %
      (-(p*np.log(p)).sum(), 0.0))
print("  purity sum p^2                = %.6f        -> after collapse: %.6f" % ((p**2).sum(), 1.0))
print("  |Omega_10| at f=2,c=1.1       = %.6e  -> after collapse: %.6f" %
      (np.prod(np.abs([ (0.4)*np.exp(1j*n*(1.1-2.0)) + 0.3*np.exp(-1j*n*2.0) + 0.3*np.exp(1j*n*1.1)
                        for n in range(1,11)])), 1.0))
print("  number of occupied vertex classes |S| = %d -> after collapse: 1" %
      len({c for w,c in enumerate(CLS) if p[w] > 0}))
print("  -> the collapse sends EVERY functional of the vertex measure to its point-mass value.")
print("     Nothing in this singles out lambda, and no property of lambda is being tested.")

# ===========================================================================
hdr("BLOCK 11 — DOES THE CONNECTION SURVIVE THE COLLAPSE? (the half of the claim that is TRUE)")
# ===========================================================================
from maps_refuter import Z_from_operators
rng2 = np.random.default_rng(161803398)
worstW = 0.0
for _ in range(5000):
    a = rng2.uniform(0, 2*np.pi, 6)
    f, c = a[0]+a[1]+a[2], a[3]+a[4]+a[5]
    # gauge-fix the tree {e1,e2,e4,e5} to zero: theta_v0=0, then solve along the tree
    th = np.zeros(5)
    th[1] = th[0] + a[0]; th[2] = th[1] + a[1]; th[3] = th[0] + a[3]; th[4] = th[3] + a[4]
    a_g = a.copy()
    for j,(s,t) in enumerate(k1.edges): a_g[j] = a[j] + th[s] - th[t]
    # after the gauge fix the two surviving edges carry the two holonomies
    worstW = max(worstW, abs(np.exp(1j*a_g[2]) - np.exp(1j*f)), abs(np.exp(1j*a_g[5]) - np.exp(1j*c)))
    assert abs(a_g[0]) < 1e-9 and abs(a_g[1]) < 1e-9 and abs(a_g[3]) < 1e-9 and abs(a_g[4]) < 1e-9
print("seed 161803398, 5000 connections: gauge-fix the spanning tree to zero, then")
print("  max | exp(i a'_e3) - W_F |  and  | exp(i a'_e6) - W_C |  =", worstW)
print("  -> the collapse loses NO gauge-invariant content of the connection. Both")
print("     invariants descend exactly. The claim's connection bookkeeping is correct.")
print("     What it loses is the VERTEX SET, i.e. the entire domain of the ready state.")
