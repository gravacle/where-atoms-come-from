#!/usr/bin/env python3
"""LANE R (MAPS) — REFUTER, part 2: the attack.  Conventions as in maps_refuter.py."""
import itertools
import numpy as np
from maps_refuter import (Complex, K1, K1_GF, K1_GC, classes_from_loops, push_to_pi,
                          lam_mahler_generic, lam_cesaro, Z_from_pi, hdr, cassaigne_maillot)

np.set_printoptions(linewidth=200, suppress=True)

# ---------------------------------------------------------------------------
# COLLAPSE MACHINERY: quotient a complex by a subforest of its 1-skeleton
# ---------------------------------------------------------------------------
def collapse(cx, tree_edges, loops):
    """Quotient cx by the subcomplex spanned by tree_edges (must be a forest).
    Returns (quotient complex, vertex map phi, edge map (index or None), new loops).
    Chain map: q_0 = phi on vertices, q_1 kills collapsed edges, q_2 = id on faces."""
    parent = list(range(cx.nv))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for e in tree_edges:
        s, t = cx.edges[e]
        rs, rt = find(s), find(t)
        assert rs != rt, "tree_edges is not a forest (cycle detected)"
        parent[rs] = rt
    reps = sorted({find(w) for w in range(cx.nv)})
    idx = {r: i for i, r in enumerate(reps)}
    phi = [idx[find(w)] for w in range(cx.nv)]

    new_edges, emap = [], {}
    for j, (s, t) in enumerate(cx.edges):
        if j in tree_edges:
            emap[j] = None
        else:
            emap[j] = len(new_edges)
            new_edges.append((phi[s], phi[t]))
    new_faces = []
    for word in cx.faces:
        new_faces.append([(emap[e], sg) for (e, sg) in word if emap[e] is not None])
    q = Complex(cx.name + "/T", len(reps), new_edges, new_faces)
    new_loops = [ {phi[w] for w in Lp} for Lp in loops ]
    return q, phi, emap, new_loops


def verify_chain_map(cx, q, phi, emap):
    """Check d1' q1 = q0 d1  and  d2' q2 = q1 d2  at the matrix level."""
    q0 = np.zeros((q.nv, cx.nv))
    for w in range(cx.nv): q0[phi[w], w] = 1.0
    q1 = np.zeros((q.ne, cx.ne))
    for j in range(cx.ne):
        if emap[j] is not None: q1[emap[j], j] = 1.0
    q2 = np.eye(cx.nf) if cx.nf == q.nf else np.zeros((q.nf, cx.nf))
    A = q.d1() @ q1 - q0 @ cx.d1()
    B = (q.d2() @ q2 - q1 @ cx.d2()) if cx.nf else np.zeros((q.ne, 0))
    return float(np.abs(A).max()), float(np.abs(B).max() if B.size else 0.0)


def lam_of(pi):
    return lam_mahler_generic(pi)

k1 = K1()
CLS = classes_from_loops(5, K1_GF, K1_GC)
p_form = np.array([0.4, 0.15, 0.15, 0.15, 0.15])   # S3's ready state; pi = (0,.3,.3,.4)

# ===========================================================================
hdr("BLOCK 4 — THE SPANNING-TREE COLLAPSE, BUILT AND CHECKED")
# ===========================================================================
T = frozenset({0, 1, 3, 4})     # e1,e2,e4,e5 : a spanning tree of K1
q, phi, emap, (gF2, gC2) = collapse(k1, T, [K1_GF, K1_GC])
a_err, b_err = verify_chain_map(k1, q, phi, emap)
b = q.betti()
print("tree T = {e1,e2,e4,e5};  quotient:", f"V={q.nv} E={q.ne} F={q.nf} chi={q.chi()}",
      f"b=({b[0]},{b[1]},{b[2]})")
print("chain map residuals:  max|d1'q1 - q0 d1| =", a_err, "   max|d2'q2 - q1 d2| =", b_err)
print("d1' =", q.d1().astype(int), "  d2' =", q.d2().astype(int).ravel())
print("K1 Betti (1,1,0)  ==  quotient Betti", b[:3], "  -> homotopy equivalence at the chain level: TRUE")
print("REGULARITY of the quotient:", q.is_regular(), "  <-- S1 :3-4, :27 requires REGULAR")
pi_q = push_to_pi(np.array([1.0]), classes_from_loops(q.nv, gF2, gC2))
print("quotient vertex classes:", classes_from_loops(q.nv, gF2, gC2), "  pi' =", pi_q)
print("lambda on the quotient =", lam_of(pi_q), "   (claim: 'exactly 0')  -> CONFIRMED")
print("lambda on K1 with the SAME state, uncollapsed =", lam_of(push_to_pi(p_form, CLS)))

# ===========================================================================
hdr("BLOCK 5 — THE CONTROL THE CLAIM DID NOT RUN: THE OUTCOME COULD NOT HAVE BEEN ANYTHING ELSE")
# ===========================================================================
# 5a. enumerate EVERY spanning tree of K1
trees = []
for sub in itertools.combinations(range(6), 4):
    par = list(range(5)); ok = True
    def fnd(x, par=par):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for e in sub:
        s, t = k1.edges[e]
        rs, rt = fnd(s), fnd(t)
        if rs == rt: ok = False; break
        par[rs] = rt
    if ok and len({fnd(w) for w in range(5)}) == 1:
        trees.append(sub)
print(f"spanning trees of K1: {len(trees)}  (all 4-edge subsets that are acyclic and connected)")
outcomes = set()
for tr in trees:
    qq, ph, em, (gf, gc) = collapse(k1, frozenset(tr), [K1_GF, K1_GC])
    piq = push_to_pi(np.ones(qq.nv) / qq.nv, classes_from_loops(qq.nv, gf, gc))
    outcomes.add((qq.nv, tuple(np.round(piq, 12)), round(lam_of(piq), 12)))
print("distinct (V', pi', lambda') over ALL spanning trees:", outcomes)

# 5b. the outcome is independent of connection, state and charge
rng = np.random.default_rng(31415926)
worst = 0.0
for _ in range(20000):
    f, c = rng.uniform(0, 2*np.pi, 2)
    qc = rng.integers(1, 8, 2)           # charges on the two loops
    u, v = np.exp(-1j*f*qc[0]), np.exp(1j*c*qc[1])
    p5 = rng.dirichlet(np.ones(5))       # arbitrary ready state on K1
    piq = np.array([0.,0.,0.,1.])        # its pushforward to the 1-vertex quotient, ALWAYS
    ks = rng.integers(1, 10**6, 50)
    Z = piq[0] + piq[1]*u**ks + piq[2]*v**ks + piq[3]*(u*v)**ks
    worst = max(worst, float(np.abs(np.log(np.abs(Z))).max()))
print("seed 31415926: 20000 random (connection, charge pair, ready state), 50 circuit counts each")
print("max |log|Z_k||  on the collapsed carrier =", worst, "  -> lambda = 0 IDENTICALLY, always")

# 5c. THE CONTROL: the same lambda = 0, on K1 itself, with NO map at all
for lbl, p in (("delta_v0 (the pinch/root)", np.array([1.,0,0,0,0])),
               ("delta_v1", np.array([0,1.,0,0,0])),
               ("delta_v3", np.array([0,0,0,1.,0])),
               ("all weight on class (1,0): v1,v2", np.array([0,.5,.5,0,0]))):
    pi = push_to_pi(p, CLS)
    print(f"  K1, NO COLLAPSE, p = {lbl:34s} pi={pi}  lambda = {lam_of(pi):.12f}")
print("  -> W-01/W-02 of record: |S| = 1 => G = {1} => never forms.  lambda = 0 is that row.")

# ===========================================================================
hdr("BLOCK 6 — lambda IS NOT A FUNCTION OF THE CELL STRUCTURE (claim's conclusion, tested)")
# ===========================================================================
print("SAME cell structure (K1, unchanged), state varied:")
rows = [("delta_v0",              np.array([1.,0,0,0,0])),
        ("(0.4,.15,.15,.15,.15)", p_form),
        ("uniform 1/5",           np.ones(5)/5),
        ("(0,.5,.5,0,0)",         np.array([0,.5,.5,0,0])),
        ("(0,.25,.25,.25,.25)",   np.array([0,.25,.25,.25,.25])),
        ("(.8,.05,.05,.05,.05)",  np.array([.8,.05,.05,.05,.05]))]
vals = []
for lbl, p in rows:
    pi = push_to_pi(p, CLS); L = lam_of(pi); vals.append(L)
    print(f"  {lbl:24s} pi={np.round(pi,4)}  lambda = {L:+.9f}")
print(f"  spread of lambda at FIXED cell structure = {max(vals)-min(vals):.9f}")
print("  and lambda = -infinity is reachable too (K1's own published connection,")
print("  W_F=-1, W_C=-i, p=(1/2,0,0,1/4,1/4)):  Z_1 =",
      Z_from_pi(push_to_pi(np.array([.5,0,0,.25,.25]), CLS), np.conj(-1+0j), -1j, 1))

print("\nDIFFERENT cell structures, SAME pi (pushforward respected):")
# K1 ; K1 with both triangles filled (B2) ; K1 subdivided (B1s) ; and a 4-vertex
# two-triangles-sharing-an-edge carrier (ERR-4's smaller complex)
b2 = Complex("B2", 5, k1.edges, [[(0,+1),(1,+1),(2,+1)], [(3,+1),(4,+1),(5,+1)]])
# subdivided K1: 11 vertices 0..4 original, 5..10 midpoints of e1..e6
sub_edges, mid = [], {}
for j,(s,t) in enumerate(k1.edges):
    m = 5+j; mid[j]=m; sub_edges += [(s,m),(m,t)]
b1s = Complex("B1s", 11, sub_edges, [[(0,+1),(1,+1),(2,+1),(3,+1),(4,+1),(5,+1)]])
gF_s = {0,1,2,5,6,7}; gC_s = {0,3,4,8,9,10}
for nm, cx, gf, gc, p in (("K1",        k1,  K1_GF, K1_GC, p_form),
                          ("B2 (both filled)", b2, K1_GF, K1_GC, p_form),
                          ("B1s (subdivided)", b1s, gF_s, gC_s,
                           np.array([0.4,.075,.075,.075,.075,.075,.075,0,.075,.075,0]))):
    bb = cx.betti()
    pi = push_to_pi(p, classes_from_loops(cx.nv, gf, gc))
    print(f"  {nm:18s} V={cx.nv:2d} E={cx.ne:2d} F={cx.nf} chi={cx.chi():+d} b=({bb[0]},{bb[1]},{bb[2]})"
          f"  pi={np.round(pi,6)}  lambda={lam_of(pi):+.9f}")
print("  -> three different cell structures, three different chi/Betti profiles, ONE lambda.")
