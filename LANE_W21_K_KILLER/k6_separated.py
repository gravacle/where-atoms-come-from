"""
k6_separated.py -- LANE W21-K, BLOCK 6.

Casini-Huerta-Rosabal's own resolution of the ambiguity is a CONTINUUM statement: the ambiguity
becomes local on the boundary, so relative entropy and mutual information of SEPARATED regions
are finite, universal and gauge independent.  k2b block 2.3 measured mutual information between
COMPLEMENTARY (touching) regions and found it as ambiguous as the entropy -- spread 3.9989 bits.
That is not a refutation of CHR; complementary regions share their boundary and there is nothing
for a boundary-local term to cancel against.

This block moves ONE VARIABLE: SEPARATION.  Two disjoint regions with a corridor between them,
same carrier, same states, same enumeration rule.  If the spread of I collapses with separation,
CHR's mechanism is visible finitely and the ambiguity is a boundary artefact.  If it does not,
the finite lattice cannot see the mechanism and that is worth knowing before anything is built
on top of mutual information here.
"""
import numpy as np, itertools
from klib import (Z2Gauge, rref_basis, span_elements, in_span, pack, unpack, sympl_perp,
                  intersect, gauge_invariant_subspace)
from klog import Sector, LogAlgebra, _reduce

def line(c="="): print(c * 100)


def all_subspaces(q):
    seen = {frozenset([0])}; out = [[]]; cur = [[]]
    while cur:
        nxt = []
        for b in cur:
            for v in range(1, 1 << q):
                if in_span(v, b):
                    continue
                nb = rref_basis(list(b) + [v])
                key = frozenset(span_elements(nb))
                if key in seen:
                    continue
                seen.add(key); out.append(nb); nxt.append(nb)
        cur = nxt
    return out


def grid(R, C):
    idx = lambda i, j: i * C + j
    E = []
    for i in range(R):
        for j in range(C - 1):
            E.append((idx(i, j), idx(i, j + 1)))
    for i in range(R - 1):
        for j in range(C):
            E.append((idx(i, j), idx(i + 1, j)))
    return R * C, E


nv, E = grid(3, 4)
G = Z2Gauge(nv, E, "grid3x4")
n = G.L
S = Sector(G)
gp = rref_basis([pack(g, n) for g in G.gauss])
GI = gauge_invariant_subspace(G)
GIred = [v for v in rref_basis([_reduce(x, gp) for x in GI]) if v]

line()
print(f"carrier {G.name}: V={G.V} L={G.L} cycle_dim={G.cycle_dim} physical dim D={S.D}")
print(f"edges {list(enumerate(G.edges))}")


def vminmax(vs):
    internal, boundary, external = G.region(vs)
    imask = sum(1 << l for l in internal)
    cyc = [c for c in span_elements(G.cycles) if c and (c & ~imask) == 0]
    Vmin = rref_basis([pack(G.X(l), n) for l in internal] + [pack(G.W(c), n) for c in cyc])
    return Vmin, internal, boundary


def choices(vsA, vsOther):
    """All algebras assignable to region vsA when the *other* region is vsOther:
    V_min(A) <= B <= (V_min(other))' cap GI."""
    VminA, intA, bndA = vminmax(vsA)
    VminO, intO, bndO = vminmax(vsOther)
    VmaxA = intersect(sympl_perp(VminO, n), GI, n)
    rmin = [v for v in rref_basis([_reduce(x, gp) for x in VminA]) if v]
    rmax = [v for v in rref_basis([_reduce(x, gp) for x in VmaxA]) if v]
    assert all(in_span(v, rmax) for v in rmin)
    q = len(rmax) - len(rmin)
    reps = []
    cur = list(rmin)
    for v in rmax:
        if not in_span(v, cur):
            reps.append(v); cur = rref_basis(cur + [v])
    out = []
    for b in all_subspaces(q):
        gens = list(rmin)
        for bv in b:
            acc = 0
            for i in range(q):
                if (bv >> i) & 1:
                    acc ^= reps[i]
            gens.append(acc)
        out.append(LogAlgebra(gens, S, gp))
    return out, q, intA, bndA


CASES = [
    ("COMPLEMENTARY (touching)", [0, 1, 4, 5], [2, 3, 6, 7, 8, 9, 10, 11]),
    ("SEPARATED  (corridor >= 1 vertex)", [0, 1, 4, 5], [10, 11]),
    ("SEPARATED  (corridor, both 2x1)", [0, 1], [10, 11]),
]

for label, VA, VB in CASES:
    line("-")
    linksAB = [l for l, (a, b) in enumerate(G.edges)
               if (a in VA and b in VB) or (b in VA and a in VB)]
    algA, qA, intA, bndA = choices(VA, VB)
    algB, qB, intB, bndB = choices(VB, VA)
    print(f"{label}")
    print(f"  A = vertices {VA} (internal links {intA}), B = vertices {VB} (internal links {intB})")
    print(f"  links directly joining A and B: {linksAB}  -> {'ADJACENT' if linksAB else 'SEPARATED'}")
    print(f"  q_A = {qA} -> {len(algA)} algebras for A ;  q_B = {qB} -> {len(algB)} algebras for B ;"
          f"  {len(algA)*len(algB)} pairs")
    CAP = 40
    if len(algA) > CAP or len(algB) > CAP:
        rr0 = np.random.default_rng(66)
        if len(algA) > CAP:
            algA = [algA[i] for i in rr0.choice(len(algA), CAP, replace=False)]
        if len(algB) > CAP:
            algB = [algB[i] for i in rr0.choice(len(algB), CAP, replace=False)]
        print(f"  (CONFOUND, RECORDED: exhaustive enumeration is out of reach here; "
              f"{CAP} algebras sampled uniformly per side, seed 66.  A SAMPLED spread is a LOWER")
        print(f"   BOUND on the true spread, so a large spread is safe and a small one is not.)")
    for g2 in (0.10, 0.50, 1.00, 3.00):
        psi, e0, gap = S.ground_state(g2)
        rho = np.outer(psi, psi.conj()); cache = {}
        sA = [A.entropy(rho, cache) for A in algA]
        sB = [B.entropy(rho, cache) for B in algB]
        pairs = list(itertools.product(range(len(algA)), range(len(algB))))
        if len(pairs) > 600:
            rr = np.random.default_rng(6)
            pairs = [pairs[i] for i in rr.choice(len(pairs), 600, replace=False)]
        Iv = []
        joincache = {}
        for i, j in pairs:
            key = (tuple(algA[i].basis), tuple(algB[j].basis))
            sj = joincache.get(key)
            if sj is None:
                J = LogAlgebra(algA[i].basis + algB[j].basis, S, gp)
                sj = J.entropy(rho, cache); joincache[key] = sj
            Iv.append(round(sA[i] + sB[j] - sj, 9))
        u = sorted(set(Iv))
        print(f"    g^2={g2:.2f}  {len(pairs)} pairs -> {len(u)} distinct I ; "
              f"min={min(u):.9f} max={max(u):.9f} SPREAD={max(u)-min(u):.9f} bits")

line()
print("""READ IT BOTH WAYS.  A spread that shrinks with separation says CHR's boundary-local
mechanism is visible at this size and mutual information is the robust quantity.  A spread that
does not shrink says this lattice is too small for the mechanism -- the 'corridor' here is one
vertex wide, which is not a separation in any continuum sense -- and mutual information carries
the same ambiguity as entropy for every purpose this program has.  The numbers above decide it
for THIS carrier only; nothing here extrapolates to the continuum, where CHR's statement lives.""")
line()
print("END BLOCK 6")
