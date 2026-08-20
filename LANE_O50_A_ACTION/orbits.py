"""Generic, assumption-free orbit / invariant-functional analyser.
   Input: a list of CONFIGURATIONS and a list of PERMUTATIONS of them (as index maps).
   Output: the generated group order (by closure), orbits, stabiliser orders, Burnside count,
   and the EXACT dimension of the space of invariant real functionals, obtained by exact
   rational Gaussian elimination on the simultaneous fixed-point equations f(g.c) = f(c).
   Nothing here assumes the action is a translation."""
from fractions import Fraction

def close_group(gens, m):
    idp = tuple(range(m)); G = {idp}; frontier = [idp]
    while frontier:
        nf = []
        for g in frontier:
            for h in gens:
                p = tuple(h[g[i]] for i in range(m))
                if p not in G: G.add(p); nf.append(p)
        frontier = nf
        if len(G) > 400000: raise RuntimeError("group too large to close")
    return G

def orbits_of(G, m):
    seen = [False] * m; orbs = []
    for c in range(m):
        if seen[c]: continue
        o = sorted({g[c] for g in G})
        for x in o: seen[x] = True
        orbs.append(o)
    return orbs

def invariant_dim_exact(G, m):
    """dim { f : R^m | f(g.c) = f(c) for all g in G } by exact rational elimination."""
    rows = []
    for g in G:
        for c in range(m):
            if g[c] == c: continue
            r = [Fraction(0)] * m
            r[g[c]] += 1; r[c] -= 1
            rows.append(r)
    if not rows: return m, 0
    # gaussian elimination over Q
    piv = 0; pivots = []
    for col in range(m):
        p = next((i for i in range(piv, len(rows)) if rows[i][col] != 0), None)
        if p is None: continue
        rows[piv], rows[p] = rows[p], rows[piv]
        pv = rows[piv][col]
        rows[piv] = [x / pv for x in rows[piv]]
        for i in range(len(rows)):
            if i != piv and rows[i][col] != 0:
                f = rows[i][col]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[piv])]
        pivots.append(col); piv += 1
        if piv == len(rows): break
    return m - len(pivots), len(pivots)

def analyse(gens, m, labels=None):
    G = close_group(gens, m)
    orbs = orbits_of(G, m)
    stab = []
    for o in orbs:
        c = o[0]
        stab.append(sum(1 for g in G if g[c] == c))
    burnside = Fraction(sum(sum(1 for c in range(m) if g[c] == c) for g in G), len(G))
    dim, rk = invariant_dim_exact(G, m)
    return dict(group_order=len(G), n_configs=m, n_orbits=len(orbs),
                orbit_sizes=sorted(len(o) for o in orbs), stabiliser_orders=stab,
                burnside_orbit_count=str(burnside),
                transitive=(len(orbs) == 1),
                simply_transitive=(len(orbs) == 1 and len(G) == m and all(s == 1 for s in stab)),
                invariant_dim_exact=dim, orbits=orbs)
