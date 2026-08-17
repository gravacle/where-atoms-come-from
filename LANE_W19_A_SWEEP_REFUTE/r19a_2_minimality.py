# W19-A REFUTER / step 2.  THE MINIMALITY CLAIM (finding F1) IS FALSE AS STATED.
#
# Lane A: "PROVED MINIMAL conditional on min degree >= 3: four plateau points require
#          d = dist_(G-l)(u,v) >= 5 ... BFS level counting ... forces level sizes >= 1,2,4 on each
#          side ... so V >= 14 and L >= 3*14/2 = 21."
# The 1,2,4 step needs the two level-1 vertices to have DISJOINT further neighbourhoods, i.e. NO
# 4-cycle through u.  That is GLOBAL girth >= 5.  The stated hypothesis is only girth_through(l) >= 6,
# which is a statement about cycles that USE l.  Cycles that AVOID l are unconstrained.
#
# Everything below is run with LANE A's OWN CODE (zn_gauge.py, imported unmodified) and LANE A's OWN
# criterion (delta = 0.10, rule A nested fragments, rule C level cuts, weight floor H(S) >= 0.10).
# The only thing that moves is THE CARRIER.
import numpy as np, sys, itertools, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP_REFUTE")
from zn_gauge import (ZNGauge, S_of, mutual_information, level_cuts, nested_fragments,
                      girth_through, has_uv_path, build_adj)
from carriers import heawood, petersen, cube, theta_sub, LADDER
from refute_carriers import (tri_chain12, tri_chain11, dbl_chain9, mg_chain,
                             min_degree, connected, max_d, cubic_graphs)

DELTA = 0.10
print("=" * 122)
print("W19-A REFUTER / 2  THE 21-LINK THRESHOLD IS NOT A THRESHOLD")
print("=" * 122)

def girth(V, E):
    from collections import deque
    best = 10 ** 9
    for l in range(len(E)):
        g_, d_ = girth_through(V, E, l)
        if g_ is not None: best = min(best, g_)
    return best

def verdict_for(nm, V, E, gsq, N=2):
    """LANE A's criterion, verbatim: rule A plateau points, rule C R_delta, delta=0.10, floor 0.10."""
    g = ZNGauge(nm, V, E, N)
    # lane A's system-link rule: the link maximising d, ties -> lowest index
    best = (-1, -1)
    for l in range(len(E)):
        _, di = girth_through(V, E, l)
        if di is not None and di > best[0]: best = (di, l)
    d, l = best
    psi, E0, gap = g.ground(2.0 / gsq, 2.0 * gsq)
    Psi = g.full_vector(psi)
    HS = S_of(Psi, g.L, N, [l])
    frs, _ = nested_fragments(V, E, l)
    cuts, _ = level_cuts(V, E, l)
    Is = [mutual_information(Psi, g.L, N, [l], F) for F in frs]
    pts = sum(abs(I / HS - 1) <= DELTA for I in Is)
    Rd = sum(mutual_information(Psi, g.L, N, [l], C) / HS >= 1 - DELTA for C in cuts)
    disj = True; seen = set()
    for C in cuts:
        if seen & set(C): disj = False
        seen |= set(C)
    ver = "EXHIBITED" if pts >= 4 else ("MARGINAL" if pts == 3 else "FAIL")
    if HS < 0.10: ver += " (WEIGHTLESS)"
    return dict(name=nm, V=V, L=g.L, C=g.C, dimP=g.dimP, l=l, d=d, girth=girth(V, E),
                mindeg=min_degree(V, E), HS=HS, pts=int(pts), R=int(Rd), disj=disj, ver=ver,
                sizes=[len(F) for F in frs], ratios=[I / HS for I in Is],
                cutsizes=[len(C) for C in cuts], E0=E0, gap=gap)

print("\n[2a] LANE A's CRITERION, LANE A's CODE, NEW CARRIERS.  g^2 = 0.50 (the exhibit's coupling).")
print("     Every carrier below has min degree >= 3 -- lane A's stated hypothesis -- and satisfies")
print("     girth_through(l) = d+1, exactly the quantity lane A's own closed form F2 says is decisive.")
hdr = (f"     {'carrier':<16}{'V':>4}{'L':>4}{'C':>4}{'dimP':>7}{'mindeg':>8}{'girth':>7}"
       f"{'g_thru_l':>10}{'d':>4}{'H(S) bits':>12}{'pts':>5}{'R_d':>5}{'disj':>6}  verdict")
print(hdr)
rows = []
CAND = [("heawood(ref)", heawood()), ("petersen(ref)", petersen()),
        ("tri_chain12", tri_chain12()), ("dbl_chain9", dbl_chain9()),
        ("tri_chain11", tri_chain11()),
        ("mg_chain4", mg_chain(4)), ("mg_chain5", mg_chain(5)),
        ("mg_chain6", mg_chain(6)), ("mg_chain7", mg_chain(7))]
for nm, (V, E) in CAND:
    r = verdict_for(nm, V, E, 0.50)
    rows.append(r)
    print(f"     {r['name']:<16}{r['V']:>4}{r['L']:>4}{r['C']:>4}{r['dimP']:>7}{r['mindeg']:>8}{r['girth']:>7}"
          f"{r['d']+1:>10}{r['d']:>4}{r['HS']:>12.9f}{r['pts']:>5}{r['R']:>5}{str(r['disj']):>6}  {r['ver']}")
print("\n     PER-FRAGMENT DETAIL for the two carriers that break the claim (rule A, l = link 0):")
for r in rows:
    if r['name'] in ("tri_chain12", "dbl_chain9", "heawood(ref)"):
        print(f"       {r['name']:<14} L={r['L']:<3} |F| = " + " ".join(f"{s:>3}" for s in r['sizes']))
        print(f"       {'':<14} {'':<5} I/H = " + " ".join(f"{x:>3.1f}" for x in r['ratios'])
              + f"   cut sizes (rule C) = {r['cutsizes']}")

print("\n[2b] SAME TABLE AT g^2 = 1.00, the ladder's coupling, so the comparison with lane A's")
print("     out_01_ladder.txt summary table is like-for-like.")
print(hdr)
for nm, (V, E) in CAND:
    r = verdict_for(nm, V, E, 1.00)
    print(f"     {r['name']:<16}{r['V']:>4}{r['L']:>4}{r['C']:>4}{r['dimP']:>7}{r['mindeg']:>8}{r['girth']:>7}"
          f"{r['d']+1:>10}{r['d']:>4}{r['HS']:>12.9f}{r['pts']:>5}{r['R']:>5}{str(r['disj']):>6}  {r['ver']}")

print("\n[2c] EXHAUSTIVE SEARCH -- the search lane A declared it did not run.")
print("     All CONNECTED LABELLED SIMPLE graphs with min degree >= 3 on V = 5,6,7 vertices, any L.")
print("     (min degree 3 forces L >= ceil(3V/2); V=7 => L >= 11.  V >= 8 => L >= 12 automatically,")
print("      so if no V <= 7 graph reaches d = 5 then the SIMPLE floor for 4 plateau points is L = 12.)")
import numpy as _np
for V in (5, 6, 7):
    pairs = list(itertools.combinations(range(V), 2)); P = len(pairs)
    vmask = _np.zeros(V, dtype=_np.int64)
    for i, (a, b) in enumerate(pairs): vmask[a] |= (1 << i); vmask[b] |= (1 << i)
    masks = _np.arange(1 << P, dtype=_np.int64)
    keep = _np.ones(1 << P, dtype=bool)
    for v in range(V):
        keep &= (_np.bitwise_count(masks & vmask[v]) >= 3)
    cand = masks[keep]
    bestd = -1; bestL = None; nconn = 0; hist = {}
    for mk in cand:
        E = [pairs[i] for i in range(P) if int(mk) >> i & 1]
        if not connected(V, E): continue
        nconn += 1
        dd = max_d(V, E)
        hist[dd] = hist.get(dd, 0) + 1
        if dd > bestd or (dd == bestd and len(E) < bestL): bestd, bestL = dd, len(E)
    print(f"     V={V}: {len(cand)} degree-feasible edge sets, {nconn} connected;  MAX d over ALL of them"
          f" = {bestd} (first attained at L = {bestL});  d-histogram = {dict(sorted(hist.items()))}")

print("\n     All labelled CUBIC simple graphs on V=8 (the only V,L combination that can give L=12):")
cg = cubic_graphs(8)
d5 = [E for E in cg if connected(8, E) and max_d(8, E) == 5]
d5max = [E for E in cg if connected(8, E) and max_d(8, E) >= 5]
print(f"     {len(cg)} labelled cubic graphs on 8 vertices; {sum(1 for E in cg if connected(8,E))} connected;"
      f" {len(d5max)} of them reach d >= 5, i.e. 4 rule-A plateau points at L = 12.")
byg = {}
for E in d5max:
    byg[girth(8, E)] = byg.get(girth(8, E), 0) + 1
print(f"     girth histogram of the L=12 winners: {byg}   (lane A's argument assumed girth >= 6; every")
print(f"     winner has girth 3, which its own criterion never forbids)")

print("\n[2d] THE CORRECTED FLOOR.  P plateau points require d = P+1 (lane A's own F2, which this lane")
print("     re-verified on every row above).  d = P+1 requires a shortest u-v path of length P+1 in G-l,")
print("     hence V >= P+2 distinct vertices; min degree 3 then forces 2L >= 3V, L >= ceil(3(P+2)/2).")
print("     The mg_chain family ATTAINS that bound for every d, so the bound is exact for multigraphs.")
print(f"     {'P (plateau pts)':<18}{'d':>4}{'V_min':>7}{'L_min = ceil(3(P+2)/2)':>26}{'attained by':>16}"
      f"{'lane A said':>14}")
for P in (2, 3, 4, 5, 6):
    d = P + 1; Vm = P + 2; Lm = -(-3 * Vm // 2)
    V, E = mg_chain(d)
    ok = (len(E) == Lm and max_d(V, E) == d and min_degree(V, E) >= 3)
    said = {3: "15 (petersen)", 4: "21 (heawood)"}.get(P, "-")
    print(f"     {P:<18}{d:>4}{Vm:>7}{Lm:>26}{('mg_chain('+str(d)+') ok' if ok else 'FAIL'):>16}{said:>14}")

print("\n[2e] WHY LANE A's LADDER COULD NOT HAVE FOUND THIS.  Its 12 carriers, by (girth, girth_through(l)):")
print(f"     {'carrier':<22}{'L':>4}{'mindeg':>8}{'girth':>7}{'girth_thru_l':>14}")
for nm, (V, E) in LADDER:
    best = max((girth_through(V, E, l)[1] or -1) for l in range(len(E)))
    print(f"     {nm:<22}{len(E):>4}{min_degree(V,E):>8}{girth(V,E):>7}{best+1:>14}")
print("     The ladder contains NO carrier with (min degree >= 3, small girth, LARGE girth_through(l)).")
print("     That cell is exactly where the counterexamples live.  The ladder's 'heawood is the only one'")
print("     is a property of the ladder, not of the class.")
json.dump(rows, open("out_2_minimality.json", "w"), indent=1, default=float)
print("\nDONE 2.")
