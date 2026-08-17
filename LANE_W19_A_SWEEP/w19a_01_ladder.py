# W19-A step 1.  THE LADDER.  One variable moved per row: THE CARRIER.  Everything else fixed:
#   N=2, g^2=1 (mag=2/g^2=2, elec=2g^2=2), plaquette set = minimum weight cycle basis,
#   system-link rule = the link maximising d = dist_{G-l}(tail,head) (ties -> lowest index),
#   fragment rule A = nested-by-BFS-distance-from-tail(l), fragment rule C = the d level cuts.
#
# PLATEAU CRITERION, FIXED BEFORE ANY MEASUREMENT (delta = 0.10):
#   H_S = S(rho_l) in bits.  WEIGHT FLOOR: the plateau is WEIGHTLESS unless H_S >= 0.10*log2(N).
#   plateau point  = a fragment F_k of rule A with | I(S:F_k)/H_S - 1 | <= delta.
#   EXHIBITED (strict)   iff #plateau points >= 4      <- the brief's "three or four points"
#   EXHIBITED (marginal) iff #plateau points == 3
#   R_delta = # pairwise DISJOINT fragments with I(S:F)/H_S >= 1 - delta  (lower bound from rule C).
import numpy as np, sys, json
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP")
from zn_gauge import *
from carriers import *
from collections import deque

DELTA = 0.10
GSQ   = 1.0
MAG, ELEC = 2.0/GSQ, 2.0*GSQ

def pick_link(V, edges):
    best=(-1,-1)
    for l in range(len(edges)):
        gi,di = girth_through(V,edges,l)
        if di is not None and di>best[0]: best=(di,l)
    return best[1], best[0]

print("="*118)
print(f"W19-A / 01 LADDER.  Z_2, g^2={GSQ} (mag={MAG}, elec={ELEC}), delta={DELTA}.  ONE VARIABLE MOVED PER ROW: the carrier.")
print("="*118)
rows=[]
for nm,(V,E) in LADDER:
    L=len(E); C=L-V+1
    if 2**L > 4_500_000:
        print(f"\n### {nm}: L={L} -> 2^{L} full vector EXCEEDS CEILING 4.5e6.  SKIPPED (ceiling recorded).")
        continue
    g = ZNGauge(nm,V,E,2)
    l,d = pick_link(V,E)
    psi,E0,gap = g.ground(MAG,ELEC)
    Psi = g.full_vector(psi)
    HS  = S_of(Psi,L,2,[l])
    frs,_  = nested_fragments(V,E,l)
    cuts,_ = level_cuts(V,E,l)
    print(f"\n### {nm}:  V={V} L={L} C={C} dim_phys={g.dimP}  |plaq|={len(g.plaq)} (lengths {[int(abs(p).sum()) for p in g.plaq]})")
    print(f"    system link l={l}={E[l]}   girth_through={d+1}   d=dist_(G-l)(u,v)={d}   E0={E0:.6f} gap={gap:.6f}")
    print(f"    H(S)=S(rho_l) = {HS:.6f} bits   (weight floor 0.10; {'WEIGHTED' if HS>=0.10 else 'WEIGHTLESS'})")
    print(f"    RULE A  nested fragments F_k = links of G-l touching BFS-ball of radius k-1 around tail(l):")
    pts=0
    for k,F in enumerate(frs,1):
        I = mutual_information(Psi,L,2,[l],F)
        r = I/HS if HS>1e-9 else float('nan')
        cyc = has_uv_path(V,E,l,F)
        ok = abs(r-1)<=DELTA
        pts += ok
        print(f"      k={k}  |F|={len(F):<3} I(S:F)={I:.6f}  I/H(S)={r:.4f}  encloses_cycle_through_l={str(cyc):<5} "
              f"{'PLATEAU' if ok else ''}")
    print(f"    RULE C  the {len(cuts)} disjoint level cuts (each satisfies X_l = X(C_i)^dag exactly):")
    Rd=0
    for i,Ci in enumerate(cuts,1):
        I = mutual_information(Psi,L,2,[l],Ci)
        r = I/HS if HS>1e-9 else float('nan')
        Rd += (r >= 1-DELTA)
        print(f"      C_{i}  |C|={len(Ci):<3} links={Ci}  I(S:C)={I:.6f}  I/H(S)={r:.4f}")
    # disjointness audit
    allp=set(); dis=True
    for Ci in cuts:
        if allp & set(Ci): dis=False
        allp |= set(Ci)
    verdict = "EXHIBITED" if pts>=4 else ("MARGINAL" if pts==3 else "FAIL")
    if HS < 0.10: verdict += " (WEIGHTLESS: H(S) below floor)"
    print(f"    >>> plateau points = {pts}   R_delta = {Rd}   cuts pairwise disjoint = {dis}   VERDICT = {verdict}")
    rows.append(dict(name=nm,V=V,L=L,C=C,dimP=g.dimP,link=l,d=d,girth=d+1,HS=HS,pts=int(pts),R=int(Rd),verdict=verdict))

print("\n"+"="*118)
print("SUMMARY TABLE  (Z_2, g^2=1)")
print(f"{'carrier':<22}{'V':>4}{'L':>4}{'C':>4}{'dimP':>7}{'link':>6}{'girth_l':>9}{'d':>4}{'H(S) bits':>12}{'plateau pts':>13}{'R_delta':>9}   verdict")
for r in rows:
    print(f"{r['name']:<22}{r['V']:>4}{r['L']:>4}{r['C']:>4}{r['dimP']:>7}{r['link']:>6}{r['girth']:>9}{r['d']:>4}{r['HS']:>12.6f}{r['pts']:>13}{r['R']:>9}   {r['verdict']}")
print("\nPREDICTION UNDER TEST (stated in the module docstring, before running):")
print("  plateau points = d - 1 = girth_through(l) - 2   and   R_delta = d = girth_through(l) - 1.")
ok = all(r['pts']==r['d']-1 and r['R']==r['d'] for r in rows)
print(f"  holds on every row above: {ok}")
for r in rows:
    if not (r['pts']==r['d']-1 and r['R']==r['d']):
        print(f"    DEVIATION: {r['name']}  pts={r['pts']} expected {r['d']-1};  R={r['R']} expected {r['d']}")
json.dump(rows, open("out_01_ladder.json","w"), indent=1)
print("\nDONE 01.")
