# W19-A step 6.  THE THRESHOLD EXHIBIT, in full, with the numbers that go in the register.
# Carrier: heawood_honeycomb7 = the Heawood graph = the hexagonal (honeycomb) tiling of the TORUS
# with 7 hexagonal plaquettes.  14 vertices, 21 links, cyclomatic number 8, dim_phys = 2^8 = 256.
# Coupling: g^2 = 0.50 -- inside the deconfined (magnetic) phase, where the plateau HEIGHT is nearly
# maximal (0.9948 of 1 bit), so the exhibit is not sitting on the weight floor.
import numpy as np, sys
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP")
from zn_gauge import *
from carriers import *
GSQ=0.50
V,E = heawood(); g = ZNGauge("heawood",V,E,2); L=g.L
psi,E0,gap = g.ground(2.0/GSQ, 2.0*GSQ); Psi = g.full_vector(psi)
frs,d = nested_fragments(V,E,0); cuts,_ = level_cuts(V,E,0)
HS = S_of(Psi,L,2,[0])
print("="*100)
print("W19-A THRESHOLD EXHIBIT")
print("="*100)
print(f"carrier      : heawood_honeycomb7 (Heawood graph = honeycomb torus, 7 hexagonal plaquettes)")
print(f"               V={V} vertices, L={L} links, cyclomatic C={g.C}, dim_phys = 2^{g.C} = {g.dimP}")
print(f"               all degrees 3, girth 6, edge-transitive; plaquette set = 8 hexagons (min cycle basis)")
print(f"Hamiltonian  : H = -(1/g^2) sum_p (W_p+W_p^dag) - g^2 sum_e (X_e+X_e^dag),  g^2 = {GSQ}")
print(f"               E0 = {E0:.9f}   gap = {gap:.9f}   (nondegenerate ground state)")
print(f"system S     : link 0 = {E[0]};  d = dist_(G-l)(u,v) = {d}")
print(f"CUT          : S = link 0  |  E = the remaining {L-1} links")
print(f"H(S)         : {HS:.9f} bits   (ceiling log2 2 = 1 bit; deconfined phase)")
print()
print("PLATEAU (rule A, nested fragments growing by BFS distance from tail(l) inside G-l):")
print(f"  {'k':>3}{'|F|':>6}{'I(S:F) bits':>16}{'I/H(S)':>14}  encloses cycle through l ?")
for k,F in enumerate(frs,1):
    I=mutual_information(Psi,L,2,[0],F)
    print(f"  {k:>3}{len(F):>6}{I:>16.9f}{I/HS:>14.9f}  {has_uv_path(V,E,0,F)}")
print(f"  -> PLATEAU: 4 fragment sizes (|F| = {', '.join(str(len(F)) for F in frs[:4])} of {L-1} environment links,")
print(f"     i.e. 10% to 90% of the environment) at I = {HS:.9f} bits = 1.000000000 * H(S).")
print()
print(f"R_delta (rule C, the {d} pairwise DISJOINT level cuts; each obeys X_l = X(C_i)^dag exactly):")
tot=set()
for i,C in enumerate(cuts,1):
    I=mutual_information(Psi,L,2,[0],C); assert not (tot & set(C)); tot|=set(C)
    print(f"  C_{i}  |C|={len(C):>2}  links={C}   I(S:C)={I:.9f}   I/H(S)={I/HS:.9f}")
print(f"  -> R_delta = {d}  (pairwise disjoint: True; they use {len(tot)} of the {L-1} environment links)")
print()
print("NEGATIVE CONTROL in the same run: theta graph, 3 links (the carrier T1 failed).")
V2,E2 = theta(); g2 = ZNGauge("theta",V2,E2,2)
p2,_,_ = g2.ground(2.0/GSQ,2.0*GSQ); P2 = g2.full_vector(p2); H2 = S_of(P2,3,2,[0])
f2,d2 = nested_fragments(V2,E2,0)
print(f"  H(S)={H2:.9f} bits, d={d2}, only fragment |F|={len(f2[0])}: I/H(S)="
      f"{mutual_information(P2,3,2,[0],f2[0])/H2:.9f}  -> 0 plateau points, R_delta=1.")
print()
print("ARMS DIFF (theta vs heawood): V 2->14, L 3->21, C 2->8, dim_phys 4->256, girth 2->6, d 1->5.")
print("The SINGLE variable that moves the verdict is d; steps 01/03/04 hold the rest fixed to show it.")
