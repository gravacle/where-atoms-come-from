"""Chamon model in the Bravyi-Leemhuis-Terhal (arXiv:1006.4871) convention:
qubits on ONE parity sublattice (fcc), the six-body XXYYZZ operator centred on the
sites of the OTHER sublattice.  Their result: for all L_i = 2 p_i even,
GSD = 2^{4 gcd(p_x,p_y,p_z)}, i.e. k = 4 gcd(p_x,p_y,p_z).  We test that."""
import sys
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_O78_FRACTON_NONABELIAN")
from gf2 import *
from math import gcd

def chamon_fcc(Lx,Ly,Lz):
    site = {}
    for x in range(Lx):
        for y in range(Ly):
            for z in range(Lz):
                if (x+y+z) % 2 == 0:
                    site[(x,y,z)] = len(site)
    n = len(site)
    q = lambda x,y,z: site[(x%Lx, y%Ly, z%Lz)]
    g = []
    for x in range(Lx):
        for y in range(Ly):
            for z in range(Lz):
                if (x+y+z) % 2 == 0: continue      # centres on the odd sublattice
                X = 0; Z = 0
                for s in (1,-1):
                    X ^= (1 << q(x+s,y,z))                                  # X
                    X ^= (1 << q(x,y+s,z)); Z ^= (1 << q(x,y+s,z))          # Y
                    Z ^= (1 << q(x,y,z+s))                                  # Z
                g.append(mk(X,Z,n))
    return g, n

print("Chamon on the fcc sublattice, cubic torus L x L x L, L = 2p:")
print(f"{'L':>3s} {'p':>3s} {'n':>6s} {'gens':>6s} {'rank':>6s} {'k':>5s} {'4*gcd(p,p,p)=4p':>16s} {'match':>7s} {'commute':>8s}")
for L in (2,4,6,8,10):
    p = L//2
    g,n = chamon_fcc(L,L,L)
    k,r = code_k(g,n)
    bad = all_commute(g,n)
    pred = 4*p
    print(f"{L:>3d} {p:>3d} {n:>6d} {len(g):>6d} {r:>6d} {k:>5d} {pred:>16d} "
          f"{'YES' if k==pred else 'no':>7s} {bad:>8d}")

print("\nAnisotropic tori, where the gcd formula has real content (k = 4 gcd(px,py,pz)):")
print(f"{'Lx,Ly,Lz':>12s} {'px,py,pz':>10s} {'n':>6s} {'k':>5s} {'4*gcd':>7s} {'match':>7s} {'commute':>8s}")
for (Lx,Ly,Lz) in [(4,4,6),(4,6,8),(4,6,10),(6,8,10),(4,8,12),(6,6,10)]:
    px,py,pz = Lx//2, Ly//2, Lz//2
    g,n = chamon_fcc(Lx,Ly,Lz)
    k,r = code_k(g,n)
    bad = all_commute(g,n)
    pred = 4*gcd(gcd(px,py),pz)
    print(f"{str((Lx,Ly,Lz)):>12s} {str((px,py,pz)):>10s} {n:>6d} {k:>5d} {pred:>7d} "
          f"{'YES' if k==pred else 'no':>7s} {bad:>8d}")

print("\nCSS test on the fcc Chamon code:")
g,n = chamon_fcc(6,6,6)
rx,rz,rt,css = css_split(g,n)
print(f"  L=6: rank(pure-X part)={rx}  rank(pure-Z part)={rz}  rank(S)={rt}  CSS = {css}")
print("  (a CSS code would need rank_X + rank_Z = rank_S)")
