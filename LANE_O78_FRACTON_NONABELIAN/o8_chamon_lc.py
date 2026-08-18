"""Is the (non-CSS) Chamon code local-Clifford equivalent to a CSS code?
Exhaustive over sublattice-uniform single-qubit Cliffords on the fcc Chamon code.
Positive control: the XZZX code, which IS local-Clifford equivalent to CSS, must be
found by the SAME sublattice-uniform machinery."""
import sys, itertools
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_O78_FRACTON_NONABELIAN")
from gf2 import *
from o7_chamon_fcc import chamon_fcc

GL2 = [(a,b,c,d) for a in range(2) for b in range(2) for c in range(2) for d in range(2)
       if (a*d ^ b*c) == 1]
IDENT = GL2.index((1,0,0,1)); HAD = GL2.index((0,1,1,0))

def apply_lc(gens, n, choice):
    out = []
    for g in gens:
        X = xpart(g,n); Z = zpart(g,n); nX = 0; nZ = 0
        for i in range(n):
            xi = (X>>i)&1; zi = (Z>>i)&1
            if not (xi or zi): continue
            a,b,c,d = GL2[choice[i]]
            if (a*xi ^ b*zi): nX |= (1<<i)
            if (c*xi ^ d*zi): nZ |= (1<<i)
        out.append(mk(nX,nZ,n))
    return out

def is_css(gens, n): return css_split(gens, n)[3]

print("POSITIVE CONTROL -- XZZX toric code, known Hadamard-equivalent to CSS.")
def xzzx(L):
    n = L*L; idx = lambda x,y: (x%L)*L + (y%L); g = []
    for x in range(L):
        for y in range(L):
            g.append(mk((1<<idx(x,y))|(1<<idx(x+1,y+1)), (1<<idx(x+1,y))|(1<<idx(x,y+1)), n))
    return g, n
xg, xn = xzzx(4)
lab = [ (i//4 + i%4) % 2 for i in range(xn) ]
hit = None
for combo in itertools.product(range(6), repeat=2):
    if is_css(apply_lc(xg, xn, [combo[lab[i]] for i in range(xn)]), xn): hit = combo; break
print(f"  XZZX L=4: 2-class sublattice-uniform search over 36 bases -> "
      f"{'CSS FOUND '+str(hit) if hit else 'NOT FOUND'}   SELF-CHECK: {'PASS' if hit else 'FAIL'}\n")

print("CHAMON (fcc convention, matches Bravyi-Leemhuis-Terhal k = 4 gcd):")
for L in (4, 6):
    gens, n = chamon_fcc(L,L,L)
    sites = []
    for x in range(L):
        for y in range(L):
            for z in range(L):
                if (x+y+z)%2==0: sites.append((x,y,z))
    for nc, cls, label in [
        (2, lambda p: (p[0])%2, "x mod 2"),
        (4, lambda p: (p[0]%2) + 2*(p[1]%2), "(x,y) mod 2"),
        (4, lambda p: (p[0]+p[1]+p[2])//1 % 4, "(x+y+z) mod 4"),
        (8, lambda p: (p[0]%2)+2*(p[1]%2)+4*(p[2]%2), "(x,y,z) mod 2 octant"),
    ]:
        lab = [cls(p) % nc for p in sites]
        ncl = max(lab)+1
        hit = None; tried = 0
        for combo in itertools.product(range(6), repeat=ncl):
            tried += 1
            if is_css(apply_lc(gens, n, [combo[lab[i]] for i in range(n)]), n): hit = combo; break
        print(f"  L={L} n={n:3d}  classes={ncl} ({label:20s}) EXHAUSTIVE {tried:7d}/{6**ncl:7d} -> "
              f"{'CSS FOUND '+str(hit) if hit else 'NO CSS FORM in this family'}")
print("""
  SCOPE.  This is exhaustive over sublattice-uniform local Cliffords, not over all 6^n
  assignments (6^32 at L=4) and not over finite-depth Clifford CIRCUITS.  It is evidence,
  not proof, that Chamon is outside the CSS class.  What IS proved unconditionally here
  is the [[5,1,3]] case (o8_no_css_5.py): no CSS code with those parameters exists at all.""")
