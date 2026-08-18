"""PHASE A / T3.  CAPACITY = m - 1, ON ARBITRARY COMPLEXES.

W-42 proved this for planar lattice patches. The proof uses no planarity and no lattice, so it should
hold for ANY 2-complex. Test it on random complexes, on non-planar ones, and on the degenerate cases
the proof's hypothesis singles out.

THE PROOF. Regions are subsets of the plaquette set P, |P| = m, so they form GF(2)^m. The boundary
map S -> bd(S) (links lying in an ODD number of the plaquettes of S) is LINEAR over GF(2). A set of
records is simultaneously protectable iff some link lies on NONE of their boundaries, and for a fixed
link L the condition bd(S)_L = 0 is ONE linear functional f_L on GF(2)^m. So the largest protectable
independent set via L is dim ker f_L = m - rank(f_L), and
        capacity = max over L of dim ker f_L = m - min over L of rank(f_L).
rank(f_L) is 0 exactly when NO plaquette contains L, and 1 otherwise. Hence:
        capacity = m - 1   if every link lies on at least one plaquette
        capacity = m       if some link lies on none
Nothing in that argument mentions planarity, dimension, a lattice, or a gauge group.
"""
import itertools, numpy as np
rng=np.random.default_rng(11)
def rank_gf2(vecs):
    piv=[]; r=0
    for v in vecs:
        for p in piv: v=min(v,v^p)
        if v: piv.append(v); piv.sort(reverse=True); r+=1
    return r
def bd(S,PL):
    c={}
    for p in S:
        for lk in PL[p]: c[lk]=c.get(lk,0)+1
    return set(lk for lk,v in c.items() if v%2)
def capacity(L,PL):
    m=len(PL); REG=[S for r in range(1,m+1) for S in itertools.combinations(range(m),r)]
    best=0
    for lk in range(L):
        free=[S for S in REG if lk not in bd(S,PL)]
        best=max(best, rank_gf2([sum(1<<i for i in S) for S in free]))
    return best
def predict(L,PL):
    m=len(PL)
    orphan=any(all(lk not in PL[p] for p in range(m)) for lk in range(L))
    return m if orphan else m-1

cases=[]
# planar 3x3 patch (the program's carrier)
cases.append(("3x3 planar patch",12,[[0,2,6,7],[1,3,7,8],[2,4,9,10],[3,5,10,11]]))
# four disconnected squares
cases.append(("4 disconnected squares",16,[[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15]]))
# a TETRAHEDRON: 6 links, 4 triangular faces -- NOT planar-lattice, not 2D grid
cases.append(("tetrahedron (4 faces, 6 links)",6,[[0,1,2],[0,3,4],[1,3,5],[2,4,5]]))
# a CUBE surface: 12 links, 6 square faces -- closed surface, genus 0
cube=[[0,1,2,3],[4,5,6,7],[0,8,4,9],[1,9,5,10],[2,10,6,11],[3,11,7,8]]
cases.append(("cube surface (6 faces, 12 links)",12,cube))
# random complexes, including some with orphan links
for t in range(6):
    m=rng.integers(3,6); L=int(rng.integers(6,14))
    PL=[sorted(rng.choice(L,size=int(rng.integers(3,5)),replace=False).tolist()) for _ in range(m)]
    cases.append((f"random complex #{t} (m={m}, L={L})",L,PL))

print("PHASE A / T3.  capacity = m - 1 on arbitrary complexes (planarity plays no role in the proof)")
print(f"\n  {'complex':>34s} {'m':>3s} {'L':>3s} {'orphan link':>12s} {'predicted':>10s} {'measured':>9s}  {'':>6s}")
print("  "+"-"*86)
ok=True
for name,L,PL in cases:
    m=len(PL)
    orph=any(all(lk not in PL[p] for p in range(m)) for lk in range(L))
    pr=predict(L,PL); me=capacity(L,PL)
    agree = (pr==me)
    ok = ok and agree
    print(f"  {name:>34s} {m:3d} {L:3d} {str(orph):>12s} {pr:10d} {me:9d}  {'ok' if agree else 'MISMATCH'}")
print(f"\n  all cases agree with the theorem: {ok}")
print("\n  THE HYPOTHESIS DOING THE WORK -- a link on no plaquette raises capacity to m.")
PL=[[0,1,2,3],[1,2,4,5]]
for L in (6,7):
    orph=any(all(lk not in PL[p] for p in range(2)) for lk in range(L))
    print(f"    m=2, L={L}: orphan link present = {str(orph):5s}  predicted {predict(L,PL)}  "
          f"measured {capacity(L,PL)}")
