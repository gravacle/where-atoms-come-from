"""PHASE C / C1.  DOES THE STORAGE/LEGIBILITY SPLIT REDUCE TO LOCALITY?

W-48 killed the decoherence half of the composed law (it is exp(-Gamma T)). The one surviving
candidate for a distinguishing prediction is the storage/legibility split: a region holds |A|-1
independent records but exports only ~perimeter of them.

BEFORE TAKING THAT ANYWHERE, THE OBLIGATION FROM THE PLAN: establish that it predicts something the
existing account does NOT. There is an obvious deflation to rule out first.

THE DEFLATION. An operator supported strictly INSIDE a region commutes with every operator outside
it. That is plain locality -- microcausality, disjoint support -- and it is not news. If legibility
is exactly "the number of independent records whose support touches the region's boundary", then the
split is locality wearing gauge-theoretic clothing and there is no prediction here.

THE TEST. For many shapes compare
    legibility  = rank over GF(2) of the outside-inside incidence
    locality    = number of independent region-combinations whose boundary touches the outside
If they agree everywhere, the split IS locality. Where they differ, that difference is the only place
any non-trivial content can live, and it must be characterised rather than hand-waved.
"""
import itertools, numpy as np
def rank_gf2(rows):
    piv=[]; r=0
    for v in rows:
        for p in piv: v=min(v,v^p)
        if v: piv.append(v); piv.sort(reverse=True); r+=1
    return r
def grid(NX,NY):
    def hid(i,j): return j*NX+i
    def vx(i,j):  return NX*(NY+1)+ j*(NX+1)+i
    return {(i,j):[hid(i,j),hid(i,j+1),vx(i,j),vx(i+1,j)] for j in range(NY) for i in range(NX)}
def analyse(shape,NX,NY):
    PL=grid(NX,NY); A=list(shape); S=set(A)
    out=[c for c in PL if c not in S]
    idx={c:k for k,c in enumerate(A)}
    # legibility: rank of the outside-inside incidence
    rows=[]
    for q in out:
        bq=set(PL[q]); v=0
        for c in A:
            if len(bq & set(PL[c]))%2: v|=(1<<idx[c])
        rows.append(v)
    leg=rank_gf2(rows)
    # locality: the subspace of region-combinations whose boundary touches NO outside link.
    # outside links = links of A that also belong to some outside plaquette.
    shared=set()
    for q in out:
        for lk in PL[q]:
            for c in A:
                if lk in PL[c]: shared.add(lk)
    inv=[]
    for r in range(1,len(A)+1):
        for comb in itertools.combinations(range(len(A)),r):
            cnt={}
            for t in comb:
                for lk in PL[A[t]]: cnt[lk]=cnt.get(lk,0)+1
            bd=set(lk for lk,v in cnt.items() if v%2)
            if not (bd & shared): inv.append(sum(1<<t for t in comb))
    loc=len(A)-rank_gf2(inv)
    return leg,loc,len(A)
def blk(nx,ny,ox=2,oy=2): return [(ox+i,oy+j) for i in range(nx) for j in range(ny)]
cases=[]
for n in (2,3,4): cases.append((f"square {n}x{n}",blk(n,n),n+6,n+6))
for (a,b) in ((1,4),(2,4),(2,5),(3,4)): cases.append((f"rect {a}x{b}",blk(a,b),a+6,b+6))
cases.append(("L pentomino",[(2,2),(3,2),(4,2),(2,3),(2,4)],9,9))
cases.append(("plus pentomino",[(3,2),(2,3),(3,3),(4,3),(3,4)],9,9))
cases.append(("T tetromino",[(2,2),(3,2),(4,2),(3,3)],9,9))
cases.append(("S tetromino",[(2,2),(3,2),(3,3),(4,3)],9,9))
cases.append(("3x3 minus centre",[(i,j) for i in range(2,5) for j in range(2,5) if (i,j)!=(3,3)],9,9))
cases.append(("U shape",[(2,2),(3,2),(4,2),(2,3),(4,3)],9,9))
print("PHASE C / C1.  IS THE STORAGE/LEGIBILITY SPLIT JUST LOCALITY?")
print(f"\n  {'shape':>20s} {'|A|':>4s} {'legibility':>11s} {'locality bound':>15s} {'differ?':>8s}")
print("  "+"-"*64)
diffs=[]
for name,shape,NX,NY in cases:
    leg,loc,na=analyse(shape,NX,NY)
    d=leg-loc; diffs.append(d)
    print(f"  {name:>20s} {na:4d} {leg:11d} {loc:15d} {('' if d==0 else f'{d:+d}'):>8s}")
print()
print(f"  cases where legibility != the locality bound: {sum(1 for d in diffs if d!=0)} of {len(diffs)}")
print(f"  differences observed: {sorted(set(diffs))}")
print()
if all(d==0 for d in diffs):
    print("  VERDICT: legibility IS the locality bound in every case tested. The storage/legibility")
    print("  split is microcausality expressed in gauge-theoretic language. It predicts nothing the")
    print("  existing account does not, and MUST NOT be taken to data as though it did.")
else:
    print("  VERDICT: they differ somewhere. That difference is the ONLY place non-trivial content")
    print("  can live and must be characterised before anything is claimed.")
