"""ADVERSARIAL VERIFY 4 -- the two exhibits the finding leans on, recomputed from scratch."""
import numpy as np, itertools
from fractions import Fraction as Fr

def table(dims, sgns):
    """dims: dimensions of the minimal projections; sgns: list of sign vectors (records)."""
    k=len(sgns); m=len(dims)
    d={}
    for s in itertools.product((1,-1),repeat=k): d[s]=0
    for j in range(m):
        cfgj=tuple(sgns[i][j] for i in range(k))
        d[cfgj]+=dims[j]
    return d

def GW(d,k):
    out=[]
    for msk in range(1<<k):
        eps=tuple((msk>>i)&1 for i in range(k))
        if all(d[s]==d[tuple(-x if e else x for e,x in zip(eps,s))] for s in d): out.append(eps)
    return out

print("="*100)
print("EXHIBIT A -- the finding's 'worked exhibit', d(++,+-,-+,--) = (3,1,1,3)")
print("="*100)
# realise it: single shell dim 8, 4 minimal projections of dims 3,1,1,3
dims=[3,1,1,3]; sg=[[1,1,-1,-1],[1,-1,1,-1]]
d=table(dims,sg)
print("   table:",{s:d[s] for s in [(1,1),(1,-1),(-1,1),(-1,-1)]})
for i,S in enumerate([(0,),(1,),(0,1)]):
    tr=sum(d[s]*np.prod([s[j] for j in S]) for s in d)
    print("   Tr(P_E R_%s) = %d   clause (iv) satisfied: %s"%("".join(str(j+1) for j in S),tr,tr==0))
G=GW(d,2); print("   G_W =",G,"  order",len(G))
seen=set(); orbs=[]
for s in d:
    if s in seen: continue
    o=sorted({tuple(-x if e else x for e,x in zip(g,s)) for g in G})
    for x in o: seen.add(x)
    orbs.append(o)
print("   orbits:",orbs,"  #orbits =",len(orbs),"  dim(invariants) = ",len(orbs))
print("   >>> MATCHES the finding.  BUT clause (v) is NEVER checked on this carrier: it is a")
print("       diagonal H on C^8 with NO regions, NO locality, NO tensor factorisation.  The")
print("       objects called 'records' here satisfy (i)-(iv) only -- 4 of the 5 clauses.")

print()
print("="*100)
print("EXHIBIT B -- the finding's smallest k=3 'open edge' witness")
print("   shell dim 6, all minimal projections rank 1,")
print("   R_1=(+,+,+,-,-,-)  R_2=(+,+,-,+,-,-)  R_3=(+,+,-,-,+,-)")
print("="*100)
dims=[1]*6
sg=[[1,1,1,-1,-1,-1],[1,1,-1,1,-1,-1],[1,1,-1,-1,1,-1]]
d=table(dims,sg)
print("   table:",{s:d[s] for s in sorted(d)})
for S in [(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)]:
    tr=sum(d[s]*np.prod([s[j] for j in S]) for s in d)
    print("     Tr(P_E R_%s) = %3d   balanced/(iv): %s"%("".join(str(j+1) for j in S),tr,tr==0))
G=GW(d,3); print("   G_W =",G,"  order",len(G),"  => EVERY chi_S is G_W-invariant")
print("   >>> CONFIRMED: chi_1 is a NON-CONSTANT G_W-invariant whose operator R_1 is BALANCED,")
print("       so clause (iv) guarantees an admissible U flips R_1 -- it just cannot fix R_2,R_3.")
print("       The lane's 'open edge' is REAL and correctly identified.")

print()
print("="*100)
print("EXHIBIT C -- is the 486,706-family search doing any work?")
print("="*100)
print("   Part 4.3's own argument: 'indep = k' MEANS each coordinate flip is realisable;")
print("   the coordinate flips GENERATE (Z_2)^k; (Z_2)^k acting on itself by translation is")
print("   regular; a regular action has one orbit.  Every step is a definition or Lagrange.")
print("   So the forbidden cell (indep=k AND orbits>1) is empty BY DEFINITION, before any")
print("   carrier is built.  Enumerating 486,706 families cannot populate a cell whose")
print("   emptiness is definitional -- the search has no power to refute and no control can")
print("   give it any.  D-15's 'populated control cells' show the code runs; they do not show")
print("   the null is informative.")
