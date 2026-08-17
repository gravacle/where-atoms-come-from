# W-08 / M4-REFUTE-1  leg B — M4's F3 "CONTINGENCY" TEST, AUDITED.
# LENS: STEELMAN W-07.  M4-5 says W-07 sec6 committed IMP-1 by voiding a THEOREM.  Its ONE piece of
# empirical support is m4_f F3: run the same sweep on C5 and "the count 2 is achieved by TEN
# 4-block partitions and NOT the discrete one", concluding "the sweep DISCRIMINATES ... its K1
# outcome is a fact ABOUT K1, not an artefact of the sweep's design."
#
# The sweep's target is not the numeral 2.  It is S1 sec4's INVARIANT PARAMETER COUNT, which on K1
# equals 2 and is DERIVED FROM THE COMPLEX: E - rank(d1) = E - (V-1) for a connected complex.
# On C5 that number is 5 - 4 = 1, not 2.  M4 held the numeral fixed while moving the complex.
#
# ISOLATION LEDGER.  B1: the complex is held fixed at K1; the PARTITION moves (all 52).
# B2: the sweep procedure is held fixed, byte-identical function; the COMPLEX moves, and the
# target moves WITH IT because the target is a function of the complex (that is the correction).
# B3: the complex moves over 200 random connected graphs; nothing else moves.
# PRECISION: exact.  All ranks by fraction-free elimination over Q on integer matrices.  No floats.
from fractions import Fraction
from collections import Counter
import itertools, random

def partitions(c):
    if len(c)==1: yield [c]; return
    first,rest=c[0],c[1:]
    for p in partitions(rest):
        for i in range(len(p)): yield p[:i]+[[first]+p[i]]+p[i+1:]
        yield [[first]]+p

def exact_rank(M):
    A=[[Fraction(x) for x in row] for row in M]
    rows=len(A); cols=len(A[0]) if A else 0; r=0
    for c in range(cols):
        piv=next((i for i in range(r,rows) if A[i][c]!=0), None)
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]
        for i in range(rows):
            if i!=r and A[i][c]!=0:
                f=A[i][c]/A[r][c]
                A[i]=[A[i][j]-f*A[r][j] for j in range(cols)]
        r+=1
    return r

def d1_of(EDGES,V,E):
    d1=[[0]*E for _ in range(V)]
    for j,(s,t) in enumerate(EDGES): d1[t][j]+=1; d1[s][j]-=1
    return d1

def sweep(EDGES,V,E,target):
    d1=d1_of(EDGES,V,E)
    cnt=Counter(); winners=[]
    for P in partitions(list(range(V))):
        B=[[1 if v in blk else 0 for blk in P] for v in range(V)]
        M=[[sum(d1[v][j]*B[v][i] for v in range(V)) for i in range(len(P))] for j in range(E)]
        inv=E-exact_rank(M)
        cnt[inv]+=1
        if inv==target: winners.append(sorted([sorted(b) for b in P]))
    return cnt,winners

EDGES_K1=[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]
EDGES_C5=[(0,1),(1,2),(2,3),(3,4),(4,0)]

print("== B1  K1, EXACT RANK.  Reproduces W-07 and M4. ==")
tgt_K1 = 6 - exact_rank([[d1_of(EDGES_K1,5,6)[v][j] for v in range(5)] for j in range(6)])
cnt,win=sweep(EDGES_K1,5,6,2)
print(f"  distribution {dict(sorted(cnt.items()))}   winners at 2: {win}")
print(f"  K1's OWN invariant count E - rank(d1) = {tgt_K1}   (S1 sec4: 6 - 4 = 2)  -> target 2 is correct FOR K1\n")

print("== B2  M4's F3, WITH THE TARGET COMPUTED FROM THE COMPLEX INSTEAD OF COPIED FROM K1 ==")
for nm,EDGES,V,E in [("K1  (V=5,E=6)",EDGES_K1,5,6), ("C5  (V=5,E=5)",EDGES_C5,5,5)]:
    d1=d1_of(EDGES,V,E)
    r=exact_rank([[d1[v][j] for v in range(V)] for j in range(E)])
    tgt=E-r
    cnt,win=sweep(EDGES,V,E,tgt)
    cnt2,win2=sweep(EDGES,V,E,2)
    disc=[[ [v] for v in range(V)]]
    print(f"  {nm}: rank(d1)={r}  ITS OWN invariant count = E - rank = {tgt}")
    print(f"      distribution {dict(sorted(cnt.items()))}")
    print(f"      winners at ITS OWN target {tgt}: {len(win)}  -> {win}")
    print(f"      winners at M4's copied numeral 2 : {len(win2)}")
print()
print("  M4's F3 compares C5 at target 2.  C5's own invariant count is 1.  At its OWN target the")
print("  sweep on C5 returns EXACTLY ONE winner and it IS the discrete partition -- the same answer")
print("  as on K1.  M4 manufactured contingency by holding a numeral fixed while moving the object")
print("  that DEFINES the numeral.  That is this program's cardinal defect, committed by the lane")
print("  whose headline convicts W-07 of it.\n")

print("== B3  AND THE SWEEP IS FORCED ON EVERY CONNECTED GRAPH, WHICH IS STRONGER THAN W-07 SAID ==")
print("  For a connected complex and a k-block partition the quotient is connected, so")
print("  rank = k-1 exactly and invariants = E-k+1.  The target E-rank(d1) = E-V+1 is therefore")
print("  attained iff k = V: the DISCRETE partition, always, uniquely, on ANY connected graph.")
random.seed(20260816)
bad=0; tested=0
for trial in range(200):
    V=random.randint(3,6)
    all_e=[(i,j) for i in range(V) for j in range(V) if i!=j]
    E=random.randint(V-1,min(9,len(all_e)))
    EDGES=random.sample(all_e,E)
    # connectivity check
    adj={v:set() for v in range(V)}
    for a_,b_ in EDGES: adj[a_].add(b_); adj[b_].add(a_)
    seen={0}; st=[0]
    while st:
        x=st.pop()
        for y in adj[x]:
            if y not in seen: seen.add(y); st.append(y)
    if len(seen)!=V: continue
    tested+=1
    d1=d1_of(EDGES,V,E)
    tgt=E-exact_rank([[d1[v][j] for v in range(V)] for j in range(E)])
    cnt,win=sweep(EDGES,V,E,tgt)
    ok = (len(win)==1 and win[0]==[[v] for v in range(V)])
    if not ok:
        bad+=1; print(f"    COUNTEREXAMPLE V={V} E={E} EDGES={EDGES} winners={win}")
print(f"  connected random graphs tested: {tested}   graphs where the unique winner is NOT the")
print(f"  discrete partition: {bad}")
print("  ==> the 52-sweep 'could not have come out otherwise' ON ANY CONNECTED INPUT -- which is")
print("      M4's OWN stated criterion for voiding a control (m4_f F3, first line).  Applied with")
print("      the target the sweep actually uses, M4's F3 REFUTES M4-5 instead of supporting it.")
