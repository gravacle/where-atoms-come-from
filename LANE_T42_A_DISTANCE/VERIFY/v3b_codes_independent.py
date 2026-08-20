"""VERIFY check 3b (independent, corrected): [[n,n-2,2]] CONFIGURATION venue.
Writers = N(S) enumerated brute force; config class of u = pairings against the k Zlogs
(computed by record_model.symplectic_logicals, same as the lane -- the pipeline downstream
of it is fully independent: no affine solver, no span_all).
Checks: n=6 all 15 nonzero classes at 2; n=8 spectrum {0,2,4} with 28 at 2 and 35 at 4.
EXTRA PROBE: swap the roles of X and Z logicals (a DIFFERENT Lagrangian) at n=8 and report
whether the class-distance multiset changes (is the venue basis-dependent beyond relabeling?)"""
import sys
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals
def pcnt(x): return bin(x).count("1")
def v2m(vec):
    m=0
    for i,b in enumerate(vec):
        if b: m|=1<<i
    return m
ok=True
def gate(lbl,b,detail=""):
    global ok
    print(("PASS " if b else "FAIL "),lbl,detail); ok&=b
def sp(u,v,n):
    mask=(1<<n)-1
    return (pcnt((u&mask)&(v>>n))+pcnt((u>>n)&(v&mask)))%2

for n in (6,8):
    full=(1<<n)-1
    k=n-2
    stab=[[1]*n+[0]*n,[0]*n+[1]*n]
    pairs=symplectic_logicals(stab,n)
    assert len(pairs)==k
    Xl=[v2m(X) for X,Z in pairs]; Zl=[v2m(Z) for X,Z in pairs]
    evens=[x for x in range(1<<n) if pcnt(x)%2==0]
    best={}
    bestX={}
    for x in evens:
        for z in evens:
            u=x|(z<<n)
            w=pcnt(x|z)
            t=tuple(sp(u,zl,n) for zl in Zl)
            if t not in best or w<best[t]: best[t]=w
            tx=tuple(sp(u,xl,n) for xl in Xl)
            if tx not in bestX or w<bestX[tx]: bestX[tx]=w
    from collections import Counter
    cnt=Counter(best.values())
    print("n=%d config-class distance histogram (Z-side labels): %s"%(n,dict(sorted(cnt.items()))))
    if n==6:
        gate("n=6 all 15 nonzero classes at 2", cnt=={0:1,2:15})
    else:
        gate("n=8 spectrum {0,2,4} with 28 at 2, 35 at 4", cnt=={0:1,2:28,4:35})
    cntX=Counter(bestX.values())
    print("n=%d SWAPPED-Lagrangian (X-side labels) histogram: %s"%(n,dict(sorted(cntX.items()))))
    print("n=%d multiset same under Lagrangian swap: %s"%(n,cnt==cntX))
print("V3B OVERALL:", "PASS" if ok else "FAIL")
