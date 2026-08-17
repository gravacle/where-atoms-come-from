# W19-A step 5.  RULE DEPENDENCE, AND THE FAILURE MODE.
#  5a  RULE A' -- nested BALLS around the whole link l (grow from BOTH endpoints).  This is a third
#      natural rule and it is HARSHER than rule A (grow from one endpoint).  Reported because the
#      verdict is rule-dependent and hiding that would be a confound.
#  5b  RULE B across the FULL ladder, including the subdivided theta graphs: does the unbiased rule
#      also get fooled by subdivision?
#  5c  The FAILURE MODE vs g^2.  The plateau itself is coupling-independent (step 02).  Is the SIZE
#      OF THE BREAK -- what I(S:F) jumps to once F encloses a cycle through l -- coupling-dependent?
import numpy as np, sys, json, itertools, time
from collections import deque
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP")
from zn_gauge import *
from carriers import *
DELTA=0.10; np.random.seed(20260817)
print("="*118); print("W19-A / 05 RULE DEPENDENCE AND FAILURE MODE"); print("="*118)

def ball_fragments(V,E,l):
    """RULE A'. F_k = links of G-l within link-distance k-1 of EITHER endpoint of l."""
    adj=build_adj(V,E); a,b=E[l]; dist={a:0,b:0}; dq=deque([a,b])
    while dq:
        x=dq.popleft()
        for (y,i) in adj[x]:
            if i==l: continue
            if y not in dist: dist[y]=dist[x]+1; dq.append(y)
    out=[]
    for k in range(1,V+2):
        F=[e for e,(x,y) in enumerate(E) if e!=l and (dist.get(x,10**9)<=k-1 or dist.get(y,10**9)<=k-1)]
        if out and F==out[-1]: continue
        out.append(sorted(F))
        if len(F)==len(E)-1: break
    return out

FULL = LADDER
print("\n[5a] RULE A' (nested balls around l, growing from BOTH endpoints).  Z_2, g^2=1.")
print(f"     {'carrier':<22}{'L':>4}{'d':>4}{'H(S)':>11}{'A pts':>7}{'A-prime pts':>13}   A-prime (|F|, I/H(S))")
rowsA=[]
for nm,(V,E) in FULL:
    if 2**len(E)>4_500_000: continue
    g=ZNGauge(nm,V,E,2); L=g.L; psi,_,_=g.ground(2.0,2.0); Psi=g.full_vector(psi); HS=S_of(Psi,L,2,[0])
    frA,d=nested_fragments(V,E,0); ptsA=sum(abs(mutual_information(Psi,L,2,[0],F)/HS-1)<=DELTA for F in frA)
    frB=ball_fragments(V,E,0); rs=[(len(F),mutual_information(Psi,L,2,[0],F)/HS) for F in frB]
    ptsB=sum(1 for _,r in rs if abs(r-1)<=DELTA)
    print(f"     {nm:<22}{L:>4}{d:>4}{HS:>11.6f}{ptsA:>7}{ptsB:>13}   " + " ".join(f"({n},{r:.4f})" for n,r in rs))
    rowsA.append(dict(name=nm,L=L,d=d,ptsA=int(ptsA),ptsAp=int(ptsB)))
print("     READ: rule A' is strictly harsher.  NO carrier in the ladder reaches 4 points under rule A'.")
print("     Under rule A' the plateau length is about ceil((d-1)/2), so girth ~10 would be needed;")
print("     the smallest min-degree-3 girth-10 graph is the (3,10)-cage, 70 vertices / 105 links -- FAR")
print("     over the exact-diagonalisation ceiling (2^105).  Recorded as a ROUTE CLOSED BY RULE CHOICE.")

print("\n[5b] RULE B (uniformly random fragments) across the FULL ladder, 48 samples/size.")
print(f"     {'carrier':<22}{'L':>4}{'C':>4}{'dimP':>7}{'d':>4}{'H(S)':>11}{'ruleB pts':>11}  verdict")
NS=48; out5b={}
for nm,(V,E) in FULL:
    if 2**len(E)>4_500_000: continue
    g=ZNGauge(nm,V,E,2); L=g.L; psi,_,_=g.ground(2.0,2.0); Psi=g.full_vector(psi); HS=S_of(Psi,L,2,[0])
    env=[e for e in range(L) if e!=0]; _,d=nested_fragments(V,E,0); means=[]
    for m in range(1,len(env)+1):
        combos=list(itertools.combinations(env,m))
        picks=combos if len(combos)<=NS else [tuple(np.random.choice(env,m,replace=False)) for _ in range(NS)]
        means.append(float(np.mean([mutual_information(Psi,L,2,[0],sorted(F))/HS for F in picks])))
    pts=sum(1 for x in means if abs(x-1)<=DELTA)
    ver="EXHIBITED" if pts>=4 else ("MARGINAL" if pts==3 else "FAIL")
    if HS<0.10: ver+=" (WEIGHTLESS)"
    print(f"     {nm:<22}{L:>4}{g.C:>4}{g.dimP:>7}{d:>4}{HS:>11.6f}{pts:>11}  {ver}")
    out5b[nm]=dict(L=L,d=d,HS=HS,pts=pts,means=means)
print("     CONFOUND ON THIS BLOCK: rule-B counts are SAMPLING NOISY at the +/-1 level, because a size")
print("     is counted whenever the sampled MEAN falls inside [0.9,1.1] and the mean crosses 1 smoothly.")
print("     Step 04a, an independent draw at 32 samples/size, gave heawood 5 rather than 4.  Only the")
print("     ORDER of the counts is load-bearing, not their exact value.")

print("\n[5c] FAILURE MODE vs g^2.  carrier=heawood, l=0.  F_5 (the first rule-A fragment that encloses")
print("     a cycle through l) and the C1 enclosure arm B.  ONE VARIABLE MOVED: g^2.")
V,E=heawood(); g=ZNGauge("heawood",V,E,2); L=g.L
frA,d=nested_fragments(V,E,0)
def shortest_path_links(V,E,l):
    a,b=E[l]; adj=build_adj(V,E); prev={a:None}; dq=deque([a])
    while dq:
        x=dq.popleft()
        for (y,i) in adj[x]:
            if i==l or y in prev: continue
            prev[y]=(x,i); dq.append(y)
    out=[]; cur=b
    while prev[cur] is not None:
        x,i=prev[cur]; out.append(i); cur=x
    return sorted(out)
sp=shortest_path_links(V,E,0); rest=[e for e in range(L) if e!=0 and e not in sp]
armB=sorted(sp+rest[:14-len(sp)])
print(f"     {'g^2':>7}{'H(S)':>12}{'I/H(F4) no cycle':>19}{'I/H(F5) cycle':>16}{'I/H(armB) cycle |F|=14':>25}")
for gsq in (0.20,0.30,0.50,0.70,1.00,1.50,3.00):
    psi,_,_=g.ground(2.0/gsq,2.0*gsq); Psi=g.full_vector(psi); HS=S_of(Psi,L,2,[0])
    r4=mutual_information(Psi,L,2,[0],frA[3])/HS
    r5=mutual_information(Psi,L,2,[0],frA[4])/HS
    rb=mutual_information(Psi,L,2,[0],armB)/HS
    print(f"     {gsq:>7.2f}{HS:>12.8f}{r4:>19.9f}{r5:>16.9f}{rb:>25.9f}")
print("     READ: the plateau value (no enclosed cycle) is 1.000000000 at every coupling -- ALGEBRAIC.")
print("     The BREAK value (enclosed cycle) moves with g^2 -- DYNAMICAL.  The only dynamical content")
print("     of this criterion sits in the height H(S) and in the size of the break, not in the plateau.")

print("\n[5d] LOCATING THE CROSSOVER.  carrier=heawood, l=0.  ONE VARIABLE MOVED: g^2 (fine scan).")
print("     H(S) is the plateau HEIGHT -- the only part of the criterion with dynamical content.")
print(f"     {'g^2':>7}{'gap':>11}{'<W_p>':>10}{'<X_e>':>10}{'H(S) bits':>12}   phase read")
prev=None; fine=[]
for gsq in [0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.90,1.00,1.10,1.30]:
    psi,E0,gap=g.ground(2.0/gsq,2.0*gsq); Psi=g.full_vector(psi); HS=S_of(Psi,L,2,[0])
    m=g._digits(); W=0.0
    for a in g.avec:
        ph=(m@np.array(a,dtype=np.int64))%2; W+=float((psi**2*np.cos(np.pi*ph)).sum())
    W/=len(g.avec)
    pw=np.array([2**c for c in range(g.C)],dtype=np.int64); X=0.0
    for t in g.tvec:
        t=np.array(t,dtype=np.int64)
        X += 1.0 if not t.any() else float((psi*psi[((m+t)%2)@pw]).sum())
    X/=len(g.tvec)
    tag = "deconfined (magnetic)" if W>X else "confined (electric)"
    print(f"     {gsq:>7.2f}{gap:>11.4f}{W:>10.5f}{X:>10.5f}{HS:>12.8f}   {tag}")
    fine.append(dict(gsq=gsq,gap=gap,W=W,X=X,HS=HS))
gmin=min(fine,key=lambda r:r["gap"]); 
cross=[r for r in fine if r["W"]<=r["X"]][0]
half =[r for r in fine if r["HS"]<0.5]
print(f"     gap minimum at g^2 = {gmin['gsq']:.2f} (gap {gmin['gap']:.4f}) -- the finite-carrier crossover.")
print(f"     <W_p> falls below <X_e> first at g^2 = {cross['gsq']:.2f}.")
print(f"     H(S) first drops below 0.5 bits at g^2 = {half[0]['gsq']:.2f}; below the 0.10 weight floor")
print(f"     at g^2 = {[r['gsq'] for r in fine if r['HS']<0.10][0] if any(r['HS']<0.10 for r in fine) else float('nan')}.")
print("     READ: the RECORD'S CONTENT switches off across the confinement crossover; the record's")
print("     SHAPE (the plateau) does not notice the crossover at all.")

json.dump(dict(ruleAp=rowsA, ruleB=out5b, fine=fine), open("out_05_ruledep.json","w"), indent=1)
print("\nDONE 05.")
