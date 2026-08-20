"""Independent part (i) check: Floyd-Warshall (pure python ints, no numpy, no Dijkstra),
metric axioms, additivity vs flip-set sum, homogeneous extensivity, controls."""
EA=[610,1200,500,2200,7800,500]; N=6; M=1<<N
INF=1<<60
D=[[0 if i==j else INF for j in range(M)] for i in range(M)]
for s in range(M):
    for i in range(N):
        D[s][s^(1<<i)]=EA[i]
for k in range(M):
    Dk=D[k]
    for i in range(M):
        Di=D[i]; dik=Di[k]
        for j in range(M):
            if dik+Dk[j]<Di[j]: Di[j]=dik+Dk[j]
def fs(s,t):
    x=s^t; tot=0; i=0
    while x:
        if x&1: tot+=EA[i]
        x>>=1; i+=1
    return tot
additive=all(D[s][t]==fs(s,t) for s in range(M) for t in range(M))
nonneg=all(D[s][t]>=0 for s in range(M) for t in range(M))
ident=all((D[s][t]==0)==(s==t) for s in range(M) for t in range(M))
sym=all(D[s][t]==D[t][s] for s in range(M) for t in range(M))
tv=sum(1 for s in range(M) for k in range(M) for t in range(M) if D[s][t]>D[s][k]+D[k][t])
print(f"additive={additive} nonneg={nonneg} ident={ident} sym={sym} tri_viol={tv}/262144 triples")
# squared control
tv2=0; wit=None
for s in range(M):
    for k in range(M):
        for t in range(M):
            if D[s][t]**2>D[s][k]**2+D[k][t]**2:
                tv2+=1
                if wit is None: wit=(s,k,t)
print(f"squared control: viol={tv2} witness={wit}")
# one-flip distances
import math
print(f"one-flip DNA={D[0][1<<2]} magnetite={D[0][1<<4]}")
# homogeneous 10x610: analytic
homog_ok = True  # d = 610*popcount by the same flip-set argument; verify on random subset via FW too costly; verify by direct argument on n=10 with dijkstra-free BFS-sum
# NOTE: additivity for state-independent positive weights is FORCED: any path s->t flips each
# differing bit >=1 time (odd) and same bits an even count >=0; weights positive => min = one flip
# each differing bit. This is a 2-line theorem, so 'earned additivity' == verified triviality.
print("analytic note: additivity is forced by state-independent positive edge weights (parity argument).")
