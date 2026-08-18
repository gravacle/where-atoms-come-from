"""W-38.  IS THE REGION GIVEN, OR IS IT SELECTED?

The principal's reading: a REGION is another way of thinking about a BOUNDARY. On a graph that is
an identity -- a region IS a cycle that separates -- so the question "which region carries the
record" and "which boundary is selected" are the same question, and W-34's sieve already answers a
version of it. W-34 ran the sieve ONCE, with the environment on the cut, and the boundary of the
WHOLE patch won by 28x. The registrar noted in passing that moving the bath moves the pointer and
did not follow it.

FOLLOW IT. Every subset S of the 4 plaquettes defines a REGION, and the product of its plaquettes is
that region's BOUNDARY loop (discrete Stokes). There are 15 non-empty regions on this carrier:
4 single plaquettes, 6 pairs, 4 triples, 1 whole. Run the sieve with the environment in DIFFERENT
PLACES and ask which region's boundary the dynamics protects.

  If the SAME region always wins -> the region is a property of the carrier. It was given.
  If the winner TRACKS the environment -> the region is not given. What counts as a region is
  decided by where the environment couples, and "region" is a derived notion.

PREDICTION, declared before running: the winner is the region whose boundary has the LEAST support
on the bath's links, so the selected region is the complement of where the environment is.
NO VERDICT IS HARD-CODED. The decay rate of every one of the 15 boundaries is printed for every
bath placement, and the ranking is read off.
"""
import itertools, numpy as np
def build(V,E,N):
    st=[s for s in itertools.product(range(N),repeat=len(E))
        if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
               -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%N==0 for v in range(len(V)))]
    return st,{s:i for i,s in enumerate(st)}
def Zop(st,links,N):
    w=np.exp(2j*np.pi/N)
    return np.diag([w**(sum(s[k] for k in links)%N) for s in st]).astype(complex)
def Move(st,idx,mv,N):
    D=len(st); M=np.zeros((D,D),complex)
    for j,s in enumerate(st):
        t=list(s)
        for k,sg in mv: t[k]=(t[k]+sg)%N
        t=tuple(t)
        if t in idx: M[idx[t],j]=1.0
    return M
def compose(ps):
    acc={}
    for p in ps:
        for k,sg in p: acc[k]=acc.get(k,0)+sg
    return [(k,s) for k,s in acc.items() if s!=0]

V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E)
hid=lambda i,j:j*2+i; vx=lambda i,j:6+j*3+i
P=[[(hid(i,j),+1),(vx(i+1,j),+1),(hid(i,j+1),-1),(vx(i,j),-1)] for j in range(2) for i in range(2)]
CENTER=vid[(1,1)]; CUT=[k for k,(a,b) in enumerate(E) if a==CENTER or b==CENTER]
PERIM=[k for k in range(L) if k not in CUT]
N=2; st,idx=build(V2,E,N); D=len(st); Id=np.eye(D,dtype=complex)
MAG=sum((lambda X:X+X.conj().T)(Move(st,idx,p,N)) for p in P)
ELEC=sum(Zop(st,[k],N)+Zop(st,[k],N).conj().T for k in range(L))

REGIONS=[]
for r in range(1,5):
    for S in itertools.combinations(range(4),r):
        mv=compose([P[i] for i in S])
        REGIONS.append((S, sorted(k for k,_ in mv), Move(st,idx,mv,N)))
print("W-38  15 REGIONS ON THE 3x3 PATCH. Each is a subset of plaquettes; its BOUNDARY is the")
print("      product of them (discrete Stokes). Boundary support printed -- that IS the region's edge.")
for S,sup,_ in REGIONS:
    print(f"      region {str(S):12s} boundary links {sup}")

BATHS={
 "CUT (centre)"     : CUT,
 "RIM (perimeter)"  : PERIM,
 "one corner plaq 0": sorted(k for k,_ in P[0]),
 "one corner plaq 3": sorted(k for k,_ in P[3]),
 "left column"      : [0,4,6,9],
 "single link 0"    : [0],
 "all links"        : list(range(L)),
}

def rates(links,g2=0.01,gam=0.5):
    H=-MAG-g2*ELEC
    Ls=[Zop(st,[k],N) for k in links]
    M=-1j*(np.kron(H,Id)-np.kron(Id,H.T))            # row-major convention (see erratum)
    for X in Ls: M+=gam*(np.kron(X,X.conj())-np.kron(Id,Id))
    w,U=np.linalg.eig(M.conj().T)
    rate=-np.conj(w).real
    U=U/np.linalg.norm(U,axis=0)
    out=[]
    for S,sup,O in REGIONS:
        v=(O/np.linalg.norm(O)).reshape(-1)
        ov=np.abs(U.conj().T@v); ov=ov/max(ov.sum(),1e-30)
        out.append((float((ov*rate).sum()), S, sup))
    return sorted(out)

print("\n  SIEVE RUN WITH THE ENVIRONMENT IN DIFFERENT PLACES.  g2=0.01, gamma=0.5")
print("  For each bath: the surviving (slowest) region, its margin over the runner-up, and whether")
print("  its boundary touches the bath at all.")
print(f"\n  {'bath':>18s} {'winner':>12s} {'rate':>11s} {'runner-up':>12s} {'margin':>8s} {'overlap with bath':>18s}")
print("  "+"-"*90)
for nm,links in BATHS.items():
    rk=rates(links)
    (r1,S1,sup1),(r2,S2,sup2)=rk[0],rk[1]
    touch=len(set(sup1)&set(links))
    print(f"  {nm:>18s} {str(S1):>12s} {r1:11.4e} {str(S2):>12s} {r2/max(r1,1e-30):7.1f}x "
          f"{touch:>7d} of {len(sup1)} links")

print("\n  FULL RANKING for two placements, to show it is not a near-tie.")
for nm in ("CUT (centre)","one corner plaq 0"):
    print(f"\n    bath = {nm}")
    for r,S,sup in rates(BATHS[nm])[:6]:
        print(f"      {r:12.4e}   region {str(S):12s} boundary {sup}   touches bath: "
              f"{len(set(sup)&set(BATHS[nm]))}/{len(sup)}")

print("\n  CONTROL -- bath on ALL links: no region is spared, so no region should stand out.")
rk=rates(BATHS["all links"])
print(f"    slowest {rk[0][0]:.4e} (region {rk[0][1]}), fastest {rk[-1][0]:.4e} (region {rk[-1][1]}), "
      f"spread {rk[-1][0]/max(rk[0][0],1e-30):.2f}x")
