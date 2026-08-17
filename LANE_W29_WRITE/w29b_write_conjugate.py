# W-29b — WRITE WITH THE CONJUGATE OPERATOR. W-29 used Z on a cut link that does not cross the
# perimeter, so [Z, W] = 0 and nothing could be written. W-27 already identified the conjugate of a
# loop: a cut that PIERCES it an ODD number of times. Use that.
import numpy as np, itertools
V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
L=len(E); NV=len(V2); N=2; CENTER=vid[(1,1)]
CUT=[k for k,(a,b) in enumerate(E) if a==CENTER or b==CENTER]; PERIM=[k for k in range(L) if k not in CUT]
st=[s for s in itertools.product(range(N),repeat=L)
    if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
           -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%N==0 for v in range(NV))]
idx={s:j for j,s in enumerate(st)}; D=len(st)
def hid(i,j): return j*2+i
def vidx(i,j): return 6+j*3+i
PLQ=[[(hid(i,j),+1),(vidx(i+1,j),+1),(hid(i,j+1),-1),(vidx(i,j),-1)] for j in range(2) for i in range(2)]
def loop(moves):
    M=np.zeros((D,D),dtype=complex)
    for j,s in enumerate(st):
        t=list(s)
        for l,sg in moves: t[l]=(t[l]+sg)%N
        M[idx[tuple(t)],j]=1.0
    return M
P=[loop(m) for m in PLQ]; W=loop([(k,+1) for k in PERIM])
Zk=[np.diag([(-1.0)**s[k] for s in st]).astype(complex) for k in range(L)]
Zset=lambda S: np.diag([(-1.0)**sum(s[k] for k in S) for s in st]).astype(complex)

print(f"  find an operator that DOES NOT commute with the record (so it can write to it):")
cands=[]
for r in (1,2,3):
    for S in itertools.combinations(PERIM,r):
        c=np.linalg.norm(Zset(S)@W-W@Zset(S))
        if c>1e-9: cands.append((S,c)); break
    if cands and len(cands[-1][0])==r: pass
for S,c in cands[:3]: print(f"    Z on {str(S):<14} || [Z_S, W] || = {c:.3e}")
WRITER=cands[0][0]
print(f"  WRITER = Z on {WRITER} -- an ODD-piercing cut, the conjugate W-27 identified.\n")

def evolve(rho,H,Ls,gamma,T,steps=3000):
    dt=T/steps
    def d(r):
        x=-1j*(H@r-r@H)
        for Lo in Ls: x+=gamma*(Lo@r@Lo.conj().T-r)
        return x
    for _ in range(steps):
        k1=d(rho); rho=rho+dt*d(rho+0.5*dt*k1); rho=(rho+rho.conj().T)/2; rho/=np.trace(rho).real
    return rho
def proj(Op,val):
    ev,U=np.linalg.eigh(Op.real); sel=U[:,np.abs(ev-val)<1e-9]; return sel@sel.conj().T
Zw=Zset(WRITER)
print(f"  Record starts UNDETERMINED. The WRITER is prepared +1 or -1. Evolve, then read the record.")
print(f"  {'g2':>5}{'gamma':>7}{'T':>6}{'<W> | writer +1':>18}{'<W> | writer -1':>18}{'SEPARATION':>13}")
for g2,gamma,T in ((0.0,0.5,10.0),(1.0,0.5,10.0),(1.0,0.0,10.0),(1.0,0.5,2.0),(2.0,0.5,10.0)):
    H=-sum(P)-g2*sum(Zk)
    out=[]
    for val in (+1.0,-1.0):
        Pw=proj(Zw,val); rho=Pw/np.trace(Pw).real
        rho=evolve(rho,H,P,gamma,T)
        out.append(float(np.real(np.trace(rho@W))))
    print(f"  {g2:>5.1f}{gamma:>7.2f}{T:>6.1f}{out[0]:>18.6f}{out[1]:>18.6f}{abs(out[0]-out[1]):>13.6f}")
print("\n  SEPARATION > 0 means the record ACQUIRED a value carrying what the writer did.")
