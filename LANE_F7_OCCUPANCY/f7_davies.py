"""F-7, first test.  Given a carrier whose systole PERMITS a record, does energy-based
dynamics OCCUPY a definite record, or only a mixture?

REGISTERED PREDICTION, made before this run (RECORD_FORMATION_V001; ledger F-6):
  a Davies bath must reach an EQUAL MIXTURE of all 2^2g record states -- no record --
  because clause (iii) makes the record distinguish states of the SAME energy, while
  detailed balance fixes populations by ENERGY alone.

Carrier: toric code, 2x2 torus, L=8 qubits, dim 256, ground space 4.
All Davies machinery is built in the ENERGY EIGENBASIS, where the spectral projections
are index masks and A(w) needs no matrix products."""
import sys, numpy as np
def say(*a): print(*a); sys.stdout.flush()

I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def op(p,L):
    M=np.array([[1]],dtype=complex)
    for l in range(L): M=np.kron(M,p.get(l,I2))
    return M

nx=ny=2; L=2*nx*ny; ind={}; k=0
for j in range(ny):
    for i in range(nx): ind[('h',i,j)]=k; k+=1; ind[('v',i,j)]=k; k+=1
STAR=[[ind[('h',i,j)],ind[('h',(i-1)%nx,j)],ind[('v',i,j)],ind[('v',i,(j-1)%ny)]] for j in range(ny) for i in range(nx)]
PLAQ=[[ind[('h',i,j)],ind[('v',(i+1)%nx,j)],ind[('h',i,(j+1)%ny)],ind[('v',i,j)]] for j in range(ny) for i in range(nx)]
H0=-sum(op({l:X for l in st},L) for st in STAR)-sum(op({l:Z for l in p},L) for p in PLAQ)

# ---- LOGICALS COMPUTED FROM STRUCTURE, NOT NOMINATED (W-62 convention) ----
def rref(rows,n):
    piv={}; 
    for r in rows:
        c=r
        for j in sorted(piv,reverse=True):
            if (c>>j)&1: c^=piv[j]
        if c: piv[c.bit_length()-1]=c
    return piv
def inspan(v,piv):
    for j in sorted(piv,reverse=True):
        if (v>>j)&1: v^=piv[j]
    return v==0
def tovec(bits,n): return [(bits>>(n-1-i))&1 for i in range(n)]

EDGES=[None]*L
for j in range(ny):
    for i in range(nx):
        EDGES[ind[('h',i,j)]]=(j*nx+i, j*nx+(i+1)%nx)
        EDGES[ind[('v',i,j)]]=(j*nx+i, ((j+1)%ny)*nx+i)
NV=nx*ny
d1=np.zeros((NV,L),dtype=np.int8)
for k_,(a,b) in enumerate(EDGES): d1[a,k_]^=1; d1[b,k_]^=1
d2=np.zeros((L,len(PLAQ)),dtype=np.int8)
for k_,pl in enumerate(PLAQ):
    for e in pl: d2[e,k_]^=1
assert not ((d1@d2)%2).any(), "d1 d2 != 0"

def nullspace(M):
    M=M.copy()%2; rows,cols=M.shape; pc=[]; r=0
    for c in range(cols):
        pv=next((i for i in range(r,rows) if M[i,c]),None)
        if pv is None: continue
        M[[r,pv]]=M[[pv,r]]
        for i in range(rows):
            if i!=r and M[i,c]: M[i]^=M[r]
        pc.append(c); r+=1
    out=[]
    for fc in [c for c in range(cols) if c not in pc]:
        v=np.zeros(cols,dtype=np.int8); v[fc]=1
        for i,pcc in enumerate(pc): v[pcc]=M[i,fc]
        out.append(v)
    return out
def toint(v): return int(''.join(map(str,v)),2)

Z1=[int(toint(v)) for v in nullspace(d1)]                 # cycles      -> Z-logicals live here
B1=rref([int(toint(d2[:,c])) for c in range(d2.shape[1])],L)
Zcand=[c for c in [x for m in range(1,1<<len(Z1)) for x in [int(np.bitwise_xor.reduce([Z1[i] for i in range(len(Z1)) if (m>>i)&1]))]] if not inspan(c,B1)]
Zc1=min(Zcand,key=lambda c:bin(c).count('1'))
Zsp=rref([Zc1]+[int(toint(d2[:,c])) for c in range(d2.shape[1])],L)
Zc2=min([c for c in Zcand if not inspan(c,Zsp)],key=lambda c:bin(c).count('1'))

Zdual=[int(toint(v)) for v in nullspace(d2.T)]            # cocycles    -> X-logicals live here
Bd=rref([int(toint(d1[r,:])) for r in range(d1.shape[0])],L)
Xcand=[c for c in [x for m in range(1,1<<len(Zdual)) for x in [int(np.bitwise_xor.reduce([Zdual[i] for i in range(len(Zdual)) if (m>>i)&1]))]] if not inspan(c,Bd)]
def ov(a,b): return bin(a&b).count('1')%2
Xc1=min([c for c in Xcand if ov(c,Zc1)==1 and ov(c,Zc2)==0],key=lambda c:bin(c).count('1'))

def zop(c): return op({l:Z for l in range(L) if (c>>(L-1-l))&1},L)
def xop(c): return op({l:X for l in range(L) if (c>>(L-1-l))&1},L)
Zbar,Zbar2,Xbar=zop(Zc1),zop(Zc2),xop(Xc1)
LOGINFO=(bin(Zc1).count('1'),bin(Zc2).count('1'),bin(Xc1).count('1'),len(Z1),len(B1),len(Z1)-len(B1))

COUP=[op({l:X},L) for l in range(L)]+[op({l:Z},L) for l in range(L)]

def kms(w,beta,gam=1.0): return gam*np.exp(beta*w/2)/(2*np.cosh(beta*w/2))

def davies_eig(Hm,beta,tol=1e-8):
    """returns (E, V, jumps) with jumps already in the ENERGY EIGENBASIS."""
    E,V=np.linalg.eigh(Hm)
    dE=E[None,:]-E[:,None]                       # dE[a,b] = E_b - E_a  (emission when >0)
    ws=np.unique(np.round(dE,6))
    jumps=[]
    for C in COUP:
        Cp=V.conj().T@C@V
        for w in ws:
            g=kms(w,beta)
            if g<1e-14: continue
            M=np.where(np.abs(dE-w)<tol,Cp,0.0)
            if np.abs(M).max()<1e-12: continue
            jumps.append(np.sqrt(g)*M)
    return E,V,jumps

def lind(rho,E,jumps):
    d=-1j*(E[:,None]*rho-rho*E[None,:])
    for J in jumps:
        Jd=J.conj().T
        d+=J@rho@Jd-0.5*(Jd@J@rho+rho@Jd@J)
    return d

say("="*104); say("0.  CARRIER AND SELF-CHECKS"); say("="*104)
E0,V0=np.linalg.eigh(H0); gs=int(np.sum(np.abs(E0-E0[0])<1e-9))
say(f"  dim Z_1={LOGINFO[3]}  dim B_1={LOGINFO[4]}  dim H_1={LOGINFO[5]}   logical weights: Zbar={LOGINFO[0]}, Zbar2={LOGINFO[1]}, Xbar={LOGINFO[2]}")
say(f"  L={L}  dim={2**L}  E_min={E0[0]:.4f}   GROUND DEGENERACY = {gs}   {'PASS (expect 4 = 2^2g)' if gs==4 else 'FAIL'}")
for nm,O in (("Zbar",Zbar),("Zbar2",Zbar2),("Xbar",Xbar)):
    n=np.linalg.norm(O@H0-H0@O); say(f"  ||[{nm},H0]||  = {n:.3e}   {'PASS' if n<1e-9 else 'FAIL'}")
ac=np.linalg.norm(Zbar@Xbar+Xbar@Zbar)
say(f"  ||{{Zbar,Xbar}}|| = {ac:.3e}   {'PASS -- Xbar WRITES the record' if ac<1e-9 else 'FAIL'}")
Pg=V0[:,:gs]@V0[:,:gs].conj().T
say(f"  Zbar on the ground space: eigenvalues {np.round(np.linalg.eigvalsh(Pg@Zbar@Pg),6)[-4:]}  -> a genuine +-1 label")

say(""); say("="*104); say("1.  POSITIVE CONTROL -- can this bath thermalise at all?  (the W-49 trap)"); say("="*104)
for beta in (0.5,2.0):
    Hc=np.diag([0.,1.]); Ec=np.array([0.,1.]); dEc=Ec[None,:]-Ec[:,None]
    js=[]
    for C in (X,Z):
        for w in np.unique(np.round(dEc,6)):
            g=kms(w,beta); M=np.where(np.abs(dEc-w)<1e-8,C,0.0)
            if np.abs(M).max()>1e-12 and g>1e-14: js.append(np.sqrt(g)*M)
    S=np.zeros((4,4),dtype=complex)
    for a in range(2):
        for b in range(2):
            B=np.zeros((2,2),dtype=complex); B[a,b]=1
            S[:,2*a+b]=lind(B,Ec,js).reshape(4)
    w_,v_=np.linalg.eig(S); ss=v_[:,np.argmin(np.abs(w_))].reshape(2,2)
    ss=ss/np.trace(ss); pops=np.real(np.diag(ss))
    bz=np.exp(-beta*Ec); bz/=bz.sum()
    say(f"  beta={beta}:  steady populations {np.round(pops,6)}   Boltzmann {np.round(bz,6)}   "
        f"{'PASS' if np.allclose(pops,bz,atol=1e-6) else 'FAIL'}")

say(""); say("="*104); say("2.  THE TEST -- the Davies steady state on the code space"); say("="*104)
say(f"  {'beta':>6}{'||L(rho_Gibbs)||':>19}{'<Zbar>':>11}{'<Zbar2>':>11}{'<Xbar>':>11}{'purity on code':>17}{'S_code/ln4':>13}")
store={}
for beta in (0.5,1.0,2.0,5.0):
    E,V,J=davies_eig(H0,beta)
    p=np.exp(-beta*(E-E.min())); p/=p.sum(); rg=np.diag(p).astype(complex)
    res=np.linalg.norm(lind(rg,E,J))
    Zb=V.conj().T@Zbar@V; Zb2=V.conj().T@Zbar2@V; Xb=V.conj().T@Xbar@V; Pge=V.conj().T@Pg@V
    rc=Pge@rg@Pge; tr=np.real(np.trace(rc)); rcn=rc/tr
    ev=np.linalg.eigvalsh(rcn); ev=ev[ev>1e-12]; S=-np.sum(ev*np.log(ev))
    say(f"  {beta:>6}{res:>19.3e}{np.real(np.trace(rg@Zb)):>11.6f}{np.real(np.trace(rg@Zb2)):>11.6f}"
        f"{np.real(np.trace(rg@Xb)):>11.6f}{np.real(np.trace(rcn@rcn)):>17.6f}{S/np.log(4):>13.6f}")
    store[beta]=(E,V,J)

say(""); say("="*104); say("3.  DOES A DEFINITE RECORD DECAY?  d<Zbar>/dt from a pure Zbar=+1 state"); say("="*104)
beta=2.0; E,V,J=store[beta]
Pge=V.conj().T@Pg@V; Zb=V.conj().T@Zbar@V
M=Pge@Zb@Pge; w_,u_=np.linalg.eigh(M); v=u_[:,-1]; v/=np.linalg.norm(v)
rho=np.outer(v,v.conj())
z0=np.real(np.trace(rho@Zb)); dz=np.real(np.trace(lind(rho,E,J)@Zb))
say(f"  initial <Zbar> = {z0:+.6f}   (a DEFINITE record)")
say(f"  d<Zbar>/dt     = {dz:+.6f}   -> {'DECAYS toward the mixture' if dz<-1e-6 else 'does not decay'}")
say(f"  steady <Zbar>  = {np.real(np.trace(np.diag(np.exp(-beta*(E-E.min()))/np.sum(np.exp(-beta*(E-E.min()))))@Zb)):+.6f}")

say(""); say("="*104); say("4.  THE CONTROL THAT MATTERS -- can a LOCAL perturbation let the bath select?"); say("="*104)
say(f"  {'eps':>7}{'ground splitting dE':>22}{'<Zbar>_Gibbs (beta=5)':>24}{'still protected?':>26}")
pert=sum(op({l:X},L) for l in range(L))
for eps in (0.0,0.02,0.05,0.10,0.20):
    Hm=H0+eps*pert
    Em,Vm=np.linalg.eigh(Hm); dEs=Em[3]-Em[0]
    p=np.exp(-5.0*(Em-Em.min())); p/=p.sum(); rg=np.diag(p).astype(complex)
    Zb=Vm.conj().T@Zbar@Vm
    say(f"  {eps:>7.3f}{dEs:>22.3e}{np.real(np.trace(rg@Zb)):>24.6f}"
        f"{('YES -- exactly degenerate' if dEs<1e-9 else 'NO -- degeneracy lifted'):>26}")

say(""); say("="*104); say("5.  WHY IS PART 4 ZERO?  a symmetry claim, and the positive control that tests it"); say("="*104)
pX=sum(op({l:X},L) for l in range(L)); pZ=sum(op({l:Z},L) for l in range(L))
say(f"  ||[sum X_l, Xbar]|| = {np.linalg.norm(pX@Xbar-Xbar@pX):.3e}   <- if 0, H+eps*sumX commutes with the WRITER Xbar")
say(f"  ||[sum Z_l, Xbar]|| = {np.linalg.norm(pZ@Xbar-Xbar@pZ):.3e}   <- if nonzero, that symmetry is broken")
say(f"  ||[sum Z_l, Zbar]|| = {np.linalg.norm(pZ@Zbar-Zbar@pZ):.3e}   <- if 0, Zbar is still a good record label")
say("")
say("  CLAIM: Xbar anticommutes with Zbar, so any rho commuting with Xbar has <Zbar> = 0 EXACTLY.")
say("  So a perturbation that PRESERVES the writer's symmetry can never bias the record,")
say("  however much it lifts the degeneracy. Test: use a perturbation that BREAKS it.")
say("")
say(f"  {'eps':>7}{'perturbation':>16}{'[pert,Xbar]':>14}{'splitting dE':>15}{'<Zbar>_Gibbs b=5':>19}{'selects?':>12}")
for nm,P in (("sum X_l",pX),("sum Z_l",pZ)):
    cX=np.linalg.norm(P@Xbar-Xbar@P)
    for eps in (0.05,0.20):
        Hm=H0+eps*P
        Em,Vm=np.linalg.eigh(Hm); dEs=Em[3]-Em[0]
        p=np.exp(-5.0*(Em-Em.min())); p/=p.sum()
        zb=np.real(np.trace(np.diag(p).astype(complex)@(Vm.conj().T@Zbar@Vm)))
        say(f"  {eps:>7.3f}{nm:>16}{cX:>14.2e}{dEs:>15.3e}{zb:>19.6f}{('YES' if abs(zb)>1e-6 else 'no'):>12}")
say("")
say("  INTERPRETATION TEST -- is the zero a symmetry, or a broken measurement?")
say("  If 'sum Z_l' also gives 0, the measurement cannot see selection and every zero above is void.")

say(""); say("="*104); say("6.  THE EXCHANGE RATE -- selection bought with protection"); say("="*104)
say(f"  {'eps':>8}{'splitting dE':>15}{'<Zbar>_Gibbs':>15}{'tanh(b dE/2)':>15}{'<Zbar>/dE':>12}   (beta=5)")
beta=5.0; rows=[]
for eps in (0.01,0.02,0.05,0.10,0.15,0.20,0.30):
    Hm=H0+eps*pZ
    Em,Vm=np.linalg.eigh(Hm); dEs=Em[3]-Em[0]
    p=np.exp(-beta*(Em-Em.min())); p/=p.sum()
    zb=np.real(np.trace(np.diag(p).astype(complex)@(Vm.conj().T@Zbar@Vm)))
    rows.append((eps,dEs,zb))
    say(f"  {eps:>8.3f}{dEs:>15.4e}{zb:>15.6f}{np.tanh(beta*dEs/2):>15.6f}{(zb/dEs if dEs>0 else 0):>12.4f}")
import numpy as _np
e=_np.array([r[0] for r in rows]); dd=_np.array([r[1] for r in rows]); zz=_np.array([r[2] for r in rows])
pe=_np.polyfit(_np.log(e),_np.log(dd),1)[0]; pz=_np.polyfit(_np.log(e),_np.log(_np.abs(zz)),1)[0]
say("")
say(f"  log-log slope  d(dE)/d(eps)    = {pe:.4f}    (Thm D: splitting ~ eps^d, derived d=2)")
say(f"  log-log slope  d(<Zbar>)/d(eps) = {pz:.4f}")
say(f"  -> SELECTION AND SPLITTING SHARE THE SAME EXPONENT: {'YES' if abs(pe-pz)<0.15 else 'NO'}")
say("")
say("  READ: the bias in the record and the loss of protection are the SAME quantity.")
say("  There is no eps at which one appears without the other.")
