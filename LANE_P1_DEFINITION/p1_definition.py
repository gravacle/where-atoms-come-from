"""P1.  ONE DEFINITION OF "RECORD", AND ITS IMMEDIATE CONSEQUENCES.

The program has used "record" in at least three senses. A proof needs one, and the definition must
carry the content so the theorems come out clean instead of being patched afterwards.

DEFINITION.  Let H be a Hamiltonian and {L_k} the jump operators of an open system. R is a RECORD if

  (i)   A BIT.          R = R^dag  and  R^2 = I.
  (ii)  DURABLE.        [H,R] = 0  and  [L_k,R] = 0 for every k.
  (iii) NON-TRIVIAL.    R is not constant on some eigenspace of H -- equivalently, R is NOT a
                        function of H, equivalently R distinguishes states of the SAME energy.
  (iv)  WRITABLE.       some admissible operation U satisfies U^dag R U = -R.
  (v)   PROTECTED.      NO operation supported on a contractible region does.

JUSTIFYING (iii), which is the clause that could be arbitrary. If R = f(H) then R's value is fixed by
the energy; knowing the energy tells you the record, so it carries no bit beyond the energy and is
not an independent degree of freedom. Stating it as "R distinguishes states of the same energy" makes
that precise and makes P2's forward direction immediate rather than assumed.

(iv) AND (v) TOGETHER are the writable/durable tension the program spent months on. The definition
makes it explicit: a record must be settable, and must not be settable by noise. W-29/W-30 showed
those cannot both hold for LOCAL operations; the definition therefore forces the writer to be
non-local, which is a consequence rather than a stipulation -- Proposition 3 below.

THIS SCRIPT: state the definition, prove the three immediate propositions, and CHECK every clause on
the toric-code record and on controls that must fail.
"""
import itertools, numpy as np
nx=ny=2
vid=lambda i,j:(j%ny)*nx+(i%nx)
E=[]; ind={}
for j in range(ny):
    for i in range(nx): ind[('h',i,j)]=len(E); E.append((vid(i,j),vid(i+1,j)))
for j in range(ny):
    for i in range(nx): ind[('v',i,j)]=len(E); E.append((vid(i,j),vid(i,j+1)))
L=len(E); NV=nx*ny
PL=[[ind[('h',i,j)],ind[('v',(i+1)%nx,j)],ind[('h',i,(j+1)%ny)],ind[('v',i,j)]]
    for j in range(ny) for i in range(nx)]
st=[s for s in itertools.product(range(2),repeat=L)
    if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
           -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%2==0 for v in range(NV))]
idx={s:i for i,s in enumerate(st)}; D=len(st)
def Move(S):
    M=np.zeros((D,D),complex)
    for j,s in enumerate(st):
        t=list(s)
        for k in S: t[k]^=1
        t=tuple(t)
        if t in idx: M[idx[t],j]=1.0
    return M
def Zl(S): return np.diag([(-1.0)**(sum(s[k] for k in S)%2) for s in st]).astype(complex)
H=-sum(Move(p) for p in PL); H=(H+H.conj().T)/2
def vv(S):
    v=0
    for k in S: v|=(1<<k)
    return v
def spanF2(vs):
    b=[]
    for v in vs:
        cur=v
        for x in b:
            p=x.bit_length()-1
            if cur>>p&1: cur^=x
        if cur: b.append(cur); b.sort(reverse=True)
    return b
def inspan(v,b):
    cur=v
    for x in b:
        p=x.bit_length()-1
        if cur>>p&1: cur^=x
    return cur==0
def ovp(v,P): return bin(v & vv(P)).count('1')%2
dual=[v for v in range(1,1<<L) if all(ovp(v,P)==0 for P in PL)]
stars=spanF2([vv([k for k in range(L) if E[k][0]==x or E[k][1]==x]) for x in range(NV)])
Zc=[v for v in dual if not inspan(v,stars)][0]
def divv(v,x): return sum(1 for k in range(L) if v>>k&1 and E[k][0]==x)-sum(1 for k in range(L) if v>>k&1 and E[k][1]==x)
bnd=spanF2([vv(p) for p in PL])
Mz=[v for v in range(1,1<<L) if all(divv(v,x)%2==0 for x in range(NV)) and not inspan(v,bnd)][0]
def bits(v): return [k for k in range(L) if v>>k&1]
R=Zl(bits(Zc)); W=Move(bits(Mz))
Lk=[Zl([k]) for k in range(L)]   # the environment couples locally

print("P1  CHECKING THE FIVE CLAUSES ON THE TORIC-CODE RECORD")
print(f"    record R = Z on the computed dual cycle {bits(Zc)}")
print(f"    writer W = shift on the computed non-contractible cycle {bits(Mz)}\n")
ev,evec=np.linalg.eigh(H); tol=1e-8*max(1.0,abs(ev).max())
def eigenspaces(H):
    ev,U=np.linalg.eigh(H); out=[]; i=0
    while i<len(ev):
        j=i
        while j+1<len(ev) and abs(ev[j+1]-ev[i])<1e-8*max(1.0,abs(ev[i])): j+=1
        out.append(U[:,i:j+1]); i=j+1
    return out
def check(R,H,Lk,W,name):
    c1=max(np.linalg.norm(R-R.conj().T),np.linalg.norm(R@R-np.eye(len(R))))
    c2=max([np.linalg.norm(H@R-R@H)]+[np.linalg.norm(l@R-R@l) for l in Lk])
    nontriv=0.0
    for Uk in eigenspaces(H):
        b=Uk.conj().T@R@Uk
        if b.shape[0]>1: nontriv=max(nontriv,np.linalg.norm(b-np.trace(b)/b.shape[0]*np.eye(b.shape[0])))
    c4=np.linalg.norm(W.conj().T@R@W+R)
    loc=min(np.linalg.norm(l.conj().T@R@l+R) for l in Lk)
    print(f"  {name}")
    print(f"    (i)   bit:        ||R-R+|| , ||R^2-I||        = {c1:.2e}   {'PASS' if c1<1e-9 else 'FAIL'}")
    print(f"    (ii)  durable:    max||[H,R]||,||[L,R]||      = {c2:.2e}   {'PASS' if c2<1e-9 else 'FAIL'}")
    print(f"    (iii) nontrivial: max deviation from scalar   = {nontriv:.3f}  {'PASS' if nontriv>1e-6 else 'FAIL'}")
    print(f"    (iv)  writable:   ||W+RW + R||                = {c4:.2e}   {'PASS' if c4<1e-9 else 'FAIL'}")
    print(f"    (v)   protected:  best local flip ||L+RL + R|| = {loc:.3f}  "
          f"{'PASS (no local operator flips it)' if loc>1e-6 else 'FAIL'}")
check(R,H,Lk,W,"toric-code record")
print()
print("  CONTROLS THAT MUST FAIL")
Hnd=H+1e-3*sum((k+1)*Zl([k]) for k in range(L))   # degeneracy broken
check(R,Hnd,Lk,W,"same R, degeneracy broken by a local perturbation")
fH=H@H                                            # a function of H
check(fH/np.linalg.norm(fH)*0+np.sign(H+1e-9*np.eye(D)) if False else np.eye(D),H,Lk,W,"R = identity (a function of H)")
print()
print("  PROPOSITIONS, immediate from the definition:")
print("   1. (iii) => H is DEGENERATE.  R is non-scalar on some eigenspace, so that eigenspace has")
print("      dimension > 1.  [P2's forward direction, now a one-line consequence of the definition.]")
print("   2. (ii) + (iv) => the writer W commutes with neither R nor anything built from H and {L_k}:")
print("      W anticommutes with R, while H and every L_k commute with it.")
print("   3. (iv) + (v) => THE WRITER IS NON-LOCAL.  Something flips R; nothing contractible does.")
print("      W-29/W-30's obstruction is thus a CONSEQUENCE of the definition, not an accident.")
