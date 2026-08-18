"""W-43b.  ARE GAUGE RECORDS AND METRIC RECORDS THE SAME KIND OF THING?

CONTAMINATION NOTICE: designed after reading Gravacle v337. W-27..W-42 were independent; this is not.

W-43's chirality test failed twice and both failures are instructive:
  (i) the MEASURE was broken -- the antipodal site was weighted -n/2 when it is equally +n/2, so a
      perfectly symmetric distribution scored -0.419, which is exactly the "nothing added" baseline;
  (ii) the CONCEPT was wrong -- a Z_2 flux cannot break time reversal, because 0 and pi are BOTH
      time-reversal invariant. Chirality separates a generic U(1) flux from a potential; it can
      never separate a Z_2 gauge record from a metric one.

THE REAL DIFFERENCE IS LOCALITY, and this program already measured one half of it.
  a GAUGE field is gauge-invariant only around a CLOSED LOOP  -> W-36: fragments learn almost
      nothing about R (I(|F|=1)/I(all) = 0.047); W-37: cut the ring and the reading is exactly 0.
  a POTENTIAL is gauge-invariant POINTWISE                     -> should be locally readable.
So: put BOTH records on the same carrier under the SAME local environment, and measure how
redundantly each is copied. Same machinery as W-36, so the two numbers are comparable.

  METRIC RECORD: a two-level "mass" m at one site, sourcing a potential V*m there. It commutes with
  everything (nothing flips it), so it is durable by construction -- and that is the point: it is
  durable AND local, which is precisely what the gauge record could not be.
"""
import itertools, numpy as np
def expm(A):
    nr=np.linalg.norm(A,np.inf)
    k=max(0,int(np.ceil(np.log2(nr)))+1) if nr>0 else 0
    B=A/(2.0**k); X=np.eye(A.shape[0],dtype=complex); T=X.copy()
    for m in range(1,60):
        T=T@B/m; X=X+T
        if np.linalg.norm(T,np.inf)<1e-18*max(1.0,np.linalg.norm(X,np.inf)): break
    for _ in range(k): X=X@X
    return X
def vn(r):
    ev=np.linalg.eigvalsh((r+r.conj().T)/2); ev=ev[ev>1e-12]
    return float(-(ev*np.log2(ev)).sum())

# ---- gauge carrier: the 3x3 patch, exactly as in W-36 ----
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
def build(N=2):
    st=[s for s in itertools.product(range(N),repeat=L)
        if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
               -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%N==0 for v in range(9))]
    return st,{s:i for i,s in enumerate(st)}
st,idx=build(); DS=len(st)
def Zop(links):
    return np.diag([(-1.0)**(sum(s[k] for k in links)%2) for s in st]).astype(complex)
def Move(mv):
    M=np.zeros((DS,DS),complex)
    for j,s in enumerate(st):
        t=list(s)
        for k,sg in mv: t[k]=(t[k]+sg)%2
        t=tuple(t)
        if t in idx: M[idx[t],j]=1.0
    return M
acc={}
for p in P:
    for k,sg in p: acc[k]=acc.get(k,0)+sg
RG=Move([(k,s) for k,s in acc.items() if s%2])          # rim Wilson loop: the GAUGE record
MAG=sum((lambda X:X+X.conj().T)(Move(p)) for p in P)
sz=np.array([[1,0],[0,-1]],complex); I2=np.eye(2,dtype=complex)
def kron_list(ops):
    o=np.array([[1]],complex)
    for x in ops: o=np.kron(o,x)
    return o

def profile(record, nq=6, kappa=8.0, T=16.0, seed=3, tag=""):
    """record: 'gauge' (rim loop on the patch) or 'metric' (a 2-level mass at one site).
       Environment: nq qubits, each coupling to a LOCAL operator -- Z on a cut link. Identical
       coupling in both cases, so the two numbers are directly comparable."""
    if record=='gauge':
        Dsys=DS; Rop=RG; Hs=-MAG
        locals_=[Zop([k]) for k in CUT]
    else:
        Dsys=2; Rop=sz.copy(); Hs=np.zeros((2,2),complex)   # mass is conserved: nothing flips it
        locals_=[sz.copy()]                                  # the environment sees the mass LOCALLY
    DE=2**nq; 
    H=np.kron(Hs,np.eye(DE,dtype=complex))
    for k in range(nq):
        O=locals_[k%len(locals_)]
        ops=[I2]*nq; ops[k]=sz
        H=H+kappa*np.kron(O,kron_list(ops))
    U=expm(-1j*H*T)
    g=np.random.default_rng(seed)
    Pp=(np.eye(Dsys)+Rop)/2; Pm=(np.eye(Dsys)-Rop)/2
    w=g.normal(size=Dsys)+1j*g.normal(size=Dsys)
    a=Pp@w; b=Pm@w; a/=np.linalg.norm(a); b/=np.linalg.norm(b)
    psiS=(a+b)/np.sqrt(2.0); psiS/=np.linalg.norm(psiS)
    plus=np.ones(2,complex)/np.sqrt(2.0)
    psi=(U@np.kron(psiS,kron_list([plus.reshape(2,1)]*nq).reshape(-1))).reshape((Dsys,)+(2,)*nq)
    def holevo(frag):
        br=[]
        for Proj in (Pp,Pm):
            v=np.tensordot(Proj,psi,axes=([1],[0])); p=float(np.vdot(v,v).real)
            if p<1e-14: br.append((0.0,None)); continue
            v=v/np.sqrt(p)
            keep=[0]+[1+i for i in frag]; tr=[ax for ax in range(1+nq) if ax not in keep]
            M=np.transpose(v,keep+tr).reshape(Dsys*2**len(frag),-1)
            rho=(M@M.conj().T).reshape(Dsys,2**len(frag),Dsys,2**len(frag))
            br.append((p,np.einsum('ijik->jk',rho)))
        avg=sum(p*m for p,m in br if m is not None)
        return vn(avg)-sum(p*vn(m) for p,m in br if m is not None)
    vals=[np.mean([holevo(c) for c in list(itertools.combinations(range(nq),f))[:20]]) for f in range(nq+1)]
    print(f"\n  {tag}  (system dim {Dsys}, {nq} local environment qubits)")
    print("    |F| : "+"  ".join(f"{i}" for i in range(nq+1)))
    print("    I   : "+"  ".join(f"{v:.3f}" for v in vals))
    return vals

print("W-43b  IS A METRIC RECORD THE SAME KIND OF THING AS A GAUGE RECORD?")
print("       Holevo information a fragment of a LOCAL environment holds about the record. Ceiling 1 bit.")
gauge = profile('gauge', tag="GAUGE record: the rim Wilson loop (nonlocal: a product over 8 links)")
metric= profile('metric', tag="METRIC record: a mass at one site (local: gauge-invariant pointwise)")
print()
print(f"  {'record':>8s} {'I(|F|=1)':>10s} {'I(all)':>9s} {'ratio':>7s}  verdict")
print("  "+"-"*74)
for nm,v in (("gauge",gauge),("metric",metric)):
    tot=v[-1]; r=v[1]/tot if tot>1e-9 else 0.0
    print(f"  {nm:>8s} {v[1]:10.4f} {tot:9.4f} {r:7.3f}  "
          f"{'REDUNDANT -> objective for free' if r>0.5 else 'NOT redundant: only the whole environment knows'}")
print()
print("  CONTROL -- kappa = 0 for both. Must be exactly 0 everywhere.")
for nm,rec in (("gauge",'gauge'),("metric",'metric')):
    v=profile(rec,kappa=0.0,tag=f"CONTROL kappa=0, {nm}")
    print(f"    {nm}: max over all fragment sizes = {max(abs(x) for x in v):.2e}")
