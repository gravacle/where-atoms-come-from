"""W-36b.  WHY IS THERE NO REDUNDANCY?  THE OBSTRUCTION IS LOCALITY.

W-36: with the environment coupled to the CUT, no fragment learns anything about R -- information
rises only near the full environment. I(|F|=1)/I(all) = 0.047. No redundancy, so no objectivity.

HYPOTHESIS. Redundancy needs the environment to couple DIRECTLY to the record observable: many
modes each scattering off the same pointer variable, the way many photons each carry the position
of a dust grain. Coupling to R itself is NOT smuggling -- L = R is symmetric between the +1 and -1
sectors, so it cannot decide which one wins, and W-34 already selected R independently by the
predictability sieve. What it IS, is a QND measurement of R.

TEST. Same construction, three couplings, everything else identical:
  A  env couples to Z on cut links   (local, gauge-invariant, [L,R]=0)          -- W-36's case
  B  env couples to R itself         (nonlocal: a product over all 8 rim links) -- QND on R
  C  env couples to a SINGLE plaquette W_p (local, gauge-invariant)             -- local magnetic
If B gives redundancy and A and C do not, the obstruction is that the protected observable is
NONLOCAL and a local environment cannot make copies of it.
"""
import itertools, numpy as np
exec(open('w36_darwinism.py').read().split('print("W-36  IS THE RECORD')[0])

def run2(nq, sysops, kappa, T, seed=3):
    DE=2**nq
    H=np.kron(-MAG,np.eye(DE,dtype=complex))
    for k in range(nq):
        O=sysops[k%len(sysops)]
        ops=[np.eye(2,dtype=complex)]*nq; ops[k]=sz
        H=H+kappa*np.kron(O,kron_list(ops))
    U=expm(-1j*H*T)
    g=np.random.default_rng(seed)
    w=g.normal(size=DS)+1j*g.normal(size=DS)
    a=Pp@w; b=Pm@w; a/=np.linalg.norm(a); b/=np.linalg.norm(b)
    psiS=(a+b)/np.sqrt(2.0); psiS/=np.linalg.norm(psiS)
    plus=np.ones(2,complex)/np.sqrt(2.0)
    psiE=kron_list([plus.reshape(2,1)]*nq).reshape(-1)
    return (U@np.kron(psiS,psiE)).reshape((DS,)+(2,)*nq)

def prof(nq,sysops,kappa,T,tag):
    psiT=run2(nq,sysops,kappa,T)
    vals=[]
    for f in range(nq+1):
        combos=list(itertools.combinations(range(nq),f))[:20]
        vals.append(float(np.mean([holevo(psiT,nq,c)[0] for c in combos])))
    print(f"\n  {tag}")
    print("    |F| : "+"  ".join(f"{i}" for i in range(nq+1)))
    print("    I   : "+"  ".join(f"{v:.3f}" for v in vals))
    return vals

Wp=Move(st,idx,P[0],N); Wp=(Wp+Wp.conj().T)/2
Zc=[Zop(st,[k],N) for k in CUT]
print("W-36b  WHAT MUST THE ENVIRONMENT TOUCH FOR THE RECORD TO BE COPIED?")
print(f"       R is a product over {len(PERIM)} rim links: NONLOCAL. Z_cut and W_p are LOCAL.")
print(f"       support check -- ||[R, Z_cut]|| = {max(np.linalg.norm(R@z-z@R) for z in Zc):.1e}"
      f"   ||[R, W_p]|| = {np.linalg.norm(R@Wp-Wp@R):.1e}")
K,T=8.0,16.0
a=prof(6,Zc,K,T,      "A  env couples to Z on cut links  (LOCAL, [L,R]=0)")
b=prof(6,[R],K,T,     "B  env couples to R itself        (NONLOCAL, QND on the record)")
c=prof(6,[Wp],K,T,    "C  env couples to one plaquette   (LOCAL magnetic)")
z=prof(6,[R],0.0,T,   "CONTROL kappa=0                   (must be 0 everywhere)")

print()
print(f"  {'case':>36s} {'I(|F|=1)':>9s} {'I(all)':>8s} {'ratio':>7s}  verdict")
print("  "+"-"*88)
for nm,v in [("A local, cut links",a),("B nonlocal, the record itself",b),("C local, one plaquette",c),("control kappa=0",z)]:
    tot=v[-1]
    if tot<1e-9:
        print(f"  {nm:>36s} {'--':>9s} {tot:8.4f} {'--':>7s}  no information anywhere"); continue
    r=v[1]/tot
    vd=("REDUNDANT -> objective: one fragment already knows" if r>0.5 else
        "partial" if r>0.2 else "NOT redundant: only the whole environment knows")
    print(f"  {nm:>36s} {v[1]:9.4f} {tot:8.4f} {r:7.3f}  {vd}")
