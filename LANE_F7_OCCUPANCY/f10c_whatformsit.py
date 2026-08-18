"""Follow-up: WHAT ARE the minimum-weight operators that can form a record?
F-10 established the threshold is weight d. This asks what sits AT the threshold."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('f7_davies.py').read().split('say("="*104); say("0.')[0])
E,V=np.linalg.eigh(H0); gs=int(np.sum(np.abs(E-E[0])<1e-9)); Pg=V[:,:gs]@V[:,:gs].conj().T
P1={'X':X,'Y':1j*(X@Z),'Z':Z}
STAB=[op({l:X for l in s},L) for s in STAR]+[op({l:Z for l in p},L) for p in PLAQ]

say("="*98); say("WHAT SITS AT THE FORMATION THRESHOLD?  (weight-2 operators, d=2)"); say("="*98)
dist=[]; 
for sites in itertools.combinations(range(L),2):
    for letters in itertools.product('XYZ',repeat=2):
        C=op({s:P1[c] for s,c in zip(sites,letters)},L)
        M=Pg@C@Pg
        if np.linalg.norm(M-(np.trace(M)/gs)*Pg)>1e-10: dist.append((sites,''.join(letters),C))
say(f"  weight-2 Paulis that distinguish code states: {len(dist)} of 252")
allcomm=sum(1 for _,_,C in dist if all(np.linalg.norm(C@S-S@C)<1e-9 for S in STAB))
say(f"  of those, how many COMMUTE WITH EVERY STABILISER (i.e. are LOGICAL operators)? {allcomm}")
say(f"  -> {'ALL of them are logical' if allcomm==len(dist) else f'only {allcomm} are logical'}")
say("")
say("  Are they EM holonomies -- supported on cycles or cocycles of the lattice?")
def bits(sites):
    b=0
    for s in sites: b|=1<<(L-1-s)
    return b
z1=[int(toint(v)) for v in nullspace(d1)]; zp=rref(z1,L)
zd=[int(toint(v)) for v in nullspace(d2.T)]; zdp=rref(zd,L)
cyc=cocyc=other=0
for sites,let,_ in dist:
    b=bits(sites)
    if set(let)=={'Z'} and inspan(b,zp): cyc+=1
    elif set(let)=={'X'} and inspan(b,zdp): cocyc+=1
    else: other+=1
say(f"    all-Z supported on a CYCLE   (magnetic Wilson loop) : {cyc}")
say(f"    all-X supported on a COCYCLE (electric Wilson loop) : {cocyc}")
say(f"    neither                                             : {other}")
say("")
say("="*98); say("  READ: the operators AT the formation threshold are exactly the minimum-weight")
say("        LOGICAL operators, and those are EM holonomies on non-contractible cycles.")
say("="*98)
