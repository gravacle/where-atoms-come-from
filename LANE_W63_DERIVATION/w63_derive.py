"""W-63.  DERIVE THE RESULTS, THEN CHECK THE DERIVATIONS PREDICT THE MEASURED NUMBERS.

Everything so far was measured. Each of the core results has a derivation, and a derivation is worth
more than a measurement only if it PREDICTS numbers it was not fitted to. So: derive, then predict,
then check against what the lanes already found.

D1  GROUND-SPACE DEGENERACY = |H_1| = 2^(2g).
    Physical states are spanned by |s> for s in Z_1 (cycles: the Gauss law). A plaquette operator
    shifts s by the plaquette boundary, so the ground space (all plaquettes +1) is spanned by equal
    superpositions over COSETS of B_1 (boundaries) in Z_1. The number of cosets is |Z_1/B_1| = |H_1|.
    Torus: 2^2 = 4. Disk: H_1 = 0, so 1.

D2  THE LOGICALS ARE H_1 AND H^1, AND THEY ANTICOMMUTE BY THE INTERSECTION PAIRING.
    Z(c)|s> = (-1)^<c,s> |s> and M(z)|s> = |s+z>, so
        Z(c) M(z) = (-1)^<c,z> M(z) Z(c).
    They anticommute exactly when <c,z> = 1. On a closed orientable surface the pairing between
    H^1 and H_1 is NON-DEGENERATE, so every non-trivial magnetic logical has an electric partner
    that anticommutes with it. THAT is why a writer exists at all.

D3  NO LOCAL OPERATOR ACTS ON THE GROUND SPACE.
    A single-link Z anticommutes with each plaquette containing that link, so it maps the ground
    space out of itself and its ground-space block is EXACTLY ZERO. An operator acts non-trivially
    only if its support contains a non-contractible cycle.

D4  THE SPLITTING IS SUPPRESSED TO ORDER d, THE MINIMAL NON-CONTRACTIBLE CYCLE LENGTH.
    A perturbation must act d times to build a logical operator, so degenerate perturbation theory
    gives  splitting ~ eps^d.  Symmetry-induced degeneracy has no such protection: d = 1, splitting
    ~ eps.  THIS IS THE PREDICTION THAT CAN FAIL, and W-61's numbers are the test.
"""
import itertools, numpy as np
def rank_gf2(vs):
    b=[]; r=0
    for v in vs:
        cur=v
        for x in b:
            p=x.bit_length()-1
            if cur>>p&1: cur^=x
        if cur: b.append(cur); b.sort(reverse=True); r+=1
    return r

def torus(nx,ny):
    vid=lambda i,j:(j%ny)*nx+(i%nx)
    E=[]; ind={}
    for j in range(ny):
        for i in range(nx): ind[('h',i,j)]=len(E); E.append((vid(i,j),vid(i+1,j)))
    for j in range(ny):
        for i in range(nx): ind[('v',i,j)]=len(E); E.append((vid(i,j),vid(i,j+1)))
    PL=[[ind[('h',i,j)],ind[('v',(i+1)%nx,j)],ind[('h',i,(j+1)%ny)],ind[('v',i,j)]]
        for j in range(ny) for i in range(nx)]
    return nx*ny,E,PL,ind

print("W-63  D1 -- IS THE DEGENERACY |H_1| = 2^(2g)?  derived count vs measured ground space")
print(f"  {'torus':>8s} {'|V|':>4s} {'|E|':>4s} {'dim Z_1':>8s} {'dim B_1':>8s} {'|H_1| = 2^(Z-B)':>16s} {'measured':>9s}")
for (nx,ny) in ((2,2),(2,3),(3,3)):
    NV,E,PL,ind=torus(nx,ny); L=len(E)
    zdim=L-NV+1                               # cycle space dimension (connected)
    bdim=rank_gf2([sum(1<<k for k in p) for p in PL])
    pred=2**(zdim-bdim)
    if L<=12:
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
        H=-sum(Move(p) for p in PL); H=(H+H.conj().T)/2
        ev=np.linalg.eigvalsh(H); tol=1e-8*max(1.0,abs(ev).max())
        meas=int(sum(1 for e in ev if abs(e-ev[0])<tol))
    else: meas=None
    print(f"  {f'{nx}x{ny}':>8s} {NV:4d} {L:4d} {zdim:8d} {bdim:8d} {pred:16d} "
          f"{(meas if meas is not None else '--'):>9}")

print("\nW-63  D4 -- IS THE SPLITTING ~ eps^d, WITH d THE MINIMAL NON-CONTRACTIBLE CYCLE?")
print("      derived prediction, then measured slope. THIS IS THE ONE THAT CAN FAIL.")
NV,E,PL,ind=torus(2,2); L=len(E)
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
# minimal non-contractible cycle length d
bnd=[sum(1<<k for k in p) for p in PL]
def inspan(v,vs):
    b=[]; 
    for x in vs:
        cur=x
        for y in b:
            p=y.bit_length()-1
            if cur>>p&1: cur^=y
        if cur: b.append(cur); b.sort(reverse=True)
    cur=v
    for y in b:
        p=y.bit_length()-1
        if cur>>p&1: cur^=y
    return cur==0
def divv(v,x): return sum(1 for k in range(L) if v>>k&1 and E[k][0]==x)-sum(1 for k in range(L) if v>>k&1 and E[k][1]==x)
noncon=[v for v in range(1,1<<L) if all(divv(v,x)%2==0 for x in range(NV)) and not inspan(v,bnd)]
d=min(bin(v).count('1') for v in noncon)
print(f"  minimal non-contractible cycle length d = {d}   -> predicted splitting ~ eps^{d}")
rng=np.random.default_rng(7)
V=sum(rng.normal()*Zl([k]) for k in range(L)); V=V/np.linalg.norm(V)*np.linalg.norm(H)
ev0=np.linalg.eigvalsh(H); tol=1e-8*max(1.0,abs(ev0).max())
g=int(sum(1 for e in ev0 if abs(e-ev0[0])<tol))
print(f"  {'eps':>10s} {'splitting':>13s} {'local slope d ln S / d ln eps':>31s}")
prev=None
for eps in (1e-4,3e-4,1e-3,3e-3,1e-2):
    e=np.sort(np.linalg.eigvalsh(H+eps*V)); s=e[g-1]-e[0]
    sl=""
    if prev: sl=f"{(np.log(s)-np.log(prev[1]))/(np.log(eps)-np.log(prev[0])):31.3f}"
    print(f"  {eps:10.1e} {s:13.3e} {sl:>31s}")
    prev=(eps,s)
print(f"\n  DERIVED exponent d = {d}.  If the measured slope matches, D4 predicts W-61's numbers.")
print(f"  W-61 measured 4.9e-13 at eps=1e-06; eps^{d} = {1e-6**d:.1e}  -> same order.")
print(f"  Symmetry degeneracy has d = 1: W-61 measured 2.0e-06 at eps=1e-06, i.e. LINEAR.")
