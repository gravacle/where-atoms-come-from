"""W-62b.  RECORD CREATION -- with the logical operators COMPUTED, not guessed.

w62 guessed the cycle operators and got ||[W_x,H]|| = 22.6. They must come out of the structure:

  MAGNETIC logical (a shift along S): commutes with every plaquette shift automatically, and
      preserves the physical sector iff div(S) = 0, i.e. S is a CYCLE. It acts non-trivially iff S is
      NOT a boundary (not a product of plaquettes). So: cycles modulo boundaries = H_1 of the torus.
  ELECTRIC logical (Z on S): commutes with every plaquette shift iff |S INTERSECT P| is EVEN for
      every plaquette P. It acts non-trivially iff S is not a product of vertex stars.

BOTH ARE COMPUTED BELOW BY GF(2) LINEAR ALGEBRA. Nothing is nominated.

THE THREE TERMS MEET HERE:
  EM       the gauge field carrying the holonomy that IS the record
  GRAVITY  H_1 of the torus is what makes non-trivial logicals exist at all. On a disk H_1 = 0 and
           there are none -- the same theory, no records. THE GENUS DOES THE WORK.
  ALPHA    whether the write completes before decoherence (measured after the structure is fixed)
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
def vecs(S): 
    v=0
    for k in S: v|=(1<<k)
    return v
def spanF2(vs):
    b=[]
    for v in vs:
        cur=v
        for x in b:
            p=x.bit_length()-1
            if cur>>p & 1: cur^=x
        if cur: b.append(cur); b.sort(reverse=True)
    return b
def inspan(v,b):
    cur=v
    for x in b:
        p=x.bit_length()-1
        if cur>>p & 1: cur^=x
    return cur==0
def allsets(L):
    for m in range(1,1<<L): yield m
def bits(v,L): return [k for k in range(L) if v>>k & 1]

# --- cycles: div(S)=0 ; boundaries: span of plaquettes ---
def div(v,vtx):
    return sum(1 for k in bits(v,L) if E[k][0]==vtx)-sum(1 for k in bits(v,L) if E[k][1]==vtx)
cycles=[v for v in allsets(L) if all(div(v,x)%2==0 for x in range(NV))]
bnd=spanF2([vecs(p) for p in PL])
noncontract=[v for v in cycles if not inspan(v,bnd)]
# --- dual: even overlap with every plaquette ; vertex stars ---
def ov(v,P): return bin(v & vecs(P)).count('1')%2
dual=[v for v in allsets(L) if all(ov(v,P)==0 for P in PL)]
stars=spanF2([vecs([k for k in range(L) if E[k][0]==x or E[k][1]==x]) for x in range(NV)])
nontriv_dual=[v for v in dual if not inspan(v,stars)]
print(f"W-62b  TORUS {nx}x{ny}: links {L}, vertices {NV}, plaquettes {len(PL)}")
print(f"  cycles {len(cycles)}   boundaries span dim {len(bnd)}   NON-CONTRACTIBLE cycles {len(noncontract)}")
print(f"  even-overlap sets {len(dual)}   vertex stars span dim {len(stars)}   non-trivial duals {len(nontriv_dual)}")

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
ev,evec=np.linalg.eigh(H); tol=1e-8*max(1.0,abs(ev).max())
g=int(sum(1 for e in ev if abs(e-ev[0])<tol)); G=evec[:,:g]
def blk(O): return G.conj().T@O@G
print(f"\n  ground-space degeneracy {g}, gap {ev[g]-ev[0]:.4f}")

M1=Move(bits(noncontract[0],L))
Zc=[v for v in nontriv_dual]
Z1=Zl(bits(Zc[0],L))
print(f"\n  STEP 1 -- the computed logicals commute with H")
print(f"    magnetic (shift on a non-contractible cycle): ||[M,H]|| = {np.linalg.norm(M1@H-H@M1):.2e}")
print(f"    electric (Z on an even-overlap set)         : ||[Z,H]|| = {np.linalg.norm(Z1@H-H@Z1):.2e}")
bm,bz=blk(M1),blk(Z1)
print(f"\n  STEP 2 -- do they ANTICOMMUTE on the ground space? (that is what makes one write the other)")
print(f"    ||{{M,Z}}|| on ground space = {np.linalg.norm(bm@bz+bz@bm):.2e}   "
      f"||[M,Z]|| = {np.linalg.norm(bm@bz-bz@bm):.2e}")
print(f"    eigenvalues of Z on ground space: {np.round(np.linalg.eigvalsh(bz),6)}")
print(f"\n  STEP 3 -- CONTROL: a CONTRACTIBLE loop must not write")
for p,pl in enumerate(PL[:2]):
    bc=blk(Move(pl))
    print(f"    plaquette {p}: ||[C,Z]|| on ground space = {np.linalg.norm(bc@bz-bz@bc):.2e}")
print(f"\n  STEP 4 -- CONTROL: no LOCAL operator writes")
w=0.0
for k in range(L): w=max(w,np.linalg.norm(blk(Zl([k]))@bz-bz@blk(Zl([k]))))
print(f"    worst single-link Z commutator with the record label: {w:.2e}")
print(f"\n  STEP 5 -- THE GENUS DOES THE WORK: the same theory on a DISK")
def diskcount(nxx,nyy):
    vv={(i,j):j*nxx+i for j in range(nyy) for i in range(nxx)}
    EE=[]; ii={}
    for j in range(nyy):
        for i in range(nxx-1): ii[('h',i,j)]=len(EE); EE.append((vv[(i,j)],vv[(i+1,j)]))
    for j in range(nyy-1):
        for i in range(nxx): ii[('v',i,j)]=len(EE); EE.append((vv[(i,j)],vv[(i,j+1)]))
    P2=[[ii[('h',i,j)],ii[('v',i+1,j)],ii[('h',i,j+1)],ii[('v',i,j)]] for j in range(nyy-1) for i in range(nxx-1)]
    LL=len(EE)
    def dv(v,x): return sum(1 for k in range(LL) if v>>k&1 and EE[k][0]==x)-sum(1 for k in range(LL) if v>>k&1 and EE[k][1]==x)
    cyc=[v for v in range(1,1<<LL) if all(dv(v,x)%2==0 for x in range(nxx*nyy))]
    bb=spanF2([vecs(p) for p in P2])
    return sum(1 for v in cyc if not inspan(v,bb)), len(cyc), len(bb)
nc,ncyc,nb=diskcount(3,3)
print(f"    disk 3x3: cycles {ncyc}, boundary span dim {nb}, NON-CONTRACTIBLE cycles = {nc}")
print(f"    -> on genus 0 there are none. Same theory, no record. THE GENUS IS THE DIFFERENCE.")

print()
print("="*80)
print("  STEP 6 -- ALPHA'S ROLE IN CREATION: the electric term breaks the protection")
ELEC=sum(Zl([k]) for k in range(L))
print(f"  {'g^2':>7s} {'ground degeneracy':>18s} {'gap':>9s} {'record label eigenvalues':>28s} {'|<M,Z> anticomm|':>17s}")
for g2 in (0.0,1e-4,1e-2,0.1,0.5):
    Hg=-sum(Move(p) for p in PL)-g2*ELEC; Hg=(Hg+Hg.conj().T)/2
    e2,v2=np.linalg.eigh(Hg); t2=1e-8*max(1.0,abs(e2).max())
    gg=int(sum(1 for x in e2 if abs(x-e2[0])<t2)); GG=v2[:,:max(gg,1)]
    bz2=GG.conj().T@Z1@GG; bm2=GG.conj().T@M1@GG
    lab=np.round(np.linalg.eigvalsh((bz2+bz2.conj().T)/2),4)
    anti=np.linalg.norm(bm2@bz2+bz2@bm2)
    print(f"  {g2:7.4f} {gg:18d} {e2[gg]-e2[0] if gg<len(e2) else float('nan'):9.4f} "
          f"{str(lab):>28s} {anti:17.2e}")
print()
print("  READING: at g^2 = 0 the ground space is 4-fold, the record label has clean +-1 eigenvalues,")
print("  and the non-contractible transport anticommutes with it -- it WRITES. As alpha rises the")
print("  degeneracy is broken and the record's home shrinks. ALPHA IS WHAT COSTS THE RECORD.")
