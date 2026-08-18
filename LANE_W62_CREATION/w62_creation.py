"""W-62.  RECORD CREATION: THE STEP WHERE ALL THREE TERMS MEET.

THE AIM: how EM, gravity and alpha COMBINE to create a record. Every prior lane measured one of them.
This is the operation in which all three appear at once.

  EM       the gauge field transported through, and the holonomy that IS the record
  GRAVITY  the NON-CONTRACTIBLE cycle. Transport around a contractible loop must change nothing;
           around a non-contractible one it must change the sector. THE GENUS IS WHAT MAKES THE
           OPERATION WRITE ANYTHING AT ALL.
  ALPHA    the coupling deciding whether the transport completes before decoherence destroys it

AND IT CLOSES THE PROGRAM'S OLDEST OBSTRUCTION. W-29/W-30: writable and durable are conjugate --
anything able to write the record can also destroy it. That was a search for a LOCAL writer for a
record only a NON-LOCAL operation can set. T1's no-go applies to operations COMMUTING with R;
transport around a cycle does not commute with it, and is not local, so noise cannot mimic it.

WHAT MUST BE SHOWN, in order:
 1. the torus ground space is 4-fold degenerate and the sectors are labelled by two Wilson loops
    on the two non-contractible cycles                                     (the record's labels)
 2. transport around a NON-CONTRACTIBLE cycle CHANGES the sector           (the write)
 3. transport around a CONTRACTIBLE loop does NOT                          (the control that decides it)
 4. no LOCAL operator changes the sector                                   (durability, retained)
Steps 3 and 4 are what make step 2 mean something.
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
def Move(links):
    M=np.zeros((D,D),complex)
    for j,s in enumerate(st):
        t=list(s)
        for k in links: t[k]^=1
        t=tuple(t)
        if t in idx: M[idx[t],j]=1.0
    return M
def Zl(links): return np.diag([(-1.0)**(sum(s[k] for k in links)%2) for s in st]).astype(complex)
H=-sum(Move(p) for p in PL); H=(H+H.conj().T)/2
print(f"W-62  TORUS {nx}x{ny}: {NV} vertices, {L} links, {len(PL)} plaquettes, physical dim {D}")

ev,evec=np.linalg.eigh(H)
tol=1e-8*max(1.0,abs(ev).max())
g=int(sum(1 for e in ev if abs(e-ev[0])<tol)); G=evec[:,:g]
print(f"\n  STEP 1 -- ground space")
print(f"    degeneracy = {g}   gap to next level = {ev[g]-ev[0]:.4f}")
# the two non-contractible cycles, as electric (Z) loops on the dual -- these label the sectors
CYC_X=[ind[('h',i,0)] for i in range(nx)]          # wraps in x
CYC_Y=[ind[('v',0,j)] for j in range(ny)]          # wraps in y
WX=Zl(CYC_X); WY=Zl(CYC_Y)
def blk(O): return G.conj().T@O@G
bx,by=blk(WX),blk(WY)
print(f"    ||[W_x, H]|| = {np.linalg.norm(WX@H-H@WX):.2e}   ||[W_y, H]|| = {np.linalg.norm(WY@H-H@WY):.2e}")
print(f"    on the ground space: eigenvalues of W_x = {np.round(np.linalg.eigvalsh(bx),6)}")
print(f"                         eigenvalues of W_y = {np.round(np.linalg.eigvalsh(by),6)}")
print(f"    ||[W_x,W_y]|| on ground space = {np.linalg.norm(bx@by-by@bx):.2e}  -> they label sectors jointly")

print(f"\n  STEP 2 -- transport around a NON-CONTRACTIBLE cycle (the write)")
# the conjugate operation: a magnetic string wrapping the OTHER cycle
STR_X=[ind[('v',i,0)] for i in range(nx)]          # magnetic string wrapping x
STR_Y=[ind[('h',0,j)] for j in range(ny)]          # magnetic string wrapping y
TX,TY=Move(STR_X),Move(STR_Y)
for nm,T in (("transport wrapping x",TX),("transport wrapping y",TY)):
    b=blk(T)
    print(f"    {nm}: ||[T,H]|| = {np.linalg.norm(T@H-H@T):.2e}   "
          f"anticommutes with W_x? {np.linalg.norm(bx@blk(T)+blk(T)@bx):.2e}   "
          f"with W_y? {np.linalg.norm(by@blk(T)+blk(T)@by):.2e}")

print(f"\n  STEP 3 -- CONTROL: transport around a CONTRACTIBLE loop must change nothing")
for p,pl in enumerate(PL):
    C=Move(pl); bc=blk(C)
    print(f"    plaquette {p} (contractible): commutes with W_x? "
          f"{np.linalg.norm(bx@bc-bc@bx):.2e}   with W_y? {np.linalg.norm(by@bc-bc@by):.2e}")

print(f"\n  STEP 4 -- CONTROL: no LOCAL operator changes the sector")
worst_z=worst_x=0.0
for k in range(L):
    z=blk(Zl([k])); x=blk(Move([k]))
    worst_z=max(worst_z,np.linalg.norm(bx@z-z@bx),np.linalg.norm(by@z-z@by))
    worst_x=max(worst_x,np.linalg.norm(bx@x-x@bx),np.linalg.norm(by@x-x@by))
print(f"    single-link Z, worst commutator with the sector labels: {worst_z:.2e}")
print(f"    single-link shift, worst commutator with the sector labels: {worst_x:.2e}")
print(f"    (a single-link shift leaves the physical sector, so its block is not a valid operation)")
