"""W-61.  DOES THE GEOMETRY EMIT THE DEGENERACY THAT RECORDS NEED?

W-60: records require exact degeneracy, and SYMMETRY-induced degeneracy is destroyed by a 1e-06
perturbation -- instantly. That is a fatal fragility if symmetry is the only source.

THE PRINCIPAL'S HYPOTHESIS: maybe gravity at the record level EMITS symmetry rather than requiring
it. There is a known kind of degeneracy that behaves exactly that way. TOPOLOGICAL degeneracy comes
from the TOPOLOGY OF THE SPACE, not from any symmetry, and its defining property is robustness to
ARBITRARY LOCAL perturbation.

THE TEST, which W-60 makes decisive because it supplies the contrast:
    same gauge theory, same perturbation strength, two SPACES --
        a TORUS (genus 1)  -> topological ground-state degeneracy expected
        a DISK  (genus 0)  -> none
    perturb with a random LOCAL operator and watch the ground-space splitting.

  if the torus degeneracy SURVIVES what killed the symmetry degeneracy, then the degeneracy records
  need is emitted by the SPACE and not supplied by a symmetry.
  if it dies the same way, topology is no better and the fragility is universal.

CONTROL: at zero perturbation the torus must show the expected degeneracy and the disk must not.
"""
import itertools, numpy as np
rng=np.random.default_rng(53)

def torus(nx,ny):
    """periodic square lattice: vertices (i,j), links right/up from each vertex"""
    NV=nx*ny
    vid=lambda i,j:(j%ny)*nx+(i%nx)
    E=[]; ind={}
    for j in range(ny):
        for i in range(nx):
            ind[('h',i,j)]=len(E); E.append((vid(i,j),vid(i+1,j)))
    for j in range(ny):
        for i in range(nx):
            ind[('v',i,j)]=len(E); E.append((vid(i,j),vid(i,j+1)))
    PL=[[ind[('h',i,j)],ind[('v',(i+1)%nx,j)],ind[('h',i,(j+1)%ny)],ind[('v',i,j)]] for j in range(ny) for i in range(nx)]
    return NV,E,PL
def disk(nx,ny):
    """open square patch"""
    vid={(i,j):j*nx+i for j in range(ny) for i in range(nx)}
    E=[]; ind={}
    for j in range(ny):
        for i in range(nx-1): ind[('h',i,j)]=len(E); E.append((vid[(i,j)],vid[(i+1,j)]))
    for j in range(ny-1):
        for i in range(nx): ind[('v',i,j)]=len(E); E.append((vid[(i,j)],vid[(i,j+1)]))
    PL=[[ind[('h',i,j)],ind[('v',i+1,j)],ind[('h',i,j+1)],ind[('v',i,j)]] for j in range(ny-1) for i in range(nx-1)]
    return nx*ny,E,PL

def build(NV,E,PL):
    L=len(E)
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
    def Zl(k): return np.diag([(-1.0)**s[k] for s in st]).astype(complex)
    H=-sum(Move(p) for p in PL); H=(H+H.conj().T)/2
    return D,H,Move,Zl,L

def gap_profile(name,NV,E,PL,eps_list):
    D,H,Move,Zl,L=build(NV,E,PL)
    ev=np.sort(np.linalg.eigvalsh(H))
    tol=1e-8*max(1.0,abs(ev).max())
    g0=int(sum(1 for e in ev if abs(e-ev[0])<tol))
    print(f"\n  {name}: physical dim {D}, links {L}, plaquettes {len(PL)}")
    print(f"    unperturbed ground-space degeneracy = {g0}")
    print(f"    {'eps':>9s} {'ground splitting':>18s} {'gap to next level':>19s} {'splitting/gap':>14s}")
    for eps in eps_list:
        V=sum(rng.normal()*Zl(k) for k in range(L))         # a random LOCAL perturbation
        V=V/np.linalg.norm(V)*np.linalg.norm(H)
        e=np.sort(np.linalg.eigvalsh(H+eps*V))
        split=e[g0-1]-e[0] if g0>1 else 0.0
        gap=e[g0]-e[g0-1] if g0<len(e) else float('nan')
        print(f"    {eps:9.1e} {split:18.3e} {gap:19.3e} {split/max(gap,1e-300):14.3e}")
    return g0

print("W-61  DOES THE SPACE EMIT THE DEGENERACY?  same theory, same perturbation, two topologies.")
eps=[0.0,1e-6,1e-3,1e-2,1e-1]
gap_profile("TORUS 2x2 (genus 1)",*torus(2,2),eps_list=eps)
gap_profile("DISK 3x3 (genus 0)",*disk(3,3),eps_list=eps)
print()
print("  CONTRAST WITH W-60 -- symmetry-induced degeneracy under the same perturbation size:")
Q,_=np.linalg.qr(rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)))
H0=Q@np.diag(np.array([0,0,1,1,2,2,3,3],float))@Q.conj().T
P=(lambda X:(X+X.conj().T)/2)(rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)))
P=P/np.linalg.norm(P)*np.linalg.norm(H0)
for e_ in (0.0,1e-6,1e-3):
    ev=np.sort(np.linalg.eigvalsh(H0+e_*P))
    print(f"    symmetry degeneracy, eps={e_:.0e}: ground splitting = {ev[1]-ev[0]:.3e}")
