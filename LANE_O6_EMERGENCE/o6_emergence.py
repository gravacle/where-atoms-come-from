"""O6.  IS THE 1-FORM SYMMETRY EMERGENT, OR EXACT-MICROSCOPIC?

This carries the whole gravity connection, so it must not rest on a citation.

WHY IT MATTERS. Harlow-Ooguri forbids EXACT GLOBAL symmetries in quantum gravity, and the modern
statement covers HIGHER-FORM ones. If records require an exact microscopic 1-form symmetry, gravity
forbids records -- absurd. If the symmetry is EMERGENT, nothing is forbidden and the claim stands.

THE DISCRIMINATOR, and it is a pair of exponents.
  EXACT symmetry      : ||[W, H]|| = 0 for the full H, at every perturbation strength.
  EMERGENT symmetry   : the microscopic H BREAKS it -- ||[W, H+eps V]|| = O(eps) -- while its
                        LOW-ENERGY consequence, the ground-space degeneracy, survives to O(eps^d).
  THE GAP BETWEEN THOSE TWO EXPONENTS (1 versus d) IS THE EMERGENCE.

CONTROL: a 0-form symmetry must show NO gap -- symmetry broken at order 1 AND degeneracy split at
order 1. If the control shows a gap too, the discriminator is worthless.
"""
import itertools, numpy as np
rng=np.random.default_rng(97)
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
H0=-sum(Move(p) for p in PL); H0=(H0+H0.conj().T)/2
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
dual=[v for v in range(1,1<<L) if all(bin(v & vv(P)).count('1')%2==0 for P in PL)]
stars=spanF2([vv([k for k in range(L) if E[k][0]==x or E[k][1]==x]) for x in range(NV)])
Zc=[v for v in dual if not inspan(v,stars)][0]
W1=Zl([k for k in range(L) if Zc>>k&1])          # the 1-FORM symmetry generator (extended)
V=sum(rng.normal()*Zl([k]) for k in range(L))    # a generic LOCAL perturbation
V=V/np.linalg.norm(V)*np.linalg.norm(H0)

def degen(H):
    ev=np.linalg.eigvalsh(H); tol=1e-8*max(1.0,abs(ev).max())
    return int(sum(1 for e in ev if abs(e-ev[0])<tol)), ev
g0,_=degen(H0)
print("O6  IS THE 1-FORM SYMMETRY EMERGENT?  toric code, 2x2 torus")
print(f"    ground degeneracy at eps=0: {g0}   ||[W_1form, H_0]|| = {np.linalg.norm(W1@H0-H0@W1):.2e}"
      f"   <- EXACT for the unperturbed H\n")
print("  ADD A GENERIC LOCAL PERTURBATION AND WATCH TWO THINGS SEPARATELY:")
print(f"  {'eps':>9s} {'||[W,H+epsV]||':>16s} {'slope':>7s} {'ground splitting':>18s} {'slope':>7s}")
print("  "+"-"*66)
prev=None
for eps in (1e-4,3e-4,1e-3,3e-3,1e-2):
    H=H0+eps*V
    br=np.linalg.norm(W1@H-H@W1)
    e=np.sort(np.linalg.eigvalsh(H)); sp=e[g0-1]-e[0]
    s1=s2=""
    if prev:
        s1=f"{(np.log(br)-np.log(prev[1]))/(np.log(eps)-np.log(prev[0])):7.3f}"
        s2=f"{(np.log(sp)-np.log(prev[2]))/(np.log(eps)-np.log(prev[0])):7.3f}"
    print(f"  {eps:9.1e} {br:16.3e} {s1:>7s} {sp:18.3e} {s2:>7s}")
    prev=(eps,br,sp)
print()
print("  CONTROL -- a 0-FORM symmetry must show NO GAP between the two exponents.")
Q,_=np.linalg.qr(rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)))
Hs=Q@np.diag(np.array([0,0,1,1,2,2,3,3],float))@Q.conj().T
S0=Q@np.diag(np.array([1,-1,1,-1,1,-1,1,-1],float))@Q.conj().T    # the 0-form generator
Vs=(lambda X:(X+X.conj().T)/2)(rng.normal(size=(8,8))+1j*rng.normal(size=(8,8)))
Vs=Vs/np.linalg.norm(Vs)*np.linalg.norm(Hs)
print(f"    ||[S_0form, H_s]|| unperturbed = {np.linalg.norm(S0@Hs-Hs@S0):.2e}")
print(f"  {'eps':>9s} {'||[S,H+epsV]||':>16s} {'slope':>7s} {'ground splitting':>18s} {'slope':>7s}")
print("  "+"-"*66)
prev=None
for eps in (1e-4,3e-4,1e-3,3e-3,1e-2):
    H=Hs+eps*Vs
    br=np.linalg.norm(S0@H-H@S0)
    e=np.sort(np.linalg.eigvalsh(H)); sp=e[1]-e[0]
    s1=s2=""
    if prev:
        s1=f"{(np.log(br)-np.log(prev[1]))/(np.log(eps)-np.log(prev[0])):7.3f}"
        s2=f"{(np.log(sp)-np.log(prev[2]))/(np.log(eps)-np.log(prev[0])):7.3f}"
    print(f"  {eps:9.1e} {br:16.3e} {s1:>7s} {sp:18.3e} {s2:>7s}")
    prev=(eps,br,sp)
print()
print("  READING: if the 1-form case shows symmetry broken at slope ~1 while the degeneracy survives")
print("  to slope ~d>1, the symmetry is BROKEN MICROSCOPICALLY AND RESTORED IN THE INFRARED --")
print("  i.e. EMERGENT, which is exactly what Harlow-Ooguri does NOT forbid.")
print("  If the 0-form control shows the SAME gap, the discriminator is worthless.")
