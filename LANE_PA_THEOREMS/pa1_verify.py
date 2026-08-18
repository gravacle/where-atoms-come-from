"""PHASE A / T1.  DURABILITY IMPLIES UNWRITABLE -- hypotheses tested at their edges.

THE PROOF (full statement and derivation in THEOREMS_V001.md). For the GKSL generator
    L[rho] = -i[H,rho] + sum_k ( L_k rho L_k^dag - (1/2){L_k^dag L_k, rho} )
    d<R>/dt = Tr(R L[rho])
             = -i Tr([R,H] rho)                                     -> 0 if [H,R]=0
               + sum_k [ Tr(L_k^dag R L_k rho) - (1/2) Tr(R {L_k^dag L_k, rho}) ]
If [L_k,R]=0 then L_k^dag R L_k = L_k^dag L_k R, and if ALSO [L_k^dag,R]=0 then
[L_k^dag L_k, R]=0 and the two remaining terms cancel exactly.

THE HYPOTHESIS THAT DOES THE WORK IS THEREFORE [L^dag,R]=0, NOT UNITARITY OF R.
And [L,R]=0 implies [L^dag,R]=0 whenever R is NORMAL: a normal R has orthogonal eigenspaces, any L
commuting with R preserves each of them, so L is block diagonal in that basis and so is L^dag.
Unitary R is a special case. Non-normal R is where it must fail.

CONSTRUCTION NOTE. The first version built commuting operators from Lagrange spectral projectors,
which accumulates error and produced residuals of 1e-9 that were mistaken for a signal. Operators are
now built by EXACT simultaneous diagonalisation, so any residual is real.
"""
import numpy as np
rng=np.random.default_rng(4)
def gksl(rho,H,Ls):
    out=-1j*(H@rho-rho@H)
    for L in Ls: out=out+L@rho@L.conj().T-0.5*(L.conj().T@L@rho+rho@L.conj().T@L)
    return out
def rand_rho(D):
    A=rng.normal(size=(D,D))+1j*rng.normal(size=(D,D)); r=A@A.conj().T
    return r/np.trace(r).real
def U(D):
    Q,_=np.linalg.qr(rng.normal(size=(D,D))+1j*rng.normal(size=(D,D))); return Q
def build(D,kind):
    """R of the requested kind, plus H and L_k built to commute with it EXACTLY."""
    Q=U(D)
    if kind=='unitary':  r=np.exp(1j*rng.uniform(0,2*np.pi,D))
    elif kind=='normal': r=np.exp(1j*rng.uniform(0,2*np.pi,D))*rng.uniform(0.4,2.0,D)
    R=Q@np.diag(r)@Q.conj().T
    H=Q@np.diag(rng.normal(size=D))@Q.conj().T                       # Hermitian, [H,R]=0 exactly
    Ls=[Q@np.diag(rng.normal(size=D)+1j*rng.normal(size=D))@Q.conj().T for _ in range(3)]
    return R,(H+H.conj().T)/2,Ls
def worst(R,H,Ls,n=40):
    D=R.shape[0]; return max(abs(np.trace(R@gksl(rand_rho(D),H,Ls))) for _ in range(n))
def cc(R,H,Ls):
    return (np.linalg.norm(H@R-R@H),
            max(np.linalg.norm(L@R-R@L) for L in Ls),
            max(np.linalg.norm(L.conj().T@R-R@L.conj().T) for L in Ls))

print("PHASE A / T1.  random operators only -- no lattice, no gauge group, no physical structure.")
print(f"\n  {'case':>44s} {'||[H,R]||':>10s} {'||[L,R]||':>10s} {'||[Ld,R]||':>11s} {'max|d<R>/dt|':>13s}")
print("  "+"-"*92)
for D in (4,6,9,12):
    for kind in ('unitary','normal'):
        R,H,Ls=build(D,kind)
        a,b,c=cc(R,H,Ls)
        lab=f"D={D:<3d} R {kind}, jumps NON-unitary"
        print(f"  {lab:>44s} {a:10.1e} {b:10.1e} {c:11.1e} {worst(R,H,Ls):13.3e}")

print("\n  NON-NORMAL R IS WHERE IT MUST FAIL: there [L,R]=0 does NOT give [L^dag,R]=0.")
for D in (5,8):
    S=rng.normal(size=(D,D))+1j*rng.normal(size=(D,D))            # similarity, not unitary
    d=rng.normal(size=D)+1j*rng.normal(size=D)
    R=S@np.diag(d)@np.linalg.inv(S)                                # diagonalisable, NOT normal
    H=S@np.diag(rng.normal(size=D))@np.linalg.inv(S); H=(H+H.conj().T)/2
    Ls=[S@np.diag(rng.normal(size=D)+1j*rng.normal(size=D))@np.linalg.inv(S) for _ in range(3)]
    a,b,c=cc(R,H,Ls)
    print(f"    D={D}  ||[H,R]||={a:.1e}  ||[L,R]||={b:.1e}  ||[Ld,R]||={c:.1e}"
          f"   max|d<R>/dt| = {worst(R,H,Ls):.3e}")

print("\n  EACH HYPOTHESIS DROPPED IN TURN -- each must break the conclusion.")
D=8; R,H,Ls=build(D,'unitary')
Hb=(lambda X:(X+X.conj().T)/2)(rng.normal(size=(D,D))+1j*rng.normal(size=(D,D)))
Lb=[rng.normal(size=(D,D))+1j*rng.normal(size=(D,D)) for _ in range(3)]
for name,Ht,Lt in (("all hypotheses hold",H,Ls),("drop [H,R]=0",Hb,Ls),("drop [L,R]=0",H,Lb)):
    w=worst(R,Ht,Lt)
    print(f"    {name:>22s}  max|d<R>/dt| = {w:.3e}   {'HOLDS' if w<1e-12 else 'BREAKS'}")
