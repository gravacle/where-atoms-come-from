"""W-55.  ITEM 24. CAN A SET OF BOUNDARIES BE AN OUTPUT RATHER THAN AN INPUT?

Every carrier in this program had its adjacency installed by the registrar, which makes every
reconstruction question circular (W-53). The only way out is a system with NO installed structure:
a bare Hilbert space and a Hamiltonian, with no tensor factorisation given.

THE QUESTION, IN THIS PROGRAM'S OWN TERMS. A "boundary" is a division of the system into a part and
its complement. In a bare Hilbert space that is a TENSOR FACTORISATION. So:

    given H on a space of dimension D with NO preferred factorisation, is there a factorisation
    D = d_A x d_B under which H is LOCAL -- H ~ H_A (x) I + I (x) H_B with small interaction?

If a generic H admits one, boundaries are an output. If it does not, then the boundaries this program
has been using were never derivable and were always an input -- which settles item 24 negatively and
says why.

FORCED-OR-NOT, FIRST. The "local" operators A(x)I and I(x)B span a subspace of dimension
2(d^2 - 1) + 1 inside the d^4-dimensional Hermitian space. For d=4 that is 31 of 256, so a RANDOM
Hermitian matrix already has a guaranteed nonzero projection onto it -- roughly sqrt(31/256) = 0.35
of its norm, by dimension counting alone. ANY result must be read against that floor, not against
zero. The floor is computed here before anything else.
"""
import numpy as np
rng=np.random.default_rng(17)

d=4; D=d*d
def herm(n):
    A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); return (A+A.conj().T)/2
def basis_local(d):
    """orthonormal basis of the LOCAL subspace: I, A(x)I, I(x)B (traceless A,B)"""
    def gens(d):
        G=[]
        for i in range(d):
            for j in range(i+1,d):
                E=np.zeros((d,d),complex); E[i,j]=1; E[j,i]=1; G.append(E/np.sqrt(2))
                F=np.zeros((d,d),complex); F[i,j]=-1j; F[j,i]=1j; G.append(F/np.sqrt(2))
        for k in range(1,d):
            v=np.zeros(d); v[:k]=1; v[k]=-k
            G.append(np.diag(v).astype(complex)/np.sqrt(k*(k+1)))
        return G
    I=np.eye(d,dtype=complex)
    B=[np.kron(I,I)/np.sqrt(D)]
    for g in gens(d): B.append(np.kron(g,I)/np.sqrt(d))
    for g in gens(d): B.append(np.kron(I,g)/np.sqrt(d))
    return B
BL=basis_local(d)
print(f"W-55  D={D} as {d}x{d}.  local subspace dimension {len(BL)} of {D*D} Hermitian dimensions")

def locality(H,V=None):
    """fraction of ||H||^2 lying in the local subspace, after the rotation V (V=None: identity)"""
    M=H if V is None else V.conj().T@H@V
    M=M-np.trace(M)/D*np.eye(D)                     # identity part is trivially local; remove it
    n2=np.real(np.vdot(M,M))
    if n2<1e-14: return 0.0
    p=0.0
    for b in BL[1:]:
        c=np.vdot(b,M); p+=np.real(c*np.conj(c))
    return p/n2

def randU(n):
    Q,R=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    return Q@np.diag(np.diag(R)/np.abs(np.diag(R)))

def optimise(H,iters=4000,step=0.35):
    """maximise locality over factorisations, by random-walk on U(D). Simple, and reported as such."""
    V=np.eye(D,dtype=complex); best=locality(H,V); s=step
    for t in range(iters):
        G=herm(D); W=V@ (np.eye(D)+1j*s*G/np.linalg.norm(G))
        Q,R=np.linalg.qr(W); W=Q@np.diag(np.diag(R)/np.abs(np.diag(R)))
        v=locality(H,W)
        if v>best: best, V = v, W
        if t%800==799: s*=0.6
    return best,V

print("\n  FLOOR FIRST -- what does a RANDOM rotation already give? Any claim is read against this.")
Hr=herm(D)
rand_vals=[locality(Hr,randU(D)) for _ in range(300)]
print(f"    random factorisations of a random H:  mean {np.mean(rand_vals):.4f}  "
      f"max of 300 {np.max(rand_vals):.4f}   (dimension-count floor ~ {31/255:.4f})")

print("\n  CASE 1 -- GENERIC H, no structure installed. Can any factorisation make it local?")
for trial in range(3):
    H=herm(D)
    b,_=optimise(H)
    print(f"    trial {trial}:  best locality found = {b:.4f}   "
          f"(random baseline {np.mean([locality(H,randU(D)) for _ in range(200)]):.4f})")

print("\n  CASE 2 -- H BUILT LOCAL, then hidden by a random rotation. The search must recover it.")
for eps in (0.0,0.15,0.4):
    A=herm(d); B=herm(d); C=herm(D)
    H0=np.kron(A,np.eye(d))+np.kron(np.eye(d),B)
    H0=H0/np.linalg.norm(H0)
    H0=H0+eps*C/np.linalg.norm(C)
    true=locality(H0)
    U=randU(D); H=U@H0@U.conj().T                    # hide the factorisation
    b,_=optimise(H)
    print(f"    interaction eps={eps:4.2f}:  true locality {true:.4f}   "
          f"recovered by search {b:.4f}   hidden-basis value {locality(H):.4f}")

print("\n  READING")
print("    If CASE 1 stays near the random floor, a generic Hamiltonian admits NO factorisation that")
print("    makes it local: boundaries are not derivable from the dynamics and must be an input.")
print("    If CASE 2 recovers its planted structure, the search works and CASE 1's null is real.")
