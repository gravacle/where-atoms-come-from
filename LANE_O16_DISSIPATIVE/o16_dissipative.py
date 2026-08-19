"""O-16: the dissipative arm. Pre-registered at commit 68ff832.

Every prior arm is CLOSED UNITARY evolution with a small bath, where information measures are
the square of an amplitude by construction. A Lindblad bath is the one place amplitudes and
rates need not stand in that relation.

  1. record DECAY exponent under single-site Lindblad jumps: predicted 2 in lambda
     (first order in gamma = lambda^2), INDEPENDENT of d
  2. STATIC SPLITTING exponent under the same operator set: predicted n*
  3. therefore the two DIFFER and the rival's squared rule does not carry over

Carriers: [[5,1,3]] (dim 32, non-CSS, d=3) and the toric 2x2 code (dim 256, d=2)."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); Y=1j*(X@Z)
P1={'I':I2,'X':X,'Y':Y,'Z':Z}
def pauli(s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,P1[c])
    return M

def five_qubit():
    gens=['XZZXI','IXZZX','XIXZZ','ZXIXZ']
    S=[pauli(g) for g in gens]
    Zb=pauli('ZZZZZ'); Xb=pauli('XXXXX')
    return 5,S,Zb,Xb
def steane():
    """[[7,1,3]] CSS. Logicals are the all-Z and all-X strings."""
    gens=['IIIXXXX','IXXIIXX','XIXIXIX','IIIZZZZ','IZZIIZZ','ZIZIZIZ']
    return 7,[pauli(g) for g in gens],pauli('ZZZZZZZ'),pauli('XXXXXXX')
def four_two_two():
    """[[4,2,2]] CSS, d=2."""
    return 4,[pauli('XXXX'),pauli('ZZZZ')],pauli('ZZII'),pauli('XXII')
def toric22():
    nq=8; ind={}; k=0
    for j in range(2):
        for i in range(2): ind[('h',i,j)]=k; k+=1; ind[('v',i,j)]=k; k+=1
    STAR=[[ind[('h',i,j)],ind[('h',(i-1)%2,j)],ind[('v',i,j)],ind[('v',i,(j-1)%2)]] for j in range(2) for i in range(2)]
    PLQ=[[ind[('h',i,j)],ind[('v',(i+1)%2,j)],ind[('h',i,(j+1)%2)],ind[('v',i,j)]] for j in range(2) for i in range(2)]
    def op(d):
        M=np.array([[1]],dtype=complex)
        for l in range(nq): M=np.kron(M,d.get(l,I2))
        return M
    S=[op({l:X for l in s}) for s in STAR[:3]]+[op({l:Z for l in p}) for p in PLQ[:3]]
    Zb=op({ind[('h',0,0)]:Z, ind[('h',0,1)]:Z}); Xb=op({ind[('v',0,0)]:X, ind[('h',0,0)]:X})
    return nq,S,Zb,Xb

def single_site(nq, letters):
    out=[]
    for q in range(nq):
        for c in letters:
            out.append(pauli(''.join(c if j==q else 'I' for j in range(nq))))
    return out

def code_projector(S,n):
    P=np.eye(2**n,dtype=complex)
    for s in S: P=P@((np.eye(2**n)+s)/2)
    return P

def decay_rate(Zb,jumps,rho):
    """d<Zbar>/dt at t=0 under the Lindblad dissipator (no Hamiltonian part: Zbar commutes with H)."""
    d=np.zeros_like(rho)
    for L in jumps:
        Ld=L.conj().T; d+=L@rho@Ld-0.5*(Ld@L@rho+rho@Ld@L)
    return np.real(np.trace(d@Zb))

def splitting(H0,V,eps,k):
    # the splitting WITHIN the code space is w[k-1] - w[0].  v1 used w[k], which is the GAP to
    # the first excited state -- it barely moves with eps and gave a slope of -0.05.
    w=np.linalg.eigvalsh(H0+eps*V); return w[k-1]-w[0]

for name,builder in (("[[4,2,2]]  (d=2, CSS)",four_two_two),("[[5,1,3]] (d=3, non-CSS)",five_qubit),("[[7,1,3]] Steane (d=3, CSS)",steane)):
    n,S,Zb,Xb=builder()
    N=2**n; P=code_projector(S,n); k=int(round(np.real(np.trace(P))))
    H0=-sum(S)
    say("="*96); say(f"  {name}   qubits {n}   dim {N}   code-space dim {k}"); say("="*96)
    c1=np.linalg.norm(Zb@H0-H0@Zb); c2=np.linalg.norm(Zb@Xb+Xb@Zb)
    say(f"  SELF-CHECK ||[Zbar,H0]|| = {c1:.2e}   ||{{Zbar,Xbar}}|| = {c2:.2e}   "
        f"{'PASS' if max(c1,c2)<1e-9 else 'FAIL -- ARM VOID'}")
    if max(c1,c2)>1e-9: say("  skipping: a nominated logical is not a logical\n"); continue
    # a DEFINITE record inside the code space
    Pp=P@((np.eye(N)+Zb)/2); tr=np.real(np.trace(Pp))
    rho=Pp/tr
    say(f"  initial <Zbar> = {np.real(np.trace(rho@Zb)):.6f}   code weight = {np.real(np.trace(P@rho)):.6f}")
    for letters in ('Z','XYZ'):
        say(f"\n  --- single-site coupling set {{{','.join(letters)}}} ---")
        ops=single_site(n,letters)
        say(f"  {'lambda':>9}{'|d<Zbar>/dt| (Lindblad)':>26}{'static splitting dE':>22}")
        rows=[]
        for lam in (0.02,0.05,0.1,0.2):
            g=lam*lam
            jumps=[np.sqrt(g)*o for o in ops]
            r=abs(decay_rate(Zb,jumps,rho))
            V=sum(o for o in ops)
            dE=splitting(H0,V,lam,k)
            rows.append((lam,r,dE)); say(f"  {lam:>9.3f}{r:>26.6e}{dE:>22.6e}")
        lg=np.log([x[0] for x in rows])
        sd=np.polyfit(lg,np.log([x[1] for x in rows]),1)[0]
        ss=np.polyfit(lg[:3],np.log([x[2] for x in rows][:3]),1)[0]
        say(f"  DECAY exponent in lambda    = {sd:.4f}   (predicted 2, independent of d)")
        say(f"  SPLITTING exponent in lambda = {ss:.4f}")
        say(f"  -> exponents {'DIFFER' if abs(sd-ss)>0.3 else 'MATCH'}")
    say("")
