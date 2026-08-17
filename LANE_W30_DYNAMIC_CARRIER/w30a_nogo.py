"""W-30a. IS THE WRITABLE/DURABLE CONJUGACY A FEATURE OF OUR PATCH, OR IS IT GENERAL?

W-29 found on the 3x3 patch: a record that survives dissipation cannot be written.
Before treating that as a fact about records, test it on OTHER carriers and other N.

CLAIM UNDER TEST.  For Lindblad  drho/dt = -i[H,rho] + sum_k (L_k rho L_k^dag - rho)
with UNITARY jump operators, if [H,R]=0 and [L_k,R]=0 for all k then d<R>/dt = 0
IDENTICALLY -- every state, every time.  Then a perfectly durable record is one whose
value never changes: it was never written, it always was.

DISCIPLINE (this lane's first version got this wrong).  A loop built as all-(+1) shifts
is only a genuine loop at N=2, where sign is irrelevant mod 2.  At N=3 it is not closed,
it leaves the physical sector, and the builder silently drops those elements -- so R can
come out as the ZERO MATRIX and score [H,R]=0 for no physical reason.  So: build every
loop as a PRODUCT OF PLAQUETTES (correctly oriented by construction, discrete Stokes),
and PRINT ||R|| and the unitarity defect. An operator proves it is real before its
commutator is allowed to mean anything.
"""
import itertools, numpy as np

def build(V, E, N):
    NV = len(V)
    st = [s for s in itertools.product(range(N), repeat=len(E))
          if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
                 -sum(s[k] for k,(a,b) in enumerate(E) if b==v)) % N == 0 for v in range(NV))]
    return st, {s:i for i,s in enumerate(st)}

def Zop(st, links, N):
    w = np.exp(2j*np.pi/N)
    return np.diag([w**(sum(s[k] for k in links) % N) for s in st])

def Move(st, idx, moves, N):
    """ONE simultaneous shift. moves = list of (link, signed amount)."""
    D = len(st); M = np.zeros((D,D), complex)
    for j,s in enumerate(st):
        t = list(s)
        for k,sg in moves: t[k] = (t[k] + sg) % N
        t = tuple(t)
        if t in idx: M[idx[t], j] = 1.0
    return M

def compose(plqs):
    """Product of plaquettes = sum of their signed moves. Interior links cancel: Stokes."""
    acc = {}
    for p in plqs:
        for k,sg in p: acc[k] = acc.get(k,0) + sg
    return [(k,sg) for k,sg in acc.items() if sg != 0]

def carriers():
    out = []
    V2=[(i,j) for j in range(3) for i in range(3)]; vid={v:k for k,v in enumerate(V2)}
    E=[]
    for j in range(3):
        for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))
    for j in range(2):
        for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))
    hid=lambda i,j: j*2+i; vx=lambda i,j: 6+j*3+i
    P=[[(hid(i,j),+1),(vx(i+1,j),+1),(hid(i,j+1),-1),(vx(i,j),-1)] for j in range(2) for i in range(2)]
    for N in (2,3): out.append((f"3x3 patch", V2, E, P, N))
    V=[(i,j) for j in range(2) for i in range(2)]; w={v:k for k,v in enumerate(V)}
    E2=[(w[(0,0)],w[(1,0)]),(w[(0,1)],w[(1,1)]),(w[(0,0)],w[(0,1)]),(w[(1,0)],w[(1,1)])]
    P2=[[(0,+1),(3,+1),(1,-1),(2,-1)]]
    for N in (2,3): out.append((f"2x2 plaquette", V, E2, P2, N))
    E3=[(0,1),(1,2),(0,2)]; P3=[[(0,+1),(1,+1),(2,-1)]]
    for N in (2,3): out.append((f"triangle", [0,1,2], E3, P3, N))
    Et=[(0,1),(0,1),(0,1)]; Pt=[[(0,+1),(1,-1)],[(1,+1),(2,-1)]]
    for N in (2,3): out.append((f"theta graph", [0,1], Et, Pt, N))
    return out

rng = np.random.default_rng(20260817)
def rand_rho(D):
    A = rng.normal(size=(D,D)) + 1j*rng.normal(size=(D,D))
    r = A @ A.conj().T
    return r / np.trace(r).real

def run(g2, header):
    print(header)
    print(f"  {'carrier':14s} {'N':>2s} {'dim':>4s} {'||R||':>7s} {'unitarity':>10s} "
          f"{'||[H,R]||':>10s} {'||[L,R]||':>10s} {'max|d<R>/dt|':>13s}  verdict")
    print("  " + "-"*99)
    for name, V, E, P, N in carriers():
        st, idx = build(V, E, N); D = len(st)
        R  = Move(st, idx, compose(P), N)               # rim loop = product of ALL plaquettes
        Ls = [Move(st, idx, p, N) for p in P]           # jump operators = the plaquettes
        nR = np.linalg.norm(R)
        uni = np.linalg.norm(R@R.conj().T - np.eye(D))  # 0 => genuinely unitary
        H  = -sum(L + L.conj().T for L in Ls) \
             - g2*sum(Zop(st,[k],N) + Zop(st,[k],N).conj().T for k in range(len(E)))
        cHR = np.linalg.norm(H@R - R@H)
        cLR = max(np.linalg.norm(L@R - R@L) for L in Ls)
        worst = 0.0
        for _ in range(20):
            rho = rand_rho(D)
            drho = -1j*(H@rho - rho@H) + sum(L@rho@L.conj().T - rho for L in Ls)
            worst = max(worst, abs(np.trace(R @ drho)))
        if nR < 1e-12:            v = "R IS ZERO -- meaningless"
        elif uni > 1e-9:          v = "R NOT UNITARY -- meaningless"
        elif worst < 1e-9:        v = "FROZEN: <R> can never change"
        else:                     v = "R can move"
        print(f"  {name:14s} {N:2d} {D:4d} {nR:7.3f} {uni:10.3e} {cHR:10.3e} {cLR:10.3e} "
              f"{worst:13.3e}  {v}")
    print()

run(0.0, "PURE MAGNETIC  g2 = 0  -- the point where W-28 measured ||[H,W]|| = 0")
run(1.0, "MIXED          g2 = 1  -- an electric term is switched on")

print("INTERPRETATION CHECK -- is 'frozen' forced by R being trivial on the physical space?")
for name, V, E, P, N in carriers():
    st, idx = build(V, E, N); D = len(st)
    R = Move(st, idx, compose(P), N)
    ev = np.linalg.eigvals(R)
    print(f"  {name:14s} N={N} dim={D:4d}  distinct eigenvalues of R: "
          f"{len(np.unique(np.round(ev,6)))}  -> "
          f"{'TRIVIAL (identity-like), cannot label anything' if len(np.unique(np.round(ev,6)))<2 else 'nontrivial: R separates states'}")
