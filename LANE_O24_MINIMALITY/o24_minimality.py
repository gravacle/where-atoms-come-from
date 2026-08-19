"""O-24: DOES MINIMALITY PICK THE CONNECTION, AND WHAT DOES CURVATURE COST?

O-23 showed the connection is not determined by the records: the same family gives holonomy 0 for
one admissible writer and 5.169 for another. T-20's flat connection was the MINIMAL-WEIGHT writer.
If the physics uses the cheapest operation that flips a record, the flat connection is SELECTED --
and curvature has a price in weight, which is a number.

Carrier: three qubits, H = 0, records Z_1 Z_2 Z_3. Weight is defined (sites acted on non-trivially)
and every admissible writer can be enumerated. A writer for record i must flip Z_i and fix the
others; that is the whole constraint, and it is checked on every candidate."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); Xm=np.array([[0,1],[1,0]],dtype=complex); Zm=np.array([[1,0],[0,-1]],dtype=complex)
def op(d,n=3):
    M=np.array([[1]],dtype=complex)
    for k in range(n): M=np.kron(M,d.get(k,I2))
    return M
n=3; N=8
R=[op({i:Zm}) for i in range(n)]
def is_writer(U,idx):
    if np.linalg.norm(U.conj().T@U-np.eye(N))>1e-9: return False
    for i in range(n):
        want = -R[i] if i==idx else R[i]
        if np.linalg.norm(U.conj().T@R[i]@U - want)>1e-9: return False
    return True
def weight(U,tol=1e-9):
    """sites U acts on non-trivially: it acts trivially on site k iff it commutes with BOTH X_k and Z_k"""
    w=0
    for k in range(n):
        if (np.linalg.norm(U@op({k:Xm})-op({k:Xm})@U)>tol or
            np.linalg.norm(U@op({k:Zm})-op({k:Zm})@U)>tol): w+=1
    return w
say("="*100); say("O-24   DOES MINIMALITY PICK THE CONNECTION?  WHAT DOES CURVATURE COST?"); say("="*100)
say(f"  carrier: {n} qubits, H = 0, records Z_1 Z_2 Z_3, dim {N}")
# candidate writers for record i: X_i times any DIAGONAL unitary (diagonal preserves every Z label)
def phase_op(sites, angs):
    D=np.ones(N,dtype=complex)
    for b in range(N):
        bits=[(b>>(n-1-k))&1 for k in range(n)]
        ph=0.0
        for s,a in zip(sites,angs): ph += a*bits[s]
        D[b]=np.exp(1j*ph)
    return np.diag(D)
say("")
say("  ENUMERATING ADMISSIBLE WRITERS BY WEIGHT")
say(f"  {'writer for R_0':<34}{'weight':>8}{'admissible':>12}")
cands=[]
base=op({0:Xm})
cands.append(("X_0  (the minimal writer)", base))
for s in range(n):
    U=base@phase_op([s],[np.pi/2])
    cands.append((f"X_0 * S_{s}   (phase on site {s})", U))
U=base@phase_op([1,2],[np.pi/2,np.pi/3])
cands.append(("X_0 * S_1 * S_2", U))
for lbl,U in cands:
    say(f"  {lbl:<34}{weight(U):>8}{str(is_writer(U,0)):>12}")
say("")
say("  HOLONOMY OF EACH CHOICE, against the minimal writers for R_1 and R_2")
W1=op({1:Xm}); W2=op({2:Xm})
say(f"  {'writer for R_0':<34}{'weight':>8}{'||[U0,W1]||':>13}{'loop ||H-I||':>14}{'verdict':>10}")
best_flat=None; best_curved=None
for lbl,U0 in cands:
    if not is_writer(U0,0): continue
    w=weight(U0)
    c=float(np.linalg.norm(U0@W1-W1@U0))
    Hol=U0@W1@np.linalg.inv(U0)@np.linalg.inv(W1)
    dev=float(np.linalg.norm(Hol-np.eye(N)))
    curved = dev>1e-9
    if curved and (best_curved is None or w<best_curved): best_curved=w
    if not curved and (best_flat is None or w<best_flat): best_flat=w
    say(f"  {lbl:<34}{w:>8}{c:>13.3e}{dev:>14.3e}{('CURVED' if curved else 'FLAT'):>10}")
say("")
say(f"  minimum weight of a FLAT   writer for R_0 : {best_flat}")
say(f"  minimum weight of a CURVED writer for R_0 : {best_curved}")
if best_flat is not None and best_curved is not None:
    say(f"  -> CURVATURE COSTS {best_curved-best_flat} EXTRA SITE(S) OF WEIGHT")
    say("")
    say("  READ: the cheapest operation that flips a record is FLAT. Curvature requires reaching")
    say("  a site the minimal writer does not touch. If the physics uses the cheapest writer, the")
    say("  flat connection is SELECTED rather than assumed -- and curvature has a price.")
