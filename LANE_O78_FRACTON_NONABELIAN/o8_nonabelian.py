"""O-8 PART B: NON-ABELIAN topological order (Levin-Wen / Fibonacci) against R1-R3.
Everything is computed from modular data (S-matrix, topological spins).  The Levin-Wen
string-net model with input category C realises the DOUBLE Z(C) = C (x) Cbar, so the
lattice-model numbers are the doubled ones."""
import numpy as np, itertools
np.set_printoptions(precision=4, suppress=True, linewidth=140)

def hr(t=""): print("\n" + "="*78); print(t); print("="*78)

PHI = (1+5**0.5)/2

def mtc(name, S, theta, labels):
    S = np.array(S, dtype=complex)
    return dict(name=name, S=S, theta=np.array(theta, dtype=complex), labels=list(labels))

Z2 = mtc("Z2 gauge theory = toric code = D(Z_2)",
         0.5*np.array([[1,1,1,1],[1,1,-1,-1],[1,-1,1,-1],[1,-1,-1,1]]),
         [1,1,1,-1], ["1","e","m","f"])
FIB = mtc("Fibonacci (chiral)",
          (1/np.sqrt(1+PHI**2))*np.array([[1,PHI],[PHI,-1]]),
          [1, np.exp(4j*np.pi/5)], ["1","t"])
ISING = mtc("Ising (chiral)",
            0.5*np.array([[1,np.sqrt(2),1],[np.sqrt(2),0,-np.sqrt(2)],[1,-np.sqrt(2),1]]),
            [1, np.exp(1j*np.pi/8), -1], ["1","s","p"])

def double(M):
    S = np.kron(M["S"], M["S"].conj())
    th = np.kron(M["theta"], M["theta"].conj())
    lab = [f"({a},{b})" for a in M["labels"] for b in M["labels"]]
    return mtc("DOUBLE of " + M["name"], S, th, lab)

def qdims(M):
    S = M["S"]; d = (S[0,:]/S[0,0]).real; D = (1/S[0,0]).real
    return d, D

def verlinde(M):
    S = M["S"]; n = S.shape[0]
    N = np.zeros((n,n,n))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                N[a,b,c] = np.sum(S[a,:]*S[b,:]*S[c,:].conj()/S[0,:]).real
    return N

def gsd(M, g):
    d, D = qdims(M)
    return float(np.sum((D/d)**(2*g-2)))

def is_pow2(x, tol=1e-6):
    if x < 1: return False
    r = round(np.log2(x))
    return abs(2.0**r - x) < tol

hr("PART B0 -- SELF-CHECKS ON THE MODULAR DATA (must all PASS)")
ALL = [Z2, FIB, ISING, double(FIB), double(ISING), double(Z2)]
allok = True
for M in ALL:
    S = M["S"]; n = S.shape[0]
    u = np.linalg.norm(S@S.conj().T - np.eye(n))
    sym = np.linalg.norm(S - S.T)
    d, D = qdims(M)
    dcheck = abs(np.sum(d**2) - D**2)
    N = verlinde(M)
    intres = float(np.max(np.abs(N - np.round(N))))
    nonneg = bool(np.all(np.round(N) >= -1e-9))
    g1 = gsd(M,1)
    ok = (u<1e-9 and sym<1e-9 and dcheck<1e-8 and intres<1e-8 and nonneg and abs(g1-n)<1e-8)
    allok = allok and ok
    print(f"{M['name'][:44]:46s} n={n}  ||SS*-I||={u:.2e}  ||S-S^T||={sym:.2e}")
    print(f"    sum d^2 == D^2 residual {dcheck:.2e} | Verlinde N integrality {intres:.2e} | "
          f"N>=0 {nonneg} | GSD(g=1)={g1:.6f} == n : {'PASS' if abs(g1-n)<1e-8 else 'FAIL'}")
print(f"\n  KNOWN FUSION RULES (independent answers):")
Nf = verlinde(FIB); Ni = verlinde(ISING); Nz = verlinde(Z2)
tests = [("Fib  N_tt^1 = 1", Nf[1,1,0], 1), ("Fib  N_tt^t = 1", Nf[1,1,1], 1),
         ("Ising N_ss^1 = 1", Ni[1,1,0], 1), ("Ising N_ss^p = 1", Ni[1,1,2], 1),
         ("Ising N_ss^s = 0", Ni[1,1,1], 0), ("Ising N_pp^1 = 1", Ni[2,2,0], 1),
         ("Z2   N_ee^1 = 1", Nz[1,1,0], 1), ("Z2   N_em^f = 1", Nz[1,2,3], 1)]
for lab, got, want in tests:
    p = abs(got-want) < 1e-8; allok = allok and p
    print(f"    {lab:20s} got {got: .6f}  {'PASS' if p else 'FAIL'}")
print(f"\n  POSITIVE CONTROL -- toric code GSD must be 4^g:")
for g in (1,2,3):
    got = gsd(Z2,g); want = 4**g; p = abs(got-want)<1e-6; allok = allok and p
    print(f"    g={g}: got {got:.6f}  want {want}  {'PASS' if p else 'FAIL'}")
print(f"\n  ALL SELF-CHECKS: {'PASS' if allok else 'FAIL'}")

hr("PART B1 -- R1: IS THE RECORD SPACE 2^k ?  (the F_2 chain-complex signature)")
print("dim H_k of ANY F_2 chain complex is 2^(dim H_k) as a Hilbert-space dimension:")
print("it is ALWAYS a power of 2.  A ground space whose dimension is not a power of 2")
print("cannot be the homology of an F_2 complex, by dimension count alone.\n")
print(f"{'theory':40s} " + "".join(f"g={g:<9d}" for g in (1,2,3,4)))
for M in [Z2, double(Z2), FIB, double(FIB), ISING, double(ISING)]:
    row = f"{M['name'][:39]:40s} "
    for g in (1,2,3,4):
        v = gsd(M,g)
        tag = "2^k" if is_pow2(v) else "NOT2^k"
        row += f"{round(v):<5d}{tag:<6s}"
    print(row)

hr("PART B2 -- R2: THE PARITY OBSTRUCTION  (new, and it is one line)")
print("""CLAUSE (i)+(iv): R = R-dagger, R^2 = I, and some admissible U has U-dagger R U = -R.
If the record is carried by the degenerate space E and the writer is a LOGICAL operation
-- i.e. U maps E to itself, which is forced if writing must not create excitations --
then R|E and -R|E are unitarily equivalent, hence

        tr(R|E) = tr(-R|E) = -tr(R|E)  =>  tr(R|E) = 0  =>  dim E IS EVEN.

An F_2 chain complex can never violate this: dim = 2^k, even for every k >= 1.  That is
why this program never met the obstruction.  Non-abelian theories violate it routinely.
""")
print(f"{'theory':40s} " + "".join(f"g={g:<8d}" for g in (1,2,3,4)))
for M in [Z2, double(Z2), FIB, double(FIB), ISING, double(ISING)]:
    row = f"{M['name'][:39]:40s} "
    for g in (1,2,3,4):
        v = round(gsd(M,g))
        row += f"{v:<4d}{'OK ' if v%2==0 else 'ODD':<5s}"
    print(row)
print("""
  ODD = clause (iv) is UNSATISFIABLE by any code-space-preserving unitary whatsoever.
  Not 'we could not find a writer' -- no writer can exist.  A trace argument, not a search.""")

print("\n  POSITIVE CONTROL for the trace argument -- explicit witnesses in even dimension:")
for dim in (2,4,5,9,16,25):
    R = np.diag([1.0]*(dim//2) + [-1.0]*(dim - dim//2))
    tr = np.trace(R)
    # a U with U^dag R U = -R exists iff R is traceless; construct it when dim even
    if dim % 2 == 0:
        U = np.zeros((dim,dim))
        for i in range(dim//2):
            U[i, dim//2+i] = 1; U[dim//2+i, i] = 1
        res = np.linalg.norm(U.conj().T@R@U + R)
        print(f"    dim={dim:3d}  tr R={tr:+.0f}  explicit U built, ||U*RU + R|| = {res:.2e}  "
              f"{'PASS' if res<1e-12 else 'FAIL'}")
    else:
        # brute force over all Hermitian involutions diagonal in some basis: trace can never be 0
        traces = set(abs(dim - 2*m) for m in range(dim+1))
        print(f"    dim={dim:3d}  possible |tr R| over ALL involutions = {sorted(traces)}  "
              f"-> 0 attainable? {0 in traces}  {'PASS (no writer)' if 0 not in traces else 'FAIL'}")

hr("PART B3 -- R2: ARE THE WILSON LOOPS EVEN UNITARY?")
print("""The writer O-8 nominated is a Wilson loop.  Two loop operators act on the torus
ground space in the flux basis |b>:
   MEASURING loop (parallel to the flux):  M_a|b> = (S_ab/S_1b) |b>      -- diagonal
   CREATING  loop (transverse):            W_a|b> = sum_c N_ab^c |c>     -- fusion matrix
Clause (i) wants an INVOLUTION; clause (iv) wants a UNITARY.""")
for M in [Z2, FIB, ISING, double(FIB)]:
    S = M["S"]; n = S.shape[0]; d, D = qdims(M); N = verlinde(M)
    print(f"\n  {M['name']}   (anyons {M['labels']})")
    for a in range(n):
        Ma = np.diag(S[a,:]/S[0,:])
        Wa = np.round(N[a]).astype(float)              # W_a in the flux basis
        inv = np.linalg.norm(Ma@Ma - np.eye(n))
        uni = np.linalg.norm(Wa@Wa.conj().T - np.eye(n))
        ev = np.sort_complex(np.diag(Ma))
        print(f"    a={M['labels'][a]:8s} d_a={d[a]:.4f}  M_a eigenvalues={np.round(ev,4)}")
        print(f"        M_a^2 = I ? residual {inv:.3e} {'INVOLUTION' if inv<1e-9 else 'NOT AN INVOLUTION'}"
              f"   |  W_a unitary ? residual {uni:.3e} {'UNITARY' if uni<1e-9 else 'NOT UNITARY'}")
print("""
  The pattern is exact and provable: W_a W_a^dag = I iff N_a is a permutation matrix
  iff d_a = 1 iff a is ABELIAN.  Every non-abelian charge has a NON-UNITARY Wilson loop,
  and its measuring loop has eigenvalues S_ab/S_1b which are quantum-dimension ratios,
  not signs.  So the nominated writer is not an admissible U and the nominated bit is
  not an involution.""")

hr("PART B4 -- IS THERE ANY WRITER AT ALL WHERE THE PARITY TEST PASSES?")
print("""Doubled Fibonacci on the TORUS has dim 4 (even) so the parity test is silent.
Search: over all 2^4 sign patterns R diagonal in the anyon basis, and over the unitary
group generated by the modular S and T matrices (the mapping class group of the torus,
which IS implementable by adiabatic surface deformation), does any U give U*RU = -R?""")
def mcg_group(M, cap=20000, with_loops=True):
    """Group generated by the modular S and T, plus every Wilson loop that is actually
       UNITARY (only the abelian charges are).  Non-unitary loops are not admissible U."""
    S = M["S"]; th = M["theta"]; T = np.diag(th); n = S.shape[0]
    gens = [S, T, S.conj().T, T.conj().T]
    if with_loops:
        N = verlinde(M)
        for a in range(n):
            Wa = np.round(N[a]).astype(complex)
            Ma = np.diag(S[a,:]/S[0,:])
            if np.linalg.norm(Wa@Wa.conj().T - np.eye(n)) < 1e-9:
                gens += [Wa, Wa.conj().T]
            if np.linalg.norm(Ma@Ma.conj().T - np.eye(n)) < 1e-9:
                gens += [Ma, Ma.conj().T]
    seen = {}; frontier = [np.eye(S.shape[0], dtype=complex)]
    def key(A):
        # projective: quotient by global phase
        i = np.unravel_index(np.argmax(np.abs(A)), A.shape)
        B = A / (A[i]/abs(A[i]))
        return np.round(B, 4).tobytes()
    seen[key(frontier[0])] = frontier[0]
    while frontier and len(seen) < cap:
        nf = []
        for A in frontier:
            for g in gens:
                B = g@A
                k = key(B)
                if k not in seen:
                    seen[k] = B; nf.append(B)
        frontier = nf
    return list(seen.values())

print("  POSITIVE CONTROL: the toric code MUST return a nonzero count -- its writer is the")
print("  Wilson loop W_e.  A zero there would mean the search is broken, not that no writer exists.")
for M in [Z2, double(FIB), double(ISING)]:
    n = M["S"].shape[0]
    G = mcg_group(M)
    print(f"\n  {M['name']}: |group generated by S,T (mod phase, capped)| = {len(G)}")
    hits = 0; examples = []
    for signs in itertools.product([1,-1], repeat=n):
        R = np.diag(np.array(signs, dtype=complex))
        if abs(np.trace(R)) > 1e-9: continue
        for U in G:
            if np.linalg.norm(U.conj().T@R@U + R) < 1e-6:
                hits += 1; examples.append((signs, U)); break
    print(f"    traceless diagonal involutions with a mapping-class-group writer: {hits}")
    if examples:
        s0, U0 = examples[0]
        Rm = np.diag(np.array(s0,dtype=complex))
        print(f"    example R = diag{s0};  ||U*RU + R|| = "
              f"{np.linalg.norm(U0.conj().T@Rm@U0 + Rm):.2e}")
print("\nDONE.")
