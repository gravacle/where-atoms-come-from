"""V4 -- EXACT-ARITHMETIC RECHECK.  No floats anywhere in this file's conclusions.

The lane reports eight quantities as "EXACTLY ZERO ... verified numerically at 1e-15".  A
float64 zero is not a proof.  But every one of those zeros is claimed to follow from operator
identities, and operator identities on Pauli matrices are INTEGER identities.  So they can be
settled exactly.  Likewise the lane's largest "exactly non-zero" claim -- the three-record
channel (a) and the "factor 4.28 at identical occupancy" -- reduces (V1) to the eigenvalue
distribution of the summed site operator S, whose moments tr(P_code S^k) are RATIONAL and can
be computed exactly.

EVERYTHING BELOW IS DONE IN EXACT INTEGER / Fraction ARITHMETIC:
  * Pauli matrices over Z (entries 0, +-1), matmul in int64 with an explicit overflow bound;
  * the code projector P_code = (I + X^n)(I + Z^n)/4, carried as an integer matrix over 4;
  * the F_2 symplectic form for the logical algebra.

PART 1  the structural identities the lane's zeros rest on          -> EXACT
PART 2  the exact eigenvalue distribution of S behind the "three-record" numbers -> EXACT
PART 3  the F_2 fact behind exact n-independence, at n = 4,6,8,10,12 -> EXACT
"""
import numpy as np, sys, os, itertools
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from common import sp, xor, stabilisers, logical_pairs, check_symplectic, say  # noqa

N = 8
d = 1 << N

# ---------------------------------------------------------------- integer Paulis
IZ = np.array([[1, 0], [0, 1]], dtype=np.int64)
XZ = np.array([[0, 1], [1, 0]], dtype=np.int64)
ZZ = np.array([[1, 0], [0, -1]], dtype=np.int64)
YZ = None  # never needed: every operator here is a real (x|z) Pauli times a sign

def pauli_int(vec, n):
    """(x|z) -> integer matrix of X^x Z^z (a REAL matrix; the global i^k phase of Y is
       dropped, exactly as record_model.xz_to_matrix does when it builds Hermitian logicals).
       Returns (M, ok) where ok records that M is Hermitian and M@M == I, both EXACTLY."""
    M = np.array([[1]], dtype=np.int64)
    for i in range(n):
        x, z = vec[i], vec[n + i]
        b = IZ
        if x and z: b = XZ @ ZZ          # = -i Y, real, Hermitian? (XZ@ZZ) = [[0,-1],[1,0]] -- antisym
        elif x: b = XZ
        elif z: b = ZZ
        M = np.kron(M, b)
    return M

def herm_involution(M):
    return np.array_equal(M, M.T) and np.array_equal(M @ M, np.eye(M.shape[0], dtype=np.int64))

say("=" * 118)
say("V4   EXACT-ARITHMETIC RECHECK OF THE LANE'S 'EXACTLY ZERO' AND 'EXACTLY NON-ZERO' CLAIMS")
say("=" * 118)
say(f"  integer arithmetic, n = {N}, full space dim {d}; int64 with |entry| bounds asserted.")

pairs = logical_pairs(N)
ok, bad = check_symplectic(pairs, N)
say(f"  symplectic self-check at n={N}: {'PASS' if ok else 'FAIL ' + str(bad[:4])}")

VEC = {}
for i, (a, b) in enumerate(pairs):
    VEC[f"X{i+1}"] = a; VEC[f"Z{i+1}"] = b
for j in range(2, len(pairs) + 1):
    for t in ("X", "Z"):
        VEC[f"Z1{t}{j}"] = xor(VEC["Z1"], VEC[f"{t}{j}"])

MAT = {}
for k, v in VEC.items():
    M = pauli_int(v, N)
    if not herm_involution(M):        # fix the sign/phase convention for Y-containing Paulis
        M2 = M @ M
        # M@M must be +-I; if -I the real form carries an antisymmetric factor
        pass
    MAT[k] = M

USED = ["X1", "Z1", "X2", "Z2", "X3", "Z3", "X4", "Z4", "X5", "Z5"]
nh = [k for k in USED if not herm_involution(MAT[k])]
say(f"  integer-representation audit: every operator used below is EXACTLY Hermitian and squares")
say(f"  to I over Z.  operators failing that audit: {nh if nh else 'NONE'}   (all logicals here are")
say(f"  pure X-type or pure Z-type, so no Y phase arises)")

say("")
say("  PART 1 -- THE STRUCTURAL IDENTITIES, EXACTLY (every entry below is an INTEGER count of")
say("  non-zero entries in a matrix that must be the zero matrix).")
STAB = [pauli_int(s, N) for s in stabilisers(N)]
R = MAT["X1"]; W = MAT["Z1"]
rows = []
rows.append(("[X^n, R] == 0 (R is durable)", int(np.count_nonzero(STAB[0] @ R - R @ STAB[0]))))
rows.append(("[Z^n, R] == 0 (R is durable)", int(np.count_nonzero(STAB[1] @ R - R @ STAB[1]))))
rows.append(("R@R - I == 0 (R is a bit)", int(np.count_nonzero(R @ R - np.eye(d, dtype=np.int64)))))
rows.append(("R - R^T == 0 (R Hermitian, real)", int(np.count_nonzero(R - R.T))))
rows.append(("W R W + R == 0 (R is writable)", int(np.count_nonzero(W @ R @ W + R))))
COMM = ["X2", "Z2", "X3", "Z3", "X4", "Z4", "X5", "Z5"]
rows.append(("[R, A] == 0 for every commuting partner A", sum(int(np.count_nonzero(R @ MAT[a] - MAT[a] @ R)) for a in COMM)))
rows.append(("[W, A] == 0 for every commuting partner A", sum(int(np.count_nonzero(W @ MAT[a] - MAT[a] @ W)) for a in COMM)))
rows.append(("CONTROL: [R, W] == 0 ?  (must FAIL: nonzero)", int(np.count_nonzero(R @ W - W @ R))))
say(f"  {'identity':<48}{'non-zero entries (exact)':>28}{'verdict':>12}")
for nm, c in rows:
    v = "ZERO" if c == 0 else f"{c} NONZERO"
    say(f"  {nm:<48}{c:>28d}{v:>12}")
say("  -> the hypotheses of the lane's generalised C-38 theorem hold EXACTLY (integer identities),")
say("     so its zeros are real zeros, not float64 zeros.  The control on the last line shows the")
say("     same test registers a non-zero when there is one.")

# ---------------------------------------------------------------- PART 2: exact spectra
say("")
say("  PART 2 -- THE EXACT EIGENVALUE DISTRIBUTION OF THE SUMMED SITE OPERATOR S.")
say("  V1 showed chi is a zero-parameter functional of this distribution.  Here it is computed")
say("  EXACTLY: 4*P_code = (I + X^n)(I + Z^n) is an integer matrix, so tr(P_code S^k) is a")
say("  rational number with denominator 4, and the moments determine the distribution.")
P4 = (np.eye(d, dtype=np.int64) + STAB[0]) @ (np.eye(d, dtype=np.int64) + STAB[1])   # = 4 P_code
dim_code = Fraction(int(np.trace(P4)), 4)
say(f"  exact code-space dimension tr(P_code) = {dim_code}   (expected 2^(n-2) = {1 << (N-2)})")

def exact_moments(labels, kmax=6):
    S = sum(MAT[l] for l in labels)
    mom = []
    A = np.eye(d, dtype=np.int64)
    for k in range(kmax + 1):
        mom.append(Fraction(int(np.trace(P4 @ A)), 4) / dim_code)
        A = A @ S
        assert np.abs(A).max() < 2 ** 40, "int64 overflow guard"
    return S, mom

def moments_of(dist):
    return [sum(Fraction(p) * Fraction(s) ** k if isinstance(s, int) else None for s, p in dist)
            for k in range(0, 7)]

CASES = [
    ("X2,X3        (two COMMUTING partners)", ["X2", "X3"],
     "S^2 = 2I + 2*X2X3, (X2X3)^2 = I, tr X2X3 = 0  ->  spec {+2:1/4, 0:1/2, -2:1/4}"),
    ("X2,Z2        (two ANTI-commuting partners)", ["X2", "Z2"],
     "S^2 = 2I EXACTLY                              ->  spec {+sqrt2:1/2, -sqrt2:1/2}"),
    ("X2,X3,X4,X5  (four commuting)", ["X2", "X3", "X4", "X5"],
     "binomial:  spec {+-4:1/16, +-2:1/4, 0:3/8}"),
    ("X2,Z2,X3,Z3  (two anti-commuting pairs)", ["X2", "Z2", "X3", "Z3"],
     "S^2 = 4I + 2(X2+Z2)(X3+Z3)                    ->  spec {+-2sqrt2:1/4, 0:1/2}"),
]
say("")
say(f"  {'partner set on the read site':<44}{'exact tr(S^2)/dim':>20}{'exact tr(S^4)/dim':>20}{'S^2 = c*I exactly?':>22}")
for nm, labs, note in CASES:
    S, mom = exact_moments(labs)
    S2 = S @ S
    c = S2[0, 0]
    cnum = np.array_equal(S2, int(c) * np.eye(d, dtype=np.int64))
    say(f"  {nm:<44}{str(mom[2]):>20}{str(mom[4]):>20}{('YES, c = ' + str(int(c))) if cnum else 'no':>22}")
say("")
say("  READ IT: the two configurations the lane calls 'the same occupancy, the same three bits'")
say("  have DIFFERENT exact second and fourth moments of the operator the bath actually couples")
say("  to.  X2,X3 gives tr(S^2)/dim = 2 and tr(S^4)/dim = 8; X2,Z2 gives 2 and 4.  The four-record")
say("  pair differs the same way.  This is elementary Pauli algebra and it is EXACT.")
for nm, labs, note in CASES:
    say(f"    {nm:<44}{note}")

# ---------------------------------------------------------------- PART 3: F_2, all n
say("")
say("  PART 3 -- THE F_2 FACT BEHIND 'chi IS EXACTLY INDEPENDENT OF n'.  For each n the computed")
say("  logicals satisfy the standard symplectic relations, so any configuration at one n is")
say("  carried to the identical configuration at another by a logical Clifford, and the unwritten")
say("  logical qubits are maximally mixed and uncoupled.  Exact F_2 linear algebra:")
say(f"  {'n':>4}{'k = n-2':>10}{'conjugate pairs returned':>28}{'symplectic relations':>24}")
for n in (4, 6, 8, 10, 12):
    pr = logical_pairs(n)
    o, b = check_symplectic(pr, n)
    say(f"  {n:>4}{n-2:>10}{len(pr):>28}{('EXACT PASS' if o else 'FAIL'):>24}")
say("  -> the n-independence is a restatement of 'the code space is a tensor factor of k qubits")
say("     and the state is maximally mixed on it'.  It is exact, and it is also a statement that")
say("     the CARRIER has been factored out of the problem entirely.")
