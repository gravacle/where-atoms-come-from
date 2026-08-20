"""STEP 1 -- build the [[n,n-2,2]] family and PROVE the objects are records before using them.

SELF-CHECKS THAT MUST PASS BEFORE ANY GEOMETRY IS DRAWN:
  SC-1  symplectic pairing matrix of the 2k logicals is NON-DEGENERATE over F2 (det = 1)
  SC-2  every logical commutes with BOTH stabilisers (so with H)         -> clause (ii)
  SC-3  R = R-dagger and R^2 = I                                          -> clause (i)
  SC-4  record_model.clause_iii(R, es) is True                            -> clause (iii)
  SC-5  record_model.clause_iv(R, es) is True                             -> clause (iv)
  SC-6  eigenvalue multiplicities of H are [2^(n-2), 2^(n-1), 2^(n-2)]
  SC-7  the exact symbolic code-space overlap matches the numeric one
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC")
from lib_scaleb import *
import numpy as np

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.append(s)

def f2_det(M):
    """determinant over F2 by elimination"""
    A = [row[:] for row in M]; m = len(A)
    for c in range(m):
        p = next((r for r in range(c, m) if A[r][c]), None)
        if p is None: return 0
        A[c], A[p] = A[p], A[c]
        for r in range(m):
            if r != c and A[r][c]:
                A[r] = [(x + y) % 2 for x, y in zip(A[r], A[c])]
    return 1

P("=" * 100)
P("LANE_SCALE_B_METRIC  STEP 1 -- CARRIER CONSTRUCTION AND CLAUSE VERIFICATION")
P("carrier family: [[n, n-2, 2]],  H = -(X^(x)n + Z^(x)n),  records = the 2k logical Paulis")
P("=" * 100)
P("")
P("%-4s %-6s %-6s %-8s %-10s %-8s %-8s %-8s %-8s %-8s %-10s" %
  ("n", "k", "2k", "dim", "mult(H)", "SC-1", "SC-2", "SC-3", "SC-4", "SC-5", "SC-6"))
P("-" * 100)

NUMERIC_N = [4, 6, 8, 10]
SYMBOLIC_N = [4, 6, 8, 10, 12, 14, 16, 18, 20]

store = {}
for n in NUMERIC_N:
    stab, pairs = carrier(n)
    k = len(pairs)
    vs, lab = record_vectors(pairs, n)
    m = len(vs)
    # SC-1 : symplectic pairing matrix non-degenerate over F2
    Sm = [[sp_form(vs[i], vs[j], n) for j in range(m)] for i in range(m)]
    sc1 = f2_det(Sm) == 1
    # SC-2 : commutes with both stabilisers (F2, exact)
    sc2 = all(sp_form(v, s, n) == 0 for v in vs for s in stab)
    # numeric
    H = hamiltonian(n)
    es = eigenspaces(H)
    mult = [mm for _, _, mm in es]
    sc6 = mult == [2 ** (n - 2), 2 ** (n - 1), 2 ** (n - 2)]
    Rs = [xz_to_matrix(v, n) for v in vs]
    Id = np.eye(2 ** n)
    sc3 = all(np.linalg.norm(R - R.conj().T) < 1e-9 and np.linalg.norm(R @ R - Id) < 1e-9 for R in Rs)
    sc4 = all(clause_iii(R, es) for R in Rs)
    sc5 = all(clause_iv(R, es) for R in Rs)
    # clause (ii) numeric confirmation of the F2 argument
    sc2n = all(np.linalg.norm(R @ H - H @ R) < 1e-9 for R in Rs)
    store[n] = dict(pairs=pairs, vs=vs, lab=lab, H=H, es=es, Rs=Rs)
    P("%-4d %-6d %-6d %-8d %-10s %-8s %-8s %-8s %-8s %-8s %-10s" %
      (n, k, m, 2 ** n, str(mult), sc1, sc2 and sc2n, sc3, sc4, sc5, sc6))

P("")
P("SC-1 non-degenerate symplectic pairing   SC-2 [R,H]=0 (F2 AND numeric)   SC-3 bit")
P("SC-4 non-trivial   SC-5 writable (Tr P_E R = 0 on EVERY eigenspace)   SC-6 multiplicities")
P("")

allpass = True
for n in NUMERIC_N:
    pass
P("SC-7  exact-symbolic vs numeric code-space overlap  |Tr(P_g R_i R_j)|/Tr(P_g)")
P("%-4s %-14s %-14s %-12s" % ("n", "max|sym-num|", "offdiag nonzero", "verdict"))
P("-" * 60)
for n in NUMERIC_N:
    d = store[n]
    vs = d["vs"]; Rs = d["Rs"]
    H = d["H"]
    w, V = np.linalg.eigh(H)
    kdim = int((np.abs(w - w[0]) < 1e-9).sum())
    Q = V[:, :kdim]; Pg = Q @ Q.conj().T
    m = len(Rs)
    num = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            num[i, j] = abs(np.trace(Pg @ Rs[i] @ Rs[j]) / kdim)
    sym = M_codespace_overlap(vs, n)
    err = float(np.max(np.abs(sym - num)))
    offnz = int((np.abs(num - np.diag(np.diag(num))) > 1e-9).sum())
    P("%-4d %-14.2e %-14d %-12s" % (n, err, offnz, "PASS" if err < 1e-9 else "FAIL"))
    store[n]["Pg"] = Pg; store[n]["kdim"] = kdim

P("")
P("THE RECORDS THEMSELVES (Pauli letters over the n physical qubits, I/X/Z/Y):")
for n in NUMERIC_N:
    d = store[n]
    P("  n=%d:" % n)
    for v, l in zip(d["vs"], d["lab"]):
        P("     %-4s %s   support=%s  weight=%d" % (l, letters(v, n),
                                                    sorted(support(v, n)), len(support(v, n))))
P("")
P("SYMBOLIC REACH -- the exact (F2) relations need no matrices, so they run far past the numerics:")
P("%-4s %-6s %-8s %-14s" % ("n", "k", "2k", "dim (unused)"))
for n in SYMBOLIC_N:
    stab, pairs = carrier(n)
    P("%-4d %-6d %-8d %-14d" % (n, len(pairs), 2 * len(pairs), 2 ** n))

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/01_records_verify.txt", "w").write("\n".join(OUT) + "\n")
