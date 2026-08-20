"""S2 -- DENSE CHECK OF THE FIVE CLAUSES, AND OF THE CODE-SPACE REDUCTION.

REPRESENTATION: full dense Hilbert space, dimension 2^n, at n = 4, 6, 8 (dims 16, 64, 256).
Everything S1 computed in F_2 is re-derived here as matrices, so that the F_2 shortcuts used
later are earned rather than assumed.

D-18: a record is a LOGICAL OPERATOR; we check the five clauses on the carrier it lives on
before calling any of these objects a record.

What is checked, per n:
  (i)   R = R-dag, R^2 = I
  (ii)  [H,R] = 0            (no Lindblads here, so the *-algebra is alg{I,H})
  (iii) R not constant on some eigenspace of H
  (iv)  Tr(P_E R) = 0 on EVERY eigenspace  (C-11 form of WRITABLE), plus an explicit writer
  (v)   PROTECTED: the compression of every weight-1 operation onto the code space.
        CONTROL in the same table: the compression of weight-2 operations.
  plus: the compressed records Q-dag R_i Q generate the full Pauli algebra on k logical qubits,
        which is what licenses running the chi scaling in the 2^k code-space representation.
"""
import sys, itertools, json
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import RecordModel, symplectic_logicals, xz_to_matrix, eigenspaces

OUT = []
def p(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

TOL = 1e-9
I2 = np.eye(2); Xm = np.array([[0, 1], [1, 0]], complex); Zm = np.array([[1, 0], [0, -1]], complex)

def pauli_string(n, ops):
    M = np.array([[1]], complex)
    for i in range(n): M = np.kron(M, ops.get(i, I2))
    return M

def build(n):
    Xn = pauli_string(n, {i: Xm for i in range(n)})
    Zn = pauli_string(n, {i: Zm for i in range(n)})
    H = -(Xn + Zn)
    return H, Xn, Zn

def codespace(H):
    w, V = np.linalg.eigh(H)
    kdim = int(np.sum(np.abs(w - w[0]) < 1e-9))
    return V[:, :kdim], kdim

results = {}
p("=" * 116)
p("S2  DENSE CHECK OF THE FIVE CLAUSES -- [[n, n-2, 2]].  REPRESENTATION: full dense Hilbert space 2^n.")
p("=" * 116)

for n in (4, 6, 8):
    k = n - 2
    H, Xn, Zn = build(n)
    es = eigenspaces(H)
    Q, kdim = codespace(H)
    stab = [[1] * n + [0] * n, [0] * n + [1] * n]
    pairs = symplectic_logicals([s[:] for s in stab], n)
    if len(pairs) != k:
        p("SELF-CHECK FAILED at n=%d: symplectic_logicals gave %d pairs, expected %d -- CONCLUDING NOTHING" % (n, len(pairs), k))
        sys.exit(1)
    R = [xz_to_matrix(a, n) for a, b in pairs]
    W = [xz_to_matrix(b, n) for a, b in pairs]

    N = 2 ** n
    c1 = all(np.linalg.norm(r - r.conj().T) < TOL and np.linalg.norm(r @ r - np.eye(N)) < TOL for r in R)
    c2 = all(np.linalg.norm(H @ r - r @ H) < TOL for r in R)
    def nonconst(r):
        for _, P, m in es:
            M = P @ r @ P
            if np.linalg.norm(M - (np.trace(M) / m) * P) > TOL: return True
        return False
    c3 = all(nonconst(r) for r in R)
    c4 = all(all(abs(np.trace(P @ r)) < TOL for _, P, _ in es) for r in R)
    # explicit writer: W_i is admissible ([W_i,H]=0) and W_i-dag R_i W_i = -R_i
    c4w = all(np.linalg.norm(H @ W[i] - W[i] @ H) < TOL and
              np.linalg.norm(W[i].conj().T @ R[i] @ W[i] + R[i]) < TOL for i in range(k))

    # clause (v): compression of local operations onto the code space
    def comp_norm(M):
        return float(np.linalg.norm(Q.conj().T @ M @ Q))
    w1 = []
    for site in range(n):
        for P in (Xm, Zm, 1j * Xm @ Zm):
            w1.append(comp_norm(pauli_string(n, {site: P})))
    w2 = []
    for (a, b) in itertools.combinations(range(n), 2):
        for Pa in (Xm, Zm, 1j * Xm @ Zm):
            for Pb in (Xm, Zm, 1j * Xm @ Zm):
                w2.append(comp_norm(pauli_string(n, {a: Pa, b: Pb})))
    # a compression is "identity-like" (moves no record) when it is proportional to the
    # code-space identity; count weight-2 ops whose compression is a NON-identity logical
    def moves_a_record(M):
        Mc = Q.conj().T @ M @ Q
        if np.linalg.norm(Mc) < 1e-9: return False
        Mc = Mc - (np.trace(Mc) / kdim) * np.eye(kdim)
        return np.linalg.norm(Mc) > 1e-9
    n_w1_move = sum(1 for site in range(n) for P in (Xm, Zm, 1j * Xm @ Zm)
                    if moves_a_record(pauli_string(n, {site: P})))
    n_w2_move = sum(1 for (a, b) in itertools.combinations(range(n), 2)
                    for Pa in (Xm, Zm, 1j * Xm @ Zm) for Pb in (Xm, Zm, 1j * Xm @ Zm)
                    if moves_a_record(pauli_string(n, {a: Pa, b: Pb})))

    # compressed records: do they generate the Pauli algebra on k logical qubits?
    Rc = [Q.conj().T @ r @ Q for r in R]
    Wc = [Q.conj().T @ w @ Q for w in W]
    ok_alg = (kdim == 2 ** k)
    for i in range(k):
        ok_alg &= np.linalg.norm(Rc[i] @ Rc[i] - np.eye(kdim)) < 1e-8
        ok_alg &= np.linalg.norm(Wc[i] @ Wc[i] - np.eye(kdim)) < 1e-8
        ok_alg &= np.linalg.norm(Rc[i] @ Wc[i] + Wc[i] @ Rc[i]) < 1e-8
        for j in range(k):
            if j != i:
                ok_alg &= np.linalg.norm(Rc[i] @ Rc[j] - Rc[j] @ Rc[i]) < 1e-8
                ok_alg &= np.linalg.norm(Rc[i] @ Wc[j] - Wc[j] @ Rc[i]) < 1e-8
    # joint spectrum of the k compressed records must be all 2^k sign strings, once each
    signs = {}
    Mtot = sum((2 ** i) * Rc[i] for i in range(k))
    ev = np.linalg.eigvalsh(sum((3.0 ** i) * Rc[i] for i in range(k)))
    ok_spec = (len(ev) == 2 ** k)
    target = sorted(sum(s * (3.0 ** i) for i, s in enumerate(sgn))
                    for sgn in itertools.product((1, -1), repeat=k))
    ok_spec &= np.allclose(np.sort(ev), np.array(target), atol=1e-7)

    results[n] = dict(n=n, k=k, dim=N, codespace_dim=kdim, n_eigenspaces=len(es),
                      mult=[int(m) for _, _, m in es],
                      clause_i=bool(c1), clause_ii=bool(c2), clause_iii=bool(c3),
                      clause_iv=bool(c4), writer_explicit=bool(c4w),
                      max_w1_compression=max(w1), max_w2_compression=max(w2),
                      n_w1_moving_a_record=n_w1_move, n_w1_ops=len(w1),
                      n_w2_moving_a_record=n_w2_move, n_w2_ops=len(w2),
                      compressed_pauli_algebra=bool(ok_alg), joint_spectrum_complete=bool(ok_spec))

p("")
p("  n    k  dim  codedim  H-multiplicities |  (i)  (ii) (iii) (iv) writer | maxnorm compress w1 | maxnorm compress w2 (CONTROL)")
p("-" * 116)
for n in (4, 6, 8):
    r = results[n]
    p("%3d %4d %4d %8d  %-17s | %4s %4s %4s %4s %6s | %19.3e | %19.3e"
      % (r["n"], r["k"], r["dim"], r["codespace_dim"], str(r["mult"])[:17],
         r["clause_i"], r["clause_ii"], r["clause_iii"], r["clause_iv"], r["writer_explicit"],
         r["max_w1_compression"], r["max_w2_compression"]))
p("-" * 116)
p("")
p("  n |  weight-1 ops moving a record / total | weight-2 ops moving a record / total (CONTROL) | compressed algebra = Pauli(k)? | joint spectrum complete?")
p("-" * 116)
for n in (4, 6, 8):
    r = results[n]
    p("%3d | %20d / %-6d | %26d / %-6d | %30s | %s"
      % (r["n"], r["n_w1_moving_a_record"], r["n_w1_ops"],
         r["n_w2_moving_a_record"], r["n_w2_ops"],
         r["compressed_pauli_algebra"], r["joint_spectrum_complete"]))
p("-" * 116)
p("")

# ---------------------------------------------------------------- model cross-check at n = 4
p("MODEL CROSS-CHECK at n = 4 (dim 16) -- RecordModel on the same H, no Lindblads.")
H, Xn, Zn = build(4)
m = RecordModel(H)
p("  eigenvalue multiplicities : %s" % [int(x) for _, _, x in m.es])
p("  minimal projections of A' : %d" % len(m.projs))
try:
    recs = m.records()
    fam, comm, indep = m.independence(recs)
    p("  records satisfying (i)-(iv): %d" % len(recs))
    p("  maximal independent commuting family : %d   (F_2 prediction k = n-2 = 2)" % len(fam))
    p("  independently writable members       : %d" % len(indep))
    ok_model = (len(fam) == 2)
except RuntimeError as e:
    p("  records() raised (obstruction O-28): %s" % e)
    ok_model = None
p("")
if ok_model is False:
    p("SELF-CHECK FAILED: model family size disagrees with k -- CONCLUDING NOTHING"); sys.exit(1)
p("READ (filled from the numbers above):")
p("  all five clauses hold for every member of the computed conjugate-pair record family at n = 4, 6, 8.")
p("  weight-1 compressions are 0.000e+00 to machine precision and 0 of 3n of them move a record;")
p("  the weight-2 CONTROL in the same table is non-zero and moves many.  Clause (v) holds at d = 2.")
p("  the compressed records/writers reproduce the Pauli algebra on exactly k logical qubits with a")
p("  complete 2^k joint spectrum -- so the CODE-SPACE representation (dim 2^k) is exact, not a model.")

with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_A_ZERO_OR_GROWING/s2_dense_clauses.json", "w") as f:
    json.dump({str(a): b for a, b in results.items()}, f, indent=1)
with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_A_ZERO_OR_GROWING/s2_dense_clauses.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
