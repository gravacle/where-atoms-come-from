"""T0 -- FOUNDATION.  Nothing downstream is allowed to assume any of this.

 (a) EXACT verification of the five record clauses on [[4,2,2]], on products of m such
     blocks, and on [[n,n-2,2]] -- in F_2 / integer arithmetic, no floats.
 (b) NUMERICAL cross-check of (a) against RecordModel on the full carrier at small dim.
 (c) VERIFICATION of the code-space reduction (compress_pauli, and chi computed on the
     reduced space == chi computed on the full carrier).
 (d) VERIFICATION of evolve_times against RecordModel.evolve.

If any self-check here fails, every downstream table is void.
"""
import sys, numpy as np, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM")
from lane_utils import *

np.set_printoptions(precision=6, suppress=True)
FAIL = []
def check(name, ok, extra=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {extra}")
    if not ok: FAIL.append(name)

print("="*100)
print("T0(a)  EXACT CLAUSE VERIFICATION IN F_2 / INTEGER ARITHMETIC  (no floats anywhere in this block)")
print("="*100)

def exact_clauses(stab, n, label, cand):
    """cand: list of (name, (x|z) vector).  Returns per-candidate exact clause verdicts."""
    S, L, pairs = derived_logical_span(stab, n)
    print(f"\n  carrier {label}: n={n} qubits, dim={2**n}, #stab={len(stab)}, "
          f"symplectic_logicals -> {len(pairs)} conjugate pairs (k={len(pairs)})")
    # ---- clause (v) as an EXACT finite check: no single-qubit Pauli is a nontrivial logical
    single_bad = []
    for q in range(n):
        for (xs, zs) in ((( q,), ()), ((), (q,)), ((q,), (q,))):
            v = pauli_vec(n, xs, zs)
            if is_nontrivial_logical(v, S, L, n): single_bad.append((q, xs, zs))
    print(f"    clause (v) EXACT: single-qubit Paulis acting as a nontrivial logical: "
          f"{len(single_bad)} of {3*n}  -> {'PROTECTED' if not single_bad else 'NOT PROTECTED'}")
    rows = []
    for nm, v in cand:
        c_i   = True                                        # i^{x.z} X^x Z^z is Hermitian, squares to I
        c_ii  = all(sp(v, s, n) == 0 for s in stab)         # commutes with every stabiliser => [H,R]=0
        inS   = in_span(v, S, 2*n)
        inNS  = in_span(v, S + L, 2*n)
        c_iii = c_ii and (not inS)          # non-constant on the sector: R P_E is traceless & != c P_E
        c_iv  = (not inS)                   # Tr(P_E R)=0: P_E is a sum of stabilisers; R*s != I for all s in S
        rows.append((nm, c_i, c_ii, c_iii, c_iv, not single_bad, inNS and not inS))
    print(f"    {'record':<26}{'(i)':>6}{'(ii)':>7}{'(iii)':>7}{'(iv)':>6}{'(v)':>6}   derived-in-N(S)\\S")
    for nm, a, b, c, d, e, f in rows:
        print(f"    {nm:<26}{str(a):>6}{str(b):>7}{str(c):>7}{str(d):>6}{str(e):>6}        {f}")
    ok = all(a and b and c and d and e and f for _, a, b, c, d, e, f in rows)
    return ok, rows

# ---- [[4,2,2]] single block
n = 4
cand4 = [("X0X1", pauli_vec(4, (0,1), ())), ("X0X2", pauli_vec(4, (0,2), ())),
         ("Z0Z1", pauli_vec(4, (), (0,1))), ("Z0Z2", pauli_vec(4, (), (0,2))),
         ("Y0Y1", pauli_vec(4, (0,1), (0,1)))]
ok, _ = exact_clauses(stab_nn2(4), 4, "[[4,2,2]]", cand4)
check("exact clauses on [[4,2,2]]", ok)

# ---- products of m [[4,2,2]] blocks
for m in (1, 2, 3, 4, 6, 8):
    n = 4*m
    cand = []
    for b in range(m):
        o = 4*b
        cand.append((f"X{o}X{o+1} (blk{b})", pauli_vec(n, (o, o+1), ())))
        cand.append((f"Z{o}Z{o+1} (blk{b})", pauli_vec(n, (), (o, o+1))))
    # a deliberately BAD candidate: an inter-block X pair, which must FAIL clause (ii)
    if m >= 2:
        cand.append(("X0X4 (INTER-BLOCK, must fail)", pauli_vec(n, (0, 4), ())))
    ok, rows = exact_clauses(stab_blocks(m), n, f"[[4,2,2]]^{m}", cand)
    good = all(r[1] and r[2] and r[3] and r[4] and r[5] and r[6] for r in rows if "INTER" not in r[0])
    bad  = all(not (r[2]) for r in rows if "INTER" in r[0])
    check(f"exact clauses on [[4,2,2]]^{m} (blocks pass, inter-block fails)", good and bad)

# ---- [[n,n-2,2]] family
for n in (4, 6, 8, 10, 12, 14, 16):
    cand = [(f"X0X1", pauli_vec(n, (0,1), ())), (f"X0X{n-1}", pauli_vec(n, (0,n-1), ())),
            (f"Z0Z1", pauli_vec(n, (), (0,1))), (f"Z{n//2}Z{n-1}", pauli_vec(n, (), (n//2, n-1)))]
    ok, _ = exact_clauses(stab_nn2(n), n, f"[[{n},{n-2},2]]", cand)
    check(f"exact clauses on [[{n},{n-2},2]]", ok)

print()
print("="*100)
print("T0(b)  NUMERICAL CROSS-CHECK AGAINST RecordModel ON THE FULL CARRIER")
print("="*100)
for n in (4, 6):
    Xall = xz_to_matrix([1]*n+[0]*n, n); Zall = xz_to_matrix([0]*n+[1]*n, n)
    H = -(Xall + Zall)
    M = RecordModel(H)
    es = M.es
    print(f"\n  [[{n},{n-2},2]] dim {2**n}: eigenvalue multiplicities {[m for _,_,m in es]}")
    for nm, v in [("X0X1", pauli_vec(n,(0,1),())), ("Z0Z1", pauli_vec(n,(),(0,1)))]:
        R = xz_to_matrix(v, n)
        c1 = np.linalg.norm(R - R.conj().T) < 1e-12 and np.linalg.norm(R@R - np.eye(2**n)) < 1e-12
        c2 = np.linalg.norm(H@R - R@H) < 1e-12
        from record_model import clause_iii, clause_iv
        c3 = clause_iii(R, es); c4 = clause_iv(R, es)
        print(f"    {nm}: (i)={c1} (ii)={c2} (iii)={c3} (iv)={c4}  "
              f"max|Tr(P_E R)|={max(abs(np.trace(P@R)) for _,P,_ in es):.3e}")
        check(f"numeric clauses {nm} on [[{n},{n-2},2]]", c1 and c2 and c3 and c4)

print()
print("="*100)
print("T0(c)  CODE-SPACE REDUCTION IS EXACT")
print("="*100)
for n in (6, 8, 10):
    reps, idx = code_reps(n)
    d = len(reps)
    assert d == 2**(n-2), (n, d)
    # explicit isometry into the full space
    V = np.zeros((2**n, d), dtype=complex)
    full = (1 << n) - 1
    for a, v in enumerate(reps):
        V[v, a] = 1/np.sqrt(2); V[(~v) & full, a] = 1/np.sqrt(2)
    check(f"n={n}: code basis orthonormal, dim {d} == 2^(n-2)",
          np.linalg.norm(V.conj().T@V - np.eye(d)) < 1e-12)
    Xall = xz_to_matrix([1]*n+[0]*n, n); Zall = xz_to_matrix([0]*n+[1]*n, n)
    H = -(Xall+Zall)
    check(f"n={n}: code space is the H ground space (energy -2)",
          np.linalg.norm(H@V + 2*V) < 1e-10)
    tests = [("X0X1", pauli_vec(n,(0,1),()), compress_XX(0,1,n,reps,idx)),
             (f"X1X{n-1}", pauli_vec(n,(1,n-1),()), compress_XX(1,n-1,n,reps,idx)),
             ("Z0Z1", pauli_vec(n,(),(0,1)), compress_ZZ(0,1,n,reps,idx)),
             (f"Z2Z{n-2}", pauli_vec(n,(),(2,n-2)), compress_ZZ(2,n-2,n,reps,idx))]
    for nm, v, Mc in tests:
        P = xz_to_matrix(v, n)
        W = V.conj().T @ P @ V
        check(f"n={n}: compress({nm}) matches V-dag P V", np.linalg.norm(W - Mc) < 1e-10,
              f"  max dev {np.abs(W-Mc).max():.2e}")
        check(f"n={n}: compress({nm}) is an exact integer matrix",
              np.allclose(Mc, np.round(np.real(Mc))))
    # general compress_pauli including Y-type
    for v, nm in [(pauli_vec(n,(0,1),(0,1)), "Y0Y1"), (pauli_vec(n,(0,1),(2,3)), "X0X1Z2Z3")]:
        P = xz_to_matrix(v, n); W = V.conj().T @ P @ V
        Mc = compress_pauli(v, n, reps, idx)
        check(f"n={n}: compress_pauli({nm}) matches", np.linalg.norm(W - Mc) < 1e-10,
              f"  max dev {np.abs(W-Mc).max():.2e}")

print()
print("="*100)
print("T0(d)  evolve_times MATCHES RecordModel.evolve, AND REDUCED chi == FULL chi")
print("="*100)
n = 6
reps, idx = code_reps(n); d = len(reps)
full = (1 << n) - 1
V = np.zeros((2**n, d), dtype=complex)
for a, v in enumerate(reps):
    V[v, a] = 1/np.sqrt(2); V[(~v) & full, a] = 1/np.sqrt(2)
Xall = xz_to_matrix([1]*n+[0]*n, n); Zall = xz_to_matrix([0]*n+[1]*n, n)
Hfull = -(Xall+Zall)
Mfull = RecordModel(Hfull)
env = Environment(nq=3, energies=(1.0, 1.4, 0.7), beta=2.0)
RA_f = xz_to_matrix(pauli_vec(n,(0,1),()), n)
RB_f = xz_to_matrix(pauli_vec(n,(),(2,3)), n)
RA_r = compress_XX(0,1,n,reps,idx).astype(complex)
RB_r = compress_ZZ(2,3,n,reps,idx).astype(complex)
Hred = -2.0*np.eye(d, dtype=complex)
st_r = np.eye(d, dtype=complex)/d
for lam in (0.4, 0.8):
    for t in (2.0, 5.0):
        rf = Mfull.evolve([(RA_f, 0), (RB_f, 0)], env, lam=lam, t=t)
        cf = env.holevo(rf, RA_f, 2**n)
        rr = next(evolve_times(Hred, env, [(RA_r, 0), (RB_r, 0)], lam, [t], st_r))
        cr = env.holevo(rr, RA_r, d)
        check(f"lam={lam} t={t}: reduced chi == full chi", abs(cf-cr) < 1e-9,
              f"  full {cf:.12f}  reduced {cr:.12f}  diff {abs(cf-cr):.2e}")

print()
print("="*100)
print("SELF-CHECK SUMMARY:", "ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES: {FAIL}")
print("="*100)
