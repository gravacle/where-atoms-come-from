"""S1 -- BUILD THE CARRIER FAMILY AND VERIFY THE FIVE CLAUSES ON EVERY RECORD.

D-18: a record is a logical operator; nothing is called a record until the clauses are checked
on the carrier it lives on.  Clauses (i)-(iv) are checked TWO WAYS:
   MATRIX  -- record_model.eigenspaces / clause_iii / clause_iv on the explicit 2^n matrices
              (feasible to n = 8, dim 256; n = 10 -> dim 1024 also done)
   F2      -- the same statements done in the symplectic algebra, valid at every n
The F2 route is validated by agreeing with the MATRIX route wherever both run.
Clause (v) is CARRIER DATA and is checked against single-qubit regions using distance 2.
"""
import sys, io, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY")
from common import *

OUT = io.StringIO()
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); OUT.write(s + "\n")

P("=" * 108)
P("S1  CARRIER FAMILY [[n, n-2, 2]] AND CLAUSE VERIFICATION")
P("=" * 108)

NS = [4, 6, 8, 10, 12, 14, 16]
MATRIX_MAX_N = 10          # dim 1024; n=12 -> dim 4096 matrices are 268 MB each, structural only

rows = []
for n in NS:
    car = carrier(n)
    k = car['k']
    ch = car['checks']
    # ---------------- F2 clause checks (valid at every n) ----------------
    stab = car['stab']
    f2 = dict(
        i_bit=True,                                   # kron of I/X/Y/Z: Hermitian, squares to I
        ii_durable=ch['commute_with_stab'],
        iii_nontrivial=ch['not_a_stabiliser'],
        iv_writable=ch['not_a_stabiliser'],           # Tr(P_E R) = 0 unless R in {I,Sx,Sz,SxSz}
    )
    # ---------------- MATRIX clause checks ----------------
    mat = dict(i=None, ii=None, iii=None, iv=None, gs_dim=None, sector_dims_ok=None)
    if n <= MATRIX_MAX_N:
        H = code_hamiltonian(n)
        Rs = record_matrices(car)
        es = eigenspaces(H)
        mat['mults'] = [m for _, _, m in es]
        ok_i = all(np.linalg.norm(R - R.conj().T) < 1e-9 and
                   np.linalg.norm(R @ R - np.eye(2 ** n)) < 1e-9 for R in Rs)
        ok_ii = all(np.linalg.norm(H @ R - R @ H) < 1e-9 for R in Rs)
        ok_iii = all(clause_iii(R, es) for R in Rs)
        ok_iv = all(clause_iv(R, es) for R in Rs)
        mat.update(i=ok_i, ii=ok_ii, iii=ok_iii, iv=ok_iv)
        # SELF-CHECK: the k records must split the ground space into 2^k sectors of EQUAL
        # dimension 1.  sum_r dim_r must equal dim(ground space).  If this fails the uniform
        # p_r = 2^-k assumed by the fast path is wrong and NO conclusion may be drawn.
        w, V = np.linalg.eigh(H)
        g = V[:, np.abs(w - w[0]) < 1e-9]
        mat['gs_dim'] = g.shape[1]
        dims = {}
        blocks = {(): g}
        for R in Rs:
            nb = {}
            for lab, C in blocks.items():
                Rs_ = C.conj().T @ R @ C
                ws, Vs = np.linalg.eigh(Rs_)
                for s in (+1, -1):
                    idx = [q for q in range(len(ws)) if (ws[q] > 0) == (s > 0)]
                    if idx: nb[lab + (s,)] = C @ Vs[:, idx]
            blocks = nb
        dims = sorted(set(C.shape[1] for C in blocks.values()))
        mat['n_sectors'] = len(blocks)
        mat['sector_dims'] = dims
        mat['sector_dims_ok'] = (len(blocks) == 2 ** k and dims == [1] and
                                 sum(C.shape[1] for C in blocks.values()) == g.shape[1])
        del H, Rs, g, blocks
    rows.append((n, k, car, f2, mat))

P("")
P(f"{'n':>3} {'k':>3} {'dim':>6} | {'pairs':>5} {'symp=I':>7} {'Zi commute':>10} | "
  f"{'F2 (i)':>7} {'F2 (ii)':>8} {'F2(iii)':>8} {'F2 (iv)':>8} | "
  f"{'M (i)':>6} {'M (ii)':>7} {'M(iii)':>7} {'M (iv)':>7} | {'gs dim':>7} {'sectors':>8} {'dims':>6} {'OK':>4}")
P("-" * 148)
for n, k, car, f2, mat in rows:
    ch = car['checks']
    m = lambda v: ("--" if v is None else ("yes" if v else "NO"))
    P(f"{n:>3} {k:>3} {2**n:>6} | {ch['n_pairs']:>5} {str(ch['pairing_is_identity']):>7} "
      f"{str(ch['records_commute']):>10} | "
      f"{str(f2['i_bit']):>7} {str(f2['ii_durable']):>8} {str(f2['iii_nontrivial']):>8} {str(f2['iv_writable']):>8} | "
      f"{m(mat['i']):>6} {m(mat['ii']):>7} {m(mat['iii']):>7} {m(mat['iv']):>7} | "
      f"{str(mat['gs_dim']):>7} {str(mat.get('n_sectors')):>8} {str(mat.get('sector_dims')):>6} "
      f"{m(mat['sector_dims_ok']):>4}")

P("")
P("EIGENVALUE MULTIPLICITIES OF H = -(X^n + Z^n)  (matrix route only):")
for n, k, car, f2, mat in rows:
    if 'mults' in mat:
        P(f"   n={n:>2}  {mat['mults']}   expected [2^k, 2^(k+1), 2^k] = "
          f"[{2**(n-2)}, {2**(n-1)}, {2**(n-2)}]")

# ---------------- clause (v) : PROTECTED against a single-qubit region ----------------
P("")
P("CLAUSE (v) PROTECTED -- carrier data, checked explicitly against SINGLE-QUBIT regions.")
P("A single-qubit operation is a combination of I,X,Y,Z on one site.  It can flip a record only")
P("if some single-qubit Pauli both (a) commutes with BOTH stabilisers, so it is admissible, and")
P("(b) anticommutes with that record.  Enumerated exhaustively:")
P("")
P(f"{'n':>3} {'k':>3} | {'single-qubit Paulis':>19} {'admissible ones':>16} {'that flip a record':>19} | READ")
P("-" * 90)
for n, k, car, f2, mat in rows:
    stab = car['stab']
    tot = adm = flip = 0
    for site in range(n):
        for (x, z) in ((1, 0), (0, 1), (1, 1)):
            v = [0] * (2 * n); v[site] = x; v[n + site] = z
            tot += 1
            if all(sp_form(v, s, n) == 0 for s in stab):
                adm += 1
                if any(sp_form(v, r, n) == 1 for r in car['recs_xz']): flip += 1
    P(f"{n:>3} {k:>3} | {tot:>19} {adm:>16} {flip:>19} | "
      f"{'PROTECTED' if flip == 0 else 'NOT PROTECTED'}")

# POSITIVE CONTROL for clause (v) (D-15): the same enumeration over TWO-qubit regions, where
# the effect CAN occur -- distance is 2, so weight-2 logicals exist and must be found.
P("")
P("POSITIVE CONTROL for the clause (v) test (D-15) -- the identical enumeration over TWO-qubit")
P("regions, where a flip CAN occur because the code distance is exactly 2:")
P("")
P(f"{'n':>3} {'k':>3} | {'two-qubit Paulis':>17} {'admissible ones':>16} {'that flip a record':>19} | READ")
P("-" * 90)
import itertools
for n, k, car, f2, mat in rows:
    stab = car['stab']
    tot = adm = flip = 0
    for s1, s2 in itertools.combinations(range(n), 2):
        for p1 in ((1, 0), (0, 1), (1, 1)):
            for p2 in ((1, 0), (0, 1), (1, 1)):
                v = [0] * (2 * n)
                v[s1], v[n + s1] = p1
                v[s2], v[n + s2] = p2
                tot += 1
                if all(sp_form(v, s, n) == 0 for s in stab):
                    adm += 1
                    if any(sp_form(v, r, n) == 1 for r in car['recs_xz']): flip += 1
    P(f"{n:>3} {k:>3} | {tot:>17} {adm:>16} {flip:>19} | "
      f"{'flips found' if flip else 'NONE -- control DEAD'}")

P("")
allok = all(r[3]['ii_durable'] and r[3]['iii_nontrivial'] and r[3]['iv_writable'] for r in rows) \
        and all(r[2]['checks']['pairing_is_identity'] and r[2]['checks']['records_commute'] for r in rows) \
        and all((r[4]['sector_dims_ok'] is None) or r[4]['sector_dims_ok'] for r in rows) \
        and all((r[4]['i'] is None) or (r[4]['i'] and r[4]['ii'] and r[4]['iii'] and r[4]['iv']) for r in rows)
P("READ: every self-check above " + ("PASSED" if allok else "DID NOT PASS -- DRAW NO CONCLUSION"))
P("      The k = n-2 operators Z_1..Z_k are records on [[n,n-2,2]] by clauses (i)-(iv), are")
P("      protected against every single-qubit region, and split the ground space into 2^k")
P("      one-dimensional joint sectors -- so the fast path's uniform prior p_r = 2^-k is exact.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY/s1_records.txt", "w").write(OUT.getvalue())
