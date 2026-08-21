"""S1 -- D-18: NOTHING IS CALLED A RECORD UNTIL THE FIVE CLAUSES ARE CHECKED ON ITS CARRIER.

Verifies clauses (i)-(iv) NUMERICALLY, on the full Hilbert space, for every record of both
families up to the largest dimension that fits, and clause (v) EXACTLY over F_2 with a
POSITIVE CONTROL in the same table (region size 1 must give ZERO flips; region size 2 must
give a NON-ZERO count -- a zero with no positive control beside it discriminates nothing).

It also cross-validates the SYMBOLIC claim used later at large k -- that the k records split
every eigenspace of H into 2^k blocks of EQUAL dimension, so every record is independently
writable -- against RecordModel.independently_writable() wherever that is affordable.
"""
import sys, time
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD")
from record_model import (RecordModel, eigenspaces, clause_iii, clause_iv, xz_to_matrix)
import carriers as C

OUT = []
def P(s=""):
    print(s); OUT.append(s)

def H_of(car):
    n = car["n"]
    return -sum(xz_to_matrix(s, n) for s in car["stabs"])

def local_logical_dim(car, region):
    """dim_F2 of the group of LOGICAL CLASSES implementable by a Pauli supported in `region`.
       = dim{Paulis in region that are in N(S)}  -  dim{stabilisers supported in region}.
       Basis-INDEPENDENT: it is a property of the region and the code, not of any chosen
       generating set of records."""
    n = car["n"]
    idx = sorted(region)
    # coordinates: for each qubit in region, an x bit and a z bit
    cols = []
    for q in idx:
        e = [0] * (2 * n); e[q] = 1; cols.append(e)
        e = [0] * (2 * n); e[n + q] = 1; cols.append(e)
    # in N(S): symplectic form with every stabiliser vanishes
    A = [[C.sp(c, s, n) for c in cols] for s in car["stabs"]]
    nul = _nullity(A, len(cols))
    # stabilisers supported inside the region
    ns = _rank([s for s in car["stabs"] if C.support(s, n) <= set(idx)], 2 * n)
    return nul - ns

def _rank(rows, ncol):
    rows = [r[:] for r in rows]; r = 0
    for c in range(ncol):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x + y) % 2 for x, y in zip(rows[i], rows[r])]
        r += 1
    return r

def _nullity(A, ncol):
    if not A: return ncol
    return ncol - _rank(A, ncol)

def clause_v_scan(car, fam, max_r=3):
    """CLAUSE (v).  For each region size r, the number of records in `fam` that some
       ADMISSIBLE (= commutes with H, i.e. lies in N(S)) operation supported on a SINGLE
       contiguous region of size r can flip.  r = 1 is the claim; r = 2 is the CONTROL."""
    n = car["n"]
    rows = []
    for r in range(1, max_r + 1):
        best = 0; total_reach = 0
        for start in range(n - r + 1):
            reg = set(range(start, start + r))
            d = local_logical_dim(car, reg)
            total_reach = max(total_reach, d)
            # how many records of fam can be flipped by SOME operator in that region
            flipped = 0
            for R in fam:
                # exists P supported in reg, P in N(S), with sp(P,R)=1 ?
                cols = []
                for q in sorted(reg):
                    e = [0] * (2 * n); e[q] = 1; cols.append(e)
                    e = [0] * (2 * n); e[n + q] = 1; cols.append(e)
                Anorm = [[C.sp(c, s, n) for c in cols] for s in car["stabs"]]
                # solvable iff sp(.,R) is NOT in the row space of the N(S) constraints
                target = [C.sp(c, R, n) for c in cols]
                if _rank(Anorm, len(cols)) < _rank(Anorm + [target], len(cols)):
                    flipped += 1
            best = max(best, flipped)
        rows.append((r, best, total_reach))
    return rows

def run(car, do_numeric=True, do_model=False):
    n, dim = car["n"], car["dim"]
    fam, part, checks = C.records_of(car)
    k = len(fam)
    P("-" * 100)
    P("CARRIER %-14s n=%-3d dim=%-6d k=%d" % (car["label"], n, dim, k))
    P("  symbolic self-checks: " + "  ".join("%s=%s" % (a, b) for a, b in checks.items()))
    if not C.all_checks_pass(checks):
        P("  SELF-CHECK FAILED -- NO CONCLUSION DRAWN FOR THIS CARRIER"); return None
    res = dict(car=car["label"], n=n, dim=dim, k=k, checks=checks)
    if do_numeric:
        t0 = time.time()
        H = H_of(car)
        es = eigenspaces(H)
        mults = [m for _, _, m in es]
        P("  H = -(sum of stabilisers); eigenvalue multiplicities %s ; sum=%d (= dim %s)"
          % (mults, sum(mults), "OK" if sum(mults) == dim else "MISMATCH"))
        rows = []
        for i, v in enumerate(fam):
            # T-35: stderr progress marker removed -- reproduce.sh captures 2>&1, so it appeared as extra lines vs the sealed output
            R = xz_to_matrix(v, n)
            c_i = (np.linalg.norm(R - R.conj().T) < 1e-9 and
                   np.linalg.norm(R @ R - np.eye(dim)) < 1e-9)
            c_ii = np.linalg.norm(H @ R - R @ H) < 1e-8
            c_iii = clause_iii(R, es)
            c_iv = clause_iv(R, es)
            rows.append((i, C.weight(v, n), c_i, c_ii, c_iii, c_iv))
        allok = all(r[2] and r[3] and r[4] and r[5] for r in rows)
        P("  clauses (i)-(iv) verified on the FULL Hilbert space for all %d records: %s"
          % (k, "ALL PASS" if allok else "FAILURE"))
        for r in rows[:4] + ([("...",)] if len(rows) > 4 else []):
            P("     " + str(r))
        res["numeric_clauses_ok"] = bool(allok)
        # joint-eigenspace equality (the symbolic claim used at large k)
        pred = dim // (2 ** (len(car["stabs"]) + k))
        P("  predicted joint eigenspace dim 2^n/2^(rank S + k) = %d" % pred)
        if do_model:
            t1 = time.time()
            M = RecordModel(H)
            Rs = [xz_to_matrix(v, n) for v in fam]
            iw = M.independently_writable(Rs)
            P("  RecordModel.independently_writable -> %d of %d records  (%.1fs)"
              % (len(iw), k, time.time() - t1))
            res["indep_writable"] = len(iw)
            res["n_minimal_projections"] = len(M.projs)
        P("  numeric block: %.1fs" % (time.time() - t0))
    # clause (v) with its positive control IN THE SAME TABLE
    # T-35: stderr progress marker removed -- reproduce.sh captures 2>&1, so it appeared as an extra line vs the sealed output
    P("  CLAUSE (v) region scan   [r=1 is the claim, r>=2 is the POSITIVE CONTROL]")
    P("     region_size  records_flippable_by_one_region   dim_of_local_logical_group")
    for r, f, d in clause_v_scan(car, fam, max_r=min(4, n - 1)):
        P("     %-12d %-34d %d" % (r, f, d))
    res["clause_v"] = clause_v_scan(car, fam, max_r=min(4, n - 1))
    return res

P("=" * 100)
P("S1  CLAUSE VERIFICATION -- every object called a record below is checked here (D-18)")
P("=" * 100)
results = []
for n in (4, 6, 8, 10):
    results.append(run(C.family_A(n), do_numeric=True, do_model=(2 ** n <= 256)))
for m in (1, 2):
    results.append(run(C.family_B(m), do_numeric=True, do_model=(4 ** (2 * m) <= 256)))
# symbolic-only reach check (no Hilbert space built)
P("-" * 100)
P("SYMBOLIC-ONLY reach (no Hilbert space built; F_2 only):")
for n in (16, 24, 32, 48, 64):
    car = C.family_A(n); fam, part, ch = C.records_of(car)
    P("  A n=%-3d k=%-3d dim=2^%-3d self-checks %s"
      % (n, len(fam), n, "PASS" if C.all_checks_pass(ch) else "FAIL"))
for m in (8, 16, 24, 32):
    car = C.family_B(m); fam, part, ch = C.records_of(car)
    P("  B m=%-3d k=%-3d dim=2^%-3d self-checks %s"
      % (m, len(fam), car["n"], "PASS" if C.all_checks_pass(ch) else "FAIL"))

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_D_THRESHOLD/s1_verify_clauses.txt",
     "w").write("\n".join(OUT) + "\n")
