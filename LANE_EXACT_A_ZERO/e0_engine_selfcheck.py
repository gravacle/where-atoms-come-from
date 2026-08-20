"""E-0  SELF-CHECK OF THE EXACT PAULI ENGINE.

Nothing downstream is trustworthy unless the symbolic Z_4 algebra reproduces, EXACTLY, the matrices
the program's own model builds.  Three independent checks, all exact:

  (1) W(a) built symbolically as a Gaussian-integer matrix == record_model.xz_to_matrix(a) (which is
      float, so this comparison is the ONE place a float appears; the tolerance is 0, we require
      bit-identical integers after rounding and we check the rounding was exact).
  (2) W(a)W(b) == i^{phi(a,b)} W(a xor b)  as exact Gaussian-integer matrices, over ALL a,b at n=2,3.
  (3) [W(a),W(b)] == 0 as an exact matrix  <=>  sp(a,b) == 0.   Exhaustive at n=2,3.
  (4) phi(a,a) == 0 (Hermitian involution), W(a)^dag == W(a), exhaustive.

If ANY of these fails, the run reports FAIL and concludes nothing.
"""
import sys
from itertools import product as iproduct

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_A_ZERO"
sys.path.insert(0, LANE)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")

from exact_pauli import (phi, sp, pmul, to_gaussian_matrix, gmatmul, gsub, gis_zero, gfrob2)
import numpy as np
from record_model import xz_to_matrix

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s)
    print(s); sys.stdout.flush()


def allvecs(n):
    return [list(v) for v in iproduct([0, 1], repeat=2 * n)]


say("=" * 110)
say("E-0  EXACT PAULI ENGINE SELF-CHECK   (integers only; the single float appears in check 1)")
say("=" * 110)

ok_all = True

# ---------------- check 1: symbolic matrix == model's xz_to_matrix
for n in (1, 2, 3):
    worst = 0.0
    nonint = 0
    for a in allvecs(n):
        G = to_gaussian_matrix((0, a), n)
        M = xz_to_matrix(a, n)
        d = len(M)
        for i in range(d):
            for j in range(d):
                re, im = G[i][j]
                worst = max(worst, abs(M[i][j].real - re), abs(M[i][j].imag - im))
                if abs(M[i][j].real - round(M[i][j].real)) > 0 or abs(M[i][j].imag - round(M[i][j].imag)) > 0:
                    nonint += 1
    ok = (worst == 0.0 and nonint == 0)
    ok_all &= ok
    say("  check 1  n=%d  vectors=%4d   max |symbolic - model| = %r   non-integer model entries = %d   -> %s"
        % (n, 4 ** n, worst, nonint, "PASS" if ok else "FAIL"))

# ---------------- check 2 & 3 & 4: exhaustive exact algebra
for n in (1, 2, 3):
    V = allvecs(n)
    bad_mul = bad_comm = bad_herm = 0
    for a in V:
        # check 4
        if phi(a, a, n) != 0:
            bad_herm += 1
        for b in V:
            A = to_gaussian_matrix((0, a), n)
            B = to_gaussian_matrix((0, b), n)
            AB = gmatmul(A, B)
            m, c = pmul((0, a), (0, b), n)
            C = to_gaussian_matrix((m, c), n)
            if not gis_zero(gsub(AB, C)):
                bad_mul += 1
            BA = gmatmul(B, A)
            comm_zero = gis_zero(gsub(AB, BA))
            if comm_zero != (sp(a, b, n) == 0):
                bad_comm += 1
    ok = (bad_mul == 0 and bad_comm == 0 and bad_herm == 0)
    ok_all &= ok
    say("  checks 2-4  n=%d  pairs=%6d   product mismatches=%d   commutator/sp mismatches=%d   phi(a,a)!=0 count=%d   -> %s"
        % (n, len(V) ** 2, bad_mul, bad_comm, bad_herm, "PASS" if ok else "FAIL"))

# ---------------- check 5: exact Frobenius norm of a non-zero commutator
say("")
say("  check 5  EXACT commutator magnitude (no floats): [W(a),W(b)] = 2 i^phi(a,b) W(a xor b) when sp=1")
for n in (1, 2, 3):
    V = allvecs(n)
    vals = set()
    for a in V:
        for b in V:
            if sp(a, b, n) == 1:
                A = to_gaussian_matrix((0, a), n)
                B = to_gaussian_matrix((0, b), n)
                vals.add(gfrob2(gsub(gmatmul(A, B), gmatmul(B, A))))
    say("           n=%d  ||[W(a),W(b)]||_F^2 over ALL anticommuting pairs = %s   (predicted 4*2^n = %d)"
        % (n, sorted(vals), 4 * 2 ** n))
    ok_all &= (vals == {4 * 2 ** n})

say("")
say("  ENGINE SELF-CHECK: %s" % ("ALL PASS" if ok_all else "**FAIL -- NOTHING DOWNSTREAM IS VALID**"))
say("=" * 110)

with open(LANE + "/e0_engine_selfcheck.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
sys.exit(0 if ok_all else 1)
