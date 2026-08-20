"""ADVERSARIAL CHECKS of T-28.

1. EXACT commutant dimensions by integer arithmetic (Fraction Gaussian elimination), fully
   independent of numpy's eig and of any tolerance: multiplicity of each exact integer eigenvalue
   E via rank(H - E I) over Q, then sum m_E^2. Must give 96 / 1536 / 24.
2. Does part 1 actually DISCRIMINATE the |lambda| filter from |Re lambda|? Recompute the t_m=inf
   slow-mode count under a |Re lambda| filter: with no dissipation it must keep everything (n^2),
   which part 1 would have caught as a mismatch. Confirms the test is discriminating.
3. INCONSISTENCY probe: grounded.clause_ii computes rate = |Re <R, L_ad R>|, so a ROTATING record
   (the very case C-75 says is NOT durable) gets rate 0 and durable=True, contradicting
   slow_modes' |lambda| filter and C-75's own sentence. Demonstrate on H = sz, R = sx."""
import sys, os, numpy as np
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'model'))
import grounded as G

I2 = np.eye(2); X = np.array([[0,1],[1,0]], dtype=complex); Z = np.array([[1,0],[0,-1]], dtype=complex)
def word(n, s):
    M = np.array([[1]], dtype=complex)
    for c in s: M = np.kron(M, {'I':I2,'X':X,'Z':Z}[c])
    return M

def rank_exact(A):
    """Rank over Q of an integer matrix, exact Gaussian elimination."""
    A = [[Fraction(int(round(x.real))) for x in row] for row in A]
    m, n = len(A), len(A[0]); r = 0
    for c in range(n):
        piv = next((i for i in range(r, m) if A[i][c] != 0), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        for i in range(m):
            if i != r and A[i][c] != 0:
                f = A[i][c]/pv
                A[i] = [a - f*b for a, b in zip(A[i], A[r])]
        r += 1
    return r

CARRIERS = [("[[4,2,2]]", -(word(4,'XXXX')+word(4,'ZZZZ')), [-2,0,2]),
            ("[[6,4,2]]", -(word(6,'X'*6)+word(6,'Z'*6)), [-2,0,2]),
            ("3-qubit Ising", -(word(3,'ZZI')+word(3,'IZZ')), [-2,0,2])]
print("1. exact commutant dimension sum m_E^2 over Q (no floating point, no tolerance):")
for name, H, evs in CARRIERS:
    d = H.shape[0]
    Hi = np.real(H)
    assert np.allclose(Hi, np.round(Hi)) and np.allclose(np.imag(H), 0)
    total = 0; mults = []
    for E in evs:
        m = d - rank_exact(Hi - E*np.eye(d))
        mults.append(m); total += m*m
    assert sum(mults) == d, (name, mults)
    print(f"   {name:<16} multiplicities {mults}  sum m^2 = {total}")
print()
print("2. would a |Re lambda| filter have been caught by part 1? (no dissipation, t_m -> inf)")
for name, H, _ in [CARRIERS[0], CARRIERS[2]]:
    n2 = H.shape[0]**2
    L = G.liouvillian(H, []).conj().T
    w = np.linalg.eigvals(L)
    keep_re = sum(1 for x in w if abs(x.real) <= 1e-9*max(1, np.max(np.abs(w))))
    keep_abs = sum(1 for x in w if abs(x) <= 1e-9*max(1, np.max(np.abs(w))))
    print(f"   {name:<16} dim {n2}: |Re| filter keeps {keep_re} (would MISMATCH), |.| filter keeps {keep_abs}")
print()
print("3. clause_ii vs slow_modes on a ROTATING record (H = sz, R = sx, no dissipation, t_m = 10):")
sz = np.array([[1,0],[0,-1]], dtype=complex); sx = X
c = G.clause_ii(sz, [], sx, 10.0)
rates, obs = G.slow_modes(sz, [], 10.0)
print(f"   clause_ii: rate = {c['rate']:.3e}, durable = {c['durable']}   (sx rotates at omega = 2)")
print(f"   slow_modes at t_m=10 keeps {len(rates)} modes (the commutant of sz is 2-dim: sx excluded)")
print("   -> clause_ii calls the rotating record DURABLE; slow_modes and C-75 say it is not.")
