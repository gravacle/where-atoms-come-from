"""ADDENDUM: which ingredient produces the continuum? Candidates: (a) metastability / finite
lifetime of the record (T-30's stated mechanism: 'a slow mode, not an exact Pauli'), (b) a
splitting dE between the record's two states, (c) a thermal state instead of the uniform trace.

Test: a CLOSED system -- no dissipation, no Lindblad operators, INFINITE barrier, the record is an
EXACT Pauli with an exact infinite lifetime -- with a splitting dE, in a Gibbs state at T.
If the trace ratio is continuous here, metastability contributes NOTHING; the continuum is
(b)+(c), i.e. leaving the degeneracy (W-asymmetry) and T->0 axes of the corner, and T-30's
'slow mode, not an exact Pauli' mechanism is not the operative one even in its own numbers
(its record operator sz IS an exact Pauli)."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'model'))
import grounded as G

I2 = np.eye(2); X = np.array([[0,1],[1,0]], dtype=complex); Z = np.array([[1,0],[0,-1]], dtype=complex)
def word(n, s):
    M = np.array([[1]], dtype=complex)
    for c in s: M = np.kron(M, {'I':I2,'X':X,'Z':Z}[c])
    return M

n = 4
Hs = -(word(n,'X'*n) + word(n,'Z'*n))
R = word(n,'ZZII')                       # exact signed Pauli, commutes with Hs: [R,Hs]=0
print("closed [[4,2,2]] carrier + splitting eps*R (R the record itself, still exact, [R,H]=0,")
print("infinite lifetime, no bath). Gibbs expectation <R> at beta:")
print(f"{'eps':>10}{'beta':>8}{'<R>_Gibbs':>16}{'in Z[i]?':>10}")
for eps in (0.0, 1e-3, 1e-1):
    H = Hs + eps*R
    w, V = np.linalg.eigh(H)
    for beta in (1.0, 10.0, 100.0):
        p = np.exp(-beta*(w - w.min())); p = p/p.sum()
        rho = (V * p) @ V.conj().T
        r = float(np.real(np.trace(rho @ R)))
        print(f"{eps:>10.0e}{beta:>8.0f}{r:>16.9f}{str(abs(r-round(r))<1e-9):>10}")
print()
print("commutator check ||[R,Hs]|| =", np.linalg.norm(R@Hs - Hs@R))
print()
print("-> non-integer values here come from an EXACT record on an EXACT carrier with NO dissipation:")
print("   the continuum is produced by (splitting + thermal state), not by metastability or by the")
print("   record failing to be an exact Pauli. T-30's control (uniform code-space trace, degenerate")
print("   shells) differs from its measurement in the STATE, and the state does all the work.")
