"""O5 ADVERSARIAL D.
 (d) is the dressed record R_c itself 'a function of H at resolution w'?  (the mirror of the
     lane's own P-2 withdrawal argument, applied to the RECORD instead of the WRITER)
 (e) the ledger claims 'local flip amplitude <= O(p^(d-1))'.  Tested on d=2 only (gives O(p)).
     Test it on d=3 (Steane) where it predicts O(p^2).  If it comes out O(p) the claim is wrong.
"""
import numpy as np, itertools
from o5_common import Zop, Xop, Yop, toric_H, sym_H, local_perturbation, Z_A_SUP, X_A_SUP, NQ

V = local_perturbation(seed=2026)

print("=" * 100)
print("(d) IS THE DRESSED RECORD ITSELF A FUNCTION OF H AT THE TOLERANCE THAT MAKES IT APPROXIMATE?")
print("    The lane withdrew P-2 with: 'on a cluster of width delta EVERY operator commutes with H")
print("    to within delta'.  The same sentence applies to R_c.  Measure how far R_c is from the")
print("    nearest exact function of H on the cluster (a +-1 diagonal in the energy eigenbasis).")
print("=" * 100)


def block(H0, R0, g, p):
    e, U = np.linalg.eigh(H0 + p * V)
    Uc = U[:, :g]
    M = Uc.conj().T @ R0 @ Uc
    M = (M + M.conj().T) / 2
    w, Q = np.linalg.eigh(M)
    Rc = Q @ np.diag(np.sign(w)) @ Q.conj().T
    return np.diag(e[:g]).astype(complex), Rc, e[g - 1] - e[0]


for nm, H0, R0, g in (("TORIC 2x2", toric_H(), Zop(Z_A_SUP), 4), ("SYMMETRY", sym_H(), Zop([0]), 2)):
    for p in (1e-3, 1e-2):
        Hc, Rc, w_ = block(H0, R0, g, p)
        best = min(np.linalg.norm(Rc - np.diag(np.array(s, float)).astype(complex), 2)
                   for s in itertools.product([1.0, -1.0], repeat=g))
        offd = np.linalg.norm(Rc - np.diag(np.diag(Rc)), 2)
        print(f"   {nm:>10s} p={p:.0e}  dist(R_c, nearest fn of H) = {best:.4f}   "
              f"||offdiag(R_c)|| = {offd:.4f}   cluster width {w_:.3e}")
print("""   R_c is O(1) away from every function of H -- so in NORM it is genuinely not f(H).  Good for
   the lane.  BUT the definitional test on offer is ||[H,R]|| <= eps, and by THAT test the exact
   function of H scores 0 and beats R_c.  The test cannot express what R_c has and f(H) lacks.""")

print()
print("=" * 100)
print("(e) 'LOCAL FLIP AMPLITUDE <= O(p^(d-1))' -- registered as a P-3 amendment on the strength of")
print("    ONE distance (d=2, giving O(p)).  Prediction at d=3 (Steane) is O(p^2).  Test it.")
print("=" * 100)
D7 = 128


def op7(kind, S):
    m = 0
    for k in S:
        m |= (1 << k)
    if kind == 'Z':
        par = np.array([bin(s & m).count('1') & 1 for s in range(D7)])
        return np.diag(np.where(par == 0, 1.0, -1.0)).astype(complex)
    M = np.zeros((D7, D7), complex)
    b = np.arange(D7)
    M[b ^ m, b] = 1.0
    return M


SX = [[3, 4, 5, 6], [1, 2, 5, 6], [0, 2, 4, 6]]
H7 = -sum(op7('X', s) for s in SX) - sum(op7('Z', s) for s in SX)
R7 = op7('Z', list(range(7)))
W7 = op7('X', list(range(7)))
rng7 = np.random.default_rng(4242)
V7 = np.zeros((D7, D7), complex)
for l in range(7):
    c = rng7.normal(size=3)
    V7 = V7 + c[0] * op7('X', [l]) + c[1] * (1j * op7('X', [l]) @ op7('Z', [l])) + c[2] * op7('Z', [l])
V7 = (V7 + V7.conj().T) / 2
V7 = V7 / np.linalg.norm(V7, 2)


def flip7(p):
    e, U = np.linalg.eigh(H7 + p * V7)
    Uc = U[:, :2]
    M = Uc.conj().T @ R7 @ Uc
    M = (M + M.conj().T) / 2
    w, Q = np.linalg.eigh(M)
    Rc = Q @ np.diag(np.sign(w)) @ Q.conj().T
    best = 0.0
    for l in range(7):
        for O in (op7('X', [l]), op7('Z', [l]), 1j * op7('X', [l]) @ op7('Z', [l])):
            Oc = Uc.conj().T @ O @ Uc
            Om = (Oc - Rc @ Oc @ Rc) / 2
            best = max(best, float(np.linalg.norm(Om, 2)))
    Wl = Uc.conj().T @ W7 @ Uc
    Wm = (Wl - Rc @ Wl @ Rc) / 2
    return best, float(np.linalg.norm(Wm, 2))


PS = [3e-2, 1e-1, 2e-1, 3e-1]
vals = []
print(f"   {'p':>9s} {'STEANE best wt-1 |O_-|':>24s} {'/p':>12s} {'/p^2':>12s} {'logical Xbar':>14s}")
for p in PS:
    b, wl = flip7(p)
    vals.append(b)
    print(f"   {p:9.2e} {b:24.6e} {b/p:12.5f} {b/p**2:12.5f} {wl:14.6f}")
k = np.polyfit(np.log(PS), np.log(vals), 1)[0]
print(f"\n   fitted exponent of the local flip amplitude on d=3 = {k:.4f}")
print(f"   ledger predicts d-1 = 2   -> {'PASS' if abs(k - 2) < 0.15 else 'FAIL'}")
print(f"   naive alternative 'always O(p)' -> {'consistent' if abs(k - 1) < 0.15 else 'REFUTED'}")
