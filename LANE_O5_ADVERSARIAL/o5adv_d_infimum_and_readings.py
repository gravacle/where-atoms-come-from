"""O5 ADVERSARIAL C.
 (a) exponent of eps_dressed (the quantity that actually sets T) on both carriers
 (b) is the STATED approximate definition self-selecting, or does it import R_0?
 (c) does the lifetime derivation adjudicate between Reading 1 and Reading 2?
"""
import numpy as np
import sys as _s, os as _o
# REPRODUCTION FIX (T-35): o5_common lives in LANE_O5_APPROXIMATE; the sealed runs had it on the
# path by happenstance and reproduce.sh could not run this lane standalone.
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..', 'LANE_O5_APPROXIMATE'))
from o5_common import Zop, Xop, toric_H, sym_H, local_perturbation, Z_A_SUP

V = local_perturbation(seed=2026)


def block(H0, R0, g, p):
    e, U = np.linalg.eigh(H0 + p * V)
    Uc = U[:, :g]
    M = Uc.conj().T @ R0 @ Uc
    M = (M + M.conj().T) / 2
    w, Q = np.linalg.eigh(M)
    Rc = Q @ np.diag(np.sign(w)) @ Q.conj().T
    Hc = np.diag(e[:g]).astype(complex)
    return e[:g], M, Rc, Hc, e[g - 1] - e[0], float(np.linalg.norm(Hc @ Rc - Rc @ Hc, 2))


print("=" * 100)
print("(a) EXPONENT OF eps_dressed -- the quantity that actually sets T -- not of the cluster width")
print("=" * 100)
PS = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1])
for nm, H0, R0, g, dd in (("SYMMETRY d=1", sym_H(), Zop([0]), 2, 1),
                          ("TORIC d=2", toric_H(), Zop(Z_A_SUP), 4, 2)):
    ed = [block(H0, R0, g, p)[5] for p in PS]
    k = np.polyfit(np.log(PS), np.log(ed), 1)[0]
    print(f"   {nm:>14s}  fitted exponent of eps_dressed = {k:.4f}   code distance {dd}   "
          f"{'PASS -- exponent claim survives on the right quantity' if abs(k - dd) < 0.02 else 'FAIL'}")

print()
print("=" * 100)
print("(b) IS THE STATED APPROXIMATE DEFINITION SELF-SELECTING?")
print("    Ledger: (ii-w) = width bound on the cluster carrying the record; P-1 amended to")
print("    '(iii-w) => H has a cluster of width <= w and dimension > 1'.  Take those at face value")
print("    and ask for the INFIMUM of eps over operators satisfying (i) + (iii-w) on that cluster.")
print("=" * 100)
for nm, H0, R0, g in (("TORIC 2x2", toric_H(), Zop(Z_A_SUP), 4), ("SYMMETRY", sym_H(), Zop([0]), 2)):
    p = 1e-3
    ee, M, Rc, Hc, w_, ed = block(H0, R0, g, p)
    Ralt = np.diag(np.array([1.0] * (g // 2) + [-1.0] * (g - g // 2))).astype(complex)
    eps_alt = float(np.linalg.norm(Hc @ Ralt - Ralt @ Hc, 2))
    nonconst = float(np.linalg.norm(Ralt - np.trace(Ralt) / g * np.eye(g), 2))
    inv = float(np.linalg.norm(Ralt @ Ralt - np.eye(g), 2))
    print(f"   {nm:>10s}  p = {p:.0e}   cluster dim {g}, width {w_:.4e}")
    print(f"   {'':>10s}    lane's R = sign(P R_0 P)          : eps = {ed:.6e}")
    print(f"   {'':>10s}    ALTERNATIVE R = diag(+1,..,-1,..) in the PERTURBED energy basis:")
    print(f"   {'':>10s}       (i)   R=Rdag, R^2=I  residual  = {inv:.2e}")
    print(f"   {'':>10s}       (iii-w) non-constant on cluster = {nonconst:.4f}")
    print(f"   {'':>10s}       eps = ||[H_c,R]||               = {eps_alt:.6e}   <-- ZERO")
print("""
   The infimum of eps over the STATED relaxed clauses (i)+(ii-w)+(iii-w) is 0 and is ATTAINED, on
   BOTH carriers, inside the SAME cluster.  This is not the lane's where_it_breaks (2) ('a chance
   near-degeneracy elsewhere in the 256 levels') -- it is the generic case in the very cluster the
   lane measures.  The 2.04e+05 separation exists only because R was FIXED to sign(P R_0 P), which
   imports the EXACT record of the UNPERTURBED carrier.  The approximate definition as stated does
   not select that operator; something else must, and only (iv)+(v) are left to do it.""")

print()
print("=" * 100)
print("(c) DOES THE LIFETIME DERIVATION ADJUDICATE BETWEEN THE TWO READINGS?")
print("    For ANY R:  ||[H,R(t)]|| = ||[H,R]|| for all t (unitary conjugation), so")
print("        ||R(t) - R||  <=  t * ||[H,R]||   EXACTLY, with no cluster and no projector.")
print("    So 'epsilon is an inverse time' is equally true of READING 1.  The Ehrenfest argument")
print("    does not choose between them.  What must be shown is that Reading 1's motion is BOUNDED")
print("    (never accumulates) while Reading 2's is SECULAR.  The lane never tests this.  Do it:")
print("=" * 100)


def fulldev_max(H, R0, tmax, n):
    e, U = np.linalg.eigh(H)
    Rb = U.conj().T @ R0 @ U
    W = e[:, None] - e[None, :]
    best = 0.0
    for t in np.linspace(0, tmax, n):
        D = Rb * (np.exp(1j * W * t) - 1.0)
        s = np.linalg.norm(D, 2)
        best = max(best, s)
    return best


for nm, H0, R0 in (("TOPOLOGICAL", toric_H(), Zop(Z_A_SUP)), ("SYMMETRY", sym_H(), Zop([0]))):
    for p in (1e-3, 1e-2):
        H = H0 + p * V
        epsf = float(np.linalg.norm(H @ R0 - R0 @ H, 2))
        mx = fulldev_max(H, R0, 200.0, 161)
        Tb = 0.01 / epsf
        print(f"   {nm:>12s} p={p:.0e}  eps_fixed = {epsf:.4e}   Reading-1 bound says eta=0.01 by "
              f"T = {Tb:.3e}")
        print(f"   {'':>12s}   ACTUAL max_t ||R_0(t)-R_0|| over t in [0,200] "
              f"(= {200/Tb:.1e} x that bound) = {mx:.4e}")
        print(f"   {'':>12s}   -> {'SECULAR' if mx > 0.5 else 'BOUNDED: saturates, never accumulates'}")
print("""
   Reading 1's epsilon is a VALID but hopelessly loose lifetime bound: the fixed record's full-space
   deviation saturates at O(p/Delta) and oscillates for ever.  THIS is the reason Reading 2 is the
   right one -- and it is an argument the lane does not make.  The lane's stated reason for choosing
   Reading 2 is that Reading 2 preserves W-61's separation, i.e. the reading is selected by the
   answer it returns.  The correct reason is available and is measured above.""")
