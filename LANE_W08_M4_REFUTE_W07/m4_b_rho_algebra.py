# W-08 / M4 leg B — IS W-07's rho THE RATIO THAT ACTUALLY GOVERNS ITS OWN OBSERVABLE?
#
# W-07 names its operative variable "ord(rho), rho the branch ratio in U(1)".  The formula for rho
# appears in exactly one place: w07_e_isolation.py, the leg that carries the five-for-five table.
#     rho_E = conj(W_F)**dF * W_C**(-dC)
# Derive the ratio that governs the dressed separation from scratch and compare.
#
# ISOLATION LEDGER (this file).  Held fixed: carrier K1, dressing tree, ready state (seed 20260816
# as in W-07 leg E), observable A_23, k-range 1..4000, code path (W-07's own dress()).
# Moved: NOTHING.  This leg is pure algebra + one arithmetic identity check.  No comparison is drawn.
# Double precision (numpy float64) throughout; the one claim that is precision-sensitive
# (exact zero vs float noise) is checked separately in m4_g_exact.py with integer arithmetic.
import numpy as np

FACE_V = {0, 1, 2}; CYC_V = {0, 3, 4}
TREE = {1: (0,), 2: (0, 1), 3: (3,), 4: (3, 4)}

def dress(s, a):                      # W-07's dress(), verbatim
    u = np.exp(1j*np.asarray(a)); t = np.array(s, dtype=complex)
    for v, p in TREE.items():
        w = 1.0+0j
        for e in p: w *= u[e]
        t[v] = s[v]/w
    return t

print("== B0  THE DERIVATION, DONE FROM THE DEFINITION ==")
print("  A_uv[s] = conj(t_u) t_v,  t_v = W(tree path)^-1 s_v.  Dressing is diagonal with UNIT-MODULUS")
print("  entries, so |t_v| = |s_v| exactly, and M_dF/M_c pass through it unchanged:")
print("    A_uv[M_dF^k s] = W_F^(k*dF) A_uv[s],   dF = [v in F] - [u in F]")
print("    A_uv[M_c^k  s] = W_C^(k*dC) A_uv[s],   dC = [v in C] - [u in C]")
print("    D_k = |A[M_dF^k s] - A[M_c^k s]| = amp * |W_F^(dF k) - W_C^(dC k)|")
print("        = amp * | (W_F^dF * W_C^-dC)^k - 1 |            <-- TRUE ratio  rho_true = W_F^dF W_C^-dC")
print("  W-07 leg E instead sets      rho_E = conj(W_F)^dF * W_C^-dC = W_F^(-dF) W_C^(-dC).")
print("  For the published observable A_23:  dF = -1, dC = +1")
print("    rho_true = W_F^-1 W_C^-1 = 1/(W_F W_C)        rho_E = W_F^+1 W_C^-1 = W_F/W_C")
print("  THESE ARE DIFFERENT ELEMENTS OF U(1) unless W_F^2 = 1.")
print("  On S1's published connection W_F = -1, so W_F^2 = 1 and they coincide up to conjugation:")
print("    rho_true = 1/((-1)(-i)) = -i ,  rho_E = (-1)/(-i) = -i .   The slip is INVISIBLE there.\n")

def D_true(a, s, u=2, v=3, K=4000):
    WF = np.exp(1j*sum(a[:3])); WC = np.exp(1j*sum(a[3:]))
    dF = (v in FACE_V) - (u in FACE_V); dC = (v in CYC_V) - (u in CYC_V)
    t = dress(s, a); amp = abs(np.conj(t[u])*t[v]); k = np.arange(1, K+1)
    return amp*np.abs(WF**(dF*k) - WC**(dC*k)), WF, WC, amp, dF, dC

def D_w07legE(a, s, u=2, v=3, K=4000):     # w07_e_isolation.py run(), verbatim
    WF = np.exp(1j*sum(a[:3])); WC = np.exp(1j*sum(a[3:]))
    dF = (v in FACE_V) - (u in FACE_V); dC = (v in CYC_V) - (u in CYC_V)
    t = dress(s, a); amp = abs(np.conj(t[u])*t[v]); k = np.arange(1, K+1)
    rho = np.conj(WF)**dF*WC**(-dC)
    D = amp*np.abs(np.exp(1j*np.angle(rho)*k) - 1)
    ordr = next((n for n in range(1, 10001) if abs(rho**n-1) < 1e-12), None)
    return D, rho, ordr

rng = np.random.default_rng(20260816)
s = rng.normal(size=5) + 1j*rng.normal(size=5); s /= np.linalg.norm(s)   # W-07 leg E's state, same seed

print("== B1  ON S1's PUBLISHED CONNECTION THE TWO AGREE (so W-07 could not have seen the slip) ==")
a_pub = np.array([np.pi/3]*3 + [np.pi/2]*3)
Dt, WF, WC, amp, dF, dC = D_true(a_pub, s); De, rhoE, ordE = D_w07legE(a_pub, s)
print(f"  rho_true = {WF**dF*WC**(-dC):+.9f}   rho_E = {rhoE:+.9f}")
print(f"  max_k |D_true - D_legE| = {np.abs(Dt-De).max():.3e}     (agree to float noise)\n")

print("== B2  A CONNECTION WHERE THEY DISAGREE, AND W-07's LEG E GETS ITS OWN OBSERVABLE WRONG ==")
print("  Construction: pick theta1 irrational; set arg(W_F) = 2pi*theta1, arg(W_C) = 2pi*(1/4 - theta1).")
print("  Then W_F*W_C = exp(2pi i/4) = i has ORDER 4  ->  rho_true has order 4  ->  EXACT ZEROS.")
print("  And W_F/W_C = exp(2pi i (2 theta1 - 1/4)) is irrational  ->  rho_E has INFINITE order.")
th1 = np.sqrt(2) % 1.0                                  # irrational
aF = 2*np.pi*th1; aC = 2*np.pi*(0.25 - th1)
a_adv = np.array([aF/3]*3 + [aC/3]*3)
Dt, WF, WC, amp, dF, dC = D_true(a_adv, s); De, rhoE, ordE = D_w07legE(a_adv, s)
rho_true = WF**dF*WC**(-dC)
print(f"  W_F = {WF:+.9f}   W_C = {WC:+.9f}   amp = {amp:.9f}")
print(f"  rho_true = {rho_true:+.9f}  arg/2pi = {np.angle(rho_true)/(2*np.pi):+.9f}   |rho_true^4 - 1| = {abs(rho_true**4-1):.3e}")
print(f"  rho_E    = {rhoE:+.9f}  arg/2pi = {np.angle(rhoE)/(2*np.pi):+.9f}   ORDER as W-07 computes it = {ordE}")
print()
print(f"  TRUE dressed separation   : min over k<=4000 = {Dt.min():.3e}   cells < 1e-9 : {int((Dt<1e-9).sum())} of 4000")
print(f"  W-07 leg E would report   : min over k<=4000 = {De.min():.3e}   cells < 1e-9 : {int((De<1e-9).sum())} of 4000")
print(f"  cells where D_true is exactly-zero-to-float (k = 0 mod 4): {int((np.arange(1,4001)%4==0).sum())}")
print()
print("  W-07's leg E, run on this connection, reports ORDER = None and NO recurrence.")
print("  The dressed observable it claims to be measuring recurs on 1000 of 4000 cells.")
print("  ==> the quantity W-07 names 'rho' is not the quantity that governs W-07's own observable.\n")

print("== B3  AND 'THE BRANCH RATIO' IS NOT ONE ELEMENT: IT DEPENDS ON THE VERTEX PAIR ==")
print("  (dF,dC) per ordered pair, and the group element that governs A_uv:")
lab = {0: "0 (F and C)", 1: "1 (F only)", 2: "2 (F only)", 3: "3 (C only)", 4: "4 (C only)"}
seen = {}
for u in range(5):
    for v in range(5):
        if u == v: continue
        dF = (v in FACE_V)-(u in FACE_V); dC = (v in CYC_V)-(u in CYC_V)
        seen.setdefault((dF, dC), []).append(f"({u},{v})")
for (dF, dC), pairs in sorted(seen.items()):
    gov = "trivial (D=0 for every k, every connection)" if (dF, dC) == (0, 0) else f"W_F^{dF} W_C^{-dC}"
    print(f"    (dF,dC) = ({dF:+d},{dC:+d})  governs by {gov:<44} pairs: {' '.join(pairs)}")
print("  On S1's published connection W_F=-1 (order 2), W_C=-i (order 4):")
for (dF, dC), pairs in sorted(seen.items()):
    r = (-1.0+0j)**dF * (-1j)**(-dC)
    o = next((n for n in range(1, 100) if abs(r**n-1) < 1e-12), None)
    print(f"    ({dF:+d},{dC:+d}) -> element {r:+.3f}, order {o}, zeros among k<=4000 : {4000 if o==1 else 4000//o}")
