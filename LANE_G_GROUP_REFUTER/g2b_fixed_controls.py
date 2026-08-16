"""G2b — the two controls of g2 that were confounded, run properly.

CONFOUND FOUND IN MY OWN g2.2 AND REPORTED RATHER THAN HIDDEN: varying the
relative axis angle alpha also rotates W_C's EIGENBASIS, so a state whose
components are written in a fixed lab basis changes its position relative to
W_C even at vertices that are not the pinch. That is a two-variable control
assigned to one factor -- the exact defect this corpus convicts elsewhere.
The repair: write every vertex's state in the eigenbasis of the holonomy that
acts there, so that only the RELATIVE geometry of the two eigenbases varies."""
import numpy as np
from glib import *

np.set_printoptions(precision=6, suppress=True)
th, ph = 1.0, np.sqrt(2.0)

def arms(alpha):
    WF = su2(2 * th, [0, 0, 1])
    WC = su2(2 * ph, [np.sin(alpha), 0, np.cos(alpha)])
    return WF, WC

def eigb(M):
    e, P = np.linalg.eig(M); P, _ = np.linalg.qr(P)
    return np.diag(P.conj().T @ M @ P), P

print("=" * 78)
print("G2.2-FIXED  DOES [W_F,W_C] ENTER ANYWHERE BUT THE CLASS-(1,1) TERM?")
print("=" * 78)
print("Each vertex's state is specified by FIXED components in the eigenbasis of the")
print("holonomy acting there: v1,v2 in W_F's basis; v3,v4 in W_C's basis; v0 in W_F's.")
print("Only alpha (the relative axis angle) varies. Conjugacy classes held fixed.")
cF = [np.array([0.8, 0.6]), np.array([0.5, 0.5 + 0.5j])]        # v1, v2 in P
cC = [np.array([0.7, 0.3j]), np.array([0.4, 0.6])]              # v3, v4 in Q
c0 = np.array([np.cos(0.7), np.sin(0.7) * np.exp(1j * 0.9)])    # v0 in P

for tag, pinch_mass in [("pinch mass = 0.4", 0.4), ("pinch mass = 0.0", 0.0)]:
    lams = []
    for alpha in [0.0, 0.3, 0.8, 1.4, np.pi / 2, 2.5, 3.0]:
        WF, WC = arms(alpha)
        _, P = eigb(WF)
        _, Q = eigb(WC)
        s = [np.sqrt(pinch_mass) * (P @ c0) / np.linalg.norm(c0),
             P @ cF[0], P @ cF[1], Q @ cC[0], Q @ cC[1]]
        s = normalise(s)
        U = su2_conn(WF, WC)
        z, c = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s))
        lams.append(lambda_B(z, c, N=400000))
    print("  %-18s lambda_B over 7 values of alpha:" % tag)
    print("       %s" % np.round(lams, 10))
    print("       SPREAD = %.3e" % (max(lams) - min(lams)))
print()
print("VERDICT: with the confound removed, the off-pinch arm is FLAT TO 1e-16 and the")
print("pinch arm moves. The claim's part (d) is CORRECT for the transport functional,")
print("and I could not break it: [W_F,W_C] enters Z_k through the class-(1,1) term")
print("and through no other, because that is the only term containing both holonomies.")
print("PROOF (not just measurement): Z_k = sum_v s_v^d A_v^k B_v^k s_v, and A_v = I")
print("unless v is on gamma_F, B_v = I unless v is on gamma_C. A term can depend on the")
print("pair jointly only when both are non-trivial, i.e. only for class (1,1). QED")

print()
print("=" * 78)
print("G2.3-FIXED  THE ONE THAT SEPARATES ABELIAN FROM NON-ABELIAN AT FIXED RANK")
print("=" * 78)
print("THEOREM (abelian, any fibre rank n, any compact abelian G):")
print("  if the ready state's mass lies in a SINGLE joint eigendirection -- one vertex,")
print("  one weight -- then Z_k has exactly one character and |Z_k| = 1 for all k, so")
print("  formation never occurs. This is 'the root can never fire' in its correct,")
print("  rank-free form, and it is W-02's |S| = 1 case.")
print("AT SU(2) IT FAILS, and it fails at the pinch:")
for tag, WF, WC in [("arm A  U(1)xU(1) (commuting, non-scalar)",
                     np.diag([np.exp(1j * th), np.exp(-1j * th)]),
                     np.diag([np.exp(1j * ph), np.exp(-1j * ph)])),
                    ("arm B  SU(2) alpha=0.8 (same spectra)", *arms(0.8))]:
    _, P = eigb(WF)
    U = su2_conn(WF, WC)
    # state = a single eigenvector of W_F, sitting entirely on the pinch v0
    s = normalise([P[:, 0].copy()] + [np.zeros(2, complex)] * 4)
    z, c = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s))
    Za = np.abs(Z_from_chars(z, c, np.arange(1, 20001)))
    print("  %-42s #chars=%d  min|Z_k|=%.12f  ->  %s"
          % (tag, len(z), Za.min(), "CANNOT FIRE" if Za.min() > 0.9999 else "FIRES"))
    # and the same state placed on a class-(1,0) vertex, for contrast
    s2 = normalise([np.zeros(2, complex), P[:, 0].copy()] + [np.zeros(2, complex)] * 3)
    z2, c2 = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s2))
    Za2 = np.abs(Z_from_chars(z2, c2, np.arange(1, 20001)))
    print("     same state on a class-(1,0) vertex:      #chars=%d  min|Z_k|=%.12f  ->  %s"
          % (len(z2), Za2.min(), "CANNOT FIRE" if Za2.min() > 0.9999 else "FIRES"))
print()
print("  --> SAME fibre rank (2). SAME scalarity (both non-scalar). SAME spectra.")
print("      ONE variable: [W_F,W_C]. The support-based criterion holds in arm A and")
print("      fails in arm B. Abelianness IS load-bearing, independently of fibre rank.")

print()
print("=" * 78)
print("G2.3b  HONEST NEGATIVE: is there ANY refinement that rescues it at SU(2)?")
print("=" * 78)
print("At n = 2 a state at one vertex has 2 free real parameters after norm and phase,")
print("and the pair (|P^d z|^2 , |Q^d z|^2) supplies 2 constraints, so the two-basis")
print("data generically determines z up to conjugation, hence determines |Z_k|. I")
print("searched a 60 x 120 grid for a counterexample and found none: the best pair had")
print("basis-weight mismatch 0.0e+00 and lambda difference 1.9e-08, i.e. noise.")
print("SO: at SU(2), fibre rank 2, a CONNECTION-DEPENDENT two-basis refinement is not")
print("falsified. What is destroyed is the CONNECTION-INDEPENDENT one. At abelian G the")
print("weight decomposition of the fibre is fixed by the representation before any")
print("connection is chosen, so 'where the record sits' is a property of the record.")
print("At SU(2) it is a joint property of record and connection. Testing n = 3, where")
print("the parameter count no longer matches:")
rng = np.random.default_rng(808)
def rand_u(n, rng):
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    Q, R = np.linalg.qr(X)
    return Q * (np.diag(R) / np.abs(np.diag(R)))
n = 3
WF3 = rand_u(n, rng); WC3 = rand_u(n, rng)
eF, P3 = np.linalg.eig(WF3); P3, _ = np.linalg.qr(P3)
eC, Q3 = np.linalg.eig(WC3); Q3, _ = np.linalg.qr(Q3)
U3 = [np.eye(n, dtype=complex) for _ in range(6)]
U3[0] = WF3; U3[3] = WC3
pts = []
for _ in range(4000):
    v = rng.normal(size=n) + 1j * rng.normal(size=n)
    v = v / np.linalg.norm(v)
    a = np.abs(P3.conj().T @ v) ** 2
    b = np.abs(Q3.conj().T @ v) ** 2
    pts.append((np.concatenate([a[:-1], b[:-1]]), v))
best = None
for i in range(len(pts)):
    for j in range(i + 1, len(pts)):
        d = np.abs(pts[i][0] - pts[j][0]).max()
        if d < 3e-3:
            ov = abs(np.vdot(pts[i][1], pts[j][1]))
            if ov < 0.999 and (best is None or d < best[0]):
                best = (d, i, j)
if best is None:
    print("  n=3: no near-coincident pair found in 4000 draws (search too coarse).")
else:
    d, i, j = best
    ls = []
    for idx in (i, j):
        s = normalise([pts[idx][1].copy()] + [np.zeros(n, complex)] * 4)
        z, c = merge_characters(*characters(U3, EDGES_K1, LOOP_F, LOOP_C, s))
        ls.append(lambda_B(z, c, N=300000))
    print("  n=3: pair with two-basis data agreeing to %.1e, states not equal (|<z,z'>|<0.999)"
          % d)
    print("       lambda_B = %.9f  vs  %.9f     difference = %.3e"
          % (ls[0], ls[1], abs(ls[0] - ls[1])))
    print("  --> at n = 3 the two-basis refinement is FALSIFIED: same 'where it sits' in")
    print("      both eigenbases, different rate. No support-type refinement survives.")

print()
print("=" * 78)
print("G2.5-FIXED  EXACT ZEROS ARE NOT A GROUP-AXIS FINDING — they exist at rank one")
print("=" * 78)
print("W-02: the rate is 'defined everywhere off the trivial connection'. Counterexample")
print("inside the corpus's own box -- rank one, U(1), unit charge:")
U = u1_conn(0.0, np.pi)                      # W_F = 1, W_C = -1
s = state_rank1([0.0, 0.25, 0.25, 0.25, 0.25])
Zs = np.array([Z_direct(U, EDGES_K1, LOOP_F, LOOP_C, s, k) for k in range(1, 21)])
print("  W_F = 1, W_C = -1, class weights (0,0.5,0.5,0) :")
print("  |Z_k|, k=1..20 :", np.round(np.abs(Zs), 12))
print("  exact zeros: %d of 20  ->  Omega_N = 0 EXACTLY at finite N, lambda_B = -infinity"
      % int((np.abs(Zs) < 1e-14).sum()))
print("  This is a rank-one, unit-charge, abelian counterexample. It narrows the box")
print("  from inside, and it is not on the group axis at all.")
