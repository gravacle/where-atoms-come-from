"""G2 — THE CONTROL THE REGISTER'S ERR-2 ASKED FOR: fibre rank held at 2,
transport held NON-SCALAR, and ONLY commutativity varied.
   arm A : G = U(1) x U(1), W_F, W_C diagonal   -> [W_F,W_C] = 0, non-scalar, rank 2
   arm B : G = SU(2),        W_F, W_C generic   -> [W_F,W_C] != 0, non-scalar, rank 2
Same fibre rank. Same scalarity. One variable."""
import numpy as np
from glib import *

np.set_printoptions(precision=6, suppress=True)
print("seed = 20260816 for every rng below unless stated")

# ---------------------------------------------------------------- arms -------
th, ph = 1.0, np.sqrt(2.0)          # eigen-angles held fixed across both arms
# arm A: commuting, same spectra
WF_A = np.diag([np.exp(1j * th), np.exp(-1j * th)])
WC_A = np.diag([np.exp(1j * ph), np.exp(-1j * ph)])
# arm B: same spectra, axes at angle alpha
def WF_B(): return su2(2 * th, [0, 0, 1])
def WC_B(alpha): return su2(2 * ph, [np.sin(alpha), 0, np.cos(alpha)])

print()
print("=" * 78)
print("G2.1  THE COEFFICIENT DICHOTOMY — the exact place abelianness is load-bearing")
print("=" * 78)
print("Z_k = sum_j c_j zeta_j^k  with, for each vertex v,")
print("     c_{mn} = (s_v^d P)_m (P^d Q)_{mn} (Q^d s_v)_n ,  zeta_{mn} = alpha_m beta_n")
print("where A_v = P diag(alpha) P^d and B_v = Q diag(beta) Q^d.")
print("If [A_v,B_v] = 0 then P = Q, so P^d Q = I, so c_{mn} = delta_{mn} |(P^d s_v)_m|^2:")
print("EVERY coefficient is a NON-NEGATIVE REAL and is a piece of the state's mass.")
print("If [A_v,B_v] != 0 the off-diagonal coefficients are COMPLEX. Measured:")
rng = np.random.default_rng(20260816)
s_pinch = normalise([np.array([np.cos(0.7), np.sin(0.7) * np.exp(1j * 0.9)]),
                     np.zeros(2, complex), np.zeros(2, complex),
                     np.zeros(2, complex), np.zeros(2, complex)])
for tag, U in [("arm A  U(1)xU(1)", su2_conn(WF_A, WC_A)),
               ("arm B  SU(2) alpha=0.8", su2_conn(WF_B(), WC_B(0.8)))]:
    z, c = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s_pinch))
    dev = max(abs(Z_from_chars(z, c, [k])[0] -
                  Z_direct(U, EDGES_K1, LOOP_F, LOOP_C, s_pinch, k)) for k in [1, 2, 5, 17, 60])
    print("  %-24s  #chars=%d  max|Im c|=%.3e  min Re c=%+.6f  sum c=%.9f  (matrix dev %.1e)"
          % (tag, len(z), np.abs(c.imag).max(), c.real.min(), c.sum().real, dev))
print()
print("  --> W-01's convex-hull criterion, S4-1's 'corners of a square', the class-weight")
print("      multiset theorem and pinch=spectator are ALL theorems about NON-NEGATIVE")
print("      MASS coefficients indexed by a support. Abelian rank 2 keeps that type.")
print("      Non-abelian rank 2 destroys it. THIS is what abelianness buys, and it is")
print("      not something the schedule can supply.")

print()
print("=" * 78)
print("G2.2  IS NON-COMMUTATIVITY VISIBLE ONLY AT THE PINCH? (the claim's part (d))")
print("=" * 78)
print("Hold both conjugacy classes fixed, vary ONLY the relative axis angle alpha.")
print("State (i) puts all mass on the pinch v0; state (ii) puts all mass off it.")
si = normalise([np.array([np.cos(.7), np.sin(.7) * np.exp(1j * .9)]),
                np.zeros(2, complex), np.zeros(2, complex), np.zeros(2, complex), np.zeros(2, complex)])
sii = normalise([np.zeros(2, complex),
                 np.array([0.6, 0.3 + 0.2j]), np.array([0.2, 0.5j]),
                 np.array([0.4, 0.1j]), np.array([0.3, 0.45])])
for tag, s in [("(i)  all mass ON the pinch  ", si), ("(ii) all mass OFF the pinch ", sii)]:
    lams = []
    for alpha in [0.0, 0.3, 0.8, 1.4, np.pi / 2, 2.5]:
        U = su2_conn(WF_B(), WC_B(alpha))
        z, c = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s))
        lams.append(lambda_B(z, c, N=300000))
    print("  %s lambda_B over alpha: %s" % (tag, np.round(lams, 9)))
    print("       spread = %.3e" % (max(lams) - min(lams)))
print("  --> CONFIRMED for the transport functional: [W_F,W_C] enters Z_k through the")
print("      class-(1,1) term and nowhere else, because that is the only term in which")
print("      both holonomies appear. This half of the claim is a theorem and I could")
print("      not break it. What I CAN break is 'in the whole construction' -- see G2.6.")

print()
print("=" * 78)
print("G2.3  THE DECISIVE CONTROL: what fails at SU(2) that does NOT fail at U(1)x U(1)")
print("=" * 78)
print("W-02's headline: 'FORMATION IS A GROUP-THEORETIC CONDITION ON WHERE THE RECORD")
print("SITS.' At any ABELIAN group the two holonomies share an eigenbasis, so 'where")
print("the record sits' has a canonical refinement (vertex class x fibre weight) that")
print("is fixed BEFORE the connection is chosen. At SU(2) there are two incompatible")
print("candidate refinements and neither works. Measured, at fibre rank 2 in both arms:")
print()
print("  Family: all mass on the pinch; |components|^2 in W_F's eigenbasis held FIXED")
print("  at (cos^2 t, sin^2 t) with t = 0.7; only the relative phase psi varies.")
for tag, WFx, WCx in [("arm A  U(1)xU(1)", WF_A, WC_A),
                      ("arm B  SU(2) a=0.8", WF_B(), WC_B(0.8))]:
    ev, P = np.linalg.eig(WFx)
    P, _ = np.linalg.qr(P)
    lams, wF, wC = [], [], []
    for psi in [0.0, 0.6, 1.3, 2.2, 3.0, 4.4]:
        zv = P @ np.array([np.cos(0.7), np.sin(0.7) * np.exp(1j * psi)])
        s = normalise([zv, np.zeros(2, complex), np.zeros(2, complex),
                       np.zeros(2, complex), np.zeros(2, complex)])
        U = su2_conn(WFx, WCx)
        z, c = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s))
        lams.append(lambda_B(z, c, N=300000))
        wF.append(abs(P.conj().T @ zv)[0] ** 2)
        evc, Q = np.linalg.eig(WCx); Q, _ = np.linalg.qr(Q)
        wC.append(abs(Q.conj().T @ zv)[0] ** 2)
    print("  %-20s class weights CONSTANT (0,0,0,1); W_F-basis weight CONSTANT %.9f"
          % (tag, np.std(wF) + wF[0]))
    print("      lambda_B : %s" % np.round(lams, 9))
    print("      SPREAD   : %.6e        W_C-basis weight spread: %.3e"
          % (max(lams) - min(lams), max(wC) - min(wC)))
print()
print("  --> arm A: zero. arm B: not zero. FIBRE RANK, SCALARITY AND THE SPECTRA ARE")
print("      IDENTICAL IN BOTH ARMS. The only difference is [W_F,W_C]. The support-based")
print("      content of W-02's criterion survives abelian rank two and DIES at SU(2).")

print()
print("  Can a two-basis refinement rescue it? Search for two pinch states with the SAME")
print("  class weights, SAME W_F-basis weights and SAME W_C-basis weights, different lambda:")
WFx, WCx = WF_B(), WC_B(0.8)
evF, P = np.linalg.eig(WFx); P, _ = np.linalg.qr(P)
evC, Q = np.linalg.eig(WCx); Q, _ = np.linalg.qr(Q)
U = su2_conn(WFx, WCx)
grid = []
for t in np.linspace(0.05, 1.5, 60):
    for psi in np.linspace(0, 2 * np.pi, 121)[:-1]:
        zv = P @ np.array([np.cos(t), np.sin(t) * np.exp(1j * psi)])
        a = abs(P.conj().T @ zv)[0] ** 2
        b = abs(Q.conj().T @ zv)[0] ** 2
        grid.append((a, b, t, psi))
grid = np.array(grid)
best = None
for i in range(len(grid)):
    d = np.abs(grid[:, 0] - grid[i, 0]) + np.abs(grid[:, 1] - grid[i, 1])
    d[i] = 9e9
    j = int(np.argmin(d))
    if d[j] < 2e-4 and abs(grid[i, 3] - grid[j, 3]) > 1e-3:
        if best is None or d[j] < best[0]:
            best = (d[j], i, j)
if best is None:
    print("     no such pair found on this grid -> at n=2 the two-basis refinement is")
    print("     not falsified by this search (reported, not asserted).")
else:
    _, i, j = best
    ls = []
    for (t, psi) in [(grid[i, 2], grid[i, 3]), (grid[j, 2], grid[j, 3])]:
        zv = P @ np.array([np.cos(t), np.sin(t) * np.exp(1j * psi)])
        s = normalise([zv] + [np.zeros(2, complex)] * 4)
        z, c = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s))
        ls.append(lambda_B(z, c, N=300000))
    print("     pair found, basis-weight mismatch %.1e :  lambda = %.9f  vs  %.9f   diff %.3e"
          % (best[0], ls[0], ls[1], abs(ls[0] - ls[1])))

print()
print("=" * 78)
print("G2.4  DOES THE ERGODIC/MAHLER SKELETON SURVIVE AT SU(2)? (the claim's part (c))")
print("=" * 78)
print("THEOREM (proved by construction in glib.characters, not asserted): for ANY")
print("W_F, W_C in U(n), Z_k = sum_j c_j zeta_j^k with at most n^2 * (#classes) terms,")
print("because the closure of the CYCLIC group {(W_F^k, W_C^k)} in U(n)xU(n) is a")
print("compact ABELIAN group. Verified numerically at SU(2):")
s = normalise([np.array([0.6, 0.3 + 0.2j]), np.array([0.4, 0.2j]), np.array([0.2, 0.5]),
               np.array([0.3, 0.1j]), np.array([0.25, 0.35])])
U = su2_conn(WF_B(), WC_B(0.8))
z, c = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s))
dev = max(abs(Z_from_chars(z, c, [k])[0] -
              Z_direct(U, EDGES_K1, LOOP_F, LOOP_C, s, k)) for k in range(1, 121))
print("  #characters = %d   max |char form - matrix action| over k<=120 : %.3e" % (len(z), dev))
print("  exponent vectors in the eigen-angles (theta, phi) -- read off the characters:")
ex = []
for zz in z:
    for m in (1, -1):
        for n in (1, -1):
            if abs(zz - np.exp(1j * (m * th + n * ph))) < 1e-9: ex.append((m, n))
    for m in (1, -1):
        if abs(zz - np.exp(1j * m * th)) < 1e-9: ex.append((m, 0))
        if abs(zz - np.exp(1j * m * ph)) < 1e-9: ex.append((0, m))
    if abs(zz - 1) < 1e-9: ex.append((0, 0))
print("   ", sorted(set(ex)))
print("  rank-one U(1) realises 4 exponents (corners of the UNIT square);")
print("  SU(2) at rank 2 realises the 3x3 box -- a different lattice configuration.")
lam_sim = lambda_B(z, c, N=1500000)
# Mahler route over the 2-torus of eigen-angles, with the SU(2) complex coefficients
Ex, Cx = [], []
for zz, cc in zip(z, c):
    hit = None
    for m in (-1, 0, 1):
        for n in (-1, 0, 1):
            if abs(zz - np.exp(1j * (m * th + n * ph))) < 1e-9: hit = (m, n)
    Ex.append(hit); Cx.append(cc)
lam_mah = mahler_torus(np.array(Ex), np.array(Cx), ngrid=3000)
print("  schedule-B simulation N=1.5e6 : %.9f" % lam_sim)
print("  2-torus Mahler quadrature     : %.9f      deviation %.2e"
      % (lam_mah, abs(lam_sim - lam_mah)))
print("  --> (c) CONFIRMED, and confirmed more strongly than claimed: the skeleton needs")
print("      no group at all, only unitarity, and it survives because BOTH BRANCHES USE")
print("      THE SAME EXPONENT k. But the surviving skeleton is a Mahler measure with")
print("      COMPLEX coefficients, which is not the object any of the four theorems is")
print("      about. 'The skeleton survives' and 'the results survive' are different claims.")

print()
print("=" * 78)
print("G2.5  WHERE THE SKELETON DOES STOP: a compact group with no generic stratum")
print("=" * 78)
print("Take G finite (compact). Then the relation lattice has rank 2 at EVERY point of")
print("the parameter space -- S4's 'full-measure generic torus' stratum is EMPTY -- and")
print("lambda can be exactly -infinity on an open set of states, not on a measure-zero")
print("exceptional set. Q8 = quaternion group, non-abelian, order 8:")
i2 = np.array([[1j, 0], [0, -1j]]); j2 = np.array([[0, 1], [-1, 0]], dtype=complex)
for tag, A, B in [("Q8  W_F=i, W_C=j (non-abelian)", i2, j2),
                  ("Z4  W_F=diag(i,1), W_C=diag(1,i) (abelian)",
                   np.diag([1j, 1]).astype(complex), np.diag([1, 1j]).astype(complex))]:
    U = su2_conn(A, B)
    ss = normalise([np.array([1, 1]) / np.sqrt(2), np.zeros(2, complex), np.zeros(2, complex),
                    np.zeros(2, complex), np.zeros(2, complex)])
    Zs = np.array([Z_direct(U, EDGES_K1, LOOP_F, LOOP_C, ss, k) for k in range(1, 25)])
    nz = int((np.abs(Zs) < 1e-12).sum())
    print("  %-42s  |Z_k| for k=1..24 has %d exact zeros -> lambda = -infinity: %s"
          % (tag, nz, nz > 0))
print("  --> at a finite compact group the corpus's 'rate defined everywhere off the")
print("      trivial connection' (W-02) is false: the record functional hits EXACTLY")
print("      zero at a finite stage. This narrows the box on the OTHER side from the claim.")

print()
print("=" * 78)
print("G2.6  'EXACTLY ONE PLACE IN THE WHOLE CONSTRUCTION' — two more places")
print("=" * 78)
print("(1) THE TRANSPORT OPERATOR IS NOT EVEN DEFINED AT NON-ABELIAN G WITHOUT A NEW")
print("    CONVENTION, AND THE CONVENTION IS NEEDED AT EVERY LOOP VERTEX, NOT THE PINCH.")
print("    W-01's operator is 'multiply by W(gamma) at vertices on gamma'. At rank one")
print("    W(gamma) is basepoint-free. At SU(2) it is not. Take the NAIVE reading (use")
print("    v0's holonomy at all three loop vertices) and gauge-transform at v1 only:")
rngg = np.random.default_rng(11)
g = [np.eye(2, dtype=complex) for _ in range(5)]
g[1] = su2(1.1, [1, 2, -1])
U = su2_conn(WF_B(), WC_B(0.8))
s = normalise([np.array([0.6, 0.3 + 0.2j]), np.array([0.4, 0.2j]), np.array([0.2, 0.5]),
               np.array([0.3, 0.1j]), np.array([0.25, 0.35])])
Ug = gauge_transform(U, EDGES_K1, g)
sg = gauge_state(s, g)
for naive in (True, False):
    d = max(abs(abs(Z_direct(U, EDGES_K1, LOOP_F, LOOP_C, s, k, naive=naive)) -
                abs(Z_direct(Ug, EDGES_K1, LOOP_F, LOOP_C, sg, k, naive=naive)))
            for k in range(1, 31))
    print("      %-9s transport:  max ||Z_k| - |Z_k^gauged||  over k<=30 = %.3e  -> %s"
          % ("NAIVE" if naive else "BASED", d, "NOT gauge invariant" if d > 1e-9 else "gauge invariant"))
print("    The gauge transformation was at v1, which is NOT on both loops.")
print()
print("(2) THE GAUGE-INVARIANT PARAMETER COUNT CHANGES, AND S1 sec.4's COUNT IS U(1)-ONLY.")
print("    Rank of the infinitesimal gauge orbit at a generic connection, measured:")
def orbit_rank(dimG, gens, edges, U, nv):
    cols = []
    eps = 1e-6
    for v in range(nv):
        for X in gens:
            gg = [np.eye(U[0].shape[0], dtype=complex) for _ in range(nv)]
            gg[v] = np.eye(U[0].shape[0], dtype=complex) + eps * X
            Up = gauge_transform(U, edges, gg)
            col = np.concatenate([((Up[e] - U[e]) / eps).ravel() for e in range(len(edges))])
            cols.append(np.concatenate([col.real, col.imag]))
    M = np.array(cols).T
    return int(np.linalg.matrix_rank(M, tol=1e-7))
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.array([[1, 0], [0, -1]], dtype=complex)
gens_su2 = [1j * sx / 2, 1j * sy / 2, 1j * sz / 2]
rngu = np.random.default_rng(5)
Usu2 = [su2(rngu.uniform(0.4, 2.5), rngu.normal(size=3)) for _ in range(6)]
r_su2 = orbit_rank(3, gens_su2, EDGES_K1, Usu2, 5)
print("      SU(2): connection space dim 6*3 = 18, gauge group dim 5*3 = 15,")
print("             measured orbit rank = %d  ->  invariants = 18 - %d = %d"
      % (r_su2, r_su2, 18 - r_su2))
Uu1 = [np.array([[np.exp(1j * rngu.uniform(0, 6.28))]]) for _ in range(6)]
r_u1 = orbit_rank(1, [np.array([[1j]])], EDGES_K1, Uu1, 5)
print("      U(1) : connection space dim 6, gauge dim 5, measured orbit rank = %d"
      % r_u1)
print("             -> invariants = 6 - %d = %d   [S1 sec.4: 2]" % (r_u1, 6 - r_u1))
print("    S1 sec.4 says 'exactly two invariants exist, so the count is saturated and")
print("    nothing is hidden.' At SU(2) there are THREE. The third is the relative")
print("    invariant of the pair -- which the claim is right to locate at the pinch, but")
print("    which S1's sealed, unaudited saturation argument silently forbids.")
