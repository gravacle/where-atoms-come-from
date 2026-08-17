# LANE W11-R-CROSS  LEG E -- THE LOAD-BEARING JOINT.  THE STEELMAN NAMED IT ITSELF:
#   "MY LOAD-BEARING MOVE HAS A SEAM AND I NAME IT FIRST ... THAT MAPPING IS MINE.  Attack it there
#    first."   This leg attacks it there.
#
# THE ARGUMENT UNDER ATTACK (steelman leg D3/D4, verbatim):
#   "SELECTION.  The pre-registered trivial-connection contact point (FOUNDING_DESIGN:117-118,
#    S2 CHOICE LEDGER C4) requires U(a=0) = I.  { V unitary : V^L = I } is a disjoint union of
#    conjugacy-class manifolds; the one containing I is the single point {I} ... THEREFORE U = D."
#
# E1.  "REQUIRES U(a=0) = I" IS FALSE.  The contact point is a condition on the FORMATION VERDICT,
#      not on the operator.  What it requires is NO FORMATION at a = 0, i.e. |Z_n| = 1 there, i.e.
#      U_F(0) = mu U_C(0) -- the two branch operators must AGREE at zero field, not equal I.
#      Exhibited: a continuous family, unitary, U^L = M_gamma EXACTLY, identity off nothing it
#      should not touch, with U_F(0) = U_C(0) at distance sqrt(3) from I, |Z_n| = 1 EXACTLY at the
#      trivial connection for every state -- and pi-BROKEN at the generic connection.
#      The component argument does not reach it, and the conclusion U = D does not follow.
# E2.  WHAT DOES EXCLUDE IT IS GAUGE COVARIANCE -- i.e. COR-J, the sealed UNDECLARED premise
#      "the record must be gauge-invariant" (S3_THE_CROSSING_AUDIT_V001.md:2.5, correction table),
#      which COR-J says is load-bearing, applied asymmetrically, and belongs in the CHOICE LEDGER,
#      and which the corpus never entered.  So "selected, not stipulated" trades one unledgered
#      premise for another.
# E3.  IMPOSE COVARIANCE AND THE CONCLUSION SURVIVES -- on a ground the steelman did not give.
#      Under covariance the root variety on K1 is NOT a continuum: it is 27 points per loop
#      (3 branch labels on each of T's three eigenvectors), T and D among them.  EXHAUSTIVE test of
#      all 27 x 27 pairs: passes-trivial-limit and is-pi-only coincide exactly.
# E4.  AND THE CRITERION ADMITS EXACTLY ONE CATEGORY, NOT ONE CANDIDATE.  Any two branch dynamics
#      that MOVE amplitude along two distinct loops disagree at zero field, so they fire there.
#      Checked on COR-F's tick, on random motions, and on the corpus's OWN ledgered alternative
#      (S2 CHOICE LEDGER A2's "a real parameter t with a Hamiltonian" = S2's addition (i), the
#      magnetic Laplacian).  The contact point does not choose among transports.  It excludes
#      TRANSPORT AS MOTION, leaving fibre-wise-ness -- which is W-06's corrected N4 mechanism and
#      Reading B's mechanism.
import numpy as np
import xlib as X

rng = np.random.default_rng(20260817)
lf, lc, NV, ne = X.K1_LOOP_F, X.K1_LOOP_C, 5, 6
a = np.array([1.0, 0.37, 0.91, 2 ** 0.5, 0.23, 1.77])
a0 = np.zeros(ne)
LF = LC = 3
MF, MC = X.M_circuit(lf, a, NV), X.M_circuit(lc, a, NV)
CLS = X.classes(lf, lc, NV)
w_base = np.array([0.40, 0.15, 0.15, 0.15, 0.15])
STATES = [np.sqrt(w_base) + 0j] + X.random_pi_identical(rng, lf, lc, NV, w_base, k=20)
assert X.arms_differ(*STATES[:6]), "STATE ARMS BYTE-IDENTICAL -- leg void"
RAND_STATES = []
for _ in range(50):
    s = rng.normal(size=NV) + 1j * rng.normal(size=NV)
    RAND_STATES.append(s / np.linalg.norm(s))


def pi_spread(uF, uC, nmax=9):
    w = 0.0
    for n in range(1, nmax + 1):
        v = [abs(X.Z(uF, uC, s, n, n)) for s in STATES]
        w = max(w, max(v) - min(v))
    return w


def forms_at(uF, uC, nmax=9, states=None):
    """the trivial-limit test, run on the VERDICT: does |Z_n| drop below 1 for any state?"""
    states = states if states is not None else RAND_STATES
    return 1.0 - min(abs(X.Z(uF, uC, s, n, n)) for s in states for n in range(1, nmax + 1))


def lam(loop, aa):
    Lg = len(loop)
    w = np.exp(1j * np.angle(X.holonomy(loop, aa)) / Lg)
    L_ = np.eye(NV, dtype=complex)
    for v in X.loop_vertices(loop):
        L_[v, v] = w
    return L_


print("== E1  THE CONTACT POINT DOES NOT REQUIRE U(0) = I, AND WITH THE TRUE REQUIREMENT THE")
print("       CONCLUSION FAILS.  AN EXPLICIT ADMITTED, INCIDENCE-VISIBLE FAMILY ==")
om = np.exp(2j * np.pi / 3)
R = np.diag([om, 1.0, 1.0, 1.0, 1.0]).astype(complex)      # order 3, class-block diagonal, != I
# G(a): a continuous unitary mixing v0 and v3 (both on gamma_C), G(0) = I
K = np.zeros((NV, NV), dtype=complex)
K[0, 3] = K[3, 0] = 1.0


def family(aa, kappa=1.0):
    t = kappa * float(np.linalg.norm(aa))
    ev, EV = np.linalg.eigh(K)
    G = EV @ np.diag(np.exp(1j * t * ev)) @ EV.conj().T
    UF = lam(lf, aa) @ R
    UC = lam(lc, aa) @ (G @ R @ G.conj().T)
    return UF, UC


UF1, UC1 = family(a)
UF0, UC0 = family(a0)
print(f"  U_F^3 = M_F ?  ||.|| = {np.linalg.norm(np.linalg.matrix_power(UF1,3)-MF):.2e}"
      f"     U_C^3 = M_C ?  ||.|| = {np.linalg.norm(np.linalg.matrix_power(UC1,3)-MC):.2e}")
print(f"  unitary ?  ||U_F*U_F - I|| = {np.linalg.norm(UF1.conj().T@UF1-np.eye(NV)):.2e}"
      f"   ||U_C*U_C - I|| = {np.linalg.norm(UC1.conj().T@UC1-np.eye(NV)):.2e}")
print(f"  AT THE TRIVIAL CONNECTION:  ||U_F(0) - I|| = {np.linalg.norm(UF0-np.eye(NV)):.6f}"
      f"   ||U_C(0) - I|| = {np.linalg.norm(UC0-np.eye(NV)):.6f}"
      f"   ||U_F(0) - U_C(0)|| = {np.linalg.norm(UF0-UC0):.2e}")
print(f"     -> NOT the identity, and yet:  max_(50 random states, n<=9) (1 - |Z_n|) at a = 0 = "
      f"{forms_at(UF0,UC0):.2e}   NO FORMATION, EXACTLY.")
print(f"  AT THE GENERIC CONNECTION:  pi-spread over 21 pi-identical states, n <= 9 = "
      f"{pi_spread(UF1,UC1):.3e}   INCIDENCE VISIBLE.")
print(f"  continuity: ||U_F(t.a) - U_F(0)|| for t = 1, 0.1, 0.01, 0.001:  " +
      " ".join(f"{np.linalg.norm(family(t*a)[0]-UF0):.4f}" for t in (1, .1, .01, .001)))
print("  -> THE STEELMAN'S SELECTION PREMISE IS FALSE AS STATED.  The pre-registered contact point")
print("     admits this family, and this family sees the incidence.  Nothing in FOUNDING_DESIGN")
print("     :117-118 or S2:583 says the operator must be the identity at zero field; they say the")
print("     ANSWER must be the known trivial answer, and it is.")

print("\n== E2  WHAT ACTUALLY EXCLUDES IT IS COR-J's UNDECLARED PREMISE, NOT THE CONTACT POINT ==")
TF, TC = X.T_edge(lf, a, NV), X.T_edge(lc, a, NV)
DF, DC = X.D_uniform(lf, a, NV), X.D_uniform(lc, a, NV)
s = RAND_STATES[0]
rows = []
for nm, mk in (("COR-F edge tick T", lambda aa: (X.T_edge(lf, aa, NV), X.T_edge(lc, aa, NV))),
               ("uniform root D", lambda aa: (X.D_uniform(lf, aa, NV), X.D_uniform(lc, aa, NV))),
               ("corpus M_gamma", lambda aa: (X.M_circuit(lf, aa, NV), X.M_circuit(lc, aa, NV))),
               ("E1's admitted family", lambda aa: family(aa))):
    worst = 0.0
    for _ in range(400):
        th = rng.uniform(0, 2 * np.pi, NV)
        g = np.diag(np.exp(1j * th))
        a2 = X.gauge_apply(a, th, X.K1_EDGES)
        uF, uC = mk(a)
        vF, vC = mk(a2)
        for n in (1, 2, 3):
            worst = max(worst, abs(X.Z(uF, uC, s, n, n) - X.Z(vF, vC, g @ s, n, n)))
    rows.append((nm, worst))
    print(f"  {nm:<24} max |Z_n(a^g, g.s) - Z_n(a,s)| over 400 gauge draws, n<=3 = {worst:.2e}")
print("  -> the admitted family is NOT gauge-invariant.  Excluding it needs the premise")
print("     'the record must be gauge-invariant' -- COR-J, sealed as UNDECLARED, LOAD-BEARING and")
print("     'part theorem and part bookkeeping convention', with the instruction 'Add it to the")
print("     CHOICE LEDGER'.  It was never added.  The steelman's selection therefore rests on an")
print("     unledgered premise of the corpus's own, which is exactly Reading B's complaint one")
print("     level up.  I do not claim the premise is wrong -- I claim it is the thing doing the work.")

print("\n== E3  IMPOSE COVARIANCE AND RE-RUN THE SELECTION, EXHAUSTIVELY ==")
print("  Under gauge covariance every operator supported on a 3-cycle is a POLYNOMIAL IN T (each")
print("  entry must carry the Wilson line of the path joining its two vertices), so the covariant")
print("  roots of M_gamma are the 3^3 = 27 branch choices u_j = W^{1/3} omega^{m_j} on T's three")
print("  eigenvectors.  T is the branch (m_0,m_1,m_2) = (0,1,2) and D is (0,0,0).  DISCRETE, not a")
print("  continuum.  All 27 x 27 pairs tested at once for both properties.")


def cov_roots(loop, aa):
    """all 27 gauge-covariant unitaries U with U^L = M_gamma(aa), identity off the loop."""
    on = sorted(X.loop_vertices(loop))
    Tl = X.T_edge(loop, aa, NV)
    sub = Tl[np.ix_(on, on)]
    ev, EV = np.linalg.eig(sub)
    w = np.exp(1j * np.angle(X.holonomy(loop, aa)) / len(loop))
    # branch index j of each eigenvector: eigenvalue = w * omega^j
    idx = [int(np.round(np.angle(e / w) * 3 / (2 * np.pi))) % 3 for e in ev]
    out = []
    for m in range(27):
        ms = [(m // 3 ** i) % 3 for i in range(3)]
        lamvec = np.array([w * om ** ms[idx[k]] for k in range(3)])
        B = EV @ np.diag(lamvec) @ np.linalg.inv(EV)
        U = np.eye(NV, dtype=complex)
        for i, v in enumerate(on):
            for j, u in enumerate(on):
                U[v, u] = B[i, j]
        out.append((tuple(ms), U))
    return out


covF, covC = cov_roots(lf, a), cov_roots(lc, a)
covF0, covC0 = cov_roots(lf, a0), cov_roots(lc, a0)
# sanity: T and D are in the list
dmin_T = min(np.linalg.norm(U - TF) for _, U in covF)
dmin_D = min(np.linalg.norm(U - DF) for _, U in covF)
print(f"  is COR-F's T among the 27?  min||U - T_F|| = {dmin_T:.2e}      "
      f"is D among the 27?  min||U - D_F|| = {dmin_D:.2e}")
bad_root = max(np.linalg.norm(np.linalg.matrix_power(U, 3) - MF) for _, U in covF)
print(f"  max ||U^3 - M_F|| over all 27 = {bad_root:.2e}")
tab = {}
for mF, UF in covF:
    UF0_ = dict(covF0)[mF]
    for mC, UC in covC:
        UC0_ = dict(covC0)[mC]
        triv = forms_at(UF0_, UC0_, nmax=6, states=RAND_STATES[:12]) < 1e-9   # passes contact point
        vis = pi_spread(UF, UC, nmax=6) > 1e-9                                # sees incidence
        tab[(triv, vis)] = tab.get((triv, vis), 0) + 1
print(f"  all {len(covF)*len(covC)} covariant pairs classified:")
for k in sorted(tab):
    print(f"     passes trivial-limit = {str(k[0]):<5}   incidence visible = {str(k[1]):<5}"
          f"   count = {tab[k]}")
survivors = [(mF, mC) for mF, UF in covF for mC, UC in covC
             if forms_at(dict(covF0)[mF], dict(covC0)[mC], 6, RAND_STATES[:12]) < 1e-9]
print(f"  the pairs that pass the contact point: {survivors}   "
      f"(branch (0,0,0) is D; COR-F's T is (0,1,2))")
print("  -> UNDER COVARIANCE THE STEELMAN'S CONCLUSION HOLDS, AND IS EXACT RATHER THAN SAMPLED:")
print("     the contact point admits exactly the uniform root and nothing else, and the admitted")
print("     cell is pi-only.  Its stated ARGUMENT (U(0) = I, plus a component argument over the")
print("     FULL unitary root variety) is not the reason -- E1 refutes that argument outright.")

print("\n== E4  AND WHAT THE CONTACT POINT ACTUALLY DOES IS EXCLUDE MOTION, NOT SELECT A TRANSPORT ==")
print("  LEMMA (checked): |Z_n| = 1 for every state at a = 0  <=>  U_F(0) = mu U_C(0).  Two")
print("  dynamics that MOVE amplitude along two DIFFERENT loops cannot agree at zero field.")
print(f"  {'candidate branch dynamics':<44}{'||U_F(0)-U_C(0)||':>20}{'1-min|Z| at a=0':>18}")
cands = [("COR-F edge tick T", X.T_edge(lf, a0, NV), X.T_edge(lc, a0, NV))]
for i in range(3):
    RF = X.random_root(lf, a0, NV, rng, "generic")
    RC = X.random_root(lc, a0, NV, rng, "generic")
    cands.append((f"random unitary root #{i+1} at a=0", RF, RC))
# the corpus's OWN ledgered alternative: S2 CHOICE LEDGER A2, "a real parameter t with a
# Hamiltonian" = S2's addition (i), the magnetic Laplacian Delta_A = D - H, restricted to each loop
def maglap(loop, aa):
    H = np.zeros((NV, NV), dtype=complex)
    deg = np.zeros(NV)
    for (src, dst, e, sg) in loop:
        ph = np.exp(1j * aa[e]) if sg > 0 else np.exp(-1j * aa[e])
        H[dst, src] += ph
        H[src, dst] += np.conj(ph)
        deg[src] += 1
        deg[dst] += 1
    return np.diag(deg) - H
def expm_h(A, t):
    ev, EV = np.linalg.eigh(A)
    return EV @ np.diag(np.exp(-1j * t * ev)) @ EV.conj().T
HF0, HC0 = maglap(lf, a0), maglap(lc, a0)
cands.append(("S2 addition (i): magnetic Laplacian, t=0.7",
              expm_h(HF0, 0.7), expm_h(HC0, 0.7)))
cands.append(("S2 addition (i): magnetic Laplacian, t=1.3",
              expm_h(HF0, 1.3), expm_h(HC0, 1.3)))
cands.append(("uniform root D (a fibre-wise phase, not a motion)",
              X.D_uniform(lf, a0, NV), X.D_uniform(lc, a0, NV)))
for nm, uF, uC in cands:
    ph = np.vdot(uC.reshape(-1), uF.reshape(-1))
    ph = ph / abs(ph) if abs(ph) > 1e-12 else 1.0
    print(f"  {nm:<44}{np.linalg.norm(uF-ph*uC):>20.4f}{forms_at(uF,uC,6,RAND_STATES[:20]):>18.3e}")
print("  -> every MOTION fires at zero field, including the corpus's own ledgered alternative.")
print("     Only the non-motion passes.  So the pre-registered contact point does not adjudicate")
print("     BETWEEN transport conventions; it eliminates the whole category of dynamics that move")
print("     amplitude, and what it leaves standing is FIBRE-WISE-NESS -- W-06's corrected N4, i.e.")
print("     precisely the stipulation Reading B names.  A criterion whose admitted set is exactly")
print("     the stipulated class does not convert the stipulation into a finding; it justifies it.")
