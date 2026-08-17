# LANE W11-R-CROSS  LEG B -- THE OPERATIVE VARIABLE, NAMED FOR THE EIGHTH TIME AND CHECKED.
#
# The registrar names it THE TRANSPORT CONVENTION.  The steelman names it THE COMPARISON TIME.
# Both names are testable and only one survives.  ONE VARIABLE MOVES IN THIS LEG: the transport
# operator, over the WHOLE root variety, with the clock held at each branch's own loop closure.
#
# B1.  For EVERY unitary U with U^L = M_gamma -- not COR-F's T, not the uniform root D, but every
#      point of the variety -- reading each branch at ITS OWN loop closure returns M_gamma^k
#      IDENTICALLY, because U^L = M is the defining equation.  So the corpus's functional is
#      INVARIANT UNDER THE ENTIRE TRANSPORT AXIS.  Sampled, not asserted.
# B2.  And the invisibility at loop closure is an IDENTITY, not a measurement: the diagonal of
#      M_F is W_F^{1[v in gamma_F]} and of M_C is W_C^{1[v in gamma_C]}, so the FOUR CLASSES ARE
#      THE JOINT LEVEL SETS OF THE TWO OPERATORS' DIAGONALS.  A function of the diagonals depends
#      only on the level sets of the diagonals.  Checked to machine zero against the class sum.
# B3.  The circuit clock is ABSORBING on the transport axis; the edge clock is not.  So the design
#      is not a symmetric 2x2, and "visibility needs BOTH" is true but mis-weighted.
import numpy as np
import xlib as X

rng = np.random.default_rng(20260817)
CASES = [("K1 ", X.K1_LOOP_F, X.K1_LOOP_C, 5, np.array([1.0, 0.37, 0.91, 2 ** 0.5, 0.23, 1.77]),
          np.array([0.40, 0.15, 0.15, 0.15, 0.15])),
         ("B0b", X.B0B_LOOP_F, X.B0B_LOOP_C, 9, rng.uniform(0, 2 * np.pi, 18),
          np.array([.10, .12, .09, .14, .11, .11, .11, .11, .11]) / 1.0)]

print("== B1  THE WHOLE TRANSPORT AXIS MOVES AND THE CORPUS'S FUNCTIONAL DOES NOT ==")
print("   400 random roots per carrier, sampled from the full variety (generic, mostly NON-diagonal),")
print("   read at each branch's OWN loop closure and compared with the corpus's own M_gamma.")
for nm, lf, lc, NV, aa, w in CASES:
    LF, LC = len(lf), len(lc)
    MF, MC = X.M_circuit(lf, aa, NV), X.M_circuit(lc, aa, NV)
    w = w / w.sum()
    states = [np.sqrt(w) + 0j] + X.random_pi_identical(rng, lf, lc, NV, w, k=5)
    pis = [X.pi_of(s, lf, lc, NV) for s in states]
    assert all(np.allclose(pis[0], p, atol=1e-12) for p in pis), "pi not held fixed"
    assert X.arms_differ(*states), "STATE ARMS BYTE-IDENTICAL -- leg void"
    worst_root = worst_fn = worst_sp = 0.0
    ndiag = 0
    for _ in range(400):
        UF = X.random_root(lf, aa, NV, rng, "generic")
        UC = X.random_root(lc, aa, NV, rng, "generic")
        assert X.arms_differ(UF, MF), "ARM COLLAPSE: sampled root equals M_gamma"
        if not np.allclose(UF, np.diag(np.diag(UF)), atol=1e-9):
            ndiag += 1
        worst_root = max(worst_root,
                         np.linalg.norm(np.linalg.matrix_power(UF, LF) - MF),
                         np.linalg.norm(np.linalg.matrix_power(UC, LC) - MC))
        for k in (1, 2, 3):
            v = [abs(X.Z(UF, UC, s, LF * k, LC * k)) for s in states]
            worst_sp = max(worst_sp, max(v) - min(v))
            worst_fn = max(worst_fn, max(abs(X.Z(UF, UC, s, LF * k, LC * k)
                                             - X.Z(MF, MC, s, k, k)) for s in states))
    print(f"  {nm} 400 roots, {ndiag} of them non-diagonal, ||U_F - M_F|| ranges over the variety:")
    print(f"       max ||U^L - M||                                    = {worst_root:.2e}")
    print(f"       max | Z[U at its own loop closure] - Z[corpus M] |  = {worst_fn:.2e}")
    print(f"       max pi-spread over 6 pi-identical states, k<=3      = {worst_sp:.2e}")
print("  -> THE TRANSPORT IS NOT THE OPERATIVE VARIABLE.  It cannot be: U^L = M is the defining")
print("     equation of the whole rival class, and the corpus reads at multiples of L.  COR-F's T")
print("     and the uniform root D are two points of a set on which the corpus's functional is")
print("     CONSTANT.  The registrar's leg B is headed 'ONE VARIABLE MOVED: THE CONVENTION' and")
print("     the variable it moved is the CLOCK.  The steelman is right about the name.")

print("\n== B2  AND AT LOOP CLOSURE THE INVISIBILITY IS AN IDENTITY, NOT A RESULT ==")
for nm, lf, lc, NV, aa, w in CASES:
    MF, MC = X.M_circuit(lf, aa, NV), X.M_circuit(lc, aa, NV)
    dF, dC = np.diag(MF), np.diag(MC)
    cls = X.classes(lf, lc, NV)
    WF, WC = X.holonomy(lf, aa), X.holonomy(lc, aa)
    bad = max(abs(dF[v] - WF ** cls[v][0]) + abs(dC[v] - WC ** cls[v][1]) for v in range(NV))
    w = w / w.sum()
    s = np.sqrt(w) * np.exp(1j * rng.uniform(0, 2 * np.pi, NV))
    p = X.pi_of(s, lf, lc, NV)
    err = 0.0
    for k in range(1, 8):
        direct = X.Z(MF, MC, s, k, k)
        byclass = (p[0] + p[1] * np.conj(WF) ** k + p[2] * WC ** k
                   + p[3] * (np.conj(WF) * WC) ** k)
        err = max(err, abs(direct - byclass))
    print(f"  {nm}  max_v | diag(M_F)_v - W_F^(1[v in gamma_F]) | + | diag(M_C)_v - W_C^(1[v in gamma_C]) |"
          f" = {bad:.2e}")
    print(f"       max_k | <M_F^k s, M_C^k s> - (p00 + p10 conj(W_F)^k + p01 W_C^k + p11 (conj(W_F)W_C)^k) |"
          f" = {err:.2e}")
print("  -> THE FOUR CLASSES ARE THE JOINT LEVEL SETS OF THE TWO OPERATORS' DIAGONALS, BY THE")
print("     DEFINITION OF BOTH.  'The functional depends only on pi' is therefore the statement")
print("     that a function of two diagonals depends only on the level sets of those diagonals.")
print("     It is an identity in the operator, true before any carrier is chosen, and it is what")
print("     the registrar's leg B1 (spread 1e-16 under the circuit convention) measures.")
print("     That leg is a CONTROL THAT COULD NOT HAVE FAILED -- the registrar says so itself.")

print("\n== B3  THE CLOCK AXIS IS ABSORBING; THE TRANSPORT AXIS IS NOT.  THE 2x2 IS NOT SQUARE ==")
print(f"  {'carrier':<6}{'clock':<10}{'transport':<28}{'worst pi-spread':>18}")
for nm, lf, lc, NV, aa, w in CASES:
    LF, LC = len(lf), len(lc)
    w = w / w.sum()
    states = [np.sqrt(w) + 0j] + X.random_pi_identical(rng, lf, lc, NV, w, k=3)
    MF, MC = X.M_circuit(lf, aa, NV), X.M_circuit(lc, aa, NV)
    TF, TC = X.T_edge(lf, aa, NV), X.T_edge(lc, aa, NV)
    DF, DC = X.D_uniform(lf, aa, NV), X.D_uniform(lc, aa, NV)
    roots = [("COR-F edge tick T", TF, TC), ("uniform root D", DF, DC)]
    for i in range(3):
        roots.append((f"random root #{i+1} of the variety",
                      X.random_root(lf, aa, NV, rng, "generic"),
                      X.random_root(lc, aa, NV, rng, "generic")))
    for rn, uF, uC in roots:
        for cn, (kF, kC) in (("EDGE   ", (1, 1)), ("CIRCUIT", (LF, LC))):
            wsp = 0.0
            for n in range(1, 7):
                v = [abs(X.Z(uF, uC, s, kF * n, kC * n)) for s in states]
                wsp = max(wsp, max(v) - min(v))
            print(f"  {nm:<6}{cn:<10}{rn:<28}{wsp:>18.2e}")
print("  -> EVERY row with the CIRCUIT clock is invisible, for EVERY transport, by the identity of")
print("     B1.  Only the EDGE clock has a transport axis at all.  So the steelman's 2x2 conclusion")
print("     -- 'visibility needs BOTH a non-class-constant root AND an edge clock' -- is true as a")
print("     conjunction and MIS-WEIGHTED as an attribution: one factor is absorbing and the other")
print("     is not.  The clock alone is SUFFICIENT for invisibility.  Nothing about the transport is.")

print("\n== B4  AND THE FIBRE-WISE-NESS MECHANISM (W-06's CORRECTED N4), QUANTIFIED ==")
print("  The steelman's item 6 rebuts 'M_gamma is in the gauge group' with a JOINT (connection,state)")
print("  invariant that M does move.  That is true and it answers a different question.  The claim")
print("  that carries Reading B is about the CARRIER's own observables: M_gamma is diagonal, so the")
print("  two branch density matrices have IDENTICAL diagonals and no diagonal observable separates")
print("  them.  COR-F's T does separate them -- which is exactly what COR-F's sealed exhibit says.")
for nm, lf, lc, NV, aa, w in CASES:
    w = w / w.sum()
    s = np.sqrt(w) * np.exp(1j * rng.uniform(0, 2 * np.pi, NV))
    MF, MC = X.M_circuit(lf, aa, NV), X.M_circuit(lc, aa, NV)
    TF, TC = X.T_edge(lf, aa, NV), X.T_edge(lc, aa, NV)
    dM = max(np.max(np.abs(np.abs(np.linalg.matrix_power(MF, k) @ s) ** 2
                           - np.abs(np.linalg.matrix_power(MC, k) @ s) ** 2)) for k in range(1, 5))
    dT = max(np.max(np.abs(np.abs(np.linalg.matrix_power(TF, n) @ s) ** 2
                           - np.abs(np.linalg.matrix_power(TC, n) @ s) ** 2)) for n in range(1, 5))
    print(f"  {nm}  max_k || diag(rho_F) - diag(rho_C) ||_inf  under M_gamma = {dM:.2e}"
          f"     under COR-F's T = {dT:.3f}")
print("  -> under the corpus's convention the carrier's own diagonal record is EMPTY, identically,")
print("     at every k, on both carriers.  That is Reading B's mechanism in one number, and it is")
print("     untouched by the steelman's item 6.")
