# LANE W11-R-CROSS  LEG C -- DIFF THE STEELMAN'S ARMS.
#
# W-08's isolation audit: "the commonest FATAL defect is not 'two variables moved' -- it is ZERO
# variables moved: a control whose two arms are byte-identical, reported as a confirmation.
# A ledger cannot catch that.  DIFF YOUR ARMS."  The steelman diffs its STATE arms and its
# TRANSPORT arms (T vs D) and reports both differ.  Both do.  It does NOT diff the four CELLS of
# its 2x2, and three of them are the same object.
#
# C1.  (T, CIRCUIT clock) and (D, CIRCUIT clock) are THE SAME EVOLUTION, exactly: both are M_gamma
#      raised to the circuit count.  ||T^L - D^L|| is machine zero.  Two cells, one object.
# C2.  (D, EDGE clock) is the CORPUS'S OWN FUNCTIONAL with the connection relabelled
#      (W_F, W_C) -> (W_F^{1/L_F}, W_C^{1/L_C}).  Verified against the closed class formula.
#      So it is the corpus's convention at another point of the corpus's own parameter space --
#      a third cell, same object.
# C3.  What is left is ONE cell.  The 2x2 has one non-degenerate entry and three restatements.
import numpy as np
import xlib as X

rng = np.random.default_rng(20260817)
CASES = [("K1 ", X.K1_LOOP_F, X.K1_LOOP_C, 5, np.array([1.0, 0.37, 0.91, 2 ** 0.5, 0.23, 1.77]),
          np.array([0.40, 0.15, 0.15, 0.15, 0.15])),
         ("B0b", X.B0B_LOOP_F, X.B0B_LOOP_C, 9, rng.uniform(0, 2 * np.pi, 18),
          np.array([.10, .12, .09, .14, .11, .11, .11, .11, .11]))]

print("== C1  THE TWO 'CIRCUIT CLOCK' CELLS OF THE 2x2 ARE THE SAME OPERATOR ==")
print(f"  {'carrier':<6}{'||T_F-D_F||':>14}{'arms differ':>13}{'||T_F^L - D_F^L||':>20}"
      f"{'||T_F^L - M_F||':>18}{'||D_F^L - M_F||':>18}")
for nm, lf, lc, NV, aa, w in CASES:
    LF = len(lf)
    TF, DF, MF = X.T_edge(lf, aa, NV), X.D_uniform(lf, aa, NV), X.M_circuit(lf, aa, NV)
    tl = np.linalg.matrix_power(TF, LF)
    dl = np.linalg.matrix_power(DF, LF)
    print(f"  {nm:<6}{np.linalg.norm(TF-DF):>14.4f}{str(X.arms_differ(TF,DF)):>13}"
          f"{np.linalg.norm(tl-dl):>20.2e}{np.linalg.norm(tl-MF):>18.2e}{np.linalg.norm(dl-MF):>18.2e}")
print("  -> the TRANSPORT arms differ in bytes (2.4495 / 2.8284 apart, as the steelman reports) and")
print("     their L-th POWERS -- the only thing the circuit clock ever evaluates -- are identical to")
print("     machine zero.  Under the circuit clock the transport axis of the 2x2 MOVES NOTHING.")

print("\n== C2  THE (D, EDGE) CELL IS THE CORPUS'S FUNCTIONAL AT RELABELLED HOLONOMIES ==")
for nm, lf, lc, NV, aa, w in CASES:
    LF, LC = len(lf), len(lc)
    DF, DC = X.D_uniform(lf, aa, NV), X.D_uniform(lc, aa, NV)
    wF = np.exp(1j * np.angle(X.holonomy(lf, aa)) / LF)
    wC = np.exp(1j * np.angle(X.holonomy(lc, aa)) / LC)
    ww = w / w.sum()
    err = 0.0
    for s in [np.sqrt(ww) + 0j] + X.random_pi_identical(rng, lf, lc, NV, ww, k=4):
        p = X.pi_of(s, lf, lc, NV)
        for n in range(1, 13):
            direct = X.Z(DF, DC, s, n, n)
            corpus = (p[0] + p[1] * np.conj(wF) ** n + p[2] * wC ** n
                      + p[3] * (np.conj(wF) * wC) ** n)
            err = max(err, abs(direct - corpus))
    print(f"  {nm}  max_n,s | Z^D_n(a) - [p00 + p10 conj(w_F)^n + p01 w_C^n + p11 (conj(w_F)w_C)^n] |"
          f" = {err:.2e}")
    print(f"       with (w_F, w_C) = (W_F^(1/{LF}), W_C^(1/{LC})) = ({wF:.6f}, {wC:.6f})")
print("  -> the (D, EDGE) cell IS W-01's functional, written at the L-th-root connection.  It is")
print("     not an independent arm: it is the corpus's own object with its argument relabelled.")

print("\n== C3  THE SCORECARD OF THE 2x2 ==")
print("  cell                 what it actually is                                     independent?")
print("  (T, CIRCUIT)         M_gamma^k, the corpus's own operator (T^L = M)           NO")
print("  (D, CIRCUIT)         M_gamma^k, the SAME matrix as the cell above             NO -- same object")
print("  (D, EDGE)            the corpus's functional at (W^{1/L_F}, W^{1/L_C})        NO -- relabelled")
print("  (T, EDGE)            COR-F's rival                                            YES -- the only one")
print("  Three of the four cells are the corpus's convention.  The steelman's conclusion from this")
print("  design -- 'neither the root alone nor the clock alone is the convention' -- rests on three")
print("  cells that could not have shown anything else.  ZERO VARIABLES MOVED in two of them.")

print("\n== C4  AND THE SAME COLLAPSE ONE LEVEL UP: 'D IS AN EQUALLY NATURAL EDGE TICK' IS NOT A ")
print("        MOTION AT ALL, AS THE STEELMAN'S OWN FLAG SAYS.  QUANTIFIED HERE ==")
for nm, lf, lc, NV, aa, w in CASES:
    NVv = NV
    TF, DF = X.T_edge(lf, aa, NV), X.D_uniform(lf, aa, NV)
    e = np.zeros(NV, dtype=complex)
    v0 = sorted(X.loop_vertices(lf))[0]
    e[v0] = 1.0
    tv = TF @ e
    dv = DF @ e
    print(f"  {nm}  T_F moves the excitation from v{v0} to v{int(np.argmax(np.abs(tv)))} "
          f"(|overlap with start| = {abs(tv[v0]):.3f});  "
          f"D_F leaves it at v{int(np.argmax(np.abs(dv)))} (|overlap| = {abs(dv[v0]):.3f})")
print("  -> D transports nothing.  The steelman's rival to COR-F's transport is a phase clock, and")
print("     the brief's decisive question -- 'if a different, EQUALLY NATURAL EDGE TICK restores")
print("     invisibility, Reading B falls' -- is therefore NOT answered by D.  D is not a tick.")
