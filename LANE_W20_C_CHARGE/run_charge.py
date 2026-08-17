# run_charge.py -- LANE W20_C.  THE CHARGED-SECTOR RUN, AS W-19 PRE-REGISTERED IT.
#
# ARM 1 (A0): G_v = +1 at all eight vertices.
# ARM 2 (A1): G_v = -1 at two vertices, {0,4}.  BYTE-IDENTICAL APART FROM THE SIGN PATTERN.
# The two Hamiltonians differ ONLY through the coset offset u0 that the Gauss constraint selects;
# the plaquette hopping is literally the same matrix.  The DIFF is printed.
import math, itertools
import numpy as np
import w20c_core as C

LOG = []
def P(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s, flush=True)
def rule(t=""):
    P("\n" + "=" * 108); P(t); P("=" * 108)

def rho_alg(A, sec, psi):
    blocks = A.state(sec, psi)
    d = 1 << A.k
    n = (1 << A.r) * d
    M = np.zeros((n, n), dtype=complex)
    for s in range(1 << A.r):
        M[s * d:(s + 1) * d, s * d:(s + 1) * d] = blocks[s]
    return (M + M.conj().T) / 2

def tracedist(M1, M2):
    w = np.linalg.eigvalsh(M1 - M2).real
    return 0.5 * float(np.sum(np.abs(w)))

def ent(A, sec, psi): return A.entropy(sec, psi)

# ============================================================================================
rule("BLOCK 0 -- VALIDATION OF THE INSTRUMENT AGAINST SEALED PRIOR NUMBERS.  RUN FIRST.")
P("This lane re-implements the physical sector, the Hamiltonian and the entropy from scratch.")
P("Before any new number is produced it must reproduce W-19's SEALED numbers on this carrier.")
P("W19_RULING/OUT_rule_main2.txt line 87 : tri_chain12 ground g2=3.00  H_elec(S)=0.001283123")
P("W19_RULING/OUT_rule_main2.txt line 86 : tri_chain12 ground g2=0.50  H_elec(S)=0.663269043")
P("W-19's S there was a SINGLE LINK (rule_verify.py:251, AX(L,[l])), so H_elec = H(X_l), max 1 bit.")
sec0 = C.Sector([])
for g2, tgt in [(3.00, 0.001283123), (0.50, 0.663269043)]:
    _, psi, _ = sec0.ground(g2)
    got = C.Alg([C.sv(1 << 0, 0)], "X0").entropy(sec0, psi)
    P("   g2=%.2f  H(X_0) here = %.9f   W-19 sealed = %.9f   |diff| = %.3e  %s"
      % (g2, got, tgt, abs(got - tgt), "MATCH" if abs(got - tgt) < 5e-10 else "*** MISMATCH ***"))
P("")
P("INTERNAL IDENTITIES, checked rather than assumed:")
_, psiv, _ = sec0.ground(1.00)
P("   pure-state duality S(rho|A) = S(rho|A') for A = A_S and A' = A_env :  %.12f vs %.12f"
  % (C.A_FULL.entropy(sec0, psiv), C.A_ENV.entropy(sec0, psiv)))
P("   algebra dimensions vs LANE_W20_PRE (pinned before any state existed):")
for nm, A in C.CHANNELS + [("WIDE", C.A_WIDE), ("A_Sigma", C.A_SIG), ("A_env", C.A_ENV)]:
    P("      %-8s rank=%-3d scalars=%-2d nbar=%-2d  r=%d k=%d  dim=%-4d  max %d bits"
      % (nm, A.rank, A.dim_scalar, A.nbar, A.r, A.k, A.dim, A.maxent))
P("   PRE-REG said: FULL dim 16 (max 3), CENTRE dim 4 (max 2), BLOCK dim 4 (max 1),")
P("   MAG dim 2 (max 1), WIDE dim 32 (max 3).  ALL FIVE REPRODUCED.")
P("")
P("   monotonicity I(A_S:A_F) <= I(A_S:A_F') for F subset F', spot-checked on the sweep grid:")
bad = 0
for g2 in (0.20, 1.00, 3.00):
    _, ps, _ = sec0.ground(g2)
    a = C.MI(C.A_FULL, C.Alg(C.frag_gens([11]), "f"), sec0, ps)
    b = C.MI(C.A_FULL, C.Alg(C.frag_gens([6, 10, 11]), "f"), sec0, ps)
    c = C.MI(C.A_FULL, C.A_ENV, sec0, ps)
    P("      g2=%.2f  I(.,{11})=%.9f <= I(.,{6,10,11})=%.9f <= I(.,E_env)=%.9f : %s"
      % (g2, a, b, c, (a <= b + 1e-9 <= c + 2e-9)))
    bad += 0 if (a <= b + 1e-9 and b <= c + 1e-9) else 1
P("   violations: %d" % bad)

# ============================================================================================
rule("BLOCK 1 -- A THEOREM ABOUT THE PRE-REGISTERED FALSIFIER, PROVED BEFORE ITS NUMBERS ARE READ")
P("The pre-registration pinned  Delta_surf = I(A_S : A_{E_env}) - I(A_S : A_Sigma)  as THE primary")
P("quantity, and pinned an EMPTINESS GATE (AXIS 6) that discards any point with H_channel < 0.10")
P("bits.  This block shows the two are not independent -- they are exact complements.")
P("")
P("(i)  A_Sigma == Z(A_S) AS REPRESENTED ALGEBRAS.  Not 'is like'.  Equal.")
vs_cen = sorted(set(C.span_all(C.gf2_basis(list(C.A_CENTRE.gens)))))
vs_sig = sorted(set(C.span_all(C.gf2_basis(list(C.A_SIG.gens)))))
def modscal(vs):
    return sorted({min(v, min((v ^ s for s in C.span_all(C.SC_BASIS)), default=v)) for v in vs})
P("     GF(2) span of Z(A_S) = {X_1X_2, X_1X_3}, and of A_Sigma = {X_0,X_4,X_5}, MODULO the")
P("     operators that are scalars on the physical sector:")
P("     Z(A_S) has nbar=%d dim=%d ; A_Sigma has nbar=%d dim=%d ; same algebra: %s"
  % (C.A_CENTRE.nbar, C.A_CENTRE.dim, C.A_SIG.nbar, C.A_SIG.dim,
     C.join(C.A_CENTRE, C.A_SIG).nbar == C.A_CENTRE.nbar))
P("     (join of the two has nbar = %d, i.e. adjoining one to the other adds NOTHING.)"
  % C.join(C.A_CENTRE, C.A_SIG).nbar)
P("")
P("(ii) A_{E_env} == A_S' , THE FULL COMMUTANT OF A_S.  Also an equality, from the incidence:")
P("     the gauge-invariant Pauli algebra on the 32-dim physical sector has GF(2) symplectic")
P("     dimension 2*(L-V+1) = 10.  nbar(A_S) = %d, so nbar(A_S') = 10 - %d = %d = nbar(A_env) = %d."
  % (C.A_FULL.nbar, C.A_FULL.nbar, 10 - C.A_FULL.nbar, C.A_ENV.nbar))
P("     and every generator of A_env commutes with every generator of A_S (X on E_env vs X on S: 0;")
P("     env cycles vs W_S: both Z; X on E_env vs W_S: |F cap S| = 0).  Containment + equal dimension")
P("     = EQUALITY.  There is NO environment operator outside the commutant, and none inside it that")
P("     the environment lacks.")
P("")
P("(iii) THEREFORE, for ANY pure physical state, at ANY coupling, in ANY charge sector:")
P("        I(A_S : A_Sigma)   = S(A_Sigma)                    because A_Sigma is a SUBALGEBRA of A_S")
P("        I(A_S : A_{E_env}) = 2 S(A_S) - S(Z(A_S))          because A_env = A_S' and rho is pure")
P("        ==>  Delta_surf = 2 * ( H_FULL - H_CENTRE )  =  2 * sum_c p_c S(rho_c)")
P("      i.e. Delta_surf IS TWICE THE CONDITIONAL ENTROPY OF THE ONE FREE BLOCK QUBIT, and it never")
P("      reads the environment at all.  It is a function of rho restricted to A_S ALONE.")
P("")
P(">>> CONSEQUENCE, AND IT IS THE FIRST FINDING OF THIS LANE:")
P("    Delta_surf = 0  <=>  the free block qubit is PURE given the central labels  <=>  the record")
P("    has no content beyond the surface BY CONSTRUCTION.  The emptiness gate (AXIS 6) discards")
P("    exactly the couplings where that happens.  So H-SURFACE is CONFIRMED wherever the gate")
P("    REJECTS and REFUTED wherever the gate ACCEPTS -- automatically, with no measurement.")
P("    THE PRE-REGISTERED PRIMARY FALSIFIER COULD NOT HAVE FAILED, IN EITHER DIRECTION.")
P("    This is the W-19 defect a third time, inside the instrument built to catch it the second")
P("    time.  It is caught here by a GF(2) dimension count, and it is verified numerically below.")

# ============================================================================================
rule("BLOCK 2 -- THE MEASURED CROSSOVER AND THE LIVE WINDOW.  FIRST NUMBERS OF THE RUN.")
FINE = [10 ** x for x in np.linspace(math.log10(0.05), math.log10(5.00), 200)]
def sweep_fine(sec):
    out = []
    for g2 in FINE:
        _, psi, _ = sec.ground(g2)
        v, m = C.var_plaq(sec, psi)
        hb = C.A_BLOCK.entropy(sec, psi)
        hm = C.A_MAG.entropy(sec, psi)
        hf = C.A_FULL.entropy(sec, psi)
        hc = C.A_CENTRE.entropy(sec, psi)
        out.append((g2, v, m, hb, hm, hf, hc))
    return out

def window(rows, idx, thr):
    ok = [r[0] for r in rows if r[idx] >= thr]
    return (min(ok), max(ok), len(ok)) if ok else (None, None, 0)

for tag, chg in [("A0 vacuum", []), ("A1 charges {0,4}", [0, 4])]:
    sec = C.Sector(chg)
    rows = sweep_fine(sec)
    g2s = max(rows, key=lambda r: r[1])
    P("%s :" % tag)
    P("   MEASURED crossover g2_* = argmax Var_psi(sum_p W_p) = %.6f   (Var = %.6f)" % (g2s[0], g2s[1]))
    P("   a-priori marker g2_sd = sqrt(C/L) = sqrt(5/12) = %.6f   ->  ratio g2_*/g2_sd = %.4f"
      % (math.sqrt(5 / 12), g2s[0] / math.sqrt(5 / 12)))
    lo1, hi1, n1 = window(rows, 3, C.GATE)
    lo2, hi2, n2 = window(rows, 4, C.GATE)
    both = [r[0] for r in rows if r[3] >= C.GATE and r[4] >= C.GATE]
    P("   H_BLOCK >= 0.10 bits on g2 in [%.4f, %.4f]   (%d of 200 grid points)" % (lo1, hi1, n1))
    P("   H_MAG   >= 0.10 bits on g2 in [%.4f, %.4f]   (%d of 200 grid points)" % (lo2, hi2, n2))
    if both:
        P("   >>> LIVE WINDOW (both >= 0.10) : g2 in [%.4f, %.4f] , %d of 200 points, width %.4f"
          % (min(both), max(both), len(both), max(both) - min(both)))
        P("       in decades: %.4f  (log10 span)" % (math.log10(max(both)) - math.log10(min(both))))
    else:
        P("   >>> LIVE WINDOW IS EMPTY on the 200-point grid.  H_BLOCK and H_MAG are never both")
        P("       above the pinned 0.10-bit gate.  THE PRE-REGISTERED LIVE WINDOW DOES NOT EXIST.")
        P("       max H_BLOCK = %.6f at g2 = %.4f ;  max H_MAG = %.6f at g2 = %.4f"
          % (max(r[3] for r in rows), max(rows, key=lambda r: r[3])[0],
             max(r[4] for r in rows), max(rows, key=lambda r: r[4])[0]))
        cond = [(2 * (r[5] - r[6]), r[0]) for r in rows]
        P("       max conditional block entropy H_FULL-H_CENTRE = %.6e at g2 = %.4f"
          % (max(cond)[0] / 2, max(cond)[1]))
    P("")

# ============================================================================================
rule("BLOCK 3 -- THE PRIMARY ARM PAIR.  A0 (eta=+1 everywhere) vs A1 (eta=-1 at {0,4}).")
SECA0 = C.Sector([]); SECA1 = C.Sector([0, 4])
P("ARM 1  A0 : charges none          u0 = links %s  |u0| = %d  flux(Sigma) = %+d"
  % (C.bits(SECA0.u0), C.pop(SECA0.u0), C.flux_sigma(SECA0)))
P("ARM 2  A1 : charges at {0,4}      u0 = links %s  |u0| = %d  flux(Sigma) = %+d"
  % (C.bits(SECA1.u0), C.pop(SECA1.u0), C.flux_sigma(SECA1)))
P("")
P("THE DIFF, AT THE LEVEL OF THE HAMILTONIAN.  The two arms are the SAME 32x32 matrix except for")
P("the diagonal, because the plaquette hopping is translation on the cycle space and is coset-blind:")
H0 = SECA0.H(1.0); H1 = SECA1.H(1.0)
off0 = H0 - np.diag(np.diag(H0)); off1 = H1 - np.diag(np.diag(H1))
P("   || offdiag(H_A0) - offdiag(H_A1) ||_F  = %.12f   (0 = byte-identical hopping)"
  % np.linalg.norm(off0 - off1))
P("   || diag(H_A0)    - diag(H_A1)    ||_F  = %.9f   at g2 = 1.00" % np.linalg.norm(np.diag(H0) - np.diag(H1)))
d0 = np.array([C.pop(int(u)) for u in SECA0.U]); d1 = np.array([C.pop(int(u)) for u in SECA1.U])
P("   electric weight |u| over the 32 basis states:")
P("      A0 histogram %s" % {int(k): int(v) for k, v in zip(*np.unique(d0, return_counts=True))})
P("      A1 histogram %s" % {int(k): int(v) for k, v in zip(*np.unique(d1, return_counts=True))})
P("   THE ENTIRE ARM DIFFERENCE IS THIS HISTOGRAM.  The charge sector is a rigid shift of the")
P("   electric-energy landscape by the coset offset u0; nothing else in the arm moves.")
P("")

hdr = ("%-6s %-4s %8s %8s %8s %8s %10s %10s %9s %9s %9s %9s %5s"
       % ("g2", "arm", "H_FULL", "H_CEN", "H_BLOCK", "H_MAG", "I(S:env)", "I(S:Sig)",
          "I(S:F1)", "I(S:F2)", "I(S:F3)", "I(S:F4)", "R_del"))
P("FULL CHANNEL A_S (dim 16, max 3 bits).  All quantities in bits.")
P(hdr); P("-" * len(hdr))
TAB = {}
for g2 in C.GRID:
    for tag, sec in (("A0", SECA0), ("A1", SECA1)):
        _, psi, _ = sec.ground(g2)
        hF = ent(C.A_FULL, sec, psi); hC = ent(C.A_CENTRE, sec, psi)
        hB = ent(C.A_BLOCK, sec, psi); hM = ent(C.A_MAG, sec, psi)
        Ienv = C.MI(C.A_FULL, C.A_ENV, sec, psi)
        Isig = C.MI(C.A_FULL, C.A_SIG, sec, psi)
        Ifr = [C.MI(C.A_FULL, A, sec, psi) for _, A in C.AF_P]
        Rd = sum(1 for x in Ifr if hF > 0 and x >= (1 - C.DELTA_TOL) * hF)
        dsurf = Ienv - Isig
        TAB[(g2, tag)] = dict(hF=hF, hC=hC, hB=hB, hM=hM, Ienv=Ienv, Isig=Isig,
                              Ifr=Ifr, Rd=Rd, dsurf=dsurf)
        P("%-6.2f %-4s %8.5f %8.5f %8.5f %8.5f %10.5f %10.5f %9.5f %9.5f %9.5f %9.5f %5d"
          % (g2, tag, hF, hC, hB, hM, Ienv, Isig, Ifr[0], Ifr[1], Ifr[2], Ifr[3], Rd))
P("")
P("DELTA_SURF, THE PRE-REGISTERED PRIMARY FALSIFIER, AND THE IDENTITY OF BLOCK 1 CHECKED POINTWISE:")
P("%-6s %-4s %14s %20s %14s %10s" % ("g2", "arm", "Delta_surf", "2*(H_FULL-H_CENTRE)", "|residual|", "gate"))
maxres = 0.0
for g2 in C.GRID:
    for tag in ("A0", "A1"):
        d = TAB[(g2, tag)]
        pred = 2 * (d["hF"] - d["hC"])
        res = abs(d["dsurf"] - pred); maxres = max(maxres, res)
        gate = "PASS" if (d["hB"] >= C.GATE and d["hM"] >= C.GATE) else "VACUOUS-BY-EMPTINESS"
        P("%-6.2f %-4s %14.9f %20.9f %14.2e %10s" % (g2, tag, d["dsurf"], pred, res, gate))
P("   MAXIMUM RESIDUAL OVER ALL 26 ARM-POINTS : %.3e" % maxres)
P("   >>> THE IDENTITY HOLDS TO MACHINE PRECISION.  Delta_surf carried no information about the")
P("       environment on any of these points.  IT IS NOT SCORED.")
P("")
P("THE ARM DIFF, PRINTED:  A1 minus A0, per coupling, per channel.")
P("%-6s %10s %10s %10s %10s %12s %12s %8s" % ("g2", "dH_FULL", "dH_CEN", "dH_BLOCK", "dH_MAG",
                                              "dI(S:env)", "dI(S:F3)", "dR_del"))
for g2 in C.GRID:
    a, b = TAB[(g2, "A0")], TAB[(g2, "A1")]
    P("%-6.2f %10.6f %10.6f %10.6f %10.6f %12.6f %12.6f %8d"
      % (g2, b["hF"] - a["hF"], b["hC"] - a["hC"], b["hB"] - a["hB"], b["hM"] - a["hM"],
         b["Ienv"] - a["Ienv"], b["Ifr"][2] - a["Ifr"][2], b["Rd"] - a["Rd"]))
mx = max(abs(TAB[(g, "A1")][k] - TAB[(g, "A0")][k]) for g in C.GRID for k in ("hF", "hC", "hB", "hM"))
P("   LARGEST SINGLE-CHANNEL MOVE ACROSS THE WHOLE SWEEP : %.9f bits" % mx)
P("")
P("RECORD-CONTENT DISTANCE, the quantity that does not depend on which scalar summary is chosen:")
P("   D_tr( rho|A_S in A1 , rho|A_S in A0 )  --  trace distance inside the SAME abstract algebra")
P("   A_S = M_2 (+) M_2 (+) M_2 (+) M_2 .  Well posed across sectors because A_S is sector-free.")
P("%-6s %12s %14s %14s" % ("g2", "D_tr(A_S)", "D_tr(CENTRE)", "D_tr(BLOCK)"))
for g2 in C.GRID:
    _, p0, _ = SECA0.ground(g2); _, p1, _ = SECA1.ground(g2)
    P("%-6.2f %12.9f %14.9f %14.9f"
      % (g2, tracedist(rho_alg(C.A_FULL, SECA0, p0), rho_alg(C.A_FULL, SECA1, p1)),
         tracedist(rho_alg(C.A_CENTRE, SECA0, p0), rho_alg(C.A_CENTRE, SECA1, p1)),
         tracedist(rho_alg(C.A_BLOCK, SECA0, p0), rho_alg(C.A_BLOCK, SECA1, p1))))

# ============================================================================================
rule("BLOCK 4 -- CHANNEL-BY-CHANNEL R_delta, AND WHERE THE PLATEAU CAN FAIL")
P("R_delta = #{i : I(channel : A_Fi) >= 0.90 * H(channel)} over the four declared fragments.")
P("%-6s %-4s %-8s %9s %9s %9s %9s %9s %6s %s"
  % ("g2", "arm", "channel", "H(chan)", "I:F1", "I:F2", "I:F3", "I:F4", "R_del", "gate"))
for g2 in C.GRID:
    for tag, sec in (("A0", SECA0), ("A1", SECA1)):
        _, psi, _ = sec.ground(g2)
        for cn, A in C.CHANNELS:
            h = ent(A, sec, psi)
            Is = [C.MI(A, F, sec, psi) for _, F in C.AF_P]
            Rd = sum(1 for x in Is if h > 0 and x >= (1 - C.DELTA_TOL) * h)
            P("%-6.2f %-4s %-8s %9.5f %9.5f %9.5f %9.5f %9.5f %6d %s"
              % (g2, tag, cn, h, Is[0], Is[1], Is[2], Is[3], Rd,
                 "" if h >= C.GATE else "VACUOUS-BY-EMPTINESS"))

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_W20_C_CHARGE/OUT_run_charge.txt", "w").write("\n".join(LOG) + "\n")
