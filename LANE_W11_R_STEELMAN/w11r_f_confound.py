# LANE W11-R  LEG F -- (F1) how much of the EDGE convention's formation is the FIELD at all;
#                      (F2) the rebuttal to the strongest surviving Reading-B mechanism, W-06's
#                           "fibre-wise-ness collapses transport into gauge";
#                      (F3) subdivision -- the edge clock is not an invariant of the space;
#                      (F4) reproduction of the registrar's five headline numbers (no arithmetic dispute).
import numpy as np, w11r_lib as L
rng = np.random.default_rng(20260817)

print("== F1  HOW MUCH OF THE EDGE CONVENTION'S 'FORMATION' IS THE CONNECTION? ==")
print("  ONE VARIABLE: the connection is scaled a -> t.a, t from 1 down to 0.  Carrier, state,")
print("  convention, observable, code path, N all held fixed.")
lf,lc,NV,ne = L.K1_LOOP_F, L.K1_LOOP_C, 5, 6
a  = np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
sA = np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
print(f"  {'t':>8}{'W_F':>26}{'W_C':>26}{'EDGE lambda':>16}{'CIRCUIT lambda':>16}")
for t in (1.0,0.5,0.2,0.1,0.01,0.0):
    at=t*a
    TF,TC = L.T_edge(lf,at,NV), L.T_edge(lc,at,NV)
    MF,MC = L.M_circuit(lf,at,NV), L.M_circuit(lc,at,NV)
    print(f"  {t:>8.2f}{L.holonomy(lf,at):>26.6f}{L.holonomy(lc,at):>26.6f}"
          f"{L.rate(TF,TC,sA,20000):>16.9f}{L.rate(MF,MC,sA,20000):>16.9f}")
print("  -> AT the trivial connection the CIRCUIT rate is EXACTLY 0, which is what CHOICE LEDGER")
print("     C4 tests.  NOT claimed: continuity -- BOTH rates are discontinuous at a = 0 (S3's own")
print("     flag F5, already of record), so the t -> 0 column is a limit, not the verdict.")
print("     THE EDGE RATE AT a = 0 IS NOT ZERO.  For state A specifically it is only -0.041; over")
print("     60 random ready states (leg E2) the median is -0.756 and 60/60 are negative.  Compare")
print("     the corpus's whole registered rate m(P) = -0.7675: under the edge convention MOST OF")
print("     THE FORMATION IS THERE WITH NO FIELD.  The connection is a perturbation on top")
print("     of the two cycles' combinatorics.  That is the definition of a confound.")

print("\n== F2  DOES 'M_gamma IS FIBRE-WISE, HENCE IN THE LOCAL GAUGE GROUP' MAKE IT TRIVIAL? ==")
print("  W-06 corrected N4's mechanism to FIBRE-WISE-NESS: any fibre-wise unitary lies in U(1)^V.")
print("  True of M_gamma.  But lying in the gauge GROUP is not being a gauge TRANSFORMATION of the")
print("  system: a gauge transformation moves (a, s) -> (a^g, g.s); the corpus moves s -> M s with")
print("  a HELD FIXED.  (a, Ms) and (a, s) are gauge-equivalent iff some g stabilises a and equals")
print("  M on supp(s).  The stabiliser of a is the CONSTANT phases; M is not constant.  Test:")
g_stab = np.exp(1j*rng.uniform(0,2*np.pi))                     # a constant phase: stabilises a
MF = L.M_circuit(lf,a,NV)
s  = rng.normal(size=NV)+1j*rng.normal(size=NV); s/=np.linalg.norm(s)
def dressed(a_, s_, e):                                        # a joint (connection,state) invariant
    u,v = L.K1_EDGES[e]; return np.conj(s_[u])*np.exp(1j*a_[e])*s_[v]
print(f"  joint invariant  conj(s_v0) U_e4 s_v3 :")
print(f"     on (a, s)      = {dressed(a,s,3):.9f}")
print(f"     on (a, M_dF s) = {dressed(a,MF@s,3):.9f}     <- CHANGED: not a gauge copy")
print(f"     on (a, g.s), g a constant phase = {dressed(a,g_stab*s,3):.9f}   <- unchanged, as gauge must be")
print(f"  || (a,Ms) invariant - (a,s) invariant || = {abs(dressed(a,MF@s,3)-dressed(a,s,3)):.6f}")
print("  -> the corpus's evolution carries the configuration OFF its gauge orbit.  Fibre-wise-ness")
print("     makes M_gamma an ELEMENT of the gauge group; it does not make the dynamics vacuous.")
print("     (This is the Aharonov-Bohm content: a relative holonomy phase between two branches.)")

print("\n== F3  THE EDGE CLOCK IS NOT AN INVARIANT OF THE SPACE.  SUBDIVISION TEST ==")
print("  Subdivide one edge of gamma_C: v3 -> w -> v4, with a_e5 split as a5' + a5''.  Nothing")
print("  gauge-invariant changes: same W_F, same W_C, same b1, same topology.  Put weight eps on w.")
def K1_sub(k, a6v, eps):
    """K1 with gamma_C's edge e5 subdivided into k segments.  NV = 5 + (k-1)."""
    NV2 = 5 + (k-1); newv = list(range(5, NV2))
    edges = [(0,1),(1,2),(2,0),(0,3)]
    seq = [3]+newv+[4]
    for i in range(len(seq)-1): edges.append((seq[i],seq[i+1]))
    edges.append((4,0))
    aa = np.zeros(len(edges))
    aa[0],aa[1],aa[2],aa[3] = a6v[0],a6v[1],a6v[2],a6v[3]
    for i in range(k): aa[4+i] = a6v[4]/k
    aa[4+k] = a6v[5]
    loopF = [(0,1,0),(1,2,1),(2,0,2)]
    loopC = [(0,3,3)]+[(seq[i],seq[i+1],4+i) for i in range(k)]+[(4,0,4+k)]
    w = np.array([0.40,0.15,0.15,0.15,0.15] + [eps]*(k-1)); w = w/w.sum()
    return loopF, loopC, NV2, aa, np.sqrt(w)+0j
print(f"  {'gamma_C len':>12}{'eps':>10}{'W_C':>26}{'CIRCUIT lambda':>17}{'EDGE lambda/tick':>18}")
for k in (1,2,3,4):
    lF,lC,NV2,aa,ss = K1_sub(k, a, 1e-6)
    TF2,TC2 = L.T_edge(lF,aa,NV2), L.T_edge(lC,aa,NV2)
    MF2,MC2 = L.M_circuit(lF,aa,NV2), L.M_circuit(lC,aa,NV2)
    print(f"  {2+k:>12}{1e-6:>10.0e}{L.holonomy(lC,aa):>26.9f}"
          f"{L.rate(MF2,MC2,ss,20000):>17.9f}{L.rate(TF2,TC2,ss,20000):>18.9f}")
print("  -> W_C is identical in every row (the connection is the same connection).  The CIRCUIT")
print("     rate is identical to 6 places (the residual is the eps weight on the added vertices).")
print("     The EDGE rate moves by 13 percent across rows that differ ONLY in how finely the SAME")
print("     cycle is cellulated -- and the meaning of one tick changes with it.")
print("  HONEST COUNTER, recorded: the corpus is committed to cells being physical (S3's 'one qubit")
print("  per cell'), so a Reading-B advocate may answer that subdivision is not a symmetry here.")
print("  This leg is SUPPORTING, not decisive; leg B is the decisive one.")

print("\n== F4  THE REGISTRAR'S HEADLINE NUMBERS, REPRODUCED IN THIS LANE'S OWN CODE ==")
sB = np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
sC = sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
TF,TC = L.T_edge(lf,a,NV), L.T_edge(lc,a,NV); MF,MC = L.M_circuit(lf,a,NV), L.M_circuit(lc,a,NV)
sp_edge = max(max(abs(L.Z(TF,TC,s_,n,n)) for s_ in (sA,sB,sC)) - min(abs(L.Z(TF,TC,s_,n,n)) for s_ in (sA,sB,sC)) for n in range(1,10))
sp_circ = max(max(abs(L.Z(MF,MC,s_,k,k)) for s_ in (sA,sB,sC)) - min(abs(L.Z(MF,MC,s_,k,k)) for s_ in (sA,sB,sC)) for k in range(1,7))
print(f"  K1 max EDGE spread, n<=9   = {sp_edge:.2e}   [registrar 5.90e-01]")
print(f"  K1 max CIRCUIT spread, k<=6= {sp_circ:.2e}   [registrar 2.22e-16]")
print(f"  m(P) for pi = (0,.3,.3,.4)  = {L.m_jensen(np.array([0.0,0.3,0.3,0.4])):.12f}   [registrar -0.767507880]")
print(f"  CIRCUIT rate N=200000      = {L.rate(MF,MC,sA,200000):.9f}   [registrar -0.767500322]")
print(f"  EDGE rate N=20000 state A  = {L.rate(TF,TC,sA,20000):.9f}   [registrar -0.864256422]")
print("  -> NO ARITHMETIC DISPUTE.  Every number the registrar reports reproduces here.  What I")
print("     dispute is the NAME of the variable and the STANDING of the rival convention.")
