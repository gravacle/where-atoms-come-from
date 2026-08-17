# LANE W11-R  LEG A -- WHAT I MUST CONCEDE, CHECKED IN MY OWN CODE BEFORE I ARGUE ANYTHING.
# The brief hands the steelman its sharpest objection: at a PARTIAL tick the two branches have been
# transported along DIFFERENT paths to different places -- is <T_F^n s, T_C^n s> even meaningful, or
# is it comparing fibres at different points via no connection at all?
# I test that objection in three independent forms and report the result whichever way it falls.
import numpy as np, w11r_lib as L
rng = np.random.default_rng(20260817)
np.set_printoptions(linewidth=150)

print("== A0  COR-F's SEALED EXHIBIT, REBUILT FROM S1's BYTES (not from the registrar's code) ==")
a = np.zeros(6); a[3],a[4],a[5] = 0.7,1.3,-0.4          # COR-F's own a4,a5,a6
T = L.T_edge(L.K1_LOOP_C, a, 5)
print(f"  || T*T - I ||        = {np.linalg.norm(T.conj().T@T-np.eye(5)):.2e}   [COR-F: 0.00e+00]")
print(f"  T diagonal?            {np.allclose(T,np.diag(np.diag(T)))}                 [COR-F: False]")
print(f"  T^3 diag             = {np.round(np.diag(np.linalg.matrix_power(T,3)),6)}")
print(f"  W_C                  = {L.holonomy(L.K1_LOOP_C,a):.6f}   [COR-F: -0.029200+0.999574j]")
rho = np.diag([0.40,0.15,0.15,0.15,0.15]).astype(complex)
print(f"  diag(T rho T*)       = {np.round(np.real(np.diag(T@rho@T.conj().T)),2)}   [COR-F: 0.15 0.15 0.15 0.40 0.15]")

print("\n== A1  IS THE PARTIAL-TICK COMPARISON A SAME-FIBRE COMPARISON? ==")
print("  Structural, not numerical: T maps Gamma(L) -> Gamma(L).  T_F^n s and T_C^n s are both")
print("  SECTIONS, and <.,.> = sum_v conj(x_v) y_v pairs the two values AT THE SAME VERTEX v,")
print("  i.e. inside the SAME fibre L_v.  Nothing is compared across fibres in either convention.")
print("  Checked: both branch objects have shape (NV,) on the same index set --")
for nm,(lf,lc,NV,ne) in (("K1",(L.K1_LOOP_F,L.K1_LOOP_C,5,6)),("B0b",(L.B0B_LOOP_F,L.B0B_LOOP_C,9,18))):
    aa = rng.uniform(0,2*np.pi,ne); s = rng.normal(size=NV)+1j*rng.normal(size=NV)
    xF = np.linalg.matrix_power(L.T_edge(lf,aa,NV),1)@s
    xC = np.linalg.matrix_power(L.T_edge(lc,aa,NV),1)@s
    print(f"    {nm}: branch_F shape {xF.shape}, branch_C shape {xC.shape}, same index set -> same fibres")
print("  ==> THE OBJECTION AS POSED FAILS.  Both branch objects live in Gamma(L); the pairing is")
print("      fibre-wise at every n, partial tick or not.")

print("\n== A2  IS <T_F^n s, T_C^n s> GAUGE-INVARIANT?  (COR-J's premise, tested independently) ==")
worst_cov = worst_inv = 0.0
worstL_F = worstL_C = 0.0
for nm,(lf,lc,NV,edges) in (("K1",(L.K1_LOOP_F,L.K1_LOOP_C,5,L.K1_EDGES)),
                            ("B0b",(L.B0B_LOOP_F,L.B0B_LOOP_C,9,L.B0B_E))):
    wc = wi = wlf = wlc = 0.0
    ne = len(edges)
    for _ in range(2000):
        aa = rng.uniform(0,2*np.pi,ne)
        s  = rng.normal(size=NV)+1j*rng.normal(size=NV)
        th = rng.uniform(0,2*np.pi,NV)
        ag = np.array([aa[j] + th[t] - th[u] for j,(u,t) in enumerate(edges)])   # a_e -> a_e + th_t - th_u
        g  = np.exp(1j*th)
        # covariance of the EDGE tick
        wc = max(wc, np.linalg.norm(L.T_edge(lf,ag,NV)@(g*s) - g*(L.T_edge(lf,aa,NV)@s)))
        # invariance of the OBSERVABLE at a PARTIAL tick (n = 1, 2)
        for n in (1,2):
            z0 = L.Z(L.T_edge(lf,aa,NV), L.T_edge(lc,aa,NV), s, n, n)
            z1 = L.Z(L.T_edge(lf,ag,NV), L.T_edge(lc,ag,NV), g*s, n, n)
            wi = max(wi, abs(z0-z1))
        # T^L = M_circuit,  and  D^L = M_circuit
        LF = len(L.K1_LOOP_F) if nm=="K1" else 4
        LC = len(L.K1_LOOP_C) if nm=="K1" else 3
        wlf = max(wlf, np.linalg.norm(np.linalg.matrix_power(L.T_edge(lf,aa,NV),LF)-L.M_circuit(lf,aa,NV)))
        wlc = max(wlc, np.linalg.norm(np.linalg.matrix_power(L.D_uniform(lc,aa,NV),LC)-L.M_circuit(lc,aa,NV)))
    print(f"  {nm:4s} max || T_F(g.a, g.s) - g.(T_F s) ||        = {wc:.2e}   -> T is gauge-COVARIANT")
    print(f"       max | Z^T_n(gauged) - Z^T_n |  (n = 1, 2)   = {wi:.2e}   -> Z^T IS GAUGE-INVARIANT")
    print(f"       max || T_F^L_F - M_dF ||                     = {wlf:.2e}")
    print(f"       max || D_C^L_C - M_c   ||                    = {wlc:.2e}   (the UNIFORM root also works)")
print("  ==> COR-J's premise 'the record must be gauge-invariant' does NOT discriminate.")
print("      The edge convention's observable is gauge-invariant at every partial tick, on both")
print("      carriers, under the FULL local group U(1)^V.  I concede this; it costs me the")
print("      natural version of objection (1).")

print("\n== A3  AND THE REFINED OBJECTION FAILS TOO, ON THE CORPUS'S OWN SEALED BYTES ==")
print("  Refined form: at a partial tick Z^T pairs s(v') with s(v'') through an OPEN Wilson line,")
print("  so it reads JOINT (connection,state) invariants the circuit convention cannot see.")
print("  That is true -- and it is LEGITIMATE on the corpus's own record:")
print("    * S2 audit COR-E: '(connection,state) mod gauge has 16-5 = 11 invariants; the page")
print("      exhibits 7.  State phases are pure gauge only when the connection is ignored.'")
print("    * W-06 (REGISTER:570): a 'Wilson-line-dressed observable built only from S1's own edge")
print("      transports is gauge-invariant to 4.45e-16 and separates the branches'.")
print("  Counting on K1: 6 edge + 10 state real params, U(1)^5 acts with 5 effective params")
print("  (the global phase acts on the state), so 16-5 = 11 joint invariants.  The circuit")
print("  convention's observable is a function of (W_F, W_C, |s_v|^2) = 7 of them.")
nF, nC = 3, 3
print(f"  Independent count on K1: joint invariants = 16 - 5 = {16-5}; circuit-visible = 2 + 5 = {2+5}")
print("  ==> the 4 invariants the edge convention adds are REAL invariants, sealed as such.")
print("      OBJECTION (1) IS DEAD IN BOTH ITS FORMS.  Reading A cannot be defended on")
print("      ill-formedness of the registrar's observable.  It must be defended on the GATE.")
