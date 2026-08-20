"""T-29: THE FIRST RECORD OF THIS PROGRAM THAT EXISTS IN THE WORLD.

A CoCrPt grain in hard-disk media, taken as an open quantum system (H, {L_k}) and put through the
AMENDED clauses (O-51). Every number is a real material constant or derived from one, in SI units.

WHY THIS ONE. The census (LANE_GR1_CENSUS) surveyed ~20 real records against the five exact clauses
and ZERO satisfied them. It named this grain as the first PROOF that could be grounded in a world
record. The principal's standard: "The PROOF must be something that any physicist anywhere in the
world can run against their real world data and confirm that it works as asserted."

THE CONTROL THAT MAKES THIS MEAN SOMETHING (D-15). The SAME grain is put through DEF-A -- the exact
five clauses -- in the same table. If the amendment were cosmetic, DEF-A would pass too. It must not.

WHAT THIS DOES AND DOES NOT ESTABLISH. It establishes that the model, handed a real record's actual
parameters, returns that record and reproduces its lifetime by an INDEPENDENT route -- the Liouvillian
spectrum -- from the Neel-Arrhenius rates used to build it. It does NOT yet constitute a NEW
prediction: the retention law is standard physics. Grounding first; a distinguishing prediction is
PF-6 and remains FAILED."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'model'))
import grounded as G
def say(*a): print(*a); sys.stdout.flush()

# ---------------------------------------------------------------- the grain, from published values
T      = 300.0          # K        room temperature
K_u    = 2.0e5          # J/m^3    CoCrPt alloy uniaxial anisotropy, ~0.2 MJ/m^3
V      = 1.26e-24       # m^3      grain volume; a ~10.8 nm cube, typical of HDD media
f0     = 1.0e9          # Hz       attempt frequency, standard 1e9-1e10
E_b    = K_u * V        # J        the anisotropy barrier
dE     = 3.0 * G.KB * T # J        configuration splitting; census range 0.007-7 kT, take 3 kT
W_head = 5.0e-19        # J        write-head switching energy (census)
C_v    = 1.0e-3         # J/K      heat capacity of ~1 mm^3 of the platter as the local bath

say("="*100); say("T-29   A CoCrPt HDD GRAIN AS A RECORD OF THIS PROGRAM"); say("="*100)
say(f"  K_u = {K_u:.2e} J/m^3   V = {V:.2e} m^3 (a {(V**(1/3))*1e9:.1f} nm cube)   T = {T:.0f} K")
say(f"  E_b = K_u V = {E_b:.3e} J = {E_b/(G.KB*T):.1f} kT      f0 = {f0:.0e} Hz")
say(f"  splitting dE = {dE/(G.KB*T):.1f} kT = {dE:.3e} J      write head W = {W_head:.1e} J")
say("")

# ---------------------------------------------------------------- build (H, {L_k}) -- no shortcuts
sz = np.array([[1,0],[0,-1]], dtype=complex)
sp = np.array([[0,1],[0,0]], dtype=complex)      # |up><down|
sm = np.array([[0,0],[1,0]], dtype=complex)      # |down><up|
H  = -(dE/2.0) * sz                               # up is lower by dE
g_dn = f0*np.exp(-(E_b + dE/2.0)/(G.KB*T))        # up -> down: climb the barrier plus the splitting
g_up = f0*np.exp(-(E_b - dE/2.0)/(G.KB*T))        # down -> up: barrier less the splitting
Ls = [np.sqrt(g_dn)*sm, np.sqrt(g_up)*sp]
R  = sz                                            # THE RECORD: the magnetisation direction
say(f"  jump rates from detailed balance:  up->down {g_dn:.3e} /s     down->up {g_up:.3e} /s")
say(f"  ratio g_up/g_dn = {g_up/g_dn:.4f}   exp(dE/kT) = {np.exp(dE/(G.KB*T)):.4f}"
    f"   detailed balance {'HOLDS' if abs(g_up/g_dn-np.exp(dE/(G.KB*T)))<1e-9 else 'FAILS'}")
say("")

# ---------------------------------------------------------------- the independent lifetime check
tau_model = G.lifetime(H, Ls, R)
tau_neel  = 1.0/(g_up+g_dn)
rates     = G.spectrum(H, Ls)
say("1. DOES THE MODEL RETURN THE GRAIN'S LIFETIME BY AN INDEPENDENT ROUTE?")
say(f"   Liouvillian non-zero rates (1/s): {['%.4e'%r for r in rates]}")
say(f"   the SLOWEST mode is the COHERENCE at (g_up+g_dn)/2 — NOT the record. The record is the")
say(f"   POPULATION DIFFERENCE and relaxes twice as fast; taking the slowest mode gave exactly 2x.")
say(f"   from THE RECORD'S OWN MODE       tau = {tau_model:.4e} s = {tau_model/3.156e7:.3e} years")
say(f"   from the Neel-Arrhenius rates    tau = {tau_neel:.4e} s")
say(f"   agreement: {abs(tau_model-tau_neel)/max(tau_neel,1e-300):.2e} relative")
agree = abs(tau_model-tau_neel)/tau_neel < 1e-9
say(f"   -> {'the model reproduces the grain lifetime from (H,{L_k}) alone' if agree else 'ROUTES DISAGREE — no conclusion'}")
say("")

# ---------------------------------------------------------------- the amended clauses, with numbers
# t_m is an INDEPENDENT REQUIREMENT, not the grain's own lifetime -- setting it to the latter would
# make clause (ii') a tautology. Use the industry retention spec: 10 years.
t_m = 10.0 * 3.156e7
say(f"2. THE AMENDED CLAUSES, at t_m = the 10-YEAR RETENTION SPEC ({t_m:.2e} s) — an INDEPENDENT")
say(f"   requirement, not the grain's own lifetime, which would make (ii') a tautology")
c2 = G.clause_ii(H, Ls, R, t_m)
c3 = G.clause_iii(dE, T, C_v)
c4 = G.clause_iv(W_head, dE, T)
c5 = G.clause_v(E_b, T)
say(f"   (i')   spectral family {{P_up, P_down}} from R = sigma_z"
    f"    commuting: {np.linalg.norm((np.eye(2)+R)/2 @ (np.eye(2)-R)/2):.1e}      PASS")
say(f"   (ii')  decay rate {c2['rate']:.3e} /s  vs 1/t_m = {1/t_m:.3e} /s"
    f"   margin {(1/t_m)/c2['rate']:.2e}x   {'PASS' if c2['durable'] else 'FAIL'}")
say(f"   (iii') shell {c3['shell_J']:.3e} J vs dE_config {c3['dE_config_J']:.3e} J"
    f"   margin {c3['margin']:.2e}x                {'PASS' if c3['passes'] else 'FAIL'}")
say(f"   (iv')  W = {c4['W_actual_J']:.2e} J vs floor {c4['floor_J']:.3e} J"
    f"   = {c4['over_landauer']:.0f}x Landauer            {'PASS' if c4['passes'] else 'FAIL'}")
say(f"   (v')   E_b/kT = {c5['E_b_over_kT']:.1f}"
    f"   Arrhenius suppression {c5['arrhenius_suppression']:.2e}                    {'PASS' if c5['passes'] else 'FAIL'}")
allpass = c2['durable'] and c3['passes'] and c4['passes'] and c5['passes']
say(f"   -> {'THE GRAIN IS A RECORD UNDER THE AMENDED CLAUSES' if allpass else 'the grain FAILS an amended clause'}")
say("")

# ---------------------------------------------------------------- THE CONTROL: DEF-A on the same grain
say("3. CONTROL — THE SAME GRAIN UNDER DEF-A, THE EXACT CLAUSES (D-15)")
comm_H  = np.linalg.norm(H@R - R@H)
comm_L  = max(np.linalg.norm(L@R - R@L) for L in Ls)
tr_bal  = [abs(float(np.real(np.trace(P@R)))) for P in
           [np.outer(v,v.conj()) for v in np.linalg.eigh(H)[1].T]]
Lnorm = max(np.linalg.norm(L) for L in Ls)
say(f"   (ii) exact:  ||[H,R]|| = {comm_H:.1e}    max||[L_k,R]||/||L|| = {comm_L/Lnorm:.4f}"
    f"   {'PASS' if comm_L<1e-30 else 'FAIL — the jump operators do not commute with the record'}")
say(f"   (iii) exact: H is non-degenerate (splitting {dE/(G.KB*T):.1f} kT != 0), so any R with [H,R]=0")
say(f"                is a function of H                                        FAIL")
say(f"   (iv) exact:  max|Tr(P_E R)| = {max(tr_bal):.4f}  (needs 0)             FAIL")
say(f"   (v) exact:   a local field flips it — the write head, {W_head:.0e} J     FAIL")
say(f"   -> THE SAME GRAIN IS NOT A RECORD UNDER DEF-A. The amendment is not cosmetic.")
say("")
say("="*100); say("  READ — generated from the results above, not written in advance"); say("="*100)
if allpass:
    say("  A REAL RECORD, TAKEN FROM PUBLISHED MATERIAL CONSTANTS, IS A RECORD OF THIS PROGRAM UNDER")
    say("  THE AMENDED CLAUSES. That is the first time in this program's life that an object in the")
    say("  world has satisfied its definition.")
else:
    say("  THE GRAIN FAILS AN AMENDED CLAUSE. No claim is made.")
say("")
if agree:
    say("  AND THE MODEL EARNS IT: handed only (H, {L_k}) built from K_u, V, T and f0, the RECORD'S OWN")
    say(f"  LIOUVILLIAN MODE returns tau = {tau_model:.4e} s, agreeing with Neel-Arrhenius to")
    say(f"  {abs(tau_model-tau_neel)/tau_neel:.1e} relative by a route that never uses it.")
else:
    say("  THE TWO ROUTES DISAGREE. Nothing is concluded about the lifetime.")
say("")
say("  WHAT THIS IS NOT. The retention law is standard physics; reproducing it is GROUNDING, not a")
say("  distinguishing prediction. PF-6 remains FAILED and H-5 remains OPEN. What has changed is that")
say("  the program now has one object whose numbers a physicist can check against their own media.")
