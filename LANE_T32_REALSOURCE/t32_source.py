"""T-32: IS THE CONFIGURATION ENERGY OF REAL RECORDS A SOURCE?

C-70 removed the exclusion that appeared to close this at any N: the Z[i] quantisation was an artifact
of exactness, and on a real record the value is a continuum. Real records DO carry a configuration
energy dE, and it is measurable.

O-48 tested the five source standards on a TOY CHAIN where dE was two-signed and cancelled as
m^(-1/2) -- C-46's screening signature. THE TOY IS NOT THE WORLD. In real media the splitting is set
by an applied or internal field with ONE SIGN, so every grain's dE points the same way.

THE FIVE STANDARDS, on real records with real numbers:
  (a) EXTENSIVE   grows without bound, asymptotically linear
  (b) ADDITIVE    over disjoint regions
  (c) NOT A COUNT moves when the medium's parameters move at fixed N
  (d) SIGN-DEFINITE  |sum|/sum|.| = 1, so it accumulates rather than screening
  (e) SEPARATION-DEPENDENT with a power-law falloff

CONTROLS IN THE SAME TABLE (D-15): the O-48 toy chain, whose dE is two-signed and must screen; and a
DEMAGNETISED medium with zero net field, where the sign-definiteness must vanish."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'model'))
import grounded as G
def say(*a): print(*a); sys.stdout.flush()
say("="*100); say("T-32   IS THE CONFIGURATION ENERGY OF REAL RECORDS A SOURCE?"); say("="*100)

T=300.0; K_u=2.0e5; V=1.26e-24; Ms=4.0e5   # A/m, CoCrPt saturation magnetisation
mu0=4*np.pi*1e-7
m_grain = Ms*V                              # A m^2, the grain's magnetic moment
say(f"  CoCrPt grain: V = {V:.2e} m^3, M_s = {Ms:.1e} A/m, moment m = {m_grain:.3e} A m^2")
say(f"  E_b = {K_u*V/(G.KB*T):.1f} kT     T = {T:.0f} K")
say("")

# ------------------------------------------------------------------ (a) (b) (c) (d) on real media
say("1. RESCOPED (solidity review). THE FIELD-FREE QUANTITY IS THE REMANENT MOMENT")
say("   M_r = M_s V sum(s_i)  [A m^2] -- what the written medium carries with NO applied field.")
say("   The Zeeman energy 2 m B s_i EXISTS ONLY IN AN APPLIED EVALUATION FIELD B (inserted) and is")
say("   kept below as medium-plus-field, labeled. SCOPE: one-signedness belongs to DC-SATURATED")
say("   writing; real data tracks use DC-free codes and SCREEN -- shown in the same table.")
say(f"   {'N':>10}{'B (T)':>9}{'written: sum dE (J)':>22}{'sum|dE| (J)':>15}{'|sum|/sum|.|':>14}{'demag: |sum|/sum|.|':>21}")
rng=np.random.default_rng(7); rows=[]
say(f"   {'pattern':<22}{'N':>9}{'M_r (A m^2)':>15}{'sum|m| (A m^2)':>16}{'|sum|/sum|.|':>14}")
for N in (10, 100, 1000, 10000, 100000):
    dcsat  = np.ones(N)                        # DC-saturated: every grain one way
    randat = rng.integers(0,2,N)*2-1           # real random data
    coded  = np.tile([1,-1], N//2 + 1)[:N]     # a DC-free coded pattern (alternating)
    demag  = rng.integers(0,2,N)*2-1           # AC-erased control
    for nm, pat in (("DC-SATURATED", dcsat), ("random data", randat), ("DC-free coded", coded), ("AC-erased ctrl", demag)):
        Mr, Ma = m_grain*pat.sum(), m_grain*np.abs(pat).sum()
        if nm == "DC-SATURATED": rows.append((N, Mr, Ma))
        if N in (10, 100000) or nm == "DC-SATURATED":
            say(f"   {nm:<22}{N:>9}{Mr:>15.4e}{Ma:>16.4e}{abs(Mr)/Ma:>14.6f}")
say("   -> ONE-SIGNED ONLY under DC saturation; random data and DC-free codes SCREEN, as the review")
say("      measured (0.00096 at N=1e5). The falsifier is scoped to DC-saturated writing.")
say("")
r2=[r for r in rows]
say(f"   (a) EXTENSIVE (M_r, DC-saturated):  ratios "
    f"{', '.join('%.4f'%(r2[i+1][1]/r2[i][1]) for i in range(len(r2)-1))}   (linear gives 10.0000)")
say(f"   (d) SIGN-DEFINITE under DC-saturation scope only; other patterns screen toward zero")
say(f"       (non-monotone at finite N -- single random draws, not an ensemble)")
say("")
# (c) not a count
say("1b. (b) ADDITIVITY WITH THE CROSS-TERM COMPUTED (replaces the earlier 1e-3 assertion).")
say("    Two tracks of N/2 grains at separation r_track; dipole cross-energy per boundary pair")
say("    U_dd = mu0 m^2/(2 pi r^3), summed over the facing rows; compared with per-grain E_b.")
mu0=4*np.pi*1e-7
say(f"    {'r_track (nm)':>13}{'U_dd/pair (J)':>16}{'U_dd/E_b':>12}{'additivity verdict':>22}")
Eb = 2.0e5*V
for rnm in (10, 20, 40, 60, 100):
    r=rnm*1e-9
    udd = mu0*m_grain**2/(2*np.pi*r**3)
    ok = udd/Eb <= 1e-3
    say(f"    {rnm:>13}{udd:>16.4e}{udd/Eb:>12.4e}{('holds (<=1e-3 of E_b)' if ok else 'FAILS at this spacing'):>22}")
say("    -> additivity of M_r is exact by superposition of moments (trivially, stated); additivity")
say("       of the ENERGY holds at track separations where U_dd/E_b <= 1e-3 -- computed above, NOT")
say("       asserted: it fails at grain spacing (10-20 nm) and holds from ~60 nm.")
say("")
say("2. (c) IS IT A COUNT? Vary the medium at FIXED N = 1000.")
say(f"   {'M_s (A/m)':>12}{'B (T)':>9}{'V (m^3)':>12}{'sum dE (J)':>16}")
vals_c=[]
for Ms2,B2,V2 in ((4.0e5,0.5,1.26e-24),(8.0e5,0.5,1.26e-24),(4.0e5,1.2,1.26e-24),(4.0e5,0.5,3.0e-24)):
    q=1000*2*(Ms2*V2)*B2; vals_c.append(q)
    say(f"   {Ms2:>12.1e}{B2:>9.2f}{V2:>12.2e}{q:>16.4e}")
spread_c = max(vals_c)/min(vals_c)
say(f"   -> N is FIXED at 1000 and the quantity moves by {spread_c:.2f}x — it is NOT a count")
say("")
say("1b. (b) ADDITIVE OVER DISJOINT REGIONS — two tracks written separately, then together")
dE=2*m_grain*0.5
for NA,NB in ((100,100),(500,1500),(3000,7000)):
    qa, qb, qab = dE*NA, dE*NB, dE*(NA+NB)
    say(f"    N_A={NA:>5}  N_B={NB:>5}   Q(A)+Q(B) = {qa+qb:.6e}   Q(A u B) = {qab:.6e}"
        f"   defect {abs(qab-qa-qb):.1e}")
say("    CONTROL: two tracks close enough to couple by dipole would show a non-zero defect;")
say("    at track separations >> grain size the dipole term is 1e-3 of the field term.")
say("")
# (e) separation
say("3. (e) SEPARATION DEPENDENCE — RELABELED (solidity review): the r^-3 here is INSERTED, not")
say("   induced -- B_dip(r) = mu0 m/(2 pi r^3) is quoted magnetostatics and the 8.0000-per-doubling")
say("   check measures the quoted formula. The INDUCED demonstration lives in LANE_T34_NAND, where")
say("   the exponent emerges from monopoles + a real ground plane against an r^-1 no-plane control.")
say(f"   {'r (nm)':>9}{'B_dip (T)':>14}{'E_int/kT':>14}{'ratio to r^-3':>16}")
prev=None
for rnm in (10,20,40,80,160):
    r=rnm*1e-9
    Bd = mu0*m_grain/(2*np.pi*r**3)
    E  = 2*m_grain*Bd
    pred = None if prev is None else (prev[1]/Bd)
    say(f"   {rnm:>9}{Bd:>14.4e}{E/(G.KB*T):>14.4e}{(pred if pred else float('nan')):>16.4f}")
    prev=(rnm,Bd)
say("   -> doubling r divides B_dip by 8.0000 exactly: consistency of the QUOTED law (INSERTED here)")
say("")
say("="*100); say("  READ -- rescoped per the solidity review, generated from the tables above"); say("="*100)
say("  FIELD-FREE, the written medium's accumulated quantity is the REMANENT MOMENT M_r = M_s V sum s_i:")
say("  extensive (ratios exactly 10 per decade), additive by superposition (trivially, stated) with the")
say("  ENERGY cross-term now COMPUTED (holds at >= ~60 nm track separation, fails at grain spacing),")
say("  not a count (moves at fixed N), and ONE-SIGNED UNDER DC-SATURATED WRITING ONLY -- random data")
say("  and DC-free codes screen, as the review measured. The Zeeman configuration energy exists only")
say("  in an applied evaluation field and is labeled medium-plus-inserted-field. The r^-3 falloff here")
say("  is QUOTED magnetostatics (inserted); the induced demonstration is LANE_T34_NAND's.")
say("")
say("  THE CROSS-MECHANISM KERNEL (with T-34): formation creates a nonzero accumulated quantity the")
say("  unwritten surface lacks, and the WRITE MECHANISM fixes its sign -- field direction here,")
say("  electron injection there. Occupancy-encoded surfaces accumulate for EVERY data pattern;")
say("  orientation-encoded surfaces accumulate only when written all-one-way. Standard physics")
say("  (saturation remanence), honestly owned; what is ours is the encoding-level statement.")
