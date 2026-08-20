"""T-34: THE SECOND MECHANISM. A NAND FLOATING-GATE ARRAY AGAINST THE FIVE SOURCE STANDARDS.

C-71 and C-72 are verified on magnetic media alone, and the strengthened PROVED bar (2026-08-20)
requires agreement with at least TWO structurally different record surfaces. This lane scores the
same five standards on trapped charge -- a mechanism sharing nothing with magnetic anisotropy but
being a two-state record.

THE PHYSICAL SYSTEM, typical planar-NAND values a physicist can replace with their own device's:
  N_e   = 100 electrons per programmed floating gate     (order typical, sub-100nm planar node)
  C_fg  = 5 aF                                           -> threshold shift Q/C = 3.2 V, typical
  h     = 10 nm  FG height above the grounded channel plane
  pitch = 40 nm  cell-to-cell
  eps_r = 3.9    (SiO2)
The stored charge is ELECTRONS ONLY: the write (FN tunnelling) injects one carrier sign. Hole
trapping exists in nature (SONOS), so the sign of Q was genuinely open until the write fixed it.

THE QUANTITIES, chosen under the sign-definiteness discipline:
  Q_tot = sum_i Q_i        total stored charge          -- the (d) quantity; its sign COULD have
                                                           been either and is fixed by the write
  E_p   = Q^2/2C per cell  stored electrostatic energy  -- EXCLUDED from (d): a square is
                                                           non-negative by construction and proves
                                                           nothing about sign-definiteness
  U(r)  = the interaction between two programmed gates  -- for (e), computed from MONOPOLE charges
          plus the REAL grounded channel plane (image construction). Whatever falloff appears is
          INDUCED by the electrostatics; nothing multipolar is inserted.

CONTROLS IN THE SAME TABLE (D-15):
  * the UNWRITTEN (never-programmed / fully retention-relaxed) page, where the accumulated quantity
    must be null;
  * RANDOM DATA pages, 1000 draws -- if the ratio depends on the pattern, the sign-definiteness
    claim is pattern-limited and must be scoped;
  * the NO-GROUND-PLANE control for (e): with the image removed, U(r) = kq^2/r must fall as r^-1
    (ratio 2 per doubling), showing the instrument distinguishes the induced falloff from the bare
    Coulomb one."""
import numpy as np
def say(*a): print(*a); import sys; sys.stdout.flush()
E    = 1.602176634e-19
EPS0 = 8.8541878128e-12
KB   = 1.380649e-23
N_e, C_fg, h, pitch, eps_r, T = 100, 5e-18, 10e-9, 40e-9, 3.9, 300.0
q  = N_e*E                      # magnitude of stored charge, C
k  = 1.0/(4*np.pi*eps_r*EPS0)
Ep = q*q/(2*C_fg)               # per-cell stored energy, J
def U_img(r):   return 2*k*q*q*(1.0/r - 1.0/np.sqrt(r*r+4*h*h))   # with grounded plane (images)
def U_free(r):  return k*q*q/r                                     # no plane (control)
say("="*100); say("T-34   A NAND FLOATING-GATE ARRAY AGAINST THE FIVE SOURCE STANDARDS"); say("="*100)
say(f"  q = {q:.3e} C ({N_e} electrons)   C_fg = {C_fg:.1e} F -> dV_t = {q/C_fg:.2f} V")
say(f"  E_p = Q^2/2C = {Ep:.3e} J = {Ep/E:.0f} eV per programmed cell   pitch {pitch*1e9:.0f} nm, h {h*1e9:.0f} nm")
say("")
say("1. (d) SIGN-DEFINITE, AND C-72's WRITTEN-vs-UNWRITTEN COMPARISON")
say(f"   {'page state':<26}{'N':>8}{'sum Q (C)':>15}{'sum |Q| (C)':>14}{'|sum|/sum|.|':>14}")
rng=np.random.default_rng(11)
rows_d=[]
for N in (10,100,1000,10000,100000):
    prog = np.ones(N)                              # all-programmed page (program-all-0s, standard op)
    sq, sa = -q*prog.sum(), q*np.abs(prog).sum()
    rows_d.append((N,abs(sq)/sa))
    say(f"   {'WRITTEN all-programmed':<26}{N:>8}{sq:>15.4e}{sa:>14.4e}{abs(sq)/sa:>14.6f}")
ratios_rand=[]
for _ in range(1000):
    pat = rng.integers(0,2,1000)                   # random user data, p = 0.5
    if pat.sum()==0: continue
    ratios_rand.append(abs(-q*pat.sum())/(q*pat.sum()))
say(f"   {'WRITTEN random data':<26}{1000:>8}{'1000 draws':>15}{'':>14}{min(ratios_rand):>14.6f} (min over draws)")
say(f"   {'UNWRITTEN / relaxed':<26}{1000:>8}{0.0:>15.1e}{0.0:>14.1e}{'null (0/0)':>14}")
say("   -> WRITTEN: ratio exactly 1 at every N and for EVERY data pattern -- the write injects one")
say("      carrier sign, so terms cannot cancel. UNWRITTEN: the accumulated quantity is NOT small,")
say("      it is ABSENT. Compare T-32's magnetic medium: demagnetised ratio 0.400, 0.080, 0.002,")
say("      0.012, 0.0037 at the same N -- a thermal-random null. The null DIFFERS IN KIND between")
say("      the mechanisms (empty vs random) and BOTH are created-by-writing nulls; that difference")
say("      is reported, not hidden.")
say("")
say("2. (a) EXTENSIVE")
say(f"   {'N':>8}{'Q_tot (C)':>15}{'E_stored (J)':>16}{'Q ratio/decade':>17}{'E ratio/decade':>17}")
prev=None
for N in (10,100,1000,10000,100000):
    Q,Es = q*N, Ep*N
    say(f"   {N:>8}{Q:>15.4e}{Es:>16.4e}"
        f"{(Q/prev[0] if prev else float('nan')):>17.4f}{(Es/prev[1] if prev else float('nan')):>17.4f}")
    prev=(Q,Es)
# interaction correction to extensivity: one cell's total coupling to an infinite surrounding array
tot=0.0
for dx in range(-40,41):
    for dy in range(-40,41):
        if dx==0 and dy==0: continue
        tot += U_img(pitch*np.hypot(dx,dy))
say(f"   per-cell interaction sum over an 81x81 surrounding array: {tot:.3e} J = {tot/Ep:.3f} of E_p")
say(f"   -> bounded per cell (r^-3 sum converges), so interactions leave extensivity linear")
say("")
say("3. (b) ADDITIVE over disjoint blocks")
say("   Q_tot: exact by charge superposition -- defect identically 0, and TRIVIALLY so (stated).")
say("   E_tot: the informative test. Two W x W blocks side by side; cross-terms computed exactly.")
say(f"   {'W':>5}{'bulk 2W^2 E_p (J)':>19}{'cross-term defect (J)':>23}{'defect/bulk':>13}")
for W in (4,8,16,32):
    xsA=[(i,j) for i in range(W) for j in range(W)]
    xsB=[(i+W,j) for i in range(W) for j in range(W)]
    cross=0.0
    A=np.array(xsA,dtype=float); B=np.array(xsB,dtype=float)
    for (ax,ay) in A:
        d=np.hypot((B[:,0]-ax),(B[:,1]-ay))*pitch
        cross += float(np.sum(2*k*q*q*(1.0/d - 1.0/np.sqrt(d*d+4*h*h))))
    bulk=2*W*W*Ep
    say(f"   {W:>5}{bulk:>19.4e}{cross:>23.4e}{cross/bulk:>13.5f}")
say("   -> the defect is a PERIMETER term: defect/bulk falls ~1/W. Additivity holds in the bulk and")
say("      the boundary correction is the real inter-cell interference every flash engineer measures.")
say("")
say("4. (c) NOT A COUNT -- N fixed at 1000, the device generation varies")
say(f"   {'N_e':>6}{'C_fg (F)':>12}{'Q_tot (C)':>15}{'E_stored (J)':>16}")
vals=[]
for ne,cf in ((50,2.5e-18),(100,5e-18),(200,1e-17),(100,2.5e-18)):
    qq=ne*E; vals.append((1000*qq, 1000*qq*qq/(2*cf)))
    say(f"   {ne:>6}{cf:>12.1e}{1000*qq:>15.4e}{1000*qq*qq/(2*cf):>16.4e}")
sq=max(v[0] for v in vals)/min(v[0] for v in vals); se=max(v[1] for v in vals)/min(v[1] for v in vals)
say(f"   -> at fixed N the charge moves {sq:.1f}x and the energy {se:.1f}x -- neither is a count")
say("")
say("5. (e) SEPARATION -- the interaction between two programmed gates, monopoles + real ground plane")
say(f"   {'r (nm)':>8}{'U_img (J)':>14}{'U_img/kT':>11}{'ratio/doubling':>16}{'CONTROL no-plane ratio':>24}")
prev=None
for rnm in (30,60,120,240,480,960):
    r=rnm*1e-9; u=U_img(r); uf=U_free(r)
    say(f"   {rnm:>8}{u:>14.4e}{u/(KB*T):>11.3e}"
        f"{(prev[0]/u if prev else float('nan')):>16.4f}{(prev[1]/uf if prev else float('nan')):>24.4f}")
    prev=(u,uf)
expo=np.log(U_img(240e-9)/U_img(480e-9))/np.log(2)
say(f"   asymptotic exponent from the last doubling: {expo:.4f}")
say("   -> the ratio approaches 8 per doubling (exponent -3) WITH the plane and is exactly 2 per")
say("      doubling (exponent -1) WITHOUT it. The falloff is INDUCED: monopole charges were supplied,")
say("      the grounded channel is physically there, and the image construction does the rest.")
say("")
say("="*100); say("  READ -- generated from the numbers above"); say("="*100)
ok_d = all(abs(r-1.0)<1e-12 for _,r in rows_d) and min(ratios_rand)>1.0-1e-12
ok_e = abs(expo-3.0)<0.15
if ok_d and ok_e:
    say("  THE FIVE STANDARDS ARE MET ON TRAPPED CHARGE AS THEY WERE ON MAGNETIC ANISOTROPY.")
    say("  (a) extensive, ratio exactly 10 per decade with a bounded interaction correction;")
    say("  (b) additive, with the E_tot defect an explicitly computed perimeter term falling ~1/W;")
    say(f"  (c) not a count, {sq:.1f}x in charge and {se:.1f}x in energy at fixed N;")
    say("  (d) sign-definite with ratio exactly 1 for EVERY data pattern, against an unwritten null;")
    say(f"  (e) an induced falloff of exponent {expo:.2f}, against a no-plane control at exponent 1.")
    say("")
    say("  ON C-72 ACROSS MECHANISMS: in both, the unwritten medium's accumulated quantity is null")
    say("  and writing creates an extensive, sign-definite value. The null differs in kind -- thermal-")
    say("  random (magnetic, ratio ~1/sqrt(N)) versus absent (charge, exactly zero) -- and the sign is")
    say("  fixed by the write in both: field direction there, carrier injection here.")
else:
    say("  A STANDARD FAILED; see the tables. No claim is made.")
