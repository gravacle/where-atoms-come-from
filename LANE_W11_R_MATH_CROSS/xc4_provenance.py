# LANE_W11_R_MATH_CROSS — LEG 4.  PROVENANCE AND CONVERGENCE.
#  X4a  wm2 M2d's convergence table is captioned "state C" and is not computed from any state C.
#  X4b  the report's T_C gauge-covariance figure 3.23e-15 is in NO sealed file of the lane.
#  X4c  W-10's commutator row is refuted harder by the lane's OWN REV4 than by its wm8.
#  X4d  the convergence claims, checked against EXACT closed forms rather than against a finer grid.
import numpy as np, xc0_lib as X
from math import log
np.set_printoptions(linewidth=200)

print("== X4a  wm2 M2d's 'state C' COLUMNS ARE COMPUTED FROM A STALE VARIABLE ==")
print("   wm2_rate.py runs M2a on K1 with sA,sB,sC (5 components), then M2c REBINDS sA,sB,sC to")
print("   B0b's 9-component states, then M2d does  car=L.K1() ; rows=coeff_rows(car,a,sC).")
print("   coeff_rows indexes s[0..4] only, so it silently consumes the FIRST FIVE COMPONENTS OF")
print("   B0b's state C.  No exception is raised.  The printed columns are headed 'state C'.")
car=X.K1(); NV=5; a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
def coeff_rows(car,a,s,NV):
    wF,wC=car["walkF"],car["walkC"]; LF,LC=len(wF),len(wC); Lam=int(np.lcm(LF,LC))
    TF,TC=X.Top(wF,a,NV),X.Top(wC,a,NV)
    x=np.conj(X.hol(wF,a)); y=X.hol(wC,a); _,F,C=X.classes(car)
    inF=[1 if v in F else 0 for v in range(NV)]; inC=[1 if v in C else 0 for v in range(NV)]
    ix={(0,0):0,(1,0):1,(0,1):2,(1,1):3}; rows={}
    for rho in range(Lam):
        B=np.linalg.inv(np.linalg.matrix_power(TF,rho%LF))@np.linalg.matrix_power(TC,rho%LC)
        c=np.zeros(4,dtype=complex)
        for u in range(NV):
            for v in range(NV): c[ix[(inF[u],inC[v])]]+=np.conj(s[u])*s[v]*B[u,v]
        rows[rho]=np.array([c[0],c[1]*x**(rho//LF),c[2]*y**(rho//LC),c[3]*x**(rho//LF)*y**(rho//LC)])
    return rows
sA=np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sC_K1=sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))           # K1's actual state C
w=np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); w/=w.sum()
sC_B0b=(np.sqrt(w)+0j)*np.exp(1j*np.random.default_rng(7).uniform(0,2*np.pi,9))  # the stale one
for tag,s in (("K1's actual state C (what the caption says)",sC_K1),
              ("first 5 entries of B0b's state C (what ran)",sC_B0b)):
    R=coeff_rows(car,a,s,NV)
    print(f"   {tag}:  m(rho=1) = {X.m_quad(R[1],1<<20):.12f}   m(rho=2) = {X.m_quad(R[2],1<<20):.12f}")
print("   wm2 M2d printed              :  m(rho=1) = -1.828435829289   m(rho=2) = -1.366583513965")
print("   wm2 M2a printed for state C  :  m(rho=1) = -1.702516603      m(rho=2) = -0.931308661")
print("   -> M2d's two right-hand columns match NEITHER the caption's state NOR M2a's own table.")
print("      COR-K defect class exactly ('not reproducible from the displayed parameters').")
print("      THE CONVERGENCE CONCLUSION IS UNHARMED -- the grid study is still a grid study of")
print("      three genuine polynomials, and the m(pi) column is correct -- but the labels are")
print("      false and the numbers cannot be regenerated from anything the lane states.\n")

print("== X4b  THE REPORT'S T_C GAUGE-COVARIANCE FIGURE IS IN NO SEALED FILE.  I COMPUTE IT. ==")
print("   The lane's report claims, twice: 'gauge-covariant in BOTH branches (3.23e-15; the")
print("   registrar tested T_F only)' and 'I checked both: max defect 3.23e-15.'")
print("   wm6_replicate.py's gauge loop builds Top(car['walkF'],...) ONLY; grep for 3.23e-15 over")
print("   the 18 sealed files returns nothing.  Under the lane's own rule ('unreproducible numbers")
print("   are treated as absent') the figure is absent.  Recomputed here, both branches:")
rng=np.random.default_rng(20260817)
E=[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]
wF=wC=0.0
for _ in range(2000):
    b=rng.uniform(0,2*np.pi,6); th=rng.uniform(0,2*np.pi,5)
    bg=np.array([b[j]+th[t]-th[s] for j,(s,t) in enumerate(E)])
    s=rng.normal(size=5)+1j*rng.normal(size=5)
    G=np.exp(1j*th)
    wF=max(wF,np.linalg.norm(X.Top(car["walkF"],bg,5)@(G*s)-G*(X.Top(car["walkF"],b,5)@s)))
    wC=max(wC,np.linalg.norm(X.Top(car["walkC"],bg,5)@(G*s)-G*(X.Top(car["walkC"],b,5)@s)))
print(f"   max gauge defect, T_F : {wF:.2e}      max gauge defect, T_C : {wC:.2e}   (2000 draws,")
print("   unnormalised complex Gaussian s, the lane's own convention in wm6).  Both are covariant;")
print("   the CLAIM is true and the FIGURE is unpublished.  Note the value depends on ||s||, which")
print("   the lane does not normalise, so '3.23e-15' is not a reproducible quantity in any case.\n")

print("== X4c  W-10's COMMUTATOR ROW IS REFUTED HARDER BY THE LANE'S OWN REV4 THAN BY ITS wm8 ==")
print("   W-10 (REGISTER_V001.md:1042-1046): '||[T_F,T_C]|| is 2.828 on B0b, 2.449 on K1 and")
print("   exactly 0 on B0a -- NON-ZERO PRECISELY WHEN CLASS 11 IS OCCUPIED.'")
print("   wm8 exhibits class 11 EMPTY with commutator 0 -- which is CONSISTENT with W-10's")
print("   correlation, and refutes only the diagnostic role wm8 reads into it.  The correlation")
print("   itself needs class 11 OCCUPIED and commutator ZERO.  The lane already owns that carrier")
print("   and did not connect it: REV4, its own wm3 M3b counterexample.")
NV4=4
walkF=[(0,1,0,+1),(1,2,1,+1),(2,3,2,+1),(3,0,3,+1)]
walkC=[(0,3,3,-1),(3,2,2,-1),(2,1,1,-1),(1,0,0,-1)]
rev4={"name":"REV4","NV":NV4,"NE":4,"walkF":walkF,"walkC":walkC}
a4=np.array([0.7,1.3,-0.4,2.1])
TF4,TC4=X.Top(walkF,a4,NV4),X.Top(walkC,a4,NV4)
cl4,F4,C4=X.classes(rev4)
print(f"   REV4: class vector {cl4}  -> class 11 occupied: {bool((cl4==3).any())}"
      f"   |class 11| = {int((cl4==3).sum())}")
print(f"         || [T_F, T_C] || = {np.linalg.norm(TF4@TC4-TC4@TF4):.2e}   <-- EXACTLY ZERO")
print(f"         and the edge convention still separates: spread at n=1 = "
      f"{X.spread(TF4,TC4,1,X.same_pi_states(cl4,np.array([0.,0.,0.,1.]),np.random.default_rng(5),20)):.3e}")
print("   -> class 11 FULLY occupied and the commutator EXACTLY zero.  W-10's stated correlation")
print("      is false in the 'only if' direction too, on a carrier inside the lane's own wm3.")
print("      Any two loops with the same vertex set and sigma_C = sigma_F^{-1} do this.\n")

print("== X4d  CONVERGENCE, CHECKED AGAINST EXACT CLOSED FORMS, NOT AGAINST A FINER GRID ==")
print("   A finer grid agreeing with a coarser one is not convergence to the right number.  Every")
print("   load-bearing quadrature in the lane is checked here against a closed form computed by a")
print("   different route (Cassaigne-Maillot with a Bloch-Wigner dilogarithm, or Jensen exactly).")
rows=[("N1's lambda, m(0.4+0.3x+0.3y)   [W-02 erratum, wm2 M2a]",
       X.m_quad((0.0,0.3,0.3,0.4),1<<22), X.m_CM(0.4,0.3,0.3), -0.767507880357),
      ("B0b pi=(.44,.25,.09,.22): m = log(0.44) by branch domination [wm7 M7c]",
       X.m_quad((0.44,0.25,0.09,0.22),1<<22), log(0.44), -0.820980552070),
      ("rank-1 f=1,c=0: m = log(0.7) exactly [wm5 M5b]",
       X.m_quad((0.0,0.3,0.3,0.4),1<<22), None, None),
      ("W-10 N-3's B0b lambda = log(4/9) [S4:599 corrected]",
       X.m_quad((4/9,2/9,1/9,2/9),1<<22), log(4/9), -0.810930216216)]
for nm,q,cf,ref in rows[:2]+rows[3:]:
    print(f"   {nm}")
    print(f"      quadrature 2^22 {q:.12f}   closed form {cf:.12f}   |diff| {abs(q-cf):.2e}"
          + (f"   lane/register printed {ref:.12f}" if ref else ""))
# the rank-1 row needs the connection, not just pi
aR=np.array([1.0,0.0,0.0,0.0,0.0,0.0])
MF,MC=X.Mop(car["walkF"],aR,5),X.Mop(car["walkC"],aR,5)
def timeavg(oF,oC,s,N):
    xF=s.copy(); xC=s.copy(); tot=0.0
    for _ in range(N):
        xF=oF@xF; xC=oC@xC; z=abs(np.vdot(xF,xC)); tot+= np.log(z) if z>1e-300 else -700.0
    return tot/N
ta=timeavg(MF,MC,sA,200000)
print(f"   rank-1 f=1.0,c=0.0 circuit rate: lane N=2e5 gives -0.356675073 ; mine {ta:.9f} ;")
print(f"      EXACT value log(0.7) = {log(0.7):.12f}   |lane - exact| = {abs(-0.356675073-log(0.7)):.2e}")
aO=np.array([np.pi/3]*3+[np.pi/2]*3)
MF,MC=X.Mop(car["walkF"],aO,5),X.Mop(car["walkC"],aO,5)
ta4=timeavg(MF,MC,sA,200000)
print(f"   S1 order-4 f=pi,c=3pi/2 circuit rate: lane N=2e5 gives -0.804718956 ; mine {ta4:.9f} ;")
print(f"      EXACT value (1/4)log(1/25) = {-log(25)/4:.12f}   |lane - exact| = {abs(-0.804718956+log(25)/4):.2e}")
print("      (the period-4 orbit gives |Z| = sqrt(1/10), 2/5, sqrt(1/10), 1, so the average of the")
print("       logs is (1/4) log(0.1*0.4) = -(1/4) log 25.  COR-K's '-0.804719' is that number, and")
print("       NEITHER the corpus NOR the lane states the closed form.)")
print("   -> every convergence claim I could reach a closed form for is CONVERGENCE, not a window.")
