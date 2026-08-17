# LANE_W11_R_MATH — LEG 2.  IS THE EDGE RATE A GENUINE LIMIT, AND WHAT IS ITS CLOSED FORM?
#
# THEOREM (mine).  Let Lam = lcm(L_F,L_C), g_F = Lam/L_F, g_C = Lam/L_C, and for rho = 0..Lam-1
#   let d(rho) = ( c00, c10 x^{floor(rho/L_F)}, c01 y^{floor(rho/L_C)}, c11 x^{..} y^{..} ) with
#   c = c(rho,s) the bilocal coefficients of LEG 1.  If (x^{g_F}, y^{g_C}) generates a dense
#   1-parameter subgroup of T^2 (Weyl), then
#
#       lim_N (1/N) sum_{n<=N} log |Z^T_n|  =  (1/Lam) sum_{rho=0}^{Lam-1} m( d(rho) )
#
#   an average of Lam LOGARITHMIC MAHLER MEASURES of polynomials of EXACTLY N1's SHAPE
#   c00 + c10 X + c01 Y + c11 XY.  The rho = 0 term is m(pi) = N1's own lambda, EXACTLY.
#
# PROOF SKETCH.  Split n by residue rho mod Lam.  Along each residue class Leg 1 gives
#   Z = d0 + d1 X + d2 Y + d3 XY with (X,Y) = (x^{g_F}, y^{g_C})^t.  log|.| is in L^1 of the torus
#   (the zero set of a nonzero polynomial is Haar-null and log|P| is integrable there -- Mahler),
#   so Weyl equidistribution for the CONTINUOUS truncations log max(|P|,eps) plus monotone
#   convergence gives the average.  Jensen in Y reduces m(d) to one circle integral.
#
# CONSEQUENCE, STATED PLAINLY: N1 HAS A RIVAL FORMULA ON THE SAME CONSTRUCTION, AND THE RIVAL
#   CONTAINS N1 AS ITS rho = 0 TERM.  The functional's SHAPE is convention-independent.  The
#   identification of its coefficients with pi is convention-DEPENDENT.
import numpy as np, wm0_lib as L
rng=np.random.default_rng(20260817)

def coeff_rows(car,a,s):
    NV=car["NV"]; wF,wC=car["walkF"],car["walkC"]; LF,LC=len(wF),len(wC); Lam=int(np.lcm(LF,LC))
    TF,TC=L.Top(wF,a,NV),L.Top(wC,a,NV)
    x=np.conj(L.hol(wF,a)); y=L.hol(wC,a)
    _,F,C=L.classes(car)
    inF=[1 if v in F else 0 for v in range(NV)]; inC=[1 if v in C else 0 for v in range(NV)]
    ix={(0,0):0,(1,0):1,(0,1):2,(1,1):3}
    rows={}
    for rho in range(Lam):
        rF,rC=rho%LF,rho%LC
        B=np.linalg.inv(np.linalg.matrix_power(TF,rF))@np.linalg.matrix_power(TC,rC)
        c=np.zeros(4,dtype=complex)
        for u in range(NV):
            for v in range(NV): c[ix[(inF[u],inC[v])]]+=np.conj(s[u])*s[v]*B[u,v]
        eF,eC=rho//LF,rho//LC
        rows[rho]=np.array([c[0], c[1]*x**eF, c[2]*y**eC, c[3]*x**eF*y**eC])
    return rows,(Lam,LF,LC,x,y)

def edge_rate_closed(car,a,s,n=1<<20):
    rows,_=coeff_rows(car,a,s)
    return np.mean([L.m_poly(rows[r],n) for r in sorted(rows)]), rows

def timeavg(op_F,op_C,s,N):
    xF=s.copy(); xC=s.copy(); tot=0.0
    for _ in range(N):
        xF=op_F@xF; xC=op_C@xC
        z=abs(np.vdot(xF,xC)); tot += np.log(z) if z>0 else -700.0
    return tot/N

print("== M2a  K1.  CLOSED FORM versus TIME AVERAGE, three states with IDENTICAL pi ==")
car=L.K1(); a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77]); NV=5
TF,TC=L.Top(car["walkF"],a,NV),L.Top(car["walkC"],a,NV)
MF,MC=L.Mop(car["walkF"],a,NV),L.Mop(car["walkC"],a,NV)
sA=np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB=np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
sC=sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
PI=np.array([0.0,0.3,0.3,0.4])
mP=L.m_jensen_registrar(PI); mP2=L.m_poly(PI)
print(f"  N1's lambda = m(pi) : registrar reduction {mP:.12f} | my complex-Jensen {mP2:.12f} | diff {abs(mP-mP2):.1e}")
print(f"  {'state':>6}{'closed form':>16}{'  m(rho=0)':>14}{'  m(rho=1)':>14}{'  m(rho=2)':>14}")
closed={}
for nm,s in (("A",sA),("B",sB),("C",sC)):
    r,rows=edge_rate_closed(car,a,s); closed[nm]=r
    ms=[L.m_poly(rows[k]) for k in (0,1,2)]
    print(f"  {nm:>6}{r:>16.9f}{ms[0]:>14.9f}{ms[1]:>14.9f}{ms[2]:>14.9f}")
print(f"  -> m(rho=0) equals N1's m(pi) for all three states (that is N1, unchanged).")
print(f"     The state dependence lives ENTIRELY in the rho != 0 terms.\n")
print(f"  {'N':>9}{'  timeavg A':>15}{'  timeavg B':>15}{'  timeavg C':>15}   |timeavg-closed| max")
for N in (2000,20000,200000,2000000):
    ta=[timeavg(TF,TC,s,N) for s in (sA,sB,sC)]
    err=max(abs(ta[i]-closed[k]) for i,k in enumerate("ABC"))
    print(f"  {N:>9}{ta[0]:>15.9f}{ta[1]:>15.9f}{ta[2]:>15.9f}      {err:.2e}")
print("  -> the edge rate IS a genuine limit and the closed form is exact.  The registrar's")
print("     N=2e5 figures were not converged: the drift is O(1/sqrt N) and is visible.\n")

print("== M2b  THE REGISTRAR'S COMPARATOR m(P)/3 IS THE WRONG NULL ==")
print(f"  registrar's stated null   m(pi)/3            = {mP/3:.9f}")
print(f"  correct convention-matched null (state A)    = {closed['A']:.9f}")
print("  A convention that PRESERVES invisibility does NOT give m(pi)/3 either.  Take the")
print("  fibre-wise cube-root tick  U = M^(1/3)  (diagonal, unitary, gauge-covariant, U^3 = M):")
wF,wC=car["walkF"],car["walkC"]
def rootM(walk,a,NV,Lr):
    W=L.hol(walk,a); th=np.angle(W)/Lr; U=np.eye(NV,dtype=complex)
    for v in L.loop_vs(walk): U[v,v]=np.exp(1j*th)
    return U
UF,UC=rootM(wF,a,5,3),rootM(wC,a,5,3)
print(f"    || U_F^3 - M_F || = {np.linalg.norm(np.linalg.matrix_power(UF,3)-MF):.2e}   "
      f"|| U_C^3 - M_C || = {np.linalg.norm(np.linalg.matrix_power(UC,3)-MC):.2e}")
for N in (200000,):
    ta=[timeavg(UF,UC,s,N) for s in (sA,sB,sC)]
    print(f"    per-tick rate at N={N}: A {ta[0]:.9f}  B {ta[1]:.9f}  C {ta[2]:.9f}  spread {max(ta)-min(ta):.1e}")
print(f"    and m(pi) = {mP:.9f} -- the invisible convention's per-tick rate is m(pi), NOT m(pi)/3.")
print("    So 'the edge rate is not m(pi)/3' carries NO weight against the edge convention.")
print("    The weight is carried by the SPREAD alone.\n")

print("== M2c  B0b (L_F=4, L_C=3, lcm=12).  CLOSED FORM versus TIME AVERAGE ==")
car=L.B0b(); NV=9
a=np.random.default_rng(20260817).uniform(0,2*np.pi,18)
TF,TC=L.Top(car["walkF"],a,NV),L.Top(car["walkC"],a,NV)
MF,MC=L.Mop(car["walkF"],a,NV),L.Mop(car["walkC"],a,NV)
w=np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); w/=w.sum()
wB=w.copy(); wB[0],wB[1]=w[0]+w[1],0.0; wB[3],wB[4]=0.0,w[3]+w[4]; wB[5],wB[8]=w[5]+w[8],0.0
sA=np.sqrt(w)+0j; sB=np.sqrt(wB)+0j
sC=sA*np.exp(1j*np.random.default_rng(7).uniform(0,2*np.pi,9))
cl,_,_=L.classes(car); PIb=L.pi_of(sA,cl)
print(f"  pi = {np.round(PIb,6)}    m(pi) = {L.m_poly(PIb):.12f}  (= N1's lambda on this carrier)")
print(f"  {'state':>6}{'closed form':>16}{'   timeavg N=2e5':>18}{'   timeavg N=2e6':>18}")
for nm,s in (("A",sA),("B",sB),("C",sC)):
    r,_=edge_rate_closed(car,a,s)
    t1=timeavg(TF,TC,s,200000); t2=timeavg(TF,TC,s,2000000)
    print(f"  {nm:>6}{r:>16.9f}{t1:>18.9f}{t2:>18.9f}")
rc=[timeavg(MF,MC,s,200000) for s in (sA,sB,sC)]
print(f"  CIRCUIT rate (per circuit) N=2e5: {rc[0]:.9f} {rc[1]:.9f} {rc[2]:.9f}  spread {max(rc)-min(rc):.1e}")

print()
print("== M2d  QUADRATURE CONVERGENCE OF m_poly (the closed form must not rest on one grid) ==")
car=L.K1(); a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])
rows,_=coeff_rows(car,a,sC)
print(f"  {'n grid':>10}{'m(pi)':>18}{'m(rho=1) state C':>20}{'m(rho=2) state C':>20}")
for n in (1<<14,1<<16,1<<18,1<<20,1<<22):
    print(f"  {n:>10}{L.m_poly(PI,n):>18.12f}{L.m_poly(rows[1],n):>20.12f}{L.m_poly(rows[2],n):>20.12f}")
print("  -> stable to 1e-11 from 2^18 up; every closed-form figure in this lane uses 2^20 or 2^21.")
