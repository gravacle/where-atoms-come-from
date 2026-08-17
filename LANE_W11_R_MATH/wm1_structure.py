# LANE_W11_R_MATH — LEG 1.  THE STRUCTURE THEOREM FOR THE EDGE SYSTEM.
#
# CLAIM (mine, proved below and then verified numerically to machine precision):
#   Let n = q_F L_F + r_F = q_C L_C + r_C, x = conj(W_F), y = W_C, and
#       B(n) := T_F^{-r_F} T_C^{r_C}      (depends on n only through rho = n mod lcm(L_F,L_C))
#   Then, EXACTLY,
#       Z^T_n = <T_F^n s, T_C^n s> = c00 + c10 x^{q_F} + c01 y^{q_C} + c11 x^{q_F} y^{q_C}
#   with     c_ab(rho,s) = sum_{u: [u in gF]=a} sum_{v: [v in gC]=b} conj(s_u) s_v B(n)_{uv}.
#
# PROOF.  T_F^{L_F} = M_F and M_F commutes with T_F, so T_F^n = T_F^{r_F} M_F^{q_F} and
#   T_F^{-n} = M_F^{-q_F} T_F^{-r_F};  likewise T_C^n = T_C^{r_C} M_C^{q_C}.  Hence
#   Z_n = s^* M_F^{-q_F} T_F^{-r_F} T_C^{r_C} M_C^{q_C} s.  M_F^{-q_F} is diagonal with entry
#   x^{q_F} on gamma_F and 1 off it;  M_C^{q_C} is diagonal with y^{q_C} on gamma_C and 1 off it.
#   Expanding the two diagonals and collecting the four (a,b) blocks gives the display.  QED
#
# WHY THIS IS THE WHOLE QUESTION.  At rho = 0, B = I, the double sum collapses to u = v and
#   c_ab = sum_{v in class (a,b)} |s_v|^2 = pi_ab.  N1's polynomial IS the rho = 0 block.
#   At rho != 0 the index (a,b) is carried by a PAIR (u,v) of DIFFERENT vertices: the coefficients
#   are a BILOCAL pushforward, sensitive to within-class distribution and to phase.
#   The SHAPE of N1 survives the convention exactly.  The IDENTIFICATION of its coefficients
#   with pi does not.
import numpy as np, wm0_lib as L
np.set_printoptions(linewidth=200)
rng=np.random.default_rng(20260817)

def legs(car, a, s, NMAX):
    NV=car["NV"]; wF,wC=car["walkF"],car["walkC"]; LF,LC=len(wF),len(wC)
    TF,TC=L.Top(wF,a,NV),L.Top(wC,a,NV)
    MF,MC=L.Mop(wF,a,NV),L.Mop(wC,a,NV)
    x=np.conj(L.hol(wF,a)); y=L.hol(wC,a)
    cl,F,C=L.classes(car)
    inF=np.array([1 if v in F else 0 for v in range(NV)])
    inC=np.array([1 if v in C else 0 for v in range(NV)])
    worst=0.0; worstpi=0.0; coeffs={}
    xF=s.copy(); xC=s.copy()
    for n in range(1,NMAX+1):
        xF=TF@xF; xC=TC@xC
        Z=np.vdot(xF,xC)                                   # direct, by repeated application
        qF,rF=divmod(n,LF); qC,rC=divmod(n,LC)
        B=np.linalg.inv(np.linalg.matrix_power(TF,rF))@np.linalg.matrix_power(TC,rC)
        c=np.zeros(4,dtype=complex)
        for u in range(NV):
            for v in range(NV):
                c[{(0,0):0,(1,0):1,(0,1):2,(1,1):3}[(inF[u],inC[v])]] += np.conj(s[u])*s[v]*B[u,v]
        Zp = c[0] + c[1]*x**qF + c[2]*y**qC + c[3]*x**qF*y**qC
        worst=max(worst,abs(Z-Zp))
        coeffs.setdefault(n % np.lcm(LF,LC), c)
    return worst, coeffs, L.pi_of(s,cl)

print("== M1  THE STRUCTURE THEOREM, VERIFIED ==")
print("   max | Z^T_n  -  [ c00 + c10 x^qF + c01 y^qC + c11 x^qF y^qC ] |   over n = 1..240")
for car,a,s in (
    (L.K1(),  np.array([1.0,0.37,0.91,2**0.5,0.23,1.77]),
              np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j),
    (L.B0b(), np.random.default_rng(20260817).uniform(0,2*np.pi,18),
              None),
):
    if s is None:
        w=np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); w/=w.sum(); s=np.sqrt(w)+0j
    w,coeffs,pi=legs(car,a,s,240)
    LF,LC=len(car["walkF"]),len(car["walkC"]); Lam=int(np.lcm(LF,LC))
    print(f"  {car['name']:>4}  L_F={LF} L_C={LC} lcm={Lam:>2}   worst residual = {w:.2e}")
    c0=coeffs[0]
    print(f"        rho=0 coefficients   c = {np.round(c0,12)}")
    print(f"        pushforward          pi = {np.round(pi,12)}")
    print(f"        || c(rho=0) - pi ||  = {np.linalg.norm(c0-pi):.2e}      <-- N1's polynomial IS the rho=0 block")
    for r in sorted(coeffs)[1:4]:
        print(f"        rho={r}  c = {np.round(coeffs[r],6)}   ||c-pi|| = {np.linalg.norm(coeffs[r]-pi):.3e}")
    print()

print("== M1b  THE COEFFICIENTS ARE BILOCAL: they see WITHIN-CLASS SPLIT and PHASE, pi does not ==")
car=L.K1(); a=np.array([1.0,0.37,0.91,2**0.5,0.23,1.77]); cl,F,C=L.classes(car)
sA=np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB=np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
sC=sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
for nm,s in (("A",sA),("B",sB),("C",sC)):
    _,co,pi=legs(car,a,s,12)
    print(f"  state {nm}:  pi = {np.round(pi,6)}")
    for r in (0,1,2):
        print(f"            rho={r}  c = {np.round(co[r],9)}")
print("  -> rho=0 row identical across A,B,C (= pi).  rho=1,2 rows differ.  That is the whole effect.")
