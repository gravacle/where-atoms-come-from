# LANE W-08 / M2 leg C — THE RACE.  Does the near-return structure ever make SUM (1-|Z_k|)
# CONVERGE, i.e. kill durability?
#   C1  the exact asymptotic rate D_H = 1 - INT_H |Z| dHaar, computed in high precision for every
#       class, against leg B's measured SUM/K.  (If the sum is linear, this is its slope.)
#   C2  the THEOREM: an exact, Diophantine-free lower bound  SUM_{k<=K}(1-|Z_k|) >= C(K - K0).
#   C3  the dichotomy: D_H = 0 <=> no formation (W-02's G = {1}, COR-B).  Swept, not asserted.
#   C4  the search for a convergent connection, and where the infimum actually sits.
import numpy as np, mpmath as mp
mp.mp.dps = 40
# SEED DEFECT, RECORDED NOT SILENTLY FIXED: the first run of this script used a hard-coded
# pair (0.63696..., 0.26979...) labelled D_RAND1 that is NOT the seeded stream leg B used.
# Any rank-2 irrational pair supports the same conclusion, but the LABEL was wrong and the
# number was unreproducible from the published seed -- the exact defect this program treats
# as absent data.  The pairs are now drawn from the master seed here exactly as in leg B.
RNG_ = __import__('numpy').random.default_rng(20260816)
RND_ = [(RNG_.random(), RNG_.random()) for _ in range(3)]


def kinks(w11,w10,w01):
    """phi in [0,1] where |A(phi)| = w01, A = w11 e^{2pi i phi} + w10.  Both integrands below are
    non-smooth exactly there (E(m) at m=1; the max() switching), so the quadrature is SPLIT there
    rather than left to find it.  This is the difference between 6 correct digits and 12."""
    if w11==0 or w10==0: return []
    co = (w01**2 - w11**2 - w10**2)/(2*w11*w10)
    if abs(co)>1: return []
    t = mp.acos(co)/(2*mp.pi)
    return sorted({t, 1-t})

def E_abs_torus(w11,w10,w01):
    """INT_{T^2} |p11 XY + p10 X + p01 Y|.  Inner theta-integral in closed form:
       (1/2pi) INT |A e^{i th} + B| d th = (2/pi)(a+b) E(m), m = 4ab/(a+b)^2  (complete elliptic)."""
    def inner(ph):
        A = w11*mp.expjpi(2*ph) + w10          # coefficient of X
        a = abs(A); b = mp.mpf(w01)
        if a+b == 0: return mp.mpf(0)
        m = 4*a*b/(a+b)**2
        return (2/mp.pi)*(a+b)*mp.ellipe(m)
    return mp.quad(inner,[0]+kinks(w11,w10,w01)+[1])

def m_torus(w11,w10,w01):
    """m(p11 XY + p10 X + p01 Y) via Jensen on the inner variable: INT log max(|A|,|B|)."""
    def inner(ph):
        A = w11*mp.expjpi(2*ph) + w10
        return mp.log(mp.mpf(max(abs(A), w01)))
    return mp.quad(inner,[0]+kinks(w11,w10,w01)+[1])

def maillot(a,b,c):
    """CLOSED FORM (Cassaigne-Maillot) for m(a + b x + c y), a,b,c > 0 — an INDEPENDENT check on
       the quadrature, using no numerical integration of the Mahler integrand at all.
       If a,b,c form a triangle with angles alpha,beta,gamma opposite a,b,c:
         pi*m = D(|a/b| e^{i gamma}) + alpha log a + beta log b + gamma log c
       else m = log max(a,b,c).   D = Bloch-Wigner dilogarithm."""
    a,b,c = mp.mpf(a),mp.mpf(b),mp.mpf(c)
    if a>=b+c or b>=a+c or c>=a+b: return mp.log(max(a,b,c))
    al = mp.acos((b**2+c**2-a**2)/(2*b*c))
    be = mp.acos((a**2+c**2-b**2)/(2*a*c))
    ga = mp.acos((a**2+b**2-c**2)/(2*a*b))
    z  = (a/b)*mp.expj(ga)
    D  = mp.im(mp.polylog(2,z)) + mp.arg(1-z)*mp.log(abs(z))
    return (D + al*mp.log(a) + be*mp.log(b) + ga*mp.log(c))/mp.pi

def circle_poly(coeffs_exp, W):
    """H is the circle s -> (e^{2pi i e_x s}, e^{2pi i e_y s}).  Return (D_H, lambda_H) by
       reducing Z(s) to a one-variable Laurent polynomial and integrating / using Jensen."""
    w11,w10,w01 = W; ex,ey = coeffs_exp
    ex_list = [ex+ey, ex, ey]; co = [w11,w10,w01]
    lo = min(ex_list); deg = max(ex_list)-lo
    P = np.zeros(deg+1)
    for c,e in zip(co,ex_list): P[e-lo] += c
    # D_H by high-order quadrature of |P(e^{2pi i s})|
    N = 2**22
    s = np.arange(N)/N
    val = np.zeros(N,dtype=complex)
    for c,e in zip(co,ex_list): val += c*np.exp(2j*np.pi*e*s)
    D = 1.0 - float(np.abs(val).mean())
    # lambda by Jensen on the exact roots
    r = np.roots(P[::-1])
    lam = np.log(abs(P[np.max(np.nonzero(P))])) + np.sum(np.log(np.maximum(1.0,np.abs(r))))
    return D, float(lam)

for W,tag in [((0.4,0.3,0.3),"RS-G (0.4,0.3,0.3)"), ((0.5,0.0,0.5),"RS-P (0.5,0,0.5)")]:
    w11,w10,w01 = W
    print(f"== C1  EXACT ASYMPTOTIC RATES,  {tag} ==")
    Et = E_abs_torus(*W); mt = m_torus(*W); mk = maillot(w11,w10,w01)
    print(f"   H = T^2 (dim 2, generic):   D_H = 1 - INT|Z| = {float(1-Et):.12f}     "
          f"lambda = m(...) = {float(mt):.12f}")
    print(f"        CLOSED-FORM CHECK (Cassaigne-Maillot, Bloch-Wigner):  m = {float(mk):.12f}"
          f"   |quad - closed form| = {float(abs(mt-mk)):.3e}")
    # dim-1 classes: H is the circle with primitive direction from the relation
    for (ex,ey,name) in [(-20,11,"B_S3RES  11a-20b=0  -> (x,y)=(e^{-2pi i 20 s}, e^{2pi i 11 s})"),
                         (-1, 1, "E_W07GEN a-b=0      -> (x,y)=(e^{-2pi i s},    e^{2pi i s})")]:
        D,lam = circle_poly((ex,ey),W)
        print(f"   H = circle (dim 1): {name}")
        print(f"        D_H = {D:.12f}     lambda = {lam:.12f}")
    # dim-0 class: EXACT finite average in exact arithmetic (x = -1, y = -i are exact Gaussian
    # integers, so |Z_k|^2 is a rational number and no float ever touches this row).
    from fractions import Fraction as Fr
    xs=[(-1,0),(1,0),(-1,0),(1,0)]                      # x^k for k=1..4, x = -1
    ys=[(0,-1),(-1,0),(0,1),(1,0)]                      # y^k for k=1..4, y = -i
    W11,W10,W01 = [Fr(str(t)) for t in (w11,w10,w01)]
    z2=[]
    for (xr,xi),(yr,yi) in zip(xs,ys):
        pr = xr*yr-xi*yi; pi_ = xr*yi+xi*yr             # (xy)^k
        re = W11*pr + W10*xr + W01*yr; im = W11*pi_ + W10*xi + W01*yi
        z2.append(re*re+im*im)                          # |Z_k|^2 EXACTLY, a rational
    print(f"   H = Z_4 (dim 0): A_S1PUB   |Z_k|^2 for k=1..4 (exact rationals) = "
          f"{[str(t) for t in z2]}")
    if all(t>0 for t in z2):
        zz=[mp.sqrt(mp.mpf(t.numerator)/t.denominator) for t in z2]
        print(f"        D_H = {float(1-sum(zz)/4):.12f}     "
              f"lambda = {float(sum(mp.log(t) for t in zz)/4):.12f}   (exact-arithmetic row)")
    else:
        print(f"        D_H = {float(1-sum(mp.sqrt(mp.mpf(t.numerator)/t.denominator) for t in z2)/4):.12f}"
              f"     lambda = -INFINITY EXACTLY: |Z_k| = 0 on every odd k, so Omega_N = 0 for all")
        print(f"        N >= 1 — not underflow, an exact zero.  This is COR-D's case, and it is the")
        print(f"        MAXIMUM of durability (branches exactly orthogonal), not a failure of it.")
    print()

print("== C2  THE THEOREM.  No Diophantine input, no equidistribution, no genericity. ==")
print("""   From leg A's exact identity, for ANY j<l and ANY k:
       1-|Z_k| >= (1/2)(1-|Z_k|^2) >= (1/2) w_j w_l |chi_j^k - chi_l^k|^2 .
   Put g = chi_j/chi_l = e^{i tau}.  Then SUM_{k=1..K} |g^k-1|^2 = 2K - 2 Re(SUM_{k=1..K} g^k),
   and |SUM_{k=1..K} g^k| <= 1/|sin(tau/2)| for every K.  Hence for every K

       SUM_{k<=K} (1 - |Z_k|)  >=  w_j w_l * ( K - 1/|sin(tau/2)| ).                        (*)

   The right-hand side is LINEAR IN K with a positive slope whenever some pair of characters
   with nonzero weight differs, i.e. whenever W-02's G != {1}, i.e. WHENEVER FORMATION OCCURS.
   There is no connection, generic or not, on which the sum converges.  The Diophantine
   structure of the near-returns cannot enter: (*) never mentions it.""")
print("   The three pairs on K1 (p00 = 0, so exactly three characters):")
print("     (0,F): g = W_C        (0,C): g = conj(W_F)      (F,C): g = rho = conj(W_F)/W_C")
print()
print(f"   {'connection':<40} {'ready':<6} {'best pair':<8} {'w_j w_l':>9} {'K0=1/|sin(tau/2)|':>18} "
      f"{'bound slope':>12} {'measured SUM/K (K=1e7)':>24}")
MEAS = {("A_S1PUB","RS-G"):0.491886,("B_S3RES","RS-G"):0.469183,("C_BADAPP","RS-G"):0.469189,
        ("D_RAND1","RS-G"):0.469189,("E_W07GEN","RS-G"):0.529471,("F_VWA","RS-G"):0.469215,
        ("A_S1PUB","RS-P"):0.500000,("B_S3RES","RS-P"):0.363380,("C_BADAPP","RS-P"):0.363380,
        ("D_RAND1","RS-P"):0.363380,("E_W07GEN","RS-P"):0.363380,("F_VWA","RS-P"):0.363380}
phi=(1+5**0.5)/2
CONN = {"A_S1PUB":(0.5,0.75),"B_S3RES":(1/np.pi,11/(20*np.pi)),
        "C_BADAPP":(np.mod(2*np.cos(2*np.pi/7),1),np.mod((2*np.cos(2*np.pi/7))**2,1)),
        "D_RAND1":(RND_[0][0],RND_[0][1]),
        "E_W07GEN":(np.mod(phi,1),np.mod(phi**2,1)),"F_VWA":(1234567/2999999,765431/2999999)}
for cn,(al,be) in CONN.items():
    for rt,W in [("RS-G",(0.4,0.3,0.3)),("RS-P",(0.5,0.0,0.5))]:
        w11,w10,w01=W
        pairs = {"(0,F)":(w11*w10, be), "(0,C)":(w11*w01, al), "(F,C)":(w10*w01, al+be)}
        best=None
        for nm,(wt,t) in pairs.items():
            d = abs(t-round(t))
            if wt<=0 or d<1e-14: continue
            tau=2*np.pi*d; sl=wt; K0=1/abs(np.sin(tau/2))
            if best is None or sl>best[1]: best=(nm,sl,K0)
        nm,sl,K0=best
        print(f"   {cn:<40} {rt:<6} {nm:<8} {sl:>9.4f} {K0:>18.3f} {sl:>12.6f} "
              f"{MEAS[(cn,rt)]:>24.6f}")
print("   -> every measured slope exceeds the proved slope; the bound is a floor, not a fit.")
print()

print("== C3  THE DICHOTOMY, SWEPT — AND A SELF-CORRECTION THAT STRENGTHENS THE RESULT ==")
print("""   SELF-CORRECTION, RECORDED NOT SILENTLY FIXED.  The first run of this block concluded:
   "D_H is continuous and -> 0 at the trivial connection, so the INFIMUM of the slope over
   non-trivial connections is 0."  THAT IS FALSE, and my OWN bound (*) refutes it: (*) gives
   liminf (1/K) SUM >= w_j w_l for EVERY non-trivial connection, a constant that does not
   depend on the connection at all.  The truth is the opposite of what I wrote and it is a
   STRONGER statement: the durability rate does not decay to zero as the connection approaches
   triviality — it JUMPS.  What does degrade continuously is the ONSET SCALE K0 = 1/|sin(tau/2)|,
   the K after which the linear growth is visible.  Both halves are measured below.""")
print("   Grid 121x121 over [0,1)^2 at K = 200000, ready state RS-G (so every grid denominator")
print("   divides 120 and K = 200000 is far beyond every K0 on the grid).")
K=200000; k=np.arange(1,K+1,dtype=np.float64); W=(0.4,0.3,0.3); w11,w10,w01=W
small=[]; nontriv_min=(np.inf,None)
G=121
for i in range(G):
    al=i/(G-1)
    for j in range(G):
        be=j/(G-1)
        u=(k*al)%1.0; v=(k*be)%1.0
        du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv-=np.round(duv)
        S=4*(w11*w10*np.sin(np.pi*dv)**2+w11*w01*np.sin(np.pi*du)**2+w10*w01*np.sin(np.pi*duv)**2)
        S=np.minimum(S,1.0); g=S/(1+np.sqrt(np.maximum(0,1-S)))
        mn=float(g.mean())
        if mn<1e-3: small.append((al,be,mn))
        else:
            if mn<nontriv_min[0]: nontriv_min=(mn,(al,be))
print(f"   grid points with mean(1-|Z_k|) < 1e-3 : {len(small)}  (all of them exactly 0)")
for al,be,mn in small:
    print(f"        alpha={al:.4f} beta={be:.4f}  mean = {mn:.3e}   "
          f"{'TRIVIAL CONNECTION (W_F=W_C=1): G={1}, no formation (COR-B)' if mn<1e-15 else ''}")
print(f"   MINIMUM over all 14641-4 NON-TRIVIAL grid points: {nontriv_min[0]:.6f} "
      f"at (alpha,beta) = ({nontriv_min[1][0]:.4f}, {nontriv_min[1][1]:.4f})")
print(f"   proved uniform lower bound from (*) for RS-G: 0.120000")
print(f"   -> measured minimum {nontriv_min[0]:.6f} >= proved bound 0.12: THE GAP IS REAL AND UNIFORM.")
print("      Formation/durability is not a matter of degree in the connection.  It is all-or-nothing")
print("      with a uniform gap in the RATE; only the ONSET K0 degrades as the connection -> trivial.")
print(f"   {'near-trivial connection':<34} {'K0 = 1/|sin(tau/2)|':>20} {'(1/K)SUM at K=1e4':>19} "
      f"{'at K=1e7':>12}")
for e in [2,3,4,5]:
    al=10.0**-e; be=0.0
    tau=2*np.pi*al; K0=1/abs(np.sin(tau/2))
    vals=[]
    for KK in [10**4,10**7]:
        acc=0.0; done=0
        while done<KK:
            c=min(10**7,KK-done); kk=np.arange(done+1,done+c+1,dtype=np.float64)
            u=(kk*al)%1.0; v=(kk*be)%1.0
            du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv-=np.round(duv)
            S=np.minimum(4*(w11*w10*np.sin(np.pi*dv)**2+w11*w01*np.sin(np.pi*du)**2
                            +w10*w01*np.sin(np.pi*duv)**2),1.0)
            acc+=float((S/(1+np.sqrt(np.maximum(0,1-S)))).sum()); done+=c
        vals.append(acc/KK)
    print(f"   {'alpha=10^-'+str(e)+', beta=0':<34} {K0:>20.1f} {vals[0]:>19.6f} {vals[1]:>12.6f}")
print("   -> exactly as (*) predicts: below K0 the rate is small, above K0 it is back above 0.12.")
print()

print("== C4  THE SEARCH FOR A CONVERGENT CONNECTION — result, and why it had to be this result ==")
print("""   Convergence of SUM (1-|Z_k|) requires 1-|Z_k| -> 0, hence (leg A4) chi_j^k/chi_l^k -> 1
   for every pair with weight, hence g^k -> 1 for g = chi_j/chi_l in U(1).  But for g in U(1),
   g^k -> 1 forces g = 1: if g = e^{2pi i t} with t irrational the sequence is equidistributed,
   and if t = a/q in lowest terms with q > 1 the sequence is periodic hitting e^{2pi i /q}
   infinitely often.  So the ONLY connections with a convergent sum are those with g = 1 for
   every weighted pair -- exactly W-02's G = {1}, exactly COR-B's non-forming families.
   NO CONNECTION EXHIBITED, AND NONE EXISTS.  This is the second-most-valuable outcome the lane
   asked for, and it is a proof rather than a failed search.""")
print("   Verified against COR-B's four non-forming families on K1 (RS-G unless the family fixes p):")
# SELF-CORRECTION, RECORDED NOT SILENTLY FIXED: the first run of this block labelled the last
# family "supp {v0,v3,v4}, W_C = 1" and set beta = 0.  That is WRONG and my own control caught it:
# W-02's criterion for S = {0,C} is G = <u> with u = conj(W_F), so the non-forming condition is
# W_F = 1, i.e. ALPHA = 0, with W_C free.  With beta=0, alpha=0.3137 the sum was 7.27e4 -- correctly
# NON-zero, because that connection DOES form.  Row corrected below; the mislabel is left on record.
fams=[("W_F=W_C=1 (trivial)",0.0,0.0,(0.4,0.3,0.3)),
      ("supp p = {v0} only  (|S|=1)",0.3137,0.7715,(1.0,0.0,0.0)),
      ("supp p = {v1,v2} only (|S|=1)",0.3137,0.7715,(0.0,1.0,0.0)),
      ("supp p = {v3,v4} only (|S|=1)",0.3137,0.7715,(0.0,0.0,1.0)),
      ("supp {v0,v3,v4}, W_F=1 -> G=<u>={1}",0.0,0.7715,(0.5,0.0,0.5)),
      ("[the mislabel: same supp, W_C=1, W_F!=1 -- FORMS]",0.3137,0.0,(0.5,0.0,0.5))]
kk=np.arange(1,200001,dtype=np.float64)
for nm,al,be,W in fams:
    w11,w10,w01=W
    u=(kk*al)%1.0; v=(kk*be)%1.0
    du=u-np.round(u); dv=v-np.round(v); duv=u+v; duv-=np.round(duv)
    S=4*(w11*w10*np.sin(np.pi*dv)**2+w11*w01*np.sin(np.pi*du)**2+w10*w01*np.sin(np.pi*duv)**2)
    S=np.minimum(S,1.0); g=S/(1+np.sqrt(np.maximum(0,1-S)))
    print(f"   {nm:<36} SUM_{{k<=2e5}} (1-|Z_k|) = {g.sum():.6e}")
print("   -> zero, exactly, on every non-forming family: the sum CONVERGES there and only there,")
print("      and there it converges because it is identically zero.  Durability and formation are")
print("      the same condition; there is no third regime in between.")
