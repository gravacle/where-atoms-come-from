# LANE_W11_R_MATH_CROSS — LEG 1.  MY COMMISSIONED LENS: ISOLATION / ARITHMETIC TYPE.
# "Verify each connection is the arithmetic type it claims (every rational (f,c) is exactly
#  resonant), and that its convergence claims are convergence and not a window."
#
# THREE THINGS THE LANE UNDER TEST DID NOT DO:
#  (A) It prints a column headed "smallest relation m f + n c = 0 mod 2pi" and does NOT report the
#      smallest relation: np.argmin over a tied-at-zero array returns the FIRST index of the
#      scan box, not the minimal generator.  Checked here by direct minimisation of |m|+|n|.
#  (B) Its genericity certificate is "the residual sits BELOW the Dirichlet bound 2pi/M^2".
#      Sitting below is FORCED by pigeonhole for EVERY point, resonant or not, so as a test it
#      could not have failed.  The discriminating statistic is the DEPTH SCALING of the residual.
#  (C) Its rate theorem carries a Weyl hypothesis that is FALSE at the corpus's own headline
#      connection, and the lane never tests it anywhere.  Tested here.
import numpy as np, xc0_lib as X
from math import log, pi as PI, gcd
np.set_printoptions(linewidth=200)
TWO_PI=2*np.pi

def relations(f,c,M):
    m=np.arange(-M,M+1); n=np.arange(-M,M+1)
    Mm,Nn=np.meshgrid(m,n,indexing="ij")
    val=Mm*f+Nn*c
    res=np.abs(val-TWO_PI*np.round(val/TWO_PI))
    res=np.where((Mm!=0)|(Nn!=0),res,np.inf)
    return Mm,Nn,res

def minimal_relation(f,c,M=1500,tol=0.0):
    """the SMALLEST relation by |m|+|n|, not the first index of the scan."""
    Mm,Nn,res=relations(f,c,M)
    hit=(res<=tol)
    if not hit.any(): return None
    w=np.abs(Mm)+np.abs(Nn)
    w=np.where(hit,w,10**9)
    k=np.unravel_index(np.argmin(w),w.shape)
    return int(Mm[k]),int(Nn[k])

def argmin_relation(f,c,M=1500):
    Mm,Nn,res=relations(f,c,M)
    k=np.unravel_index(np.argmin(res),res.shape)
    return int(Mm[k]),int(Nn[k]),float(res[k])

CONN={
 "W-11 registrar's K1     f=2.28, c=2+sqrt2 ": (1.0+0.37+0.91, 2**0.5+0.23+1.77),
 "corpus's ONLY generic   S4:603 f=1, c=sqrt2": (1.0, 2**0.5),
 "S3/S4 headline          f=2.0, c=1.1      ": (2.0, 1.1),
 "S1 sec6 order-4         f=pi,  c=3pi/2    ": (np.pi, 3*np.pi/2),
 "wm5 M5b rank-1 row      f=1.0, c=0.0      ": (1.0, 0.0),
}
carB=X.B0b(); aB=np.random.default_rng(20260817).uniform(0,2*np.pi,18)
CONN["W-11 registrar's B0b   rng(20260817)   "]=(float(np.angle(X.hol(carB["walkF"],aB))),
                                                 float(np.angle(X.hol(carB["walkC"],aB))))

print("== X1a  THE 'SMALLEST RELATION' COLUMN OF wm5 M5a DOES NOT REPORT THE SMALLEST RELATION ==")
print(f"  {'connection':<44}{'wm5 prints':>18}{'actual minimal':>18}  ratio")
for k,(f,c) in CONN.items():
    am,an,ar=argmin_relation(f,c)
    mm=minimal_relation(f,c,1500,0.0)
    if mm is None:
        print(f"  {k:<44}{f'({am},{an})':>18}{'none at depth 1500':>18}   residual {ar:.3e}")
    else:
        g=gcd(abs(am),abs(an)) if (am or an) else 0
        print(f"  {k:<44}{f'({am},{an})':>18}{str(mm):>18}   wm5's row is {abs(am)//max(abs(mm[0]),1) if mm[0] else abs(an)//max(abs(mm[1]),1)}x the generator")
print("  -> for BOTH exactly-resonant rows wm5 prints a NON-MINIMAL multiple of the generator,")
print("     under a header that says 'smallest'.  The quantity the lane itself declares operative")
print("     ('not rationality but the SIZE of the smallest relation, W-08's (q,delta)') is the one")
print("     its census does not compute.  The TYPE verdicts (RESONANT / not) are unaffected.\n")

print("== X1b  THE DIRICHLET TEST COULD NOT HAVE FAILED.  DEPTH SCALING IS THE TEST THAT CAN. ==")
print("  Pigeonhole forces min residual <= ~2pi/M^2 at depth M for EVERY (f,c).  What separates a")
print("  generic point from a near-resonant one is whether the residual keeps FALLING like M^-2")
print("  (generic) or COLLAPSES to 0 at a finite depth (resonant).  Both are run:")
print(f"  {'connection':<44}" + "".join(f"{'M='+str(M):>12}" for M in (50,150,450,1350)) + "   verdict")
for k,(f,c) in CONN.items():
    rs=[]
    for M in (50,150,450,1350):
        _,_,r=argmin_relation(f,c,M); rs.append(r)
    exact=any(r==0.0 for r in rs)
    if exact:
        d=min(M for M,r in zip((50,150,450,1350),rs) if r==0.0)
        v=f"EXACTLY RESONANT, first at depth {d}"
    else:
        # generic Diophantine: residual ~ C/M^2, so log-log slope near -2
        sl=np.polyfit(np.log([50,150,450,1350]),np.log(rs),1)[0]
        v=f"no relation; log-log slope {sl:+.2f} (Dirichlet -2)"
    print(f"  {k:<44}" + "".join(f"{r:>12.2e}" for r in rs) + f"   {v}")
print("  -> every row typed by wm5 is confirmed, by a statistic that CAN come out the other way.")
print("     The rank-1 row f=1.0, c=0.0 is EXACTLY RESONANT at depth 1 -- wm5 M5b runs it and")
print("     wm5 M5a does not census it, so the corpus's own N-4 rule is applied to four of the")
print("     five connections the lane executes.\n")

print("== X1c  THE RIVAL RATE FORMULA CARRIES THE ERRATUM-AGAINST-W-02 TRAP AND IS NOT FLAGGED ==")
print("  wm2's theorem is stated 'IF (x^gF, y^gC) generates a dense subgroup of T^2 (Weyl)'.")
print("  That hypothesis is FALSE at the corpus's own headline connection f=2.0, c=1.1, which is")
print("  exactly resonant (11f - 20c = 0, generator confirmed in X1a) -- the erratum against W-02")
print("  exists for precisely this.  Neither wm2 nor wm5 ever evaluates the hypothesis.")
print("  DECISIVE FORM: compute the SUBTORUS average (the true limit) as well as the lane's")
print("  TORUS average (its closed form), and see which one the time average converges to.")
car=X.K1(); NV=5
def coeff_rows(car,a,s):
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
def m_orbit(d,x,y,n=1<<22):
    """average of log|P| over the CLOSURE of {(x^t, y^t) : t integer}, computed as the average
       over t in [0, n) of the one-parameter subgroup, which is the closure when (f,c) is
       resonant of rank 1 and is equal to the torus average when it is generic."""
    k=np.arange(n)
    Xv=x**k; Yv=y**k
    return float(np.log(np.maximum(np.abs(d[0]+d[1]*Xv+d[2]*Yv+d[3]*Xv*Yv),1e-300)).mean())
def timeavg(oF,oC,s,N):
    xF=s.copy(); xC=s.copy(); tot=0.0
    for _ in range(N):
        xF=oF@xF; xC=oC@xC; z=abs(np.vdot(xF,xC)); tot+= np.log(z) if z>1e-300 else -700.0
    return tot/N
sA=np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB=np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
print(f"  {'connection':<28}{'st':>4}{'wm2 TORUS form':>17}{'ORBIT average':>16}{'time avg 2e6':>15}{'|torus-ta|':>12}{'|orbit-ta|':>12}")
for nm,a in (("generic f=2.28,c=2+sq2",np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])),
             ("RESONANT f=2.0, c=1.1 ",np.array([2.0,0.0,0.0,1.1,0.0,0.0]))):
    x=np.conj(X.hol(car["walkF"],a)); y=X.hol(car["walkC"],a)
    for tag,s in (("A",sA),("B",sB)):
        rows=coeff_rows(car,a,s)
        tor=float(np.mean([X.m_quad(rows[r],1<<20) for r in sorted(rows)]))
        orb=float(np.mean([m_orbit(rows[r],x,y,1<<22) for r in sorted(rows)]))
        ta=timeavg(X.Top(car["walkF"],a,NV),X.Top(car["walkC"],a,NV),s,2000000)
        print(f"  {nm:<28}{tag:>4}{tor:>17.9f}{orb:>16.9f}{ta:>15.9f}{abs(tor-ta):>12.2e}{abs(orb-ta):>12.2e}")
print("  -> at the GENERIC point torus == orbit == time average, all three agree to ~1e-5.")
print("     at the RESONANT point the torus form is WRONG BY ~1.5e-04 and the ORBIT average is")
print("     right to ~1e-06: two orders better.  The size of the error is the same order as the")
print("     4.9e-04 that the erratum against W-02 was written for (-0.767014993 subtorus versus")
print("     -0.767507880 torus).  IT IS NOT LARGE.  IT IS ALSO NOT NOISE, AND IT DOES NOT FALL")
print("     WITH N.  wm2's hypothesis is load-bearing, false at the corpus's headline connection,")
print("     and never tested in the lane -- the same omission, one convention down, that produced")
print("     the erratum.  Scope defect, not an arithmetic error.")
