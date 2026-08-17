# REFUTER 1 of LANE W-08 / M2 — LENS: MATHEMATICS.  Leg R4: the RATES, by routes M2 did not use.
#   R4.1  lambda on T^2 by JENSEN on the inner variable (an exact reduction to a 1-d integral with
#         the kink split by hand) — independent of both of M2's routes (its own quadrature and
#         Cassaigne-Maillot).  Target: -0.767507880358 (RS-G), -0.693147180560 (RS-P).
#   R4.2  the subtorus and circle rates, by direct 1-d quadrature instead of np.roots+Jensen.
#         Targets: -0.767014992998 (B_S3RES) and -1.203972804326 (E_W07GEN, claimed = log 0.3).
#   R4.3  the |H| = 4 rate in EXACT rational arithmetic.  Target: -0.804718956217.
#   R4.4  D_H.  M2-9 states "D_H = 1-2/pi = 0.363380227639".  1-2/pi is checked to 30 digits.
# PRECISION: mpmath at 50 dps throughout; the |H|=4 row in exact Fractions.
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps=50

print("== R4.1  lambda ON T^2 BY JENSEN — A THIRD, INDEPENDENT ROUTE ==")
print("   m(p11 XY + p10 X + p01 Y): fix X = e^{2 pi i phi}; the Y-integral is Jensen's formula for")
print("   the linear polynomial (p11 X) Y + p10 X, giving log max(|p11 X + p10|, p01) -- wait: the")
print("   Y-coefficient is p11 X + p01 ... so group as  Y*(p11 X) + (p10 X + 0) is WRONG.  Group")
print("   correctly:  Z = X*(p11 Y + p10) + p01 Y  is not linear in one variable alone either.")
print("   USE INSTEAD:  Z = Y*(p11 X + p01) + p10 X.  Jensen in Y:")
print("       int log|Z| dY = log max( |p11 X + p01| , |p10 X| ) = log max(|p11 X + p01|, p10).")
for tag,(p11,p10,p01) in [("RS-G",(mp.mpf('0.4'),mp.mpf('0.3'),mp.mpf('0.3'))),
                          ("RS-P",(mp.mpf('0.5'),mp.mpf(0),mp.mpf('0.5')))]:
    def inner(ph):
        A=abs(p11*mp.expjpi(2*ph)+p01)
        return mp.log(A if A>p10 else p10)
    # kink where |p11 X + p01| = p10  ->  p11^2+p01^2+2 p11 p01 cos = p10^2
    ks=[]
    if p10>0:
        co=(p10**2-p11**2-p01**2)/(2*p11*p01)
        if abs(co)<=1:
            t=mp.acos(co)/(2*mp.pi); ks=sorted({t,1-t})
    # log-singularities of the integrand (|p11 X + p01| = 0) must be split too, or mp.quad
    # evaluates AT them and returns -inf.  MY defect on the first run, recorded not hidden.
    sing=[mp.mpf('0.5')] if abs(p11-p01)<mp.mpf('1e-30') else []
    val=mp.quad(inner,[0]+sorted(ks+sing)+[1])
    print(f"   {tag}: lambda = {mp.nstr(val,15)}   kinks at phi = {[mp.nstr(k,8) for k in ks]}"
          f"   log-singularities split at {[mp.nstr(s,4) for s in sing]}")
print(f"   M2 reports  RS-G -0.767507880358   RS-P -0.693147180560  (= -log 2: {mp.nstr(-mp.log(2),15)})")
print()

print("== R4.2  THE dim-1 RATES BY DIRECT QUADRATURE (M2 used np.roots + Jensen) ==")
for tag,(ex,ey),(p11,p10,p01),claim in [
    ("B_S3RES 11a-20b=0",(-20,11),(mp.mpf('0.4'),mp.mpf('0.3'),mp.mpf('0.3')),'-0.767014992998'),
    ("E_W07GEN a-b=0    ",(-1,1), (mp.mpf('0.4'),mp.mpf('0.3'),mp.mpf('0.3')),'-1.203972804326'),
    ("B_S3RES RS-P      ",(-20,11),(mp.mpf('0.5'),mp.mpf(0),mp.mpf('0.5')),'-0.693147180560')]:
    def f(s):
        z=lambda e: mp.expjpi(2*e*s)
        return mp.log(abs(p11*z(ex+ey)+p10*z(ex)+p01*z(ey)))
    pts=list(mp.linspace(0,1,4*(abs(ex)+abs(ey))+1))
    # SPLIT AT THE LOG-SINGULARITIES OF THE INTEGRAND, or the quadrature is silently wrong.
    # E_W07GEN: Z(s) = 0.4 + 0.6 cos(2 pi s) vanishes at s = +-acos(-2/3)/2pi.
    # RS-P rows: Z(s) = 0.5(1 + z^{n}) vanishes at the n odd half-integers.
    if p10==0:
        n=abs(ey-ex); pts+= [mp.mpf(2*j+1)/(2*n) for j in range(n)]
    if (ex,ey)==(-1,1) and p10>0:
        t=mp.acos(-(p11)/(2*p01))/(2*mp.pi); pts+=[t,1-t]
    v=mp.quad(f,sorted(set(pts)))
    print(f"   {tag}  lambda(quadrature) = {mp.nstr(v,13):<18} M2: {claim}")
print("   E_W07GEN's polynomial is 0.3 + 0.4 z + 0.3 z^2 (after the monomial shift).  Its roots:")
r=mp.polyroots([mp.mpf('0.3'),mp.mpf('0.4'),mp.mpf('0.3')])
print(f"     roots {[mp.nstr(x,10) for x in r]}   |roots| = {[mp.nstr(abs(x),20) for x in r]}")
print(f"     both EXACTLY on the unit circle (self-reciprocal, discriminant 0.16-0.36 < 0), so")
print(f"     Jensen gives lambda = log(0.3) = {mp.nstr(mp.log(mp.mpf('0.3')),15)} EXACTLY.  M2's")
print(f"     -1.203972804326 is right, and 'exactly log 0.3' is right.")
print()

print("== R4.3  THE |H| = 4 RATE, EXACT RATIONALS ==")
p11,p10,p01=Fr(2,5),Fr(3,10),Fr(3,10)
# x = conj(W_F) = -1, y = W_C = -i  (S1 sec6: f = pi, c = 3pi/2)
xs=[(-1,0),(1,0),(-1,0),(1,0)]; ys=[(0,-1),(-1,0),(0,1),(1,0)]
z2=[]
for (xr,xi),(yr,yi) in zip(xs,ys):
    pr,pj=xr*yr-xi*yi, xr*yi+xi*yr
    re=p11*pr+p10*xr+p01*yr; im=p11*pj+p10*xi+p01*yi
    z2.append(re*re+im*im)
print(f"   |Z_k|^2, k=1..4 (exact) = {[str(t) for t in z2]}   product = {z2[0]*z2[1]*z2[2]*z2[3]}")
prod=z2[0]*z2[1]*z2[2]*z2[3]
lam=mp.log(mp.mpf(prod.numerator)/prod.denominator)/8
print(f"   lambda = (1/8) log(prod) = {mp.nstr(lam,15)}   M2: -0.804718956217")
print(f"   and (1/4) log(1/25) = {mp.nstr(mp.log(mp.mpf(1)/25)/4,15)}  (M2 writes (1/4)log(0.04))")
print()

print("== R4.4  D_H ON RS-P: M2-9 SAYS 'D_H = 1-2/pi = 0.363380227639' ==")
print(f"   1 - 2/pi to 30 digits = {mp.nstr(1-2/mp.pi,30)}")
print(f"   i.e. 1-2/pi = 0.363380227632...  M2's C1 prints 0.363380227639 for the T^2 row and")
print(f"   0.363380227633 / 0.363380227632 for its two circle rows.  The T^2 row's last two")
print(f"   digits are wrong (its mp.quad of the elliptic-integral form is not split at m=1);")
print(f"   the value quoted in finding M2-9 as an EQUALITY with 1-2/pi is that wrong one.")
print(f"   Size of the error: {mp.nstr(mp.mpf('0.363380227639')-(1-2/mp.pi),4)}.  Immaterial to every")
print(f"   claim, but this corpus has two corrections of record (COR-E, COR-K) for exactly this")
print(f"   class of defect: a computed figure printed as an identity to digits it does not have.")
# reproduce M2's own route to show where it loses the digits
def E_abs_torus(w11,w10,w01):
    def inner(ph):
        A=w11*mp.expjpi(2*ph)+w10; a=abs(A); b=mp.mpf(w01)
        if a+b==0: return mp.mpf(0)
        m=4*a*b/(a+b)**2
        return (2/mp.pi)*(a+b)*mp.ellipe(m)
    return mp.quad(inner,[0,1])
v=E_abs_torus(0.5,0.0,0.5)
print(f"   reproduction of M2's route on RS-P: 1 - E|Z| = {mp.nstr(1-v,15)}  (integrand is the")
print(f"   CONSTANT 2/pi there, so the whole error is mp.quad/ellipe at m = 1)")
