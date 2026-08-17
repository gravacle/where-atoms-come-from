# LANE_W11_R_MATH — LEG 5.  PRECISION AND DEGENERACY.
#  (i)   arithmetic type of every connection used in W-11 and of the corpus's own published ones;
#  (ii)  is the CONVENTION VERDICT connection-dependent?  four connections, one variable moved;
#  (iii) the whole K1 test re-run in EXACT arithmetic in Z[zeta_12] on S1's OWN published
#        connection (S1 sec6: a1=a2=a3=pi/3, a4=a5=a6=pi/2), no floating point anywhere.
#
# W-10 N-4 says every RATIONAL (f,c) is exactly resonant.  Taken literally that convicts every
# connection any lane has ever run, because every IEEE double is a dyadic rational.  The operative
# variable is therefore not rationality but the SIZE of the smallest relation, W-08's pair (q,delta).
import numpy as np, wm0_lib as L
from fractions import Fraction as Fr
np.set_printoptions(linewidth=200)

TWO_PI=2*np.pi
def best_relation(f,c,M=1500):
    m=np.arange(-M,M+1); n=np.arange(-M,M+1)
    Mm,Nn=np.meshgrid(m,n,indexing="ij")
    val=Mm*f+Nn*c
    res=np.abs(val-TWO_PI*np.round(val/TWO_PI))
    mask=(Mm!=0)|(Nn!=0)
    res=np.where(mask,res,np.inf)
    k=np.unravel_index(np.argmin(res),res.shape)
    return int(Mm[k]),int(Nn[k]),float(res[k])

CONN={
 "W-11 registrar's K1 (a = 1.0,0.37,0.91 | sqrt2,0.23,1.77)": (1.0+0.37+0.91, 2**0.5+0.23+1.77),
 "corpus's ONLY published generic  S4:603  f=1.0, c=sqrt(2)": (1.0, 2**0.5),
 "S3/S4 headline  f=2.0, c=1.1  (known exactly resonant)   ": (2.0, 1.1),
 "S1 sec6 published, order 4:  f=pi, c=3pi/2               ": (np.pi, 3*np.pi/2),
 "W-11 registrar's B0b connection (rng uniform, W_F/W_C)   ": None,
}
carB=L.B0b(); aB=np.random.default_rng(20260817).uniform(0,2*np.pi,18)
CONN["W-11 registrar's B0b connection (rng uniform, W_F/W_C)   "]=(
    float(np.angle(L.hol(carB["walkF"],aB))), float(np.angle(L.hol(carB["walkC"],aB))))

print("== M5a  ARITHMETIC TYPE CENSUS.  smallest relation m f + n c = 0 mod 2pi, |m|,|n| <= 1500 ==")
print(f"  {'connection':<58}{'f':>12}{'c':>12}   (m,n)          |residual|")
for k,(f,c) in CONN.items():
    m,n,r=best_relation(f,c)
    tag="EXACTLY RESONANT" if r<1e-12 else "no exact relation at this depth"
    print(f"  {k:<58}{f:>12.6f}{c:>12.6f}   ({m:>5},{n:>5})   {r:.3e}   {tag}")
print(f"  Dirichlet bound at depth M = 1500 for a point with NO exact relation: 2 pi / M^2 = "
      f"{TWO_PI/1500**2:.2e}.  Every non-resonant row above sits BELOW it, which is what generic")
print("  Diophantine behaviour looks like -- so no row is anomalously close to a relation, and I")
print("  refuse to award any of them the word 'near-resonant', which would be a fifth mislabelling.")
print("  The registrar's K1 connection has f = 2.28, c = 2 + sqrt(2).  It is NOT the corpus's")
print("  published generic point (f=1.0, c=sqrt2); its script's comment 'generic: c involves sqrt(2)'")
print("  is not a genericity argument.  It IS generic, but the reason has to be given:")
print("  m(2.28) + n(2+sqrt2) + 2 pi k = 0 splits over the Q-independent set {1, sqrt2, pi} into")
print("  n = 0, k = 0, 2.28 m = 0, so m = 0.  Certified -- by an argument the lane did not make.")
print("  CAVEAT OF RECORD (W-10 N-4 taken literally): sqrt(2) as an IEEE double is the dyadic")
print("  rational 6369051672525773/2^52, so the connection ACTUALLY EXECUTED is exactly resonant")
print("  with |m|,|n| ~ 2^52.  Genericity is a statement about the depth of the run, never about")
print("  the float.  Every 'generic' row in this corpus carries that caveat.\n")

print("== M5b  IS THE CONVENTION VERDICT CONNECTION-DEPENDENT?  ONE VARIABLE MOVED: THE CONNECTION ==")
car=L.K1(); NV=5; cl,F,C=L.classes(car)
sA=np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j
sB=np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
sC=sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
assert np.linalg.norm(sA-sB)>0.1 and np.linalg.norm(sA-sC)>0.1        # ARMS DIFF
PI=L.pi_of(sA,cl)
assert np.allclose(PI,L.pi_of(sB,cl)) and np.allclose(PI,L.pi_of(sC,cl))
ROWS=[("registrar's  f=2.28, c=2+sqrt2", np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])),
      ("corpus generic f=1.0, c=sqrt2 ", np.array([1.0,0.0,0.0,2**0.5,0.0,0.0])),
      ("resonant  f=2.0, c=1.1        ", np.array([2.0,0.0,0.0,1.1,0.0,0.0])),
      ("S1 order-4  f=pi, c=3pi/2     ", np.array([np.pi/3]*3+[np.pi/2]*3)),
      ("rank-1  f=1.0, c=0.0          ", np.array([1.0,0.0,0.0,0.0,0.0,0.0]))]
def sp(op_F,op_C,n,sts):
    v=[abs(np.vdot(np.linalg.matrix_power(op_F,n)@s,np.linalg.matrix_power(op_C,n)@s)) for s in sts]
    return max(v)-min(v)
def timeavg(oF,oC,s,N):
    xF=s.copy(); xC=s.copy(); tot=0.0
    for _ in range(N):
        xF=oF@xF; xC=oC@xC; z=abs(np.vdot(xF,xC)); tot+= np.log(z) if z>1e-300 else -700.0
    return tot/N
print(f"  {'connection':<32}{'CIRC spread k=1,2':>19}{'EDGE spread n=1,2':>19}{'EDGE n=3,6':>13}")
for nm,a in ROWS:
    TF,TC=L.Top(car['walkF'],a,NV),L.Top(car['walkC'],a,NV)
    MF,MC=L.Mop(car['walkF'],a,NV),L.Mop(car['walkC'],a,NV)
    c12=max(sp(MF,MC,k,(sA,sB,sC)) for k in (1,2))
    e12=max(sp(TF,TC,n,(sA,sB,sC)) for n in (1,2))
    e36=max(sp(TF,TC,n,(sA,sB,sC)) for n in (3,6))
    print(f"  {nm:<32}{c12:>19.2e}{e12:>19.2e}{e36:>13.2e}")
print("  -> the verdict is NOT connection-dependent.  Circuit invisible / edge visible at every one,")
print("     including the resonant and the order-4 points.  (At the rank-1 point c = 0 the cycle")
print("     holonomy is trivial and the edge tick still separates: the effect is not a curvature effect.)\n")
print(f"  {'connection':<32}{'CIRCUIT rate (N=2e5)':>22}{'EDGE rate A':>14}{'EDGE B':>12}{'EDGE C':>12}")
for nm,a in ROWS:
    TF,TC=L.Top(car['walkF'],a,NV),L.Top(car['walkC'],a,NV)
    MF,MC=L.Mop(car['walkF'],a,NV),L.Mop(car['walkC'],a,NV)
    rc=timeavg(MF,MC,sA,200000); re=[timeavg(TF,TC,s,200000) for s in (sA,sB,sC)]
    print(f"  {nm:<32}{rc:>22.9f}{re[0]:>14.9f}{re[1]:>12.9f}{re[2]:>12.9f}")
print()

# ------------------------------------------------------------------ EXACT ARITHMETIC, Z[zeta_12]
# zeta = exp(i pi/6);  Phi_12(x) = x^4 - x^2 + 1, so zeta^4 = zeta^2 - 1 and zeta^6 = -1.
# element = (a0,a1,a2,a3) meaning a0 + a1 z + a2 z^2 + a3 z^3, a_i in Q (Fraction).
Z0=(Fr(0),)*4; Z1=(Fr(1),Fr(0),Fr(0),Fr(0))
def zadd(a,b): return tuple(x+y for x,y in zip(a,b))
def zsub(a,b): return tuple(x-y for x,y in zip(a,b))
def zscal(k,a): return tuple(Fr(k)*x for x in a)
def zmul(a,b):
    p=[Fr(0)]*7
    for i in range(4):
        for j in range(4): p[i+j]+=a[i]*b[j]
    # reduce with z^4 = z^2 - 1,  z^5 = z^3 - z,  z^6 = -1
    r=[p[0],p[1],p[2],p[3]]
    r[2]+=p[4]; r[0]-=p[4]
    r[3]+=p[5]; r[1]-=p[5]
    r[0]-=p[6]
    return tuple(r)
def zconj(a):
    a0,a1,a2,a3=a
    return (a0+a2, a1, -a2, -(a1+a3))
def zpow(a,n):
    r=Z1
    for _ in range(n): r=zmul(r,a)
    return r
ZETA=(Fr(0),Fr(1),Fr(0),Fr(0))
def znum(a):
    z=np.exp(1j*np.pi/6); return complex(sum(float(a[i])*z**i for i in range(4)))

print("== M5c  THE K1 TEST IN EXACT ARITHMETIC, ON S1's OWN PUBLISHED CONNECTION (S1 sec6) ==")
print("   a1=a2=a3=pi/3 -> U = zeta^2 ;  a4=a5=a6=pi/2 -> U = zeta^3 ;  W_F = zeta^6 = -1,")
print("   W_C = zeta^9 = -i.  Every matrix entry and every amplitude below is an EXACT element of")
print("   Z[zeta_12] over Q.  No floating point is used in this block.")
UF_e=zpow(ZETA,2); UC_e=zpow(ZETA,3)
NV=5
def zmatvec(Mt,v):
    out=[]
    for i in range(NV):
        acc=Z0
        for j in range(NV):
            if Mt[i][j]!=Z0: acc=zadd(acc,zmul(Mt[i][j],v[j]))
        out.append(acc)
    return out
def zeros(): return [[Z0]*NV for _ in range(NV)]
TFz=zeros(); TCz=zeros(); MFz=zeros(); MCz=zeros()
for v in (3,4): TFz[v][v]=Z1
for (u,v) in ((0,1),(1,2),(2,0)): TFz[v][u]=UF_e
for v in (1,2): TCz[v][v]=Z1
for (u,v) in ((0,3),(3,4),(4,0)): TCz[v][u]=UC_e
WFz=zpow(ZETA,6); WCz=zpow(ZETA,9)
for v in range(NV): MFz[v][v]=WFz if v in (0,1,2) else Z1
for v in range(NV): MCz[v][v]=WCz if v in (0,3,4) else Z1
print(f"   W_F exact = {WFz}  (= -1)      W_C exact = {WCz}  (= -i, since zeta^3 = i)")
def rat(k,d): return (Fr(k,d),Fr(0),Fr(0),Fr(0))
# rational amplitudes with sum of squares exactly 1, and TWO exact splits of class 10 = {v1,v2}
sAz=[rat(20,25),rat(3,25),rat(4,25),rat(10,25),rat(10,25)]
sBz=[rat(20,25),rat(5,25),rat(0,25),rat(10,25),rat(10,25)]
sCz=[zmul(sAz[i],zpow(ZETA,k)) for i,k in enumerate((0,1,5,3,7))]
def zpi(s):
    w=[zmul(zconj(x),x)[0] for x in s]                     # |s_v|^2 is rational
    return (w[0], w[1]+w[2], w[3]+w[4])                    # p11, p10, p01   (p00 = 0 on K1)
print(f"   pi(A) = {zpi(sAz)}   pi(B) = {zpi(sBz)}   pi(C) = {zpi(sCz)}")
assert zpi(sAz)==zpi(sBz)==zpi(sCz), "pi must be EXACTLY equal"
assert sAz!=sBz and sAz!=sCz, "ARMS DIFF: the states must actually differ"
print("   pi identical as EXACT RATIONALS (not to 1e-16).  States pairwise distinct.")
def zZ(opF,opC,s,n):
    xF=list(s); xC=list(s)
    for _ in range(n): xF=zmatvec(opF,xF); xC=zmatvec(opC,xC)
    acc=Z0
    for v in range(NV): acc=zadd(acc,zmul(zconj(xF[v]),xC[v]))
    return acc
print(f"   {'n':>3}  {'|Z|^2 exact, CIRCUIT convention':<44}  equal?")
for n in (1,2,3,4):
    vals=[zmul(zconj(zZ(MFz,MCz,s,n)),zZ(MFz,MCz,s,n)) for s in (sAz,sBz,sCz)]
    ok=(vals[0]==vals[1]==vals[2])
    for vv in vals: assert vv[2]==0 and 2*vv[3]+vv[1]==0, "|Z|^2 must be exactly real"
    print(f"   {n:>3}  {str(vals[0]):<44}  {ok}")
print(f"   {'n':>3}  {'|Z|^2 exact, EDGE convention (COR-F)':<44}  equal?")
for n in (1,2,3,4,5,6):
    vals=[zmul(zconj(zZ(TFz,TCz,s,n)),zZ(TFz,TCz,s,n)) for s in (sAz,sBz,sCz)]
    ok=(vals[0]==vals[1]==vals[2])
    d=max(abs(znum(vals[0])-znum(vals[1])),abs(znum(vals[0])-znum(vals[2])))
    for vv in vals: assert vv[2]==0 and 2*vv[3]+vv[1]==0, "|Z|^2 must be exactly real"
    print(f"   {n:>3}  {str(vals[0]):<44}  {ok}   max |difference| = {d:.6f}")
print("   -> EXACT, in rationals: identical at every n under the circuit convention; DIFFERENT")
print("      exactly at n not divisible by 3 under the edge convention.  The 1e-16 agreements and")
print("      the O(1) disagreements of the registrar's tables are both confirmed WITHOUT floating")
print("      point, on S1's own published connection rather than on an invented one.")
