# Registrar's own verification of the three W-08 findings that correct W-07. Nothing taken on trust.
import numpy as np
from fractions import Fraction as F

print("== V1  WHERE WERE W-01's RECURRENCE FIGURES ACTUALLY MEASURED? ==")
# S2_AUDIT:439-440 and S3:423 both state the test point verbatim: f=2.0, c=1.1, p=(0.4,.15,.15,.15,.15)
# Z_k = P0 conj(W_F)^k W_C^k + PF conj(W_F)^k + PC W_C^k   (root / F-only / C-only)
p=[0.4,0.15,0.15,0.15,0.15]; P0,PF,PC=p[0],p[1]+p[2],p[3]+p[4]
def Zabs(aF,aC,K):
    k=np.arange(1,K+1)
    return np.abs(P0*np.exp(1j*k*(aC-aF)) + PF*np.exp(-1j*k*aF) + PC*np.exp(1j*k*aC))
z=Zabs(2.0,1.1,400)
print(f"  f=2.0 c=1.1 :  min |Z_k| over k<=400 = {z.min():.6f} at k = {int(z.argmin())+1}"
      f"   [register/S2-audit: 0.024654 at n=42]")
z4=Zabs(2.0,1.1,4000)
print(f"                 sup |Z_k| over k<=4000 = {z4.max():.6f} at k = {int(z4.argmax())+1}"
      f"   [register: 0.99994 / S3: 0.999941 at k=377]")
# and the arithmetic type of that connection
print(f"  -11f + 20c = {-11*2.0 + 20*1.1:.1f}   -> EXACTLY RESONANT (erratum against W-02)")
rho=np.exp(1j*(1.1-2.0))
ordr=next((n for n in range(1,10**7) if abs(rho**n-1)<1e-12), None)
print(f"  branch ratio order at f=2.0,c=1.1 : {ordr}  (None = infinite)")
# same figures at S1's published order-4 connection?
zp=Zabs(np.pi,3*np.pi/2,400)
print(f"  S1 PUBLISHED (order 4): |Z_k| takes {len(set(np.round(zp,9)))} distinct values: {sorted(set(np.round(zp,6)))}")
print(f"  distance from 0.024654 to nearest attained value there = {min(abs(v-0.024654) for v in set(np.round(zp,9))):.6f}")
print("  ==> W-01's recurrence figures were measured at the RESONANT connection, NOT at order 4.")
print("      MY W-07 HEADLINE IS FALSE OF THEM. It is true of W-06's 1000-of-4000 figure only.\n")

print("== V2  IS 3*sqrt(3)/10 REACHABLE ON S1's ORDER-4 CONNECTION?  (exact) ==")
# W-08 claims |s|^2 = (3/4, 4/25, 0, 9/100, 0), pair (v0,v3): (dF,dC)=(-1,0) so D = amp*|(-1)^k - 1|
w=[F(3,4),F(4,25),F(0),F(9,100),F(0)]
print(f"  sum of |s_v|^2 = {sum(w)}   (must be 1: {sum(w)==1})")
amp2=w[0]*w[3]                                  # |t_0||t_3| squared; dressing is unimodular
print(f"  amp^2 = |s_0|^2 |s_3|^2 = {amp2}      (2*amp)^2 = {4*amp2}")
print(f"  (3*sqrt(3)/10)^2 = 27/100 = {F(27,100)}   EQUAL: {4*amp2 == F(27,100)}")
print(f"  numerically: 2*amp = {float(2*(amp2**0.5)):.16f}   3*sqrt(3)/10 = {3*3**0.5/10:.16f}")
print("  ==> REACHABLE, exactly. The sqrt(3) is in the READY STATE's amplitudes, not the group.")
print("      MY W-07 REGISTER SENTENCE 'a factor sqrt3 needs an element of order 3' IS FALSE.")
print("      (The W-07 PAGE listed 'a different normalisation' as an escape; the ROW did not.)\n")

print("== V3  IS |Omega_N| MONOTONE?  i.e. does a near-return UN-WRITE earlier cells? ==")
for tag,aF,aC in [("resonant f=2.0,c=1.1",2.0,1.1),("S1 published (order 4)",np.pi,3*np.pi/2),
                  ("random",1.7231,4.9014)]:
    zz=Zabs(aF,aC,10**6)
    print(f"  {tag:<24} max_k(|Z_k| - 1) = {(zz-1).max():+.3e}   exceedances of 1: {int((zz>1+1e-12).sum())}")
print("  |Z_k| <= 1 always (triangle inequality on non-negative weights summing to 1),")
print("  so |Omega_N| = prod |Z_k| is MONOTONE NON-INCREASING. A cell that writes nothing")
print("  does not un-write the previous N-1.  The founding obstruction reads the SINGLE-CELL")
print("  observable |Z_k| and infers a property of the RECORD. That inference does not hold.")
