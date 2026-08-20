"""VERIFY 6: does the lane's C-46 conclusion (criterion (d) FAILED, screening) survive on a
BOUNDARY-FREE instrument?  The lane measured on an OPEN chain, where V5 shows the odd and even
sublattice branches are contaminated by the boundary Friedel modulation by up to a factor 7.8.
Here: the same ratio on a RING, where there is no boundary, plus the analytic expectation
ln2/ln(R) that follows from J(r) = (-1)^{r+1}/(pi r)."""
import numpy as np
OUT=[]
def P(*a):
    s=" ".join(str(x) for x in a); print(s); OUT.append(s)
def T_pbc(N, epsF=0.0, t=1.0):
    n=np.arange(N); k=2*np.pi*n/N; eps=-2.0*t*np.cos(k)
    occ=np.where(eps<epsF-1e-12)[0]; emp=np.where(eps>epsF+1e-12)[0]
    G=np.zeros(N)
    for p in occ:
        G+=np.bincount((p-emp)%N, weights=1.0/(eps[emp]-eps[p]), minlength=N)
    return np.real(np.fft.fft(G))/(N*N)
P("="*112)
P("V6  C-46 RATIO |sum J|/sum|J| ON A BOUNDARY-FREE RING, AND THE ANALYTIC VALUE")
P("="*112)
P("")
P(f"{'N':>8} {'R':>6} {'sum J':>13} {'sum |J|':>13} {'|sum|/sum|.|':>14} "
  f"{'analytic ln2/lnR':>18} {'LANE (open m=8192)':>20}")
lane={8:0.203847087,32:0.148083504,128:0.122186423,512:0.130662126}
for N in (8194,32770):
    J=-8.0*T_pbc(N)
    for R in (8,32,128,512,2048):
        if R>N//8: continue
        v=J[1:R+1]
        P(f"{N:>8} {R:>6} {v.sum():>13.9f} {np.abs(v).sum():>13.9f} "
          f"{abs(v.sum())/np.abs(v).sum():>14.9f} {np.log(2)/np.log(R):>18.9f} "
          f"{(f'{lane[R]:.9f}' if R in lane and N==8194 else '-'):>20}")
P("")
P("D-15 CONTROL: a deliberately sign-definite A/r^3 through the identical statistic.")
for R in (8,32,128,512):
    v=np.array([0.30*r**-3.0 for r in range(1,R+1)])
    P(f"{'-':>8} {R:>6} {v.sum():>13.9f} {np.abs(v).sum():>13.9f} "
      f"{abs(v.sum())/np.abs(v).sum():>14.9f} {'-':>18} {'-':>20}")
P("")
P("D-15 CONTROL: sum|J| growth on the ring -- must be logarithmic and unbounded if J ~ 1/(pi r).")
J=-8.0*T_pbc(32770)
for R in (8,32,128,512,2048):
    P(f"   R={R:>5}  sum|J| = {np.abs(J[1:R+1]).sum():>10.6f}   "
      f"(1/pi)*(ln R + gamma) = {(np.log(R)+0.5772157)/np.pi:>10.6f}")
open("v6_c46_clean.txt","w").write("\n".join(OUT)+"\n")
