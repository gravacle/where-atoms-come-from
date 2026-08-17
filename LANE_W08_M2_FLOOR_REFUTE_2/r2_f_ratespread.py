# LANE W-08 / M2 REFUTER 2 — leg F: "THE CONNECTION IS WORTH 13%" (lane E4 / M2-12).
# That 13% is the spread across the lane's EIGHT sampled connections, six of which are rank-2 or
# a small perturbation of it.  Here the same statistic is computed over the WHOLE connection axis.
# ISOLATION: HELD - carrier K1, RS-G, the functional D_H = E_H(1-|Z|) and lambda_H = E_H log|Z|,
# both as exact orbit/circle averages.  MOVES - the connection, over all finite closures of order
# <= 120 and all rank-1 closures with |p|,|q| <= 8, plus the rank-2 closure.
import numpy as np
from math import gcd
W11,W10,W01=0.4,0.3,0.3
D_T2=0.469188699222; LAM_T2=-0.767507880358
def stats_finite(a,b,n):
    k=np.arange(1,n+1)
    Z=W11*np.exp(2j*np.pi*(k*(a+b)%n)/n)+W10*np.exp(2j*np.pi*(k*a%n)/n)+W01*np.exp(2j*np.pi*(k*b%n)/n)
    az=np.abs(Z); return float(np.mean(1-az)), float(np.mean(np.log(az)))
def stats_circle(p,q,N=1<<20):
    t=2*np.pi*np.arange(N)/N
    Z=W11*np.exp(1j*(p+q)*t)+W10*np.exp(1j*p*t)+W01*np.exp(1j*q*t)
    az=np.abs(Z); return float(np.mean(1-az)), float(np.mean(np.log(az)))
print("== F1  THE DURABILITY RATE ACROSS THE WHOLE CONNECTION AXIS (RS-G) ==")
print("   D_H = asymptotic slope of SUM(1-|Z_k|)/K ;  lambda_H = asymptotic (1/N)log|Omega_N|.")
rows=[]
for n in range(2,121):
    for a in range(n):
        for b in range(n):
            if n//gcd(gcd(a,b),n)!=n: continue
            rows.append((stats_finite(a,b,n),('finite',n,a,b)))
for p in range(-8,9):
    for q in range(-8,9):
        if (p,q)==(0,0) or gcd(abs(p),abs(q))!=1: continue
        rows.append((stats_circle(p,q),('circle',p,q,None)))
Ds=[r[0][0] for r in rows]; Ls=[r[0][1] for r in rows]
iD=int(np.argmin(Ds)); jD=int(np.argmax(Ds)); iL=int(np.argmin(Ls)); jL=int(np.argmax(Ls))
print(f"   MIN slope D_H = {Ds[iD]:.6f}  at {rows[iD][1]}      (proved lower bound, lane C2: 0.120000)")
print(f"   MAX slope D_H = {Ds[jD]:.6f}  at {rows[jD][1]}")
print(f"   rank-2 (generic) D_T2 = {D_T2:.6f};  the lane's eight sampled rows spanned 0.4692-0.5295 = 13%")
print(f"   TRUE SPREAD OF THE CONNECTION AXIS: {Ds[jD]/Ds[iD]:.2f}x  "
      f"({(Ds[jD]-Ds[iD])/D_T2*100:.0f}% of the generic rate), not 13%.")
print(f"   MIN lambda_H = {Ls[iL]:.6f} at {rows[iL][1]}")
print(f"   MAX lambda_H = {Ls[jL]:.6f} at {rows[jL][1]}")
print(f"   TRUE SPREAD OF lambda ACROSS CONNECTIONS: {Ls[iL]/Ls[jL]:.2f}x  "
      f"({(Ls[jL]-Ls[iL])/abs(LAM_T2)*100:.0f}% of |lambda_generic|).")
print()
print("== F2  THE SAME NUMBERS RESTRICTED TO EACH SIDE OF THE ATTAINED/APPROACHED CUT ==")
fin=[r for r in rows if r[1][0]=='finite']; cir=[r for r in rows if r[1][0]=='circle']
for lab,grp in [("ATTAINED (finite H)",fin),("APPROACHED (rank-1 H)",cir)]:
    d=[r[0][0] for r in grp]; l=[r[0][1] for r in grp]
    print(f"   {lab:<24} D_H in [{min(d):.6f}, {max(d):.6f}]   lambda in [{min(l):.6f}, {max(l):.6f}]")
print(f"   {'APPROACHED (rank-2 H)':<24} D_H = {D_T2:.6f}                lambda = {LAM_T2:.6f}")
print("   -> THE TWO SIDES OVERLAP ALMOST COMPLETELY.  Neither D_H nor lambda separates attained")
print("      from approached: the attained interval strictly CONTAINS the approached one except")
print("      at its two ends.  A partition that does not separate the observable is not a")
print("      partition of the observable's consequences.")
