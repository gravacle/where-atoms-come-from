# S_03 — (a) the record can never recover: |Omega_N| is monotone.  (b) the ADVERSARIAL schedule,
# which is the only place the founding obstruction survives.  (c) the share of the decay budget
# carried by the near-return band, as a function of the threshold.
import numpy as np
RSG=(0.4,0.3,0.3)
def prof(f,c,w,K,chunk=10**6):
    out=[]
    for lo in range(1,K+1,chunk):
        hi=min(lo+chunk-1,K); k=np.arange(lo,hi+1,dtype=np.float64)
        z=np.abs(w[0]*np.exp(1j*k*(c-f))+w[1]*np.exp(-1j*k*f)+w[2]*np.exp(1j*k*c))
        out.append(1.0-np.minimum(z,1.0))
    return np.concatenate(out)

print("== S3a  |Omega_N| IS MONOTONE NON-INCREASING.  |Z_k| <= p11+p10+p01 = 1 by the triangle")
print("        inequality, so no circuit can ever RESTORE the record.  Checked, not assumed:")
for tag,f,c in [("RESONANT",2.0,1.1),("DIOPH",-2*np.pi*2**(1/3),2*np.pi*4**(1/3))]:
    g=prof(f,c,RSG,10**6)
    print(f"   {tag:<10}  max_k |Z_k| - 1 = {(-g).max():+.3e}   (must be <= 0)   #(|Z_k|>1) = {(g<0).sum()}")

print()
print("== S3b  THE HONEST SCHEDULE k_n = n, and the ADVERSARIAL schedule of M2-refuter-2 ==")
print("   Adversary writes only the delta_K*K cells of SMALLEST 1-|Z_k| with delta_K = K^{-1/2}.")
print(f"   {'K':>9} {'honest SUM(1-|Z|)':>19} {'#writes adv':>12} {'adversarial SUM':>16} {'|Omega| adv':>12}")
f,c = -2*np.pi*2**(1/3), 2*np.pi*4**(1/3)
for K in (10**4,10**5,10**6,10**7):
    g=prof(f,c,RSG,K); J=int(round(K**0.5))
    part=np.partition(g,J)[:J]
    print(f"   {K:>9} {g.sum():>19.4f} {J:>12} {part.sum():>16.4f} {np.exp(-part.sum()):>12.4f}")
print("   The adversary writes J = sqrt(K) -> infinity cells and the record stays bounded away")
print("   from 0.  Under the corpus's OWN registered criterion (W-02: divergence of SUM(1-z_n))")
print("   durability FAILS on an approached connection too.  Two-way, and said as such.")

print()
print("== S3c  SHARE OF THE DECAY BUDGET CARRIED BY THE NEAR-RETURN BAND (DIOPH, K=1e7) ==")
g=prof(f,c,RSG,10**7); tot=g.sum()
print(f"   total decay budget SUM(1-|Z_k|) = {tot:.4f}   density c = {tot/1e7:.6f}")
print(f"   {'threshold eps':>14} {'#cells':>9} {'#/K':>11} {'0.8388*eps':>11} {'band SUM':>10} {'band SHARE':>12}")
for eps in (1e-1,1e-2,1e-3,1e-4,1e-5,1e-6):
    m=g<eps
    print(f"   {eps:>14.0e} {int(m.sum()):>9d} {m.sum()/1e7:>11.3e} {0.838820*eps:>11.3e} "
          f"{g[m].sum():>10.4f} {g[m].sum()/tot:>12.3e}")
print("   #cells/K reproduces M2's UNFITTED prediction N(eps)/K = 0.838820*eps (Haar measure of")
print("   the quadratic sublevel set, det M = 0.036).  The band's SHARE of the budget is")
print("   ~ 0.8388*eps^2/2c, i.e. it FALLS quadratically in eps and is INDEPENDENT of K.")
