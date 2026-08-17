# S_09 — THE SCHEDULE HALF, CLOSED QUANTITATIVELY.  M2 self-flagged that its schedule finding has
# no theorem.  Here it does.  DEFINE:  J_max(K) = the largest number of cells an adversary can
# write out of k<=K while keeping the TOTAL record SUM(1-|Z_k|) <= 1 nat (i.e. |Omega| >= 1/e).
# ONE VARIABLE MOVES per pair of rows: d_eff (the dimension of the orbit closure H).  Ready state
# RS-G, budget B=1, evaluator and code path identical in every row.  Double precision.
import numpy as np, math
w=(0.4,0.3,0.3)
def gaps(f,c,K,chunk=10**6):
    out=[]
    for lo in range(1,K+1,chunk):
        hi=min(lo+chunk-1,K); k=np.arange(lo,hi+1,dtype=np.float64)
        z=np.abs(w[0]*np.exp(1j*k*(c-f))+w[1]*np.exp(-1j*k*f)+w[2]*np.exp(1j*k*c))
        out.append(1.0-np.minimum(z,1.0))
    return np.concatenate(out)
def Jmax(g,B=1.0):
    s=np.sort(g); cs=np.cumsum(s); return int(np.searchsorted(cs,B))
phi=(1+5**0.5)/2
ROWS=[("d_eff=2  DIOPH  -2^(1/3),4^(1/3) ",-2*np.pi*2**(1/3),2*np.pi*4**(1/3),2),
      ("d_eff=2  RANDOM seed 20260816    ",None,None,2),
      ("d_eff=1  RESONANT f=2.0,c=1.1    ",2.0,1.1,1),
      ("d_eff=1  W07GEN 2pi.phi,2pi.phi^2",2*np.pi*phi,2*np.pi*phi**2,1)]
rng=np.random.default_rng(20260816); a,b=rng.uniform(0,1,2)
ROWS[1]=(ROWS[1][0],2*np.pi*a,2*np.pi*b,2)
Ks=(10**4,10**5,10**6,10**7)
print("PREDICTION, derived from M2-4's own N(eps) law and nothing else:")
print("  N(eps) ~ C_d K eps^{d/2}  =>  the J-th smallest gap over k<=K is eps_J ~ (J/(C_d K))^{2/d}")
print("  =>  SUM of the J smallest ~ (d/(d+2)) J (J/(C_d K))^{2/d}  =>  at fixed budget")
print("        d_eff = 2 :  J_max ~ K^{1/2},  write density ~ K^{-1/2}")
print("        d_eff = 1 :  J_max ~ K^{2/3},  write density ~ K^{-1/3}")
print("  i.e. a RANK-1 connection is MORE exposed to the adversary than a rank-2 one, not less.")
print()
print(f"   {'connection':<35} " + " ".join(f"{'J(K=%.0e)'%K:>12}" for K in Ks) + f" {'fitted exp':>11} {'theory':>7}")
for tag,f,c,d in ROWS:
    Js=[]
    for K in Ks:
        Js.append(Jmax(gaps(f,c,K)))
    lg=np.polyfit(np.log([float(K) for K in Ks]), np.log([max(j,1) for j in Js]),1)[0]
    print(f"   {tag:<35} " + " ".join(f"{j:>12d}" for j in Js) + f" {lg:>11.3f} {2/(d+2)*2 if False else (0.5 if d==2 else 2/3):>7.3f}")
print()
print("   And the ATTAINED case, for contrast: on S1's order-4 connection 1-|Z_k| = 0 EXACTLY on")
g=gaps(np.pi,3*np.pi/2,10**6)
print(f"   a fraction 1/4 of all cells -- {int((g<=1e-15).sum())} of 10^6 -- so J_max = infinity at")
print("   ZERO budget.  THAT, and only that, is what 'ATTAINED' buys the adversary: a FIXED")
print("   positive density instead of a density that must fall like K^{-1/2} or K^{-1/3}.")
