# Registrar's verification of the two W-13 findings that correct register rows.
import numpy as np, itertools
def lam(pi,f,c,N):
    """(1/N) sum log|Z_k|, Z_k = p00 + p10 u^k + p01 v^k + p11 (uv)^k, u=conj(W_F)=e^{-if}, v=e^{ic}"""
    k=np.arange(1,N+1); u=np.exp(-1j*f*k); v=np.exp(1j*c*k)
    Z=pi[0]+pi[1]*u+pi[2]*v+pi[3]*u*v
    return float(np.log(np.abs(Z)+1e-300).mean())
def m_jensen(p,n=1<<20):
    a,b,c_,d=p; t=2*np.pi*np.arange(n)/n; ct=np.cos(t)
    A=np.sqrt(np.maximum(a*a+b*b+2*a*b*ct,0)); B=np.sqrt(np.maximum(c_*c_+d*d+2*c_*d*ct,0))
    return float(np.log(np.maximum(A,B)+1e-300).mean())

print("== V1  N2's 24-of-24 PERMUTATION INVARIANCE: does it hold at EVERY connection? ==")
print("   REGISTER:196-197 and W10_SCOPE_TABLE 3.2 state it with NO connection qualification.")
print("   One variable moving: the permutation. Same connection, same estimator, same N.")
base=(0.0,0.3,0.3,0.4); N=2_000_000
print(f"   weights {base}   m(P) = {m_jensen(base):.9f}   N = {N}")
perms=sorted(set(itertools.permutations(base)))
print(f"   {'connection':<34}{'#distinct pi':>13}{'spread of lambda':>19}{'distinct values':>17}")
for tag,f,c in [("f=1, c=sqrt2   GENERIC",1.0,2**0.5),
                ("f=2.0, c=1.1   RESONANT",2.0,1.1),
                ("f=pi, c=pi/2   ORDER 4",np.pi,np.pi/2)]:
    vals=[lam(p,f,c,N) for p in perms]
    print(f"   {tag:<34}{len(perms):>13}{max(vals)-min(vals):>19.3e}{len(set(np.round(vals,6))):>17}")
print("   -> the invariance is a theorem about m(P). It holds for lambda ONLY where lambda = m(P).")
print("      N2 inherits N1's hypothesis exactly. My W-10 scope table states it unqualified.\n")

print("== V2  W-03's INVOLUTION (00<->11, 10<->01) AT THE SAME THREE CONNECTIONS ==")
inv=lambda p:(p[3],p[2],p[1],p[0])
for tag,f,c in [("GENERIC",1.0,2**0.5),("RESONANT",2.0,1.1),("ORDER 4",np.pi,np.pi/2)]:
    d=abs(lam(base,f,c,N)-lam(inv(base),f,c,N))
    print(f"   {tag:<12} |lambda(pi) - lambda(involution pi)| = {d:.3e}")
print("   -> the involution survives at EVERY connection. Only the full S_4 needs the hypothesis.\n")

print("== V3  S1's OWN REGISTERED READY STATE IS A CLASSICAL OBJECT ==")
print("   S1 sec6 publishes p = (1/2,0,0,1/4,1/4); v0 is class 11, v1,v2 class 10, v3,v4 class 01,")
print("   so pi = (0, 0, 1/2, 1/2)  ->  P = (1/2) y (1 + x),  which FACTORS.")
pi1=(0.0,0.0,0.5,0.5)
al=(2**0.5-1)
k=np.arange(1,200001); u=np.exp(-1j*2*np.pi*al*k)
Z=0.5*np.abs(1+u)                                   # |v^k| = 1
print(f"   max | |Z_k| - |cos(pi k alpha)| | over k<=2e5 = {np.abs(Z-np.abs(np.cos(np.pi*al*k))).max():.2e}")
print(f"   m(P) = log(1/2) = {np.log(0.5):.9f}     [Jensen: m(1+x) = 0, m(y) = 0]")
print(f"   {'N':>8}{'(1/N) sum log|Z_k|':>24}")
for NN in (100,1000,10000,100000,200000):
    print(f"   {NN:>8}{float(np.log(np.abs(np.cos(np.pi*al*np.arange(1,NN+1)))+1e-300).mean()):>24.9f}")
print("   -> sum log|Z_k| is a SUDLER-TYPE PRODUCT. Sudler 1964, Erdos-Szekeres 1959, Lubinsky 1999.")
print("      The corpus's own published ready state is a 65-year-old named problem and no lane")
print("      looked it up until now.")
