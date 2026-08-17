# W-14 leg A — the state side is NOT the connection side, and the difference is structural.
# W-12 closed the connection side by a theorem: a |-> (W_F,W_C) is onto T^2 with Haar pushforward,
# so no absolutely continuous measure on connections moves the rate (N3). Does the same argument
# close the STATE side? It cannot even be posed the same way, and this leg shows why.
import numpy as np
rng=np.random.default_rng(20260821)

def m_jensen(p,n=1<<20):
    """m(p00 + p10 x + p01 y + p11 xy) by Jensen in y; continuous integrand, no torus-zero noise."""
    a,b,c,d=p; t=2*np.pi*np.arange(n)/n; ct=np.cos(t)
    A=np.sqrt(np.maximum(a*a+b*b+2*a*b*ct,0)); B=np.sqrt(np.maximum(c*c+d*d+2*c*d*ct,0))
    return float(np.log(np.maximum(A,B)+1e-300).mean())

print("== A1  THE CONNECTION SIDE: lambda is CONSTANT almost everywhere (N3, now proved at W-12) ==")
pi=np.array([0.0,0.3,0.3,0.4])                       # K1's registered ready state
print(f"   pi held fixed = {pi};  m(P) = {m_jensen(pi):.12f}")
print("   the rate is m(P) for EVERY connection off a Haar-null set -- so no a.c. measure on")
print("   connections can move it. That is the wall W-12 proved holds on every carrier.\n")

print("== A2  THE STATE SIDE: lambda is CONTINUOUS AND NON-CONSTANT. There is no analogue of N3. ==")
S=rng.dirichlet(np.ones(4),4000)                     # uniform on the 3-simplex
vals=np.array([m_jensen(q, 1<<14) for q in S])
print(f"   4000 uniform draws on the simplex:  lambda ranges {vals.min():.6f} .. {vals.max():.6f}")
print(f"   spread {vals.max()-vals.min():.6f} nats, versus 0 across the whole connection torus.")
print("   -> ANY constraint on the state moves the rate. The state side cannot be closed the way")
print("      the connection side was, because there is no generic value to swallow the constraint.\n")

print("== A3  AND AVERAGING OVER STATES DOES NOT COLLAPSE, BECAUSE m IS NOT AFFINE ==")
mean_pi=S.mean(axis=0)
lhs=vals.mean(); rhs=m_jensen(mean_pi)
print(f"   E[ m(P(pi)) ]      = {lhs:.9f}")
print(f"   m(P( E[pi] ))      = {rhs:.9f}")
print(f"   difference         = {abs(lhs-rhs):.6f}   -> NOT equal: m is strictly concave-ish here,")
print("      so a prior over states is NOT equivalent to its mean state. On the connection side the")
print("      corresponding statement was an equality and that is exactly what made N3 bite.")
