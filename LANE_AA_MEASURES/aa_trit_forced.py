# IS THE RESIDUAL TRIT FORCED? The spokes all leave the hub, so the hub's Gauss law fixes their
# SUM: one constraint on n links. Max entropy for the spokes is then n-1, not n -- a deficit of
# exactly 1, for free, at any temperature. The rim links are not jointly constrained, so max is n.
import numpy as np, itertools
def wheel(n): return [(0,k+1) for k in range(n)]+[(k+1,(k+1)%n+1) for k in range(n)], n+1
N=3; n=5; E,V=wheel(n); L=len(E)
st=[s for s in itertools.product(range(N),repeat=L)
    if all((sum(s[i] for i,(a,b) in enumerate(E) if a==v)
           -sum(s[i] for i,(a,b) in enumerate(E) if b==v))%N==0 for v in range(V))]
RIM=[n+k for k in range(n)]; SPK=list(range(n))
distinct=lambda keep: len({tuple(s[l] for l in keep) for s in st})
print(f"  distinct SPOKE configurations reachable in the physical sector: {distinct(SPK)} = {N}^{np.log(distinct(SPK))/np.log(N):.0f}")
print(f"  distinct RIM   configurations reachable in the physical sector: {distinct(RIM)} = {N}^{np.log(distinct(RIM))/np.log(N):.0f}")
print(f"  so max S(spokes) = {np.log(distinct(SPK))/np.log(N):.4f} and max S(rim) = {np.log(distinct(RIM))/np.log(N):.4f}  (log base {N})")
print(f"  measured at T=50: S(spokes) = 3.998077 -> deficit {n-3.998077:.4f};  S(rim) = 4.996157 -> deficit {n-4.996157:.4f}")
print()
print("  ==> THE TRIT IS FORCED. The hub's Gauss law fixes the sum of the spokes, so they can never")
print("      exceed n-1 at any temperature, in any state, under any Hamiltonian. SIXTH INSTANCE.")
print("      The 'residual order surviving thermalisation' is the constraint, not a record.")
