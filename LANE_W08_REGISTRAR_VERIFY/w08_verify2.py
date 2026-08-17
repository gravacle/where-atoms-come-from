import numpy as np
rng=np.random.default_rng(20260816)

print("== V4  THE ADVERSARIAL SCHEDULE — does the record survive a chosen cell schedule? ==")
# RS-G generic ready state; honest schedule k_n = n vs adversary picking the J smallest (1-|Z_k|).
P0,PF,PC=0.4,0.3,0.3
def z(aF,aC,K):
    k=np.arange(1,K+1)
    return np.abs(P0*np.exp(1j*k*(aC-aF)) + PF*np.exp(-1j*k*aF) + PC*np.exp(1j*k*aC))
aF,aC=1.7231,4.9014                     # generic
print(f"  {'K':>9} {'honest -log|Omega_K|':>22} {'adversary J=sqrt(K)':>22} {'|Omega| adversary':>19}")
for K in [10**4,10**5,10**6,10**7]:
    zz=z(aF,aC,K); d=-np.log(zz)
    J=int(K**0.5); adv=np.sort(d)[:J].sum()
    print(f"  {K:>9} {d.sum():>22.3f} {adv:>22.4f} {np.exp(-adv):>19.4f}")
print("  -> honest schedule: the record decays linearly and irreversibly.")
print("     adversary writing only sqrt(K) best-chosen cells: |Omega| ~ 0.55 FOREVER,")
print("     with unboundedly many writes. Durability is a (connection, SCHEDULE) property.")
print("     The corpus never stated a schedule stipulation. This is where the obstruction survives.\n")

print("== V5  W-01's CONVEX HULL OFF K1 — one variable: the number of occupied characters ==")
# On K1 (p00 = 0) the coefficients are the THREE characters {u, v, uv}. With a spectator
# vertex (p00 > 0) there are FOUR: {1, u, v, uv}.  Same hull test, same grid, both arms.
def in_hull(pts):
    # 0 in conv of unit-modulus points  <=>  they do not all lie in an open half-plane
    a=np.sort(np.angle(pts))
    gaps=np.diff(np.concatenate([a,[a[0]+2*np.pi]]))
    return gaps.max() <= np.pi + 1e-12
N=200000
f=rng.uniform(-np.pi,np.pi,N); c=rng.uniform(-np.pi,np.pi,N)
u=np.exp(-1j*f); v=np.exp(1j*c)
three=np.array([in_hull(np.array([u[i],v[i],u[i]*v[i]])) for i in range(N)])
four =np.array([in_hull(np.array([1+0j,u[i],v[i],u[i]*v[i]])) for i in range(N)])
closed=(np.cos(f)+np.cos(c) <= 0)
print(f"  K1 (three characters, p00=0)      firing region = {three.mean():.6f}   [exact 1/4 = 0.25]")
print(f"  spectator carrier (four chars)    firing region = {four.mean():.6f}   [exact 1/2 = 0.50]")
print(f"  closed form  cos f + cos c <= 0  agrees with the four-character hull on "
      f"{int((four==closed).sum())} of {N}")
print("  -> the firing region DOUBLES on one added vertex, and the criterion acquires a closed form.\n")
print("  W-01's advertised virtue: 'it distinguishes curvature from flat holonomy'. Test it:")
uu=np.exp(+1j*f)                         # f -> -f alone, c untouched
three_neg=np.array([in_hull(np.array([uu[i],v[i],uu[i]*v[i]])) for i in range(N)])
four_neg =np.array([in_hull(np.array([1+0j,uu[i],v[i],uu[i]*v[i]])) for i in range(N)])
print(f"    f -> -f changes the verdict on K1               : {int((three!=three_neg).sum())} of {N}")
print(f"    f -> -f changes the verdict on spectator carrier: {int((four !=four_neg ).sum())} of {N}")
print("  -> the property the register advertises is a coincidence of p00 = 0.")
