# W-09 leg C — the two measured regions are EXACT, and the closed forms name themselves.
# Reporting 0.2507 and 0.4989 as measurements when they are 1/4 and 1/2 exactly is the COR-E
# defect class (a window figure stated as a value). Both are derived here.
import numpy as np
from math import comb
rng=np.random.default_rng(20260816)

print("== C1  THREE CLASSES -> EXACTLY 1/4, BY WENDEL'S THEOREM (1962) ==")
print("  Any three occupied classes reduce, after dividing by one character (a rotation, which")
print("  does not move 0 relative to the hull), to {1, e^{i.theta1}, e^{i.theta2}}. For K1's set")
print("  {v,u,uv}, dividing by v gives {1, e^{-i(f+c)}, e^{-if}}, and (f+c, f) is jointly uniform")
print("  on T^2 whenever (f,c) is. So the two free angles are INDEPENDENT UNIFORM.")
print("  Wendel: P(N iid symmetric points in R^d all lie in some half-space) = 2^{-N+1} sum_{k<d} C(N-1,k).")
for N in (3,4):
    p_half = 2.0**(-N+1)*sum(comb(N-1,k) for k in range(2))
    print(f"    N={N}, d=2 :  P(all in a half-plane) = {p_half}   ->  P(0 in hull) = {1-p_half}")
print("  MEASURED (leg B, four distinct three-class sets): 0.250700, 0.250700, 0.248225, 0.248225")
print("  against EXACTLY 1/4. Monte-Carlo sigma at N=200000 is 9.7e-04; the spread is under 2 sigma.\n")

print("== C2  FOUR CLASSES -> EXACTLY 1/2, AND NOT BY WENDEL ==")
print("  With all four occupied the points are {1, u, v, uv} and uv is DETERMINED by u and v,")
print("  so Wendel does not apply. The closed form does: 0 in conv{1,u,v,uv} <=> cos f + cos c <= 0.")
N=400000
f=rng.uniform(-np.pi,np.pi,N); c=rng.uniform(-np.pi,np.pi,N)
def zih(pts):
    # DEFECT RECORDED, NOT SILENTLY FIXED: this function first read
    #   np.diff(np.concatenate([A, A[:1]+2*np.pi], 0), 0)
    # where np.diff's second POSITIONAL argument is n (difference order), not axis. n=0 returns
    # the array untouched, so C2's first run printed hull.mean() = 0.000000 and a meaningless
    # agreement count. Legs A and B pass axis=0 by keyword and were never affected; their numbers
    # stand. Exactly the class of defect this program requires be recorded rather than patched.
    A = np.sort(np.angle(np.stack(pts, axis=0)), axis=0)
    g = np.diff(np.concatenate([A, A[:1] + 2*np.pi], axis=0), axis=0)
    return g.max(axis=0) <= np.pi + 1e-12
hull=zih([np.ones_like(f,dtype=complex),np.exp(-1j*f),np.exp(1j*c),np.exp(1j*(c-f))])
cf=(np.cos(f)+np.cos(c)<=0)
print(f"    closed form agrees with the hull on {int((hull==cf).sum())} of {N}")
print("    and the map (f,c) -> (pi-f, pi-c) preserves uniform measure and sends")
print("    cos f + cos c -> -(cos f + cos c), so P(cos f + cos c <= 0) = 1/2 EXACTLY.")
print(f"    measured: {hull.mean():.6f}  against EXACTLY 1/2.\n")
print("  ==> THE FIRING REGION DOUBLES FROM EXACTLY 1/4 TO EXACTLY 1/2 WHEN THE INCIDENCE")
print("      OCCUPIES ALL FOUR CLASSES. Neither number is a measurement.")
