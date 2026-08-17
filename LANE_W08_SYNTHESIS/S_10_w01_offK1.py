# S_10 — THE WEAKEST LOAD-BEARING CLAIM, CHECKED.  W-01's criterion as the register states it
# (REGISTER_V001.md:43, "0 in the convex hull of THREE unit-modulus coefficients") off K1.
# ONE VARIABLE MOVES: the number of OCCUPIED CHARACTERS (3 on K1; 4 with a spectator vertex).
# Everything else -- the criterion, the grid, the evaluator -- is identical between the two arms.
import numpy as np
def zero_in_hull(pts):
    """0 in conv{unit-modulus pts} <=> every consecutive angular gap <= pi."""
    a = np.sort(np.mod(np.angle(pts), 2*np.pi))
    g = np.diff(np.concatenate([a, a[:1]+2*np.pi]))
    return g.max() <= np.pi + 1e-12
rng = np.random.default_rng(20260816); N = 200000
F, C = rng.uniform(0, 2*np.pi, N), rng.uniform(0, 2*np.pi, N)
u, v = np.exp(-1j*F), np.exp(1j*C)
h3 = np.array([zero_in_hull(np.array([u[i]*v[i], u[i], v[i]]))            for i in range(N)])
h4 = np.array([zero_in_hull(np.array([1+0j, u[i], v[i], u[i]*v[i]]))      for i in range(N)])
cl = (np.cos(F) + np.cos(C)) <= 0
print("== S10a  MEASURE OF THE FIRING REGION, ONE VARIABLE (occupied characters) ==")
print(f"   K1, THREE characters {{uv,u,v}}          fires on {h3.mean():.6f} of T^2   (exact 1/4)")
print(f"   spectator carrier, FOUR {{1,u,v,uv}}     fires on {h4.mean():.6f} of T^2   (exact 1/2)")
print(f"   closed form  cos f + cos c <= 0 :        agrees with the 4-char hull on "
      f"{int((h4==cl).sum())} of {N}  (mismatches {int((h4!=cl).sum())})")
print("   -> THE FIRING REGION DOUBLES.  The register's word 'three' is p00 = 0, i.e. K1's")
print("      incidence, and it is not flagged as a hypothesis anywhere in the W-01 row.")
print()
print("== S10b  W-01's ADVERTISED VIRTUE, ONE VARIABLE (the same) ==")
print("   'it DISTINGUISHES CURVATURE FROM FLAT HOLONOMY, which K1 exists to separate':")
u2 = np.exp(+1j*F)      # f -> -f, c held
h3m = np.array([zero_in_hull(np.array([u2[i]*v[i], u2[i], v[i]]))        for i in range(N)])
h4m = np.array([zero_in_hull(np.array([1+0j, u2[i], v[i], u2[i]*v[i]]))  for i in range(N)])
print(f"   K1  (3 chars): sending f -> -f alone CHANGES the verdict at {int((h3!=h3m).sum())} of {N} points")
print(f"   spectator (4): sending f -> -f alone CHANGES the verdict at {int((h4!=h4m).sum())} of {N} points")
print("   -> the criterion becomes SEPARABLE in (f,c) the moment a fourth character appears.")
print("      The property the register advertises as the criterion's virtue is a coincidence")
print("      of p00 = 0.  Witness: (pi/2, pi/2) fires on K1, (-pi/2, pi/2) does not; both fire")
print("      on the spectator carrier.")
for f0, c0 in [(np.pi/2, np.pi/2), (-np.pi/2, np.pi/2)]:
    U, V = np.exp(-1j*f0), np.exp(1j*c0)
    print(f"      f={f0:+.4f} c={c0:+.4f}   K1 fires: {zero_in_hull(np.array([U*V,U,V]))}"
          f"   spectator fires: {zero_in_hull(np.array([1+0j,U,V,U*V]))}")
