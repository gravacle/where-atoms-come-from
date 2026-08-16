"""RC-5  THE DECISIVE CASE, IN CLOSED FORM, HAND-CHECKABLE.

The claim's wall (1) says S4's and W-03's results HOLD for class-homogeneous charge --
"any single integer q != 0 shared by every vertex of a class".  The strictest reading of
that is a GLOBALLY CONSTANT q.  Here is an exact, four-term, hand-checkable computation
in which globally constant q = 2 moves lambda_B.
"""
import numpy as np
from math import pi, log, sqrt
import rclib as R

p = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
f = c = pi / 2                       # u = e^{-i pi/2} = -i ,  v = e^{i pi/2} = i
print("=" * 78)
print("RC-5  GLOBALLY CONSTANT CHARGE MOVES lambda_B ON THE EXCEPTIONAL SET")
print("=" * 78)
print(f"\n   connection  f = c = pi/2  ->  u = -i , v = +i")
rk, gens = R.relation_lattice_rank(f, c, mmax=40)
print(f"   relation lattice L: rank {rk}, generators include {gens[:4]}  -- rank L = 2,")
print("   so this is S4 TIER 1, where every exact value of record lives.")
print("   L is IDENTICAL for every charge: it depends on (f,c) only.")

for q in (1, 2, 3, 4):
    chi = []
    for e in R.support_points([q]*5, p)[0]:
        chi.append(np.exp(1j * (-f * e[0] + c * e[1])))
    pts, w = R.support_points([q]*5, p)
    Zs = [R.Z_closed(n, f, c, [q]*5, p) for n in range(1, 13)]
    mags = [abs(z) for z in Zs]
    # period of the orbit
    lam, T = R.lambda_B_finite_orbit(f, c, [q]*5, p, None)
    print(f"\n   q = {q} (globally constant):")
    print(f"      exponent points {[tuple(map(int,x)) for x in pts]}  weights {[round(float(x),3) for x in w]}")
    print(f"      characters chi = {[complex(round(z.real,9), round(z.imag,9)) for z in chi]}")
    print(f"      |Z_n|, n=1..12 : {[round(m,9) for m in mags]}")
    print(f"      orbit period T = {T}")
    print(f"      lambda_B (exact finite average) = {lam:.15f}")

print("\n   HAND CHECK, q = 1:  Z_n = 0.4 + 0.3(-i)^n + 0.3(i)^n")
print("      n=1: 0.4   n=2: -0.2   n=3: 0.4   n=4: 1.0")
print("      lambda_B = (1/4) log(0.4 * 0.2 * 0.4 * 1.0) = (1/4) log(0.032)")
print(f"               = {0.25*log(0.032):.15f}")
print("\n   HAND CHECK, q = 2:  exponents (2,2),(2,0),(0,2); chi = (1,-1,-1)")
print("      Z_n = 0.4 + 0.6(-1)^n :  n odd -> -0.2 , n even -> 1.0")
print("      lambda_B = (1/2) log(0.2 * 1.0) = (1/2) log(0.2)")
print(f"               = {0.5*log(0.2):.15f}")
print(f"\n   DIFFERENCE = {abs(0.25*log(0.032) - 0.5*log(0.2)):.15f}")
print("   Same carrier.  Same ready state.  Same weights.  Same relation lattice L.")
print("   Charge globally constant -- ZERO inhomogeneity of any kind, within a class or")
print("   between classes.  |S| = 3, rank Delta = 2 throughout.  lambda_B MOVES.")
print("\n   ==> 'S4's and W-03's results hold for CLASS-HOMOGENEOUS charge' is FALSE at")
print("       the very first non-unit charge, on the set where S4's Tier-1 and Tier-2")
print("       exact values -- including B0b = log(4/9) -- are defined.")
print("       S4's Q3 headline 'lambda_B is a function of L alone' does not survive q=2.")

print("\n--- and the TAXONOMY's content moves too, at globally constant charge:")
for q in (1, 2, 3):
    p0C = np.array([0.5, 0.0, 0.0, 0.25, 0.25])       # support {0, C}
    pts, w = R.support_points([q]*5, p0C)
    D = pts[1:] - pts[0]
    print(f"   q={q}  support {{0,C}}  points {[tuple(map(int,x)) for x in pts]}  Delta gen {tuple(int(x) for x in D[0])}")
    print(f"        S4-1 says 'lambda sees W_F only'; with charge it is W_F^{q} :")
    print(f"        non-formation locus = {{ f in (2pi/{q}) Z }} = {q} circle(s), not 1.")
