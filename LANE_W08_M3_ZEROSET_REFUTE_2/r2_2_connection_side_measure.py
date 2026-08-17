#!/usr/bin/env python3
# LANE W08 / M3 — REFUTER 2 — LENS: SCOPE.
# script 2: THE CONNECTION SIDE DOES NOT "SURVIVE UNCHANGED".  ITS FIRING REGION DOUBLES.
#
# The lane's operative-variable block says, as ruled-out alternative (i):
#     "'four classes vs three' -- refuted, the CONNECTION-side criterion (Theorem A) holds
#      unchanged at four classes, 0 disagreements on 361 grid points [m3_5 (B)]"
# and m3_5 (B) prints "=> THE CONNECTION-SIDE CRITERION SURVIVES AT FOUR CLASSES, UNCHANGED."
#
# The FORM of the criterion is unchanged and that is trivially true: for ANY finite set of unit
# vectors, {exists p in the simplex with sum p_i z_i = 0} <=> 0 in conv{z_i} <=> every angular
# gap <= pi.  That is a tautology about convex hulls; it can never depend on the class count.
# What the criterion CUTS OUT of the torus is not unchanged.  On K1's OWN PAGE (m3_5) the two
# counts are printed four lines apart and never compared:
#         three classes:  342 of 1369  = 0.2498
#         four  classes:  180 of  361  = 0.4986
# This script derives BOTH measures in closed form, EXACTLY.
#
#   THREE CLASSES {uv, u, v}   : measure exactly 1/4.
#   FOUR  CLASSES {1, u, v, uv}: 0 in conv  <=>  cos f + cos c <= 0.  Measure exactly 1/2.
#
# The four-class criterion has a closed form the lane never found, and it is as simple as the
# state-side one.  PRECISION: the closed forms are proved on paper below; the grid checks are
# double, the counts are exact integers.
import numpy as np

L = []
def out(s=""):
    print(s); L.append(s)

out("=" * 100)
out("R2-2  THE CONNECTION-SIDE FIRING REGION DOUBLES WHEN A SPECTATOR VERTEX IS ADDED")
out("=" * 100)
out("numpy %s ; IEEE double for grid checks; counts are exact integers; measures proved on paper."
    % np.__version__)
out()

out("PROOF, FOUR CLASSES.  The coefficients are z = (1, u, v, uv) with u = e^{-if}, v = e^{ic}.")
out("Because 1*(uv) = u*v, BOTH PAIRS {1,uv} and {u,v} have the same product, so both pairs are")
out("exchanged by the SAME reflection R(z) = (uv) conj(z), i.e. reflection in the line at angle")
out("m = arg(uv)/2.  Rotate by e^{-im}: the four points sit at angles  +-m  and  +-e  with")
out("      m = (c-f)/2 ... i.e. arg(uv)/2      e = arg(u/v)/2 = -(f+c)/2 ,")
out("a MIRROR-SYMMETRIC configuration about the real axis.  If an OPEN half-plane H contains all")
out("four, so does conj(H), hence so does H n conj(H) -- and that intersection is an open wedge")
out("symmetric about the real axis (on the + side or the - side), of half-angle < pi/2.  Hence")
out("      0 NOT in conv   <=>   cos m and cos e are BOTH > 0, or BOTH < 0   <=>   cos m cos e > 0.")
out("And  cos m cos e = (cos(m+e) + cos(m-e))/2 = (cos f + cos c)/2  after substituting")
out("      m+e = (c-f)/2 - (f+c)/2 = -f ,     m-e = (c-f)/2 + (f+c)/2 = c .")
out("THEREFORE   0 in conv{1,u,v,uv}   <=>   cos f + cos c <= 0.   [R2-2A]")
out("Its Haar measure is EXACTLY 1/2: (f,c) -> (pi-f, pi-c) is measure preserving and sends")
out("cos f + cos c to its negative, so the two closed regions have equal measure and overlap in")
out("the null set {cos f + cos c = 0}.")
out()

out("PROOF, THREE CLASSES (K1).  Coefficients {uv, u, v} at angles s+t, s, t with s=-f, t=c.")
out("Rotate by e^{-i(s+t)/2}: angles become  m = (s+t)/2  and  +-e  with e = (s-t)/2.  Only ONE")
out("mirror pair now, so the above collapse does not happen.  With eps = |e| reduced to [0,pi]:")
out("  eps <= pi/2 : the gap NOT containing 0 has width 2pi-2eps >= pi and m must split it:")
out("                m in [pi-eps, pi+eps], an arc of length 2eps.")
out("  eps >  pi/2 : the gap containing 0 has width 2eps > pi and m must split it:")
out("                m in [eps-pi, pi-eps], an arc of length 2(pi-eps).")
out("m is uniform on [0,2pi) and eps has density 1/pi on [0,pi], independently, so")
out("  P = (1/pi) int_0^{pi/2} (2 eps)/(2pi) d eps + (1/pi) int_{pi/2}^{pi} 2(pi-eps)/(2pi) d eps")
out("    = (1/pi^2)(pi^2/8 + pi^2/8) = 1/4   EXACTLY.   [R2-2B]")
out()


def zero_in_hull(zs, tol=1e-12):
    a = np.sort(np.mod(np.angle(np.asarray(zs)), 2 * np.pi))
    gaps = np.diff(np.concatenate([a, [a[0] + 2 * np.pi]]))
    return bool(gaps.max() <= np.pi + tol)


# ------------------------------------------------------------------ (a) verify [R2-2A]
out("(a) VERIFY [R2-2A] AGAINST THE ANGULAR-GAP CRITERION.  Random (f,c), and the lane's own grid.")
rng = np.random.default_rng(20260816)
bad = 0
for _ in range(200000):
    f, c = rng.uniform(0, 2 * np.pi, 2)
    u, v = np.exp(-1j * f), np.exp(1j * c)
    if zero_in_hull([1 + 0j, u, v, u * v]) != (np.cos(f) + np.cos(c) <= 1e-12):
        bad += 1
out("    200000 random (f,c): #{gap criterion != (cos f + cos c <= 0)} = %d" % bad)
G4 = 19
n_gap = n_cf = 0
for i in range(G4):
    for j in range(G4):
        f, c = 2 * np.pi * i / G4, 2 * np.pi * j / G4
        u, v = np.exp(-1j * f), np.exp(1j * c)
        n_gap += zero_in_hull([1 + 0j, u, v, u * v])
        n_cf += (np.cos(f) + np.cos(c) <= 1e-12)
out("    the lane's own 19x19 grid: gap criterion fires %d of %d ; [R2-2A] fires %d of %d"
    % (n_gap, G4 * G4, n_cf, G4 * G4))
out("    (the lane printed 180 of 361 = %.4f and read it as 'unchanged'.)" % (180 / 361))
out()

# ------------------------------------------------------------------ (b) the two measures
out("(b) THE TWO MEASURES, side by side.  ONE THING MOVES: whether a fourth (trivial) character")
out("    is present.  The carrier, the observable Z_1, the quantifier (EXISTS a ready state),")
out("    the criterion (0 in conv) and the grid are all held fixed.")
for G in (37, 181, 721, 1441):
    n3 = n4 = 0
    for i in range(G):
        f = 2 * np.pi * i / G
        for j in range(G):
            c = 2 * np.pi * j / G
            u, v = np.exp(-1j * f), np.exp(1j * c)
            n3 += zero_in_hull([u * v, u, v])
            n4 += zero_in_hull([1 + 0j, u, v, u * v])
    out("    grid %4dx%-4d : THREE classes (K1) %8d/%-8d = %.6f   FOUR classes (K1S) %8d/%-8d = %.6f"
        % (G, G, n3, G * G, n3 / G / G, n4, G * G, n4 / G / G))
out("    -> 1/4 and 1/2, the two closed forms [R2-2B] and [R2-2A].")
out()
out("*** THE CORRECTION.  'The connection-side criterion holds unchanged at four classes' is true")
out("*** only of the criterion's FORM, and its form could not have been otherwise: '0 in conv of")
out("*** the coefficients' is the definition of the range of a simplex-weighted sum, for any")
out("*** number of coefficients.  m3_5 (B) is therefore a CONTROL THAT COULD NOT HAVE FAILED,")
out("*** and by this program's own rule that voids the control -- not the theorem.  The CONTENT")
out("*** of the criterion doubles: half the torus fires on K1S, a quarter on K1.  The lane used")
out("*** (B) to rule out 'four classes' as the operative variable; (B) cannot do that work,")
out("*** because it is insensitive to the class count BY CONSTRUCTION.")
out()

# ------------------------------------------------------------------ (c) what this costs W-01
out("(c) WHAT THIS COSTS THE REGISTER.  W-01's row reports its criterion verified on '1369 grid")
out("    points, 0 mismatches' and reports firing on S1's published connection.  Both are K1")
out("    facts.  On K1S the same sentence, same quantifier, same proof, cuts out TWICE the torus.")
f0, c0 = np.pi, 3 * np.pi / 2
u0, v0 = np.exp(-1j * f0), np.exp(1j * c0)
out("    S1's published connection (f,c) = (pi, 3pi/2): cos f + cos c = %.6f <= 0, so it fires on"
    % (np.cos(f0) + np.cos(c0)))
out("    K1S too -- the register's headline exhibit is NOT what breaks.")
out("    What breaks is the COMPLEMENT: the set of connections at which NO ready state can fire")
out("    shrinks from 3/4 of the torus to 1/2.  A carrier's spectator vertex makes formation")
out("    STRICTLY EASIER on the connection side and STRICTLY HARDER on the state side.")
out("    (state side: firing volume is 1/4 of the simplex at BOTH class counts -- m3_1 and m3_2 --")
out("     so it is the CONNECTION side, the one the lane called unchanged, that moves.)")
out()

# ------------------------------------------------------------------ (d) three loops, out of lens
out("(d) OUT OF LENS, RECORDED NOT SCORED.  With a THIRD loop the coefficients are the 8 vertices")
out("    of a cube in the character lattice and the same reflection argument fails (the four")
out("    pairs do not share one axis).  P = A(x,y) + z B(x,y) and a zero exists iff |A| = |B|")
out("    somewhere on T^2 -- a codimension-1 condition on a 2-torus, not a finite inequality on")
out("    the weights.  No closed-form weight criterion is claimed here; only that neither M3-2")
out("    nor the polygon reading can be it.  This lane did not run three loops.")
out()
out("DONE.")
open("r2_2_connection_side_measure.OUT.txt", "w").write("\n".join(L) + "\n")
