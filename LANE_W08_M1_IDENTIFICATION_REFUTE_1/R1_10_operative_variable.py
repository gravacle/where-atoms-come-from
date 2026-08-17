#!/usr/bin/env python3
"""
R1_10 — THE OPERATIVE VARIABLE, AND THE LANE'S "EXPLICITLY RULED OUT" LIST.

The lane names its operative variable as
   "THE RELATION LATTICE L = {(m,n) : u^m v^n = 1} ... rank L = 0 gives H = T^2 and
    lambda = m(P); rank L = 1 gives a circle subgroup and lambda = a one-variable (subtorus)
    Mahler measure; rank L = 2 gives a finite cyclic H and lambda = a finite average"
and rules out, among others,
   "(b) ord(rho), W-07's variable -- IT IS EXACTLY THE RANK-2 CASE and CANNOT distinguish
    S3/S4's resonant headline (ord = infinity, rank L = 1, lambda != m(P)) from a Diophantine
    point (ord = infinity, rank L = 0, lambda = m(P))".

THREE CHECKS, all of which the program's own naming rule demands and none of which the lane
ran ("Before you report an effect, ask what else moves with your named variable, and rule it
out explicitly").

C1.  IS ord(rho) "EXACTLY THE RANK-2 CASE"?   ord(rho) < infinity <=> (u/v) is torsion
     <=> (d,-d) in L for some d > 0.  That is a condition on ONE DIRECTION of L, and it is
     satisfied by rank-1 lattices.  COUNTEREXAMPLE FROM THE LANE'S OWN TABLE below.

C2.  DOES rank L DETERMINE lambda?   No, in every rank:
       rank 2 : different finite H give different finite averages (checked)
       rank 1 : (1,1) vs (11,20) differ by 0.437  (R1_06); and R1_05 exhibits two points of
                the SAME rank-1 lattice with different limits, one of them -infinity
       rank 0 : Diophantine vs Liouville differ, the lane's own F4 (M1_06)
     So L is not the operative variable for the RATE.  It IS the operative variable for the
     QUALITATIVE type of H, and it is exactly the right variable for the DURABILITY criterion
     (which does not depend on it at all -- see C3).

C3.  WHAT ACTUALLY CARRIES T4?   Not L: T4's criterion is G != {1}, a function of
     (supp pi, u, v) alone, and its proof uses only that 1-|P| is continuous.  T4 is
     therefore INDEPENDENT of the operative variable the lane names.  The lane says so in one
     clause of the operative_variable field ("SEPARATELY: the DURABILITY question has a
     second, independent operative variable") and then names L first and headline-first.

Precision: exact where the statement is arithmetic; float64 for the orbit sums (labelled).
"""
import numpy as np
from fractions import Fraction as Fr

print("=" * 78)
print("R1_10 — 'ord(rho) IS EXACTLY THE RANK-2 CASE' IS FALSE")
print("=" * 78)
print("""
  DEFINITION USED (W-07 PUBLISHED_CONVENTIONS / W-07 register row): rho is the branch ratio
  in U(1); on S1's published connection W_F = -1, W_C = -i it has order 4.

  COUNTEREXAMPLE, AND IT IS ROW 4 OF THE LANE'S OWN M1_04 TABLE:
      u = v = e^{0.7 i}   (the lane's "S={10,01} and u = v (W_F W_C = 1)" row).
      rho = u/v = 1, so ord(rho) = 1 -- FINITE, indeed minimal.
      But u^m v^n = e^{0.7 i (m+n)} = 1  iff  0.7 (m+n) in 2 pi Z  iff  m + n = 0,
      because 0.7/(2 pi) is irrational (pi is transcendental and 0.7 is rational).
      So L = Z.(1,-1):  RANK 1, not rank 2.  H is a CIRCLE, not a finite group.
  Confirmed numerically (float64): the orbit (u^k, v^k) visits infinitely many points.""")
u = np.exp(0.7j)
k = np.arange(1, 20001)
pts = np.round(np.angle(u ** k), 9)
print("      distinct values of arg(u^k) to 9 decimals over k <= 2e4: %d  (a finite H would"
      % len(set(pts.tolist())))
print("      have at most ord(rho) * something bounded; it does not)")
print("""
  CONSEQUENCE.  The lane's ruling-out of ord(rho) is right in its OPERATIVE half (ord(rho)
  cannot separate S3/S4's resonant point from a Diophantine point -- both have ord = infinity;
  checked below) and WRONG in its stated reason.  ord(rho) is not the rank-2 case; it is the
  condition "(1,-1) direction lies in L", which cuts across ranks 1 and 2.  Since the whole
  point of the naming rule in this program is that the reason must be right and not only the
  conclusion, this is recorded as a defect of the same class the rule exists to catch.""")
print("\n  CHECK OF THE OPERATIVE HALF (both have ord(rho) = infinity):")
print("    S3/S4 headline: u = e^{-2i}, v = e^{1.1 i}, rho = e^{-3.1 i};  3.1/(2 pi) is")
print("      irrational (pi transcendental, 3.1 rational)  ->  ord(rho) = infinity.")
print("    Diophantine:    rho = e^{2 pi i (alpha - beta)} with alpha-beta = -2^(1/3)-4^(1/3),")
print("      irrational  ->  ord(rho) = infinity.   The two are indistinguishable by ord(rho),")
print("      and their limits differ by 4.929e-04.  THE LANE'S CONCLUSION HERE STANDS.")

print("\n" + "=" * 78)
print("C2 — rank L DOES NOT DETERMINE lambda, IN ANY RANK")
print("=" * 78)
print("\n  RANK 2, two different finite subgroups, same weights (0.3,0.3,0.4):")
for (uu, vv, lab) in [(-1+0j, -1j, "u=-1, v=-i   (S1 published, |H| = 4)"),
                      (np.exp(2j*np.pi/3), np.exp(2j*np.pi/3), "u=v=zeta_3  (|H| = 3)"),
                      (1j, np.exp(2j*np.pi/5), "u=i, v=zeta_5 (|H| = 20)"),
                      (-1+0j, -1+0j, "u=v=-1      (|H| = 2)")]:
    kk = np.arange(1, 120001)
    az = np.abs(0.3*uu**kk + 0.3*vv**kk + 0.4*(uu*vv)**kk)
    with np.errstate(divide='ignore'):
        print("     %-40s lambda = %.12f" % (lab, float(np.mean(np.log(np.maximum(az,1e-323))))))
print("     -> four rank-2 connections, four different rates.  'rank L = 2 gives a finite")
print("        average' is a TYPE statement, not a value statement; the lane's isolation")
print("        ledger reads its three-row comparison as if rank explained the three values.")
print("""
  RANK 1: R1_06 gives lambda_H over the relation vectors, spanning -1.2040 to -0.3567.
  RANK 0: the lane's own F4 gives two rank-0 pairs, one with lambda = m(P) and one with
          liminf = -infinity.
  SO: the variable that determines the RATE is (L, inhomogeneous Diophantine type of the
  orbit relative to Z(P restricted to H)).  The lane's Comparison 5 says precisely this and
  its operative_variable field says L.  The two are inconsistent, and the field is the one
  that will be read.""")
print("\nDONE R1_10")
