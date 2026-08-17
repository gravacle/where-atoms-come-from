# LANE W11-R  LEG B -- WHAT THE FORMATION GATE ACTUALLY REQUIRES.
#
# THE CRITERION IS NOT MINE AND IT IS NOT NEW.  It is PRE-REGISTERED:
#   FOUNDING_DESIGN_V001.md:117-118, sec8 "EXTERNAL CONTACT, NAMED IN ADVANCE":
#       "S2 the trivial-connection limit must give the known trivial answer"
#   S2_FORMATION_CONDITION_ON_K1_V001.md:227-231, CHOICE LEDGER C4 (status: CLOSED):
#       "A formation condition must be a NON-CONSTANT function of the connection's gauge
#        invariants (W_F, W_C).  A condition that returns the same verdict at the trivial
#        connection as at a generic one is not a formation condition -- it is a fact about the
#        carrier's shape, true before any connection exists."
#   It DID WORK: it disqualified readings (i) and (ii) of S2 sec3.2, one of them named in S2's brief.
#   REGISTER:49 carries "correct trivial limit" as a checked property of the surviving construction.
#   W10_SCOPE_TABLE row 1.6 marks it CARRIER_INDEPENDENT, basis [T].
#
# ISOLATION LEDGER FOR THIS LEG.
#   Held fixed: the carrier, the ready state, the observable <branch_F, branch_C>, the code path,
#   the seed, the tick indices.
#   ARM 1 vs ARM 2:  ONE variable -- the CONNECTION (generic -> trivial).  Convention fixed.
#   ARM 3 vs ARM 4:  ONE variable -- the CONVENTION (circuit -> edge).  Connection fixed at TRIVIAL.
#   Arms are byte-diffed below before any number is read.
import numpy as np, w11r_lib as L
from fractions import Fraction as Fr
rng = np.random.default_rng(20260817)

GEN_K1  = np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])     # f=1.0,c=sqrt(2) idiom: the ONLY generic
TRIV_K1 = np.zeros(6)                                    #   connection the corpus publishes (W-10 N-4)
assert L.arms_differ(GEN_K1, TRIV_K1), "ARMS IDENTICAL -- leg void"

print("== B0  THE EXACT WITNESS.  K1, TRIVIAL CONNECTION, EDGE CONVENTION, GATE FIRES AT ZERO ==")
print("  At a_e = 0 the edge tick is a pure permutation:  T_F s = (s2,s0,s1,s3,s4),")
print("  T_C s = (s4,s1,s2,s0,s3), so  Z^T_1 = conj(s2)s4 + conj(s0)s1 + conj(s1)s2 + conj(s3)s0")
print("  + conj(s4)s3.  Take the RATIONAL normalised ready state s = (6/11, 7/11, -6/11, 0, 0):")
p,q,r,u,w = Fr(6,11), Fr(7,11), Fr(-6,11), Fr(0), Fr(0)
nrm = p*p+q*q+r*r+u*u+w*w
Zex = r*w + p*q + q*r + u*p + w*u
print(f"    ||s||^2 = {nrm}   (EXACT, in Fraction)")
print(f"    Z^T_1   = {Zex}   (EXACT, in Fraction)   <-- THE GATE FIRES.  W_F = W_C = 1.")
s_ex = np.array([6/11,7/11,-6/11,0,0], dtype=complex)
TF0, TC0 = L.T_edge(L.K1_LOOP_F,TRIV_K1,5), L.T_edge(L.K1_LOOP_C,TRIV_K1,5)
MF0, MC0 = L.M_circuit(L.K1_LOOP_F,TRIV_K1,5), L.M_circuit(L.K1_LOOP_C,TRIV_K1,5)
print(f"    double-precision check |Z^T_1| = {abs(L.Z(TF0,TC0,s_ex,1,1)):.3e}")
print(f"    SAME STATE, CIRCUIT convention:  |Z_1| = {abs(L.Z(MF0,MC0,s_ex,1,1)):.12f}   (identically 1)")
print("  This state has support {11,10}, so W-02's registered CARRIER_INDEPENDENT criterion")
print("  FORMATION <=> G != {1} says NO FORMATION here (G = <W_C> = {1}).  The edge convention")
print("  fires anyway, exactly, with no field.")

print("\n== B1  AND NOT ONLY ON A DEGENERATE SUPPORT.  FULL-SUPPORT ZERO AT THE TRIVIAL CONNECTION ==")
# Z^T_1 at a=0 is (1/2) x^T A x for A the adjacency of the 5-cycle v0-v1-v2-v4-v3-v0 in the
# variable order (s0,s1,s2,s3,s4).  0 lies strictly inside its numerical range, so a FULL-support
# real zero exists.  Build it from two eigenvectors.
A5 = np.zeros((5,5))
for (i,j) in [(0,1),(1,2),(2,4),(4,3),(3,0)]: A5[i,j]=A5[j,i]=1.0
ev, EV = np.linalg.eigh(A5)
vp, vm = EV[:,3], EV[:,0]                      # eigenvalues ~ +0.618 and ~ -1.618
lp, lm = ev[3], ev[0]
c2 = 1.0; c1 = np.sqrt(-lm/lp)*c2
x = c1*vp + c2*vm; x /= np.linalg.norm(x)
s_full = x.astype(complex)
print(f"  s = {np.round(np.real(s_full),9)}   ||s||^2 = {np.vdot(s_full,s_full).real:.12f}")
print(f"  pi(s) = {np.round(L.pi_of(s_full,L.K1_LOOP_F,L.K1_LOOP_C,5),9)}   (all three K1 classes occupied)")
print(f"  EDGE    convention at the TRIVIAL connection:  |Z^T_1| = {abs(L.Z(TF0,TC0,s_full,1,1)):.3e}  <-- FIRES")
print(f"  CIRCUIT convention at the TRIVIAL connection:  |Z_1|   = {abs(L.Z(MF0,MC0,s_full,1,1)):.12f}  <-- cannot fire")

print("\n== B2  IT IS NOT A KNIFE-EDGE: THE WHOLE RATE IS THERE WITH NO FIELD ==")
print("  ARM 1 (generic connection f=1.0,c=sqrt(2))  vs  ARM 2 (trivial connection).  ONE VARIABLE.")
sA = np.sqrt(np.array([0.40,0.15,0.15,0.15,0.15]))+0j          # the registrar's own three states
sB = np.sqrt(np.array([0.40,0.30,0.00,0.05,0.25]))+0j
sC = sA*np.exp(1j*np.array([0.0,1.3,-0.7,2.2,0.4]))
assert L.arms_differ(sA,sB,sC)
print(f"  {'':<34}{'state A':>15}{'state B':>15}{'state C':>15}")
for tag, aa in (("GENERIC  f=1.0,c=sqrt(2)",GEN_K1), ("TRIVIAL  a_e = 0",TRIV_K1)):
    TF,TC = L.T_edge(L.K1_LOOP_F,aa,5), L.T_edge(L.K1_LOOP_C,aa,5)
    MF,MC = L.M_circuit(L.K1_LOOP_F,aa,5), L.M_circuit(L.K1_LOOP_C,aa,5)
    re_ = [L.rate(TF,TC,s,20000) for s in (sA,sB,sC)]
    rc_ = [L.rate(MF,MC,s,20000) for s in (sA,sB,sC)]
    print(f"  EDGE    lambda/tick  {tag:<14}{re_[0]:>15.9f}{re_[1]:>15.9f}{re_[2]:>15.9f}")
    print(f"  CIRCUIT lambda/circ  {tag:<14}{rc_[0]:>15.9f}{rc_[1]:>15.9f}{rc_[2]:>15.9f}")
print("  -> at the TRIVIAL connection the CIRCUIT rate is exactly 0 (no formation, as the founding")
print("     design demands) and the EDGE rate is strongly negative (formation, with no field).")

print("\n== B3  THE SAME ON B0b, WHERE THE TWO LOOP LENGTHS DIFFER (4 and 3) ==")
GEN_B, TRIV_B = rng.uniform(0,2*np.pi,18), np.zeros(18)
assert L.arms_differ(GEN_B, TRIV_B)
wA = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wA/=wA.sum()
tB = np.sqrt(wA)+0j
for tag, aa in (("GENERIC",GEN_B), ("TRIVIAL",TRIV_B)):
    TF,TC = L.T_edge(L.B0B_LOOP_F,aa,9), L.T_edge(L.B0B_LOOP_C,aa,9)
    MF,MC = L.M_circuit(L.B0B_LOOP_F,aa,9), L.M_circuit(L.B0B_LOOP_C,aa,9)
    print(f"  {tag:8s}  EDGE lambda/tick = {L.rate(TF,TC,tB,20000):>14.9f}    "
          f"CIRCUIT lambda/circ = {L.rate(MF,MC,tB,20000):>14.9f}")
print("  -> same verdict.  The edge convention reports a durable record forming on a carrier with")
print("     ZERO connection, zero curvature and zero holonomy.")

print("\n== B4  THE RULING, IN THE CORPUS'S OWN WORDS ==")
print("  S2 sec3.2 disqualified reading (i) -- the pair (W_F z, W_C z) in C^2 -- with exactly this")
print("  test: 'it fires identically at the trivial connection, failing the criterion ... This")
print("  reading is bookkeeping wearing the costume of a result.'  And reading (ii), the")
print("  orthogonality of the two loops as chains: 'it holds for every connection and indeed")
print("  before any connection exists.  DISQUALIFIED BY THE CRITERION.  It is geometry, not")
print("  formation.'")
print("  COR-F's edge tick under the edge clock fails the SAME test in the SAME way: at a_e = 0")
print("  its Z is the overlap of two PERMUTATIONS of the ready state -- pure carrier combinatorics,")
print("  true before any connection exists.  On CHOICE LEDGER C4, which is CLOSED and entailed by")
print("  a contact point NAMED IN ADVANCE in FOUNDING_DESIGN sec8, it is not a formation condition.")
