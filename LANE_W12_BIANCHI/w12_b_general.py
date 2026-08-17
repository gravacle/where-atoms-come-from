# W-12 leg B — the general statement the failure of leg A points at.
#
# THEOREM (mine).  Let X be any finite complex, gamma_F and gamma_C any two designated 1-cycles
# that are R-independent.  The map  phi : (R/2piZ)^E -> T^2,  a |-> (<gamma_F,a>, <gamma_C,a>)
# is a continuous surjective homomorphism of compact connected groups.  Its image is a CONNECTED
# CLOSED subgroup of T^2, hence {1}, a circle, or T^2 -- and its dimension is the R-rank of the
# 2 x E incidence matrix, which is 2.  SO THE IMAGE IS ALL OF T^2, and phi pushes Haar to Haar.
#
# THREE COROLLARIES, none of them in the register:
#   (1) NO carrier and NO designation can constrain (W_F, W_C).  The Bianchi route is closed for
#       every complex, not just for K1 -- so leg A's null is structural, not a property of B0b.
#   (2) N1's HYPOTHESIS H2 IS PURELY ARITHMETIC.  It is a condition on the individual connection
#       and can never be violated by the carrier or by the loop designation.
#   (3) N3 UPGRADES FROM ARGUED TO PROVED, AND FROM K1-SCOPED TO CARRIER-INDEPENDENT: any
#       absolutely continuous measure on connections pushes to an absolutely continuous measure
#       on T^2, whose resonant set is Haar-null.  So no local action moves the rate ON ANY CARRIER.
import numpy as np, itertools
rng=np.random.default_rng(20260820)

def image_dim(gF,gC,E):
    M=np.zeros((2,E))
    for (e,s) in gF: M[0,e]+=s
    for (e,s) in gC: M[1,e]+=s
    return int(np.linalg.matrix_rank(M)), M

print("== B1  THE RANK, ON EVERY CARRIER AND DESIGNATION THE CORPUS OWNS ==")
CASES=[
 ("K1        gF=filled tri, gC=unfilled tri", [(0,1),(1,1),(2,1)], [(3,1),(4,1),(5,1)], 6),
 ("B0b       gF=face(0,0),  gC=row j=0",      [(0,1),(10,1),(3,-1),(9,-1)], [(0,1),(1,1),(2,1)], 18),
 ("B0b alt   gF=row j=0,    gC=row j=1",      [(0,1),(1,1),(2,1)], [(3,1),(4,1),(5,1)], 18),
 ("adversarial: loops SHARE two of three edges", [(0,1),(1,1),(2,1)], [(0,1),(1,1),(3,1)], 6),
 ("DEGENERATE: gC = gF exactly",              [(0,1),(1,1),(2,1)], [(0,1),(1,1),(2,1)], 6),
 ("DEGENERATE: gC = gF reversed",             [(0,1),(1,1),(2,1)], [(0,-1),(1,-1),(2,-1)], 6),
]
for tag,gF,gC,E in CASES:
    r,M=image_dim(gF,gC,E)
    print(f"  {tag:<44} rank = {r}   image = {'T^2' if r==2 else ('a circle' if r==1 else '{1}')}")
print()
print("  -> rank 2 for every pair of DISTINCT loops; rank 1 only when gamma_C = +/- gamma_F,")
print("     because a simple cycle has 0/+-1 coefficients, so two simple cycles are R-dependent")
print("     ONLY IF they have the same support.\n")

print("== B2  AND THE ONLY DEGENERATE DESIGNATION IS EXCLUDED BY W-02's OWN CRITERION ==")
print("  If gamma_C = +/- gamma_F the two loops have the SAME VERTEX SET, so only classes 00 and 11")
print("  are occupied and P = p00 + p11*x*y.  Two sub-cases:")
print("    gC = +gF  : (u,v) = (conj W, W), so uv = 1 and Z_k = p00 + p11 = 1 for all k.")
print("                G = {1}. NO FORMATION -- W-02's criterion excludes it.")
print("    gC = -gF  : (u,v) = (u,u), confined to the diagonal circle. H2 DOES fail.")
print("                But P = p00 + p11*x*y depends only on the PRODUCT, and on the diagonal")
print("                the product is equidistributed, so the subtorus average EQUALS m(P).")
p00,p11=0.35,0.65
# DEFECT RECORDED, NOT SILENTLY PATCHED: this first used 400000 Monte Carlo draws and printed a
# difference of 1.10e-03, which is ~2 sigma of MC noise on a bounded smooth integrand and NOT a
# real gap. Replaced with deterministic trapezoid, which converges spectrally on a periodic
# analytic integrand. The conclusion is unchanged; the number was not evidence as first computed.
th=2*np.pi*np.arange(1<<20)/(1<<20)
sub=np.log(np.abs(p00+p11*np.exp(2j*th))).mean()            # diagonal circle, u^2 equidistributed
print(f"    measured: diagonal-circle average = {sub:.12f}")
print(f"              m(p00 + p11 xy)         = {np.log(max(p00,p11)):.12f}   (Jensen, two-term)")
print(f"              difference              = {abs(sub-np.log(max(p00,p11))):.2e}")
print()
print("  ==> WHENEVER FORMATION OCCURS AT ALL, H2 CANNOT FAIL BY TOPOLOGY OR DESIGNATION.")
print("      The one designation that confines the pair also collapses P onto the confined")
print("      direction, so the average is unchanged. N1's H2 is arithmetic, and only arithmetic.\n")

print("== B3  THE PUSHFORWARD IS HAAR -- so N3 holds on EVERY carrier, not just K1 ==")
for tag,gF,gC,E in CASES[:3]:
    a=rng.uniform(0,2*np.pi,(200000,E))
    x=np.array([sum(s*a[:,e] for (e,s) in gF)]).ravel()%(2*np.pi)
    y=np.array([sum(s*a[:,e] for (e,s) in gC)]).ravel()%(2*np.pi)
    Hg,_,_=np.histogram2d(x,y,bins=24,range=[[0,2*np.pi]]*2)
    exp=200000/576
    print(f"  {tag:<44} cells {int((Hg>0).sum())}/576  chi2/dof = {float(((Hg-exp)**2/exp).sum()/575):.3f}")
print("  -> uniform on T^2 to sampling error. Any absolutely continuous measure on connections")
print("     pushes forward absolutely continuous, and the resonant set is Haar-null there.")
print("     N3 IS THEREFORE A THEOREM ON EVERY CARRIER, not an argued K1 fact.")
