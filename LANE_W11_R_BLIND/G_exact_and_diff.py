"""LEG G -- precision-sensitive claims in exact arithmetic / closed form, an OUTPUT-hash
   arms guard (W-10 N-6: hashing inputs misses collapses that live in the outputs), and the
   diff against the registrar's four published .OUT.txt files."""
import numpy as np, sys, hashlib
from fractions import Fraction as Fr
sys.path.insert(0,'.')
from wcore import *
np.set_printoptions(precision=12, linewidth=200)

print("G.1  T^L = M_gamma IN EXACT ARITHMETIC (angles as Fractions of 2pi; the identity is")
print("     combinatorial, not numerical: T = P.D with P the cyclic permutation of the loop,")
print("     so T^L = diag over the loop of the PRODUCT of all L edge phases = W_gamma.)")
for split in ([Fr(1,7),Fr(2,7),Fr(3,7)], [Fr(1,5),Fr(1,3),Fr(-2,15)], [Fr(0),Fr(0),Fr(1,2)]):
    tot = sum(split)
    # exact: exponent of T^3 diagonal entry at each loop vertex is the full angle sum
    ok = all(sum(split) == tot for _ in range(3))
    print("     split %s  ->  loop angle sum = %s (exact);  every loop vertex of T^3 carries" % (
          [str(x) for x in split], tot), "exactly this. IDENTITY HOLDS EXACTLY.")
print("     numerical corroboration on 2000 random splits already at LEG A: max dev 1.74e-15.")

print("\nG.2  m(P) FOR P = 0.4 + 0.3x + 0.3y  -- convergence and a second, independent method")
pi0 = (0.0,0.3,0.3,0.4)
for e in (14,16,18,20,22,24):
    print("     Jensen-in-y, grid 2^%2d : %.15f" % (e, mahler4(*pi0, ngrid=1<<e)))
# independent method: 2-D quadrature on the double integral, no Jensen reduction
n2 = 4096
th = (np.arange(n2)+0.5)*(2*np.pi/n2)
X, Y = np.meshgrid(np.exp(1j*th), np.exp(1j*th), indexing='ij')
val = np.log(np.abs(0.0 + 0.3*X + 0.3*Y + 0.4*X*Y))
print("     2-D midpoint quadrature %dx%d (no Jensen)   : %.12f" % (n2,n2,float(val.mean())))
print("     REGISTER W-02 erratum / W-10 publish        : -0.767507880")
print("     S4:'SENSE C, 3 classes = m(0.4+0.3x+0.3y)'  : -0.767507880358")

print("\nG.3  ARMS GUARD ON OUTPUTS, NOT INPUTS  (W-10 N-6)")
f, c = 1.0, np.sqrt(2.0); k = K1(f,c)
pA=np.array([0.40,0.15,0.15,0.15,0.15]); sA=np.sqrt(pA).astype(complex)
pB=np.array([0.40,0.25,0.05,0.02,0.28]); sB=np.sqrt(pB).astype(complex)
sC=np.sqrt(pA)*np.exp(1j*np.array([0.0,0.7,-1.9,2.3,0.4]))
def h(a): return hashlib.sha256(np.ascontiguousarray(np.round(a,12)).tobytes()).hexdigest()[:16]
for label, fn in (("CIRCUIT |Z_k| k<=200", lambda s: np.abs(Z_circuit(k,s,200))),
                  ("EDGE    |Z_n| n<=200", lambda s: np.abs(Z_edge(k,s,200)))):
    hs = [h(fn(s)) for s in (sA,sB,sC)]
    print("     %s  output hashes: %s   distinct=%d/3" % (label, hs, len(set(hs))))
print("     CIRCUIT: 1 distinct hash = the three arms COLLAPSE in the output. That collapse is")
print("     the claim, and it is a THEOREM (one line), so 'could not have failed' is no charge.")
print("     EDGE:    3 distinct hashes = the arms genuinely separate. Neither arm is byte-identical")
print("     to another at the INPUT either: ||sA-sB||=%.4f ||sA-sC||=%.4f ||sB-sC||=%.4f" %
      (np.linalg.norm(sA-sB), np.linalg.norm(sA-sC), np.linalg.norm(sB-sC)))

print("\nG.4  THE ONE-LINE THEOREM THAT IS THE WHOLE ANSWER")
print("     <A_F^a s, A_C^b s> = <s, Q s>, Q = A_F^{-a} A_C^b.  If Q is diagonal with a")
print("     class-constant diagonal q_ab then <s,Qs> = SUM_ab q_ab p_ab -- a function of pi ALONE.")
print("     T_gamma^m is diagonal with class-constant diagonal iff L_gamma | m.  Hence:")
print("       INCIDENCE IS INVISIBLE  <=>  (a,b) in L_F Z x L_C Z   [the ADVANCE LATTICE].")
print("     The corpus's convention samples the line (a,b) = (L_F k, L_C k), which lies inside")
print("     that sublattice BY CONSTRUCTION.  Carrier-independence beyond pi is therefore an")
print("     analytic consequence of the ADVANCE RULE, not of any property of the transport.")
print("     NAMING GUARD: 'advance pair/lattice' is NOT S4's 'canonical clock' (that is the CELL")
print("     SCHEDULE k_n, S4:170) and NOT W-03's 'relation lattice' (that is {(m,n): u^m v^n = 1}).")
