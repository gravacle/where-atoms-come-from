"""S5-CHARGE  step 0 : rebuild the corpus from the sealed text and reproduce its numbers
BEFORE departing from it.  Nothing below is inherited numerically."""
import numpy as np
from fractions import Fraction
from s5lib import *

print("="*94)
print("0.  K1 REBUILT FROM S1_CARRIER_K1_V001.md sections 1 and 3, AS INCIDENCE")
print("="*94)

# S1 :19-21   e1: v0->v1  e2: v1->v2  e3: v2->v0  e4: v0->v3  e5: v3->v4  e6: v4->v0
# S1 :24      face F attached along e1.e2.e3
K1 = CW("B1  K1", 5,
        [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)],
        [ {0:1,1:1,2:1} ])
gF = {0:1,1:1,2:1}
gC = {3:1,4:1,5:1}

print("d1 (rows v0..v4, cols e1..e6):")
print(K1.d1)
print("d2 (rows e1..e6, col F):")
print(K1.d2.T, " (transposed for display)")
print("d1 @ d2 =", (K1.d1 @ K1.d2).ravel(), " max|entry| =", np.abs(K1.d1@K1.d2).max())
print("betti:", K1.betti())
print("gamma_F is a cycle:", K1.is_cycle(gF), " bounds:", K1.bounds(gF))
print("gamma_C is a cycle:", K1.is_cycle(gC), " bounds:", K1.bounds(gC))
a,b = K1.classes(gF,gC)
print("a_v =", a, "  b_v =", b, "   -> classes (a,b) per vertex:", list(zip(a,b)))

print()
print("="*94)
print("1.  UNIT CHARGE: S3/S4 CORPUS NUMBERS REPRODUCED FROM MY OWN CODE")
print("="*94)
q1 = np.ones(5, dtype=np.int64)
E1 = exponents_from_charge(a,b,q1)
p_corpus = np.array([0.4,0.15,0.15,0.15,0.15])
f0, c0 = 2.0, 1.1

ks = np.arange(1,4001)
Zk = Z_of_k(ks, E1, p_corpus, f0, c0)
az = np.abs(Zk)
print(f"min |Z_k| over k<=400  = {az[:400].min():.6f} at k = {int(ks[:400][np.argmin(az[:400])])}   (S3/S4: 0.024654 at 42)")
print(f"sup |Z_k| over k<=4000 = {az.max():.6f} at k = {int(ks[np.argmax(az)])}   (S3/S4: 0.999941 at 377)")
print(f"#(|Z_k|>0.99, k<=4000) = {(az>0.99).sum()}   (S3/S4: 37)")
for N in (1,2,5,10,20,42,50,100,200,400,1000,2000,4000):
    lam = np.log(az[:N]).mean()
    print(f"   N={N:5d}  |Z_N| = {az[N-1]:.6f}   (1/N)log|Omega_N| = {lam:.6f}")
lam200k = lambda_B_direct(E1, p_corpus, f0, c0, N=200000)
print(f"lambda_B direct N=200000 = {lam200k:.9f}    (S3/S4: -0.767026)")
Z1e5 = np.abs(Z_of_k(np.arange(1,100001), E1, p_corpus, f0, c0))
print(f"sum(1-|Z_n|), n<=100000  = {(1-Z1e5).sum():.3f}    (S3/S4: 46918.264)")

# closed form vs DIRECT MATRIX ACTION on C^5 with random section phases
rng = np.random.default_rng(20260816)          # SEED PUBLISHED
ph = rng.uniform(0, 2*np.pi, 5)
s = np.sqrt(p_corpus) * np.exp(1j*ph)
dev = max(abs(Z_of_k_matrix(k, a,b,q1, s, f0,c0) - Z_of_k(k,E1,p_corpus,f0,c0)[0]) for k in range(1,201))
print(f"closed form vs DIRECT MATRIX ACTION on C^5, max dev over k<=200 : {dev:.3e}")

# gauge invariance
gauge_dev = 0.0
for t in range(8):
    g = np.exp(1j*rng.uniform(0,2*np.pi,5))
    s2 = g*s
    gauge_dev = max(gauge_dev, max(abs(Z_of_k_matrix(k,a,b,q1,s2,f0,c0)-Z_of_k_matrix(k,a,b,q1,s,f0,c0)) for k in range(1,20)))
print(f"8 random gauge transformations, spread : {gauge_dev:.3e}")

# S1's own published connection (S1 section 6):  W_F = -1, W_C = -i
p_pub = np.array([0.5,0,0,0.25,0.25])
zz = [Z_of_k(k, E1, p_pub, np.pi, 3*np.pi/2)[0] for k in range(1,7)]
print("K1's own connection Z1..Z6 =", [f"{z.real:+.6f}{z.imag:+.6f}j" for z in zz], " (corpus: 0,-1,0,+1,0,-1)")

print()
print("="*94)
print("2.  EXACT VALUES OF RECORD, RE-DERIVED WITH MY OWN MAHLER CODE")
print("="*94)
pc = np.array([0.4,0.3,0.3]);  Ec = np.array([[1,1],[1,0],[0,1]])
print(f"generic torus  m(0.4+0.3x+0.3y)  Nx=4096  = {mahler_generic(Ec,pc,Nx=4096):.9f}   (record: -0.767507880)")
print(f"                                  Nx=16384 = {mahler_generic(Ec,pc,Nx=16384):.9f}")
print(f"L=<(11,20)>  exact                        = {lambda_B_exact(Ec,pc,((11,20),)):.9f}   (record: -0.767014993)")
print(f"L=<(1,1)>    exact                        = {lambda_B_exact(Ec,pc,((1,1),)):.9f}   (record: -1.203972804)")
print(f"L=<(3,0)>=3<(1,0)> exact                  = {lambda_B_exact(Ec,pc,((3,0),)):.9f}   (record: -0.798965257)")
print(f"B0b 4-class (4/9,1/9,2/9,2/9) SENSE U     = {mahler_generic(np.array([[0,0],[0,1],[1,0],[1,1]]),np.array([4,1,2,2])/9.,Nx=16384):.10f}   (W-03: log(4/9) = {np.log(4/9):.10f})")
print(f"SENSE C 4 classes (1/4 each)              = {mahler_generic(np.array([[0,0],[0,1],[1,0],[1,1]]),np.array([.25]*4),Nx=4096):.9f}   (record: log(1/4) = {np.log(0.25):.9f})")
print(f"m(1+x+y)/3 route: log(1/3)+m(1+x+y)       = {mahler_generic(Ec,np.array([1,1,1])/3.,Nx=16384):.9f}   (record: -0.775546341)")

# S4's exceptional table, three rows
for mn in [(1,0),(1,1),(1,-1),(2,1),(3,2),(11,20),(13,8)]:
    print(f"   L=<{mn}>  exact = {lambda_B_exact(Ec,pc,(mn,)):.9f}")
