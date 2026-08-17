"""LEG C -- WHAT THE EDGE CLOCK ACTUALLY COMPUTES.
   Claim (derived here, then verified): with T^L = M and M diagonal,
       Z^T_{L q + r} = sum_{classes ab} c^{(r)}_{ab} u^{a q} v^{b q},
       c^{(r)}_{ab} = sum_{w in class ab} conj( (T_F^r s)(w) ) * (T_C^r s)(w)   [COMPLEX]
   so the edge-clock rate is  lambda_T = (1/L) sum_{r=0}^{L-1} m(P_r),  and  P_0 = pi EXACTLY.
   i.e. N1's Mahler-measure FORM is convention-robust; only the coefficient vector generalises
   from the real probability vector pi to L complex class-pairings."""
import numpy as np, sys
sys.path.insert(0,'.')
from wcore import *
np.set_printoptions(precision=12, linewidth=200)

f, c = 1.0, np.sqrt(2.0)
k = K1(f, c)
pA = np.array([0.40,0.15,0.15,0.15,0.15]); sA = np.sqrt(pA).astype(complex)
pB = np.array([0.40,0.25,0.05,0.02,0.28]); sB = np.sqrt(pB).astype(complex)
sC = np.sqrt(pA)*np.exp(1j*np.array([0.0,0.7,-1.9,2.3,0.4]))
states = {'sA':sA, 'sB':sB, 'sC':sC}
cls = k.classes(); order = [(0,0),(1,0),(0,1),(1,1)]

def coeffs(car, s, r):
    tr = np.linalg.matrix_power(car.T('F'), r) @ s
    wr = np.linalg.matrix_power(car.T('C'), r) @ s
    prod = np.conj(tr)*wr
    return np.array([prod[cls[o]].sum() for o in order])

print("LEG C -- THE EDGE CLOCK COMPUTES AN AVERAGE OF MAHLER MEASURES.  K1, L=3, f=1.0, c=sqrt(2)\n")
N = 200000
for nm, s in states.items():
    print("  state %s" % nm)
    C = [coeffs(k, s, r) for r in range(3)]
    print("     P_0 coeffs = %s" % np.array2string(C[0], precision=9))
    print("     pi(s)      = %s      ||P_0 - pi|| = %.2e" %
          (np.array2string(k.pi(s).astype(complex), precision=9), np.linalg.norm(C[0]-k.pi(s))))
    ms = [mahler4(*C[r]) for r in range(3)]
    print("     m(P_0)=%.12f  m(P_1)=%.12f  m(P_2)=%.12f" % tuple(ms))
    pred = sum(ms)/3.0
    Ze = Z_edge(k, s, N); meas = rate(Ze)
    print("     PREDICTED edge rate (1/3)sum m(P_r) = %.12f" % pred)
    print("     MEASURED  edge rate  N=%d          = %.12f    dev = %.2e" % (N, meas, abs(pred-meas)))
    # residue-resolved check
    for r in range(3):
        idx = np.arange(N); sel = (idx+1) % 3 == r
        print("        residue r=%d: measured %.9f   m(P_%d) %.9f   dev %.1e" %
              (r, np.mean(np.log(np.abs(Ze[sel]))), r, ms[r], abs(np.mean(np.log(np.abs(Ze[sel])))-ms[r])))
    print()

print("CONSEQUENCE, STATED PLAINLY:")
print("  * the r=0 term of the edge-clock rate is EXACTLY N1's polynomial m(p00+p10x+p01y+p11xy);")
print("  * the whole edge-clock rate is still a Mahler measure average of 4-term polynomials in")
print("    the SAME two characters u,v.  The IDENTIFICATION N1 publishes survives the convention.")
print("  * what does NOT survive is 'the coefficients are the real pushforward pi'.")
