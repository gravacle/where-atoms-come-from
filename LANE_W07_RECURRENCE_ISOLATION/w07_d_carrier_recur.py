# W-07 leg D — the SAME question one level down: the FOUNDING obstruction, on the CARRIER's own Z_k.
# FOUNDING_DESIGN sec4 promotes "a finite discrete spectrum is recurrent" to the obstruction the
# construction must answer AT THE START. W-01/W-02/COR-E measure it. On WHICH connection?
# ISOLATION: p held fixed, carrier held fixed, k-range held fixed. Only the connection's arithmetic moves.
import numpy as np
# Z_k = P0*conj(W_F)^k W_C^k + PF*conj(W_F)^k + PC*W_C^k       (S2 Theorem A, three unit-modulus coeffs)
def Zsup(aF, aC, P0, PF, PC, K):
    k=np.arange(1,K+1)
    z = P0*np.exp(1j*k*(aC-aF)) + PF*np.exp(-1j*k*aF) + PC*np.exp(1j*k*aC)
    return np.abs(z)

print("== D  sup_k |Z_k| : ATTAINED, or only APPROACHED? ==")
print("   K1's published ready state p = (1/2, 0, 0, 1/4, 1/4)  ->  P0=1/2, PF=0, PC=1/2")
P0,PF,PC = 0.5, 0.0, 0.5
phi=(1+5**0.5)/2
cases = [
  ("S1 PUBLISHED  a=(pi/3 x3, pi/2 x3)   W_F=-1, W_C=-i", np.pi,        3*np.pi/2),
  ("GENERIC       (badly approximable)                  ", 2*np.pi*phi,  2*np.pi*phi**2),
  ("S3/S4 HEADLINE f=2.0, c=1.1  (exactly resonant)     ", 2.0,          1.1),
]
print(f"   {'connection':<54} {'max|Z_k|':>18} {'k<=1e6 with |Z_k|>1-1e-12':>28}")
for tag,aF,aC in cases:
    d=Zsup(aF,aC,P0,PF,PC,10**6)
    print(f"   {tag:<54} {d.max():>18.15f} {int((d>1-1e-12).sum()):>28}")
print()
print("   Same, with a GENERIC ready state P0=0.4, PF=0.3, PC=0.3 (all three coefficients live):")
P0,PF,PC = 0.4,0.3,0.3
for tag,aF,aC in cases:
    d=Zsup(aF,aC,P0,PF,PC,10**6)
    print(f"   {tag:<54} {d.max():>18.15f} {int((d>1-1e-12).sum()):>28}")
print()
print("   sup|Z_k| = 1 is ATTAINED, exactly and periodically, on S1's published connection.")
print("   Off it, |Z_k| < 1 for every k: approached by Kronecker, never reached.")
print("   COR-E already said 'the true supremum is 1' and labelled it a LOWER BOUND.")
print("   Nothing in the record distinguishes ATTAINED from APPROACHED. They are different obstructions.")
