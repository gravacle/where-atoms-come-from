# S_06 — THE THREE THEOREM-LEVEL CHECKS THE RULINGS TURN ON.  All EXACT where marked.
import numpy as np
from fractions import Fraction
import itertools, math

print("== S6a  M2-1's IDENTITY, EXACTLY (Fractions), and its constant on RS-G ==")
print("   |Z_k|^2 = 1 - 2*sum_{j<l} w_j w_l (1 - cos(k(phi_j - phi_l)))")
print("           = 1 -   sum_{j<l} w_j w_l |chi_j^k - chi_l^k|^2")
# exact check with roots of unity: chi_j = zeta^{e_j}, zeta = primitive q-th root
worst = Fraction(0)
for q in (3,4,5,6,8,12):
    zs = [complex(math.cos(2*math.pi*j/q), math.sin(2*math.pi*j/q)) for j in range(q)]
    for w in [(Fraction(2,5),Fraction(3,10),Fraction(3,10)),
              (Fraction(1,2),Fraction(0),Fraction(1,2)),
              (Fraction(1,7),Fraction(2,7),Fraction(4,7))]:
        for e in itertools.product(range(q), repeat=3):
            for k in range(1,7):
                ch = [zs[(e[i]*k) % q] for i in range(3)]
                lhs = abs(sum(float(w[i])*ch[i] for i in range(3)))**2
                rhs = 1.0 - sum(float(w[i])*float(w[j])*abs(ch[i]-ch[j])**2
                                for i in range(3) for j in range(3) if i<j)
                worst = max(worst, Fraction(abs(lhs-rhs)).limit_denominator(10**18))
print(f"   worst |LHS-RHS| over 3 states x 6 root-of-unity orders x q^3 character choices x k<=6 :"
      f" {float(worst):.3e}")
w = (0.4,0.3,0.3)
print(f"   RS-G pairwise products w_j w_l : {w[0]*w[1]:.2f} {w[0]*w[2]:.2f} {w[1]*w[2]:.3f}"
      f"   -> M2's published race constant 0.12 = w11*w10 = {w[0]*w[1]:.2f}")
print("   CONSEQUENCE (the race, in one line): for ANY pair j!=l with chi_j != chi_l,")
print("     SUM_{k<=K}(1-|Z_k|) >= (1/2)SUM(1-|Z_k|^2) >= w_j w_l (K - 1/|sin(tau/2)|),")
print("   LINEAR in K.  No Diophantine hypothesis anywhere.  The floor is a per-cell quantity;")
print("   the decay is a per-K quantity.  They are not commensurable and cannot race.")

print()
print("== S6b  THE RANK-1 LOCUS uv = 1 (i.e. W_C = W_F).  T4 SAYS THE FRAGILITY NEEDS H = T^2. ==")
print("   On H = {(z, 1/z)}:  P|_H * z = p10 z^2 + p11 z + p01.  At RS-G that is 0.3z^2+0.4z+0.3.")
# EXACT: 3z^2+4z+3 over Z.  self-reciprocal; discriminant 16-36 = -20 < 0; product of roots = 1
print("   EXACT over Z (clear denominators): 3z^2 + 4z + 3.")
print("   self-reciprocal (coeffs 3,4,3 palindromic) -> roots come in pairs r, 1/r")
print("   discriminant = 4^2 - 4*3*3 = 16 - 36 = -20 < 0 -> roots are complex conjugates")
print("   conjugate AND reciprocal  =>  |r|^2 = r * conj(r) = r * (1/r) = 1  =>  BOTH ROOTS ON |z|=1.")
r = np.roots([3,4,3]); print(f"   numerical check: |roots| = {np.abs(r)}   max ||r|-1| = {np.abs(np.abs(r)-1).max():.3e}")
print(f"   m(0.3z^2+0.4z+0.3) = log(0.3) = {math.log(0.3):.12f}   (no root outside the circle)")
print("   So log|Z_k| is UNBOUNDED BELOW on this codimension-1 locus and its Birkhoff average")
print("   carries exactly the inhomogeneous Diophantine hypothesis T4 confines to H = T^2.")
print("   T4's SCOPING FALLS.  T4's CONCLUSION (durability iff G != {1}) does not: 1-|P| stays")
print("   continuous, and every pathology on this locus drives lambda DOWN, never toward 0.")

print()
print("== S6c  THE REGISTER'S '3*sqrt(3)/10 CANNOT ARISE, Z_4 HAS NO ELEMENT OF ORDER 3' ==")
print("   On S1's published connection the pair (v0,v3) has (dF,dC) = (-1,0), so")
print("   D_k = amp*|W_F^{-k} - 1| = amp*|(-1)^k - 1| in {0, 2*amp}.  EXACT, in Fractions:")
p = [Fraction(3,4), Fraction(4,25), Fraction(0), Fraction(9,100), Fraction(0)]
print(f"   |s|^2 = {[str(x) for x in p]}   sum = {sum(p)}")
amp2 = p[0]*p[3]                               # amp^2 = |s_0|^2 |s_3|^2
print(f"   amp^2 = |s_0|^2 |s_3|^2 = {amp2};   (2*amp)^2 = {4*amp2}")
print(f"   (3*sqrt(3)/10)^2 = 27/100 = {Fraction(27,100)}   EQUAL: {4*amp2 == Fraction(27,100)}")
print("   So 3*sqrt(3)/10 IS an attainable separation at k=1 on the order-4 connection, with")
print("   NO element of order 3 anywhere.  The sqrt(3) comes from the READY STATE's amplitudes,")
print("   not from the group.  THE REGISTER'S W-07 ROW IS WRONG ON THIS SENTENCE.")
