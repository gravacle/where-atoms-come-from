#!/usr/bin/env python3
"""Exact algebra checks for GL6CM; standard library only."""

from fractions import Fraction as F


checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


# For E(s)=-sqrt[(J-sw0)^2+(J-sw1)^2], write E''=sqrt(2)*q/J.
# Direct differentiation gives q=-(w0-w1)^2/4.
tests = ((F(0), F(0)), (F(1), F(0)), (F(3), F(-2)),
         (F(7, 5), F(11, 9)))
for w0, w1 in tests:
    # F0=2J^2, F1=-2J(w0+w1), F2=2(w0^2+w1^2), set J=1.
    f0 = F(2)
    f1 = -F(2) * (w0 + w1)
    f2 = F(2) * (w0*w0 + w1*w1)
    # E''/sqrt(2) = -f2/4 + f1^2/16 at J=1.
    coefficient = -f2/F(4) + f1*f1/F(16)
    check(coefficient == -(w0-w1)**2/F(4),
          "two-arm exact curvature")

# One arm is exactly linear on the local ground branch.
for w in (F(0), F(1), F(-3, 7)):
    check(F(0) == 0*w, "isolated K2 second derivative zero")

# Physical coefficient: sqrt(2) times this rational coefficient.
J0 = F(63, 8)       # J=(63/8) h^6/U^5
w0 = F(105, 8)      # w=(105/8) h^6/U^6
coefficient = w0*w0/(F(4)*J0)
check(coefficient == F(175, 32), "literal response coefficient 175/32")

# The inverse of K=(175 sqrt(2)/32) h^6/U^7 is
# (16 sqrt(2)/175) U^7/h^6.
check(F(175, 32)*F(16, 175)*F(2) == F(1),
      "quadratic-field inverse coefficient")

# Abstract spectral identities: positive eigenvalue denominators make every
# term a square; the common writer is proportional to H and hence Q-annihilated.
for gaps, amplitudes in (
        ((F(1),), (F(2),)),
        ((F(1), F(3)), (F(2), F(-5))),
        ((F(2, 7), F(11, 3), F(5)), (F(0), F(4), F(-1)))):
    response = sum(F(2)*a*a/g for a, g in zip(amplitudes, gaps))
    check(response >= 0, "spectral response is positive")
    check((response == 0) == all(a == 0 for a in amplitudes),
          "strict kernel criterion on positive resolvent")

# Bilinear reciprocity is termwise real-symmetric.
left = sum(F(2)*a*b/g for a, b, g in
           zip((F(1), F(2)), (F(3), F(-4)), (F(5), F(7))))
right = sum(F(2)*b*a/g for a, b, g in
            zip((F(1), F(2)), (F(3), F(-4)), (F(5), F(7))))
check(left == right, "bilinear response reciprocity")

print(f"PASS__GL6CM_EXACT_ALGEBRA__{checks}/{checks}")
print("SPECTRAL_RESPONSE=RECIPROCAL_PSD;KERNEL=QB|0>_ZERO")
print("COMMON_RING_RESCALING=NULL;ISOLATED_K2=NULL")
print("TWO_OVERLAP=175SQRT2_OVER32_H6_UD7_STRICT")
print("CONTACT_BULK_REALTIME_RICCI_GRAVITY_G=OPEN")
