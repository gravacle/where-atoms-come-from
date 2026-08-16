# C2 REDONE EXACTLY. rho = -i has order 4, so rho^k depends only on k mod 4 — compute it that way,
# not by complex exponentiation (which accumulates error and undercounted the exact zeros).
import numpy as np
phi=(1+5**0.5)/2
print(f"  {'K':>10} | {'PUBLISHED rho=-i: min':>22} {'exact zeros':>12} | {'GENERIC: min':>14} {'zeros':>7}")
for K in [10**3,10**4,10**5,10**6,10**7]:
    k=np.arange(1,K+1)
    zeros_pub=int((k%4==0).sum())                      # rho^k = 1 exactly iff 4 | k
    dpub_min=0.0
    dgen=np.abs(np.exp(2j*np.pi*((1/phi**2)*k % 1.0))-1)
    print(f"  {K:>10} | {dpub_min:>22.3e} {zeros_pub:>12} | {dgen.min():>14.3e} {int((dgen<1e-12).sum()):>7}")
print()
print("  PUBLISHED: exactly K/4 cells annihilate the dressed record TO ZERO, at every K, forever.")
print("  GENERIC  : never zero at any k; the worst near-return floor falls like ~2pi/K.")
print()
print("  W-06 registered '1000 of 4000'.  1000 = 4000/4.  That is ord(-i) = 4, not a recurrence.")
