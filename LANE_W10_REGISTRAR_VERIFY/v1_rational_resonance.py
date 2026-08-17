# Registrar's verification 1 — "EVERY RATIONAL CONNECTION IS EXACTLY RESONANT", and the census.
from fractions import Fraction as F
print("== PROOF ==")
print("  Resonance means: EXISTS integers (m,n) != (0,0) with m*f + n*c = 2*pi*j for some integer j.")
print("  If f and c are RATIONAL reals, m*f + n*c is rational for all integer m,n.")
print("  2*pi*j is irrational for every j != 0 (pi transcendental), so j = 0 and we need m*f+n*c = 0.")
print("  Write f = p/q, c = r/s in lowest terms. Then (m,n) = (r*q, -p*s) is a nonzero integer")
print("  solution: r*q*(p/q) - p*s*(r/s) = r*p - p*r = 0.  EVERY RATIONAL (f,c) IS EXACTLY RESONANT.\n")
def primitive_relation(f, c):
    f, c = F(f).limit_denominator(10**9), F(c).limit_denominator(10**9)
    m, n = c.numerator*f.denominator, -f.numerator*c.denominator
    from math import gcd
    g = gcd(abs(m), abs(n)) or 1
    m, n = m//g, n//g
    assert m*f + n*c == 0
    return m, n, m*f + n*c
print("== CENSUS OF EVERY CONNECTION THE CORPUS PUBLISHES OR USES ==")
CASES = [
 ("S3/S4 headline           f=2.0,     c=1.1",      "2.0", "1.1",      "erratum v W-02 already convicts this"),
 ("S4:973 'they are GENERIC' 3.14159,  1.57080",    "3.14159","1.57080","S4 asserts GENERIC in terms"),
 ("S4 second row            f=2.0,     c=2.0",      "2.0","2.0",       ""),
 ("W-10 lane D hard-coded   f=1.3,     c=2.0",      "1.3","2.0",       ""),
]
for tag, f, c, note in CASES:
    m, n, chk = primitive_relation(f, c)
    print(f"  {tag:<45} primitive relation ({m}, {n})   m*f+n*c = {chk}   {note}")
print()
print("  S4:603's  f = 1.0, c = sqrt(2)  : c is IRRATIONAL, so the argument does not apply.")
print("  sqrt(2) is irrational and 1, sqrt(2) are Q-linearly independent, so m*1 + n*sqrt(2) = 0")
print("  forces m = n = 0.  THIS IS THE ONLY GENERIC CONNECTION THE CORPUS PUBLISHES,")
print("  and S4 used it to verify the entire lambda column.")
