"""ADVERSARY CHECKS on T-32 (magnetic) and T-34 (NAND).

A. T-32: what happens to sign-definiteness (d) when the WRITTEN track carries
   real DATA (both bit values), as every real recorded track does (DC-free
   RLL/PRML codes force near-zero running digital sum)?
B. T-32: is the additivity justification "dipole term is 1e-3 of the field term"
   true at the separations in the lane's own tables?
C. T-32: field-free evaluation. B_applied = 0 on a stored platter. What is the
   lane's quantity dE = 2 m B s_i then?
D. T-34: image-construction interaction energy -- check the factor of 2.
   Standard grounded-plane result for two charges q at height h, lateral r:
   U_int = k q^2 (1/r - 1/sqrt(r^2+4h^2)).  The lane uses 2 * that.
E. T-34: the random-data "control" -- show the coded expression is identically 1
   for every pattern with at least one programmed cell (tautology, measures its
   own construction).
"""
import numpy as np

E    = 1.602176634e-19
EPS0 = 8.8541878128e-12
KB   = 1.380649e-23
mu0  = 4*np.pi*1e-7

print("="*90)
print("A. T-32 (d) UNDER A REAL DATA PATTERN (the lane only ran all-ones and demag)")
print("="*90)
Ms, V, B = 4.0e5, 1.26e-24, 0.5
m = Ms*V
rng = np.random.default_rng(123)
print(f"   {'N':>8}{'all-ones ratio':>16}{'random-data ratio':>19}{'DC-free-coded ratio':>21}")
for N in (10,100,1000,10000,100000):
    ones  = np.ones(N)
    data  = rng.integers(0,2,N)*2-1          # random user data, p=0.5
    # DC-free coded track: running digital sum bounded (model: balanced blocks of 8)
    blocks = N//8
    coded = np.tile(np.array([1,-1,1,-1,-1,1,-1,1]), blocks+1)[:N]
    def ratio(s):
        num, den = abs((2*m*B*s).sum()), np.abs(2*m*B*s).sum()
        return num/den
    print(f"   {N:>8}{ratio(ones):>16.6f}{ratio(data):>19.6f}{ratio(coded):>21.6f}")
print("   -> a WRITTEN track holding actual data screens like the demagnetised control;")
print("      the ratio-1 result belongs to the ALL-SAME-BIT (DC-saturated) pattern only.")
print()

print("="*90)
print("B. T-32 additivity: is the dipole term 1e-3 of the field term? (lane asserts, never computes)")
print("="*90)
for rnm in (10, 13, 20, 40, 59, 80, 160):
    r = rnm*1e-9
    Bd = mu0*m/(2*np.pi*r**3)
    print(f"   r = {rnm:>4} nm   B_dip = {Bd:.3e} T   B_dip/B_applied = {Bd/B:.3e}")
print("   -> at grain spacing (~10-13 nm, the lane's own (e) table) the ratio is 0.1-0.2,")
print("      i.e. 100-200x the claimed 1e-3. The 1e-3 figure only holds at r >~ 59 nm.")
print("      The additivity table itself contains NO interaction term at all: the computed")
print("      defect-0 is the linearity of dE*N, a model identity, not a property tested.")
print()

print("="*90)
print("C. T-32 field-free evaluation (stored platter, no applied field)")
print("="*90)
for Bv in (0.5, 0.0):
    print(f"   B = {Bv} T:  sum dE over N=1000 written grains = {2*m*Bv*1000:.3e} J")
print("   -> the lane's (a)(b)(c)(d) quantity is IDENTICALLY ZERO field-free. What survives")
print("      B=0 is the remanent moment N*m (an A m^2, not an energy) and the grain-grain")
print("      dipole energy (pattern-dependent, two-signed across pairs).")
print()

print("="*90)
print("D. T-34 image energy: factor check at r = 40 nm (pitch), h = 10 nm, eps_r = 3.9")
print("="*90)
q  = 100*E; h = 10e-9; k = 1.0/(4*np.pi*3.9*EPS0)
r  = 40e-9
U_lane    = 2*k*q*q*(1.0/r - 1.0/np.sqrt(r*r+4*h*h))
U_correct =   k*q*q*(1.0/r - 1.0/np.sqrt(r*r+4*h*h))
print(f"   lane formula   : {U_lane:.4e} J = {U_lane/E:.1f} eV")
print(f"   standard image : {U_correct:.4e} J = {U_correct/E:.1f} eV")
print("   Standard grounded-plane interaction between two charges q at height h:")
print("   U = k q^2 (1/r - 1/sqrt(r^2+4h^2)).  The lane doubles it. Exponent and all")
print("   ratio columns are unaffected (multiplicative constant); absolute J/eV and the")
print("   81x81 interaction sum (1.164 E_p -> 0.582 E_p) are 2x too large.")
print()

print("="*90)
print("E. T-34 random-data 'control' is a tautology as coded")
print("="*90)
rng = np.random.default_rng(11)
vals = set()
for _ in range(1000):
    pat = rng.integers(0,2,1000)
    if pat.sum()==0: continue
    vals.add(abs(-q*pat.sum())/(q*pat.sum()))
print(f"   distinct values of abs(-q*S)/(q*S) over draws: {vals}")
print("   -> abs(-q*S)/(q*S) == 1 identically for any S>0. The 1000 draws cannot fail by")
print("      construction; the physical content (one carrier sign) is in the premise -q,")
print("      not in the draw. ok_d's min(ratios_rand) gate is vacuous.")
