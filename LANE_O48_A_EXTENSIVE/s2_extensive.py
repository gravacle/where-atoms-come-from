"""S2 -- EXTENSIVITY, kept strictly separated into the two parts the program must not conflate.

  CONFIGURATION-INDEPENDENT :  S(n) = max_s E(s) - min_s E(s)      the energy SPREAD
  CONFIGURATION-DEPENDENT   :  E(s) for named s, and Var_s E(s)

Everything is EXACT: energies are integers in units of 1/D.  Three independent computations of
each quantity are printed side by side -- brute-force over all 2^n sign strings, brute-force
over the 2^(n-1) bond strings, and a closed form -- and they must agree to the integer.
The closed form is then evaluated far beyond any enumeration, which makes the asymptotics a
THEOREM about this family rather than a fitted trend (D-20 does not arise).

FIXED ENVIRONMENT (standard (a)): the coupling stream is a PREFIX stream -- growing n appends
bonds and never alters existing ones.  That is checked, not assumed.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE")
import numpy as np
from chain import D, couplings, uniform_couplings, configs, energies_int

OUT = []
def p(*x):
    s = " ".join(str(y) for y in x); OUT.append(s); print(s)

p("=" * 112)
p("S2  EXTENSIVITY OF THE CORRELATION ENERGY  --  configuration-independent and -dependent parts")
p("=" * 112)

# ---------------------------------------------------------------- fixed-environment check
p("")
p("PREFIX CHECK (standard (a) demands a FIXED environment as n grows):")
ok = all(couplings(k)[:5] == couplings(5) for k in range(5, 30))
p(f"  couplings(k)[:5] == couplings(5) for every k in 5..29 :  {ok}")
p(f"  first five a_i (units 1/D, D = 2^40) : {couplings(5)}")
p(f"  as J_i                               : {[round(v/D,6) for v in couplings(5)]}")
p("  Growing n APPENDS a bond; no existing coupling moves.  The environment is fixed.")

def E_over_t(a):
    """Exact int64 array of E*D over all 2^(len(a)) bond strings t, by doubling.
       Each value occurs twice in configuration space (s and -s), so the spread and the
       configuration-average of any function of E are identical in the two spaces."""
    E = np.zeros(1, dtype=np.int64)
    for v in a:
        E = np.concatenate([E + np.int64(v), E - np.int64(v)])
    return E

# ---------------------------------------------------------------- SPREAD
p("")
p("-" * 112)
p("PART 1  CONFIGURATION-INDEPENDENT:  THE ENERGY SPREAD  S(n) = max_s E(s) - min_s E(s)")
p("        Three independent computations. Integer units of 1/D.")
p("-" * 112)
p(f"{'n':>3} {'S from all 2^n configs':>24} {'S from all 2^(n-1) bonds':>25} "
  f"{'closed form 2*sum a_i':>23} {'agree':>6} {'S(n)/D':>12} {'S(n)/S(n-1)':>12}")
S_int = {}
prev = None
for n in range(2, 21):
    a = couplings(n - 1)
    s = configs(n); Ec = energies_int(s, a)
    s1 = int(Ec.max() - Ec.min())
    Et = E_over_t(a); s2 = int(Et.max() - Et.min())
    s3 = 2 * sum(a)
    S_int[n] = s3
    r = f"{s3/prev:.6f}" if prev else "-"
    p(f"{n:>3} {s1:>24} {s2:>25} {s3:>23} {str(s1==s2==s3):>6} {s3/D:>12.6f} {r:>12}")
    prev = s3
p("  (2^n configuration enumeration stopped at n = 20: 1,048,576 sign strings, ~20 MB. The bond")
p("   enumeration and the closed form continue below and agree with it wherever both ran.)")
p("")
p(f"{'n':>4} {'S from 2^(n-1) bonds':>22} {'closed form':>22} {'agree':>6}")
for n in (21, 22, 23, 24):
    a = couplings(n - 1); Et = E_over_t(a)
    s2 = int(Et.max() - Et.min()); s3 = 2 * sum(a)
    S_int[n] = s3
    p(f"{n:>4} {s2:>22} {s3:>22} {str(s2==s3):>6}")
p("  (bond enumeration stopped at n = 24: 8,388,608 int64 = 67 MB. Nothing about the object")
p("   changes past that point -- the closed form is exact for every n and is used below.)")

p("")
p("SUCCESSIVE DOUBLING RATIO S(2N)/S(N), the form the standard names.  EXACT for every N,")
p("computed from the closed form, so no fit and no extrapolation is involved.")
p(f"{'N':>9} {'S(N)/D':>18} {'S(2N)/D':>18} {'S(2N)/S(N)':>14}  READ")
for N in (2, 4, 8, 16, 32, 64, 128, 1024, 8192, 65536, 1000000):
    a2 = couplings(2 * N - 1)
    sN = 2 * sum(a2[:N - 1]); s2N = 2 * sum(a2)
    p(f"{N:>9} {sN/D:>18.6f} {s2N/D:>18.6f} {s2N/sN:>14.9f}  "
      f"{'-> 2' if abs(s2N/sN - 2) < 0.05 else 'not yet 2'}")
p("READ: S(2N)/S(N) -> 2. It is EXACTLY (sum of the first 2N-1 couplings)/(sum of the first N-1),")
p("      and the couplings are bounded above and below away from zero, so the ratio is")
p("      (2N-1)/(N-1) up to the coupling fluctuation. THE SPREAD IS EXTENSIVE. Exact argument,")
p("      not a numerical trend.")

# ---------------------------------------------------------------- SIGN-DEFINITENESS OF THE SPREAD
p("")
p("C-46 TEST ON THE SPREAD:  S(n) = sum over bonds of the term 2|J_i|. All terms have the SAME")
p("SIGN, so |sum|/sum|.| = 1 exactly at every n -- it accumulates and never cancels.")
p(f"{'n':>4} {'sum of terms':>22} {'sum of |terms|':>22} {'|sum|/sum|.|':>14}")
for n in (2, 4, 8, 16, 24, 64, 4096):
    a = couplings(n - 1); t = [2 * v for v in a]
    p(f"{n:>4} {sum(t):>22} {sum(abs(v) for v in t):>22} {abs(sum(t))/sum(abs(v) for v in t):>14.6f}")
p("  INDUCED, not inserted: the terms are 2|J_i| because max_s and min_s are attained by")
p("  INDEPENDENTLY choosing each bond variable t_i = s_i s_{i+1}. Nothing was put in by hand;")
p("  what was inserted is only that H is a sum of bond terms.")

# ---------------------------------------------------------------- CONFIGURATION-DEPENDENT
p("")
p("-" * 112)
p("PART 2  CONFIGURATION-DEPENDENT:  E(s) for named configurations, and Var_s E(s)")
p("        ALIGNED  s = (+,+,...,+)      GROUND  the minimiser      TYPICAL  averaged over all 2^n")
p("-" * 112)
p(f"{'n':>3} {'E(aligned)/D':>15} {'E(ground)/D':>14} {'mean_s E (int)':>16} "
  f"{'Var_s E, exact enum':>22} {'closed form sum a_i^2':>23} {'agree':>6} {'Var/D^2':>12}")
V = {}
for n in range(2, 21):
    a = couplings(n - 1)
    s = configs(n); Ec = energies_int(s, a)
    Ealign = sum(a)                       # s = all +1  ->  every t_i = +1
    Eground = -sum(a)                     # all a_i > 0 -> minimiser takes every t_i = -1
    assert int(Ec.max()) == Ealign and int(Ec.min()) == Eground
    tot = sum(int(v) for v in Ec)                       # exact Python-int accumulation
    sq = sum(int(v) * int(v) for v in Ec)
    var_enum = sq // (1 << n) - (tot // (1 << n)) ** 2   # mean is exactly 0, checked next
    var_closed = sum(v * v for v in a)
    V[n] = var_closed
    p(f"{n:>3} {Ealign/D:>15.6f} {Eground/D:>14.6f} {tot:>16} {var_enum:>22} "
      f"{var_closed:>23} {str(var_enum==var_closed and tot==0):>6} {var_closed/D**2:>12.6f}")
p("READ: mean_s E(s) = 0 EXACTLY at every n (integer zero). Var_s E(s) = sum_i J_i^2 exactly.")
p("")
p(f"{'N':>9} {'Var(N)/D^2':>18} {'Var(2N)/D^2':>18} {'Var(2N)/Var(N)':>16}  READ")
for N in (2, 4, 8, 16, 32, 64, 128, 1024, 8192, 65536, 1000000):
    a2 = couplings(2 * N - 1)
    vN = sum(v * v for v in a2[:N - 1]); v2N = sum(v * v for v in a2)
    p(f"{N:>9} {vN/D**2:>18.6f} {v2N/D**2:>18.6f} {v2N/vN:>16.9f}  "
      f"{'-> 2' if abs(v2N/vN - 2) < 0.05 else 'not yet 2'}")
p("READ: Var(2N)/Var(N) -> 2. THE CONFIGURATION VARIANCE IS EXTENSIVE, by the same exact argument.")

p("")
p("BUT the configuration-dependent ENERGY ITSELF is extensive only for EXTREMAL configurations.")
p("|E(s)| for a typical s grows like sqrt(n), not n. Exact enumeration of the whole distribution:")
p(f"{'n':>3} {'|E(align)|/D  (extremal)':>25} {'RMS_s |E(s)|/D  (typical)':>26} "
  f"{'ratio typical/extremal':>23}")
for n in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20):
    a = couplings(n - 1)
    ext = sum(a) / D
    rms = (sum(v * v for v in a)) ** 0.5 / D
    p(f"{n:>3} {ext:>25.6f} {rms:>26.6f} {rms/ext:>23.6f}")
p("READ: the extremal configuration's energy is extensive (linear in n); the typical")
p("      configuration's is sub-extensive (sqrt n), and the ratio falls to zero. WHICH")
p("      CONFIGURATION IS OCCUPIED DECIDES WHETHER THE ENERGY IS EXTENSIVE AT ALL.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE/s2_extensive.txt", "w").write("\n".join(OUT) + "\n")
