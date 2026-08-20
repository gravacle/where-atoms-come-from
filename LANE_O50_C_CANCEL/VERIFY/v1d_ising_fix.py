"""VERIFY 1-D (CORRECTED).  My first pass used a RING correlation; the lane's transfer matrix is an
OPEN chain (free boundaries).  Redone with the open chain, which is what the lane actually used."""
from fractions import Fraction
import math
def line(s=""): print(s, flush=True)
def var_open(m, t):
    return m + 2 * sum((m - d) * t ** d for d in range(1, m))
line("D (CORRECTED).  OPEN-CHAIN Ising, <s_i s_j> = t^|i-j|.  Var = m + 2*sum_d (m-d) t^d.")
line(f"  {'m':>5} {'t':>10} {'Var exact':>22} {'lane Var':>12} {'Var/m':>9} {'(1+t)/(1-t)':>12} {'match?':>7}")
lane = {(16, Fraction(1,3)): "30.50000", (16, Fraction(3,5)): "56.50212",
        (16, Fraction(9,11)): "112.49619", (16, Fraction(99,101)): "230.84338",
        (18, Fraction(1,3)): "34.50000", (18, Fraction(9,11)): "131.83629",
        (18, Fraction(99,101)): "288.49093"}
for (m, t), lv in sorted(lane.items()):
    v = var_open(m, t)
    line(f"  {m:>5} {float(t):>10.6f} {float(v):>22.5f} {lv:>12} {float(v)/m:>9.4f} "
         f"{float((1+t)/(1-t)):>12.4f} {str(abs(float(v)-float(lv))<5e-5):>7}")
line()
line("  LARGE m, coh*sqrt(m) = sqrt(2/pi)*sqrt(Var/m), open chain, against the lane's I-3 table:")
lane_i3 = {(64, .5): 1.367505, (256, .5): 1.378373, (1024, .5): 1.381077, (4096, .5): 1.381752,
           (16384, .5): 1.381920, (262144, .5): 1.381973,
           (64, .818182): 2.423595, (4096, .818182): 2.521607, (262144, .818182): 2.523109,
           (64, .998002): 6.249793, (4096, .998002): 23.641679, (262144, .998002): 25.207251}
line(f"  {'m':>8} {'t':>10} {'mine':>12} {'lane':>12} {'|diff|':>10}")
for (m, t), lv in sorted(lane_i3.items()):
    s = 1.0 + 2.0 * sum((1 - d / m) * t ** d for d in range(1, min(m, 20000)))
    val = math.sqrt(2 / math.pi) * math.sqrt(s)
    line(f"  {m:>8} {t:>10.6f} {val:>12.6f} {lv:>12.6f} {abs(val-lv):>10.2e}")
line("  READ: the lane's I-3 numbers reproduce exactly on the OPEN-chain Ising measure.  Var/m ->")
line("  (1+t)/(1-t) is the textbook 1D susceptibility sum; nothing here involves the carrier.")
