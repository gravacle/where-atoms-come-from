"""V1 -- ADVERSARIAL CHECK: does the lane's ENTIRE float layer collapse to ONE BATH QUBIT?

Claim under test (mine, not theirs): every chi number in LANE_EXACT_C_FORM is produced by a
model that contains NO carrier, NO code, NO records, NO n and NO geometry.  Specifically:

  * on the code space the lane sets Hs = -2*I (a multiple of the identity) and rho0 = I/d,
  * the m coupling operators are mutually commuting +-1 operators all coupled to bath site 0,
  * Environment.HB is a SUM OF SINGLE-QUBIT terms with NO bath-bath coupling and the probe is
    a single bath qubit's X.

Therefore the joint evolution block-diagonalises into sectors s in {+-1}^m and inside each
sector the bath qubit that is coupled evolves under   h_K = Z + lam*K*X ,  K = sum_a s_a ,
while every other bath qubit is untouched and cancels exactly out of the Holevo difference.

If that is right, I can reproduce the lane's headline numbers -- 0.521527300760,
0.136408688972, the whole m = 1..8 saturation table -- from a 2x2 matrix, with no reference to
codes, records, separation or N.  If I can, then the lane's "exact" nulls over separation,
orientation and carrier size are properties of a single qubit, not discoveries about records.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM")
from lane_utils import *
from itertools import product

I2 = np.eye(2); Xp = np.array([[0,1],[1,0]], complex); Zp = np.array([[1,0],[0,-1]], complex)
TIMES_ = np.linspace(1.0, 13.0, 25)

def vn(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
    return float(-(e*np.log2(e)).sum())

def toy_chi(m, lam, beta=2.0, times=TIMES_):
    """ONE BATH QUBIT.  m mutually commuting records sharing that qubit.  Readout = record 1."""
    w = np.array([1.0, -1.0])                      # H_b = Z, energies 1.0
    p = np.exp(-beta*w); p /= p.sum()
    tau = np.diag(p).astype(complex)
    # branch weights: K = s1 + sum_{a>=2} s_a
    from math import comb
    acc = 0.0
    for t in times:
        rp = np.zeros((2,2), complex); rm = np.zeros((2,2), complex)
        for j in range(m):                          # j = number of +1 among the other m-1... 
            pass
        for j in range(m):                          # placeholder
            break
        for j in range(m):                          # (real loop below)
            break
        tot = 0.0
        for j in range(m):                          # unused
            break
        for j in range(0, m):                       # unused
            break
        for j in range(0, m):                       # unused
            break
        # enumerate the other m-1 signs by their count of +1
        for j in range(0, m):                       # j = #(+1) among the other m-1, 0..m-1
            if j > m-1: break
            wgt = comb(m-1, j) / 2**(m-1)
            Ko = 2*j - (m-1)
            for s1, target in ((+1, 'p'), (-1, 'm')):
                K = s1 + Ko
                h = Zp + lam*K*Xp
                ev, U = np.linalg.eigh(h)
                ph = np.exp(-1j*ev*t)
                Ut = (U*ph) @ U.conj().T
                sig = Ut @ tau @ Ut.conj().T
                if s1 == +1: rp += wgt*sig
                else:        rm += wgt*sig
        av = (rp + rm)/2
        acc += max(vn(av) - 0.5*vn(rp) - 0.5*vn(rm), 0.0)
    return acc/len(times)

print("="*100)
print("V1(a)  ONE-BATH-QUBIT TOY vs THE LANE'S OWN PUBLISHED NUMBERS")
print("="*100)
ref = {(1,0.4):0.276635159834, (1,0.8):0.521527300760, (1,1.2):0.599566691971,
       (2,0.4):0.125625482447, (2,0.8):0.136408688972, (2,1.2):0.127539682306}
print(f"  {'m':>3}{'lam':>6}{'TOY (1 bath qubit, 2x2)':>26}{'LANE (code dim 64..1024)':>26}{'|diff|':>12}")
worst = 0.0
for (m, lam), v in sorted(ref.items()):
    t = toy_chi(m, lam)
    worst = max(worst, abs(t-v))
    print(f"  {m:>3}{lam:>6}{t:>26.12f}{v:>26.12f}{abs(t-v):>12.2e}")
print(f"  WORST DISAGREEMENT: {worst:.3e}")

print()
print("="*100)
print("V1(b)  THE WHOLE m = 1..8 'SATURATION' TABLE FROM THE SAME 2x2 MODEL  (lam = 0.8)")
print("="*100)
lane = {1:0.521527300760, 2:0.272817377944, 3:0.215349659532, 4:0.180650299961,
        5:0.167652293211, 6:0.145842259357, 7:0.141549402657, 8:0.124984198966}
print(f"  {'m':>3}{'TOY  Qchi = m*chi':>22}{'LANE Qchi':>18}{'|diff|':>12}"
      f"{'  <- Qchi is the TOTAL Holevo of ONE bath qubit; bound = 1 bit'}")
w2 = 0.0
for m in range(1, 9):
    t = m*toy_chi(m, 0.8)
    w2 = max(w2, abs(t-lane[m]))
    print(f"  {m:>3}{t:>22.12f}{lane[m]:>18.12f}{abs(t-lane[m]):>12.2e}")
print(f"  WORST DISAGREEMENT: {w2:.3e}")
print("  ORDINARY EXPLANATION: Qchi(m) is the total Holevo information m records write into a")
print("  SINGLE QUBIT.  A qubit holds at most 1 bit.  Qchi(8) = %.6f." % (8*toy_chi(8,0.8)))
print("  Hence defect = m*chi(1) - Qchi(m) -> chi(1)*m - (bounded)  ==>  LINEAR IN m, NECESSARILY.")
print("  That is C-36, already registered, not a new form result.")
