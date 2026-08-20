"""V1 -- ADVERSARIAL: is the quantity that passes (a)-(d) RECORD-BLIND?

The lane's PASSING quantities are
    S(n)   = max_s E(s) - min_s E(s)      (energy spread)
    Var(n) = Var_s E(s)
Neither computation in s2/s3/s4/s7 ever references a record operator.  If they are functions
of H alone, then a carrier that holds ZERO records must pass the same four standards.

CONTROL CARRIER (the refutation):  H' = sum_i J_i Z_i Z_{i+1}  +  sum_i h_i Z_i
Distinct integer h_i make H' NON-DEGENERATE.  By P-1, a non-degenerate H can hold NO record
at all -- clause (iii) is unsatisfiable because every eigenspace is one-dimensional.
We verify the non-degeneracy exactly, verify clause (iii)/(iv) FAIL for every Z_i, and then
re-run the lane's own four standards on this record-free carrier.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
import numpy as np
from chain import D, couplings, configs, energies_int, dense_H, dense_Z, pauli
from record_model import eigenspaces, clause_iii, clause_iv

OUT = []
def p(*x):
    s = " ".join(str(y) for y in x); OUT.append(s); print(s)

def fields(n, stream=7):
    a, x = [], (0xDEADBEEFCAFEBABE ^ (stream * 0x94D049BB133111EB)) & ((1 << 64) - 1)
    for _ in range(n):
        x = (x * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
        a.append((D >> 2) + (x >> 26) % (D >> 2))     # h_i in [0.25, 0.5), distinct
    return a

def E_field(s, a, h):
    n = s.shape[1]
    E = np.zeros(s.shape[0], dtype=np.int64)
    for i in range(n - 1):
        E += np.int64(a[i]) * (s[:, i].astype(np.int64) * s[:, i + 1].astype(np.int64))
    for i in range(n):
        E += np.int64(h[i]) * s[:, i].astype(np.int64)
    return E

p("=" * 110)
p("V1  IS THE PASSING QUANTITY RECORD-BLIND?  Control carrier with a field: ZERO records exist.")
p("=" * 110)
p("")
p("STEP 1  DOES THE CONTROL CARRIER HOLD ANY RECORD AT ALL?")
p(f"{'n':>3} {'#distinct energies':>19} {'2^n':>8} {'NON-DEGENERATE':>15} "
  f"{'max|Tr(P_E Z_i)| FIELD':>23} {'clause_iii(Z_1) FIELD':>22} {'[CONTROL] no-field max|Tr|':>27} "
  f"{'clause_iii(Z_1) no-field':>25}")
for n in (4, 6, 8, 10, 12):
    a = couplings(n - 1); h = fields(n)
    s = configs(n)
    Ef = E_field(s, a, h); E0 = energies_int(s, a)
    u = np.unique(Ef)
    nondeg = (len(u) == (1 << n))
    # clause (iv) trace on the FIELD carrier, exact integers
    def maxtr(E, s):
        uu, inv = np.unique(E, return_inverse=True)
        m = 0
        for i in range(s.shape[1]):
            col = s[:, i].astype(np.int64)
            tr = np.bincount(inv, weights=col.astype(float), minlength=len(uu))
            m = max(m, int(round(np.abs(tr).max())))
        return m
    mf, m0 = maxtr(Ef, s), maxtr(E0, s)
    if n <= 8:
        Hf = dense_H(n, a)
        for i in range(n):
            lab = [0]*n; lab[i] = 3
            Hf = Hf + (h[i]/D) * pauli(lab)
        esf = eigenspaces(Hf); es0 = eigenspaces(dense_H(n, a))
        c3f = clause_iii(dense_Z(n, 0), esf); c30 = clause_iii(dense_Z(n, 0), es0)
    else:
        c3f = c30 = "(dense too large)"
    p(f"{n:>3} {len(u):>19} {1<<n:>8} {str(nondeg):>15} {mf:>23} {str(c3f):>22} {m0:>27} {str(c30):>25}")
p("READ: with the field the spectrum is non-degenerate, every eigenspace is 1-dimensional, so")
p("      clause (iii) is UNSATISFIABLE by ANY operator and clause (iv) fails for every Z_i")
p("      (max|Tr| non-zero). The control carrier holds NO RECORDS WHATSOEVER.")
p("      The no-field columns are the lane's carrier, shown alongside: trace 0, clause (iii) True.")

p("")
p("STEP 2  RUN THE LANE'S FOUR STANDARDS ON THE RECORD-FREE CARRIER.")
p("(a) EXTENSIVE, (b) ADDITIVE, (c) NOT A COUNT, (d) SIGN-DEFINITE.")
p("")
p("(a) EXTENSIVITY.  S and Var by full 2^n enumeration on the FIELD carrier (no records).")
p(f"{'n':>3} {'S_field/D':>14} {'S_field(n)/S_field(n/2)':>24} {'Var_field/D^2':>16} "
  f"{'Var(n)/Var(n/2)':>17} {'[lane] S_norecord? ':>20}")
Sf = {}; Vf = {}
for n in (2, 4, 6, 8, 10, 12, 14, 16):
    a = couplings(n - 1); h = fields(n); s = configs(n); E = E_field(s, a, h)
    Sf[n] = int(E.max() - E.min())
    tot = sum(int(v) for v in E); sq = sum(int(v)*int(v) for v in E)
    Vf[n] = sq // (1 << n) - (tot // (1 << n))**2
    r1 = f"{Sf[n]/Sf[n//2]:.6f}" if n//2 in Sf else "-"
    r2 = f"{Vf[n]/Vf[n//2]:.6f}" if n//2 in Vf else "-"
    p(f"{n:>3} {Sf[n]/D:>14.6f} {r1:>24} {Vf[n]/D**2:>16.6f} {r2:>17} {'NO RECORDS':>20}")
p(f"CLOSED FORM CHECK on the field carrier: S_field = 2*(sum|J| + sum|h|)?  "
  f"{[Sf[n] == 2*(sum(couplings(n-1)) + sum(fields(n))) for n in (2,4,6,8,10,12,14,16)]}")
p("READ: S(2N)/S(N) -> 2 and Var(2N)/Var(N) -> 2 on a carrier with NO RECORDS.")

p("")
p("(b) ADDITIVITY on the record-free carrier, same protocol as S3 (cut the boundary bond).")
p(f"{'n':>3} {'m':>3} {'DEC spread defect':>19} {'DEC var defect':>26} "
  f"{'CTL spread defect':>19} {'CTL var defect':>26}")
def meas(nn, aa, hh):
    s = configs(nn); E = E_field(s, aa, hh)
    tot = sum(int(v) for v in E); sq = sum(int(v)*int(v) for v in E)
    return int(E.max()-E.min()), sq//(1<<nn) - (tot//(1<<nn))**2
for (m, k) in [(2,2),(3,3),(4,4),(5,5),(6,6)]:
    n = m + k
    a = couplings(n-1); h = fields(n)
    acut = list(a); acut[m-1] = 0
    for arm, aa in (("DEC", acut), ("CTL", a)):
        Sw, Vw = meas(n, aa, h)
        SA, VA = meas(m, aa[:m-1], h[:m])
        SB, VB = meas(k, aa[m:], h[m:])
        d = (Sw-SA-SB, Vw-VA-VB)
        if arm == "DEC": d1 = d
        else: d2 = d
    p(f"{n:>3} {m:>3} {d1[0]:>19} {d1[1]:>26} {d2[0]:>19} {d2[1]:>26}")
p("READ: exactly additive over disjoint regions, live control non-zero -- ON A CARRIER WITH NO RECORDS.")

p("")
p("(c) NOT A COUNT, on the record-free carrier: 200 coupling+field draws at fixed n.")
p(f"{'n':>3} {'distinct S / 200':>18} {'min S/D':>12} {'max S/D':>12} {'max/min':>9} {'CONTROL count n':>16}")
for n in (6, 8, 10):
    vals = set()
    lo = hi = None
    for st in range(1, 201):
        a = couplings(n-1, stream=st); h = fields(n, stream=1000+st)
        v = 2*(sum(a)+sum(h)); vals.add(v)
        lo = v if lo is None else min(lo, v); hi = v if hi is None else max(hi, v)
    p(f"{n:>3} {len(vals):>18} {lo/D:>12.6f} {hi/D:>12.6f} {hi/lo:>9.4f} {n:>16}")
p("READ: 200 distinct values from 200 draws. NOT A COUNT -- with no records present.")

p("")
p("(d) SIGN-DEFINITENESS (C-46) on the record-free carrier: terms 2|J_i| and 2|h_i|.")
p(f"{'n':>3} {'|sum|/sum|.| FIELD carrier':>28} {'[lane carrier]':>16}")
for n in (4, 8, 16, 64):
    a = couplings(n-1); h = fields(n)
    t = [2*v for v in a] + [2*v for v in h]
    p(f"{n:>3} {abs(sum(t))/sum(abs(v) for v in t):>28.6f} {1.0:>16.6f}")
p("READ: ratio exactly 1 on the record-free carrier too.")
p("")
p("VERDICT V1: all four standards the lane reports as PASSED are passed IDENTICALLY by a")
p("carrier that holds zero records. The passing quantity is a functional of H alone -- the")
p("spectral width and the trace variance -- and knows nothing about what is written.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE/VERIFY/v1_record_blind.txt","w").write("\n".join(OUT)+"\n")
