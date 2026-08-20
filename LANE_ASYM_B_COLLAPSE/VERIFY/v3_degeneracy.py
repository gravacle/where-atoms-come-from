"""V3 -- THE VENUE-SCALE VARIATION THE LANE DID NOT RUN (D-17, applied to the COUPLING PATTERN).

S3's fixed-bath venue couples EVERY record to its bath site with the SAME strength lam.  With
equal couplings the bath does not merely SHARE a site between records (C-36) -- it is literally
BLIND to which record is which: records in the same group enter H_B only through the integer
c_j = SUM s_i.  That degeneracy is what the "counting compression" exploits, and it is exactly
the kind of structure that manufactures a 1/sqrt(m) dilution.

The lane's D-17 sweep varied lam, nq, beta and the time window -- the SCALE of the coupling.
It never varied the PATTERN.  This script does.

  EQUAL     lam_i = lam for all i           (reproduces S3 -- the control that must match)
  SPREAD    lam_i uniform in [0.4, 1.2]     (records distinguishable, same mean coupling)
  GEOMETRIC lam_i = lam * 1.35^(i within group)  (maximally non-degenerate)

Bath: nq = 3 qubits, H_B = SUM_j e_j Z_j, coupling operator per site = X_j.  H_B + coupling is
a SUM OF SINGLE-QUBIT TERMS, so the bath factorises EXACTLY and chi(R_i : B) = chi(R_i : qubit
j(i)).  Verified against the lane's own chi_lib in the EQUAL case.

TIME-AVERAGED over the same 25 times in [1,13].
"""
import sys, itertools, numpy as np
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE"
sys.path.insert(0, LANE)
from chi_lib import total_chi_fixed, TIMES

OUT = []
def say(s=""):
    print(s); OUT.append(s)

Z = np.array([[1, 0], [0, -1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)

def vN(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
    return float(-(e*np.log2(e)).sum())

def qubit_states(e, beta, lams, t):
    """rho(g) for every sign pattern of the group's records; returns dict pattern -> rho."""
    h = e*Z
    w0, V0 = np.linalg.eigh(h)
    p = np.exp(-beta*w0); p /= p.sum()
    r0 = (V0*p) @ V0.conj().T
    m = len(lams)
    out = {}
    for s in itertools.product((1, -1), repeat=m):
        g = float(np.dot(lams, s))
        H = h + g*X
        w, U = np.linalg.eigh(H); ph = np.exp(-1j*w*t)
        Uc = U.conj().T @ r0 @ U
        out[s] = U @ (ph[:, None]*Uc*ph.conj()[None, :]) @ U.conj().T
    return out

def group_chis(e, beta, lams, t):
    """per-record chi and the group's joint chi, on that group's single bath qubit."""
    m = len(lams)
    if m == 0: return [], 0.0
    st = qubit_states(e, beta, lams, t)
    avg = sum(st.values())/2**m
    cond = sum(vN(r) for r in st.values())/2**m
    cj = max(vN(avg)-cond, 0.0)
    per = []
    for i in range(m):
        P = sum(r for s, r in st.items() if s[i] > 0)/2**(m-1)
        M = sum(r for s, r in st.items() if s[i] < 0)/2**(m-1)
        per.append(max(vN(0.5*(P+M)) - 0.5*(vN(P)+vN(M)), 0.0))
    return per, cj

def totals(k, pattern, nq=3, lam=0.8, beta=2.0, times=TIMES, seed=11):
    en = [(1.0, 1.4, 0.7)[j % 3] for j in range(nq)]
    rng = np.random.default_rng(seed)
    groups = [[i for i in range(k) if i % nq == j] for j in range(nq)]
    LAM = []
    for j in range(nq):
        m = len(groups[j])
        if pattern == "EQUAL":      LAM.append(np.full(m, lam))
        elif pattern == "SPREAD":   LAM.append(rng.uniform(0.4, 1.2, m))
        elif pattern == "GEOMETRIC":LAM.append(lam*1.35**np.arange(m)/np.mean(1.35**np.arange(m)))
        else: raise ValueError(pattern)
    tot, jnt = [], []
    for t in times:
        s = 0.0; cj = 0.0
        for j in range(nq):
            per, c = group_chis(en[j], beta, LAM[j], t)
            s += sum(per); cj += c
        tot.append(s); jnt.append(cj)
    tot, jnt = np.array(tot), np.array(jnt)
    return tot.mean(), tot.std(ddof=1)/np.sqrt(len(tot)), jnt.mean(), jnt.std(ddof=1)/np.sqrt(len(jnt))

say("="*116)
say("V3  D-17 ON THE COUPLING **PATTERN**, not just its scale.  Fixed 3-qubit bath, nq = 3.")
say("="*116)
say()

# ---------------- control: my factorised engine must reproduce the lane's chi_lib exactly
say("-"*116)
say("CONTROL FIRST (D-15).  My factorised engine vs the lane's chi_lib, EQUAL couplings.")
say("   if these disagree nothing below means anything.")
say()
say("      k   lane SUM chi_i   mine SUM chi_i     |diff|   |  lane chi_joint   mine chi_joint    |diff|")
ok = True
for k in (2, 4, 6, 9, 12):
    tm, tse, tsd, jm, jse = total_chi_fixed(k)
    a, ae, b, be = totals(k, "EQUAL")
    d1 = abs(tm-a); d2 = abs(jm-b); ok &= (d1 < 1e-9 and d2 < 1e-9)
    say("   %4d %16.9f %16.9f %10.2e   | %16.9f %16.9f %10.2e" % (k, tm, a, d1, jm, b, d2))
say()
say("   CONTROL %s" % ("PASS -- same numbers, so the pattern sweep below is on the same footing"
                       if ok else "FAIL -- STOP, conclude nothing"))
say()
if not ok:
    open(LANE+"/VERIFY/v3_degeneracy.txt", "w").write("\n".join(OUT)+"\n"); sys.exit(1)

# ---------------- the sweep
KS = [2, 4, 6, 8, 10, 12, 15, 18, 21, 24, 27, 30]
say("-"*116)
say("TABLE.  SUM_i chi(R_i : fixed 3q bath), time-averaged over 25 t in [1,13], nq = 3.")
say("        EQUAL is the lane's venue.  SPREAD and GEOMETRIC break the record-record degeneracy")
say("        while keeping the bath, the temperature and the mean coupling the same.")
say()
say("      N=k |  EQUAL (lane)   |   SPREAD        |  GEOMETRIC      ||  chi_joint EQUAL  SPREAD   GEOM   (cap = 3)")
res = {p: [] for p in ("EQUAL", "SPREAD", "GEOMETRIC")}
resj = {p: [] for p in ("EQUAL", "SPREAD", "GEOMETRIC")}
for k in KS:
    row = {}
    for p in ("EQUAL", "SPREAD", "GEOMETRIC"):
        a, ae, b, be = totals(k, p)
        res[p].append(a); resj[p].append(b); row[p] = (a, ae, b)
    say("   %6d | %8.5f+-%.4f | %8.5f+-%.4f | %8.5f+-%.4f || %8.5f %8.5f %8.5f"
        % (k, row["EQUAL"][0], row["EQUAL"][1], row["SPREAD"][0], row["SPREAD"][1],
           row["GEOMETRIC"][0], row["GEOMETRIC"][1],
           row["EQUAL"][2], row["SPREAD"][2], row["GEOMETRIC"][2]))
say()
say("   RATIOS (gravity's requirement (a) needs Q(2N)/Q(N) -> 2):")
for p in ("EQUAL", "SPREAD", "GEOMETRIC"):
    v = dict(zip(KS, res[p]))
    dbl = [(a, v[2*a]/v[a]) for a in KS if 2*a in v]
    say("     %-10s  %s" % (p, "  ".join("Q(%d)/Q(%d)=%.3f" % (2*a, a, r) for a, r in dbl)))
say()
say("   MONOTONICITY: does SUM_i chi_i still DECAY once the degeneracy is broken?")
for p in ("EQUAL", "SPREAD", "GEOMETRIC"):
    v = np.array(res[p])
    say("     %-10s  first %.5f  max %.5f (at N=%d)  last %.5f   decaying over N>=6? %s"
        % (p, v[0], v.max(), KS[int(v.argmax())], v[-1],
           bool(np.all(np.diff(v[2:]) < 0))))
say()
say("   CAP CHECK (exact chain SUM_i chi_i <= chi_joint <= S(rho_B) <= nq = 3 bits):")
for p in ("EQUAL", "SPREAD", "GEOMETRIC"):
    v = np.array(res[p]); j = np.array(resj[p])
    say("     %-10s  max SUM = %.5f   max joint = %.5f   SUM<=joint everywhere? %s   joint<=3? %s"
        % (p, v.max(), j.max(), bool(np.all(v <= j+1e-9)), bool(np.all(j <= 3+1e-9))))
say()
say("   POWER-LAW EXPONENT of SUM_i chi_i over N >= 6, per pattern:")
for p in ("EQUAL", "SPREAD", "GEOMETRIC"):
    v = np.array(res[p]); n = np.array(KS, float)
    m = n >= 6
    A = np.c_[np.log(n[m]), np.ones(m.sum())]
    beta, *_ = np.linalg.lstsq(A, np.log(v[m]), rcond=None)
    pred = A@beta; r = np.log(v[m])-pred
    se = float(np.sqrt(np.sum(r**2)/max(len(r)-2, 1)/np.sum((np.log(n[m])-np.log(n[m]).mean())**2)))
    say("     %-10s  exponent %+.3f +- %.3f   (rms log-residual %.4f)" % (p, beta[0], se, float(np.sqrt(np.mean(r**2)))))
say("-"*116)
open(LANE+"/VERIFY/v3_degeneracy.txt", "w").write("\n".join(OUT)+"\n")
