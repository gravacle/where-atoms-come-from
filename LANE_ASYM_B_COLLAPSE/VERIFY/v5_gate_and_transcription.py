"""V5 -- (a) TRANSCRIPTION CHECK: does the reported finding match the lane's own .npy data?
        (b) GATE FRAGILITY THE OTHER WAY: add the POWER LAW to the model set and re-run the
            engine's extensivity gate on the KNOWN-EXACTLY-LINEAR chi control.
"""
import sys, numpy as np
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE"
sys.path.insert(0, LANE)
from fss_lib import fit_form, aicc, _wls

OUT = []
def say(s=""):
    print(s); OUT.append(s)

s2 = np.load(LANE+"/s2_data.npy")
s3 = np.load(LANE+"/s3_data.npy")

say("="*112)
say("V5a TRANSCRIPTION CHECK -- the numbers in the reported finding against the lane's .npy files")
say("="*112)
say()
# finding's RAW chi table rows: (n, N, SUM fixed, joint fixed, SUM grown)
CLAIM = [(4, 2, 0.992506, 0.992506, 0.992506), (6, 4, 1.214476, 1.443735, 1.984713),
         (8, 6, 0.772985, 1.403562, 2.926371), (10, 8, 0.644711, 1.679272, 3.918878),
         (12, 10, 0.577766, 1.787757, 4.911084), (14, 12, 0.512064, 1.731021, 5.852743),
         (22, 20, 0.400387, 1.987584, 9.771620), (32, 30, 0.320029, 2.049119, 14.631857),
         (42, 40, 0.280845, 2.164087, 19.542941)]
lut = {int(r[0]): r for r in s3}
say("      N |   claimed SUM   npy SUM      d  |  claimed joint  npy joint    d  | claimed grown  npy grown   d")
allok = True
for n, N, a, b, c in CLAIM:
    r = lut[N]
    d1, d2, d3 = abs(r[1]-a), abs(r[4]-b), abs(r[6]-c)
    ok = max(d1, d2, d3) < 1e-6; allok &= ok
    say("   %4d | %12.6f %10.6f %7.0e | %12.6f %10.6f %7.0e | %12.6f %10.6f %7.0e  %s"
        % (N, a, r[1], d1, b, r[4], d2, c, r[6], d3, "ok" if ok else "MISMATCH"))
say()
say("   chi transcription faithful: %s" % allok)
say()
# combinatorial claims
CL2 = [(4, 2, 0, 0.0, 1, 1.00, 4, 2, 1.00, 6, 4.00, 8),
       (8, 6, 0, 0.0, 7, 6.12, 12, 6, 1.00, 34, 15.54, 24),
       (14, 12, 0, 0.0, 30, 14.44, 24, 12, 1.00, 132, 33.57, 48),
       (24, 22, 0, 0.0, 111, 58.84, 44, 22, 1.00, 466, 126.60, 88),
       (40, 38, 0, 0.0, 343, 163.44, 76, 38, 1.00, 1410, 341.45, 152),
       (60, 58, 0, 0.0, 813, 367.16, 116, 58, 1.00, 3310, 755.91, 232)]
l2 = {int(r[0]): r for r in s2}
say("      n |  Ps  lam_s   Po    lam_o   W  |  Bps  Blam_s  Bpo   Blam_o   BW  | matches npy")
ok2 = True
for row in CL2:
    n = row[0]; r = l2[n]
    got = (r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11])
    want = row[2:]
    d = max(abs(g-w) for g, w in zip(got, want))
    good = d < 0.011; ok2 &= good
    say("   %4d | %3.0f %6.2f %5.0f %7.2f %4.0f | %4.0f %6.2f %6.0f %8.2f %5.0f | %s (maxdiff %.3f)"
        % (n, r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11],
           "ok" if good else "MISMATCH", d))
say()
say("   combinatorial transcription faithful: %s" % ok2)
say()

say("="*112)
say("V5b GATE FRAGILITY IN THE OTHER DIRECTION -- add POW = a N^b to the model set")
say("="*112)
say()
say("   The engine's model set is {LIN, LOG, SAT1, SAT2}.  A power law is not in it.  For any")
say("   series that is nearly a power law -- which a LINEAR series trivially is, with exponent 1 --")
say("   POW and LIN are near-degenerate, so the engine's stopping rule 'dAICc < 4 => CANNOT")
say("   DISTINGUISH, stop' would fire on its own known-exactly-linear control.")
say()

def fit_pow(N, Q, sig):
    X = np.c_[np.log(N), np.ones_like(N)]
    beta, _, _, _ = _wls(X, np.log(Q), np.maximum(sig/Q, 1e-12))
    pred = np.exp(beta[1])*N**beta[0]
    return dict(name="POW", p=2, beta=list(beta), chi2=float(np.sum(((Q-pred)/sig)**2)))

N = s3[:, 0]
for lab, Q, sg in (("CONTROL-LIN2  SUM chi_i, GROWN bath (KNOWN EXACTLY LINEAR)", s3[:, 6], s3[:, 7]),
                   ("CONTROL-SAT   chi_joint, FIXED bath (KNOWN BOUNDED <= 3 bits)", s3[:, 4], s3[:, 5]),
                   ("Q1            SUM chi_i, FIXED bath (probe)", s3[:, 1], s3[:, 2])):
    f4 = {nm: fit_form(nm, N, Q, sg) for nm in ("LIN", "LOG", "SAT1", "SAT2")}
    a4 = sorted((aicc(f["chi2"], len(N), f["p"]), nm) for nm, f in f4.items())
    f5 = dict(f4); f5["POW"] = fit_pow(N, Q, sg)
    a5 = sorted((aicc(f["chi2"], len(N), f["p"]), nm) for nm, f in f5.items())
    say("   %s" % lab)
    say("     model set {LIN,LOG,SAT1,SAT2}      best %-5s dAICc %6.2f   -> engine says %s"
        % (a4[0][1], a4[1][0]-a4[0][0], "CANNOT DISTINGUISH" if a4[1][0]-a4[0][0] < 4 else a4[0][1]))
    say("     model set + POW = a N^b            best %-5s dAICc %6.2f   -> engine says %s"
        % (a5[0][1], a5[1][0]-a5[0][0], "CANNOT DISTINGUISH" if a5[1][0]-a5[0][0] < 4 else a5[0][1]))
    if "POW" in [nm for _, nm in a5[:2]]:
        say("     POW fitted exponent b = %.3f" % f5["POW"]["beta"][0])
    say()
say("   READ: adding one obvious competitor form flips the KNOWN-EXACTLY-LINEAR control from")
say("   'GROWING (linear), EXT? YES' to 'CANNOT DISTINGUISH'.  The flip is benign in substance --")
say("   POW's exponent there is 0.996, i.e. POW IS the linear law -- but it shows the dAICc >= 4")
say("   gate is carrying almost no weight on the chi series.  The chi verdicts in this lane rest")
say("   on the RAW numbers and the EXACT cap, not on the model selection.")
say("="*112)
open(LANE+"/VERIFY/v5_gate_and_transcription.txt", "w").write("\n".join(OUT)+"\n")
