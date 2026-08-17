# W-08 / M4-REFUTE-1  leg F — ON WHICH CONNECTION WAS THE CORPUS'S RECURRENCE ACTUALLY MEASURED?
#
# W-07 sec1 (its headline, one sentence): "S1's published connection has finite order 4 in U(1),
# and EVERY RECURRENCE FIGURE IN THIS CORPUS WAS MEASURED THERE."
# W-07 sec4: "W-01 measures it (0.0247 at n=42, recurring to 0.99994); COR-H corrects the figure;
#             COR-E corrects it again ... On which connection was it measured?  ON S1's."
# M4-9 ratifies this: "Confirmed true: ... every recurrence figure in the corpus was measured at
# that point."  BOTH lanes assert it.  Neither checked the source lines.
#
# The corpus states its own test point, verbatim, in three places:
#   S3_THE_CROSSING_V001.md:421   "At the audit's own test point f = 2.0, c = 1.1,
#                                  p = (0.4, 0.15, 0.15, 0.15, 0.15)"  -> 0.024654 at k=42,
#                                                                         0.999941 at k=377
#   S3_THE_CROSSING_V001.md:565   record-breakers to k=200000 "of the frozen carrier at f=2.0,c=1.1"
#                                  -> COR-E's 0.999999981 at k=106123
#   S2_..._AUDIT_V001.md:440      "at f = 2.0, c = 1.1, p = (0.4,0.15,0.15,0.15,0.15)"
#                                  -> min 0.024654 at n=42, max 0.999941   (W-01's own figures)
#
# ISOLATION LEDGER.  Held fixed: carrier K1, the closed form Z_k = P0 e^{ik(c-f)} + PF e^{-ikf}
# + PC e^{ikc} (S2 Theorem A / S3 sec4.1, and W-07's own leg D uses it verbatim), the k-ranges the
# corpus used.  MOVED: the connection alone, between S1's published (f=pi, c=3pi/2) and the corpus's
# stated test point (f=2.0, c=1.1), at the corpus's own stated ready state.
# PRECISION: float64.  The exact-zero / exact-one claims are decided by the closed form, stated.
import numpy as np

def Z(f, c, p, K):
    P0 = p[0]; PF = p[1]+p[2]; PC = p[3]+p[4]
    k = np.arange(1, K+1)
    return np.abs(P0*np.exp(1j*k*(c-f)) + PF*np.exp(-1j*k*f) + PC*np.exp(1j*k*c)), k

print("== F1  REPRODUCE THE CORPUS'S OWN RECURRENCE FIGURES AT ITS OWN STATED TEST POINT ==")
p_test = [0.4, 0.15, 0.15, 0.15, 0.15]
d, k = Z(2.0, 1.1, p_test, 400)
i = d.argmin()
print(f"  f=2.0 c=1.1 p=(0.4,0.15,0.15,0.15,0.15):  min |Z_k| over k<=400  = {d[i]:.6f} at k = {k[i]}")
print( "                                            corpus records          0.024654 at k = 42")
d, k = Z(2.0, 1.1, p_test, 4000)
i = d.argmax()
print(f"                                            sup |Z_k| over k<=4000 = {d[i]:.6f} at k = {k[i]}")
print( "                                            corpus records          0.999941 at k = 377")
d, k = Z(2.0, 1.1, p_test, 200000)
i = d.argmax()
print(f"                                            sup |Z_k| over k<=2e5  = {d[i]:.9f} at k = {k[i]}")
print( "                                            COR-E records           0.999999981 at k = 106123")
print()

print("== F2  THE SAME FIGURES AT S1's PUBLISHED CONNECTION — WHERE W-07 SAYS THEY WERE MEASURED ==")
d, k = Z(np.pi, 3*np.pi/2, p_test, 400)
i = d.argmin(); j = d.argmax()
print(f"  f=pi c=3pi/2 (ord 4), same p:  min |Z_k| over k<=400 = {d[i]:.6f} at k = {k[i]}")
print(f"                                 sup |Z_k| over k<=400 = {d[j]:.6f} at k = {k[j]}")
print( "  0.024654 and 0.999941 DO NOT OCCUR at S1's published connection at any k.  On the order-4")
print( "  connection |Z_k| takes only FOUR values, exactly and periodically -- there is no k=42, no")
print( "  k=377, and no near-miss to five decimal places, because there are no near-misses at all:")
d4, k4 = Z(np.pi, 3*np.pi/2, p_test, 40)
print(f"    |Z_k| for k = 1..8 : {np.round(d4[:8],9)}")
print(f"    distinct values of |Z_k| over k<=4000: {sorted(set(np.round(Z(np.pi,3*np.pi/2,p_test,4000)[0],12)))}")
print()

print("== F3  W-07's OWN sec4 TABLE ALREADY CONTRADICTS ITS OWN sec1/sec4 PROSE ==")
print("  W-07 sec4 lists 'S3/S4 HEADLINE f=2.0,c=1.1' as a SEPARATE ROW from 'S1 PUBLISHED', and")
print("  scores it 0 cells > 1-1e-12 (APPROACHED).  Recomputed here at K=1e6, W-07's own grid:")
for tag, (f, c) in [("S1 PUBLISHED  f=pi, c=3pi/2", (np.pi, 3*np.pi/2)),
                    ("S3/S4 HEADLINE f=2.0, c=1.1", (2.0, 1.1))]:
    for pnm, pp in [("published p=(1/2,0,0,1/4,1/4)", [0.5,0,0,0.25,0.25]),
                    ("corpus test p=(0.4,.15,.15,.15,.15)", p_test)]:
        d, k = Z(f, c, pp, 10**6)
        print(f"    {tag:<28} {pnm:<36} max|Z_k| = {d.max():.15f}  cells>1-1e-12 : {int((d>1-1e-12).sum())}")
print()
print("  CONCLUSION.  The corpus's recurrence figures that the REGISTER cites — W-01's 0.0247/0.99994,")
print("  S3 sec4.1's 0.024654 at k=42 and 0.999941 at k=377, COR-E's 0.999999981 at k=106123, and")
print("  S2's 0.994373 (S2 build :448-455, a CONTINUOUS-time overlap at f=3.663319822, c=2.194049746;")
print("  the S2 audit's replacement figure :440 uses f=2.0, c=1.1) — NONE were measured at S1's")
print("  published connection.  The k-indexed ones were all at f=2.0, c=1.1, the")
print("  EXACTLY RESONANT connection of the erratum against W-02, whose ord is INFINITE and at which")
print("  W-07's own leg D scores APPROACHED, 0 cells.  W-06's '1000 of 4000' is the ONLY recurrence")
print("  figure in the corpus measured at S1's order-4 point.")
print()
print("  W-07 sec1's headline sentence is FALSE as written, and W-07 sec4's rhetorical answer")
print("  ('On which connection was it measured? On S1's.') is FALSE.  M4-9 ratified both.")
print("  This does not touch W-07's degeneracy FACT (ord = 4 at S1's published connection) nor its")
print("  second, correct observation that the corpus has TWO distinguished connections and both are")
print("  arithmetically degenerate.  It refutes the UNIVERSAL QUANTIFIER W-07 built its sec1")
print("  headline on, and it is the one claim in the exchange that neither lane tested.")
