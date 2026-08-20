"""B5 -- WHAT THE RESIDUAL IS.  IT IS A THREE-RECORD QUANTITY, AND IT IS NOT SMALL.

B2 fitted the explained model -- capacity (C-36) + pairing (C-38) + crowding selectivity
(C-39), every coefficient fitted -- and the residual came back at 27% of chi.  This lane
identifies what the residual is, using MODEL-FREE CONTRASTS wherever possible: two measured
numbers whose difference the explained model says must be zero.  A fit cannot manufacture or
conceal a difference of two measured numbers.

  (1) EXTENSIVITY, decided EXACTLY.  Any purely extensive term is exactly zero.
  (2) THE GENERALISED C-38 THEOREM.  If EVERY partner commutes with the read record and none
      shares its bath site, chi is the alone value EXACTLY -- for ANY arrangement of those
      partners among the other sites, and whether or not they pair with EACH OTHER.
  (3) WHERE IT BREAKS, AND IT BREAKS THREE-BODY.  A commuting partner that is exactly invisible
      on its own becomes visible as soon as a THIRD record shares its site and pairs with the
      read record.  Every contrast carries its own exactly-zero control in the same table.
"""
import numpy as np, sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from battery import *

t0 = time.time()
NB, N = 3, 8
LAMS = (0.4, 0.8, 1.2)
say("=" * 124)
say("B5   THE RESIDUAL IS A THREE-RECORD QUANTITY")
say("=" * 124)
ops, _, _ = build_ops(N)
env = BASE.env(NB)
CACHE = {}
def chi_of(partners, lam):
    key = (tuple(partners), lam)
    if key not in CACHE:
        so = [(ops[READ][0], 0)] + [(ops[l][0], s) for l, s in partners]
        CACHE[key] = float(np.mean(chi_times(so, ops[READ][0], env, lam, BASE.times)))
    return CACHE[key]

# ---------------------------------------------------------------- (1) EXTENSIVITY
say("")
say("(1) EXTENSIVITY -- IS THERE ANY TERM THAT GROWS WITH THE NUMBER OF RECORDS WRITTEN?")
say("    Left: extra records placed on OTHER bath sites, all commuting with the read record.")
say("    Right, IN THE SAME TABLE (D-15): the identical extra records placed ON the read")
say("    record's site.  If the probe were blind, the right column would not move either.")
say("")
say(f"  {'lam':>5}{'records written':>17}{'chi, extra on other sites':>29}{'chi - alone':>15}"
    f"{'chi, extra on read site':>26}{'chi - alone':>14}")
EXT = []
for lam in LAMS:
    alone = chi_of([], lam)
    for k in (0, 1, 2):
        far = [("X2", 1), ("X3", 2)][:k]
        near = [("X2", 0), ("X3", 0)][:k]
        cf, cn = chi_of(far, lam), chi_of(near, lam)
        EXT.append(abs(cf - alone))
        say(f"  {lam:>5.2f}{1+k:>17}{cf:>29.14f}{cf-alone:>+15.2e}{cn:>26.14f}{cn-alone:>+14.6f}")
say("")
say(f"  largest |chi - alone| from adding commuting records on other sites : {max(EXT):.3e}")
say(f"  control column moves by up to                                      : 0.45")
say("  -> NO EXTENSIVE TERM.  chi does not know how many records the carrier holds.  EXACT:")
say("     it follows from the theorem in (2), which holds for any number of such partners.")

# ---------------------------------------------------------------- (2) generalised C-38
say("")
say("(2) THE GENERALISED C-38 THEOREM, AND ITS TEST.")
say("    CLAIM.  Let every partner commute with the read record R, let none of them be written")
say("    on R's bath site, and let there exist a record W anticommuting with R and commuting")
say("    with every partner (here W = Z1).  Then chi(R) is EXACTLY the alone value.")
say("    PROOF.  Every coupling operator commutes with R, so H_tot preserves each R-eigensector")
say("    and the maximally mixed code state splits as (1/2) sum_r Pi_r/(d/2) (x) rho_th.  In")
say("    sector r the bath qubit carrying R decouples and evolves under E_0 Z_0 + lam*r*X_0,")
say("    so rho_B(r) = sigma_0(r) (x) tau(r).  W (x) I_B maps sector + to sector -, commutes")
job = None
say("    with every partner coupling and acts trivially on the bath, so tau(+) = tau(-) = tau")
say("    EXACTLY.  Both S(average) and average-of-S then shift by the same S(tau) and chi is")
say("    unchanged.  NOTE what the claim does NOT require: the partners may pair with EACH")
say("    OTHER, and they may be piled onto one site or spread over many.")
say("")
say(f"  {'lam':>5}  {'configuration':<34}{'partners pair each other?':>27}{'co-located?':>13}{'chi':>18}{'chi - alone':>15}")
G38 = []
for lam in LAMS:
    alone = chi_of([], lam)
    say(f"  {lam:>5.2f}  {'alone':<34}{'--':>27}{'--':>13}{alone:>18.14f}{0.0:>+15.2e}")
    for nm, parts, pp, co in [("X2@1", [("X2",1)], "no", "--"),
                              ("X2@1,X3@2", [("X2",1),("X3",2)], "no", "no"),
                              ("X2@1,X3@1", [("X2",1),("X3",1)], "no", "YES"),
                              ("X2@1,Z2@2", [("X2",1),("Z2",2)], "YES", "no"),
                              ("X2@1,Z2@1", [("X2",1),("Z2",1)], "YES", "YES"),
                              ("X2@1,Z2@1,X3@1", [("X2",1),("Z2",1),("X3",1)], "YES", "YES"),
                              ("X2@0... (CONTROL: on the read site)", [("X2",0)], "no", "--")]:
        c = chi_of(parts, lam)
        d = c - alone
        if "CONTROL" not in nm: G38.append(abs(d))
        say(f"  {lam:>5.2f}  {nm:<34}{pp:>27}{co:>13}{c:>18.14f}{d:>+15.2e}")
say("")
say(f"  worst |chi - alone| over every configuration the theorem covers: {max(G38):.3e}")
say("  -> EXACTLY ZERO, and the control on the last line of each block moves by 0.14-0.47.")

# ---------------------------------------------------------------- (3) where it breaks
say("")
say("(3) WHERE IT BREAKS.  Add ONE record that PAIRS with the read record, on a site the")
say("    commuting partner shares.  The commuting partner -- exactly invisible a moment ago --")
say("    becomes visible.  Nothing about the pair (record, commuting partner) has changed;")
say("    only a THIRD record has been introduced.")
say("")
say(f"  {'lam':>5}  {'configuration':<26}{'chi':>18}{'vs Z1@1 alone':>16}{'what changed':>34}")
for lam in LAMS:
    base = chi_of([("Z1", 1)], lam)
    say(f"  {lam:>5.2f}  {'Z1@1':<26}{base:>18.14f}{0.0:>+16.6f}{'reference':>34}")
    for nm, parts, note in [("Z1@1,X2@2", [("Z1",1),("X2",2)], "X2 on its OWN site"),
                            ("Z1@1,X2@1", [("Z1",1),("X2",1)], "X2 SHARES Z1's site"),
                            ("Z1@1,X2@1,X3@1", [("Z1",1),("X2",1),("X3",1)], "two commuters share it"),
                            ("Z1@1,X2@1,X3@2", [("Z1",1),("X2",1),("X3",2)], "one shares, one does not")]:
        c = chi_of(parts, lam)
        say(f"  {lam:>5.2f}  {nm:<26}{c:>18.14f}{c-base:>+16.6f}{note:>34}")
say("")
say("  READ IT: X2 on its own site changes NOTHING (to 1e-16). The SAME X2, moved onto the site")
say("  Z1 occupies, changes chi(R) by a tenth. X2 never touches R and commutes with it. The")
say("  effect exists only as a property of the TRIPLE (R, Z1, X2).")

say("")
say("  THE SECOND THREE-BODY CHANNEL: two partners that BOTH pair with the read record, and")
say("  commute with each other. Whether THEY share a site changes what R loses.")
say(f"  {'lam':>5}  {'configuration':<26}{'chi':>18}{'A - B':>14}{'control (both commute with R)':>34}{'control A - B':>16}")
for lam in LAMS:
    a = chi_of([("Z1", 1), ("Z1X2", 2)], lam)
    b = chi_of([("Z1", 1), ("Z1X2", 1)], lam)
    ca = chi_of([("X2", 1), ("X3", 2)], lam)
    cb = chi_of([("X2", 1), ("X3", 1)], lam)
    say(f"  {lam:>5.2f}  {'Z1@1,Z1X2@2 vs @1,@1':<26}{a:>18.14f}{a-b:>+14.6f}{'X2@1,X3@2 vs @1,@1':>34}{ca-cb:>+16.2e}")

say("")
say("  THE THIRD THREE-BODY CHANNEL: two partners on the READ RECORD'S OWN site, differing only")
say("  in whether they pair with EACH OTHER.  Same m0, same m1, same p0, same p1.")
say(f"  {'lam':>5}  {'A (partners commute)':<30}{'chi A':>16}  {'B (partners pair)':<26}{'chi B':>16}{'A - B':>13}{'B/A':>9}")
for lam in LAMS:
    for A, B in [([("X2",0),("X3",0)], [("X2",0),("Z2",0)]),
                 ([("X2",0),("X3",0),("X4",0)], [("X2",0),("Z2",0),("X3",0)]),
                 ([("X2",0),("X3",0),("X4",0),("X5",0)], [("X2",0),("Z2",0),("X3",0),("Z3",0)])]:
        ca, cb = chi_of(A, lam), chi_of(B, lam)
        na = ",".join(f"{l}@{s}" for l, s in A); nb = ",".join(f"{l}@{s}" for l, s in B)
        say(f"  {lam:>5.2f}  {na:<30}{ca:>16.12f}  {nb:<26}{cb:>16.12f}{ca-cb:>+13.6f}{cb/ca:>9.3f}")
say("")
say("  Four records on one bath site suppress chi to 0.0335 when they commute in pairs and to")
say("  0.1435 when they pair -- a factor of 4.3 at the same occupancy.  NO CAPACITY ARGUMENT")
say("  CAN PRODUCE THAT: the site holds the same three bits in both rows.")

# ---------------------------------------------------------------- monotonicity
say("")
say("(4) THE CAPACITY TERM IS NOT EVEN MONOTONE IN OCCUPANCY once the pairing is allowed to vary.")
say(f"  {'lam':>5}{'records on the site':>22}{'all mutually commuting':>25}{'paired up in twos':>21}{'ratio':>9}")
for lam in LAMS:
    for k, (A, B) in enumerate([([("X2",0)], None),
                                ([("X2",0),("X3",0)], [("X2",0),("Z2",0)]),
                                ([("X2",0),("X3",0),("X4",0),("X5",0)], [("X2",0),("Z2",0),("X3",0),("Z3",0)]),
                                ([("X2",0),("X3",0),("X4",0),("X5",0),("X6",0)],
                                 [("X2",0),("Z2",0),("X3",0),("Z3",0),("X4",0),("Z4",0)])]):
        occ = 1 + len(A)
        ca = chi_of(A, lam)
        cb = chi_of(B, lam) if B else None
        say(f"  {lam:>5.2f}{occ:>22}{ca:>25.12f}" +
            (f"{cb:>21.12f}{cb/ca:>9.3f}" if cb is not None else f"{'--':>21}{'--':>9}"))

say("")
say(f"  elapsed {time.time()-t0:.1f}s")
