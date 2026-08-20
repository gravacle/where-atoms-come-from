"""V3 -- ADVERSARIAL CHECK OF THE 'COMPACT SUPPORT' / 'INFLUENCE RADIUS = BATH COUPLING RANGE'
CLAIM (lane B7/B8).

THE ORDINARY EXPLANATION UNDER TEST.  The lane's ring bath is
        H_B = E * sum_j Z_j  +  J * sum_j Z_j Z_{j+r}
which is COMPLETELY DIAGONAL.  The only non-diagonal operators in the whole problem are the
two coupling terms lam*R (x) X_0 and lam*P (x) X_s.  Therefore Z_j is EXACTLY CONSERVED on
every bath qubit j not in {0, s}.  Conditioning on those conserved classical values, the bath
splits into the block containing qubit 0 and the block containing qubit s, and when qubit s
carries no H_B term that touches qubit 0 the two blocks are exactly independent: the partner's
entropy contribution is identical in both R-sectors and cancels in chi.

That predicts something much narrower than "compact support with radius = the coupling range":
it predicts the zero set is exactly

        {s : H_B contains NO term coupling qubit 0 to qubit s}

i.e. ADJACENCY IN H_B, not DISTANCE, and the exactness of the zeros should be destroyed by
making the bath coupling NON-DIAGONAL, which breaks the spectator-Z conservation.

TWO DISCRIMINATING VENUES, neither of which the lane ran:
  (B) A SINGLE-BOND bath: the only bond is 0--3, the MAXIMUM ring distance.  "Radius = coupling
      range" predicts nothing at distance 3.  Adjacency predicts the effect appears at s = 3 and
      is exactly zero at s = 1, 2, 4, 5.
  (C) A NON-DIAGONAL (XX) nearest-neighbour bath.  "Radius = coupling range" predicts exact
      zeros beyond distance 1.  Adjacency-through-Z-conservation predicts the exact zeros are
      GONE, because nothing is conserved any more.

POSITIVE CONTROL (D-15): venue (A) reproduces the lane's own diagonal range-1 ring, whose
zeros must come back, and the read site s = 0 entry, which must be large in every venue.

METHOD.  The partner commutes with R, so H_tot preserves each R-eigensector and the whole
computation reduces EXACTLY to a bath-only calculation: in sector r the bath evolves under
H_B + lam*r*X_0 + lam*p*X_s with p = +-1 equally weighted.  This reduction is validated
against the lane's own full pipeline (common.chi_times) in the first table.
"""
import numpy as np, sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from common import *                      # noqa
from battery import build_ops, READ       # noqa

t0 = time.time()
L, LAM, BETA, E = 6, 0.8, 2.0, 1.0
TIMES = np.linspace(1.0, 13.0, 25)

I2 = np.eye(2); Xb = np.array([[0, 1], [1, 0]], complex); Zb = np.array([[1, 0], [0, -1]], complex)
def bop(j, P, L=L):
    M = np.array([[1]], complex)
    for k in range(L): M = np.kron(M, P if k == j else I2)
    return M
XS = [bop(j, Xb) for j in range(L)]
ZS = [bop(j, Zb) for j in range(L)]

def HB_ring(J, rng):
    H = E * sum(ZS)
    for r in range(1, rng + 1):
        H = H + J * sum(ZS[j] @ ZS[(j + r) % L] for j in range(L))
    return H
def HB_singlebond(J, a=0, b=3):
    return E * sum(ZS) + J * (ZS[a] @ ZS[b])
def HB_XX(J):
    return E * sum(ZS) + J * sum(XS[j] @ XS[(j + 1) % L] for j in range(L))

def thermal(HB):
    w, V = np.linalg.eigh(HB)
    p = np.exp(-BETA * w); p /= p.sum()
    return (V * p) @ V.conj().T

def _vn(rho):
    e = np.linalg.eigvalsh(rho); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

def chi_bath(HB, partner_site, lam=LAM, times=TIMES):
    """EXACT reduction for a partner COMMUTING with the read record.
       partner_site = None -> no partner."""
    rho0 = thermal(HB)
    branches = {}
    for r in (+1, -1):
        ps = [(+1, 0.5), (-1, 0.5)] if partner_site is not None else [(0, 1.0)]
        bl = []
        for p, wgt in ps:
            H = HB + lam * r * XS[0]
            if partner_site is not None: H = H + lam * p * XS[partner_site]
            w, V = np.linalg.eigh(H)
            bl.append((wgt, w, V))
        branches[r] = bl
    acc = []
    for t in times:
        rr = {}
        for r in (+1, -1):
            tot = 0
            for wgt, w, V in branches[r]:
                U = (V * np.exp(-1j * w * t)[None, :]) @ V.conj().T
                tot = tot + wgt * (U @ rho0 @ U.conj().T)
            rr[r] = tot
        av = 0.5 * (rr[+1] + rr[-1])
        acc.append(max(_vn(av) - 0.5 * (_vn(rr[+1]) + _vn(rr[-1])), 0.0))
    return float(np.mean(acc))

say("=" * 126)
say("V3   IS IT AN INFLUENCE RADIUS SET BY THE BATH'S COUPLING RANGE, OR JUST ADJACENCY IN H_B?")
say("=" * 126)

# ------------------------------------------------------- validation of the reduction (D-15)
say("")
say("  VALIDATION OF THE BATH-ONLY REDUCTION against the lane's own full pipeline (n = 6, ring,")
say("  range 1, J = 0.6).  If these disagree, nothing below means anything.")
ops, _, _ = build_ops(6)
R = ops[READ][0]
class _B:
    def __init__(s, HB):
        s.nq, s.dim, s.beta = L, 2 ** L, BETA
        s.HB = HB; s.site = XS; s.probe = sum(XS)
    def thermal(s): return thermal(s.HB)
BR = _B(HB_ring(0.6, 1))
def chi_full(parts, bath, lam=LAM):
    so = [(R, 0)] + [(ops[l][0], s) for l, s in parts]
    return float(np.mean(chi_times(so, R, bath, lam, TIMES)))
say(f"  {'configuration':<22}{'full pipeline':>20}{'bath-only reduction':>24}{'|diff|':>12}")
for lab, parts, ps in [("alone", [], None), ("X2 at s=1", [("X2", 1)], 1),
                       ("X2 at s=2", [("X2", 2)], 2), ("X2 at s=3", [("X2", 3)], 3)]:
    a = chi_full(parts, BR); b = chi_bath(BR.HB, ps)
    say(f"  {lab:<22}{a:>20.14f}{b:>24.14f}{abs(a-b):>12.2e}")

# ------------------------------------------------------- the three venues
def sweep(name, HB, note=""):
    al = chi_bath(HB, None)
    vals = [chi_bath(HB, s) - al for s in range(1, L)]
    say(f"  {name:<34}{al:>13.9f}" + "".join(f"{v:>+15.6e}" for v in vals) + ("   " + note if note else ""))
    return vals

def dist(a, b): return min((a - b) % L, (b - a) % L)
say("")
say("  chi(with a COMMUTING partner at ring site s) - chi(alone).  Read record on site 0.")
say(f"  {'bath H_B':<34}{'chi(alone)':>13}" + "".join(f"{'s='+str(s)+' d='+str(dist(0,s)):>15}" for s in range(1, L)))
A = sweep("(A) diagonal ZZ, range 1, J=0.6", HB_ring(0.6, 1), "POSITIVE CONTROL: reproduces the lane")
B1 = sweep("(B) diagonal, ONE bond 0--3, J=0.6", HB_singlebond(0.6), "bond spans the MAXIMUM distance")
B2 = sweep("(B) diagonal, ONE bond 0--3, J=0.2", HB_singlebond(0.2))
B3 = sweep("(B) diagonal, ONE bond 0--2, J=0.6", HB_singlebond(0.6, 0, 2))
C1 = sweep("(C) NON-diagonal XX, range 1, J=0.6", HB_XX(0.6), "spectator Z no longer conserved")
C2 = sweep("(C) NON-diagonal XX, range 1, J=0.2", HB_XX(0.2))
C3 = sweep("(C) NON-diagonal XX, range 1, J=0.05", HB_XX(0.05))

say("")
say("  THE FLOOR IN EACH ROW, from the venue's own exact reflection symmetry s <-> L-s where it")
say("  exists (rows A, C: reflection about site 0 is a symmetry; row B 0--3 is also reflection")
say("  symmetric).  Reported beside the largest |effect| at the sites the lane calls OUT OF RANGE.")
say(f"  {'bath H_B':<34}{'floor':>12}{'|effect| d=1 (s=1,5)':>22}{'|effect| d=2 (s=2,4)':>22}{'|effect| d=3 (s=3)':>20}")
for nm, v, sym in [("(A) ZZ range1 J=0.6", A, True), ("(B) one bond 0--3 J=0.6", B1, True),
              ("(B) one bond 0--3 J=0.2", B2, True), ("(B) one bond 0--2 J=0.6", B3, False),
              ("(C) XX range1 J=0.6", C1, True), ("(C) XX range1 J=0.2", C2, True),
              ("(C) XX range1 J=0.05", C3, True)]:
    fl = f"{max(abs(v[0]-v[4]), abs(v[1]-v[3])):.2e}" if sym else "n/a (asym)"
    say(f"  {nm:<34}{fl:>12}{max(abs(v[0]),abs(v[4])):>22.4e}{max(abs(v[1]),abs(v[3])):>22.4e}{abs(v[2]):>20.4e}")
say("  NOTE: the 0--2 bond venue is NOT reflection-symmetric about site 0, so it has no")
say("  self-floor column; its own zeros (1e-16 at s=1,3,4,5) are read against the other rows.")

say("")
say(f"  elapsed {time.time()-t0:.1f}s")
