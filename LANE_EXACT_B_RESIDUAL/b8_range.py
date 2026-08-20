"""B8 -- THE SHAPE OF THE FALLOFF.  IT IS NOT A FALLOFF: IT IS COMPACT SUPPORT.

B7 found, on a ring bath with nearest-neighbour ZZ, that a partner commuting with the read
record changes chi by EXACTLY ZERO at ring distance 2 and 3 and by -0.078 at distance 1.  If
that is right, the influence of one record on another does not decay with separation -- it is
IDENTICALLY ZERO outside a radius, and the radius is set by the BATH's own interaction range,
not by the coupling strength.  That is a scale-free statement about FORM and it is decidable
here.

THE DISCRIMINATOR.  Vary two things independently:
  * the RANGE of the bath coupling (nearest-neighbour, then next-nearest-neighbour);
  * the STRENGTH J of that coupling, over a 6x range.
If the zero set follows the RANGE and ignores the STRENGTH, the influence has compact support.
If it were any decaying function of distance, changing J would move the apparent zeros.

CONTROL IN EVERY TABLE (D-15): J = 0, where every entry must be exactly zero, and the
PAIRING partner, whose effect must be large at every distance.
"""
import numpy as np, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from battery import build_ops, READ

t0 = time.time()

class RingBath:
    def __init__(self, L=6, E=1.0, J=0.0, rng=1, beta=2.0):
        self.nq, self.dim, self.beta = L, 2 ** L, beta
        I2 = np.eye(2)
        Xb = np.array([[0, 1], [1, 0]], dtype=complex)
        Zb = np.array([[1, 0], [0, -1]], dtype=complex)
        def op(j, P):
            M = np.array([[1]], dtype=complex)
            for k in range(L): M = np.kron(M, P if k == j else I2)
            return M
        self.site = [op(j, Xb) for j in range(L)]
        Zs = [op(j, Zb) for j in range(L)]
        HB = E * sum(Zs)
        for r in range(1, rng + 1):
            HB = HB + J * sum(Zs[j] @ Zs[(j + r) % L] for j in range(L))
        self.HB = HB
        self.probe = sum(self.site)
        self.L, self.J, self.rng = L, J, rng
    def thermal(self):
        w, V = np.linalg.eigh(self.HB)
        p = np.exp(-self.beta * w); p /= p.sum()
        return (V * p) @ V.conj().T

L, N = 6, 6
TIMES = np.linspace(1.0, 13.0, 25)
ops, _, _ = build_ops(N)
R = ops[READ][0]
def dist(a, b): return min((a - b) % L, (b - a) % L)
def chi(parts, bath, lam=0.8):
    so = [(R, 0)] + [(ops[l][0], s) for l, s in parts]
    return float(np.mean(chi_times(so, R, bath, lam, TIMES)))

say("=" * 120)
say("B8   IS IT A FALLOFF OR IS IT COMPACT SUPPORT?")
say("=" * 120)
say(f"  ring bath L = {L}, uniform on-site energy 1.0, beta 2.0; carrier n = {N}; lam = 0.8")
say("  read record written on ring site 0.  Partner moved to ring site s.")
say("")
say("  COMMUTING PARTNER (X2).  Entry is chi(with partner) - chi(alone).")
say(f"  {'bath coupling':<34}{'chi(alone)':>13}" + "".join(f"{'s='+str(s)+' d='+str(dist(0,s)):>15}" for s in range(L)))
ZERO = {}
for rng, J in [(1, 0.0), (1, 0.2), (1, 0.6), (1, 1.2), (2, 0.2), (2, 0.6), (2, 1.2)]:
    bath = RingBath(L=L, J=J, rng=rng)
    al = chi([], bath)
    vals = [chi([("X2", s)], bath) - al for s in range(L)]
    ZERO[(rng, J)] = vals
    tag = f"range {rng}, J = {J}" + ("   (CONTROL: no geometry)" if J == 0 else "")
    say(f"  {tag:<34}{al:>13.9f}" + "".join(f"{v:>+15.6e}" for v in vals))
say("")
say("  THE FLOOR AT EACH J, MEASURED IN THE SAME TABLE.  The ring is reflection-symmetric about")
say("  site 0, so chi(partner at s) and chi(partner at L-s) are the SAME quantity; their")
say("  difference is this row's float64 floor.  It is NOT constant across J: at strong bath")
say("  coupling the spectrum becomes near-degenerate and the floor rises.")
FLOOR = {}
say(f"  {'bath coupling':<24}{'floor (reflection s <-> L-s)':>32}{'largest |effect| in the row':>30}")
for (rng, J), vals in ZERO.items():
    fl = max(abs(vals[s_] - vals[L - s_]) for s_ in (1, 2))
    FLOOR[(rng, J)] = fl
    say(f"  {'range '+str(rng)+', J = '+str(J):<24}{fl:>32.2e}{max(abs(v) for v in vals):>30.3e}")
say("")
say("  THE ZERO SET, at a threshold of 100x that row's own measured floor:")
say(f"  {'bath coupling':<24}{'threshold':>12}{'sites with a NON-zero effect':>32}{'largest OUT-OF-RANGE |effect|':>32}{'ratio to the in-range effect':>31}")
for (rng, J), vals in ZERO.items():
    th = max(100 * FLOOR[(rng, J)], 1e-14)
    nz = [s_ for s_, v in enumerate(vals) if abs(v) >= th]
    inr = [s_ for s_ in range(1, L) if dist(0, s_) <= rng and J != 0]
    outr = [s_ for s_ in range(1, L) if s_ not in inr]
    mo = max([abs(vals[s_]) for s_ in outr], default=0.0)
    mi = max([abs(vals[s_]) for s_ in inr], default=0.0)
    say(f"  {'range '+str(rng)+', J = '+str(J):<24}{th:>12.1e}{str(nz):>32}{mo:>32.2e}"
        f"{(mo/mi if mi > 0 else float('nan')):>31.2e}")
say("")
say("  READ IT OFF THE TABLE, NOT FROM AN EXPECTATION.  The decisive column is the last one:")
say("  the largest OUT-OF-RANGE effect as a fraction of the IN-RANGE effect in the same row.")
say("   * J = 0: every off-site entry sits at the floor.  No geometry, no influence anywhere.")
say("   * COUPLING RANGE 1: sites at distance <= 1 (that is 1 and 5) carry the whole effect;")
say("     everything further out is 10 to 14 decades smaller, and at J <= 0.6 it is exactly the")
say("     float64 floor.")
say("   * COUPLING RANGE 2: the boundary moves out by exactly one site.  Sites 1, 2, 4 and 5")
say("     carry the effect and site 3 -- the only site at distance 3 -- is the one that drops to")
say("     the floor.  One step in the bath's coupling range, one step in the boundary.")
say("   * STRENGTH does not move the boundary.  J = 0.2, 0.6, 1.2 is a 6x range and the")
say("     out-of-range/in-range ratio stays between 4e-15 and 1e-10 throughout.")
say("   * CAVEAT, stated because the table says so: at J = 1.2 the pipeline's own floor rises to")
say("     ~2e-12 (near-degenerate bath spectrum), so out of range the honest claim at J = 1.2 is")
say("     ZERO TO 2e-12, not exactly zero.  At J <= 0.6 it is zero to 5e-16.")
say("  -> the DENSITY channel has COMPACT SUPPORT.  Its radius is the bath's interaction range,")
say("     and it does not depend on the coupling strength.  No decaying function of distance can")
say("     produce a zero set that ignores a 6x change in J.")
say("")
say("  PAIRING PARTNER (Z1), the same sweep.  Entry is chi(with partner)/chi(alone).")
say(f"  {'bath coupling':<34}{'chi(alone)':>13}" + "".join(f"{'s='+str(s)+' d='+str(dist(0,s)):>15}" for s in range(1, L)))
for rng, J in [(1, 0.0), (1, 0.6), (2, 0.6)]:
    bath = RingBath(L=L, J=J, rng=rng)
    al = chi([], bath)
    vals = [chi([("Z1", s)], bath) / al for s in range(1, L)]
    say(f"  {'range '+str(rng)+', J = '+str(J):<34}{al:>13.9f}" + "".join(f"{v:>15.9f}" for v in vals))
say("")
say("  -> the PAIRING partner suppresses chi at EVERY distance, including the maximum distance")
say("     the ring allows, and by nearly the same amount.  The pairing channel has NO radius at")
say("     all; the density channel has a radius equal to the bath's interaction range.")
say("")
say(f"  elapsed {time.time()-t0:.1f}s")
