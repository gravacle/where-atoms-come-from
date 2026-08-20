"""B7 -- DOES ANY OF IT DEPEND ON SEPARATION?  FORM IS SCALE-FREE, SO THIS IS DECIDABLE HERE.

Two separations exist in this construction and they are different questions.

  SEPARATION IN THE CODE.  B1/E7 settled it EXACTLY: chi depends on the symplectic pairing and
  on nothing else about a partner -- not its weight, not its physical support, not its overlap
  with the read record.  There is no falloff because there is no dependence at all.

  SEPARATION IN THE BATH.  The baseline bath has NO GEOMETRY: H_B is a sum of independent
  qubits, so "distance" between bath sites is undefined and the only meaningful variable is
  COINCIDENCE.  Answering the falloff question therefore requires GIVING the bath a geometry,
  which is a change of venue and is done here explicitly (D-17):

      RING BATH:  H_B = E * sum_j Z_j  +  J * sum_j Z_j Z_{j+1 mod L}

  with UNIFORM E, so that the only thing distinguishing one site from another is its position
  on the ring.  J = 0 is the CONTROL in the same table: at J = 0 the ring is the baseline bath
  and every site is exactly equivalent, so any distance dependence there is numerical noise.
"""
import numpy as np, sys, os, time, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from battery import build_ops, READ

t0 = time.time()

class RingBath:
    """A bath of L qubits on a ring, uniform on-site energy E, nearest-neighbour ZZ coupling J.
       Same interface as record_model.Environment as far as chi_times needs it."""
    def __init__(self, L=6, E=1.0, J=0.0, beta=2.0):
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
        self.HB = E * sum(Zs) + J * sum(Zs[j] @ Zs[(j + 1) % L] for j in range(L))
        self.probe = sum(self.site)
        self.E, self.J, self.L = E, J, L
    def thermal(self):
        w, V = np.linalg.eigh(self.HB)
        p = np.exp(-self.beta * w); p /= p.sum()
        return (V * p) @ V.conj().T

L = 6
N = 6
TIMES = np.linspace(1.0, 13.0, 25)
ops, _, _ = build_ops(N)
R = ops[READ][0]
say("=" * 122)
say("B7   SEPARATION.  RING BATH OF 6 QUBITS, UNIFORM ON-SITE ENERGY, NEAREST-NEIGHBOUR ZZ.")
say("=" * 122)
say(f"  carrier n = {N} (code space {2**(N-2)}), bath ring L = {L}, joint dimension {2**(N-2)*2**L}")
say("  distances on the ring: d(0,1) = d(0,5) = 1, d(0,2) = d(0,4) = 2, d(0,3) = 3.")
say("  The read record is always written on ring site 0.")

def chi(parts, bath, lam=0.8):
    so = [(R, 0)] + [(ops[l][0], s) for l, s in parts]
    return float(np.mean(chi_times(so, R, bath, lam, TIMES)))

for J in (0.0, 0.6):
    bath = RingBath(L=L, E=1.0, J=J, beta=2.0)
    tag = "CONTROL  J = 0 (no geometry: every site exactly equivalent)" if J == 0 else f"J = {J} (the ring now has a geometry)"
    say("")
    say("-" * 122)
    say(f"  {tag}")
    say("-" * 122)
    alone = chi([], bath)
    say(f"  chi(alone) = {alone:.14f}")
    say("")
    say("  (a) TWO-BODY.  One PAIRING partner (Z1) at ring distance d from the read record.")
    say(f"      {'d':>4}{'chi':>20}{'chi/alone':>14}{'spread vs d=1':>16}")
    v = []
    for d in (1, 2, 3):
        c = chi([("Z1", d)], bath); v.append(c)
        say(f"      {d:>4}{c:>20.14f}{c/alone:>14.9f}{c-v[0]:>+16.2e}")
    say(f"      -> range over d: {max(v)-min(v):.3e}")
    say("")
    say("  (b) TWO-BODY CONTROL.  One COMMUTING partner (X2) at ring distance d.")
    say(f"      {'d':>4}{'chi':>20}{'chi - alone':>16}")
    w = []
    for d in (1, 2, 3):
        c = chi([("X2", d)], bath); w.append(c)
        say(f"      {d:>4}{c:>20.14f}{c-alone:>+16.2e}")
    say(f"      -> worst |chi - alone|: {max(abs(x-alone) for x in w):.3e}")
    say("")
    say("  (c) THREE-BODY.  Pairing partner Z1 fixed on ring site 1; commuting partner X2 moved")
    say("      to ring site s.  s = 1 is COINCIDENT with Z1; the others are separated from it.")
    say(f"      {'s':>4}{'d(s, Z1)':>10}{'d(s, R)':>9}{'chi':>20}{'chi - (Z1 only)':>18}")
    base = chi([("Z1", 1)], bath)
    tb = []
    for s in range(L):
        c = chi([("Z1", 1), ("X2", s)], bath)
        d1 = min((s - 1) % L, (1 - s) % L); d0 = min(s % L, (-s) % L)
        tb.append((s, d1, c - base))
        say(f"      {s:>4}{d1:>10}{d0:>9}{c:>20.14f}{c-base:>+18.6e}")
    off = [abs(x[2]) for x in tb if x[1] != 0]
    coi = [abs(x[2]) for x in tb if x[1] == 0]
    say(f"      coincident (d=0): |effect| = {coi[0]:.6e}")
    say(f"      separated  (d>0): largest |effect| = {max(off):.3e}, smallest = {min(off):.3e}")
    if J == 0:
        say("      -> at J = 0 the effect is EXACTLY zero at every non-zero separation and O(0.1)")
        say("         at coincidence.  Not a falloff: a step function of coincidence.")
    else:
        say("      -> with a geometry present, look at whether the separated entries have left zero,")
        say("         and by how much relative to the coincident one.")
say("")
say(f"  elapsed {time.time()-t0:.1f}s")
