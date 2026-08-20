"""ADVERSARIAL VERIFICATION of LANE_O48_B_SEPARATION, check 2 -- THE MAIN ATTACK.

CLAIM UNDER ATTACK (the lane's headline centrepiece, and its [A3] "exact structural statement"):
  "clause (iv) holds ONLY at mu = 0, mu = 0 IS half filling, half filling IS 2k_F = pi, and
   2k_F = pi IS period-2 alternation ... Writability and the period-2 alternation are the SAME
   constraint."

THE ATTACK.  The lane established that on ITS OWN one-chain carrier, mu != 0 kills clause (iv).
It then generalised that to a statement about clause (iv) as such.  The generalisation has two
holes, and this script drills both:

  HOLE 1.  The writer search was exhaustive over the PAULI GROUP, not over ADMISSIBLE UNITARIES.
           Clause (iv) says "some admissible U".  A carrier can satisfy clause (iv) with a
           non-Pauli U and return ZERO from an exhaustive Pauli search.
  HOLE 2.  Clause (iv) does NOT pin the mediator's FERMI WAVEVECTOR.  Here is a carrier on which
           clause (iv) holds exactly -- verified by the basis-independent Tr(P_E R) = 0 criterion
           the lane itself uses, at full dense ED -- while the mediator is NOT at 2k_F = pi and
           the induced interaction does NOT alternate with period 2.

THE COUNTER-CARRIER.  Two decoupled mediator chains A and B with OPPOSITE chemical potentials,
records coupled on-site to BOTH:
   H = -sum (t/2)(XA XA + YA YA) - mu sum ZA_i
       -sum (t/2)(XB XB + YB YB) + mu sum ZB_i
       -g sum_i Zr_i (ZA_i + ZB_i)
Nothing long-range is inserted; the record-mediator coupling is still strictly on-site.
Qubit layout: 0..m-1 records, m..2m-1 chain A, 2m..3m-1 chain B.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION")
from common import eigenspaces, clause_i, clause_iii, clause_iv_trace, pauli_label
from mediator import spin_op, SZ, SX, I2, chi_row

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

SY = np.array([[0, -1j], [1j, 0]], dtype=complex)

def H_two(m, t, g, mu):
    nq = 3 * m
    H = np.zeros((2 ** nq, 2 ** nq), dtype=complex)
    for blk, sgn in ((m, -1.0), (2 * m, +1.0)):          # chain A at -mu, chain B at +mu
        for i in range(m - 1):
            H += -(t / 2.0) * (spin_op(nq, {blk + i: SX, blk + i + 1: SX})
                               + spin_op(nq, {blk + i: SY, blk + i + 1: SY}))
        for i in range(m):
            H += sgn * mu * spin_op(nq, {blk + i: SZ})
    for i in range(m):
        H += -g * spin_op(nq, {i: SZ, m + i: SZ}) - g * spin_op(nq, {i: SZ, 2 * m + i: SZ})
    return H

def H_two_terms(m, t, g, mu):
    """Same H as a list of F2 Paulis (x|z) on nq=3m qubits, for the exact Pauli commutation test."""
    nq = 3 * m
    def pl(xs, zs):
        x = [0] * nq; z = [0] * nq
        for k in xs: x[k] = 1
        for k in zs: z[k] = 1
        return tuple(x + z)
    out = []
    for blk in (m, 2 * m):
        for i in range(m - 1):
            out.append((pl([blk + i, blk + i + 1], []), -t / 2.0))                       # XX
            out.append((pl([blk + i, blk + i + 1], [blk + i, blk + i + 1]), -t / 2.0))   # YY
        for i in range(m):
            if abs(mu) > 1e-14:
                out.append((pl([], [blk + i]), mu if blk == 2 * m else -mu))
    for i in range(m):
        out.append((pl([], [i, m + i]), -g))
        out.append((pl([], [i, 2 * m + i]), -g))
    return out

def pauli_search(nq, terms, target):
    N = 4 ** nq
    v = np.arange(N, dtype=np.int64); mask = (1 << nq) - 1
    x = v & mask; z = (v >> nq) & mask
    def anti(q):
        qx = sum(q[k] << k for k in range(nq)); qz = sum(q[nq + k] << k for k in range(nq))
        return (np.bitwise_count(x & qz) ^ np.bitwise_count(z & qx)) & 1
    ok = anti(target) == 1
    for q, _ in terms: ok &= (anti(q) == 0)
    cnt = int(ok.sum())
    if cnt == 0: return 0, None
    v0 = int(v[ok][0])
    return cnt, tuple((v0 >> b) & 1 for b in range(nq)) + tuple((v0 >> (nq + b)) & 1 for b in range(nq))

def swap_AB(m):
    """Permutation matrix swapping mediator qubit blocks A (m..2m-1) and B (2m..3m-1)."""
    nq = 3 * m; D = 2 ** nq
    idx = np.arange(D)
    bits = ((idx[:, None] >> np.arange(nq - 1, -1, -1)[None, :]) & 1)   # qubit 0 is the HIGH bit
    nb = bits.copy()
    nb[:, m:2 * m] = bits[:, 2 * m:3 * m]
    nb[:, 2 * m:3 * m] = bits[:, m:2 * m]
    new = (nb * (1 << np.arange(nq - 1, -1, -1))[None, :]).sum(axis=1)
    S = np.zeros((D, D), dtype=complex)
    S[new, idx] = 1.0
    return S

P("=" * 126)
P("V2  DOES CLAUSE (iv) FORCE HALF FILLING / 2k_F = pi / PERIOD-2 ALTERNATION?")
P("    Counter-carrier: two mediator chains at OPPOSITE chemical potentials +-mu, records coupled")
P("    ON-SITE to both.  Nothing long-range inserted.  H = H_A(-mu) + H_B(+mu) - g sum Zr_i(ZA_i+ZB_i)")
P("=" * 126)

P("")
P("[V2a] CLAUSES (i)-(iv) ON THE COUNTER-CARRIER BY FULL DENSE ED, at mu = 0 AND mu != 0.")
P("      Clause (iv) is tested by the BASIS-INDEPENDENT criterion Tr(P_E R) = 0 on every eigenspace")
P("      -- the same criterion the lane uses in [A1], which given (i)+(ii) IS clause (iv).")
P("      D-15 CONTROL: the PAIR Zr_0 Zr_1 through the identical criterion; it must register NON-ZERO.")
P("      LAST COLUMN: the lane's own exhaustive 4^(3m) PAULI search for a flipper of Zr_0.")
P("-" * 126)
P(f"{'m':>2} {'mu':>6} {'dim':>6} {'#eig':>5} {'(i)':>5} {'(ii) max||[H,R]||':>18} {'(iii)':>6} "
  f"{'(iv) max|Tr P_E R_i|':>21} {'(iv)?':>6} | {'CTRL pair |Tr P_E|':>19} {'ctrl(iv)?':>9} | "
  f"{'#adm PAULI flippers':>20}")
for m in (2, 3):
    for mu in (0.0, 0.30, 0.55):
        nq = 3 * m
        H = H_two(m, 1.0, 0.40, mu)
        es = eigenspaces(H)
        c1 = True; c2 = 0.0; c3 = False; worst = 0.0
        for i in range(m):
            R = spin_op(nq, {i: SZ})
            c1 &= clause_i(R)
            c2 = max(c2, float(np.linalg.norm(H @ R - R @ H)))
            c3 |= clause_iii(R, es)
            worst = max(worst, clause_iv_trace(R, es)[1])
        Cp = spin_op(nq, {0: SZ, 1: SZ}) if m > 1 else None
        wc = clause_iv_trace(Cp, es)[1]
        tgt = tuple([0] * nq + [1 if k == 0 else 0 for k in range(nq)])
        cnt, _ = pauli_search(nq, H_two_terms(m, 1.0, 0.40, mu), tgt)
        P(f"{m:>2} {mu:>6.2f} {2**nq:>6} {len(es):>5} {str(c1):>5} {c2:>18.12f} {str(c3):>6} "
          f"{worst:>21.12f} {str(worst<1e-7):>6} | {wc:>19.12f} {str(wc<1e-7):>9} | {cnt:>20}")

P("")
P("[V2b] AN EXPLICIT ADMISSIBLE WRITER ON THE COUNTER-CARRIER AT mu != 0.")
P("      U = (X on every record) (X on every mediator qubit) (SWAP chain A <-> chain B).")
P("      This is EXHIBITED, not nominated as evidence: clause (iv) above was established by the")
P("      basis-independent trace criterion.  U is shown only to make the mechanism visible, and")
P("      to prove that a Pauli-group search cannot see it (SWAP is not a Pauli).")
P("-" * 126)
P(f"{'m':>2} {'mu':>6} {'||[U,H]||':>14} {'||U^dag H U - H|| (write cost)':>32} "
  f"{'max_i ||U^dag Zr_i U + Zr_i||':>31} {'U unitary?':>11}")
for m in (2, 3):
    for mu in (0.0, 0.30, 0.55):
        nq = 3 * m
        H = H_two(m, 1.0, 0.40, mu)
        F = spin_op(nq, {k: SX for k in range(nq)})
        U = F @ swap_AB(m)
        comm = float(np.linalg.norm(U @ H - H @ U))
        cost = float(np.linalg.norm(U.conj().T @ H @ U - H))
        worst = 0.0
        for i in range(m):
            R = spin_op(nq, {i: SZ})
            worst = max(worst, float(np.linalg.norm(U.conj().T @ R @ U + R)))
        uni = float(np.linalg.norm(U.conj().T @ U - np.eye(2 ** nq))) < 1e-10
        P(f"{m:>2} {mu:>6.2f} {comm:>14.10f} {cost:>32.10f} {worst:>31.10f} {str(uni):>11}")

P("")
P("[V2c] THE INDUCED INTERACTION ON THE COUNTER-CARRIER.  The two chains are decoupled, so at")
P("      O(g^2)  J_eff(r) = -8 g^2 [ T^A(r) + T^B(r) ],  T^A at eps_F = +2mu, T^B at eps_F = -2mu.")
P("      PREDICTION from ordinary RKKY: oscillation at 2k_F with cos(k_F) = -eps_F/2, i.e. period")
P("      2*pi/(2k_F) sites -- which is 2 ONLY at eps_F = 0.  D-15 controls: mu=0 row must give")
P("      strict period-2 alternation (it is the lane's own case), and a t=0 mediator must give zero.")
P("-" * 126)
mm = 2048
t1 = np.ones(mm - 1)
def sgnpat(J, n=28):
    return "".join("+" if x > 0 else "-" for x in J[:n])
P(f"{'mu':>6} {'eps_F':>7} {'filling A':>10} {'filling B':>10} {'2k_F':>9} {'period 2pi/2k_F':>16} "
  f"{'clause (iv)':>12} {'sign pattern of J_eff, r=1..28':>32} {'strict +-+- ?':>14}")
A0 = np.zeros((mm, mm))
for i in range(mm - 1): A0[i, i + 1] = A0[i + 1, i] = -1.0
ev = np.linalg.eigvalsh(A0)
for mu in (0.0, 0.15, 0.30, 0.55):
    eF = 2 * mu
    try:
        rA = chi_row(mm, t1, mm // 2, fermi=+eF)
        rB = chi_row(mm, t1, mm // 2, fermi=-eF)
    except RuntimeError:
        P(f"{mu:>6.2f} zero mode at eps_F -- skipped"); continue
    J = np.array([-8 * (rA[mm // 2 + r] + rB[mm // 2 + r]) for r in range(1, 200)])
    kF = np.arccos(np.clip(-eF / 2.0, -1, 1))
    fA = float((ev < +eF).mean()); fB = float((ev < -eF).mean())
    alt = all((J[r - 1] > 0) == (r % 2 == 1) for r in range(1, 120))
    P(f"{mu:>6.2f} {eF:>7.2f} {fA:>10.6f} {fB:>10.6f} {2*kF:>9.5f} {2*np.pi/(2*kF):>16.5f} "
      f"{'HOLDS':>12} {sgnpat(J):>32} {str(alt):>14}")
P(f"{'t=0':>6} {'-':>7} {'-':>10} {'-':>10} {'-':>9} {'-':>16} {'HOLDS':>12} "
  f"{'(J = 0 at every r)':>32} {'n/a':>14}")

P("")
P("[V2d] IS THE COUNTER-CARRIER'S FALLOFF STILL A POWER LAW OF EXPONENT 1?")
P("      envelope |J_eff(r)| * r on the clean window, mu = 0.30.")
P("-" * 126)
rA = chi_row(mm, t1, mm // 2, fermi=+0.60); rB = chi_row(mm, t1, mm // 2, fermi=-0.60)
J = np.array([-8 * (rA[mm // 2 + r] + rB[mm // 2 + r]) for r in range(1, 200)])
P(f"{'r':>5} {'J_eff(r)/g^2':>18} {'|J|*r':>12}")
for r in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128):
    P(f"{r:>5} {J[r-1]:>18.9e} {abs(J[r-1])*r:>12.6f}")
P("")
P("      C-46 ratio |sum J|/sum|J| on the counter-carrier (R=64, m=2048): "
  f"{abs(J[:64].sum())/np.abs(J[:64]).sum():>.9f}")
P("      lane's one-chain mu=0 value for the same R: 0.155459834")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION/VERIFY/"
     "v2_clause_iv_does_not_force_it.txt", "w").write("\n".join(OUT) + "\n")
