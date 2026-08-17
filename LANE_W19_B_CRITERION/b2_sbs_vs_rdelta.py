"""
b2 -- KORBICZ SBS  versus  QUANTUM DARWINISM R_delta.  THEY ARE NOT THE SAME CONDITION.

THE TWO CONDITIONS, STATED PRECISELY.

  (QD)  QUANTUM DARWINISM / REDUNDANCY.  Fix a system S and a partition of the environment into
        DISJOINT fragments F_1..F_m.  Let H(S) = S(rho_S).  Then
             R_delta  :=  # { k : I(S : F_k) >= (1 - delta) H(S) },
        and a PLATEAU is the statement that the curve |F| -> I(S:F) is flat at (1-delta)H(S) over
        a range of fragment sizes.  QD IS A CONDITION ON MUTUAL INFORMATIONS ALONE, evaluated one
        fragment at a time.  It never looks at the joint state of two fragments.

  (SBS) KORBICZ SPECTRUM BROADCAST STRUCTURE.  rho_{S:F_1..F_m} has SBS form iff there is a basis
        {|i>} of S with
             rho  =  sum_i p_i |i><i| (x) rho_i^{F_1} (x) ... (x) rho_i^{F_m},
        with  rho_i^{F_k} rho_j^{F_k} = 0  for every k and every i != j.
        THREE INDEPENDENT CLAUSES:
          (a) DECOHERENCE   rho is block-diagonal in the pointer basis of S;
          (b) STRONG INDEPENDENCE  the conditional environment state factorises ACROSS FRAGMENTS;
          (c) PERFECT DISTINGUISHABILITY  the conditional fragment states have orthogonal support.
        SBS looks at the JOINT state of the fragments.  Clause (b) has no counterpart in QD.

  DEFECTS MEASURED HERE (all zero iff SBS holds), in trace distance:
        D_coh  = (1/2)|| rho - Delta_S(rho) ||_1
        D_prod = max_i (1/2)|| sigma_i - (x)_k sigma_i^{(k)} ||_1
        D_orth = max_{k, i!=j} || sqrt(sigma_i^{(k)}) sqrt(sigma_j^{(k)}) ||_1     (0 = orthogonal)

ISOLATION LEDGER
  held fixed : carrier theta_9 (9 links, Gauss law prod_l X_l, dim H_ext = 512, dim H_phys = 256),
               S = {link 0}, fragments F1={1,2} F2={3,4} F3={5,6}, unobserved U={7,8},
               delta = 0.1, both criteria evaluated on the SAME rho by the SAME functions.
  moved      : the state, and (in the last block only, marked) whether U is traced out.
"""
import numpy as np, itertools
from lib_b import *

np.set_printoptions(precision=6, suppress=True)
L = 9
car = theta(L)
S = [0]; FR = [[1, 2], [3, 4], [5, 6]]; U = [7, 8]
DELTA = 0.1

print("=" * 104)
print("b2  SBS versus R_delta")
print("=" * 104)
print(f"carrier {car['name']}  L={L}  indep Gauss = {indep_gauss(car)}  dim H_phys = {phys_dim(car)}")
print(f"S = {S}   fragments = {FR}   unobserved U = {U}   delta = {DELTA}")
print()

# --------------------------------------------------------------------------- states
def ghz():
    return sym_basis_state(car, 0)

def qd_yes_sbs_no():
    """CONSTRUCTED so that every single-fragment mutual information saturates H(S) while the
       conditional environment state is CORRELATED ACROSS FRAGMENTS (SBS clause (b) fails).
         |psi> = ( |0>_0 (x) |chi>  +  |1>_0 (x) X^{(x)8} |chi> ) / sqrt2 ,
         |chi> = ( |00,00,00> + |01,01,01> )/sqrt2  (x) |00>_U .
       X^{(x)8} on the environment makes it gauge invariant; U's two branches are orthogonal so
       tracing U decoheres S exactly."""
    psi = np.zeros(1 << L, dtype=complex)
    c1 = 0
    c2 = (1 << 2) | (1 << 4) | (1 << 6)
    full = (1 << L) - 1
    for c in (c1, c2):
        psi[c] += 0.5
        psi[c ^ full] += 0.5                      # the gauge image: link 0 flipped too
    return project_physical(psi, car)

def qd_no_sbs_no():
    return haar_physical(car, 777)

def sbs_and_qd():
    """A genuine SBS state that is NOT the GHZ: the conditional fragment states are MIXED but
       still orthogonal and still a product across fragments."""
    psi = np.zeros(1 << L, dtype=complex)
    full = (1 << L) - 1
    # link0 = 0 branch: each fragment independently in {|00>,|01>}, U = |00>
    for b1 in (0, 1):
        for b2 in (0, 1):
            for b3 in (0, 1):
                c = (b1 << 2) | (b2 << 4) | (b3 << 6)
                psi[c] += 1.0
                psi[c ^ full] += 1.0
    psi /= np.linalg.norm(psi)
    return project_physical(psi, car)

# --------------------------------------------------------------------------- the two criteria
def qd_report(psi, frags):
    HS = vn_entropy(reduce_links(psi, L, S))
    per = [mi_ext(psi, L, S, f) for f in frags]
    cum = [mi_ext(psi, L, S, sum(frags[:k + 1], [])) for k in range(len(frags))]
    R = int(sum(1 for v in per if v >= (1 - DELTA) * HS - 1e-9))
    return HS, per, cum, R

def sbs_report(psi, frags, trace_out_U=True):
    order = S + sum(frags, [])
    if not trace_out_U:
        order = order + U
    rho = rdm_ordered(psi, L, order)
    dims = [2] + [1 << len(f) for f in frags] + ([1 << len(U)] if not trace_out_U else [])
    nS = 2
    # pointer basis: eigenbasis of rho_S; report degeneracy, and also use the Z basis of the link
    rhoS = ptrace_blocks(rho, dims, [0])
    ev = np.linalg.eigvalsh(rhoS)
    degenerate = abs(ev[0] - ev[1]) < 1e-8
    # dephase S in the computational (Z) basis of the link -- the einselected basis of this carrier
    d = int(np.prod(dims))
    dephased = np.zeros_like(rho)
    step = d // nS
    for i in range(nS):
        dephased[i * step:(i + 1) * step, i * step:(i + 1) * step] = \
            rho[i * step:(i + 1) * step, i * step:(i + 1) * step]
    D_coh = 0.5 * trace_norm(rho - dephased)
    ps = []; sig = []
    for i in range(nS):
        blk = rho[i * step:(i + 1) * step, i * step:(i + 1) * step]
        p = float(np.trace(blk).real); ps.append(p)
        sig.append(blk / p if p > 1e-12 else blk)
    sub = dims[1:]
    D_prod = 0.0
    marg = []
    for i in range(nS):
        m = [ptrace_blocks(sig[i], sub, [k]) for k in range(len(sub))]
        marg.append(m)
        prod = np.array([[1.0 + 0j]])
        for mm in m: prod = np.kron(prod, mm)
        D_prod = max(D_prod, 0.5 * trace_norm(sig[i] - prod))
    D_orth = 0.0
    for k in range(len(sub)):
        for i in range(nS):
            for j in range(i + 1, nS):
                D_orth = max(D_orth, trace_norm(msqrt(marg[i][k]) @ msqrt(marg[j][k])))
    return dict(p=ps, D_coh=D_coh, D_prod=D_prod, D_orth=D_orth, degenerate=degenerate,
                sbs=(D_coh < 1e-8 and D_prod < 1e-8 and D_orth < 1e-8))

# --------------------------------------------------------------------------- run
ARMS = [
    ("GHZ    magnetic GHZ (the textbook broadcast state)", ghz()),
    ("SBSm   SBS with MIXED orthogonal conditional fragment states", sbs_and_qd()),
    ("QDNS   constructed: QD plateau saturated, fragments CORRELATED given S", qd_yes_sbs_no()),
    ("HAAR   Haar-random physical state, seed 777", qd_no_sbs_no()),
]

print("=" * 104)
print(f"{'ARM':<56}{'H(S)':>8}{'Rdelta':>8}{'D_coh':>10}{'D_prod':>10}{'D_orth':>10}{'SBS':>6}")
print("-" * 104)
rows = []
for name, psi in ARMS:
    HS, per, cum, R = qd_report(psi, FR)
    sb = sbs_report(psi, FR, trace_out_U=True)
    rows.append((name, HS, per, cum, R, sb))
    print(f"{name:<56}{HS:8.5f}{R:8d}{sb['D_coh']:10.2e}{sb['D_prod']:10.2e}"
          f"{sb['D_orth']:10.2e}{('YES' if sb['sbs'] else 'no'):>6}")

print()
print("PER-FRAGMENT AND CUMULATIVE MUTUAL INFORMATIONS (the QD plateau axis has 3 points here)")
print("-" * 104)
for name, HS, per, cum, R, sb in rows:
    print(f"{name}")
    print(f"     H(S) = {HS:.6f}   p(pointer) = {np.array(sb['p'])}"
          f"   rho_S degenerate: {sb['degenerate']}")
    print(f"     I(S:F_k) for k=1,2,3      = {['%.6f' % v for v in per]}")
    print(f"     I(S:F_1..F_j) for j=1,2,3 = {['%.6f' % v for v in cum]}")

print()
print("=" * 104)
print("THE DISAGREEMENT, IN ONE ROW")
print("=" * 104)
n, HS, per, cum, R, sb = rows[2]
print(f"  QDNS:  R_delta = {R} of 3   (every disjoint fragment saturates H(S) = {HS:.6f}),")
print(f"         plateau across |F| = 1,2,3 fragments at I = {['%.6f' % v for v in cum]},")
print(f"         and yet SBS FAILS: D_prod = {sb['D_prod']:.6f} while D_coh = {sb['D_coh']:.2e}")
print(f"         and D_orth = {sb['D_orth']:.2e}.")
print("  READING: the ONLY SBS clause that fails is STRONG INDEPENDENCE -- the environment")
print("  fragments are correlated with each other GIVEN the system's pointer value.  QD's")
print("  mutual-information test cannot see this because it never evaluates a joint fragment")
print("  state.  QD IS STRICTLY WEAKER THAN SBS, and the gap is exactly clause (b).")

print()
print("=" * 104)
print("THE SECOND DISAGREEMENT: SBS NEEDS A TRACE-OUT, QD DOES NOT")
print("=" * 104)
print("  Evaluated on rho_{S:F1F2F3} with U = {7,8} traced out, GHZ is SBS.  Evaluated on the")
print("  SAME state with NOTHING traced out (all 8 environment links observed), rho is pure and")
print("  entangled, hence not separable, hence not of SBS form -- while R_delta is unchanged.")
for name, psi in ARMS[:1]:
    sb_t = sbs_report(psi, FR, trace_out_U=True)
    sb_f = sbs_report(psi, [[1, 2], [3, 4], [5, 6], [7, 8]], trace_out_U=True)
    # the second call observes ALL of the environment: order = S + all env links, nothing left over
    print(f"  GHZ, U traced out   : D_coh = {sb_t['D_coh']:.2e}  D_prod = {sb_t['D_prod']:.2e}  "
          f"D_orth = {sb_t['D_orth']:.2e}  SBS = {sb_t['sbs']}")
    print(f"  GHZ, U OBSERVED     : D_coh = {sb_f['D_coh']:.2e}  D_prod = {sb_f['D_prod']:.2e}  "
          f"D_orth = {sb_f['D_orth']:.2e}  SBS = {sb_f['sbs']}")
    HS, per, cum, R = qd_report(psi, [[1, 2], [3, 4], [5, 6], [7, 8]])
    print(f"  GHZ, 4 fragments, QD: H(S) = {HS:.6f}  R_delta = {R} of 4  "
          f"I(S:F_k) = {['%.6f' % v for v in per]}")
print("  => the two criteria do not even take the same argument.  R_delta is a function of the")
print("     global state and a fragment; SBS is a function of a STATE AFTER A PARTIAL TRACE, and")
print("     its verdict flips when the observer's cut of the environment changes with the state")
print("     held fixed.  A sweep that reports 'a record' must say WHICH criterion and WHICH cut.")

print()
print("=" * 104)
print("THE THIRD DISAGREEMENT, AND IT IS GAUGE-SPECIFIC:")
print("SBS TAKES ITS POINTER BASIS FROM rho_S, AND THE GAUSS LAW FIXES THAT BASIS TO BE ELECTRIC")
print("=" * 104)
print("  CLAIM.  On theta_L the single-link Gauss operator G = prod_l X_l anticommutes with Z_0 and")
print("  with Y_0.  For any PHYSICAL rho (G rho G = rho),  <Z_0> = Tr(rho Z_0) = Tr(rho G Z_0 G)")
print("  = -<Z_0>, so <Z_0> = 0, and likewise <Y_0> = 0.  Hence rho_{link 0} = (I + <X_0> X)/2 for")
print("  EVERY physical state: its eigenbasis is the ELECTRIC (X) basis, or is degenerate.")
print("  CONSEQUENCE.  Korbicz's SBS reads the pointer basis off rho_S.  On a gauged link that")
print("  basis is never the magnetic one.  A magnetic (holonomy) record therefore cannot be")
print("  certified by SBS on its own terms -- the pointer basis must be imported from outside.")
print()
print("  MEASURED, over every state used in this file plus 40 Haar-random physical states:")
worst_z = worst_y = 0.0
Z0 = SP(0, 1); Y0 = SP(1, 1); X0 = SP(1, 0)
allpsi = [p for _, p in ARMS] + [haar_physical(car, 9000 + s) for s in range(40)]
xs = []
for p in allpsi:
    Ep = pauli_table(p, L)
    worst_z = max(worst_z, abs(sp_expect(Ep, Z0)))
    worst_y = max(worst_y, abs(sp_expect(Ep, Y0)))
    xs.append(abs(sp_expect(Ep, X0).real))
print(f"     max |<Z_0>| over {len(allpsi)} physical states = {worst_z:.3e}")
print(f"     max |<Y_0>| over {len(allpsi)} physical states = {worst_y:.3e}")
print(f"     |<X_0>| ranges over [{min(xs):.6f}, {max(xs):.6f}] -- the only surviving direction.")
print("  CONTROL (this could have failed): the same three expectations on NON-physical states.")
rng = np.random.default_rng(4242)
nz = []
for s in range(10):
    v = rng.normal(size=1 << L) + 1j * rng.normal(size=1 << L); v /= np.linalg.norm(v)
    Ev = pauli_table(v, L)
    nz.append(abs(sp_expect(Ev, Z0)))
print(f"     max |<Z_0>| over 10 UNPROJECTED random states = {max(nz):.6f}  (nonzero, as it must be)")
