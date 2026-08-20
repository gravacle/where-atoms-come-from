"""LANE_O48_B_SEPARATION -- PART 2, STEP C: THE SAME QUESTION ASKED FOUR MORE WAYS.

p2b used ONE mediator (a hopping chain) and ONE gap knob (dimerisation).  If the answer there is
a property of the mediator rather than of that particular knob, it must survive:
  [C1] a SECOND, INDEPENDENT gap knob -- pairing anisotropy w, a different term in H entirely;
  [C2] a DISORDERED mediator with all-distinct couplings (D-22's own prescription);
  [C3] DILUTE records, occupying only a sublattice of the mediator;
  [C4] mediators of DIMENSION 2 and 3, where the exponent has somewhere else to go.
Every falloff reported here is INDUCED: H still contains no record-record term (p2b [B0]).
The record PLACEMENT in [C3] is INSERTED and is labelled as such.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION")
from mediator import (hop_A, pair_B, j_eff_fd, chi_row, chi_row_general, square_lattice_A,
                      E0_fields, H_full_dense, H_full_terms, spin_op, SZ)
from common import fit_power_vs_exp, eigenspaces, clause_iv_trace, clause_iii, pauli_label

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

def gap_of(A0, B):
    return float(np.linalg.svd(A0 + B, compute_uv=False).min())

P("=" * 128)
P("PART 2 STEP C -- THE SAME QUESTION ASKED FOUR MORE WAYS")
P("=" * 128)

# =============================================================== [C0] instrument
P("")
P("-" * 128)
P("[C0] THE INSTRUMENT USED HERE.  d^2 E0/dh_i dh_j by a CENTRAL SECOND DIFFERENCE of the exact")
P("     block energy in the on-site fields h_i = g z_i, taken at ZERO background so the mediator is")
P("     never doped away from the half filling that clause (iv) demands.  It equals J_eff(i,j)/g^2.")
P("     Validated below against the independent closed-form orbital sum, which only exists at w=0.")
P('     Sizes: m=256 in [C0]; m=1024 in [C1]; m=2048 in [C2] and [C3]; 4032 (2D) and 4080 (3D) in [C4].')
P("-" * 128)
m = 256; t = np.ones(m - 1); i0 = 128
A0 = hop_A(m, t); Bz = np.zeros((m, m))
js = [i0 + r for r in (1, 2, 4, 8, 16, 32)]
fd = j_eff_fd(A0, Bz, i0, js, delta=1e-2)
ex = np.array([-8 * chi_row(m, t, i0)[j] for j in js])
P(f"{'r':>4} {'second difference':>22} {'closed-form orbital sum':>26} {'relative difference':>21}")
for k, j in enumerate(js):
    P(f"{j-i0:>4} {fd[k]:>22.12e} {ex[k]:>26.12e} {abs(fd[k]-ex[k])/abs(ex[k]):>21.3e}")
P("")
P("READ: the two independent instruments agree to ~2e-4 relative, which is the O(delta^2) systematic")
P("      of the difference stencil.  The second-difference route is therefore trusted where the")
P("      orbital sum does not apply (pairing, and arbitrary lattices).")

# =============================================================== [C1] second gap knob
P("")
P("-" * 128)
P("[C1] A SECOND, INDEPENDENT GAP KNOB: PAIRING ANISOTROPY w, the (XaXa - YaYa) term.  A DIFFERENT")
P("     term in H from dimerisation, gapping the mediator by a different mechanism.  p2a verified by")
P('     exhaustive Pauli search that clause (iv) and the free writer survive w != 0 (for w < t).')
P("     m = 1024, i0 = 512, envelope fitted on r in [6,40], well inside the clean window r <= m/16.")
P("     THE GAP IS THE BULK GAP.  The open XY chain carries ONE near-zero MAJORANA EDGE MODE whose")
P("     energy is 1.06e-23 at w=0.05 and 6.90e-245 at w=0.50; it is a boundary state, not the gap,")
P("     and reporting it as the gap would be an artefact.  The bulk gap is the next singular value")
P("     and is printed beside the analytic 2|w| so the identification can be checked.")
P("-" * 128)
mC = 1024; iC = mC // 2
A0C = hop_A(mC, np.ones(mC - 1)); BzC = np.zeros((mC, mC))
rsC = list(range(6, 41, 2))
P(f"{'w':>6} {'edge mode':>12} {'bulk gap':>10} {'2|w|':>8} {'J(r=8)':>13} {'J(r=16)':>13} "
  f"{'J(r=32)':>13} {'xi':>8} {'xi*bulkgap':>11} {'p':>8} {'OUT pow':>9} {'OUT exp':>9} {'VERDICT':>14}")
for wv in (0.0, 0.02, 0.05, 0.10, 0.20, 0.50):
    B = pair_B(mC, np.full(mC - 1, wv))
    sig = np.sort(np.linalg.svd(A0C + B, compute_uv=False))
    edge, bulk = float(sig[0]), float(sig[1])
    vals = j_eff_fd(A0C, B, iC, [iC + r for r in rsC], delta=1e-2)
    env = np.array([((-1) ** (r + 1)) * vals[k] for k, r in enumerate(rsC)])
    f = fit_power_vs_exp(rsC, env, split=24, floor=1e-12)
    pick = lambda R: vals[rsC.index(R)]
    if f.get("ok"):
        if f["pow_out"] < f["exp_out"] / 3.0:   v = "POWER"
        elif f["exp_out"] < f["pow_out"] / 3.0: v = "EXPONENTIAL"
        else:                                    v = "CANNOT DECIDE"
        P(f"{wv:>6.2f} {edge:>12.2e} {bulk:>10.6f} {2*wv:>8.4f} {pick(8):>13.5e} {pick(16):>13.5e} "
          f"{pick(32):>13.5e} {f['xi']:>8.4f} {f['xi']*bulk:>11.4f} {f['p_exponent']:>8.4f} "
          f"{f['pow_out']:>9.4f} {f['exp_out']:>9.4f} {v:>14}")
    else:
        P(f"{wv:>6.2f} {edge:>12.2e} {bulk:>10.6f} {2*wv:>8.4f} "
          f"only {f['n']} points survived the declared floor 1e-12 -- NOT CLASSIFIED")
P("")
P('READ (filled from the numbers above): the bulk gap tracks 2|w| to five figures (0.040425 vs')
P('      0.0400, 0.100180 vs 0.1000, 0.200091 vs 0.2000, 0.400045 vs 0.4000, 1.000014 vs 1.0000), so')
P('      the Majorana edge mode was correctly separated out and the gap column really is the bulk gap.')
P("      THE VERDICT COLUMN IS MOSTLY 'CANNOT DECIDE' AND THAT IS THE HONEST ANSWER FOR THOSE ROWS.")
P('      Only w=0.00 (POWER) and w=0.20 (EXPONENTIAL) beat the other form by the required factor of 3.')
P('      Two reasons, both stated rather than worked around: (1) the true gapped form is an exponential')
P('      TIMES a power prefactor, so a two-parameter pure exponential is misspecified and cannot win by')
P('      3x at the smaller gaps; (2) at w=0.50 the signal has already reached the declared floor --')
P('      J(r=16) prints as 0.00000e+00 and J(r=32) as -2.8e-10, which is numerical noise.')
P('      WHAT THE TABLE DOES SETTLE, WITH NO MODEL SELECTION AT ALL: the fitted power exponent p climbs')
P('      MONOTONICALLY with the bulk gap -- 1.0238, 1.6527, 2.5202, 3.9012, 6.6349, 7.8115 -- and a')
P("      genuine power law must have an exponent that does NOT depend on the mediator's gap.  So the")
P('      gapped rows are not power laws whatever else they are, while the w=0 row sits at p=1.0238 and')
P('      agrees with the structurally independent dimerisation knob of p2b [B4] (p=1.023).  xi falls')
P('      monotonically 13.21, 8.08, 5.26, 3.38, 1.98, 1.88 as the gap rises, while xi*bulkgap CLIMBS')
P('      0.04, 0.33, 0.53, 0.68, 0.79, 1.88 rather than holding constant -- so xi ~ 1/gap is a')
P('      direction over this range, not a law, and it is reported as a direction.')

# =============================================================== [C2] disorder
P("")
P("-" * 128)
P("[C2] A DISORDERED MEDIATOR.  All hoppings distinct -- D-22's own prescription for putting")
P("     geometry into a carrier -- with t_i = 1 + s*(uniform on [-1,1]).  There is NO gap opened:")
P("     the band is untouched, only made inhomogeneous.  So if the falloff still turns exponential,")
P("     GAPLESS IS NOT SUFFICIENT for a power law.  w = 0 here, so the exact closed-form orbital sum")
P("     applies and m = 2048 is affordable; median over 16 disorder realisations, i0 = 1024.")
P("-" * 128)
mD = 2048; iD = mD // 2
rsD = list(range(6, 129, 2))
rngD = np.random.default_rng(2024)
P(f"{'s':>6} {'#distinct t_i':>14} {'med |J| r=8':>14} {'r=16':>13} {'r=32':>13} {'r=64':>13} "
  f"{'r=128':>13} {'xi':>9} {'p':>8} {'OUT pow':>9} {'OUT exp':>9} {'VERDICT':>14}")
for sd in (0.0, 0.1, 0.2, 0.4, 0.7):
    acc = []; ndist = 1
    for trial in range(16 if sd > 0 else 1):
        tt = 1.0 + sd * (2 * rngD.random(mD - 1) - 1.0)
        ndist = len(set(np.round(tt, 9)))
        try:
            rr = chi_row(mD, tt, iD)
        except RuntimeError:
            continue
        acc.append(np.array([abs(-8 * rr[iD + r]) for r in rsD]))
    med = np.median(np.array(acc), axis=0)
    f = fit_power_vs_exp(rsD, med, split=64, floor=1e-13)
    pick = lambda R: med[rsD.index(R)]
    if f.get("ok"):
        if f["pow_out"] < f["exp_out"] / 3.0:   v = "POWER"
        elif f["exp_out"] < f["pow_out"] / 3.0: v = "EXPONENTIAL"
        else:                                    v = "CANNOT DECIDE"
        P(f"{sd:>6.2f} {ndist:>14} {pick(8):>14.5e} {pick(16):>13.5e} {pick(32):>13.5e} "
          f"{pick(64):>13.5e} {pick(128):>13.5e} {f['xi']:>9.4f} {f['p_exponent']:>8.4f} "
          f"{f['pow_out']:>9.4f} {f['exp_out']:>9.4f} {v:>14}")
    else:
        P(f"{sd:>6.2f} {ndist:>14} not classifiable above the declared floor 1e-13")
P("")
P("READ (filled from the numbers above): the '#distinct t_i' column is the D-22 check -- 2047")
P('      distinct couplings on 2048 sites at every s>0, so the disordered carriers really do carry')
P('      geometry, against the s=0 reference with one repeated value.  WEAK DISORDER DOES NOT DESTROY')
P('      THE POWER LAW at the ranges reached: p = 1.0263, 1.0001, 1.1693 at s = 0.0, 0.1, 0.2, all')
P('      classified POWER out of sample.  STRONG DISORDER SUPPRESSES THE TAIL HARD, but the fits cannot')
P('      name the replacement form: at s=0.7 the median |J| at r=128 is 6.80e-08 against 2.14e-03 in')
P('      the clean case, a factor of 3.1e4, while at r=8 it is down by only a factor of 12 -- a')
P('      suppression that GROWS WITH SEPARATION, which is what localisation looks like.  Both s=0.4 and')
P('      s=0.7 return CANNOT DECIDE and are left there (D-20).  THE CONCLUSION THAT SURVIVES: A GAPLESS')
P('      MEDIATOR IS NOT SUFFICIENT FOR A POWER LAW -- the band must also be clean enough -- and how')
P('      much disorder is too much is not settled at these sizes.')

# =============================================================== [C3] dilute records
P("")
P("-" * 128)
P("[C3] DILUTE RECORDS: records occupy only every k-th mediator site.  The MEDIATOR is unchanged;")
P("     only WHERE the records sit changes.  *** THE PLACEMENT IS INSERTED, NOT INDUCED. ***  The")
P("     falloff is induced; the sublattice choice is a free choice the five clauses do not fix.")
P("     Separations below are in MEDIATOR SITES.  m=2048, clean gapless mediator.")
P("-" * 128)
mL = 2048; aL = mL // 2; row = chi_row(mL, np.ones(mL - 1), aL)
P(f"{'k':>3} {'record separations used':>26} {'sign pattern (first 12)':>25} {'R':>5} "
  f"{'sum J':>13} {'sum |J|':>13} {'|sum|/sum|.|':>14}")
for k in (1, 2, 3, 4):
    seps = [k * n for n in range(1, 1 + 128 // k)]
    seps = [r for r in seps if aL + r < mL - 2]
    J = np.array([-8 * row[aL + r] for r in seps])
    pat = "".join("+" if x > 0 else "-" for x in J[:12])
    for R in (16, 64, min(128, len(J))):
        if R > len(J): continue
        v = J[:R]
        P(f"{k:>3} {('every %d-th site' % k):>26} {pat:>25} {R:>5} {v.sum():>13.8f} "
          f"{np.abs(v).sum():>13.8f} {abs(v.sum())/np.abs(v).sum():>14.9f}")
P("")
P("READ (filled from the numbers above): with records on EVERY site (k=1) the induced coupling")
P("      alternates and the cancellation ratio is small.  With records on every SECOND site (k=2)")
P("      every separation is even, the sign never changes, and the ratio is exactly 1.000000000 --")
P("      the same induced 1/r interaction now ACCUMULATES instead of cancelling.  k=3 alternates")
P("      again; k=4 does not.  *** WHETHER THIS INTERACTION SCREENS OR ACCUMULATES IS DECIDED BY")
P("      THE SUBLATTICE THE RECORDS OCCUPY, WHICH IS PUT IN BY HAND AND IS NOT FIXED BY ANY OF THE")
P("      FIVE CLAUSES. ***  That is a finding about the clauses, not a source.")

# clause check for dilute records
P("")
P("     Clause re-check for the DILUTE carrier (records on even mediator sites only), full ED:")
P(f"{'mediator sites':>16} {'#records':>9} {'dim':>7} {'(iii)':>7} {'(iv) max|Tr P_E R|':>21} "
  f"{'(iv)?':>7} {'#adm. Pauli flippers':>21}")
for mm in (4, 6):
    nrec = mm // 2
    nq = nrec + mm
    # build H: mediator XY chain on qubits nrec..nrec+mm-1, records on even mediator sites
    import numpy as _np
    from mediator import SX, SY, I2, _kron
    H = _np.zeros((2 ** nq, 2 ** nq), dtype=complex)
    for i in range(mm - 1):
        H += -(0.5) * (spin_op(nq, {nrec + i: SX, nrec + i + 1: SX})
                       + spin_op(nq, {nrec + i: SY, nrec + i + 1: SY}))
    for a in range(nrec):
        H += -0.4 * spin_op(nq, {a: SZ, nrec + 2 * a: SZ})
    es = eigenspaces(H)
    worst = max(clause_iv_trace(spin_op(nq, {a: SZ}), es)[1] for a in range(nrec))
    c3 = all(clause_iii(spin_op(nq, {a: SZ}), es) for a in range(nrec))
    # exhaustive Pauli search
    N = 4 ** nq; v = _np.arange(N, dtype=_np.int64); msk = (1 << nq) - 1
    x = v & msk; z = (v >> nq) & msk
    def anti(qx, qz):
        return (_np.bitwise_count(x & qz) ^ _np.bitwise_count(z & qx)) & 1
    terms = []
    for i in range(mm - 1):
        terms.append(((1 << (nrec + i)) | (1 << (nrec + i + 1)), 0))
        terms.append(((1 << (nrec + i)) | (1 << (nrec + i + 1)),
                      (1 << (nrec + i)) | (1 << (nrec + i + 1))))
    for a in range(nrec):
        terms.append((0, (1 << a) | (1 << (nrec + 2 * a))))
    ok = anti(0, 1) == 1
    for qx, qz in terms: ok &= (anti(qx, qz) == 0)
    P(f"{mm:>16} {nrec:>9} {2**nq:>7} {str(c3):>7} {worst:>21.12f} {str(worst<1e-7):>7} "
      f"{int(ok.sum()):>21}")
P("     READ: the dilute carrier still satisfies (iii) and (iv) and still admits admissible")
P("           flippers, so the k=2 placement is a legitimate record set, not a cheat.")

# =============================================================== [C4] dimension
P("")
P("-" * 128)
P("[C4] DOES THE EXPONENT TRACK THE MEDIATOR'S DIMENSION?  Same construction, mediator on a 1D")
P("     chain / 2D rectangle / 3D box, all bipartite (so half filling stays particle-hole symmetric")
P("     and clause (iv)'s flip survives), all open, all nearest-neighbour, all at half filling.")
P("     NOTE: in 1D the mediator is a genuinely LOCAL spin chain under Jordan-Wigner; in 2D and 3D")
P("     the mediator is a FERMIONIC bath and its spin representation is non-local.  Clauses (i)-(iv)")
P("     are unaffected (none of them needs locality); only clause (v) would be, and (v) is NOT")
P("     claimed anywhere in this lane.")
P("-" * 128)
def local_slope(rs, vals):
    rs = np.asarray(rs, float); v = np.abs(np.asarray(vals, float))
    return [(rs[i], float(-(np.log(v[i+1]) - np.log(v[i])) / (np.log(rs[i+1]) - np.log(rs[i]))))
            for i in range(len(rs) - 1)]
P(f"{'mediator':>22} {'sites':>7} {'min |eps| (no zero mode)':>26} {'r':>4} {'|J|/g^2':>14} "
  f"{'|J| r^1':>10} {'|J| r^2':>10} {'|J| r^3':>10} {'local slope':>12}")
# 1D
m1 = 4096; A1 = hop_A(m1, np.ones(m1 - 1)); c1 = chi_row(m1, np.ones(m1 - 1), m1 // 2)
e1 = np.linalg.eigvalsh(A1); rs1 = [2, 4, 6, 8, 12, 16, 24, 32]
v1 = [abs(-8 * c1[m1 // 2 + r]) for r in rs1]
sl1 = dict(local_slope(rs1, v1))
for k, r in enumerate(rs1):
    P(f"{'1D chain 4096':>22} {m1:>7} {abs(e1).min():>26.6e} {r:>4} {v1[k]:>14.6e} "
      f"{v1[k]*r:>10.5f} {v1[k]*r**2:>10.4f} {v1[k]*r**3:>10.3f} {sl1.get(r, float('nan')):>12.5f}")
# 2D
Lx, Ly = 63, 64
A2, idx2 = square_lattice_A(Lx, Ly)
e2 = np.linalg.eigvalsh(A2)
i2 = idx2(Lx // 2, Ly // 2)
c2 = chi_row_general(A2, i2)
rs2 = [2, 4, 6, 8, 10, 12, 16, 20]
v2 = [abs(-8 * c2[idx2(Lx // 2 + r, Ly // 2)]) for r in rs2]
sl2 = dict(local_slope(rs2, v2))
for k, r in enumerate(rs2):
    P(f"{'2D 63x64':>22} {Lx*Ly:>7} {abs(e2).min():>26.6e} {r:>4} {v2[k]:>14.6e} "
      f"{v2[k]*r:>10.5f} {v2[k]*r**2:>10.4f} {v2[k]*r**3:>10.3f} {sl2.get(r, float('nan')):>12.5f}")
# 3D
Lx3, Ly3, Lz3 = 15, 16, 17
m3 = Lx3 * Ly3 * Lz3
A3 = np.zeros((m3, m3)); idx3 = lambda a, b, c: (a * Ly3 + b) * Lz3 + c
for a in range(Lx3):
    for b in range(Ly3):
        for c in range(Lz3):
            if a + 1 < Lx3: A3[idx3(a,b,c), idx3(a+1,b,c)] = A3[idx3(a+1,b,c), idx3(a,b,c)] = -1
            if b + 1 < Ly3: A3[idx3(a,b,c), idx3(a,b+1,c)] = A3[idx3(a,b+1,c), idx3(a,b,c)] = -1
            if c + 1 < Lz3: A3[idx3(a,b,c), idx3(a,b,c+1)] = A3[idx3(a,b,c+1), idx3(a,b,c)] = -1
e3 = np.linalg.eigvalsh(A3)
i3 = idx3(Lx3 // 2, Ly3 // 2, Lz3 // 2)
c3r = chi_row_general(A3, i3)
rs3 = [2, 3, 4, 5, 6, 7]
v3 = [abs(-8 * c3r[idx3(Lx3 // 2 + r, Ly3 // 2, Lz3 // 2)]) for r in rs3]
sl3 = dict(local_slope(rs3, v3))
for k, r in enumerate(rs3):
    P(f"{'3D 15x16x17':>22} {m3:>7} {abs(e3).min():>26.6e} {r:>4} {v3[k]:>14.6e} "
      f"{v3[k]*r:>10.5f} {v3[k]*r**2:>10.4f} {v3[k]*r**3:>10.3f} {sl3.get(r, float('nan')):>12.5f}")
P("")
P('READ (filled from the numbers above): ONLY THE 1D ANSWER IS ESTABLISHED.  In 1D |J|*r is flat at')
P('      0.3130-0.3159 from r=4 to r=24 and the local slope runs 0.94835, 0.98281, 0.99333, 0.99966,')
P('      1.00433, 1.00861, 1.01354 -- it passes through 1.000 at r=8 and drifts up only as the finite')
P('      chain begins to bite.  EXPONENT 1, measured rather than fitted.')
P('      IN 2D NO EXPONENT IS EXTRACTED.  The local slope runs 0.66053, 0.68634, 0.81630, 0.99857,')
P('      1.22867, 1.68647, 2.65558 and is still CLIMBING at the largest separation reached; it settles')
P('      nowhere.  The reachable range is r <= 20 from the centre of a 63x64 open lattice whose nearest')
P('      edge is only 31 sites away, which is too short (D-20).  All that can be said is that by r ~ 12')
P('      the 2D slope has passed 1 and is steeper than the 1D value at the same r.')
P('      IN 3D THE DATA IS UNUSABLE FOR AN EXPONENT.  The local slope is 1.93826, -0.02981, 0.49544,')
P('      7.28170, -8.35469 over r = 2..7: not monotone, and it changes sign.  A 15x16x17 box gives only')
P('      seven sites of separation from the centre.  NO 3D CLAIM IS MADE HERE.')
P('      SO THE DIMENSION QUESTION IS OPEN, AND IT MATTERS: whether the induced exponent tracks the')
P("      mediator's dimension decides whether this family could ever carry a gravity-shaped falloff in")
P('      three dimensions.  Settling it needs a mediator of order 1e5-1e6 sites, which needs a sparse')
P('      or momentum-space method rather than the dense O(m^3) eigendecomposition used here.')

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION/p2c_construction_b.txt","w").write("\n".join(OUT)+"\n")
