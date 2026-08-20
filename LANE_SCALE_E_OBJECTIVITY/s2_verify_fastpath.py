"""S2 -- VERIFY THE FAST PATH AGAINST record_model ITSELF.

Nothing downstream is interpretable unless the factorised broadcast state reproduces
record_model.evolve + Environment.holevo exactly.  Three carriers, single-qubit fragments,
two-qubit fragments and the whole bath, several times.

Also verified here:
  * the 4-outcome (pair) Holevo, against a reference built directly from the model's full
    evolved state with the model's own partial traces;
  * the PARITY observable R_i R_j, which is a +-1 operator, against Environment.holevo directly;
  * the required CONTROL -- k literally independent carriers with literally disjoint baths,
    built as a tensor product of two [[4,2,2]] carriers and run through the FULL model, where
    no relation can exist by construction.  Its synergy must be 0 and its single-record chi
    must not be (D-15: the zero is only worth printing beside a positive control).
"""
import sys, io, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY")
from common import *

OUT = io.StringIO()
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.write(s + "\n")

LAM = 0.8

# ---------------------------------------------------------------- reference readouts
def build_HINT(Rs, W, env):
    """sum_ij W[i,j] R_i (x) X_j, assembled as sum_j (sum_i W[i,j] R_i) (x) X_j -- nq krons,
       not k*nq, which matters at dim 4096."""
    out = None
    for j in range(env.nq):
        Aj = sum(W[i, j] * Rs[i] for i in range(len(Rs)))
        t = np.kron(Aj, env.site[j])
        out = t if out is None else out + t
    return out


def ref_chi_multi(env, r, projs, nS, fragment=None):
    """Holevo of a multi-outcome observable, assembled exactly as Environment.holevo does."""
    out = []
    for pr in projs:
        Pk = np.kron(pr, np.eye(env.dim))
        blk = Pk @ r @ Pk
        p = float(np.real(np.trace(blk)))
        if p < 1e-12: continue
        rB = env._trace_system(blk / p, nS)
        out.append((p, env._fragment(rB, fragment) if fragment is not None else rB))
    if len(out) < 2: return 0.0
    def vn(m):
        e = np.linalg.eigvalsh(m); e = e[e > 1e-13]
        return float(-(e * np.log2(e)).sum())
    av = sum(p * rb for p, rb in out)
    return max(vn(av) - sum(p * vn(rb) for p, rb in out), 0.0)


P("=" * 108)
P("S2  FAST PATH vs record_model -- exactness check")
P("=" * 108)
P("")
P("H_int = lam * sum_ij W[i,j] R_i (x) X_j handed to record_model.evolve as a full operator.")
P(f"lam = {LAM}, bath energies {ENERGIES[:4]}, beta = {BETA}, weights = 'crowded' (generic, column-normalised)")
P("")

CASES = [(4, 4, (2.0, 5.0, 9.0)), (4, 6, (3.0, 8.0)), (6, 4, (2.0, 5.0, 9.0)),
         (6, 5, (4.0,)), (8, 4, (3.0,))]
maxdiff_single = maxdiff_pair = maxdiff_par = 0.0
P(f"{'n':>3} {'k':>3} {'nq':>3} {'dim_tot':>8} {'t':>5} {'fragment':>10} {'rec':>4} | "
  f"{'chi model':>11} {'chi fast':>11} {'|diff|':>10} | {'pair model':>11} {'pair fast':>11} {'|diff|':>10} | "
  f"{'par model':>10} {'par fast':>10} {'|diff|':>10}")
P("-" * 160)
for n, nq, ts in CASES:
    car = carrier(n); k = car['k']
    H = code_hamiltonian(n); Rs = record_matrices(car)
    W = weights('crowded', k, nq)
    env = Environment(nq=nq, energies=ENERGIES[:nq], beta=BETA)
    B = Broadcast(k, nq, W, LAM, times=ts)
    m = RecordModel(H); HINT = build_HINT(Rs, W, env); nS = H.shape[0]
    fraglist = ([0], [1], [0, 1], list(range(nq))) if nS * env.dim <= 2048 else ([0], [0, 1], list(range(nq)))
    for ti, t in enumerate(ts):
        r = m.evolve(HINT, env, lam=LAM, t=t)
        for frag in fraglist:
            fr = None if len(frag) == nq else frag
            i, j = 0, 1
            cm = env.holevo(r, Rs[i], nS, fragment=fr)
            cf = B.chi_single(frag, i, ti)
            pi = [(np.eye(nS) + a * Rs[i]) @ (np.eye(nS) + b * Rs[j]) / 4
                  for a in (+1, -1) for b in (+1, -1)]
            pm = ref_chi_multi(env, r, pi, nS, fragment=fr)
            pf = B.chi_pair(frag, i, j, ti)
            am = env.holevo(r, Rs[i] @ Rs[j], nS, fragment=fr)
            af = B.chi_parity(frag, i, j, ti)
            maxdiff_single = max(maxdiff_single, abs(cm - cf))
            maxdiff_pair = max(maxdiff_pair, abs(pm - pf))
            maxdiff_par = max(maxdiff_par, abs(am - af))
            P(f"{n:>3} {k:>3} {nq:>3} {nS*env.dim:>8} {t:>5.1f} {str(frag):>10} {i:>4} | "
              f"{cm:>11.7f} {cf:>11.7f} {abs(cm-cf):>10.2e} | "
              f"{pm:>11.7f} {pf:>11.7f} {abs(pm-pf):>10.2e} | "
              f"{am:>10.7f} {af:>10.7f} {abs(am-af):>10.2e}")
        del r
    del m, HINT
P("-" * 160)
P(f"MAX |diff| single-record chi : {maxdiff_single:.3e}")
P(f"MAX |diff| pair chi          : {maxdiff_pair:.3e}")
P(f"MAX |diff| parity chi        : {maxdiff_par:.3e}")
ok_fast = max(maxdiff_single, maxdiff_pair, maxdiff_par) < 1e-9
P("SELF-CHECK: " + ("PASSED -- the fast path IS record_model on this coupling class."
                    if ok_fast else "FAILED -- setup broken, no conclusion may be drawn."))

# ------------------------------------------------------------------ parity is itself a record
P("")
P("=" * 108)
P("IS THE PARITY R_i R_j ITSELF A RECORD?  (D-18: nothing is called a record unchecked)")
P("=" * 108)
P(f"{'n':>3} {'k':>3} | {'pairs (i,j)':>12} {'(i) bit':>8} {'(ii) durable':>13} {'(iii) nontriv':>14} {'(iv) writable':>14}")
P("-" * 78)
for n in (4, 6, 8):
    car = carrier(n); k = car['k']
    H = code_hamiltonian(n); Rs = record_matrices(car); es = eigenspaces(H)
    npair = ci = cii = ciii = civ = 0
    for i, j in itertools.combinations(range(k), 2):
        Rp = Rs[i] @ Rs[j]; npair += 1
        ci += (np.linalg.norm(Rp - Rp.conj().T) < 1e-9 and
               np.linalg.norm(Rp @ Rp - np.eye(2 ** n)) < 1e-9)
        cii += (np.linalg.norm(H @ Rp - Rp @ H) < 1e-9)
        ciii += clause_iii(Rp, es)
        civ += clause_iv(Rp, es)
    P(f"{n:>3} {k:>3} | {npair:>12} {ci:>8} {cii:>13} {ciii:>14} {civ:>14}")
P("READ: every pair product R_i R_j satisfies (i)-(iv) too -- the parity of two records is a")
P("      record in its own right, so 'the environment holds the relation' is a statement about")
P("      a bona fide record, not about a derived statistic.")

# ------------------------------------------------------------------ THE INDEPENDENT-CARRIER CONTROL
P("")
P("=" * 108)
P("CONTROL, RUN THROUGH THE FULL MODEL: TWO LITERALLY INDEPENDENT CARRIERS, DISJOINT BATHS")
P("=" * 108)
P("System = [[4,2,2]] (x) [[4,2,2]], H = H_A (x) I + I (x) H_B, dim 256.  Record 1 lives on")
P("carrier A and couples ONLY to bath sites {0,1}; record 2 lives on carrier B and couples ONLY")
P("to bath sites {2,3}.  No relation between the two records can exist in any fragment.")
P("")
c4 = carrier(4)
H4 = code_hamiltonian(4); R4 = record_matrices(c4)
I16 = np.eye(16, dtype=complex)
Hc = np.kron(H4, I16) + np.kron(I16, H4)
RA = np.kron(R4[0], I16)
RB = np.kron(I16, R4[0])
esc = eigenspaces(Hc)
P(f"  composite eigenvalue multiplicities: {[m for _,_,m in esc]}")
P(f"  clause (ii) [Hc,RA]=0 : {np.linalg.norm(Hc@RA-RA@Hc):.2e}   [Hc,RB]=0 : {np.linalg.norm(Hc@RB-RB@Hc):.2e}")
P(f"  clause (iii) RA {clause_iii(RA, esc)}  RB {clause_iii(RB, esc)}   "
  f"clause (iv) RA {clause_iv(RA, esc)}  RB {clause_iv(RB, esc)}")
envc = Environment(nq=4, energies=ENERGIES[:4], beta=BETA)
mc = RecordModel(Hc)
Wc = np.array([[1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0]])
HINTc = build_HINT([RA, RB], Wc, envc)
P("")
P(f"{'t':>5} {'fragment':>10} | {'chi(R1)':>9} {'chi(R2)':>9} {'chi(R1,R2)':>11} {'SYNERGY':>10} "
  f"{'chi(parity)':>12} | READ")
P("-" * 92)
maxsyn = 0.0
poscontrol = 0.0
for t in (3.0, 9.0):
    rc = mc.evolve(HINTc, envc, lam=LAM, t=t)
    for frag in ([0], [2], [0, 2], [0, 1], [0, 1, 2, 3]):
        fr = None if len(frag) == 4 else frag
        c1 = envc.holevo(rc, RA, 256, fragment=fr)
        c2 = envc.holevo(rc, RB, 256, fragment=fr)
        pj = [(np.eye(256) + a * RA) @ (np.eye(256) + b * RB) / 4 for a in (+1, -1) for b in (+1, -1)]
        c12 = ref_chi_multi(envc, rc, pj, 256, fragment=fr)
        syn = c12 - c1 - c2
        cpar = envc.holevo(rc, RA @ RB, 256, fragment=fr)
        maxsyn = max(maxsyn, abs(syn))
        poscontrol = max(poscontrol, c1, c2)
        P(f"{t:>5.1f} {str(frag):>10} | {c1:>9.6f} {c2:>9.6f} {c12:>11.6f} {syn:>10.2e} "
          f"{cpar:>12.6f} | {'relation' if abs(syn)>1e-8 else 'no relation'}")
    del rc
P("-" * 92)
P(f"MAX |synergy| over independent carriers : {maxsyn:.3e}")
P(f"POSITIVE CONTROL, max single-record chi  : {poscontrol:.6f}   (non-zero: the instrument fires)")
P("READ: with k independent carriers on disjoint baths the synergy is " +
  ("ZERO to machine precision, while the same instrument registers a large single-record chi."
   if maxsyn < 1e-8 else "NOT zero -- the control is broken."))
P("      Note chi(parity) is NOT zero here: two independent bits still have a parity, and a")
P("      fragment straddling both baths knows it BECAUSE it knows both bits separately.  The")
P("      quantity that must vanish for 'no relation' is the SYNERGY, and it does.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY/s2_verify_fastpath.txt", "w").write(OUT.getvalue())
