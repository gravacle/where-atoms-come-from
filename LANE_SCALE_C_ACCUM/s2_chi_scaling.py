"""LANE_SCALE_C_ACCUM  --  script 2: does the ENVIRONMENT'S information accumulate with k?

Measured on [[n,n-2,2]] (k = n-2 records) against the CONTROL, a product of m independent
[[4,2,2]] blocks (k = 2m records).

THE EXACT REDUCTION, and why it is not a shortcut.
Records are logical operators; the formation coupling is  sum_i R_i (x) X_site(i), i.e. the bath
reads the records without moving them (the model's own remark: "under a commuting coupling the
record's value does not change").  H commutes with every R_i, so H and R_1..R_k share an
eigenbasis; in that basis EVERY system-side operator in H_tot is diagonal, so H_tot is block
diagonal in the system index and the joint state never leaves its block:

    H_tot = sum_a |a><a| (x) [ h_a I + H_B + lam * sum_i r_i(a) X_site(i) ]

The code space carries exactly 2^k joint sectors, each ONE-dimensional, labelled by the sign
vector r(a) in {+-1}^k.  So the whole chi computation collapses to bath-sized (2^nq x 2^nq)
evolutions, one per sign vector -- and the bath operator depends on the sign vector ONLY through
the per-site sums, which cuts it further to prod_j (k_j + 1) evolutions.

That reduction is EXACT, and it is CHECKED against RecordModel.evolve + Environment.holevo on
the real 2^n matrices at n = 4, 6, 8 and on the control at dim 256 before any of it is believed.

IT ALSO CARRIES A FINDING.  The reduced problem depends only on (k, the site map, lam, beta,
the bath) -- the CARRIER HAS DROPPED OUT.  The [[n,n-2,2]] family and the product-of-blocks
control must therefore give IDENTICAL chi at equal k.  That is verified on the full matrices,
not merely asserted, and it means chi cannot discriminate a collective carrier from k
independent ones.
"""
import sys, itertools, time, json, os
from math import comb, log2
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_C_ACCUM")
import numpy as np
from record_model import (RecordModel, Environment, symplectic_logicals, xz_to_matrix,
                          eigenspaces, clause_iii, clause_iv)
from s1_combinatorics import carrier_nn2, carrier_product, embed, sp2, f2_rank

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.append(s)

TIMES = np.linspace(1.0, 13.0, 25)          # TIME-AVERAGE chi: unitary evolution recurs
BETA = 2.0

def energies_for(nq):
    """deterministic, non-degenerate bath energies; nq<=3 reproduces the model's own default"""
    base = [1.0, 1.4, 0.7]
    if nq <= 3: return tuple(base[:nq])
    out = list(base)
    for j in range(3, nq):
        out.append(round(0.7 + 0.9 * ((j * 0.6180339887) % 1.0), 6))
    return tuple(out)

def vN(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

# ================================================================= REDUCED ENGINE
def reduced_factorised(k, sites, nq, lam, beta=BETA, times=TIMES):
    """The ONE-RECORD-PER-SITE case, done site by site.  H_B is a sum of single-qubit terms and
       each record couples to its own qubit, so the bath state is a product across qubits for
       EVERY sign vector; entropies add and the sign bits are independent, hence
       chi_joint = sum_j chi_j exactly.  This is not an approximation -- it is checked against
       the general engine below -- and it is what makes a bath of nq = k qubits reachable at all
       (the general engine would need a 2^k x 2^k bath)."""
    assert len(set(sites)) == len(sites) == k, "factorised path needs one record per site"
    E = energies_for(nq)
    chi = []
    for j in sites:
        env1 = Environment(nq=1, energies=(E[j],), beta=beta)   # THIS site's own energy
        rth = env1.thermal()
        vals = []
        for t in times:
            sig = {}
            for s in (+1, -1):
                K = env1.HB + lam * s * env1.site[0]
                w_, V_ = np.linalg.eigh(K)
                ph = np.exp(-1j * w_ * t)
                M = V_.conj().T @ rth @ V_
                sig[s] = V_ @ (ph[:, None] * M * ph.conj()[None, :]) @ V_.conj().T
            av = 0.5 * sig[+1] + 0.5 * sig[-1]
            vals.append(max(vN(av) - 0.5 * vN(sig[+1]) - 0.5 * vN(sig[-1]), 0.0))
        chi.append(float(np.mean(vals)))
    return dict(chi_joint=float(sum(chi)), chi_each=chi, chi_sum=float(sum(chi)),
                chi_site={j: chi[i] for i, j in enumerate(sites)},
                kj=[1] * nq, cap=float(nq), factorised=True)

def reduced(k, sites, nq, lam, beta=BETA, times=TIMES):
    """EXACT chi for k records read by a bath of nq qubits, record i coupling to site sites[i].
       Returns time-averaged chi_joint, per-record chi, and their sum."""
    env = Environment(nq=nq, energies=energies_for(nq), beta=beta)
    rth = env.thermal()
    kj = [sum(1 for s in sites if s == j) for j in range(nq)]
    ranges = [range(kj[j] + 1) for j in range(nq)]
    uvecs = list(itertools.product(*ranges))
    # bath generator for each per-site up-count vector, and its evolution
    W = {}; SIG = {}
    for u in uvecs:
        K = env.HB + lam * sum((2 * u[j] - kj[j]) * env.site[j] for j in range(nq))
        w_, V_ = np.linalg.eigh(K)
        SIG[u] = (w_, V_)
        W[u] = np.prod([comb(kj[j], u[j]) / 2.0 ** kj[j] for j in range(nq)])
    tot = sum(W.values())
    assert abs(tot - 1.0) < 1e-12, f"weights must sum to 1, got {tot}"
    # conditional weights given record i (at site j0) is +1 / -1
    def cond_w(j0, sign):
        out = {}
        kk = kj[j0]
        for u in uvecs:
            u0 = u[j0]
            if sign > 0:
                c = comb(kk - 1, u0 - 1) / 2.0 ** (kk - 1) if 0 <= u0 - 1 <= kk - 1 else 0.0
            else:
                c = comb(kk - 1, u0) / 2.0 ** (kk - 1) if 0 <= u0 <= kk - 1 else 0.0
            out[u] = c * np.prod([comb(kj[j], u[j]) / 2.0 ** kj[j] for j in range(nq) if j != j0])
        return out
    condW = {}
    for j0 in range(nq):
        if kj[j0] == 0: continue
        condW[j0] = (cond_w(j0, +1), cond_w(j0, -1))
        for cw in condW[j0]:
            s_ = sum(cw.values())
            assert abs(s_ - 1.0) < 1e-10, f"conditional weights must sum to 1, got {s_}"
    cj_t, cs_t = [], {j: [] for j in condW}
    for t in times:
        sig = {}
        for u in uvecs:
            w_, V_ = SIG[u]
            ph = np.exp(-1j * w_ * t)
            M = V_.conj().T @ rth @ V_
            sig[u] = V_ @ (ph[:, None] * M * ph.conj()[None, :]) @ V_.conj().T
        rbar = sum(W[u] * sig[u] for u in uvecs)
        Sbar = vN(rbar)
        cj_t.append(max(Sbar - sum(W[u] * vN(sig[u]) for u in uvecs), 0.0))
        for j0, (wp, wm) in condW.items():
            rp = sum(wp[u] * sig[u] for u in uvecs); rm = sum(wm[u] * sig[u] for u in uvecs)
            cs_t[j0].append(max(Sbar - 0.5 * vN(rp) - 0.5 * vN(rm), 0.0))
    chi_joint = float(np.mean(cj_t))
    chi_site = {j: float(np.mean(v)) for j, v in cs_t.items()}
    chi_each = [chi_site[s] for s in sites]
    return dict(chi_joint=chi_joint, chi_each=chi_each, chi_sum=float(sum(chi_each)),
                chi_site=chi_site, kj=kj, cap=float(nq))

# ================================================================= FULL-MATRIX REFERENCE
def build_carrier(kind, par):
    """returns (n, S, records as F2 vectors, partners as F2 vectors)"""
    if kind == "family":
        n, S = carrier_nn2(par)
        pr = symplectic_logicals([s[:] for s in S], n)
        return n, S, [p[1] for p in pr], [p[0] for p in pr]
    n, S = carrier_product(par)
    pr4 = symplectic_logicals([s[:] for s in carrier_product(1)[1]], 4)
    R = [embed(p[1], 4, 4 * b, n) for b in range(par) for p in pr4]
    W = [embed(p[0], 4, 4 * b, n) for b in range(par) for p in pr4]
    return n, S, R, W

def full_reference(kind, par, nq, lam, times, checks):
    """chi computed the long way: real 2^n matrices, RecordModel.evolve, Environment.holevo."""
    n, S, Rv, Wv = build_carrier(kind, par)
    k = len(Rv)
    Smat = [xz_to_matrix(s, n) for s in S]
    H = -sum(Smat)
    R = [xz_to_matrix(v, n) for v in Rv]
    Wc = [xz_to_matrix(v, n) for v in Wv]
    dim = 2 ** n
    es = eigenspaces(H)
    tag = f"{kind}:{par}"
    # ---- CLAUSES (i)-(iv) on the actual matrices  (D-18)
    for i, Rm in enumerate(R):
        checks.append((tag, f"R{i} clause (i)  R=R+ and R^2=I",
                       np.linalg.norm(Rm - Rm.conj().T) < 1e-9
                       and np.linalg.norm(Rm @ Rm - np.eye(dim)) < 1e-9, ""))
        checks.append((tag, f"R{i} clause (ii) [H,R]=0",
                       np.linalg.norm(H @ Rm - Rm @ H) < 1e-9, ""))
        checks.append((tag, f"R{i} clause (iii) non-constant on some eigenspace",
                       clause_iii(Rm, es), ""))
        checks.append((tag, f"R{i} clause (iv) Tr(P_E R)=0 on every eigenspace",
                       clause_iv(Rm, es), ""))
    Pg, gk = None, None
    w_, V_ = np.linalg.eigh(H)
    gk = int(np.sum(np.abs(w_ - w_[0]) < 1e-9))
    Q = V_[:, :gk]; Pg = Q @ Q.conj().T
    checks.append((tag, "ground (code) space has dimension 2^k", gk == 2 ** k, f"{gk} vs {2**k}"))
    # ---- joint sectors: must be 2^k of them, each of dimension dim/2^k, summing to dim
    groups = {(): np.eye(dim, dtype=complex)}
    cols = {(): V_}
    basis = {(): np.eye(dim, dtype=complex)}
    cur = {(): np.eye(dim, dtype=complex)}
    for Rm in R:
        nxt = {}
        for lab, C in cur.items():
            Rs = C.conj().T @ Rm @ C
            ws, Vs = np.linalg.eigh(Rs)
            for s in (+1, -1):
                idx = [i for i in range(len(ws)) if (ws[i] > 0) == (s > 0)]
                if idx: nxt[lab + (s,)] = C @ Vs[:, idx]
        cur = nxt
    dims = [C.shape[1] for C in cur.values()]
    checks.append((tag, "joint record sectors: 2^k of them, dims sum to dim",
                   len(cur) == 2 ** k and sum(dims) == dim and len(set(dims)) == 1,
                   f"{len(cur)} sectors, dims {sorted(set(dims))}, sum {sum(dims)} vs {dim}"))
    # ---- evolve on the FULL space and read chi the model's way
    env = Environment(nq=nq, energies=energies_for(nq), beta=BETA)
    sites = [i % nq for i in range(k)]
    coupling = [(R[i], sites[i]) for i in range(k)]
    nB = env.dim
    HINT = sum(np.kron(A, env.site[j]) for A, j in coupling)
    Ht = np.kron(H, np.eye(nB)) + np.kron(np.eye(dim), env.HB) + lam * HINT
    wt, Ut = np.linalg.eigh(Ht)
    r0 = np.kron(Pg / gk, env.thermal())
    Uc = Ut.conj().T @ r0 @ Ut
    # For the JOINT chi, change to the joint record basis ONCE instead of applying 2^k
    # projectors per time step: the sectors then become contiguous index ranges.
    labs = list(cur.keys())
    Vall = np.concatenate([cur[l] for l in labs], axis=1)
    CHECKS_local = np.linalg.norm(Vall.conj().T @ Vall - np.eye(dim))
    checks.append((tag, "joint record basis is unitary", CHECKS_local < 1e-8,
                   f"||V+V - I|| = {CHECKS_local:.2e}"))
    bounds = np.cumsum([0] + [cur[l].shape[1] for l in labs])
    # CACHE only the expensive full-space evolution; every CHECK above is recomputed every run.
    ckey = f"{kind}|{par}|{nq}|{lam}|{len(times)}|{times[0]}|{times[-1]}"
    cpath = "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_C_ACCUM/s2_fullref_cache.json"
    cache = {}
    if os.path.exists(cpath):
        try: cache = json.load(open(cpath))
        except Exception: cache = {}
    if ckey in cache:
        c = cache[ckey]
        return dict(k=k, n=n, dim=dim, chi_joint=c["cj"], chi_each=c["ce"],
                    chi_partner=c["cw"], chi_sum=float(sum(c["ce"])), sites=sites, cached=True)
    Utj = np.kron(Vall, np.eye(nB)).conj().T @ Ut          # done once
    cj, ce, cw = [], [[] for _ in range(k)], [[] for _ in range(k)]
    for t in times:
        ph = np.exp(-1j * wt * t)
        r = Ut @ (ph[:, None] * Uc * ph.conj()[None, :]) @ Ut.conj().T
        for i in range(k):
            ce[i].append(env.holevo(r, R[i], dim))
            cw[i].append(env.holevo(r, Wc[i], dim))          # D-15 POSITIVE/NEGATIVE CONTROL
        rj = Utj @ (ph[:, None] * Uc * ph.conj()[None, :]) @ Utj.conj().T
        T = rj.reshape(dim, nB, dim, nB)
        blocks = []
        for s in range(len(labs)):
            a0, a1 = bounds[s], bounds[s + 1]
            rb = T[a0:a1, :, a0:a1, :].trace(axis1=0, axis2=2)
            p = float(np.real(np.trace(rb)))
            if p < 1e-12: continue
            blocks.append((p, rb / p))
        ptot = sum(p for p, _ in blocks)
        assert abs(ptot - 1.0) < 1e-8, f"sector probabilities must sum to 1, got {ptot}"
        av = sum(p * rb for p, rb in blocks)
        cj.append(max(vN(av) - sum(p * vN(rb) for p, rb in blocks), 0.0))
    res = dict(k=k, n=n, dim=dim, chi_joint=float(np.mean(cj)),
               chi_each=[float(np.mean(c)) for c in ce],
               chi_partner=[float(np.mean(c)) for c in cw],
               chi_sum=float(sum(np.mean(c) for c in ce)), sites=sites, cached=False)
    cache[ckey] = dict(cj=res["chi_joint"], ce=res["chi_each"], cw=res["chi_partner"])
    json.dump(cache, open(cpath, "w"), indent=1)
    return res

# ================================================================= RUN
def run():
    CHECKS = []
    P("=" * 118)
    P("LANE_SCALE_C_ACCUM  script 2  --  does the environment's information ACCUMULATE with k?")
    P("=" * 118)

    # ---------------------------------------------------------------- VALIDATION
    P("")
    P("VALIDATION 1  --  the reduced engine against the FULL 2^n matrices (RecordModel.evolve +")
    P("                  Environment.holevo).  Nothing below is believed until this agrees.")
    P("-" * 118)
    TV = np.linspace(1.0, 13.0, 9)          # fewer times: the full route is the expensive one
    P(f"{'carrier':>14} {'dim':>6} {'k':>3} {'nq':>3} | {'chi_joint FULL':>14} {'chi_joint RED':>14}"
      f" {'|diff|':>10} | {'chi_sum FULL':>13} {'chi_sum RED':>12} {'|diff|':>10} |"
      f" {'chi(partner) FULL':>18}")
    ref = {}
    for kind, par in [("family", 4), ("family", 6), ("family", 8), ("control", 1), ("control", 2)]:
        t0 = time.time()
        f = full_reference(kind, par, nq=3, lam=0.8, times=TV, checks=CHECKS)
        rd = reduced(f["k"], f["sites"], 3, 0.8, times=TV)
        d1 = abs(f["chi_joint"] - rd["chi_joint"]); d2 = abs(f["chi_sum"] - rd["chi_sum"])
        CHECKS.append((f"{kind}:{par}", "reduced engine reproduces the full-matrix chi_joint",
                       d1 < 1e-8, f"|diff| = {d1:.3e}"))
        CHECKS.append((f"{kind}:{par}", "reduced engine reproduces the full-matrix per-record chi",
                       d2 < 1e-8, f"|diff| = {d2:.3e}"))
        mp = max(f["chi_partner"])
        CHECKS.append((f"{kind}:{par}", "D-15 control: chi about the CONJUGATE partner is zero "
                       "while chi about the record is not",
                       mp < 1e-8 and min(f["chi_each"]) > 1e-3,
                       f"partner max {mp:.3e}, record min {min(f['chi_each']):.4f}"))
        name = f"[[{par},{par-2},2]]" if kind == "family" else f"[[4,2,2]]^{par}"
        P(f"{name:>14} {f['dim']:>6} {f['k']:>3} {3:>3} | {f['chi_joint']:>14.9f} "
          f"{rd['chi_joint']:>14.9f} {d1:>10.2e} | {f['chi_sum']:>13.9f} {rd['chi_sum']:>12.9f} "
          f"{d2:>10.2e} | {max(f['chi_partner']):>18.2e}   ({time.time()-t0:.0f}s)")
        ref[(kind, par)] = f
    P("")
    P("VALIDATION 2  --  CARRIER-BLINDNESS.  At equal k the collective carrier [[6,4,2]] and the")
    P("                  product control [[4,2,2]]^2 are DIFFERENT carriers (dim 64 vs 256, 2 vs 4")
    P("                  stabilisers).  Full-matrix chi at k=4, same bath, same site map:")
    P("-" * 118)
    a = ref[("family", 6)]; b = ref[("control", 2)]
    P(f"   [[6,4,2]]    dim {a['dim']:>4}  chi_joint {a['chi_joint']:.9f}  chi_sum {a['chi_sum']:.9f}")
    P(f"   [[4,2,2]]^2  dim {b['dim']:>4}  chi_joint {b['chi_joint']:.9f}  chi_sum {b['chi_sum']:.9f}")
    dj = abs(a["chi_joint"] - b["chi_joint"]); ds = abs(a["chi_sum"] - b["chi_sum"])
    P(f"   |difference| chi_joint {dj:.3e}   chi_sum {ds:.3e}")
    CHECKS.append(("k=4", "collective carrier and product control give the SAME chi",
                   dj < 1e-8 and ds < 1e-8, f"{dj:.2e}, {ds:.2e}"))
    P("   => at equal k the CONTROL COLUMN IS THE FAMILY COLUMN.  chi does not see the carrier;")
    P("      it sees only k, the site map and the bath.  Every chi table below therefore reports")
    P("      one column that IS both, and the discriminating control is the BATH-SCALING contrast.")

    P("")
    P("VALIDATION 3  --  the FACTORISED path (one record per bath qubit) against the general")
    P("                  reduced engine.  The growing-bath control needs a bath of k qubits, which")
    P("                  the general engine cannot build past k ~ 10; the factorised path can.")
    P("-" * 118)
    P(f"{'k':>3} | {'general chi_joint':>18} {'factorised chi_joint':>21} {'|diff|':>10}")
    for k in (2, 4, 6, 8):
        gg = reduced(k, list(range(k)), k, 0.8)
        ff = reduced_factorised(k, list(range(k)), k, 0.8)
        d = abs(gg["chi_joint"] - ff["chi_joint"])
        CHECKS.append((f"k={k}", "factorised path reproduces the general reduced engine", d < 1e-9,
                       f"|diff| = {d:.3e}"))
        P(f"{k:>3} | {gg['chi_joint']:>18.10f} {ff['chi_joint']:>21.10f} {d:>10.2e}")

    # ---------------------------------------------------------------- SELF-CHECKS
    P("")
    P("SELF-CHECKS  (a FAILING check voids every conclusion below it)")
    P("-" * 118)
    bad = [c for c in CHECKS if not c[2]]
    for c in bad: P(f"   FAIL  {c[0]:12s} {c[1]}   {c[3]}")
    P(f"   {len(CHECKS)-len(bad)} / {len(CHECKS)} checks pass"
      + ("   -- ALL PASS" if not bad else "   -- SOME FAILED"))
    if bad:
        P("   CONCLUSIONS VOID -- stopping.")
        return CHECKS

    # ---------------------------------------------------------------- TABLE A
    P("")
    P("TABLE A   FIXED ENVIRONMENT (nq bath qubits, held constant while k grows) versus the")
    P("          POSITIVE CONTROL of a GROWING environment (one bath qubit per record).")
    P("          Same engine, same lam, same beta, same time average -- only the venue's scale differs.")
    P("          lam = 0.8, beta = 2.0, chi averaged over 25 times in [1,13].")
    P("-" * 118)
    KS = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
    P(f"{'k':>3} | {'nq=1':>19} | {'nq=2':>19} | {'nq=3':>19} | {'nq=4':>19} || {'GROWING nq=k':>21}")
    P(f"{'':>3} | {'joint':>9}{'sum':>10} | {'joint':>9}{'sum':>10} | {'joint':>9}{'sum':>10} |"
      f" {'joint':>9}{'sum':>10} || {'joint':>10}{'sum':>11}")
    P("-" * 118)
    A = {}
    for k in KS:
        row = f"{k:>3} |"
        for nq in (1, 2, 3, 4):
            r = reduced(k, [i % nq for i in range(k)], nq, 0.8)
            A[(k, nq)] = r
            row += f" {r['chi_joint']:>9.5f}{r['chi_sum']:>10.5f} |"
        g = reduced_factorised(k, list(range(k)), k, 0.8)
        A[(k, "grow")] = g
        row += f"| {g['chi_joint']:>10.5f}{g['chi_sum']:>11.5f}"
        P(row)
    P("-" * 118)
    P("joint = chi the bath holds about the FULL k-bit record string;  sum = sum of the k individual chi")
    P(f"bath capacity ceiling log2(dim bath) = nq bits: 1, 2, 3, 4 respectively (and k for the growing bath)")

    # ---------------------------------------------------------------- TABLE B
    P("")
    P("TABLE B   IS THE WHOLE MORE OR LESS THAN THE SUM?  gap = chi_joint - sum_i chi_i,")
    P("          and the SPREAD/CROWDED contrast (D-16: divide by the SPREAD control, never by 'alone').")
    P("-" * 118)
    P(f"{'k':>3} | {'nq=3 joint':>11} {'nq=3 sum':>9} {'gap':>9} {'joint/cap':>10} {'chi/record':>11}"
      f" | {'crowded(nq=3,all site0)':>24} {'crowded/spread':>15} || {'grow joint':>11} {'grow gap':>9}")
    P("-" * 118)
    for k in KS:
        r = A[(k, 3)]; g = A[(k, "grow")]
        cr = reduced(k, [0] * k, 3, 0.8)
        ratio = cr["chi_sum"] / r["chi_sum"] if r["chi_sum"] > 1e-12 else float("nan")
        P(f"{k:>3} | {r['chi_joint']:>11.5f} {r['chi_sum']:>9.5f} "
          f"{r['chi_joint']-r['chi_sum']:>9.5f} {r['chi_joint']/3.0:>10.5f} "
          f"{r['chi_sum']/k:>11.5f} | {cr['chi_sum']:>24.5f} {ratio:>15.5f} || "
          f"{g['chi_joint']:>11.5f} {g['chi_joint']-g['chi_sum']:>9.2e}")
    P("-" * 118)

    # ---------------------------------------------------------------- TABLE C
    P("")
    P("TABLE C   D-17: VARY THE VENUE'S OWN SCALE.  Same measurement at three couplings.")
    P("-" * 118)
    P(f"{'k':>3} |" + "".join(f" {'lam='+str(l):>26} |" for l in (0.4, 0.8, 1.2)))
    P(f"{'':>3} |" + "".join(f" {'joint(nq3)':>12}{'sum(nq3)':>13} |" for _ in range(3)))
    for k in KS:
        row = f"{k:>3} |"
        for lam in (0.4, 0.8, 1.2):
            r = reduced(k, [i % 3 for i in range(k)], 3, lam)
            row += f" {r['chi_joint']:>12.5f}{r['chi_sum']:>13.5f} |"
        P(row)

    # ---------------------------------------------------------------- SCALING READ
    P("")
    P("SCALING READ   (computed from Table A, not written in advance)")
    P("-" * 118)
    def classify(ks, vals, name):
        K = np.array(ks, float); V = np.array(vals, float)
        if V.max() - V.min() < 1e-9:
            return f"{name}: CONSTANT at {V[0]:.5f}"
        A1 = np.vstack([K, np.ones_like(K)]).T
        s1 = np.linalg.lstsq(A1, V, rcond=None)[0]
        r1 = float(np.sqrt(np.mean((V - A1 @ s1) ** 2)))
        A2 = np.vstack([np.log2(K), np.ones_like(K)]).T
        s2 = np.linalg.lstsq(A2, V, rcond=None)[0]
        r2 = float(np.sqrt(np.mean((V - A2 @ s2) ** 2)))
        # saturating test: last-decade increment against first
        inc_first = V[1] - V[0]; inc_last = V[-1] - V[-2]
        sat = inc_last < 0.15 * abs(inc_first) if abs(inc_first) > 1e-9 else True
        lab = "SATURATING" if sat else ("LINEAR" if r1 < 0.02 * max(1e-9, V.max()) else "GROWING(sublinear)")
        return (f"{name}: {lab}  linear-slope={s1[0]:+.5f} rms_lin={r1:.5f} rms_log={r2:.5f} "
                f"first-step={inc_first:+.5f} last-step={inc_last:+.5f} "
                f"[k x{K[-1]/K[0]:.0f} -> value x{V[-1]/V[0]:.2f}]")
    for nq in (1, 2, 3, 4):
        P("  " + classify(KS, [A[(k, nq)]["chi_joint"] for k in KS], f"chi_joint, FIXED bath nq={nq}"))
        P("  " + classify(KS, [A[(k, nq)]["chi_sum"] for k in KS], f"chi_sum  , FIXED bath nq={nq}"))
    P("  " + classify(KS, [A[(k, "grow")]["chi_joint"] for k in KS], "chi_joint, GROWING bath nq=k"))
    P("  " + classify(KS, [A[(k, "grow")]["chi_sum"] for k in KS], "chi_sum  , GROWING bath nq=k"))
    P("  " + classify(KS, [A[(k, 3)]["chi_sum"] / k for k in KS], "chi PER RECORD, fixed bath nq=3"))

    P("")
    P("NOISE FLOOR.  The chi values are von Neumann entropies of matrices assembled from exact")
    P("eigendecompositions; the reduced engine reproduces the full 2^n route to better than 1e-9")
    P("(validation table).  Anything below 1e-9 is not distinguishable from zero here.")
    P(f"Largest k reached for chi: k = {KS[-1]} records.  What stopped it: nothing in this engine --")
    P("the reduction removed the 2^n cost.  What stopped the FULL-MATRIX route is stated in script 2's")
    P("validation table: n = 8 (dim 256, joint dim 2048) was the largest full evolution run here.")
    return CHECKS

if __name__ == "__main__":
    run()
    with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_C_ACCUM/s2_chi_scaling.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
