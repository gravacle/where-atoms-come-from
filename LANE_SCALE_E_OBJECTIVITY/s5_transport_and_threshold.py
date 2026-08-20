"""S5 -- (a) the k = 1 baseline, (b) how the relational effect survives disorder, (c) TRANSPORT.

(a) DOES OBJECTIVITY REQUIRE MANY RECORDS?  k = 1 is not in the [[n,n-2,2]] family (n odd), but a
    one-record broadcast is exactly the k = 1 case of the same construction: couple ONE of the
    [[4,2,2]] records and leave the other uncoupled.  It is the baseline every k > 1 row must be
    read against.

(b) IS THE RELATIONAL EFFECT AN ARTIFACT OF AN EXACTLY DEGENERATE VENUE (D-17)?  The effect shows
    up when every record couples to a site with the SAME MAGNITUDE.  Two probes:
      SIGN disorder   W[i,j] = +-1/sqrt(k), random signs.  R and -R are the same record, so this
                      must NOT destroy the effect -- if it does, the effect is a phase artifact.
      MAGNITUDE disorder  W[i,j] = +-(1 + eta * xi_ij), xi ~ N(0,1).  eta = 0 is the degenerate
                      venue, eta -> large is the generic venue.  Sweeping eta says whether the
                      effect is a measure-zero coincidence or has a finite basin.

(c) TRANSPORT.  C-43's gauge transport A_h lives on D(G) and has no counterpart on a stabiliser
    carrier.  The honest counterpart HERE is the carrier's own admissible symmetry group: every
    qubit permutation pi commutes with X^(x)n and Z^(x)n, hence with H, so U_pi is admissible, and
    it MOVES the records -- it permutes the record group.  Transport is therefore applied exactly
    as in C-43 (an admissible operation that moves R) and the question is whether the environment
    still holds the transported record.
"""
import sys, io, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY")
from common import *
from s4_helpers import parity_matrix, mask_rows, entropy_stack, group_scan_chi

OUT = io.StringIO()
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.write(s + "\n")

LAM = 0.8
NQ = 6


def disordered_weights(k, nq, eta, seed=7, signs=True):
    rng = np.random.default_rng(seed)
    mag = 1.0 + eta * rng.normal(size=(k, nq))
    sg = rng.choice([-1.0, 1.0], size=(k, nq)) if signs else np.ones((k, nq))
    W = sg * mag
    return W / np.sqrt((W ** 2).sum(axis=0, keepdims=True))


def profile(k, nq, W, lam, frags, times=TIMES):
    """group-scan chi plus the redundancy readouts, from one broadcast."""
    chi = group_scan_chi(k, nq, W, lam, frags, times=times)
    depth = np.array([bin(m).count('1') for m in range(2 ** k)])
    return chi, depth


P("=" * 130)
P("S5  BASELINE, DISORDER, TRANSPORT")
P("=" * 130)

# ------------------------------------------------------------------ (a) k = 1 baseline
P("")
P("=" * 130)
P("TABLE G  THE k = 1 BASELINE.  One coupled record, same bath, same total coupling per site.")
P("          N1(0.1) = # single sites holding >= 0.9 of the whole-bath chi;  R(0.1) = nq/f(0.1).")
P("=" * 130)
P(f"{'k_coupled':>9} {'nq':>3} {'lam':>5} | {'chi_whole':>9} {'chi_site':>9} {'N1(.1)':>7} "
  f"{'f(.1)':>6} {'R(.1)':>6} | {'N1(.25)':>8} {'R(.25)':>7}")
P("-" * 84)
for lam in (0.4, 0.8, 1.2):
    for k in (1, 2, 4, 6, 8, 10):
        nq = NQ
        W = weights('crowded', k, nq) if k > 1 else np.ones((1, nq))
        B = Broadcast(k, nq, W, lam)
        allsub = [tuple(c) for f in range(1, nq + 1) for c in itertools.combinations(range(nq), f)]
        specs = [(f"R{i}", [(0.5, {i: +1}), (0.5, {i: -1})]) for i in range(k)]
        chi = {s: np.zeros(k) for s in allsub}
        for ti in range(len(TIMES)):
            for s in allsub:
                r = chi_batch(B, list(s), ti, specs)
                chi[s] += np.array([r[f"R{i}"] for i in range(k)])
        for s in allsub: chi[s] /= len(TIMES)
        whole = tuple(range(nq)); cw = chi[whole]
        csite = np.array([chi[(j,)] for j in range(nq)])
        line = [f"{k:>9} {nq:>3} {lam:>5.1f} | {cw.mean():>9.5f} {csite.mean():>9.5f}"]
        for d in (0.10, 0.25):
            thr = (1 - d) * cw
            N1 = (csite >= thr[None, :]).sum(axis=0).mean()
            fd = []
            for i in range(k):
                got = nq
                for f in range(1, nq + 1):
                    av = float(np.mean([chi[s][i] for s in allsub if len(s) == f]))
                    if av >= thr[i]: got = f; break
                fd.append(got)
            fd = np.array(fd, dtype=float)
            if d == 0.10:
                line.append(f" {N1:>7.3f} {fd.mean():>6.3f} {float(np.mean(nq/fd)):>6.3f} |")
            else:
                line.append(f" {N1:>8.3f} {float(np.mean(nq/fd)):>7.3f}")
        P("".join(line))
    P("-" * 84)
P("READ: k = 1 is the baseline.  Objectivity (many fragments each holding the record) is either")
P("      present already at k = 1 or it is not; whether it needs MANY records is read off this")
P("      column against the k > 1 rows.")

# ------------------------------------------------------------------ (b) disorder sweep
P("")
P("=" * 130)
P("TABLE H  DISORDER SWEEP (D-17).  EXCESS = chi_star - chi_depth1 and depth* of the best-known")
P("          record, as the coupling magnitudes go from exactly equal (eta = 0) to generic.")
P("          Sign-randomised throughout, because R and -R are the SAME record.")
P("          Control column: the 'separate' geometry at the same k and bath.")
P("=" * 130)
P(f"{'k':>3} {'eta':>6} {'seed':>5} | {'site: chi_d1':>12} {'chi_star':>9} {'EXCESS':>9} {'depth*':>7} "
  f"{'Nbetter':>8} | {'whole: chi_d1':>13} {'chi_star':>9} {'EXCESS':>9} {'depth*':>7} | "
  f"{'sep EXCESS':>10} {'sep depth*':>10}")
P("-" * 150)
ETAS = [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, None]     # None = fully generic (record_model 'crowded')
for k in (2, 4, 6, 8):
    frags = [[j] for j in range(NQ)] + [list(range(NQ))]
    Ws = weights('separate', k, NQ)
    chs, depth = profile(k, NQ, Ws, LAM, frags)
    sd1 = chs[:, depth == 1].max(axis=1); sst = chs[:, 1:].max(axis=1)
    sds = depth[1 + chs[:, 1:].argmax(axis=1)]
    for eta in ETAS:
        for seed in (7, 11):
            W = weights('crowded', k, NQ, seed=seed) if eta is None else disordered_weights(k, NQ, eta, seed=seed)
            chi, _ = profile(k, NQ, W, LAM, frags)
            site = slice(0, NQ); wh = slice(NQ, NQ + 1)
            d1 = chi[:, depth == 1].max(axis=1); st = chi[:, 1:].max(axis=1)
            ds = depth[1 + chi[:, 1:].argmax(axis=1)]
            nbet = np.array([int((chi[q, 1:] > d1[q] + 1e-12).sum()) for q in range(chi.shape[0])])
            P(f"{k:>3} {('gen' if eta is None else f'{eta:.2f}'):>6} {seed:>5} | "
              f"{d1[site].mean():>12.5f} {st[site].mean():>9.5f} {(st-d1)[site].mean():>9.6f} "
              f"{ds[site].mean():>7.3f} {nbet[site].mean():>8.2f} | "
              f"{d1[wh].mean():>13.5f} {st[wh].mean():>9.5f} {(st-d1)[wh].mean():>9.6f} "
              f"{ds[wh].mean():>7.3f} | {(sst-sd1)[site].mean():>10.2e} {sds[site].mean():>10.3f}")
    P("-" * 150)
P("Nbetter = how many of the 2^k - 1 records in the group the fragment knows BETTER than it knows")
P("any single coupled record.  It is 0 exactly when the environment's information is record-shaped.")

# ------------------------------------------------------------------ (c) transport
P("")
P("=" * 130)
P("TABLE I  TRANSPORT BY AN ADMISSIBLE CARRIER SYMMETRY.  Every qubit permutation pi commutes")
P("          with X^(x)n and Z^(x)n, hence with H, so U_pi is admissible (O-4) and it MOVES the")
P("          records.  R'_i = U_pi^dagger R_i U_pi is again a record; we locate it in the record")
P("          group and read its chi straight off the group scan.")
P("=" * 130)


def induced_map(n, perm):
    """Express each permuted record in the record basis (mod the stabiliser Sz, which acts as the
       identity on the code space).  Returns the list of subset masks S(i) with R'_i = g_{S(i)}."""
    car = carrier(n); k = car['k']
    Z = [np.array(v[n:], dtype=int) for v in car['recs_xz']]      # Z-parts (records are Z-type)
    one = np.ones(n, dtype=int)
    basis = Z + [one]                                             # spans the even-weight subspace
    Mb = np.array(basis, dtype=int).T % 2                         # n x (k+1)
    out = []
    for i in range(k):
        tgt = np.zeros(n, dtype=int)
        for q in range(n): tgt[perm[q]] = Z[i][q]
        A = np.concatenate([Mb, tgt[:, None]], axis=1) % 2
        # gaussian elimination over F2
        A = A.copy(); r = 0; piv = []
        for c in range(k + 1):
            p = next((q for q in range(r, n) if A[q, c]), None)
            if p is None: continue
            A[[r, p]] = A[[p, r]]
            for q in range(n):
                if q != r and A[q, c]: A[q] = (A[q] + A[r]) % 2
            piv.append(c); r += 1
        if A[r:, -1].any(): out.append(None); continue            # not in the span: setup broken
        coef = np.zeros(k + 1, dtype=int)
        for q, c in enumerate(piv): coef[c] = A[q, -1]
        out.append(int(sum((1 << i2) for i2 in range(k) if coef[i2])))
    return out


P(f"{'n':>3} {'k':>3} {'geom':>9} {'permutation':>26} | {'depth of R_i after':>18} | "
  f"{'chi_site before':>15} {'chi_site after':>14} | {'chi_whole before':>16} {'chi_whole after':>15} | "
  f"{'N1(.1) before':>13} {'N1(.1) after':>12}")
P("-" * 168)
rng = np.random.default_rng(3)
for n in (6, 8, 10):
    k = n - 2
    perms = [tuple(range(1, n)) + (0,)]                                  # cyclic shift
    perms.append(tuple(int(q) for q in rng.permutation(n)))
    perms.append(tuple([1, 0] + list(range(2, n))))                      # a transposition
    frags = [[j] for j in range(NQ)] + [list(range(NQ))]
    for kind in ('crowded', 'separate'):
        W = weights(kind, k, NQ)
        chi, depth = profile(k, NQ, W, LAM, frags)
        cw = chi[NQ]                       # whole bath
        cs = chi[:NQ]                      # sites
        base = [1 << i for i in range(k)]
        thr_b = 0.9 * np.array([cw[m] for m in base])
        N1b = np.array([(cs[:, m] >= thr_b[i]).sum() for i, m in enumerate(base)], dtype=float)
        for perm in perms:
            img = induced_map(n, perm)
            if any(v is None or v == 0 for v in img):
                P(f"{n:>3} {k:>3} {kind:>9} {str(perm):>26} | SELF-CHECK FAILED: a transported "
                  f"record left the record group -- no conclusion")
                continue
            thr_a = 0.9 * np.array([cw[m] for m in img])
            N1a = np.array([(cs[:, m] >= thr_a[i]).sum() for i, m in enumerate(img)], dtype=float)
            dep = np.array([bin(m).count('1') for m in img], dtype=float)
            P(f"{n:>3} {k:>3} {kind:>9} {str(perm):>26} | {dep.mean():>18.3f} | "
              f"{np.mean([cs[:, m].mean() for m in base]):>15.5f} "
              f"{np.mean([cs[:, m].mean() for m in img]):>14.5f} | "
              f"{np.mean([cw[m] for m in base]):>16.5f} {np.mean([cw[m] for m in img]):>15.5f} | "
              f"{N1b.mean():>13.3f} {N1a.mean():>12.3f}")
    P("-" * 168)
P("READ: 'depth of R_i after' is the depth in the record group of the transported record.  If it")
P("      is 1 the permutation only relabels the records; if it is > 1 the transported record is a")
P("      relation among the coupled ones, and whether the environment still holds it is exactly")
P("      what the chi columns report.  The 'separate' rows are the control: there a single site")
P("      knows one coupled record and nothing of depth > 1, so transport out of depth 1 must")
P("      take the single-site chi to zero.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY/s5_transport_and_threshold.txt",
     "w").write(OUT.getvalue())
