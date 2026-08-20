"""STEP 4 -- VALIDATE the code-space reduction and the fast chi against the model's own
RecordModel.evolve + Environment.holevo in the FULL 2^n space.  Nothing in step 5 is reported
unless every check here passes.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC")
from lib_scaleb import *
from lib_operational import *
import numpy as np

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

P("=" * 108)
P("LANE_SCALE_B_METRIC  STEP 4 -- VALIDATION OF THE CODE-SPACE REDUCTION")
P("=" * 108)
P("")
P("%-4s %-6s %-12s %-12s %-12s %-12s %-8s" %
  ("n", "k", "SC-8a |W'W-I|", "SC-8b |WW'-Pg|", "SC-9 Paulis", "SC-10 |W'HW+2|", "verdict"))
P("-" * 80)
iso = {}
for n in [4, 6, 8, 10]:
    stab, pairs = carrier(n)
    W, Xb, Zb = code_isometry(n, pairs)
    c = isometry_checks(W, n, pairs)
    ok = all(v < 1e-8 for v in c.values())
    iso[n] = (W, pairs, Xb, Zb)
    P("%-4d %-6d %-12.2e %-12.2e %-12.2e %-12.2e %-8s" %
      (n, len(pairs), c["SC-8a"], c["SC-8b"], c["SC-9"], c["SC-10"], "PASS" if ok else "FAIL"))

P("")
P("SC-11  reduced-space chi  vs  FULL-space chi from RecordModel.evolve + Environment.holevo")
P("       coupling = the record itself distributed to the bath sites its support touches")
P("")
P("%-4s %-6s %-6s %-8s %-14s %-14s %-12s %-8s" %
  ("n", "drive", "read", "t", "chi FULL", "chi REDUCED", "|diff|", "verdict"))
P("-" * 88)

ENERGIES3 = (1.0, 1.4, 0.7)
env3 = Environment(nq=3, energies=ENERGIES3, beta=2.0)
worst = 0.0
for n in [4, 6]:
    W, pairs, Xb, Zb = iso[n]
    k = len(pairs)
    vs, lab = record_vectors(pairs, n)
    Rfull = [xz_to_matrix(v, n) for v in vs]
    Heff = -2 * np.eye(2 ** k, dtype=complex)
    mdl_full = RecordModel(hamiltonian(n))
    for di in [0, k]:                                    # one X-record and one Z-record
        supp = sorted(support(vs[di], n))
        coup_full = [(Rfull[di], q) for q in supp]
        HINT_red = sum(np.kron(W.conj().T @ Rfull[di] @ W, env3.site[q % env3.nq]) for q in supp)
        state0 = np.eye(2 ** k, dtype=complex) / 2 ** k
        for t in [2.0, 5.5]:
            rho_full = mdl_full.evolve(coup_full, env3, lam=0.8, t=t)
            rho_red = next(evolve_cached(Heff, env3, HINT_red, 0.8, [t], state0))
            for ri in [di, (di + 1) % k + (k if di >= k else 0)]:
                q = ri % k; letter = 'X' if ri < k else 'Z'
                cf = env3.holevo(rho_full, Rfull[ri], 2 ** n)
                cr = chi_fast(rho_red, 2 ** k, env3.dim, k, q, letter)
                d = abs(cf - cr); worst = max(worst, d)
                P("%-4d %-6s %-6s %-8.1f %-14.9f %-14.9f %-12.2e %-8s" %
                  (n, lab[di], lab[ri], t, cf, cr, d, "PASS" if d < 1e-9 else "FAIL"))
P("")
P("worst |chi_FULL - chi_REDUCED| over all validation points = %.3e" % worst)
P("SC-11 %s" % ("PASS -- the reduction and the fast chi are exact" if worst < 1e-9 else "FAIL"))

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/04_operational_validate.txt", "w").write("\n".join(OUT) + "\n")
