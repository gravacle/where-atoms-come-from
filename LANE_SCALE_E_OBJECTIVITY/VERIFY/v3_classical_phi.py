"""V3 -- IS Phi AN ORDINARY FACT?  (attack axis 5)

Phi = 1 - (sum_i chi_i)/chi_ALL is the lane's ONLY quantity that grows monotonically with k.

THE ORDINARY EXPLANATION.  In the lane's model every bath site j is driven by the single scalar
    h_j(r) = lam * sum_i W[i,j] r_i
-- a LINEAR FUNCTIONAL of the record vector.  Site j therefore measures one coordinate of the
linear map r -> W^T r, not any record.  When k > nq that map is many-to-one, so no individual
r_i is recoverable while the collective is; that is precisely what Phi > 0 says.  Nothing in it
requires quantum mechanics, a bath, a Holevo quantity, or a record.

THE TEST.  Replace the whole quantum apparatus with the most boring classical channel that has
the same coupling matrix:
        y = W^T r + noise,   r uniform on {+-1}^k,  noise ~ N(0, sigma^2 I_nq)
and compute the SAME functional
        Phi_cl = 1 - (sum_i I(y ; r_i)) / I(y ; r).
If Phi_cl reproduces the lane's Phi curve -- 0 at k = 1, rising monotonically, exactly 0 in the
'separate' geometry -- then Phi measures the rank/overlap structure of W and nothing else.

CONTROL IN THE SAME TABLE (D-15): the 'separate' geometry, where W has disjoint support, must
give Phi_cl = 0 identically -- and the SAME code must give a large Phi_cl for 'crowded', or the
estimator is dead.
"""
import sys
import numpy as np

sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY")
from common import weights, sign_patterns

RNG = np.random.default_rng(101)
NS = 20000          # Monte-Carlo samples for the differential entropies


def mi_gaussian_mixture(W, sigma, k, nq, groups):
    """r uniform on {+-1}^k; y = W^T r + N(0,sigma^2 I).
       Returns I(y ; g) for each grouping g in `groups`, where a grouping is a list of record
       indices and the 'message' is the tuple of those signs (so groups=[[0],[1],..] gives the
       individual informations and groups=[[0..k-1]] gives the collective one).

       I(y;g) = h(y) - h(y|g).  Both are estimated with the SAME samples, so the shared
       differential-entropy offset cancels the way it does in the Holevo quantity."""
    S = sign_patterns(k).astype(float)            # (nP,k)
    nP = S.shape[0]
    means = S @ W                                  # (nP,nq)
    idx = RNG.integers(0, nP, size=NS)
    y = means[idx] + sigma * RNG.normal(size=(NS, nq))

    def logmix(sel_means, sel_w):
        # log p(y) for a uniform mixture of the listed component means
        d2 = ((y[:, None, :] - sel_means[None, :, :]) ** 2).sum(axis=2)
        a = -d2 / (2 * sigma ** 2) + np.log(sel_w)[None, :]
        m = a.max(axis=1, keepdims=True)
        return (m[:, 0] + np.log(np.exp(a - m).sum(axis=1)))

    w_all = np.full(nP, 1.0 / nP)
    lp_y = logmix(means, w_all)
    out = []
    for g in groups:
        g = list(g)
        # conditional on the signs of the records in g being those of the drawn pattern
        key = ((S[:, g] + 1) / 2).astype(int) @ (2 ** np.arange(len(g)))
        lp_cond = np.zeros(NS)
        for kv in np.unique(key):
            sel = np.where(key == kv)[0]
            rows = np.where(key[idx] == kv)[0]
            if len(rows) == 0:
                continue
            d2 = ((y[rows][:, None, :] - means[sel][None, :, :]) ** 2).sum(axis=2)
            a = -d2 / (2 * sigma ** 2) - np.log(len(sel))
            m = a.max(axis=1, keepdims=True)
            lp_cond[rows] = m[:, 0] + np.log(np.exp(a - m).sum(axis=1))
        out.append(float(np.mean(lp_cond - lp_y) / np.log(2)))     # bits
    return out


P = print
P("=" * 112)
P("V3  Phi FROM A PURELY CLASSICAL LINEAR-MIXING CHANNEL -- no bath, no Holevo, no records.")
P("    y = W^T r + N(0, sigma^2),  sigma = 1.0,  same W as the lane (col-normalised), nq = 6.")
P("=" * 112)
P(f"{'k':>3} | {'I(y;r) crd':>11} {'sum_i I crd':>12} {'Phi_cl crd':>11} | "
  f"{'I(y;r) SEP':>11} {'sum_i I SEP':>12} {'Phi_cl SEP':>11} | {'lane Phi crd (quantum)':>23}")
P("-" * 112)
LANE = {1: 0.000000, 2: 0.010593, 3: 0.097487, 4: 0.197543, 5: 0.305917,
        6: 0.470363, 7: 0.534806, 8: 0.608622, 9: 0.663847, 10: 0.694720}
NQ = 6
SIG = 1.0
for k in range(1, 11):
    row = []
    for kind in ('crowded', 'separate'):
        W = weights(kind, k, NQ, seed=7)
        groups = [[i] for i in range(k)] + [list(range(k))]
        vals = mi_gaussian_mixture(W, SIG, k, NQ, groups)
        ind = float(sum(vals[:k]))
        tot = float(vals[k])
        row.append((tot, ind, 1.0 - ind / tot if tot > 1e-9 else 0.0))
    P(f"{k:>3} | {row[0][0]:>11.5f} {row[0][1]:>12.5f} {row[0][2]:>11.5f} | "
      f"{row[1][0]:>11.5f} {row[1][1]:>12.5f} {row[1][2]:>11.5f} | {LANE[k]:>23.5f}")
P("-" * 112)
P("READ: the last two columns are the comparison.  If Phi_cl(crowded) tracks the lane's quantum")
P("      Phi and Phi_cl(separate) is ~0 (Monte-Carlo noise only), then Phi is a property of the")
P("      coupling matrix W -- a linear map from k inputs to nq outputs -- and carries no")
P("      information about records, quantum mechanics, or objectivity.")

# ---------------------------------------------------------------- sigma sweep
P("")
P("=" * 112)
P("SIGMA SWEEP -- the MAGNITUDE of Phi_cl is set by the channel noise, a free parameter.")
P("Control in the same table: 'separate' at each sigma, which must stay ~0.")
P("=" * 112)
P(f"{'sigma':>6} | " + " ".join(f"{'k='+str(k):>8}" for k in (2, 4, 6, 8, 10))
  + " | " + " ".join(f"{'SEPk='+str(k):>9}" for k in (2, 6, 10)))
P("-" * 112)
for sig in (0.35, 0.5, 0.75, 1.0, 1.5):
    a, b = [], []
    for k in (2, 4, 6, 8, 10):
        W = weights('crowded', k, NQ, seed=7)
        v = mi_gaussian_mixture(W, sig, k, NQ, [[i] for i in range(k)] + [list(range(k))])
        a.append(1.0 - sum(v[:k]) / v[k])
    for k in (2, 6, 10):
        W = weights('separate', k, NQ, seed=7)
        v = mi_gaussian_mixture(W, sig, k, NQ, [[i] for i in range(k)] + [list(range(k))])
        b.append(1.0 - sum(v[:k]) / v[k])
    P(f"{sig:>6.2f} | " + " ".join(f"{x:>8.4f}" for x in a) + " | " + " ".join(f"{x:>9.2e}" for x in b))
P("-" * 112)
P("READ (from the numbers above, not in advance): Phi_cl is positive and RISING WITH k at every")
P("      sigma, and is machine-zero in the disjoint control at every sigma.  Its MAGNITUDE tops out")
P("      near 0.21 here, BELOW the lane's 0.69 -- so the additive-Gaussian surrogate reproduces the")
P("      SHAPE and the exact control zero but NOT the size.  The claim supported is qualitative:")
P("      Phi > 0 requires only a many-to-one linear coupling map, not records and not quantum mechanics.")
