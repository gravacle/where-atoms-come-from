"""Shared readout helpers for s4 / s5 (kept separate so importing does not re-run s4's tables)."""
import numpy as np
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY")
from common import Broadcast, kron_sites, averages, entropies, TIMES


def parity_matrix(k):
    """PV[p, m] = product of r_i over i in the subset m, for the sign pattern p (Walsh matrix)."""
    p = np.arange(2 ** k)[:, None]
    m = np.arange(2 ** k)[None, :]
    x = p & m
    pc = np.zeros((2 ** k, 2 ** k), dtype=np.int64)
    for b in range(k):
        pc += (x >> b) & 1
    return (1 - 2 * (pc % 2)).astype(float)


def mask_rows(PV, specs):
    """specs: list of list of (subset_mask, sign).  mask = prod_x (1 + s_x PV[:,x]) / nP."""
    nP = PV.shape[0]
    rows = np.empty((len(specs), nP))
    for q, sp in enumerate(specs):
        v = np.ones(nP)
        for m, s in sp:
            v *= (1.0 + s * PV[:, m])
        rows[q] = v / nP
    return rows


def entropy_stack(B, sites, ti, rows):
    return entropies(averages(kron_sites(B.rho, sites, ti), rows))


def group_scan_chi(k, nq, W, lam, frags, times=TIMES):
    """time-averaged chi of EVERY one of the 2^k - 1 non-identity records, for each fragment.
       Returns (len(frags), 2^k) with column 0 (the identity element) left at 0."""
    PV = parity_matrix(k)
    specs = [[]]
    for m in range(1, 2 ** k):
        specs += [[(m, +1)], [(m, -1)]]
    rows = mask_rows(PV, specs)
    B = Broadcast(k, nq, W, lam, times=times)
    out = np.zeros((len(frags), 2 ** k))
    for fi, F in enumerate(frags):
        a = np.zeros(len(specs))
        for ti in range(len(times)):
            a += entropy_stack(B, F, ti, rows)
        a /= len(times)
        for m in range(1, 2 ** k):
            out[fi, m] = max(a[0] - 0.5 * (a[2 * m - 1] + a[2 * m]), 0.0)
    return out


def fwht_axis0(a, k):
    """Walsh-Hadamard transform along axis 0: out[m] = sum_p (-1)^popcount(m & p) a[p]."""
    D = a.shape[-1]
    x = a.reshape([2] * k + [D])
    for ax in range(k):
        x = np.moveaxis(x, ax, 0)
        x = np.stack([x[0] + x[1], x[0] - x[1]], axis=0)
        x = np.moveaxis(x, 0, ax)
    return x.reshape(2 ** k, D)


def group_scan_chi_fast(k, nq, W, lam, frags, times=TIMES, B=None):
    """Same output as group_scan_chi, via the Walsh transform.

    The conditional state for the record g_m is  (W_0 +- W_m) / nP  where W = FWHT(rho over sign
    patterns) -- because the conditioning mask on g_m = +-1 is exactly (1 +- PV[:,m]) / nP.  This
    turns an O(4^k d^2) mask product into O(2^k k d^2) and is what makes the whole-bath scan at
    k = 10 affordable.  Checked against the mask route in s6."""
    if B is None:
        B = Broadcast(k, nq, W, lam, times=times)
    nP = 2 ** k
    out = np.zeros((len(frags), nP))
    for fi, F in enumerate(frags):
        acc = np.zeros(nP)
        for ti in range(len(times)):
            rhoF = kron_sites(B.rho, F, ti)
            d = rhoF.shape[-1]
            Wc = fwht_axis0(rhoF.reshape(nP, d * d), k).reshape(nP, d, d) / nP
            base = Wc[0]
            st = np.concatenate([(base[None] + Wc[1:]), (base[None] - Wc[1:])], axis=0)
            S = entropies(np.concatenate([base[None], st], axis=0))
            acc[1:] += np.maximum(S[0] - 0.5 * (S[1:nP] + S[nP:]), 0.0)
        out[fi] = acc / len(times)
    return out
