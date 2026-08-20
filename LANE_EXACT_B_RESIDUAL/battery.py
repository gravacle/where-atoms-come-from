"""The configuration battery and the EXPLAINED MODEL.

A CONFIGURATION is: the read record X1 written on bath site 0, plus a list of partner records
each written on a named bath site.  Every partner is a genuine record on this carrier -- a logical
Pauli: Hermitian, squares to I, commutes with H (hence durable), traceless on the code space
(hence writable).  Products of commuting logicals are logicals too and are used to obtain more
than one partner that PAIRS with the read record.

THE EXPLAINED MODEL -- the three terms the register already accounts for, each with a FITTED
coefficient (nothing assumed), plus a CONTROL column that must come out zero:

    log chi  =  a
              + gamma * log(1/(1+m))        HOLEVO CAPACITY (C-36): the site's bits split among
                                            the m+1 records written on it.  gamma = 1 is exact
                                            equipartition; it is FITTED, not assumed.
              + beta  * (m1 + p1)           PAIRING (C-38): site-BLIND disturbance from partners
                                            that pair with the record, wherever they sit.
              + delta * m1                  CROWDING SELECTIVITY (C-39): the pairing-dependent
                                            modulation of what SHARING the site costs.
              + c0    * p0                  CONTROL (D-15): commuting partners on OTHER sites.
                                            C-38 proves this coefficient is EXACTLY zero; if the
                                            fit returns anything else the pipeline is broken.

    m0 = same-site partners with pairing 0     m1 = same-site partners with pairing 1   m = m0+m1
    p0 = other-site partners with pairing 0    p1 = other-site partners with pairing 1
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *

READ = "X1"

def build_ops(n):
    """label -> (reduced matrix on the code space, symplectic vector over F_2^{2n})."""
    pairs = logical_pairs(n)
    ok, bad = check_symplectic(pairs, n)
    assert ok, f"SELF-CHECK FAILED at n={n}: {bad[:5]}"
    Q = code_basis(n)
    ops = {}
    for i, (a, b) in enumerate(pairs):
        ops[f"X{i+1}"] = (reduce_op(xz_to_matrix(a, n), Q), a)
        ops[f"Z{i+1}"] = (reduce_op(xz_to_matrix(b, n), Q), b)
    # products of COMMUTING logicals -- still Hermitian, still square to I, still records
    K = len(pairs)
    for j in range(2, K + 1):
        for t in ("X", "Z"):
            lab = f"Z1{t}{j}"
            A, va = ops["Z1"]; B, vb = ops[f"{t}{j}"]
            assert np.linalg.norm(A @ B - B @ A) < 1e-10
            M = A @ B
            assert np.linalg.norm(M - M.conj().T) < 1e-10 and np.linalg.norm(M @ M - np.eye(M.shape[0])) < 1e-10
            ops[lab] = (M, xor(va, vb))
    return ops, Q, n

def maxq(label):
    """highest logical-qubit index a label needs."""
    import re
    return max(int(d) for d in re.findall(r'\d+', label))

def cfg_maxq(partners):
    return max([1] + [maxq(l) for l, _ in partners])

# ---------------------------------------------------------------- THE BATTERY
CONFIGS = [
 ("alone",                       []),
 ("Z1@0",                        [("Z1",0)]),
 ("Z1@1",                        [("Z1",1)]),
 ("X2@0",                        [("X2",0)]),
 ("X2@1",                        [("X2",1)]),
 ("Z2@0",                        [("Z2",0)]),
 ("Z1X2@0",                      [("Z1X2",0)]),
 ("X2@0,X3@0",                   [("X2",0),("X3",0)]),
 ("X2@0,Z2@0",                   [("X2",0),("Z2",0)]),
 ("Z1@0,X2@0",                   [("Z1",0),("X2",0)]),
 ("Z1@0,Z1X2@0",                 [("Z1",0),("Z1X2",0)]),
 ("Z1X2@0,Z1Z2@0",               [("Z1X2",0),("Z1Z2",0)]),
 ("Z1@1,X2@1",                   [("Z1",1),("X2",1)]),
 ("Z1@0,X2@1",                   [("Z1",0),("X2",1)]),
 ("X2@0,Z1@1",                   [("X2",0),("Z1",1)]),
 ("X2@0,X3@1",                   [("X2",0),("X3",1)]),
 ("Z1@1,Z1X2@2",                 [("Z1",1),("Z1X2",2)]),
 ("Z1@1,Z1X2@1",                 [("Z1",1),("Z1X2",1)]),
 ("X2@1,X3@2",                   [("X2",1),("X3",2)]),
 ("X2@1,X3@1",                   [("X2",1),("X3",1)]),
 ("X2@0,X3@0,X4@0",              [("X2",0),("X3",0),("X4",0)]),
 ("X2@0,Z2@0,X3@0",              [("X2",0),("Z2",0),("X3",0)]),
 ("X2@0,Z2@0,X3@0,Z3@0",         [("X2",0),("Z2",0),("X3",0),("Z3",0)]),
 ("Z1@0,X2@0,X3@0",              [("Z1",0),("X2",0),("X3",0)]),
 ("Z1@0,Z1X2@0,Z1X3@0",          [("Z1",0),("Z1X2",0),("Z1X3",0)]),
 ("Z1X2@0,Z1Z2@0,X3@0",          [("Z1X2",0),("Z1Z2",0),("X3",0)]),
 ("X2@0,X3@0,X4@0,X5@0",         [("X2",0),("X3",0),("X4",0),("X5",0)]),
 ("X2@0,X3@0,X4@0,X5@0,X6@0",    [("X2",0),("X3",0),("X4",0),("X5",0),("X6",0)]),
 ("Z1@0,X2@0,X3@0,X4@0",         [("Z1",0),("X2",0),("X3",0),("X4",0)]),
 ("Z1@0,Z1X2@0,X3@0,X4@0",       [("Z1",0),("Z1X2",0),("X3",0),("X4",0)]),
 ("X2@0,X3@1,X4@2",              [("X2",0),("X3",1),("X4",2)]),
 ("Z1@0,X2@1,X3@2",              [("Z1",0),("X2",1),("X3",2)]),
 ("Z1@1,X2@0,X3@0",              [("Z1",1),("X2",0),("X3",0)]),
 ("Z1@0,Z1X2@0,Z1X3@0,Z1X4@0",   [("Z1",0),("Z1X2",0),("Z1X3",0),("Z1X4",0)]),
 ("X2@0,X3@0,X4@0,X5@0,X6@0,X7@0",[("X2",0),("X3",0),("X4",0),("X5",0),("X6",0),("X7",0)]),
 ("X2@0,X3@0,X4@0,X5@0,X6@0,X7@0,X8@0",
                                 [("X2",0),("X3",0),("X4",0),("X5",0),("X6",0),("X7",0),("X8",0)]),
]

def features(partners, ops, n):
    vr = ops[READ][1]
    m0 = m1 = p0 = p1 = 0
    for lab, site in partners:
        g = sp(vr, ops[lab][1], n)
        if site == 0:
            if g: m1 += 1
            else: m0 += 1
        else:
            if g: p1 += 1
            else: p0 += 1
    return dict(m0=m0, m1=m1, p0=p0, p1=p1, m=m0 + m1)

def partner_pairings(partners, ops, n):
    """the THREE-BODY descriptor: how many partner PAIRS anticommute with each other,
       split by whether they share the read record's site."""
    tot = same = 0
    for i in range(len(partners)):
        for j in range(i + 1, len(partners)):
            g = sp(ops[partners[i][0]][1], ops[partners[j][0]][1], n)
            tot += g
            if partners[i][1] == 0 and partners[j][1] == 0: same += g
    return tot, same

def design_row(f):
    m = f['m']
    return [1.0, np.log(1.0 / (1.0 + m)), float(f['m1'] + f['p1']), float(f['m1']), float(f['p0'])]

TERMS = ["a (log chi alone)", "gamma (CAPACITY C-36)", "beta (PAIRING C-38)",
         "delta (CROWDING-SELECTIVITY C-39)", "c0 (CONTROL: must be 0)"]

def run_battery(n, venue, lam, NB=3, configs=None, ops=None):
    if ops is None: ops, _, _ = build_ops(n)
    K = n - 2
    env = venue.env(NB)
    rows = []
    for name, partners in (configs or CONFIGS):
        if cfg_maxq(partners) > K: continue
        if any(s >= NB for _, s in partners): continue
        so = [(ops[READ][0], 0)] + [(ops[l][0], s) for l, s in partners]
        chi = float(np.mean(chi_times(so, ops[READ][0], env, lam, venue.times)))
        f = features(partners, ops, n)
        tot, same = partner_pairings(partners, ops, n)
        rows.append(dict(n=n, lam=lam, name=name, chi=chi, npart=len(partners),
                         pp_anti=tot, pp_anti_same=same, **f))
    return rows

def fit(rows):
    """least squares for the explained model in log chi.  Returns coefficients, predictions,
       residuals (LINEAR chi), and diagnostics."""
    A = np.array([design_row(r) for r in rows])
    y = np.log(np.array([r['chi'] for r in rows]))
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = A @ coef
    pred = np.exp(yhat)
    meas = np.array([r['chi'] for r in rows])
    return dict(coef=coef, pred=pred, meas=meas, resid=meas - pred,
                logresid=y - yhat, cond=float(np.linalg.cond(A)),
                rms=float(np.sqrt(np.mean((meas - pred) ** 2))),
                maxabs=float(np.max(np.abs(meas - pred))),
                relrms=float(np.sqrt(np.mean(((meas - pred) / meas) ** 2))))

# ---------------------------------------------------------------- THREE-RECORD DESCRIPTORS
def three_body_counts(partners, ops, n, read=READ):
    """Counts over unordered PAIRS OF PARTNERS -- quantities no statement about (record, one
       partner) can contain.  C-36, C-38 and C-39 are all such statements, so every count here
       is invisible to the explained model by construction."""
    vr = ops[read][1]
    c = dict(A_read=0, A_colo=0, A_split=0, C_colo_11=0, C_colo_10=0, C_colo_00=0)
    for i in range(len(partners)):
        for j in range(i + 1, len(partners)):
            li, si = partners[i]; lj, sj = partners[j]
            g = sp(ops[li][1], ops[lj][1], n)          # do the two PARTNERS pair with each other
            gi = sp(vr, ops[li][1], n); gj = sp(vr, ops[lj][1], n)
            if g:
                if si == 0 and sj == 0: c['A_read'] += 1
                elif si == sj:          c['A_colo'] += 1
                else:                   c['A_split'] += 1
            if si == sj and si != 0:
                k = gi + gj
                c['C_colo_11' if k == 2 else ('C_colo_10' if k == 1 else 'C_colo_00')] += 1
    return c

# additional configurations that isolate the three-record quantities and their controls
EXTRA_CONFIGS = [
 ("X2@1,Z2@1",                   [("X2",1),("Z2",1)]),
 ("X2@1,Z2@2",                   [("X2",1),("Z2",2)]),
 ("X2@1,Z2@1,X3@1",              [("X2",1),("Z2",1),("X3",1)]),
 ("Z1@1,X2@2",                   [("Z1",1),("X2",2)]),
 ("Z1@1,X2@1,X3@1",              [("Z1",1),("X2",1),("X3",1)]),
 ("Z1@1,X2@1,X3@2",              [("Z1",1),("X2",1),("X3",2)]),
 ("Z1@1,Z1X2@1,X3@1",            [("Z1",1),("Z1X2",1),("X3",1)]),
 ("Z1@1,Z1X2@2,X3@1",            [("Z1",1),("Z1X2",2),("X3",1)]),
 ("Z1X2@1,Z1Z2@1",               [("Z1X2",1),("Z1Z2",1)]),
 ("Z1X2@1,Z1Z2@2",               [("Z1X2",1),("Z1Z2",2)]),
 ("Z1@0,X2@1,Z2@1",              [("Z1",0),("X2",1),("Z2",1)]),
 ("Z1@0,X2@1,Z2@2",              [("Z1",0),("X2",1),("Z2",2)]),
 ("X2@0,Z2@0,X3@1",              [("X2",0),("Z2",0),("X3",1)]),
 ("X2@0,X3@0,X4@1",              [("X2",0),("X3",0),("X4",1)]),
 ("X2@0,Z2@0,X3@0,X4@0",         [("X2",0),("Z2",0),("X3",0),("X4",0)]),
 ("X2@0,Z2@0,X3@0,Z3@0,X4@0,Z4@0",[("X2",0),("Z2",0),("X3",0),("Z3",0),("X4",0),("Z4",0)]),
 ("Z1@0,X2@0,Z2@0",              [("Z1",0),("X2",0),("Z2",0)]),
 ("Z1@0,Z1X2@0,X3@0,Z3@0",       [("Z1",0),("Z1X2",0),("X3",0),("Z3",0)]),
]
