"""C-71/C-72 FORMATION GATE -- THE HARDENED CONTRAST BLOCK (T-50 part b) -- 2026-08-21.

DROP-IN replacement for the C-71/C-72 formation section of model/validate_geometry.py
(everything between the "C-71/C-72: formation" marker and the "summary + chain" marker).
Same check() idiom.  The registrar installs; this lane writes only under LANE_T50/.

THE RULING (the principal, 2026-08-20): stop trying to rescue an absolute single-page
ratio.  The observable is the BASELINE-FREE DIFFERENTIAL -- the contrast between written
and unwritten sectors of the SAME part in the SAME read -- and its EXPONENT in block
size N.  This block gates the model's formation layer in exactly that vocabulary, with
NO LITERAL EXPECTED VALUE ON ANY DECISION PATH: every decision compares one quantity
measured in this run against another quantity measured in this run.  The claim's own
algebra supplies the only fixed numbers on structural decision paths -- the factor 2
(E|sum| grows as N for an accumulating contrast, as N^(1/2) for a screening one; the
threshold sits at the MIDPOINT of the two measured endpoints, beta_null and 2*beta_null,
i.e. excess > beta_null/2).

THE THREE COMPUTED CONSTRAINTS, honoured where each binds:
  1. The raw exponent is NOT baseline-free.  No raw single-sector statistic sits on any
     decision path.  The raw unwritten-sector exponent is computed as a SHADOW, printed,
     and never gated; the mutation harness asserts that the shadow FLIPS under an
     in-spec non-zero-mean residual while every decision check is invariant -- the
     invariance is measured, not assumed.
  2. Centred variance discriminates nothing.  The decision statistic is mean|D| --
     the first absolute moment E|Q| literally.  The uncentred second moment (slope = 2*beta) and the centred
     variance (slope ~ 1 for BOTH encodings) are printed as displays, never gated.
  3. The programmed density is held fixed in the observable's definition: blocks
     subsample fixed-density sectors, and check S6 gates the read's own measured
     per-block programmed fraction for an N-trend, with a failing branch demonstrated
     at rate 1.000 (mutation M3, refuter A's counterexample).

REPLACEMENT MAP (old check -> what stands in its place, with the reason):
  * "floor RHO(f) bounds every treated page"      -> DROPPED.  Reason: INST-14 -- the
    check is the triangle inequality given the generator's own premises; 200,000
    adversarial in-domain draws produced zero violations; it has no failing branch
    inside the declared physics, and the ruling removes the absolute ratio from every
    decision path.  Its structural content (occupancy's accumulation is N-free in
    kind, orientation's screening falls as N^(-1/2)) is gated by S1 and S4 in exponent
    form, where the two-signed-write mutation flips them at rate 1.000.
  * "DISCRIMINATOR at matched N"                  -> S4 (cross-encoding excess gap,
    exponent form; no matched-N requirement -- exponents are dimensionless).
  * "DISCRIMINATOR WIDENS with N"                 -> S1.  The old widens check had
    0.18 power (9/50 seeds) under the two-signed write; growth of the separation IS
    the excess exponent xi > beta_null/2, and S1 fails at 1.000 under that mutation.
  * control A (all-programmed page == RHO(1))     -> DROPPED, reason stated: INST-15 --
    RHO(1) == 1 identically for every Delta, and the all-programmed page has zero
    erased cells for the treatment to act on; verified incapable of failing at any
    tolerance (Delta = 1..99).  A comparison against the literal 1.0 with a function
    call interposed.  SLOT REPAIRED as A2: the f=1 endpoint read through the SAME
    contrast pipeline must ACCUMULATE (xi_f1 > beta_null/2) -- a decision boolean
    with a failing branch measured at rate 1.000 (mutation M1).
  * control B (unwritten page screens)            -> B2: the unwritten pool's
    cross-family pseudo-contrast must show NO excess (xi_pseudo < beta_null/2).
    Failing branch measured at rate 1.000 (mutation M2b: in-spec sector-differential
    residual -- the within-part null premise broken from inside spec).
  * control C (DC-saturated accumulates)          -> S3 (same content, exponent form,
    failing branch at 1.000 under mutation M5).
  * C-72 orientation literal checks (COMP-11: 0.4, 0.0027, +270, -2/+262, literal 1.0)
                                                  -> S2, S3, S5.  No sealed number is
    compared on any path; screening and accumulation are gated as measured exponents
    against the same read's own null.
  * C-71 "unwritten sum == -38 e within tolerance" -> DROPPED as a decision (the -38 is
    a literal; the tolerance half is |sum r| <= N*Delta for |r_i| <= Delta, an identity
    with no failing branch).  The null claim is gated by B2 and the S5 sign-consistency
    half, both with measured failing branches.

R8 REPOINTING MAP (for the registrar; instrument defect 18): gate cells naming the old
checks verbatim repoint as follows -- P-FORM-9 (the floor) -> S1 + S4 (the floor's
N-free-vs-N^(-1/2) content in exponent form); P-ROLES-2 -> S3 + S5.

DECLARED CONSTANTS ON DECISION PATHS (tolerances and analysis constants -- none is an
expected value of any measurement):
  GRID, SECTOR_LEN, K_PAIRS  -- the measurement geometry (8 grid points, 2.1 decades).
  DELTA_E = 5                -- the declared over-erase tolerance of the T-34 surface,
                                carried unchanged from the sealed block.
  FRAC = DELTA_E/N_E         -- the declared orientation read tolerance, as sealed.
  DENSITY_TOL = 0.05         -- S6's trend tolerance; measured into place by the
                                harness (see MARGINS in the sealed matrix: baseline
                                worst |slope| and mutated best |slope| bracket it by
                                more than an order of magnitude on each side).
  The structural factor 2 (midpoint threshold beta_null/2) is the claim's own algebra,
  not a tuned constant.  Sign checks use strict inequalities and the model's OWN
  declared write law (electron injection, sign(-N_E)) -- a model constant, not a
  sealed output.

T-50 ITEM 3 (design booleans): NOT ADDED.  Both designs of this lane -- Design One
(same-read contrast exponent) and Design Two (self-referenced doubling contrast) --
were REFUTED by the independent instrument refuter (LANE_T50/VERIFY_B).  No design
survives, so no design decision boolean enters this block; the gate is built
regardless, as the ruling requires.

MODES: (a) installed inside validate_geometry.py -- the trailing driver finds the
surrounding check() and runs; (b) imported by the mutation harness -- the harness
drives run_c72_contrast_checks() with its own check() and mutation; (c) standalone
(python3 c72_check_block.py) -- self-drives with a local check() and exits 0/1.
"""
import sys, os
import numpy as np

try:
    import geometry as GE            # installed inside model/ (validate_geometry.py)
except ImportError:                  # LANE_T50/GATE/ standalone or harness import
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "..", "model"))
    import geometry as GE

# ---------------------------------------------------------------- declared constants
GRID = (64, 128, 256, 512, 1024, 2048, 4096, 8192)   # block sizes, 2.1 decades
SECTOR_LEN = GRID[-1]                                # every block inside ONE sector
K_PAIRS = 64                                         # disjoint pairs per grid point
K_F1 = 32                                            # pairs for the f=1 endpoint (A2)
K_DISP = 16                                          # display-only families
DELTA_E = 5                                          # declared tolerance, e/cell (T-34)
FRAC = DELTA_E / GE.N_E                              # orientation read tolerance (sealed)
DENSITY_TOL = 0.05                                   # S6 trend tolerance (margins: harness)
DECLARED_SEED = 20260821
_LOG2N = np.log2(np.asarray(GRID, dtype=float))
_LOG10N = np.log10(np.asarray(GRID, dtype=float))
_R0_FIXED_RECORD = 32                                # M3's fixed record length (harness)


def _c72g_slope(y, x):
    x = np.asarray(x, float); y = np.asarray(y, float)
    xc = x - x.mean()
    return float((xc * (y - y.mean())).sum() / (xc * xc).sum())


def _c72g_beta(A):
    """Log-log slope of a per-grid-point statistic; None if degenerate (treated as a
    loud failure by every check that consumes it)."""
    A = np.asarray(A, float)
    if (not np.all(np.isfinite(A))) or np.any(A <= 0):
        return None
    return _c72g_slope(np.log2(A), _LOG2N)


def _c72g_mean_abs(D):
    """mean_k |D_k| per grid point: the first absolute moment E|Q| (constraint 2)."""
    return np.abs(D).mean(axis=1)


def _c72g_block_sums(values):
    """(n_sectors, SECTOR_LEN) value arrays -> (len(GRID), n_sectors) block sums
    (each block is the first N cells of its sector; blocks at one N are disjoint)."""
    c = np.cumsum(values, axis=1)
    return np.stack([c[:, N - 1] for N in GRID], axis=0)


# ---------------------------------------------------------------- the occupancy read
def _c72g_occupancy_read(rng, mutation=None):
    kind = (mutation or {}).get("kind")
    L = SECTOR_LEN
    prog_model = -float(GE.N_E)          # the model's declared write law (electron injection)

    def residuals(n, parity):
        if kind == "residual_mean_inspec":            # M2a-ii: mean +4, every |r| <= Delta
            return rng.integers(3, 6, n).astype(float)
        if kind == "residual_diff_inspec":            # M2b: +-4 by sector parity, in-spec
            return (rng.integers(3, 6, n) if parity == 0
                    else rng.integers(-5, -2, n)).astype(float)
        return rng.integers(-DELTA_E, DELTA_E + 1, n).astype(float)

    def written_sector(f, sector_idx):
        prog = rng.random(L) < f
        v = residuals(L, sector_idx % 2)
        level = prog_model
        if kind == "carrier_split" and sector_idx % 2 == 1:   # M6
            level = -prog_model
        if kind == "two_signed_write":                        # M1
            signs = rng.integers(0, 2, L) * 2 - 1
            v[prog] = (abs(prog_model) * signs[prog])
        else:
            v[prog] = level
        return v

    fams = [("Wf05", 0.5, K_PAIRS), ("Wf1", 1.0, K_F1),
            ("Wf025", 0.25, K_DISP), ("Wf075", 0.75, K_DISP)]
    read = {}
    order = []                       # physical sector order (drift acts on it)
    for name, f, n in fams:
        W = np.stack([written_sector(f, j) for j in range(n)])
        U = np.stack([residuals(L, j % 2) for j in range(n)])
        read[name] = W
        read[name + "_U"] = U
        for j in range(n):
            order.append((name, j)); order.append((name + "_U", j))
    UU = np.stack([residuals(L, j % 2) for j in range(2 * K_PAIRS)])
    read["UU"] = UU
    for j in range(2 * K_PAIRS):
        order.append(("UU", j))

    if kind == "offset_common":                               # M2a-i (constraint 1's case)
        for k in read:
            read[k] = read[k] + mutation.get("v_e", 0.5)
    if kind == "drift_inspec":                                # M7 (supplementary)
        total = len(order) * L
        pos = 0
        ramp_rate = 2.0 * DELTA_E / total                     # -Delta -> +Delta, in-spec
        for (nm, j) in order:
            read[nm][j] += (-DELTA_E + ramp_rate * (pos + np.arange(L)))
            pos += L

    BS = {k: _c72g_block_sums(v) for k, v in read.items()}

    # M3: the fixed-length record in a growing block -- per grid point N the written
    # block carries exactly _R0_FIXED_RECORD programmed cells (a DIFFERENT preparation
    # per N, which is precisely refuter A's counterexample).  Applied to the f=0.5
    # structural family only.
    if kind == "density_falling":
        R0 = _R0_FIXED_RECORD
        res = np.stack([residuals(L, j % 2) for j in range(K_PAIRS)])  # one per Wf05 pair
        c = np.cumsum(res, axis=1)
        BS["Wf05"] = np.stack([prog_model * R0 + c[:, (N - R0) - 1] for N in GRID], axis=0)

    # per-block programmed-cell counts, classified from the READ's own values: the
    # threshold is the midpoint of the read's own extreme populations (self-referenced,
    # no datasheet number).
    gmin = min(float(read[k].min()) for k in read)
    gmax = max(float(read[k].max()) for k in read)
    thr = 0.5 * (gmin + gmax)
    cls = np.cumsum(read["Wf05"] < thr, axis=1)
    counts = np.stack([cls[:, N - 1] for N in GRID], axis=0).astype(float)
    if kind == "density_falling":
        counts = np.stack([np.full(K_PAIRS, float(min(_R0_FIXED_RECORD, N)))
                           for N in GRID], axis=0)           # the record's cells all read
                                                             # below thr; erased cells not
    return BS, counts, prog_model


# ---------------------------------------------------------------- the orientation read
def _c72g_orientation_read(rng, mutation=None):
    kind = (mutation or {}).get("kind")
    L = SECTOR_LEN

    def grains(style, parity):
        if style == "data":
            g = (rng.integers(0, 2, L) * 2 - 1).astype(float)
            if kind == "one_way_data":                        # M4
                g = np.ones(L)
            return g
        if style == "dc":
            if kind == "ac_erased_dc":                        # M5
                return (rng.integers(0, 2, L) * 2 - 1).astype(float)
            return np.ones(L)
        g = (rng.integers(0, 2, L) * 2 - 1).astype(float)     # AC-erased by procedure
        if kind == "ori_diff_inspec":                         # M8 (supplementary)
            g = g + (FRAC if parity == 0 else -FRAC)
        return g

    def noisy(g):
        return g + rng.uniform(-FRAC, FRAC, g.shape)

    DATA = np.stack([noisy(grains("data", j % 2)) for j in range(K_PAIRS)])
    DATA_U = np.stack([noisy(grains("u", j % 2)) for j in range(K_PAIRS)])
    DC = np.stack([noisy(grains("dc", j % 2)) for j in range(K_PAIRS)])
    DC_U = np.stack([noisy(grains("u", j % 2)) for j in range(K_PAIRS)])
    UU = np.stack([noisy(grains("u", j % 2)) for j in range(2 * K_PAIRS)])
    read = {"DATA": DATA, "DATA_U": DATA_U, "DC": DC, "DC_U": DC_U, "UU": UU}
    if kind == "offset_common":
        for k in read:
            read[k] = read[k] + mutation.get("v_g", FRAC)
    return {k: _c72g_block_sums(v) for k, v in read.items()}


# ---------------------------------------------------------------- measurement assembly
def _c72g_measure(rng, mutation=None):
    BS, counts, prog_model = _c72g_occupancy_read(rng, mutation)
    OBS = _c72g_orientation_read(rng, mutation)
    K = K_PAIRS
    st = {}

    # occupancy contrasts
    D_wu05 = BS["Wf05"] - BS["Wf05_U"]
    D_wu1 = BS["Wf1"] - BS["Wf1_U"]
    D_wu025 = BS["Wf025"] - BS["Wf025_U"]
    D_wu075 = BS["Wf075"] - BS["Wf075_U"]
    uu = BS["UU"]
    D_cross = uu[:, 0::2] - uu[:, 1::2]                  # adjacent, opposite parity
    ev, od = uu[:, 0::2], uu[:, 1::2]                    # parity families
    D_null = np.concatenate([ev[:, 0::2] - ev[:, 1::2],  # within-parity: the null
                             od[:, 0::2] - od[:, 1::2]], axis=1)

    st["b_wu05"] = _c72g_beta(_c72g_mean_abs(D_wu05))
    st["b_wu1"] = _c72g_beta(_c72g_mean_abs(D_wu1))
    st["b_wu025"] = _c72g_beta(_c72g_mean_abs(D_wu025))
    st["b_wu075"] = _c72g_beta(_c72g_mean_abs(D_wu075))
    st["b_uu_null"] = _c72g_beta(_c72g_mean_abs(D_null))
    st["b_uu_cross"] = _c72g_beta(_c72g_mean_abs(D_cross))
    st["b_shadow"] = _c72g_beta(np.abs(uu[:, 0::2]).mean(axis=1))        # RAW: no pairing

    # signs (C-71): the write law's carrier against every structural pair, every N
    carrier = float(np.sign(prog_model))
    st["sign_all_carrier"] = bool(np.all(np.sign(D_wu05) == carrier))
    npos = (D_cross > 0).sum(axis=1); nneg = (D_cross < 0).sum(axis=1)
    st["uu_cross_maxrun"] = int(np.max(np.maximum(npos, nneg)))
    st["uu_cross_n"] = int(D_cross.shape[1])
    st["uu_cross_split_ok"] = bool(st["uu_cross_maxrun"] < D_cross.shape[1])

    # density (constraint 3): the read's own per-block programmed fraction vs N
    fhat = counts / np.asarray(GRID, float)[:, None]
    st["dens_slope"] = _c72g_slope(fhat.mean(axis=1), _LOG10N)

    # orientation contrasts
    D_data = OBS["DATA"] - OBS["DATA_U"]
    D_dc = OBS["DC"] - OBS["DC_U"]
    ouu = OBS["UU"]
    oev, ood = ouu[:, 0::2], ouu[:, 1::2]
    D_onull = np.concatenate([oev[:, 0::2] - oev[:, 1::2],
                              ood[:, 0::2] - ood[:, 1::2]], axis=1)
    st["b_data_u"] = _c72g_beta(_c72g_mean_abs(D_data))
    st["b_dc_u"] = _c72g_beta(_c72g_mean_abs(D_dc))
    st["b_uu_null_ori"] = _c72g_beta(_c72g_mean_abs(D_onull))

    # constraint-2 displays (never gated): uncentred second moment and centred variance
    st["disp_q2_wu"] = _c72g_beta((D_wu05 ** 2).mean(axis=1))
    st["disp_cvar_wu"] = _c72g_beta(np.var(D_wu05, axis=1))
    st["disp_q2_data"] = _c72g_beta((D_data ** 2).mean(axis=1))
    st["disp_cvar_data"] = _c72g_beta(np.var(D_data, axis=1))

    def sub(a, b):
        return None if (st[a] is None or st[b] is None) else st[a] - st[b]
    st["xi_occ"] = sub("b_wu05", "b_uu_null")
    st["xi_f1"] = sub("b_wu1", "b_uu_null")
    st["xi_pseudo"] = sub("b_uu_cross", "b_uu_null")
    st["xi_shadow"] = sub("b_shadow", "b_uu_null")
    st["xi_data"] = sub("b_data_u", "b_uu_null_ori")
    st["xi_dc"] = sub("b_dc_u", "b_uu_null_ori")
    return st


def _f(x):
    return "None" if x is None else f"{x:+.4f}"


def run_c72_contrast_checks(check, seed=DECLARED_SEED, mutation=None, verbose=True):
    """Run the hardened C-71/C-72 formation checks through the supplied check().
    Returns the measured-statistics dict (the mutation harness asserts on it)."""
    rng = np.random.default_rng(seed)
    st = _c72g_measure(rng, mutation)

    def gt(a, b):
        return (a is not None) and (b is not None) and (a > b)

    half = None if st["b_uu_null"] is None else st["b_uu_null"] / 2.0
    ohalf = None if st["b_uu_null_ori"] is None else st["b_uu_null_ori"] / 2.0
    gap = (None if (st["b_uu_null"] is None or st["b_uu_null_ori"] is None)
           else (st["b_uu_null"] + st["b_uu_null_ori"]) / 4.0)
    xgap = (None if (st["xi_occ"] is None or st["xi_data"] is None)
            else st["xi_occ"] - st["xi_data"])

    check("C-72 S1 ACCUMULATION (widens successor): the written-vs-unwritten contrast's "
          "exponent excess clears the midpoint of the read's own null and doubled null",
          gt(st["xi_occ"], half),
          f"xi_occ={_f(st['xi_occ'])} vs beta_null/2={_f(half)} "
          f"(b_WU={_f(st['b_wu05'])}, b_null={_f(st['b_uu_null'])})")

    check("C-72 S2 SCREENING: the orientation data contrast shows NO excess over the "
          "same read's null (below the same self-referenced midpoint)",
          (st["xi_data"] is not None and ohalf is not None and st["xi_data"] < ohalf),
          f"xi_data={_f(st['xi_data'])} vs beta_null/2={_f(ohalf)}")

    check("C-72 S3 POSITIVE CONTROL (control C successor): the same read's DC-saturated "
          "contrast ACCUMULATES -- the instrument demonstrably sees accumulation",
          gt(st["xi_dc"], ohalf),
          f"xi_dc={_f(st['xi_dc'])} vs beta_null/2={_f(ohalf)}")

    check("C-72 S4 DISCRIMINATOR: the occupancy excess exceeds the orientation data "
          "excess by more than half the mean measured null exponent",
          gt(xgap, gap),
          f"xi_occ-xi_data={_f(xgap)} vs (b_null+b_null_ori)/4={_f(gap)}")

    check("C-71 S5 SIGN LAW: every structural pair at every N carries the model's own "
          "write-law carrier sign, and the unwritten cross pairs are strictly "
          "sign-split (no within-part consistency absent writing)",
          st["sign_all_carrier"] and st["uu_cross_split_ok"],
          f"all_carrier={st['sign_all_carrier']} "
          f"uu_cross_maxrun={st['uu_cross_maxrun']}/{st['uu_cross_n']}")

    check("C-72 S6 DENSITY GUARD (constraint 3): the read's own per-block programmed "
          "fraction shows no trend in N",
          abs(st["dens_slope"]) <= DENSITY_TOL,
          f"|slope|={abs(st['dens_slope']):.4f} vs tol {DENSITY_TOL}")

    check("C-72 B2 UNWRITTEN NULL (control B repaired): the unwritten pool's "
          "cross-family pseudo-contrast shows NO excess over the within-family null",
          (st["xi_pseudo"] is not None and half is not None and st["xi_pseudo"] < half),
          f"xi_pseudo={_f(st['xi_pseudo'])} vs beta_null/2={_f(half)}")

    check("C-72 A2 PATTERN ENDPOINT (control A dropped per INST-15, slot repaired): "
          "the f=1 all-programmed read through the SAME pipeline accumulates",
          gt(st["xi_f1"], half),
          f"xi_f1={_f(st['xi_f1'])} vs beta_null/2={_f(half)}")

    if verbose:
        print("  ----- non-decision displays (constraints 1 and 2; never gated) -----")
        print(f"  NOTE shadow RAW unwritten exponent b_shadow={_f(st['b_shadow'])} "
              f"xi_shadow={_f(st['xi_shadow'])} -- constraint 1: the raw exponent is "
              "not baseline-free; asserted only in the mutation harness")
        print(f"  NOTE constraint-2 table (occupancy W-U): mean|D| {_f(st['b_wu05'])} "
              f"| uncentred D^2 {_f(st['disp_q2_wu'])} (= 2*beta) | centred var "
              f"{_f(st['disp_cvar_wu'])} (non-discriminating)")
        print(f"  NOTE constraint-2 table (orientation DATA-U): mean|D| "
              f"{_f(st['b_data_u'])} | uncentred D^2 {_f(st['disp_q2_data'])} | "
              f"centred var {_f(st['disp_cvar_data'])} (~1 for BOTH encodings)")
        print(f"  NOTE pattern sweep (displayed, not gated -- no in-suite mutation "
              f"gives exponent pattern-dependence a failing branch, and a check that "
              f"cannot fail is the defect this block removes): "
              f"b_WU(f=0.25)={_f(st['b_wu025'])} b_WU(f=0.75)={_f(st['b_wu075'])}")
    return st


# ---------------------------------------------------------------- mode dispatch
def _selfdrive():
    tally = {"pass": 0, "fail": 0}

    def _check(name, cond, detail=""):
        tally["pass" if cond else "fail"] += 1
        print(f"  {'PASS' if cond else 'FAIL'}  {name}  {detail}")

    print("C-71/C-72 HARDENED CONTRAST BLOCK -- standalone run "
          f"(declared seed {DECLARED_SEED})")
    print("=" * 78)
    run_c72_contrast_checks(_check)
    print("=" * 78)
    print(f"  BLOCK: {tally['pass']} PASS, {tally['fail']} FAIL")
    return tally["fail"]


try:
    check                      # exists when installed inside validate_geometry.py
    _C72G_INSTALLED = True
except NameError:
    _C72G_INSTALLED = False

if _C72G_INSTALLED:
    run_c72_contrast_checks(check)
elif __name__ == "__main__":
    sys.exit(1 if _selfdrive() else 0)
# else: imported by the mutation harness, which drives run_c72_contrast_checks itself
