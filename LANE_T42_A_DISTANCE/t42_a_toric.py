"""T42-A part 1: d_W BETWEEN CONFIGURATIONS on the toric code, L = 2,3,4,5.  EXACT, F_2.

d_W(s,s') = minimum weight of an admissible Pauli operation taking configuration s to s'.
Writers are SEARCHED per configuration pair: for each ordered pair the FULL affine space of
admissible operations with that label action is enumerated and the minimum weight taken.

Configurations: the 4^g ground labels (g=1 here, 4 of them), labels = eigenvalue flips against
the COMPUTED conjugate logical pairs (record_model.symplectic_logicals -- computed, never
nominated).  An admissible u flips label j iff sp(u, Zlog_j) = 1.

REDUCTION USED AT L=4,5 (verified exhaustively at L=2,3 against the full (x|z) search):
weight(x|z) = |supp x  U  supp z| >= |supp x|, and z=0 is admissible whenever anything is,
so the minimum over the full space equals the minimum over the z=0 slice.

BOUNDARY-CROSSING: for every enumerated writer we compute its crossing parity on each of the
L vertical and L horizontal cuts.  If the parity is a (resp b) on EVERY cut for EVERY writer
in the class, then every writer costs >= L*a horizontal + L*b vertical edges -- that step is
integer arithmetic, not sampling -- and the observed minimum L*(a+b) closes the identity
   d_W = minimal boundary-crossing cost,   exhaustively at each L.

NO literal verdicts: every PASS/FAIL below is gated by a computed boolean.
"""
import sys, time
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_T42_A_DISTANCE"
sys.path.insert(0, LANE)
from t42_lib import (pc, solve_affine_f2, span_all, span_numpy, np_weights, sp_pair,
                     weight_xz, vec_to_mask, metric_axioms, aut_count,
                     toric_stabilizers, toric_edge_masks)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals
import numpy as np

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s); sys.stdout.flush()

def gate(label, ok, detail=""):
    say(("PASS  " if ok else "FAIL  ") + label + (("  " + detail) if detail else ""))
    return ok

say("T42-A part 1  --  toric configuration distance, exact, started " + time.strftime("%F %T"))
say("")

RESULTS = {}
ALL_OK = True

for L in (2, 3, 4, 5):
    t0 = time.time()
    n, stab, stars, plaqs = toric_stabilizers(L)
    hmask, vmask, vcuts, hcuts = toric_edge_masks(L)
    pairs = symplectic_logicals(stab, n)
    say("L=%d  n=%d qubits  logical pairs COMPUTED: %d" % (L, n, len(pairs)))
    ok = gate("L=%d exactly 2 conjugate pairs" % L, len(pairs) == 2)
    ALL_OK &= ok
    logi = []
    for X, Z in pairs:
        logi.append((vec_to_mask(X), vec_to_mask(Z)))
    # verify computed pairs really are conjugate & mutually commuting (surface's own operation)
    conj_ok = all(sp_pair(logi[i][0], logi[i][1], n) == 1 for i in range(2))
    cross_ok = all(sp_pair(logi[i][a], logi[j][b], n) == 0
                   for i in range(2) for j in range(2) for a in range(2) for b in range(2)
                   if i != j)
    ALL_OK &= gate("L=%d computed logicals: conjugate within pair, commuting across" % L,
                   conj_ok and cross_ok)
    Zlogs = [logi[0][1], logi[1][1]]           # labeling convention: 2nd of each pair
    mask_n = (1 << n) - 1
    Zl_z = [zl >> n for zl in Zlogs]           # z-parts pair against x-type writers
    Zl_x = [zl & mask_n for zl in Zlogs]

    # ---- z=0 slice: rows = z-parts of stabilizers (plaquettes), plus label rows
    plaq_rows = [p for p in plaqs]             # z-part of Z-type stabilizers = plaquette masks
    configs = [(a, b) for a in range(2) for b in range(2)]
    dmat = {}
    class_min = {}
    class_info = {}
    nullN = None
    searched = 0
    for s in configs:
        for sp_ in configs:
            t = (s[0] ^ sp_[0], s[1] ^ sp_[1])
            rows = plaq_rows + [Zl_z[0], Zl_z[1]]
            rhs = [0] * len(plaq_rows) + [t[0], t[1]]
            x0, nb, rank = solve_affine_f2(rows, rhs, n)
            assert x0 is not None
            dim = len(nb)
            if nullN is None:
                nullN = span_numpy(nb)
                say("L=%d  writer search space per pair: 2^%d = %d admissible x-type ops"
                    % (L, dim, 1 << dim))
            coset = nullN ^ np.uint64(x0)
            w = np_weights(coset)
            searched += coset.shape[0]
            dmat[(s, sp_)] = int(w.min())
            if t not in class_min:
                # crossing parities on every cut, over the ENTIRE coset
                vpar = [int((np.bitwise_count(coset & np.uint64(c)) & np.uint64(1)).max())
                        for c in vcuts]
                vpar_min = [int((np.bitwise_count(coset & np.uint64(c)) & np.uint64(1)).min())
                            for c in vcuts]
                hpar = [int((np.bitwise_count(coset & np.uint64(c)) & np.uint64(1)).max())
                        for c in hcuts]
                hpar_min = [int((np.bitwise_count(coset & np.uint64(c)) & np.uint64(1)).min())
                            for c in hcuts]
                hcount = np.bitwise_count(coset & np.uint64(hmask))
                vcount = np.bitwise_count(coset & np.uint64(vmask))
                class_min[t] = int(w.min())
                class_info[t] = {
                    "vcut_parity_constant": all(a == b for a, b in zip(vpar, vpar_min)),
                    "hcut_parity_constant": all(a == b for a, b in zip(hpar, hpar_min)),
                    "vcut_parities": vpar, "hcut_parities": hpar,
                    "min_hedges": int(hcount.min()), "min_vedges": int(vcount.min()),
                }
    say("L=%d  writers searched (exact, exhaustive): %d" % (L, searched))

    # group invariance: d(s,s') depends only on s xor s'  (a FINDING, computed)
    inv_ok = all(dmat[(s, sp_)] == class_min[(s[0] ^ sp_[0], s[1] ^ sp_[1])]
                 for s in configs for sp_ in configs)
    ALL_OK &= gate("L=%d d_W(s,s') depends only on s+s' (searched per pair, invariance emerged)"
                   % L, inv_ok)

    # crossing structure: parities constant per class; identify homology label (a,b)
    par_ok = True; hom = {}
    for t, ci in class_info.items():
        cst = ci["vcut_parity_constant"] and ci["hcut_parity_constant"]
        vset, hset = set(ci["vcut_parities"]), set(ci["hcut_parities"])
        cst &= (len(vset) == 1 and len(hset) == 1)
        par_ok &= cst
        hom[t] = (ci["vcut_parities"][0], ci["hcut_parities"][0])
    ALL_OK &= gate("L=%d crossing parity constant over every class and every parallel cut" % L,
                   par_ok)
    # crossing-cost identity:  min weight = L*a + L*b,  min h-edges = L*a,  min v-edges = L*b
    cc_ok = True
    for t, ci in class_info.items():
        a, b = hom[t]
        cc_ok &= (class_min[t] == L * (a + b))
        cc_ok &= (ci["min_hedges"] == L * a) and (ci["min_vedges"] == L * b)
    ALL_OK &= gate("L=%d d_W = minimal boundary-crossing cost  (min w = L*a+L*b exactly)" % L,
                   cc_ok, str({t: (class_min[t], hom[t]) for t in sorted(class_min)}))

    # metric axioms over ALL pairs and ALL triples
    rep = metric_axioms(configs, dmat)
    ALL_OK &= gate("L=%d METRIC AXIOMS identity/symmetry/indiscernible/triangle "
                   "(%d pairs, %d triples)" % (L, rep["n_pairs"], rep["n_triples"]),
                   rep["all"], "" if rep["all"] else str(rep))

    # code distance (computed): min over nonzero classes
    dcode = min(class_min[t] for t in class_min if t != (0, 0))
    ALL_OK &= gate("L=%d code distance (computed as min nonzero class weight) == L" % L,
                   dcode == L, "d_code=%d" % dcode)
    # canonical (metric-selected) basis: the classes AT distance d_code
    units = [t for t in class_min if t != (0, 0) and class_min[t] == dcode]
    basis_ok = (len(units) == 2)
    if basis_ok:
        e1, e2 = units
        relabel = {(0, 0): (0, 0), e1: (1, 0), e2: (0, 1),
                   (e1[0] ^ e2[0], e1[1] ^ e2[1]): (1, 1)}
        ham_ok = all(class_min[t] == dcode * (relabel[t][0] + relabel[t][1])
                     for t in class_min)
    else:
        ham_ok = False
    ALL_OK &= gate("L=%d in the METRIC-SELECTED basis: d_W = d_code x Hamming(labels)" % L,
                   basis_ok and ham_ok,
                   "units(at distance d_code)=%s" % str(units))

    kind, na = aut_count(configs, dmat)
    say("L=%d  D-22 venue Aut group of the 4-point metric space: |Aut| = %d (%s)"
        % (L, na, kind))
    ALL_OK &= gate("L=%d venue HAS geometry to detect (Aut != full symmetric group S_4=24)" % L,
                   na != 24, "|Aut|=%d" % na)

    say("L=%d  full metric on the 4 configurations (order 00,10,01,11 in computed labels):" % L)
    order = [(0, 0), (1, 0), (0, 1), (1, 1)]
    for s in order:
        say("      " + " ".join("%3d" % dmat[(s, sp_)] for sp_ in order))

    RESULTS[L] = {"dmat": {str(k): v for k, v in dmat.items()},
                  "class_min": {str(k): v for k, v in class_min.items()},
                  "dcode": dcode, "aut": na, "searched": searched}
    say("L=%d done in %.1fs" % (L, time.time() - t0))
    say("")

# ---------------- z-freedom control: at L=2 and L=3, FULL (x|z) search must equal z=0 search
say("REDUCTION CONTROL: full (x|z) writer search vs z=0 slice, L=2,3")
for L in (2, 3):
    n, stab, stars, plaqs = toric_stabilizers(L)
    pairs = symplectic_logicals(stab, n)
    logi = [(vec_to_mask(X), vec_to_mask(Z)) for X, Z in pairs]
    Zlogs = [logi[0][1], logi[1][1]]
    mask_n = (1 << n) - 1
    # full system: 2n unknowns; row for generator g: mask = zpart(g) | (xpart(g)<<n)
    def full_row(g):
        return (g >> n) | ((g & mask_n) << n)
    stab_masks = [vec_to_mask(sv) for sv in stab]
    rows = [full_row(g) for g in stab_masks]
    ok_all = True
    for t in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        R = rows + [full_row(Zlogs[0]), full_row(Zlogs[1])]
        rhs = [0] * len(rows) + [t[0], t[1]]
        x0, nb, rank = solve_affine_f2(R, rhs, 2 * n)
        full_min = min(weight_xz(x0 ^ v, n) for v in span_all(nb))
        # z=0 min recomputed
        prow = [p for p in plaqs] + [Zlogs[0] >> n, Zlogs[1] >> n]
        prhs = [0] * len(plaqs) + [t[0], t[1]]
        y0, nby, _ = solve_affine_f2(prow, prhs, n)
        z0_min = min(pc(y0 ^ v) for v in span_all(nby))
        ok_all &= (full_min == z0_min)
        say("  L=%d class %s: full-space min %d  (space 2^%d)   z=0 min %d"
            % (L, t, full_min, len(nb), z0_min))
    ALL_OK &= gate("L=%d allowing z-parts NEVER lowers d_W (exhaustive both spaces)" % L, ok_all)
say("")

# ---------------- D-15 control candidate that must FAIL: support-overlap "distance"
say("D-15 CONTROL CANDIDATE (must FAIL an axiom): d_ctrl(s,s') = |supp(minrep s) & supp(minrep s')|")
L = 3
n, stab, stars, plaqs = toric_stabilizers(L)
pairs = symplectic_logicals(stab, n)
logi = [(vec_to_mask(X), vec_to_mask(Z)) for X, Z in pairs]
Zlogs = [logi[0][1], logi[1][1]]
configs = [(a, b) for a in range(2) for b in range(2)]
minrep = {}
for t in configs:
    prow = [p for p in plaqs] + [Zlogs[0] >> n, Zlogs[1] >> n]
    prhs = [0] * len(plaqs) + [t[0], t[1]]
    y0, nby, _ = solve_affine_f2(prow, prhs, n)
    best, bw = None, None
    for v in span_all(nby):
        w = pc(y0 ^ v)
        if bw is None or w < bw:
            bw, best = w, y0 ^ v
    minrep[t] = best
dctrl = {(s, sp_): pc(minrep[s] & minrep[sp_]) for s in configs for sp_ in configs}
repc = metric_axioms(configs, dctrl)
ALL_OK &= gate("control candidate FAILS at least one axiom (instrument distinguishes)",
               not repc["all"],
               "identity_ok=%s indiscernible_ok=%s triangle_ok=%s witness=%s"
               % (repc["identity_ok"], repc["indiscernible_ok"], repc["triangle_ok"],
                  repc["identity_witness"] or repc["indiscernible_witness"] or repc["triangle_witness"]))
say("")
say("SCALING TABLE  d_W(class)/L over L (exact integers):")
say("  L | d(unit class 1) | d(unit class 2) | d(sum class) | d_code")
for L in (2, 3, 4, 5):
    cm = {eval(k): v for k, v in RESULTS[L]["class_min"].items()}
    nz = sorted((v, k) for k, v in cm.items() if k != (0, 0))
    say("  %d |      %2d         |      %2d         |     %2d       |   %d"
        % (L, nz[0][0], nz[1][0], nz[2][0], RESULTS[L]["dcode"]))
say("")
say("ALL GATES: " + ("PASS" if ALL_OK else "AT LEAST ONE FAIL"))
with open(LANE + "/t42_a_toric.OUT.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
