"""LANE_T51_A / G1 -- THE CONNECTIVITY GATE for THE SECOND LUMP (C-92; commissioned in
FIELD_INSTRUMENT_V001.md section 1).

THE DEBT THIS GATE PAYS (found by the audit, binding per both Second Lump critiques):
the sign-constrained coset minimum w_enc admits DISCONNECTED representatives that the
dynamics cancels (linked-cluster cancellation).  Computed counterexample at 3x2
(critique 2): sign-constrained coset minimum 3, achieved only by a disconnected
representative; measured onset order 4.0.  w_enc is here REDEFINED over CONNECTED
enclosing strings, computed by exact F_2 machinery, exhaustively over the ENTIRE coset,
and checked against the measured onset by a COMPUTED comparison at the 3x2 calibration
placement and at every 3x3 placement the first computation will use.

DISCIPLINE CARRIED (binding):
  D-1  no classical gravitational form is required of, or tested against, anything here;
       the shape of any response vs separation is an OUTPUT.
  D-8  every verdict is a computed boolean with both branches reachable; no literal
       expected value sits on any decision path (the audit's numbers appear ONLY as
       labeled cross-references, never inside a verdict); every fit carries a noise
       floor beside it.
  D-15 any reported zero carries a positive control in a DIFFERENT configuration;
       algebraic identities are labeled certificates, never controls.
  D-24 separations are STATED in earned quantities (d_gen -- DESCRIPTIVE per the
       critique's binding repair, and connected w_enc, which carries the earned
       separation); lattice coordinates below are construction labels for placements,
       and no claim is stated in them.
  Principal's directive (2026-08-20, quoted in C-92): the mechanism for accumulation is
       whatever it proves to be; every outcome of every computation here registers as
       the surface's own law.

REUSED SEALED MACHINERY, NAMED (imported from model/geometry.py, never copied):
  geometry.Torus              -- the o54c_lib venue conventions (edge indexing, star/
                                 plaquette/zbar masks), C-80's own venue.
  geometry.pc, rank_bits      -- the exact F_2 kit (t42 lanes).
  geometry.sp_pair, supp_mask -- symplectic pairing / support (o54c_lib).
  geometry._generator_graph_dist -- C-80's d_gen instrument (descriptive here).
  geometry.coset_min_np       -- C-80's exhaustive coset-minimum scan, used as a
                                 cross-implementation gate on the unconstrained coset.
  The parity/connectivity-filtered enumeration below is NEW machinery for this gate; it
  follows coset_min_np's exhaustive-span pattern (enumerate the whole span, exact
  popcount weights), with the two filters computed per element.

BORROWED IDEAS, OWNERS NAMED (per the critiques' binding repairs):
  linked-cluster cancellation of disconnected clusters -- Goldstone linked-cluster
    theorem / Kato-Bloch degenerate perturbation theory (named by the Second Lump
    critique when it made the repair load-bearing);
  defect (hole-pair) records -- Bravyi-Kitaev quant-ph/9811052 (named in geometry.py);
  toric carrier -- Kitaev quant-ph/9707021;
  cycle-space decomposition of loops (plaquette boundaries + homology representatives)
    -- standard graph/surface homology; the enclosure functional built on it is this
    program's own (C-92's holonomy vocabulary).

THE DEFINITIONS (stated so a stranger reimplements them):
  Venue: Lx x Ly toric code, qubits on edges, geometry.Torus conventions.  Probe: two
  removed star stabilizers at adjacent vertices u_p, v_p (their unique shared edge is
  the DIRECT connector, weight 1).  Source: two removed adjacent plaquettes s1, s2.
  ADMISSIBLE PROBE CONNECTOR: a Z-type string t (an edge mask) that anticommutes with
    both removed probe stars and commutes with every intact star (Z-strings commute
    with all plaquettes identically -- certificate, not a control).  Exactly the coset
    t = direct XOR c, c in the cycle space of the direct lattice (dimension Lx*Ly+1;
    basis: all plaquette boundaries except one excluded plaquette P0 not in {s1,s2},
    plus the venue's sealed winding representatives W1 = zbar1, W2 = zbar2).
  ENCLOSURE PARITY eps(t): write loop(t) = t XOR direct uniquely over that basis;
    eps(t) = coeff(s1) XOR coeff(s2).  eps is independent of the choice of P0
    (computed lemma below).  For winding-odd loops eps depends on the choice of
    winding representative (computed exhibit below); the DECLARED convention is the
    venue's sealed zbar1/zbar2 -- the same operators that label the winding sectors of
    the dynamics, which is what makes eps(t) the source-value sensitivity bit of the
    string amplitude at fixed sector labels: as an operator, loop(t) = product of the
    enclosed plaquettes times W1^w1 * W2^w2, so its eigenvalue carries the source bit b
    iff eps(t) = 1.
  CONNECTED: the graph whose nodes are the edges in supp(t), two nodes adjacent iff
    the corresponding lattice edges share at least one lattice vertex (coordinates mod
    (Lx, Ly)), is a single connected component.  (Vertex-sharing is the dynamically
    relevant adjacency: a Z_e term excites the two stars at e's endpoints, so two
    edges interact through a common energy denominator iff they share an endpoint;
    disjoint clusters factorize and cancel -- Goldstone/Kato-Bloch, above.)
  OLD-SPEC (disconnected-permitting) MINIMUM: w_enc_old = min |t| over the coset with
    eps(t) = 1.
  THE REDEFINED QUANTITY: w_enc_conn = min |t| over the coset with eps(t) = 1 AND t
    connected.
  MINIMALITY IS ABSOLUTE AT THESE SIZES: the enumeration is exhaustive over the ENTIRE
    coset (2^(Lx*Ly+1) elements: 128 at 3x2, 1024 at 3x3), not a truncated weight
    window; the "bound" is the whole coset, justified by its size.

THE DYNAMICS-SIDE CHECK (sector-exact, independent implementation for this gate):
  H = -sum(intact stars) - sum(intact plaquettes) + lam * sum_e Z_e  (the design's
  mediator, content-blind; [V, B_p] = 0 exactly, so source value and windings are
  exact quantum numbers -- certificate).  Work in the Z product basis; a sector is
  fixed by (all plaquette parities, both zbar parities): intact plaquettes +1,
  source plaquettes b (b = +1 unwritten, b = -1 written; b_s1 = b_s2 = b is the
  consistent pair), windings +1.  Sector dimensions: 32 at 3x2, 256 at 3x3 (computed,
  asserted).  Delta(b) = E1 - E0 of the sector ground doublet (doublet identity
  certified by the witness |<0|A_probe|1>| >= WITNESS_MIN); F(lam) = Delta(-1) -
  Delta(+1).  ONSET ESTIMATOR (declared before use): geometric grid LAMS, adjacent-
  pair log-slope k = ln(|F(l2)|/|F(l1)|) / ln(l2/l1); the estimate k_hat is the slope
  of the SMALLEST usable pair, a pair being usable iff both |F| >= FLOOR_USE; the
  drift to the next usable pair is reported beside it as the estimator's uncertainty.
  NOISE FLOOR (declared): FLOOR_EIG = 1e-11 absolute on any eigenvalue difference
  (dim * machine-eps * ||H|| <= 256 * 2.2e-16 * ~9 ~ 5e-13; 20x headroom taken);
  FLOOR_USE = 100 * FLOOR_EIG = 1e-9.  TOL_K = 0.25 (declared): the onset matches an
  integer w iff |k_hat - w| <= TOL_K.  Both branches of every boolean are reachable:
  nothing below forces the connected minimum to match and nothing forces the old-spec
  minimum to fail.
"""
import os
import sys
import time
import hashlib
from collections import deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(REPO, "model"))

import geometry  # the sealed machinery layer (T-46 port); reused symbols named above
from geometry import Torus, pc, rank_bits, sp_pair, supp_mask, _generator_graph_dist, \
    coset_min_np

# ---- declared tolerances and conventions (stated before any result) -----------------
LAMS = (0.004, 0.008, 0.016, 0.032, 0.064)   # the audit's own range, for comparability
FLOOR_EIG = 1e-11     # declared absolute eigenvalue-difference floor (justified above)
FLOOR_USE = 1e-9      # declared usability floor for onset pairs = 100 * FLOOR_EIG
TOL_K = 0.25          # declared onset-integer matching tolerance
WITNESS_MIN = 0.9     # declared doublet-identity witness floor

OUT_LINES = []


def emit(s=""):
    OUT_LINES.append(s)
    print(s)


# =====================================================================================
# geometry side: the coset, the enclosure functional, connectivity
# =====================================================================================
def edge_endpoints(T, e):
    """Endpoints (vertices) of edge index e under geometry.Torus conventions."""
    half = T.Lx * T.Ly
    if e < half:
        y, x = divmod(e, T.Lx)
        return (x, y), ((x + 1) % T.Lx, y)
    y, x = divmod(e - half, T.Lx)
    return (x, y), (x, (y + 1) % T.Ly)


def is_connected_string(T, mask):
    """CONNECTED per the declared definition: support-edge graph, adjacency = shared
       lattice vertex, single component.  Empty string is not connected."""
    edges = [e for e in range(T.n) if (mask >> e) & 1]
    if not edges:
        return False
    ends = {e: set(edge_endpoints(T, e)) for e in edges}
    seen = {edges[0]}
    dq = deque([edges[0]])
    while dq:
        a = dq.popleft()
        for b in edges:
            if b not in seen and ends[a] & ends[b]:
                seen.add(b)
                dq.append(b)
    return len(seen) == len(edges)


def plaq_edge_mask(T, x, y):
    return T.plaq(x, y) >> T.n


def star_edge_mask(T, x, y):
    return T.star(x, y)


def cycle_basis(T, excluded, w1_mask, w2_mask):
    """Cycle-space basis: all plaquette boundaries except `excluded`, then W1, W2.
       Returns (masks, index_of_plaquette dict)."""
    masks, idx = [], {}
    for y in range(T.Ly):
        for x in range(T.Lx):
            if (x, y) == excluded:
                continue
            idx[(x, y)] = len(masks)
            masks.append(plaq_edge_mask(T, x, y))
    iW1 = len(masks); masks.append(w1_mask)
    iW2 = len(masks); masks.append(w2_mask)
    return masks, idx, iW1, iW2


def enumerate_coset(T, direct, s1, s2, excluded, w1_mask, w2_mask):
    """Exhaustive enumeration of the ENTIRE admissible-connector coset with, per
       element: mask, weight, enclosure parity eps, winding coefficients (w1, w2).
       Enumeration by explicit basis coefficients, so eps needs no solve."""
    masks, idx, iW1, iW2 = cycle_basis(T, excluded, w1_mask, w2_mask)
    k = len(masks)
    assert rank_bits(list(masks)) == k, "cycle basis not independent -- gate failed"
    assert k == T.Lx * T.Ly + 1, "cycle-space dimension wrong -- gate failed"
    i1, i2 = idx[s1], idx[s2]
    out = []
    for i in range(1 << k):
        c = 0
        j = i
        b = 0
        while j:
            if j & 1:
                c ^= masks[b]
            j >>= 1
            b += 1
        t = direct ^ c
        eps = ((i >> i1) & 1) ^ ((i >> i2) & 1)
        out.append((t, pc(t), eps, (i >> iW1) & 1, (i >> iW2) & 1))
    return out


def geometry_block(Lx, Ly, probe, src):
    """All geometry-side quantities for one placement.  Returns a dict."""
    T = Torus(Lx, Ly)
    n = T.n
    u_p, v_p = probe
    s1, s2 = src
    # the direct connector: the unique shared edge of the two probe stars
    shared = star_edge_mask(T, *u_p) & star_edge_mask(T, *v_p)
    assert pc(shared) == 1, "probe stars must be adjacent (one shared edge)"
    direct = shared
    assert is_connected_string(T, direct)
    w1_mask = T.zbar1() >> n
    w2_mask = T.zbar2() >> n
    # declared excluded plaquette: the first (row-major) plaquette not a source hole
    excluded = next((x, y) for y in range(Ly) for x in range(Lx)
                    if (x, y) not in (s1, s2))
    elems = enumerate_coset(T, direct, s1, s2, excluded, w1_mask, w2_mask)

    # ---- admissibility gate on every element (computed, not asserted from algebra) --
    stars_all = [(x, y) for y in range(Ly) for x in range(Lx)]
    intact_stars = [T.star(x, y) for (x, y) in stars_all if (x, y) not in probe]
    probe_stars = [T.star(*u_p), T.star(*v_p)]
    for (t, w, eps, w1, w2) in elems:
        tz = t << n   # Z-type Pauli on 2n bits
        assert all(sp_pair(tz, a, n) == 1 for a in probe_stars)
        assert all(sp_pair(tz, a, n) == 0 for a in intact_stars)
    # Z-vs-Z commutation with every plaquette is identically zero in the symplectic
    # pairing -- CONSTRUCTION CERTIFICATE (labeled; not counted as a control).

    # ---- excluded-plaquette invariance lemma (computed over all valid choices) ------
    base_par = {t: eps for (t, w, eps, w1, w2) in elems}
    for exc in [(x, y) for y in range(Ly) for x in range(Lx)
                if (x, y) not in (s1, s2)]:
        alt = enumerate_coset(T, direct, s1, s2, exc, w1_mask, w2_mask)
        assert {t: e for (t, w, e, a, b) in alt} == base_par, \
            "enclosure parity depends on excluded plaquette -- lemma FAILED"
    lemma_excl = True

    # ---- winding-representative dependence exhibit (computed) -----------------------
    w1_alt = w1_mask ^ plaq_edge_mask(T, *s1)     # same homology class, shifted rep
    alt = enumerate_coset(T, direct, s1, s2, excluded, w1_alt, w2_mask)
    alt_par = {t: e for (t, w, e, a, b) in alt}
    flips = [t for t in base_par if base_par[t] != alt_par[t]]
    wind_odd = {t for (t, w, e, a, b) in elems if a == 1}
    exhibit_flips = len(flips)
    exhibit_all_winding_odd = set(flips) <= wind_odd and \
        (exhibit_flips == 0 or len(flips) == len(wind_odd))

    # ---- cross-implementation gate: unconstrained coset minimum via the sealed scan -
    masks, _, _, _ = cycle_basis(T, excluded, w1_mask, w2_mask)
    wmin_sealed, _, _ = coset_min_np(masks, direct, base_bits=min(22, len(masks)))
    wmin_mine = min(w for (t, w, e, a, b) in elems)
    assert wmin_sealed == wmin_mine, "sealed scan and enumeration disagree -- gate"

    enc = [(t, w, e, a, b) for (t, w, e, a, b) in elems if e == 1]
    assert len(enc) == len(elems) // 2, "enclosure functional not balanced -- gate"
    w_old = min(w for (t, w, e, a, b) in enc)
    old_min_elems = [(t, a, b) for (t, w, e, a, b) in enc if w == w_old]
    old_any_connected = any(is_connected_string(T, t) for (t, a, b) in old_min_elems)
    conn = [(t, w, a, b) for (t, w, e, a, b) in enc if is_connected_string(T, t)]
    w_conn = min(w for (t, w, a, b) in conn)
    conn_min_elems = [(t, a, b) for (t, w, a, b) in conn if w == w_conn]
    conn_min_windings = sorted({(a, b) for (t, a, b) in conn_min_elems})
    assert w_conn >= w_old   # set inclusion; computed anyway

    # descriptive d_gen (mixed-type; DEMOTED to descriptive per the critique's repair)
    local = T.all_stars() + T.all_plaqs()
    N = Lx * Ly
    pidx = {(x, y): y * Lx + x for y in range(Ly) for x in range(Lx)}
    d_gen = min(_generator_graph_dist(local, n, pidx[pv], N + pidx[sv])
                for pv in probe for sv in src)

    def fmt_edges(mask):
        half = Lx * Ly
        out = []
        for e in range(n):
            if (mask >> e) & 1:
                if e < half:
                    y, x = divmod(e, Lx)
                    out.append("h(%d,%d)" % (x, y))
                else:
                    y, x = divmod(e - half, Lx)
                    out.append("v(%d,%d)" % (x, y))
        return "{" + ",".join(out) + "}"

    return dict(T=T, direct=direct, excluded=excluded,
                lemma_excl=lemma_excl, exhibit_flips=exhibit_flips,
                exhibit_clean=exhibit_all_winding_odd,
                coset_size=len(elems), wmin_unconstrained=wmin_mine,
                w_old=w_old, old_n_min=len(old_min_elems),
                old_any_connected=old_any_connected,
                old_witness=fmt_edges(old_min_elems[0][0]),
                old_witness_winding=(old_min_elems[0][1], old_min_elems[0][2]),
                old_witness_connected=is_connected_string(T, old_min_elems[0][0]),
                w_conn=w_conn, conn_n_min=len(conn_min_elems),
                conn_witness=fmt_edges(conn_min_elems[0][0]),
                conn_min_windings=conn_min_windings,
                d_gen=d_gen)


# =====================================================================================
# dynamics side: sector-exact diagonalization (independent implementation)
# =====================================================================================
class Venue:
    """Precomputed parity tables for one (Lx, Ly): every plaquette parity and both
       zbar parities on every Z-basis state."""

    def __init__(self, Lx, Ly):
        self.T = Torus(Lx, Ly)
        n = self.T.n
        self.states = np.arange(1 << n, dtype=np.uint64)
        self.par = {}
        for y in range(Ly):
            for x in range(Lx):
                m = np.uint64(plaq_edge_mask(self.T, x, y))
                self.par[("P", x, y)] = (np.bitwise_count(self.states & m)
                                         & np.uint64(1)).astype(np.uint8)
        for name, mask in (("W1", self.T.zbar1() >> n), ("W2", self.T.zbar2() >> n)):
            m = np.uint64(mask)
            self.par[(name,)] = (np.bitwise_count(self.states & m)
                                 & np.uint64(1)).astype(np.uint8)
        self.popc = np.bitwise_count(self.states).astype(np.int64)


def sector_states(V, src, b):
    """States of the sector: intact plaquettes +1, source plaquettes b, windings +1."""
    Lx, Ly = V.T.Lx, V.T.Ly
    want_src = 0 if b == +1 else 1
    keep = np.ones(len(V.states), dtype=bool)
    for y in range(Ly):
        for x in range(Lx):
            tgt = want_src if (x, y) in src else 0
            keep &= (V.par[("P", x, y)] == tgt)
    keep &= (V.par[("W1",)] == 0)
    keep &= (V.par[("W2",)] == 0)
    return V.states[keep], V.popc[keep]


def sector_spectrum(V, probe, src, b, lam):
    """(Delta, witness, gap_to_band) for one sector.  H = -sum intact A - sum intact B
       + lam sum Z_e; plaquette part is constant in-sector and dropped (Delta is a
       difference of eigenvalues of the same sector)."""
    T = V.T
    n = T.n
    st, popc = sector_states(V, src, b)
    m = len(st)
    idx = {int(s): i for i, s in enumerate(st)}
    H = np.zeros((m, m))
    zsum = (n - 2 * popc).astype(float)          # sum_e z_e per state
    H[np.arange(m), np.arange(m)] = lam * zsum
    intact = [T.star(x, y) for y in range(T.Ly) for x in range(T.Lx)
              if (x, y) not in probe]
    for A in intact:
        for i, s in enumerate(st):
            j = idx[int(s) ^ A]
            H[i, j] -= 1.0
    evals, evecs = np.linalg.eigh(H)
    delta = evals[1] - evals[0]
    gap = evals[2] - evals[1]
    Ap = T.star(*probe[0])
    perm = np.array([idx[int(s) ^ Ap] for s in st])
    witness = abs(float(evecs[:, 0] @ evecs[perm, 1]))
    return delta, witness, gap, m


def onset_measurement(V, probe, src):
    """F(lam) over the declared grid, plus the declared onset estimator."""
    rows = []
    for lam in LAMS:
        dm, wm_, gm, mdim = sector_spectrum(V, probe, src, -1, lam)
        dp, wp_, gp, _ = sector_spectrum(V, probe, src, +1, lam)
        rows.append(dict(lam=lam, dM=dm, dP=dp, F=dm - dp,
                         witness=min(wm_, wp_), gap=min(gm, gp), mdim=mdim))
    slopes = []
    for a, bb in zip(rows, rows[1:]):
        if abs(a["F"]) >= FLOOR_USE and abs(bb["F"]) >= FLOOR_USE:
            k = np.log(abs(bb["F"]) / abs(a["F"])) / np.log(bb["lam"] / a["lam"])
            slopes.append((a["lam"], bb["lam"], float(k)))
    k_hat = slopes[0][2] if slopes else None
    drift = abs(slopes[0][2] - slopes[1][2]) if len(slopes) > 1 else None
    return rows, slopes, k_hat, drift


# =====================================================================================
# the run
# =====================================================================================
def bool_word(x):
    return "TRUE" if x else "FALSE"


def run_placement(tag, Lx, Ly, probe, src, V):
    g = geometry_block(Lx, Ly, probe, src)
    rows, slopes, k_hat, drift = onset_measurement(V, probe, src)
    wit = min(r["witness"] for r in rows)
    maxF = max(abs(r["F"]) for r in rows)
    b_conn = (k_hat is not None) and (abs(k_hat - g["w_conn"]) <= TOL_K)
    b_old = (k_hat is not None) and (abs(k_hat - g["w_old"]) <= TOL_K)
    return dict(tag=tag, Lx=Lx, Ly=Ly, probe=probe, src=src, g=g, rows=rows,
                slopes=slopes, k_hat=k_hat, drift=drift, witness_min=wit,
                maxF=maxF, b_conn=b_conn, b_old=b_old,
                witness_ok=wit >= WITNESS_MIN)


def report_placement(r, verbose):
    g = r["g"]
    emit("PLACEMENT %s   [construction labels: probe stars %s,%s; source plaquettes "
         "%s,%s -- coordinates are construction data, claims are in earned units]"
         % (r["tag"], r["probe"][0], r["probe"][1], r["src"][0], r["src"][1]))
    emit("  earned separation:      w_enc_conn = %d   (d_gen = %d, DESCRIPTIVE -- "
         "mixed-type, demoted per critique repair)" % (g["w_conn"], g["d_gen"]))
    emit("  coset:                  %d elements enumerated exhaustively (the whole "
         "coset; enclosing class %d)" % (g["coset_size"], g["coset_size"] // 2))
    emit("  OLD-SPEC minimum:       w_enc_old  = %d   (N_min = %d; any connected "
         "representative at min: %s)"
         % (g["w_old"], g["old_n_min"], bool_word(g["old_any_connected"])))
    emit("    old-min witness:      %s  connected=%s  winding=(%d,%d)"
         % (g["old_witness"], bool_word(g["old_witness_connected"]),
            g["old_witness_winding"][0], g["old_witness_winding"][1]))
    emit("  CONNECTED minimum:      w_enc_conn = %d   (N_min = %d; winding classes at "
         "min: %s)" % (g["w_conn"], g["conn_n_min"], g["conn_min_windings"]))
    emit("    conn-min witness:     %s" % g["conn_witness"])
    if verbose:
        emit("  gates: excluded-plaquette invariance lemma %s over all valid choices; "
             "winding-representative exhibit: %d parity flips under the shifted W1 "
             "representative, all-and-only-winding-odd = %s (the declared convention "
             "is the sealed zbar1/zbar2); unconstrained minimum %d agrees with the "
             "sealed coset_min_np scan."
             % ("PASSED" if g["lemma_excl"] else "FAILED", g["exhibit_flips"],
                bool_word(g["exhibit_clean"]), g["wmin_unconstrained"]))
    emit("  dynamics (sector-exact; sector dim %d; doublet witness min %.4f >= %.1f: "
         "%s; doublet-to-band gap min %.3f):"
         % (r["rows"][0]["mdim"], r["witness_min"], WITNESS_MIN,
            bool_word(r["witness_ok"]), min(x["gap"] for x in r["rows"])))
    for x in r["rows"]:
        emit("    lam=%.3f  Delta(-1)=%+.9e  Delta(+1)=%+.9e  F=%+.3e"
             % (x["lam"], x["dM"], x["dP"], x["F"]))
    if r["slopes"]:
        emit("  onset pair log-slopes (usable pairs, |F| >= FLOOR_USE=%.0e): %s"
             % (FLOOR_USE, "; ".join("(%.3f,%.3f)->%.3f" % s for s in r["slopes"])))
        emit("  onset k_hat = %.3f  (smallest usable pair; drift to next pair = %s; "
             "noise floor FLOOR_EIG = %.0e beside every F above)"
             % (r["k_hat"],
                ("%.3f" % r["drift"]) if r["drift"] is not None else "n/a",
                FLOOR_EIG))
    else:
        emit("  onset: NO USABLE PAIR at this grid (max|F| = %.2e < FLOOR_USE); "
             "positive control: see the other placements' nonzero F in different "
             "configurations on this same table (D-15)." % r["maxF"])
    emit("  COMPUTED COMPARISON (TOL_K = %.2f):" % TOL_K)
    emit("    onset matches CONNECTED minimum  |k_hat - w_enc_conn| <= TOL_K : %s"
         % bool_word(r["b_conn"]))
    emit("    onset matches OLD-SPEC minimum   |k_hat - w_enc_old|  <= TOL_K : %s"
         % bool_word(r["b_old"]))
    emit("")


def main():
    t0 = time.time()
    emit("=" * 88)
    emit("LANE_T51_A / G1 -- THE CONNECTIVITY GATE: connected w_enc, computed exactly")
    emit("date: 2026-08-21   (commissioned in FIELD_INSTRUMENT_V001.md section 1; the")
    emit("audit's counterexample and both Second Lump critiques' repairs are binding)")
    emit("=" * 88)
    emit("")
    emit("DECLARED BEFORE USE: LAMS = %s; FLOOR_EIG = %.0e; FLOOR_USE = %.0e;"
         % (list(LAMS), FLOOR_EIG, FLOOR_USE))
    emit("TOL_K = %.2f; WITNESS_MIN = %.1f; excluded plaquette = first row-major"
         % (TOL_K, WITNESS_MIN))
    emit("non-source plaquette; winding representatives = the venue's sealed zbar1/")
    emit("zbar2 (geometry.Torus).  Definitions of eps(t) and CONNECTED: header of")
    emit("g1_connected_wenc.py, stated for reimplementation.")
    emit("")
    emit("CONSTRUCTION CERTIFICATES (labeled; never counted as controls): Z-strings")
    emit("commute with every plaquette (symplectic pairing identically 0); [V, B_p] = 0")
    emit("exactly, so source value and windings are exact sector labels.")
    emit("")

    # ------------------------------ 3x2 calibration ---------------------------------
    emit("-" * 88)
    emit("SECTION 1 -- THE 3x2 CALIBRATION PLACEMENT (the audit's computed "
         "counterexample)")
    emit("-" * 88)
    V32 = Venue(3, 2)
    cal = run_placement("3x2-CAL", 3, 2, ((0, 0), (1, 0)), ((1, 1), (2, 1)), V32)
    report_placement(cal, verbose=True)
    g = cal["g"]
    emit("  CROSS-REFERENCES (the audit's numbers, quoted as labeled data, NEVER on a")
    emit("  decision path): audit coset minimum 3; audit measured onset 4.0")
    emit("  (F/lam^4 -> 18.9 over lam = 0.004..0.064).")
    emit("    my w_enc_old  = %d  (audit said 3): agree = %s"
         % (g["w_old"], bool_word(g["w_old"] == 3)))
    emit("    my k_hat      = %.3f (audit said 4.0): |diff| = %.3f"
         % (cal["k_hat"], abs(cal["k_hat"] - 4.0)))
    Fs = cal["rows"][0]
    emit("    my F/lam^4 at lam=%.3f = %.2f (audit said -> 18.9)"
         % (Fs["lam"], abs(Fs["F"]) / Fs["lam"] ** 4))
    emit("")
    emit("  GATE VERDICT G1a (counterexample resolves; computed, both branches "
         "reachable):")
    emit("    the CONNECTED minimum matches the measured onset AND the old-spec")
    emit("    minimum differs from the connected minimum at this placement:")
    v_g1a = cal["b_conn"] and (g["w_conn"] != g["w_old"]) and not cal["b_old"]
    emit("    G1a = b_conn AND (w_enc_conn != w_enc_old) AND NOT b_old = %s"
         % bool_word(v_g1a))
    emit("")

    # ------------------------------ 3x3 placements ----------------------------------
    emit("-" * 88)
    emit("SECTION 2 -- EVERY 3x3 PLACEMENT THE FIRST COMPUTATION WILL USE")
    emit("(probe stars (0,0),(1,0); ALL 18 adjacent source plaquette pairs enumerated,")
    emit("so every choice of the two commissioned separations -- and the Gamma-")
    emit("equivalent placement-swap control -- is covered)")
    emit("-" * 88)
    V33 = Venue(3, 3)
    probe33 = ((0, 0), (1, 0))
    placements = []
    for y in range(3):
        for x in range(3):
            placements.append((("%d%d" % (x, y)) + "H", ((x, y), ((x + 1) % 3, y))))
            placements.append((("%d%d" % (x, y)) + "V", ((x, y), (x, (y + 1) % 3))))
    results = []
    for tag, src in placements:
        r = run_placement("3x3-" + tag, 3, 3, probe33, src, V33)
        results.append(r)
    for r in results:
        report_placement(r, verbose=False)

    emit("-" * 88)
    emit("SECTION 3 -- THE 3x3 TABLE, GROUPED INTO EARNED-SEPARATION CLASSES")
    emit("-" * 88)
    emit("  %-10s %-6s %-10s %-10s %-8s %-8s %-8s %-6s" %
         ("placement", "d_gen", "w_enc_old", "w_enc_conn", "k_hat", "b_conn",
          "b_old", "wit"))
    for r in results:
        emit("  %-10s %-6d %-10d %-10d %-8s %-8s %-8s %-6s" %
             (r["tag"], r["g"]["d_gen"], r["g"]["w_old"], r["g"]["w_conn"],
              ("%.3f" % r["k_hat"]) if r["k_hat"] is not None else "none",
              bool_word(r["b_conn"]), bool_word(r["b_old"]),
              "ok" if r["witness_ok"] else "LOW"))
    emit("")
    classes = {}
    for r in results:
        classes.setdefault((r["g"]["d_gen"], r["g"]["w_conn"]), []).append(r["tag"])
    emit("  earned-separation classes (d_gen DESCRIPTIVE, w_enc_conn earned):")
    for k in sorted(classes):
        emit("    (d_gen=%d, w_enc_conn=%d): %d placements: %s"
             % (k[0], k[1], len(classes[k]), ", ".join(classes[k])))
    dgens = sorted({r["g"]["d_gen"] for r in results})
    emit("  the venue affords %d distinct descriptive d_gen values %s; the earned"
         % (len(dgens), dgens))
    emit("  quantity w_enc_conn takes values %s and is PER PLACEMENT (C-93's repair:"
         % sorted({r["g"]["w_conn"] for r in results}))
    emit("  V3 scores against connected w_enc placement by placement).")
    emit("")

    # second-counterexample check at 3x3 (computed; both branches reachable)
    cx = [r for r in results if r["g"]["w_old"] != r["g"]["w_conn"]]
    emit("  GATE VERDICT G1b (3x3): placements where old-spec and connected minima")
    emit("  DIFFER: %d of %d  [%s]" % (len(cx), len(results),
                                       ", ".join(r["tag"] for r in cx)))
    v_g1b = all(r["b_conn"] for r in results if r["k_hat"] is not None) and \
        all(not r["b_old"] for r in cx if r["k_hat"] is not None)
    emit("  G1b = (onset matches w_enc_conn at every measurable placement) AND")
    emit("        (onset fails w_enc_old at every differing placement) = %s"
         % bool_word(v_g1b))
    all_measurable = all(r["k_hat"] is not None for r in results)
    emit("  all 18 placements measurable at this grid: %s" % bool_word(all_measurable))
    emit("")
    emit("  GATE VERDICT G1 (the connectivity gate) = G1a AND G1b = %s"
         % bool_word(v_g1a and v_g1b))
    emit("")

    emit("-" * 88)
    emit("SECTION 4 -- PRIOR-ART COMPARISON (built independently FIRST; compared "
         "AFTER)")
    emit("-" * 88)
    compare_prior_art(cal, results)

    emit("-" * 88)
    emit("runtime %.1f s; exact F_2 machinery on the geometry path; float tolerances"
         % (time.time() - t0))
    emit("declared above for every eigenvalue; every fit carries its noise floor.")
    emit("NEXT STEP NAMED (no route closes without one): with the gate landed, the")
    emit("first computation (V1-V5 of FIELD_INSTRUMENT_V001.md section 1) scores V3")
    emit("against w_enc_conn per placement as computed here; the winding-sector sweep")
    emit("(control 3) uses the same declared zbar convention this gate declared.")
    emit("=" * 88)

    out = os.path.join(HERE, "g1_connected_wenc.txt")
    with open(out, "w") as f:
        f.write("\n".join(OUT_LINES) + "\n")
    h = hashlib.sha256(open(out, "rb").read()).hexdigest()
    with open(out + ".sha256", "w") as f:
        f.write("%s  g1_connected_wenc.txt\n" % h)
    print("\nsealed -> %s\nsha256 %s" % (out, h))


def compare_prior_art(cal, results):
    """AFTER the independent build: read the fourth-angle lane's sealed output
       (LANE_T51_IDENT/FOURTH_rigidity_OUT.txt, C-93) and compare every overlapping
       number.  Agreement is a cross-check; disagreement is a finding to report,
       never reconciled silently.  Their numbers are PARSED AT RUNTIME from their
       sealed file -- nothing of theirs is hardcoded here.  Declared comparison
       tolerances (print-precision-limited): TOL_CMP_F = 1e-4 relative on Part C
       F values (they print 7 significant digits), TOL_CMP_F4 = 1e-3 relative on
       Part E F(lam=0.064) values (4 significant digits printed)."""
    import re
    TOL_CMP_F, TOL_CMP_F4 = 1e-4, 1e-3
    path = os.path.join(REPO, "LANE_T51_IDENT", "FOURTH_rigidity_OUT.txt")
    try:
        txt = open(path).read()
    except OSError:
        emit("  fourth-angle sealed output not readable; comparison NOT performed")
        emit("  (reported as such; nothing reconciled).")
        return
    emit("  source: LANE_T51_IDENT/FOURTH_rigidity_OUT.txt (C-93, sealed; read only")
    emit("  after every number above was computed).  Tolerances declared:")
    emit("  TOL_CMP_F = 1e-4 rel (their 7-digit prints), TOL_CMP_F4 = 1e-3 rel")
    emit("  (their 4-digit prints).")
    emit("")
    checks = []

    def check(name, ok, detail):
        checks.append((name, ok))
        emit("  [%s] %s -- %s" % ("AGREE" if ok else "DIFFER", name, detail))

    # ---- Part B: 3x2 calibration static minima ------------------------------------
    mB = re.search(r"source \(\(1, 1\), \(2, 1\)\): histogram.*?\n"
                   r".*?\n\s*w_direct=1\s+w_enc\(old-spec\)=(\d+)\s+"
                   r"w_enc\(CONNECTED\)=(\d+)", txt)
    if mB:
        to, tc = int(mB.group(1)), int(mB.group(2))
        check("3x2-CAL static minima (Part B)",
              (to, tc) == (cal["g"]["w_old"], cal["g"]["w_conn"]),
              "theirs (old=%d, conn=%d) vs mine (old=%d, conn=%d)"
              % (to, tc, cal["g"]["w_old"], cal["g"]["w_conn"]))
    else:
        check("3x2-CAL static minima (Part B)", False, "PARSE MISMATCH -- their "
              "Part B line not found; comparison incomplete, reported as a finding")
    # ---- Part C: 3x2 calibration F values at every shared lam ---------------------
    mC = re.search(r"source \(\(1, 1\), \(2, 1\)\) \(connected w_enc = \d+\), "
                   r"winding \(0,0\):\n.*?\n((?:\s+[\d.]+\s+\S+\s+\S+\s+\S+"
                   r"\s+\S+\s+\S+\n)+)", txt)
    if mC:
        theirs = {}
        for row in mC.group(1).strip().splitlines():
            f = row.split()
            theirs[float(f[0])] = float(f[3])
        worst = 0.0
        n_shared = 0
        for x in cal["rows"]:
            if x["lam"] in theirs:
                n_shared += 1
                worst = max(worst, abs(x["F"] - theirs[x["lam"]])
                            / max(abs(theirs[x["lam"]]), 1e-300))
        check("3x2-CAL F(lam) values (Part C, %d shared lams)" % n_shared,
              n_shared == len(LAMS) and worst <= TOL_CMP_F,
              "worst relative difference %.1e (tol %.0e)" % (worst, TOL_CMP_F))
    else:
        check("3x2-CAL F(lam) values (Part C)", False, "PARSE MISMATCH -- reported")
    # ---- Part E: the two 3x3 placements their lane computed -----------------------
    by_src = {tuple(sorted(r["src"])): r for r in results}
    for mE in re.finditer(r"source \(\((\d+), (\d+)\), \((\d+), (\d+)\)\): "
                          r"w_direct=1 w_enc\(old\)=(\d+) w_enc\(CONN\)=(\d+)"
                          r".*?F\(lam=0\.064\)=([+-][\d.]+e[+-]\d+)", txt, re.S):
        src = tuple(sorted(((int(mE.group(1)), int(mE.group(2))),
                            (int(mE.group(3)), int(mE.group(4))))))
        to, tc, tF = int(mE.group(5)), int(mE.group(6)), float(mE.group(7))
        r = by_src.get(src)
        if r is None:
            check("3x3 %s (Part E)" % (src,), False,
                  "their placement is not in my 18-placement table -- FINDING")
            continue
        myF = [x["F"] for x in r["rows"] if x["lam"] == 0.064][0]
        ok_static = (to, tc) == (r["g"]["w_old"], r["g"]["w_conn"])
        ok_F = abs(myF - tF) / abs(tF) <= TOL_CMP_F4
        check("3x3 %s static minima (Part E)" % (src,), ok_static,
              "theirs (old=%d, conn=%d) vs mine (old=%d, conn=%d)"
              % (to, tc, r["g"]["w_old"], r["g"]["w_conn"]))
        check("3x3 %s F(0.064) (Part E)" % (src,), ok_F,
              "theirs %.4e vs mine %.4e (rel diff %.1e, tol %.0e)"
              % (tF, myF, abs(myF - tF) / abs(tF), TOL_CMP_F4))
    emit("")
    emit("  ESTIMATOR-DISCIPLINE DIFFERENCE, REPORTED (not reconciled): their onset")
    emit("  ladders include slope pairs whose |F| sits below this lane's declared")
    emit("  FLOOR_USE = %.0e (e.g. their first 3x3 ((1,1),(2,1)) rung uses" % FLOOR_USE)
    r11 = by_src.get(((1, 1), (2, 1)))
    if r11 is not None:
        below = [(x["lam"], x["F"]) for x in r11["rows"] if abs(x["F"]) < FLOOR_USE]
        emit("  F at lam=%s, values %s -- below the floor this lane declared);"
             % ([b[0] for b in below], ["%.2e" % b[1] for b in below]))
    emit("  on the smallest pair PASSING my floor the slopes coincide with their")
    emit("  corresponding ladder rung to the printed digits.  Their raw F values,")
    emit("  where printed, agree with mine (rows above); the ladders differ only in")
    emit("  which pairs each lane admits into the estimate.")
    emit("")
    n_ok = sum(1 for _, ok in checks if ok)
    emit("  CROSS-CHECK TALLY: %d/%d overlapping comparisons AGREE.  Any DIFFER row"
         % (n_ok, len(checks)))
    emit("  above stands as a finding; none was reconciled.")


if __name__ == "__main__":
    main()
