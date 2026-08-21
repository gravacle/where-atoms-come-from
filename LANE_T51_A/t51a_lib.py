"""LANE_T51_A shared machinery for THE MEASUREMENT (the Second Lump first computation,
commissioned in FIELD_INSTRUMENT_V001.md section 1; C-92).

REUSED MACHINERY, NAMED: everything geometric and sector-exact is imported from this
lane's own sealed gate script g1_connected_wenc.py (which itself imports the sealed
model/geometry.py by name: Torus venue conventions, exact F_2 kit, coset machinery).
The generalizations here (arbitrary winding-sector labels, source-absent and
probe-absent configurations, explicit plaquette-constant accounting so absolute
sector energies and the port certificates are real computations) are NEW for the
measurement and follow the gate's construction line by line.

DISCIPLINE CARRIED (binding; stated once here, referenced by every script):
  D-1  no classical gravitational form is required of, or tested against, anything in
       this lane; the shape of F vs separation is an OUTPUT.
  D-8  verdicts are computed booleans with both branches reachable; no literal
       expected value sits on any physics decision path.  (The a1 calibration section
       is a REPRODUCTION comparison against the design exploration's quoted numbers,
       which are prior data parsed at runtime and labeled as such -- a port check in
       the gate's prior-art pattern, not a physics verdict.)
  D-15 every reported zero carries a positive control in a DIFFERENT configuration in
       the same table; construction certificates are labeled certificates and are
       never counted as controls.
  D-24 separation claims are stated in earned quantities only: connected w_enc per
       placement (from the sealed gate table, seal-verified at parse time) with d_gen
       descriptive.  Lattice coordinates are construction labels.
  Principal's directive (2026-08-20, binding, quoted in C-92): the mechanism for
       accumulation is whatever it proves to be; every outcome registers as the
       surface's own law.

BORROWED IDEAS, OWNERS NAMED: hole-pair (defect) records -- Bravyi-Kitaev
quant-ph/9811052 (named in geometry.py); toric carrier -- Kitaev quant-ph/9707021;
linked-cluster cancellation behind connected w_enc -- Goldstone linked-cluster /
Kato-Bloch degenerate perturbation theory (named binding by the Second Lump
critiques); induced two-lump energy as the mechanism-class -- degenerate perturbation
theory / Schrieffer-Wolff, RKKY and Casimir-Polder as comparison-only owners (named
in the design, comparison never enters any verdict).
"""
import ast
import hashlib
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from g1_connected_wenc import (  # noqa: E402  (the gate's sealed machinery, by name)
    Venue, Torus, plaq_edge_mask, star_edge_mask, geometry_block,
    is_connected_string, FLOOR_EIG, FLOOR_USE, TOL_K, WITNESS_MIN,
    LAMS as ONSET_GRID, bool_word,
)
from geometry import sp_pair  # noqa: E402  (sealed symplectic pairing, o54c_lib)

# ---- declared tolerances for THE MEASUREMENT (stated before any result) -------------
COMM_LAMS = (0.02, 0.05, 0.10)  # the commissioned mediator grid (FIELD_INSTRUMENT s.1)
TOL_CAL_REL = 5e-3   # calibration reproduction: quoted values carry 3 significant
#                      figures, so half-ulp quantization is <= 2.5e-3 relative; 2x taken
TOL_PORT_D = 1e-8    # port certificate, Delta vs the gate's sealed prints (10 digits)
TOL_PORT_F = 1e-3    # port certificate, F vs the gate's sealed prints (4 digits)
TOL_SWAP = 0.25      # C2: onset agreement inside a Gamma-equivalent swap pair
#                      (= TOL_K, the same integer-resolution scale the gate declared)
TOL_BACK = 1e-3      # C5: back-action bound, energy units of the unit stabilizer term
#                      (basis: three orders below the sector gap ~2 and one order
#                      above the largest commissioned reading at the far placement;
#                      scale fixed from a scratchpad magnitude survey BEFORE this
#                      sealed run -- logged in D24_AUDIT.txt; both branches reachable:
#                      a stronger probe or coupling fails it)
BOUND_LOG2 = 11      # the declared sector-dimension bound ~2^11 (FIELD_INSTRUMENT s.1)

GATE_TXT = os.path.join(HERE, "g1_connected_wenc.txt")


def sha256_file(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def verify_gate_seal():
    """Recompute the gate output's sha256 and compare against its sealed sidecar.
       Returns the hash; raises if the seal does not verify."""
    h = sha256_file(GATE_TXT)
    recorded = open(GATE_TXT + ".sha256").read().split()[0]
    if h != recorded:
        raise RuntimeError("gate seal FAILED to verify: %s != %s" % (h, recorded))
    return h


_BLOCK_RE = re.compile(
    r"PLACEMENT (3x3-\d\d[HV])\s+\[construction labels: probe stars "
    r"\((\d), (\d)\),\((\d), (\d)\); source plaquettes \((\d), (\d)\),\((\d), (\d)\)"
    r".*?earned separation:\s+w_enc_conn = (\d)\s+\(d_gen = (\d)"
    r".*?OLD-SPEC minimum:\s+w_enc_old  = (\d)"
    r".*?CONNECTED minimum:\s+w_enc_conn = (\d)\s+\(N_min = (\d+); "
    r"winding classes at min: (\[[^\]]*\])"
    r".*?onset k_hat = ([\d.]+)", re.S)

_DYN_RE = re.compile(
    r"lam=([\d.]+)\s+Delta\(-1\)=([+-][\d.]+e[+-]\d+)\s+"
    r"Delta\(\+1\)=([+-][\d.]+e[+-]\d+)\s+F=([+-][\d.]+e[+-]\d+)")


def parse_gate_table():
    """Parse the seal-verified gate output: per 3x3 placement, the earned quantities
       (connected w_enc; d_gen descriptive), the gate's onset, the winding classes of
       the connected minima, and the gate's dynamics rows (for the port certificate).
       The gate is INPUT here, never re-run."""
    seal = verify_gate_seal()
    txt = open(GATE_TXT).read()
    out = {}
    for m in _BLOCK_RE.finditer(txt):
        tag = m.group(1)
        block_start = m.start()
        block_end = txt.find("PLACEMENT", m.end())
        block = txt[block_start:block_end if block_end > 0 else len(txt)]
        dyn = {float(d.group(1)): (float(d.group(2)), float(d.group(3)),
                                   float(d.group(4)))
               for d in _DYN_RE.finditer(block)}
        out[tag] = dict(
            probe=((int(m.group(2)), int(m.group(3))),
                   (int(m.group(4)), int(m.group(5)))),
            src=((int(m.group(6)), int(m.group(7))),
                 (int(m.group(8)), int(m.group(9)))),
            w_conn=int(m.group(13)), d_gen=int(m.group(11)),
            w_old=int(m.group(12)), n_min=int(m.group(14)),
            conn_windings=ast.literal_eval(m.group(15)),  # list of int 2-tuples
            k_hat_gate=float(m.group(16)), dyn=dyn)
    if len(out) != 18:
        raise RuntimeError("gate table parse found %d placements, not 18" % len(out))
    return out, seal


# =====================================================================================
# sector-exact machinery, generalized from the gate (plaquette constants explicit)
# =====================================================================================
def sector_states_g(V, src, b, w1, w2):
    """Z-basis states of the sector: intact plaquette parities 0 (+1 eigenvalue),
       source plaquette parities (0 if b=+1 else 1), winding parities (w1, w2)
       [parity p <-> eigenvalue (-1)^p of the sealed zbar1/zbar2]."""
    Lx, Ly = V.T.Lx, V.T.Ly
    want_src = 0 if b == +1 else 1
    keep = np.ones(len(V.states), dtype=bool)
    for y in range(Ly):
        for x in range(Lx):
            tgt = want_src if (x, y) in src else 0
            keep &= (V.par[("P", x, y)] == tgt)
    keep &= (V.par[("W1",)] == w1)
    keep &= (V.par[("W2",)] == w2)
    return V.states[keep], V.popc[keep]


def sector_spectrum_g(V, probe, src, b, lam, w1=0, w2=0, probe_present=True,
                      src_present=True):
    """Sector-exact spectrum of H = -sum(stars in H) - sum(plaquettes in H)
       + lam * sum_e Z_e restricted to the sector.

       Configurations: src_present=True -> the two source plaquettes are HOLES
       (terms absent; parities are the exact labels b).  src_present=False -> the
       source plaquettes are RESTORED as terms (requires b=+1); this is the
       source-absent configuration.  probe_present likewise for the two probe stars.
       Plaquette terms in H are exact constants in-sector and are INCLUDED, so
       absolute energies are comparable across configurations (the certificates are
       then real computed identities, not comparisons of one call with itself).
       Returns dict(E0, E1, delta, gap, m, witness[, ...]); witness only when the
       probe is present."""
    T = V.T
    n = T.n
    if not src_present:
        assert b == +1, "source-absent configuration has no written value"
        assert not src, "src must be () when src_present=False"
    st, popc = sector_states_g(V, src, b, w1, w2)
    m = len(st)
    assert m == (1 << (n - (T.Lx * T.Ly - 1) - 2)), \
        "sector dimension %d differs from the computed constraint count" % m
    assert m <= (1 << BOUND_LOG2), \
        "sector dimension %d exceeds the declared 2^%d bound" % (m, BOUND_LOG2)
    idx = {int(s): i for i, s in enumerate(st)}
    H = np.zeros((m, m))
    zsum = (n - 2 * popc).astype(float)
    H[np.arange(m), np.arange(m)] = lam * zsum
    # plaquette terms present in H: all plaquettes if src absent, intact ones if the
    # source holes exist.  In-sector each contributes -(eigenvalue) = -(+1) for
    # intact/all-restored (parity 0); source holes contribute no term.
    n_plaq_terms = T.Lx * T.Ly if not src_present else T.Lx * T.Ly - len(src)
    H[np.arange(m), np.arange(m)] += -1.0 * n_plaq_terms
    holes = probe if probe_present else ()
    stars = [T.star(x, y) for y in range(T.Ly) for x in range(T.Lx)
             if (x, y) not in holes]
    for A in stars:
        for i, s in enumerate(st):
            H[i, idx[int(s) ^ A]] -= 1.0
    evals, evecs = np.linalg.eigh(H)
    out = dict(E0=float(evals[0]), E1=float(evals[1]),
               delta=float(evals[1] - evals[0]),
               gap=float(evals[2] - evals[1]), m=m)
    if probe_present:
        Ap = T.star(*probe[0])
        perm = np.array([idx[int(s) ^ Ap] for s in st])
        out["witness"] = abs(float(evecs[:, 0] @ evecs[perm, 1]))
    return out


def F_reading(V, probe, src, lam, w1=0, w2=0):
    """The field reading in one winding sector: F = Delta(b=-1) - Delta(b=+1),
       with the doublet witness and band gap minima alongside."""
    rm = sector_spectrum_g(V, probe, src, -1, lam, w1, w2)
    rp = sector_spectrum_g(V, probe, src, +1, lam, w1, w2)
    return dict(lam=lam, w1=w1, w2=w2, dM=rm["delta"], dP=rp["delta"],
                F=rm["delta"] - rp["delta"],
                witness=min(rm["witness"], rp["witness"]),
                gap=min(rm["gap"], rp["gap"]), m=rm["m"],
                E0M=rm["E0"], E0P=rp["E0"])


def onset_measure(V, probe, src, grid=ONSET_GRID):
    """The declared onset estimator (identical discipline to the gate): adjacent-pair
       log-slopes on the geometric grid, pairs usable iff both |F| >= FLOOR_USE;
       k_hat = slope of the smallest usable pair; drift to the next usable pair
       reported beside it; FLOOR_EIG beside every F."""
    rows = [F_reading(V, probe, src, lam) for lam in grid]
    slopes = []
    for a, b in zip(rows, rows[1:]):
        if abs(a["F"]) >= FLOOR_USE and abs(b["F"]) >= FLOOR_USE:
            k = np.log(abs(b["F"]) / abs(a["F"])) / np.log(b["lam"] / a["lam"])
            slopes.append((a["lam"], b["lam"], float(k)))
    k_hat = slopes[0][2] if slopes else None
    drift = abs(slopes[0][2] - slopes[1][2]) if len(slopes) > 1 else None
    return rows, slopes, k_hat, drift


# =====================================================================================
# sealing helpers
# =====================================================================================
class Emitter:
    def __init__(self):
        self.lines = []

    def __call__(self, s=""):
        self.lines.append(s)
        print(s)

    def seal(self, path):
        with open(path, "w") as f:
            f.write("\n".join(self.lines) + "\n")
        h = sha256_file(path)
        with open(path + ".sha256", "w") as f:
            f.write("%s  %s\n" % (h, os.path.basename(path)))
        print("\nsealed -> %s\nsha256 %s" % (path, h))
        return h


def write_json(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=float)
    h = sha256_file(path)
    with open(path + ".sha256", "w") as f:
        f.write("%s  %s\n" % (h, os.path.basename(path)))
    return h


def load_json_verified(path):
    h = sha256_file(path)
    recorded = open(path + ".sha256").read().split()[0]
    if h != recorded:
        raise RuntimeError("sidecar seal FAILED: %s" % path)
    return json.load(open(path))
