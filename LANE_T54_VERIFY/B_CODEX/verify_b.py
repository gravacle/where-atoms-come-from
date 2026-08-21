#!/usr/bin/env python3
"""Independent, default-refuted T-54 verification for classes and writing.

Sealed T-44/T-48 artifacts are inputs only. This verifier does not execute a sealed
lane program. It parses sealed outputs, calls repaired model APIs afresh, and uses
independent walk/counting machinery for off-sealed controls. Expected refusal probes
catch only ValueError and accept it only when its message identifies the intended
domain gate; all positive controls are called without an exception wrapper.
"""

import ast
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import re
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "model"
sys.path.insert(0, str(MODEL))

import classes as CL  # noqa: E402
import writing as WR  # noqa: E402
from project_model import URM  # noqa: E402


REQUIRED = {
    "shared": {"S-SEALS", "S-EXCEPTION-HYGIENE"},
    "classes": {
        "CL-SUITE", "CL-SEALED", "CL-DP", "CL-MUC", "CL-KERNEL",
        "CL-EXACT", "CL-OFF", "CL-REFUSE-PRICE", "CL-REFUSE-VENUE",
        "CL-D8", "CL-D15",
    },
    "writing": {
        "WR-SUITE", "WR-SEALED-A", "WR-SEALED-B", "WR-SEALED-C",
        "WR-OFF", "WR-D15", "WR-REFUSE-E1", "WR-REFUSE-E2",
        "WR-REFUSE-E3", "WR-SIGNED", "WR-SURFACE", "WR-EXACT", "WR-D8",
    },
}
seen = {key: set() for key in REQUIRED}
failures = {key: 0 for key in REQUIRED}
passes = {key: 0 for key in REQUIRED}


def check(group, ident, cond, detail=""):
    """Record one required predicate; every family starts REFUTED until all pass."""
    if ident in seen[group]:
        cond = False
        detail = (detail + "; " if detail else "") + "duplicate predicate id"
    seen[group].add(ident)
    if bool(cond):
        passes[group] += 1
        print(f"OK      {ident}  {detail}")
    else:
        failures[group] += 1
        print(f"REFUTE  {ident}  {detail}")


def expected_value_error(call, message_fragment):
    """Return success only for the intended ValueError; other exceptions propagate."""
    try:
        call()
    except ValueError as exc:
        message = str(exc)
        return message_fragment in message, f"ValueError={message}"
    return False, "call returned instead of refusing"


def all_fractions(rows):
    return all(isinstance(x, F) for row in rows for x in row.values())


def independent_dp3(K):
    """Fresh dictionary DP on Z^3; no classes.py counting function is used."""
    steps = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
             (0, 0, 1), (0, 0, -1))
    grid = {(0, 0, 0): 1}
    out = [dict(grid)]
    for _ in range(K):
        nxt = defaultdict(int)
        for (x, y, z), n in grid.items():
            for dx, dy, dz in steps:
                nxt[(x + dx, y + dy, z + dz)] += n
        grid = dict(nxt)
        out.append(grid)
    return out


def independent_walk_counts(adj, src, K):
    """Fresh exact transfer action on (neighbor,multiplicity) adjacency rows."""
    v = [0] * len(adj)
    v[src] = 1
    out = [list(v)]
    for _ in range(K):
        nxt = [0] * len(adj)
        for i, n in enumerate(v):
            for j, mult in adj[i]:
                nxt[j] += n * mult
        v = nxt
        out.append(list(v))
    return out


def manifest_is_intact(lane):
    manifest = lane / "SEALS.sha256"
    if not manifest.is_file():
        return False, 0
    count = 0
    for raw in manifest.read_text().splitlines():
        if not raw.strip():
            continue
        digest, rel = raw.split(None, 1)
        target = lane / rel.strip().lstrip("*")
        if not target.is_file() or sha256(target.read_bytes()).hexdigest() != digest:
            return False, count
        count += 1
    return count > 0, count


def capture_groups(pattern, text):
    match = re.search(pattern, text, re.MULTILINE)
    return match.groups() if match else None


def suite_result(filename, label):
    run = subprocess.run(
        [sys.executable, str(MODEL / filename)], cwd=MODEL,
        text=True, capture_output=True, check=False,
    )
    parsed = capture_groups(rf"{label}:\s+(\d+) PASS,\s+(\d+) FAIL", run.stdout)
    return run.returncode, parsed


def literal_check_conditions(path):
    """D-8 static tell: check(..., True/False) and assert True/False."""
    tree = ast.parse(path.read_text(), filename=str(path))
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = node.func.id if isinstance(node.func, ast.Name) else None
            if name == "check" and len(node.args) >= 2:
                arg = node.args[1]
                if isinstance(arg, ast.Constant) and isinstance(arg.value, bool):
                    bad.append(node.lineno)
        if (isinstance(node, ast.Assert) and isinstance(node.test, ast.Constant)
                and isinstance(node.test.value, bool)):
            bad.append(node.lineno)
    return bad


def float_literals_in_functions(path, function_names):
    tree = ast.parse(path.read_text(), filename=str(path))
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in function_names:
            for child in ast.walk(node):
                if isinstance(child, ast.Constant) and isinstance(child.value, float):
                    bad.append((node.name, child.lineno))
    return bad


print("T-54 VERIFIER B RE-RUN -- classes + writing (default REFUTED)")
print("=" * 96)

# The five source lanes are read, never re-executed.
sealed_lanes = [
    ROOT / "LANE_T44_A_CORNER", ROOT / "LANE_T44_B_WORLD",
    ROOT / "LANE_T48_A_DERIVATION", ROOT / "LANE_T48_B_CORNER",
    ROOT / "LANE_T48_C_WORLD",
]
seal_results = [manifest_is_intact(lane) for lane in sealed_lanes]
check("shared", "S-SEALS", all(ok for ok, _ in seal_results),
      "five input manifests intact; %s entries checked" %
      "/".join(str(count) for _, count in seal_results))

self_tree = ast.parse(Path(__file__).read_text(), filename=__file__)
handlers = [node for node in ast.walk(self_tree) if isinstance(node, ast.ExceptHandler)]
handler_names = [node.type.id if isinstance(node.type, ast.Name) else None for node in handlers]
check("shared", "S-EXCEPTION-HYGIENE", handler_names == ["ValueError"],
      f"exception handlers={handler_names}; positive controls are unwrapped")

# Full builder suites, captured rather than echoed so output is deterministic.
cl_rc, cl_suite = suite_result("checks_classes.py", "CLASSES")
wr_rc, wr_suite = suite_result("checks_writing.py", "WRITING")
check("classes", "CL-SUITE", cl_rc == 0 and cl_suite == ("52", "0"),
      f"exit={cl_rc} summary={cl_suite}")
check("writing", "WR-SUITE", wr_rc == 0 and wr_suite == ("57", "0"),
      f"exit={wr_rc} summary={wr_suite}")

# ---------------------------------------------------------------- sealed inputs, parsed directly
t44 = (ROOT / "LANE_T44_B_WORLD" / "t44b_world.OUT.txt").read_text()
t44_adv = (ROOT / "LANE_T44_B_WORLD" / "VERIFY" / "adv_verify.OUT.txt").read_text()
t48a = (ROOT / "LANE_T48_A_DERIVATION" / "t48a_derivation.OUT.txt").read_text()
t48b = (ROOT / "LANE_T48_B_CORNER" / "t48b_corner.OUT.txt").read_text()
t48c = (ROOT / "LANE_T48_C_WORLD" / "t48c_world.OUT.txt").read_text()

g0_anchor = capture_groups(
    r"S4 G\(0\).*bracket \[([0-9.]+), ([0-9.]+)\]", t44)
r24_anchor = capture_groups(
    r"axis doubling pair d=2->4.*\[([0-9.]+), ([0-9.]+)\]", t44)
r48_anchor = capture_groups(
    r"axis doubling pair d=4->8.*\[([0-9.]+), ([0-9.]+)\]", t44)
deep_ratio_anchor = capture_groups(
    r"REFUTE  E6[^\n]*\('axis', 4, 8\): \('([0-9.]+)', '([0-9.]+)'\)",
    t44_adv)
deep_coef_anchor = capture_groups(
    r"c_ax\(M=2800\)=\[([0-9.]+),([0-9.]+)\]", t44_adv)

# ---------------------------------------------------------------- classes: independent D=3 and fresh APIs
dp = independent_dp3(14)
targets = ((0, 0, 0), (2, 0, 0), (4, 0, 0), (2, 2, 0), (2, 2, 2))
small_exact = all(
    CL.n3_even_row(m, *target) == dp[2 * m].get(target, 0)
    for m in range(8) for target in targets
)
check("classes", "CL-DP", small_exact,
      "fresh Z3 dictionary DP == ported counts; five targets, m=0..7")

w4, _cells4, _idx4 = CL.world_venue(4)
loc4 = CL.mu_c_of(w4, certify="full")
check("classes", "CL-MUC",
      loc4.get("located") and loc4["mu_c"] == F(1, 6)
      and loc4["resolvent_singular_at_mu_c"]
      and loc4["solvable_below"] and loc4["solvable_above"],
      f"fresh exact resolvent route mu_c={loc4.get('mu_c')}")

A2, A4, A8, A16 = (2, 0, 0), (4, 0, 0), (8, 0, 0), (16, 0, 0)
tc = CL.tail_constants()
B5 = tc["B5"]
kp = CL.kernel_pass({A2: 2800, A4: 2800, A8: 2800, A16: 2800}, 2800,
                    snapshots=(350, 700, 1400))
ker1400 = kp["ker_at"][1400]
tails1400 = {t: CL.diff_tail_bound(t, 1400, B5) for t in (A2, A4, A8, A16)}
g0_1400 = (kp["s0"][1400],
           kp["s0"][1400] + CL.abs_tail_bound(1400, B5, kp["p2m0"][1400]))
r24 = CL.doubling_ratio(ker1400, tails1400, (A2, A4), (A4, A8))
r48 = CL.doubling_ratio(ker1400, tails1400, (A4, A8), (A8, A16))
tails2800 = {t: CL.diff_tail_bound(t, 2800, B5) for t in (A2, A4, A8, A16)}
r48_deep = CL.doubling_ratio(kp["ker"], tails2800, (A4, A8), (A8, A16))
h816 = CL.increment_interval(kp["ker"], tails2800, A8, A16)
coef2800 = (16 * h816[0], 16 * h816[1])

fresh_t44 = ((CL.ff(g0_1400[0]), CL.ff(g0_1400[1])),
             (CL.ff(r24[0]), CL.ff(r24[1])),
             (CL.ff(r48[0]), CL.ff(r48[1])),
             (CL.ff(r48_deep[0]), CL.ff(r48_deep[1])),
             (CL.ff(coef2800[0]), CL.ff(coef2800[1])))
check("classes", "CL-SEALED",
      None not in (g0_anchor, r24_anchor, r48_anchor, deep_ratio_anchor, deep_coef_anchor)
      and fresh_t44 == (g0_anchor, r24_anchor, r48_anchor,
                        deep_ratio_anchor, deep_coef_anchor),
      "fresh=" + repr(fresh_t44))
check("classes", "CL-KERNEL",
      r48_deep[1] < F(1, 2) and coef2800[0] <= CL.C_3D[0] <= coef2800[1],
      "registered deep caveat retained: r48_hi<1/2; owner coefficient inside")
check("classes", "CL-EXACT",
      all(isinstance(x, F) for x in
          (tc["Q3"], B5, *g0_1400, *r24, *r48, *r48_deep, *coef2800))
      and not float_literals_in_functions(ROOT / "model" / "classes.py", {
          "_coupling_price", "venue", "walk_counts", "venue_series",
          "resolvent_exact", "annihilates_constant", "mu_c_of", "series_3d",
          "subcritical_row", "divergence_witness", "class_verdict",
      }),
      "ints/Fractions on sampled and static core measurement paths")

# Off-sealed K4: independent transfer action provides the expected series.
k4_adj = [[(j, 1) for j in range(4) if j != i] for i in range(4)]
k4 = CL.venue("off-sealed K4", k4_adj, provenance="verifier synthetic graph")
k4_loc = CL.mu_c_of(k4, certify="full")
k4_api = CL.venue_series(k4, F(1, 7), 0, range(4), 12)
k4_counts = independent_walk_counts(k4_adj, 0, 12)
k4_direct = {t: sum(F(k4_counts[k][t], 7 ** k) for k in range(13))
             for t in range(4)}
check("classes", "CL-OFF",
      k4_loc.get("located") and k4_loc["mu_c"] == F(1, 3)
      and all(k4_api[t][0] == k4_direct[t] for t in range(4))
      and CL.class_verdict(k4, F(1, 7))["exponential"],
      f"off-sealed K4 mu_c={k4_loc.get('mu_c')}; exact independent series")

neg_series_ok, neg_series_detail = expected_value_error(
    lambda: CL.venue_series(CL.chain_venue(10), -F(1, 4), 0, [1], 4),
    "coupling price mu must be nonnegative")
neg_class_ok, _ = expected_value_error(
    lambda: CL.class_verdict(CL.chain_venue(10), -F(1, 4)),
    "coupling price mu must be nonnegative")
chain = CL.chain_venue(10)
zero_price = CL.venue_series(chain, F(0), 0, [1], 4)[1]
inside_price = CL.venue_series(chain, F(1, 4), 0, [1], 4)[1]
check("classes", "CL-REFUSE-PRICE",
      neg_series_ok and neg_class_ok and zero_price == (F(0), F(0))
      and inside_price[0] >= 0 and inside_price[1] >= 0,
      neg_series_detail + "; valid neighbors mu=0,1/4 returned exact bounds")

bad_adj_ok, bad_adj_detail = expected_value_error(
    lambda: CL.venue("negative multiplicity", [[(1, -1)], [(0, 1)]],
                     provenance="verifier invalid graph"),
    "edge multiplicities must be strictly positive integers")
good_graph = CL.venue("positive multiplicity", [[(1, 1)], [(0, 1)]],
                      provenance="verifier valid neighboring graph")
good_loc = CL.mu_c_of(good_graph, certify="full")
check("classes", "CL-REFUSE-VENUE",
      bad_adj_ok and good_graph.row_sums() == [1, 1]
      and good_loc.get("located") and good_loc["mu_c"] == F(1),
      bad_adj_detail + "; +1 neighbor accepted with mu_c=1")

check("classes", "CL-D8",
      not literal_check_conditions(ROOT / "model" / "checks_classes.py"),
      "AST scan: no check(..., True/False) or assert True/False")
lead_zero = all(dp[k].get((2, 0, 0), 0) == 0 for k in (0, 1))
lead_positive = dp[2].get((2, 0, 0), 0)
check("classes", "CL-D15",
      loc4["resolvent_singular_at_mu_c"] and loc4["solvable_below"]
      and lead_zero and lead_positive == 1 and r48_deep[1] < F(1, 2),
      "singular/solvable, below-distance/leading-count, and deep-caveat controls computed")

# ---------------------------------------------------------------- writing: sealed values through fresh APIs
move1_anchor = capture_groups(r"S1\.4 moving .*: (\d+) of (\d+)", t48a)
move2_anchor = capture_groups(r"S1\.11 moving .*: (\d+) of (\d+)", t48a)
det_anchor = capture_groups(r"S3\.4 C8: det\(I - \(9/10\)K\) = ([0-9/]+)", t48a)
tuple_anchor = capture_groups(r"invariant tuple .* = \(([^)]+)\) on all (\d+) links", t48b)
e1_anchor = capture_groups(r"induced mu set = \{([^}]+)\}", t48c)
e2_anchor = capture_groups(r"b=1/2: \(m_fresh, m_back\) = \(([^,]+), ([^)]+)\)", t48c)
e3_anchor = capture_groups(r"'1/2': '([^']+)'", t48c)

ring8 = WR.ring_venue(8)
t3 = WR.plaquette_venue(3)
cen1 = WR.moving_census(ring8["deltas"], 8, even_only=False)
cen2 = WR.moving_census(t3["deltas"], 9, even_only=True)
leak_det = WR.crit_det(WR.leak_kernel(ring8["adj"], F(9, 10)))
check("writing", "WR-SEALED-A",
      move1_anchor is not None and move2_anchor is not None and det_anchor is not None
      and (str(cen1["moving"]), str(cen1["total"])) == move1_anchor
      and (str(cen2["moving"]), str(cen2["total"])) == move2_anchor
      and str(leak_det) == det_anchor[0],
      f"moving={cen1['moving']}/{cen1['total']},{cen2['moving']}/{cen2['total']} det={leak_det}")

cv46 = WR.corner_venue(4, 6)
inv46 = WR.writer_invariants(cv46)
tuple_api = ", ".join(str(x) for x in inv46[0])
check("writing", "WR-SEALED-B",
      tuple_anchor is not None and tuple_api == tuple_anchor[0]
      and str(len(inv46)) == tuple_anchor[1]
      and all(row == inv46[0] for row in inv46),
      f"tuple=({tuple_api}) on {len(inv46)} links")

_cells8, _idx8, nbr8 = WR.torus3(8)
v1_sealed = WR.transport_verdict(nbr8, WR.ensemble_transport(nbr8, F(1, 40)))
v2_sealed = WR.retreat_verdict(
    nbr8, WR.ensemble_trail_retreat(nbr8, F(1, 20), F(1, 2)))
v3_sealed = WR.decay_verdict(
    nbr8, WR.ensemble_trail_decay(nbr8, F(1, 20), F(1, 2), "H1"))
check("writing", "WR-SEALED-C",
      e1_anchor is not None and e2_anchor is not None and e3_anchor is not None
      and str(v1_sealed["mu"]) == e1_anchor[0]
      and (str(v2_sealed["m_fresh"]), str(v2_sealed["m_back"])) == e2_anchor
      and str(v3_sealed["mu"]) == e3_anchor[0]
      and v1_sealed["at_criticality"] and v2_sealed["doubly_stochastic"]
      and v3_sealed["below_criticality"],
      f"E1 mu={v1_sealed['mu']} E2={v2_sealed['m_fresh']},{v2_sealed['m_back']} E3 mu={v3_sealed['mu']}")

# Off every sealed world grid: 5^3, with fresh u,b and both counting rules.
_cells5, _idx5, nbr5 = WR.torus3(5)
u, b = F(3, 80), F(5, 11)
W1 = WR.ensemble_transport(nbr5, u * b)
W2 = WR.ensemble_trail_retreat(nbr5, u, b)
W3 = WR.ensemble_trail_decay(nbr5, u, b, "H1")
W3nb = WR.ensemble_trail_decay(nbr5, u, b, "NB")
v1, v2 = WR.transport_verdict(nbr5, W1), WR.retreat_verdict(nbr5, W2)
v3, v3nb = WR.decay_verdict(nbr5, W3), WR.decay_verdict_nb(nbr5, W3nb)
check("writing", "WR-OFF",
      all_fractions(W1) and all_fractions(W2) and all_fractions(W3) and all_fractions(W3nb)
      and v1["at_criticality"] and v1["mu"] == F(1, 6)
      and v2["doubly_stochastic"]
      and (v2["m_fresh"], v2["m_back"]) == (F(5, 36), F(11, 36))
      and v3["mu"] == F(5, 41) and v3["mass_ratio"] == F(41, 30)
      and v3nb["mu"] == F(5, 36) and v3nb["mass_ratio"] == F(36, 25),
      "off-sealed 5^3 exact E1/E2/E3-H1/E3-NB identities")

# D-15: the exact zero determinant must sit beside a nonzero leak on the same ring.
ring10 = WR.ring_venue(10)
Kc = WR.kernel_uniform(ring10["adj"], F(2, 9))
Kl = WR.leak_kernel(ring10["adj"], F(17, 19))
te = WR.trail_energetics(nbr5)
check("writing", "WR-D15",
      WR.crit_det(Kc) == 0 and WR.crit_det(Kl) != 0
      and te["transport"] == [0] * 6
      and te["write"] == [1] * 5 and te["erase"] == [-1]
      and v1["conserving"] and v3["loss"] == u,
      "det zero/leak, transport zero/write-erase, and no-loss/loss controls computed")

# Exact intended refusals plus unwrapped boundary and just-inside valid controls.
e1_neg, _ = expected_value_error(
    lambda: WR.ensemble_transport(nbr5, -F(1, 100)), "writing E1 REFUSES")
e1_over, e1_detail = expected_value_error(
    lambda: WR.ensemble_transport(nbr5, F(1, 5)), "writing E1 REFUSES")
e1_boundary = WR.transport_verdict(nbr5, WR.ensemble_transport(nbr5, F(1, 6)))
e1_inside = WR.transport_verdict(nbr5, WR.ensemble_transport(nbr5, F(1, 7)))
check("writing", "WR-REFUSE-E1",
      e1_neg and e1_over and e1_boundary["nonnegative"]
      and e1_boundary["at_criticality"] and e1_inside["at_criticality"],
      e1_detail + "; valid a=1/6 and 1/7 controls returned")

e2_negu, _ = expected_value_error(
    lambda: WR.ensemble_trail_retreat(nbr5, -F(1, 20), F(1)), "writing E2 REFUSES")
e2_negb, _ = expected_value_error(
    lambda: WR.ensemble_trail_retreat(nbr5, F(1, 20), -F(1)), "writing E2 REFUSES")
e2_over, e2_detail = expected_value_error(
    lambda: WR.ensemble_trail_retreat(nbr5, F(1, 2), F(1)), "writing E2 REFUSES")
e2_boundary = WR.retreat_verdict(
    nbr5, WR.ensemble_trail_retreat(nbr5, F(1, 6), F(1)))
e2_inside = WR.retreat_verdict(
    nbr5, WR.ensemble_trail_retreat(nbr5, F(1, 7), F(1)))
check("writing", "WR-REFUSE-E2",
      e2_negu and e2_negb and e2_over
      and e2_boundary["nonnegative"] and e2_boundary["doubly_stochastic"]
      and e2_inside["doubly_stochastic"],
      e2_detail + "; valid u=1/6 and 1/7 controls returned")

e3_bad_name, _ = expected_value_error(
    lambda: WR.ensemble_trail_decay(nbr5, F(1, 20), F(1), "mystery"),
    "counting must be exactly 'H1' or 'NB'")
e3_over_h1, _ = expected_value_error(
    lambda: WR.ensemble_trail_decay(nbr5, F(1, 2), F(1), "H1"), "writing E3 REFUSES")
e3_over_nb, e3_detail = expected_value_error(
    lambda: WR.ensemble_trail_decay(nbr5, F(1, 2), F(1), "NB"), "writing E3 REFUSES")
e3_h_boundary = WR.decay_verdict(
    nbr5, WR.ensemble_trail_decay(nbr5, F(1, 7), F(1), "H1"))
e3_n_boundary = WR.decay_verdict_nb(
    nbr5, WR.ensemble_trail_decay(nbr5, F(1, 6), F(1), "NB"))
check("writing", "WR-REFUSE-E3",
      e3_bad_name and e3_over_h1 and e3_over_nb
      and e3_h_boundary["nonnegative"] and e3_h_boundary["substochastic"]
      and e3_n_boundary["nonnegative"] and e3_n_boundary["substochastic"],
      e3_detail + "; valid H1/NB zero-stay boundaries returned")

signed_e1 = WR.kernel_pos(nbr5, -F(1, 5), F(1, 5))
signed_e2 = WR.kernel_edge(nbr5, -F(2), F(1, 2), F(1, 2))
signed_v1 = WR.transport_verdict(nbr5, signed_e1)
signed_v2 = WR.retreat_verdict(nbr5, signed_e2)
check("writing", "WR-SIGNED",
      set(WR.srow_sums(signed_e1)) == {F(1)}
      and set(WR.srow_sums(signed_e2)) == {F(1)}
      and not signed_v1["nonnegative"] and not signed_v1["conserving"]
      and not signed_v1["at_criticality"]
      and not signed_v2["nonnegative"] and not signed_v2["doubly_stochastic"],
      "unit row sums cannot hide signed entries")

mystery = URM.surface("verifier provenance probe", "synthetic", 1e-20, 1e-19,
                      300.0, 1e9, provenance="verifier source removed below")
del mystery.provenance
missing_prov, prov_detail = expected_value_error(
    lambda: WR.surface_gap(mystery),
    "world-tier surface must enter through the D-25 provenance gate")
eV = 1.602176634e-19
nand = URM.surface("NAND floating gate", "trapped charge", 0.30 * eV, 1.60 * eV,
                   300.0, 1e13)
nand_gap = WR.surface_gap(nand)
extreme = URM.surface("extreme verifier surface", "synthetic", 1e-18, 1e-21,
                      1.0, 1e9, provenance="verifier underflow control")
extreme_dial = WR.surface_boltzmann(extreme)
extreme_gap = WR.surface_gap(extreme)
check("writing", "WR-SURFACE",
      missing_prov and nand_gap is not None and nand_gap["contained"]
      and nand_gap["closed_form_agrees"] and nand_gap["u_independent"]
      and extreme_dial is not None and extreme_dial["b_underflow"]
      and extreme_dial["b"] == 0.0 and extreme_gap is None,
      prov_detail + "; resolved NAND returned; underflow surface declined")

check("writing", "WR-EXACT",
      all(isinstance(x, F) for x in (
          v1["mu"], v2["m_fresh"], v2["m_back"], v3["mu"],
          v3["mass_ratio"], v3nb["mu"], v3nb["mass_ratio"], leak_det))
      and not float_literals_in_functions(ROOT / "model" / "writing.py", {
          "is_stochastic", "is_doubly_stochastic", "dict_doubly_stochastic",
          "ensemble_transport", "ensemble_trail_retreat", "ensemble_trail_decay",
          "transport_verdict", "retreat_verdict", "decay_verdict", "decay_verdict_nb",
          "closed_form_gap_ratio",
      }),
      "Fractions dynamically and no float literal in static core measurement functions")
check("writing", "WR-D8",
      not literal_check_conditions(ROOT / "model" / "checks_writing.py"),
      "AST scan: no check(..., True/False) or assert True/False")

# ---------------------------------------------------------------- terminal default-refuted decision
for group, required in REQUIRED.items():
    missing = sorted(required - seen[group])
    if missing:
        failures[group] += len(missing)
        for ident in missing:
            print(f"REFUTE  {ident}  required predicate was not executed")

classes_ok = failures["shared"] == 0 and failures["classes"] == 0
writing_ok = failures["shared"] == 0 and failures["writing"] == 0
classes_verdict = "CONFIRMED" if classes_ok else "REFUTED"
writing_verdict = "CONFIRMED" if writing_ok else "REFUTED"
total_pass = sum(passes.values())
total_refute = sum(failures.values())

print("=" * 96)
print(f"VERIFIER B: {total_pass} OK, {total_refute} REFUTE, "
      f"{total_pass + total_refute} required predicates")
print(f"FAMILY VERDICTS: classes={classes_verdict}; writing={writing_verdict}")
sys.exit(0 if classes_ok and writing_ok else 1)
