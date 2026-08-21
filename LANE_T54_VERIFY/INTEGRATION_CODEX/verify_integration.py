#!/usr/bin/env python3
"""Independent, default-refuted integration verifier for T-54/T-55.

This program does not run any sealed lane program and does not rerun validate_urm.py.
The one umbrella execution is tied to its source hashes in UMBRELLA.OUT.txt.  These
predicates exercise the integrated live API off the sealed grids, its refusal paths,
delegate wiring, validator topology, and the already-produced verifier artifacts.
"""
from __future__ import annotations

import ast
import hashlib
import inspect
import sys
from fractions import Fraction as Fr
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "model"
sys.path.insert(0, str(MODEL))

import arrow as AW
import classes as CC
import countlaw as CL
import project_model as PM
import writing as WW
from checks_d25 import scan_direct_constructions
from record_model import Environment, RecordModel


DEFAULT_VERDICT = "REFUTED"
REQUIRED = 25
results: list[tuple[str, bool, str]] = []


def predicate(name, fn):
    try:
        detail = fn()
        if detail is False:
            raise AssertionError("predicate returned False")
        results.append((name, True, str(detail or "computed predicate passed")))
    except Exception as exc:  # an unexpected exception is evidence against integration
        results.append((name, False, f"{type(exc).__name__}: {exc}"))


def expect_value_error(call, marker):
    try:
        call()
    except ValueError as exc:
        assert marker in str(exc), f"wrong ValueError contract: {exc}"
        return str(exc)
    raise AssertionError("required ValueError was not raised")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def p01_contract():
    doc = inspect.getdoc(PM) or ""
    gate = inspect.getdoc(PM.URM) or ""
    for token in ("new surface", "new law", "new external number", "validator gate",
                  "claim row", "stated tolerance", "power control"):
        assert token in doc, token
    assert "surface() REFUSES" in gate and "corner() REFUSES" in gate
    return "three doors and both public refusal constructors are explicit"


def p02_api_surface():
    required = {
        "arrow_threshold", "arrow_ledger", "arrow_invariance", "arrow_history",
        "arrow_redundancy", "arrow_observation", "census", "count_widths",
        "coupling_venue", "world_coupling_venue", "corner_coupling_venue",
        "chain_coupling_venue", "critical_price", "reachable_class", "coupling",
        "critical_kernel", "class_discriminator", "writing_kernel_verdict",
        "writing_uniformity", "writing_transport", "writing_trail_retreat",
        "writing_trail_decay", "writing_gap",
    }
    missing = sorted(n for n in required if not callable(getattr(PM.ProjectModel, n, None)))
    assert not missing, missing
    assert issubclass(PM.URM, PM.ProjectModel)
    return f"{len(required)} integrated layer methods reachable through URM"


def p03_surface_refusal():
    m = PM.URM()
    msgs = [
        expect_value_error(lambda: m.surface("fresh", "x", 1., 2., 300., 1e9),
                           "PROVENANCE"),
        expect_value_error(lambda: m.surface("fresh", "x", 1., 2., 300., 1e9,
                                             provenance="   "), "PROVENANCE"),
        expect_value_error(lambda: m.surface("corner", "x", 0., 0., 0., 0.,
                                             tier="corner"), "DEF-A"),
        expect_value_error(lambda: m.surface("corner", "x", 0., 0., 0., 0.,
                                             tier="corner", provenance="def-a"), "DEF-A"),
    ]
    assert all(msgs)
    return "blank/whitespace world and absent/lowercase corner declarations refused"


def p04_surface_acceptance():
    m = PM.URM()
    reg = m.surface("NAND floating gate", "trapped charge", .03 * CL.EV, .9 * CL.EV,
                    321., 2e9)
    ext = m.surface("fresh device", "thermal", .04 * CL.EV, .8 * CL.EV, 333., 3e9,
                    provenance="pinned off-sealed integration specimen")
    cor = m.surface("fresh corner", "thermal", .02 * CL.EV, .7 * CL.EV, 311., 7e8,
                    tier="corner", provenance="DEF-A")
    assert reg.tier == ext.tier == "world" and reg.provenance and ext.provenance
    assert cor.tier == "corner" and cor.provenance == "DEF-A"
    return "registry, explicit world, and exact DEF-A entries accepted"


def p05_corner_refusal_and_positive():
    m = PM.URM()
    H = np.diag([-1.3, .2, 2.7])
    expect_value_error(lambda: m.corner(H), "DEF-A")
    expect_value_error(lambda: m.corner(H, provenance="def-a"), "DEF-A")
    got = m.corner(H, provenance="DEF-A")
    assert got == {"slow_dim": 3, "commutant_dim": 3}
    return "URM.corner refuses two bad declarations beside a computed DEF-A result"


def p06_d25_model_scan():
    off = scan_direct_constructions(root=str(MODEL))
    assert off == [], off
    return "AST scan found no RecordSurface binding/reference outside project_model.py"


def p07_countlaw_raw_bypass():
    m = PM.URM()
    raw = PM.RecordSurface("raw", "thermal", .04 * CL.EV, .8 * CL.EV, 333., 3e9)
    expect_value_error(lambda: m.census([raw], 10.), "COUNTLAW REFUSES")
    expect_value_error(lambda: m.count_widths(raw, 10.), "COUNTLAW REFUSES")
    raw.tier, raw.provenance = "corner", "not DEF-A"
    expect_value_error(lambda: CL.census([raw], 10.), "DEF-A")
    return "raw and false-corner constructor bypasses refused at consumption"


def p08_countlaw_mutation_bypass():
    m = PM.URM()
    s = m.surface("NAND floating gate", "trapped charge", .04 * CL.EV, .8 * CL.EV,
                  333., 3e9)
    s.provenance = None
    for fn in (lambda: CL.record_rate(s), lambda: CL.drop_time(s),
               lambda: CL.drop_time_formula(s), lambda: CL.census([s], 10.)):
        expect_value_error(fn, "COUNTLAW REFUSES")
    return "all four public surface-consuming paths recheck removed provenance"


def fresh_countlaw_pair():
    m = PM.URM()
    s = m.surface("fresh census grain", "thermal", .07 * CL.EV, .87 * CL.EV,
                  333., 3e10, provenance="pinned off-sealed census specimen")
    return s, m.census([s], 10.), m.census([s], 1e6)


def p09_countlaw_offsealed():
    s, short, long = fresh_countlaw_pair()
    assert short["k"] == short["k_formula"] == 1
    assert long["k"] == long["k_formula"] == 0
    row = short["schedule"][0]
    assert row["provenance"] == s.provenance and row["tier"] == "world"
    assert abs(row["t_star"] / row["t_star_formula"] - 1) < 1e-12
    return "fresh 333 K record computes 1->0 staircase with both routes equal"


def p10_countlaw_wrapper_fidelity():
    m = PM.URM()
    s, short, _ = fresh_countlaw_pair()
    assert m.census([s], 10.) == CL.census([s], 10.) == short
    got = m.count_widths(s, 10.)
    want = dict(delta_pop=CL.delta_pop(s.E_b + s.dE, s.T, s.f0, 10.),
                delta_coh=CL.delta_coh(10.), t_star=CL.drop_time(s),
                t_star_formula=CL.drop_time_formula(s))
    assert got == want
    return "census and both-width wrapper equal the family API off sealed grids"


def tiny_arrow_rows():
    m = PM.URM()
    X = np.array([[0., 1.], [1., 0.]])
    Z = np.diag([1., -1.])
    rm = RecordModel(np.zeros((2, 2)), [])
    env = Environment(nq=1, energies=(.9,), beta=1.7)
    args = dict(model=rm, lam=.37, t=.83, tier="corner", provenance="DEF-A")
    no_record = m.arrow_observation(env, X, record=Z, **args)
    holds = m.arrow_observation(env, X, record=X, **args)
    return m, env, rm, X, Z, args, no_record, holds


def p11_arrow_observation_offsealed():
    m, env, rm, X, Z, args, no_record, holds = tiny_arrow_rows()
    expect_value_error(lambda: m.arrow_observation(env, X, record=Z, model=rm),
                       "PROVENANCE")
    assert no_record["entangled_without_record"] and not no_record["holds_record_bits"]
    assert holds["holds_record_bits"] and not holds["entangled_without_record"]
    assert no_record["I_SB"] > 0 and no_record["chi_whole"] == 0
    direct = AW.score_bath_observation(env, X, record=Z, **args)
    for key in ("I_SB", "chi_whole", "holds_record_bits", "entangled_without_record",
                "redundant_fragments"):
        assert no_record[key] == direct[key]
    return "fresh one-qubit venue reaches record/no-record branches and wrapper equals family"


def p12_arrow_delegate_wiring():
    m = PM.URM()
    coupling, env = object(), object()
    sentinels = [object() for _ in range(5)]
    seen = {}

    def stub(tag, value):
        def f(*args, **kwargs):
            seen[tag] = (args, kwargs)
            return value
        return f

    with patch.object(AW, "arrow_threshold", stub("threshold", sentinels[0])), \
         patch.object(AW, "arrow_ledger", stub("ledger", sentinels[1])), \
         patch.object(AW, "arrow_invariance", stub("invariance", sentinels[2])), \
         patch.object(AW, "arrow_history", stub("history", sentinels[3])), \
         patch.object(AW, "arrow_redundancy", stub("redundancy", sentinels[4])):
        assert m.arrow_threshold(lam=.37, weights=(3,), coupling=coupling) is sentinels[0]
        assert m.arrow_ledger(lam=.41) is sentinels[1]
        assert m.arrow_invariance(n_unitaries=3, seed=17, lam=.43) is sentinels[2]
        assert m.arrow_history((-.2, .7), coupling=coupling, lam=.47, env=env,
                               keep_states=(.7,)) is sentinels[3]
        assert m.arrow_redundancy(coupling=coupling, lam=.53, t=.91,
                                  env=env) is sentinels[4]
    assert seen["threshold"] == ((), {"lam": .37, "weights": (3,), "coupling": coupling})
    assert seen["ledger"] == ((), {"lam": .41})
    assert seen["invariance"] == ((), {"n_unitaries": 3, "seed": 17, "lam": .43})
    assert seen["history"][1]["keep_states"] == (.7,) and seen["history"][1]["env"] is env
    assert seen["redundancy"][1] == {"coupling": coupling, "lam": .53, "t": .91, "env": env}
    return "five non-observation arrow delegates preserve off-default arguments exactly"


def p13_classes_entry_refusal():
    m = PM.URM()
    adj = [[((i + 1) % 5, 1), ((i - 1) % 5, 1)] for i in range(5)]
    expect_value_error(lambda: m.coupling_venue("fresh C5", adj), "provenance")
    expect_value_error(lambda: m.coupling_venue("fresh C5", adj, tier="corner",
                                                provenance="def-a"), "DEF-A")
    v = m.coupling_venue("fresh C5", adj, provenance="off-sealed graph source")
    assert v.n == 5 and v.tier == "world" and v.provenance
    return "class venue gate refuses two undeclared cases beside a valid fresh C5"


def p14_classes_offsealed():
    m = PM.URM()
    v = m.chain_coupling_venue(11)
    loc = m.critical_price(v)
    rows = [m.reachable_class(v, mu) for mu in (Fr(1, 3), Fr(1, 2), Fr(2, 3))]
    assert v.n == 11 and loc["located"] and loc["mu_c"] == Fr(1, 2)
    assert [(r["exponential"], r["critical"], r["divergent"]) for r in rows] == [
        (True, False, False), (False, True, False), (False, False, True)]
    assert m.coupling(Fr(1, 9), (2, 1, 0), 13) == CC.series_3d(Fr(1, 9), (2, 1, 0), 13)
    targets = [(2, 0, 0), (4, 0, 0)]
    assert m.critical_kernel(targets, 17) == CC.crit_kernel_3d(targets, 17)
    return "fresh C11 reaches all three computed classes; two layer wrappers equal family"


def p15_classes_constructor_delegates():
    m = PM.URM()
    w = m.world_coupling_venue(5)
    c, ccells, cidx = m.corner_coupling_venue(4, 5)
    wd, cells, idx = w
    wd0, cells0, idx0 = CC.world_venue(5)
    c0, ccells0, cidx0 = CC.corner_venue(4, 5)
    assert cells == cells0 and idx == idx0 and wd.adj == wd0.adj
    assert ccells == ccells0 and cidx == cidx0 and c.adj == c0.adj
    assert c.provenance == c0.provenance == "DEF-A"
    sentinel, seen = object(), {}

    def disc(**kwargs):
        seen.update(kwargs)
        return sentinel

    with patch.object(CC, "discriminator", disc):
        assert m.class_discriminator(K2=731, K1=1913) is sentinel
    assert seen == {"K2": 731, "K1": 1913}
    return "fresh world/corner constructors and discriminator arguments preserve family API"


def p16_writing_kernel_domain():
    m = PM.URM()
    expect_value_error(lambda: m.writing_kernel_verdict("C8", Fr(-1, 10)), "0 <= c <= 1")
    expect_value_error(lambda: m.writing_kernel_verdict("C8", Fr(11, 10)), "0 <= c <= 1")
    end = m.writing_kernel_verdict("C8", Fr(1))
    mid = m.writing_kernel_verdict("C8", Fr(3, 7))
    assert end["conserving"] and end["critical"] and end["per_crossing"] is None
    assert end["link_amplitudes"] == {Fr(0)} and end["leak_det"] != 0
    assert mid["conserving"] and mid["critical"] and mid["per_crossing"] == Fr(1, 2)
    return "outside [0,1] refused; c=1 is critical with no crossing; interior computed"


def p17_writing_offsealed():
    m = PM.URM()
    _, _, nbr = WW.torus3(5)
    a, u, b = Fr(1, 13), Fr(1, 17), Fr(3, 7)
    assert m.writing_transport(5, a) == WW.transport_verdict(nbr, WW.ensemble_transport(nbr, a))
    assert m.writing_trail_retreat(5, u, b) == WW.retreat_verdict(
        nbr, WW.ensemble_trail_retreat(nbr, u, b))
    assert m.writing_trail_decay(5, u, b, "H1") == WW.decay_verdict(
        nbr, WW.ensemble_trail_decay(nbr, u, b, "H1"))
    assert m.writing_trail_decay(5, u, b, "NB") == WW.decay_verdict_nb(
        nbr, WW.ensemble_trail_decay(nbr, u, b, "NB"))
    return "all four world writer delegates equal family machinery on fresh 5^3 venue"


def p18_writing_corner_and_gap():
    m = PM.URM()
    uni = m.writing_uniformity(3, 5)
    assert uni["identical"] and uni["n_links"] == 30 and uni["mu_c"]["mu_c"] == Fr(1, 4)
    s = m.surface("fresh gap surface", "thermal", .04 * CL.EV, .83 * CL.EV, 329., 4e9,
                  provenance="pinned off-sealed writing specimen")
    got = m.writing_gap(s, n=5, den=10 ** 6)
    assert got == WW.surface_gap(s, n=5, den=10 ** 6) and got is not None
    raw = PM.RecordSurface("raw gap", "thermal", .04 * CL.EV, .83 * CL.EV, 329., 4e9)
    expect_value_error(lambda: m.writing_gap(raw, n=5), "D-25")
    return "fresh (3,5) corner and surface-gap wrappers computed; raw gap refused"


def p19_validator_topology():
    vu = (MODEL / "validate_urm.py").read_text()
    vg = (MODEL / "validate_geometry.py").read_text()
    vp = (MODEL / "validate_project.py").read_text()
    assert all(f'("{name}", {n},' in vu for name, n in
               (("ARROW", 27), ("COUNTLAW", 40), ("CLASSES", 52), ("WRITING", 57)))
    assert "gates {n_pass + n_fail}/176" in vu and '"validate_geometry.py"' in vu
    assert '"validate_project.py"' in vg and "all 24 project/D-25 gates" in vg
    assert "from project_model import URM" in vp and "run_d25_checks(check)" in vp
    assert "RecordSurface" not in ast.dump(ast.parse(vp))
    return "176-family -> geometry -> 24-project/D25 chain is explicit and bypass-free"


def p20_umbrella_capture():
    text = (Path(__file__).with_name("UMBRELLA.OUT.txt")).read_text()
    required = (
        "Exit: 0", "ARROW: 27 PASS, 0 FAIL", "COUNTLAW: 40 PASS, 0 FAIL",
        "CLASSES: 52 PASS, 0 FAIL", "WRITING: 57 PASS, 0 FAIL",
        "T-54 FAMILIES: 176 PASS, 0 FAIL", "GEOMETRY: 33 PASS, 0 FAIL",
        "PROJECT/D-25: 24 PASS, 0 FAIL", "URM OVERALL: PASS",
    )
    assert all(x in text for x in required)
    for line in text.splitlines():
        if "  model/" in line:
            digest, rel = line.split("  ", 1)
            assert sha256(ROOT / rel) == digest, f"umbrella snapshot changed: {rel}"
    proof = (Path(__file__).with_name("PROOF_GATE.OUT.txt")).read_text()
    assert "GATE PASSED" in proof and "R11 PASS  83 blocks; 115 distinct rows cited" in proof
    assert sha256(ROOT / "PROOF_V002.md") == "f5496567630ddaeb983fa070f43e6f1c3b546cc9b4c070054c39de2eff2105c9"
    assert sha256(ROOT / "replicate/check_proof.py") == "00b72d58bb9ff896277d437ac87384b29bc67bee5f23d0065942e31f13ed1bd3"
    return ("one umbrella passes exact 27/40/52/57, 33, 24 counts; proof R1-R11 passes "
            "at 83 blocks/115 rows; all captures match live hashes")


def p21_d8_static_scan():
    offenders = []
    paths = [MODEL / f for f in ("checks_arrow.py", "checks_countlaw.py",
                                  "checks_classes.py", "checks_writing.py", "checks_d25.py")]
    for path in paths:
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) \
                    and node.func.id == "check" and len(node.args) >= 2:
                if isinstance(node.args[1], ast.Constant) and isinstance(node.args[1].value, bool):
                    offenders.append((path.name, node.lineno, node.args[1].value))
    assert offenders == [], offenders
    return "five integrated check blocks contain no literal boolean decision predicate"


def p22_d15_computed_controls():
    _, _, _, _, _, _, no_record, holds = tiny_arrow_rows()
    _, short, long = fresh_countlaw_pair()
    m = PM.URM()
    v = m.chain_coupling_venue(11)
    cls = [m.reachable_class(v, x) for x in (Fr(1, 3), Fr(1, 2), Fr(2, 3))]
    wk = m.writing_kernel_verdict("C8", Fr(1))
    assert no_record["chi_whole"] == 0 < holds["chi_whole"]
    assert short["k"] == 1 and long["k"] == 0
    assert [sum((r["exponential"], r["critical"], r["divergent"])) for r in cls] == [1, 1, 1]
    assert wk["critical"] and wk["leak_det"] != 0
    return "arrow zero/positive, count 1/0, three class branches, critical/leak computed"


def p23_builder_artifacts():
    expected = [
        ROOT / "model/arrow.py", ROOT / "model/countlaw.py", ROOT / "model/classes.py",
        ROOT / "model/writing.py", ROOT / "model/checks_arrow.py",
        ROOT / "model/checks_countlaw.py", ROOT / "model/checks_classes.py",
        ROOT / "model/checks_writing.py", ROOT / "model/checks_d25.py",
    ]
    assert all(p.is_file() and p.stat().st_size > 0 for p in expected)
    b = (ROOT / "LANE_T54_VERIFY/B_CODEX/VERIFY_B.OUT.txt").read_text()
    bd = (ROOT / "LANE_T54_VERIFY/B_CODEX/D24_AUDIT.txt").read_text()
    assert "VERIFIER B: 26 OK, 0 REFUTE, 26 required predicates" in b
    assert "classes=CONFIRMED; writing=CONFIRMED" in b
    assert "CORRECTION LOG" in bd and "25 OK / 1 REFUTE" in bd
    return "nine builder files present; verifier B final and initial-refutation custody present"


def p24_verifier_a_artifacts():
    out = ROOT / "LANE_T54_VERIFY/A_CODEX/verifier_t54_a.OUT.txt"
    audit = ROOT / "LANE_T54_VERIFY/A_CODEX/D24_AUDIT.txt"
    assert out.is_file(), "verifier A final transcript not present"
    text = out.read_text()
    assert "FINAL REPAIRED OVERALL: CONFIRMED" in text
    assert "FINAL VERDICT ARROW: CONFIRMED" in text
    assert "FINAL VERDICT COUNTLAW: CONFIRMED" in text
    assert "FINAL VERDICT D25: CONFIRMED" in text
    at = audit.read_text() if audit.is_file() else ""
    assert "Default verdict: REFUTED" in at
    assert "INITIAL BUILDER SNAPSHOT" in at and "initial builder overall: REFUTED" in at
    assert "math domain error" in at and "D25 integration landed" in at
    return "verifier A final confirms arrow/countlaw/D25 and preserves default-refuted audit"


def p25_proof_semantic_cutover():
    text = (ROOT / "PROOF_V002.md").read_text()
    stale = (
        "URM.surface is called nowhere in the repository outside the three validator gates",
        "C-86 and C-87 appear nowhere in model/ at all",
        "C-86 appears nowhere in model/",
        "it is bypassable, and is bypassed inside the proof's own validator",
        "`python3 model/validate_geometry.py` runs 31 checks with 0 fail",
        "geometry 36 PASS / 0 FAIL and, chained, project 14 PASS / 0 FAIL",
        "the count law, the emergence material and ten of the roles statement's seventeen "
        "load-bearing rows are outside it",
        "the two PROVED rows named in \u00a70.1",
    )
    found = [phrase for phrase in stale if phrase in text]
    assert found == [], found
    for fn in ("ProjectModel.arrow_threshold", "ProjectModel.census",
               "ProjectModel.reachable_class", "ProjectModel.writing_kernel_verdict"):
        assert f"`{fn}`" in text, fn
    assert "not a tamper-proof provenance certificate" in text
    assert "model-integrity refusal, not an authentication or security boundary" in text
    assert "33 PASS" in text and "24 PASS" in text
    return "pre-T54 absence/bypass/count claims removed; four families and D25 residual explicit"


for name, fn in (
    ("I01 THREE-DOOR CONTRACT", p01_contract),
    ("I02 INTEGRATED API SURFACE", p02_api_surface),
    ("I03 URM.SURFACE REFUSAL", p03_surface_refusal),
    ("I04 URM.SURFACE POSITIVE", p04_surface_acceptance),
    ("I05 URM.CORNER REFUSAL/POSITIVE", p05_corner_refusal_and_positive),
    ("I06 NO DIRECT RECORDSURFACE", p06_d25_model_scan),
    ("I07 COUNTLAW RAW BYPASS", p07_countlaw_raw_bypass),
    ("I08 COUNTLAW MUTATION BYPASS", p08_countlaw_mutation_bypass),
    ("I09 COUNTLAW OFF-SEALED", p09_countlaw_offsealed),
    ("I10 COUNTLAW WRAPPER", p10_countlaw_wrapper_fidelity),
    ("I11 ARROW OFF-SEALED", p11_arrow_observation_offsealed),
    ("I12 ARROW DELEGATES", p12_arrow_delegate_wiring),
    ("I13 CLASSES ENTRY REFUSAL", p13_classes_entry_refusal),
    ("I14 CLASSES OFF-SEALED", p14_classes_offsealed),
    ("I15 CLASSES DELEGATES", p15_classes_constructor_delegates),
    ("I16 WRITING KERNEL DOMAIN", p16_writing_kernel_domain),
    ("I17 WRITING WORLD DELEGATES", p17_writing_offsealed),
    ("I18 WRITING CORNER/GAP", p18_writing_corner_and_gap),
    ("I19 VALIDATOR TOPOLOGY", p19_validator_topology),
    ("I20 UMBRELLA CAPTURE", p20_umbrella_capture),
    ("I21 D-8 STATIC", p21_d8_static_scan),
    ("I22 D-15 COMPUTED CONTROLS", p22_d15_computed_controls),
    ("I23 BUILDER/VERIFIER B ARTIFACTS", p23_builder_artifacts),
    ("I24 VERIFIER A ARTIFACTS", p24_verifier_a_artifacts),
    ("I25 PROOF SEMANTIC CUTOVER", p25_proof_semantic_cutover),
):
    predicate(name, fn)

print("T-54/T-55 FINAL INTEGRATION VERIFIER -- DEFAULT REFUTED")
print("=" * 96)
for name, ok, detail in results:
    print(f"{'PASS' if ok else 'REFUTE':7s} {name} :: {detail}")
print("=" * 96)
npass = sum(ok for _, ok, _ in results)
nref = len(results) - npass
confirmed = len(results) == REQUIRED and nref == 0
print(f"INTEGRATION PREDICATES: {npass} PASS, {nref} REFUTE, {REQUIRED} required")
print(f"FINAL INTEGRATION VERDICT: {'CONFIRMED' if confirmed else DEFAULT_VERDICT}")
raise SystemExit(0 if confirmed else 1)
