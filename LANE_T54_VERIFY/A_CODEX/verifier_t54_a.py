#!/usr/bin/env python3
"""Independent T-54 verifier A: arrow, countlaw, and d25.

Default verdict is REFUTED.  A family becomes CONFIRMED only when every mandatory
test passes, including every sealed-input manifest and every expected refusal.
No builder, model-integration, register, ledger, manifest, or shared file is written.
"""

from __future__ import annotations

import ast
import contextlib
import hashlib
import io
import inspect
import math
import os
from pathlib import Path
import re
import subprocess
import sys
from dataclasses import dataclass
import types
from typing import Callable

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "model"
BASE_COMMIT = "d3b181d41719dd4e605c43cd26c2a261ba5bbf5e"
sys.path.insert(0, str(MODEL))

import checks_d25 as D25  # noqa: E402
import grounded as G  # noqa: E402
from project_model import PROVENANCE, ProjectModel, RecordSurface, URM  # noqa: E402
from record_model import Environment, RecordModel  # noqa: E402


def git_blob(rel: str) -> str:
    """Read the exact pre-repair builder snapshot, pinned by commit, without checkout."""
    cp = subprocess.run(
        ["git", "show", f"{BASE_COMMIT}:{rel}"], cwd=ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if cp.returncode:
        raise RuntimeError(f"cannot read {BASE_COMMIT}:{rel}: {cp.stderr.strip()}")
    return cp.stdout


BASE_SOURCES = {
    rel: git_blob(rel) for rel in (
        "model/arrow.py", "model/checks_arrow.py", "INTEGRATION_arrow.md",
        "model/countlaw.py", "model/checks_countlaw.py", "INTEGRATION_countlaw.md",
    )
}


def source_module(name: str, rel: str) -> types.ModuleType:
    module = types.ModuleType(name)
    module.__file__ = str(ROOT / rel)
    sys.modules[name] = module
    exec(compile(BASE_SOURCES[rel], f"{BASE_COMMIT}:{rel}", "exec"), module.__dict__)
    return module


# The first verifier pass found the defects before concurrent repairs landed.  All arrow and
# countlaw judgments below are therefore tied to the immutable pre-repair builder snapshot.
AR = source_module("t54_arrow_builder_snapshot", "model/arrow.py")
CL = source_module("t54_countlaw_builder_snapshot", "model/countlaw.py")


@dataclass(frozen=True)
class Result:
    family: str
    name: str
    passed: bool
    detail: str


RESULTS: list[Result] = []


def add(family: str, name: str, passed: bool, detail: str = "") -> None:
    RESULTS.append(Result(family, name, bool(passed), str(detail).replace("\n", " | ")))


def close(a: float, b: float, *, rtol: float = 1e-9, atol: float = 0.0) -> bool:
    return bool(abs(float(a) - float(b)) <= atol + rtol * abs(float(b)))


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_manifest(name: str) -> tuple[bool, int, str]:
    manifest = ROOT / name
    if not manifest.is_file():
        return False, 0, f"missing {name}"
    checked = 0
    bad: list[str] = []
    for raw in manifest.read_text(encoding="utf-8").splitlines():
        if not raw.strip():
            continue
        parts = raw.split(None, 1)
        if len(parts) != 2:
            bad.append(f"malformed:{raw}")
            continue
        want, rel = parts
        target = ROOT / rel.strip()
        if not target.is_file():
            bad.append(f"missing:{rel.strip()}")
            continue
        checked += 1
        got = digest(target)
        if got != want:
            bad.append(f"hash:{rel.strip()}")
    return bool(checked and not bad), checked, "ok" if not bad else ",".join(bad)


MANIFESTS = {
    "arrow": ["LANE_F1_ARROW.sha256", "LANE_PF2_DYNAMICAL.sha256", "LANE_T9_AUDIT.sha256"],
    "countlaw": ["LANE_T47_A_WIDTH.sha256", "LANE_T47_B_STAIRCASE.sha256",
                 "LANE_T47_D_REGISTER.sha256", "LANE_T31_ASYMMETRY.sha256"],
    "d25": ["LANE_T41_EXTERNAL.sha256"],
}


for family, manifests in MANIFESTS.items():
    for manifest in manifests:
        ok, count, detail = verify_manifest(manifest)
        add(family, f"sealed input {manifest}", ok, f"{count} files; {detail}")


def run_live_suite(script: str, summary: str, timeout: int = 420) -> tuple[subprocess.CompletedProcess[str], int, int]:
    cp = subprocess.run(
        [sys.executable, "-B", str(MODEL / script)],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    m = re.search(summary, cp.stdout)
    if not m:
        return cp, -1, -1
    return cp, int(m.group(1)), int(m.group(2))


# Builder check files are executed with __name__ == '__main__' against their exact pinned
# modules, with stdout captured.  This is the standalone path without materialising or editing
# a shared-tree copy.  The populated module caches are then reused by the independent probes.
def run_pinned_suite(check_rel: str, imported_name: str, imported_module: types.ModuleType,
                     summary: str) -> tuple[int, int, int]:
    old = sys.modules.get(imported_name)
    sys.modules[imported_name] = imported_module
    namespace = {"__name__": "__main__", "__file__": str(ROOT / check_rel)}
    output = io.StringIO()
    code = 0
    try:
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
            try:
                exec(compile(BASE_SOURCES[check_rel], f"{BASE_COMMIT}:{check_rel}", "exec"),
                     namespace)
            except SystemExit as exc:
                code = int(exc.code or 0)
    finally:
        if old is None:
            del sys.modules[imported_name]
        else:
            sys.modules[imported_name] = old
    m = re.search(summary, output.getvalue())
    return code, int(m.group(1)) if m else -1, int(m.group(2)) if m else -1


try:
    code, n_pass, n_fail = run_pinned_suite(
        "model/checks_arrow.py", "arrow", AR, r"ARROW:\s+(\d+) PASS,\s+(\d+) FAIL")
    add("arrow", "standalone checks_arrow.py at pinned builder snapshot",
        code == 0 and (n_pass, n_fail) == (26, 0),
        f"exit={code}; {n_pass} PASS {n_fail} FAIL")
except Exception as exc:  # default-refuted on any harness failure
    add("arrow", "standalone checks_arrow.py at pinned builder snapshot", False,
        f"{type(exc).__name__}: {exc}")

try:
    code, n_pass, n_fail = run_pinned_suite(
        "model/checks_countlaw.py", "countlaw", CL,
        r"COUNTLAW:\s+(\d+) PASS,\s+(\d+) FAIL")
    add("countlaw", "standalone checks_countlaw.py at pinned builder snapshot",
        code == 0 and (n_pass, n_fail) == (39, 0),
        f"exit={code}; {n_pass} PASS {n_fail} FAIL")
except Exception as exc:
    add("countlaw", "standalone checks_countlaw.py at pinned builder snapshot", False,
        f"{type(exc).__name__}: {exc}")

try:
    cp, n_pass, n_fail = run_live_suite("checks_d25.py", r"D25:\s+(\d+) PASS,\s+(\d+) FAIL")
    no_offenders = "OFFENDERS" not in cp.stdout
    add("d25", "integrated standalone checks_d25.py",
        cp.returncode == 0 and (n_pass, n_fail) == (10, 0) and no_offenders,
        f"exit={cp.returncode}; {n_pass} PASS {n_fail} FAIL; zero offenders={no_offenders}")
except Exception as exc:
    add("d25", "integrated standalone checks_d25.py", False,
        f"{type(exc).__name__}: {exc}")


def must(pattern: str, text: str, flags: int = 0) -> re.Match[str]:
    m = re.search(pattern, text, flags)
    if not m:
        raise AssertionError(f"sealed parse failed: {pattern}")
    return m


def explicit_refusal(call: Callable[[], object], family_word: str, required: str) -> tuple[bool, str]:
    try:
        value = call()
    except Exception as exc:
        msg = str(exc)
        intentional = (isinstance(exc, ValueError) and "REFUSES" in msg
                       and family_word in msg.upper() and required.lower() in msg.lower())
        return intentional, f"{type(exc).__name__}: {msg}"
    rendered = repr(value)
    return False, f"ACCEPTED: {rendered[:500]}{'...' if len(rendered) > 500 else ''}"


# =====================================================================================
# ARROW: parse the sealed outputs directly, then call the API fresh.
# =====================================================================================
try:
    f1 = (ROOT / "LANE_F1_ARROW/f1_arrow.txt").read_text(encoding="utf-8")
    f1b = (ROOT / "LANE_F1_ARROW/f1b_invariance.txt").read_text(encoding="utf-8")
    pf2 = (ROOT / "LANE_PF2_DYNAMICAL/pf2_history.txt").read_text(encoding="utf-8")

    sealed_chi = float(must(r"Zbar\s+\(logical\)\s+2\s+32\.000\s+([0-9.]+)\s+([0-9.]+)", f1).group(1))
    w1 = must(r"^\s*1\s+24\s+([0-9.]+)\s+\S+\s*$", f1, re.MULTILINE)
    w2 = must(r"^\s*2\s+252\s+([0-9.]+)\s+.+$", f1, re.MULTILINE)
    sealed_w1, sealed_w2 = float(w1.group(1)), float(w2.group(1))
    inv_i = float(must(r"max \|I\(S:B\) change\|[^=]*=\s*([0-9.eE+-]+)", f1b).group(1))
    inv_fixed = float(must(r"FIXED label Zbar\| change\s*=\s*([0-9.eE+-]+)", f1b).group(1))
    ze = must(r"Z_e\s+\(single site\)\s+1\s+([0-9.]+)\s+([0-9.]+)", f1b)
    sealed_ze_i, sealed_ze_chi = float(ze.group(1)), float(ze.group(2))
    hist_rows = re.findall(
        r"^\s+([+-]?\d+\.\d+)\s+[-+0-9.]+\s+([0-9.]+)\s+[-+0-9.]+\s*$",
        pf2, re.MULTILINE)
    sealed_hist = {float(t): float(v) for t, v in hist_rows}
    frag = must(r"Zbar \(logical\)\s+2\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)", pf2)
    sealed_red = [float(frag.group(i)) for i in range(1, 5)]

    C = AR.carrier()
    th = AR.arrow_threshold()
    ok = (th[1]["n_swept"] == 24 and th[2]["n_swept"] == 252
          and close(th[1]["max_chi"], sealed_w1, atol=1e-10)
          and close(th[2]["max_chi"], sealed_w2, atol=5e-9)
          and close(th["logical"]["chi_Zbar"], sealed_chi, atol=5e-9)
          and close(th["logical"]["chi_Zbar"], th["logical"]["closed_form"], atol=1e-9))
    add("arrow", "sealed F-17 threshold re-derived then API-fresh", ok,
        f"w1={th[1]['max_chi']:.3e}/24; w2={th[2]['max_chi']:.8f}/252; sealed={sealed_chi:.8f}")

    led = {row["coupling"]: row for row in AR.arrow_ledger()}
    ok = (close(led["Ze"]["I_SB"], sealed_ze_i, atol=5e-9)
          and close(led["Ze"]["chi_record"], sealed_ze_chi, atol=1e-10)
          and close(led["Zbar"]["I_SB"], sealed_chi, atol=5e-9)
          and led["identity"]["I_SB"] < 1e-10)
    add("arrow", "sealed F-18 ledger re-derived then API-fresh", ok,
        f"Ze I={led['Ze']['I_SB']:.8f} chi={led['Ze']['chi_record']:.3e}; identity I={led['identity']['I_SB']:.3e}")

    inv = AR.arrow_invariance()
    ok = (f"{inv['mutual_worst']:.3e}" == f"{inv_i:.3e}"
          and f"{inv['fixed_label_worst']:.3e}" == f"{inv_fixed:.3e}"
          and inv["covariance_worst"] < 1e-8)
    add("arrow", "sealed F-19 invariance re-derived then API-fresh", ok,
        f"mutual={inv['mutual_worst']:.3e}; fixed={inv['fixed_label_worst']:.3e}; covariance={inv['covariance_worst']:.3e}")

    hist, _ = AR.arrow_history(sorted(sealed_hist))
    fresh_hist = {float(row["t"]): float(row["chi"]) for row in hist}
    ok = all(close(fresh_hist[t], v, atol=5e-9) for t, v in sealed_hist.items())
    add("arrow", "sealed PF-2 history re-derived then API-fresh", ok,
        " ".join(f"t={t:g}:{fresh_hist[t]:.8f}" for t in sorted(fresh_hist)))

    red = AR.arrow_redundancy()
    got_red = [red["whole"], *red["fragments"]]
    ok = (close(got_red[0], sealed_red[0], atol=5e-9)
          and all(close(a, b, atol=5e-7) for a, b in zip(got_red[1:], sealed_red[1:])))
    add("arrow", "sealed F-21 redundancy re-derived then API-fresh", ok,
        f"whole={got_red[0]:.8f}; fragments={[f'{x:.6f}' for x in got_red[1:]]}")

    # Independent off-grid closed form: this verifier owns the entropy calculation.
    lam = 0.37
    eb = np.array([0.0, 0.7, 1.3, 2.1])
    bb = np.array([1.0, 0.3, -0.2, -0.9])
    p_plus = np.exp(-2.0 * (eb + lam * bb))
    p_minus = np.exp(-2.0 * (eb - lam * bb))
    z_plus, z_minus = p_plus.sum(), p_minus.sum()
    q_plus, q_minus = z_plus / (z_plus + z_minus), z_minus / (z_plus + z_minus)
    p_plus, p_minus = p_plus / z_plus, p_minus / z_minus

    def entropy(p: np.ndarray) -> float:
        p = p[p > 1e-13]
        return float(-(p * np.log2(p)).sum())

    independent = entropy(q_plus * p_plus + q_minus * p_minus) \
        - q_plus * entropy(p_plus) - q_minus * entropy(p_minus)
    r = AR.mean_force_state(C["Zbar"], lam=lam)
    api = AR.chi(r, C["Zbar"], C["nS"], AR.NB4)
    add("arrow", "off-grid lambda=0.37 independent API-fidelity",
        close(api, independent, atol=1e-10) and abs(api - sealed_chi) > 1e-3,
        f"API={api:.12f}; independent={independent:.12f}; not sealed={abs(api-sealed_chi):.3e}")

    ok_ref, detail = explicit_refusal(
        lambda: AR.score_bath_observation(Environment(nq=1, energies=(0.9,), beta=2.0),
                                          np.array([[1, 0], [0, -1]], complex)),
        "ARROW", "provenance")
    add("arrow", "world bath without provenance explicitly refused", ok_ref, detail)

    z2 = np.array([[1, 0], [0, -1]], dtype=complex)
    custom = RecordModel(np.zeros((2, 2), dtype=complex))
    custom_env = Environment(nq=1, energies=(0.9,), beta=2.0)
    scored = AR.score_bath_observation(custom_env, z2, record=z2, model=custom,
                                       tier="corner", provenance="DEF-A")
    add("arrow", "custom-model positive control with explicit record",
        scored["holds_record_bits"] and close(scored["chi_whole"], scored["I_SB"], atol=1e-10),
        f"I={scored['I_SB']:.12f}; chi={scored['chi_whole']:.12f}")

    ok_ref, detail = explicit_refusal(
        lambda: AR.score_bath_observation(custom_env, z2, model=custom,
                                          tier="corner", provenance="DEF-A"),
        "ARROW", "record")
    add("arrow", "custom model with omitted record must be intentionally refused", ok_ref, detail)

    integ_arrow = BASE_SOURCES["INTEGRATION_arrow.md"]
    sig = must(r"def arrow_observation\((.*?)\):", integ_arrow, re.DOTALL).group(1)
    add("arrow", "integration method exposes and forwards custom model",
        "model" in sig and "model=model" in integ_arrow,
        "signature=" + " ".join(sig.split()))
except Exception as exc:
    add("arrow", "independent arrow audit completed", False, f"{type(exc).__name__}: {exc}")


# =====================================================================================
# COUNTLAW: parse sealed outputs directly, then call the API fresh.
# =====================================================================================
try:
    atext = (ROOT / "LANE_T47_A_WIDTH/t47_a_width.txt").read_text(encoding="utf-8")
    btext = (ROOT / "LANE_T47_B_STAIRCASE/t47b_staircase.txt").read_text(encoding="utf-8")
    rate_rows = re.findall(
        r"^\s+(0\.000|0\.050|0\.158|0\.300)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)",
        atext, re.MULTILINE)
    sealed_rates = {float(de): float(rate) for de, rate, _ in rate_rows}
    sealed_delta = float(must(r"EXACT derived threshold[^:]*:\s+([0-9.eE+-]+) J", atext).group(1))
    drop_rows = re.findall(r"^\s+r([0-5])\s+[0-9.]+\s+([0-9.eE+-]+)\s+", atext, re.MULTILINE)
    sealed_drops = [float(v) for _, v in sorted(drop_rows, key=lambda pair: int(pair[0]))]
    stair_rows = re.findall(
        r"^\s+(~9 h|0\.1 y|1 y|10 y|100 y|1 ky|10 ky)\s+(\d+)\s+(\d+)\s+ok$",
        atext, re.MULTILINE)
    sealed_k = [int(k1) for _, k1, _ in stair_rows]

    EV = CL.EV

    def corner(name: str, B_eV: float, dE_eV: float, T: float = 300.0, f0: float = 1e9):
        return URM.surface(name, "thermal two-well verifier corner", dE_eV * EV,
                           (B_eV - dE_eV) * EV, T, f0,
                           tier="corner", provenance="DEF-A")

    rate_api = {de: CL.record_rate(corner(f"rate-{de}", 1.2, de)) for de in sealed_rates}
    rate_formula = {
        de: 1e9 * math.exp(-((1.2 - de) * EV) / (G.KB * 300.0))
            + 1e9 * math.exp(-(1.2 * EV) / (G.KB * 300.0))
        for de in sealed_rates
    }
    ok = all(close(rate_api[de], sealed_rates[de], rtol=1e-9)
             and close(rate_api[de], rate_formula[de], rtol=1e-12)
             for de in sealed_rates)
    add("countlaw", "sealed instrument rates re-derived then API-fresh", ok,
        "; ".join(f"dE={de:.3f}:{rate_api[de]:.12e}" for de in sorted(rate_api)))

    d_api = CL.delta_pop(1.2 * EV, 300.0, 1e9, 10.0 * CL.YEAR)
    y = (1.2 * EV) / (G.KB * 300.0) - math.log(1e9 * 10.0 * CL.YEAR)
    d_independent = G.KB * 300.0 * math.log(math.expm1(y))
    add("countlaw", "sealed exact width re-derived then API-fresh",
        close(d_api, sealed_delta, rtol=1e-12) and close(d_api, d_independent, rtol=1e-14),
        f"API={d_api:.15e}; sealed={sealed_delta:.15e}; independent={d_independent:.15e}")

    ens6 = [corner(f"r{i}", 1.2, de) for i, de in enumerate((0, .05, .10, .15, .20, .25))]
    c = CL.census(ens6, CL.YEAR)
    got_drops = [row["t_star"] for row in c["schedule"]]
    ok = len(sealed_drops) == 6 and all(close(a, b, rtol=1e-6)
                                       for a, b in zip(got_drops, sealed_drops))
    add("countlaw", "sealed six-record drop schedule re-derived then API-fresh", ok,
        " ".join(f"r{i}={v:.6e}" for i, v in enumerate(got_drops)))

    t_grid = (3.156e4, 3.156e6, 3.156e7, 3.156e8, 3.156e9, 3.156e10, 3.156e11)
    got_k = [CL.census(ens6, t)["k"] for t in t_grid]
    add("countlaw", "sealed staircase re-derived then API-fresh",
        got_k == sealed_k == [6, 6, 5, 4, 2, 1, 0], f"k={got_k}")

    t31_rows = re.findall(
        r"^\s+(0\.00|0\.05|0\.10|0\.16)\s+([0-9.]+)\s+([0-9.]+)\s+"
        r"([0-9.eE+-]+)\s+([0-9.eE+-]+)", btext, re.MULTILINE)
    t31_ok = len(t31_rows) == 4
    t31_detail: list[str] = []
    for eps_s, de1_s, de2_s, tau1_s, tau2_s in t31_rows:
        eps = float(eps_s)
        got = CL.t31_basin(eps)
        t31_ok &= (abs(got["dE"][0] - float(de1_s)) < 5e-4
                   and abs(got["dE"][1] - float(de2_s)) < 5e-4
                   and close(got["taus"][0], float(tau1_s), rtol=5e-4)
                   and close(got["taus"][1], float(tau2_s), rtol=5e-4))
        t31_detail.append(f"eps={eps:g}:dE={got['dE'][0]:.4f}/{got['dE'][1]:.4f},"
                          f"tau={got['taus'][0]:.4e}/{got['taus'][1]:.4e},v2={got['v2_exact']}")
    add("countlaw", "sealed T-31 control re-derived then API-fresh", t31_ok, "; ".join(t31_detail))

    sealed_departure = float(must(r"magnetisation-weighted departure = .*? = ([0-9.]+) \(of N", btext).group(1))
    kTg = G.KB * 350.0
    xs = [0.2, 0.5, 0.8, 1.2, 1.8, 2.5, 3.2, 4.0, 5.0, 6.5]
    grains = [corner(f"g{i}", 1.1, x * kTg / EV, 350.0, 1e9) for i, x in enumerate(xs)]
    cg = CL.census(grains, 1e6)
    add("countlaw", "sealed departure term re-derived then API-fresh",
        cg["k"] == 4 and abs(cg["departure"] - sealed_departure) < 5e-4,
        f"k={cg['k']}; departure={cg['departure']:.4f}; sealed={sealed_departure:.4f}")

    rejected: list[str] = []
    for keyword in ("width", "tol", "delta", "cluster_width", "margin"):
        try:
            CL.census(ens6, 1.0, **{keyword: 1e-2})
        except TypeError:
            rejected.append(keyword)
    add("countlaw", "C-76 chosen-width API is unreachable",
        rejected == ["width", "tol", "delta", "cluster_width", "margin"]
        and list(inspect.signature(CL.census).parameters) == ["surfaces", "t_m"],
        f"rejected={rejected}; signature={inspect.signature(CL.census)}")

    # Off-grid 500 K check uses verifier-owned bisection against the record-mode API.
    B = 1.2 * EV
    target_t = 1.0
    lo, hi = 0.0, 0.999 * B
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        rmid = CL.record_rate(URM.surface(
            "500K probe", "thermal two-well verifier corner", mid, B - mid, 500.0, 1e9,
            tier="corner", provenance="DEF-A"))
        if rmid <= 1.0 / target_t:
            lo = mid
        else:
            hi = mid
        if hi - lo <= np.spacing(hi):
            break
    bisected = 0.5 * (lo + hi)
    formula_500 = CL.delta_pop(B, 500.0, 1e9, target_t)
    add("countlaw", "off-grid 500 K API-fidelity bisection",
        close(bisected, formula_500, rtol=1e-11),
        f"bisected={bisected/EV:.12f} eV; formula={formula_500/EV:.12f} eV")

    valid_world = URM.surface("NAND floating gate", "trapped charge", .05 * EV,
                              1.0 * EV, 358.0, 1e9)
    valid_corner = corner("valid-corner", 1.2, 0.05)
    valid = CL.census([valid_world, valid_corner], 1e3)
    add("countlaw", "valid world and DEF-A corner entries accepted",
        bool(valid["schedule"][0]["provenance"])
        and valid["schedule"][1]["provenance"] == "DEF-A",
        f"world={valid['schedule'][0]['provenance']!r}; corner={valid['schedule'][1]['provenance']!r}")

    raw_world = RecordSurface("raw world", "thermal", .05 * EV, 1.0 * EV, 300.0, 1e9)
    ok_ref, detail = explicit_refusal(lambda: CL.census([raw_world], 1e3),
                                      "COUNTLAW", "provenance")
    add("countlaw", "raw unprovenanced world surface must be refused at census", ok_ref, detail)

    raw_corner = RecordSurface("raw corner", "thermal", .05 * EV, 1.0 * EV, 300.0, 1e9)
    raw_corner.tier = "corner"
    raw_corner.provenance = "not DEF-A"
    ok_ref, detail = explicit_refusal(lambda: CL.census([raw_corner], 1e3),
                                      "COUNTLAW", "DEF-A")
    add("countlaw", "corner census must require exact DEF-A self-declaration", ok_ref, detail)
except Exception as exc:
    add("countlaw", "independent countlaw audit completed", False, f"{type(exc).__name__}: {exc}")


# =====================================================================================
# D25: independent AST scan, exact integration-footprint proof, and fresh anchor calculation.
# =====================================================================================
def scan_source(source: str, relpath: str) -> list[tuple[str, int, str]]:
    tree = ast.parse(source)
    out: list[tuple[str, int, str]] = []
    for node in ast.walk(tree):
        kind = None
        if isinstance(node, ast.ImportFrom) and any(a.name == "RecordSurface" for a in node.names):
            kind = "import"
        elif isinstance(node, ast.Name) and node.id == "RecordSurface":
            kind = "reference"
        elif isinstance(node, ast.Attribute) and node.attr == "RecordSurface":
            kind = "attribute"
        if kind:
            out.append((relpath, int(getattr(node, "lineno", 0)), kind))
    return sorted(out, key=lambda row: (row[0], row[1], row[2]))


def scan_tree(root: Path) -> list[tuple[str, int, str]]:
    """Verifier-owned D-25 scan; it does not call the builder's tree walker."""
    out: list[tuple[str, int, str]] = []
    for path in sorted(root.rglob("*.py")):
        if path.name == "project_model.py" or "__pycache__" in path.parts:
            continue
        rel = str(path.relative_to(root))
        try:
            out.extend(scan_source(path.read_text(encoding="utf-8"), rel))
        except SyntaxError as exc:
            out.append((rel, int(exc.lineno or 0), "unparseable"))
    return sorted(out, key=lambda row: (row[0], row[1], row[2]))


try:
    validate_path = MODEL / "validate_project.py"
    live_source = validate_path.read_text(encoding="utf-8")
    independent_off = scan_tree(MODEL)
    builder_off = sorted(D25.scan_direct_constructions(root=str(MODEL)),
                         key=lambda row: (row[0], row[1], row[2]))
    add("d25", "integrated model tree has zero D-25 construction offenders",
        independent_off == builder_off == [],
        f"independent={independent_off}; builder={builder_off}")

    # Rebuild validate_project.py from the immutable pre-integration blob using only the
    # observed D-25/URM routing changes. Exact equality rules out unrelated edits there.
    base_validate = git_blob("model/validate_project.py")
    live_edits = [
        ("from project_model import RecordSurface, ProjectModel",
         "from project_model import URM"),
        ("M = ProjectModel(); n_pass = 0; n_fail = 0",
         "M = URM(); n_pass = 0; n_fail = 0"),
        ('grain = RecordSurface("CoCrPt grain", "magnetic anisotropy", 3.0 * G.KB * 300, 2.0e5 * 1.26e-24, 300.0, 1e9)',
         'grain = M.surface("CoCrPt grain", "magnetic anisotropy", 3.0 * G.KB * 300,\n'
         '                  2.0e5 * 1.26e-24, 300.0, 1e9)'),
        ("    s = RecordSurface(nm, mech, dE, Eb, T, f0)",
         "    s = M.surface(nm, mech, dE, Eb, T, f0)"),
        ('zirc = RecordSurface("Zircon U-238", "nuclear decay", 4.27e6 * eV, 4.27e6 * eV, 300.0, 1e21)\n'
         'cmb = RecordSurface("CMB photon", "free flight", 0.0, 0.0, 2.7, 0.0, thermal=False)',
         'zirc = M.surface("Zircon U-238", "nuclear decay", 4.27e6 * eV, 4.27e6 * eV,\n'
         '                 300.0, 1e21, provenance="control surface, census GR1 entry 3 (zircon U-Pb): "\n'
         '                 "decay is temperature-independent; the model must DECLINE")\n'
         'cmb = M.surface("CMB photon", "free flight", 0.0, 0.0, 2.7, 0.0, thermal=False,\n'
         '                provenance="control surface, census GR1 entry 4 (CMB photon polarisation): "\n'
         '                "no bath, free flight; the model must DECLINE")'),
        ('    r = M.corner(H)', '    r = M.corner(H, provenance="DEF-A")'),
        ('# D-25: the URM\'s public gate must REFUSE undeclared surfaces — the guard is tested, not assumed\n'
         'from project_model import URM\n'
         'try:',
         '# D-25: the URM\'s public gate must REFUSE undeclared surfaces — the guard is tested, not assumed\n'
         'try:'),
        ('    URM.surface("mystery device", "unknown", 1e-20, 1e-19, 300.0, 1e9)',
         '    M.surface("mystery device", "unknown", 1e-20, 1e-19, 300.0, 1e9)'),
        ('    URM.surface("toy torus", "stabiliser", 0.0, 0.0, 0.0, 0.0, tier="corner")',
         '    M.surface("toy torus", "stabiliser", 0.0, 0.0, 0.0, 0.0, tier="corner")'),
        ('ok = URM.surface("CoCrPt grain", "magnetic anisotropy", 3.0*G.KB*300, 2.0e5*1.26e-24, 300.0, 1e9)',
         'ok = M.surface("CoCrPt grain", "magnetic anisotropy", 3.0*G.KB*300,\n'
         '               2.0e5*1.26e-24, 300.0, 1e9)'),
        ('check("D-25 registry supplies pinned provenance", "Weller" in ok.provenance, ok.provenance[:60])\n'
         'print("=" * 78)',
         'check("D-25 registry supplies pinned provenance", "Weller" in ok.provenance, ok.provenance[:60])\n'
         '# D-25 (T-55): construction scan + stage-(2) external-anchor gates — see model/checks_d25.py\n'
         'from checks_d25 import run_d25_checks\n'
         'run_d25_checks(check)\n'
         'print("=" * 78)'),
    ]
    counts = [base_validate.count(old) for old, _ in live_edits]
    candidate = base_validate
    for old, new in live_edits:
        candidate = candidate.replace(old, new, 1)
    add("d25", "validate_project integration is exactly the minimal URM/D-25 routing patch",
        counts == [1] * len(live_edits) and candidate == live_source,
        f"old-occurrences={counts}; exact-live={candidate == live_source}; edits={len(live_edits)}")

    expected_registry = {"CoCrPt grain", "CoCrPt HDD grain", "NAND floating gate",
                         "DNA base tautomer", "Fe(II) spin crossover", "Azobenzene",
                         "Azobenzene cis/trans", "Alanine enantiomer"}
    aliases = {"CoCrPt HDD grain": "CoCrPt grain",
               "Azobenzene cis/trans": "Azobenzene"}
    alias_ok = all(name in PROVENANCE and f"alias of {target}" in PROVENANCE[name]
                   for name, target in aliases.items())
    add("d25", "display rows resolve through two explicit provenance aliases",
        expected_registry <= set(PROVENANCE) and alias_ok,
        f"resolved={sorted(expected_registry & set(PROVENANCE))}; aliases={alias_ok}")

    idoc = (ROOT / "INTEGRATION_d25.md").read_text(encoding="utf-8")
    prospective_pairs = [
        ("from project_model import RecordSurface, ProjectModel",
         "from project_model import ProjectModel, URM"),
        ('grain = RecordSurface("CoCrPt grain", "magnetic anisotropy", 3.0 * G.KB * 300, 2.0e5 * 1.26e-24, 300.0, 1e9)',
         'grain = URM.surface("CoCrPt grain", "magnetic anisotropy", 3.0 * G.KB * 300, 2.0e5 * 1.26e-24, 300.0, 1e9)'),
        ('SURFACES = [("CoCrPt HDD grain", "magnetic anisotropy", 3.0 * G.KB * 300, 60.8 * G.KB * 300, 300.0, 1e9),',
         'SURFACES = [("CoCrPt grain", "magnetic anisotropy", 3.0 * G.KB * 300, 60.8 * G.KB * 300, 300.0, 1e9),'),
        ('            ("Azobenzene cis/trans", "photoisomerisation", 0.60 * eV, 1.05 * eV, 300.0, 1e13),',
         '            ("Azobenzene", "photoisomerisation", 0.60 * eV, 1.05 * eV, 300.0, 1e13),'),
        ("    s = RecordSurface(nm, mech, dE, Eb, T, f0)",
         "    s = URM.surface(nm, mech, dE, Eb, T, f0)"),
        ('zirc = RecordSurface("Zircon U-238", "nuclear decay", 4.27e6 * eV, 4.27e6 * eV, 300.0, 1e21)',
         'zirc = URM.surface("Zircon U-238", "nuclear decay", 4.27e6 * eV, 4.27e6 * eV, 300.0, 1e21, provenance="control surface, census GR1 entry 3 (zircon U-Pb): decay is temperature-independent; the model must DECLINE")'),
        ('cmb = RecordSurface("CMB photon", "free flight", 0.0, 0.0, 2.7, 0.0, thermal=False)',
         'cmb = URM.surface("CMB photon", "free flight", 0.0, 0.0, 2.7, 0.0, thermal=False, provenance="control surface, census GR1 entry 4 (CMB photon polarisation): no bath, free flight; the model must DECLINE")'),
    ]
    doc_exact = all(f"old: {old}" in idoc and f"new: {new}" in idoc
                    for old, new in prospective_pairs)
    chain_exact = ("from checks_d25 import run_d25_checks" in idoc
                   and "run_d25_checks(check)" in idoc)
    alias_option = ("keep both display names and add alias entries" in idoc
                    and '"CoCrPt HDD grain" / "Azobenzene cis/trans"' in idoc)
    add("d25", "integration brief carries exact prospective edits, alias option, and chain hook",
        doc_exact and alias_option and chain_exact,
        f"seven old/new pairs={doc_exact}; alias-option={alias_option}; chain={chain_exact}")

    citations = (ROOT / "LANE_T41_EXTERNAL/CITATIONS.md").read_text(encoding="utf-8")
    pinned = ("1.4 days in benzene at 35 °C" in citations and "0.915 eV" in citations
              and "−50.2 J/mol/K" in citations)
    EV = 1.602176634e-19
    T = 308.15
    dE = 0.60 * EV
    dH = 0.915 * EV
    dS = -50.2
    f0_eff = (G.KB * T / D25.H_PLANCK) * math.exp(dS / D25.R_GAS)
    surface = URM.surface("Azobenzene", "photoisomerisation", dE, dH, T, f0_eff)
    tau_api = ProjectModel().lifetime(surface)
    gu = f0_eff * math.exp(-dH / (G.KB * T))
    gd = f0_eff * math.exp(-(dH + dE) / (G.KB * T))
    tau_independent = 1.0 / (gu + gd)
    half_api = tau_api * math.log(2.0)
    dec = math.log10(half_api / 1.2e5)
    midpoint_half = (1.0 / (
        f0_eff * math.exp(-(dH - dE / 2.0) / (G.KB * T))
        + f0_eff * math.exp(-((dH - dE / 2.0) + dE) / (G.KB * T)))) * math.log(2.0)
    midpoint_dec = math.log10(midpoint_half / 1.2e5)
    add("d25", "pinned azobenzene anchor computed independently through fresh API",
        pinned and close(tau_api, tau_independent, rtol=1e-9)
        and abs(dec) <= 1.0 and abs(midpoint_dec) > 1.0,
        f"t1/2 API={half_api:.3e}s; measured=1.2e5s; delta={dec:+.2f} decades; midpoint={midpoint_dec:+.2f}")

    ok_ref, detail = explicit_refusal(
        lambda: URM.surface("mystery", "unknown", 1e-20, 1e-19, 300.0, 1e9,
                            provenance="   "),
        "URM", "provenance")
    add("d25", "whitespace world provenance explicitly refused", ok_ref, detail)
    ok_ref, detail = explicit_refusal(
        lambda: URM.surface("toy", "stabiliser", 0.0, 0.0, 0.0, 0.0,
                            tier="corner", provenance="def-a"),
        "URM", "DEF-A")
    add("d25", "lowercase corner declaration explicitly refused", ok_ref, detail)
except Exception as exc:
    add("d25", "independent d25 audit completed", False, f"{type(exc).__name__}: {exc}")


# =====================================================================================
# D-8 / D-15 source scan.  No fit appears in these blocks; every condition must be an
# expression rather than a literal truth value, and each block must carry control markers.
# =====================================================================================
def static_check_scan(source: str) -> tuple[bool, str]:
    tree = ast.parse(source)
    calls = [node for node in ast.walk(tree)
             if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
             and node.func.id == "check" and len(node.args) >= 2]
    literal = [int(getattr(node, "lineno", 0)) for node in calls
               if isinstance(node.args[1], ast.Constant)
               and isinstance(node.args[1].value, bool)]
    fits = [token for token in ("polyfit", "curve_fit", "lstsq") if token in source]
    controls = ("positive control" in source.lower() or "d-15" in source.lower())
    return bool(calls and not literal and not fits and controls), \
        f"{len(calls)} check call sites; literal-bool lines={literal}; fits={fits}; controls={controls}"


baseline_check_sources = {
    "arrow": BASE_SOURCES["model/checks_arrow.py"],
    "countlaw": BASE_SOURCES["model/checks_countlaw.py"],
    "d25": (MODEL / "checks_d25.py").read_text(encoding="utf-8"),
}
for family, source in baseline_check_sources.items():
    try:
        ok, detail = static_check_scan(source)
        add(family, "D-8/D-15 static scan", ok, detail)
    except Exception as exc:
        add(family, "D-8/D-15 static scan", False, f"{type(exc).__name__}: {exc}")


# Preserve the first-pass judgment as its own immutable stage before testing the repairs.
RESULTS[:] = [Result(f"{row.family}_initial" if row.family in ("arrow", "countlaw") else row.family,
                     row.name, row.passed, row.detail) for row in RESULTS]


# =====================================================================================
# REPAIRED LIVE TREE: full standalone suites plus the exact adversarial attacks that
# refuted the builder snapshot.  These are a new default-refuted judgment, not a rewrite
# of the initial evidence above.
# =====================================================================================
import importlib  # noqa: E402

AR_LIVE = importlib.import_module("arrow")
CL_LIVE = importlib.import_module("countlaw")

for family, manifests in (("arrow", MANIFESTS["arrow"]),
                          ("countlaw", MANIFESTS["countlaw"])):
    for manifest in manifests:
        ok, count, detail = verify_manifest(manifest)
        add(family, f"sealed input {manifest} after repair", ok, f"{count} files; {detail}")

try:
    cp, n_pass, n_fail = run_live_suite("checks_arrow.py", r"ARROW:\s+(\d+) PASS,\s+(\d+) FAIL")
    add("arrow", "repaired standalone checks_arrow.py",
        cp.returncode == 0 and (n_pass, n_fail) == (27, 0),
        f"exit={cp.returncode}; {n_pass} PASS {n_fail} FAIL")
except Exception as exc:
    add("arrow", "repaired standalone checks_arrow.py", False, f"{type(exc).__name__}: {exc}")

try:
    z2 = np.diag([1.0, -1.0])
    custom = RecordModel(np.zeros((2, 2)), [])
    env2 = Environment(nq=2, energies=(0.9, 1.6), beta=2.0)
    ok_ref, detail = explicit_refusal(
        lambda: AR_LIVE.score_bath_observation(env2, z2, model=custom,
                                                tier="corner", provenance="DEF-A"),
        "ARROW", "record")
    add("arrow", "repaired custom model omitting record explicitly refuses", ok_ref, detail)
    ok_shape, shape_detail = explicit_refusal(
        lambda: AR_LIVE.score_bath_observation(env2, z2, record=AR_LIVE.carrier()["Zbar"],
                                                model=custom, tier="corner", provenance="DEF-A"),
        "ARROW", "dimension")
    add("arrow", "repaired custom model rejects mismatched record shape", ok_shape, shape_detail)
    scored = AR_LIVE.score_bath_observation(env2, z2, record=z2, model=custom,
                                            tier="corner", provenance="DEF-A")
    add("arrow", "repaired custom-model explicit-record positive control",
        scored["chi_whole"] > 0.0 and scored["I_SB"] + 1e-9 >= scored["chi_whole"],
        f"I={scored['I_SB']:.12f}; chi={scored['chi_whole']:.12f}; fragments={len(scored['fragments'])}")
    live_integration = (ROOT / "INTEGRATION_arrow.md").read_text(encoding="utf-8")
    live_sig = must(r"def arrow_observation\((.*?)\):", live_integration, re.DOTALL).group(1)
    add("arrow", "repaired integration exposes and forwards custom model",
        "model" in live_sig and "model=model" in live_integration,
        "signature=" + " ".join(live_sig.split()))
    ok_static, static_detail = static_check_scan(
        (MODEL / "checks_arrow.py").read_text(encoding="utf-8"))
    add("arrow", "repaired D-8/D-15 static scan", ok_static, static_detail)
except Exception as exc:
    add("arrow", "repaired arrow adversarial audit completed", False,
        f"{type(exc).__name__}: {exc}")

try:
    cp, n_pass, n_fail = run_live_suite("checks_countlaw.py",
                                        r"COUNTLAW:\s+(\d+) PASS,\s+(\d+) FAIL")
    add("countlaw", "repaired standalone checks_countlaw.py",
        cp.returncode == 0 and (n_pass, n_fail) == (40, 0),
        f"exit={cp.returncode}; {n_pass} PASS {n_fail} FAIL")
except Exception as exc:
    add("countlaw", "repaired standalone checks_countlaw.py", False,
        f"{type(exc).__name__}: {exc}")

try:
    ev = CL_LIVE.EV
    valid_world = URM.surface("NAND floating gate", "trapped charge", .05 * ev,
                              1.0 * ev, 358.0, 1e9)
    valid_corner = URM.surface("valid repaired corner", "thermal", .05 * ev, 1.0 * ev,
                               300.0, 1e9, tier="corner", provenance="DEF-A")
    valid = CL_LIVE.census([valid_world, valid_corner], 1e3)
    add("countlaw", "repaired valid world and corner positive control",
        len(valid["schedule"]) == 2 and valid["schedule"][0]["tier"] == "world"
        and valid["schedule"][1]["provenance"] == "DEF-A",
        f"k={valid['k']}; tiers={[row['tier'] for row in valid['schedule']]}")

    raw_world = RecordSurface("raw repaired world", "thermal", .05 * ev, 1.0 * ev,
                              300.0, 1e9)
    ok_ref, detail = explicit_refusal(lambda: CL_LIVE.census([raw_world], 1e3),
                                      "COUNTLAW", "provenance")
    add("countlaw", "repaired census refuses raw unprovenanced world", ok_ref, detail)

    false_corner = RecordSurface("false repaired corner", "thermal", .05 * ev, 1.0 * ev,
                                 300.0, 1e9)
    false_corner.tier = "corner"
    false_corner.provenance = "not DEF-A"
    ok_ref, detail = explicit_refusal(lambda: CL_LIVE.census([false_corner], 1e3),
                                      "COUNTLAW", "DEF-A")
    add("countlaw", "repaired census refuses false corner declaration", ok_ref, detail)

    mutated = URM.surface("NAND floating gate", "trapped charge", .05 * ev,
                          1.0 * ev, 358.0, 1e9)
    mutated.provenance = None
    direct_results = []
    direct_details = []
    for label, call in (("record_rate", lambda: CL_LIVE.record_rate(mutated)),
                        ("drop_time", lambda: CL_LIVE.drop_time(mutated)),
                        ("drop_time_formula", lambda: CL_LIVE.drop_time_formula(mutated))):
        ok_direct, direct_detail = explicit_refusal(call, "COUNTLAW", "provenance")
        direct_results.append(ok_direct)
        direct_details.append(f"{label}={direct_detail}")
    add("countlaw", "repaired every public surface consumer refuses post-entry mutation",
        all(direct_results), "; ".join(direct_details))

    # Fresh off-grid definition check after the guard repair.
    B = 1.17 * ev
    t_m = 2.5
    f0_probe = 3e9
    d_api = CL_LIVE.delta_pop(B, 425.0, f0_probe, t_m)
    y = B / (G.KB * 425.0) - math.log(f0_probe * t_m)
    d_ind = G.KB * 425.0 * math.log(math.expm1(y))
    add("countlaw", "repaired off-grid width still equals independent definition",
        close(d_api, d_ind, rtol=1e-14),
        f"API={d_api/ev:.12f} eV; independent={d_ind/ev:.12f} eV")
    ok_static, static_detail = static_check_scan(
        (MODEL / "checks_countlaw.py").read_text(encoding="utf-8"))
    add("countlaw", "repaired D-8/D-15 static scan", ok_static, static_detail)
except Exception as exc:
    add("countlaw", "repaired countlaw adversarial audit completed", False,
        f"{type(exc).__name__}: {exc}")


PINNED_FILES = [
    "model/arrow.py", "model/checks_arrow.py", "INTEGRATION_arrow.md",
    "model/countlaw.py", "model/checks_countlaw.py", "INTEGRATION_countlaw.md",
]
LIVE_FILES = PINNED_FILES + [
    "model/checks_d25.py", "INTEGRATION_d25.md", "model/project_model.py",
    "model/validate_project.py",
]


print("T-54 VERIFIER A (CODEX) -- ARROW / COUNTLAW / D25")
print("DEFAULT: REFUTED.  CONFIRMED requires EVERY mandatory test, sealed input, and expected refusal.")
print("=" * 100)
print(f"INITIAL BUILDER SNAPSHOT {BASE_COMMIT} (sha256 of git blobs)")
for rel in PINNED_FILES:
    print(f"  {hashlib.sha256(BASE_SOURCES[rel].encode('utf-8')).hexdigest()}  {rel}")
print("LIVE REPAIRED SNAPSHOT (sha256)")
for rel in LIVE_FILES:
    print(f"  {digest(ROOT / rel)}  {rel}")

initial_overall = True
for family in ("arrow_initial", "countlaw_initial"):
    rows = [row for row in RESULTS if row.family == family]
    print("-" * 100)
    print(f"INITIAL FAMILY {family.removesuffix('_initial').upper()}")
    for row in rows:
        print(f"  {'PASS' if row.passed else 'FAIL'}  {row.name} :: {row.detail}")
    confirmed = bool(rows) and all(row.passed for row in rows)
    initial_overall &= confirmed
    print(f"  INITIAL VERDICT {family.removesuffix('_initial').upper()}: "
          f"{'CONFIRMED' if confirmed else 'REFUTED'} "
          f"({sum(row.passed for row in rows)}/{len(rows)} mandatory tests pass)")

overall = True
for family in ("arrow", "countlaw", "d25"):
    rows = [row for row in RESULTS if row.family == family]
    print("-" * 100)
    print(f"FINAL FAMILY {family.upper()}")
    for row in rows:
        print(f"  {'PASS' if row.passed else 'FAIL'}  {row.name} :: {row.detail}")
    confirmed = bool(rows) and all(row.passed for row in rows)
    overall &= confirmed
    print(f"  FINAL VERDICT {family.upper()}: {'CONFIRMED' if confirmed else 'REFUTED'} "
          f"({sum(row.passed for row in rows)}/{len(rows)} mandatory tests pass)")

print("=" * 100)
print(f"INITIAL BUILDER OVERALL: {'CONFIRMED' if initial_overall else 'REFUTED'}")
print(f"FINAL REPAIRED OVERALL: {'CONFIRMED' if overall else 'REFUTED'}")
sys.exit(0 if overall else 1)
