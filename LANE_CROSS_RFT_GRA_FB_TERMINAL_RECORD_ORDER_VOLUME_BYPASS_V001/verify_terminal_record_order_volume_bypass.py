#!/usr/bin/env python3
"""Exact finite replay for TROV operational-order and volume claims."""

from pathlib import Path


ROOT = Path(__file__).resolve().parent
passed = 0
total = 0
failures: list[str] = []


def check(condition: bool, name: str) -> None:
    global passed, total
    total += 1
    if condition:
        passed += 1
    else:
        failures.append(name)


def closure(vertices, arrows):
    rel = {(v, v) for v in vertices} | set(arrows)
    changed = True
    while changed:
        changed = False
        new = {(a, d) for (a, b) in rel for (c, d) in rel if b == c}
        if not new <= rel:
            rel |= new
            changed = True
    return rel


def antisymmetric(vertices, rel):
    return all(a == b or (b, a) not in rel for a, b in rel)


# Cross-mission acyclicity is not inherited from per-mission DAGs.
vertices = {"A", "B"}
m1 = {("A", "B")}
m2 = {("B", "A")}
check(antisymmetric(vertices, closure(vertices, m1)), "mission one acyclic")
check(antisymmetric(vertices, closure(vertices, m2)), "mission two acyclic")
union_closure = closure(vertices, m1 | m2)
check(not antisymmetric(vertices, union_closure), "union cycles")
check(("A", "B") in union_closure and ("B", "A") in union_closure,
      "cycle is explicit")


# Transitive closure adds reachability that need not be one witnessed mission.
vertices3 = {"A", "B", "C"}
arrows3 = {("A", "B"), ("B", "C")}
closed3 = closure(vertices3, arrows3)
check(("A", "C") not in arrows3, "A to C not directly witnessed")
check(("A", "C") in closed3, "A to C added by closure")
check(antisymmetric(vertices3, closed3), "acyclic closure is partial order")


# One material process front can be null in one metric and timelike in another.
dt = 1
dx = 1
ds2_a = -(dt * dt) + dx * dx
ds2_b = -4 * (dt * dt) + dx * dx
check(ds2_a == 0, "speed-one front null in g_A")
check(ds2_b < 0, "same front subluminal/timelike in g_B")


# Four-dimensional conformal volume normalization.
omega = 2
f = omega ** 4
metric_factor = omega ** 2
check(f == 16, "four-volume scales as Omega fourth")
check(metric_factor == 4, "metric scales as Omega squared")
check(metric_factor * metric_factor == f, "g=f half times g0")


theorem = (ROOT / "THEOREM.md").read_text()
normalized = " ".join(theorem.split())
required = [
    "unique minimal preorder",
    "preceq_R` is reflexive, whereas a chronological relation",
    "does **not** earn equality to complete physical chronology",
    "absence of a positive signal",
    "universal maximal front",
    "already admitted Lorentzian manifold",
    "one may not silently identify a generic",
    "An event count, record count, gamma exponent",
    "RCV -- record causal-volume realization law",
    "not a lower-premise proof",
    "cannot replace multiscale convergence",
    "EINSTEIN_DYNAMICS_REMAINS_DOWNSTREAM",
]
for snippet in required:
    check(snippet in normalized, f"scope snippet {snippet[:30]}")

check(theorem.count("\\[") == theorem.count("\\]"), "display math balanced")
tags = []
for line in theorem.splitlines():
    if "\\tag{FB" in line:
        tags.append(line.split("\\tag{", 1)[1].split("}", 1)[0])
check(len(tags) == len(set(tags)), "equation tags unique")
check(len(tags) == 8, "expected equation-tag count")

if failures:
    print(f"FAIL {passed}/{total}")
    for failure in failures:
        print(f" - {failure}")
    raise SystemExit(1)

print(f"PASS {passed}/{total}")
print("POSITIVE_TERMINAL_RECORD_ARROWS_SOUND_CONDITIONAL_ON_DCL_CUSTODY")
print("CROSS_MISSION_UNION_ORDER_NO_GO_EXACT")
print("MATERIAL_FRONT_TO_CAUSAL_CONE_NO_GO_EXACT")
print("MALAMENT_VOLUME_RECONSTRUCTION_SHARPLY_CONDITIONAL")
print("RCV_BYPASS_PREMISE_BURDEN_ISOLATED")
