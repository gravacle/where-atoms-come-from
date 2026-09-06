#!/usr/bin/env python3
"""Exact GL6CY development composition from the frozen GL6CU V002 engine.

This is a deliberately bounded development calculation.  It imports the
author-frozen GL6CU V002 engine by pinned byte hash, reconstructs its complete
selected-Q4 diagonal source Jets and all 64 active alternating-cycle Jets,
and composes both terms of CY02 on one explicitly declared Q4
plane-wave/literal-anchor probe path.

It does *not* compute a connected response, CTP kernel, 1PI kernel, physical
Ward residual, Einstein kernel, gravity, C_R, or G.  GL6CU V002 is used only
as development evidence until its distinct hostile audit passes.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CU_DIR = ROOT / "LANE_CROSS_RFT_GRA_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V002"
CU_ENGINE = CU_DIR / "derive_complete_six_pair_h6_source_jet.py"
CU_LEDGER = CU_DIR / "EXACT_LEDGER.json"
CU_MANIFEST = CU_DIR / "MANIFEST.sha256"
CU_SEAL = CU_DIR / "SEAL.sha256"

PINNED = {
    CU_ENGINE: "a62141731da025bdad51fde9223f9fa7d57f5a78bfff5d3375c4eeb1fa8ec84f",
    CU_LEDGER: "75c08daa934ab2de4a3300a6abb2a9fe4be72f99ce9790e7a73287d322c16de0",
    CU_MANIFEST: "58f313b1053756473fa8bfd1ff23b25e2f1740b1903bb366771e26cd3b3c89ec",
    CU_SEAL: "27cb7ad002e2539f756493697a1d90c9c259a1711c2c5d12a10970818b002810",
}

PAIR_ORDER = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
TENSOR_ORDER = ("xx", "yy", "zz", "xy", "xz", "yz")
CENTER_ORDER = ("x", "y", "z")
T = (
    (F(1), F(1), F(1)),
    (F(1), F(-1), F(-1)),
    (F(-1), F(1), F(-1)),
    (F(-1), F(-1), F(1)),
)
V = tuple(tuple(x / 2 for x in row) for row in T)
AXES = ((F(1), F(0), F(0)),
        (F(0), F(1), F(0)),
        (F(0), F(0), F(1)))
CHECKS = 0


def check(condition, label):
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def file_digest(path):
    return sha256(path.read_bytes()).hexdigest()


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def payload_digest(value):
    return sha256(canonical_json(value).encode()).hexdigest()


def qtext(value):
    value = F(value)
    return (str(value.numerator) if value.denominator == 1
            else f"{value.numerator}/{value.denominator}")


@dataclass(frozen=True)
class GQ:
    """Exact Gaussian rational a+i b."""

    re: F = F(0)
    im: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, "re", F(self.re))
        object.__setattr__(self, "im", F(self.im))

    def __add__(self, other):
        other = gq(other)
        return GQ(self.re + other.re, self.im + other.im)

    __radd__ = __add__

    def __neg__(self):
        return GQ(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-gq(other))

    def __rsub__(self, other):
        return gq(other) - self

    def __mul__(self, other):
        other = gq(other)
        return GQ(self.re * other.re - self.im * other.im,
                  self.re * other.im + self.im * other.re)

    __rmul__ = __mul__

    def conjugate(self):
        return GQ(self.re, -self.im)

    def __bool__(self):
        return bool(self.re or self.im)


def gq(value):
    return value if isinstance(value, GQ) else GQ(F(value))


I_POWERS = (GQ(1), GQ(0, 1), GQ(-1), GQ(0, -1))


def gtext(value):
    value = gq(value)
    if not value.im:
        return qtext(value.re)
    if not value.re:
        return f"{qtext(value.im)}i"
    sign = "+" if value.im > 0 else "-"
    return f"{qtext(value.re)}{sign}{qtext(abs(value.im))}i"


def dot(left, right):
    return sum((F(a) * F(b) for a, b in zip(left, right)), F(0))


def vadd(left, right):
    return tuple(F(a) + F(b) for a, b in zip(left, right))


def matvec(matrix, vector):
    return tuple(sum((F(matrix[i][j]) * F(vector[j])
                      for j in range(len(vector))), F(0))
                 for i in range(len(matrix)))


def symmetric_basis():
    """A_s with h_s=A_s+A_s^T in xx,yy,zz,xy,xz,yz order."""
    rows = []
    for i in range(3):
        matrix = [[F(0) for _ in range(3)] for _ in range(3)]
        matrix[i][i] = F(1, 2)
        rows.append(tuple(tuple(row) for row in matrix))
    for i, j in ((0, 1), (0, 2), (1, 2)):
        matrix = [[F(0) for _ in range(3)] for _ in range(3)]
        matrix[i][j] = matrix[j][i] = F(1, 2)
        rows.append(tuple(tuple(row) for row in matrix))
    return tuple(rows)


def dual_source_rows():
    """Columns of D_C^* in the standard six tensor-coordinate basis."""
    tensor_slots = ((0, 0), (1, 1), (2, 2),
                    (0, 1), (0, 2), (1, 2))
    weights = (F(1), F(1), F(1), F(2), F(2), F(2))
    rows = []
    for a, b in PAIR_ORDER:
        tensor = tuple(tuple(V[a][i] * V[b][j] + V[b][i] * V[a][j]
                             for j in range(3)) for i in range(3))
        rows.append(tuple(tensor[i][j] * weights[s]
                          for s, (i, j) in enumerate(tensor_slots)))
    return tuple(rows)


def t_direction(vector):
    return tuple(dot(vector, vadd(V[a], V[b])) for a, b in PAIR_ORDER)


def load_cu():
    for path, expected in PINNED.items():
        check(path.is_file() and file_digest(path) == expected,
              f"pinned CU V002 dependency {path.name}")
    seal = CU_SEAL.read_text().strip()
    check(seal == f"{PINNED[CU_MANIFEST]}  MANIFEST.sha256",
          "CU V002 seal points to pinned manifest")
    spec = importlib.util.spec_from_file_location("gl6cu_v002_development", CU_ENGINE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    check(tuple(module.PAIR_ORDER) == PAIR_ORDER, "CU pair order")
    check(len(module.NODES) == 128 and len(module.ACTIVE_CYCLES) == 64,
          "CU Q4 node and active-cycle census")
    return module


def build_cu_components(cu):
    """Reconstruct the exact CU Jets without calling its ledger writer."""
    engine = cu.LinkedOwners(cu.BASE)
    accumulators = {2: cu.JetAccumulator(), 4: cu.JetAccumulator(),
                    6: cu.JetAccumulator()}
    for edge in range(len(cu.EDGES)):
        for order in (2, 4, 6):
            accumulators[order].add(engine.coefficient((edge,), order))
    for support in sorted(cu.WEDGES, key=lambda row: tuple(sorted(row))):
        owner = tuple(support)
        for order in (4, 6):
            accumulators[order].add(engine.coefficient(owner, order))
    for support in sorted(cu.STARS | cu.PATHS,
                          key=lambda row: tuple(sorted(row))):
        owner = tuple(support)
        accumulators[6].add(engine.coefficient(owner, 6))
        engine.release(owner, 6)
    diagonal = {order: accumulator.freeze()
                for order, accumulator in accumulators.items()}

    bare_gradient = {}
    for node in cu.NODES:
        base = 6 * cu.NODE_INDEX[node]
        for pair, value in enumerate(cu.pair_word(cu.BASE, node)):
            bare_gradient[base + pair] = value
    bare = cu.Jet.linear(bare_gradient)

    cycles = []
    for index, cycle in enumerate(cu.ACTIVE_CYCLES):
        canonical, cycle_nodes = cu.cycle_canonical_jet(cycle, cu.BASE)
        cycles.append((index, tuple(cycle), tuple(cycle_nodes), canonical[6][0][1]))
    check(tuple(diagonal[n].value for n in (2, 4, 6)) ==
          (F(-128), F(-224, 3), F(-28576, 135)),
          "CU diagonal source-off regression")
    check(len(cycles) == 64 and all(row[3].value == F(-63, 8) for row in cycles),
          "CU 64 active-cycle source-off regression")
    return bare, diagonal, cycles


def contract_hessian(jet, left, right):
    answer = GQ()
    for (a, b), value in jet.hess.items():
        if a == b:
            answer += value * left.get(a, GQ()) * right.get(b, GQ())
        else:
            answer += value * (
                left.get(a, GQ()) * right.get(b, GQ())
                + left.get(b, GQ()) * right.get(a, GQ())
            )
    return answer


def contract_gradient(jet, field):
    return sum((value * field.get(key, GQ())
                for key, value in jet.grad.items()), GQ())


def matrix_payload(matrix):
    return [[gtext(value) for value in row] for row in matrix]


def matrix_add(left, right):
    return tuple(tuple(left[i][j] + right[i][j]
                       for j in range(len(left[0])))
                 for i in range(len(left)))


def zero_matrix():
    return tuple(tuple(GQ() for _ in CENTER_ORDER) for _ in TENSOR_ORDER)


def character_phase(character, cell):
    exponent = sum(int(character[i]) * int(cell[i]) for i in range(3)) % 4
    return I_POWERS[exponent]


def source_fields(cu, character):
    """Return eta_s, eta_r, eta_sr for the declared CY02 path.

    At every node v with cell character chi_v and local orientation
    sigma_v=+1(P),-1(C), the path is

      F_v = I + u_s chi_v^* A_s,
      b_v = u_r chi_v e_r,
      eta_total = u_s chi_v^* (D_C^* h_s) + eta_phys(F_v,b_v).

    Thus A_r=0, b_s=0, b_r=chi_v e_r, and b_sr=0 by the explicitly
    bilinear-free literal-anchor definition.  CY14 nevertheless gives the
    nonzero mixed source eta_sr=2 sigma_v t_(A_s e_r); it is evaluated rather
    than set to zero.
    """
    qrows = dual_source_rows()
    abasis = symmetric_basis()
    tensor = [dict() for _ in TENSOR_ORDER]
    center = [dict() for _ in CENTER_ORDER]
    second = [[dict() for _ in CENTER_ORDER] for _ in TENSOR_ORDER]
    for node in cu.NODES:
        base = 6 * cu.NODE_INDEX[node]
        sigma = F(1) if node[0] == "P" else F(-1)
        phase = character_phase(character, node[1])
        inverse_phase = phase.conjugate()
        for s in range(6):
            for pair in range(6):
                value = inverse_phase * qrows[pair][s]
                if value:
                    tensor[s][base + pair] = value
        for r, axis in enumerate(AXES):
            first = t_direction(axis)
            for pair in range(6):
                value = phase * (-2 * sigma * first[pair])
                if value:
                    center[r][base + pair] = value
            for s, a_s in enumerate(abasis):
                mixed_vector = matvec(a_s, axis)
                mixed = t_direction(mixed_vector)
                for pair in range(6):
                    # CY14 with A_r=0, b_s=0, b_sr=0:
                    # eta_sr=+2 sigma (A_s b_r).(v_a+v_b).
                    value = GQ(2 * sigma * mixed[pair])
                    if value:
                        second[s][r][base + pair] = value
    check(any(second[s][r] for s in range(6) for r in range(3)),
          "CY02 second source field is explicitly nonzero")
    return tensor, center, second


def compose(jet, fields):
    tensor, center, second = fields
    first = tuple(tuple(contract_hessian(jet, tensor[s], center[r])
                        for r in range(3)) for s in range(6))
    contact = tuple(tuple(contract_gradient(jet, second[s][r])
                          for r in range(3)) for s in range(6))
    return {
        "H_AB_eta_s_eta_r": first,
        "H_A_eta_sr": contact,
        "sum": matrix_add(first, contact),
    }


def coordinate_data(cu, key):
    node = cu.NODES[key // 6]
    return node, key % 6


def add_kernel_value(kernel, delta, s, r, value):
    if not value:
        return
    block = kernel.setdefault(delta, [[F(0) for _ in CENTER_ORDER]
                                      for _ in TENSOR_ORDER])
    block[s][r] += F(value)


def compile_hessian_character_kernel(cu, jet):
    """Compile H_AB eta_A,s eta_B,r by exact Q4 translation difference.

    The expensive arbitrary-nonuniform Hessian is traversed once.  Its
    contraction on every Q4 character is then the finite Fourier polynomial

        sum_delta K_sr(delta) i^(m dot delta).

    This is also the finite-range stencil export needed by an overlap-family
    calculation; no cell-level Ward interpretation is made.
    """
    qrows = dual_source_rows()
    tdirs = tuple(t_direction(axis) for axis in AXES)
    kernel = {}
    for (left, right), hvalue in jet.hess.items():
        left_node, left_pair = coordinate_data(cu, left)
        right_node, right_pair = coordinate_data(cu, right)
        left_cell, right_cell = left_node[1], right_node[1]
        left_sigma = F(1) if left_node[0] == "P" else F(-1)
        right_sigma = F(1) if right_node[0] == "P" else F(-1)
        forward_delta = tuple((right_cell[i] - left_cell[i]) % 4 for i in range(3))
        reverse_delta = tuple((-value) % 4 for value in forward_delta)
        for s in range(6):
            qleft = qrows[left_pair][s]
            if qleft:
                for r in range(3):
                    add_kernel_value(
                        kernel, forward_delta, s, r,
                        hvalue * qleft * (-2 * right_sigma * tdirs[r][right_pair]),
                    )
            if left != right:
                qright = qrows[right_pair][s]
                if qright:
                    for r in range(3):
                        add_kernel_value(
                            kernel, reverse_delta, s, r,
                            hvalue * qright * (-2 * left_sigma * tdirs[r][left_pair]),
                        )
    return {
        delta: tuple(tuple(value for value in row) for row in block)
        for delta, block in kernel.items()
        if any(value for row in block for value in row)
    }


def evaluate_character_kernel(kernel, character):
    matrix = [[GQ() for _ in CENTER_ORDER] for _ in TENSOR_ORDER]
    for delta, block in kernel.items():
        phase = I_POWERS[sum(character[i] * delta[i] for i in range(3)) % 4]
        for s in range(6):
            for r in range(3):
                if block[s][r]:
                    matrix[s][r] += phase * block[s][r]
    return tuple(tuple(value for value in row) for row in matrix)


def add_character_kernels(*kernels):
    result = {}
    for kernel in kernels:
        for delta, block in kernel.items():
            target = result.setdefault(delta, [[F(0) for _ in CENTER_ORDER]
                                                for _ in TENSOR_ORDER])
            for s in range(6):
                for r in range(3):
                    target[s][r] += block[s][r]
    return {
        delta: tuple(tuple(value for value in row) for row in block)
        for delta, block in result.items()
        if any(value for row in block for value in row)
    }


def character_kernel_payload(kernel):
    return [
        {
            "delta_mod_4": list(delta),
            "matrix": [[qtext(value) for value in row] for row in kernel[delta]],
        }
        for delta in sorted(kernel)
    ]


def compile_second_chain_matrix(cu, jet):
    """Compute H_A eta_A,sr once; eta_sr has zero net Q4 character."""
    abasis = symmetric_basis()
    result = [[GQ() for _ in CENTER_ORDER] for _ in TENSOR_ORDER]
    for key, value in jet.grad.items():
        node, pair = coordinate_data(cu, key)
        sigma = F(1) if node[0] == "P" else F(-1)
        for s, a_s in enumerate(abasis):
            for r, axis in enumerate(AXES):
                mixed = t_direction(matvec(a_s, axis))[pair]
                if mixed:
                    result[s][r] += value * (2 * sigma * mixed)
    return tuple(tuple(value for value in row) for row in result)


def compile_component(cu, jet):
    return {
        "hessian_kernel": compile_hessian_character_kernel(cu, jet),
        "second_matrix": compile_second_chain_matrix(cu, jet),
    }


def evaluate_component(compiled, character):
    first = evaluate_character_kernel(compiled["hessian_kernel"], character)
    second = compiled["second_matrix"]
    return {
        "H_AB_eta_s_eta_r": first,
        "H_A_eta_sr": second,
        "sum": matrix_add(first, second),
    }


def raw_pair_unit_center_matrix(cu, jet):
    """Fast exact replay of the preliminary raw-pair/unit-center matrix."""
    tdirs = tuple(t_direction(axis) for axis in AXES)
    result = [[F(0) for _ in CENTER_ORDER] for _ in TENSOR_ORDER]
    for (left, right), value in jet.hess.items():
        left_node, left_pair = coordinate_data(cu, left)
        right_node, right_pair = coordinate_data(cu, right)
        left_sigma = F(1) if left_node[0] == "P" else F(-1)
        right_sigma = F(1) if right_node[0] == "P" else F(-1)
        for r in range(3):
            result[left_pair][r] += value * right_sigma * tdirs[r][right_pair]
            if left != right:
                result[right_pair][r] += value * left_sigma * tdirs[r][left_pair]
    return tuple(tuple(GQ(value) for value in row) for row in result)


def classify_character_function(rows):
    groups = {}
    for character, matrix in rows:
        payload = matrix_payload(matrix)
        key = canonical_json(payload)
        bucket = groups.setdefault(key, {"matrix": payload, "characters": []})
        bucket["characters"].append(list(character))
    result = list(groups.values())
    result.sort(key=lambda row: (len(row["characters"]), canonical_json(row["matrix"])))
    return {
        "class_count": len(result),
        "classes": result,
        "all_64_values_sha256": payload_digest([
            {"character": list(character), "matrix": matrix_payload(matrix)}
            for character, matrix in rows
        ]),
    }


def run():
    started = time.monotonic()
    cu = load_cu()
    bare, diagonal, cycles = build_cu_components(cu)
    characters = tuple((x, y, z) for x in range(4)
                       for y in range(4) for z in range(4))

    component_rows = {
        "h0_bare": {term: [] for term in
                    ("H_AB_eta_s_eta_r", "H_A_eta_sr", "sum")},
        "h2_diagonal": {term: [] for term in
                        ("H_AB_eta_s_eta_r", "H_A_eta_sr", "sum")},
        "h4_diagonal": {term: [] for term in
                        ("H_AB_eta_s_eta_r", "H_A_eta_sr", "sum")},
        "h6_diagonal": {term: [] for term in
                        ("H_AB_eta_s_eta_r", "H_A_eta_sr", "sum")},
        "h6_active_cycles_64": {term: [] for term in
                                ("H_AB_eta_s_eta_r", "H_A_eta_sr", "sum")},
        "h6_diagonal_plus_cycles": {term: [] for term in
                                    ("H_AB_eta_s_eta_r", "H_A_eta_sr", "sum")},
    }
    cycle_functions = [[] for _ in cycles]

    # Compile each sparse nonuniform Jet once into an exact translation-
    # difference stencil.  The 64-character evaluation below never reruns an
    # owner or rescans the 52,608-entry h6 diagonal Hessian.
    compiled_components = {
        "h0_bare": compile_component(cu, bare),
        "h2_diagonal": compile_component(cu, diagonal[2]),
        "h4_diagonal": compile_component(cu, diagonal[4]),
        "h6_diagonal": compile_component(cu, diagonal[6]),
    }
    compiled_cycles = [compile_component(cu, cycle[3]) for cycle in cycles]
    cycle_kernel_sum = add_character_kernels(*[
        row["hessian_kernel"] for row in compiled_cycles
    ])
    total_h6_kernel = add_character_kernels(
        compiled_components["h6_diagonal"]["hessian_kernel"],
        cycle_kernel_sum,
    )

    for character in characters:
        pieces = {
            name: evaluate_component(compiled, character)
            for name, compiled in compiled_components.items()
        }
        cycle_parts = [evaluate_component(compiled, character)
                       for compiled in compiled_cycles]
        cycle_sum = {term: zero_matrix() for term in
                     ("H_AB_eta_s_eta_r", "H_A_eta_sr", "sum")}
        for index, part in enumerate(cycle_parts):
            for term in cycle_sum:
                cycle_sum[term] = matrix_add(cycle_sum[term], part[term])
            cycle_functions[index].append({
                "character": list(character),
                "H_AB_eta_s_eta_r": matrix_payload(part["H_AB_eta_s_eta_r"]),
                "H_A_eta_sr": matrix_payload(part["H_A_eta_sr"]),
                "sum": matrix_payload(part["sum"]),
            })
        pieces["h6_active_cycles_64"] = cycle_sum
        pieces["h6_diagonal_plus_cycles"] = {
            term: matrix_add(pieces["h6_diagonal"][term], cycle_sum[term])
            for term in cycle_sum
        }
        for name, piece in pieces.items():
            for term, matrix in piece.items():
                component_rows[name][term].append((character, matrix))

    # The full second chain term is nonzero as a source field but its
    # contraction vanishes on this declared path for every frozen component.
    # This is an executed result, not an imposed eta_sr=0 assumption.
    for name in component_rows:
        check(all(matrix == zero_matrix() for _character, matrix in
                  component_rows[name]["H_A_eta_sr"]),
              f"executed second-chain contraction vanishes on path: {name}")

    # Reproduce the old raw-pair/unit-center k=0 diagnostic directly from
    # freshly rebuilt Jets before changing to D_C^* tensor rows and the -2
    # literal physical center normalization.
    # Directly exercise the source-field builder once so the declared
    # nonzero eta_sr path is checked independently of the compiled formula.
    source_fields(cu, (0, 0, 0))

    def old_matrix(jet):
        return raw_pair_unit_center_matrix(cu, jet)

    p_matrix = (
        (GQ(), GQ(), GQ()),
        (GQ(), GQ(1), GQ()),
        (GQ(), GQ(), GQ(-1)),
        (GQ(), GQ(), GQ(1)),
        (GQ(), GQ(-1), GQ()),
        (GQ(), GQ(), GQ()),
    )
    check(old_matrix(diagonal[2]) == zero_matrix(), "preliminary h2 matrix")
    check(old_matrix(diagonal[4]) ==
          tuple(tuple(F(8032, 27) * x for x in row) for row in p_matrix),
          "preliminary h4 matrix")
    check(old_matrix(diagonal[6]) ==
          tuple(tuple(F(5010712, 6075) * x for x in row) for row in p_matrix),
          "preliminary h6 diagonal matrix")
    old_cycle = zero_matrix()
    for cycle in cycles:
        old_cycle = matrix_add(old_cycle, old_matrix(cycle[3]))
    check(old_cycle == tuple(tuple(F(-528) * x for x in row)
                             for row in p_matrix),
          "preliminary h6 64-cycle matrix")

    function_result = {}
    for component, terms in component_rows.items():
        function_result[component] = {
            term: classify_character_function(rows)
            for term, rows in terms.items()
        }

    cycle_records = []
    for (index, edges, nodes, _jet), compiled, values in zip(
            cycles, compiled_cycles, cycle_functions):
        stencil_payload = character_kernel_payload(compiled["hessian_kernel"])
        cycle_records.append({
            "index": index,
            "edges": list(edges),
            "nodes": [cu.NODE_INDEX[node] for node in nodes],
            "translation_stencil_term_count": len(stencil_payload),
            "translation_stencil_sha256": payload_digest(stencil_payload),
            "all_64_two_term_matrices_sha256": payload_digest(values),
        })

    result = {
        "schema": "GL6CY_DEVELOPMENT_CU_V002_RAW_SOURCE_COMPOSITION_V001",
        "status": "PASS_DEVELOPMENT_ONLY_PENDING_GL6CU_V002_DISTINCT_HOSTILE_AUDIT",
        "checks": CHECKS,
        "runtime_seconds": round(time.monotonic() - started, 3),
        "scope": {
            "object": "raw Hamiltonian source-family CY02 composition",
            "not": ["connected response", "CTP kernel", "1PI kernel",
                    "physical Ward residual", "gravity", "C_R", "G"],
            "cu_status": "development evidence pending distinct hostile audit",
            "q4": "128 nodes; 256 links; 64 active alternating cycles",
        },
        "orders": {
            "diagonal": [0, 2, 4, 6],
            "configuration_changing": "64 active h6 cycle owners",
        },
        "basis": {
            "pair_order": [list(pair) for pair in PAIR_ORDER],
            "tensor_test_order": list(TENSOR_ORDER),
            "center_order": list(CENTER_ORDER),
            "D_C_star_columns_pair_by_tensor": [
                [qtext(value) for value in row] for row in dual_source_rows()
            ],
        },
        "declared_path": {
            "character": "chi_m(v)=i^(m dot cell(v)), m in Z4^3; same cell phase on P and C",
            "orientation": "sigma(P)=+1, sigma(C)=-1",
            "F": "F_v=I+u_s chi_m(v)^* A_s",
            "A_s": "symmetric half-basis with h_s=A_s+A_s^T",
            "A_r": "0",
            "b_s": "0",
            "b_r": "chi_m(v) e_r",
            "b_sr": "0 by declared bilinear-free literal-anchor path",
            "eta_s": "chi_m(v)^* (D_C^* h_s)",
            "eta_r": "-2 sigma(v) chi_m(v) t_r",
            "eta_sr": "+2 sigma(v) t_(A_s e_r), nonzero and executed",
            "warning": "this probe-plus-anchor path is an algebraic development diagnostic, not an authenticated physical source law",
        },
        "character_functions": function_result,
        "cycle_owner_functions": {
            "count": 64,
            "all_records_sha256": payload_digest(cycle_records),
            "records": cycle_records,
        },
        "translation_difference_stencils": {
            "meaning": "exact Q4 finite-range coefficient K_sr(delta); character function is sum_delta K_sr(delta) i^(m dot delta)",
            "use_boundary": "local raw-source input for G_L overlap/seam gluing; not a cell Ward inference",
            "h2_diagonal": character_kernel_payload(
                compiled_components["h2_diagonal"]["hessian_kernel"]),
            "h4_diagonal": character_kernel_payload(
                compiled_components["h4_diagonal"]["hessian_kernel"]),
            "h6_diagonal": character_kernel_payload(
                compiled_components["h6_diagonal"]["hessian_kernel"]),
            "h6_active_cycles_64": character_kernel_payload(cycle_kernel_sum),
            "h6_diagonal_plus_cycles": character_kernel_payload(total_h6_kernel),
        },
        "k0_matrices": {
            component: {
                term: matrix_payload(dict(rows)[(0, 0, 0)])
                for term, rows in terms.items()
            }
            for component, terms in component_rows.items()
        },
        "preliminary_raw_pair_unit_center_regression": {
            "h2_diagonal": matrix_payload(old_matrix(diagonal[2])),
            "h4_diagonal": matrix_payload(old_matrix(diagonal[4])),
            "h6_diagonal": matrix_payload(old_matrix(diagonal[6])),
            "h6_active_cycles_64": matrix_payload(old_cycle),
        },
        "dependencies": {str(path.relative_to(ROOT)): digest
                         for path, digest in PINNED.items()},
    }
    check(result["cycle_owner_functions"]["count"] == 64,
          "all 64 cycle functions exposed")
    result["checks"] = CHECKS
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(f"PASS__GL6CY_DEVELOPMENT_RAW_SOURCE_COMPOSITION__{result['checks']}/{result['checks']}")
        print(f"RUNTIME_SECONDS={result['runtime_seconds']}")
        print("OBJECT=RAW_HAMILTONIAN_CY02__NOT_CONNECTED_CTP_1PI_OR_WARD")
        print("PIECES=H0_BARE__H2_H4_H6_DIAGONAL__64_ACTIVE_H6_CYCLES")
        print("SECOND_TERM=ETA_SR_NONZERO_AND_EXPLICITLY_CONTRACTED")
        print("STATUS=DEVELOPMENT_ONLY_PENDING_DISTINCT_GL6CU_V002_HOSTILE_AUDIT")


if __name__ == "__main__":
    main()
