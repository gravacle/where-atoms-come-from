#!/usr/bin/env python3
"""Independent exact audit of GL6CY's raw six-output source composition.

This reconstruction imports only the pinned, previously independently audited
GL6CU V002 source-Jet engine.  It does not import any code from the GL6CY
development target.  Its character compiler and source directions are stated
here independently, then compared with a fresh full replay of the target.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "GL6CY_DEVELOPMENT_CU_V002_RAW_SOURCE_COMPOSITION"
TARGET_PROGRAM = TARGET / "derive_raw_source_composition.py"
CU_DIR = ROOT / "LANE_CROSS_RFT_GRA_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V002"
CU_PROGRAM = CU_DIR / "derive_complete_six_pair_h6_source_jet.py"

TARGET_PINS = {
    "ADDENDUM.md": "a5a9d4d99db6da268182fc9884d56cd8caf9af3bac122c7d39a6eeeda77dd528",
    "README.md": "f63fd2e2cb4fde1ce451a98a14dce50c368fe28e95a8213a47e500cca940d87a",
    "VERIFICATION.txt": "f23fd4307f8ac52afd18284a6f15129e33add699cbea2a68f229472367f0e655",
    "derive_raw_source_composition.py": "e494e2ba3612f68e4be8283c6e4a1c1c5861d2465c668d898717e9e02ed82141",
    "verify_development.py": "312ceb388693b7342b90c2498713f190f66a885cdf937ac347ab033804f28191",
    "MANIFEST.sha256": "c2ec2e52b684b39c6a2cf8efd79c48efd11de4121c858adde671f6f69b10e3dd",
}
CU_PINS = {
    "derive_complete_six_pair_h6_source_jet.py": "a62141731da025bdad51fde9223f9fa7d57f5a78bfff5d3375c4eeb1fa8ec84f",
    "EXACT_LEDGER.json": "75c08daa934ab2de4a3300a6abb2a9fe4be72f99ce9790e7a73287d322c16de0",
    "MANIFEST.sha256": "58f313b1053756473fa8bfd1ff23b25e2f1740b1903bb366771e26cd3b3c89ec",
    "SEAL.sha256": "27cb7ad002e2539f756493697a1d90c9c259a1711c2c5d12a10970818b002810",
}

PAIR_ORDER = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
TENSOR_ORDER = ("xx", "yy", "zz", "xy", "xz", "yz")
CHARACTERS = tuple((x, y, z) for x in range(4) for y in range(4)
                   for z in range(4))
T = ((F(1), F(1), F(1)), (F(1), F(-1), F(-1)),
     (F(-1), F(1), F(-1)), (F(-1), F(-1), F(1)))
V = tuple(tuple(value / 2 for value in row) for row in T)
AXES = ((F(1), F(0), F(0)), (F(0), F(1), F(0)),
        (F(0), F(0), F(1)))
I_POWERS = ((F(1), F(0)), (F(0), F(1)), (F(-1), F(0)), (F(0), F(-1)))
CHECKS = 0


def check(condition, label):
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def qtext(value):
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def gadd(left, right):
    return left[0] + right[0], left[1] + right[1]


def gmul(left, right):
    return (left[0] * right[0] - left[1] * right[1],
            left[0] * right[1] + left[1] * right[0])


def gscale(value, scalar):
    return value[0] * scalar, value[1] * scalar


def gtext(value):
    real, imag = value
    if not imag:
        return qtext(real)
    if not real:
        return f"{qtext(imag)}i"
    return f"{qtext(real)}{'+' if imag > 0 else '-'}{qtext(abs(imag))}i"


def dot(left, right):
    return sum((F(a) * F(b) for a, b in zip(left, right)), F(0))


def matvec(matrix, vector):
    return tuple(sum((F(matrix[i][j]) * F(vector[j]) for j in range(3)), F(0))
                 for i in range(3))


def qrows():
    slots = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
    weights = (F(1), F(1), F(1), F(2), F(2), F(2))
    rows = []
    for a, b in PAIR_ORDER:
        dyad = tuple(tuple(V[a][i] * V[b][j] + V[b][i] * V[a][j]
                           for j in range(3)) for i in range(3))
        rows.append(tuple(dyad[i][j] * weights[s]
                          for s, (i, j) in enumerate(slots)))
    return tuple(rows)


def symmetric_half_basis():
    output = []
    for i in range(3):
        matrix = [[F(0) for _ in range(3)] for _ in range(3)]
        matrix[i][i] = F(1, 2)
        output.append(tuple(tuple(row) for row in matrix))
    for i, j in ((0, 1), (0, 2), (1, 2)):
        matrix = [[F(0) for _ in range(3)] for _ in range(3)]
        matrix[i][j] = matrix[j][i] = F(1, 2)
        output.append(tuple(tuple(row) for row in matrix))
    return tuple(output)


def tdir(vector):
    return tuple(dot(vector, tuple(V[a][i] + V[b][i] for i in range(3)))
                 for a, b in PAIR_ORDER)


def load_cu():
    for name, expected in CU_PINS.items():
        check((CU_DIR / name).is_file() and digest(CU_DIR / name) == expected,
              f"pinned GL6CU byte {name}")
    check((CU_DIR / "SEAL.sha256").read_text().strip() ==
          "58f313b1053756473fa8bfd1ff23b25e2f1740b1903bb366771e26cd3b3c89ec  MANIFEST.sha256",
          "GL6CU seal")
    spec = importlib.util.spec_from_file_location("independent_gl6cu_v002", CU_PROGRAM)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    check(tuple(module.PAIR_ORDER) == PAIR_ORDER, "pair order")
    check(len(module.NODES) == 128 and len(module.EDGES) == 256 and
          len(module.ACTIVE_CYCLES) == 64, "Q4 census")
    return module


def build_source_jets(cu):
    """Independent owner aggregation: each selected owner is included once."""
    owners = cu.LinkedOwners(cu.BASE)
    sums = {order: cu.JetAccumulator() for order in (2, 4, 6)}
    for edge in range(len(cu.EDGES)):
        for order in (2, 4, 6):
            sums[order].add(owners.coefficient((edge,), order))
    for support in sorted(cu.WEDGES, key=lambda row: tuple(sorted(row))):
        for order in (4, 6):
            sums[order].add(owners.coefficient(tuple(support), order))
    for support in sorted(cu.STARS | cu.PATHS, key=lambda row: tuple(sorted(row))):
        owner = tuple(support)
        sums[6].add(owners.coefficient(owner, 6))
        owners.release(owner, 6)
    diagonal = {order: sums[order].freeze() for order in sums}
    gradient = {}
    for node in cu.NODES:
        offset = 6 * cu.NODE_INDEX[node]
        for pair, value in enumerate(cu.pair_word(cu.BASE, node)):
            gradient[offset + pair] = value
    bare = cu.Jet.linear(gradient)
    cycles = []
    for cycle in cu.ACTIVE_CYCLES:
        canonical, _nodes = cu.cycle_canonical_jet(cycle, cu.BASE)
        cycles.append(canonical[6][0][1])
    check(tuple(diagonal[order].value for order in (2, 4, 6)) ==
          (F(-128), F(-224, 3), F(-28576, 135)), "diagonal source-off values")
    check(len(cycles) == 64 and all(jet.value == F(-63, 8) for jet in cycles),
          "active alternating-cycle source-off values")
    return {"h0_bare": bare, "h2_diagonal": diagonal[2],
            "h4_diagonal": diagonal[4], "h6_diagonal": diagonal[6],
            "h6_active_cycles_64": sum_jets(cu, cycles)}, cycles


def sum_jets(cu, jets):
    accumulator = cu.JetAccumulator()
    for jet in jets:
        accumulator.add(jet)
    return accumulator.freeze()


def node_pair(cu, coordinate):
    return cu.NODES[coordinate // 6], coordinate % 6


def zero_matrix():
    return tuple(tuple((F(0), F(0)) for _ in range(3)) for _ in range(6))


def add_matrix(left, right):
    return tuple(tuple(gadd(left[s][r], right[s][r]) for r in range(3))
                 for s in range(6))


def compile_first_kernel(cu, jet):
    """Exact Fourier coefficient of H_AB eta_s eta_r, independently coded."""
    q = qrows()
    directions = tuple(tdir(axis) for axis in AXES)
    coefficients = {}

    def add(delta, s, r, value):
        if not value:
            return
        block = coefficients.setdefault(delta, [[F(0) for _ in range(3)]
                                                 for _ in range(6)])
        block[s][r] += value

    for (left, right), hessian in jet.hess.items():
        left_node, left_pair = node_pair(cu, left)
        right_node, right_pair = node_pair(cu, right)
        ls = F(1) if left_node[0] == "P" else F(-1)
        rs = F(1) if right_node[0] == "P" else F(-1)
        delta = tuple((right_node[1][i] - left_node[1][i]) % 4 for i in range(3))
        reverse = tuple((-item) % 4 for item in delta)
        for s in range(6):
            if q[left_pair][s]:
                for r in range(3):
                    add(delta, s, r, hessian * q[left_pair][s] *
                        (-2 * rs * directions[r][right_pair]))
            if left != right and q[right_pair][s]:
                for r in range(3):
                    add(reverse, s, r, hessian * q[right_pair][s] *
                        (-2 * ls * directions[r][left_pair]))
    return {delta: tuple(tuple(value for value in row) for row in block)
            for delta, block in coefficients.items()
            if any(value for row in block for value in row)}


def evaluate_kernel(kernel, character):
    answer = [[(F(0), F(0)) for _ in range(3)] for _ in range(6)]
    for delta, block in kernel.items():
        power = sum(character[i] * delta[i] for i in range(3)) % 4
        phase = I_POWERS[power]
        for s in range(6):
            for r in range(3):
                if block[s][r]:
                    answer[s][r] = gadd(answer[s][r], gscale(phase, block[s][r]))
    return tuple(tuple(value for value in row) for row in answer)


def second_chain(cu, jet):
    answer = [[F(0) for _ in range(3)] for _ in range(6)]
    for coordinate, derivative in jet.grad.items():
        node, pair = node_pair(cu, coordinate)
        sigma = F(1) if node[0] == "P" else F(-1)
        for s, matrix in enumerate(symmetric_half_basis()):
            for r, axis in enumerate(AXES):
                term = tdir(matvec(matrix, axis))[pair]
                answer[s][r] += derivative * 2 * sigma * term
    return tuple(tuple((value, F(0)) for value in row) for row in answer)


def matrix_payload(matrix):
    return [[gtext(value) for value in row] for row in matrix]


def stencil_payload(kernel):
    return [{"delta_mod_4": list(delta),
             "matrix": [[qtext(value) for value in row] for row in kernel[delta]]}
            for delta in sorted(kernel)]


def target_rows(target_data, component):
    output = {}
    for group in target_data["character_functions"][component]["sum"]["classes"]:
        for character in group["characters"]:
            key = tuple(character)
            check(key not in output, f"unique target character {component} {key}")
            output[key] = group["matrix"]
    check(set(output) == set(CHARACTERS), f"64 target characters {component}")
    return output


def main():
    for name, expected in TARGET_PINS.items():
        check((TARGET / name).is_file() and digest(TARGET / name) == expected,
              f"pinned target byte {name}")
    check((TARGET / "SEAL.sha256").read_text().strip() ==
          "c2ec2e52b684b39c6a2cf8efd79c48efd11de4121c858adde671f6f69b10e3dd  MANIFEST.sha256",
          "target seal")

    cu = load_cu()
    jets, individual_cycles = build_source_jets(cu)
    kernels = {name: compile_first_kernel(cu, jet) for name, jet in jets.items()}
    cycle_kernel_sum = compile_first_kernel(cu, jets["h6_active_cycles_64"])
    diag_plus_cycle = {}
    for delta in set(kernels["h6_diagonal"]) | set(cycle_kernel_sum):
        left = kernels["h6_diagonal"].get(delta, tuple(tuple(F(0) for _ in range(3)) for _ in range(6)))
        right = cycle_kernel_sum.get(delta, tuple(tuple(F(0) for _ in range(3)) for _ in range(6)))
        block = tuple(tuple(left[s][r] + right[s][r] for r in range(3)) for s in range(6))
        if any(value for row in block for value in row):
            diag_plus_cycle[delta] = block
    kernels["h6_diagonal_plus_cycles"] = diag_plus_cycle

    first_values = {name: {character: evaluate_kernel(kernel, character)
                           for character in CHARACTERS}
                    for name, kernel in kernels.items()}
    second_values = {name: second_chain(cu, jet) for name, jet in jets.items()}
    second_values["h6_diagonal_plus_cycles"] = add_matrix(
        second_values["h6_diagonal"],
        second_values["h6_active_cycles_64"],
    )
    for name, value in second_values.items():
        check(value == zero_matrix(), f"executed nonzero eta_sr contracts to zero {name}")
    for index, cycle in enumerate(individual_cycles):
        check(second_chain(cu, cycle) == zero_matrix(),
              f"executed eta_sr cycle cancellation {index}")

    expected_k0 = {
        "h0_bare": F(0), "h2_diagonal": F(0),
        "h4_diagonal": F(32128, 27),
        "h6_diagonal": F(20042848, 6075),
        "h6_active_cycles_64": F(-2112),
        "h6_diagonal_plus_cycles": F(7212448, 6075),
    }
    support = {(3, 2): F(-1), (4, 1): F(1)}
    for name, coefficient in expected_k0.items():
        matrix = first_values[name][(0, 0, 0)]
        for s in range(6):
            for r in range(3):
                wanted = (coefficient * support.get((s, r), F(0)), F(0))
                check(matrix[s][r] == wanted, f"exact k0 {name} {s} {r}")

    # Fresh full target replay.  Its structured 64-character output is
    # compared component by component to the independent compilation above.
    replay = subprocess.run([sys.executable, "-B", str(TARGET_PROGRAM), "--json"],
                            cwd=ROOT, capture_output=True, text=True, check=False)
    check(replay.returncode == 0, "fresh target replay exits zero")
    target_data = json.loads(replay.stdout)
    check(target_data["checks"] == 21, "fresh target replay 21 checks")
    check(target_data["scope"]["object"] == "raw Hamiltonian source-family CY02 composition",
          "target scope exact")
    for name in ("h0_bare", "h2_diagonal", "h4_diagonal", "h6_diagonal",
                 "h6_active_cycles_64", "h6_diagonal_plus_cycles"):
        target = target_rows(target_data, name)
        for character in CHARACTERS:
            check(target[character] == matrix_payload(first_values[name][character]),
                  f"all-character composition agreement {name} {character}")
        # The target's all-character table includes the combined h6 block,
        # which is already compared above.  Its k0 convenience table omits
        # only that redundant combined entry, so do not turn a presentation
        # omission into a false arithmetic failure.
        k0_entry = target_data["k0_matrices"].get(name)
        if k0_entry is not None and "H_A_eta_sr" in k0_entry:
            check(k0_entry["H_A_eta_sr"] == matrix_payload(second_values[name]),
                  f"second chain target agreement {name}")
        else:
            check(name == "h6_diagonal_plus_cycles",
                  f"only combined h6 k0 contact convenience entry may be absent: {name}")

    stencils = target_data["translation_difference_stencils"]
    for name in ("h2_diagonal", "h4_diagonal", "h6_diagonal",
                 "h6_active_cycles_64", "h6_diagonal_plus_cycles"):
        check(stencils[name] == stencil_payload(kernels[name]),
              f"translation stencil agreement {name}")

    check(target_data["cycle_owner_functions"]["count"] == 64,
          "target exposes all cycle owner functions")
    result = {
        "schema": "AUDIT_G_GL6CY_RAW_SOURCE_COMPOSITION_V001",
        "target": "GL6CY_DEVELOPMENT_CU_V002_RAW_SOURCE_COMPOSITION",
        "target_manifest_sha256": TARGET_PINS["MANIFEST.sha256"],
        "disposition": "PASS",
        "independent_checks": CHECKS,
        "independent_scope": "all 64 Q4 characters; six source rows; three centre directions; raw CY02 composition",
        "input_custody": "GL6CU V002 independently audited source-Jet",
        "raw_result": {
            "second_source": "eta_sr nonzero; H_A eta_sr is exactly zero only on the declared probe-plus-literal-anchor path",
            "k0": "H0=H2=0; H4=(32128/27)S; H6diag+cycles=(7212448/6075)S",
            "character_comparison": "fresh full target replay agrees at all 64 characters",
            "cycle_owners": 64,
        },
        "not_earned": ["connected response", "CTP kernel", "1PI quotient",
                       "physical Ward null", "GL6CR", "Einstein gravity", "C_R", "G"],
    }
    (HERE / "INDEPENDENT_RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"PASS__INDEPENDENT_GL6CY_RAW_SOURCE_COMPOSITION__{CHECKS}/{CHECKS}")
    print("REPLAY=FRESH_FULL_TARGET_64_CHARACTER_COMPARISON")
    print("SOURCE=H_AB_ETA_S_ETA_R_PLUS_H_A_ETA_SR__SECOND_SOURCE_NONZERO_AND_CONTRACTED")
    print("K0=H0_H2_ZERO__H4_32128_OVER_27_S__H6_DIAG_PLUS_CYCLE_7212448_OVER_6075_S")
    print("CEILING=RAW_HAMILTONIAN_ONLY__NO_CONNECTED_CTP_1PI_WARD_GL6CR_GRAVITY_CR_OR_G")


if __name__ == "__main__":
    main()
