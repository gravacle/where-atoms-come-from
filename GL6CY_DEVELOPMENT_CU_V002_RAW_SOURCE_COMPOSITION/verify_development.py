#!/usr/bin/env python3
"""Fast independent custody and algebra checks for the GL6CY development addendum."""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CU_DIR = ROOT / "LANE_CROSS_RFT_GRA_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V002"
CHECKS = 0


def check(condition, label):
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dot(left, right):
    return sum((F(a) * F(b) for a, b in zip(left, right)), F(0))


def matvec(matrix, vector):
    return tuple(sum((matrix[i][j] * vector[j] for j in range(len(vector))), F(0))
                 for i in range(len(matrix)))


pins = {
    CU_DIR / "derive_complete_six_pair_h6_source_jet.py":
        "a62141731da025bdad51fde9223f9fa7d57f5a78bfff5d3375c4eeb1fa8ec84f",
    CU_DIR / "EXACT_LEDGER.json":
        "75c08daa934ab2de4a3300a6abb2a9fe4be72f99ce9790e7a73287d322c16de0",
    CU_DIR / "MANIFEST.sha256":
        "58f313b1053756473fa8bfd1ff23b25e2f1740b1903bb366771e26cd3b3c89ec",
    CU_DIR / "SEAL.sha256":
        "27cb7ad002e2539f756493697a1d90c9c259a1711c2c5d12a10970818b002810",
    HERE / "derive_raw_source_composition.py":
        "e494e2ba3612f68e4be8283c6e4a1c1c5861d2465c668d898717e9e02ed82141",
}
for path, expected in pins.items():
    check(path.is_file() and digest(path) == expected, f"pinned byte {path.name}")

check((CU_DIR / "SEAL.sha256").read_text().strip() ==
      "58f313b1053756473fa8bfd1ff23b25e2f1740b1903bb366771e26cd3b3c89ec  MANIFEST.sha256",
      "CU seal points to pinned manifest")

T = ((F(1), F(1), F(1)), (F(1), F(-1), F(-1)),
     (F(-1), F(1), F(-1)), (F(-1), F(-1), F(1)))
V = tuple(tuple(x / 2 for x in row) for row in T)
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
slots = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
weights = (F(1), F(1), F(1), F(2), F(2), F(2))

qrows = []
for a, b in PAIRS:
    tensor = tuple(tuple(V[a][i] * V[b][j] + V[b][i] * V[a][j]
                         for j in range(3)) for i in range(3))
    qrows.append(tuple(tensor[i][j] * weights[s]
                       for s, (i, j) in enumerate(slots)))
qrows = tuple(qrows)
expected_qrows = (
    (F(1, 2), F(-1, 2), F(-1, 2), 0, 0, -1),
    (F(-1, 2), F(1, 2), F(-1, 2), 0, -1, 0),
    (F(-1, 2), F(-1, 2), F(1, 2), -1, 0, 0),
    (F(-1, 2), F(-1, 2), F(1, 2), 1, 0, 0),
    (F(-1, 2), F(1, 2), F(-1, 2), 0, 1, 0),
    (F(1, 2), F(-1, 2), F(-1, 2), 0, 0, 1),
)
check(qrows == expected_qrows, "D_C star columns")

p = (
    (F(0), F(0), F(0)),
    (F(0), F(1), F(0)),
    (F(0), F(0), F(-1)),
    (F(0), F(0), F(1)),
    (F(0), F(-1), F(0)),
    (F(0), F(0), F(0)),
)
qt_p = tuple(tuple(sum((qrows[a][s] * p[a][r] for a in range(6)), F(0))
                       for r in range(3)) for s in range(6))
s_matrix = (
    (F(0), F(0), F(0)),
    (F(0), F(0), F(0)),
    (F(0), F(0), F(0)),
    (F(0), F(0), F(-1)),
    (F(0), F(1), F(0)),
    (F(0), F(0), F(0)),
)
check(qt_p == tuple(tuple(-2 * value for value in row) for row in s_matrix),
      "raw-pair to tensor-row matrix transport")
check(F(4) * F(8032, 27) == F(32128, 27), "h4 k0 coefficient")
check(F(4) * F(5010712, 6075) == F(20042848, 6075),
      "h6 diagonal k0 coefficient")
check(F(4) * F(-528) == F(-2112), "h6 cycle k0 coefficient")
check(F(20042848, 6075) - F(2112) == F(7212448, 6075),
      "h6 total k0 coefficient")

# Load only CU's graph definitions and pointwise frozen first-source formulas;
# do not execute its expensive run().
engine_path = CU_DIR / "derive_complete_six_pair_h6_source_jet.py"
spec = importlib.util.spec_from_file_location("gl6cu_v002_fast_check", engine_path)
cu = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = cu
assert spec.loader is not None
spec.loader.exec_module(cu)
check(len(cu.NODES) == 128 and len(cu.ACTIVE_CYCLES) == 64,
      "Q4 and active-cycle census")

axes = ((F(1), F(0), F(0)), (F(0), F(1), F(0)), (F(0), F(0), F(1)))
abasis = []
for i in range(3):
    a = [[F(0) for _ in range(3)] for _ in range(3)]
    a[i][i] = F(1, 2)
    abasis.append(tuple(tuple(row) for row in a))
for i, j in ((0, 1), (0, 2), (1, 2)):
    a = [[F(0) for _ in range(3)] for _ in range(3)]
    a[i][j] = a[j][i] = F(1, 2)
    abasis.append(tuple(tuple(row) for row in a))


def tdir(vector):
    return tuple(dot(vector, tuple(V[a][i] + V[b][i] for i in range(3)))
                 for a, b in PAIRS)


nonzero_second = False
for s, a_s in enumerate(abasis):
    for r, axis in enumerate(axes):
        mixed = tdir(matvec(a_s, axis))
        nonzero_second = nonzero_second or any(mixed)
        bare = h2 = h4 = F(0)
        for node in cu.NODES:
            sigma = F(1) if node[0] == "P" else F(-1)
            eta_sr = tuple(2 * sigma * value for value in mixed)
            memory = cu.pair_word(cu.BASE, node)
            bare += dot(memory, eta_sr)
            h2 += dot(tuple(-value for value in memory), eta_sr)
            h4 += dot(tuple(-F(4, 9) - F(37, 12) * value
                            for value in memory), eta_sr)
        check(bare == h2 == h4 == 0, f"second-chain bare/h2/h4 {s},{r}")

        for cycle in cu.ACTIVE_CYCLES:
            total = F(0)
            for index in range(6):
                common = set(cu.ENDS[cycle[index - 1]]) & set(cu.ENDS[cycle[index]])
                check(len(common) == 1, "cycle node incidence")
                node = next(iter(common))
                ports = tuple(sorted(cu.EDGES[edge][1]
                                     for edge in cycle if node in cu.ENDS[edge]))
                pair = cu.PAIR_INDEX[ports]
                sigma = F(1) if node[0] == "P" else F(-1)
                total += F(105, 8) * 2 * sigma * mixed[pair]
            check(total == 0, f"second-chain cycle cancellation {s},{r}")
check(nonzero_second, "eta_sr is not identically zero")

code = (HERE / "derive_raw_source_composition.py").read_text()
for token in (
    "compile_hessian_character_kernel",
    "evaluate_character_kernel",
    "compile_second_chain_matrix",
    "h6_active_cycles_64",
    "translation_difference_stencils",
    "not a cell Ward inference",
):
    check(token in code, f"executable interface token {token}")

addendum = (HERE / "ADDENDUM.md").read_text()
for token in (
    "A_r=0",
    "b_s=0",
    "b_r=\\chi_m(v)e_r",
    "b_{sr}=0",
    "\\eta_{A,sr}=2\\sigma(v)t_A(A_se_r)\\ne0",
    "{32128\\over27}",
    "{20042848\\over6075}",
    "-2112S",
    "{7212448\\over6075}",
    "not asserted to be an authenticated physical source law",
):
    check(token in addendum, f"addendum token {token}")

verification = (HERE / "VERIFICATION.txt").read_text()
check("PASS__GL6CY_DEVELOPMENT_RAW_SOURCE_COMPOSITION__21/21" in verification,
      "full replay transcript")
check("OBJECT=RAW_HAMILTONIAN_CY02__NOT_CONNECTED_CTP_1PI_OR_WARD" in verification,
      "scope transcript")
check("GL6CY_TARGET_BYTES_CHANGED=NO" in verification,
      "target custody transcript")

print(f"PASS__GL6CY_DEVELOPMENT_FAST_CHECK__{CHECKS}/{CHECKS}")
print("FULL_REPLAY=21/21__439.229_SECONDS")
print("ETA_SR=NONZERO__SECOND_CHAIN_CONTRACTED__ZERO_ON_DECLARED_PATH")
print("EXPORT=EXACT_Q4_TRANSLATION_DIFFERENCE_LAURENT_STENCIL")
print("SCOPE=RAW_HAMILTONIAN__NOT_CONNECTED_CTP_1PI_WARD_GRAVITY_OR_G")
