#!/usr/bin/env python3
"""Exact replay for EU Alpha-to-q4 dimension lock.

This verifier proves the finite projector and engineering-dimension algebra,
checks the declared counterexamples and boundary tokens, and authenticates the
mutable packet.  It does not certify the physical QFRONT-DIM, same-front
Maxwell, or MARGINAL-ALPHA premises.
"""

from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
SCIENCE_ONLY = "--science-only" in sys.argv[1:]
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def eye(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def mm(left, right):
    right_t = transpose(right)
    return [
        [sum((x * y for x, y in zip(row, col)), F(0)) for col in right_t]
        for row in left
    ]


def mv(matrix, vector):
    return [
        sum((x * y for x, y in zip(row, vector)), F(0))
        for row in matrix
    ]


def rank(matrix):
    work = [row[:] for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    pivot_row = 0
    for col in range(cols):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][col]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][col]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row:
                continue
            value = work[row][col]
            if value:
                work[row] = [
                    entry - value * pivot_entry
                    for entry, pivot_entry in zip(work[row], work[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


# EU02--EU05: exact general-q count/contrast algebra.
for q in range(2, 13):
    identity = eye(q)
    projector = [
        [identity[i][j] - F(1, q) for j in range(q)]
        for i in range(q)
    ]
    one = [F(1)] * q
    check(mm(projector, projector) == projector, f"P_{q} idempotent")
    check(transpose(projector) == projector, f"P_{q} symmetric")
    check(mv(projector, one) == [F(0)] * q, f"P_{q} kills one")
    check(rank(projector) == q - 1, f"rank P_{q}=q-1")
    check(sum(projector[i][i] for i in range(q)) == q - 1,
          f"trace P_{q}=q-1")
    spatial = rank(projector)
    spacetime = 1 + spatial
    check(spacetime == q, f"QFRONT-DIM gives D=q for q={q}")


# EU06--EU14: exact engineering dimensions in general D.
marginal_dimensions = []
for dimension in range(2, 13):
    a_dim = F(dimension - 2, 2)
    psi_dim = F(dimension - 1, 2)
    e_dim = F(4 - dimension, 2)
    alpha_dim = F(4 - dimension)
    g_squared_dim = F(4 - dimension)

    check(2 + 2 * a_dim == dimension,
          f"Maxwell kinetic dimension D={dimension}")
    check(1 + 2 * psi_dim == dimension,
          f"Dirac kinetic dimension D={dimension}")
    check(e_dim + a_dim == 1,
          f"covariant derivative dimension D={dimension}")
    check(e_dim + a_dim + 2 * psi_dim == dimension,
          f"charged vertex dimension D={dimension}")
    check(alpha_dim == 2 * e_dim,
          f"alpha engineering dimension D={dimension}")
    check(g_squared_dim == alpha_dim,
          f"alternative gauge normalization D={dimension}")
    check(alpha_dim + F(dimension - 4) == 0,
          f"scale-dressed alpha is dimensionless D={dimension}")
    check(F(dimension - 4, 2) == -e_dim,
          f"canonical e beta coefficient D={dimension}")
    if e_dim == 0:
        marginal_dimensions.append(dimension)

check(marginal_dimensions == [4], "classical Maxwell marginality selects D=4")


# EU12: composition with QFRONT-DIM selects q=4 and only q=4.
admitted_q = [q for q in range(2, 13) if F(4 - q, 2) == 0]
check(admitted_q == [4], "D=q plus marginal Maxwell charge selects q=4")


# Explicit controls: these evade the conclusion only by failing a premise.
for dimension in (2, 3, 5, 6, 7, 8):
    e_dim = F(4 - dimension, 2)
    alpha_dim = 2 * e_dim
    check(e_dim != 0, f"D={dimension} charge is not classically marginal")
    check(alpha_dim + F(dimension - 4) == 0,
          f"D={dimension} admits explicit mu-dressed dimensionless alpha")
check(F(4 - 5, 2) + F(1, 2) == 0,
      "D=5 scale compensator is load-bearing")
check(F(4 - 3, 2) - F(1, 2) == 0,
      "D=3 inverse scale compensator is load-bearing")


if not SCIENCE_ONLY:
    required = {
        "THEOREM.md": (
            "`QFRONT-DIM`",
            "`MARGINAL-ALPHA`",
            "\\boxed{D=1+(q-1)=q.}",
            "\\boxed{[e_D]=1-[A_\\mu]={4-D\\over2}}",
            "\\widehat\\alpha(\\mu)=C_D e_D^2\\mu^{D-4}",
            "An interacting fixed point can cancel the canonical term",
            "same record front earns `QFRONT-DIM`",
            "It neither produces the four\nphysical operations nor calculates",
            "Relabel the q=4 selector as a metric, curvature,",
        ),
        "README.md": (
            "This is a narrow REQUIRE-side selector.",
            "The same-front requirement is load-bearing.",
            "This packet is mutable and requests independent hostile prescreen.",
        ),
        "RESULT.md": (
            "MUTABLE_PRESCREEN_READY",
            "ALPHA_MARGINALITY_SELECTS_D4_AND_Q4_CONDITIONALLY",
        ),
        "AUDIT.md": (
            "not independent\nhostile review",
            "cannot freeze or seal the lane.",
        ),
        "PRESCREEN_REQUEST.md": (
            "REQUESTED__MUTABLE_SOURCES__NO_FREEZE_AUTHORIZED",
            "Do not edit builder files.",
        ),
    }
    for name, tokens in required.items():
        text = (HERE / name).read_text()
        for token in tokens:
            check(token in text, f"{name} required token: {token}")

    dependency_lines = [
        line for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
        if line.strip()
    ]
    check(len(dependency_lines) == 6, "dependency ledger has six pins")
    for line in dependency_lines:
        expected, relative = line.split(None, 1)
        target = (HERE / relative.strip()).resolve()
        check(target.is_file(), f"dependency exists: {relative.strip()}")
        check(digest(target) == expected,
              f"dependency hash: {relative.strip()}")

    manifest_path = HERE / "MANIFEST.sha256"
    if manifest_path.exists():
        manifest_lines = [
            line for line in manifest_path.read_text().splitlines()
            if line.strip()
        ]
        check(len(manifest_lines) == 8, "builder manifest has eight entries")
        expected_names = {
            "README.md", "THEOREM.md", "RESULT.md", "AUDIT.md",
            "PRESCREEN_REQUEST.md", "DEPENDENCIES.sha256",
            "VERIFICATION.txt", "verify_alpha_q4_dimension_lock.py",
        }
        actual_names = set()
        for line in manifest_lines:
            expected, name = line.split(None, 1)
            name = name.strip()
            actual_names.add(name)
            check(digest(HERE / name) == expected, f"manifest hash: {name}")
        check(actual_names == expected_names, "builder manifest exact file census")

    for path in sorted(HERE.iterdir()):
        if not path.is_file():
            continue
        data = path.read_bytes()
        check(b"\x00" not in data, f"no NUL bytes: {path.name}")
        check(b"\r" not in data, f"no CR bytes: {path.name}")
        check(data.endswith(b"\n"), f"terminal newline: {path.name}")


print("GENERAL_Q_COUNT_CONTRAST_D_EQUALS_Q_EXACT_CONDITIONAL")
print("CANONICAL_MAXWELL_ENGINEERING_DIMENSIONS_EXACT")
print("MARGINAL_ALPHA_REQUIRES_D4_AND_Q4_EXACT_CONDITIONAL")
print("D_NE_4_SCALE_DRESSED_AND_FIXED_POINT_BOUNDARY_EXPLICIT")
print(f"PASS {checks}/{checks}")
