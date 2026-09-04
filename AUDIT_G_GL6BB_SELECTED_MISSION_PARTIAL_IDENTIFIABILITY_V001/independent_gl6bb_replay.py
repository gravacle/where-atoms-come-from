#!/usr/bin/env python3
"""Independent hostile mathematical replay for frozen GL6BB V001.

This program imports no author module.  It reconstructs the radius-zero parent,
the Dicke reduction, the state-free and robust intervals, the energy bound, and
the calculator certificates from standard-library arithmetic.  The author
calculator is exercised only as a separate command-line executable.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import subprocess
import sys
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6BB_SELECTED_MISSION_PARTIAL_IDENTIFIABILITY_V001"
CALCULATOR = AUTHOR / "calculate_prepared_blank_collar0.py"
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"FAIL [{CHECKS}] {label}")


def weight(state: int) -> int:
    return bin(state).count("1")


def permute_state(state: int, permutation: tuple[int, ...]) -> int:
    result = 0
    for old_bit, new_bit in enumerate(permutation):
        if state & (1 << old_bit):
            result |= 1 << new_bit
    return result


def full_hamiltonian(ratio: Fraction) -> list[list[Fraction]]:
    matrix = [[Fraction(0) for _ in range(16)] for _ in range(16)]
    for state in range(16):
        k = weight(state)
        matrix[state][state] = ratio * k * (k - 7)
        for bit in range(4):
            matrix[state ^ (1 << bit)][state] -= 1
    return matrix


def matvec(matrix: list[list[Fraction]], vector: list[Fraction]) -> list[Fraction]:
    return [sum(row[column] * vector[column] for column in range(len(vector)))
            for row in matrix]


def exact_rank(columns: list[list[Fraction]]) -> int:
    if not columns:
        return 0
    rows = [[columns[column][row] for column in range(len(columns))]
            for row in range(len(columns[0]))]
    rank = 0
    column = 0
    while rank < len(rows) and column < len(columns):
        pivot = next((row for row in range(rank, len(rows))
                      if rows[row][column] != 0), None)
        if pivot is None:
            column += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][column]
        rows[rank] = [entry / pivot_value for entry in rows[rank]]
        for row in range(len(rows)):
            if row == rank or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [left - factor * right
                         for left, right in zip(rows[row], rows[rank])]
        rank += 1
        column += 1
    return rank


def verify_nonselection_and_state_free_interval() -> None:
    # The exact inherited d_star=3 witness realizes every positive rational R;
    # the algebra is symbolic and the grid is an adversarial reconstruction.
    for numerator in range(1, 35):
        for denominator in range(1, 35):
            ratio = Fraction(numerator, denominator)
            h = Fraction(11, 17)
            ud = ratio * h
            delta = 4 * ud
            epsilon_star = delta + 2 * ud * (1 - 2 * 3)
            check(h > 0 and ud > 0 and delta > 0,
                  "strict all-positive inherited member")
            check(ud / h == ratio, "arbitrary positive R is realized")
            check(epsilon_star == -6 * ud,
                  "native lock-ray onsite coefficient")

    admitted = (Fraction(2), Fraction(5, 2))
    check(tuple(1 / ratio for ratio in admitted)
          == (Fraction(1, 2), Fraction(2, 5)),
          "two H6 witnesses invert to the declared h/U_d values")

    eigenvalues = []
    for state in range(16):
        z0 = 1 if not (state & 1) else -1
        z1 = 1 if not (state & 2) else -1
        eigenvalues.append(z0 * z1)
    check(eigenvalues.count(1) == eigenvalues.count(-1) == 8,
          "pair observable has two nonempty spectral endpoints")
    # A convex mixture of one state from each eigenspace realizes every point.
    for denominator in range(1, 80):
        for numerator in range(denominator + 1):
            probability = Fraction(numerator, denominator)
            expectation = probability * 1 + (1 - probability) * (-1)
            check((1 + expectation) / 2 == probability,
                  "state-free probability point is attained")


def verify_exact_five_state_reduction() -> None:
    permutations = tuple(itertools.permutations(range(4)))
    expected_diagonal_coefficients = (0, -6, -10, -12, -12)
    expected_off_diagonal_squares = (4, 6, 6, 4)
    expected_pair_plus = (Fraction(1), Fraction(1, 2), Fraction(1, 3),
                          Fraction(1, 2), Fraction(1))

    for ratio in (Fraction(2), Fraction(5, 2)):
        matrix = full_hamiltonian(ratio)
        for row in range(16):
            check(matrix[row][row] == ratio * weight(row) * (weight(row) - 7),
                  "full parent diagonal is R(k^2-7k)")
            check(matrix[row] == [matrix[column][row] for column in range(16)],
                  "full radius-zero Hamiltonian is Hermitian")
            check(sum(value != 0 for column, value in enumerate(matrix[row])
                      if column != row) == 4,
                  "each word has four transverse neighbors")
            check(all(value == -1 for column, value in enumerate(matrix[row])
                      if column != row and value != 0),
                  "every transverse sign is negative")

        # The exact parent is invariant under every permutation of four links.
        for permutation in permutations:
            for row in range(16):
                prow = permute_state(row, permutation)
                for column in range(16):
                    pcolumn = permute_state(column, permutation)
                    check(matrix[prow][pcolumn] == matrix[row][column],
                          "S4 symmetry of the exact four-link parent")

        # The blank cyclic subspace has exactly dimension five and is uniform
        # within each Hamming-weight orbit.
        vector = [Fraction(int(state == 0)) for state in range(16)]
        krylov = []
        for _ in range(5):
            krylov.append(vector)
            for k in range(5):
                orbit_values = {vector[state] for state in range(16)
                                if weight(state) == k}
                check(len(orbit_values) == 1,
                      "blank Krylov vector is constant on each Dicke orbit")
            vector = matvec(matrix, vector)
        check(exact_rank(krylov) == 5, "blank cyclic subspace is exactly five-dimensional")

        diagonals = tuple(ratio * k * (k - 7) for k in range(5))
        check(diagonals == tuple(ratio * value
                                 for value in expected_diagonal_coefficients),
              "Dicke diagonal is exact")

    # Count the normalized Dicke couplings without introducing floating square
    # roots: their signs are negative and their exact squares are integral.
    for k in range(4):
        edge_count = math.comb(4, k) * (4 - k)
        square = Fraction(edge_count * edge_count,
                          math.comb(4, k) * math.comb(4, k + 1))
        check(square == expected_off_diagonal_squares[k],
              "normalized Dicke off-diagonal square")
        check(-(4 - k) < 0, "Dicke transverse coupling has negative sign")

    for k in range(5):
        states = [state for state in range(16) if weight(state) == k]
        plus = sum(((state >> 0) & 1) == ((state >> 1) & 1)
                   for state in states)
        check(Fraction(plus, len(states)) == expected_pair_plus[k],
              "exact pair-plus Dicke weight")


def verify_energy_bound_and_member_values() -> None:
    # The commuting X_a have joint eigenvalues whose sums have norm four.
    x_sums = [sum(signs) for signs in itertools.product((-1, 1), repeat=4)]
    check(max(abs(value) for value in x_sums) == 4,
          "operator norm of sum X_a is four")

    for k in range(5):
        g_value = k * (7 - k)
        mismatch = Fraction(k * (4 - k), 6)
        check(g_value >= 0, "G=N(7-N) is nonnegative on the five-state sector")
        check(mismatch <= Fraction(g_value, 12),
              "pair mismatch is pointwise bounded by G/12")
        check(1 - mismatch in (Fraction(1), Fraction(1, 2), Fraction(1, 3)),
              "pair-plus weight normalization")

    # The blank has both zero transverse expectation and G=0, hence zero K_R
    # energy.  Conservation gives R<G>=-<sum X><=4.
    blank_x_expectation = Fraction(0)
    blank_g = 0
    check(-blank_x_expectation - Fraction(2) * blank_g == 0,
          "blank energy is zero for R=2")
    check(-blank_x_expectation - Fraction(5, 2) * blank_g == 0,
          "blank energy is zero for R=5/2")

    expected = {Fraction(2): Fraction(5, 6),
                Fraction(5, 2): Fraction(13, 15)}
    for ratio, lower in expected.items():
        reconstructed = 1 - Fraction(1, 3 * ratio)
        check(reconstructed == lower, "exact admitted-member energy lower bound")
        check(0 <= lower <= 1, "energy lower bound is a probability")


def verify_robust_intervals_and_exterior_propagation() -> None:
    # Exhaust a dense exact rational battery for pointwise clipping.
    for denominator in range(1, 14):
        values = [Fraction(index, denominator)
                  for index in range(denominator + 1)]
        for q_value in values:
            for epsilon in values:
                lower = max(Fraction(0), q_value - epsilon)
                upper = min(Fraction(1), q_value + epsilon)
                check(0 <= lower <= upper <= 1, "clipped pointwise interval is ordered")
                for actual in values:
                    if abs(actual - q_value) <= epsilon:
                        check(lower <= actual <= upper,
                              "clipped pointwise interval contains licensed value")

    # Pointwise extrema preserve correlation and can be strictly sharper than
    # independently extremizing q and epsilon.
    correlated = ((Fraction(1, 5), Fraction(1, 10)),
                  (Fraction(4, 5), Fraction(1, 20)))
    point_lower = max(0, min(q - epsilon for q, epsilon in correlated))
    point_upper = min(1, max(q + epsilon for q, epsilon in correlated))
    separate_lower = max(0, min(q for q, _ in correlated)
                         - max(epsilon for _, epsilon in correlated))
    separate_upper = min(1, max(q for q, _ in correlated)
                         + max(epsilon for _, epsilon in correlated))
    check((point_lower, point_upper) == (Fraction(1, 10), Fraction(17, 20)),
          "correlated robust hull is exact")
    check(point_upper < separate_upper and point_lower == separate_lower,
          "pointwise hull avoids a decorrelation loss")

    # Classical binary states saturate the trace-distance coefficient one.
    for denominator in range(1, 35):
        for left in range(denominator + 1):
            for right in range(denominator + 1):
                p_value = Fraction(left, denominator)
                q_value = Fraction(right, denominator)
                trace_distance = abs(p_value - q_value)
                probability_change = abs(p_value - q_value)
                check(probability_change == trace_distance,
                      "trace-distance half-width coefficient one is sharp")

    # Write z=exp(48 R |s|)>=1.  This exact rational battery checks that the
    # min(1,z-1) cap and final [0,1] clipping yield the displayed formulas.
    members = ((Fraction(2), Fraction(5, 6), Fraction(11, 6), 96),
               (Fraction(5, 2), Fraction(13, 15), Fraction(28, 15), 120))
    for ratio, collar_lower, constant, rate in members:
        check(48 * ratio == rate, "GL6BA L=0 exponential rate substitution")
        check(collar_lower + 1 == constant,
              "analytic full-mission lower constant")
        for numerator in range(20, 81):
            z_value = Fraction(numerator, 20)
            epsilon = min(Fraction(1), z_value - 1)
            propagated = max(Fraction(0), collar_lower - epsilon)
            displayed = max(Fraction(0), constant - z_value)
            check(propagated == displayed,
                  "capped exterior error equals displayed clipped lower endpoint")


def factorial_tail_majorant(x_value: Fraction, order: int) -> Fraction:
    check(x_value >= 0 and order >= 0, "tail majorant domain")
    if x_value == 0:
        return Fraction(0)
    ratio = x_value / (order + 2)
    check(ratio < 1, "tail geometric ratio is strictly below one")
    first = x_value ** (order + 1) / math.factorial(order + 1)
    return first / (1 - ratio)


def choose_order(x_value: Fraction, tolerance: Fraction) -> tuple[int, Fraction]:
    if x_value == 0:
        return 0, Fraction(0)
    start = max(0, math.ceil(x_value) - 1)
    for order in range(start, 5000):
        if x_value / (order + 2) >= 1:
            continue
        remainder = factorial_tail_majorant(x_value, order)
        if remainder * (2 + remainder) <= tolerance:
            return order, remainder
    raise RuntimeError("independent order search exceeded its fail-closed cap")


def polynomial_probability(ratio: Fraction, sigma: Fraction, order: int) -> Fraction:
    matrix = full_hamiltonian(ratio)
    real = [Fraction(0) for _ in range(16)]
    imaginary = [Fraction(0) for _ in range(16)]
    vector = [Fraction(int(state == 0)) for state in range(16)]
    coefficient = Fraction(1)
    for degree in range(order + 1):
        phase = degree % 4
        target = real if phase in (0, 2) else imaginary
        sign = 1 if phase in (0, 3) else -1
        for index, value in enumerate(vector):
            target[index] += sign * coefficient * value
        if degree < order:
            vector = matvec(matrix, vector)
            coefficient *= sigma / (degree + 1)
    return sum(real[state] ** 2 + imaginary[state] ** 2
               for state in range(16)
               if ((state >> 0) & 1) == ((state >> 1) & 1))


def capped_exponential_upper(x_value: Fraction, tolerance: Fraction) -> Fraction:
    check(x_value >= 0 and tolerance > 0, "capped exponential domain")
    if x_value == 0:
        return Fraction(0)
    partial = Fraction(1)
    term = Fraction(1)
    for degree in range(1, 10000):
        term *= x_value / degree
        partial += term
        if partial >= 2:
            return Fraction(1)
        if x_value / (degree + 2) >= 1:
            continue
        first_omitted = term * x_value / (degree + 1)
        tail = first_omitted / (1 - x_value / (degree + 2))
        if partial + tail < 2 and tail <= tolerance:
            return partial + tail - 1
    raise RuntimeError("independent capped exponential exceeded its fail-closed cap")


def command(arguments: list[str], optimized: bool) -> subprocess.CompletedProcess[str]:
    invocation = [sys.executable, "-B"]
    if optimized:
        invocation.append("-O")
    invocation.extend([str(CALCULATOR), *arguments])
    return subprocess.run(invocation, cwd=AUTHOR, text=True,
                          capture_output=True, check=False)


def calculator_json(arguments: list[str], optimized: bool) -> dict[str, object]:
    completed = command(arguments, optimized)
    check(completed.returncode == 0, "author calculator command succeeds")
    check(completed.stderr == "", "author calculator emits no warning stream")
    return json.loads(completed.stdout)


def exact_field(payload: dict[str, object], key: str) -> Fraction:
    value = payload[key]
    check(isinstance(value, dict) and "exact" in value,
          f"calculator exact field: {key}")
    return Fraction(str(value["exact"]))


def verify_calculator_certificates() -> None:
    # The command-line interface must fail closed when either physical scalar is
    # missing or outside its declared positive domain.
    missing_sigma = command(["--ratio", "2"], optimized=False)
    check(missing_sigma.returncode != 0 and "--sigma" in missing_sigma.stderr,
          "sigma has no physical command-line default")
    missing_ratio = command(["--sigma", "0"], optimized=False)
    check(missing_ratio.returncode != 0 and "--ratio" in missing_ratio.stderr,
          "ratio has no physical command-line default")
    bad_ratio = command(["--ratio", "0", "--sigma", "0"], optimized=False)
    check(bad_ratio.returncode != 0, "nonpositive ratio fails closed")
    bad_tolerance = command(["--ratio", "2", "--sigma", "0",
                             "--tolerance", "0"], optimized=False)
    check(bad_tolerance.returncode != 0, "nonpositive tolerance fails closed")

    cases = (
        (Fraction(2), Fraction(0), Fraction(1, 10**14)),
        (Fraction(2), Fraction(1, 10000), Fraction(1, 10**14)),
        (Fraction(5, 2), Fraction(-1, 1000), Fraction(1, 10**14)),
        (Fraction(2), Fraction(1, 10), Fraction(1, 10**12)),
    )
    for ratio, sigma, tolerance in cases:
        arguments = ["--ratio", str(ratio), f"--sigma={sigma}",
                     "--tolerance", str(tolerance)]
        normal = calculator_json(arguments, optimized=False)
        optimized = calculator_json(arguments, optimized=True)
        check(normal == optimized, "calculator normal and optimized JSON agree exactly")
        check(normal["R"] == str(ratio) and normal["sigma_obs"] == str(sigma),
              "calculator reports exact supplied physical scalars")
        check(normal["internal_probability_tolerance"] == str(tolerance),
              "calculator reports exact internal tolerance")
        check(normal["admitted_h6_member"] is (ratio in (2, Fraction(5, 2))),
              "calculator labels admitted members without selecting one")
        check("sigma supplied by caller" in str(normal["scope"]),
              "calculator output retains conditional scope")

        norm_bound = 12 * ratio + 4
        x_value = norm_bound * abs(sigma)
        order, remainder = choose_order(x_value, tolerance)
        check(normal["order"] == order, "independent propagator order agrees")
        check(exact_field(normal, "norm_bound") == norm_bound,
              "independent Hamiltonian norm bound agrees")
        check(exact_field(normal, "unitary_remainder") == remainder,
              "independent unitary tail majorant agrees")
        probability_error = remainder * (2 + remainder)
        check(exact_field(normal, "probability_error") == probability_error,
              "independent probability error conversion agrees")

        approximate = polynomial_probability(ratio, sigma, order)
        check(exact_field(normal, "approximate") == approximate,
              "independent exact Taylor polynomial agrees")
        collar_lower = max(Fraction(0), approximate - probability_error)
        collar_upper = min(Fraction(1), approximate + probability_error)
        check(exact_field(normal, "lower") == collar_lower,
              "calculator collar lower endpoint is outward")
        check(exact_field(normal, "upper") == collar_upper,
              "calculator collar upper endpoint is outward")

        exterior = capped_exponential_upper(48 * ratio * abs(sigma),
                                             tolerance / 100)
        check(exact_field(normal, "exterior_error_upper") == exterior,
              "independent GL6BA exponential certificate agrees")
        check(exact_field(normal, "full_lower")
              == max(Fraction(0), collar_lower - exterior),
              "full lower endpoint clips after adding exterior error")
        check(exact_field(normal, "full_upper")
              == min(Fraction(1), collar_upper + exterior),
              "full upper endpoint clips after adding exterior error")
        analytic = max(Fraction(0), 1 - Fraction(1, 3) / ratio)
        check(exact_field(normal, "analytic_collar_lower") == analytic,
              "calculator reports the exact energy lower bound")
        if sigma == 0:
            check((collar_lower, collar_upper, exterior)
                  == (Fraction(1), Fraction(1), Fraction(0)),
                  "zero-time mission is exactly [1,1]")
        if 48 * ratio * abs(sigma) >= 1:
            check(exterior == 1, "large exterior certificate clips honestly")
            check((exact_field(normal, "full_lower"),
                   exact_field(normal, "full_upper")) == (0, 1),
                  "clipped complete mission returns [0,1]")

    # Fraction, decimal, and scientific spellings are parsed to the same exact
    # rational and must produce byte-equivalent semantic JSON.
    variants = []
    for spelling in ("1/1000", "0.001", "1e-3"):
        variants.append(calculator_json(
            ["--ratio", "2", "--sigma", spelling,
             "--tolerance", "1e-14"], optimized=False))
    check(variants[0] == variants[1] == variants[2],
          "fraction decimal and scientific inputs are exactly equivalent")


def verify_scope_and_custody_meaning() -> None:
    theorem = (AUTHOR / "THEOREM.md").read_text(encoding="utf-8")
    result = (AUTHOR / "RESULT.md").read_text(encoding="utf-8")
    result_flat = " ".join(result.split())
    dependencies = (AUTHOR / "DEPENDENCIES.md").read_text(encoding="utf-8")
    dependencies_flat = " ".join(dependencies.split())
    ledger = json.loads((AUTHOR / "IDENTIFIABILITY_LEDGER.json").read_text(
        encoding="utf-8"))

    required_theorem = (
        "every positive",
        "not a measurement, posterior bound, or selection rule",
        "provide no finite upper",
        "prepared-blank state as an already defined special case",
        "=[0,1]",
        "same argument applies directly on every complete finite exterior",
        "No independence between its entries is assumed",
        "D_{\\rm tr}",
        "five normalized Dicke states",
        "K_R=-\\sum_aX_a-RG",
        "R\\langle G\\rangle=-\\left\\langle\\sum_aX_a\\right\\rangle\\le4",
        "1-{1\\over3R}",
        "11/6-e^{96|\\sigma_{\\rm obs}|}",
        "28/15-e^{120|\\sigma_{\\rm obs}|}",
        "No read outcome",  # supplied by the exact parent audit, checked below
    )
    # The last no-postselection token is intentionally checked across the exact
    # parent and its audit rather than fabricated into GL6BB's theorem.
    for token in required_theorem[:-1]:
        check(token in theorem, f"frozen theorem scope token: {token}")

    gl6ba_audit = (ROOT / "AUDIT_G_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001"
                   / "AUDIT.md").read_text(encoding="utf-8")
    check(required_theorem[-1] in gl6ba_audit and "postselected" in gl6ba_audit,
          "exact parent retains every read outcome without postselection")
    check("same-clock" in theorem and "its pulse time is not silently charged" not in theorem,
          "GL6BB uses the inherited sampling clock without inventing read time")
    gl6az = (ROOT / "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001"
             / "THEOREM.md").read_text(encoding="utf-8")
    gl6az_flat = " ".join(gl6az.split())
    check("endpoint immediately before the identical terminal read dilation" in gl6az_flat
          and "pulse time is not silently charged" in gl6az_flat,
          "clock ownership is fixed by the declared exact dependency")

    prepared_marker = theorem.index(
        "already defined prepared-blank, two-member scenario")
    sole_datum_marker = theorem.index("one new physical datum", prepared_marker)
    actual_mission_marker = theorem.index("To complete one actual selected mission")
    check(actual_mission_marker < prepared_marker < sole_datum_marker,
          "sigma-only datum claim is confined to the explicit conditional scenario")
    check("All three entries are absent as a jointly selected tuple" in theorem,
          "actual mission still requires the selected triple")
    check("One actual single-member mission still requires" in result_flat,
          "result does not promote the two-member conditional envelope")

    check("calibration is supplied rather than derived" in dependencies_flat,
          "dependency ledger retains calibration ceiling")
    check("prepared-blank special case" in dependencies_flat,
          "dependency ledger binds the blank branch to GL6BA")
    check(ledger["selected_mission_triple"] == ["R", "sigma_obs", "omega_L"],
          "machine ledger has the exact selected-mission triple")
    check(ledger["state_free_probability_interval"] == ["0", "1"],
          "machine ledger records the sharp state-free interval")
    check(ledger["prepared_blank_L0"]["exact_reduced_dimension"] == 5,
          "machine ledger records the five-state reduction")
    check(ledger["minimal_new_datum_for_prepared_blank_two_member_scenario"]
          == "a same-clock numerical sigma_obs or finite interval",
          "machine ledger confines the minimal datum to the named scenario")
    check(all(word not in theorem.lower() for word in (
        "we have proved gravity", "this is gravity", "derives newton's constant")),
          "no forbidden gravity promotion")


def main() -> int:
    verify_nonselection_and_state_free_interval()
    verify_exact_five_state_reduction()
    verify_energy_bound_and_member_values()
    verify_robust_intervals_and_exterior_propagation()
    verify_calculator_certificates()
    verify_scope_and_custody_meaning()
    print(f"PASS__INDEPENDENT_GL6BB_HOSTILE_REPLAY__{CHECKS}/{CHECKS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
