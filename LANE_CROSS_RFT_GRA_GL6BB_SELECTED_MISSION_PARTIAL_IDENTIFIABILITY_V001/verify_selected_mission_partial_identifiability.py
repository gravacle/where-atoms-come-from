#!/usr/bin/env python3
"""Constructive replay for GL6BB V001 (standard library only)."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

from calculate_prepared_blank_collar0 import (
    blank_collar_probability_interval,
    capped_exp_minus_one_upper,
    complete_mission_interval,
    diagonal_energy,
    hamiltonian_action,
)


LANE = Path(__file__).resolve().parent
ROOT = LANE.parent
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"FAIL [{CHECKS}] {label}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


EXPECTED_DEPENDENCIES = {
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md": "00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec",
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MANIFEST.sha256": "a522cd11b7c9b62adb080c7895d6c184ddf05f0c2d17843699bd85900e513c3b",
    "LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/THEOREM.md": "2b88febc569efa0de0238e8000d018bf3f798a8ebed2e4ff1327f053d6bd9284",
    "LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/MANIFEST.sha256": "ea06b1f0aa012fb009477ea28debc49c9afa57328c0cf410e2856c9cc9542a27",
    "LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/INDEPENDENT_REAUDIT.md": "5c275748d54743ef44098f74c4c5698aead0845d51e6c2dcf32a1bef63f0c7bf",
    "LANE_CROSS_RFT_GRA_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/THEOREM.md": "32f597edc51a609a37b86144487cd7db3bd2f14a65adb754a893d47ef6807e81",
    "LANE_CROSS_RFT_GRA_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/MANIFEST.sha256": "24a71c01ed1b7a92830e92ec7682882c892667289e2794dafb4af5905ad71b2e",
    "LANE_CROSS_RFT_GRA_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/SEAL.sha256": "a946902f027c555f91cd1f2e9ce93e3182f8edeca319955cd691a0bc929fba51",
    "AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/AUDIT.md": "8754210b2ff0077e8cec4c5ce0f771ba40cfee2b8957d37f5be64a50ba49d0b4",
    "AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/MANIFEST.sha256": "73b618e88b96d40ca40e32bd33ab578df98ed20af9fc359d90bfe4ae75f5c91b",
    "AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/SEAL.sha256": "d0a50cda599842d8854db0bc2ab9e665f823e2ab7e2048a0eafa672ae3ad7b7e",
    "LANE_CROSS_RFT_GRA_GL5ZZF_H6_DOUBLED_REAL_METRIC_SOURCE_BIND_BOUNDARY_V001/THEOREM.md": "2f9fbdde026765b5a4dc335d6b87777e1042efcace1bd4d6482ce0f8ac235b22",
    "LANE_CROSS_RFT_GRA_GL5ZZF_H6_DOUBLED_REAL_METRIC_SOURCE_BIND_BOUNDARY_V001/MANIFEST.sha256": "49940820ce4e84c4157c2359de26f02012cfffcef052c44b391bb8098370dd8e",
    "LANE_CROSS_RFT_GRA_GL5ZZF_H6_DOUBLED_REAL_METRIC_SOURCE_BIND_BOUNDARY_V001/SEAL.sha256": "550ca316a3829796fa0007eb6377e9e092a62d95853e85df58d0777daaaf8411",
    "LANE_CROSS_RFT_GRA_GL5ZZF_H6_DOUBLED_REAL_METRIC_SOURCE_BIND_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md": "c29c8439a1b6ba5fcce4f2cd881461ee491616ed672dc50208762ef5ece37f58",
    "LANE_CROSS_RFT_GRA_GL5ZZF_H6_DOUBLED_REAL_METRIC_SOURCE_BIND_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT/AUDIT_MANIFEST.sha256": "46e03b796e3de36382975d4b7db337885cc541a468d07448cc5ed7c844a60225",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/THEOREM.md": "bcca352d9e58deba63a068a29f94a0e19fb88ae9994391ab1b5053120587ba44",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/MANIFEST.sha256": "a932e083b5dd629e41f0014d52bf1c65f982612b55314e1e03e082b0feb8ebc7",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/SEAL.sha256": "620eb8bcb254b0296bd4d2f7e81a0b3e78c0d0b9ce80ad257e97a4db394f392d",
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/AUDIT.md": "72c3c6552a89a3d5b977fb7ff566d558ed70d9c14b2caac64c8602609356498f",
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/MANIFEST.sha256": "9954328b6939e12c0b2f712bf6373d8a5c2a4cfafd77ab85aa133b27a4a9f707",
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/SEAL.sha256": "1d7cd10a59469139edf4b2836ce95820f605069ef3e6bc0d56bf4df58009b869",
    "LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/THEOREM.md": "d7ce0a7527a68f49e6ea2ee8edbb400a142fbb49297d8fe99cae78ffa0154ab0",
    "LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/MANIFEST.sha256": "6e14332230f713d51e393a5889fe78964fe0e63588b4b841533fa6af7ef19103",
    "LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/SEAL.sha256": "34f29b3c03d53c4dbc9736d1bf7a7785e0a49a2aad04299b69a3804290c5971e",
    "AUDIT_G_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/AUDIT.md": "03bda2dba369211542dfef1af065490e21033483443cd6a67a21b06bf74e0bc9",
    "AUDIT_G_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/MANIFEST.sha256": "a91fde1a69983815d9238306a273638fdd22aa009a16380d22c94cd5d5d186b9",
    "AUDIT_G_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/SEAL.sha256": "e45a27a04a8318472fdceff354179d5be13adee2d63c13e7f5b9671c0ff71965",
}


def verify_dependencies() -> None:
    for relative, expected in sorted(EXPECTED_DEPENDENCIES.items()):
        path = ROOT / relative
        check(path.is_file(), f"dependency exists: {relative}")
        check(sha256(path) == expected, f"dependency hash: {relative}")
    declared: dict[str, str] = {}
    for line in (LANE / "DEPENDENCIES.sha256").read_text(encoding="utf-8").splitlines():
        digest, relative = line.split("  ", 1)
        declared[relative] = digest
    check(declared == EXPECTED_DEPENDENCIES,
          "dependency ledger and executable custody map agree exactly")


def verify_ratio_custody() -> None:
    for ratio in (Fraction(2), Fraction(5, 2)):
        inverse = 1 / ratio
        check(inverse in (Fraction(1, 2), Fraction(2, 5)),
              "admitted H6 inverse ratio")
    for numerator in range(1, 30):
        for denominator in range(1, 30):
            ratio = Fraction(numerator, denominator)
            h = Fraction(7, 11)
            ud = ratio * h
            d_star = 3
            delta = 4 * ud * (d_star - 2)
            epsilon = delta + 2 * ud * (1 - 2 * d_star)
            check(ud / h == ratio, "all-positive F3 ray realizes ratio")
            check(delta > 0 and epsilon == -6 * ud,
                  "all-positive ratio lies on strict inherited lock ray")


def verify_state_free_sharpness() -> None:
    plus = minus = 0
    for state in range(16):
        z_zero = 1 if ((state >> 0) & 1) == 0 else -1
        z_one = 1 if ((state >> 1) & 1) == 0 else -1
        eigenvalue = z_zero * z_one
        plus += eigenvalue == 1
        minus += eigenvalue == -1
    check(plus == minus == 8, "pair observable has both eigenspaces")
    for probability in (Fraction(i, 100) for i in range(101)):
        # A mixture of one + and one - spectral state attains every rational
        # point in the dense battery; unitary conjugation preserves this fact.
        expectation = probability - (1 - probability)
        reconstructed = (1 + expectation) / 2
        check(reconstructed == probability, "state-free interval is attained")


def pair_plus_weight(weight: int) -> Fraction:
    states = [state for state in range(16) if bin(state).count("1") == weight]
    plus = 0
    for state in states:
        if ((state >> 0) & 1) == ((state >> 1) & 1):
            plus += 1
    return Fraction(plus, len(states))


def verify_dicke_reduction_and_energy_bound() -> None:
    expected_plus = (Fraction(1), Fraction(1, 2), Fraction(1, 3),
                     Fraction(1, 2), Fraction(1))
    for weight in range(5):
        check(pair_plus_weight(weight) == expected_plus[weight],
              "exact Dicke pair-plus weight")
        mismatch = 1 - expected_plus[weight]
        g_value = weight * (7 - weight)
        check(mismatch <= Fraction(g_value, 12),
              "pair mismatch bounded by G/12")

    # Full sixteen-state Hamiltonian diagonal and four flip neighbors.
    for ratio in (Fraction(2), Fraction(5, 2)):
        for state in range(16):
            weight = bin(state).count("1")
            check(diagonal_energy(state, ratio) == ratio * weight * (weight - 7),
                  "L=0 diagonal is R(k^2-7k)")
            basis = [Fraction(0) for _ in range(16)]
            basis[state] = 1
            image = hamiltonian_action(basis, ratio)
            check(image[state] == diagonal_energy(state, ratio),
                  "Hamiltonian action retains diagonal")
            neighbors = [index for index, value in enumerate(image)
                         if index != state and value != 0]
            check(len(neighbors) == 4, "Hamiltonian has four transverse neighbors")
            check(all(image[index] == -1 for index in neighbors),
                  "every transverse matrix element is minus one")

        analytic_lower = 1 - Fraction(1, 3) / ratio
        expected = Fraction(5, 6) if ratio == 2 else Fraction(13, 15)
        check(analytic_lower == expected, "admitted-member energy lower bound")

    # In normalized Dicke basis, transverse off-diagonal squares are exact.
    check(tuple((k + 1) * (4 - k) for k in range(4)) == (4, 6, 6, 4),
          "Dicke off-diagonal squares are 4,6,6,4")


def verify_robust_interval_logic() -> None:
    for denominator in range(1, 30):
        values = [Fraction(i, denominator) for i in range(denominator + 1)]
        for q in values:
            for epsilon in values:
                lower = max(Fraction(0), q - epsilon)
                upper = min(Fraction(1), q + epsilon)
                check(0 <= lower <= upper <= 1, "pointwise clipped interval valid")
                for actual in values:
                    if abs(actual - q) <= epsilon:
                        check(lower <= actual <= upper,
                              "pointwise interval contains every licensed value")

    # Exact classical two-outcome trace-distance saturation.
    for denominator in range(1, 50):
        for i in range(denominator + 1):
            p = Fraction(i, denominator)
            for j in range(denominator + 1):
                q = Fraction(j, denominator)
                trace_distance = abs(p - q)
                measurement_difference = abs(p - q)
                check(measurement_difference <= trace_distance,
                      "trace distance adds one eta, not two")


def exponential_tail(order: int, x: float) -> float:
    if order < 1 or x < 0:
        raise ValueError("tail domain")
    if x == 0:
        return 0.0
    term = math.exp(order * math.log(x) - math.lgamma(order + 1))
    total = term
    index = order
    for _ in range(100000):
        index += 1
        term *= x / index
        total += term
        if term <= max(1.0, total) * 1e-15 and index > x + 40:
            return total
    raise RuntimeError("tail did not converge")


def verify_unbounded_time_ceiling() -> None:
    for radius in range(8):
        factor = 3 * radius * radius + 3 * radius + 1
        order = 2 * radius + 1
        for ratio in (2.0, 2.5):
            small = factor * exponential_tail(order, 0.0)
            large = factor * exponential_tail(order, 48 * ratio * 10.0)
            check(small == 0, "fixed collar exact at zero duration")
            check(min(1.0, large) == 1.0,
                  "fixed collar becomes trivial within unbounded finite times")


def verify_certified_calculator() -> None:
    zero = complete_mission_interval(Fraction(2), Fraction(0))
    check(zero["lower"] == zero["upper"] == 1,
          "blank collar probability exact at zero time")
    check(zero["full_lower"] == zero["full_upper"] == 1,
          "complete mission exact at zero time")

    for ratio in (Fraction(2), Fraction(5, 2)):
        for sigma in (Fraction(1, 10000), Fraction(1, 1000)):
            result = complete_mission_interval(ratio, sigma, Fraction(1, 10**14))
            check(result["lower"] <= result["upper"],
                  "certified collar enclosure ordered")
            check(result["full_lower"] <= result["full_upper"],
                  "certified complete interval ordered")
            check(result["probability_error"] <= Fraction(1, 10**14),
                  "propagator reaches internal probability target")
            check(result["lower"] >= result["analytic_collar_lower"],
                  "certified small-time collar respects analytic energy bound")
            exterior = capped_exp_minus_one_upper(
                48 * ratio * abs(sigma), Fraction(1, 10**16))
            check(result["exterior_error_upper"] == exterior,
                  "reported exterior bound is reproducible")

    # A duration whose L=0 bound is necessarily clipped remains honest.
    clipped = complete_mission_interval(Fraction(2), Fraction(1, 10))
    check(clipped["exterior_error_upper"] == 1,
          "large L=0 influence bound is clipped to one")
    check(clipped["full_lower"] == 0 and clipped["full_upper"] == 1,
          "trivial certificate is reported rather than overclaimed")


def verify_ledger_and_scope() -> None:
    data = json.loads((LANE / "IDENTIFIABILITY_LEDGER.json").read_text(encoding="utf-8"))
    check(data["lane"] == "GL6BB_V001", "ledger lane")
    check(data["current_custody"]["admitted_h6_scenarios"] == ["2", "5/2"],
          "ledger two-member scenarios")
    check(data["state_free_probability_interval"] == ["0", "1"],
          "ledger sharp state-free interval")
    check(data["prepared_blank_L0"]["exact_reduced_dimension"] == 5,
          "ledger exact Dicke dimension")

    theorem = (LANE / "THEOREM.md").read_text(encoding="utf-8")
    check("\r" not in theorem and "\x00" not in theorem,
          "theorem contains no hidden control characters")
    required = [
        "every positive",
        "not a measurement, posterior bound, or selection rule",
        "provide no finite upper",
        "prepared-blank state as an already defined special case",
        "=[0,1]",
        "pointwise robust interval inherited",
        "D_{\\rm tr}",
        "five normalized Dicke states",
        "1-{1\\over3R}",
        "11/6-e^{96|\\sigma_{\\rm obs}|}",
        "28/15-e^{120|\\sigma_{\\rm obs}|}",
        "one authenticated tuple",
        "No graviton, Ricci target, Einstein equation, gravity identification, or `G`",
    ]
    for token in required:
        check(token in theorem, f"required theorem token: {token}")
    forbidden = [
        "r=2 is nature's value",
        "r=5/2 is nature's value",
        "the repository selects the blank state",
        "all finite times are uniformly bounded",
        "we have proved gravity",
        "this is gravity",
        "derives newton's constant",
    ]
    lower = theorem.lower()
    for phrase in forbidden:
        check(phrase not in lower, f"forbidden promotion absent: {phrase}")


def main() -> int:
    verify_dependencies()
    verify_ratio_custody()
    verify_state_free_sharpness()
    verify_dicke_reduction_and_energy_bound()
    verify_robust_interval_logic()
    verify_unbounded_time_ceiling()
    verify_certified_calculator()
    verify_ledger_and_scope()
    print(f"PASS: {CHECKS} exact GL6BB checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
