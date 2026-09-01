#!/usr/bin/env python3
"""Independent-style constructive replay for GL6AZ V001.

Standard library only.  Every check is explicit so ``python -O`` has the same
semantics and count as an ordinary run.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import random
from fractions import Fraction
from pathlib import Path


LANE = Path(__file__).resolve().parent
ROOT = LANE.parent
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"FAIL [{CHECKS}] {label}")


def close(a: float, b: float, label: str, tol: float = 1e-12) -> None:
    check(abs(a - b) <= tol * max(1.0, abs(a), abs(b)), label)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


EXPECTED_DEPENDENCIES = {
    "LANE_CROSS_RFT_GRA_GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE_V001/THEOREM.md": "5b86ab5eb2998eb719dffd09e05add131863fd2a3290d87fb749dc8aebc1891c",
    "LANE_CROSS_RFT_GRA_GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE_V001/MANIFEST.sha256": "e81ec1cfd4bdcdc43b4709b8f90f9eceac3dfba82be80701dd4a2a7e08de089b",
    "LANE_CROSS_RFT_GRA_GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE_V001/SEAL.sha256": "740f051b3347d7387e481a9991f536bd61a7e47ad51d80b680475dff394e5cbb",
    "LANE_CROSS_RFT_GRA_GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE_V001/PRIMARY_SOURCES.md": "4d45c19b4ae742e1d27b397e03cf1679491a335bb3b736a5b4c439198842b325",
    "AUDIT_G_GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE_POST_REPAIR_V001/AUDIT.md": "9828df0f2a376ed1615a617772b87a667a6e31ed77cc2237b72cbdc196f2556e",
    "AUDIT_G_GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE_POST_REPAIR_V001/MANIFEST.sha256": "77af0925b66f5c854285f25d574d4a995ba173337c77e634d8da94f2a86c2975",
    "AUDIT_G_GL6AY_FINITE_COUPLING_PRETHERMAL_LOCKED_BRIDGE_POST_REPAIR_V001/SEAL.sha256": "3f787b28c44bb896aaff2e98c84691a7fa86bcddb3e5b2b4b6934c8a6b398bf1",
    "LANE_CROSS_RFT_GRA_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001/THEOREM.md": "8407cee5196bfa4240f02159a5f59f941903dcf7a10e2baa18cf52a01ac8f743",
    "LANE_CROSS_RFT_GRA_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001/MANIFEST.sha256": "8b4ac6f6ceda2acb117480201ee96ce22be97fe0a99d8c097d8267100efa8c44",
    "AUDIT_G_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001/AUDIT.md": "2659392e6f2fe3c0062068426faf5d516cceaf3a106017742b8ea88c21517b00",
    "AUDIT_G_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001/MANIFEST.sha256": "d35a2a9b2e581db963ca0513a26ffbcbbbbae28efd7c3a019dfa4a9f0db50301",
    "AUDIT_G_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001/SEAL.sha256": "d1d70e294e4be73f8efdc26301a37bf11808f9f248b433b3ec5aaf7f87e324c7",
    "LANE_CROSS_RFT_GRA_GL6V_Q4_PAIR_PHASE_SOURCE_READ_DILATION_V001/THEOREM.md": "3f5b52aa066d4d6f56f75a06a1f6623d49d988c531d0b9ef82e590dd92aec51d",
    "LANE_CROSS_RFT_GRA_GL6V_Q4_PAIR_PHASE_SOURCE_READ_DILATION_V001/MANIFEST.sha256": "b23723bc85ec2083399fa981a0b0100306afe4470f30196706dfa7beb2e377d9",
    "AUDIT_G_GL6V_Q4_PAIR_PHASE_SOURCE_READ_DILATION_V001/README.md": "de7aa00b7ee99b4919086e5e0313cef587dcd24b236faca3c1d121fd234164d2",
    "AUDIT_G_GL6V_Q4_PAIR_PHASE_SOURCE_READ_DILATION_V001/MANIFEST.sha256": "e602349fd2d584a06fb5cefafb7cf29e431814faee9ae22cf303cf7dfe10db38",
    "AUDIT_G_GL6V_Q4_PAIR_PHASE_SOURCE_READ_DILATION_V001/SEAL.sha256": "b4308fb9f7af68f7738ab555ee9fe38b371ee1cda588a21042f4e28ef4a95c66",
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
    "LANE_RFT_UDCL_CONDITIONAL_UNIVERSAL_COVERAGE_V001/THEOREM.md": "30cdd93998b556306e72ee7cee9ad434e84d2e7e8ddc08a820c01394895ebc32",
    "LANE_RFT_UDCL_CONDITIONAL_UNIVERSAL_COVERAGE_V001/TYPED_CLARIFICATION_V002.md": "3570f3a2f5851e9620935b235439b67bc3a4c45420d2f2149c8eeeae4f90b34c",
    "LANE_RFT_UDCL_CONDITIONAL_UNIVERSAL_COVERAGE_V001/MANIFEST.sha256": "3bfe5cf2614ae0df3775e36084afe4870197c425c01f9de0ada83f34a04de1d3",
    "LANE_RFT_UDCL_CONDITIONAL_UNIVERSAL_COVERAGE_V001/AUDIT.md": "bd81c3270fcd8e29ee6ec230becebdfbc4a9366e943f72180b819b8e908f608c",
    "LANE_RFT_UDCL_CONDITIONAL_UNIVERSAL_COVERAGE_V001/INDEPENDENT_AUDIT.md": "bbc4b2fffd9a76cb2026ec6988f60d8102cbad5dc6f69cc64bd0924cf1b9c0da",
}


def verify_dependencies() -> None:
    for rel, expected in sorted(EXPECTED_DEPENDENCIES.items()):
        path = ROOT / rel
        check(path.is_file(), f"dependency exists: {rel}")
        check(sha256(path) == expected, f"dependency hash: {rel}")


def verify_one_flip_native_identity() -> None:
    # Target link plus the three other links at each degree-four endpoint.
    # Enumerate every locked local endpoint pattern compatible with the target
    # being absent or occupied.
    sample_scales = [
        (Fraction(2, 3), Fraction(5, 7)),
        (Fraction(7, 5), Fraction(11, 13)),
        (Fraction(19, 4), Fraction(23, 9)),
    ]
    for target in (0, 1):
        required_other = 2 - target
        triples = [bits for bits in itertools.product((0, 1), repeat=3)
                   if sum(bits) == required_other]
        check(len(triples) == 3, "three locked endpoint completions")
        for left in triples:
            for right in triples:
                k_left = target + sum(left)
                k_right = target + sum(right)
                check(k_left == 2 and k_right == 2, "initial endpoints locked")
                flipped = 1 - target
                k_left_f = flipped + sum(left)
                k_right_f = flipped + sum(right)
                ndef_initial = (k_left - 2) ** 2 + (k_right - 2) ** 2
                ndef_final = (k_left_f - 2) ** 2 + (k_right_f - 2) ** 2
                check(ndef_initial == 0, "initial local defect number zero")
                check(ndef_final == 2, "one flip creates two unit defects")
                for ud, h in sample_scales:
                    diagonal_gap = ud * (ndef_final - ndef_initial)
                    offdiag = abs(-h)
                    check(diagonal_gap == 2 * ud,
                          "native diagonal defect separation is 2 U_d")
                    check(offdiag == h, "native flip amplitude is h")
                    check(diagonal_gap / (2 * offdiag) == ud / h,
                          "native ratio identity")


def verify_complete_pair_coarsening() -> None:
    outcomes = list(itertools.product((-1, 1), repeat=4))
    check(len(outcomes) == 16, "complete four-link read has sixteen outcomes")
    for a in range(4):
        for b in range(a + 1, 4):
            plus = [q for q in outcomes if q[a] * q[b] == 1]
            minus = [q for q in outcomes if q[a] * q[b] == -1]
            check(len(plus) == 8 and len(minus) == 8,
                  "pair coarsening partitions complete read 8+8")
            for q in outcomes:
                m = q[a] * q[b]
                p_plus = Fraction(1 + m, 2)
                p_minus = Fraction(1 - m, 2)
                check(p_plus in (0, 1) and p_minus in (0, 1),
                      "pair spectral projectors are Boolean on read basis")
                check(p_plus + p_minus == 1, "pair PVM resolves identity")


def verify_binary_total_variation() -> None:
    # Dense exact rational battery.  A binary M expectation is 2 p_plus - 1.
    for den in range(1, 81):
        for i in range(den + 1):
            p = Fraction(i, den)
            m_p = 2 * p - 1
            for j in range(den + 1):
                q = Fraction(j, den)
                m_q = 2 * q - 1
                tv = (abs(p - q) + abs((1 - p) - (1 - q))) / 2
                check(tv == abs(p - q), "binary TV equals plus-probability contrast")
                check(tv == abs(m_p - m_q) / 2,
                      "binary TV equals half M-expectation contrast")


def a_of_x(x: Fraction) -> Fraction:
    return Fraction(8, 63) / (x ** 6) * (
        1 - x ** 2 - Fraction(37, 12) * x ** 4
        - Fraction(16247, 900) * x ** 6
    )


def verify_two_member_witness() -> None:
    x1 = Fraction(2, 5)
    x2 = Fraction(1, 2)
    check(1 / x1 == Fraction(5, 2), "x=2/5 gives R=5/2")
    check(1 / x2 == 2, "x=1/2 gives R=2")
    a1 = a_of_x(x1)
    a2 = a_of_x(x2)
    check(a1 == Fraction(2415673, 113400), "first exact GL5ZZF coefficient")
    check(a2 == Fraction(31706, 14175), "second exact GL5ZZF coefficient")
    check(a1 - a2 == Fraction(3203, 168), "exact source-owner difference")
    for j6 in (Fraction(1, 7), Fraction(3, 2), Fraction(29, 11)):
        descendants = []
        for x in (x1, x2):
            ud = Fraction(8, 63) * j6 / (x ** 6)
            h = x * ud
            reconstructed = Fraction(63, 8) * ud * (h / ud) ** 6
            check(reconstructed == j6, "same source-free H6 coefficient")
            descendants.append((ud, h))
        check(descendants[0] != descendants[1], "microscopic members remain distinct")


def verify_all_positive_ratios_admitted() -> None:
    # GL6AN inherited ray: epsilon_*=-6 U_d iff
    # Delta=4 U_d(d_*-2), with h,Delta,U_d>0 and d_*>2.
    ratios = [Fraction(n, d) for d in range(1, 41) for n in range(1, 41)]
    h_values = [Fraction(1, 7), Fraction(5, 3), Fraction(31, 11)]
    d_stars = [Fraction(3, 1)]
    for ratio in ratios:
        for h in h_values:
            ud = ratio * h
            for d_star in d_stars:
                delta = 4 * ud * (d_star - 2)
                epsilon = delta + 2 * ud * (1 - 2 * d_star)
                check(h > 0 and ud > 0 and delta > 0,
                      "constructed strict inherited member is positive")
                check(epsilon == -6 * ud, "constructed member lies on lock ray")
                check(ud / h == ratio, "constructed member realizes requested R")


def verify_sufficient_domain_boundary() -> None:
    threshold_small = 36.0 * math.pi * math.e
    threshold_sep = 432.0 * math.pi * math.e ** 2
    close(threshold_small, 307.4304320162484, "36 pi e numerical value")
    close(threshold_sep, 10028.190682380982, "432 pi e^2 numerical value")
    check(2.0 < threshold_small and 2.5 < threshold_small,
          "both admitted R below first-smallness floor")
    close(2.0 / threshold_small, 0.006505536836035462,
          "R=2 first-smallness fraction")
    close(2.5 / threshold_small, 0.008131921045044328,
          "R=5/2 first-smallness fraction")

    # f'(kappa) has sign(4*kappa-1); check the exact sign identity densely and
    # the analytic minimizer numerically over many decades.
    kappas = [10 ** (-4 + 6 * i / 20000) for i in range(20001)]
    minimum = float("inf")
    argmin = None
    for kappa in kappas:
        f = 9.0 * math.pi * math.exp(4.0 * kappa) / kappa
        derivative = 9.0 * math.pi * math.exp(4.0 * kappa) * (4.0 * kappa - 1.0) / (kappa ** 2)
        check((derivative < 0) == (kappa < 0.25) or abs(kappa - 0.25) < 1e-14,
              "derivative sign below minimizer")
        check((derivative > 0) == (kappa > 0.25) or abs(kappa - 0.25) < 1e-14,
              "derivative sign above minimizer")
        g = 108.0 * math.pi * math.exp(4.0 * kappa) / (kappa ** 2)
        g_derivative_sign = 4.0 - 2.0 / kappa
        check(f + 1e-10 >= threshold_small,
              "sampled first-smallness threshold above analytic floor")
        check(g + 1e-8 >= threshold_sep,
              "sampled scale-separation threshold above analytic floor")
        check((g_derivative_sign < 0) == (kappa < 0.5)
              or abs(kappa - 0.5) < 1e-14,
              "scale-separation derivative sign below minimizer")
        check((g_derivative_sign > 0) == (kappa > 0.5)
              or abs(kappa - 0.5) < 1e-14,
              "scale-separation derivative sign above minimizer")
        if f < minimum:
            minimum, argmin = f, kappa
    check(argmin is not None, "threshold grid has minimizer")
    check(abs(math.log(argmin / 0.25)) < 5e-4, "grid minimizer near kappa=1/4")

    for power in range(-120, 121):
        kappa = 0.25 * math.exp(power / 20)
        f = 9.0 * math.pi * math.exp(min(4.0 * kappa, 700.0)) / kappa
        if 4.0 * kappa < 700.0:
            check(f >= threshold_small * (1 - 1e-13),
                  "log-grid first-smallness floor")

    # The compact theorem statement's n_*>=1 does not select the high branch.
    # This exact low-branch counterexample motivates restoring nu>=nu_0.
    low_x = math.exp(-0.9)
    low_value = low_x / (1.0 + math.log(low_x)) ** 3
    close(low_value, 406.5696597405994, "low-branch n-star quotient")
    check(math.floor(low_value) - 2 == 404, "low branch has n_star>=1")
    check(low_x < 1.0, "low branch violates scale separation")

    # On x>=1, n_*>=1 iff floor(f(x))>=3, hence f(x)>=3.  The
    # unique high root lies above e^2, where f is increasing.
    def high_branch_f(x: float) -> float:
        return x / (1.0 + math.log(x)) ** 3

    check(high_branch_f(1.0) < 3.0, "high branch starts below n-star threshold")
    check(high_branch_f(math.e ** 2) < 3.0, "high-branch minimum below threshold")
    lo, hi = math.e ** 2, 1_000_000.0
    check(high_branch_f(hi) > 3.0, "high-root bracket closes")
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if high_branch_f(mid) < 3.0:
            lo = mid
        else:
            hi = mid
    x_star = (lo + hi) / 2.0
    close(x_star, 1861.32559690908, "unique high-branch x_star", tol=1e-13)
    close(high_branch_f(x_star), 3.0, "x_star solves n-star boundary")
    full_floor = x_star * threshold_sep
    close(full_floor, 18665728.00780086, "complete universal proof-domain floor",
          tol=1e-12)
    check(2.0 < full_floor and 2.5 < full_floor,
          "both admitted members below complete floor")


def verify_grouped_strong_support_witness() -> None:
    displacements = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
    for a, da in enumerate(displacements):
        support = {
            tuple(da[i] - db[i] for i in range(3))
            for db in displacements
        }
        check(len(support) == 4, f"native strong support has four cells, port {a}")

    # Even if several edge potentials are grouped under one identical support,
    # only the target flip connects n to n^target.
    for width in range(2, 9):
        for n in itertools.product((0, 1), repeat=width):
            for target in range(width):
                target_state = list(n)
                target_state[target] ^= 1
                target_state = tuple(target_state)
                contributing = []
                for edge in range(width):
                    out = list(n)
                    out[edge] ^= 1
                    if tuple(out) == target_state:
                        contributing.append(edge)
                check(contributing == [target], "grouped edge sum isolates target flip")


def verify_scale_covariance() -> None:
    rng = random.Random(602061)
    for _ in range(10000):
        ud = Fraction(rng.randint(1, 1000), rng.randint(1, 1000))
        h = Fraction(rng.randint(1, 1000), rng.randint(1, 1000))
        t = Fraction(rng.randint(1, 1000), rng.randint(1, 1000))
        c = Fraction(rng.randint(1, 1000), rng.randint(1, 1000))
        r_before = ud / h
        s_before = h * t
        r_after = (c * ud) / (c * h)
        s_after = (c * h) * (t / c)
        check(r_before == r_after, "common energy scaling preserves R")
        check(s_before == s_after, "inverse time scaling preserves dimensionless clock")


def verify_theorem_scope() -> None:
    theorem = (LANE / "THEOREM.md").read_text(encoding="utf-8")
    required = [
        "R={U_d\\over h}",
        "\\sigma_{\\rm obs}:={h(t_Q-t_F)\\over\\hbar}",
        "D_{\\rm TV}(p^H,p^{eff})",
        "{\\Delta_{\\rm def}\\over2A_X}",
        "36\\pi e\\approx307.430432",
        "432\\pi e^2\\approx10028.190682",
        "18665728.0078",
        "R\\ge\\bar\\nu_0",
        "selected ideal ready/`MATCH` active factor",
        "sharper/different finite-mission local-observable theorem",
        "not physical phase failure",
        "No exact infinite-volume locked projector",
        "No graviton, Ricci target, Einstein equation, gravity identification, or",
    ]
    for token in required:
        check(token in theorem, f"required theorem scope token: {token}")
    forbidden = [
        "we have proved gravity",
        "this is gravity",
        "derives Newton's constant",
        "the admitted members are unstable",
        "necessary and sufficient prethermal",
        "exact all-time finite-coupling locked phase is proved",
    ]
    lower = theorem.lower()
    for phrase in forbidden:
        check(phrase not in lower, f"forbidden promotion absent: {phrase}")


def verify_ledger() -> None:
    data = json.loads((LANE / "IDENTIFIABILITY_LEDGER.json").read_text(encoding="utf-8"))
    check(data["lane"] == "GL6AZ_V001", "ledger lane")
    check(data["dimensionless_inputs"]["R"] == "U_d/h", "ledger R")
    check(len(data["admitted_members"]) == 2, "ledger two admitted members")
    check(data["application_gate"][1] == "R >= nubar_0", "ledger restored scale separation")
    check("0 < r_1 < ln(3/2)/4" in data["application_gate"],
          "ledger local-dynamics exponent range")
    check(data["proof_domain_floors"]["complete_high_branch_universal"]
          == "R >= 18665728.0078", "ledger complete floor")
    check("outside_domain" in data["decision_fork"], "ledger explicit fallback branch")
    check(data["interpretation"] == "theorem-domain failure, not physical-phase failure",
          "ledger ceiling wording")
    for item in ("no graviton", "no Ricci", "no gravity", "no G"):
        check(item in data["ceilings"], f"ledger ceiling {item}")


def main() -> None:
    verify_dependencies()
    verify_one_flip_native_identity()
    verify_complete_pair_coarsening()
    verify_binary_total_variation()
    verify_grouped_strong_support_witness()
    verify_all_positive_ratios_admitted()
    verify_two_member_witness()
    verify_sufficient_domain_boundary()
    verify_scale_covariance()
    verify_theorem_scope()
    verify_ledger()
    print(f"PASS: {CHECKS}/{CHECKS} GL6AZ constructive checks")


if __name__ == "__main__":
    main()
