#!/usr/bin/env python3
"""Certified GL6BB prepared-blank L=0 collar calculation.

The physical dimensionless observation time has no default: callers must
supply it.  Decimal inputs are parsed as exact rational numbers.  The collar
propagator is summed on the equivalent sixteen-state basis with exact rational
real/imaginary parts; an induced-norm exponential-tail bound encloses the
omitted series.  The sealed GL6BA L=0 exterior bound is then added.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, localcontext
from fractions import Fraction
from typing import Iterable


Vector = list[Fraction]


def parse_fraction(text: str) -> Fraction:
    """Parse an integer, fraction, decimal, or scientific-notation string."""
    try:
        value = Fraction(text)
    except (ValueError, ZeroDivisionError) as exc:
        raise argparse.ArgumentTypeError(f"not an exact number: {text}") from exc
    return value


def fraction_decimal(value: Fraction, digits: int = 18) -> str:
    with localcontext() as context:
        context.prec = digits + 8
        result = Decimal(value.numerator) / Decimal(value.denominator)
        return format(result, f".{digits}g")


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def exponential_remainder_upper(x: Fraction, order: int) -> Fraction:
    """Upper-bound sum_{n=order+1}^infinity x^n/n! for x>=0.

    Once order+2>x, every ratio after the first omitted term is bounded by
    x/(order+2), so a geometric majorant is exact rational arithmetic.
    """
    if x < 0 or order < 0:
        raise ValueError("tail requires x>=0 and order>=0")
    if x == 0:
        return Fraction(0)
    ratio = x / (order + 2)
    if ratio >= 1:
        raise ValueError("order too small for geometric tail majorant")
    first = x ** (order + 1)
    for integer in range(2, order + 2):
        first /= integer
    return first / (1 - ratio)


def choose_propagator_order(x: Fraction, tolerance: Fraction) -> tuple[int, Fraction]:
    if x < 0 or tolerance <= 0:
        raise ValueError("x must be nonnegative and tolerance positive")
    if x == 0:
        return 0, Fraction(0)
    start = max(0, ceil_fraction(x) - 1)
    for order in range(start, 5000):
        if x / (order + 2) >= 1:
            continue
        remainder = exponential_remainder_upper(x, order)
        probability_error = remainder * (2 + remainder)
        if probability_error <= tolerance:
            return order, remainder
    raise RuntimeError("failed to certify propagator within 5000 terms")


def diagonal_energy(state: int, ratio: Fraction) -> Fraction:
    weight = bin(state).count("1")
    return ratio * weight * (weight - 7)


def hamiltonian_action(vector: Vector, ratio: Fraction) -> Vector:
    if len(vector) != 16:
        raise ValueError("L=0 active vector must have length 16")
    result = [Fraction(0) for _ in range(16)]
    for state in range(16):
        result[state] += diagonal_energy(state, ratio) * vector[state]
        for bit in range(4):
            result[state] -= vector[state ^ (1 << bit)]
    return result


def blank_collar_probability_interval(
    ratio: Fraction,
    sigma: Fraction,
    tolerance: Fraction = Fraction(1, 10**12),
) -> dict[str, Fraction | int]:
    """Certified exact-series enclosure for q_0^blank(R,sigma)."""
    if ratio <= 0:
        raise ValueError("R must be positive")
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")

    absolute_sigma = abs(sigma)
    # Maximum absolute row sum of the exact sixteen-state Hamiltonian.
    norm_bound = 12 * ratio + 4
    series_argument = norm_bound * absolute_sigma
    order, unitary_remainder = choose_propagator_order(series_argument, tolerance)

    real = [Fraction(0) for _ in range(16)]
    imaginary = [Fraction(0) for _ in range(16)]
    power = [Fraction(0) for _ in range(16)]
    power[0] = Fraction(1)
    scalar = Fraction(1)

    for degree in range(order + 1):
        phase = degree % 4
        if phase == 0:
            for index, value in enumerate(power):
                real[index] += scalar * value
        elif phase == 1:  # (-i)^1
            for index, value in enumerate(power):
                imaginary[index] -= scalar * value
        elif phase == 2:
            for index, value in enumerate(power):
                real[index] -= scalar * value
        else:  # (-i)^3 = +i
            for index, value in enumerate(power):
                imaginary[index] += scalar * value

        if degree != order:
            power = hamiltonian_action(power, ratio)
            scalar *= sigma / (degree + 1)

    approximate = Fraction(0)
    for state in range(16):
        bit_zero = (state >> 0) & 1
        bit_one = (state >> 1) & 1
        if bit_zero == bit_one:
            approximate += real[state] ** 2 + imaginary[state] ** 2

    probability_error = unitary_remainder * (2 + unitary_remainder)
    lower = max(Fraction(0), approximate - probability_error)
    upper = min(Fraction(1), approximate + probability_error)
    return {
        "order": order,
        "norm_bound": norm_bound,
        "unitary_remainder": unitary_remainder,
        "probability_error": probability_error,
        "approximate": approximate,
        "lower": lower,
        "upper": upper,
    }


def capped_exp_minus_one_upper(
    x: Fraction,
    tolerance: Fraction = Fraction(1, 10**15),
) -> Fraction:
    """Rational upper bound on min(1, exp(x)-1), x>=0."""
    if x < 0 or tolerance <= 0:
        raise ValueError("x must be nonnegative and tolerance positive")
    if x == 0:
        return Fraction(0)

    partial = Fraction(1)
    term = Fraction(1)
    for degree in range(1, 10000):
        term *= x / degree
        partial += term
        if partial >= 2:
            return Fraction(1)
        ratio = x / (degree + 2)
        if ratio >= 1:
            continue
        first_omitted = term * x / (degree + 1)
        tail = first_omitted / (1 - ratio)
        upper = partial + tail
        if upper < 2 and tail <= tolerance:
            return upper - 1
    raise RuntimeError("failed to enclose capped exponential within 10000 terms")


def complete_mission_interval(
    ratio: Fraction,
    sigma: Fraction,
    tolerance: Fraction = Fraction(1, 10**12),
) -> dict[str, Fraction | int]:
    collar = blank_collar_probability_interval(ratio, sigma, tolerance)
    exterior = capped_exp_minus_one_upper(48 * ratio * abs(sigma), tolerance / 100)
    full_lower = max(Fraction(0), collar["lower"] - exterior)  # type: ignore[operator]
    full_upper = min(Fraction(1), collar["upper"] + exterior)  # type: ignore[operator]
    return {
        **collar,
        "exterior_error_upper": exterior,
        "full_lower": full_lower,
        "full_upper": full_upper,
        "analytic_collar_lower": max(Fraction(0), 1 - Fraction(1, 3) / ratio),
    }


def result_as_json(
    ratio: Fraction,
    sigma: Fraction,
    tolerance: Fraction,
    result: dict[str, Fraction | int],
) -> str:
    encoded: dict[str, object] = {
        "scope": "conditional prepared-blank L=0 collar; sigma supplied by caller",
        "R": str(ratio),
        "sigma_obs": str(sigma),
        "internal_probability_tolerance": str(tolerance),
        "admitted_h6_member": ratio in (Fraction(2), Fraction(5, 2)),
    }
    for key, value in result.items():
        if isinstance(value, Fraction):
            encoded[key] = {
                "exact": str(value),
                "decimal": fraction_decimal(value),
            }
        else:
            encoded[key] = value
    return json.dumps(encoded, indent=2, sort_keys=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ratio", required=True, type=parse_fraction,
                        help="selected or scenario value R=U_d/h, e.g. 2 or 5/2")
    parser.add_argument("--sigma", required=True, type=parse_fraction,
                        help="same-clock dimensionless observation time")
    parser.add_argument("--tolerance", type=parse_fraction,
                        default=Fraction(1, 10**12),
                        help="internal rational probability enclosure target")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.ratio <= 0:
        raise SystemExit("--ratio must be positive")
    if args.tolerance <= 0:
        raise SystemExit("--tolerance must be positive")
    result = complete_mission_interval(args.ratio, args.sigma, args.tolerance)
    print(result_as_json(args.ratio, args.sigma, args.tolerance, result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
