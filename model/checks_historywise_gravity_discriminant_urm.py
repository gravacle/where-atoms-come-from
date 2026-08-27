#!/usr/bin/env python3
"""Twenty-four URM delegation, isolation, ceiling, and conjunction checks."""

from __future__ import annotations

import ast
import inspect
from pathlib import Path
from unittest import mock

import historywise_gravity_discriminant as hgd
from project_model import URM


HERE = Path(__file__).resolve().parent


def refused(callable_object) -> bool:
    try:
        callable_object()
    except hgd.HistorywiseGravityDiscriminantRefusal:
        return True
    return False


def type_error(callable_object) -> bool:
    try:
        callable_object()
    except TypeError:
        return True
    return False


def main() -> int:
    checks = 0

    def check(condition: bool) -> None:
        nonlocal checks
        assert condition
        checks += 1

    public_names = {
        name for name in dir(URM) if name.startswith("historywise_gravity_")
    }
    check(public_names == {"historywise_gravity_discriminant", "historywise_gravity_discriminant_certificate"})  # 1
    check(tuple(inspect.signature(URM.historywise_gravity_discriminant).parameters) == ())  # 2
    check(tuple(inspect.signature(URM.historywise_gravity_discriminant_certificate).parameters) == ())  # 3

    delegated = URM.historywise_gravity_discriminant()
    check(isinstance(delegated, hgd.HistorywiseGravityDiscriminant))  # 4
    check(delegated.manifest_sha256 == hgd.MANIFEST_SHA256)  # 5
    delegated_certificate = URM.historywise_gravity_discriminant_certificate()
    direct_certificate = hgd.historywise_gravity_discriminant_certificate()
    check(delegated_certificate == direct_certificate)  # 6
    check(delegated_certificate is not direct_certificate)  # 7
    check(delegated_certificate["statuses"] is not direct_certificate["statuses"])  # 8

    check(type_error(lambda: URM.historywise_gravity_discriminant(True)))  # 9
    check(type_error(lambda: URM.historywise_gravity_discriminant(root="replacement")))  # 10
    check(type_error(lambda: URM.historywise_gravity_discriminant_certificate({})))  # 11
    check(type_error(lambda: URM.historywise_gravity_discriminant_certificate(packet={})))  # 12

    check(delegated_certificate["claim_class"] == "FORMAL_FINITE_GROUP_DISCRIMINANT_ONLY")  # 13
    check(not any(delegated_certificate["authorizations"].values()))  # 14
    check(delegated_certificate["statuses"]["caller_input_scientific_weight"] == "ZERO")  # 15
    check(delegated_certificate["statuses"]["physical_GARH_D"] == "NOT_ADMITTED_BY_THIS_DISCRIMINANT")  # 16
    check(delegated_certificate["statuses"]["GARH_Q"] == "NOT_DERIVED_NOT_FORCED_BY_THIS_DISCRIMINANT")  # 17
    check(delegated_certificate["nonpromotion"]["positive_boundary_promotes_physical_orientation"] is False)  # 18
    check(delegated_certificate["nonpromotion"]["failed_GARH_D_promotes_GARH_Q"] is False)  # 19
    check(delegated_certificate["executable_scope"]["general_finite_group_theorem_machine_proved"] is False)  # 20

    project_source = (HERE / "project_model.py").read_text(encoding="utf-8")
    check("return historywise_gravity_discriminant()" in project_source)  # 21
    check("return historywise_gravity_discriminant_certificate()" in project_source)  # 22

    with mock.patch.object(
        hgd,
        "historywise_gravity_discriminant",
        side_effect=hgd.HistorywiseGravityDiscriminantRefusal("injected custody failure"),
    ):
        check(refused(URM.historywise_gravity_discriminant))  # 23

    validate_source = (HERE / "validate_urm.py").read_text(encoding="utf-8")
    tree = ast.parse(validate_source)
    overall_assignments = [
        node for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "overall" for target in node.targets)
    ]
    overall_uses_gate = bool(overall_assignments) and any(
        isinstance(node, ast.Name) and node.id == "historywise_gravity_ok"
        for node in ast.walk(overall_assignments[-1].value)
    )
    chained = "validate_historywise_gravity_discriminant.py" in validate_source
    zero_weight = "zero physical or empirical proof weight" in validate_source
    check(overall_uses_gate and chained and zero_weight)  # 24

    assert checks == 24
    print("HISTORYWISE_GRAVITY_DISCRIMINANT_URM_CHECKS: 24/24 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

