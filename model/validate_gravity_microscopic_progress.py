#!/usr/bin/env python3
"""Fail-closed zero-input validator for the GL6T--GL6AI URM checkpoint."""

from __future__ import annotations

import inspect
from types import MappingProxyType

import gravity_formation_theory as gft
import gravity_microscopic_progress as gmp
from project_model import URM


passed = 0


def check(condition: bool) -> None:
    global passed
    if not condition:
        raise AssertionError("gravity microscopic progress validation failed")
    passed += 1


def refused(callable_object) -> bool:
    try:
        callable_object()
    except gmp.GravityMicroscopicProgressRefusal:
        return True
    return False


def type_error(callable_object) -> bool:
    try:
        callable_object()
    except TypeError:
        return True
    return False


def mutate(mapping, key, value) -> None:
    mapping[key] = value


def main() -> int:
    check(tuple(inspect.signature(gmp.gravity_microscopic_progress).parameters) == ())
    check(tuple(inspect.signature(gmp.gravity_microscopic_progress_certificate).parameters) == ())
    check(type_error(lambda: gmp.gravity_microscopic_progress({})))
    check(type_error(lambda: gmp.gravity_microscopic_progress_certificate(data={})))

    progress = gmp.gravity_microscopic_progress()
    certificate = progress.certificate()
    check(progress.claim_class == gmp.CLAIM_CLASS)
    check(certificate["schema"] == gmp.SCHEMA)
    check(certificate["claim_class"] == gmp.CLAIM_CLASS)
    check(certificate["relationship_to_V014"] == "ADDITIVE_PROGRESS_SURFACE__V014_MEANING_UNCHANGED")
    check(isinstance(certificate, MappingProxyType))
    check(isinstance(certificate["exact_results"], MappingProxyType))
    check(isinstance(certificate["custody"]["packets"], tuple))
    check(refused(lambda: gmp._root_path("../outside")))
    check(refused(lambda: gmp._verify_seal(gmp._PACKETS[-1].author_dir, "0" * 64, "0" * 64)))
    check(type_error(lambda: mutate(certificate, "schema", "changed")))
    check(type_error(lambda: mutate(certificate["exact_results"], "gravity", True)))

    check(certificate["custody"]["packet_count"] == 7)
    check(tuple(row["gate"] for row in certificate["custody"]["packets"]) == (
        "GL6T", "GL6U", "GL6AA", "GL6AF", "GL6AG", "GL6AH", "GL6AI"
    ))
    check(certificate["custody"]["declared_hash_rows_checked"] == 404)
    check(all(row["audit_disposition"].startswith("PASS") for row in certificate["custody"]["packets"]))
    check(all(len(row["theorem_sha256"]) == 64 for row in certificate["custody"]["packets"]))
    check(all(len(row["audit_sha256"]) == 64 for row in certificate["custody"]["packets"]))

    exact = certificate["exact_results"]
    expected_exact = {
        "record_gated_local_pair_response": {
            "premise": "N=0;U_d=0;h>0;Delta>0;0<vartheta<2*pi",
            "KEEP": "D=-8*h*x*I6-4*h*x*z^2*A_L",
            "sectors": (
                "D_A1=-8*h*x*(1+2*z^2)",
                "D_E2=-8*h*x*(1-z^2)",
                "D_T2=-8*h*x",
            ),
            "open_domain_rank": 6,
            "matched_BREAK": "D_BREAK=0",
        },
        "interaction_owned_nonfactorization": (
            "y-x*z^2=-(16/3)*h^3*U_d*(tau/hbar)^4+O((tau/hbar)^6)"
        ),
        "authenticated_finite_atlas": {
            "shared_child": "J_e^1=J_f^1 iff c=c_prime",
            "edge_distance": "d_G(m,n)=one_half*norm_1(n-m)=1",
            "translation_cocycle": "sum_closed_path(delta)=0",
            "frame_cocycle": "T_ln*T_nm=T_lm",
            "premise": "INDEPENDENTLY_FORMED_NODE_ID_AND_PORT_LABEL_RECORDS",
        },
        "formation_pattern_source_threshold": {
            "premise": "N=0;U_d=0;h>0;Delta>0;0<vartheta<2*pi",
            "fixed_E_restriction": (
                "ZERO_OR_ONE_FORMED_LINK_GIVES_ZERO__EXACTLY_TWO_GIVES_RANK_ONE"
            ),
            "two_formed": "E^T*D^(ab)*E=-16*h*x*(1-z)*w_ab^T*w_ab",
            "scope": "FIXED_E_X_E_RESTRICTION_NOT_BROKEN_S4_INVARIANT_BLOCK",
        },
        "matched_neighbor_response": {
            "premise": (
                "N=1;h=U_d=E0>0;Delta=6*E0;d_star=2;delta=0;"
                "receiver_c_in_{1,2,3}"
            ),
            "reference": "MATCHED_TO_SOURCE_PATTERN_0000_WITH_RECEIVERS_FIXED_FORMED_KEEP",
            "expansion_variable": "u=E0*tau/hbar",
            "expansion": (
                "Delta_mu_c^E=-(5626/42525)*u^12*kappa_c*w_0c+O(u^14)"
            ),
        },
        "homogeneous_direct_edge_propagation": {
            "premise": "GENERAL_N_DIRECT_EDGE;delta=0",
            "raw": "Delta_Q_n6=-128*h^4*U_d^2*u_b_in",
            "normalized_witness": "Delta_mu_n=(8/45)*u^6*u_b_in+O(u^8)",
            "irrep": "u_b_in=one_half*one+one_half*s_b in A1_plus_T2",
            "direct_E": "E_TRANSPOSE*x_q=0_AT_EVERY_ORDER",
            "displayed_helper_E": (
                "NONZERO_ON_INCOMPLETE_ENDPOINT_INCIDENCE__CANCELS_ON_COMPLETE_"
                "HOMOGENEOUS_INTERIOR_ENDPOINT_INCIDENCE"
            ),
        },
        "uniform_quasilocal_influence_envelope": {
            "active_link_degree": "degree_L(e)=q(child)+2<=6",
            "lambda_F3": "48*abs(U_d)/hbar",
            "tail": (
                "norm_commutator<=2*norm(A)*norm(B)*T_distance(lambda_F3*abs(t))"
            ),
            "exponential": (
                "norm_commutator<=2*norm(A)*norm(B)*exp(-mu*(d_cell-v_mu*abs(t)))"
            ),
            "v_mu": "lambda_F3*exp(mu)/mu",
            "v_1": (
                "48*e*abs(U_d)/hbar_CELL_STEPS_PER_PARENT_TIME_UPPER_ENVELOPE_ONLY"
            ),
        },
    }
    check(exact == expected_exact)
    check(exact["record_gated_local_pair_response"]["open_domain_rank"] == 6)
    check(exact["record_gated_local_pair_response"]["matched_BREAK"] == "D_BREAK=0")
    check("U_d" in exact["interaction_owned_nonfactorization"])
    check(exact["authenticated_finite_atlas"]["premise"] == "INDEPENDENTLY_FORMED_NODE_ID_AND_PORT_LABEL_RECORDS")
    check(
        "EXACTLY_TWO_GIVES_RANK_ONE"
        in exact["formation_pattern_source_threshold"]["fixed_E_restriction"]
    )
    check(
        "5626/42525" in exact["matched_neighbor_response"]["expansion"]
        and exact["matched_neighbor_response"]["expansion_variable"]
        == "u=E0*tau/hbar"
        and "O(u^14)" in exact["matched_neighbor_response"]["expansion"]
        and "O(tau^14)" not in exact["matched_neighbor_response"]["expansion"]
    )
    check(
        exact["homogeneous_direct_edge_propagation"]["premise"]
        == "GENERAL_N_DIRECT_EDGE;delta=0"
        and exact["homogeneous_direct_edge_propagation"]["direct_E"]
        == "E_TRANSPOSE*x_q=0_AT_EVERY_ORDER"
    )
    check("CANCELS" in exact["homogeneous_direct_edge_propagation"]["displayed_helper_E"])
    check(exact["uniform_quasilocal_influence_envelope"]["lambda_F3"] == "48*abs(U_d)/hbar")
    check("UPPER_ENVELOPE_ONLY" in exact["uniform_quasilocal_influence_envelope"]["v_1"])

    check(len(certificate["open_gates"]) == 6)
    check(certificate["open_gates"][0].startswith("STATIONARY_BULK_IR_LAW"))
    check(certificate["open_gates"][-1] == "MICROSCOPIC_RICCI_COEFFICIENT_AND_G_MODEL")
    check(all(value is False for value in certificate["ceilings"].values()))
    check(certificate["executable_scope"]["caller_arguments"] == 0)
    check(certificate["executable_scope"]["physics_recalculated"] is False)
    check(certificate["executable_scope"]["gravity_solver"] is False)

    delegated = URM.gravity_microscopic_progress_certificate()
    check(delegated["schema"] == gmp.SCHEMA)
    check(URM.gravity_microscopic_progress().claim_class == gmp.CLAIM_CLASS)
    check(URM.gravity_formation_theory_certificate()["schema"] == gft.SCHEMA)
    check(gft.SCHEMA == "WAC_GRAVITY_FORMATION_THEORY_CERTIFICATE_V014")
    check(delegated["ceilings"]["gravity_derived_here"] is False)
    check(delegated["ceilings"]["G_calculated_here"] is False)

    expected = 46
    check(passed + 1 == expected)
    print(f"GRAVITY MICROSCOPIC PROGRESS: {passed}/{expected} PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
