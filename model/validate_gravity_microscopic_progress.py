#!/usr/bin/env python3
"""Fail-closed zero-input validator for the GL6T--GL6AY URM checkpoint."""

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
    check(isinstance(certificate["controlled_evidence"], MappingProxyType))
    check(isinstance(certificate["custody"]["packets"], tuple))
    check(refused(lambda: gmp._root_path("../outside")))
    check(refused(lambda: gmp._verify_seal(gmp._PACKETS[-1].author_dir, "0" * 64, "0" * 64)))
    check(type_error(lambda: mutate(certificate, "schema", "changed")))
    check(type_error(lambda: mutate(certificate["exact_results"], "gravity", True)))

    check(certificate["custody"]["packet_count"] == 21)
    check(tuple(row["gate"] for row in certificate["custody"]["packets"]) == (
        "GL6T", "GL6U", "GL6AA", "GL6AF", "GL6AG", "GL6AH", "GL6AI",
        "GL6AK", "GL6AM", "GL6AN", "GL6AO", "GL6AP", "GL6AQ", "GL6AR",
        "GL6AS", "GL6AT", "GL6AU", "GL6AV", "GL6AW", "GL6AX", "GL6AY",
    ))
    check(certificate["custody"]["declared_hash_rows_checked"] == 1012)
    check(all(row["audit_disposition"].startswith("PASS") for row in certificate["custody"]["packets"]))
    check(all(len(row["author_claim_sha256"]) == 64 for row in certificate["custody"]["packets"]))
    check(all(len(row["audit_sha256"]) == 64 for row in certificate["custody"]["packets"]))
    check(certificate["custody"]["packets"][7]["audit_claim_file"] == "POSTFREEZE_AUDIT.md")
    check(certificate["custody"]["packets"][15]["author_claim_file"] == "RESULT.md")
    check(certificate["custody"]["packets"][15]["theorem_sha256"] is None)
    check(certificate["custody"]["historical_audit_count"] == 1)
    check(len(certificate["custody"]["historical_audits"]) == 1)
    check(certificate["custody"]["historical_audits"][0]["gate"] == "GL6AY_PRE_REPAIR_FAIL")
    check(certificate["custody"]["historical_audits"][0]["audit_disposition"] == (
        "FAIL__REPAIR_REQUIRED__GLOBAL_LOCK_PROJECTOR_AND_DRESSED_SUBSPACE_SCOPE"
    ))
    check(certificate["custody"]["historical_audits"][0]["status"] == (
        "SUPERSEDED_BY_REPAIRED_AUTHOR_AND_POST_REPAIR_PASS"
    ))

    exact = certificate["exact_results"]
    expected_legacy_exact = {
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
    check(all(exact[key] == value for key, value in expected_legacy_exact.items()))
    check(tuple(exact) == (
        "record_gated_local_pair_response",
        "interaction_owned_nonfactorization",
        "authenticated_finite_atlas",
        "formation_pattern_source_threshold",
        "matched_neighbor_response",
        "homogeneous_direct_edge_propagation",
        "uniform_quasilocal_influence_envelope",
        "authenticated_a3_bulk_dynamics",
        "authenticated_finite_window_bulk_response",
        "native_degree_lock_sector",
        "complete_order6_locked_hamiltonian",
        "locked_ir_representation_and_response_boundary",
        "authenticated_E_loop_selection_boundary",
        "locked_hexagon_thermodynamic_sector",
        "native_hexagon_collective_response",
        "order6_quantum_ice_crosswalk",
        "vg0_first_character_static_closure",
        "record_conditioned_collective_clock_and_typed_atlas",
        "anisotropic_folner_twist_closure",
        "all_fixed_order_port_and_twist_stability",
        "finite_coupling_prethermal_locked_bridge",
    ))
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

    check(exact["authenticated_a3_bulk_dynamics"] == {
        "premise": "SELECTED_HOMOGENEOUS_ALL_FORMED_F3_MEMBER_ON_A3_X_FOUR_PORTS",
        "site_set": "L=A3_x_{1,2,3,4}",
        "pair_interaction_degree": 6,
        "lambda_F3": "48*abs(U_d)/hbar",
        "boundary_comparison": (
            "norm(tau_t^S(A)-tau_t^R(A))<=3*norm(A)*abs(X)*"
            "sum_{r=R}^infinity(2r+1)^3*T_{r-r_X+1}(lambda_F3*abs(t))"
        ),
        "bulk_limit": "BOUNDARY_INDEPENDENT_STRONGLY_CONTINUOUS_QUASILOCAL_DYNAMICS",
        "stationary_state": "AT_LEAST_ONE_JOINT_TIME_TRANSLATION_S4_INVARIANT_STATE_EXISTS_NOT_SELECTED",
        "spectral_measure": "mu_AB(B)=mu_A1(B)*P_A1+mu_E(B)*P_E+mu_T2(B)*P_T2",
    })
    check(exact["authenticated_finite_window_bulk_response"]["retarded_kernel"] == (
        "G^R_{beta,alpha}(t)=i*E_star^2/(2*hbar)*Theta(t)*"
        "omega([tau_t(M_beta),M_alpha])"
    ))
    check(exact["authenticated_finite_window_bulk_response"]["commutator_measure"] == (
        "NOT_POSITIVE_WITHOUT_PASSIVITY_OR_KMS"
    ))
    check(exact["native_degree_lock_sector"]["hamiltonian"] == (
        "H=-h*sum_e X_e+U_d*sum_v(k_v-2)^2+C"
    ))
    check(exact["native_degree_lock_sector"]["linear_Ward_no_go"] == (
        "[H,k_v-2]=-i*h*sum_{e_incident_v}Y_e"
    ))
    check("893M/1080" in exact["complete_order6_locked_hamiltonian"]["formula"])
    check(exact["complete_order6_locked_hamiltonian"]["hexagon_amplitude"] == "-63/8")
    check(exact["locked_ir_representation_and_response_boundary"]["mismatch"] == (
        "Hom_S4(T2,E)=0"
    ))
    check(exact["authenticated_E_loop_selection_boundary"]["uniform_locked_variance"] == (
        "(8/3)*norm(c)^2"
    ))
    check("product_{e_in_C}kappa_e" in (
        exact["authenticated_E_loop_selection_boundary"]["retained_support_gate"]
    ))
    check(exact["locked_hexagon_thermodynamic_sector"]["finite_component"] == (
        "H_C=-t*A_C;Delta_C=t*(rho_C-lambda_2(C))"
    ))
    check(exact["locked_hexagon_thermodynamic_sector"]["variance_bound"] == (
        "Delta_L<=18*t*norm(w)_infinity^2*L/Var(F_L)"
    ))
    check(exact["native_hexagon_collective_response"]["conserved_density"] == (
        "CENTERED_PORT_T2"
    ))
    check(exact["native_hexagon_collective_response"]["single_mode_bound"] == (
        "Delta_T2(chi;u)<=f_u(chi)/S_u_plus(chi)"
    ))
    check(exact["order6_quantum_ice_crosswalk"]["exact_parameter"] == "v/g=0")
    check(exact["order6_quantum_ice_crosswalk"]["distinct_RK_point"] == "v/g=1")
    check(exact["vg0_first_character_static_closure"]["oscillator_bound"] == (
        "f_u(q_L)<=6*J*sin(pi/L)^2"
    ))
    check(exact["vg0_first_character_static_closure"]["component_gap_bound"] == (
        "Delta_C(L)<=6*J*sin(pi/L)^2/S_u,L(q_L)"
    ))
    check("alpha<2" in (
        exact["vg0_first_character_static_closure"]["static_exponent_premise"]
    ))
    check("O(J/L)" in (
        exact["vg0_first_character_static_closure"]["alpha_one_consequence"]
    ))
    check("product_{e_in_c}kappa_e" in (
        exact["record_conditioned_collective_clock_and_typed_atlas"]["record_conditioned_loop"]
    ))
    check("q^6*H_hex(1)" in (
        exact["record_conditioned_collective_clock_and_typed_atlas"]["homogeneous_formal_family"]
    ))
    check(exact["record_conditioned_collective_clock_and_typed_atlas"]["tetrahedral_evaluation"] == (
        "rank(E)=4;im(E)=A1_plus_T2;ker(E)=E"
    ))
    check("NOT_ONE_PHYSICAL_METRIC_TANGENT" in (
        exact["record_conditioned_collective_clock_and_typed_atlas"]["typed_atlas"]
    ))
    check(exact["anisotropic_folner_twist_closure"]["translation_character"] == (
        "Y*U0*Y^-1=-U0"
    ))
    check("4*pi^2*J*L0*L2/L1" in (
        exact["anisotropic_folner_twist_closure"]["gap_bound"]
    ))
    check(exact["anisotropic_folner_twist_closure"]["folner_sequence"] == (
        "(L0,L1,L2)=(m,2*m^3,m);m>=5_ODD"
    ))
    check("SELECTED_GNS" in (
        exact["anisotropic_folner_twist_closure"]["closure_scope"]
    ))
    check(exact["all_fixed_order_port_and_twist_stability"]["contractible_port_law"] == (
        "Delta_N_a=0_FOR_a=0,1,2,3"
    ))
    check(exact["all_fixed_order_port_and_twist_stability"]["minimum_winding_hamming_distance"] == (
        "2*L_min"
    ))
    check("r<2*L_min" in (
        exact["all_fixed_order_port_and_twist_stability"]["fixed_order_scope"]
    ))
    check("T_L" in (
        exact["all_fixed_order_port_and_twist_stability"]["quasilocal_gap_dichotomy"]
    ))
    check("NO_UNIFORM_FINITE_COUPLING" in (
        exact["all_fixed_order_port_and_twist_stability"]["exact_ceiling"]
    ))
    check(exact["finite_coupling_prethermal_locked_bridge"]["parent"] == (
        "H=U_d*N_def-h*sum_e(X_e);N_def=sum_v(k_v-2)^2"
    ))
    check("REMAINDER_RETAINED" in (
        exact["finite_coupling_prethermal_locked_bridge"]["exact_normal_form"]
    ))
    check("P_S^0=chi(N_S=0)" in (
        exact["finite_coupling_prethermal_locked_bridge"]["local_collar"]
    ))
    check("D_2(L)<=" in (
        exact["finite_coupling_prethermal_locked_bridge"]["finite_second_twist_moment"]
    ))
    check("0<r_1<ln(3/2)/4" in (
        exact["finite_coupling_prethermal_locked_bridge"]["local_observable_horizon"]
    ))
    check("P_L_TO_Q_L_LEAKAGE_NEEDS_NO_WINDING" in (
        exact["finite_coupling_prethermal_locked_bridge"]["topology_boundary"]
    ))
    check("NOT_EXACT_ALL_TIME_LOCKED_PHASE" in (
        exact["finite_coupling_prethermal_locked_bridge"]["exact_ceiling"]
    ))

    evidence = certificate["controlled_evidence"]
    check(tuple(evidence) == ("quantum_ice_v_over_g_zero",))
    check(evidence["quantum_ice_v_over_g_zero"]["claim_class"] == (
        "NUMERICAL_AND_EFFECTIVE_EVIDENCE_NOT_A_PHASE_GAP_OR_POLE_THEOREM"
    ))
    check("0.6_plus_or_minus_0.1" in evidence["quantum_ice_v_over_g_zero"]["Benton"])

    check(len(certificate["open_gates"]) == 10)
    check(certificate["open_gates"][0].startswith("ISOTROPIC_V_OVER_G_ZERO_PHASE_CONTROL"))
    check(certificate["open_gates"][1].startswith("EXACT_ALL_TIME_FINITE_H_OVER_U_D"))
    check(certificate["open_gates"][2].startswith("PHYSICAL_U_D_OVER_H_AND_CLOCK_CALIBRATION"))
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

    expected = 101
    check(passed + 1 == expected)
    print(f"GRAVITY MICROSCOPIC PROGRESS: {passed}/{expected} PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
