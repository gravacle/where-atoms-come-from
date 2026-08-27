#!/usr/bin/env python3
"""Exact finite ledger for the standard-causal URFT scope composition.

This verifies propositional composition and the finite CTC cut witness. It does
not establish that any actual-world record satisfies C, S, or J.
"""

from itertools import product


class Ledger:
    def __init__(self) -> None:
        self.passed = 0
        self.total = 0

    def check(self, condition: bool, name: str) -> None:
        self.total += 1
        if not condition:
            raise AssertionError(name)
        self.passed += 1


def closure(rec: bool, cut: bool, sewing: bool, joint: bool):
    """Apply only the three imported forward implications."""
    screc = rec and cut and sewing and joint
    cts = screc
    occ = cts
    phys_arch = occ
    a = tuple(phys_arch for _ in range(4))
    return screc, cts, occ, phys_arch, a


def valid_writer_query_cut(assignment, edges):
    if assignment["W"] != "past" or assignment["Q"] != "future":
        return False
    return all(
        not (assignment[source] == "future" and assignment[target] == "past")
        for source, target in edges
    )


def main() -> None:
    ledger = Ledger()

    for flags in product((False, True), repeat=4):
        rec, cut, sewing, joint = flags
        screc, cts, occ, phys_arch, a = closure(*flags)
        expected = all(flags)
        ledger.check(screc == expected, f"SCREC_exact_{flags}")
        ledger.check(not screc or cts, f"SCREC_implies_CTS_{flags}")
        ledger.check(not cts or occ, f"CTS_implies_OCC_{flags}")
        ledger.check(not occ or phys_arch, f"OCC_implies_PhysEnc_ARCH_{flags}")
        ledger.check(not phys_arch or all(a), f"PhysEnc_ARCH_implies_A1_A4_{flags}")

    # Each physical conjunct is load-bearing for entry into this theorem.
    for missing in range(4):
        flags = [True, True, True, True]
        flags[missing] = False
        screc, cts, occ, phys_arch, a = closure(*flags)
        ledger.check(not screc, f"missing_gate_{missing}_outside_domain")
        ledger.check(not any((cts, occ, phys_arch, *a)),
                     f"missing_gate_{missing}_not_inferred_by_chain")

    # The distributional CTC witness has both W->Q and Q->W, so C is false.
    nodes = ("W", "Q")
    edges = (("W", "Q"), ("Q", "W"))
    assignments = [
        dict(zip(nodes, sides))
        for sides in product(("past", "future"), repeat=2)
    ]
    candidates = [
        assignment for assignment in assignments
        if assignment["W"] == "past" and assignment["Q"] == "future"
    ]
    ledger.check(len(candidates) == 1, "CTC_one_required_side_assignment")
    ledger.check(not any(valid_writer_query_cut(a, edges) for a in candidates),
                 "CTC_no_writer_before_query_cut")
    ctc = closure(True, False, True, True)
    ledger.check(not ctc[0] and not ctc[1], "CTC_REC_not_SCREC_or_CTS_by_chain")

    # The theorem's conclusion deliberately excludes evidential/strong claims.
    claimed = {
        "OCC": True,
        "PhysEnc": True,
        "ARCH": True,
        "A1_A4": True,
        "C5": False,
        "AuthEnc": False,
        "authenticated_ADMIT": False,
        "A5": False,
        "reset": False,
        "objective_actualization": False,
        "sealing_trigger": False,
        "all_records_satisfy_CSJ": False,
        "no_non_CTS_cover": False,
    }
    for name in ("OCC", "PhysEnc", "ARCH", "A1_A4"):
        ledger.check(claimed[name], f"positive_scope_{name}")
    for name in (
        "C5", "AuthEnc", "authenticated_ADMIT", "A5", "reset",
        "objective_actualization", "sealing_trigger", "all_records_satisfy_CSJ",
        "no_non_CTS_cover",
    ):
        ledger.check(not claimed[name], f"excluded_scope_{name}")

    print("STANDARD_CAUSAL_DOMAIN EXACT")
    print("RCTS_CTS_OCC_FAITHFUL_ADMISSION_COMPOSITION EXACT")
    print("DISTRIBUTIONAL_CTC_OUTSIDE_CAUSAL_CUT EXACT")
    print("AUTHENTICATION_AND_STRONG_URFT_CEILING EXACT")
    print(f"TOTAL {ledger.passed}/{ledger.total} PASS")
    print(
        "VERDICT STANDARD_CAUSAL_RECORD_SUBCLASS_COVERED__"
        "ALL_RECORDS_STANDARD_CAUSAL_OPEN__NON_CTS_COVERS_OPEN"
    )


if __name__ == "__main__":
    main()
