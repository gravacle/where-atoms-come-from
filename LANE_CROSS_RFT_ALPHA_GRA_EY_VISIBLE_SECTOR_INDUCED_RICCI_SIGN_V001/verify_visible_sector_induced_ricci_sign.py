#!/usr/bin/env python3
"""Exact arithmetic checks for VSIRS (standard library only)."""

from fractions import Fraction as Q
from pathlib import Path


ROOT = Path(__file__).resolve().parent
THEOREM = (ROOT / "THEOREM.md").read_text(encoding="utf-8")
THEOREM_FLAT = " ".join(THEOREM.split())


def require(condition, label):
    if not condition:
        raise AssertionError(label)
    print(f"PASS {label}")


# Visser-convention low-spin coefficients used by the theorem.
k_scalar_minimal = Q(1, 6)
k_weyl = Q(-1, 6)
k_vector = Q(-2, 3)

# Minimal visible gauge-basis Standard Model: 4 real Higgs components,
# 45 two-component Weyl fields, and 12 gauge vectors.  The explicit minus
# before the Weyl term is the single (-1)^F supertrace weight.
visible_minimal = 4 * k_scalar_minimal - 45 * k_weyl + 12 * k_vector
require(visible_minimal == Q(1, 6), "minimal-Higgs visible supertrace is 1/6")

# A common Higgs nonminimal coupling contributes -4 xi_H.
# At the sign boundary xi_H=1/24 the supertrace vanishes.
xi_boundary = Q(1, 24)
visible_at_boundary = visible_minimal - 4 * xi_boundary
require(visible_at_boundary == 0, "Higgs sign boundary is 1/24")
require(visible_minimal - 4 * Q(0) > 0, "minimal Higgs gives negative Ricci coefficient")
visible_conformal = visible_minimal - 4 * Q(1, 6)
require(visible_conformal == Q(-1, 2), "conformal Higgs supertrace is -1/2")

# Six additional minimal real scalar modes worsen the minimal-Higgs sign.
visible_plus_six = visible_minimal + 6 * k_scalar_minimal
require(visible_plus_six == Q(7, 6), "six scalar pair modes give 7/6 at minimal Higgs")

# From C_R=-str(k1) Delta(kappa^2)/(32 pi^2), store only the exact
# rational numerator multiplying Delta(kappa^2)/pi^2.
c_visible = -visible_minimal / 32
c_visible_plus_six = -visible_plus_six / 32
require(c_visible == Q(-1, 192), "minimal visible C_R factor is -1/192")
require(c_visible_plus_six == Q(-7, 192), "visible-plus-six C_R factor is -7/192")
require(-visible_conformal / 32 == Q(1, 64), "conformal visible C_R factor is 1/64")

# G=1/(16 pi C_R) gives the theorem's exact 12 pi numerator.
# Algebraically: 1 / (16 pi * [(24xi-1-Np)/(192 pi^2)] Delta)
# = 12 pi / [(24xi-1-Np) Delta], in the convention where
# EY02 gives the same result after the exact supertrace census.
# Retain the explicit pi comparison as a textual normalization guard.
require("{12\\pi\\over" in THEOREM, "conditional G_eff numerator recorded as 12 pi")
require("-{\\operatorname{str}k_1\\over2\\pi}" in THEOREM, "Visser 1/G normalization recorded")
require("-{\\operatorname{str}k_1\\over32\\pi^2}" in THEOREM, "Ricci coefficient normalization recorded")

# Claim-scope guards: a positive visible contribution is not the total without
# spectrum and matching custody, and Lambda is not derived.
for phrase, label in [
    ("complete ultraviolet spectrum", "incomplete-spectrum ceiling"),
    ("zero-bare-term matching", "matching ceiling"),
    ("cosmological-constant problem", "Lambda ceiling"),
    ("not promoted for the actual world", "actual-world ceiling"),
    ("R_R=-R_V", "curvature-convention dictionary"),
    ("not an exact momentum projector", "proper-time shell typing"),
    ("not an assumption of a thermal cosmological", "nonthermal ultraviolet typing"),
    ("two gauge ancestors", "joint B-W3 photon ancestry"),
    ("complete anomaly-free representation set", "Weyl anomaly premise"),
]:
    require(phrase in THEOREM_FLAT, label)

print("PASS VSIRS exact checks")
