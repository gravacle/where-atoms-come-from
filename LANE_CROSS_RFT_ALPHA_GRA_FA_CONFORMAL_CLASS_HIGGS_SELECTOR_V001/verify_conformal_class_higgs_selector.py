#!/usr/bin/env python3
"""Deterministic exact checks for CCHS (standard library only)."""

from fractions import Fraction as Q
from pathlib import Path


ROOT = Path(__file__).resolve().parent
THEOREM = (ROOT / "THEOREM.md").read_text(encoding="utf-8")
FLAT = " ".join(THEOREM.split())


def require(condition, label):
    if not condition:
        raise AssertionError(label)
    print(f"PASS {label}")


def conformal_data(D):
    """Exact CCHS coefficients for an integer dimension D>2."""
    D = Q(D)
    w = -(D - 2) / 2
    xi = (D - 2) / (4 * (D - 1))
    wp = w - 2
    cross = 2 * w + D - 2
    box_sigma = -w - 2 * xi * (D - 1)
    grad_sigma_sq = -w * (w + D - 2) - xi * (D - 1) * (D - 2)
    return w, xi, wp, cross, box_sigma, grad_sigma_sq


# Check the symbolic identities by their exact rational specialization in a
# broad set of dimensions.  THEOREM.md contains the algebraic derivation.
for D in range(3, 33):
    w, xi, wp, cross, box_sigma, grad_sigma_sq = conformal_data(D)
    require(w == Q(-(D - 2), 2), f"D={D} scalar weight")
    require(xi == Q(D - 2, 4 * (D - 1)), f"D={D} conformal xi")
    require(wp == Q(-(D + 2), 2), f"D={D} output weight")
    require(cross == 0, f"D={D} derivative cross term vanishes")
    require(box_sigma == 0, f"D={D} laplacian-sigma term vanishes")
    require(grad_sigma_sq == 0, f"D={D} gradient-sigma-square term vanishes")

# The actual q4/D4 composition.
w4, xi4, wp4, *_ = conformal_data(4)
require(w4 == -1, "D4 Higgs weight is -1")
require(wp4 == -3, "D4 equation output weight is -3")
require(xi4 == Q(1, 6), "D4 Visser-table conformal coupling is 1/6")

# Clean EY visible-only arithmetic with one statistics insertion.
str_visible = Q(1, 6) - 4 * xi4
c_visible = -str_visible / 32
require(str_visible == Q(-1, 2), "EY visible supertrace is -1/2")
require(c_visible == Q(1, 64), "EY visible Ricci factor is 1/64")

# Additional minimal scalar margin at xi=1/6.
for Np in range(0, 9):
    str_with_pairs = Q(1 + Np, 6) - 4 * xi4
    c_with_pairs = -str_with_pairs / 32
    require(c_with_pairs == Q(3 - Np, 192), f"Np={Np} scalar-margin coefficient")
    require((c_with_pairs > 0) == (Np < 3), f"Np={Np} positivity classification")

# Scope and no-go guards.
for phrase, label in [
    ("PREVOLUME-WEYL-HIGGS", "new physical premise named"),
    ("characteristic propagation alone", "characteristic-only no-go"),
    ("does not enter the principal symbol", "principal-symbol reason"),
    ("not supplied by the current", "adopted-premise ceiling"),
    ("trace/Weyl anomaly", "quantum anomaly ceiling"),
    ("No beta function is assumed here", "RG ceiling"),
    ("visible-only (`N_p=0`)", "visible-only EY typing"),
    ("no factor-of-two change", "Higgs real-component normalization"),
    ("positive only for integer `N_p=0,1,2`", "pair-scalar margin"),
    ("not an actual-world closure", "actual-world ceiling"),
    ("record-origin gravity theorem", "gravity ceiling"),
]:
    require(phrase in FLAT, label)

print("PASS CCHS exact checks")
