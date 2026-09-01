# GL6AT primary-literature quantum-ice crosswalk

**Status:** primary-source screen frozen on 2026-08-31; exact effective-model
crosswalk established; no new microscopic or infrared theorem; independent
audit required before promotion.

## Decisive result

After removal of its common scalar, the sealed GL6AO order-six locked
Hamiltonian is exactly the zero-potential fully-packed-loop quantum-ice
Hamiltonian on the diamond graph:

\[
 H_{\rm AO}^{(6)}-E_{\rm scalar}P
 =-g\sum_cT_c,
 \qquad g={63\over8}{h^6\over U_d^5}>0,
 \qquad {v\over g}=0.
\]

Here every vertex has exactly two occupied links, and `T_c` toggles an
alternating elementary hexagon and annihilates a nonalternating one.  The
comparison potential `v` counts flippable hexagons; GL6AO has no such term.
The exactly soluble Rokhsar--Kivelson point is `v/g=1`, not the GL6AO point.

The mapping is exact for the displayed order-six linked interaction.  It is
not an exact all-orders identification of the finite-`h/U_d` F3 Hamiltonian:
higher-order diagonal and longer-loop terms have not been excluded.

## What the primary papers support

- Hermele--Fisher--Balents give the exact two-dimers-per-diamond-site / ice
  constraint and alternating-hexagon operator map, including the sign-changing
  unitary and the distinct RK point.
- Shannon--Sikora--Pollmann--Penc--Fulde study this same fully-packed-loop
  Hamiltonian by zero-temperature GFMC and finite exact diagonalization.  Their
  evidence places `v/g=0` inside a `U(1)` liquid region.  This is numerical
  phase evidence, not a rigorous spectral theorem.
- Benton--Sikora--Shannon fit a lattice Gaussian field theory to equal-time
  QMC and infer a linear low-character photon branch and a velocity at the
  zero-potential point.  Their real-frequency pole formula is field-theoretic,
  not direct real-time QMC.
- Fu--Rau--Gingras--Perkins show, within an approximate XXZ/slave-particle
  Raman treatment, that a bilinear electric-field tensor can create two
  photons.  Their active cubic `T2g` operator is not the sealed local pair-`E`
  irrep, so it does not establish pair-`E` spectral overlap.

No primary result located in this screen proves the `v/g=0` phase, its gap, or
its pole structure mathematically.  No located primary result computes the
exact positive-frequency spectral measure of the sealed local pair-`E`
operator in the pure projected Hamiltonian.

## Packet map

- `RESULT.md` gives the exact graph/operator/parameter crosswalk and all
  mapping conditions.
- `PRIMARY_SOURCES.md` records only primary research papers, direct DOI and
  arXiv links, equation anchors, and use ceilings.
- `EVIDENCE_LADDER.md` separates exact algebra, exact results at the wrong
  parameter, numerical phase inference, effective dispersion, and approximate
  composite response.
- `DEPENDENCIES.md` and `DEPENDENCIES.sha256` pin the sealed GL6AO and GL6AP
  author/audit bytes used here.
- `SELF_AUDIT.md` records exclusions and promotion attacks.
- `verify_packet.py`, `VERIFICATION.txt`, `MANIFEST.sha256`, and `SEAL.sha256`
  provide fail-closed local custody.  The verifier checks the recorded source
  metadata and claim structure; it does not pretend to rederive the papers.

Nothing in this packet identifies a translation character with calibrated
physical momentum or derives a physical length, speed, cone, Ricci tensor,
gravity, or Newton's constant `G`.
