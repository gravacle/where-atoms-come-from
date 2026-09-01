# GL6AY primary-source pinning

The external theorem inputs are pinned to exact arXiv versions.  GL6AY uses
their stated hypotheses and bounds, not informal recollections or later
review-language.

## 1. Closed-system prethermal normal form

- D. Abanin, W. De Roeck, W. W. Ho, and F. Huveneers,
  *A rigorous theory of many-body prethermalization for periodically driven
  and closed quantum systems*, arXiv:`1509.05386v3` (9 July 2017), especially
  Sections 2.2 and 3.1--3.3:
  <https://arxiv.org/html/1509.05386v3>.
- Exact theorem objects used: the exponential potential norm in equation
  (2.2 of the paper's section numbering display), the closed-system
  hypotheses (3.1), Theorems 3.1--3.3, and the volume-independent constants.
- Scope used: a finite spin lattice with finite local Hilbert space,
  integer-spectrum number operator, bounded local potential, exact
  `N`-commuting normal form plus an exponentially small potential-norm
  remainder, and quasi-exponential local-dynamics control.
- Exact ceiling retained: Section 1.2.3 explicitly explains why the
  elimination is stopped at an optimal finite order rather than continued
  to an assumed convergent all-orders reduction.

## 2. Finite-range overlapping constraint extension

- D. V. Else, P. Fendley, J. Kemp, and C. Nayak,
  *Prethermal Strong Zero Modes and Topological Qubits*,
  arXiv:`1704.08703v2` (26 September 2017), Appendix A:
  <https://arxiv.org/html/1704.08703v2>.
- Exact theorem object used: the extension from single-site summands of `N`
  to commuting finite-radius summands through **strong support**.  A term
  strongly supported on `S` commutes with every constraint summand whose
  support is not contained in `S`; the exponential potential norm is then
  evaluated on those strong supports.  Evolution by `N` preserves strong
  support, and commutators take unions, so the ADHH proof carries over.

## 3. Whole-band Schrieffer--Wolff boundary

- S. Bravyi, D. P. DiVincenzo, and D. Loss,
  *Schrieffer-Wolff transformation for quantum many-body systems*,
  arXiv:`1105.0675v1` (3 May 2011), especially Section 4:
  <https://arxiv.org/html/1105.0675v1>.
- Exact scope used: in a many-body lattice the global perturbation norm is
  extensive while the unperturbed gap is not.  The global direct-rotation
  criterion therefore need not be volume-uniform, the Taylor series may
  diverge, and their controlled many-body result is a fixed-order linked
  truncation rather than an exact all-orders block diagonalization.

No phase result, particle ontology, gauge photon, graviton, Ricci ansatz,
Einstein equation, or value of `G` is imported from these papers.

