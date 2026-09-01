# GL6AY primary-source replay

The audit checked the exact source versions named in the frozen author
packet.

## ADHH closed-system normal form

Source: Abanin, De Roeck, Ho, and Huveneers,
arXiv:`1509.05386v3`, 9 July 2017,
<https://arxiv.org/pdf/1509.05386v3>.

The versioned PDF gives in Section 3:

- integer-spectrum local number operators and `N=sum_xN_x`;
- `G=nu N+H`, `D=<H>`, `V=H-D`;
- equation (3.1), including `nu>=9pi||V||_kappa0/kappa0`;
- `nu_0=(54pi/kappa_0^2)(||D||+2||V||)` and the stated `n_*`;
- Theorem 3.1, exact `YGY^*=nu N+D_hat+V_hat`,
  `[D_hat,N]=0`, remainder `(2/3)^n_*||V||`, and quasi-local
  conjugation control;
- Theorem 3.3, local-observable error `K_3(O)/nu` through
  `exp(r_1n_*)` for `r_1<ln(3/2)/(d+1)`.

The source explicitly explains that the iteration stops at finite optimal
order; it does not supply a convergent all-orders block diagonalization.

## Strong-support extension

Source: Else, Fendley, Kemp, and Nayak,
arXiv:`1704.08703v2`, 26 September 2017, Appendix A,
<https://arxiv.org/pdf/1704.08703v2>.

Appendix A replaces ordinary support with strong support for commuting
finite-radius summands of `N`.  A strongly supported term commutes with every
constraint summand whose full support is not contained in its strong support.
Evolution by `N` preserves that support, and commutators take unions, so the
ADHH proof carries over.

## Whole-band boundary

Source: Bravyi, DiVincenzo, and Loss,
arXiv:`1105.0675v1`, 3 May 2011, Section 4,
<https://arxiv.org/pdf/1105.0675v1>.

Section 4 states that the global perturbation norm in a many-body lattice is
typically extensive while the unperturbed gap is not.  The standard global
SW transformation may then be undefined and its Taylor series may diverge;
fixed-order truncations and linked-cluster structure remain controlled.  It
does not prove a universal impossibility of exact block diagonalization for
special models.
