# GL6AY finite-coupling prethermal locked bridge

**Status:** author frozen and sealed; a distinct independent hostile audit is
required before promotion.

This packet gives the shortest rigorous bridge from the sealed `GL6AX`
all-fixed-order theorem to a genuinely nonzero finite ratio `h/U_d`.  The
native F3 Hamiltonian

```text
H=U_d N_def-h sum_e X_e,
N_def=sum_v(k_v-2)^2
```

meets the primary-source hypotheses of the closed-system prethermal normal
form after the finite-range, commuting constraint terms are handled with
strong support.  At sufficiently small but finite `h/U_d`, uniformly in
volume, an exact quasi-local change of frame gives

```text
Y H Y^*=U_d N_def+D_hat+V_hat,
[D_hat,N_def]=0,
||V_hat||_(kappa_*) <= (2/3)^n_* ||V_0||_(kappa_0).
```

The identity is exact; replacing the full Hamiltonian by
`U_d N_def+D_hat` is prethermal.  For every strongly supported term define
the finite collar projector onto zero defect for all constraint stars fully
contained in its support.  The resulting local interaction

```text
Phi_S=P_S^0 D_hat(S) P_S^0
```

agrees with `D_hat(S)` on globally locked configurations, conserves all four
native port totals there term by term, and has a finite local second twist
moment.  It therefore falls exactly within the controlled quasi-local form
of `GL6AX` without introducing a nonexistent infinite-volume global lock
projector.

The remainder is not promoted away.  Local charged excursions can leave the
locked space without winding.  Winding is required only for a process that
begins and ends locked while changing a port total.  Local/potential dressing
is controlled, but no volume-uniform global dressed-subspace norm is claimed.
Raw
whole-band Feshbach/direct-rotation bounds lose their convergence radius with
volume, and an explicit finite-torus winding coefficient is nonzero at its
first allowed order.  Thus the packet proves a finite-coupling **prethermal**
locked bridge, not an exact microscopic finite-coupling phase, gravity, or
`G`.
