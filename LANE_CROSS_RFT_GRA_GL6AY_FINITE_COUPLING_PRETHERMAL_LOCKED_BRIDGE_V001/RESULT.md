# GL6AY result

## Exact finite-coupling normal form

Let

```text
N_def=sum_v(k_v-2)^2,
H=U_d N_def-h sum_e X_e.
```

For any `kappa_0>0`, decompose the flip potential by the exact `N_def`
average as `W=D_0+V_0`.  If

```text
U_d >= 9 pi ||V_0||^str_(kappa_0)/kappa_0,
n_* = floor{(U_d/nu_0)/[1+ln(U_d/nu_0)]^3}-2 >=1,
nu_0=(54 pi/kappa_0^2)
     (||D_0||^str_(kappa_0)+2||V_0||^str_(kappa_0)),
```

then, uniformly in volume, there is an exact quasi-local change of frame

```text
Y H Y^*=U_d N_def+D_hat+V_hat,
[D_hat,N_def]=0,
||V_hat||^str_(kappa_*)
 <=(2/3)^n_* ||V_0||^str_(kappa_0),
kappa_*=kappa_0/[1+log(n_*+1)].
```

For every strongly supported term, set

```text
N_S=sum_(v:supp(q_v^2) subset S)q_v^2,
P_S^0=chi(N_S=0),
Phi_S=P_S^0 D_hat(S)P_S^0.
```

This is a finite-support interaction.  It agrees with `D_hat(S)` on globally
locked configurations and conserves the four native port totals between
such locked endpoints.  Its exponential strong-support norm gives a finite
local second twist moment, so the quasi-local `GL6AX` obstruction applies
exactly to this locally typed locked interaction.  In finite volume its
locked compression has first non-scalar coefficient

```text
-(63/8)(h^6/U_d^5) sum_c T_c.
```

## Operational finite-horizon corollary

For a local observable `O`, dimension `d=3`, and any
`0<r_1<ln(3/2)/4`, the primary theorem gives

```text
||tau_t^H(O)-tau_t^(U_d N_def+D_hat)(O)||
 <=K_3(O)/U_d,
t<=t_*:=exp(r_1 n_*),
```

in the theorem's clock units, with `K_3(O)` independent of volume and
`U_d`.  Therefore any separately declared finite observation or
record-retention horizon `t_obs<=t_*` is governed locally, to the displayed
error, by the exactly defect-conserving effective dynamics, whose locally
typed locked representation `{Phi_S}` conserves the four ports.  The F3
value of `U_d/h` and the conversion of `t_*` to a
physical clock remain open.  The inequality is not gravity.

## Exact ceiling

The full microscopic Hamiltonian is not proved to preserve the locked space.
In finite volume `P_L V_hat_L P_L=0`, but a local remainder term may make a
charged excursion `P_L -> Q_L` and change bare port occupation without
winding.  Only a returning `P_L -> P_L` process that changes a port total
must wrap; its minimum torus support is `2L_min`.  The source controls
dressing of local potentials/observables, not the norm distance of a global
dressed spectral subspace.

Global whole-band Feshbach/direct-rotation bounds have a radius that shrinks
with volume.  On each fixed torus, the first canonical matrix element around
the minimal alternating winding cycle is nonzero at order `2L_j`.  Thus
GL6AY is a rigorous finite-coupling **prethermal** theorem, not exact
all-time finite-coupling locked-phase closure, selected-GNS gravity, or `G`.
