# GL6AY self-audit

## Primary-theorem attacks

1. **Single-site mismatch:** the F3 constraint terms overlap.  The theorem
   uses the explicit strong-support extension for commuting finite-radius
   summands, not the single-site theorem without repair.
2. **Integer spectrum:** every `q_v^2` has spectrum `{0,1,4}`; hence
   `N_def` has integer spectrum and `exp(2 pi i N_def)=1`.
3. **Local Hilbert and locality:** a coarse `A3` cell contains four link
   qubits and has dimension sixteen.  Each constraint star and each strongly
   supported flip term has volume-independent finite radius.
4. **Norm substitution:** the small parameter is stated in the exact
   exponential strong-support potential norm.  A global operator norm is
   never substituted for it.
5. **System size:** every quoted constant is independent of the finite
   volume.  No global norm estimate is mislabeled volume-uniform.
6. **Static frequency:** `nu=U_d` is the static large energy coefficient in
   the closed-system theorem, in its `hbar=1` convention.  It is not a
   periodic external drive.

## Effective-interaction attacks

1. **Charged resonances:** `D_0` is not set to zero globally.  Single flips
   can be resonant in charged sectors.  Only finite-volume
   `P_L D_0 P_L=0` is used.
2. **No infinite-volume global projector:** `P_L=chi(N_def=0)` is used only
   in finite volume.  Every infinite-volume term is instead restricted by the
   finite collar projector
   `P_S^0=product_(v:supp(q_v^2) subset S)chi(q_v^2=0)`.
3. **Termwise port conservation:** strong support and termwise `N_def`
   pinching imply `[D_hat(S),N_S]=0`.  Thus
   `Phi_S=P_S^0D_hat(S)P_S^0` maps globally locked configurations to globally
   locked configurations locally, and the sealed `GL6AX` endpoint theorem
   applies.
4. **Torus exception:** periodized wrapping supports are not declared port
   conserving.  They are placed in the explicit quasi-local tail.
5. **Second moment:** `P_S^0` is diagonal in link occupation and commutes
   with `A_S`, so the double commutator is the collar projection of the
   unprojected one.  Exponential support-size control beats its fourth power.
6. **Order-six coefficient:** agreement with `GL6AO` is coefficient
   matching in the normal-form expansion.  It is not a claim that the full
   power series converges.

## Remainder and topology attacks

1. Finite-volume `P_L V_hat_L P_L=0` does not make `V_hat_L` zero and is not
   converted into an infinite-volume global projection.
2. Local `P_L -> Q_L` leakage is not assigned a winding threshold.  A local
   charged excursion can change a bare port total.
3. The `2L_min` threshold is used only for locked-to-locked port-sector
   mixing.  A contractible excursion that returns locked has zero net port
   change even if its intermediate states are charged.
4. The theorem supplies an upper bound on the remainder, not a lower bound;
   it does not prove generic nonzero leakage from the bound alone.
5. The potential/local-observable dressing estimates are not used to bound
   `||Y_L^*P_LY_L-P_L||`; no volume-uniform global dressed-subspace closeness
   is claimed.

## Exact-closure attacks

1. The exact equality with `V_hat` retained is not called exact reduction to
   `D_hat`.
2. The extensive global-norm obstruction shows failure of the standard
   volume-uniform whole-band contour argument; it is not elevated into an
   impossibility theorem for every specially integrable normal form.
3. The finite-torus winding coefficient is a fixed-volume analytic fact.
   Its convergence radius is not claimed uniform in volume.
4. A prethermal horizon is not an all-time statement or an equilibrium
   phase theorem.
5. The conditional observation-horizon corollary leaves both `U_d/h` and
   physical clock calibration open.
6. No record observation, photon, graviton, Ricci tensor, Einstein equation,
   gravity law, or numerical `G` is inferred.
