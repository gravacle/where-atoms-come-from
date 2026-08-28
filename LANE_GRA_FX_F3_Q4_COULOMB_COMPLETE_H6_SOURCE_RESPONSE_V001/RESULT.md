# Result: complete homogeneous H6 source response

Under FV `S1`--`S10`/`FV-PURE`, exact enumeration of all diagonal closed
flip words and all BW/Feshbach folds through H6 gives

\[
 a_2=-60,\qquad a_4=-35,\qquad a_6=-893/9.
\]

After applying `Q=-2 dH/dj`, restoring `x=h/U_d`, and using
`J_6=(63/8)U_dx^6`, the previously open `Q_diag^(2,4,6)` terms reduce exactly
to the direct pair source plus identities.  The complete nonidentity
diagonal coefficient is

\[
 f_E(x)=1-x^2-{37\over12}x^4-{16247\over900}x^6,
 \qquad \rho_E={8f_E(x)\over63x^6}.
\]

Consequently, because `f_E(x)` is a formal power-series unit, the generic
finite homogeneous response hierarchy remains

\[
 \boxed{5\ \to\ 3\ \to\ 2\ \to\ 2}
\]

for operator rank modulo identity, `ad_H` rank, ground retarded rank, and
`M_1` rank.  The two FW poles remain at `2+2 sqrt(2)` and `4+2 sqrt(2)`;
only the first residue's `E` amplitude changes from `rho` to `rho_E`.

The evaluated finite through-H6 cancellation occurs at the unique positive
root `x=0.5398271903...`, where the truncated ranks are `4 -> 2 -> 2 -> 2`.
It is not a physical threshold: H8+ is uncontrolled there and may move or
remove it.

The exact response statement follows from the exact pair-plus-identity
reduction composed with audited FW algebra.  NumPy 180-state spectra are
independent replays.  The theorem is complete only for the homogeneous
source on the selected FO 180-state winding component, through H6, modulo
`H_id`.  Nonzero momentum, other components, CTP/Ward closure,
thermodynamic response, RGRL-B, gravity, `G`, and Newton's constant remain
outside this result.
