# GL6BB result

Current sealed custody does not select any one of the three inputs required
for a numerical `GL6BA` mission:

\[
 (R,\sigma_{\rm obs},\omega_L).
\]

The values `R=2` and `R=5/2` are two exact admitted H6 witnesses, not a
selection; `sigma_obs` has neither a numerical value nor a finite upper bound;
and the prepared-blank state is an available authenticated branch rather than
the selected postformation state.  With the state unrestricted, the exact
binary interval is sharply `[0,1]` even after `R` and `sigma_obs` are fixed.

For any supplied uncertainty set of selected-mission triples, the strongest
pointwise interval inherited directly from `GL6BA` is

\[
 \left[
 \max\{0,\inf(q_L-\varepsilon_L)\},
 \min\{1,\sup(q_L+\varepsilon_L)\}
 \right],
\]

with each collar value paired with its own exterior error.  A trace-distance
state radius `eta` adds exactly `eta` to the corresponding binary-probability
half-width.

For the already defined prepared-blank branch, the radius-zero collar reduces
exactly to a five-state Dicke Hamiltonian.  It gives

\[
 q_0^{\rm blank}(R,s)
 =|a_0|^2+\tfrac12|a_1|^2+\tfrac13|a_2|^2
  +\tfrac12|a_3|^2+|a_4|^2,
 \qquad a=e^{-isK_R}e_0,
\]

and energy conservation proves the all-time collar bound

\[
 q_0^{\rm blank}(R,s)\ge1-{1\over3R}.
\]

Thus the complete finite mission conditionally obeys

\[
\begin{aligned}
R=2:&\quad p_+^\Omega\in
[\max\{0,11/6-e^{96|\sigma_{\rm obs}|}\},1],\\
R=5/2:&\quad p_+^\Omega\in
[\max\{0,28/15-e^{120|\sigma_{\rm obs}|}\},1].
\end{aligned}
\]

No physical value of `sigma_obs` is manufactured.  A same-clock value or
finite interval for it is the one new datum needed to evaluate this
prepared-blank two-member result.  One actual single-member mission still
requires the complete authenticated `(R,sigma_obs,omega_L)` tuple.
