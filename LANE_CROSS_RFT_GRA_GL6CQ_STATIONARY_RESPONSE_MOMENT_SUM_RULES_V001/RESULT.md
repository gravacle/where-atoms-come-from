# GL6CQ Result

For the bare connected stationary cycle susceptibility

\[
 K^{\rm bare}_{dd'}(R)=2\operatorname{Re}
 \langle0|T_{0d}Q(H_0-E_0)^{-1}QT_{Rd'}|0\rangle,
\]

define

\[
 Z_{dd'}=\sum_RK^{\rm bare}_{dd'}(R),\qquad
 M_{dd'}^{mn}=\sum_RX_{R;dd'}^mX_{R;dd'}^nK^{\rm bare}_{dd'}(R).
\]

Absolute second-moment convergence and inversion give

\[
 \widehat K(k)=Z-\frac12k_mk_nM^{mn}+o(|k|^2).
\]

With the inherited orthogonal cycle solder, the six GL6CO coefficients are
exact observables:

\[
 \kappa={Z_T\over3},\quad
 \alpha=-{1\over6}\sum_mM_{AA}^{mm},\quad
 \eta=-{1\over3}(M_{AT,x}^{yz}+M_{AT,y}^{zx}+M_{AT,z}^{xy}),
\]

\[
 b=-{M_\perp\over12},\qquad
 c=-{M_\parallel\over6}+{M_\perp\over12},\qquad
 d=-{M_\times\over6}.
\]

The writer normalizations obey

\[
 \mu={105\over8}{h^6\over U_d^6}=2\lambda_T,
\]

and `K^bare` contains no writer coefficient.  The spectral pullback in the
GL6CL common coordinate is `mu^2 B_T^* K^bare B_T`; in the orthonormal
parent/child common source used by GL6BV it is
`(mu^2/2) B_T^* K^bare B_T`.  Combining the latter with the contact scale
`g_ct=h^2/(4U_d^3)` and the same-state probability
`p=<Pi_same>` turns CO29 and CO30 into

\[
 \boxed{{\mu^2\over2}[-2Z_T+M_\perp-2M_\parallel-2M_\times]
       +4g_{\rm ct}(2p-1)=0,}
\]

\[
 \boxed{-{\mu^2\over2}[Z_T+M_\perp]+2g_{\rm ct}(1-4p)=0.}
\]

These equations are exact if-and-only-if tests for the quadratic-gradient
shape of the explicitly defined contact-plus-two-writer sector.  They do
not assert that the left sides vanish, remove the zero-momentum response,
complete other source-second vertices, invert the response to a 1PI kernel,
or establish Ricci, gravity, or `G`.  If the second moment diverges at a
critical point, the analytic `k^2` test ceases to exist and the leading
nonanalytic kernel must be measured directly.
