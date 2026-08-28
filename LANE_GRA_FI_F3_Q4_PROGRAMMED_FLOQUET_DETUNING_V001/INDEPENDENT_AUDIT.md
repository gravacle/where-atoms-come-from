# Independent hostile audit -- programmed F3/q4 Floquet detuning

**Lane:** `GRA-FI-F3-Q4-PFCD-V001`

**Date:** 2026-08-27

## Verdict

`ACCEPT_AFTER_TWO_NARROW_TYPING_AND_OWNERSHIP_REPAIRS__EXACT_PROGRAMMED_RESULT_ONLY`

The exact finite Floquet spectrum, its uniform principal-branch quasienergy
gap, the dressed-parent functional calculus, and the operator remainder bound
survive independent hostile recomputation.  Two nonconclusion-changing defects
were repaired: the zero-pulse value of the spectral isometry is the canonical
inclusion `iota_P:P -> P \oplus C`, not the ill-typed identity `I_P`; and the
absence of a new interaction is now scoped to the **bulk carrier generator**.
Concrete controller couplings, matrices, work, and calibration remain supplied
physical antecedents within the inherited port architecture.

The accepted result is therefore exact but conditional: F3 already owns the
onsite and incidence-gated hopping generators, while `FPSS` admits the finite
fixed-program control architecture and complete port types.  This packet adds
the specifically supplied `PROGRAMMED-NEXT-SLAB-DETUNING` schedule.  It does
not derive an autonomous controller or the static source-off `FD05` block.

## Exact spectrum recomputation

Let `B=B_N` and reduce its injective singular channel to a parent/child pair
with singular value `sigma>0`.  With

\[
 a=\eta\sigma,\qquad c=\cos a,\qquad s=\sin a,
\]

the two pulses give

\[
 U_F=
 \begin{pmatrix}1&0\\0&-i\end{pmatrix}
 \begin{pmatrix}c&is\\is&c\end{pmatrix}
 =\begin{pmatrix}c&is\\s&-ic\end{pmatrix}.          \tag{A1}
\]

Its trace and determinant are `c(1-i)` and `-i`.  Defining

\[
 \omega=\arccos\!\left({\cos a\over\sqrt2}\right)
\]

therefore yields exactly

\[
 \lambda_P=e^{i(\omega-\pi/4)},\qquad
 \lambda_C=e^{-i(\omega+\pi/4)}.                    \tag{A2}
\]

Because q4 incidence is injective, there are no unpaired parent modes.  Each
vector in `ker B^dagger` is fixed by the hop pulse and receives `-i` from the
child pulse, giving exactly `dim(C)-dim(P)` dark child modes.

Every q4 column has sum four and every row has sum at most four, so
`||B||<=4`.  The frozen ceiling `0<|eta|<=pi/16` gives `|a|<=pi/4`, hence

\[
 \omega\in[\pi/4,\pi/3],\quad
 \arg\lambda_P\in[0,\pi/12],\quad
 \arg\lambda_C\in[-7\pi/12,-\pi/2].                \tag{A3}
\]

The dark phase is `-pi/2`.  On the displayed principal branch the distance
from every parent-connected phase to every child phase is consequently at
least `pi/2`, proving the stated uniform quasienergy gap
`pi hbar/(2 T_F)`.  Quasienergy periodicity creates no competing shorter gap
inside this branch window.

## Dressed-parent pullback and kernel bound

The normalized parent-connected eigenvectors define an isometry
`W_P(eta):P -> P \oplus C` with `W_P^dagger W_P=I_P` and continuous limit
`W_P(0)=iota_P`, where `iota_P p=(p,0)`.  Since each eigenphase depends only
on `sigma^2`, pulling the principal logarithm back through this isometry gives
the exact spectral function

\[
 H_P^F=-{\hbar\over T_F}
 \left[\arccos\!\left({\cos(\eta\sqrt K)\over\sqrt2}\right)
 -{\pi\over4}I\right],\qquad K=B^\dagger B.          \tag{A4}
\]

For `z=eta^2 kappa in [0,pi^2/16]`, set

\[
 g(z)=\arccos\!\left({\cos\sqrt z\over\sqrt2}\right)-{\pi\over4}.
\]

Direct differentiation independently reproduces

\[
 g'(z)={\sin x\over2x\sqrt{1+\sin^2x}},\qquad
 g''(z)={x\cos x-\sin x(1+\sin^2x)
 \over4x^3(1+\sin^2x)^{3/2}},\quad x=\sqrt z.       \tag{A5}
\]

The numerator of `-g''` is
`(sin x-x cos x)+sin^3 x`.  On `0<=x<=pi/4`, its first term is the integral
of `u sin u` and is at most `x^3/3`, while the second is at most `x^3`.
The denominator in (A5) is at least `4x^3`.  Thus
`-1/3<=g''<=0`, including the continuous zero limit, and

\[
 {z\over2}-{z^2\over6}\le g(z)\le {z\over2}.        \tag{A6}
\]

Applying (A6) eigenvalue by eigenvalue to positive `K` proves exactly

\[
 0\preceq H_P^F+{\hbar\eta^2\over2T_F}K
 \preceq {\hbar\eta^4\over6T_F}K^2.                \tag{A7}
\]

With the independently established q4 identity `K=4I+A_N`, the first
nonconstant term is therefore `-(hbar eta^2/(2T_F))A_N`; the remainder is
retained rather than discarded.  This is the earned common-child Floquet
kernel.  It is not the static Schur Hamiltonian `FD05` and does not establish
a massless or gravitational phase.

## Inherited-physics and ownership audit

1. **F3 generator ownership.**  BS09 contains the uniform active-slab onsite
   term and incidence-gated content-blind carrier hopping.  In the original
   `V_0--V_1` slab its onsite contribution is a one-carrier scalar.  In the
   next `V_1--V_2` slab, `P=V_0` is absent, `C=V_1` is present, and blank
   `G=V_2` carries no occupation.  Conditional on exact isolation, the same
   generator therefore restricts to `epsilon_psi Pi_C` without a label-based
   stagger.
2. **Parameter domain.**  F3 admits `epsilon_psi>=0`; this theorem lawfully
   selects the narrower `epsilon_psi>0` domain.  It also requires `t!=0` and a
   calibrated pulse satisfying the stated `eta` ceiling.  None of those
   numerical values is predicted.
3. **Support ownership.**  `FPSS` physically binds the supplied finite q4
   append support to F3 site/link factors and permits the saturated old-slab
   incidence hold.  Blank next-slab incidence makes its BS09 hopping vanish,
   provided incidence flips and formation/copy terms are isolated as the
   supplied program requires.
4. **Control ownership.**  `FPSS` admits fixed orthogonal programming,
   switching/refocusing, and controller/clock/work/failure port types, but it
   does not derive the exact physical matrices for this new cycle.  The packet
   now states that those matrices and their calibration are supplied.  Thus
   no new **bulk** carrier generator or register/port type is introduced, but
   no control interaction is conjured for free from symbolic BS12.
5. **Static boundary.**  `CLDNG-3` remains intact.  With the program removed,
   the exchange-symmetric source-off F3 bulk has no child-only carrier onsite
   offset.  The result is programmed and stroboscopic, not autonomous or a
   retroactive derivation of `DETUNED-Q4-CARRIER-LIFT` as originally stated.

## Independent numerical replay

An audit-only replay used exact q4 incidence matrices for `N=0,...,8` and
both signs of `eta`, including both pulse-ceiling endpoints.  It found a
maximum unitarity residual `1.203e-14`, maximum spectrum-phase residual
`4.505e-15`, zero pullback residual at printed precision, and no detected
violation of either side of (A7).  The smallest sampled principal phase gap
was `0.500162636191970 pi`; the analytic result above supplies the uniform
limiting bound `pi/2`.

The packet verifier was then rerun after the repairs and reported:

```text
SUMMARY 42/42 checks passed
```

The accepted payload hashes before adding this audit were:

```text
7cec16096123a0846de8327cbcc00a7948d8669084124aeed40f920b22edb513  THEOREM.md
1492ef978e7fb966f918c43367dd176c73d9b4b7fe1e348f0c428c065c28901b  verify_programmed_floquet_detuning.py
0183907837387b6e25363140f67cfa36154e0669483a6d449da7b765c05298e1  SELF_AUDIT.md
```

## Final ceiling

The packet closes the finite **programmed** positive-detuning role with an
exact Floquet substitute assembled from inherited F3 bulk generators.  It
does not close autonomous schedule formation, static source-off detuning,
collective-phase stability, physical time/refinement binding, pair-field
dynamics, tensor constraints, universal stress coupling, RGRL-B, gravity, or
the value of `G`.
