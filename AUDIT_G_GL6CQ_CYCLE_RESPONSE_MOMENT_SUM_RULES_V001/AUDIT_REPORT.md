# Independent hostile audit — GL6CQ

## Disposition

**PASS**, within the exact analytic, same-state, and sector ceilings stated
by the target.  No stale common-source normalization remains.

The target's author program was neither imported nor executed.  A separate
exact implementation rebuilt the raw orientation-space moment tensors,
projected them through an independently constructed orthogonal cycle solder,
reconstructed the contact coefficients, and derived both boxed observable
equations.

## 1. Frozen target and dependency custody

The audited target is the packet sealed by manifest digest

```text
bc0c979698a0738db9c841bfabd718023940fe3a43d397f04332bbe5009fe469
```

`TARGET.sha256` pins its ten payload bytes, manifest, and seal.  All 27
upstream dependency hashes resolve exactly, including the repaired GL6CO
bytes whose manifest digest is
`f085a73c4d7590a44ac89117f53bdb00583646153396a16552a84829ebd323b6`.
The target seal closes its manifest, and every pinned upstream seal closes
the corresponding upstream manifest.

## 2. Independent moment projection

The audit independently forms

\[
u_d={1\over2},\qquad Q_{di}={1\over2}(T_d)_i,
\qquad S=(u\;Q),\qquad S^TS=I_4.
\]

It starts from a generic invariant `A1+T2` quadratic symbol, differentiates
it twice with respect to momentum, converts the result to raw cycle-
orientation coordinates with `S`, and only then projects the raw moments
back with `u` and `Q`.  All 144 projected second-moment components are
checked for each of four independent rational coefficient assignments.  The
contractions give

\[
 Z_T=3\kappa,\qquad M_\perp=-12b,\qquad
 M_\parallel=-6(b+c),\qquad M_\times=-6d,
\]

and hence

\[
\boxed{
 \kappa={Z_T\over3},\quad
 \alpha=-{1\over6}\sum_mM_{AA}^{mm},\quad
 \eta=-{M_{AT,x}^{yz}+M_{AT,y}^{zx}+M_{AT,z}^{xy}\over3},}
\]

\[
\boxed{
 b=-{M_\perp\over12},\qquad
 c=-{M_\parallel\over6}+{M_\perp\over12},\qquad
 d=-{M_\times\over6}.}
\]

The signs follow directly from
`Khat(k)=Z-(1/2) k_m k_n M^{mn}+o(k^2)`.  In particular, the two copies of
an off-diagonal monomial in the Einstein sum are correctly accounted for;
there is no missing factor of two in `eta` or `d`.

## 3. Independent normalization challenge

For a pure `T2` pair source and a complementary pair `pbar`,

\[
 j\cdot\Theta_p=j_p-j_{\bar p}=2j_p,
\]

so the two established vertex conventions obey

\[
 \mu=2\lambda_T,
 \qquad \lambda_T(j\cdot\Theta_p)=\mu j_p.
\]

The GL6CL sublattice convention is

\[
 j_P=j_++j_-,\qquad j_C=j_+-j_-,
 \qquad j_+={j_P+j_C\over2}.
\]

Its common embedding has Gram two.  Therefore the orthonormal coordinate is
`jhat_+=(j_P+j_C)/sqrt(2)=sqrt(2)j_+`, and

\[
 \delta a={\mu\over\sqrt2}B_T\widehat j_+,
 \qquad
 \widehat{\cal H}^{H6}_T={\mu^2\over2}B_T^*K^{bare}B_T.
\]

The target uses `mu^2` only for the unnormalized GL6CL coordinate in CQ20
and CQ21; it uses `mu^2/2` for the normalized contact comparison in CQ20a,
CQ25--CQ28, the result, and the ledger.  The bare susceptibility contains no
writer factor, and the contact carries `g_ct` with no `mu`.  Thus neither a
stale factor-two error nor a hidden second writer factor remains.

## 4. Contact and the two observable equations

The audit reconstructs all four tetrahedral defect projectors in the
orthonormal pair-`T2` basis and verifies

\[
 \sum_a(k\cdot T_a)^2Q_a={4\over3}|k|^2I+{8\over3}O(k).
\]

For four independent probabilities `p`, this gives

\[
 A_{ct}={4\over3}(1-4p),\qquad B_{ct}=0,
 \qquad C_{ct}={8\over3}(2p-1).
\]

Substitution of the independently recovered moment formulas into the GL6CO
extension mismatch gives

\[
 -4\kappa+8(c+d)
 ={2\over3}[-2Z_T+M_\perp-2M_\parallel-2M_\times].
\]

Multiplying the complete normalized condition by `3/2` therefore yields

\[
\boxed{
 {\mu^2\over2}[-2Z_T+M_\perp-2M_\parallel-2M_\times]
 +4g_{ct}(2p-1)=0.}
\]

Likewise,

\[
 -2\kappa+8b=-{2\over3}(Z_T+M_\perp)
\]

gives the held-out reference-shape test

\[
\boxed{
 -{\mu^2\over2}[Z_T+M_\perp]+2g_{ct}(1-4p)=0.}
\]

The audit checks both equivalences on independent rational observable,
writer-scale, contact-scale, and probability assignments.  Both equations
are exact if-and-only-if tests for the displayed quadratic-gradient sector;
they are not claims that a physical state makes either left side vanish.

## 5. Analytic boundary and claim screen

- Entrywise absolute second-moment convergence is sufficient for the stated
  twice-differentiable expansion with `o(k^2)` remainder.
- A finite fourth absolute moment plus inversion supports `O(k^4)`; an
  exponential moment supports analyticity near zero.
- If the second moment diverges, the analytic coefficients `b,c,d` and these
  sum rules do not exist.  The leading nonanalytic kernel must instead be
  tested directly.
- The state supplying `K^bare` must be the same state supplying
  `p=<Pi_same>`.
- The zero-momentum tensor response is not removed or interpreted.
- Additional source-second blocks, if found, must be added explicitly.
- No phase, background-stationarity, masslessness, gauge null, 1PI kernel,
  Ricci/Einstein endpoint, gravity, or `G` is established here.

No material defect was found.
