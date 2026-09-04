# GL6CK theorem — the first two-overlap global-`H6` tensor response

## 1. Inputs and boundary

Let `c_0,c_1` be the selected period-four `Q4` alternating hexagons of
`GL6BY`.  They share the physical edge `((0,0,0),0)`.  With the exterior
frozen to the authenticated `BY` witness, their locked subspace is

```text
P_* = span{|B>,|A>,|C>},
A=B xor c_0,  C=B xor c_1,                         (CK.1)
```

and the double flip is excluded.  The order-six source-free operator on this
subspace, after removal of a common diagonal shift, is

```text
H_0=-J A_*,                  J=(63/8)h^6/U_d^5,
A_*=[[0,1,1],[1,0,0],[1,0,0]].                    (CK.2)
```

Only this geometry and locked-star topology are imported from `BY`.  In
particular, its finite-collar order-four `T2` vertex is boundary contaminated
and is **not** an input to this theorem.

The global order-six first-source law derived in `GL6CH` is, for an
alternating ring `c`,

```text
<sigma xor c|H_eff(j)|sigma>
 = -(63/8)h^6/U_d^5
   +(105/16)h^6/U_d^6 sum_(v in c) j_v.Theta_(v,c)
   +O(j_T^2,h^8).                                  (CK.3)
```

Here `Theta_(v,c)=e_ab-e_cd`, where `{a,b}` is the pair of cycle ports at
`v` and `{c,d}` is its complement.  Thus `Theta` is pure `T2` and
`||Theta||^2=2`.  The packet's independent 720-order replay recovers the
source-free coefficient `-63/8`, the full canonical pair gradient
`+(105/8)e_ab`, and its `T2` projection `+(105/16)Theta`.

## 2. Isolated-ring null

For one ring the locked graph is `K2`.  If `s` scales any fixed source
direction, (CK.3) gives

```text
H_K2(s)=(-J+s w)sigma_x.                           (CK.4)
```

The perturbation and source-free Hamiltonian commute.  The local ground
branch is `E_g(s)=-J+s w`; hence

```text
E_g''(0)=0.                                        (CK.5)
```

The field dependence is real, but in an isolated ring it only rescales the
same flip amplitude.  It creates no retained pole.

## 3. Two-overlap star theorem

On the ordered basis `(B,A,C)`, let `w_0,w_1` be the first-source changes of
the `B-A` and `B-C` amplitudes.  Then

```text
H_*(s)=
 [[0,-J+s w_0,-J+s w_1],
  [-J+s w_0,0,0],
  [-J+s w_1,0,0]].                                 (CK.6)
```

At `s=0`, define

```text
|g>=(sqrt(2)|B>+|A>+|C>)/2,
|m>=(|A>-|C>)/sqrt(2),
|u>=(sqrt(2)|B>-|A>-|C>)/2.                        (CK.7)
```

Their energies are `-sqrt(2)J,0,+sqrt(2)J`.  Writing

```text
wbar=(w_0+w_1)/2,   delta=(w_0-w_1)/2,             (CK.8)
```

the exact first-source matrix elements are

```text
<g|B|g>=sqrt(2)wbar,
<m|B|g>=delta,
<u|B|g>=0.                                         (CK.9)
```

The common part is again only amplitude rescaling.  The unequal part couples
the ground state to the relative arm mode across the gap `sqrt(2)J`.
Therefore the retained stationary spectral Hessian is

```text
E_g''(0)
 = 2 |<m|B|g>|^2/(-sqrt(2)J)
 = -2 delta^2/(sqrt(2)J)
 = -sqrt(2)(w_0-w_1)^2/(4J).                       (CK.10)
```

Equivalently, the exact branch

```text
E_g(s)=-sqrt[(-J+s w_0)^2+(-J+s w_1)^2]            (CK.11)
```

has precisely (CK.10) as its second derivative.  Consequently the spectral
response is nonzero exactly when `w_0 != w_1`.

## 4. Pure-`T2` witness and normalization

At the shared parent endpoint `P(0,0,0)`, the local cycle pairs are
`{0,2}` for `c_0` and `{0,1}` for `c_1`.  At the child endpoint they are
interchanged.  In either case,

```text
Theta_0.Theta_0=2,
Theta_1.Theta_1=2,
Theta_0.Theta_1=0.                                 (CK.12)
```

Put a source at one shared endpoint only and choose the literal pure source
direction `d j/ds=Theta_0`, where `s` has energy units.  Equation (CK.3)
gives the two amplitude derivatives

```text
w_0=(105/16)(h^6/U_d^6)(Theta_0.Theta_0)
    =(105/8)h^6/U_d^6,
w_1=(105/16)(h^6/U_d^6)(Theta_0.Theta_1)=0.         (CK.13)
```

Substitution into (CK.10) yields

```text
E_g''(0)=-(175sqrt(2)/32)h^6/U_d^7.                (CK.14)
```

This Hessian is with respect to the energy-valued coordinate `s` in the
literal direction.  For the unit direction
`j_hat=Theta_0/sqrt(2)`, it is

```text
E_g,unit''(0)=-(175sqrt(2)/64)h^6/U_d^7.           (CK.15)
```

For an arbitrary local `T2` source vector, the full spectral Hessian is the
rank-one tensor

```text
H_v^spectral
 = -sqrt(2)/(4J) [105h^6/(16U_d^6)]^2
   (Theta_0-Theta_1)(Theta_0-Theta_1)^T.            (CK.16)
```

## 5. Meaning and exact scope

Equations (CK.5) and (CK.10) separate mere writer existence from accumulated
response.  A single ring carries a candidate-field-dependent future writer
but cannot respond spectrally to it.  The smallest overlapping locked graph
has a relative collective mode; a source that weights its two real ring
channels differently excites that mode and produces an exact nonzero
stationary energy curvature.  This is the first exact two-overlap
composition of the universal global-`H6` writer into a retained `T2`
spectral response on this route.

The theorem is conditional on the selected finite `BY` locked star and on
the `CH` global writer law.  It does not establish a parent-selected
thermodynamic phase, an owner-once linked-cluster/bulk coefficient, the full
contact-plus-pole tensor, record authentication, a continuum solder map,
Ricci/Einstein form, gravity, or Newton's `G`.
