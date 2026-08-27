# Canonical-\(U(1)\) alternative-alpha record-world theorem

**Claim class:** exact ideal canonical-\(U(1)\) finite-mission EFT model;
active-EM `REC` and `FCLPD_W` construction; exact non-singleton
record-forming-alpha set

**Not claimed:** a second complete universe with dynamical charged matter,
atoms, chemistry, gravity, or a UV completion

## 1. The decisive sets must remain distinct

For a frozen electromagnetic comparison context \(\chi\) and record packet
\(\Pi\), define

\[
 {\cal A}_{\rm RF}^{\rm EM,toy}(\Pi;\chi)
 :=\{a\ge0:\text{the frozen canonical-}U(1)\text{ parent at }\alpha=a
 \text{ contains the constructed active-EM `REC` and `FCLPD` episode}\}.
                                                               \tag{1}
\]

This is not the stronger set

\[
 {\cal A}_{\rm RF}^{\rm full}
 :=\{a:\text{a complete self-consistent universe at }\alpha=a
 \text{ has dynamical charged matter, stable structures, and records}\}. \tag{2}
\]

Let \(I_\chi\subseteq[0,\infty)\) be the prospectively frozen alpha domain on
which this ideal parent family is asserted. No conclusion is extrapolated
beyond \(I_\chi\).

Two members of (1) exactly refute selection of alpha by present recordhood even
when electromagnetism is load-bearing.  They do not settle the cardinality of
(2).

## 2. Canonically typed electromagnetic parent

Work in rationalized natural units at one frozen renormalization context:

\[
 {\cal L}_{\rm EM}
 =-{1\over4}F_{\mu\nu}F^{\mu\nu}+eJ^\mu A_\mu,
 \qquad e=\sqrt{4\pi\alpha}.                            \tag{3}
\]

The photon kinetic term is canonically normalized.  The physical eigenmode,
unit-current normalization, scale, scheme, cavity, boundary condition, and all
non-alpha controls are fixed before comparing parents.  Changing \(e\) therefore
changes the invariant interaction strength; it is not a field-coordinate or
charge-unit relabeling.

Take a bounded ideal cavity and a normalized transverse mode
\(\mathbf u(\mathbf x)\) of frequency \(\omega>0\):

\[
 \mathbf A_T(\mathbf x)
 =\sqrt{\hbar\over2\omega}
 \left[\mathbf u(\mathbf x)a+
 \mathbf u^*(\mathbf x)a^\dagger\right].               \tag{4}
\]

Freeze a real smooth pulse \(f\in C_c^\infty(0,T_W)\) and a smooth spatial
current \(\mathbf j\) satisfying

\[
 \nabla\!\cdot\mathbf j=0,
 \qquad
 \mathbf n\!\cdot\mathbf j|_{\partial V}=0.            \tag{5}
\]

The authorized source vertex receives \(x\in\{0,1\}\) and emits a
classical-central physical current/program register carrying

\[
 J_x^\mu(t,\mathbf x)
 =\bigl(0,x f(t)\mathbf j(\mathbf x)\bigr),
 \qquad \partial_\mu J_x^\mu=0.                        \tag{6}
\]

Thus \(e\int J\cdot A\) is gauge-consistent.  Choose a nonzero overlap with the
declared mode.  In the exact one-mode toy parent, all other modes have zero
overlap by construction. The theorem below uses that one-mode parent. A
multimode extension additionally requires an integrable \(\ell^2\)-valued drive,
a finite Magnus phase, a square-summable final displacement, and every retained
mode in the census.

The current is an explicit physical source/boundary port.  Its microscopic
charged-matter production, recoil, and battery are not derived by this ideal
parent.

## 3. Exact driven write

The electromagnetic write module applies one arm-common physical rule

\[
 J\longmapsto U_e[J].                                  \tag{7a}
\]

For the input register value \(J_x\), its interaction-picture Hamiltonian is

\[
 H_I^{(e)}[J_x](t)
 =-xe\,[\ell(t)a+\ell(t)^*a^\dagger],                  \tag{7}
\]

where \(\ell\) is fixed entirely by the current pulse, cavity mode, and
frequency.  Since \([H_I(t),H_I(s)]\) is a scalar multiple of the identity, the
Magnus series terminates after its scalar phase term. With the convention

\[
 D(\beta):=\exp(\beta a^\dagger-\beta^*a),             \tag{8a}
\]

the complete current-controlled rule is

\[
 U_e^{I}=\sum_{x=0}^1
 |J_x\rangle\!\langle J_x|_{\rm cur}\otimes
 e^{i\chi_{x,e}}D(xeb),
 \qquad
 b={i\over\hbar}\int_0^{T_W}\ell(t)^*dt,
 \qquad b\ne0.                                         \tag{8}
\]

For the exact one-mode theorem define

\[
 B^2:=|b|^2>0.                                         \tag{9}
\]

Both source arms begin with the identical prospective photon blank.  Ignoring
the irrelevant scalar phase, the post-write carrier states are

\[
 x=0:\ |0\rangle,
 \qquad
 x=1:\ |\beta_e\rangle,
 \qquad
 \beta_e=eb.                                            \tag{10}
\]

The alpha-bearing field is therefore on the only event-to-carrier path.  It is
not a decoupled spectator.

## 4. Exact writer-off hold and pre-read distinction

Because \(f(t)=0\) for \(t\ge T_W\), set \(G=T_W\).  After \(G\), isolate the
source/current generator and apply only free field evolution for a fixed
\(T_H>0\).  Each coherent amplitude gains a phase, so

\[
 |\beta_e(T_H)|^2=e^2B^2.                              \tag{11}
\]

The distinction already present in the carrier before reading is

\[
 D_{\rm tr}\!\left(
 |0\rangle\!\langle0|,
 |\beta_e(T_H)\rangle
 \!\langle\beta_e(T_H)|\right)
 =\sqrt{1-e^{-e^2B^2}}.                                \tag{12}
\]

This prevents the terminal query from defining or creating recordhood.

## 5. Fixed label-blind read

Let

\[
 P_0=|0\rangle\!\langle0|,
 \qquad P_1=I-P_0.                                     \tag{13}
\]

Couple a blank pointer qubit \(A\) with the fixed unitary

\[
 U_Q=P_0\otimes I_A+P_1\otimes X_A.                    \tag{14}
\]

The complete terminal pointer instrument retains \(Y=0,1\).  It has no source
arm, current-generator, writer-control, or outcome-dependent setting port.  Its
exact output laws are

\[
 \Pr(Y=1\mid x=0,e)=0,
 \qquad
 \Pr(Y=1\mid x=1,e)=1-e^{-e^2B^2},                    \tag{15}
\]

and therefore

\[
 \boxed{D_{\rm TV}(Y_1,Y_0)=1-e^{-e^2B^2}
 =1-e^{-4\pi\alpha B^2}.}                              \tag{16}
\]

Equation (16) is strictly positive for \(e>0\), strictly increasing in \(e\)
at fixed controls, and exactly zero at \(e=0\).  Removing the electromagnetic
write therefore removes this declared photon record:

\[
 e=0\Longrightarrow
 U_{0}^{I}=I,\quad D_{\rm tr}=D_{\rm TV}=0.            \tag{17}
\]

This ablation does not say that the source controller loses its original copy
of \(x\), only that no record forms in the nominated photon lineage.

## 6. REC theorem

Let \(r_e\) be the derivative photon-record episode with the frozen ancestry

```text
joint preparation
  -> authorized source arm x
  -> conserved physical current/program register J_x
  -> arm-common canonical rule J -> U_e[J]
  -> physical closure/cut G with isolated source/current exhaust
  -> free positive hold
  -> fixed photon-pointer query
  -> complete Y outcome
```

The preparation also supplies the common photon vacuum, clocks, source/writer
program, pointer blank, and fixed query program. The module at \(G\) is a
declared physical closure/cut: it isolates and terminates the post-write
event-dependent source, current-register, and controller outputs. Arm-common
clock and hold context may feed the closure/hold. The written photon is the only
event-dependent query-reaching write-side carrier crossing \(G\); event-bearing
exhaust routes are censused but terminate without reaching the query. Assume
this is the complete source-to-query ancestry of the ideal toy parent.

Then for every \(e>0\):

1. all physical roles, controls, graph, outcomes, and the photon lineage are
   nominated independently of the contrast;
2. both arms begin in the same exact photon blank;
3. the conserved-current interaction creates (12);
4. the source is exactly off and isolated after \(G\);
5. a strictly positive common hold preserves the distinction;
6. the fixed complete query reads a pre-existing distinction and has no label
   port; and
7. every event-dependent source-to-query route crosses the electromagnetic
   write; censused event-bearing exhaust routes terminate outside the query
   ancestry.

These are the finite-mission physical record clauses, so at zero mathematical
relevance floor

\[
 \boxed{e>0\Longrightarrow\operatorname{REC}(r_e).}     \tag{18}
\]

For a frozen physical relevance floor \(0\le\delta<1\), the exact pass condition
is

\[
 1-e^{-e^2B^2}>\delta
 \quad\Longleftrightarrow\quad
 \boxed{
 \alpha>{-\ln(1-\delta)\over4\pi B^2}.}                \tag{19}
\]

## 7. FCLPD and DCL theorem

The same episode satisfies `FCLPD_W`:

- the displayed module graph is finite, acyclic, outcome-independent,
  mission-separated, and complete by the toy-parent definition;
- the preparation supplies one normalized joint source/field/control/pointer
  state;
- the smooth pulse contracts exactly to the unitary displacement (8), while
  hold and query are unitary;
- the terminal two-outcome pointer measurement is a complete CP instrument;
- only the declared source vertex receives \(x\); it emits \(J_x\), and the
  electromagnetic module applies the one common rule (7a); and
- every later rule is arm-common and receives source dependence only through
  the displayed physical edges.

The oscillator dimension is infinite, but FCLPD requires a finite exact module
graph, not a finite-dimensional Hilbert space. Use the Fock representation and
the unital algebra

\[
 C(\{J_0,J_1\})\ \bar\otimes\ B({\cal F})\ \bar\otimes\
 B(\mathbb C_A^2),                                     \tag{19a}
\]

together with the finite control/clock registers. The controlled displacement,
free hold, vacuum projector, pointer unitary, and complete query instrument all
belong to this declared algebra.

The existing FCLPD class theorem therefore gives, with the relevance condition
kept separate,

\[
 \boxed{\operatorname{FCLPD}_W(r_e)\quad\text{for every asserted }e
 \text{ with }e^2/(4\pi)\in I_\chi,}                   \tag{20a}
\]

\[
 \boxed{
 1-e^{-e^2B^2}>\delta
 \Longrightarrow\operatorname{REC}(r_e),}              \tag{20b}
\]

\[
 \boxed{
 \operatorname{FCLPD}_W(r_e)
 \Longrightarrow DCL_{\rm phys}(r_e).}                 \tag{20c}
\]

The simpler \(e>0\Rightarrow\operatorname{REC}(r_e)\) statement is reserved for
the ideal \(\delta=0\) packet in (18).

No U-DCL postulate is used in this per-record construction.

## 8. Two alpha parents with identical controls

Freeze the one pulse so that

\[
 B^2={100\over9}\ln2.                                  \tag{21}
\]

Assume the frozen validity domain \(I_\chi\) contains the two displayed values.
Compare the same parent family and apparatus at

\[
 e_1={3\over10},
 \qquad \alpha_1={9\over400\pi},                       \tag{22}
\]

and

\[
 e_2={3\over5},
 \qquad \alpha_2={9\over100\pi}.                      \tag{23}
\]

Then

\[
 D_{\rm TV}^{(1)}=1-e^{-\ln2}={1\over2},
 \qquad
 D_{\rm TV}^{(2)}=1-e^{-4\ln2}={15\over16}.            \tag{24}
\]

Thus any one common frozen floor \(\delta<1/2\) admits both parents as exact
`REC`, `FCLPD_W`, and hence `DCL_phys` episodes.  Current, pulse, cavity,
duration, hold, read, and thresholds are unchanged.  Their different response
at fixed controls also proves that \(\alpha_1\) and \(\alpha_2\) are physically
inequivalent, not conventions.

Consequently

\[
 \boxed{
 \left|{\cal A}_{\rm RF}^{\rm EM,toy}(\Pi;\chi)\right|\ge2.} \tag{25}
\]

More strongly, for the ideal packet with floor \(\delta\),

\[
 \boxed{
 {\cal A}_{\rm RF}^{\rm EM,toy}(\Pi;\chi)
 =I_\chi\cap
 \left({-\ln(1-\delta)\over4\pi B^2},\infty\right).}   \tag{26}
\]

If a vacuum-subtracted mean photon-excitation-energy ceiling \(E_{\max}\) is
part of the one-mode packet, then

\[
 \Delta E=e^2\hbar\omega B^2,
 \qquad
 \alpha\le {E_{\max}\over4\pi\hbar\omega B^2},         \tag{27}
\]

and the set in (26) is additionally intersected with that closed upper bound,
provided the resulting interval is nonempty. In the admitted multimode
extension the corresponding formulas are

\[
 \Delta E=e^2\sum_k\hbar\omega_k|b_k|^2,
 \qquad
 \alpha\le
 {E_{\max}\over4\pi\sum_k\hbar\omega_k|b_k|^2}.        \tag{27a}
\]

These are mean-energy constraints, not hard energy-support cutoffs; a coherent
state has an unbounded number and energy tail.

## 9. Decisive logical consequence

Equation (25) exactly refutes

\[
 \operatorname{REC}
 +\text{the declared canonical-}U(1)\text{ cavity interaction is
 load-bearing on the nominated photon path}
 \Longrightarrow \alpha=\alpha_*                       \tag{28}
\]

within the current finite-mission RFT model class.  Record formation does not
numerically select alpha even after the spectator loophole and control-retuning
loophole are removed.

This does not conflict with SAI or AWAI.  Across the two possible parents alpha
differs; within each parent, every same-sector record still inherits that
parent's one alpha.

## 10. Exact ceiling

The construction supplies an exact ideal active-EM record world, not a complete
cosmos.  It assumes:

- an externally supplied smooth conserved current with fixed normalization;
- an ideal one-mode cavity, or the stronger stated multimode drive/phase
  premises;
- an exact vacuum/nonvacuum pointer instrument; and
- complete toy-parent source, mode, boundary, controller, and detector custody.

It does not derive the current's charged matter, recoil, battery, or dependence
on alpha.  A coherent state has finite mean energy but unbounded number support,
so a hard occupation cutoff would make the model approximate.  A microscopic
detector built from the same varying-alpha matter is not supplied.  The theorem
does not establish charge quantization, anomaly cancellation, radiative/RG
consistency, stable atoms, chemistry, gravity, empirical realization,
actualization, U-DCL, or the viability of a complete alternate universe.

Therefore the scientifically complete fork remains

\[
 |{\cal A}_{\rm RF}^{\rm full}|>1
 \quad\text{(contingent full-universe alpha)}
\]

versus

\[
 {\cal A}_{\rm RF}^{\rm full}=\{\alpha_{\rm obs}\}
 \quad\text{(full-universe consistency selects alpha)}. \tag{29}
\]

Failure to construct another complete universe would not prove the singleton;
that conclusion requires an exhaustive no-go for every other value.
