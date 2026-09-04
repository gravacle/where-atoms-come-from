# Independent hostile audit — GL6CL global Fourier pair writer

**Target:** `LANE_CROSS_RFT_GRA_GL6CL_GLOBAL_FOURIER_PAIR_WRITER_V001`  
**Disposition:** **PASS**

The packet's exact Fourier, rank, and smooth-access results reproduce
independently.  The repaired scope is also correct: the complete
arbitrary-profile order-six writer is only the projected tensor map
`B_+ P_T` (and `B_- P_T` in the relative sector).  The unprojected
canonical-direct row is bookkeeping.  The locked diagonal read `D`, not an
unclassified off-diagonal writer, supplies `A1+E` in the six-direction
composition.

## 1. Independent method

The audit imports no target code.  Starting from the four tetrahedral link
vectors and the audited GL6CH tensor vertex, it independently:

- walks each six-link incidence cycle and derives its six node positions,
  center, and parent/child pair offsets;
- constructs the exact zero-mode incidence, locked-read, and projected
  writer matrices;
- expands the common symbol as sparse rational polynomials through the
  claimed orders;
- derives the relative leading matrix, all four `3x3` minors, and exact sine
  dependencies;
- evaluates the Brillouin-boundary rank-loss point algebraically in the
  basis `cos(pi/8), cos(3pi/8)`; and
- checks the scalar storage identity on all 16 four-bit words.

The independent result is frozen in `verify_gl6cl_independent.py` and
`INDEPENDENT_RESULT.json`.

## 2. Geometry and exact Fourier rows — PASS

The incidence walk yields four orientations and twelve child offsets
`rho_{d,p}`.  For each pair, the parent and child positions are opposite
about the derived ring center.  Every offset has

\[
 |\rho_{d,p}|^2={11\over4}.
\]

The aggregate moments are

\[
 \sum_{d,p}\rho\rho^T=11I,
 \qquad
 \sum_{d,p}\rho_i^4={83\over4},
 \qquad
 \sum_{d,p}\rho_i^2\rho_j^2={19\over4}\quad(i\ne j).
\]

Factoring the ring-center phase gives exactly

\[
 B_d^P=\sum_{p\subset\bar d}e^{-ik\cdot\rho_{d,p}}e_p^T,
 \qquad
 B_d^C=\sum_{p\subset\bar d}e^{+ik\cdot\rho_{d,p}}e_p^T,
\]

and hence the common cosine and relative sine rows.  The physical
arbitrary-profile statement inherited from GL6CH is the tensor projection

\[
 \delta a_d^T=\mu\,[B_d^+P_Tj_+ + B_d^-P_Tj_-],
 \qquad \mu={105\over8}{h^6\over U_d^6}.
\]

The audit found no basis, sign, phase, or normalization discrepancy.

## 3. Zero-mode composition and claim scope — PASS

At zero momentum the unprojected canonical-direct row satisfies

\[
 B_+(0)^*B_+(0)=24P_A+8P_T,
\]

so it has rank four on `A1+T2`.  This does **not** make its `A1` column a
complete nonuniform writer or its absent `E` column a physical null law.
GL6CH classifies only the `T2` off-diagonal completion at arbitrary profile.
The target now preserves that limitation throughout.

The complete tensor writer instead obeys

\[
 (B_+P_T)^*(B_+P_T)=8P_T,
\]

with rank three.  Independently enumerating the six locked degree-two words
gives

\[
 D^*D=4P_A+16P_E,
\]

with rank three and kernel `T2`.  Therefore the properly typed stack
`C=(D;B_+P_T)` has

\[
 C(0)^*C(0)=4P_A+16P_E+8P_T,
\]

rank six, determinant `524288`, and inverse normal

\[
 {1\over4}P_A+{1\over16}P_E+{1\over8}P_T.
\]

This confirms the advertised six-direction operator-access result without
using an unclassified `A1/E` off-diagonal writer.  For unsoldered parent and
child fields, the audit obtains rank nine at zero momentum and the general
dimension ceiling `rank<=3+3+4=10` from twelve inputs.  A common-field law is
therefore a genuine extra physical restriction, not an incidence theorem.

## 4. Smooth neighborhood and determinant — PASS

The exact canonical-direct determinant in the unnormalized `(A,t1,t2,t3)`
basis is

\[
 768-1408|k|^2+1072|k|^4
 -{416\over3}\sum_i k_i^4+O(|k|^6).
\]

This determinant is a useful bookkeeping diagnostic, not a completed
physical five-shear determinant.

Using `|2(cos x-1)|<=x^2`, the twelve equal-radius offsets give

\[
 \|B_+(k)-B_+(0)\|_2^2
 \le \|B_+(k)-B_+(0)\|_F^2
 \le {363\over4}|k|^4.
\]

The zero-mode tensor singular-value square is `8`, so the tensor rank and
the combined rank remain full on the explicit open ball

\[
 |k|^4<{32\over363}.
\]

The analytic left inverse `[C^*C]^{-1}C^*` is valid there.  This is a
finite-dimensional analytic rank theorem; the target correctly does not
call it a continuum limit or spacetime result.

## 5. Cubic tensor block and `SO(3)` ceiling — PASS

The independent quadratic expansion in an orthonormal `T2` basis is

\[
 N_T(k)=8I-2|k|^2I+12kk^T
       -28\operatorname{diag}(k_x^2,k_y^2,k_z^2)+O(|k|^4).
\]

The last term records cubic-lattice structure in the restricted block, but
cannot by itself diagnose physical rotational anisotropy.  The audit gives a
direct representation check: a 45-degree rotation about `z` sends the `xy`
off-diagonal shear into `diag(-1,1,0)`.  Thus the three-dimensional `T2`
space is not closed under `SO(3)`; it mixes with the two diagonal traceless
`E2` directions.  A consistently normalized full `E2+T2` operator is needed
before either isotropy or anisotropy can be concluded.  The target states
this boundary explicitly and makes no `SO(3)` promotion from CL23.

## 6. Relative sector — PASS

The derived leading relative map vanishes at `k=0` and has the exact four by
three coefficient matrix in CL27.  Its Cauchy--Binet minor sum is

\[
 9\sum_i k_i^6-9\sum_{i\ne j}k_i^4k_j^2
 +58k_x^2k_y^2k_z^2.
\]

It has rank three for generic direction and rank two at leading order on the
six nonzero Cartesian face-diagonal directions.  Direct sine-series algebra
also finds one exact dependent column pair along each of those six momentum
lines.  The exact dependency gives rank at most two; isolated special
momenta may lower it further.  This is consistent with the target's precise
leading-rank statement and does not create an arbitrary-support inverse.

## 7. Exact finite-momentum loss — PASS

At reciprocal coordinate `q=(pi,0,0)`, equivalently
`k=(pi/4)(1,1,1)`, the three complete common tensor columns are identical.
In the exact symbols `C=cos(pi/8)` and `S=cos(3pi/8)`, their rows are

```text
(-2C,-2C,-2C)
( 2S, 2S, 2S)
( 2S, 2S, 2S)
( 2S, 2S, 2S).
```

The tensor writer therefore has rank one at this point.  The smooth-access
theorem is correctly local in momentum and is not promoted to a full-zone
inverse.

## 8. Uniform scalar identity — PASS

All sixteen local words satisfy

\[
 \sum_{a<b}Z_aZ_b=2(n-2)^2-2.
\]

A uniform scalar source therefore sends `U_d` to `U_d+2q`, apart from the
constant `-2q`.  Differentiating the source-free ring coefficient gives

\[
 {d\over dq}\left[-{63\over8}{h^6\over(U_d+2q)^5}\right]_{q=0}
 ={315\over4}{h^6\over U_d^6},
\]

exactly equal to `6*(105/8)h^6/U_d^6`.  This establishes the claimed uniform
storage-energy/future-writer linkage.  It does not establish autonomous
source formation or reciprocal bulk dynamics.

## 9. Custody and verdict

The target science replay passes `122/122`; its packet/custody verifier
passes `129/129`.  All twelve imported GL6CH target/audit dependencies match,
the ten target payload hashes match, and the one-line seal authenticates the
manifest.  The audit separately pins those bytes and its own exact replay.

**PASS.**  GL6CL is a sound global operator-jet and smooth-access theorem.
It is not yet an autonomous field equation, response kernel, rotational
completion, continuum geometry, gravity theorem, or calculation of `G`.
