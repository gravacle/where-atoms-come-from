# Independent hostile audit -- finite programmed q4-to-F3 support solder

**Lane:** `GRA-FH-F3-Q4-FPSS-V001`

**Date:** 2026-08-27

## Verdict

`ACCEPT__EXACT_FINITE_PROGRAMMED_SOLDER_WITH_CORRECT_PHYSICAL_CEILINGS`

The final frozen theorem, verifier, verification capture, and self-audit pass
independent hostile review.  The packet proves a lawful finite fixed-program
construction, not autonomous support emergence.  No remaining material defect
was found.

## Independent result

For fixed `N` and cap `R >= max(4,N+1)`, the census is exact:

\[
 |S_N|={N+3\choose3},\qquad
 |S_{N+1}|={N+4\choose3},\qquad
 g_N={N+3\choose2}.
\]

Choosing the common F3 layer width `M=|S_(N+1)|` gives `M^2` possible
cross-layer links.  The `4|S_N|` append keys are distinct, leaving exactly
`M^2-4|S_N|` nonedges.  The supplied injective address maps therefore fit the
whole q4 slab into the already admitted equal-width F3 allocation without
identifying labels or deleting guard factors.

The formation, KEEP, and controlled-link stages are reversible on the complete
retained register set.  Edge writers use distinct endpoint-owned slots; the
FPMH link factor is reused literally as the F3 incidence factor; nonedges keep
explicit blank `L/K/G/n` factors.  The final pulse statement correctly
requires the raw link flip and every noncommuting incidence detuning,
degree-return, current, and incidence-gated hopping term to be off or exactly
refocused.  Under that isolated schedule,

\[
 e^{i\pi P_e^KX_{n_e}/2}
\]

is identity on `K_e=0` and maps `|K_e=1,n_e=0>` to
`i|K_e=1,n_e=1>`.  No simultaneous omitted term is hidden in the claimed
preparation unitary.

## Quarantine, history, and selection boundary

With formation and route operations ended, `K_e=0,n_e=0` nonedges are
invariant: the gated flip vanishes, diagonal incidence terms preserve the
blank, and carrier/current terms carry the zero incidence projector.  Guard
carriers are likewise invariant when blank because no occupied incident link
or formation coupling reaches them.  On the fixed orthogonal support program,
the carrier hold is identity on retained BQ4 and compiler history factors.

This is passive memory, not selection.  ASSC conserves every support word,
including wrong, dense, fragmented, and non-q4 words.  The address map, q4
edge list, coexisting F3 hardware, source tokens, schedule, and physical port
realization remain supplied.  The packet also correctly excludes an arbitrary
coherent isometry from one BQ4 count-front factor into many coexisting F3
sites.  It removes a finite hardware-type obstruction only.

The port claim is properly limited to logical and custody completeness.  The
theorem explicitly leaves physical energies, port matrices, calibration,
timing, work, heat, recoil, and actual apparatus realization supplied; the
symbolic F3 port slot is not promoted into a calibrated physical completion.

## FD and FE claim audit

For the FD lane, saturated eligible incidence and `lambda_J=0` give exactly
the off-diagonal one-carrier block

\[
 -t\begin{pmatrix}0&B_N^\dagger\\ B_N&0\end{pmatrix}.
\]

Blank guards decouple.  This closes the finite programmed site/edge portion
only.  The positive uniform child/parent detuning, its owned maintenance
ports, collective phase, and scaling law remain open.

For the FE lane, the raw finite slab cannot have a global `d_*=2` sector.  The
child `c_*=(N+1,0,0,0)` has only one eligible parent, so every incidence word
satisfies `d_(c_*) <= 1` and therefore

\[
 \Omega_2(E_N)=\varnothing.
\]

The earned FE result is consequently finite edge binding and deep-interior
local support/operator inheritance.  A regular boundary or periodic diamond
completion is separately supplied physics, not a BQ4 consequence.  FD
saturation has parent degree four whereas FE ice requires degree two, so the
same-`n` incompatibility remains exact.  No `K_eT_e` interaction, second
support field, visible electromagnetic identification, tensor response, or
gravity theorem is introduced.

## Deterministic and packaging checks

The independent verifier rerun reports:

```text
SUMMARY 37/37 checks passed
```

Fresh verifier output is byte-for-byte identical to `VERIFICATION.txt`.  The
forbidden-promotion check normalizes both its needles and theorem text before
comparison.  Display mathematics and Markdown delimiters pass, and all five
declared dependencies exist.

Before this audit file was added, the frozen payload hashes independently
matched their four-entry manifest:

```text
40f70b76dd1b9ab32c2c47cece371cd9bf97247f18073e283c0a898f50b947e6  THEOREM.md
6cf264db1b39cd1b99d3f6969bb127b392d9a5fe5cdfcc6bb331a8eabede877d  verify_finite_programmed_support_solder.py
52e13148201d7b900a4eda60a8ec810af60cfc962982328ee069a3833becde6d  VERIFICATION.txt
31ff1302e0d9f9e015fc78e3392cc21cac09bf2d280a36a12e15b5a979988b98  SELF_AUDIT.md
```

The independent audit was added only after those four payload bytes were
accepted; none of them was modified during this audit.
