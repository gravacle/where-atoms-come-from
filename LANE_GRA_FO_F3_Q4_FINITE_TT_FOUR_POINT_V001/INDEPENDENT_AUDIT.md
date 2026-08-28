# Independent hostile audit: finite inherited F3/q4 TT composite screen

**Lane ID:** `GRA-FO-F3-Q4-FTT4-V001`

**Audit date:** 2026-08-27

**Disposition:** `PASS_AFTER_COMPOSITE_TYPING_AND_SEAL_REPAIR`

## 1. Scope and independence

I did not accept the packet's stated sector, matrix spectrum, response
coefficients, or interpretation as premises.  I rebuilt the periodic graph
from the shifts `(0,1,5,19) mod 30`, enumerated its simple cycles, generated
the ice component from the declared seed, constructed the inherited
Hamiltonian directly from legal alternating-hexagon flips, and recomputed the
source response.  I then used independent reductions for the load-bearing
numbers: translation-orbit Fourier blocks for the spectrum, nonsingular
augmented linear solves for reduced resolvents, and a four-step Lanczos
reconstruction for the selected spectral measure.

The finite arithmetic is reproducible.  The original scientific typing was
not fully sound: it treated a four-`Q` composite cumulant too much like the
connected four-one-link kernel requested by `FM`.  That defect is material
because a quadratic composite of a Gaussian fundamental field can have
nonzero higher composite cumulants.  The theorem and self-audit have been
narrowed accordingly; no inherited matrix element or numerical checksum was
changed.

## 2. Quotient, cycle set, sector, and translation closure

The twenty unordered triple sums of the four shifts are distinct modulo 30.
Direct enumeration gives 60 vertices, 120 physical links, and exactly 120
simple six-cycles.  Every six-cycle has three q4 labels, each appearing twice,
so the finite Hamiltonian uses the plaquette-complete elementary hexagon set
rather than a hand-selected subset.

The stated seed loop contains eight distinct links.  Each successive pair
shares exactly one vertex, its occupation alternates on the frozen reference,
and its lifted winding has magnitude 60.  Closure under every legal hexagon
flip yields exactly 180 ice states and 420 undirected transitions, with graph
degree from four through six.  Every state remains two-in/two-out.  Translation
by one cyclic cell maps every state back into the same component.  More
strongly, the component decomposes into **six free translation orbits**, each
of length 30.  Global occupation complement produces a disjoint isomorphic
180-state component; the audited packet correctly does not call either one
the full periodic ice Hilbert space.

## 3. Hamiltonian and algebraic spectrum

With `J6=1`, the directly constructed matrix has one `-1` for each legal
inherited ring transition and no diagonal or fitted term.  It is real
symmetric, contains exactly 420 undirected graph edges, and commutes with the
cyclic translation permutation.

Reducing the six translation orbits at each of the 30 momenta gives thirty
Hermitian `6 x 6` blocks.  Their combined eigenvalues reproduce the full
`180 x 180` spectrum.  The `m=0` block has

`{-2-2 sqrt(2), 0, 0, -2+2 sqrt(2), 2, 2}`,

while `m=5` (and its conjugate `m=25`) has

`{-1-2 sqrt(3), -1, 1, 1, 1, -1+2 sqrt(3)}`.

The remaining-block floor is `-4.410987667370205`.  Thus the audited sector
has the unique ground energy `-2(1+sqrt(2))` and first gap
`1+2sqrt(2)-2sqrt(3) = 0.364325509608438`.  These statements are exact for the
declared connected component only; they do not identify the ground sector of
the complete periodic ice Hilbert space.

## 4. Source typing and normalization

The local tensor `Q=STF(F F^T)` is symmetric, traceless, and invariant under
global occupation complement.  Expanding it on two-in/two-out states confirms
that its one-link-square contribution is isotropic and disappears under STF,
leaving an even pair observable.  The selected reciprocal representative is
`q=(1,5,-11)/30`; it is the unique shortest alias in the checked reciprocal
neighborhood.  The fractional `B`-basis offset maps to the declared
label-zero bond vector.  Both polarization tensors are unit Frobenius norm,
mutually orthogonal, transverse, and traceless, and the complex source
transforms with cyclic momentum `+1`.

This earns a kinematic finite-volume TT projection, not a helicity-two
particle assignment.  The factors `1/sqrt(60)` and `sqrt(2)` are consistent
finite-quotient conventions.  They are not a continuum residue, a physical
cell-volume calibration, or the one-link normalization required for external
photon-leg amputation.

## 5. Independent response replay

The static two-composite susceptibility matrix and its doublet eigenvalues
reproduce the theorem.  For the selected plus-cosine source, an augmented
linear solve using `H-E0+|0><0|` instead of the spectral inverse gives

- `A = 0.565423567997862`,
- `C = 0.161972957692486`,
- `B = 0.085882282424106`,
- `E4 = -B+AC = 0.005701045233546`,
- `W2 = 1.130847135995723`,
- `W4 = -0.136825085605100`, and
- `Gamma_comp^(4) = -W4/W2^4 = 0.083666214307836`.

The low- and high-susceptibility cosine eigenchannels independently give
positive scalar composite Legendre coefficients `18.586988116578` and
`0.090544748133189`, respectively.  These are zero-frequency,
time-integrated derivatives of the ground energy.  They are not equal-time
moments.

The load-bearing correction is this: `W2` is a connected two-`Q`
susceptibility, and `W4` is a connected **four-Q** cumulant.  Since each `Q`
is bilinear, the latter is an eight-one-link object.  A scalar Legendre
transform in the composite coordinate amputates composite susceptibilities;
it does not perform four independent one-link external-leg amputations or a
two-particle-reducible subtraction.  Consequently these numbers prove
non-Gaussian statistics of the selected composite coordinate, but not a
nonzero photon four-one-link 1PI/channel-2PI interaction, attraction, or
binding.

## 6. Spectral and threshold audit

Starting only from `O|0>`, the Krylov space closes after four Lanczos steps.
Its gaps and weights are

| gap / `J6` | weight |
|---:|---:|
| 3.194109035554332 | 0.005026104004432 |
| 3.490165912028476 | 1.965864248197576 |
| 6.166688337463908 | 0.003649484732984 |
| 9.139267639373482 | 0.000000908886868 |

They reproduce the selected `W2`.  Independently scanning transverse
one-link sources gives the declared finite proxy minimum at cyclic momenta
`(9,22)`, with sum `2.059674505691458 J6`.  The lowest selected composite pole
is above that proxy and the same energy also has one-link response.

This is a useful negative finite-sector diagnostic, but **not a no-bound-state theorem**.
The proxy is not the thermodynamic two-particle
continuum; finite-volume shifts, other sectors, and nonperturbative channel
resummation remain unresolved.  Likewise, equal energies in the doubled
complement blocks do not prove equal spin or operator identity: the even and
odd parity states remain distinct.

## 7. Dependency, seal, and claim-boundary audit

The final theorem and independent-audit bytes of `CW`, `FK`, `FL`, and `FM`
were recomputed and agree with all eight hashes pinned in this lane.  Each
dependency was also subjected in memory to an appended-byte tamper, and every
tampered digest failed the expected hash.  The local manifest seals the final
theorem, self-audit, independent audit, executable verifier, and verification
transcript.  A full manifest replay is part of the executable gate.

The corrected packet no longer promotes a finite composite insertion to a
photon-leg 1PI vertex, a channel-2PI kernel, a Bethe--Salpeter result, a
thermodynamic pole, helicity two, a Ward identity, RGRL-B, or gravity.  It also
inherits `FM`'s corrected analytic boundary: finite-order locality does not
imply that every dressed or nonperturbatively resummed propagator is incapable
of shifting or acquiring a pole.

## 8. Sharp next lawful datum

The next calculation is the normalized, connected, frequency- and relative-
momentum-resolved four-one-link response with four independent transverse
sources.  The tuple must freeze all four external `(omega,p,polarization)`
labels, the one-link two-point function and residue convention, cell/field
normalization, and the declared two-particle-reducible subtraction.  Only
after external-leg amputation and channel-2PI subtraction is there a finite
kernel suitable for Bethe--Salpeter analysis.  Volume scaling, a surviving
pole and residue, a common cone, helicity-two transformation, and a rank-two
Ward test would still be required for the gravitational carrier claim.

## Final disposition

`PASS_AFTER_COMPOSITE_TYPING_AND_SEAL_REPAIR`.  The quotient, 180-state
translation-closed sector, 420 transitions, exact block spectrum, static
two-`Q`/four-`Q` responses, composite Legendre coefficients, four selected
poles, and finite threshold proxy all reproduce.  The material claim-typing
defect has been repaired.  The packet now earns an exact finite composite
precursor and a negative below-proxy screen, while leaving the actual
connected four-one-link/channel-2PI and thermodynamic gravity targets open.
