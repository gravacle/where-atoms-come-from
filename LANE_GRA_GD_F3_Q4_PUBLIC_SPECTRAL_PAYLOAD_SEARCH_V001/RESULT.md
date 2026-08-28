# Result

## Physics result

No inspected public release is directly admissible for the GC three-size
spectral screen.  In particular, none supplies the same projected pure-H6
Hamiltonian, fixed sector, complete native FY/GC tensor source, three GC
small-momentum rays, prospectively declared normalization, full `2x2` TT
imaginary-time correlator, raw covariance, and common-clock photon comparator.

The strongest new public root is Zhou et al. (2026):

- data DOI `10.25442/hku.32404548.v1` provides ten machine-readable QFI figure
  tables;
- software DOI `10.25442/hku.32412273.v1` provides pyrochlore SSE/MDL-QMC
  Fortran and a general ED/DSSF archive.

The tables are integrated/derived QFI values rather than raw `G(q,tau)` arrays.
The software is an implementation substrate, not an already matched result.  A
static inspection found that the ED builder natively stores one-, two-, and
three-body terms and its advertised ring-exchange method explicitly throws;
there is no native six-body H6 operator.  The QMC code measures sublattice spin
correlators for the full XXZ model and writes averaged correlators, but does not
own the complete FY/GC tensor insertion or a per-bin covariance packet.

## Strongest executable consequence

After the native `G_L` source ledger is complete, the narrowest no-lab route is
to reuse the public pyrochlore SSE/MDL-QMC infrastructure while making four
physics changes:

1. replace the full-XXZ operator string by the projected pure-kinetic H6 ring
   Hamiltonian on `G_L`, with a fixed zero-flux/complement-symmetric sector;
2. insert the complete six-component native H6 source before TT projection;
3. measure the `2x2` TT and one-link photon correlators at GC19 on
   `G_5,G_10,G_20` in one declared `J_6` clock and count normalization;
4. retain independent-bin outputs so the complete time/channel covariance can
   be reconstructed.

This is source adaptation and a new calculation.  It cannot be represented as
reanalyzing the released QFI tables.

## Claim ceiling

The inspected public corpus contains useful model, phase, static, dynamic-spin,
and implementation evidence, but not the matched response packet required by
GC16--GC19.  The conclusion is non-exhaustive and may be overturned by a new
repository, an author release of raw arrays, or a future simulation.  No
source-visible tensor pole, Ward identity, common cone, gravity response, or
Newton constant follows from this lane.

