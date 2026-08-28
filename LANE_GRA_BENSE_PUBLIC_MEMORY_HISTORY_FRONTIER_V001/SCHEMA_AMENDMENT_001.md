# Schema amendment 001 -- deterministic duplicate handling and custody ceiling

**Recorded:** 2026-08-27, after source acquisition and schema inspection but
before the production executable and sealed result.

## Disclosure

The four CSV headers are exactly Time, Extension, and Load, with seconds,
millimetres, and newtons. They contain no deposited hysteron-state word,
trial identifier, covariance, work, heat, or separate query channel. During
schema inspection the official force-jump rule from the publication-linked
scripts was exercised and basic row ranges and branch counts were viewed.
Therefore this is a retrospective DEVELOPMENT diagnostic, not an untouched
validation execution.

The pre-registration left two numerical conventions underspecified. This
amendment fixes them without selecting on an effect direction:

1. At one exact repeated Extension value within a monotone branch, use the
   arithmetic mean of all Load samples at that exact deposited value.
   Opposite-branch exact grids use the intersection, not the union, of
   deposited Extension values. If fewer than three values intersect, collapse
   duplicates by the same mean and interpolate on the frozen 101-point grid.
2. “Overlap-integrated difference” means the trapezoidal integral of absolute
   Load difference divided by overlap width and by the pooled Load range. This
   is dimensionless. The per-file sample-grid spacing is the median positive
   absolute Extension increment.

Zero Extension increments are assigned to the nearest nonzero increment in
sample index; ties go to the earlier increment. Local extrema are the shared
endpoints of successive opposite monotone branches.

## Script-owned transition proxy

The official scripts define candidate force events by

\[
|\Delta \mathrm{Load}|>0.0079\ {\rm N}
\]

and collapse consecutive candidates unless separated by more than 15 samples.
The production executable reproduces that rule exactly and reports every
retained event. It is a force-jump proxy. It is not promoted to a complete
hysteron state word because the deposited scripts depend on unavailable
BasicFunctions/CreateData modules and bind files through directory order rather
than immutable file names. In particular, Fig6.py addresses two experiment
objects while its Figshare item deposits only one CSV.

Accordingly D3 stops at the missing complete state word and common future-query
custody even if Extension and Load happen to match numerically.

