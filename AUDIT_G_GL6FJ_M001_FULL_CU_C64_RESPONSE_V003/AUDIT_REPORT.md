# Independent hostile audit — GL6FJ V003

## Scope

The subject is the immutable `L=4`, non-self-conjugate `m=(0,0,1)` full
24-real-channel / 300-ray launch target.  The audit treats missing production
output as a required condition, not a defect: it audits the target and its
authorization boundary before the response is collected.

## Independent checks

The audit separately verifies the target's closed-tree manifest and seal,
the exact `m001` measurement plan, source content, and the full upstream
dependency contract.  It replays the target verifier in ordinary and
optimized modes.  It then independently compiles the frozen sampler at
`-O0`, `-O3 -DNDEBUG`, and UBSan, checks its ray and character self-tests,
and samples unit and mixed directions.  At every sampled checkpoint it
recomputes the balance

\[
  \texttt{owner\_total}=\texttt{stay\_contact}
  +\texttt{writer\_contact}+\texttt{stay\_spectral}
  +\texttt{writer\_spectral}+\texttt{stay\_writer\_interference}
  +\texttt{disconnected}.
\]

The source-zero history digest and decision digests agree across all three
compilation modes.  Invalid launch authorization fails without a result,
staging tree, or lock.  The independent replay completed **1,627** checks
without launching production.

## Claim ceiling

This audit certifies only that the frozen `m001` *measurement target* is
admissible for a subsequent nonproduction-free 300-ray launch.  It does not
close the global F3 material-dual owner census.  Consequently it does not
evaluate the Gate-A descent defect, a physical 1PI quotient, a Ward null,
Einstein form, a continuum limit, or gravity.
