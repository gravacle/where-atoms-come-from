# S1 — THE CARRIER K1 — V001 — 2026-08-16

**The project's first construction.** A specific finite oriented regular CW complex of
dimension at most two, with a one-dimensional Hermitian fiber at each vertex and a discrete
unitary connection on edges. Everything below is exhibited, not asserted.

Stage S1 of `FOUNDING_DESIGN_V001.md`. Carrier chosen by the principal, 2026-08-16.

---

## 1. THE COMPLEX

**Vertices (5):** `v0 v1 v2 v3 v4`
**Root:** `r = v0` — the join vertex, the only vertex lying in both triangles.

**Edges (6), each oriented source → target:**

```
e1 : v0 -> v1        e4 : v0 -> v3
e2 : v1 -> v2        e5 : v3 -> v4
e3 : v2 -> v0        e6 : v4 -> v0
```

**Faces (1):** `F`, attached along the closed edge-path `e1 · e2 · e3`.
The triangle `v0 v1 v2` is **filled**. The triangle `v0 v3 v4` is **unfilled**.

Regular: each attaching map is injective on the boundary. Oriented: every cell carries the
orientation displayed. Connected: `v0` joins both triangles.

## 2. TOPOLOGY — COMPUTED, NOT ASSUMED

```
V = 5   E = 6   F = 1
chi = V - E + F = 5 - 6 + 1 = 0
b0 = 1                      (connected)
b2 = 0                      (F has non-empty boundary; no closed 2-cycle)
chi = b0 - b1 + b2   =>   0 = 1 - b1 + 0   =>   b1 = 1
```

**`b1 = 1`.** Exactly one independent 1-cycle survives the filling: the unfilled triangle
`e4 · e5 · e6`. The filled triangle's boundary is a boundary and contributes no class.

This is why the carrier has this shape and not a simpler one: it is the smallest complex
carrying **one face** and **one independent cycle** at the same time.

## 3. FIBERS AND CONNECTION

**Fibers.** `L_v = C` at each of the five vertices, with the standard Hermitian inner product
`<z,w> = conj(z) w`. Rank one, as the predecessor's adopted ruling types it.

**Connection.** A unitary `U_e in U(1)` on each edge, written `U_e = exp(i a_e)` with
`a_e in R / 2 pi Z`. Parallel transport along `e : u -> v` is `z |-> U_e z`. Reverse traversal
transports by `U_e^{-1}`.

Six real parameters `a_1 ... a_6`, each modulo `2 pi`.

## 4. GAUGE, AND THE PARAMETER COUNT THAT FIXES WHAT IS PHYSICAL

**Gauge group.** `g_v in U(1)` at each vertex, `g_v = exp(i theta_v)`. Action on an edge
`e : u -> v`:

```
a_e  ->  a_e + theta_v - theta_u
```

**Count.** 6 edge parameters. 5 vertex parameters. The global phase `theta_v = theta` for all
`v` acts trivially, so the gauge group acts through **4** effective parameters.

```
6 - 4 = 2 gauge-invariant real parameters.
```

**And exactly two invariants exist**, so the count is saturated and nothing is hidden:

```
W_F = exp( i ( a_1 + a_2 + a_3 ) )      the FACE holonomy   — a curvature
W_C = exp( i ( a_4 + a_5 + a_6 ) )      the CYCLE holonomy  — flat, bounded by nothing
```

Each is invariant because every vertex phase enters its exponent once positively and once
negatively around the closed path. `W_F` and `W_C` are the complete gauge-invariant content of
the connection on `K1`.

## 5. CURVATURE, AND THE IDENTITY CHECKED

The discrete curvature two-form on the single face is

```
f(F) = a_1 + a_2 + a_3   (mod 2 pi),        W_F = exp( i f(F) )
```

which is `f = da` evaluated on `F`: the sum of the connection one-form over `partial F`.

**`d^2 = 0`, verified on this complex.** With `partial F = e1 + e2 + e3` as a 1-chain,

```
partial(partial F) = (v1 - v0) + (v2 - v1) + (v0 - v2) = 0
```

so `f` is well defined on `F` and independent of where the boundary traversal starts.

**Bianchi is vacuous here and that is stated rather than passed over:** the identity constrains
sums of face curvatures over closed 2-cycles, and `b2 = 0` means `K1` has none. `W_F` is
therefore unconstrained — any value in `U(1)` is realizable.

**`W_C` is flat and unconstrained.** No face bounds `e4 · e5 · e6`, so no curvature determines
it. It is pure holonomy.

**This is the point of the carrier:** `K1` separates the two, and neither determines the other.

## 6. A WORKED INSTANCE

```
a_1 = a_2 = a_3 = pi/3      =>   f(F) = pi           W_F = exp(i pi)      = -1
a_4 = a_5 = a_6 = pi/2      =>   holonomy = 3 pi/2    W_C = exp(i 3 pi/2)  = -i
```

Non-trivial curvature and non-trivial flat holonomy, independently set. Under any gauge
transformation both values are unchanged; both are computed here by summation over the closed
paths of §4.

## 7. WHAT S1 DELIVERS, AND WHAT IT DOES NOT

**Delivers:** a fully specified carrier; its topology computed; the connection's complete
gauge-invariant content derived and counted; the holonomy/curvature identity checked; a worked
instance.

**Does not deliver, and does not claim to:** any formation condition (S2); any crossing (S3);
any durability. `K1` is a carrier and nothing has been formed on it.

**S1's falsifier was "the project stops if this cannot be written down." It is written down.**

## 8. THE EXTERNAL CONTACT POINT, DISCHARGED

The design names S1's contact as *the holonomy/curvature identity on a discrete connection*.
Discharged at §5: `f = da` on `partial F`, `d^2 = 0` verified on the complex, and the
invariant count of §4 matching the two exhibited invariants exactly.

## 9. CUSTODY

Built under `CUSTODY_V001.md`. Every claim above is exhibited on this page; nothing is cited
from the predecessor, so no transfer grade is consumed. Sealed on creation.
