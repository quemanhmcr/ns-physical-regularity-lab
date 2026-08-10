# Residual shape holonomy: kill the gauge, expose the common-strain escape

## THINK — a closed Gram loop is not automatically a physical holonomy

The productive pair state is the directed Gram shape

`G=(alpha,beta,gamma,sign T)`

with

`alpha=xi_a.n`, `beta=xi_b.n`, `gamma=xi_a.xi_b`,

`T^2=det G`,

and

`P=(-alpha beta) det G`.

Because all of these quantities are invariant under a common `SO(3)` rotation of the triad, the map from a lifted triad `(xi_a,n,xi_b)` to Gram shape has a common-rotation gauge freedom.

Therefore a closed path in Gram shape does **not** determine a unique final rotation of the lifted triad.  If one chooses an external connection and calls its accumulated rotation a holonomy, the result depends on that chosen gauge.

This must be attacked before any holonomy is promoted as a nonreusable cost.

## ATTACK I — arbitrary lifted rotation on the same Gram loop

Take any nontrivial closed path of intrinsic coordinates and one concrete lifted triad realizing it.  Multiply the whole lift by any time-dependent common rotation `Q(t)`.

All pairwise dot products and the oriented triple product are unchanged:

`alpha, beta, gamma, T, det G, P` are identical at every instant.

But `Q(T)` can be chosen arbitrarily.

Hence

**closed productive Gram history does not carry a gauge-independent rotational holonomy by itself.**

KILL any cost based on lifted net rotation unless a separate physical connection is supplied by the dynamics.

## THINK II — does Gram motion itself force ancestry replacement?

Once common spin is removed, the remaining Gram current contains symmetric deformation and relative/viscous sources.  It is tempting to identify all residual shape motion with new ancestry geometry.

The exact quadratic Navier-Stokes family attacks that idea:

`u=(a(t)x+eps yz,-a(t)y,0)`,

with

`a=-d log(lambda)/dt`.

For the material pair used in the renewal-memory microscope,

`R=(-Q/lambda,-2lambda,1)`,

`Q=L+eps integral lambda^2 dt`,

while

`omega_a=(0,eps lambda,0)`,

`omega_b=(0,-eps lambda,-eps)`.

The normalized directions remain nontrivial for every `eps>0`, and their Gram geometry can execute a finite excursion driven by `lambda`.

But the non-affine pair-cell renewal is controlled by

`Qdot=eps lambda^2`.

Numerically, the deposited memory must be observed directly as `Delta Q=eps Delta I`, not reconstructed as `Q_after-Q_before` when `eps` is tiny.  The latter needlessly subtracts a small physical tail from an order-one parent bridge coordinate.

Let `eps -> 0+` while keeping a fixed finite-amplitude stretch excursion `lambda:1->2->1`.  Then the Gram excursion remains finite while the deposited bridge memory and relative pair-cell renewal tend to zero.

Therefore

**Gram-shape motion does not universally imply ancestry replacement.**

There is a genuine second source channel: common symmetric strain.

## AUTOPSY

### lifted rotational holonomy as intrinsic reset cost

Killed.  Common `SO(3)` lift is gauge.

### every nontrivial Gram excursion must import new bridge/viscous ancestry

Killed.  Common affine symmetric strain can drive a finite excursion while non-affine renewal tends to zero.

### common symmetric strain as another gauge

Killed.  Unlike common spin, symmetric strain changes dot products and material shape.  It must be physically supported by the velocity field.

## PROMOTE

1. Productive shape space should be treated directly through Gram invariants, not through an arbitrarily chosen lifted connection.
2. The exact Gram current has at least two physically distinct source classes:
   - common symmetric deformation;
   - relative endpoint/bridge/viscous ancestry renewal.
3. A regularity mechanism must therefore account for the physical support of common strain rather than trying to identify all shape motion with ancestry replacement.

## Next frontier — ask who supports the reset strain

The Hodge strain microscope is already the natural source microscope for a symmetric strain at a physical point.  It separates

`S=S_h(r)+S_v(r)`

without a Fourier cutoff.

The next task is to project this exact Hodge split directly onto the productive Gram current itself.
