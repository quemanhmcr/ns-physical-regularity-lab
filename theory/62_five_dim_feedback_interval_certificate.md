# Five-dimensional coupled Hodge fixed-point certificate

## Physical chart

The stationary capacity flow has an orientation-reversing reflection symmetry `x<->y`.  Velocity is polar and vorticity is axial, so the degree-six `T4` feedback sector has a five-dimensional physically fixed subspace.  The coupled maintenance equation preserves this subspace.

The reduced equation is

`G(a)=0`, `a in R^5`,

obtained by embedding `a` into the axial-symmetric `T4` feedback, eliminating degree four with `K44^{-1}`, eliminating all degree-six amplitudes with `K66^{-1}`, and projecting the remaining feedback mismatch back to the five symmetry coordinates.

## Search is not proof

Midpoint Newton is allowed only to locate an approximate root.  Its Jacobian is the analytic Euler/Hodge directional derivative; no finite differences are used.  Midpoint recentering deliberately discards enclosure information, so the Newton orbit is never called a theorem.

## Krawczyk attack

For a candidate `a0`, a box `X=a0+[-r,r]^5`, and a point preconditioner `B` approximating `DG(a0)^{-1}`, evaluate

`K(X)=a0-B G(a0) + [I-B DG(X)](X-a0)`.

All evaluations of `G(a0)` and `DG(X)` in this stage retain their Arb interval enclosures.  If

`K(X) subset interior(X)`

componentwise, standard Krawczyk theory certifies a unique zero of `G` inside `X`.

Module 134 is intentionally an ATTACK rather than an automatic green theorem.  It reports the inclusion result.  If 160-bit inclusion fails because fixed Hodge inversions inflate interval widths while 256/512 succeed, the observer must be repaired; the box is not enlarged merely to force a pass.
