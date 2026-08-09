# Research turn 04: strain support, then vortex-line closure

## Physics question
Can a singular core obtain scale-critical strain from an arbitrarily short donor, thereby evading the `Gamma^2 per-unit-length` toll by shrinking coherent length together with core thickness?

## Canonical open-segment test
A straight finite vortex filament of circulation `Gamma`, length `L`, at perpendicular distance `d` induces

`u = Gamma L / (4 pi d sqrt(d^2+(L/2)^2))`.

For `alpha=L/(2d)`, the dimensionless shear-strain efficiency is

`E(alpha)=4 pi d^2 |S_xy|/Gamma = alpha(alpha^2+2)/(1+alpha^2)^(3/2)`.

Hence `E(alpha) ~ 2 alpha` for very short donors, but `E(1)=3/(2 sqrt(2))`, already order one. The maximum occurs at `alpha=sqrt(2)`, and `E -> 1` as `alpha -> infinity`.

So the naive conjecture `strong scale-critical strain requires L/d -> infinity` is false.

## The physical structure we had omitted
An isolated open vortex segment is not a complete admissible vorticity geometry: `div omega=0`, so vortex lines do not terminate in the fluid. The segment must connect onward or close.

The correct constraint is therefore not support length alone, but vortex-line closure / connected global geometry.

## Canonical closed-loop completion
For a circular vortex filament with isotropic core regularization, rotational and reflection symmetry force self-induced velocity at every material point to be the same axial vector. Therefore the closed loop translates but has zero tangential material stretching rate.

Thus the compact closed completion of the efficient open donor does not inherit its stretching action on itself. Sustained stretching requires symmetry breaking and/or interaction with other/nonlocal vorticity.

## Decision
Killed: `strong strain implies a scale-independent lower bound on donor length`.

Promoted working principle:
> Strong stretching is controlled by the non-cancelling anisotropic part of a divergence-free closed/connected vorticity geometry, not by donor length alone.

Next: quantify the anisotropy/nonlocal-interaction price required for a closed vortex geometry to sustain stretching of its high-vorticity core.
