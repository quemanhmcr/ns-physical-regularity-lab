# Hodge-lifted degree-four null servo

## AUTOPSY — vorticity control without velocity lifting is incomplete

The six-dimensional affine solve `P_null L_S V=-N_4` is a valid vorticity-level statement but not yet a physical local velocity mechanism.  Vorticity does not uniquely determine velocity without the Hodge boundary condition.

Inside the Hodge source ball the vortical velocity must satisfy

`curl U=V`, `div U=0`, `U.n=0` on the source boundary.

The generated degree-four null sector has two pieces: poloidal `l=3` and toroidal `l=4`.

## Exact Hodge lift

For a harmonic homogeneous `H_3`,

`curl[r^2 x cross grad H_3]`

is the canonical degree-four poloidal `l=3` vorticity.  Its velocity is already tangent.

For a harmonic homogeneous `H_4`, the toroidal vorticity

`V_T=x cross grad H_4`

has the unique tangent Hodge lift

`U_T=(7/22) r^2 grad H_4 -(4/11)H_4 x -(5/22)grad H_4`.

The last term is curl-free, divergence-free and harmonic.  It is not optional: it enforces `U_T.n=0` at the physical source boundary.

For `l=2` the same general formula reproduces the original tangent Hodge strain carrier exactly.  Thus this is an intrinsic continuation of the established microscope, not a new gauge convention.

## Correct degree-four linear maintenance operator

Because the `l=4` Hodge lift contains a degree-three harmonic velocity `U_3(V)`, the full degree-four linear response of a servo vorticity is

`[V,Sx]+[omega_2,U_3(V)]`.

Therefore define

`K_H V=P_null([V,Sx]+[omega_2,U_3(V)])`.

The earlier affine operator is demoted to an intermediate coordinate model.

Generate the Krylov space of `K_H` from the actual nonlinear null field `N_4`, let the physics determine its dimension, and solve

`K_H V_servo=-N_4`.

If this solve succeeds, the first null generation is maintainable even after respecting the exact Hodge velocity relation.  The next unavoidable responses occur at higher homogeneity through the degree-five vortical part of the servo velocity and nonlinear servo self-interaction.
