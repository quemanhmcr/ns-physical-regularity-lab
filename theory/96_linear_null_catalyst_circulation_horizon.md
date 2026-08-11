# Finite-energy circulation horizon of the linear null catalyst

Use the exact converter

`omega=A Bx`, `B=diag(2,-1,-1)`.

This choice is especially clean because

`(B^2)_TF=B`.

Module 168 then gives

`dS_v/dt=(A^2 r^2/7)B`.

To supply a core maintenance source-rate coefficient

`j_core=k nu/epsilon^2`,

choose

`A^2=7 k nu/epsilon^4`.

## Exact energy of the null catalyst

Its tangent Hodge velocity is

`u=-(A/3)x cross Bx`.

Direct ball integration gives

`E(B_L)=(4pi/105)A^2 L^7`.

At finite-energy horizon `E=E0`,

`L=[15 E0 epsilon^4/(4pi k nu)]^(1/7)`.

Therefore

`L~epsilon^(4/7)`,

`L/epsilon~epsilon^-3/7`.

The null catalyst can indeed be a broad reservoir compared with the maintained core.

## Actual circulation relocation

On a sphere of radius `R`,

`omega.n=A R(3 n_x^2-1)`.

The positive cap `n_x>1/sqrt(3)` has exact flux

`Gamma_cap(R)=4pi A R^3/(3sqrt(3))`.

By Stokes this is an actual circulation around the cap boundary.

At the core radius,

`Gamma_core/nu~epsilon ->0`.

At the energy-horizon halo radius,

`Gamma_halo/nu~epsilon^-2/7 ->infinity`.

So the transaction-null nonlinear catalyst does something very specific: it keeps the circulation attached to the tiny productive core arbitrarily small, but only by storing increasingly large circulation ancestry in a broader shrinking halo.

This is not a contradiction with finite kinetic energy. The velocity amplitude and occupied volume scale so that energy stays finite while vorticity flux grows.

## Next ancestry gate

The cap flux is not automatically a single closed frozen lineage. Earlier specific-volume theory gives a strong bound only for complete closed lineages contained in the source.

The next attack therefore compares `Gamma_halo` with that closed-lineage capacity. If the required halo flux exceeds what smooth frozen closed ancestry can fit in `B_L`, the only survivors are precisely the already isolated branches:

- through-going winding/closure outside the source;
- material recruitment across a moving source boundary;
- viscous ancestry renewal.
