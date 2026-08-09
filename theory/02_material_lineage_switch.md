# Canonical stress-test 02: material-lineage switching by viscous circulation leakage

The first canonical test showed that a coherent strained vortex that stays with its material lineage pays a `Gamma^2`-type toll per strain time.  A putative singular cascade can try to escape by abandoning the old material core and transferring its active circulation into a new smaller core.  This module asks whether that switch can be physically free.

Use the exact Lamb-Oseen vortex

`omega(r,t) = Gamma/(4*pi*nu*t) exp(-r^2/(4 nu t))`.

Its velocity is purely azimuthal.  Therefore every circle `r=R` is a material loop: fluid particles rotate around it but have no radial velocity.  The fraction of total circulation inside that material circle is

`f(q) = Gamma_R/Gamma = 1-exp(-q)`,

with

`q = R^2/(4 nu t)`.

As time increases, viscosity diffuses vorticity outward and `Gamma_R` falls.  For a material loop, the nonlinear advection and pressure do not supply the radial circulation leakage; the loss is exactly the viscous vorticity flux across the loop.

The exact squared-vorticity integral is

`Integral omega^2 dA = Gamma^2/(8*pi*nu*t)`.

Hence the viscosity-weighted dissipation per unit axial length is

`D(t) = nu Integral omega^2 dA = Gamma^2/(8*pi*t)`.

Suppose a fixed material core changes from containing fraction `f1` of the circulation to fraction `f2`, with `0 < f2 < f1 < 1`.  Define

`q_i = -log(1-f_i)`.

Since `t = R^2/(4 nu q)`, the episode satisfies

`t2/t1 = q1/q2`.

The integrated viscous cost is therefore

`C_switch = Integral_{t1}^{t2} D(t) dt = Gamma^2/(8*pi) log(q1/q2)`.

Crucially, `R` and `nu` cancel.  Thus in this exact NS diffusion model, transferring any prescribed nonzero fraction of circulation out of a material core has a positive scale-independent cost proportional to `Gamma^2`.

This is not yet a general 3D theorem.  Lamb-Oseen is a canonical diffusion model with non-finite total kinetic energy in the infinite plane, and it has no vortex stretching.  Its role is narrower: it stress-tests the *escape branch* in which singular activity changes material lineage.  It shows that the same candidate currency `Gamma^2` appears in an exact NS model of viscous lineage switching that appeared in the exact strained-vortex model of material-lineage persistence.

The emerging canonical dichotomy is therefore:

- stay with the material lineage -> strained-vortex `Gamma^2` toll;
- switch material lineage through viscous leakage -> Oseen `Gamma^2` toll.

The next research question is whether this common currency survives when stretching and lineage switching occur simultaneously in genuinely three-dimensional, non-axisymmetric local geometry.
