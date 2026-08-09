# Hodge transition radius: a PDE-defined shrinking source scale

The harmonic strain floor and the Hodge strain microscope combine into an exact local-to-global statement for any smooth finite-energy incompressible velocity field in `R^3`.

For every radius `r`, let `S_h(r)` be the center strain of the harmonic Hodge component in `B_r(x0)`. The energy floor gives

`E0 >= E(B_r) >= E_h(B_r) >= (2pi/15) |S_h(r)|^2 r^5`.

Hence

`|S_h(r)| <= sqrt(15 E0/(2pi)) r^(-5/2)`.

Smoothness gives

`S_h(r) -> S(x0)` as `r -> 0`,

because the degree-two boundary-flux projector reads the linear strain term of the local Taylor expansion while translations, rigid rotations, and higher orders vanish or decouple.

Let `s=|S(x0)|>0` and choose `0<theta<1`. Define

`R_theta = [15 E0/(2pi theta^2 s^2)]^(1/5)`.

At sufficiently small radius, `|S_h(r)|` is close to `s`; at `R_theta`, the energy floor forces `|S_h(R_theta)| <= theta s`. By continuity there is a transition radius `r_theta <= R_theta` where the harmonic fraction reaches `theta` (or crosses it earlier). At that scale

`|S_v(r_theta)| = |S-S_h(r_theta)| >= (1-theta)s`.

Thus a large strain necessarily acquires a genuinely vortical Hodge support scale no larger than

`r_theta <= C_theta E0^(1/5) s^(-2/5)`.

For `theta=1/2`,

`r_1/2 <= [30 E0/(pi s^2)]^(1/5)`.

This is more precise than the earlier donor-horizon heuristic. It does not assume compact vorticity support, a vortex-ring donor, a multipole truncation, or an arbitrary near/far cutoff. It uses the actual velocity flux through physical spheres and the exact Hodge decomposition.

It is not yet a regularity theorem. It says that if strain becomes arbitrarily large while total kinetic energy remains finite, the radius at which at least half of that strain must be carried by the vortical Hodge component shrinks to zero. The next regularity problem is therefore forced into the vortical branch at a PDE-defined shrinking scale.
