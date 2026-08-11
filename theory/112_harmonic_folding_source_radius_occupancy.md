# Harmonic folding source-radius occupancy gate

Module 184 shows that a degenerate zero manifold can be folded by an **irrotational** non-affine velocity.  Therefore the folding price cannot be assigned to local vorticity ancestry or null enstrophy.  The correct first question is how far the harmonic actor can reach from the vorticity that sources it.

Use the exact homogeneous harmonic potential

`phi_m = -Re[(x+i z)^(m+1)]/(m+1)`.

Its velocity has degree `m`.  On the `z` axis,

`|partial_z^2 u_x| = m(m-1)|z|^(m-2)`.

Let the productive/zero-manifold core be sampled at `z=epsilon`, and let `B_L` be a ball on which the folding actor is harmonic.  Thus `L` is a physical harmonicity radius: vorticity sourcing this actor lies outside that ball.

## Exact kinetic occupancy

Set `n=m+1` and `P_n=Re[(x+i z)^n]`.  Green's identity gives

`int_{B_L} |grad phi_m|^2 dx = c_m L^(2m+3)`,

where

`c_m = (1/(m+1)) int_{S^2} P_(m+1)^2 dS`

and explicitly

`int_{S^2} P_n^2 dS = pi 2^(2n+1) (n!)^2/(2n+1)!`.

Now multiply the mode by amplitude `a` and prescribe a physical curvature source `Q` at the core:

`Q = |a| m(m-1) epsilon^(m-2)`.

Writing

`lambda=L/epsilon`,

one obtains the exact scale-free occupancy law

`E_harm/(Q^2 epsilon^7)`

`= c_m lambda^(2m+3)/[m^2(m-1)^2]`.

The entire source-clearance penalty is therefore

`E_harm(lambda)/E_harm(1)=lambda^(2m+3)`.

## Physical consequence

For any fixed `lambda>1`, high angular degree is exponentially expensive in kinetic occupancy.  The coefficient `c_m/[m^2(m-1)^2]` is only polynomial in `m`, so bounded normalized occupancy forces

`log lambda = O(log m/m)`

as `m->infinity`.

Thus increasingly fine harmonic folding cannot be supplied from a fixed relative distance.  It is routed into a **near-contact source branch**: the vorticity generating the harmonic field must approach the shrinking core in relative distance as the required degree rises.

## Scope

This is an exact theorem for one homogeneous harmonic folding mode normalized by its physical curvature source.  It does **not** yet prove that fixed-time cap multiplicity `N` requires a single mode with `m comparable to N`; mixtures can interfere, and a clever low-degree time-dependent actor may fold by repeated action.

The next attack should therefore target the surviving near-contact escape directly: if the harmonicity gap `delta=L-epsilon` collapses, what circulation / Kelvin / turnover ancestry must live in the vorticity collar that sources the folding field?
