# Canonical stress-test 01: strained viscous vortex

Before attacking a general singular cascade, test the proposed toll mechanism on an exact Navier–Stokes balance where stretching and viscosity are both explicit.

Take an axisymmetric extensional strain with rate `a > 0`:

- radial velocity `u_r = -(a/2) r`,
- axial velocity `u_z = a z`.

A steady axial vorticity profile is

`omega(r) = Gamma * a/(4*pi*nu) * exp(-a r^2/(4 nu))`.

It has total circulation `Gamma`.  Its squared-vorticity integral per unit axial length is

`Integral omega^2 dA = Gamma^2 a/(8*pi*nu)`.

Therefore the viscosity-weighted vorticity integral is

`D_omega/L = nu Integral omega^2 dA = Gamma^2 a/(8*pi)`.

Over one strain time `tau_s = 1/a`,

`(D_omega/L) tau_s = Gamma^2/(8*pi)`.

The core radius obeys `delta^2 = 4 nu/a`, so making the coherent core thinner by increasing strain makes instantaneous viscous activity larger in exactly the amount needed to keep the one-strain-time toll independent of core scale.

This is not a regularity proof.  The background strain is non-decaying and the model is a canonical local balance, not a finite-energy global blow-up scenario.  Its role is narrower and valuable: it is an exact NS stress-test showing that a scale-independent toll is physically compatible with the nonlinear stretching/diffusion mechanism.
