# Stationary affine strain has a viscous spectral horizon

The remote collar is cheap after it exists. Can the reusable common strain of the stationary productive amplifier continually create it?

For the validated strain

`S_* = sigma [[1,0,3/(2sqrt2)],[0,1,3/(2sqrt2)],[3/(2sqrt2),3/(2sqrt2),-2]]`,

its eigenvalues are

`c_+ sigma`, `sigma`, `-c_- sigma`,

with

`c_+=(3sqrt2-1)/2`,

`c_-=(1+3sqrt2)/2`.

Consider the exact affine vorticity equation

`partial_t omega + (S_* x).grad omega = S_* omega + nu Delta omega`.

In Fourier variables,

`dot k = -S_* k`,

`dot omega_hat = S_* omega_hat - nu |k|^2 omega_hat`.

Choose `k` on the strongest compressive eigenaxis and vorticity polarization on the strongest stretching eigenaxis.  This is the best affine configuration for simultaneous frequency growth and vorticity stretching. Then

`k(t)=k0 exp(c_- sigma t)`

and

`log |omega_hat(t)/omega_hat(0)|`

`= c_+ sigma t - nu k0^2 [exp(2c_- sigma t)-1]/(2c_- sigma)`.

The instantaneous growth rate is

`c_+ sigma - nu |k|^2`.

Hence the exact affine spectral horizon is

`nu |k|^2 <= c_+ sigma`.

A collar at `|k|^2=M` therefore requires

`sigma >= nu M/c_+`

merely to avoid instantaneous decay. Maintaining it for a finite fraction of its own viscous lifetime has the same scaling `sigma~nu M` with a larger dimensionless constant.

The result is also a pure scaling identity. For any fixed smooth profile `Omega` and

`omega_epsilon(x)=A_epsilon Omega(x/epsilon)`,

the affine Euler operator

`S omega-(Sx).grad omega`

has the same epsilon scaling as `omega`, while `nu Delta omega` carries the extra factor `nu epsilon^-2`. Thus fixed affine strain loses against viscosity as `epsilon->0`; balance requires `sigma epsilon^2/nu=O(1)`.

Now combine this with the already validated stationary-amplifier harmonic occupancy

`E_h(B_R)=(7pi/5)sigma^2 R^5`.

Finite energy gives `R=O(sigma^-2/5)`. At the spectral horizon `sigma~nu M`,

`R_source=O(M^-2/5)`.

Since the collar length is `epsilon=M^-1/2`,

`R_source=O(epsilon^(4/5))`.

The same `4/5` exponent that appeared in the accelerated Hodge reset scale therefore reappears from a completely different physical requirement: replenishing a remote viscous collar with the reusable stationary affine strain.

This does not close regularity.  In fact `R_source/epsilon~epsilon^-1/5` grows, so finite energy still permits a strain source region larger than the collar core.  What is ruled out is a **fixed large-scale affine reservoir** maintaining arbitrarily high-frequency bookkeeping.  The source itself must intensify and move inward, or the maintenance must be supplied by genuinely non-affine same-scale dynamics.
