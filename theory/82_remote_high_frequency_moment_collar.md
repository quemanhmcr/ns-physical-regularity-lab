# Remote high-frequency moment collar

The fixed-band costs of modules 152-154 are sharp, but real Navier–Stokes is not obliged to keep the cancellation spectrum inside a fixed band.

This freedom produces a direct analogue of the earlier spatial closure-collar escape.

## Low carrier

Take a screened feedback density

`rho_0(s)=1`,  `1<=s<=2`.

Its moments are

`mu_q = integral_1^2 s^q ds`.

At any fixed positive time, for example `tau=1`, it carries a definite positive heat feedback.

## Remote moment collar

Fix any finite `n`.  Let `phi_j(x)`, `j=0,...,n-1`, be the exact dual moment polynomials on `[1,2]`:

`integral_1^2 x^q phi_j(x) dx = delta_qj`.

For a large spectral scale `M`, define on `M<=s<=2M`

`eta_M(Mx) = - sum_{j=0}^{n-1} mu_j M^(-j-1) phi_j(x)`.

Then, exactly,

`integral_M^(2M) s^q eta_M(s) ds = -mu_q`,  `q=0,...,n-1`.

Thus `rho_0+eta_M` hides its first `n` heat-feedback derivatives at `tau=0`.

The correction contributes essentially nothing at a positive time because

`|F_collar(tau)| <= exp(-M tau) ||eta_M||_1`.

## Physical enstrophy of the remote collar

For the angular Hodge channel `l`, effective dimension is `d=2l+3` and

`m_l(z)=(l+1)[1-Phi_d(z)]/z^2`.

A quantitative lower bound follows directly from the spherical-average representation.  For `z>=2`, restrict the first-coordinate integral to

`1/(2z) <= |theta_1| <= 1/z`.

On this set,

`1-cos(z theta_1) >= 1-cos(1/2)`

and the spherical density is bounded below. Therefore

`1-Phi_d(z) >= c_l^*/z`

with an explicit positive `c_l^*`, hence

`m_l(z) >= b_l z^-3`.

The exact Hankel/Plancherel identity from module 153 then gives on `s in [M,2M]`

`Z_collar <= C_l M^(5/2-l) ||eta_M||_2^2`.

But the dual-moment construction has

`||eta_M||_2^2 = O_n(M^-1)`.

Consequently

`Z_collar <= C_l,n M^(3/2-l)`.

For every `l>=2`,

`3/2-l < 0`.

Therefore, for every fixed finite number of hidden moments, the cancellation bookkeeping can be sent to sufficiently high viscous frequency with **arbitrarily small additional physical enstrophy**.

After normalizing the total future feedback to one, the total enstrophy approaches that of the unchanged low-band carrier.

## KILL

This kills the global inference

`many hidden initial heat moments => large total enstrophy`

when spectral support is unrestricted.

The sharp `16^n` fixed-band floor remains correct and useful, but its fixed-band hypothesis is essential.

The escape is not dynamically free.  A collar at `s~M` dies on the viscous clock

`Delta tau ~ M^-1`.

Thus the new question is no longer whether a static state can hide the feedback cheaply.  It can.  The question is whether nonlinear Navier–Stokes can **continually regenerate** the remote high-frequency cancellation collar while keeping the productive low-frequency ancestry alive.

That is a maintenance / ancestry problem, not a static norm problem.
