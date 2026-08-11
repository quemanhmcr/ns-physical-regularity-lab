# Remote collars can also have vanishing viscous palinstrophy

Module 155 used a deliberately elementary lower bound `m_l(z) >= const z^-3`.  The actual high-frequency transfer is stronger.

For effective dimension `d=2l+3`, the first-coordinate density of the spherical average is

`p_l(t)=C_l (1-t^2)^l`,  `-1<=t<=1`.

The normalized radial heat eigenfunction is

`Phi_d(z)=integral_-1^1 cos(zt) p_l(t) dt`.

Because `(1-t^2)^l` has a zero of order `l` at `t=+-1`, integration by parts `l` times has no boundary terms. Therefore

`|Phi_d(z)| <= H_l z^-l`,

where

`H_l = ||p_l^(l)||_L1(-1,1)`

and a completely explicit coefficientwise upper bound is available.

Hence for all sufficiently large `z`,

`|Phi_d(z)| <= 1/2`.

On that regime,

`1-Phi_d(z) >= 1/2`

and the screened transfer satisfies

`m_l(z) >= (l+1)/(2 z^2)`.

Let a remote moment collar live on spectral variable `s=z^2 in [M,2M]`.  The exact Hankel enstrophy representation is

`Z = 2/c_d^2 integral rho(s)^2/[m_l(sqrt(s))^2 s^((d-2)/2)] ds`.

The high-frequency bound gives

`Z_collar <= 8/[c_d^2(l+1)^2] M^(3/2-l) ||rho||_2^2`.

For the exact dual-moment collars of module 155,

`||rho||_2^2 = O_n(M^-1)`.

Therefore

`Z_collar = O_n(M^(1/2-l))`.

The radial palinstrophy has one additional factor `s`; on `[M,2M]`,

`P_collar <= 2M Z_collar = O_n(M^(3/2-l))`.

For every `l>=2`, both exponents are negative:

`1/2-l < 0`,

`3/2-l < 0`.

Thus an already-created remote moment collar can simultaneously have

- arbitrarily small vorticity enstrophy;
- arbitrarily small viscous palinstrophy / dissipation throughput;
- exact cancellation of any fixed finite collection of initial Hodge-feedback moments;
- exponentially negligible influence at any fixed later time.

So the missing obstruction cannot be a static enstrophy tax or a static viscous-dissipation tax on the collar itself.

The physical question moves again: **what nonlinear process creates the collar and replaces it as it disappears?**  Any successful regularity mechanism must charge that creation/routing process, not merely the existence of high-frequency vorticity.
