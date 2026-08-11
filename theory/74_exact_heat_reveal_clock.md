# Exact heat reveal clock for maximally hidden radial profiles

Finite derivative hiding can be misleading. Work in the source-scaled coordinate `s=r/L`, and write a fixed angular channel as

`omega = Omega b(s) X_l`,

where `X_l = n cross grad_S Y_l`. The exact viscous time is

`tau = nu t / L^2`.

For the degree-`2n` maximally delayed profile,

`C_l[D_l^j b]=0` for `j<n`,

while `D_l^(n+1)b=0`. Therefore pure heat gives the exact identity

`F_l(tau) = tau^n/n! C_l[D_l^n b_0]`.

There are no neglected higher Taylor terms.

Let

`Z = integral_0^1 s^(2l+2)|b|^2 ds`.

For the corresponding actual vorticity enstrophy `Zeta_0` and generated lower harmonic velocity `u_h`, the angular harmonic identities give

`2 E_h /(L^2 Zeta_0) = |F_l|^2 / ((l+1) Z)`.

Thus the basis-independent physical Hodge-energy reveal clock is

`Theta_l,n = [n! sqrt((l+1) Z) / |C_l[D_l^n b]|]^(1/n)`.

At

`t = Theta_l,n L^2/nu`,

the generated lower harmonic kinetic occupancy satisfies

`2 E_h = L^2 Zeta_0`.

This is a physical source-scale comparison; rescaling the spherical harmonic basis cancels between input enstrophy and output harmonic energy. A second radial-coordinate clock omitting the factor `sqrt(l+1)` is reported only as an observer diagnostic.

If `Theta_l,n` grows with hiding order, increasingly elaborate radial structure genuinely buys viscous invisibility time. If it stays bounded or decreases, arbitrarily many vanishing initial derivatives are only a coordinate illusion: the first nonzero derivative grows fast enough that physical feedback reappears on an `O(L^2/nu)` or shorter clock.

Module 147 measures this exact clock to hiding order 16 at several angular channels. No conclusion is based on polynomial degree itself.
