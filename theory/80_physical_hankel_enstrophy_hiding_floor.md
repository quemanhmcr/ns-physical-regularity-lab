# Physical Hankel-enstrophy floor for fixed-band hiding

The sharp Legendre cost from module 152 can be expressed directly in the actual vorticity enstrophy of one angular Hodge channel.

Normalize the toroidal sphere field `X_l` by

`integral_S2 |X_l|^2 dOmega = 1`.

Write

`omega(r,n) = a(r) r^l X_l(n)`.

Then

`Z_l = integral_R3 |omega|^2 dx = integral_0^infty |a(r)|^2 r^(2l+2) dr`.

With

`d=2l+3`,  `nu_d=d/2-1`,

this is exactly the radial `L^2(r^(d-1)dr)` norm in the effective dimension selected by the angular channel.

Let

`Phi_d(kr)=0F1(;d/2;-k^2 r^2/4)`

and

`c_d = [2^(d/2-1) Gamma(d/2)]^-1`.

The transform

`A(k)=c_d integral_0^infty a(r) Phi_d(kr) r^(d-1) dr`

is self-inverse and unitary:

`Z_l = integral_0^infty |A(k)|^2 k^(d-1) dk`.

For source radius `L=1`, set `s=k^2`.  The screened transfer from module 148 is

`C_l[Phi_d(kr)] = -m_l(k)`

with

`m_l(k)=(l+1)[1-Phi_d(k)]/k^2`.

Define the feedback spectral density by

`rho(s) = c_d/2 m_l(sqrt(s)) A(sqrt(s)) s^((d-2)/2)`.

Then the heat feedback is

`F(tau)=-integral_0^infty rho(s)e^(-s tau)ds`,

and Plancherel becomes the exact identity

`Z_l = 2/c_d^2 integral_0^infty rho(s)^2 / [m_l(sqrt(s))^2 s^((d-2)/2)] ds`.

Now restrict the cancellation spectrum to the fixed viscous band

`1 <= s <= 2`.

Since `Phi_d(z)` is a spherical average of `cos(z theta_1)`, the elementary physical inequality

`1-cos x <= x^2/2`

gives

`m_l(z) <= (l+1)/(2d)`.

Also `s^((d-2)/2) <= 2^((d-2)/2)`. Therefore

`Z_l >= K_l ||rho||_L2(1,2)^2`,

where

`K_l = 2/[c_d^2 m_max^2 2^((d-2)/2)]`

and equivalently

`K_l = 2^(d/2+2) d^2 Gamma(d/2)^2/(l+1)^2`.

If the feedback derivatives through order `n-1` vanish while the `n`th dimensionless derivative has magnitude `|R|`, module 152 gives

`||rho||_2^2 >= |R|^2 (2n+1) binomial(2n,n)^2`.

Hence the actual channel enstrophy obeys

`Z_l >= K_l |R|^2 (2n+1) binomial(2n,n)^2`.

As `n -> infinity`,

`Z_l >= K_l |R|^2 [2/pi + o(1)] 16^n`.

This is not a Sobolev penalty imposed before looking at the PDE.  The norm here is literally the vorticity enstrophy after the exact Hodge/viscous channel has been diagonalized by its own radial heat operator.

Restoring the source scale: if `r=L s` and the physical vorticity amplitude is `Omega`, then the channel enstrophy is `Omega^2 L^3` times the dimensionless `Z_l`, while viscous time is `tau=nu t/L^2`.

Scope remains strict.  The theorem is a fixed positive spectral-band statement.  The next attacks must let the band move and widen; otherwise one would be freezing exactly the escape that real NS may exploit.
