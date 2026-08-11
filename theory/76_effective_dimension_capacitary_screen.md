# Effective-dimensional capacitary origin of the Hodge screen

Fix one three-dimensional toroidal angular channel of degree `l`.  Its radial operator is

`D_l = d^2/dr^2 + (2l+2)/r d/dr`.

This is the radial Laplacian in the effective dimension

`d = 2l+3`.

The Hodge feedback functional at source radius `L` is

`C_l^L[a] = -(l+1)/(2l+1) integral_0^L r a(r) [1-(r/L)^(2l+1)] dr`.

Set

`psi_l,L(r) = r^(2-d) - L^(2-d)`.

Then

`r^(d-1) psi_l,L(r) = r [1-(r/L)^(d-2)]`,

and `d-2=2l+1`.  Therefore

`C_l^L[a] = -(l+1)/(d-2) integral_0^L a(r) psi_l,L(r) r^(d-1) dr`.

The screen is thus a capacitary pairing.  The function `psi_l,L` is the radial Newtonian potential of a point source in dimension `d`, grounded at the source sphere:

`D_l psi_l,L = 0` for `0<r<L`,  `psi_l,L(L)=0`.

Green identity in the measure `r^(d-1)dr` gives

`integral psi D_l a r^(d-1)dr = (d-2)[a(L)-a(0)]`

for smooth regular radial profiles.  Multiplying by the Hodge normalization yields

`C_l^L[D_l a] = (l+1)[a(0)-a(L)]`.

This explains several structures at once:

- the original quadrupole screen `1-(r/L)^5` is the `l=2`, `d=7` capacitary kernel;
- the general exponent `2l+1` is simply the Newtonian exponent `d-2` of the effective radial dimension;
- viscosity is `D_l`, while Hodge screening pairs against a `D_l`-harmonic capacitary potential;
- for a spectral mode `D_l a=-k^2 a`, Green identity necessarily produces a `k^-2` transfer factor.

So the Hodge screen and the viscous endpoint law are not two unrelated formulas.  They are the inverse-Laplacian and Laplacian sides of the same angular-channel potential theory.
