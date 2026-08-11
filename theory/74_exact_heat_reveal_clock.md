# Exact heat reveal clock for maximally hidden radial profiles

Finite derivative hiding can be misleading. Work in the source-scaled coordinate `s=r/L`, and write a fixed angular channel as

`omega = Omega b(s) X_l`.

The exact viscous time is

`tau = nu t/L^2`.

For the degree-`2n` maximally delayed profile,

`C_l[D_l^j b]=0` for `j<n`,

while `D_l^(n+1)b=0`. Therefore pure heat gives the exact identity

`F_l(tau) = tau^n/n! C_l[D_l^n b_0]`.

There are no neglected higher Taylor terms.

Let

`Z = integral_0^1 s^(2l+2)|b|^2 ds`.

A natural radial-coordinate comparison is

`A_rad = |C_l[D_l^n b]|/sqrt(Z)`,

with exact coordinate reveal clock

`Theta_rad = [n! sqrt(Z)/|C_l[D_l^n b]|]^(1/n)`.

Earlier versions also reported a rescaling by `sqrt(l+1)`.  Those numbers are retained in the experiment for historical cross-run comparison only.

## AUTOPSY — the companion is not an independent Hodge energy component

The field `C_l[a] grad H_l` is the harmonic **companion inside the tangent div-curl representation** of the vorticity channel.  It is not the Hodge-orthogonal harmonic field `h` in the physical split `u=h+v`.

The remaining radial vortical part of the tangent velocity can cancel the companion strongly away from the core.  Consequently `|C_l|^2` by itself is not a lower bound for total kinetic energy, and the previous interpretation of the rescaled clock as an independent "Hodge kinetic occupancy reveal time" is demoted.

This resolves an apparent tension with curl-Poincare: highly concentrated radial vorticity can have small total enstrophy and small tangent velocity energy while producing a large local companion/Taylor coefficient near the center; the cancellation occurs in the full tangent field, not in the coefficient `C_l` itself.

What survives unchanged is the exact heat statement:

- arbitrary finite-order initial hiding is possible;
- for the maximally delayed polynomial the first reveal is an exact monomial in `tau`;
- the first nonzero coefficient can grow very rapidly with hiding order.

Modules 148 onward therefore treat `C_l` as a **screened local feedback/jet coordinate**, not as a separately occupied harmonic-energy reservoir.
