# Canonical obstruction 07: near/far geometric gain--tax

The working question is physical rather than norm-based: how can a small active vortex core acquire the extensional strain needed to keep amplifying?

Two clean channels are isolated.

## Far channel: a closed vortex ring as donor

A circular filament of circulation `Gamma_d` and radius `R` induces on its symmetry axis

`u_z(z) = Gamma_d R^2 / (2 (R^2+z^2)^(3/2))`.

Hence its axial strain magnitude is

`|S_zz| = (Gamma_d/R^2) E(q)`, `q=z/R`,

with

`E(q)=3 q / (2 (1+q^2)^(5/2))`.

The efficiency has a unique maximum at `q=1/2`,

`E_max = 24/(25 sqrt(5))`.

If a target core of scale `ell` and active circulation `Gamma_c` needs an order-one turnover strain `s_* = Gamma_c/ell^2`, then even an optimally placed closed ring of size `R=Lambda ell` needs

`Gamma_d/Gamma_c >= Lambda^2/E_max`.

Thus remote strain is not free: moving the donor geometry `Lambda` core-scales away costs quadratically in donor circulation.  This is a physical collateral law, not yet a global budget theorem.

## Local channel: exact helical/Beltrami calibration

Consider on a periodic domain

`u = -(A/k) (cos(kz), sin(kz), 0)`.

Then

`omega = curl u = A (cos(kz), sin(kz), 0) = -k u`,

so `u x omega = 0`.  The nonlinear Euler part is a pure pressure gradient; Navier--Stokes evolves this field only by viscous decay `A(t)=A_0 exp(-nu k^2 t)`.

Writing `omega=rho xi`, we have `rho=A(t)`, `|grad xi|^2=k^2`, and exactly

`xi.S.xi = 0`,

while

`D_t rho = -nu rho k^2 = -nu rho |grad xi|^2`.

This exact solution is a calibration point: rapid directional variation by itself does not create stretching, but viscosity taxes it exactly.  Productive stretching therefore requires a more special non-Beltrami geometric misalignment.

Together these models sharpen the dilemma:

- export the strain source far away -> donor circulation collateral grows like separation squared;
- keep strong geometry local -> direction variation is viscously taxed, and geometry alone can still be nonproductive through nonlinear cancellation.
