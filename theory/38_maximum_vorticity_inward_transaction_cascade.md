# Maximum-vorticity inward productive-transaction cascade

## THINK — if blow-up exists, which physical source must keep amplifying the maximum?

Let

`M(t)=||omega(t)||_infinity`.

At a spatial maximum point `x_*(t)` with `omega=M xi_*`, the vorticity-magnitude equation is

`(partial_t+u.grad)|omega|=(xi.S.xi)|omega|+nu[Delta|omega|-|omega||grad xi|^2]`.

At the spatial maximum, `grad|omega|=0` and `Delta|omega|<=0`.  Hence viscosity cannot be a positive source of the maximum:

`D^+ log M(t) <= [xi_*.S(x_*,t).xi_*]_+`

in the usual upper-Dini/maximizer sense.  Therefore finite-time blow-up `M(t)->infinity` forces infinite positive strain action at maximum-vorticity points.

The question is then purely physical: from what spatial scales can that positive strain keep arriving?

## Hodge split at the moving maximum

Fix any physical outer radius `R>0`.  At each maximizing point use the exact Hodge split in `B_R(x_*)`:

`S(x_*)=S_h(R)+S_v(R)`.

For the instantaneous maximum direction define the STF production tensor

`A_*=xi_* tensor xi_* - I/3`,

so `|A_*|^2=2/3` and, because all strain and transaction tensors are trace free,

`xi_*.S.xi_*=A_*:S`.

For any `0<delta<R`, split the vortical Hodge transaction into

`J_inner(delta,R)`

`=integral_0^delta [1-(rho/R)^5] (A_*:Q(rho)) d rho/rho`,

and

`J_outer(delta,R)`

`=integral_delta^R [1-(rho/R)^5] (A_*:Q(rho)) d rho/rho`.

Thus

`xi_*.S.xi_*=A_*:S_h(R)+J_outer+J_inner`.

No Fourier shell or dyadic scale is introduced.  `delta` and `R` are ordinary physical radii in the exact Hodge representation.

## Fixed harmonic source has finite action

The exact harmonic Hodge energy floor gives

`E0 >= (2pi/15)|S_h(R)|^2 R^5`.

Hence at every time

`|A_*:S_h(R)| <= sqrt(2/3) sqrt(15 E0/(2pi)) R^(-5/2)`.

On every finite time interval, the harmonic contribution at fixed `R` therefore has finite `L1_t` action.

## Sharp projector makes every fixed outer annulus finite-action

The sharp shell transaction theorem gives

`integral_{S^2}|n cross omega|^2 dOmega >= (20pi/9)|Q(rho)|^2`.

Cauchy-Schwarz in the native radial measure yields

`|J_outer|^2`

`<= [9 |A_*|^2/(20pi)] I(delta,R) Z_Q(delta,R)`,

where

`Z_Q(delta,R)=(20pi/9) integral_delta^R rho^2 |Q(rho)|^2 d rho`

is the minimum productive enstrophy in that annulus, and

`I(delta,R)`

`=integral_delta^R [1-(rho/R)^5]^2 rho^(-4) d rho`

`=1/(3 delta^3)-25/(21 R^3)+delta^2/R^5-delta^7/(7 R^10)`.

Because `|A_*|^2=2/3`,

`|J_outer|^2 <= [3 I(delta,R)/(10pi)] Z_Q(delta,R)`.

The productive projector is an orthogonal part of actual vorticity enstrophy, so even though the center `x_*(t)` moves,

`Z_Q(delta,R,t) <= integral_{B_R(x_*(t))}|omega|^2 dx <= integral_R3 |omega|^2 dx`.

The Navier-Stokes energy identity therefore gives, for every `T0<T`,

`integral_0^T0 |J_outer|^2 dt`

`<= [3 I(delta,R)/(10pi)] E0/nu`.

Thus every fixed outer annulus has finite `L2_t`, hence finite `L1_t`, productive strain action on finite time intervals.

## Consequence — blow-up forces inward transaction action through every physical radius

Suppose `M(t)->infinity` as `t->T<infinity`.  The positive maximum-vorticity strain action must diverge.  But for every fixed `R` the harmonic term has finite action, and for every fixed `0<delta<R` the outer transaction term has finite action.  Since

`[a+b+c]_+ <= |a|+|b|+[c]_+`,

the only remaining possibility is

`integral_0^T [J_inner(delta,R,t)]_+ dt = infinity`

for every `0<delta<R`.

Equivalently, a finite-time maximum-vorticity blow-up necessarily drives positive Hodge transaction action into arbitrarily small physical neighborhoods of the maximum:

`for every delta>0, productive transaction radii below delta are unavoidable.`

This is an inward cascade forced by the PDE itself, not a cascade assumed by an analyst.

## AUTOPSY

### Viscosity as a positive source of the maximum

Killed.  At a spatial vorticity maximum the viscous contribution to the magnitude equation is nonpositive.

### A fixed nonlocal/harmonic source can provide infinite finite-time growth

Killed by the harmonic Hodge occupancy floor at every fixed physical radius.

### A fixed annulus can provide infinite productive action

Killed, conditional only on the already-sharp shell projector: its action is `L2_t`-controlled by the genuine viscous enstrophy budget.

### This proves regularity

No.  It proves a necessary inward productive-transaction cascade if blow-up exists.  The remaining obstruction is temporal/material reuse: the same circulation ancestry may try to regenerate productive carriers at successively smaller radii.

## PROMOTE if the sharp projector survives remote attack

1. Finite-time vorticity blow-up forces infinite positive strain action at spatial vorticity maxima.
2. At every fixed Hodge radius, harmonic production has finite finite-time action.
3. At every fixed annulus `[delta,R]`, productive Hodge transaction has finite finite-time action.
4. Therefore blow-up forces positive productive transaction action into `rho<delta` for every `delta>0`.

## Next frontier — ancestry reusability under compulsory inward cascade

Once inward productive scale migration is unavoidable, the regularity problem becomes much narrower:

**can a finite material circulation ancestry stock be reused through infinitely many shrinking productive Hodge scales without either viscous ancestry replacement, transaction-null excess, or irreversible closure/occupancy deformation?**

That is the next material-spacetime microscope.
