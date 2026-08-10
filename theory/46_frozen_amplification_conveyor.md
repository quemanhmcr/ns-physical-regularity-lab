# Frozen high-Re branch as a material amplification conveyor

## THINK — what does stationary productive amplification do to material geometry?

The pair ancestry cell is

`D=rho_a rho_b r T`.

A common incompressible affine action has zero determinant rate, so its contribution to `d log|D|` vanishes exactly.

If the productive Gram shape is held fixed (`Tdot=0`) and both endpoint vorticity magnitudes amplify at the same rate `s`, then

`sigma_a=sigma_b=s`,

and determinant preservation forces

`sigma_R=-2s`.

This is realized instantaneously by the stationary capacity amplifier derived earlier: its common strain gives endpoint rates `s,s` and bridge rate `-2s` while the Gram shape rates vanish.

Thus stationary productive amplification is not “nothing happens geometrically.”  It is a hyperbolic material conveyor: two vorticity directions stretch while the source bridge contracts twice as fast.

## Conditional frozen shape-lock scaling

Let

`dN=s dt`

be strain count.  If pair-cell renewal is negligible and shape remains locked, then

`rho_a/rho_a0 = rho_b/rho_b0 = e^N`,

`r/r0=e^(-2N)`,

and

`rho_a rho_b r=constant`.

If the bridge/source scale follows

`r~tau^alpha`, `tau=T-t`,

then

`N~(alpha/2) log(1/tau)`,

endpoint vorticity magnitudes scale like

`tau^(-alpha/2)`,

and the required gain rate is

`s~alpha/(2tau)`.

This is a kinematic calibration, not a claim that the instantaneous stationary strain tensor remains an exact finite-time Navier-Stokes solution.

## Recruitment is flux turnover, not volume turnover

The cumulative source-Reynolds theorem says a truly frozen infinite cascade requires source flux `Phi` to grow without bound.  In a representative scaling

`Phi~e^(aN)`, `L~e^(-bN)`,

the net positive recruitment per strain count is

`dPhi/dN=a Phi`,

which diverges, while

`integral L^3 dN`

can remain finite.

Moreover

`Phi/L^2~e^[(a+2b)N]`,

on the same scale as the source production rate `s` in the scalar escape family.

So an infinite high-Re source does not need infinite fresh volume.  It needs smaller and smaller material parcels carrying larger and larger vorticity flux density.

## Sharp localized carrier confirms the incoming material is already amplified

For the extremal localized profile

`q=q0 x^2(1-x^p)`

with `E=diag(2,-1,-1)`, the Hodge source strain is `s=C_p q0`.  The radial transaction amplitude reaches its maximum at

`x_m=[2/(p+2)]^(1/p)`.

The minimum productive vorticity has

`max_n |omega_prod|=(5/2) q`.

At the radial maximum, the ratio

`|omega_prod|_peak/s`

is order one and tends to `7` as `p->infinity`.

Thus the flux being recruited into a high-Re extremal source is not weak background circulation that becomes productive only after arrival.  It is already vorticity-amplified on the current source-strain scale.

The high-Re branch is therefore a **material amplification conveyor**: ancestry amplified at an earlier/larger stage becomes the input flux of a later/smaller stage.

## AUTOPSY

- Killed: fresh material volume as the source of an infinite cascade.
- Killed: the idea that a stationary Gram shape means no material geometric deformation.
- Promoted: shape-lock amplification necessarily couples to bridge contraction in the frozen pair cell.
- Promoted: high-Re source recruitment is multiplicative flux turnover of already amplified ancestry.

## Next frontier — can the conveyor self-feed indefinitely?

The remaining question is no longer how much raw material or circulation exists.  It is whether amplified material ancestry can be passed from stage to stage indefinitely while preserving enough sharp productive geometry.

There are two obvious failure modes to attack:

1. **folding / packing:** long stretched material lineages must be arranged inside or through an increasingly small source, potentially creating transaction-null vorticity beyond the sharp carrier;
2. **ancestry-current curvature:** any radial/source localization away from the smooth frozen zero mode activates `d j`.

The next microscope should measure how far a frozen conveyor can remain near the sharp extremal nested-loop geometry while its material line length grows and its source radius contracts.
