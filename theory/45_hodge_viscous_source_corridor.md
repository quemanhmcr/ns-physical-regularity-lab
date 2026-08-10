# Hodge-viscous source corridor and the failure of a material-volume toll

## THINK — where can a persistent singular source physically live?

Two independently derived scales now constrain a large-strain source.

From the half-strain Hodge transition and finite kinetic energy,

`L_E(s)=[30 E0/(pi s^2)]^(1/5)`.

Any radius at which at least half of a large center strain has become genuinely vortical lies no larger than this energy horizon.

From the localized productive spectral gap, define the unit-exposure viscous length by requiring the pure-diffusion squared-norm exponent over one production time to equal one:

`28 nu/(s L_nu^2)=1`,

so

`L_nu(s)=sqrt(28nu/s)`.

A source below this scale is strongly renewal-dominated over one strain time; a near-frozen localized source must live above it.

## Exact corridor identity

The two scales satisfy

`(L_E/L_nu)^10=s/s_c`,

where

`s_c=(28^5 pi^2/900) nu^5/E0^2`.

Also

`Re_source(L_E)=28 (s/s_c)^(1/5)`.

Thus above the critical strain `s_c`, the only remaining near-frozen source geometry lies in the intrinsic corridor

`L_nu(s) <~ L(s) <= L_E(s)`.

Below `s_c`, even the unit-exposure viscous scale is outside the finite-energy Hodge horizon, so a large localized source is necessarily renewal-dominated at that calibration threshold.

The tenth root is noteworthy but should not be mystified: it is simply the composition of the Hodge exponent `2/5` with the diffusive exponent `1/2`.

## Finite-time power-law corridor

If a candidate has the natural divergent strain clock

`s~1/tau`, `tau=T-t`,

and source scale

`L~tau^alpha`,

then:

- Hodge energy-horizon compatibility gives `alpha>=2/5`;
- high source Re gives `alpha<1/2`;
- finite cumulative localized viscous exposure gives the same strict condition `alpha<1/2`;
- the sharp productive-enstrophy spacetime integral only requires `alpha>1/3`.

Hence the scalar-admissible frozen branch occupies

`2/5 <= alpha < 1/2`.

Throughout this corridor

`Gamma_Q=sL^2~tau^(2alpha-1)->infinity`.

So the price of staying cumulatively frozen is unbounded source circulation-dimensional transaction, exactly as the source-Reynolds dichotomy predicted.

## ATTACK — material volume is not the missing finite resource

A tempting next conjecture is that a shrinking Eulerian source must consume a fixed positive amount of new material on every strain time.

This is false.

Parameterize by strain count `N` and take

`L=e^(-bN)`.

Even if the source replaces an entire ball volume `~L^3` on every increment `dN`, the total distinct material volume needed from any late `N` onward is

`integral_N^infinity L^3 dN < infinity`.

Likewise a spherical source with `|dL/dt|=b s L` sweeps boundary volume `~b s L^3 dt=~b L^3 dN`, whose infinite-cascade total is finite.

Thus an infinite sequence of shrinking source regions can be fed by a finite total material volume.  The later parcels simply become smaller.

The real requirement is much sharper: those vanishing parcels must carry an **increasing circulation flux and increasing vorticity amplitude**.  Material amount is not the ancestry currency.

## PROMOTE / KILL

PROMOTE:
1. The Hodge energy horizon and localized viscous length define an intrinsic source-scale corridor.
2. In the power-law finite-time clock, the near-frozen scalar-admissible corridor is `2/5<=alpha<1/2`.
3. Every such frozen corridor has unbounded `Gamma_Q=sL^2`.
4. Infinite source recruitment can use finite total material volume; the hard variable is the amplification/flux carried by that material, not its volume.

KILL:
1. Raw material volume per event as a scale-independent nonreusable cost.
2. A regularity proof based only on excluding high-Re scales by finite energy or enstrophy scalings.

## Next frontier — amplification conveyor of recruited ancestry

If the source recruits vanishing material parcels carrying unbounded flux density, those parcels cannot arrive as weak fresh background material.  They must already have been strongly amplified before recruitment.  The high-Re branch is therefore an **amplification conveyor**: amplified material ancestry is transported from larger/older source geometry into a smaller current source.

The next microscope should trace that conveyor backward in material labels and determine whether it can be indefinitely self-feeding without generating transaction-null folding or viscous current curvature.
