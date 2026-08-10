# Frozen-flux persistence: stretch -> thin -> leak

## THINK — follow deformation of the same circulation ancestry

The previous microscopes separated two objects that should not be conflated:

- an instantaneous vorticity-flux line carries the signed Hodge strain transaction;
- a material circulation collar carries ancestry when viscosity is present.

There is still a missing branch.  A singular cascade does not have to replace its ancestry.  It can try to keep the **same** material circulation flux and geometrically reuse it at ever larger vorticity.

That escape has an intrinsic kinematic consequence before any estimate is imposed.

In the inviscid/frozen-flux limit, let a material vortex-tube element carry circulation `Gamma`, vorticity magnitude `rho`, cross-sectional area `A_perp`, and material length `ell`.  Flux conservation and incompressibility give

`Gamma = rho A_perp = constant`,

`A_perp ell = constant`.

If the same ancestry amplifies by

`Lambda = rho/rho0`,

then necessarily

`ell/ell0 = Lambda`,

`A_perp/A_perp0 = 1/Lambda`.

Thus vorticity amplification is not merely "large magnitude".  It is simultaneously **material line stretch and transverse area collapse**.

For transverse principal widths `d1,d2`, their product scales like `1/Lambda`.  Hence at least one direction satisfies

`d_min/d_min,0 <= Lambda^(-1/2)`

up to the initial cross-sectional shape factor.  Axisymmetric stretching is the most favorable case for keeping the smallest width large: both widths shrink exactly as `Lambda^(-1/2)`.  Anisotropy can only make one transverse direction thinner sooner.

The viscous crossing clock of a collar thickness `delta` is

`tau_nu ~ delta^2/nu`.

Therefore, on a persistent frozen-flux branch,

`tau_nu(Lambda) <= const * tau_nu(1)/Lambda`.

The same stretching that amplifies the lineage automatically accelerates the viscous mechanism capable of mutating its ancestry.

This is the missing connection between the early stretch-thin-diffuse picture and the later material-circulation microscope: **thinning is not an auxiliary model assumption; it is the geometric dual of frozen vorticity-flux amplification.**

## PREDICT — exact linear Euler calibration

Use the exact spatially linear incompressible Euler family

`u = A(t) x`,

with

`A(t) = [[-a/2, -W(t)/2, 0],`
`        [ W(t)/2, -a/2, 0],`
`        [ 0,       0,      a]]`,

where `a>0` is constant and

`W'(t)=a W(t)`.

The antisymmetric part of `A' + A^2` vanishes exactly, so the remaining matrix is symmetric and can be absorbed into a quadratic pressure.  The vorticity is

`omega = W(t) e_z`,

with

`W(t)=W0 exp(a t)`.

Set

`Lambda=exp(a t)`.

A material axial line element obeys

`ell(t)=ell0 Lambda`,

while every transverse radius obeys

`r(t)=r0 Lambda^(-1/2)`.

For a material circular flux tube,

`Gamma(t)=pi r(t)^2 W(t)=pi r0^2 W0`.

Its volume is also exactly conserved:

`pi r(t)^2 ell(t)=pi r0^2 ell0`.

This is a direct calibration of the ancestry mechanism: arbitrarily large inviscid vorticity amplification of the same flux requires exactly the same factor of material-length growth and inverse factor of transverse-area collapse.

## PREDICT — isolated circulation collar has a dual law

Suppose the material tube is surrounded by an isolated annular collar with inner radius `r(t)`, outer radius `b(t)`, and both radii carried by the same axisymmetric transverse contraction.  Then `b/r` stays fixed and the harmonic circulation energy of a material segment is

`E_h(t) = Gamma^2 ell(t)/(4 pi) log(b/r)`.

Hence

`E_h(t)/E_h(0)=Lambda`.

The collar thickness

`delta(t)=b(t)-r(t)`

obeys

`delta(t)^2/delta0^2 = 1/Lambda`,

so its viscous crossing clock satisfies

`tau_nu(t)/tau_nu(0)=1/Lambda`,

where `tau_nu=delta^2/nu`.

Therefore the exact canonical duality is

`E_h(Lambda) = Lambda E_h0`,

`tau_nu(Lambda) = tau_nu0/Lambda`,

and

`E_h(Lambda) tau_nu(Lambda) = E_h0 tau_nu0`.

This product identity is not proposed as a universal invariant.  Its value is structural: **persistent amplification cannot simultaneously keep an isolated circulation collar cheap in kinetic occupancy and slow to viscously penetrate.**

Define the ancestry-survival clock relative to one strain time `1/a`:

`Chi = a delta^2/nu`.

On this frozen-flux branch,

`Chi(Lambda)=Chi0/Lambda`.

For every finite initial `Chi0`, sufficiently large amplification reaches `Chi=O(1)`, where transverse diffusion acts on the same timescale as stretching.

## ATTACK I — compare with the exact viscous Burgers balance

The steady Burgers vortex uses the same extensional rate `a` but includes viscosity exactly.  Its Gaussian core scale is

`delta_B^2 = 4 nu/a`.

Therefore

`a delta_B^2/nu = 4`.

This is precisely the order-one ancestry-survival clock predicted by the frozen-flux thinning branch.  The exact NS solution does not continue inviscid thinning indefinitely; viscosity arrests it at the stretch-diffuse balance.

If an inviscid material collar began with thickness `delta0`, its formal frozen-flux trajectory would reach the Burgers scale at

`Lambda_cross = a delta0^2/(4 nu) = Chi0/4`.

No length scale remains after this comparison.  The crossover is controlled by the ratio of the initial viscous crossing time to the strain time.

The old Burgers `Gamma^2` toll now acquires a clearer ancestry interpretation: once persistent stretching drives a circulation carrier to the critical transverse clock, viscosity is no longer a distant tax.  It acts on the same physical time as the amplification mechanism.

## ATTACK II — anisotropic cross-sections

Escape attempt: avoid diffusion by flattening the cross-section anisotropically instead of shrinking it isotropically.

For a material tube element, transverse area still falls like `1/Lambda`.  If transverse linear factors are `lambda1,lambda2`, then

`lambda1 lambda2 = 1/Lambda`.

Consequently

`min(lambda1,lambda2) <= Lambda^(-1/2)`.

So at least one transverse direction becomes no thicker than in the axisymmetric case.  Isotropic diffusion can use that shortest direction.  Anisotropy does not provide a kinematic route to keep every material thickness large.

What is **not** proved is that the relevant circulation-isolation distance in a highly folded 3D collar equals this infinitesimal singular width.  Folding and topology can change which path viscosity must cross.  That remains a geometric branch to attack, not something to hide inside a norm.

## ATTACK III — make Gamma tiny

The harmonic occupancy floor scales like `Gamma^2`, so a cascade can try to use ever smaller circulation lineages.

But on the **same material ancestry**, Euler stretching cannot shrink `Gamma`; Kelvin circulation is frozen.  A decrease of `Gamma` is therefore not a persistence maneuver.  It moves the event into the viscosity-driven ancestry-mutation branch already isolated by the Kelvin microscope.

Across different lineages, arbitrarily small `Gamma` remains an escape candidate.  The exact Oseen survival calculation shows that `Gamma/nu -> 0` donors lose strain before one strain time in that canonical NS geometry.  A fully 3D self-regenerating small-circulation donor has not been ruled out.

## ATTACK IV — reuse the same ancestry repeatedly

The collar energy is an **occupancy** resource, not automatically an irreversible expenditure.  A material lineage could stretch, later compress, and then stretch again.  It would be wrong to sum the same harmonic energy on every visit.

However, a genuine finite-time singularity requires an unbounded instantaneous amplification along some sequence.  For one fixed isolated ancestry with nondegenerate collar geometry,

`E_h ~ Lambda`

already prevents `Lambda -> infinity` under a finite instantaneous kinetic-energy budget.  Repeated bounded oscillations do not defeat that specific statement.

The remaining escape is to abandon that ancestry/collar before `Lambda` becomes large.  That is exactly where material circulation mutation, collar penetration, cancellation, or reconnection must enter.

## ATTACK V — local segment versus whole vortex line

A high-vorticity event need only stretch a local material tube segment.  The entire closed line need not lengthen by the same factor.

The persistence statement should therefore be charged to the **material segment whose flux supplies the productive transaction**, not automatically to a whole closed tube.  Its local volume conservation still gives transverse-area collapse when its length grows.  But a curved or folded collar needs a geometric lower bound replacing the straight-annulus formula before the harmonic energy growth can be made universal.

This is an important remaining gap.

## AUTOPSY — a sharper physical dichotomy

The current picture is no longer "persistence costs Gamma^2, switching costs Gamma^2" as two unrelated canonical observations.

A more natural mechanism is:

### Persistence branch

The same circulation ancestry supplies productive strain and remains isolated.

Then frozen-flux amplification forces

- line stretch `~ Lambda`,
- transverse area `~ 1/Lambda`,
- at least one material width `<= Lambda^(-1/2)`,
- isolated harmonic collar occupancy `~ Gamma^2 Lambda` in the canonical tube geometry,
- viscous crossing clock `<= const/Lambda`.

Persistence therefore drives itself toward either an energy-occupancy obstruction or a viscous-isolation breakdown.

### Renewal branch

Before that happens, the flow replaces or cancels the ancestry.

Then either material circulation changes, which only viscosity can do, or the isolation collar is penetrated by competing vorticity.  The exact Oseen modules show a localized `Gamma^2` spacetime toll for fixed-fraction circulation leakage in one canonical NS realization.

The missing theorem is the bridge saying that every general 3D escape from the persistence branch enters a quantitatively chargeable renewal event without double counting.

## PROMOTE / KILL

### PROMOTE

1. **Frozen-flux amplification has an intrinsic transverse-thinning dual.**  For a persistent material vortex-tube element, `rho/rho0=Lambda` forces cross-sectional area `1/Lambda`.
2. **At least one transverse material scale shrinks no slower than `Lambda^(-1/2)`.**  Anisotropy cannot keep both directions thick.
3. **Canonical isolated-collar dual law:** in exact axisymmetric Euler stretching, harmonic circulation energy grows like `Lambda` while the viscous crossing time falls like `1/Lambda`.
4. **Burgers criticality matches the ancestry clock:** the exact steady NS core satisfies `a delta_B^2/nu=4`, an order-one stretch/diffusion balance.
5. **A change of Gamma is not persistence.**  On the inviscid ancestry branch circulation is frozen; reducing circulation necessarily invokes viscosity/renewal.

### KILL / DO NOT PROMOTE

1. Do not sum harmonic collar occupancy repeatedly as if it were irreversible dissipation.
2. Do not claim the straight-annulus energy factor for arbitrary folded 3D collars without a geometric capacity/inductance argument.
3. Do not claim transverse singular-value thinning alone proves a quantitative circulation leak; it only exposes the shortening viscous route.
4. Do not claim the Oseen `Gamma^2` renewal toll is already universal in 3D.
5. No regularity contradiction has yet been obtained.

## Next frontier — geometric circulation inductance under deformation

The natural next object is the kinetic energy of the harmonic circulation class on a **deformed finite-thickness material collar**.  For a straight annulus it is

`E_h = Gamma^2 ell/(4 pi) log(b/a)`.

For a general curved/folded collar, define the circulation inductance `L_C` by the minimum kinetic energy among divergence-free fields carrying unit circulation through the collar:

`E_min(Gamma,C) = (Gamma^2/2) L_C`.

If `L_C` can be related from below to the material stretch and to the shortest diffusion path, then the persistence-or-renewal dichotomy would no longer depend on a straight tube.  That is the next geometric bottleneck worth attacking.
