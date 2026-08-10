# Transaction-to-closure / vortex-lineage microscope

## THINK — follow the object that actually carries vorticity

The Hodge transaction microscope identified the exact signed vorticity observable that supplies the vortical center strain,

`S_v(r) = integral_0^r [1-(rho/r)^5] Q(rho) d rho/rho`.

Its main autopsy also exposed a structural defect in any same-ball geometry-tax argument: a locally straight vortex segment can transact strain inside `B_r` while having zero local curvature/direction tax there.  Because `div omega=0`, that segment cannot end.  Its necessary turning and closure may lie elsewhere.

The physical object that links the gain region to the tax region is therefore not a shell, a norm, or an arbitrary dyadic block.  It is the **connected vorticity-flux lineage itself**.

Write `omega=varpi xi` away from zeros.  Infinitesimal vortex tubes carry circulation flux `dGamma`; instantaneously, `div omega=0` gives conservation of that flux along the tube.  In flux-line coordinates,

`varpi dV = dGamma ds`,

where `s` is arclength along the vortex line.  Substituting this identity into the exact Hodge-screened Biot-Savart representation yields a line functional rather than a magnitude bound.

For a vortex line `gamma` through `x`, with `dx = xi ds`, define

`K_r[gamma] = [3/(8 pi)] [1-(|x|/r)^5] / |x|^3`

`             * [ n tensor (n cross dx) + (n cross dx) tensor n ]`,

with `n=x/|x|`, and set the kernel to zero outside `B_r`.  Then, at the level of an infinitesimal flux-tube disintegration,

`S_v(r) = integral dGamma integral_gamma K_r[gamma]`.

The vorticity magnitude has disappeared from the geometric factor.  It survives only as circulation weight `dGamma`.  The rest is the actual oriented path taken by the vortex line.

For the singularity-relevant direction `e`, the scalar transaction of one lineage is

`T_e^r[gamma,dGamma]`

` = [3 dGamma/(4 pi)] integral_gamma [1-(|x|/r)^5]`

`       * (e.x) e.(x cross dx) / |x|^5`.

In cylindrical coordinates around `e`, where `z=e.x`, `R` is cylindrical radius and `phi` is azimuth,

`e.(x cross dx) = R^2 dphi`,

so

`T_e^r = [3 dGamma/(4 pi)] integral [1-(|x|/r)^5]`

`                  * z R^2/(R^2+z^2)^(5/2) dphi`.

This is the **productive winding one-form**.  It is parameterization-free.  It separates the two pieces nature actually uses:

- `dGamma`: how much vorticity flux the lineage carries;
- the oriented spatial winding of that lineage around the target direction.

No absolute value is inserted before the signed transaction is formed.

A useful extra feature is inherited from the Hodge screen: the one-form vanishes on `partial B_r`.  Thus a line may transact strain while crossing the Hodge ball and then carry all of its required closure geometry outside the ball.  That is not a loophole in the representation; it is precisely what the representation says must happen.

## PREDICT — a straight donor should be productive with zero local curvature

Take `e=e_z` and the straight line

`gamma(s)=(d,s,z)`,

with `h^2=d^2+z^2<r^2`.  Inside `B_r`, `|s|<=L`, `L=sqrt(r^2-h^2)`, and

`e.(x cross dx)=d ds`.

The exact Hodge-screened lineage transaction is therefore

`T_e^r = [3 Gamma d z/(4 pi)]`

` * { 2/h^4 [ L/r - L^3/(3r^3) ] - 2L/r^5 }`.

For an infinite unscreened straight filament this tends to

`T_e^infinity = Gamma d z / [pi (d^2+z^2)^2]`.

Yet the centerline curvature is exactly zero throughout the productive segment.

This gives a sharper statement than the old open-filament test: **the exact Hodge transaction itself can be nonzero on a zero-curvature segment.**  Any geometry tax that is genuinely forced by `div omega=0` must therefore be transported along the same vortex lineage to its continuation/closure.

## ATTACK I — can closure curvature be made arbitrarily cheap?

For a coherent thin flux tube, let

`kappa = |(xi.grad) xi|`

be centerline curvature.  The exact directional viscous term contains

`nu integral varpi |grad xi|^2 dV`

and therefore, using `|grad xi|^2 >= |(xi.grad)xi|^2` and `varpi dV=dGamma ds`,

`Tax_dir >= nu integral dGamma integral kappa^2 ds`.

For a single closed flux `Gamma` following a closed centerline of length `ell`, Fenchel plus Cauchy gives

`integral kappa^2 ds >= 4 pi^2/ell`,

hence

`Tax_dir >= 4 pi^2 nu Gamma/ell`.

This is real, but by itself it does **not** solve the problem.  A lineage can make `ell` arbitrarily large and drive the curvature tax toward zero.

An explicit stadium closure makes the escape physical rather than formal.  Keep one straight side passing through the Hodge ball and place the return side far away.  Join the two sides by two semicircles of radius `R`.  The productive straight chord in the ball is unchanged, while all turning occurs outside.  The two bends have

`integral_bends kappa^2 ds = 2 pi/R`,

so their directional tax is

`Tax_bend >= 2 pi nu Gamma/R -> 0` as `R->infinity`.

Therefore the claim

> productive Hodge transaction forces a scale-independent curvature tax

is **killed**.

## ATTACK II — the large-closure escape carries a circulation collar

The long-closure escape changes another physical quantity: the amount of velocity field that must remain wrapped around the circulation.

Define a **circulation-isolation collar** around a tube segment.  Let `a` be a core radius and `b>a` a collar radius.  Suppose every meridional loop `C_rho`, `a<=rho<=b`, encloses at least a fraction `theta` of the lineage circulation:

`|integral_{C_rho} u.dl| >= theta Gamma`, `0<theta<=1`.

This definition is physical: it asks what a family of actual circulation loops measures.  It does not require a Fourier support decomposition.

For a straight collar, Cauchy on each meridional circle gives

`integral_{C_rho} |u|^2 dl >= theta^2 Gamma^2/(2 pi rho)`.

Integrating across the annulus and along length `ell` gives the exact collar floor

`E_collar >= theta^2 Gamma^2 ell/(4 pi) log(b/a)`.

For a curved tube with `kappa b<1`, the tubular Jacobian contributes the conservative factor `1-kappa b`.  On a stadium bend of radius `R>b`, total bend length `2 pi R`, hence

`E_bend_collar >= [theta^2 Gamma^2 R/2] (1-b/R) log(b/a)`.

Multiplying by the bend curvature tax removes the artificial closure radius:

`Tax_bend * E_bend_collar`

` >= pi nu theta^2 Gamma^3 (1-b/R) log(b/a)`.

This is the first closure law in which making the vortex lineage longer does not create a free asymptotic parameter: the curvature tax falls as `1/R`, while the unavoidable circulation-collar energy rises as `R`.

If total kinetic energy is bounded by `E0` and the bend is safely tubular, `b/R<=epsilon<1`, then necessarily

`R <= 2 E0/[theta^2 Gamma^2 (1-epsilon) log(b/a)]`,

and therefore

`Tax_bend >= pi nu theta^2 Gamma^3 (1-epsilon) log(b/a)/E0`.

This statement is rigorous for the declared collar geometry.  It is not yet a universal theorem for arbitrary intertwined vorticity.

## ATTACK III — shrink the circulation to escape the collar law

The closure-collar lower bound scales like `Gamma^3`.  A singular cascade might therefore try to use ever smaller circulation while moving the donor ever closer to the target so that `Gamma/d^2` still produces large strain.

But `Gamma/nu` is not a bookkeeping ratio.  It is the circulation Reynolds number, and exact viscous dynamics tests whether such a donor can survive for even one strain time.

Use an exact Lamb-Oseen vortex.  At distance `d`, with

`q=d^2/(4 nu t)`, `Re_Gamma=Gamma/nu`,

its shear-strain magnitude is

`s(q) = Gamma/[2 pi d^2] * F(q)`,

`F(q)=1-(1+q) exp(-q)`.

Start at `q0` and wait one initial strain time `Delta t=1/s(q0)`.  The new similarity coordinate is exactly

`q1 = q0 / [1 + 8 pi q0/(Re_Gamma F(q0))]`,

so the surviving strain fraction is

`s1/s0 = F(q1)/F(q0)`.

As `Re_Gamma->0`,

`q1 ~ Re_Gamma F(q0)/(8 pi)`,

and

`s1/s0 ~ Re_Gamma^2 F(q0)/(128 pi^2) -> 0`.

Thus the small-circulation escape fails in this exact NS diffusion model: a donor with vanishing `Gamma/nu` loses essentially all of its strain before one strain e-fold can occur.

The same dimensionless ratio appears directly at a Hodge scale.  A strain transaction `sigma` acting across distance `r` has circulation dimension

`Gamma_H = sigma r^2`,

and

`Gamma_H/nu = (r^2/nu)/(1/sigma)`

is exactly the ratio of the viscous crossing time of scale `r` to the strain time.  The canonical Burgers balance likewise places the viscous core at `r^2 sigma/nu=O(1)`.

This strongly suggests a physical **survival gate**: productive circulation much below `nu` cannot remain a coherent strain donor for a full amplification time unless some simultaneous 3D mechanism regenerates it.  That final clause is essential; the Oseen calculation alone does not prove the gate for every 3D Navier-Stokes geometry.

## AUTOPSY — attack the new mechanism before trusting it

### 1. Huge closure

Escape attempt: `R->infinity` to make curvature tax vanish.

Result: killed **within an isolated circulation collar**.  The collar energy grows linearly with `R`; finite `E0` converts that growth back into a positive curvature-tax floor.

### 2. Collapse the collar

Escape attempt: `b/a->1` or `theta->0`.

Result: survives.  The factor

`theta^2 log(b/a)`

can vanish.  Physically this means the lineage is no longer surrounded by an annulus carrying its own circulation cleanly: opposite/competing flux, self-approach, or another branch has entered at essentially the core scale.  The correct next object is therefore not a stronger collar inequality; it is the **event by which circulation isolation is lost or renewed**.

### 3. Shrink circulation

Escape attempt: `Gamma->0` while donor distance shrinks so `Gamma/d^2` stays large.

Result: strongly obstructed but not universally killed.  Exact Oseen dynamics gives `s1/s0=O((Gamma/nu)^2)` after one strain time as `Gamma/nu->0`.  A general 3D flow could try to use stretching/compression to regenerate the donor while viscosity diffuses it.  Any such rescue must itself be part of the transaction ledger.

### 4. Opposite-sign cancellation / paired tubes

Escape attempt: pack opposite circulation close enough that large-loop circulation and collar energy cancel.

Result: open.  If the pair is separated by a genuine collar, each lineage pays its own circulation floor.  If separation collapses to the core scale, `theta` or `log(b/a)` collapses and the present law deliberately stops claiming a bound.  This is exactly where the earlier compactness/cancellation physics and the present lineage picture must be joined.

### 5. Beltrami / helical coherence

Nothing in the lineage representation promotes direction complexity by itself.  The productive winding one-form must still be nonzero.  Beltrami zero-stretch calibrations therefore remain null controls, not escapes.

### 6. Reconnection and material-lineage switching

Instantaneous vortex lines are not exact material labels when `nu>0`.  Viscosity permits flux to cross material surfaces and reconnect line geometry.  The old exact Oseen switch module already shows a scale-independent `Gamma^2` dissipation toll for fixed-fraction circulation leakage in a canonical material core.  The present result says where that module must now attach: **collar loss/renewal is the physical switching event that can defeat closure accounting.**

The remaining hard problem is to prove a non-double-counted spacetime alternative:

- either a productive circulation lineage remains isolated long enough to be charged by closure/collar geometry,
- or its circulation isolation changes by an `O(1)` fraction and that switch is charged by viscous flux leakage / reconnection,
- while repeated use of the same donor must not be charged twice unless a genuinely new amplification transaction occurs.

That is a much more constrained target than a generic norm estimate.

## PROMOTE / KILL

### PROMOTE

1. **Productive winding one-form.**  The exact Hodge vortical strain can be disintegrated along connected vorticity-flux lineages; circulation is the weight, and oriented line geometry is the transaction.
2. **Spatial displacement is structural.**  A straight zero-curvature segment can carry a nonzero Hodge strain transaction while its required turning/closure lies outside the Hodge ball.
3. **Closure alone is insufficient.**  Curvature tax can be driven to zero by increasing closure length.
4. **Circulation collar tradeoff.**  For an isolated tube collar, making closure large exchanges curvature tax for kinetic-energy occupancy, with the stadium product floor
   `Tax_bend E_bend_collar >= pi nu theta^2 Gamma^3(1-b/R) log(b/a)`.
5. **Small-circulation survival obstruction in exact NS diffusion.**  In Lamb-Oseen, a strain donor with `Gamma/nu->0` loses its strain over one initial strain time like `O((Gamma/nu)^2)`.

### KILL / DO NOT PROMOTE

1. Kill any universal claim `productive transaction => local same-ball direction tax`.
2. Kill any universal claim `closed lineage => scale-independent curvature tax` without a finite-resource/collar condition.
3. Do not promote the collar law through `b/a->1` or `theta->0`; that is exactly the unresolved cancellation/reconnection branch.
4. Do not promote `Gamma_H/nu >= c` as a general theorem yet.  It is a physically sharp survival gate supported by exact Oseen/Burgers balances, but a genuinely 3D self-regenerating donor has not been excluded.
5. No regularity contradiction has yet been obtained.  The remaining frontier is a **spacetime lineage-renewal microscope** that distinguishes persistence from replacement without double charging.

## Next frontier

The next representation should attach a material circulation collar to each productive Hodge-lineage event and observe only two physical changes:

- persistence: how long the same isolated circulation flux remains capable of the signed winding transaction;
- renewal: how much circulation crosses the material collar when the productive donor is replaced, reconnected, or cancelled.

The goal is a spacetime ancestry ledger in which every productive strain e-fold either occupies finite energy for a definite duration or irreversibly spends a fixed fraction of circulation through viscosity.  Only after that non-reuse structure is explicit should one ask whether an infinite singular cascade demands infinite total resource.
