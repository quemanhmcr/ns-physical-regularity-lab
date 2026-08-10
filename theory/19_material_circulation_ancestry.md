# Material circulation ancestry and the collar-cohomology microscope

## THINK — instantaneous vortex lines are not the ancestry object

The transaction-to-closure microscope follows the correct instantaneous carrier: connected vorticity flux.  But viscosity means vortex lines are not material labels.  They may reconnect, exchange partners, or cease to be meaningful near a vorticity zero.

There is, however, an older physical object built directly into Navier–Stokes that survives this ambiguity: **circulation around a material loop**.

Let `C(t)=X(t,C0)` be a closed loop transported by the fluid and

`Gamma_C(t)=integral_{C(t)} u.dl`.

For incompressible Navier–Stokes,

`D_t u = -grad p + nu Delta u`,

and the exact material-loop identity is

`d Gamma_C/dt = nu integral_{C(t)} Delta u.dl`

`             = -nu integral_{C(t)} (curl omega).dl`.

The pressure contribution is exact, and the entire inviscid Euler transport/stretching part cancels from the circulation derivative.  This is not a chosen conservation law: it is the physical ancestry rule of the equation.

**Stretching can amplify vorticity magnitude and alter geometry without changing material circulation.  Only viscosity mutates the circulation ancestry of a material loop.**

This suggests replacing fragile line identity by a material circulation genealogy.

## PREDICT — Lamb–Oseen should calibrate the ancestry rule exactly

For the exact Lamb–Oseen vortex,

`u_theta(R,t) = Gamma/(2 pi R) [1-exp(-q)]`,

`omega_z(R,t) = Gamma/(4 pi nu t) exp(-q)`,

`q = R^2/(4 nu t)`.

A circle of fixed radius `R` is a material loop because the velocity is purely azimuthal.  Its material circulation is

`Gamma_R(t)=Gamma[1-exp(-q)]`.

Therefore

`d Gamma_R/dt = -Gamma q exp(-q)/t`.

Meanwhile

`(curl omega)_theta = -partial_R omega_z`

`                   = omega_z R/(2 nu t)`,

so the viscous Kelvin defect gives

`-nu integral_C curl(omega).dl`

` = -nu 2 pi R [omega_z R/(2 nu t)]`

` = -Gamma q exp(-q)/t`.

The two routes are exactly identical.  This gives an exact NS calibration of material circulation ancestry, independent of any vortex-line labeling convention.

## ATTACK I — one material loop is too thin an observer

A tempting next step would be to argue that an `O(Gamma)` circulation change of one material loop must consume an `O(Gamma^2)` amount of ordinary energy dissipation.  That is false as a static statement.

Take cylindrical coordinates around the `z` axis and a reference circle `C: r=R,z=0`.  Let a thin axisymmetric toroidal swirl be

`u_delta = A_delta (r/R) Phi((r-R)/delta,z/delta) e_theta`,

where `Phi` is a smooth bump supported in `|r-R|<delta`, `|z|<delta`, with

`Phi(0,0)=1`, `grad Phi(0,0)=0`, `Delta_{X,Z} Phi(0,0) != 0`,

and `delta<R/2`.  Pure axisymmetric swirl is divergence-free, and the support stays away from the symmetry axis.

For the theta component,

`(Delta u)_theta = u_rr + r^{-1}u_r + u_zz - r^{-2}u`.

At the reference circle the `r^{-1}u_r` and `r^{-2}u` terms cancel because of the factor `r/R`, while the first derivatives of `Phi` vanish.  Hence

`(Delta u_delta)_theta|_C`

` = A_delta delta^{-2} Delta_{X,Z}Phi(0,0)`.

Thus the instantaneous Kelvin defect is

`J_delta := integral_C Delta u_delta.dl`

` = 2 pi R A_delta delta^{-2} Delta Phi(0,0)`.

Choose `A_delta ~ delta^2`.  Then `J_delta` is `O(1)` as `delta->0`.

But the vorticity is only `O(A_delta/delta)` in a toroidal volume `O(R delta^2)`, so

`integral |omega_delta|^2 dV = O(A_delta^2 R) = O(R delta^4) -> 0`.

Therefore an arbitrarily small-enstrophy smooth neighborhood can produce a finite **instantaneous** Kelvin circulation derivative on a single curve.

This kills any scale-free static inequality of the form

`|d Gamma_C/dt| <= F(nu integral |omega|^2)`

for a single material loop.

The failure is physically informative.  A loop is codimension two; a second derivative can be concentrated onto its tiny neighborhood while a volume energy observer barely sees it.  The correct ancestry observer must have thickness.

This attack does **not** show that an `O(Gamma)` circulation change over a finite physical time is free.  The thin defect layer diffuses on its own time `tau_delta ~ delta^2/nu`.  A finite instantaneous derivative that lives only for `tau_delta` produces a vanishing integrated change as `delta->0`.  The unresolved quantity is therefore intrinsically spacetime, not static.

## ATTACK II — thicken the observer into a circulation collar

Consider a straight annular collar

`A_{a,b} x [0,ell]`, `a<rho<b`,

and suppose every meridional circle has the same circulation `Gamma`.

The canonical harmonic circulation field is

`h_Gamma = Gamma/(2 pi rho) e_theta`.

It is divergence-free and curl-free in the annulus.  Write

`u = h_Gamma + v`,

where `v` has zero circulation on every meridional circle.  Then the angular average of `v_theta` vanishes at each `rho`, and the kinetic-energy cross term is exactly zero:

`integral h_Gamma.v dV = 0`.

Therefore

`E_collar = E_harm + E_zero`,

with

`E_harm = Gamma^2 ell/(4 pi) log(b/a)`

and `E_zero >= 0`.

This recovers the earlier circulation-collar floor, but now exposes its actual physical structure: it is the energy of the **harmonic cohomology mode** carried by the noncontractible circulation cycle of the annulus.

The collar is therefore more than a convenient family of Cauchy loops.  It is a topological/Hodge detector of circulation ancestry.

### What counts as a real renewal event?

As long as the material annulus remains an honest collar and its family of material loops retains their circulation, the harmonic ancestry mode persists even if the instantaneous vortex centerline bends, stretches, or reconnects elsewhere.

To erase or replace that ancestry, nature must do at least one of two physical things:

1. **viscously mutate the circulation of a material family of loops**, according to the exact Kelvin defect; or
2. **breach the collar geometry itself** by bringing vorticity/cancellation/contact through the annular region so that the noncontractible circulation class is no longer isolated.

This is a sharper renewal definition than "the vortex line changed identity".

## ATTACK III — can a local breach destroy global ancestry cheaply?

Yes, this remains the dangerous branch.

A connected annular collar can lose its topological isolation through a small spatial gate.  Instantaneous reconnection is local, and the single-loop thin-layer construction warns that trace quantities can change very rapidly in a tiny neighborhood.  Therefore the harmonic collar floor alone does not prove a nonreusable cost per renewal.

The missing structure is a **parabolic breach law**: a localized gate of transverse size `delta` has a viscous lifetime/penetration time of order `delta^2/nu`.  To change an `O(1)` fraction of a finite-thickness circulation family, a breach must either

- persist long enough for diffusion to cross a finite material thickness,
- sweep across a finite portion of the collar,
- or be continually regenerated by 3D stretching/advection.

Each possibility uses spacetime extent.  A static norm cannot see this correctly.

This is where the Oseen survival gate and the material Kelvin rule meet naturally: the former is an exact parabolic example of a donor losing coherence; the latter identifies the only term that can mutate circulation ancestry.

## AUTOPSY — what the new representation does and does not buy

### Vortex-line reconnection

Instantaneous line connectivity can change, so line labels are not promoted as material identities.  Material circulation survives as the ancestry observer.

### Single-loop concentration

Killed as a sufficient observer.  A codimension-two loop admits an `O(1)` Kelvin defect with vanishing enstrophy by concentrating second derivatives in a thin toroidal layer.

### Finite-thickness collar

Promoted.  Its circulation is represented by an exact harmonic mode with an exact orthogonal energy floor.  This is robust against zero-circulation fluctuations inside the collar.

### Collar breach

Open.  A local gate can destroy global isolation.  No universal finite energy-dissipation toll has been proved for such a spacetime event.

### Material-lineage switching

Reformulated.  "Switching" should no longer mean choosing a new instantaneous centerline.  A true ancestry renewal is either material circulation mutation or collar breach.  This removes observer-dependent relabeling from the ledger.

## PROMOTE / KILL

### PROMOTE

1. **Material circulation is the natural ancestry currency.**  Euler transport/stretching cannot change it; viscosity alone appears in its exact evolution.
2. **Circulation collar = harmonic cohomology mode.**  For a straight annulus the energy split is exact and orthogonal, with `E_harm = Gamma^2 ell log(b/a)/(4 pi)`.
3. **Renewal should be defined physically as circulation mutation or collar breach**, not as a change of instantaneous vortex-line label.
4. **The next obstruction is spacetime/parabolic.**  Thin static concentration defeats single-loop energy estimates, but its lifetime shrinks like `delta^2/nu`.

### KILL / DO NOT PROMOTE

1. Kill any scale-free static bound from a single-loop Kelvin defect to ordinary enstrophy/energy dissipation.
2. Do not identify instantaneous vortex-line connectivity with material ancestry when `nu>0`.
3. Do not claim the harmonic collar energy is irreversibly spent when a donor is reused; it is an occupancy resource until a genuine breach/mutation is shown.
4. No universal cost per collar breach has yet been obtained.

## Next frontier — spacetime collar-breach microscope

The next microscope should follow a finite-thickness material annulus through a productive strain episode and record only invariant physical events:

- the signed productive Hodge transaction supplied by circulation carried through the collar;
- the harmonic circulation coefficient of the material collar;
- the first spacetime point where a fixed fraction of collar circulation is lost or the collar is pierced by competing vorticity;
- the transverse breach thickness `delta`, its lifetime, and the viscous clock `delta^2/nu`.

The target is not a static inequality.  It is a parabolic alternative: either the ancestry collar persists for an amplification time, or a breach large/long enough to replace it must occur.  Only that spacetime alternative has a chance to make nonreuse a property of Navier–Stokes itself rather than a bookkeeping convention.
