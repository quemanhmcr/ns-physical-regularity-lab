# Productive radial diffusion and Hodge-kernel self-duality

## THINK — which part of the minimum productive carrier is actually viscously active?

The sharp carrier on each sphere is

`omega_prod=-(5/3)n cross Q(r)n`.

For fixed `r`, this is a toroidal degree-two angular field.  Each Cartesian component is a spherical harmonic of degree `l=2`.  Therefore the ordinary three-dimensional Laplacian acts on the STF radial profile through

`L_2 Q=Q''+(2/r)Q'-(6/r^2)Q`.

Equivalently,

`L_2 Q=r^-4 d_r [ r^6 d_r(Q/r^2) ]`.

This formula is componentwise and does not require a preferred eigenframe of `Q`.

Numerically, the `r^2` zero mode must be certified from the structural factor `(m-2)(m+3)=0`, not from cancellation of the three order-one radial-operator terms.  Likewise a small radial defect `Q-r^2 C` should be represented directly rather than recovered by subtracting nearly equal parent profiles.

## The old tangent carrier is the unique smooth viscosity-null productive profile

The radial homogeneous equation has powers `r^2` and `r^-3`.  Smoothness at the center excludes `r^-3`, so the unique smooth kernel is

`Q(r)=r^2 C`,

with constant STF tensor `C` in radius.

This is exactly the radial profile of the self-contained tangent Hodge strain carrier discovered much earlier.

Thus its special persistence was not accidental: the minimum transaction carrier sits on the smooth `l=2` harmonic zero mode of vorticity diffusion in the interior.

This also corrects an overly naive reading of the local transaction Reynolds clock.  A small observation-radius value `r^2|Q|/nu` does not by itself imply rapid diffusion when the carrier is part of the smooth `r^2` zero mode.  Diffusion sees **radial departure of `Q/r^2`**, hence source localization/taper, not radius alone.

## Same Hodge screen inverts the radial viscous defect

Let

`C=lim_{r->0} Q(r)/r^2`

for a smooth profile.  Integrating the factorized radial operator twice gives the exact tensor identity

`Q(r)-r^2 C`

`=(r^2/5) integral_0^r [1-(rho/r)^5] (L_2 Q)(rho) d rho/rho`.

The Green weight is **the same** `1-(rho/r)^5` that appeared in the Hodge strain microscope:

`S_v(r)=integral_0^r [1-(rho/r)^5] Q(rho) d rho/rho`.

The exponent five is therefore not two unrelated accidents.  It is generated both by the degree-two Hodge boundary problem and by the two radial `l=2` Laplacian modes `r^2,r^-3`.

## Physical interpretation

The profile splits as

`Q=r^2 C + Q_def`.

- `r^2 C` is the smooth viscosity-null productive mode;
- `Q_def` is the radial localization/deformation needed to depart from that mode;
- `Q_def` is exactly a Hodge-screened accumulation of `L_2 Q`.

Since the sharp angular projector is the toroidal `l=2` projection, the Laplacian preserves this sector.  Thus `nu L_2 Q` is precisely the transaction tensor carried by the viscous diffusion of the minimum productive component.

This is the first direct representation-level bridge between the inward transaction profile and a viscous-active ancestry channel.

## ATTACK — can localization hide in a sharper taper?

Use

`q(r)=q0 (r/L)^2 [1-(r/L)^p]`.

The inner part approaches the viscosity-null `r^2` mode as `p` grows, while the profile is forced to vanish at `r=L`.  Exact computation gives

`L_2 q=-p(p+5) q0 (r/L)^p/L^2`.

The squared productive viscous defect integrated over the ball is proportional to

`(q0^2/L) p^2(p+5)^2/(2p+3)`,

which grows like `p^3 q0^2/L`.  Moreover the fraction of this defect outside `alpha L` is

`1-alpha^(2p+3)`.

Sharper taper therefore moves the viscous-active defect into a thinner outer collar while increasing its amplitude; it does not erase the defect.

This is a geometric/viscous localization statement, not yet a finite-resource inequality because global energy does not control palinstrophy.

## Connection to compulsory inward cascade

The maximum-vorticity theorem says a blow-up candidate cannot obtain divergent positive action from any fixed outer annulus.  A single smooth zero-mode `Q=r^2 C(t)` extending through a fixed annulus would produce inner and outer transactions with the same instantaneous tensor coefficient.  If its positive coefficient tried to supply divergent inner action, the fixed outer annulus would also see divergent action unless another radial component canceled it.

Hence an inward blow-up cascade cannot remain purely in the smooth viscosity-null radial mode.  Either the zero-mode contribution stays action-finite, or cancellation/localization requires a nonzero `Q_def` of comparable action.  In both cases the divergent inward mechanism is forced toward the radial sector detected by `L_2`.

This is not yet a lower bound on total viscous ancestry replacement, but it identifies the correct physical defect variable.

## PROMOTE / KILL

PROMOTE:
1. `L_2` is the exact radial viscous operator of the minimum productive carrier.
2. `Q=r^2 C` is its unique smooth radial zero mode.
3. `Q-r^2 C` is reconstructed from `L_2 Q` by the same Hodge screen.
4. Source localization necessarily creates a viscous-active radial defect; sharpening a taper concentrates rather than removes it.

KILL / DEMOTE:
1. Do not use the observation-radius transaction Reynolds number alone as a universal diffusion theorem.
2. Do not identify `int |L_2 Q|^2` with a known finite NS resource; no such global bound has been proved here.
3. Do not claim the inward cascade contradiction is closed yet.  The remaining question is whether the required radial defect has a nonreusable material-spacetime ancestry cost.

## Next frontier — spacetime current of the productive radial defect

The next step is to project the exact material-spacetime ancestry law onto this toroidal `l=2` sector and ask whether repeated creation of `Q_def` at shrinking source scales can occur with finite viscous ancestry flux.  That is now a sharply defined operator-level question rather than a generic norm-growth problem.
