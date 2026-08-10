# Productive Kelvin ancestry curvature

## THINK — nonzero viscous current is not the same as ancestry replacement

The material circulation law uses the Kelvin current one-form

`j = nu curl omega . dl`,

with

`d/dt integral_{C(t)} u.dl = - integral_{C(t)} j`.

The material-spacetime vorticity law is

`partial_t Omega_hat + d j_hat = 0`.

The previous radial microscope identified `L_2 Q` as the viscous-active defect of the minimum productive carrier.  The next question is whether that defect is exactly the part of the Kelvin current that can change circulation ancestry.

## Exact current of the minimum productive carrier

For a fixed STF ray `Q=q(r)E`,

`omega_prod=-(5/3) q(r) n cross E n`.

Direct vector calculus gives

`j_r = 5 nu [q(r)/r] (n.E.n)`,

`j_t = (5nu/3) [q'(r)+q(r)/r] (I-nn) E n`.

The tensor formula follows componentwise for a general radial STF profile `Q(r)`.

## The smooth radial zero mode carries an exact current

For

`Q(r)=r^2 C`,

the minimum productive vorticity is

`omega_prod=-(5/3) x cross C x`.

Its Kelvin current is

`j = 5 nu C x`

`  = grad[(5nu/2) x.C.x]`.

So `j` is generally nonzero, but

`d j = 0`

and every loop period vanishes exactly.

This is an important observer correction: **the magnitude of the viscous current is not ancestry replacement**.  Only its non-exact part / curvature can change material circulation.

## Current curvature is exactly the radial productive defect

Because `div omega=0`,

`curl j = nu curl curl omega = -nu Delta omega`.

The minimum productive carrier is the toroidal `l=2` sector, and the Laplacian acts there through `L_2`.  Hence

`curl j_prod = (5nu/3) n cross (L_2 Q)n`.

Applying the Hodge transaction map to this toroidal field gives the exact tensor identity

`K_j(r) := Transaction[curl j](r) = -nu L_2 Q(r)`.

Thus `L_2 Q` is not merely a diffusion diagnostic.  Up to `-nu`, it is precisely the productive transaction tensor of `d j`, the only term that changes the pulled-back material vorticity 2-form.

## Same Hodge screen reconstructs transaction localization from ancestry-current curvature

The radial Green identity becomes

`Q(r)-r^2 C`

`= - [r^2/(5nu)] integral_0^r [1-(rho/r)^5] K_j(rho) d rho/rho`.

Therefore the departure of the transaction profile from the smooth circulation-preserving zero mode is exactly a Hodge-screened accumulation of **Kelvin ancestry-current curvature**.

The causal chain is now representation-level exact:

`transaction localization`

`<-> radial l=2 defect`

`<-> non-exact Kelvin current`

`<-> d j in the material-spacetime ancestry law`.

No auxiliary norm or preferred vortex-tube coordinate is needed.

## Connection to compulsory inward cascade

The maximum-vorticity theorem already forces positive transaction action into every shrinking physical radius if blow-up occurs.

A pure `r^2 C(t)` zero mode extending through a fixed annulus cannot be the sole divergent source: its inner and outer Hodge transactions carry the same instantaneous tensor coefficient.  If that coefficient has finite positive action, the zero mode cannot drive the blow-up; if it has infinite positive action, the fixed outer annulus would also have infinite action unless a radial defect cancels it.  Either branch forces the non-zero-mode sector to participate with unbounded cumulative action.

Since that sector is reconstructed from `K_j`, a blow-up candidate cannot remain entirely inside the circulation-preserving exact-current geometry.  It must create Kelvin ancestry-current curvature at shrinking source scales.

This still does **not** prove that the total variation of material circulation is infinite.  The curvature can move, cancel, and act on shrinking flux cells whose circulation amount tends to zero.  The missing theorem is a spacetime nonreusability statement for this current curvature.

## PROMOTE / KILL

PROMOTE:
1. The smooth productive `r^2 C` mode has nonzero but exact Kelvin current and causes no loop-circulation change.
2. `Transaction[curl j] = -nu L_2 Q` exactly.
3. `Q-r^2 C` is a Hodge-screened accumulation of Kelvin ancestry-current curvature.
4. Compulsory inward transaction cascade therefore cannot live solely in the material-circulation-preserving radial zero mode.

KILL:
1. Nonzero `j` as a proxy for ancestry replacement.
2. A regularity proof based only on total circulation amount; shrinking Clebsch flux cells can carry vanishing circulation.
3. Any claim that `K_j` already has a known finite total-variation budget.

## Next frontier — nonreusability of spacetime current curvature

The remaining question is now sharp:

**can `d j` repeatedly create the required productive radial defect on shrinking flux cells in finite time while its material-spacetime current ledger remains globally finite and telescoping?**

The next microscope should follow the flux of `j` through the extremal Clebsch cell boundaries and test whether repeated inward regeneration deposits an unavoidable spacetime boundary/current measure.
