# Source-relative ancestry flux: recruitment and viscosity are the two supply channels

## THINK — a shrinking Hodge source is not a material object

The inward-cascade theorem is Eulerian: at each time it centers a physical Hodge ball/source surface on a vorticity maximum and asks which radii supply the productive strain.  In the high source-Re branch, circulation ancestry is approximately frozen over one production time.  But a shrinking Eulerian source region cannot be a fixed material volume, because incompressibility preserves material volume.

The correct question is therefore not merely whether viscosity changes circulation.  It is:

**how does vorticity ancestry enter or leave a moving source surface that is not transported by the fluid?**

## Exact moving-surface flux law

Let `Omega` be the vorticity 2-form and

`(partial_t+L_u)Omega=-d j`,

where `j=nu curl omega . dl` is the Kelvin current one-form.

Let `Sigma(t)` be any moving oriented surface with velocity `V`.  The transport theorem gives

`d/dt integral_{Sigma(t)} Omega`

`= integral_{Sigma(t)} (partial_t+L_V)Omega`

`= integral [L_{V-u}Omega-dj]`.

Since `dOmega=0`, Cartan's formula yields

`L_{V-u}Omega=d i_{V-u}Omega`.

Hence

`d Phi_Sigma/dt`

`= integral_{partial Sigma} [ i_{V-u}Omega - j ]`.

In vector notation,

`d Phi_Sigma/dt`

`=oint_{partial Sigma} [ omega cross (V-u) - nu curl omega ].dl`.

This is an exact source-ancestry ledger.

## Two physically distinct supply channels

1. **Material recruitment / source-boundary crossing**:

   `i_{V-u}Omega`.

   It vanishes only when the source surface is material (`V=u`).  A shrinking Eulerian Hodge/Clebsch source generally cuts through material flux and can aggregate new ancestry even when viscosity is weak.

2. **Viscous ancestry replacement**:

   `j`.

   On a material surface this is the only circulation-change channel and reduces to the Kelvin law already derived.

Thus low-Re and high-Re source behavior are two terms of one exact boundary-current law, not separate mechanisms invented after the fact.

## Calibration

For `omega=B(t)e_z`, a disk of radius `R(t)`, `u=0`, and

`j=(kappa/2)(-y,x,0)`,

with `Bdot=-kappa`, the flux is

`Phi=pi R^2 B`.

The direct derivative is

`Phidot=2pi B R Rdot-pi kappa R^2`.

The two boundary currents are exactly

`recruitment=2pi B R Rdot`,

`viscous=pi kappa R^2`.

Remote Arb checks this identity over extreme amplitudes and scales.

## Consequence for inward cascade

If a high-Re source tries to preserve material circulation while its physical Hodge scale shrinks, it cannot make the source region itself material without losing the shrinking-volume geometry.  It must instead:

- recruit different material flux through the moving source boundary;
- keep only a changing portion of a larger material structure inside the source;
- or deform the material ancestry into an increasingly anisotropic/folded geometry.

The first two are measured by the recruitment term above.  The third is the frozen-flux geometry branch already exposed by pair/closure microscopes.

## PROMOTE

1. A non-material source has an exact ancestry ledger with only two boundary supply channels: relative material crossing and Kelvin viscous current.
2. High source Re does not mean ancestry inventory is fixed in an Eulerian shrinking source; material recruitment can dominate even when Kelvin circulation of each material loop is nearly frozen.
3. To follow a singular source, the natural current is source-relative, not the Kelvin current alone.

## Next frontier

Combine this ledger with the cumulative source-Re dichotomy.  If viscosity cannot supply the infinite cascade, the source must aggregate increasing material flux through its moving boundary.  The next obstruction must therefore concern ancestry aggregation/recruitment geometry rather than a scalar circulation toll.
