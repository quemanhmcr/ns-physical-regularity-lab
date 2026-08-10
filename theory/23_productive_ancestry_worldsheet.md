# Productive ancestry worldsheet: metric deformation, renewal, and the transaction-network obstruction

## THINK — let the material collar carry its own geometry

The previous frontier left three intrinsic objects:

- productive strain transaction `Gamma_Q`;
- intrinsic circulation inductance `I_C` of a material collar;
- the closed material-spacetime vorticity flux form that records ancestry renewal.

The next temptation is to force them into one scalar inequality.  Do not do that.  First ask how `I_C` itself changes under the physical deformation of the collar.

Let `C(t)=X_t(C0)` be a material collar with one circulation cycle.  Let `h(t)` be the physical harmonic 1-form on `C(t)` satisfying the tangent/absolute boundary condition and normalized by

`integral_gamma(t) h(t) = 1`.

Define

`I_C(t)=integral_C(t) |h(t)|^2 dV`.

Pull `h` back to `C0`.  The pulled-back form stays in the same cohomology class, while the pulled-back metric evolves with

`g_dot = 2 X_t^* S`,

because the antisymmetric part of `grad u` is a local rotation and incompressibility fixes the material volume form.

The harmonic representative minimizes the metric energy inside its fixed cohomology class.  Therefore its own first variation, which has zero periods and is exact in the one-cycle case, is orthogonal to `h`.  Only the metric variation survives.  This gives the exact material shape law

`I_C_dot = -2 integral_C h . S h dV`.

Define the intrinsic period-geometry strain rate

`sigma_C = [integral h.S.h dV]/I_C`.

Then

`d log I_C/dt = -2 sigma_C`.

This is not an estimate and not a norm choice.  It is the shape derivative of the physical harmonic circulation class.

For several independent circulation cycles, with unit-period harmonic basis `h_a` and inductance matrix

`I_ab = integral h_a.h_b dV`,

the same argument gives the matrix law

`I_ab_dot = -2 integral h_a.S.h_b dV`.

## PREDICT — compression of the circulation direction makes persistence expensive

If the circulation cycle is stretched by the material deformation, a unit period can be carried by a smaller velocity and `I_C` can fall.  If the circulation cycle is compressed, the unit-period harmonic field must intensify and `I_C` rises.

This immediately reconciles two earlier calibrations that looked different:

- in the affine shear family, lengthening the circulation cycle makes `I_C` fall like `1/(1+k^2)`;
- in the frozen-flux axial stretching family, the azimuthal circulation cycle is transversely compressed, so `I_C` grows like the axial amplification factor.

The sign is not arbitrary.  It is exactly the sign of the strain seen by the normalized harmonic period field.

## ATTACK I — affine material deformations

For a volume-preserving diagonal deformation

`M=diag(L,L^(-1/2),L^(-1/2))`,

with instantaneous strain rate

`S=diag(a,-a/2,-a/2)`,

a circulation cycle along the extensional `x` direction has

`I_x/I_x0 = L^(-2)`,

`sigma_x=a`,

and hence `I_x_dot/I_x=-2a`.

A cycle along either compressed transverse direction has

`I_perp/I_perp0=L`,

`sigma_perp=-a/2`,

and `I_perp_dot/I_perp=a`.

For the shear

`M=[[1,0,0],[k,1,0],[0,0,1]]`,

with `k_dot=g`, the one-cycle flat-collar law is

`I_C=I0/(1+k^2)`,

`sigma_C=g k/(1+k^2)`,

so again

`I_C_dot=-2 sigma_C I_C`.

Rigid rotation has `S=0`, hence leaves `I_C` unchanged even though the collar moves in space.

## ATTACK II — do not confuse harmonic ancestry with every material loop

Once vorticity penetrates a collar, the circulation of different homologous material loops need not agree.  Then a single scalar harmonic coefficient and the individual loop circulations are different physical observables.

An exact Navier-Stokes calibration makes this explicit.  On a periodic box take

`u=(H + A exp(-nu k^2 t) cos(k y), 0, 0)`.

The nonlinear term vanishes identically, so this is an exact viscous shear solution.  The spatially constant `H` is the harmonic circulation mode.  The cosine is a vortical mode.

Viscosity damps the cosine exponentially but cannot change `H` on the closed periodic topology.  The kinetic-energy split is exactly orthogonal:

`E = E_harm + E_vort`,

`E_harm = constant`,

`E_vort ~ exp(-2 nu k^2 t)`.

However a material `x`-cycle at fixed `y` has

`Gamma_y(t)=P[H + A exp(-nu k^2 t) cos(k y)]`,

and its circulation changes by the exact Kelvin-viscous term.

Thus:

**the finite-thickness harmonic ancestry mode can remain unchanged while individual codimension-two material-loop circulations renew.**

This is not a contradiction.  It is the observer hierarchy that the single-loop counterexample was already warning about.

## ATTACK III — a steady Eulerian vortex can be a material-ancestry conveyor belt

Use the exact steady Burgers vortex.  Its radial material motion is

`R_dot=-(a/2)R`,

so for a material circle

`q(t)=a R(t)^2/(4 nu)=q0 exp(-a t)`.

The circulation carried by that material circle is

`Gamma_m(t)=Gamma_inf [1-exp(-q(t))]`.

Although the Eulerian vortex profile is exactly steady,

`Gamma_m_dot = -a Gamma_inf q exp(-q)`.

The same quantity is obtained from the Kelvin viscous current around the moving circle.

Define the fractional material-renewal rate

`lambda_nu = -Gamma_m_dot/Gamma_m`

`          = a q exp(-q)/(1-exp(-q))`.

The transverse material-compression rate is `a/2`.  Their ratio is

`lambda_nu/(a/2)=2 q exp(-q)/(1-exp(-q))`.

There is a unique order-one crossover solving

`exp(q)-1=2q`,

numerically near `q=1.256...`, or ancestry clock

`Chi=4q about 5.03`.

For larger `q`, transverse compression outruns fractional circulation renewal.  Once contraction carries the material circle to order-one `q`, renewal catches and then exceeds that geometric rate.

This is a stronger interpretation of Burgers balance: the steady Eulerian vortex is not one material vortex persisting forever.  It is a **steady structure maintained while material circulation ancestry is continuously replaced by viscosity**.

Do not promote the numerical crossover as a universal constant.  Promote the mechanism: Eulerian persistence can hide material renewal.

## ATTACK IV — productive strain does not force the donor's own `I_C` to change

The most dangerous proposed bridge would be

`productive transaction by donor C  =>  |d log I_C/dt| is large`.

It is false even at the canonical Biot-Savart geometry level.

A circular vortex filament of radius `R` and circulation `Gamma` induces on its axis

`u_z(z)=Gamma R^2/[2(R^2+z^2)^(3/2)]`,

so

`S_zz(z)=partial_z u_z`

`       =-3 Gamma R^2 z/[2(R^2+z^2)^(5/2)]`.

On the `z<0` side this is positive stretching of the axial target direction and is `O(Gamma/R^2)` when `|z|~R`.

But rotational symmetry forces the self-induced velocity of the circular filament, with isotropic core regularization, to be the same axial translation at every point of the ring.  The donor ring therefore translates without changing its radius or its circulation-cycle geometry.  Its own instantaneous period-geometry rate is zero while its outgoing productive strain is nonzero.

This kills every universal instantaneous law that charges a donor's outgoing strain to deformation of that same donor's `I_C`.

The earlier straight-segment example already showed the local version; the closed ring removes the objection that an open segment is not globally admissible.

## AUTOPSY — the natural object is now a directed transaction network

A single lineage is not a closed accounting system.

Each persistent ancestry `C_i` has a node geometry:

- circulation/period data;
- inductance `I_i`;
- intrinsic geometry rate `sigma_i`.

But productive strain is carried on directed edges:

`T_{i -> j}` = signed strain transaction supplied by ancestry `i` to target lineage or target direction `j`.

The closed-ring attack shows that a large outgoing edge need not create a large self-node deformation rate.  The target can be stretched while the donor remains geometrically steady.

Therefore the next regularity mechanism, if it exists, cannot be a sum of independent one-lineage tolls.  It must exploit a constraint on the **whole directed interaction network**: closure, finite spatial packing, mutual Biot-Savart geometry, or the need for some lineage in an amplification chain to be deformed/renewed.

This is physically closer to vortex dynamics.  Vorticity does not stretch itself in isolation; stretching is an interaction.

## PROMOTE / KILL

### PROMOTE

1. Exact material Hodge shape law:
   `I_C_dot = -2 integral h.S.h`.
2. Intrinsic period-geometry strain rate:
   `sigma_C=(integral h.S.h)/I_C` and `d log I_C/dt=-2 sigma_C`.
3. The multi-cycle inductance matrix obeys the analogous bilinear shape derivative.
4. Harmonic ancestry and individual material-loop circulation must be distinguished after vorticity penetrates a collar.
5. Exact Burgers flow demonstrates steady Eulerian coherence with continuous material-circulation renewal.
6. Productive strain is an inter-ancestry transaction, not necessarily a self-deformation of the donor.

### KILL / DO NOT PROMOTE

1. Kill any universal implication `outgoing productive q => large |I_C_dot|` for the same donor.
2. Do not identify a harmonic collar coefficient with every material-loop circulation once the collar contains vorticity.
3. Do not treat the Burgers crossover `Chi about 5.03` as universal.
4. Do not add node occupancy energies over time as irreversible costs.
5. No regularity contradiction has yet been obtained.

## Next frontier — transaction network / amplification ancestry graph

The next microscope should stop asking whether one donor pays for its own outgoing strain.  Instead build the directed graph of productive transactions among coherent ancestry carriers over a material time window.

A blow-up lineage requires an unbounded accumulated positive incoming stretching history.  The question is whether an arbitrarily long directed amplification chain can avoid, at every stage, all three events:

- substantial period-geometry deformation of the participating ancestries;
- viscous ancestry renewal across a material collar;
- spatial/energetic crowding required to keep supplying productive `Gamma/r^2` transactions.

The physical target is a network alternative, not a scalar Gronwall estimate: **either an amplification path closes into a mutually stretching cycle, escapes to ever-new donors, or repeatedly renews ancestry.**  Each branch should then be attacked using the geometry already discovered rather than by imposing a preselected functional norm.
