# Intrinsic 3D Hodge period-flux duality

## THINK — when cross-sections fail, keep the topology and metric

The Frobenius gate kills a universal reduction of a 3D vortex collar to orthogonal planar slices.  But the circulation ancestry itself is still a one-dimensional cohomology class.  The correct replacement is therefore not a more elaborate slicing rule; it is the intrinsic Hodge pairing of that class with its dual cut flux.

Let `C` be a smooth oriented three-dimensional collar with one noncontractible circulation cycle, so the relevant first cohomology is one-dimensional.  Let `gamma` be a generator of that cycle and let `Sigma` be a dual relative two-surface intersecting it once.

Let `h` be the harmonic 1-form satisfying the physical tangent/absolute boundary condition and normalized by

`integral_gamma h = 1`.

Define its circulation inductance coefficient

`I_C = integral_C |h|^2 dV`.

A circulation `Gamma` then has minimum harmonic kinetic energy

`E_h = Gamma^2 I_C/2`.

Now let `eta` be the harmonic relative 2-form normalized by

`integral_Sigma eta = 1`.

Poincare-Lefschetz duality fixes the pairing

`integral_C h wedge eta = 1`.

Because Hodge star maps the one-dimensional harmonic absolute 1-space to the dual harmonic relative 2-space,

`star h = I_C eta`.

Indeed, pairing with `h` gives

`integral h wedge star h = I_C`.

Therefore

`integral_C |eta|^2 dV = 1/I_C`,

and the intrinsic reciprocity is

`||h||_2^2 ||eta||_2^2 = 1`.

This requires no radius, no planar cross-section, and no orthogonality of the vorticity-normal plane field.

## PREDICT — planar capacity is a special shadow

For a product collar `A x [0,ell]`, where `A` is a doubly connected planar cross-section,

`I_C = ell/Cap(A)`.

Hence the previous planar law

`E_h/ell = Gamma^2/[2 Cap(A)]`

is exactly the product-geometry specialization of the 3D period-flux duality.

The quantity that survives genuine 3D deformation is `I_C`, the L2 size of the normalized harmonic circulation class, not an effective radius or even necessarily a scalar cross-sectional capacity.

## ATTACK — affine sheared flat collar

Use a flat periodic collar with coordinates

`x mod P`, `0<y<d`, `0<z<ell`,

then apply a volume-preserving affine map `F=M(x,y,z)`.

The physical circulation cycle tangent is

`a=M e_x`.

The constant vector field

`h_vec = a/[P |a|^2]`

is divergence-free, curl-free, tangent to all transformed side faces, and has unit circulation along the transformed periodic cycle.

Since `det M=1`, the physical volume remains `P d ell`, so

`I_C = d ell/[P |a|^2]`.

The dual cut at fixed `x` has oriented area vector

`M e_y cross M e_z`.

The flux of `star h` through that cut is

`d ell/[P |a|^2] = I_C`.

Thus the normalized dual two-form has norm squared `1/I_C` exactly, even when the collar is sheared and the circulation direction is not orthogonal to the coordinate cut.

For the shear

`M=[[1,0,0],[k,1,0],[0,0,1]]`,

`|a|^2=1+k^2`, so

`I_C = d ell/[P(1+k^2)]`.

Making the circulation cycle geometrically long by shear makes the fixed-period harmonic velocity smaller and the circulation occupancy cheaper, while the normalized dual cut-flux norm grows reciprocally.

This is not an energy contradiction: the dual 2-form is a geometric test field, not automatically a physical expenditure.  Its value is representation.  It tells us exactly which geometry has been exchanged when the circulation harmonic mode becomes cheap.

## AUTOPSY

### Universal planar capacity

Demoted to a special product/slice representation.  It remains exact where applicable but is not the fundamental 3D object.

### Intrinsic circulation inductance `I_C`

Promoted.  It is defined directly by the physical metric, topology, and harmonic circulation class.

### Reciprocal dual-cut norm

Promoted as an exact Hodge-Poincare identity, not as an irreversible cost.

### Shear as a free escape

Not killed by energy alone.  Shear can reduce `I_C` by lengthening the circulation cycle.  But the geometry has not disappeared: it reappears exactly in the reciprocal dual cut field.  A future spacetime theorem must determine whether that geometric exchange accelerates viscous ancestry crossing or merely stores complexity reversibly.

## PROMOTE / KILL

### PROMOTE

1. For a one-cycle 3D collar, the normalized harmonic circulation 1-form and normalized dual-cut harmonic 2-form have reciprocal squared L2 norms.
2. `E_h=Gamma^2 I_C/2` is the intrinsic 3D circulation occupancy law.
3. Planar Dirichlet capacity satisfies `I_C=ell/Cap(A)` only in the corresponding product geometry.
4. The representation survives affine shear without inventing orthogonal cross-sections.

### KILL / DO NOT PROMOTE

1. Do not equate the reciprocal dual 2-form norm with physical dissipation.
2. Do not infer a regularity contradiction from `I_C` alone; it is still an occupancy geometry.
3. Do not force helical collars into planar capacity when the Frobenius gate fails.

## Next frontier

The remaining natural bridge is now sharply stated:

- instantaneous productive strain is a signed circulation transaction;
- persistent ancestry occupies the 3D harmonic period class `I_C`;
- renewal is sideways flux of the closed material-spacetime ancestry 2-form.

The next theorem candidate should couple the **change of the intrinsic period geometry `I_C(t)`** to the **spacetime viscous flux crossing that changes material circulation**.  That coupling, if it exists, would be genuinely three-dimensional and would no longer depend on a preferred vortex-tube coordinate system.
