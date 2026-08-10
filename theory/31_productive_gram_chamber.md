# Productive Gram chamber: the intrinsic shape of a mutual-stretching pair

## THINK — quotient the observer's rigid frame before describing cycle shape

A productive pair uses three unit directions:

`a=xi_a`, `n=R/|R|`, `b=xi_b`.

Cartesian components contain a common rigid-frame redundancy.  The intrinsic shape of the directed triad is encoded by the three pairwise dot products

`alpha=a.n`,

`beta=b.n`,

`gamma=a.b`,

plus the orientation sign of the scalar triple product

`T=a.(n cross b)`.

The normalized Gram matrix is

`G=[[1,alpha,gamma],[alpha,1,beta],[gamma,beta,1]]`.

Its determinant is exactly

`det G = T^2`

`=1+2 alpha beta gamma-alpha^2-beta^2-gamma^2`.

Thus the transverse noncoplanarity gate `T` and the longitudinal gates `alpha,beta` are not arbitrary separate coordinates.  They are different geometric faces of the same intrinsic Gram shape.

The two directed transaction factors are

`K_{b->a}=T alpha`,

`K_{a->b}=-T beta`.

Their product is therefore

`P=K_{b->a}K_{a->b}`

`=(-alpha beta) det G`.

A positive mutual-stretching chamber lies on one of the two orientation sheets where

`T alpha>0`,

`-T beta>0`.

Its unsigned capacity is `P>0`.

## PREDICT — there is a universal geometric capacity

Resolve each vorticity direction into longitudinal and transverse pieces relative to `n`:

`a=a_perp+alpha n`,

`b=b_perp+beta n`.

Then

`|a_perp|=sqrt(1-alpha^2)`,

`|b_perp|=sqrt(1-beta^2)`,

and

`T=a_perp.(n cross b_perp)`.

Hence

`T^2 <= (1-alpha^2)(1-beta^2)`.

For a positive product put

`x=|alpha|`, `y=|beta|`.

Then

`P <= [x(1-x^2)] [y(1-y^2)]`.

The one-dimensional factor has the exact factorization

`2/(3 sqrt(3)) - x(1-x^2)`

`=(x-1/sqrt(3))^2 (x+2/sqrt(3)) >=0`

for `0<=x<=1`.

Therefore

`x(1-x^2)<=2/(3 sqrt(3))`

and globally

`P <= 4/27`.

Equality requires simultaneously

`|alpha|=|beta|=1/sqrt(3)`,

`alpha beta<0`,

and perpendicular transverse projections, so

`|T|=2/3`.

On the positive orientation sheet this gives

`K_{b->a}=K_{a->b}=2/(3 sqrt(3))`,

recovering the optimum found earlier from one explicit construction, now as a global geometric theorem for every unit triad.

## ATTACK — the chamber boundary has physical death modes

The formula

`P=(-alpha beta) det G`

shows three distinct routes to the boundary:

1. **longitudinal starvation**: `alpha -> 0` or `beta -> 0`;
2. **coplanarity**: `det G=T^2 ->0` while longitudinal access can remain finite;
3. **axial/transverse starvation**: `|alpha|->1` or `|beta|->1`, which forces the corresponding transverse component and hence `T` to zero.

These are not different norms becoming large or small.  They are actual geometric ways a mutual-stretching transaction loses one of the directions needed to exist.

Moreover if `P>=kappa>0`, then automatically

`T^2>=kappa`,

`|alpha|>=kappa`, `|beta|>=kappa`,

and

`1-alpha^2>=kappa`, `1-beta^2>=kappa`.

Thus a quantitatively productive pair stays inside a compact region of Gram shape space separated from every death boundary.  The estimates are deliberately crude but intrinsic.

## Productive Gram current — common spin cancels before any budget is chosen

Let the three unit directions evolve as

`adot=W a+f_a`,

`ndot=W n+f_n`,

`bdot=W b+f_b`,

where `W` is any common skew matrix and the residuals are perpendicular to their respective unit vectors.

Then common rigid spin cancels exactly from every Gram entry:

`alphadot=f_a.n+a.f_n`,

`betadot=f_b.n+b.f_n`,

`gammadot=f_a.b+a.f_b`.

The Gram determinant current is therefore

`Fdot=2[alphadot beta gamma+alpha betadot gamma+alpha beta gammadot`

`       -alpha alphadot-beta betadot-gamma gammadot]`,

where `F=det G=T^2`.

For `T!=0`, the same current is

`Fdot=2 T Tdot`,

with the signed triple-product residual current

`Tdot=f_a.(n cross b)+a.(f_n cross b)+a.(n cross f_b)`.

Finally

`Pdot=-(alphadot beta+alpha betadot)F+(-alpha beta)Fdot`.

This is the normalized shape-current version of the earlier full cycle balance.  It is formulated directly after common rigid motion has been removed.

For Navier-Stokes, a natural choice is to take `W` as the common bridge-average spin.  Then `f_n` contains bridge-average symmetric strain, while `f_a,f_b` contain endpoint strain, spin mismatch relative to the bridge, and viscous direction renewal.  The detailed dynamics remain rich, but the observer no longer pays for a common rotation of the entire physical triad.

## AUTOPSY

### `T` and longitudinal access are unrelated gates

Killed geometrically.  They are dynamically distinct channels, but they belong to one Gram shape and obey the determinant constraint.

### the optimum from module 38 was only a special symmetric construction

Killed.  `P<=4/27` is global for every unit triad, and the module-38 geometry saturates it.

### arbitrary angular motion is a shape current

Killed.  Any common skew mode drops out of all three Gram coordinates and from `T^2` and `P`.

### the Gram chamber itself proves regularity

Do not promote.  It gives the correct intrinsic state space and capacity but no finite resource preventing infinite residual shape-current action.

## PROMOTE

1. The productive pair's intrinsic normalized state is its directed Gram shape `(alpha,beta,gamma,sign T)`.
2. `T^2=det G` and `P=(-alpha beta)det G` unify the transverse and longitudinal gates without erasing their physical distinction.
3. Every positive mutual-stretching pair obeys the universal capacity `P<=4/27`.
4. Quantitative productivity confines the triad away from all Gram-degeneracy boundaries.
5. Common rigid spin is an exact gauge mode of the Gram current; only residual deformation, spin mismatch, bridge inhomogeneity, and viscosity can move the productive shape.

## Next frontier — residual shape holonomy

The regularity question is now sharper.  A singular candidate that keeps recycling productive pairs must drive their Gram shapes through repeated excursions inside the productive chamber while avoiding secular occupancy and ancestry loss.

The next physical object is therefore the **residual shape holonomy** generated after common rigid motion is removed.  The key question is whether a finite-energy material network can accumulate infinite productive Gram-current action in finite time without either:

- depositing macroscopic material shape memory;
- crossing a transaction-death boundary;
- or importing new relative geometry through bridge inhomogeneity / viscous ancestry current.

That is a much more intrinsic target than any raw norm of the velocity gradient.
