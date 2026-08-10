# Frozen vortex specific volume: self-contained branch closes, through-going ancestry survives

## THINK — raw circulation and raw material volume failed separately

The high-Re source can use circulation flux that grows while the material volume recruited at each later stage shrinks.  Neither circulation amount nor material volume alone is therefore the right frozen-ancestry state variable.

For a material vortex-line element in the Euler/frozen limit, Cauchy transport gives

`omega(X,t)=F omega0`,

and if `ds0 xi0` is initially tangent to `omega0`,

`ds(t)=|F xi0| ds0`

while

`|omega(t)|=|omega0| |F xi0|`.

Hence pointwise

`ds/|omega| = ds0/|omega0|`.

For a closed frozen vortex lineage define

`mu(C)=oint_C ds/|omega|`.

This is materially invariant.

## Physical meaning — specific material volume per circulation flux

For an infinitesimal vortex-flux tube element,

`dGamma=|omega| dA`,

and

`dV=dA ds`.

Therefore

`dV=dGamma ds/|omega|`.

After integrating along a complete flux lineage,

`dV_lineage = mu dGamma`.

So `mu` is not an artificial line norm.  It is exactly the material volume occupied per unit frozen circulation flux of that lineage.

## Smooth initial data gives a universal lower bound for closed lineages

Let `xi=omega/|omega|` along a regular closed vortex line.  Its curvature is

`kappa=|(xi.grad)xi|`.

Since

`|(xi.grad)xi| <= |grad omega|/|omega|`,

Fenchel's theorem gives

`2pi <= oint kappa ds`

`<= ||grad omega0||_infinity oint ds/|omega0|`.

Thus every closed initial vortex lineage satisfies

`mu0 >= 2pi/K0`,

where

`K0=||grad omega0||_infinity`.

In the frozen branch this lower bound persists materially.

## Exact sharp-loop cutoff

For the sharp axisymmetric transaction ray

`E=diag(2,-1,-1)`,

an extremal vortex line on a sphere of radius `r` is a latitude `n_x=c`.  Its length and vorticity magnitude are

`ell=2pi r sqrt(1-c^2)`,

`|omega_prod|=5 q |c| sqrt(1-c^2)`.

Therefore exactly

`mu=2pi r/[5q|c|]`.

At the representative strongly productive latitude `|c|=1/sqrt(2)`, frozen smooth ancestry forces

`q <= (sqrt(2)/5) K0 r`,

and hence

`q r^2/nu <= (sqrt(2)/5)(K0/nu) r^3 ->0`.

An exactly sharp closed frozen lineage cannot stay high-Re while migrating to zero radius.

## Stronger self-contained flux-volume bound

Suppose a source ball `B_L` contains complete disjoint closed frozen flux lineages.  Since every flux element has specific volume at least `mu0`,

`V(B_L) >= mu0 Phi_closed`,

where `Phi_closed` is the total absolute circulation flux carried by those complete lineages.  Thus

`Phi_closed <= (2/3) K0 L^3`.

For the localized sharp source `Phi_ext=c_p nu Re_source` with `c_p` order one, this gives

`Re_source <= C (K0/nu)L^3`.

Therefore the self-contained frozen high-Re source is impossible below the initial-data cube-root scale

`L ~ (nu/K0)^(1/3)`

up to the explicit order-one geometry coefficient.

This is the first high-Re frozen branch closed by a genuine material ancestry quantity rather than a scalar instantaneous toll.

## ATTACK — near-sharp does not force closed ancestry

The surviving escape is closure outside the Hodge source.

To test whether the exact null-enstrophy remainder controls this topology, construct on a spherical shell a divergence-free flow-box

`omega_r=A sin^2(theta) cos(theta)/r^2`,

`omega_phi=B sin(theta) cos(theta)/r`,

`omega_theta=0`.

The azimuthal part has the same axisymmetric angular pattern as the sharp productive carrier.  The radial part contributes exactly zero to `n cross omega` and is therefore transaction-null.

At `theta=pi/4`, a vortex line obeys

`dr/dphi=A/(2B)`.

To traverse a fixed shell thickness in `N` turns, choose `A/B~1/N`.  The full-shell radial-null/tangential enstrophy ratio is then proportional to

`N^-2`.

Thus a lineage can thread from one radius to another and close outside while becoming arbitrarily close in instantaneous L2 to the sharp tangent carrier, simply by using more winding.

This kills the tempting theorem

`small transaction-null remainder => closed/self-contained productive ancestry`.

The payment has moved into **material winding/line length**, which the frozen amplification conveyor naturally increases.

## PROMOTE / KILL

PROMOTE:
1. `mu=oint ds/|omega|` is a frozen material invariant and the specific material volume per circulation flux of a closed vortex lineage.
2. Smooth initial data gives `mu>=2pi/||grad omega0||_infinity` for every closed vortex line.
3. Exact sharp self-contained frozen lineages become low-Re at small radius; complete self-contained flux inside `B_L` is bounded by `(2/3)K0L^3`.
4. Therefore a high-Re inward source must eventually use viscous ancestry renewal, lineage switching, or through-going closure outside the source.

KILL:
1. Instantaneous small transaction-null enstrophy as a topological guarantee of closed ancestry.
2. A self-contained frozen-loop cascade all the way to zero scale.

## Next frontier — winding/tether ancestry

The only frozen high-Re escape left by this microscope is increasingly long, through-going material ancestry whose closure/specific volume lies outside the shrinking Hodge source.  Such lineages can make their radial transaction-null component arbitrarily small by increasing winding.

The next question is therefore exact and geometric:

**can the required winding/line-length growth be packed through a shrinking source while each flux tube remains above its own viscous persistence scale and the global finite-energy/dissipation ledgers stay finite?**

This is a tether-packing problem, not a circulation-stock problem.
