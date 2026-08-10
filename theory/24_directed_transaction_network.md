# Directed transaction network: mutual stretching is not a conservative edge flow

## THINK — after self-toll fails, inspect the edge law itself

The closed-ring attack killed the rule

`donor creates productive strain => donor must deform its own period geometry`.

The natural replacement is a directed transaction graph.  Before using such a graph for any global bookkeeping, ask a more basic physical question:

**are the strain transactions on opposite directions of an interaction antisymmetric, like a conserved transfer, or can two ancestries stretch each other simultaneously?**

The Biot-Savart kernel answers this directly.

Let two small coherent vorticity elements have unit directions `a` and `b`, separated by unit vector `n` from the target `a` element to the donor `b` element.  Up to the common positive Biot-Savart strength/distance prefactor, the angular stretching transaction from `b` to `a` is

`K_{b->a} = (a.n) a.(n cross b)`.

When the roles are reversed, the separation direction reverses, giving

`K_{a->b} = (b.(-n)) b.((-n) cross a)`

`          = -(b.n) a.(n cross b)`.

Define the genuine-3D triple product

`T = a.(n cross b)`.

Then

`K_{b->a}=T(a.n)`,

`K_{a->b}=-T(b.n)`.

Nothing here is antisymmetric in general.

## PREDICT — one-way transactions and positive two-cycles both exist

### One-way edge

Take

`n=e_z`,

`a=(1/sqrt(2),0,1/sqrt(2))`,

`b=(0,-1,0)`.

Then

`T=1/sqrt(2)`,

`a.n=1/sqrt(2)`,

`b.n=0`,

so

`K_{b->a}=1/2`,

`K_{a->b}=0`.

The network is genuinely directed.

### Symmetric mutual-stretching pair

Let

`a=(sin theta,0,cos theta)`,

`b=(0,-sin theta,-cos theta)`,

`n=e_z`.

Then

`T=sin^2 theta`,

`a.n=cos theta`,

`b.n=-cos theta`,

and therefore

`K_{b->a}=K_{a->b}=sin^2 theta cos theta`.

For `0<theta<pi/2`, both are strictly positive.

So two vorticity elements can stretch each other at the same instant.  Vortex stretching is not a conservative scalar transfer between nodes.

The symmetric efficiency

`f(c)=c(1-c^2)`, `c=cos theta`,

is maximized at

`c=1/sqrt(3)`,

with

`f_max=2/(3 sqrt(3))`.

The corresponding angle is the familiar order-one three-dimensional tilt, not an asymptotically singular alignment.

## ATTACK — what geometry is required for a positive two-cycle?

The product of the two directed angular transactions is

`K_{b->a} K_{a->b} = -T^2 (a.n)(b.n)`.

Therefore both directions can have the same nonzero sign only if

`(a.n)(b.n)<0`.

A positive two-cycle requires the two vorticity directions to have opposite longitudinal projections along their separation axis.

It also requires

`T=a.(n cross b) != 0`.

Thus coplanar geometry cannot produce the mutual-stretching cycle.  The cycle is intrinsically three-dimensional.

This is a useful geometric restriction, but it is not yet a cost.  Straight coherent segments can realize the local direction geometry while putting their bends/closures elsewhere, as earlier modules already showed.

## AUTOPSY — transaction graph is not an energy-flow graph

The graph edges are signed *rates of stretching action*, not transfers of a conserved edge currency.

A positive edge `i->j` does not imply a compensating negative edge `j->i`.  In fact:

- one-way edges exist;
- positive two-cycles exist;
- negative two-cycles also exist by reversing the transverse orientation;
- coplanar configurations kill the edge through the triple product.

Therefore an amplification proof cannot sum edge transactions and hope for pairwise cancellation.

The finite resource remains kinetic energy/circulation geometry of the whole flow, while the transaction graph tells us how that resource is *used to amplify vorticity*, not how it is conserved.

## PROMOTE / KILL

### PROMOTE

1. Productive interactions are naturally directed edges.
2. The angular edge kernel is controlled by the triple product `T=a.(n cross b)` and longitudinal target/donor projections.
3. Positive mutual-stretching two-cycles are allowed by the exact Biot-Savart angular kernel.
4. Mutual stretching is genuinely 3D and requires opposite longitudinal tilts along the pair separation axis.

### KILL / DO NOT PROMOTE

1. Kill any antisymmetric edge-transfer ledger for vortex stretching.
2. Kill the assumption that an amplification graph must be acyclic.
3. Do not infer a cost merely from the order-one angular tilt; straight coherent pieces can displace closure geometry.
4. No regularity contradiction follows from pairwise angular structure alone.

## Next frontier — persistence of a positive transaction cycle

A singular amplification route can now try to reuse a finite set of mutually stretching ancestries rather than recruit endlessly new donors.

The next physical question is dynamic:

**can a positive transaction cycle remain geometrically productive for arbitrarily many strain e-folds while all participating circulation ancestries remain isolated and finite-energy?**

If not, the failure mode must be one of the structures already discovered:

- the pair separations/orientations deform and productive edges collapse or change sign;
- a participating period geometry becomes expensive or thin;
- circulation ancestry renews through viscosity;
- closure/packing of the mutually interacting lineages consumes nonlocal space/energy.

This is the next place to search for a nonreusable cost: not on one edge and not on one node, but on the **lifetime of a productive directed cycle**.
