# Productive cycle lifetime: an exact material pair-cell balance

## THINK — stop assigning a lifetime to an edge in isolation

A positive two-cycle is not a conservative edge transfer.  The next object should therefore not be a graph-theoretic cost imposed from outside.  Follow the physical geometry already carried by two material vorticity ancestries.

Let `X_a(t), X_b(t)` be two material particles while the Navier-Stokes solution is smooth.  Set

`R = X_b-X_a`, `r=|R|`, `n=R/r`,

`omega_a=omega(X_a,t)`, `omega_b=omega(X_b,t)`,

`rho_a=|omega_a|`, `rho_b=|omega_b|`,

`xi_a=omega_a/rho_a`, `xi_b=omega_b/rho_b`.

The directed angular transaction gate from the previous turn contains the triple product

`T_ab = xi_a . (n cross xi_b)`.

Do not normalize it prematurely.  The physical unnormalized object is

`D_ab = omega_a . (R cross omega_b)`.

Exactly,

`D_ab = rho_a rho_b r T_ab`.

`D_ab` is translation invariant, rotation invariant, and symmetric under exchanging the labels `a,b` together with reversing the bridge.  It vanishes exactly when the two endpoint vorticities and their material bridge become coplanar.  It is therefore the oriented three-dimensional **ancestry cell** of the pair.

## Exact NS balance — common incompressible deformation cancels

Write

`A(x,t)=grad u(x,t)`.

The exact material bridge velocity is

`Rdot = u(X_b)-u(X_a) = Abar_ab R`,

where the physically distinguished bridge average is

`Abar_ab = integral_0^1 A(X_a+sR,t) ds`.

Because `div u=0`,

`tr Abar_ab=0`.

At the endpoints the vorticity equation is

`omega_dot_a = A_a omega_a + nu Delta omega_a`,

`omega_dot_b = A_b omega_b + nu Delta omega_b`.

Differentiate the determinant `D_ab=det[omega_a,R,omega_b]`.  Add and subtract the common bridge matrix `Abar_ab` from the two endpoint gradients.  The three common-affine terms obey the determinant identity

`det[Abar p,q,r] + det[p,Abar q,r] + det[p,q,Abar r]`

`= tr(Abar) det[p,q,r] = 0`.

Therefore the exact material pair-cell law is

`Ddot_ab = J_bridge + J_nu`,

with

`J_bridge = det[(A_a-Abar) omega_a, R, omega_b]`

`         + det[omega_a, R, (A_b-Abar) omega_b]`,

and

`J_nu = nu det[Delta omega_a,R,omega_b]`

`     + nu det[omega_a,R,Delta omega_b]`.

This is the first exact cycle-lifetime balance in the program with no preferred tube radius, no planar slicing, and no artificial graph conservation law.

A common incompressible affine deformation may strongly stretch both endpoint vorticity vectors and the bridge, but it cannot change `D_ab`.  To change the oriented ancestry cell, the flow must use spatial variation of its velocity gradient across the bridge or viscosity.

### Bridge-current representation

The endpoint mismatch is itself a line moment of the actual Hessian along the physical bridge.  With

`G(s) = (R.grad) A(X_a+sR)`,

one has exactly

`A_a-Abar = - integral_0^1 (1-s) G(s) ds`,

`A_b-Abar =   integral_0^1 s G(s) ds`.

Thus `J_bridge` is not an abstract norm of `grad^2 u`; it is a signed bridge moment selected by the two actual vorticity ancestries and their connecting material separation.

## PREDICT — exact survival trilemma

Whenever `D_ab`, `T_ab` are nonzero, define endpoint amplification factors

`Lambda_a = rho_a(t)/rho_a(0)`,

`Lambda_b = rho_b(t)/rho_b(0)`,

and pair-cell renewal factor

`R_D = |D_ab(t)|/|D_ab(0)|`.

The identity `D=rho_a rho_b r T` gives the exact factorization

`Lambda_a Lambda_b = R_D [r(0)/r(t)] [|T(0)|/|T(t)|]`.

This exposes three and only three multiplicative ways for large pair amplification to occur:

1. **bridge compression:** `r(t)` becomes small;
2. **angular degeneration:** `|T(t)|` becomes small, which kills the genuinely-3D transaction gate rather than sustaining the cycle;
3. **ancestry-cell renewal:** `|D_ab|` grows, and the exact balance says this requires bridge-gradient inhomogeneity and/or viscosity.

For a genuinely persistent positive two-cycle suppose both directed angular weights stay above `kappa>0`.  Since each edge is `T` times a longitudinal projection of magnitude at most one, this implies `|T|>=kappa`.  If the two ancestries also remain separated by `r(t)>=delta r(0)`, then

`Lambda_a Lambda_b <= R_D |T(0)|/(delta kappa)`.

So arbitrarily large amplification in a uniformly productive, non-contacting two-cycle requires arbitrarily large pair-cell renewal.  This is an exact geometric statement.  It is **not yet a finite-resource contradiction**, because `J_bridge` can be reversible and can in principle regenerate `D_ab`.

## ATTACK 1 — divergence-free quadratic field with an actual positive pair cycle

Use the smooth divergence-free initial velocity

`u=(a x + eps y z, -a y, 0)`.

Then

`A=[[a,eps z,eps y],[0,-a,0],[0,0,0]]`,

`omega=(0,eps y,-eps z)`,

and `Delta omega=0`.

Choose material particles at the instant

`X_a=(0,1,0)`,

`X_b=(-L,-1,1)`, `L>0`.

Then

`R=(-L,-2,1)`,

`omega_a=(0,eps,0)`,

`omega_b=(0,-eps,-eps)`.

The pair has a positive angular two-cycle:

`K_{b->a}=sqrt(2) L/(L^2+5)>0`,

`K_{a->b}=L/[2(L^2+5)]>0`.

Its oriented ancestry cell is

`D_ab=-L eps^2`,

while the exact balance gives

`Ddot_ab=-eps^3`,

hence

`Ddot_ab/D_ab=eps/L`.

A numerical observer must respect this cancellation structurally: measuring the tiny renewal by subtracting separate `O(|a| eps^2)` determinant terms is ill-conditioned when `|a|/eps` is huge.  The physical bridge-average decomposition removes the common mode before observation; numerically, the projected mismatch should be represented directly rather than recovered by subtracting two nearly equal parent gradients.

The parameter `a`, which can make both endpoint vorticity magnitudes amplify arbitrarily rapidly when `a<0`, cancels completely from the pair-cell renewal rate:

`d log rho_a/dt=-a`,

`d log rho_b/dt=-a/2`,

so

`d log(rho_a rho_b)/dt=-3a/2`.

Therefore

`d log(r |T|)/dt = eps/L + 3a/2`.

For `a=-g<0` with `g >> eps/L`, the two endpoint vorticities amplify while the productive pair geometry `r|T|` must collapse at rate approximately `3g/2` unless the non-affine `eps` bridge current grows comparably.

At `L=2`, the common affine part contributes exactly zero to the instantaneous separation rate, so the collapse occurs almost entirely by loss of the triple-product angle rather than by bridge approach.  This is an exact **cycle-death by coplanarization** calibration.

For other bridge orientations the same product collapse can instead be shared with physical approach.  The invariant object is `r|T|`, not either factor separately.

## ATTACK 2 — exact Navier-Stokes heat flow isolates viscous pair-cell renewal

Take the exact nonlinear-null solution

`u(y,t)=(A exp(-nu k^2 t) cos(k y), 0, B exp(-nu m^2 t) cos(m y))`.

It is divergence free and `(u.grad)u=0`, so it solves Navier-Stokes with constant pressure.  Its vorticity is

`omega=(-B m exp(-nu m^2 t) sin(m y), 0, A k exp(-nu k^2 t) sin(k y))`.

Material particles keep their `y` coordinates.  Choose two `y` values for which the endpoint vorticities are not parallel.  Because every vorticity has zero `y` component while every velocity gradient has only a `y` column,

`A omega = 0`

at every point.  Hence `J_bridge=0` exactly and the pair cell changes only by viscosity.

The determinant is bilinear in one `m`-mode component and one `k`-mode component, so

`D_ab(t)=D_ab(0) exp[-nu(k^2+m^2)t]`.

Therefore

`Ddot_ab/D_ab=-nu(k^2+m^2)`.

This is a clean null control showing that the viscous term in the pair-cell law is a real ancestry-cell current, not a bookkeeping artifact of the bridge decomposition.

## AUTOPSY

### A universal finite lifetime for every positive cycle

Not proved.  The trilemma leaves a live escape: `D_ab` can be renewed by non-affine bridge dynamics.  A cycle may therefore recruit spatial gradient variation to avoid both contact and angular death.

### Bridge compression alone as the universal obstruction

Killed.  The quadratic attack at `L=2` shows large endpoint amplification can kill the cycle almost entirely through `T -> 0` while the separation has no common-affine compression at that instant.

### Angular tilt alone as a reusable free resource

Demoted.  `T` is part of the exact ancestry cell `D=rho_a rho_b r T`; if endpoint magnitudes grow without cell renewal, angular noncoplanarity and/or separation must be consumed.

### Pair-cell renewal as irreversible cost

Do not promote.  `J_bridge` is geometric and can be reversible.  `J_nu` is viscous, but even its sign need not universally reduce `|D|` for arbitrary data.

## PROMOTE / KILL

### PROMOTE

1. `D_ab=omega_a.(R cross omega_b)` is the natural unnormalized three-dimensional pair ancestry cell behind the transaction triple product.
2. Its exact NS material balance removes every common incompressible affine deformation and retains only bridge-gradient inhomogeneity plus viscosity.
3. The exact survival factorization is

   `Lambda_a Lambda_b = R_D (r0/r) (|T0|/|T|)`.

4. Hence sustained productive amplification has a physical trilemma: near-contact, angular cycle death, or pair-cell renewal.
5. The nonlinear bridge current is a signed line moment of the actual velocity-gradient variation along the bridge, not an imposed global norm.

### KILL / DO NOT PROMOTE

1. Do not claim positive cycles have a universal finite lifetime yet.
2. Do not treat `D_ab` itself as a conserved or finite global resource.
3. Do not call `J_bridge` dissipation; it can be reversible.
4. Do not force the trilemma into a scalar Gronwall estimate.
5. Do not identify two material sample points with complete vortex tubes; the pair cell is a microscope of ancestry geometry, not a tube decomposition.

## Next frontier — can renewal be reused indefinitely?

The problem is now sharper.  A positive cycle can escape compression/coplanarization only by repeatedly renewing its oriented ancestry cells through `J_bridge` or viscosity.

The next microscope should ask whether **bridge renewal itself can be recycled without generating new physical structure**.  Because `J_bridge` is a weighted moment of `(R.grad)grad u` along the material bridge, sustaining it while the endpoints remain strongly productive may require a hierarchy of non-affine bridge structure.  The natural next attack is therefore not a Sobolev norm of the Hessian, but the lineage of the bridge-current itself: does renewal of a productive pair cell force new smaller bridge structure, or can a finite coherent geometry circulate that current indefinitely?
