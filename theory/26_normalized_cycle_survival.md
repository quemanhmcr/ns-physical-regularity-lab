# Normalized cycle survival: a positive cycle can remain positive and still become sterile

## THINK — pair-cell renewal is not yet cycle survival

The pair ancestry cell

`D_ab = omega_a . (R cross omega_b) = rho_a rho_b r T_ab`

separates common affine deformation from true bridge/viscous renewal.  But a positive interaction cycle is controlled by the **normalized** three-dimensional fraction

`T_ab = D_ab/(rho_a rho_b r)`,

plus the two longitudinal projections.  Large `|D_ab|` alone does not mean the cycle remains productively three-dimensional: the endpoint vorticity magnitudes and the bridge can grow even faster.

Whenever all factors are nonzero, the exact logarithmic survival law is

`d log|T|/dt = lambda_D - sigma_a - sigma_b - sigma_R`,

where

`lambda_D = Ddot/D`,

`sigma_a = d log rho_a/dt`,

`sigma_b = d log rho_b/dt`,

`sigma_R = d log r/dt`.

Define the **cycle renewal deficit**

`Delta_cycle = sigma_a + sigma_b + sigma_R - lambda_D`.

Then exactly

`Delta_cycle = - d log|T|/dt`.

So a persistent noncoplanar cycle is not maintained merely by renewing `D`.  Renewal must keep pace with the combined expansion of both endpoint vorticity magnitudes **and** the material bridge.

If `Delta_cycle` stays positive for many strain times, the genuinely-three-dimensional angular cell is consumed exponentially even if `D` itself grows rapidly.

## PREDICT — an exact quadratic NS family can expose the deficit

Take constants `g>0`, `eps>0` and the steady divergence-free velocity

`u=(-g x + eps y z, g y, 0)`.

It is an exact Euler and Navier-Stokes solution on `R^3` because

`Delta u=0`

and

`(u.grad)u=(g^2 x,g^2 y,0)=-grad p`

for

`p=-(g^2/2)(x^2+y^2)`.

Its vorticity is

`omega=(0,eps y,-eps z)`.

This field has infinite global energy and is used only as an exact local-mechanism calibration, never as a finite-energy counterexample or regularity model.

Follow the same two material particles as in the instantaneous microscope:

`X_a(0)=(0,1,0)`,

`X_b(0)=(-L,-1,1)`, `L>0`.

Use strain time

`N=g t`

and the dimensionless bridge-renewal parameter

`c=eps/(2g)>0`.

Set

`x=exp(-N)`.

The exact trajectories give the scaled material bridge

`x R = (-H,-2,x)`,

where

`H=c+(L-c)x^2`.

Hence

`r = sqrt(H^2+4+x^2)/x`.

The endpoint vorticity directions are

`xi_a=(0,1,0)`,

`xi_b=(0,-1,-x)/sqrt(1+x^2)`.

The normalized three-dimensional pair cell is therefore

`T = - x H/[sqrt(1+x^2) sqrt(H^2+4+x^2)]`.

The longitudinal projections are

`alpha=xi_a.n=-2/sqrt(H^2+4+x^2)`,

`beta=xi_b.n=(2-x^2)/[sqrt(1+x^2) sqrt(H^2+4+x^2)]`.

For every finite `N>=0`, `H>0`, `alpha<0`, `beta>0`, `T<0`.  Thus both directed angular transactions remain strictly positive:

`K_{b->a}=T alpha`

`= 2 x H/[sqrt(1+x^2)(H^2+4+x^2)] > 0`,

and

`K_{a->b}=-T beta`

`= x H(2-x^2)/[(1+x^2)(H^2+4+x^2)] > 0`.

The cycle never needs to flip sign.

## ATTACK — exact full-time renewal deficit

The exact pair-cell magnitude is

`|D| = eps^2 [L+c(exp(2N)-1)]`

`    = eps^2 H/x^2`.

Therefore its renewal rate per strain time is

`lambda_D/g = d log|D|/dN = 2c/H`.

The endpoint vorticity gain rates are

`sigma_a/g = 1`,

`sigma_b/g = 1/(1+x^2)`.

The bridge rate is

`sigma_R/g`

`= 1 - x^2[2H(L-c)+1]/(H^2+4+x^2)`.

Consequently

`d log|T|/dN`

`= 2c/H - 1 - 1/(1+x^2) - sigma_R/g`.

This agrees identically with direct differentiation of the closed form for `T`.

As `N -> infinity`, `x -> 0`, `H -> c`, and

`lambda_D/g -> 2`,

`sigma_a/g -> 1`,

`sigma_b/g -> 1`,

`sigma_R/g -> 1`.

Thus

`Delta_cycle/g -> 1`,

and

`d log|T|/dN -> -1`.

So even though the physical ancestry cell `|D|` is being renewed at asymptotic rate `2g`, the two endpoint magnitudes consume `2g` and the expanding material bridge consumes one additional `g`.  The normalized three-dimensional transaction geometry loses by exactly one strain exponent.

The longitudinal projections do **not** collapse:

`alpha -> -2/sqrt(c^2+4)`,

`beta -> 2/sqrt(c^2+4)`.

Therefore the transaction edges die specifically because `T -> 0`, not because the positive-cycle sign condition fails.

Moreover

`K_{b->a} ~ K_{a->b}`

`~ [2c/(c^2+4)] exp(-N)`.

The positive two-cycle survives topologically in sign but becomes exponentially sterile in productive angular strength.

## AUTOPSY

### Pair-cell renewal alone as sufficient persistence

Killed.  `|D|` can grow exponentially and still fail to preserve the normalized productive cell because `rho_a rho_b r` grows faster.

### Positive edge signs as sufficient persistence

Killed.  Both edges remain positive for every finite strain time in the exact family while their magnitudes vanish exponentially.

### Near-contact as the necessary cycle-death mechanism

Killed again, now dynamically rather than instantaneously.  In this exact family the bridge grows like `exp(N)` asymptotically; the cycle becomes sterile by coplanarization of the normalized triad.

### Universal exponential death of every positive cycle

Do not promote.  This family has infinite global energy and very special polynomial geometry.  It is a canonical exact mechanism, not a universal theorem.

## PROMOTE

1. The normalized cycle-survival balance

   `d log|T|/dt = lambda_D-sigma_a-sigma_b-sigma_R`

   is exact whenever the factors are nonzero.
2. `Delta_cycle=sigma_a+sigma_b+sigma_R-lambda_D` is the exact **renewal deficit** of the genuinely-three-dimensional pair geometry.
3. A positive cycle can keep both directed edge signs forever and nevertheless lose all productive strength through `T -> 0`.
4. Pair-cell renewal must pay not only for endpoint vorticity amplification but also for material-bridge deformation if angular productivity is to persist.
5. In the exact quadratic NS calibration, `Delta_cycle/g -> 1`, so both edge weights decay as `exp(-N)` with explicit coefficient `2c/(c^2+4)`.

## Next frontier — renewal hierarchy, now normalized

The surviving escape is sharper than before.  A singular candidate must repeatedly arrange

`lambda_D >= sigma_a+sigma_b+sigma_R`

on enough productive pairs to prevent angular-cell depletion, while simultaneously keeping the longitudinal projections and Hodge transaction signs productive.

The next physical question is therefore not merely whether `J_bridge` can be large.  It is whether **bridge-current renewal can keep pace with the entire normalized survival demand** without creating a new, smaller non-affine structure.

That suggests a hierarchy built from the actual bridge current density, not a Hessian norm:

- level 0: endpoint vorticity amplification;
- level 1: pair ancestry-cell renewal `D`;
- level 2: normalized survival deficit of `T`;
- level 3: material lineage of the bridge-gradient current that feeds `D`.

The key future attack is whether a finite coherent geometry can circulate this normalized renewal current indefinitely, or whether sustained productivity forces a cascade of new bridge-scale structure.
