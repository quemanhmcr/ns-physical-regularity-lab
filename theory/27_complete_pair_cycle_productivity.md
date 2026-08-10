# Complete pair-cycle productivity: transverse cell and longitudinal access are distinct physical gates

## THINK — normalized noncoplanarity is only one gate

The directed Biot-Savart angular transactions for a target ancestry direction `xi_a`, donor direction `xi_b`, and separation direction `n` are

`K_{b->a}=T alpha`,

`K_{a->b}=-T beta`,

where

`T=xi_a.(n cross xi_b)`,

`alpha=xi_a.n`,

`beta=xi_b.n`.

The previous microscope followed `T`, the normalized three-dimensional ancestry cell.  But a positive two-cycle needs two independent pieces of geometry:

1. noncoplanarity: `T != 0`;
2. opposite longitudinal access: `alpha beta < 0`.

A cycle can therefore remain strongly noncoplanar and still become sterile if either vorticity direction loses longitudinal access to the pair separation.

The exact symmetric object already supplied by the transaction law is the product of the two directed edges:

`P_ab = K_{b->a} K_{a->b}`

`     = -T^2 alpha beta`.

For a positive two-cycle, `P_ab>0`.  This is not an imposed graph norm: it is the exact product of the two physical angular stretching transactions.

Define the **longitudinal gate**

`G_ab=-alpha beta`.

Then

`P_ab=T^2 G_ab`.

Thus full mutual productivity factorizes into a transverse genuinely-3D gate and a longitudinal-access gate.

## PREDICT — three unnormalized ancestry observables expose the full gate

For the material pair, write

`D = omega_a.(R cross omega_b)`,

`L_a = omega_a.R`,

`L_b = omega_b.R`,

with

`rho_a=|omega_a|`, `rho_b=|omega_b|`, `r=|R|`.

Then exactly

`T = D/(rho_a rho_b r)`,

`alpha=L_a/(rho_a r)`,

`beta=L_b/(rho_b r)`.

Therefore

`K_{b->a}=D L_a/(rho_a^2 rho_b r^2)`,

`K_{a->b}=-D L_b/(rho_a rho_b^2 r^2)`,

and

`P_ab = -D^2 L_a L_b/(rho_a^3 rho_b^3 r^4)`.

The sign condition for a positive two-cycle is simply

`D L_a>0`,

`D L_b<0`.

This decomposition is useful because `D`, `L_a`, and `L_b` respond differently to the same velocity gradient.

## ATTACK — exact balances distinguish transverse renewal from longitudinal renewal

Let

`A_a=grad u(X_a)`, `A_b=grad u(X_b)`,

and let the exact chord-average velocity gradient be

`Abar = integral_0^1 grad u(X_a+sR) ds`,

so that

`Rdot=Abar R`.

For incompressible flow, `tr Abar=0`.

The pair-cell balance is

`Ddot = J_bridge + J_nu`,

where the common affine action cancels exactly and only endpoint/chord mismatch plus viscosity can renew `D`.

The longitudinal overlaps obey a different law.  With `Sbar=(Abar+Abar^T)/2`,

`Ldot_a`

`= 2 omega_a.Sbar R`

`  + [(A_a-Abar)omega_a].R`

`  + nu (Delta omega_a).R`,

and similarly

`Ldot_b`

`= 2 omega_b.Sbar R`

`  + [(A_b-Abar)omega_b].R`

`  + nu (Delta omega_b).R`.

So common affine strain cannot create the oriented three-dimensional pair cell `D`, but it **can** create, destroy, or reverse longitudinal access.  The two gates have different physical renewal channels.

Whenever all factors stay nonzero, define

`lambda_D=Ddot/D`,

`lambda_La=Ldot_a/L_a`,

`lambda_Lb=Ldot_b/L_b`,

`sigma_a=d log rho_a/dt`,

`sigma_b=d log rho_b/dt`,

`sigma_R=d log r/dt`.

Then

`d log G/dt`

`= lambda_La + lambda_Lb - sigma_a - sigma_b - 2 sigma_R`,

and the exact full-cycle product balance is

`d log P/dt`

`= 2 lambda_D + lambda_La + lambda_Lb`

`  - 3 sigma_a - 3 sigma_b - 4 sigma_R`.

Equivalently, if

`Delta_T = -d log|T|/dt`,

`Delta_G = -d log G/dt`,

`Delta_P = -d log P/dt`,

then

`Delta_P = 2 Delta_T + Delta_G`.

A cycle may therefore refill its noncoplanar cell while simultaneously losing longitudinal access faster.

## AUTOPSY

### `T` alone as the cycle lifetime variable

Killed.  `T` detects coplanarization but is blind to loss of the opposite-longitudinal-projection gate.

### pair-cell renewal as full cycle renewal

Killed more strongly.  `D` can be renewed while either `T` or `G` is consumed.

### one common renewal channel for every gate

Killed.  Incompressible common affine deformation cancels from `Ddot`, but common strain acts directly on `L_a` and `L_b`.

## PROMOTE

1. `P_ab=K_{b->a}K_{a->b}=-T^2 alpha beta` is the intrinsic symmetric strength of a positive mutual-stretching pair.
2. Full cycle productivity factorizes into transverse noncoplanarity `T^2` and longitudinal access `G=-alpha beta`.
3. The unnormalized observables `(D,L_a,L_b)` are the natural ancestry variables behind the exact angular kernel.
4. `D` and `(L_a,L_b)` have physically different renewal laws.
5. The exact full-cycle deficit is `Delta_P=2 Delta_T+Delta_G`.

## Next attack

The immediate question is whether repeated `Delta_T<0` events can preserve the full product `P`, or whether renewal of one gate deposits geometric memory that drains the other gate.  An exact time-dependent quadratic Navier-Stokes family gives a clean attack on this question.
