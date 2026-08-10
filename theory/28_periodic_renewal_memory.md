# Periodic bridge renewal: angular surplus can recur forever while longitudinal access is spent

## THINK — can the same bridge-current mechanism refill `T` repeatedly?

A pointwise law `Delta_T>=0` is already false.  The stronger possibility is that a finite bridge-current geometry might repeatedly produce `Delta_T<0`, continually refilling the normalized three-dimensional cell without creating smaller structure.

To attack this without numerical time stepping, use the exact quadratic Navier-Stokes family

`u=(a(t)x + eps y z, -a(t)y, 0)`,

with constant `eps>0` and arbitrary smooth scalar `a(t)`.

Its vorticity is

`omega=(0,eps y,-eps z)`,

`Delta u=0`,

and direct substitution gives

`u_t+(u.grad)u=((adot+a^2)x,(a^2-adot)y,0)=-grad p`.

Hence every smooth `a(t)` gives an exact Euler/Navier-Stokes solution on `R^3`.  It has infinite global energy and is used only as a local mechanism microscope.

Choose a positive periodic material stretch factor directly:

`lambda(t)=1+eta sin^2 t`, `eta>0`,

and set

`a(t)=-lambdadot/lambda`.

Then `lambda(t)>=1` and returns exactly to one at every `t=k pi`.

## Exact material bridge

Follow

`X_a(0)=(0,1,0)`,

`X_b(0)=(-L,-1,1)`, `L>0`.

Define

`I(t)=integral_0^t lambda(s)^2 ds`,

`Q(t)=L+eps I(t)`.

Because

`(lambda x)dot=eps y_0 z lambda^2`,

the exact pair separation is

`R=(-Q/lambda,-2 lambda,1)`.

The endpoint vorticities are

`omega_a=(0,eps lambda,0)`,

`omega_b=(0,-eps lambda,-eps)`.

The pair-cell is

`D=-eps^2 Q`,

so

`Qdot=eps lambda^2>0`.

Thus the same bridge-current mechanism continuously renews `|D|`, and `Q` is an exact deposited bridge-memory coordinate in this family.

For the chosen periodic `lambda`,

`I(t)` is elementary:

`I(t)`

`= t(1+eta+3 eta^2/8)`

`  -(eta/2+eta^2/4) sin(2t)`

`  +(eta^2/32) sin(4t)`.

At every stroboscopic return `t_k=k pi`,

`lambda=1`, `lambdadot=0`,

`Q_k=L+eps k pi C_eta`,

where

`C_eta=1+eta+3 eta^2/8>0`.

So the local strain state returns, but the material bridge remembers every previous renewal period through the secular variable `Q_k`.

## PREDICT — `T` can be refilled forever

At the stroboscopic times,

`r_k=sqrt(Q_k^2+5)`,

`T_k=-Q_k/[sqrt(2) sqrt(Q_k^2+5)]`,

`alpha_k=-2/sqrt(Q_k^2+5)`,

`beta_k=1/[sqrt(2) sqrt(Q_k^2+5)]`.

Hence both edges stay positive:

`K_{b->a}=sqrt(2) Q_k/(Q_k^2+5)>0`,

`K_{a->b}=Q_k/[2(Q_k^2+5)]>0`.

Moreover, because `Qdot=eps` and all endpoint magnitude rates vanish at each strobe,

`d log|T|/dt`

`= eps/Q_k - eps Q_k/(Q_k^2+5)`

`= 5 eps/[Q_k(Q_k^2+5)] > 0`.

Thus **there are infinitely many return times with genuine angular renewal surplus**.  A theorem saying the bridge can only refill `T` finitely many times is false even in a smooth exact NS mechanism.

But the total angular reserve gained is finite:

`|T_k| -> 1/sqrt(2)`,

and

`log(|T_infinity|/|T_0|)`

`= (1/2) log(1+5/L^2)`.

Repeated surplus events become weaker and only fill a bounded angular reservoir.

## ATTACK — the longitudinal gate records the renewal history

At the same stroboscopic times,

`G_k=-alpha_k beta_k`

`= sqrt(2)/(Q_k^2+5)`.

Therefore

`G_k -> 0` like `Q_k^{-2}`.

The full cycle product is

`P_k=T_k^2 G_k`

`= Q_k^2/[sqrt(2)(Q_k^2+5)^2]`,

so

`P_k ~ 1/(sqrt(2) Q_k^2)`.

Both directed edges therefore decay as `1/Q_k`, even though `T` is being repeatedly refilled and approaches a nonzero limit.

The stroboscopic product rate is

`d log P/dt`

`= 2 eps(5-Q_k^2)/[Q_k(Q_k^2+5)]`.

After `Q_k>sqrt(5)`, every stroboscopic return is productively losing even while `T` is gaining.

This identifies the physical memory channel:

`Q_k-L=eps k pi C_eta`

is simultaneously

- accumulated pair-cell renewal, because `|D|=eps^2 Q`;
- accumulated bridge displacement, because `R_x=-Q` at each return;
- the variable that drives `alpha,beta -> 0`.

Renewal of the transverse gate is therefore not reusable for free in this family.  It leaves secular material geometry that drains longitudinal access.

## AUTOPSY

### finite number of `T`-surplus bursts

Killed.  Exact smooth NS geometry supports infinitely many stroboscopic times with `d log|T|/dt>0`.

### repeated `T` renewal as sufficient full-cycle survival

Killed.  `T` approaches `1/sqrt(2)` while both directed transactions vanish.

### every renewal burst must create a smaller spatial scale

Killed as a local statement.  Here the polynomial bridge-current structure stays at the same differential order; the payment is instead secular bridge displacement/occupancy.

### `Q` as a universal monotone NS quantity

Do not promote.  `Q` is special to this exact infinite-energy family.

## PROMOTE

1. Full cycle survival must follow both `T` and `G=-alpha beta`, equivalently `P=T^2G`.
2. Repeated angular renewal surplus can occur indefinitely without preserving mutual-stretching strength.
3. In the exact periodic family, renewal deposits a material memory `Q` that grows linearly per strain-return cycle and drains longitudinal access.
4. The cumulative angular surplus is finite even though surplus events recur forever.
5. A future nonreusable-cost theorem must allow at least two physical payment channels: creation of finer bridge structure **or** secular deformation/occupancy/closure memory.

## Next frontier — can finite-energy closure erase renewal memory?

The polynomial family escapes by letting the material bridge acquire unbounded secular displacement.  A finite-energy, closed, localized ancestry network cannot obviously use that escape indefinitely.

The next question is therefore sharper:

**when a productive closed ancestry network repeatedly resets its local strain geometry, can it also reset the bridge-memory/longitudinal-access ledger without either transporting that memory into closure, renewing material ancestry, or creating new smaller bridge structure?**

That is a closure-level reusability problem, not a norm-growth problem.
