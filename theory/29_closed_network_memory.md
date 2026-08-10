# Closed-network renewal memory: closure cancels edge sums but not material shape

## THINK — bridge memory is an exact material-graph field

For material particles `X_i(t)`, every oriented bridge is

`R_ij = X_j-X_i`.

Therefore bridge geometry is already an exact one-cochain on the material graph.  Around every closed polygon,

`R_12+R_23+...+R_m1=0`.

The same is true of the bridge velocity current:

`Rdot_ij = u(X_j)-u(X_i)`,

so

`sum_loop Rdot_ij=0`.

This is not a dissipative conservation law.  It is the kinematic fact that bridge variables are node-potential differences.  A closure theorem based only on summing signed edge memory can therefore miss arbitrarily large equal-and-opposite deposits.

Use again the exact quadratic Navier-Stokes family

`u=(a(t)x+eps y z,-a(t)y,0)`,

with positive material stretch `lambda` and

`a=-lambdadot/lambda`,

`I(t)=integral_0^t lambda^2 ds`.

A material node with initial data `(x_i0,y_i,z_i)` follows

`y_i(t)=lambda y_i`,

`z_i(t)=z_i`,

`x_i(t)=[x_i0+eps y_i z_i I]/lambda`.

Define the scaled node coordinate

`Phi_i=lambda x_i=x_i0+eps c_i I`,

where

`c_i=y_i z_i`.

Then the scaled edge memory is exactly

`q_ij=lambda R_ij,x = Phi_j-Phi_i`.

Its deposited part is

`mu_ij=eps I(c_j-c_i)=phi_j-phi_i`,

with node memory potential

`phi_i=eps I c_i`.

Hence every closed-loop memory sum vanishes identically even when every individual edge deposit is large.

## PREDICT — a closed productive triangle can still become sterile

Choose three distinct ordered parameters

`theta_1<theta_2<theta_3`

and initial transverse coordinates on the unit hyperbola

`y_i=cosh(theta_i)`,

`z_i=sinh(theta_i)`.

Then

`y_i^2-z_i^2=1`

and

`c_i=y_i z_i=(1/2)sinh(2 theta_i)`

is strictly increasing.

Choose the initial x coordinate to share one positive memory clock:

`x_i0=ell c_i`, `ell>0`.

At every stroboscopic return with `lambda=1`, define

`S=ell+eps I>0`.

Then

`X_i=(S c_i,y_i,z_i)`

and every oriented x-bridge is

`q_ij=S(c_j-c_i)`.

For every pair `i<j`, put

`delta=theta_j-theta_i>0`,

`h=cosh(delta)-1>0`,

`C_i=cosh(2 theta_i)`, `C_j=cosh(2 theta_j)`.

The endpoint vorticities are

`omega_i=eps(0,y_i,-z_i)`.

The exact pair quantities are

`D_ij=eps^2 q_ij sinh(delta) > 0`,

`L_i=omega_i.R_ij=eps h > 0`,

`L_j=omega_j.R_ij=-eps h < 0`.

Thus every unordered pair has the correct opposite longitudinal signs for a positive mutual-stretching two-cycle.

Writing

`r_ij^2=q_ij^2+(y_j-y_i)^2+(z_j-z_i)^2`,

we obtain

`T_ij=q_ij sinh(delta)/[sqrt(C_i C_j) r_ij] > 0`,

`alpha_i=h/[sqrt(C_i) r_ij] > 0`,

`beta_j=-h/[sqrt(C_j) r_ij] < 0`.

Therefore both directed transactions are positive:

`K_{j->i}=T_ij alpha_i>0`,

`K_{i->j}=-T_ij beta_j>0`.

So the closed three-node network is not merely topologically closed: every one of its three unordered edges is a genuine positive two-cycle at each stroboscopic return.

## ATTACK — closure redistributes memory into material occupancy

The three oriented memory edges obey

`q_12+q_23+q_31=0`

exactly, because they are differences of the same node potential `S c_i`.

Nevertheless, as `S -> infinity`, every edge has `|q_ij| -> infinity` because the `c_i` are distinct.

For each pair,

`T_ij -> sinh(delta)/sqrt(C_i C_j)`,

a finite nonzero three-dimensional limit, while

`G_ij=-alpha_i beta_j`

`=h^2/[sqrt(C_i C_j) r_ij^2]`

`~ h^2/[sqrt(C_i C_j) q_ij^2]`.

Hence the full cycle product satisfies

`P_ij=T_ij^2 G_ij`

`~ [sinh(delta)^2 h^2]/[(C_i C_j)^(3/2) q_ij^2]`.

All three pairwise positive cycles therefore become longitudinally sterile like `S^-2`, even though their normalized transverse cells approach nonzero limits.

The closure cancellation has not erased memory.  It has hidden it in the shape of the material packet.

Let

`cbar=(1/3) sum_i c_i`.

At a stroboscopic return the axial central second moment is

`O_x=(1/3) sum_i (x_i-xbar)^2`

`   =S^2 V_c`,

where

`V_c=(1/3) sum_i(c_i-cbar)^2>0`.

For the deposited displacement relative to the initial packet,

`M_x=(S-ell)^2 V_c`.

Equivalently, if

`mu_ij=(S-ell)(c_j-c_i)`,

then for three nodes

`sum_{i<j} mu_ij^2 = 9 M_x`.

This identity is not introduced as an abstract graph norm.  `M_x` is the physical central second moment of the deposited material displacement; the pairwise formula merely states the same occupancy in edge coordinates.

Thus a signed loop observer reports zero memory while the material packet itself records an ever-growing spatial history.

## AUTOPSY

### closure forces renewal memory to cancel physically

Killed.  Closure forces the **signed edge sum** to vanish because bridge memory is an exact graph field.  Individual edge memories and packet occupancy can still grow without bound.

### a nonzero loop memory circulation is the right cost

Killed.  In the exact family the loop circulation is identically zero for purely kinematic reasons, even while every productive edge is being driven toward longitudinal starvation.

### every closed edge can stay productively strong while memory is redistributed

Killed in this exact hyperbolic triangle.  If the common memory clock `S` diverges, all three positive pair cycles have `P_ij~S^-2`.

### material occupancy growth is already an irreversible cost theorem

Do not promote.  The polynomial exact solution has infinite global energy and allows the packet to spread arbitrarily far.  The second moment is a physical memory channel, not yet a finite-energy bound.

## PROMOTE

1. Bridge displacement and bridge velocity current are exact material-graph one-cochains: they are node-potential differences.
2. Closed-loop signed memory sums therefore vanish kinematically and cannot serve as a universal nonreusable-cost ledger.
3. Equal-and-opposite edge memory is stored in material **shape/occupancy**, not erased.
4. A hyperbolic three-node exact NS calibration supports positive mutual-stretching on every pair while a single secular memory clock drives every full cycle product to zero by longitudinal starvation.
5. The centered material second moment is a gauge-invariant physical observer of memory that signed closure sums miss.

## Next question

Closure alone is not enough.  A finite-energy candidate could try to avoid secular occupancy by resetting shape repeatedly rather than letting `S` grow.

The next attack must therefore distinguish two very different notions that an observer might confuse:

- repeated **events** in which productive geometry is refilled;
- net **deposited material memory**.

Can infinitely many apparent renewal events be packed into finite physical time while the deposited memory remains bounded?  If yes, event count and even cumulative positive variation are observer artifacts, and the real payment must lie in the deformation needed to execute the resets.
