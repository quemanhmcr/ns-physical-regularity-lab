# Recurrence is not memory: finite-time reset action after quotienting rigid motion

## THINK — event count can be an observer artifact

The periodic exact family showed infinitely many angular-surplus returns, but there the returns occurred over infinite physical time and deposited bridge memory grew secularly.

An observer might therefore try to count productive refills, or integrate the positive variation of a normalized gate, as an irreversible cost.

That is too fast.  The exact quadratic family permits arbitrary positive material stretch histories:

`u=(a(t)x+eps y z,-a(t)y,0)`,

`a=-d log(lambda)/dt`.

Choose an accelerated phase on `0<=t<1`,

`phi(t)=t/(1-t)`,

and

`lambda(t)=1+eta sin^2(phi(t))`.

Then infinitely many exact returns

`phi(t_k)=k pi`

accumulate at finite time `t=1`, with

`t_k=k pi/(1+k pi) -> 1`.

At every return `lambda=1` and `lambdadot=0`; at every intervening crest `lambda=1+eta` and again `lambdadot=0`.

The sampled local strain state is therefore completely benign at the return and crest instants even though the transitions become faster and faster.

## PREDICT — infinitely many refills can coexist with bounded deposited memory

The bridge memory clock remains

`I(t)=integral_0^t lambda(s)^2 ds`.

Because

`1 <= lambda <= 1+eta`,

for every `t<1`,

`t <= I(t) <= (1+eta)^2 t < (1+eta)^2`.

Thus finite-time accumulation of infinitely many recurrence events does **not** force unbounded deposited bridge memory.

Use the same pair as the periodic microscope:

`X_a(0)=(0,1,0)`,

`X_b(0)=(-L,-1,1)`.

Write

`Q=L+eps I`.

For `eta=1`, `L=2`, `eps=1`, all times before `t=1` obey

`2 <= Q < 6`.

At a return (`lambda=1`),

`|T_return|=Q/[sqrt(2) sqrt(Q^2+5)]`.

At a crest (`lambda=2`),

`|T_crest|=Q/[sqrt(5) sqrt(Q^2+68)]`.

Both expressions increase with `Q`.  Consequently every return-to-crest excursion obeys the uniform bounds

`|T_return| >= 2/[sqrt(2)sqrt(9)]`,

`|T_crest| <= 6/[sqrt(5)sqrt(104)]`.

Their difference is strictly larger than `0.2`.

Therefore every accelerated cycle contains a finite-amplitude decrease and refill of `|T|`, while the entire deposited memory `Q-L` remains below four.

There are infinitely many cycles before `t=1`, so the total positive variation of `|T|` is infinite even though the material bridge-memory deposit is bounded.

At each return the exact instantaneous angular-surplus rate is

`d log|T|/dt=5 eps/[Q(Q^2+5)]`.

Since `Q<6`, every return has a rate larger than

`5/(6*41)>0.02`.

Thus even infinitely many **uniformly positive instantaneous surplus samples** do not by themselves represent an accumulating physical resource.

The recurrence tail itself must be observed intrinsically as `tau_k=1/(1+k pi)`.  Computing it as `1-t_k` after `t_k` has numerically saturated near one recreates observer complexity by subtracting a tiny physical tail from a large parent state.

## ATTACK — what actually becomes nonreusable when resets accelerate?

The stretch rate is not an arbitrary analytical quantity here.  It is exactly the material logarithmic deformation rate:

`a=-d log(lambda)/dt`.

During each complete excursion

`lambda: 1 -> 1+eta -> 1`,

we have the parameterization-independent action

`integral_cycle |a| dt`

`= integral_cycle |d log(lambda)|`

`=2 log(1+eta)`.

It does not matter how rapidly the cycle is executed.

Hence after `N` complete resets,

`A_reset(N)=2N log(1+eta)`.

Packing infinitely many finite-amplitude resets into finite time while keeping the net bridge memory bounded therefore forces

`integral_0^1 |a(t)| dt = infinity`.

This is a very different object from event count.  It is the actual path length travelled by a material stretch degree of freedom in multiplicative deformation space.

The exact accelerated family therefore exposes a clean trilemma:

1. let bridge memory/occupancy grow secularly;
2. let productive geometry lose amplitude;
3. repeatedly reset shape, paying unbounded material deformation action if infinitely many finite-amplitude resets are compressed into finite time.

This is still only a calibration, not a global Navier-Stokes theorem.

## AUTOPSY — rigid motion must be quotiented out

One more observer trap remains.  A packet can execute arbitrarily large rigid rotation without changing its productive pair geometry.

For any common `Q in SO(3)`,

`p -> Qp`, `R -> QR`, `q -> Qq`

preserves all dot products and scalar triple products.  Hence it preserves

`D`, `T`, `alpha`, `beta`, `G`, and `P`.

Infinitely large rotational path length is therefore not a productive reset cost.

Infinitesimally, for a common skew matrix `W`,

`d(p.R)/dt=(Wp).R+p.(WR)=0`,

while for a general common velocity gradient `A=S+W`,

`d(p.R)/dt=2 p.S R`.

The common rigid spin `W` drops out exactly.  Only deformation, not frame rotation, changes longitudinal access.

Thus any future reset-action theorem must live in material **shape space** -- equivalently in the Cauchy-Green/symmetric-strain part after common rigid motion has been removed.

## KILL

- number of renewal-surplus events as an irreversible cost;
- infinite recurrence in finite time as proof of unbounded bridge-memory deposition;
- cumulative positive variation of `T` as a universal nonreusable resource;
- total angular/rotational motion of the packet as a productive reset cost.

The accelerated exact family gives infinitely many finite-amplitude `T` refills and infinite positive variation before finite time while its deposited bridge memory stays bounded.  Common rigid rotation can also accumulate arbitrary motion with exactly zero effect on productive geometry.

## PROMOTE

1. Net deposited bridge memory and recurrence count are physically different currencies.
2. Finite-amplitude shape resets have a parameterization-independent logarithmic stretch action in the exact family:

   `A_cycle=2 log(1+eta)`.
3. Infinite reset count in finite time can avoid secular occupancy only by accumulating infinite shape-deformation action in this calibration.
4. Common rigid rotation is a gauge mode for pair productivity; it must be removed before measuring reset action.
5. The next plausible physical object is therefore not a norm of `grad u`, but the **non-common material deformation/holonomy required to return a productive closed network to reusable shape**.

## Next frontier — deformation holonomy versus ancestry renewal

A genuine finite-energy NS candidate is free to combine strain, rotation, bridge inhomogeneity, viscosity, and ancestry replacement.  The next question is whether a closed productive network can execute infinitely much **shape-reset holonomy** in finite time while avoiding all three observable consequences:

- secular packet occupancy;
- longitudinal/transverse productive starvation;
- new ancestry entering through the viscous spacetime current.

If common rigid motion is quotiented out, any such reset must be produced by symmetric deformation or non-affine/viscous relative motion.  That is the natural place to reconnect the cycle-survival hierarchy to the material-spacetime ancestry ledger.
