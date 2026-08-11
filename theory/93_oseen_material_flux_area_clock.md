# Lamb-Oseen material flux-area clock

The relay-throughput construction is conditional on a physical premise: a lineage whose transverse flux-area diffusion clock is shorter than the transaction clock should lose an order-one fraction of **useful local material circulation** on that clock.

This premise can be attacked in an exact Navier-Stokes solution.

For Lamb-Oseen,

`Gamma_R(t)=Gamma_inf [1-exp(-q)]`,

`q=R^2/(4 nu t)`.

A circle of fixed `R` is material because the velocity is purely azimuthal.

At time `t0`, the vorticity at the loop is

`omega(R,t0)=Gamma_inf exp(-q0)/(4 pi nu t0)`.

Define the local circulation flux-area scale

`A_Gamma=Gamma_R(t0)/omega(R,t0)`.

Exactly,

`A_Gamma=4 pi nu t0 [exp(q0)-1]`.

Its viscous clock is

`tau_A=A_Gamma/nu`.

After one such clock,

`t1/t0=1+4 pi [exp(q0)-1]`,

`q1=q0/(t1/t0)`,

and the same material loop retains the fraction

`[1-exp(-q1)]/[1-exp(-q0)]`.

The material-loop half-life is also exact. If `Gamma_R(t_half)=Gamma_R(t0)/2`, then

`q_half=-log[(1+exp(-q0))/2]`,

`t_half/t0=q0/q_half`.

Therefore the number of local flux-area clocks to lose half the circulation is

`[(q0/q_half)-1] / [4 pi(exp(q0)-1)]`.

## Physical reading

For a concentrated core (`q0` order one), Oseen gives order-one material-circulation loss on an order-one-or-shorter flux-area clock. This is exactly the ancestry mutation measured by the Kelvin viscous term; Euler swirl does not change the loop circulation.

There is also an honest escape: `q0->0` describes a diffuse tail. Its local loop circulation changes only weakly on one `A_Gamma/nu` clock.

So Oseen supports, but does not universalize, the relay premise. A global Navier-Stokes proof would still need to show that a lineage producing a localized maintenance transaction cannot hide most of its circulation ancestry in a diffuse transverse tail while only a tiny core enters the source.

That diffuse-tail escape is now the next geometry to attack.
