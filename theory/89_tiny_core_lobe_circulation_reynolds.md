# Lobe circulation is the self-maintenance Reynolds number

For the explicit compact tiny core use

`H_l=Re(x+iy)^l`,

`omega=a(r) x cross grad H_l`,

with `C_l[a]=1`.

On the equatorial disk `z=0`,

`omega_z=a(R) partial_phi H_l`

and

`partial_phi H_l=-l R^l sin(l phi)`.

One positive angular lobe has width `pi/l`, and its angular integral is exactly `2`. Therefore its vorticity flux is

`Gamma_lobe = 2 |A_epsilon| epsilon^l integral_0^1 (1-s^2)^2 s^(l+1) ds`

`=16 |A_epsilon| epsilon^l /[(l+2)(l+4)(l+6)]`.

Because `curl u=omega`, this is the actual Stokes circulation around the boundary of that lobe sector.

Module 161 gives, for the same core,

`||Euler_self||_2 / ||nu Delta omega||_2 = c_l epsilon^l/nu`.

Hence

`||Euler_self||_2 / ||nu Delta omega||_2`

`= kappa_l(epsilon) Gamma_lobe/nu`,

where `kappa_l(epsilon)` tends to a finite nonzero shape constant as `epsilon->0`.

This identifies the physically relevant self-maintenance parameter:

`R_Gamma = Gamma_lobe/nu`.

A unit-screened tiny core has

`R_Gamma ->0`.

That is why its self-induced Euler dynamics loses to viscosity.  The failure is not mysterious and is not a generic norm phenomenon: the core carries vanishing circulation ancestry at the scale where it asks to maintain an `O(1)` local screened feedback coordinate.

By contrast modules 158-160 require an external maintenance rate

`sigma epsilon^2/nu = O(1)`.

A vortical strain source at that same scale has transaction circulation `Gamma_Q~sigma epsilon^2`, so its transaction Reynolds is order one.  The maintenance actor therefore needs circulation/transaction ancestry of order `nu`, while the tiny core itself carries only `o(nu)` circulation.

This is the next ancestry bottleneck.  The question is no longer whether a small core can exist cheaply; it can.  The question is where the order-`nu` circulation lineage needed to keep producing or replacing such cores comes from, and whether one material lineage can be reused through infinitely many shrinking maintenance events.
