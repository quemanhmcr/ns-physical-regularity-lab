# Conditional relay ancestry throughput

Module 164 says that in the exact throughgoing winding calibration

`tau_nu/tau_trans = (4/5)(Gamma/nu)/sqrt(1+DeltaR^2/(2 pi^2 N^2 r^2))`.

Suppose only for the purpose of an adversarial relay test that a material lineage is useful for at most one characteristic transverse viscous clock `tau_nu` before its local flux core must be replaced or renewed.

Then the number of lineage-lifetimes needed to cover one transaction clock is

`N_life=tau_trans/tau_nu`.

The aggregate circulation presented by those lineage-lifetimes is

`A_relay=Gamma N_life`.

The exact shell identity gives

`A_relay = (5/4) nu sqrt(1+DeltaR^2/(2 pi^2 N^2 r^2))`.

Hence in the high-winding limit

`A_relay -> (5/4)nu`.

This is the dual cancellation to the winding gain:

- `Gamma/nu ->0`;
- `N -> infinity` to keep `Gamma_Q/nu` fixed;
- each lineage becomes short-lived relative to the transaction clock;
- the number of lineage-lifetimes grows like `nu/Gamma`;
- circulation per lineage shrinks like `Gamma`;
- their product remains order `nu`.

## Source-relative meaning

A shrinking Hodge source is not material. Its exact ancestry flux law is

`d/dt int_Sigma Omega = int_boundary [i_(V-u) Omega - j]`.

Therefore any actual relay of ancestry into the active source has only two physical supply channels:

1. relative material crossing of the moving source boundary;
2. Kelvin viscous current.

`A_relay` should be read as a conditional throughput those channels would have to serve, not as a new conserved stock and not as a dissipation norm.

## Crucial caveat

The whole interpretation depends on the statement that one characteristic transverse diffusion clock destroys an order-one fraction of the circulation that is useful to the local transaction.

That is not a theorem for arbitrary 3D tubes.

The next attack therefore uses an exact Navier-Stokes solution, Lamb-Oseen, where circular loops are material and their circulation change is known exactly. If an Oseen core can retain essentially all useful local circulation for many transverse clocks, the relay-throughput interpretation must be demoted.
