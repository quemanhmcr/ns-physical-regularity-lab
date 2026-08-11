# Harmonic folding shell: circulation Reynolds is the angular persistence clock

Module 185 routes high-degree harmonic folding from fixed relative clearance into a near-contact source branch.  The tempting stronger claim

`fine folding => one order-nu independent circulation packet per fold`

is false in the simplest exact source geometry.  The missing variable is lifetime.

## Exact matched harmonic source

Let `n=m+1` be odd, `n>=3`.  With `Y_n` the degree-`n` spherical harmonic whose equatorial trace is `cos(n theta)`, set

`phi_in=-(a/n) r^n Y_n`

inside `r<L`, and

`phi_out=a L^(2n+1)/(n+1) r^(-n-1)Y_n`

outside.  Their radial velocities agree at `r=L`, so the piecewise field is distributionally incompressible and irrotational away from the sphere.  The source is a tangential vortex sheet.

On the equator,

`[u_theta]=-a(2n+1)/(n+1)L^(n-1) sin(n theta)`.

A Stokes loop straddling one positive jump lobe therefore has actual circulation

`Gamma_lobe=2|a|(2n+1)L^n/[n(n+1)]`.

The same harmonic mode gives boundary straight-line curvature source

`Q_L=|partial_z^2 u_x|=|a|(n-1)(n-2)L^(n-3)`.

Hence

`Gamma_lobe/[Q_L L^3]=2(2n+1)/[n(n+1)(n-1)(n-2)] ~ 4/n^3`.

So static circulation per lobe can vanish as angular degree grows.  That KILLs a fixed-packet ancestry obstruction for near-contact harmonic folding.

## Persistence restores the Reynolds gate

The tangential sheet is a degree-`n` toroidal angular pattern.  Its intrinsic tangential viscous clock is

`tau_ang=L^2/[nu n(n+1)]`.

One angular wavelength is `~L/n`; curvature of order `n/L` is needed to bend a straight lineage on that scale.  Define the one-wavelength folding clock

`tau_fold=(n/L)/Q_L`.

Then exactly

`Gamma_lobe/nu`

`=[2n(2n+1)/((n-1)(n-2))] (tau_ang/tau_fold)`.

The bracket tends to `4`.

Thus the same `n^-3` leverage that makes the instantaneous circulation cheap makes the angular source short-lived.  In this calibration,

`Gamma_lobe/nu -> 0 => tau_ang/tau_fold ->0`:

the pattern diffuses angularly before it can build one wavelength of material curvature.

This is the folding analogue of the winding/persistence cancellation and the tiny-core circulation-Reynolds gate.

## Scope and next attack

The source is an ideal vortex-sheet limit.  Smooth divergence-free thin-shell interpolants can approximate the matched inner/exterior fields; radial smoothing introduces another viscous clock and is not assumed to help persistence.  This module does not prove that fixed-time cap multiplicity `N` requires degree `n~N`, nor does it exclude mixtures or repeated low-degree folding.

The surviving escape is regeneration/reuse:

**Can one near-contact high-Re material lineage continually regenerate a rapidly diffusing angular folding pattern while servicing infinitely many shrinking maintenance events, without Kelvin/material ancestry throughput?**
