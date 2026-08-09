# Smooth impulse packing and the viscous clock

The filament packing ladder could be an artifact of choosing thin vortex tubes.  Replace the filament by a smooth localized vector potential

`A = A_z e_z`,

and define physical vorticity by

`omega = curl A`.

This vorticity is exactly divergence-free.  For localized `A`, its hydrodynamic impulse is

`I = (1/2) integral x cross omega dx = integral A dx`.

A Gaussian packet therefore provides a smooth impulse carrier.  Signed packets can be arranged so that the net impulse, then the first spatial impulse moment, cancel exactly, without introducing filament singularities.

Under physical dilation by `L`, a smooth packet impulse is proportional to amplitude times `L^3`.  If the first surviving impulse moment has spatial order `n`, then keeping that moment fixed requires packet impulse amplitude proportional to `L^-n`.  The velocity energy kernel contributes `L^-3`, while one additional vorticity derivative contributes `L^-5`.  Thus the predicted smooth-field ladder is

`E ~ M_n^2 L^-(2n+3)`,

`Omega = integral |omega|^2 ~ M_n^2 L^-(2n+5)`.

Consequently

`E/(nu Omega) ~ L^2/nu`,

independent of cancellation depth.  The surviving moment may change, but viscosity keeps the same physical diffusion clock.

The experiment also lets a single smooth impulse carrier change aspect ratio.  With

`A_z ~ exp[-r_perp^2/(2s^2)-z^2/(2h^2)]`

and fixed RMS size `2s^2+h^2=3L^2`, shape optimization is allowed to flatten or elongate the carrier.  The question is whether this can make the fixed-impulse kinetic-energy bill disappear.  The parallel diagnostic is enstrophy, because a geometry that cheats energy by creating sharp directional/spatial structure may simply move the bill into viscous dissipation.

This remains a canonical smooth family, not a universal variational theorem.  Its purpose is to remove the filament-core artifact and identify the next escape route honestly.
