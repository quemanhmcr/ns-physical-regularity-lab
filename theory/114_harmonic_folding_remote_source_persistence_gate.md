# Remote harmonic folding source: near contact minimizes the persistence threshold

Module 186 gives a thin-shell identity between source-lobe circulation Reynolds and the ratio of angular diffusion time to one-wavelength folding time.  Could the source move farther away to gain a longer diffusion clock while keeping the circulation ancestry cheap?

Let the core radius be `epsilon` and the matched degree-`n` source sphere be

`L=lambda epsilon`, `lambda>=1`.

For one homogeneous inner harmonic mode, if

`Q_epsilon=|partial_z^2 u_x|` at the core,

then homogeneity gives

`Q_L/Q_epsilon=lambda^(n-3)`.

Using module 186,

`Gamma_source/nu = c_G H_epsilon lambda^n`,

where

`H_epsilon=Q_epsilon epsilon^3/nu`,

`c_G=2(2n+1)/[n(n+1)(n-1)(n-2)]`.

The source angular viscous clock is

`tau_ang(L)=L^2/[nu n(n+1)]`,

while the core one-wavelength folding clock is

`tau_fold(epsilon)=(n/epsilon)/Q_epsilon`.

Therefore

`tau_ang(L)/tau_fold(epsilon)=H_epsilon lambda^2/[n^2(n+1)]`.

Eliminating `H_epsilon` gives the exact source-radius law

`Gamma_source/nu`

`= C_n lambda^(n-2) [tau_ang(L)/tau_fold(epsilon)]`,

with

`C_n=2n(2n+1)/[(n-1)(n-2)] ->4`.

Hence, if the degree-`n` source pattern is required to survive at least one core one-wavelength folding time,

`Gamma_source/nu >= C_n lambda^(n-2)`.

For every `lambda>=1`, the minimum threshold occurs at `lambda=1`.  Moving the source outward buys only a `lambda^2` diffusion-time gain but pays a `lambda^n` harmonic-transmission cost in source circulation.

So remote placement is not an escape from the module-186 persistence gate for one matched harmonic degree.  Near contact is the optimal source placement in circulation-Reynolds currency.

## Scope

This is still a single-mode theorem with a deliberately explicit one-wavelength persistence criterion.  It does not control nonlinear regeneration, degree mixtures, repeated low-degree folding, moving source geometry, or a fully smooth 3D source collar.  Those are now the correct reuse escapes to attack; static source radius is not.
