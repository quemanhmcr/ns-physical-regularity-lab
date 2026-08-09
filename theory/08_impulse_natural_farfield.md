# Natural far-field variable: hydrodynamic impulse

The previous closed-ring donor result should not be generalized by imposing a preferred ring shape.  Let the closure of vorticity itself decide the far-field variables.

For localized incompressible vorticity,

`u(x) = (1/4pi) integral omega(y) x (x-y)/|x-y|^3 dy`.

At large `r=|x|`, the apparent `1/r^2` term is proportional to `integral omega dy`.  For localized divergence-free vorticity this vanishes.  For a closed filament it is simply `Gamma integral dX = 0`.  Thus nature removes the monopole exactly.

The first surviving moment is the hydrodynamic impulse

`I = (1/2) integral y x omega(y) dy`.

For a thin closed filament,

`I = (Gamma/2) integral X x dX = Gamma * (vector area)`.

Using divergence-free closure, the leading far velocity is

`u_far(x) = [3 e (I.e) - I]/(4 pi r^3)`, where `e=x/r`.

Consequently the leading far strain is order `|I|/r^4`; differentiating gives

`S_ik = [3 delta_ik(I.e) + 3(e_i I_k + e_k I_i) - 15(I.e)e_i e_k]/(4 pi r^4)`.

This is a physically generated coarse-graining: the remote field forgets microscopic vortex shape and retains the first non-cancelled vortex moment.

There is a second structural fact.  For sufficiently localized unforced Navier-Stokes flow, impulse is globally conserved.  From the vorticity equation,

`dI/dt = integral (u x omega) dx + (nu/2) integral x x Delta omega dx`.

The viscous term is a boundary term and vanishes under localization.  Also

`u x omega = grad(|u|^2/2) - (u.grad)u`,

whose space integral vanishes for decaying incompressible flow.  Hence `dI/dt=0`.

This does **not** yet bound local positive/negative impulse pieces; local structures may rearrange and cancel while total vector impulse stays fixed.  It does identify the natural far-field collateral variable much more cleanly than circulation of an arbitrarily chosen ring.

If a target core at scale `ell` requires strain `s ~ Gamma_c/ell^2` from a remote compact donor at distance `d`, the dipole law suggests the universal far-field requirement

`|I_d| ~ s d^4` up to orientation/efficiency factors.

For `d = Lambda ell`, this becomes `|I_d| ~ Gamma_c Lambda^4 ell^2`.  When the donor itself has size comparable to `d`, `I_d ~ Gamma_d d^2` and the earlier ring collateral `Gamma_d ~ Lambda^2 Gamma_c` is recovered as a special case.

The next obstruction question is therefore not "how much circulation is far away?" but "how much hydrodynamic impulse can a localized donor deploy without requiring compensating impulse/energy elsewhere?"
