# Finite-energy turnover collar of the linear null catalyst

The polynomial null catalyst is a local germ, not a finite-energy global field. The physically correct way to localize it is to taper its **velocity**:

`u_chi=chi(r/L) u0`,

then define

`omega_chi=curl u_chi`.

This keeps `div omega_chi=0` automatically and avoids imposing an artificial radial cutoff on vorticity.

For `B=diag(2,-1,-1)`, the uncut field is

`u0=(0,-A x z,A x y)`.

Exact angular integration gives

`E=(4pi A^2/15) int r^6 chi^2 dr`.

For the vorticity,

`omega_chi=chi omega0+chi' n cross u0`.

The raw enstrophy contains three terms. After radial integration by parts, the bulk and cross terms cancel exactly:

`Z=(8pi A^2/15) int r^6 (chi')^2 dr`.

This identifies the physical location of viscosity. The pure linear germ has `Delta omega=0`, but any finite-energy realization must turn its velocity over, and the entire enstrophy of this localized field is carried by that turnover derivative.

The radial vorticity-flux density also has a simple structure:

`omega_chi.n=chi omega0.n`.

Hence the positive cap circulation grows like

`Gamma_cap(r) proportional to A r^3 chi(r/L)`

and returns to zero at the outer boundary. The closure/turnover collar is simultaneously the place where the through-going ancestry flux bends back and where viscosity reappears.

For a fixed taper shape,

`E/(nu Z)=c_chi L^2/nu`.

So a halo with `L>>epsilon` can survive many core clocks `epsilon^2/nu`. Localization restores a viscous clock but does not alone defeat broad reusable ancestry.

There is a sharp weighted derivative floor. If `chi(R)=1` and `chi(infinity)=0`,

`int_R^infinity r^6 (chi')^2 dr >= 5R^5`,

by Cauchy with weight `r^-6`. The relaxed equality profile is `(R/r)^5` outside `R`.

This floor is a closure-collar geometry statement, not yet a nonreusable cost. The remaining question is spacetime: as a singular cascade moves `epsilon,L` inward, can successive turnover collars have summable viscous dissipation and still carry the required high-circulation ancestry?
