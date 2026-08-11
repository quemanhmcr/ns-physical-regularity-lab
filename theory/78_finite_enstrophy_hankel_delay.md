# Finite-enstrophy Hankel states with arbitrarily deep transient hiding

Discrete radial Helmholtz modes are generalized eigenfunctions, so a finite-frequency packet by itself does not settle the finite-enstrophy question.  The spectral construction can be made genuinely `L^2`.

Use the viscous spectral variable

`s = k^2`

and choose, for any integer `N>=2`,

`eta_N(s) = (s-1)^N (2-s)^N` on `1<=s<=2`, zero outside.

This function is nonnegative and vanishes to order `N` at both endpoints.  Define

`rho_N(s) = d^(N-1) eta_N / ds^(N-1)`.

After `N-1` integrations by parts,

`F_N(tau) = integral_1^2 exp(-s tau) rho_N(s) ds`

becomes

`F_N(tau) = tau^(N-1) integral_1^2 exp(-s tau) eta_N(s) ds`.

Hence

`F_N^(q)(0)=0` for `q=0,...,N-2`,

but

`F_N(tau)>0` for every `tau>0`.

So there is arbitrarily deep finite-order hiding without any positive-time interval of exact silence.

Now connect this feedback density to an actual radial Hankel spectrum.  In effective dimension `d=2l+3`, the screened transfer multiplier from module 148 is

`m_l(sqrt(s)) > 0` for every `s>0`.

On the compact interval `[1,2]` it is continuous and therefore has a strictly positive minimum.  In a standard radial Hankel normalization, the feedback spectral density is proportional to

`m_l(sqrt(s)) ahat(sqrt(s)) s^((d-2)/2)`.

Thus choose

`ahat(sqrt(s)) = const * rho_N(s) / [m_l(sqrt(s)) s^((d-2)/2)]`.

Because `rho_N` is compactly supported and square-integrable, while the denominator is smooth and bounded away from zero on `[1,2]`, `ahat` belongs to the radial spectral `L^2` space.  Hankel Plancherel therefore gives a finite-enstrophy radial profile.

This kills a tempting shortcut:

`finite enstrophy + viscosity` does **not** imply a fixed finite number of feedback derivatives detects every hidden radial state.

The correct rigid statement is weaker but exact:

- every individual radial viscous frequency is Hodge-visible;
- finite-enstrophy signed spectra can cancel arbitrarily many initial feedback derivatives;
- these constructions cannot remain exactly silent on a positive time interval.

The remaining physical issue is quantitative rather than qualitative: after normalizing the actual finite-enstrophy spectrum, how long and how strongly can the screened feedback remain small?  That requires measuring the Hankel `L^2` cost of the cancellation rather than counting derivative order.
