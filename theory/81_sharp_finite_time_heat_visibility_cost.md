# Sharp finite-time cost after hidden initial feedback moments

A fixed revealed derivative is still a local-in-time target.  A more physical question is:

> after hiding the first `n` screened heat moments, what is the smallest spectral burden that can produce a prescribed feedback at an actual positive time `tau_*`?

Stay first on the fixed viscous band `1<=s<=2`.  Let `rho(s)` be the screened feedback spectral density and impose

`integral s^q rho(s) ds = 0`,  `q=0,...,n-1`.

Require

`F(tau_*) = integral rho(s) exp(-s tau_*) ds = 1`.

In `L^2(1,2)`, let `Pi_{n-1}` denote orthogonal projection onto polynomials of degree below `n`.  Since every admissible `rho` is orthogonal to that polynomial space,

`F(tau_*) = <rho, (I-Pi_{n-1}) exp(-s tau_*)>`.

Cauchy-Schwarz is sharp here because the equality direction is itself admissible. Therefore

`||rho||_min = 1 / ||(I-Pi_{n-1}) exp(-s tau_*)||_2`.

This is not a generic norm estimate imposed on NS.  It is the exact Riesz representer of the actual heat/Hodge observation functional after the physically required moment cancellations have been imposed.

Write `s=1+x`, `0<=x<=1`, and use shifted Legendre polynomials `L_j(x)`. Their heat coefficients are

`I_j(T)=integral_0^1 exp(-T x)L_j(x)dx`

and can be written

`I_j(T)=(-1)^j exp(-T/2) i_j(T/2)`,

where `i_j` is the modified spherical Bessel function. Hence

`||(I-Pi_{n-1}) exp(-T x)||_2^2 = sum_{j=n}^infinity (2j+1) I_j(T)^2`.

All terms are nonnegative.  This positive-tail representation is the correct numerical observer for deep hiding; subtracting the first `n` projection energies from the full heat-kernel norm can produce false zero at modest Arb precision.

Restoring the shift from `[0,1]` to `[1,2]` multiplies the squared tail norm by `exp(-2T)`.

Combining with module 153 gives, for each angular channel,

`Z_l >= K_l / ||(I-Pi_{n-1}) exp(-s tau_*)||_2^2`

when unit screened feedback is required at `tau_*`.

For small `T`, the first surviving Legendre component gives

`||(I-Pi_{n-1}) exp(-T x)||_2`

`~ T^n / [n! binomial(2n,n) sqrt(2n+1)]`.

Thus finite-time hiding at a time shorter than the band decay scale reproduces an even steeper burden than simply fixing the first revealed derivative.

The next step must release the fixed band.  A true NS mechanism may move the spectrum to higher decay rates or broaden it; both change the heat clock and must be allowed before any contradiction is claimed.
