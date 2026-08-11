# Signed viscous spectral delay: visible modes can still interfere

Module 148 shows that every real radial viscous frequency has a strictly nonzero screened Hodge gain.  That does **not** imply that a signed superposition cannot hide transiently.

Take `N` distinct viscous decay rates

`lambda_j = j+1`,  `j=0,...,N-1`,

and prescribe the screened contribution weights

`A_j = (-1)^j binomial(N-1,j)`.

Because each individual transfer gain `m_l(sqrt(lambda_j))` is positive, a generalized spectral vorticity coefficient can be chosen as

`c_j = -A_j/m_l(sqrt(lambda_j))`.

The total screened heat signal is then

`F_N(tau) = sum_j A_j exp(-lambda_j tau)`

and the binomial identity gives the exact closed form

`F_N(tau) = exp(-tau) [1-exp(-tau)]^(N-1)`.

Therefore

`F_N^(q)(0)=0` for `q=0,...,N-2`,

while

`F_N^(N-1)(0)=(-1)^(N-1)(N-1)!`.

So arbitrarily deep finite derivative hiding is possible even after viscosity has been diagonalized.  The polynomial delayed profiles were not merely artifacts of using a Taylor basis.

But the same formula gives

`F_N(tau)>0` for every `tau>0`.

The packet is never silent on a positive time interval.  Its maximum occurs when `exp(-tau)=1/N`, hence

`tau_peak = log N`,

with

`F_N(tau_peak) = N^-1 (1-N^-1)^(N-1)`.

This sharpens the spectral conclusion:

- there is no **single-frequency blind spot**;
- signed frequency interference can suppress arbitrarily many initial time derivatives;
- finite packets cannot maintain exact silence for any positive time interval.

The discrete Helmholtz modes used here are generalized radial eigenfunctions, not finite-energy states on the whole effective radial space.  The next physical question is whether finite-energy Hankel wave packets can approximate the same long transient cancellation while respecting actual enstrophy and spatial closure.  That question must be attacked directly rather than inferred from the discrete packet.
