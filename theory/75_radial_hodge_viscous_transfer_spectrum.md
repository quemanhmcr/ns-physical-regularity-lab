# Radial Hodge-viscous transfer spectrum

The radial profile in a fixed toroidal angular channel should not be organized by polynomial degree once viscosity is the mechanism under study.  In source-scaled radius `s=r/L`,

`D_l = d^2/ds^2 + (2l+2)/s d/ds`

is exactly the radial Laplacian in effective dimension

`d = 2l+3`.

Its normalized regular generalized eigenmode is

`Phi_d(z s) = 0F1(;d/2;-z^2 s^2/4)`,

with

`D_l Phi_d(z s) = -z^2 Phi_d(z s)`,  `Phi_d(0)=1`.

The screened Hodge feedback functional already satisfies

`C_l[D_l a] = (l+1)(a(0)-a(1))`.

Applying it to one spectral mode gives the exact transfer law

`C_l[Phi_d(z s)] = -(l+1)/z^2 [1-Phi_d(z)]`.

Define the positive gain

`m_l(z) = (l+1)[1-Phi_d(z)]/z^2`.

The value `Phi_d(z)` has the geometric representation

`Phi_d(z) = average_{theta in S^(d-1)} cos(z theta_1)`.

Therefore for every real `z != 0`,

`1-Phi_d(z) = average [1-cos(z theta_1)] > 0`.

There is no real radial viscous frequency hidden from the screened Hodge feedback.  At zero frequency,

`m_l(0) = (l+1)/(2(2l+3))`

by the regular limit.

A generalized frequency `z=kL` decays under viscosity as

`exp(-z^2 tau)`,  `tau=nu t/L^2`,

so its dimensionless lifetime is `tau_nu=1/z^2`.  The transfer law becomes

`m_l(z) = (l+1)[1-Phi_d(z)] tau_nu`.

Thus high radial frequency can make the instantaneous lower Hodge response small, but only on the same `z^-2` scale on which viscosity destroys that frequency.  For `z -> infinity`, `Phi_d(z) -> 0` and

`m_l(z)/tau_nu -> l+1`.

For a finite packet of distinct frequencies,

`F(tau)=sum_j A_j exp(-lambda_j tau)`,  `lambda_j=z_j^2`,

the first `N` time derivatives form a Vandermonde system.  If all `N` derivatives vanish for `N` distinct frequencies, every `A_j` vanishes.  Since `m_l(z_j)>0`, no nontrivial finite spectral packet can be permanently Hodge-silent.

For sufficiently decaying continuous Hankel spectra, the same statement is expected from uniqueness of the Laplace transform: a feedback signal that vanishes on a time interval is analytic for positive time, hence vanishes for all positive time, and the nonvanishing transfer multiplier then forces the spectral density to vanish.  This last continuous-spectrum extension should be promoted only after its functional hypotheses are stated carefully.

The polynomial delayed-visibility families from modules 145-147 are still valid local calibrations, but they are not viscosity eigenmodes.  Their high-order cancellations should therefore be understood as coordinated superpositions of radial spectral content, not as blind frequencies of the heat operator.
