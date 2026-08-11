# Continuous screened-feedback maintenance ledger

Static hiding and continuous hiding are physically different.

For one source-scaled toroidal angular channel write

`partial_tau a_l = g_l + D_l a_l`,

where `D_l` is the exact radial viscous operator. In true Navier–Stokes, `g_l` is not a free control: it is the projection of the Euler transport/stretching term and all angular-channel couplings into this radial channel.

Define

`M_q = C_l[D_l^q a_l]`,

`J_q = C_l[D_l^q g_l]`.

Linearity gives the exact hierarchy

`dot M_q = J_q + M_(q+1)`.

The capacitary endpoint law gives simultaneously

`M_(q+1) = (l+1)[D_l^q a_l(0)-D_l^q a_l(1)]`.

Therefore maintaining

`M_0=...=M_(n-1)=0`

through a positive time interval requires

`J_q=-M_(q+1)`,  `q=0,...,n-1`,

at every time. Initial moment cancellation is not enough.

For a maximally delayed degree-`2n` polynomial, initially

`M_0=...=M_(n-1)=0`,  `M_n != 0`,

so only the top maintenance moment is nonzero:

`J_0=...=J_(n-2)=0`,  `J_(n-1)=-M_n`.

The artificial calibration `g_l=-D_l a_l` freezes the entire profile and realizes the ledger exactly. It is only a microscope. Navier–Stokes has no arbitrary forcing knob of this form.

Modules 155-156 show why this distinction is now central: once a remote high-frequency moment collar exists, both its enstrophy and its viscous palinstrophy can be arbitrarily small. The missing cost cannot be assigned to static storage or static diffusion of the collar. It must be sought in the **nonlinear process that creates and continually replaces the required source moments**.

The next step is therefore to compute `J_q` for the actual Euler operator in the D2/eigenframe Hodge language and reconnect that source to material vorticity/circulation ancestry.
