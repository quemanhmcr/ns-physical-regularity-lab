# Actual affine Euler source gate for a screened tiny core

The continuous ledger from module 157 says a screened coefficient can remain constant only if the non-viscous source cancels its viscous drift at every time.  The next question is whether the reusable stationary common strain can supply that source.

Use the actual affine Euler vorticity operator

`E_S omega = S omega - (Sx).grad omega`

with the validated stationary productive strain `S_*`.

Align harmonic coordinates with the three strain eigenaxes and choose D2-invariant toroidal harmonics in each eigenplane.  For

`omega=a(r) T_l`,

the projection of `E_S omega` back into the same angular `T_l` channel has the exact form

`Pi_l E_S omega = sigma [alpha_l a(r) + beta_l r a'(r)] T_l`.

The coefficients `alpha_l,beta_l` depend only on angular geometry.  They do not depend on the radial localization scale.

Now take the compact core from module 159,

`a_epsilon=A_epsilon epsilon^-2(1-r^2/epsilon^2)^2`,

normalized so `C_l[a_epsilon]=1`.

Its same-channel affine source moment is

`J_0^aff = sigma C_l[alpha_l a_epsilon + beta_l r a_epsilon']`.

After `r=epsilon s`, this is `sigma` times a finite screened integral; it has no `epsilon^-2` amplification.

The viscous drift is instead exact:

`M_1=C_l[D_l a_epsilon]=(l+1)[a_epsilon(0)-a_epsilon(1)]`

`=(l+1)A_epsilon epsilon^-2`.

The maintenance condition

`dot C_l=J_0^aff+nu M_1=0`

therefore requires, for every nonblind orientation,

`sigma ~ nu epsilon^-2`.

So the `epsilon^-2` gate is not merely dimensional and not merely a Fourier characteristic. It is the exact screened source ledger for the actual affine Euler operator.

This leaves two escapes, both genuinely nonlinear/local:

1. use strain whose magnitude itself grows like `nu/epsilon^2`, forcing its physical Hodge source scale inward as module 158 quantified;
2. use non-affine same-scale velocity/vorticity structure and route the maintenance source through other angular channels.

The second escape is where material ancestry becomes unavoidable: same-scale non-affine source cannot be represented as one reusable large-scale strain reservoir.
