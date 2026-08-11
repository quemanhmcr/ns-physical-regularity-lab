# Physical-space tiny-core realization of the remote collar

The remote high-frequency collar should not be trusted merely because a spectral normalization allows it.  It has a direct compact realization in physical space.

Fix source radius `L=1` and a toroidal angular degree `l`.  Choose

`a_epsilon(r)=A_epsilon epsilon^-2 (1-r^2/epsilon^2)^2`,  `0<=r<=epsilon`,

and zero outside.  Let

`omega_epsilon=a_epsilon(r) x cross grad H_l`.

The radial factor and its first derivative vanish at `r=epsilon`, so the zero extension is `C^1`; the vorticity is divergence-free because every radial multiple of the toroidal angular field is divergence-free.

Choose `A_epsilon` from the exact screened Hodge companion functional so

`C_l[a_epsilon]=1`.

As `epsilon->0`, the screen kernel tends to one on the core and `A_epsilon` approaches a finite nonzero constant.

Because `x cross grad H_l` is homogeneous of degree `l`, scaling `x=epsilon y` gives exactly

`int |omega_epsilon|^2 dx = A_epsilon^2 epsilon^(2l-1) Z_l^base`,

and

`int |grad omega_epsilon|^2 dx = A_epsilon^2 epsilon^(2l-3) P_l^base`.

For every `l>=2`, both powers are positive. Hence

`Z_epsilon -> 0`,

`P_epsilon -> 0`,

while

`C_l[a_epsilon]=1`.

This does **not** mean finite enstrophy produces finite kinetic energy for free. `C_l grad H_l` is an internal harmonic companion of the tangent div-curl representation; the rest of the tangent velocity cancels it outside the tiny core.  Module 147's earlier independent-energy interpretation was explicitly demoted for this reason.

The result shows in physical space what modules 155-156 showed spectrally: a local screened/Taylor feedback coefficient can be carried by an ever smaller vorticity core with vanishing static enstrophy and palinstrophy.

The unresolved physics is therefore dynamic: how does real Euler transport/stretching repeatedly create such shrinking cores before viscosity erases them?
