# Coupled degree-four/six Hodge feedback search

## THINK — solve only the feedback that physics leaves unslaved

The complete simultaneous degree-four/degree-six maintenance problem reduces to the physical nine-dimensional fixed-point law

`F(y)=Pi_T4 K66^{-1}[-R6(V4(y),y)]-y=0`,

where

`V4(y)=-K44^{-1}(N4+K46 y)`.

The other 79 amplitudes are not free optimization variables: Hodge response invertibility slaves them to `y`.

## Exact directional Jacobian

The map is quadratic because `V4` depends linearly on `y` and the Euler bracket is bilinear.  For a direction `e_a`, let

`dV4_a=-K44^{-1}K46 e_a`.

Then the derivative of the degree-six nonlinear source is assembled directly from the Euler brackets

`[dV4,u3] + [omega2,dU5_4]`

`+ [dV4,U3_4] + [V4,dU3_4]`

`+ [dV4,U3_y] + [V4,dU3_y]`.

Applying `-K66^{-1}` and projecting to the physical `T4` domain gives `D Phi(y)e_a`; subtract `e_a` to obtain `DF(y)e_a`.

No finite-difference derivative is used.

## Observer discipline

The search starts from the `T4` component of the old sequential degree-six servo.  Because the nested linear inversions are strongly conditioned at 160 bits, the search observer recenters each Newton iterate at its Arb midpoint.  This deliberately makes module 132 a **root search/calibration**, not a root-existence certificate.

A candidate is promoted only if a later interval/Krawczyk calculation encloses a true fixed point without midpoint recentering.

## PREDICT / ATTACK

Three outcomes are informative:

1. Newton converges rapidly to a stable candidate: the sequential obstruction is an observer artifact of solving a triangular system in the wrong order, and a coupled local servo may exist.
2. The physical 9x9 Jacobian loses rank: coupled compatibility hits an intrinsic response degeneracy.
3. Newton remains nonconvergent across precision: the feedback may have no nearby branch, but nonexistence still requires a structural or interval argument.

No conclusion about full Navier-Stokes regularity is made at this stage.
