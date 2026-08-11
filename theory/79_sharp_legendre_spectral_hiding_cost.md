# Sharp fixed-band cost of hiding viscous feedback derivatives

Work with the feedback spectral density `rho(s)` on the fixed viscous decay-rate band `1<=s<=2`.

To hide the screened heat feedback derivatives through order `n-1`, require

`integral_1^2 s^q rho(s) ds = 0`,  `q=0,...,n-1`.

Fix the first revealed derivative by

`integral_1^2 s^n rho(s) ds = R`.

This is an exact Hilbert-space minimization problem in `L^2(1,2)`.  Let

`L_n(s)=P_n(2(s-1)-1)`

be the shifted Legendre polynomial on `[1,2]`.  It spans the one-dimensional component of degree `n` orthogonal to all lower polynomials.  Its exact identities are

`||L_n||_2^2 = 1/(2n+1)`,

`leading(L_n) = binomial(2n,n)`,

and therefore

`integral_1^2 s^n L_n(s) ds = 1/[(2n+1) binomial(2n,n)]`.

By orthogonal projection / Riesz representation, every admissible density satisfies the sharp bound

`||rho||_2 >= |R| sqrt(2n+1) binomial(2n,n)`.

Equality holds for the normalized shifted Legendre density.

Using the central-binomial asymptotic,

`binomial(2n,n) ~ 4^n/sqrt(pi n)`,

we obtain

`||rho||_2 / |R| ~ sqrt(2/pi) 4^n`.

Thus, **on a fixed positive viscous-frequency band and at fixed first revealed feedback derivative**, each additional hidden initial derivative costs asymptotically a factor four in spectral-amplitude `L^2` burden, or a factor sixteen in its square.

This cost is not imposed by choosing a Sobolev norm before understanding the mechanism.  It arises after the physical Hodge/viscous transfer has been diagonalized and from the exact moment cancellations required to keep that transfer hidden.

The scope is essential.  The result does not yet exclude:

- shifting the spectral support with `n`;
- broadening the support;
- using nonlinear time-dependent replenishment;
- mixing angular channels;
- changing the source radius.

Those are the next escapes to attack before this can enter a regularity contradiction.
